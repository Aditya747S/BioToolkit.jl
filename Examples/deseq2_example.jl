using CSV
using DataFrames
using Statistics
using LinearAlgebra
using SparseArrays

if !isdefined(Main, :BioToolkit)
    include(joinpath(@__DIR__, "..", "src", "BioToolkit.jl"))
end

const BT = BioToolkit

const COUNTS_PATH = joinpath(@__DIR__, "counts_matrix (2).csv")
const COLDATA_PATH = joinpath(@__DIR__, "col_data.csv")
const R_RESULTS_PATH = joinpath(@__DIR__, "r_results.csv")
const R_DISPERSIONS_PATH = joinpath(@__DIR__, "r_dispersions.csv")
const JULIA_RESULTS_PATH = joinpath(@__DIR__, "julia_results.csv")
const JULIA_DISPERSIONS_PATH = joinpath(@__DIR__, "julia_dispersions.csv")

function read_counts_csv(csv_path::AbstractString)
    counts_df = CSV.read(csv_path, DataFrame)
    gene_ids = String.(counts_df[:, 1])
    sample_ids = String.(names(counts_df)[2:end])
    counts_matrix = Matrix{Int}(counts_df[:, 2:end])
    return sparse(counts_matrix), gene_ids, sample_ids
end

function read_coldata_csv(csv_path::AbstractString, expected_sample_ids::Vector{String})
    raw_coldata = CSV.read(csv_path, DataFrame; missingstring=["NA"])
    sample_ids = String.(raw_coldata[:, 1])
    condition = Symbol.(String.(raw_coldata[:, 2]))

    lookup = Dict(sample_id => index for (index, sample_id) in enumerate(sample_ids))
    ordered_indices = [lookup[sample_id] for sample_id in expected_sample_ids]
    ordered_conditions = condition[ordered_indices]

    return DataFrame(condition = ordered_conditions), ordered_conditions
end

function read_result_table(csv_path::AbstractString)
    return CSV.read(csv_path, DataFrame; missingstring=["NA"])
end

function resolve_column_name(results_df::DataFrame, candidates::AbstractVector{<:AbstractString})
    available = Set(String.(names(results_df)))
    for candidate in candidates
        if candidate in available
            return Symbol(candidate)
        end
    end
    throw(ArgumentError("none of the expected columns were found: $(join(candidates, ", "))"))
end

function summarize_results(results_df::DataFrame, counts_matrix::SparseMatrixCSC{Int,Int}; label::AbstractString, padj_cutoff::Real=0.1, low_count_cutoff::Real=12.0)
    padj_column = resolve_column_name(results_df, ["padj"])
    lfc_column = resolve_column_name(results_df, ["log2_fold_change", "log2FoldChange"])
    base_mean_column = resolve_column_name(results_df, ["base_mean", "baseMean"])

    padj_values = results_df[!, padj_column]
    lfc_values = results_df[!, lfc_column]
    base_mean_values = results_df[!, base_mean_column]

    significant_rows = [index for index in eachindex(padj_values) if !ismissing(padj_values[index]) && Float64(padj_values[index]) <= padj_cutoff]
    up_count = count(index -> Float64(lfc_values[index]) > 0, significant_rows)
    down_count = count(index -> Float64(lfc_values[index]) < 0, significant_rows)
    low_count = count(value -> !ismissing(value) && Float64(value) < low_count_cutoff, base_mean_values)
    outlier_count = count(index -> ismissing(padj_values[index]) && !ismissing(base_mean_values[index]) && Float64(base_mean_values[index]) >= low_count_cutoff, eachindex(padj_values))
    total_nonzero = count(vec(sum(counts_matrix, dims=2)) .> 0)

    println()
    println(label)
    println("-------")
    println("out of $(total_nonzero) with nonzero total read count")
    println("adjusted p-value < $(padj_cutoff)")
    println("LFC > 0 (up)       : $(up_count), $(round(100 * up_count / max(total_nonzero, 1), digits=1))%")
    println("LFC < 0 (down)     : $(down_count), $(round(100 * down_count / max(total_nonzero, 1), digits=1))%")
    println("outliers [1]       : $(outlier_count), $(round(100 * outlier_count / max(total_nonzero, 1), digits=1))%")
    println("low counts [2]     : $(low_count), $(round(100 * low_count / max(total_nonzero, 1), digits=1))%")
    println("(mean count < $(low_count_cutoff))")
    println("[1] see 'cooksCutoff' argument of ?results")
    println("[2] see 'independentFiltering' argument of ?results")
end

function compare_results(reference_df::DataFrame, julia_df::DataFrame)
    ref_gene_column = resolve_column_name(reference_df, ["gene_id"])
    julia_gene_column = resolve_column_name(julia_df, ["gene_id"])
    ref_lfc_column = resolve_column_name(reference_df, ["log2FoldChange", "log2_fold_change"])
    julia_lfc_column = resolve_column_name(julia_df, ["log2_fold_change", "log2FoldChange"])
    ref_pvalue_column = resolve_column_name(reference_df, ["pvalue"])
    julia_pvalue_column = resolve_column_name(julia_df, ["pvalue"])
    ref_padj_column = resolve_column_name(reference_df, ["padj"])
    julia_padj_column = resolve_column_name(julia_df, ["padj"])

    ref_lookup = Dict(String(reference_df[row, ref_gene_column]) => row for row in 1:nrow(reference_df))
    julia_lookup = Dict(String(julia_df[row, julia_gene_column]) => row for row in 1:nrow(julia_df))

    common_genes = intersect(collect(keys(ref_lookup)), collect(keys(julia_lookup)))
    ref_lfc_all = Float64[]
    julia_lfc_all = Float64[]
    ref_lfc_padj = Float64[]
    julia_lfc_padj = Float64[]
    ref_pvalue_all = Float64[]
    julia_pvalue_all = Float64[]

    for gene_id in common_genes
        ref_row = reference_df[ref_lookup[gene_id], :]
        julia_row = julia_df[julia_lookup[gene_id], :]

        ref_lfc_value = ref_row[ref_lfc_column]
        julia_lfc_value = julia_row[julia_lfc_column]
        if !ismissing(ref_lfc_value) && !ismissing(julia_lfc_value)
            rv = Float64(ref_lfc_value)
            jv = Float64(julia_lfc_value)
            if isfinite(rv) && isfinite(jv)
                push!(ref_lfc_all, rv)
                push!(julia_lfc_all, jv)
                if !ismissing(ref_row[ref_padj_column]) && !ismissing(julia_row[julia_padj_column])
                    push!(ref_lfc_padj, rv)
                    push!(julia_lfc_padj, jv)
                end
            end
        end

        ref_pvalue_value = ref_row[ref_pvalue_column]
        julia_pvalue_value = julia_row[julia_pvalue_column]
        if !ismissing(ref_pvalue_value) && !ismissing(julia_pvalue_value)
            rp = Float64(ref_pvalue_value)
            jp = Float64(julia_pvalue_value)
            if isfinite(rp) && isfinite(jp)
                push!(ref_pvalue_all, rp)
                push!(julia_pvalue_all, jp)
            end
        end
    end

    if !isempty(ref_lfc_all)
        println()
        println("Reference comparison")
        println("--------------------")
        println("LFC correlation (finite common, n=$(length(ref_lfc_all)))      : $(round(cor(ref_lfc_all, julia_lfc_all), digits=4))")
        if !isempty(ref_lfc_padj)
            println("LFC correlation (both padj finite, n=$(length(ref_lfc_padj))) : $(round(cor(ref_lfc_padj, julia_lfc_padj), digits=4))")
        end
        if !isempty(ref_pvalue_all)
            println("p-value correlation (finite common, n=$(length(ref_pvalue_all)))  : $(round(cor(ref_pvalue_all, julia_pvalue_all), digits=4))")
        end
    end
end

function compare_dispersions(reference_df::DataFrame, gene_ids::Vector{String}, dds::BT.DESeqDataSet; top_n::Integer=15)
    if dds.gene_wise_dispersions === nothing || dds.dispersion_fit === nothing || dds.dispersions === nothing
        println()
        println("Dispersion comparison")
        println("---------------------")
        println("Julia dispersion intermediates are unavailable; skipping per-gene comparison.")
        return
    end

    ref_gene_column = resolve_column_name(reference_df, ["gene_id"])
    ref_dispersion_column = resolve_column_name(reference_df, ["dispersion"])
    available_columns = Set(String.(names(reference_df)))
    ref_gene_wise_column = "dispGeneEst" in available_columns ? :dispGeneEst : ref_dispersion_column
    ref_trend_column = "dispFit" in available_columns ? :dispFit : ref_dispersion_column
    reference_lookup = Dict(String(reference_df[row, ref_gene_column]) => row for row in 1:nrow(reference_df))

    comparison_gene_ids = String[]
    reference_dispersions = Float64[]
    reference_gene_wise = Float64[]
    reference_trend = Float64[]
    julia_gene_wise = Float64[]
    julia_trend = Float64[]
    julia_map = Float64[]
    map_log_error = Float64[]
    gene_wise_log_error = Float64[]
    trend_log_error = Float64[]

    for (index, gene_id) in enumerate(gene_ids)
        reference_row = get(reference_lookup, gene_id, 0)
        reference_row == 0 && continue

        reference_value = reference_df[reference_row, ref_dispersion_column]
        ismissing(reference_value) && continue
        reference_gene_wise_value = reference_df[reference_row, ref_gene_wise_column]
        reference_trend_value = reference_df[reference_row, ref_trend_column]

        julia_gene_wise_value = Float64(dds.gene_wise_dispersions[index])
        julia_trend_value = Float64(dds.dispersion_fit[index])
        julia_map_value = Float64(dds.dispersions[index])
        reference_dispersion = Float64(reference_value)
        reference_gene_wise_dispersion = ismissing(reference_gene_wise_value) ? reference_dispersion : Float64(reference_gene_wise_value)
        reference_trend_dispersion = ismissing(reference_trend_value) ? reference_dispersion : Float64(reference_trend_value)

        push!(comparison_gene_ids, gene_id)
        push!(reference_dispersions, reference_dispersion)
        push!(reference_gene_wise, reference_gene_wise_dispersion)
        push!(reference_trend, reference_trend_dispersion)
        push!(julia_gene_wise, julia_gene_wise_value)
        push!(julia_trend, julia_trend_value)
        push!(julia_map, julia_map_value)
        push!(map_log_error, abs(log(julia_map_value) - log(reference_dispersion)))
        push!(gene_wise_log_error, abs(log(julia_gene_wise_value) - log(reference_gene_wise_dispersion)))
        push!(trend_log_error, abs(log(julia_trend_value) - log(reference_trend_dispersion)))
    end

    safe_corr(x::Vector{Float64}, y::Vector{Float64}) = begin
        mask = isfinite.(x) .& isfinite.(y)
        count(mask) > 1 ? cor(x[mask], y[mask]) : NaN
    end

    comparison_df = DataFrame(
        gene_id = comparison_gene_ids,
        r_dispersion = reference_dispersions,
        r_gene_wise = reference_gene_wise,
        r_trend = reference_trend,
        julia_gene_wise = julia_gene_wise,
        julia_trend = julia_trend,
        julia_map = julia_map,
        abs_log_map_error = map_log_error,
        abs_log_gene_wise_error = gene_wise_log_error,
        abs_log_trend_error = trend_log_error,
    )
    sort!(comparison_df, :abs_log_map_error, rev=true)
    CSV.write(JULIA_DISPERSIONS_PATH, comparison_df; missingstring="NA")

    println()
    println("Dispersion comparison")
    println("---------------------")
    println("common genes      : $(nrow(comparison_df))")
    println("gene-wise correlation : $(round(safe_corr(comparison_df.r_gene_wise, comparison_df.julia_gene_wise), digits=4))")
    println("trend correlation     : $(round(safe_corr(comparison_df.r_trend, comparison_df.julia_trend), digits=4))")
    println("map correlation       : $(round(safe_corr(comparison_df.r_dispersion, comparison_df.julia_map), digits=4))")
    println("worst MAP divergences :")
    show(first(sort(comparison_df, :abs_log_map_error, rev=true), min(top_n, nrow(comparison_df))); allrows=true, allcols=true)
    println()
    println("worst gene-wise divergences :")
    show(first(sort(comparison_df, :abs_log_gene_wise_error, rev=true), min(top_n, nrow(comparison_df))); allrows=true, allcols=true)
    println()
    println("worst trend divergences :")
    show(first(sort(comparison_df, :abs_log_trend_error, rev=true), min(top_n, nrow(comparison_df))); allrows=true, allcols=true)
    println()
end

println("Saved data for Julia.")
counts_sparse, gene_ids, sample_ids = read_counts_csv(COUNTS_PATH)
coldata, design = read_coldata_csv(COLDATA_PATH, sample_ids)

println("Running Julia DESeq2 (poscounts mode)...")
println("converting counts to integer mode")

dds = BT.DESeqDataSetFromMatrix(counts_sparse, coldata, design; gene_ids=gene_ids, sample_ids=sample_ids)

println("using pre-existing size factors")
dds = BT.estimateSizeFactors(dds; type=:poscounts)

println("estimating dispersions")
# dds = BT.DESeq(dds; sfType=:poscounts, fitType=:local, reference_level=:A, target_level=:B)
dds = BT.DESeq(dds; sfType=:poscounts, fitType=:parametric, reference_level=:A, target_level=:B)

println("fitting model and testing")

r_results = isfile(R_RESULTS_PATH) ? read_result_table(R_RESULTS_PATH) : nothing
r_dispersions = isfile(R_DISPERSIONS_PATH) ? read_result_table(R_DISPERSIONS_PATH) : nothing

julia_results = BT.results(dds, [:condition, :B, :A])

CSV.write(JULIA_RESULTS_PATH, julia_results; missingstring="NA")

summarize_results(julia_results, counts_sparse; label="Julia results")

if r_results !== nothing
    summarize_results(r_results, counts_sparse; label="R reference")
    compare_results(r_results, julia_results)
end

if r_dispersions !== nothing
    compare_dispersions(r_dispersions, gene_ids, dds)
end

println()
println("Saved Julia results.")