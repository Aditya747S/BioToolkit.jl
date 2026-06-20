# ==============================================================================
# bioconductor_compat.jl — Bioconductor-style compatibility utilities
#
# This module provides scientifically motivated API equivalents for commonly
# used Bioconductor workflows, layered on top of BioToolkit native primitives.
#
# References:
#   - Lun et al. (2016) Genome Biology 17:75 (scran pooling/deconvolution)
#   - McCarthy et al. (2017) F1000Research 6:748 (scater/scuttle QC patterns)
#   - Zappia et al. (2017) Genome Biology 18:174 (Splatter simulation)
#   - Lun & Smyth (2016) NAR 44(5):e45 (csaw window-based DB)
#   - Stark & Brown (2011) DiffBind package (differential binding workflow)
#   - Schep et al. (2017) Nature Methods 14:975-978 (chromVAR deviations)
#   - Aryee et al. (2014) Bioinformatics 30(10):1363-1369 (minfi)
#   - Hansen et al. (2012) Genome Biology 13:R83 (BSmooth/BSseq context)
#   - Smith et al. (2006) Analytical Chemistry 78(3):779-787 (XCMS)
#   - Rohart et al. (2017) PLoS Comput Biol 13(11):e1005752 (mixOmics)
#   - Argelaguet et al. (2018) Mol Syst Biol 14:e8124 (MOFA)
#   - Obenchain et al. (2014) Bioinformatics 30(14):2076-2078 (VariantAnnotation)
# ==============================================================================

module BioconductorCompat

using Base.Threads
using DataFrames
using Distributions
using LinearAlgebra
using Random
using SparseArrays
using Statistics

using ..BioPlotting: clustered_heatmap
using ..BioToolkit: annotate_variants
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_record, provenance_result!, register_provenance!, with_provenance

@inline function _register_bioc_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end
using ..DifferentialExpression: CountMatrix, benjamini_hochberg, differential_expression
using ..Epigenetics: MethylationCall, PeakSet, bin_methylation, compute_motif_deviations, differential_binding, differential_methylation
using ..GWAS: ebi_lookup
using ..Metabolomics: metabolomics_differential_abundance
using ..Proteomics: MassSpecExperiment, align_samples, detect_peaks
using ..SingleCell: SingleCellExperiment, annotate_cell_types, fit_singlecell_projection_model, project_singlecell
using ..SystemsBio: multi_omics_factor_analysis

export biocparallel_map
export scater_qc_metrics, scran_pool_size_factors, scuttle_aggregate_counts, singler_annotate, scmap_project, splatter_simulate_counts
export diffbind_differential_binding, csaw_window_differential_binding, chromvar_deviations, minfi_dmp, bsseq_dmr
export xcms_peak_workflow, lipidr_differential_abundance
export mixomics_factor_analysis, mofa2_factor_analysis
export variantannotation_annotate, varianttools_filter_variants, gwascat_lookup
export complexheatmap_payload, gviz_track_table

@inline function _bh(p::Vector{Float64})
    p = clamp.(p, 0.0, 1.0)
    return benjamini_hochberg(p)
end

"""
    biocparallel_map(f, collection; nworkers=Threads.nthreads())

Parallel map utility equivalent in spirit to BiocParallel backends.
"""
function biocparallel_map(f::Function, collection; nworkers::Integer=Threads.nthreads())
    ntasks = max(1, Int(nworkers))
    return collect(asyncmap(f, collection; ntasks=ntasks))
end

"""
    scater_qc_metrics(experiment; mito_prefixes=["MT-", "MT_", "mt-", "mt_"])

Compute per-cell QC metrics analogous to common scater/scuttle summaries.
"""
function scater_qc_metrics(experiment::SingleCellExperiment; mito_prefixes::AbstractVector{<:AbstractString}=["MT-", "MT_", "mt-", "mt_"])
    _ctx = active_provenance_context()

    counts = experiment.counts
    n_cells = size(counts, 2)

    total_counts = vec(sum(counts, dims=1))
    detected_features = if counts isa SparseMatrixCSC
        [counts.colptr[j + 1] - counts.colptr[j] for j in 1:n_cells]
    else
        [count(>(0), @view counts[:, j]) for j in 1:n_cells]
    end

    mito_idx = [i for (i, gene) in enumerate(experiment.gene_ids) if any(prefix -> startswith(gene, prefix), mito_prefixes)]
    mito_counts = isempty(mito_idx) ? zeros(Float64, n_cells) : Float64.(vec(sum(counts[mito_idx, :], dims=1)))
    pct_mito = 100.0 .* mito_counts ./ max.(Float64.(total_counts), 1.0)

    result = with_provenance(DataFrame(
        cell_id=copy(experiment.cell_ids),
        total_counts=Float64.(total_counts),
        detected_features=Int.(detected_features),
        mito_counts=mito_counts,
        pct_mito=pct_mito), "BioconductorCompatTable", "BioconductorCompat/scater_qc_metrics"; parameters=(cell_count=n_cells, mito_prefix_count=length(mito_prefixes)))
    return _register_bioc_result!(_ctx, result, "scater_qc_metrics"; parents=provenance_parent_ids(experiment), parameters=(cell_count=n_cells, n_mito=length(mito_idx)))
end

function _geometric_mean_positive(values::AbstractVector{<:Real})
    positive = Float64[v for v in values if isfinite(v) && v > 0]
    isempty(positive) && return 0.0
    return exp(mean(log.(positive)))
end

"""
    scran_pool_size_factors(experiment; pool_size=20, overlap=10)

Approximate scran-style deconvolution size factors via pooled median-ratio
normalization and least-squares deconvolution of pool factors.
"""
function scran_pool_size_factors(experiment::SingleCellExperiment; pool_size::Integer=20, overlap::Integer=10)
    _ctx = active_provenance_context()

    counts = experiment.counts
    n_genes, n_cells = size(counts)
    n_cells >= 2 || throw(ArgumentError("scran_pool_size_factors requires at least two cells"))

    k = clamp(Int(pool_size), 2, n_cells)
    step = clamp(Int(overlap), 1, k)

    libsize = Float64.(vec(sum(counts, dims=1)))
    order = sortperm(libsize)

    starts = collect(1:step:max(1, n_cells - k + 1))
    terminal = n_cells - k + 1
    terminal >= 1 && (isempty(starts) || last(starts) != terminal) && push!(starts, terminal)

    pools = [order[s:(s + k - 1)] for s in starts]
    n_pools = length(pools)
    n_pools > 0 || throw(ArgumentError("unable to construct pooling design"))

    pooled = Matrix{Float64}(undef, n_genes, n_pools)
    for (j, pool) in enumerate(pools)
        pooled[:, j] = vec(sum(counts[:, pool], dims=2))
    end

    reference = [_geometric_mean_positive(@view pooled[g, :]) for g in 1:n_genes]
    valid_genes = findall(>(0.0), reference)
    isempty(valid_genes) && return libsize ./ max(median(libsize), eps(Float64))

    pool_sf = ones(Float64, n_pools)
    for j in 1:n_pools
        ratios = Float64[]
        sizehint!(ratios, length(valid_genes))
        for g in valid_genes
            v = pooled[g, j]
            v > 0 || continue
            push!(ratios, v / reference[g])
        end
        pool_sf[j] = isempty(ratios) ? 1.0 : median(ratios)
    end

    A = zeros(Float64, n_pools, n_cells)
    for (j, pool) in enumerate(pools)
        @inbounds for c in pool
            A[j, c] = 1.0
        end
    end

    rhs = log.(max.(pool_sf, eps(Float64)))
    log_sf = A \ rhs
    sf = exp.(log_sf)
    sf[.!isfinite.(sf)] .= 1.0
    sf ./= max(median(sf), eps(Float64))
    return _register_bioc_result!(_ctx, sf, "scran_pool_size_factors"; parents=provenance_parent_ids(experiment), parameters=(pool_size=Int(pool_size), overlap=Int(overlap), n_cells=n_cells, n_pools=n_pools))
end

"""
    scuttle_aggregate_counts(experiment, groups)

Aggregate cell-level counts into pseudobulk group-level counts.
"""
function scuttle_aggregate_counts(experiment::SingleCellExperiment, groups::AbstractVector)
    length(groups) == length(experiment.cell_ids) || throw(DimensionMismatch("groups must match number of cells"))
    labels = String.(groups)
    unique_labels = sort!(unique(labels))

    counts = experiment.counts
    aggregated = zeros(Int, size(counts, 1), length(unique_labels))
    for (j, label) in enumerate(unique_labels)
        idx = findall(==(label), labels)
        isempty(idx) && continue
        aggregated[:, j] = Int.(round.(vec(sum(counts[:, idx], dims=2))))
    end

    result = CountMatrix(sparse(aggregated), copy(experiment.gene_ids), unique_labels)
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, result, "scuttle_aggregate_counts"; parents=provenance_parent_ids(experiment), parameters=(n_groups=length(unique_labels), n_cells=length(groups)))
end

"""
    singler_annotate(experiment; reference="HumanPrimaryCellAtlas", labels=nothing, reference_labels=nothing)

SingleR-style annotation wrapper using BioToolkit cell annotation internals.
"""
function singler_annotate(experiment::SingleCellExperiment; reference::Union{String,SingleCellExperiment}="HumanPrimaryCellAtlas", labels=nothing, reference_labels=nothing, top_n::Int=20)
    method = reference isa SingleCellExperiment ? :reference : :marker
    encoded_labels = labels
    if labels !== nothing
        if !(eltype(labels) <: Integer)
            label_to_id = Dict{String,Int}()
            encoded_labels = Vector{Int}(undef, length(labels))
            for (i, label) in enumerate(labels)
                key = String(label)
                encoded_labels[i] = get!(label_to_id, key, length(label_to_id) + 1)
            end
        else
            encoded_labels = Int.(labels)
        end
    end
    return annotate_cell_types(experiment; labels=encoded_labels, reference=reference, reference_labels=reference_labels, method=method, top_n=top_n)
end

"""
    scmap_project(query, reference; labels=nothing, n_components=20, reduction_name="scmap")

scmap-like projection with nearest-reference assignment in projected PCA space.
"""
function scmap_project(query::SingleCellExperiment, reference::SingleCellExperiment; labels=nothing, n_components::Integer=20, reduction_name::String="scmap")
    model = fit_singlecell_projection_model(reference; n_components=n_components, use_variable_features=false)
    reference_embedding = project_singlecell(reference, model; reduction_name="_scmap_reference", n_components=n_components)
    query_embedding = project_singlecell(query, model; reduction_name=reduction_name, n_components=n_components)

    ref_labels = labels === nothing ? ["ref_$(i)" for i in 1:size(reference_embedding, 1)] : String.(labels)
    length(ref_labels) == size(reference_embedding, 1) || throw(DimensionMismatch("labels must match reference cell count"))

    nearest_index = zeros(Int, size(query_embedding, 1))
    nearest_label = Vector{String}(undef, size(query_embedding, 1))

    for i in 1:size(query_embedding, 1)
        q = @view query_embedding[i, :]
        best_j = 1
        best_d = Inf
        for j in 1:size(reference_embedding, 1)
            d = sum(abs2, q .- @view(reference_embedding[j, :]))
            if d < best_d
                best_d = d
                best_j = j
            end
        end
        nearest_index[i] = best_j
        nearest_label[i] = ref_labels[best_j]
    end
    return (embedding=query_embedding, nearest_index=nearest_index, predicted_label=nearest_label)
end

"""
    splatter_simulate_counts(; n_genes=2000, n_cells=500, n_groups=2, de_prob=0.1, de_lfc=1.0, seed=1)

Splatter-inspired Gamma-Poisson single-cell simulator with group-specific DE and
mean-dependent dropout.
"""
function splatter_simulate_counts(; n_genes::Integer=2000, n_cells::Integer=500, n_groups::Integer=2, de_prob::Real=0.1, de_lfc::Real=1.0, dropout_mid::Real=1.5, dropout_shape::Real=1.0, seed::Integer=1)
    n_genes > 0 || throw(ArgumentError("n_genes must be positive"))
    n_cells > 0 || throw(ArgumentError("n_cells must be positive"))
    n_groups > 0 || throw(ArgumentError("n_groups must be positive"))

    rng = MersenneTwister(seed)
    group_ids = [mod(i - 1, n_groups) + 1 for i in 1:n_cells]

    base_mean = rand(rng, Gamma(2.0, 2.0), n_genes)
    dispersion = rand(rng, Gamma(2.0, 0.15), n_genes)
    lib_factor = rand(rng, LogNormal(0.0, 0.35), n_cells)

    de_mask = rand(rng, n_genes) .< clamp(Float64(de_prob), 0.0, 1.0)
    de_dir = rand(rng, Bool, n_genes)

    counts = Matrix{Int}(undef, n_genes, n_cells)
    for g in 1:n_genes
        phi = max(dispersion[g], 1e-6)
        for c in 1:n_cells
            mu = base_mean[g] * lib_factor[c]
            if de_mask[g] && group_ids[c] > 1
                sign = de_dir[g] ? 1.0 : -1.0
                mu *= 2.0^(sign * Float64(de_lfc))
            end
            mu = max(mu, 1e-8)
            lambda = rand(rng, Gamma(1.0 / phi, mu * phi))
            value = rand(rng, Poisson(lambda))

            p_drop = 1.0 / (1.0 + exp((log1p(mu) - dropout_mid) / max(Float64(dropout_shape), eps(Float64))))
            counts[g, c] = rand(rng) < p_drop ? 0 : Int(value)
        end
    end

    gene_ids = ["gene_$(i)" for i in 1:n_genes]
    cell_ids = ["cell_$(i)" for i in 1:n_cells]
    sce = SingleCellExperiment(counts, gene_ids, cell_ids; metadata=Dict("group" => group_ids, "splatter_like" => true))

    truth = DataFrame(
        gene_id=gene_ids,
        is_de=Bool.(de_mask),
        direction=[de_mask[i] ? (de_dir[i] ? "up" : "down") : "none" for i in 1:n_genes],
        base_mean=Float64.(base_mean))
    result = (experiment=sce, truth=truth)
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, result, "splatter_simulate_counts"; parents=String[], parameters=(n_genes=Int(n_genes), n_cells=Int(n_cells), n_groups=Int(n_groups), de_prob=Float64(de_prob), seed=Int(seed)))
end

"""
    diffbind_differential_binding(fragments_by_sample, peaks, design; kwargs...)

DiffBind-style differential peak binding wrapper.
"""
function diffbind_differential_binding(fragments_by_sample::AbstractDict, peaks::PeakSet, design::AbstractVector, kwargs...)
    de = differential_binding(fragments_by_sample, peaks, Symbol.(design); kwargs...)
    result = DataFrame(de)
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, result, "diffbind_differential_binding"; parents=String[], parameters=(n_samples=length(fragments_by_sample)))
end

@inline function _frag_prop(fragment, name::Symbol)
    hasproperty(fragment, name) || throw(ArgumentError("fragment record must provide property :$(name)"))
    return getproperty(fragment, name)
end

"""
    csaw_window_differential_binding(fragments_by_sample, design; window_size=200, step_size=50, min_count=5)

csaw-like sliding-window differential binding analysis.
"""
function csaw_window_differential_binding(fragments_by_sample::AbstractDict, design::AbstractVector; window_size::Integer=200, step_size::Integer=50, min_count::Integer=5)
    window_size > 0 || throw(ArgumentError("window_size must be positive"))
    step_size > 0 || throw(ArgumentError("step_size must be positive"))

    sample_ids = sort!(collect(String.(keys(fragments_by_sample))))
    length(design) == length(sample_ids) || throw(DimensionMismatch("design length must match number of samples"))

    chrom_span = Dict{String,Tuple{Int,Int}}()
    for sample_id in sample_ids
        for fragment in fragments_by_sample[sample_id]
            chrom = String(_frag_prop(fragment, :chrom))
            left = Int(_frag_prop(fragment, :left))
            right = Int(_frag_prop(fragment, :right))
            lo, hi = get(chrom_span, chrom, (typemax(Int), typemin(Int)))
            chrom_span[chrom] = (min(lo, left), max(hi, right))
        end
    end
    isempty(chrom_span) && throw(ArgumentError("fragments_by_sample has no fragments"))

    window_chrom = String[]
    window_left = Int[]
    window_right = Int[]
    window_id = String[]
    chrom_rows = Dict{String,UnitRange{Int}}()

    for chrom in sort!(collect(keys(chrom_span)))
        lo, hi = chrom_span[chrom]
        start_row = length(window_id) + 1
        for left in lo:step_size:hi
            right = left + Int(window_size) - 1
            push!(window_chrom, chrom)
            push!(window_left, left)
            push!(window_right, right)
            push!(window_id, "$(chrom):$(left)-$(right)")
        end
        chrom_rows[chrom] = start_row:length(window_id)
    end

    counts = zeros(Int, length(window_id), length(sample_ids))

    for (sample_col, sample_id) in enumerate(sample_ids)
        for fragment in fragments_by_sample[sample_id]
            chrom = String(_frag_prop(fragment, :chrom))
            left = Int(_frag_prop(fragment, :left))
            right = Int(_frag_prop(fragment, :right))
            rows = get(chrom_rows, chrom, 0:-1)
            isempty(rows) && continue

            row0 = first(rows)
            lo = window_left[row0]
            idx_start = max(1, Int(floor((left - lo) / step_size)) + 1)
            idx_stop = min(length(rows), Int(floor((right - lo) / step_size)) + 1)

            for idx in idx_start:idx_stop
                row = row0 + idx - 1
                overlap_left = max(left, window_left[row])
                overlap_right = min(right, window_right[row])
                overlap_left <= overlap_right || continue
                counts[row, sample_col] += 1
            end
        end
    end

    cm = CountMatrix(sparse(counts), window_id, sample_ids)
    de = differential_expression(cm, Symbol.(design); min_total=Int(min_count), sfType=:poscounts)
    tab = DataFrame(de)

    chrom = String[]
    start = Int[]
    stop = Int[]
    for id in tab.gene_id
        value = String(id)
        parts = split(value, ':')
        coord = split(parts[2], '-')
        push!(chrom, parts[1])
        push!(start, parse(Int, coord[1]))
        push!(stop, parse(Int, coord[2]))
    end
    tab[!, :chromosome] = chrom
    tab[!, :start] = start
    tab[!, :stop] = stop
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, tab, "csaw_window_differential_binding"; parents=String[], parameters=(n_samples=length(fragments_by_sample), window_size=Int(window_size), step_size=Int(step_size)))
end

"""
    chromvar_deviations(chromatin, motif_peak_membership)

chromVAR-like motif deviation wrapper returning motif-by-cell z-scores.
"""
function chromvar_deviations(chromatin, motif_peak_membership::AbstractDict)
    deviations = compute_motif_deviations(chromatin, motif_peak_membership)
    motifs = sort!(collect(String.(keys(deviations))))
    table = DataFrame(motif=motifs)

    cell_ids = chromatin.base.cell_ids
    for (cell_index, cell_id) in enumerate(cell_ids)
        table[!, Symbol(cell_id)] = [Float64(deviations[motif][cell_index]) for motif in motifs]
    end
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, table, "chromvar_deviations"; parents=String[], parameters=(n_motifs=length(motifs), n_cells=length(cell_ids)))
end

"""
    minfi_dmp(beta_matrix, design; feature_ids=nothing, prior_df=4.0)

minfi/limma-style differential methylation at probe level using M-values and
moderated t-statistics.
"""
function minfi_dmp(beta_matrix::AbstractMatrix{<:Real}, design::AbstractVector; feature_ids=nothing, prior_df::Real=4.0, offset::Real=1e-3)
    B = clamp.(Float64.(beta_matrix), Float64(offset), 1.0 - Float64(offset))
    n_features, n_samples = size(B)
    length(design) == n_samples || throw(DimensionMismatch("design length must match sample count"))

    groups = unique(Symbol.(design))
    length(groups) == 2 || throw(ArgumentError("minfi_dmp currently supports exactly two groups"))
    g1 = findall(==(groups[1]), Symbol.(design))
    g2 = findall(==(groups[2]), Symbol.(design))
    (isempty(g1) || isempty(g2)) && throw(ArgumentError("both groups must contain at least one sample"))
    (length(g1) < 2 || length(g2) < 2) && throw(ArgumentError("each group must contain at least two samples"))

    M = log2.(B ./ (1.0 .- B))
    n1, n2 = length(g1), length(g2)
    dof = max(n1 + n2 - 2, 1)

    delta_beta = vec(mean(B[:, g2], dims=2) .- mean(B[:, g1], dims=2))
    effect = vec(mean(M[:, g2], dims=2) .- mean(M[:, g1], dims=2))

    var1 = vec(var(M[:, g1], dims=2; corrected=true))
    var2 = vec(var(M[:, g2], dims=2; corrected=true))
    s2 = ((n1 - 1) .* var1 .+ (n2 - 1) .* var2) ./ dof
    finite_s2 = s2[isfinite.(s2)]
    isempty(finite_s2) && throw(ArgumentError("unable to estimate probe variance from provided beta matrix"))
    s2[.!isfinite.(s2)] .= median(finite_s2)

    s0_sq = median(s2[isfinite.(s2)])
    d0 = max(Float64(prior_df), 0.0)
    s2_post = (d0 * s0_sq .+ dof .* s2) ./ (d0 + dof)

    se = sqrt.(s2_post .* (1.0 / n1 + 1.0 / n2))
    tstat = effect ./ max.(se, eps(Float64))
    pvalue = 2.0 .* ccdf.(TDist(dof + d0), abs.(tstat))
    padj = _bh(Float64.(pvalue))

    ids = feature_ids === nothing ? ["cg_$(i)" for i in 1:n_features] : String.(feature_ids)
    length(ids) == n_features || throw(DimensionMismatch("feature_ids length must match number of rows"))

    result = DataFrame(
        feature_id=ids,
        delta_beta=delta_beta,
        mvalue_effect=effect,
        t_stat=tstat,
        pvalue=pvalue,
        padj=padj)
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, result, "minfi_dmp"; parents=provenance_parent_ids(beta_matrix), parameters=(n_features=n_features, n_samples=n_samples, prior_df=Float64(prior_df), n_sig=count(<(0.05), padj)))
end

"""
    bsseq_dmr(calls, design; bin_width=500, min_total=10)

BSseq-style DMR wrapper using binned methylation and beta-binomial testing.
"""
function bsseq_dmr(calls::AbstractVector{<:MethylationCall}, design::AbstractVector; bin_width::Integer=500, min_total::Integer=10, sample_metadata::DataFrame=DataFrame())
    experiment = bin_methylation(calls; bin_width=Int(bin_width), sample_metadata=sample_metadata)
    length(design) == length(experiment.sample_ids) || throw(DimensionMismatch("design length must match methylation sample count"))
    result = differential_methylation(experiment, Symbol.(design); min_total=Int(min_total))
    df = DataFrame(result)
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, df, "bsseq_dmr"; parents=provenance_parent_ids(calls), parameters=(bin_width=Int(bin_width), min_total=Int(min_total), n_dmrs=nrow(df)))
end

function _match_feature_index(mz::Float64, rt::Float64, feature_mz::Vector{Float64}, feature_rt::Vector{Float64}; mz_tolerance::Real, rt_tolerance::Real)
    for i in eachindex(feature_mz)
        if abs(mz - feature_mz[i]) <= mz_tolerance && abs(rt - feature_rt[i]) <= rt_tolerance
            return i
        end
    end
    return 0
end

"""
    xcms_peak_workflow(experiment; scales=[1.0,2.0,4.0], threshold=1.25, mz_tolerance=0.01, rt_tolerance=15.0)

XCMS-like peak detection/alignment summary returning a feature-by-sample peak table.
"""
function xcms_peak_workflow(experiment::MassSpecExperiment; scales::AbstractVector{<:Real}=[1.0, 2.0, 4.0], threshold::Real=1.25, align_method::Symbol=:dtw, mz_tolerance::Real=0.01, rt_tolerance::Real=15.0)
    _ctx = active_provenance_context()

    n_samples = length(experiment.spectra)
    n_samples > 0 || throw(ArgumentError("MassSpecExperiment must contain at least one spectrum"))

    sample_names = if nrow(experiment.design) == n_samples && :sample_id in names(experiment.design)
        String.(experiment.design.sample_id)
    else
        ["sample_$(i)" for i in 1:n_samples]
    end

    peak_results = [detect_peaks(spectrum; scales=scales, threshold=threshold) for spectrum in experiment.spectra]

    reference = experiment.spectra[1].intensity
    alignment_cost = zeros(Float64, n_samples)
    alignment_cost[1] = 0.0
    for i in 2:n_samples
        aligned = align_samples(reference, experiment.spectra[i].intensity; method=align_method)
        alignment_cost[i] = aligned.cost
    end

    feature_mz = Float64[]
    feature_rt = Float64[]
    feature_n = Int[]
    intensity = zeros(Float64, 0, n_samples)

    for sample_index in 1:n_samples
        for peak in peak_results[sample_index].peaks
            idx = _match_feature_index(peak.mz, peak.rt, feature_mz, feature_rt; mz_tolerance=mz_tolerance, rt_tolerance=rt_tolerance)
            if idx == 0
                push!(feature_mz, peak.mz)
                push!(feature_rt, peak.rt)
                push!(feature_n, 1)
                intensity = vcat(intensity, zeros(1, n_samples))
                idx = size(intensity, 1)
            else
                n_prev = feature_n[idx]
                feature_mz[idx] = (feature_mz[idx] * n_prev + peak.mz) / (n_prev + 1)
                feature_rt[idx] = (feature_rt[idx] * n_prev + peak.rt) / (n_prev + 1)
                feature_n[idx] = n_prev + 1
            end
            intensity[idx, sample_index] = max(intensity[idx, sample_index], peak.intensity)
        end
    end

    table = DataFrame(
        feature_id=["feature_$(i)" for i in 1:length(feature_mz)],
        mz=feature_mz,
        rt=feature_rt)
    for (j, name) in enumerate(sample_names)
        table[!, Symbol(name)] = intensity[:, j]
    end

    result = (peak_table=table, peak_results=peak_results, alignment_cost=alignment_cost)
    return _register_bioc_result!(_ctx, result, "xcms_peak_workflow"; parents=String[], parameters=(n_samples=n_samples, n_features=nrow(table), align_method=align_method))
end

"""
    lipidr_differential_abundance(matrix, groups; covariates=nothing)

lipidr-style differential lipid abundance wrapper.
"""
function lipidr_differential_abundance(matrix::AbstractMatrix{<:Real}, groups::AbstractVector; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing)
    result = metabolomics_differential_abundance(matrix, groups; covariates=covariates)
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, result, "lipidr_differential_abundance"; parents=provenance_parent_ids(matrix), parameters=(n_features=size(matrix,1), n_samples=size(matrix,2)))
end

"""
    mixomics_factor_analysis(assays; n_components=3)

mixOmics-style latent factor integration wrapper.
"""
function mixomics_factor_analysis(assays::Vector{<:AbstractMatrix{<:Real}}; n_components::Integer=3)
    result = multi_omics_factor_analysis(assays; n_factors=Int(n_components))
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, result, "mixomics_factor_analysis"; parents=String[], parameters=(n_components=Int(n_components), n_assays=length(assays)))
end

"""
    mofa2_factor_analysis(assays; n_factors=3)

MOFA2-style latent factor integration wrapper.
"""
function mofa2_factor_analysis(assays::Vector{<:AbstractMatrix{<:Real}}; n_factors::Integer=3)
    result = multi_omics_factor_analysis(assays; n_factors=Int(n_factors))
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, result, "mofa2_factor_analysis"; parents=String[], parameters=(n_factors=Int(n_factors), n_assays=length(assays)))
end

"""
    variantannotation_annotate(variant_records, gene_features; reference_sequences=nothing)

VariantAnnotation-style convenience wrapper around BioToolkit variant effect annotation.
"""
function variantannotation_annotate(variant_records::AbstractVector, gene_features::AbstractVector; reference_sequences=nothing)
    annotations = annotate_variants(variant_records, gene_features; reference_sequences=reference_sequences)
    result = DataFrame(annotations)
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, result, "variantannotation_annotate"; parents=String[], parameters=(n_variants=length(variant_records), n_features=length(gene_features)))
end

"""
    varianttools_filter_variants(variants; min_qual=0.0, pass_only=true, snv_only=false)

VariantTools-style filtering utility for parsed variant record vectors.
"""
function varianttools_filter_variants(variants::AbstractVector; min_qual::Real=0.0, pass_only::Bool=true, snv_only::Bool=false)
    out = Any[]
    for variant in variants
        qual = hasproperty(variant, :qual) ? getproperty(variant, :qual) : missing
        filt = hasproperty(variant, :filter) ? String(getproperty(variant, :filter)) : "PASS"
        ref = hasproperty(variant, :ref) ? String(getproperty(variant, :ref)) : ""
        alt = hasproperty(variant, :alt) ? String(getproperty(variant, :alt)) : ""

        pass_only && !(filt in ("PASS", ".", "")) && continue
        if min_qual > 0
            q = qual === missing ? NaN : Float64(qual)
            (isfinite(q) && q >= min_qual) || continue
        end
        snv_only && !(ncodeunits(ref) == 1 && ncodeunits(alt) == 1) && continue

        push!(out, variant)
    end
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, out, "varianttools_filter_variants"; parents=String[], parameters=(n_input=length(variants), n_output=length(out)))
end

"""
    gwascat_lookup(snp_ids; max_requests=50, max_concurrency=16)

GWAS Catalog-like query wrapper using EBI lookup service.
"""
function gwascat_lookup(snp_ids::AbstractVector{<:AbstractString}; max_requests::Int=50, max_concurrency::Int=16)
    result = ebi_lookup(snp_ids; max_requests=max_requests, max_concurrency=max_concurrency)
    _ctx = active_provenance_context()
    return _register_bioc_result!(_ctx, result, "gwascat_lookup"; parents=String[], parameters=(n_snps=length(snp_ids), max_requests=max_requests))
end

"""
    complexheatmap_payload(matrix; row_labels=nothing, column_labels=nothing, kwargs...)

ComplexHeatmap-style clustered heatmap payload wrapper.
"""
function complexheatmap_payload(matrix::AbstractMatrix{<:Real}; row_labels=nothing, column_labels=nothing, kwargs...)
    return clustered_heatmap(matrix; row_labels=row_labels, column_labels=column_labels, kwargs...)
end

"""
    gviz_track_table(chromosomes, starts, stops; scores=nothing, labels=nothing, track="signal")

Create a tidy genomic-track table for genome-browser style plotting.
"""
function gviz_track_table(chromosomes::AbstractVector, starts::AbstractVector{<:Integer}, stops::AbstractVector{<:Integer}; scores=nothing, labels=nothing, track::AbstractString="signal")
    n = length(chromosomes)
    length(starts) == n || throw(DimensionMismatch("starts length must match chromosomes"))
    length(stops) == n || throw(DimensionMismatch("stops length must match chromosomes"))

    score_values = scores === nothing ? fill(1.0, n) : Float64.(scores)
    length(score_values) == n || throw(DimensionMismatch("scores length must match chromosomes"))

    label_values = labels === nothing ? ["feature_$(i)" for i in 1:n] : String.(labels)
    length(label_values) == n || throw(DimensionMismatch("labels length must match chromosomes"))
    return DataFrame(
        track=fill(String(track), n),
        chromosome=String.(chromosomes),
        start=Int.(starts),
        stop=Int.(stops),
        score=score_values,
        label=label_values)
end

end
