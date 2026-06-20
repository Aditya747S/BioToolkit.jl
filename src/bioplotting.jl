module BioPlotting

using Statistics
using LinearAlgebra
using Plots
using Distributions

using ..GWAS: GWASResult, MetaAnalysisResult
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id
import ..BioToolkit: ProvenanceContext, analysis_result_fields, analysis_result_summary, container_provenance_summary, ensure_provenance_id!, provenance_result!, provenance_parent_ids, register_container_provenance!, register_provenance!

export VolcanoPoint, VolcanoPlotResult, MAPoint, MAPlotResult, ClusteredHeatmapResult
export ManhattanPoint, ManhattanPlotResult, QQPoint, QQPlotResult, ForestPoint, ForestPlotResult
export volcano_data, volcano_plot, ma_data, ma_plot, clustered_heatmap, export_plot
export manhattan_data, manhattan_plot, qq_data, qq_plot, gwas_forest_plot

struct VolcanoPoint
    gene_id::String
    log2_fold_change::Float64
    pvalue::Float64
    padj::Float64
    negative_log10_pvalue::Float64
    negative_log10_padj::Float64
    category::Symbol
    labeled::Bool
end

struct MAPoint
    gene_id::String
    abundance::Float64
    log2_fold_change::Float64
    pvalue::Float64
    padj::Float64
    category::Symbol
end

struct VolcanoPlotResult <: AbstractAnalysisResult
    points::Vector{VolcanoPoint}
    summary::NamedTuple
    figure::Any
    metadata::Dict{Symbol,Any}
    provenance::ResultProvenance
end

VolcanoPlotResult(points, summary, figure, metadata) =
    VolcanoPlotResult(points, summary, figure, metadata, provenance_record("VolcanoPlotResult", "bioplotting"))

struct MAPlotResult <: AbstractAnalysisResult
    points::Vector{MAPoint}
    summary::NamedTuple
    figure::Any
    metadata::Dict{Symbol,Any}
    provenance::ResultProvenance
end

MAPlotResult(points, summary, figure, metadata) =
    MAPlotResult(points, summary, figure, metadata, provenance_record("MAPlotResult", "bioplotting"))

struct ClusteredHeatmapResult <: AbstractAnalysisResult
    matrix::Matrix{Float64}
    row_order::Vector{Int}
    column_order::Vector{Int}
    row_labels::Vector{String}
    column_labels::Vector{String}
    figure::Any
    metadata::Dict{Symbol,Any}
    provenance::ResultProvenance
end

ClusteredHeatmapResult(matrix, row_order, column_order, row_labels, column_labels, figure, metadata) =
    ClusteredHeatmapResult(matrix, row_order, column_order, row_labels, column_labels, figure, metadata, provenance_record("ClusteredHeatmapResult", "bioplotting"))

struct ManhattanPoint
    snp_id::String
    chromosome::String
    position::Int
    cumulative_position::Float64
    pvalue::Float64
    negative_log10_pvalue::Float64
    category::Symbol
    labeled::Bool
end

struct QQPoint
    expected::Float64
    observed::Float64
    pvalue::Float64
end

struct ForestPoint
    label::String
    beta::Float64
    standard_error::Float64
    lower::Float64
    upper::Float64
    pvalue::Float64
    qvalue::Float64
end

struct ManhattanPlotResult <: AbstractAnalysisResult
    points::Vector{ManhattanPoint}
    summary::NamedTuple
    figure::Any
    metadata::Dict{Symbol,Any}
    provenance::ResultProvenance
end

ManhattanPlotResult(points, summary, figure, metadata) =
    ManhattanPlotResult(points, summary, figure, metadata, provenance_record("ManhattanPlotResult", "bioplotting"))

struct QQPlotResult <: AbstractAnalysisResult
    points::Vector{QQPoint}
    summary::NamedTuple
    figure::Any
    metadata::Dict{Symbol,Any}
    provenance::ResultProvenance
end

QQPlotResult(points, summary, figure, metadata) =
    QQPlotResult(points, summary, figure, metadata, provenance_record("QQPlotResult", "bioplotting"))

struct ForestPlotResult <: AbstractAnalysisResult
    points::Vector{ForestPoint}
    summary::NamedTuple
    figure::Any
    metadata::Dict{Symbol,Any}
    provenance::ResultProvenance
end

ForestPlotResult(points, summary, figure, metadata) =
    ForestPlotResult(points, summary, figure, metadata, provenance_record("ForestPlotResult", "bioplotting"))

@inline function _plot_result_metadata(metadata::AbstractDict=Dict{Symbol,Any}())
    metadata_copy = Dict{Symbol,Any}(metadata)
    ensure_provenance_id!(metadata_copy)
    return metadata_copy
end

VolcanoPlotResult(points::Vector{VolcanoPoint}, summary::NamedTuple, figure; metadata::AbstractDict=Dict{Symbol,Any}()) = VolcanoPlotResult(points, summary, figure, _plot_result_metadata(metadata))
MAPlotResult(points::Vector{MAPoint}, summary::NamedTuple, figure; metadata::AbstractDict=Dict{Symbol,Any}()) = MAPlotResult(points, summary, figure, _plot_result_metadata(metadata))
ClusteredHeatmapResult(matrix::Matrix{Float64}, row_order::Vector{Int}, column_order::Vector{Int}, row_labels::Vector{String}, column_labels::Vector{String}, figure; metadata::AbstractDict=Dict{Symbol,Any}()) = ClusteredHeatmapResult(matrix, row_order, column_order, row_labels, column_labels, figure, _plot_result_metadata(metadata))
ManhattanPlotResult(points::Vector{ManhattanPoint}, summary::NamedTuple, figure; metadata::AbstractDict=Dict{Symbol,Any}()) = ManhattanPlotResult(points, summary, figure, _plot_result_metadata(metadata))
QQPlotResult(points::Vector{QQPoint}, summary::NamedTuple, figure; metadata::AbstractDict=Dict{Symbol,Any}()) = QQPlotResult(points, summary, figure, _plot_result_metadata(metadata))
ForestPlotResult(points::Vector{ForestPoint}, summary::NamedTuple, figure; metadata::AbstractDict=Dict{Symbol,Any}()) = ForestPlotResult(points, summary, figure, _plot_result_metadata(metadata))

analysis_result_fields(::Type{VolcanoPlotResult}) = (:points, :summary, :figure)
analysis_result_fields(::Type{MAPlotResult}) = (:points, :summary, :figure)
analysis_result_fields(::Type{ClusteredHeatmapResult}) = (:matrix, :row_order, :column_order, :row_labels, :column_labels, :figure)
analysis_result_fields(::Type{ManhattanPlotResult}) = (:points, :summary, :figure)
analysis_result_fields(::Type{QQPlotResult}) = (:points, :summary, :figure)
analysis_result_fields(::Type{ForestPlotResult}) = (:points, :summary, :figure)

function Base.show(io::IO, result::VolcanoPlotResult)
    print(io, analysis_result_summary(result), " | ", container_provenance_summary(result))
end

function Base.show(io::IO, result::MAPlotResult)
    print(io, analysis_result_summary(result), " | ", container_provenance_summary(result))
end

function Base.show(io::IO, result::ClusteredHeatmapResult)
    print(io, analysis_result_summary(result), " | ", container_provenance_summary(result))
end

function Base.show(io::IO, result::ManhattanPlotResult)
    print(io, analysis_result_summary(result), " | ", container_provenance_summary(result))
end

function Base.show(io::IO, result::QQPlotResult)
    print(io, analysis_result_summary(result), " | ", container_provenance_summary(result))
end

function Base.show(io::IO, result::ForestPlotResult)
    print(io, analysis_result_summary(result), " | ", container_provenance_summary(result))
end

@inline function _plot_result_finalize(result, _ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, parents, operation::AbstractString, parameters)
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

function _safe_probability(value::Real)
    finite = isfinite(Float64(value)) ? Float64(value) : 1.0
    return clamp(finite, eps(Float64), 1.0)
end

function _top_labels(results, limit::Integer)
    limit <= 0 && return Set{String}()
    ranked = sort(collect(results); by = result -> (isfinite(result.padj) ? result.padj : 1.0, -abs(result.log2_fold_change)))
    return Set(result.gene_id for result in Iterators.take(ranked, min(limit, length(ranked))))
end

function volcano_data(results; lfc_cutoff::Real=1.0, fdr_cutoff::Real=0.05, label_top::Integer=10, threaded::Bool=true)
    labels = _top_labels(results, label_top)
    points = Vector{VolcanoPoint}(undef, length(results))
    if threaded && length(results) > 1 && Threads.nthreads() > 1
        Threads.@threads for index in eachindex(results)
            result = results[index]
            pvalue = _safe_probability(result.pvalue)
            padj = _safe_probability(result.padj)
            category = padj <= fdr_cutoff && result.log2_fold_change >= lfc_cutoff ? :up : padj <= fdr_cutoff && result.log2_fold_change <= -lfc_cutoff ? :down : :neutral
            points[index] = VolcanoPoint(result.gene_id, Float64(result.log2_fold_change), pvalue, padj, -log10(pvalue), -log10(padj), category, result.gene_id in labels)
        end
    else
        for index in eachindex(results)
            result = results[index]
            pvalue = _safe_probability(result.pvalue)
            padj = _safe_probability(result.padj)
            category = padj <= fdr_cutoff && result.log2_fold_change >= lfc_cutoff ? :up : padj <= fdr_cutoff && result.log2_fold_change <= -lfc_cutoff ? :down : :neutral
            points[index] = VolcanoPoint(result.gene_id, Float64(result.log2_fold_change), pvalue, padj, -log10(pvalue), -log10(padj), category, result.gene_id in labels)
        end
    end
    summary = (
        total = length(points),
        significant = count(point -> point.category != :neutral, points),
        upregulated = count(point -> point.category == :up, points),
        downregulated = count(point -> point.category == :down, points))

    return points, summary
end

function ma_data(results; lfc_cutoff::Real=1.0, fdr_cutoff::Real=0.05, threaded::Bool=true)
    points = Vector{MAPoint}(undef, length(results))
    if threaded && length(results) > 1 && Threads.nthreads() > 1
        Threads.@threads for index in eachindex(results)
            result = results[index]
            pvalue = _safe_probability(result.pvalue)
            padj = _safe_probability(result.padj)
            abundance = log2(max(Float64(result.base_mean), 0.0) + 1.0)
            category = padj <= fdr_cutoff && abs(result.log2_fold_change) >= lfc_cutoff ? :significant : :background
            points[index] = MAPoint(result.gene_id, abundance, Float64(result.log2_fold_change), pvalue, padj, category)
        end
    else
        for index in eachindex(results)
            result = results[index]
            pvalue = _safe_probability(result.pvalue)
            padj = _safe_probability(result.padj)
            abundance = log2(max(Float64(result.base_mean), 0.0) + 1.0)
            category = padj <= fdr_cutoff && abs(result.log2_fold_change) >= lfc_cutoff ? :significant : :background
            points[index] = MAPoint(result.gene_id, abundance, Float64(result.log2_fold_change), pvalue, padj, category)
        end
    end
    summary = (
        total = length(points),
        significant = count(point -> point.category == :significant, points))

    return points, summary
end

function _scale_matrix(matrix::AbstractMatrix{<:Real}; scale::Symbol=:row, threaded::Bool=true)
    values = Matrix{Float64}(matrix)
    scale == :none && return values
    if scale == :row
        if threaded && size(values, 1) > 1 && Threads.nthreads() > 1
            Threads.@threads for row in axes(values, 1)
                series = values[row, :]
                centered = series .- mean(series)
                spread = std(series)
                values[row, :] = spread > 0 ? centered ./ spread : centered .* 0
            end
        else
            for row in axes(values, 1)
                series = values[row, :]
                centered = series .- mean(series)
                spread = std(series)
                values[row, :] = spread > 0 ? centered ./ spread : centered .* 0
            end
        end
    elseif scale == :column
        if threaded && size(values, 2) > 1 && Threads.nthreads() > 1
            Threads.@threads for column in axes(values, 2)
                series = values[:, column]
                centered = series .- mean(series)
                spread = std(series)
                values[:, column] = spread > 0 ? centered ./ spread : centered .* 0
            end
        else
            for column in axes(values, 2)
                series = values[:, column]
                centered = series .- mean(series)
                spread = std(series)
                values[:, column] = spread > 0 ? centered ./ spread : centered .* 0
            end
        end
    else
        throw(ArgumentError("scale must be :row, :column, or :none"))
    end
    return values
end

function _distance(a::AbstractVector{<:Real}, b::AbstractVector{<:Real}; metric::Symbol=:correlation)
    metric == :euclidean && return norm(Float64.(a) .- Float64.(b))
    metric == :correlation || throw(ArgumentError("metric must be :euclidean or :correlation"))
    ax = Float64.(a)
    bx = Float64.(b)
    da = ax .- mean(ax)
    db = bx .- mean(bx)
    sa = norm(da)
    sb = norm(db)
    (sa == 0 || sb == 0) && return 1.0
    return 1.0 - dot(da, db) / (sa * sb)
end

function _centroid_distance(matrix::AbstractMatrix{Float64}, left::Vector{Int}, right::Vector{Int}; metric::Symbol=:correlation)
    left_centroid = vec(mean(matrix[left, :], dims=1))
    right_centroid = vec(mean(matrix[right, :], dims=1))
    return _distance(left_centroid, right_centroid; metric=metric)
end

struct _ClusterNode
    members::Vector{Int}
    left::Union{Nothing,_ClusterNode}
    right::Union{Nothing,_ClusterNode}
end

function _leaf_order(node::_ClusterNode)
    node.left === nothing && node.right === nothing && return node.members
    left = node.left === nothing ? Int[] : _leaf_order(node.left)
    right = node.right === nothing ? Int[] : _leaf_order(node.right)
    return vcat(left, right)
end

function _hierarchical_order(matrix::AbstractMatrix{Float64}; metric::Symbol=:correlation)
    n = size(matrix, 1)
    n == 0 && return Int[]
    n == 1 && return [1]
    clusters = [_ClusterNode([index], nothing, nothing) for index in 1:n]
    while length(clusters) > 1
        best_left = 1
        best_right = 2
        best_distance = Inf
        for left in 1:length(clusters)-1
            for right in left+1:length(clusters)
                distance = _centroid_distance(matrix, clusters[left].members, clusters[right].members; metric=metric)
                if distance < best_distance
                    best_distance = distance
                    best_left = left
                    best_right = right
                end
            end
        end
        merged = _ClusterNode(vcat(clusters[best_left].members, clusters[best_right].members), clusters[best_left], clusters[best_right])
        deleteat!(clusters, best_right)
        deleteat!(clusters, best_left)
        push!(clusters, merged)
    end
    return _leaf_order(first(clusters))
end

function _reorder_labels(labels, order::Vector{Int})
    labels === nothing && return [string(index) for index in order]
    return [String(labels[index]) for index in order]
end

function _save_plots_figure(path::Union{Nothing,String}, figure)
    path === nothing && return nothing
    savefig(figure, path)
    return path
end

function volcano_plot(results; lfc_cutoff::Real=1.0, fdr_cutoff::Real=0.05, label_top::Integer=10, title::String="Volcano plot", size=(980, 700), save_path::Union{Nothing,String}=nothing, kwargs...)
    points, summary = volcano_data(results; lfc_cutoff=lfc_cutoff, fdr_cutoff=fdr_cutoff, label_top=label_top)
    figure = plot(; title=title, xlabel="log2 fold change", ylabel="-log10 p-value", size=size, legend=:outerright, background_color=:white, left_margin=15Plots.mm, bottom_margin=15Plots.mm, top_margin=10Plots.mm, right_margin=28Plots.mm, kwargs...)
    palette = Dict(:up => :firebrick, :down => :royalblue, :neutral => :gray)
    labels = Dict(:up => "Upregulated", :down => "Downregulated", :neutral => "Neutral")
    for category in (:neutral, :up, :down)
        selected = [point for point in points if point.category == category]
        isempty(selected) && continue
        scatter!(figure, [point.log2_fold_change for point in selected], [point.negative_log10_pvalue for point in selected]; color=palette[category], markersize=7, markerstrokewidth=0, label=labels[category])
    end
    vline!(figure, [-Float64(lfc_cutoff), Float64(lfc_cutoff)]; color=:gray, linestyle=:dash, label=false)
    hline!(figure, [-log10(_safe_probability(fdr_cutoff))]; color=:gray, linestyle=:dash, label=false)
    labeled = sort([point for point in points if point.labeled]; by = point -> point.negative_log10_pvalue, rev=true)
    for (index, point) in enumerate(Iterators.take(labeled, 8))
        offset_y = 0.28 + 0.16 * ((index - 1) % 3)
        annotate!(figure, point.log2_fold_change, point.negative_log10_pvalue + offset_y, text(point.gene_id, 7, :black, :center, :bottom; rotation=90))
    end
    _save_plots_figure(save_path, figure)
    result = VolcanoPlotResult(points, summary, figure)
    provenance_parameters = Dict{Any,Any}(kwargs)
    provenance_parameters[:lfc_cutoff] = lfc_cutoff
    provenance_parameters[:fdr_cutoff] = fdr_cutoff
    provenance_parameters[:label_top] = label_top
    provenance_parameters[:title] = title
    provenance_parameters[:size] = size
    _ctx = active_provenance_context()
    return _plot_result_finalize(result, _ctx, provenance_parent_ids(results), "volcano_plot", provenance_parameters)
end

function ma_plot(results; lfc_cutoff::Real=1.0, fdr_cutoff::Real=0.05, title::String="MA plot", size=(900, 660), save_path::Union{Nothing,String}=nothing, kwargs...)
    points, summary = ma_data(results; lfc_cutoff=lfc_cutoff, fdr_cutoff=fdr_cutoff)
    figure = plot(; title=title, xlabel="log2 abundance", ylabel="log2 fold change", size=size, legend=:topright, background_color=:white, left_margin=15Plots.mm, bottom_margin=15Plots.mm, top_margin=10Plots.mm, kwargs...)
    palette = Dict(:significant => :darkmagenta, :background => :slategray)
    labels = Dict(:significant => "Significant", :background => "Background")
    for category in (:background, :significant)
        selected = [point for point in points if point.category == category]
        isempty(selected) && continue
        scatter!(figure, [point.abundance for point in selected], [point.log2_fold_change for point in selected]; color=palette[category], markersize=7, markerstrokewidth=0, label=labels[category])
    end
    hline!(figure, [0.0]; color=:gray, linestyle=:dash, label=false)
    _save_plots_figure(save_path, figure)
    result = MAPlotResult(points, summary, figure)
    provenance_parameters = Dict{Any,Any}(kwargs)
    provenance_parameters[:lfc_cutoff] = lfc_cutoff
    provenance_parameters[:fdr_cutoff] = fdr_cutoff
    provenance_parameters[:title] = title
    provenance_parameters[:size] = size
    _ctx = active_provenance_context()
    return _plot_result_finalize(result, _ctx, provenance_parent_ids(results), "ma_plot", provenance_parameters)
end

function clustered_heatmap(matrix::AbstractMatrix{<:Real}; row_labels=nothing, column_labels=nothing, scale::Symbol=:row, metric::Symbol=:correlation, title::String="Clustered heatmap", size=(1080, 740), save_path::Union{Nothing,String}=nothing, kwargs...)
    scaled = _scale_matrix(matrix; scale=scale)
    row_order = _hierarchical_order(scaled; metric=metric)
    column_order = _hierarchical_order(permutedims(scaled); metric=metric)
    reordered = scaled[row_order, column_order]
    row_names = _reorder_labels(row_labels, row_order)
    column_names = _reorder_labels(column_labels, column_order)
    figure = heatmap(reordered; color=:balance, title=title, xlabel="Samples", ylabel="Features", size=size, colorbar=true, colorbar_title="Row-scaled z-score", clim=(-2.5, 2.5), left_margin=18Plots.mm, bottom_margin=18Plots.mm, top_margin=10Plots.mm, right_margin=28Plots.mm, kwargs...)
    xticks!(figure, (1:length(column_names), column_names))
    yticks!(figure, (1:length(row_names), row_names))
    plot!(figure; xrotation=45)
    _save_plots_figure(save_path, figure)
    result = ClusteredHeatmapResult(reordered, row_order, column_order, row_names, column_names, figure)
    provenance_parameters = Dict{Any,Any}(kwargs)
    provenance_parameters[:scale] = scale
    provenance_parameters[:metric] = metric
    provenance_parameters[:title] = title
    provenance_parameters[:size] = size
    _ctx = active_provenance_context()
    return _plot_result_finalize(result, _ctx, provenance_parent_ids(matrix), "clustered_heatmap", provenance_parameters)
end

function export_plot(result, save_path::String)
    result isa VolcanoPlotResult && result.figure !== nothing && savefig(result.figure, save_path)
    result isa MAPlotResult && result.figure !== nothing && savefig(result.figure, save_path)
    result isa ClusteredHeatmapResult && result.figure !== nothing && savefig(result.figure, save_path)

    return save_path
end

function _chromosome_key(chromosome::String)
    stripped = replace(lowercase(String(chromosome)), "chr" => "")
    parsed = tryparse(Int, stripped)
    return parsed === nothing ? typemax(Int) : parsed
end

function _chromosome_offsets(result::GWASResult)
    order = unique(result.chromosomes)
    sort!(order; by = _chromosome_key)
    offsets = Dict{String,Float64}()
    cumulative = 0.0
    for chromosome in order
        positions = [Float64(position) for (chrom, position) in zip(result.chromosomes, result.positions) if chrom == chromosome]
        max_position = isempty(positions) ? 0.0 : maximum(positions)
        offsets[chromosome] = cumulative
        cumulative += max(max_position, 1.0)
    end
    return offsets
end

function _chromosome_layout(result::GWASResult)
    order = unique(result.chromosomes)
    sort!(order; by = _chromosome_key)
    offsets = Dict{String,Float64}()
    centers = Dict{String,Float64}()
    cumulative = 0.0
    for chromosome in order
        positions = [Float64(position) for (chrom, position) in zip(result.chromosomes, result.positions) if chrom == chromosome]
        span = isempty(positions) ? 1.0 : max(maximum(positions) - minimum(positions), 1.0)
        offsets[chromosome] = cumulative
        centers[chromosome] = cumulative + span / 2
        cumulative += span
    end
    return order, offsets, centers, cumulative
end

function manhattan_data(result::GWASResult; pvalue_threshold::Real=5e-8, label_top::Integer=10, threaded::Bool=true)
    order, offsets, _, total_span = _chromosome_layout(result)
    ranked = sortperm(result.pvalue)
    labels = Set(result.snp_ids[index] for index in Iterators.take(ranked, min(label_top, length(ranked))))
    points = Vector{ManhattanPoint}(undef, length(result.snp_ids))
    if threaded && length(result.snp_ids) > 1 && Threads.nthreads() > 1
        Threads.@threads for index in eachindex(result.snp_ids)
            pvalue = clamp(isfinite(result.pvalue[index]) ? Float64(result.pvalue[index]) : 1.0, eps(Float64), 1.0)
            category = pvalue <= pvalue_threshold ? :significant : :background
            cumulative = get(offsets, result.chromosomes[index], 0.0) + Float64(result.positions[index])
            points[index] = ManhattanPoint(result.snp_ids[index], result.chromosomes[index], result.positions[index], cumulative, pvalue, -log10(pvalue), category, result.snp_ids[index] in labels)
        end
    else
        for index in eachindex(result.snp_ids)
            pvalue = clamp(isfinite(result.pvalue[index]) ? Float64(result.pvalue[index]) : 1.0, eps(Float64), 1.0)
            category = pvalue <= pvalue_threshold ? :significant : :background
            cumulative = get(offsets, result.chromosomes[index], 0.0) + Float64(result.positions[index])
            points[index] = ManhattanPoint(result.snp_ids[index], result.chromosomes[index], result.positions[index], cumulative, pvalue, -log10(pvalue), category, result.snp_ids[index] in labels)
        end
    end
    summary = (
        total = length(points),
        significant = count(point -> point.category == :significant, points),
        chromosomes = length(order),
        span = total_span)

    return points, summary
end

function manhattan_plot(result::GWASResult; pvalue_threshold::Real=5e-8, label_top::Integer=10)
    points, summary = manhattan_data(result; pvalue_threshold=pvalue_threshold, label_top=label_top)
    order, offsets, centers, total_span = _chromosome_layout(result)
    figure = plot(; size=(1200, 620), dpi=180, background_color=:white, foreground_color=:black, legend=false, grid=false, left_margin=6Plots.mm, right_margin=4Plots.mm, top_margin=5Plots.mm, bottom_margin=6Plots.mm, titlefontsize=14, guidefontsize=12, tickfontsize=10)
    xs = [point.cumulative_position for point in points]
    ys = [point.negative_log10_pvalue for point in points]
    alternating = Dict(order[index] => (isodd(index) ? "#4f6d7a" : "#9fb3c8") for index in eachindex(order))
    colors = [point.category == :significant ? "#d1495b" : alternating[point.chromosome] for point in points]
    scatter!(figure, xs, ys; color=colors, markersize=4.8, markerstrokewidth=0, xlabel="genomic position", ylabel="-log10(p)", title="Manhattan plot")
    xticks!(figure, ([get(centers, chromosome, 0.0) for chromosome in order], order))
    xlims!(figure, (0.0, max(total_span, 1.0)))
    ypeak = isempty(ys) ? 0.0 : maximum(ys)
    ylims!(figure, (0.0, max(ypeak * 1.08, -log10(clamp(Float64(pvalue_threshold), eps(Float64), 1.0)) * 1.15)))
    if any(point.labeled for point in points)
        labeled = [point for point in points if point.labeled]
        annotate!(figure, [(point.cumulative_position, point.negative_log10_pvalue, text(point.snp_id, 7, "#000000", halign=:left)) for point in labeled])
    end
    hline!(figure, [-log10(clamp(Float64(pvalue_threshold), eps(Float64), 1.0))]; linestyle=:dash, color="#7d1f2f", linewidth=1.5)
    result_plot = ManhattanPlotResult(points, summary, figure)
    _ctx = active_provenance_context()


    return _plot_result_finalize(result_plot, _ctx, provenance_parent_ids(result), "manhattan_plot", (pvalue_threshold=pvalue_threshold, label_top=label_top))
end

function qq_data(result::GWASResult; threaded::Bool=true)
    pvalues = sort(clamp.(Float64.(result.pvalue), eps(Float64), 1.0))
    n = length(pvalues)
    points = Vector{QQPoint}(undef, n)
    if threaded && n > 1 && Threads.nthreads() > 1
        Threads.@threads for index in eachindex(pvalues)
            pvalue = pvalues[index]
            expected = -log10((index - 0.5) / n)
            observed = -log10(pvalue)
            points[index] = QQPoint(expected, observed, pvalue)
        end
    else
        for (index, pvalue) in enumerate(pvalues)
            expected = -log10((index - 0.5) / n)
            observed = -log10(pvalue)
            points[index] = QQPoint(expected, observed, pvalue)
        end
    end
    chisq = quantile.(Ref(Chisq(1)), 1 .- pvalues)
    lambda_gc = isempty(chisq) ? 1.0 : median(chisq) / 0.4549364
    summary = (
        total = n,
        lambda_gc = lambda_gc)

    return points, summary
end

function qq_plot(result::GWASResult)
    points, summary = qq_data(result)
    figure = plot(; size=(840, 760), dpi=180, background_color=:white, foreground_color=:black, legend=false, grid=false, left_margin=6Plots.mm, right_margin=4Plots.mm, top_margin=5Plots.mm, bottom_margin=6Plots.mm, titlefontsize=14, guidefontsize=12, tickfontsize=10, aspect_ratio=:equal)
    xs = [point.expected for point in points]
    ys = [point.observed for point in points]
    if !isempty(points)
        n = length(points)
        confidence_x = copy(xs)
        confidence_low = [max(-log10(quantile(Beta(index, n - index + 1), 0.025)), 0.0) for index in 1:n]
        confidence_high = [max(-log10(quantile(Beta(index, n - index + 1), 0.975)), 0.0) for index in 1:n]
        plot!(figure, confidence_x, confidence_high; fillrange=confidence_low, fillalpha=0.18, linecolor=:transparent, fillcolor="#9db4c0")
    end
    scatter!(figure, xs, ys; color="#173f5f", markersize=4.2, markerstrokewidth=0, xlabel="expected -log10(p)", ylabel="observed -log10(p)", title="QQ plot")
    max_value = isempty(points) ? 1.0 : max(maximum(xs), maximum(ys))
    plot!(figure, [0.0, max_value], [0.0, max_value]; color="#7a7a7a", linestyle=:dash, linewidth=1.4)
    result_plot = QQPlotResult(points, summary, figure)
    _ctx = active_provenance_context()


    return _plot_result_finalize(result_plot, _ctx, provenance_parent_ids(result), "qq_plot", NamedTuple())
end

function gwas_forest_plot(result::MetaAnalysisResult; threaded::Bool=true)
    points = Vector{ForestPoint}(undef, length(result.snp_ids))
    if threaded && length(result.snp_ids) > 1 && Threads.nthreads() > 1
        Threads.@threads for index in eachindex(result.snp_ids)
            beta = result.beta[index]
            stderr = result.standard_error[index]
            lower = beta - 1.96 * stderr
            upper = beta + 1.96 * stderr
            points[index] = ForestPoint(result.snp_ids[index], beta, stderr, lower, upper, result.pvalue[index], result.qvalue[index])
        end
    else
        for index in eachindex(result.snp_ids)
            beta = result.beta[index]
            stderr = result.standard_error[index]
            lower = beta - 1.96 * stderr
            upper = beta + 1.96 * stderr
            points[index] = ForestPoint(result.snp_ids[index], beta, stderr, lower, upper, result.pvalue[index], result.qvalue[index])
        end
    end
    order = sortperm([point.qvalue for point in points])
    selected = points[order]
    figure = plot(; size=(920, max(320, 42 * max(length(selected), 1))), dpi=180, background_color=:white, foreground_color=:black, legend=false, grid=false, left_margin=7Plots.mm, right_margin=5Plots.mm, top_margin=6Plots.mm, bottom_margin=6Plots.mm, titlefontsize=14, guidefontsize=12, tickfontsize=10)
    if !isempty(selected)
        y_positions = collect(length(selected):-1:1)
        scatter!(figure, [point.beta for point in selected], y_positions; xerror=[1.96 * point.standard_error for point in selected], color="#173f5f", markersize=5.5, markerstrokewidth=0)
        for (index, point) in enumerate(selected)
            plot!(figure, [point.lower, point.upper], [y_positions[index], y_positions[index]]; color="#173f5f", linewidth=2.2)
        end
        yticks!(figure, (y_positions, [point.label for point in selected]))
        vline!(figure, [0.0]; color="#7d1f2f", linestyle=:dash, linewidth=1.3)
        xlims!(figure, (minimum([point.lower for point in selected]) - 0.1, maximum([point.upper for point in selected]) + 0.1))
        ylims!(figure, (0.5, length(selected) + 0.5))
    end
    summary = (
        total = length(selected),
        study_count = result.study_count,
        random_effects = result.method == "random_effects")
    result_plot = ForestPlotResult(selected, summary, figure)
    _ctx = active_provenance_context()


    return _plot_result_finalize(result_plot, _ctx, provenance_parent_ids(result), "gwas_forest_plot", (threaded=threaded,))
end

end
