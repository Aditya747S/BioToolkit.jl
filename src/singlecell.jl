module SingleCell

using CUDA
using JSON
using Mmap
using Optim
using OrdinaryDiffEq
using SparseArrays
using Statistics
using LinearAlgebra
using Random
using Graphs
using SimpleWeightedGraphs: SimpleWeightedGraph, SimpleWeightedEdge
using SpecialFunctions: erfc
using Plots: plot, scatter, scatter!, plot!, heatmap, quiver!, bar

using Serialization

using ..DifferentialExpression: CountMatrix, DEResult, benjamini_hochberg, calc_norm_factors, differential_expression, estimate_dispersions, filter_low_counts, vst
using ..BioToolkit: threaded_foreach, threaded_map_collect
using ..BioToolkit: AbstractAnalysisResult, ProvenanceContext, ResultProvenance, ThreadSafeProvenanceContext, active_provenance_context, analysis_result_summary, metadata_provenance, new_provenance_id, provenance_parent_ids, provenance_record, provenance_summary, provenance_result!, register_provenance!, stamp_provenance!, update_provenance!, with_provenance

# ---------------------------------------------------------------------------
# Provenance helper
# ---------------------------------------------------------------------------
"""
    _register_singlecell_result!(_ctx, result, operation; parameters=NamedTuple())

Zero-cost provenance registration for SingleCell operations.
When `_ctx === nothing` this is a no-op — no allocation, no SHA, no dict.
"""
@inline function _register_singlecell_result!(::Nothing, result, ::AbstractString; parameters=NamedTuple())
    return result
end

function _register_singlecell_result!(ctx::ProvenanceContext, result, operation::AbstractString;
        parents::AbstractVector{<:AbstractString}=String[],
        parameters=NamedTuple())
    register_provenance!(ctx, operation; parents=parents, parameters=parameters)
    return result
end
export SingleCellExperiment, SingleCellProjectionModel, count_matrix, normalize_counts, sctransform, find_variable_features, fit_singlecell_projection_model, project_singlecell, run_pca, find_neighbors, find_spatial_neighbors, run_umap, cluster_cells, find_clusters, calculate_pseudotime, find_cluster_markers, find_markers, summarize_clusters, cluster_marker_summary
export detect_doublets, integrate_batches, integrate_data, score_cell_cycle
export RNAVelocityResult, DynamicalRNAVelocityResult, attach_velocity_layers!, calculate_rna_velocity, calculate_dynamical_rna_velocity, plot_trajectory, plot_rna_velocity, plot_communication_network, plot_communication_pathways, plot_ligand_receptor_report, plot_rna_velocity_quiver, gpu_singlecell_experiment, materialize_singlecell_experiment, save_singlecell_experiment, load_singlecell_experiment, SingleCellArchive, save_singlecell_archive, load_singlecell_archive, archive_to_singlecell_experiment, SingleCellArchiveBrowser, browse_singlecell_archive, archive_layer, archive_embedding
export LigandReceptorPair, CellCommunicationResult, CellCommunicationNetwork, CommunicationPathwaySummary, LigandReceptorReport, SpatialMoranResult, default_ligand_receptor_pairs, find_cell_communication, communication_network, communication_pathway_summary, rank_ligand_receptor_report, attach_spatial_coords!, moran_i_test, find_spatially_variable_genes
export WNNResult, weighted_nearest_neighbors, PerturbationPredictionResult, predict_perturbation
export read_h5ad, write_h5ad, CellTypeAnnotationResult, annotate_cell_types, AmbientRNARemovalResult, remove_ambient_rna, SingleCellViewer, interactive_singlecell_viewer, lasso_select_cells, recluster_singlecell_viewer!, cell_hover_text, top_expressed_genes
export bayesspace_like_domains, spagc_like_domains, sparkx_spatial_de, spatialde_gp_de, spatial_trajectory_graph, niche_weighted_communication, cell2location_like_segmentation, mixscape_like_contrast
export spatial_markov_refine_domains, spatial_lr_permutation_test
export mixscape_multibatch_contrast, perturbseq_pseudobulk, milo_like_neighborhood_da, perturbation_synergy_scores

"""
    SingleCellExperiment

Container for single-cell count matrices, cell and gene identifiers, metadata,
dimensionality reductions, clustering assignments, selected variable features,
and cached neighbor graphs.
"""
mutable struct SingleCellExperiment{M<:AbstractMatrix{Int}}
    counts::M
    gene_ids::Vector{String}
    cell_ids::Vector{String}
    spatial_coords::Union{Nothing,Matrix{Float64}}
    metadata::Dict{String,Any}
    reductions::Dict{String,Matrix{Float64}}
    clusters::Dict{String,Vector{Int}}
    variable_features::Dict{String,Vector{String}}
    neighbors::Dict{String,Vector{Vector{Int}}}
end

@inline function _singlecell_metadata_copy(metadata::AbstractDict)
    return Dict{String,Any}(string(key) => value for (key, value) in metadata)
end

function _singlecell_metadata_with_provenance(metadata::AbstractDict; source::AbstractString, notes::AbstractVector{<:AbstractString}=String[], parameters::NamedTuple=NamedTuple(), status::Symbol=:ok)
    copied = _singlecell_metadata_copy(metadata)
    update_provenance!(copied; label="SingleCellExperiment", source=source, status=status, notes=notes, parameters=parameters)
    return copied
end

"""
    SingleCellProjectionModel

Reusable normalization, scaling, and PCA parameters for projecting query data
onto a fitted single-cell reference basis.
"""
mutable struct SingleCellProjectionModel
    feature_ids::Vector{String}
    scale_factor::Float64
    log_transform::Bool
    center_features::Bool
    scale_features::Bool
    feature_means::Vector{Float64}
    feature_scales::Vector{Float64}
    pca_loadings::Matrix{Float64}
end

"""
    SingleCellExperiment(counts, gene_ids, cell_ids; metadata=Dict())

Construct a single-cell experiment from a gene-by-cell integer count matrix.
The constructor validates dimensions, normalizes metadata keys to strings, and
initializes caches for reductions, clusters, variable features, spatial
coordinates, and neighbor graphs.
"""
function _normalize_spatial_coords(spatial_coords::AbstractMatrix{<:Real}, n_cells::Int)
    if size(spatial_coords, 1) == n_cells && size(spatial_coords, 2) >= 2
        return Matrix{Float64}(spatial_coords[:, 1:2])
    elseif size(spatial_coords, 2) == n_cells && size(spatial_coords, 1) >= 2
        return Matrix{Float64}(permutedims(spatial_coords[1:2, :]))
    else
        throw(DimensionMismatch("spatial_coords must have one dimension matching the number of cells and at least two coordinates"))
    end
end

function SingleCellExperiment(counts::AbstractMatrix{<:Integer}, gene_ids::AbstractVector{<:String}, cell_ids::AbstractVector{<:String}; metadata::AbstractDict=Dict{String,Any}(), spatial_coords::Union{Nothing,AbstractMatrix}=nothing)
    matrix = counts isa SparseMatrixCSC ? sparse(Int.(counts)) : counts isa CUDA.CuArray ? Int.(counts) : sparse(Int.(counts))
    size(matrix, 1) == length(gene_ids) || throw(ArgumentError("gene_ids must match the number of rows in counts"))
    size(matrix, 2) == length(cell_ids) || throw(ArgumentError("cell_ids must match the number of columns in counts"))
    standardized_coords = spatial_coords === nothing ? nothing : _normalize_spatial_coords(spatial_coords, length(cell_ids))
    metadata_copy = _singlecell_metadata_copy(metadata)
    if metadata_provenance(metadata_copy) === nothing
        stamp_provenance!(metadata_copy;
            label="SingleCellExperiment",
            source="constructor",
            notes=["constructed from in-memory counts"],
            parameters=(n_genes=length(gene_ids), n_cells=length(cell_ids), has_spatial_coords=spatial_coords !== nothing))
    end
    return SingleCellExperiment{typeof(matrix)}(matrix, String.(gene_ids), String.(cell_ids), standardized_coords, metadata_copy, Dict{String,Matrix{Float64}}(), Dict{String,Vector{Int}}(), Dict{String,Vector{String}}(), Dict{String,Vector{Vector{Int}}}())
end

function _show_singlecell_experiment(io::IO, experiment::SingleCellExperiment)
    print(io, "SingleCellExperiment(", length(experiment.gene_ids), " genes, ", length(experiment.cell_ids), " cells")
    if experiment.spatial_coords !== nothing
        print(io, ", spatial=", size(experiment.spatial_coords, 1), "x", size(experiment.spatial_coords, 2))
    end
    if !isempty(experiment.reductions)
        print(io, ", reductions=", join(sort!(collect(keys(experiment.reductions))), ", "))
    end
    if !isempty(experiment.clusters)
        print(io, ", clusters=", join(sort!(collect(keys(experiment.clusters))), ", "))
    end
    provenance = metadata_provenance(experiment.metadata)
    provenance !== nothing && print(io, ", provenance=", provenance_summary(provenance))
    print(io, ")")
end

function Base.show(io::IO, experiment::SingleCellExperiment)
    _show_singlecell_experiment(io, experiment)
end

function Base.show(io::IO, ::MIME"text/plain", experiment::SingleCellExperiment)
    _show_singlecell_experiment(io, experiment)
end

function attach_spatial_coords!(experiment::SingleCellExperiment, spatial_coords::AbstractMatrix{<:Real})
    experiment.spatial_coords = _normalize_spatial_coords(spatial_coords, length(experiment.cell_ids))
    return experiment
end

"""
    count_matrix(experiment)

Return the `CountMatrix` view used by the differential-expression layer.
"""
function count_matrix(experiment::SingleCellExperiment)
    counts = experiment.counts isa CUDA.CuArray ? Array(experiment.counts) : experiment.counts
    return CountMatrix(counts, experiment.gene_ids, experiment.cell_ids)
end

function _dense_counts(experiment::SingleCellExperiment)
    return Array{Float64}(experiment.counts)
end

function _matrix_triplets(matrix::SparseMatrixCSC{<:Real,Int})
    rows, cols, values = findnz(matrix)
    return rows, cols, Float64.(values)
end

function _matrix_triplets(matrix::CUDA.CuArray{<:Real,2})
    return _matrix_triplets(Array(matrix))
end

function _matrix_triplets(matrix::AbstractMatrix{<:Real})
    rows = Int[]
    cols = Int[]
    values = Float64[]
    nrows, ncols = size(matrix)
    for column in 1:ncols
        for row in 1:nrows
            value = Float64(matrix[row, column])
            value == 0.0 && continue
            push!(rows, row)
            push!(cols, column)
            push!(values, value)
        end
    end
    return rows, cols, values
end

function _sparse_row_mean_variance(matrix::SparseMatrixCSC{<:Real,Int})
    nrows, ncols = size(matrix)
    if nrows == 0
        return Float64[], Float64[]
    end

    row_sums = zeros(Float64, nrows)
    row_sumsq = zeros(Float64, nrows)
    rows, _, values = findnz(matrix)
    @inbounds for index in eachindex(values)
        row = rows[index]
        value = Float64(values[index])
        row_sums[row] += value
        row_sumsq[row] += value * value
    end

    inv_ncols = ncols == 0 ? 0.0 : 1 / ncols
    means = row_sums .* inv_ncols
    variances = max.(row_sumsq .* inv_ncols .- means .^ 2, 0.0)
    return means, variances
end

_dense_row_mean_variance(matrix::AbstractMatrix{<:Real}) = (
    vec(mean(Matrix{Float64}(matrix), dims=2)),
    vec(var(Matrix{Float64}(matrix), dims=2; corrected=false)))

_row_mean_variance(matrix::SparseMatrixCSC{<:Real,Int}) = _sparse_row_mean_variance(matrix)
_row_mean_variance(matrix::AbstractMatrix{<:Real}) = _dense_row_mean_variance(matrix)

"""
    normalize_counts(experiment; scale_factor=1e4, log_transform=true)

Library-size normalize a single-cell count matrix and optionally apply
`log1p` stabilization.
"""
function normalize_counts(experiment::SingleCellExperiment; scale_factor::Real=1e4, log_transform::Bool=true, use_cuda::Bool=false)
    _ctx = active_provenance_context()
    counts = experiment.counts
    if use_cuda
        CUDA.functional() || throw(ArgumentError("CUDA is not available"))
        gpu_counts = counts isa CUDA.CuArray ? counts : CUDA.CuArray{Float32}(Array(counts))
        library_sizes = max.(vec(sum(gpu_counts, dims=1)), eps(Float32))
        normalized = gpu_counts ./ reshape(library_sizes, 1, :) .* Float32(scale_factor)
        log_transform && (normalized = log1p.(normalized))
        return _register_singlecell_result!(_ctx, normalized, "normalize_counts";
            parameters=(scale_factor=Float64(scale_factor), log_transform=Bool(log_transform), use_cuda=true, n_genes=size(counts,1), n_cells=size(counts,2)))
    end

    library_sizes = vec(sum(counts, dims=1))
    library_sizes = max.(Float64.(library_sizes), eps(Float64))
    rows, cols, values = _matrix_triplets(counts)
    normalized_values = values ./ library_sizes[cols] .* Float64(scale_factor)
    log_transform && (normalized_values = log1p.(normalized_values))
    result = sparse(rows, cols, normalized_values, size(counts, 1), size(counts, 2))

    return _register_singlecell_result!(_ctx, result, "normalize_counts";
        parameters=(scale_factor=Float64(scale_factor), log_transform=Bool(log_transform), use_cuda=false, n_genes=size(counts,1), n_cells=size(counts,2)))
end

function _projection_feature_ids(experiment::SingleCellExperiment; use_variable_features::Bool=true)
    if use_variable_features
        feature_indices = _variable_feature_indices(experiment)
        isempty(feature_indices) || return [experiment.gene_ids[index] for index in feature_indices]
    end
    return copy(experiment.gene_ids)
end

function _normalized_projection_matrix(experiment::SingleCellExperiment, feature_ids::AbstractVector{<:String}; scale_factor::Real=1e4, log_transform::Bool=true)
    counts = experiment.counts isa CUDA.CuArray ? Array(experiment.counts) : experiment.counts
    library_sizes = vec(sum(counts, dims=1))
    library_sizes = max.(Float64.(library_sizes), eps(Float64))
    feature_lookup = Dict(feature => index for (index, feature) in pairs(feature_ids))
    rows, cols, values = _matrix_triplets(counts)
    selected_rows = Int[]
    selected_cols = Int[]
    selected_values = Float64[]
    for index in eachindex(values)
        feature_index = get(feature_lookup, experiment.gene_ids[rows[index]], 0)
        feature_index == 0 && continue
        value = values[index] / library_sizes[cols[index]] * Float64(scale_factor)
        value = log_transform ? log1p(value) : value
        push!(selected_rows, feature_index)
        push!(selected_cols, cols[index])
        push!(selected_values, value)
    end
    return sparse(selected_rows, selected_cols, selected_values, length(feature_ids), size(counts, 2))
end

function _fit_projection_basis(normalized::AbstractMatrix{<:Real}; n_components::Integer=20, center_features::Bool=true, scale_features::Bool=true)
    means, variances = _row_mean_variance(normalized)
    feature_means = center_features ? means : zeros(Float64, length(means))
    feature_scales = scale_features ? sqrt.(max.(variances, eps(Float64))) : ones(Float64, length(means))
    centered = Matrix{Float64}(permutedims(normalized))
    center_features && (centered .-= reshape(feature_means, 1, :))
    scale_features && (centered ./= reshape(feature_scales, 1, :))
    svd_result = svd(centered; full=false)
    components = min(Int(n_components), size(svd_result.U, 2), size(svd_result.Vt, 1), length(svd_result.S))
    embedding = Array(svd_result.U[:, 1:components] * Diagonal(svd_result.S[1:components]))
    loadings = Matrix(transpose(svd_result.Vt[1:components, :]))
    return embedding, feature_means, feature_scales, loadings
end

"""
    fit_singlecell_projection_model(experiment; n_components=20, use_variable_features=true, scale_factor=1e4, log_transform=true, center_features=true, scale_features=true)

Fit a reusable projection model from a reference experiment. The fitted model
stores the feature order together with normalization, centering, scaling, and
PCA loadings so query data can be projected without rebuilding the basis.
"""
function fit_singlecell_projection_model(experiment::SingleCellExperiment; n_components::Integer=20, use_variable_features::Bool=true, scale_factor::Real=1e4, log_transform::Bool=true, center_features::Bool=true, scale_features::Bool=true)
    feature_ids = _projection_feature_ids(experiment; use_variable_features=use_variable_features)
    isempty(feature_ids) && throw(ArgumentError("No features available to fit a projection model"))
    normalized = _normalized_projection_matrix(experiment, feature_ids; scale_factor=scale_factor, log_transform=log_transform)
    _, feature_means, feature_scales, loadings = _fit_projection_basis(normalized; n_components=n_components, center_features=center_features, scale_features=scale_features)
    model = SingleCellProjectionModel(String.(feature_ids), Float64(scale_factor), log_transform, center_features, scale_features, feature_means, feature_scales, loadings)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, model, "fit_singlecell_projection_model")
end

"""
    project_singlecell(experiment, model; reduction_name="pca", n_components=nothing)

Project a query experiment onto a fitted single-cell projection model using the
stored feature alignment and PCA basis.
"""
function project_singlecell(experiment::SingleCellExperiment, model::SingleCellProjectionModel; reduction_name::String="pca", n_components::Union{Nothing,Integer}=nothing)
    query_genes = Set(experiment.gene_ids)
    shared_features = count(feature -> feature in query_genes, model.feature_ids)
    shared_features > 0 || throw(ArgumentError("No genes overlap between the query experiment and the projection model"))

    normalized = _normalized_projection_matrix(experiment, model.feature_ids; scale_factor=model.scale_factor, log_transform=model.log_transform)
    component_count = n_components === nothing ? size(model.pca_loadings, 2) : min(Int(n_components), size(model.pca_loadings, 2))
    loadings = model.pca_loadings[:, 1:component_count]
    weights = model.scale_features ? loadings ./ reshape(model.feature_scales, :, 1) : loadings
    offset = model.center_features ? model.feature_means' * weights : zeros(Float64, component_count)
    embedding = permutedims(normalized) * weights
    model.center_features && (embedding .-= reshape(offset, 1, :))
    experiment.reductions[reduction_name] = Array{Float64}(embedding)
    update_provenance!(experiment.metadata; source="SingleCell/project_singlecell", notes=["projected onto fitted single-cell basis"], parameters=(reduction_name=reduction_name, n_components=component_count, shared_features=shared_features))
    return experiment.reductions[reduction_name]
end

"""
    sctransform(experiment; min_total=10)

Apply a variance-stabilizing transform to the filtered count matrix and return
the transformed matrix together with the retained gene and cell identifiers.
"""
function sctransform(experiment::SingleCellExperiment; min_total::Integer=10)
    cm = count_matrix(experiment)
    filtered = filter_low_counts(cm; min_total=min_total)
    norm_factors = calc_norm_factors(filtered)
    dispersions = estimate_dispersions(filtered, norm_factors)
    transformed = vst(filtered, dispersions; norm_factors=norm_factors)
    sc_result = (matrix=transformed, gene_ids=filtered.gene_ids, sample_ids=filtered.sample_ids, provenance=provenance_record("SingleCellTransform", "SingleCell/sctransform"; parameters=(min_total=Int(min_total), output_genes=length(filtered.gene_ids), output_cells=length(filtered.sample_ids))))
    _ctx = active_provenance_context()
    _register_singlecell_result!(_ctx, sc_result, "sctransform";
        parameters=(min_total=Int(min_total), n_genes=length(filtered.gene_ids), n_cells=length(filtered.sample_ids)))

    return sc_result
end

function _variable_feature_indices(experiment::SingleCellExperiment)
    isempty(experiment.variable_features) && return Int[]
    feature_set = if haskey(experiment.variable_features, "variance")
        experiment.variable_features["variance"]
    elseif haskey(experiment.variable_features, "dispersion")
        experiment.variable_features["dispersion"]
    elseif haskey(experiment.variable_features, "vst")
        experiment.variable_features["vst"]
    else
        first(values(experiment.variable_features))
    end
    gene_index = Dict(gene => index for (index, gene) in pairs(experiment.gene_ids))
    return [gene_index[feature] for feature in feature_set if haskey(gene_index, feature)]
end

"""
    find_variable_features(experiment; n_features=2000, method=:vst, min_mean=0.1, max_mean=8.0)

    Rank genes by variance or dispersion and store the selected feature set on the
    experiment for later PCA or clustering steps. The sparse counts are screened
    without densifying the matrix, which keeps large single-cell datasets memory-safe.

    `:vst` remains available as a sparse-safe heuristic on the log-normalized matrix.
    Use `sctransform` directly when you want the full variance-stabilized matrix.
"""
function find_variable_features(experiment::SingleCellExperiment; n_features::Integer=2000, method::Symbol=:vst, min_mean::Real=0.1, max_mean::Real=8.0)
    n_features >= 0 || throw(ArgumentError("n_features must be nonnegative"))
    method in (:vst, :variance, :dispersion) || throw(ArgumentError("method must be :vst, :variance, or :dispersion"))

    matrix = normalize_counts(experiment)
    means, variances = _row_mean_variance(matrix)
    scores = if method == :variance
        variances
    elseif method == :dispersion
        variances ./ (means .+ eps(Float64))
    else
        variances ./ sqrt.(means .+ eps(Float64))
    end

    if length(means) == 0
        experiment.variable_features[string(method)] = String[]
        return String[]
    end

    keep = [index for index in eachindex(means) if isfinite(means[index]) && isfinite(scores[index]) && means[index] >= min_mean && means[index] <= max_mean]
    isempty(keep) && (experiment.variable_features[string(method)] = String[]; return String[])

    limit = min(Int(n_features), length(keep))
    selected_indices = sort(keep; by = index -> scores[index], rev=true)[1:limit]
    selected = String[experiment.gene_ids[index] for index in selected_indices]
    experiment.variable_features[string(method)] = selected
    update_provenance!(experiment.metadata; source="SingleCell/find_variable_features", notes=["selected variable features"], parameters=(method=method, n_features=Int(n_features), selected_count=length(selected), min_mean=Float64(min_mean), max_mean=Float64(max_mean)))
    _ctx = active_provenance_context()
    _register_singlecell_result!(_ctx, selected, "find_variable_features";
        parameters=(method=method, n_features=Int(n_features), selected_count=length(selected)))

    return selected
end

function _center_rows(matrix::AbstractMatrix{<:Real})
    dense = Matrix{Float64}(matrix)
    dense .-= mean(dense, dims=1)
    return dense
end

"""
    run_pca(experiment; normalized=nothing, n_components=20, use_variable_features=true, projection_model=nothing)

Compute a PCA embedding for cells and cache it on the experiment. When variable
features have been selected, PCA is computed on that subset by default.
"""
function run_pca(experiment::SingleCellExperiment; normalized::Union{Nothing,AbstractMatrix}=nothing, n_components::Integer=20, use_variable_features::Bool=true, use_cuda::Bool=false, projection_model::Union{Nothing,SingleCellProjectionModel}=nothing)
    if projection_model !== nothing
        return project_singlecell(experiment, projection_model; n_components=n_components)
    end
    matrix = normalized === nothing ? normalize_counts(experiment; use_cuda=use_cuda) : normalized
    if normalized === nothing && use_variable_features
        feature_indices = _variable_feature_indices(experiment)
        isempty(feature_indices) || (matrix = matrix[feature_indices, :])
    end
    centered = if use_cuda
        CUDA.CuArray{Float32}(Array(permutedims(matrix)))
    else
        Matrix{Float64}(permutedims(matrix))
    end
    centered .-= mean(centered, dims=1)
    svd_result = svd(centered; full=false)
    components = min(n_components, size(svd_result.U, 2))
    embedding = Array(svd_result.U[:, 1:components] * Diagonal(svd_result.S[1:components]))
    experiment.reductions["pca"] = embedding
    update_provenance!(experiment.metadata; source="SingleCell/run_pca", notes=["computed PCA embedding"], parameters=(n_components=components, use_variable_features=Bool(use_variable_features), use_cuda=Bool(use_cuda), used_projection_model=false))
    _ctx = active_provenance_context()
    _register_singlecell_result!(_ctx, embedding, "run_pca";
        parameters=(n_components=components, use_variable_features=Bool(use_variable_features), use_cuda=Bool(use_cuda)))

    return embedding
end

function _pairwise_distance(data::Matrix{Float64}, left::Int, right::Int)
    return norm(view(data, left, :) .- view(data, right, :))
end

function _knn_neighbors(data::Matrix{Float64}, k::Int)
    n = size(data, 1)
    neighbors = Vector{Vector{Int}}(undef, n)
    for index in 1:n
        distances = [(other, _pairwise_distance(data, index, other)) for other in 1:n if other != index]
        sort!(distances; by = last)
        neighbors[index] = [neighbor for (neighbor, _) in Iterators.take(distances, min(k, length(distances)))]
    end
    return neighbors
end

"""
    find_neighbors(experiment; embedding=nothing, reduction="pca", k=15, graph_name="neighbors")

Build and cache a k-nearest-neighbor graph for cells. The graph is stored on the
experiment so downstream clustering can reuse it.
"""
function find_neighbors(experiment::SingleCellExperiment; embedding::Union{Nothing,AbstractMatrix}=nothing, reduction::String="pca", k::Int=15, graph_name::String="neighbors", use_cuda::Bool=false)
    k >= 0 || throw(ArgumentError("k must be nonnegative"))
    data = if embedding === nothing
        haskey(experiment.reductions, reduction) ? experiment.reductions[reduction] : run_pca(experiment; n_components=max(2, k), use_cuda=use_cuda)
    else
        Matrix{Float64}(embedding)
    end
    neighbors = _knn_neighbors(Matrix{Float64}(data), k)
    experiment.neighbors[graph_name] = neighbors
    update_provenance!(experiment.metadata; source="SingleCell/find_neighbors", notes=["computed k-nearest-neighbor graph"], parameters=(graph_name=graph_name, k=Int(k), reduction=reduction, use_cuda=Bool(use_cuda)))
    _ctx = active_provenance_context()
    _register_singlecell_result!(_ctx, neighbors, "find_neighbors";
        parameters=(k=Int(k), reduction=reduction, n_cells=length(neighbors)))

    return neighbors
end

"""
    SpatialMoranResult

Summary statistics from a Moran's I spatial autocorrelation test.
"""
struct SpatialMoranResult <: AbstractAnalysisResult
    gene::String
    moran_i::Float64
    expected_i::Float64
    p_value::Float64
    graph_name::String
    provenance::ResultProvenance
end

SpatialMoranResult(moran_i, expected_i, p_value, graph_name) =
    SpatialMoranResult(moran_i, expected_i, p_value, graph_name, provenance_record("SpatialMoranResult", "singlecell"))

SpatialMoranResult(gene, moran_i, expected_i, p_value, graph_name) =
    SpatialMoranResult(
        String(gene),
        Float64(moran_i),
        Float64(expected_i),
        Float64(p_value),
        String(graph_name),
        provenance_record(
            "SpatialMoranResult",
            "SingleCell/SpatialMoranResult";
            notes=["constructed spatial Moran's I summary"],
            parameters=(graph_name=String(graph_name))))

function find_spatial_neighbors(experiment::SingleCellExperiment; spatial_coords::Union{Nothing,AbstractMatrix}=nothing, k::Int=6, graph_name::String="spatial_neighbors")
    coords = spatial_coords === nothing ? experiment.spatial_coords : _normalize_spatial_coords(spatial_coords, length(experiment.cell_ids))
    coords === nothing && throw(ArgumentError("spatial coordinates are required to build a spatial neighbor graph"))
    neighbors = _knn_neighbors(Matrix{Float64}(coords), k)
    experiment.neighbors[graph_name] = neighbors
    update_provenance!(experiment.metadata; source="SingleCell/find_spatial_neighbors", notes=["computed spatial neighbor graph"], parameters=(graph_name=graph_name, k=Int(k)))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, neighbors, "find_spatial_neighbors")
end

function _spatial_neighbors(experiment::SingleCellExperiment; spatial_coords::Union{Nothing,AbstractMatrix}=nothing, neighbors::Union{Nothing,Vector{Vector{Int}}}=nothing, graph_name::String="spatial_neighbors", k::Int=6)
    if neighbors !== nothing
        return neighbors
    end
    haskey(experiment.neighbors, graph_name) && return experiment.neighbors[graph_name]
    return find_spatial_neighbors(experiment; spatial_coords=spatial_coords, k=k, graph_name=graph_name)
end

function _moran_i_statistic(values::AbstractVector{<:Real}, neighbors::Vector{Vector{Int}})
    n = length(values)
    n == length(neighbors) || throw(DimensionMismatch("values and neighbors must have the same length"))
    values_vector = Float64.(values)
    centered = values_vector .- mean(values_vector)
    denominator = sum(abs2, centered)
    denominator <= eps(Float64) && return 0.0
    numerator = 0.0
    total_weight = 0.0
    for source in 1:n
        for target in neighbors[source]
            source == target && continue
            weight = 1.0
            numerator += weight * centered[source] * centered[target]
            total_weight += weight
        end
    end
    total_weight <= 0 && return 0.0
    return (n / total_weight) * (numerator / denominator)
end

function _moran_i_pvalue(values::AbstractVector{<:Real}, neighbors::Vector{Vector{Int}}, observed::Float64; n_permutations::Int=0, rng::AbstractRNG=Random.default_rng())
    n_permutations <= 0 && return NaN
    extreme = 1
    absolute_observed = abs(observed)
    values_vector = Float64.(values)
    for _ in 1:n_permutations
        permuted = values_vector[randperm(rng, length(values_vector))]
        permuted_i = _moran_i_statistic(permuted, neighbors)
        abs(permuted_i) >= absolute_observed && (extreme += 1)
    end
    return extreme / (n_permutations + 1)
end

function moran_i_test(experiment::SingleCellExperiment, gene::String; spatial_coords::Union{Nothing,AbstractMatrix}=nothing, neighbors::Union{Nothing,Vector{Vector{Int}}}=nothing, graph_name::String="spatial_neighbors", k::Int=6, normalize::Bool=true, n_permutations::Int=0, rng::AbstractRNG=Random.default_rng())
    coords_neighbors = _spatial_neighbors(experiment; spatial_coords=spatial_coords, neighbors=neighbors, graph_name=graph_name, k=k)
    expression = normalize ? normalize_counts(experiment) : Array{Float64}(experiment.counts)
    gene_lookup = Dict(gene_name => index for (index, gene_name) in enumerate(experiment.gene_ids))
    haskey(gene_lookup, String(gene)) || throw(ArgumentError("gene '$gene' was not found"))
    values = Vector{Float64}(Array(expression[gene_lookup[String(gene)], :]))
    moran_i = _moran_i_statistic(values, coords_neighbors)
    expected_i = length(values) <= 1 ? 0.0 : -1 / (length(values) - 1)
    result = SpatialMoranResult(
        String(gene),
        moran_i,
        expected_i,
        _moran_i_pvalue(values, coords_neighbors, moran_i; n_permutations=n_permutations, rng=rng),
        graph_name,
        provenance_record(
            "SpatialMoranResult",
            "SingleCell/moran_i_test";
            notes=["computed Moran's I for a single gene"],
            parameters=(gene=String(gene), graph_name=graph_name, k=Int(k), normalize=Bool(normalize), n_permutations=Int(n_permutations))))
    update_provenance!(experiment.metadata; source="SingleCell/moran_i_test", notes=["computed Moran's I for a spatial feature"], parameters=(gene=String(gene), graph_name=graph_name, k=Int(k), normalize=Bool(normalize), n_permutations=Int(n_permutations)))
    return result
end

function find_spatially_variable_genes(experiment::SingleCellExperiment; spatial_coords::Union{Nothing,AbstractMatrix}=nothing, neighbors::Union{Nothing,Vector{Vector{Int}}}=nothing, graph_name::String="spatial_neighbors", k::Int=6, normalize::Bool=true, top_n::Int=50, min_expression::Real=0.0, n_permutations::Int=0, rng::AbstractRNG=Random.default_rng())
    coords_neighbors = _spatial_neighbors(experiment; spatial_coords=spatial_coords, neighbors=neighbors, graph_name=graph_name, k=k)
    expression = normalize ? normalize_counts(experiment) : Array{Float64}(experiment.counts)
    results = SpatialMoranResult[]
    for (gene_index, gene_name) in enumerate(experiment.gene_ids)
        values = Vector{Float64}(Array(expression[gene_index, :]))
        mean(values) < min_expression && continue
        moran_i = _moran_i_statistic(values, coords_neighbors)
        expected_i = length(values) <= 1 ? 0.0 : -1 / (length(values) - 1)
        push!(results, SpatialMoranResult(
            gene_name,
            moran_i,
            expected_i,
            _moran_i_pvalue(values, coords_neighbors, moran_i; n_permutations=n_permutations, rng=rng),
            graph_name,
            provenance_record(
                "SpatialMoranResult",
                "SingleCell/find_spatially_variable_genes";
                notes=["ranked gene by spatial Moran's I"],
                parameters=(gene=gene_name, graph_name=graph_name, k=Int(k), normalize=Bool(normalize), min_expression=Float64(min_expression), n_permutations=Int(n_permutations)))))
    end
    sort!(results; by = result -> result.moran_i, rev=true)
    selected = results[1:min(top_n, length(results))]
    update_provenance!(experiment.metadata; source="SingleCell/find_spatially_variable_genes", notes=["ranked spatially variable genes by Moran's I"], parameters=(graph_name=graph_name, k=Int(k), normalize=Bool(normalize), top_n=Int(top_n), min_expression=Float64(min_expression), n_permutations=Int(n_permutations), result_count=length(selected)))
    return selected
end

function _snn_graph(neighbors::Vector{Vector{Int}}; shared_threshold::Int=1)
    n = length(neighbors)
    adjacency = [Dict{Int,Float64}() for _ in 1:n]
    neighbor_sets = [Set(neighbor_list) for neighbor_list in neighbors]

    for left in 1:n-1
        left_set = neighbor_sets[left]
        for right in neighbors[left]
            right <= left && continue
            shared = count(in(left_set), neighbors[right])
            shared >= shared_threshold || continue
            weight = Float64(shared)
            adjacency[left][right] = weight
            adjacency[right][left] = weight
        end
    end

    return adjacency
end

function _community_totals(adjacency::Vector{Dict{Int,Float64}}, labels::Vector{Int})
    totals = Dict{Int,Float64}()
    for (node, label) in enumerate(labels)
        totals[label] = get(totals, label, 0.0) + sum(values(adjacency[node]))
    end
    return totals
end

function _local_move_communities(adjacency::Vector{Dict{Int,Float64}}, labels::Vector{Int}; rng::AbstractRNG=Random.default_rng(), max_passes::Int=20)
    n = length(labels)
    n == 0 && return labels
    degree = [sum(values(node)) for node in adjacency]
    total_weight = sum(degree) / 2
    total_weight <= 0 && return labels

    nodes = collect(1:n)
    for _ in 1:max_passes
        moved = false
        shuffle!(rng, nodes)
        totals = _community_totals(adjacency, labels)

        for node in nodes
            current = labels[node]
            node_degree = degree[node]
            node_degree == 0 && continue

            totals[current] = get(totals, current, 0.0) - node_degree
            totals[current] <= 0 && delete!(totals, current)

            weights_to_community = Dict{Int,Float64}()
            for (neighbor, weight) in adjacency[node]
                community = labels[neighbor]
                weights_to_community[community] = get(weights_to_community, community, 0.0) + weight
            end
            weights_to_community[current] = get(weights_to_community, current, 0.0)

            best_community = current
            best_score = weights_to_community[current] - get(totals, current, 0.0) * node_degree / (2 * total_weight)
            for (community, in_weight) in weights_to_community
                score = in_weight - get(totals, community, 0.0) * node_degree / (2 * total_weight)
                if score > best_score + 1e-12
                    best_score = score
                    best_community = community
                end
            end

            labels[node] = best_community
            totals[best_community] = get(totals, best_community, 0.0) + node_degree
            moved |= best_community != current
        end

        moved || break
    end

    return labels
end

function _refine_communities(adjacency::Vector{Dict{Int,Float64}}, labels::Vector{Int})
    n = length(labels)
    n == 0 && return labels

    refined = fill(0, n)
    next_label = 0
    grouped = Dict{Int,Vector{Int}}()
    for (node, label) in enumerate(labels)
        push!(get!(grouped, label, Int[]), node)
    end

    for label in sort(collect(keys(grouped)))
        members = grouped[label]
        member_set = Set(members)
        seen = Set{Int}()
        for seed in members
            seed in seen && continue
            next_label += 1
            stack = [seed]
            push!(seen, seed)
            refined[seed] = next_label

            while !isempty(stack)
                node = pop!(stack)
                for neighbor in keys(adjacency[node])
                    neighbor in member_set || continue
                    neighbor in seen && continue
                    push!(seen, neighbor)
                    refined[neighbor] = next_label
                    push!(stack, neighbor)
                end
            end
        end
    end

    return refined
end

function _cluster_neighbors(neighbors::Vector{Vector{Int}}; method::Symbol=:leiden, shared_threshold::Int=1, random_seed::Int=1)
    adjacency = _snn_graph(neighbors; shared_threshold=shared_threshold)
    labels = collect(1:length(neighbors))
    rng = MersenneTwister(random_seed)
    for _ in 1:10
        previous = copy(labels)
        _local_move_communities(adjacency, labels; rng=rng)
        if method == :leiden
            labels = _refine_communities(adjacency, labels)
            _local_move_communities(adjacency, labels; rng=rng)
        end
        labels == previous && break
    end

    remap = Dict{Int,Int}()
    next_label = 0
    clustered = similar(labels)
    for (index, label) in enumerate(labels)
        community = get!(remap, label) do
            next_label += 1
            next_label
        end
        clustered[index] = community
    end
    return clustered
end

function _mutual_knn_graph(data::Matrix{Float64}, k::Int)
    neighbors = _knn_neighbors(data, k)
    n = length(neighbors)
    adjacency = [Int[] for _ in 1:n]
    for source in 1:n
        for target in neighbors[source]
            source in neighbors[target] || continue
            push!(adjacency[source], target)
        end
    end
    for node in adjacency
        sort!(unique!(node))
    end
    return adjacency
end

function _graph_components(adjacency::Vector{Vector{Int}})
    n = length(adjacency)
    labels = fill(0, n)
    current = 0
    for seed in 1:n
        labels[seed] != 0 && continue
        current += 1
        stack = [seed]
        labels[seed] = current
        while !isempty(stack)
            node = pop!(stack)
            for neighbor in adjacency[node]
                labels[neighbor] != 0 && continue
                labels[neighbor] = current
                push!(stack, neighbor)
            end
        end
    end
    return labels
end

function _connected_components(neighbors::Vector{Vector{Int}}; shared_threshold::Int=1)
    n = length(neighbors)
    labels = fill(0, n)
    current = 0
    for seed in 1:n
        labels[seed] != 0 && continue
        current += 1
        stack = [seed]
        labels[seed] = current
        while !isempty(stack)
            node = pop!(stack)
            for candidate in 1:n
                labels[candidate] == 0 || continue
                shared = length(intersect(neighbors[node], neighbors[candidate]))
                mutual = candidate in neighbors[node] && node in neighbors[candidate]
                if mutual || shared >= shared_threshold
                    labels[candidate] = current
                    push!(stack, candidate)
                end
            end
        end
    end
    return labels
end

function _kmeans(data::Matrix{Float64}, k::Int; max_iter::Int=100, random_seed::Int=1)
    n = size(data, 1)
    n == 0 && return Int[]
    k = clamp(k, 1, n)
    rng = MersenneTwister(random_seed)
    centers = data[rand(rng, 1:n, k), :]
    labels = fill(1, n)
    for _ in 1:max_iter
        changed = false
        for index in 1:n
            distances = [norm(view(data, index, :) .- view(centers, center, :)) for center in 1:k]
            label = argmin(distances)
            if label != labels[index]
                labels[index] = label
                changed = true
            end
        end
        new_centers = similar(centers)
        for center in 1:k
            members = findall(==(center), labels)
            if isempty(members)
                new_centers[center, :] .= vec(data[rand(rng, 1:n), :])
            else
                new_centers[center, :] .= vec(mean(data[members, :], dims=1))
            end
        end
        centers .= new_centers
        changed || break
    end
    return labels
end

"""
    cluster_cells(experiment; ...)

Cluster cells from a PCA embedding using kNN, graph-connected-components, or
k-means-style assignments. This is the core clustering primitive used by the
Seurat-style wrapper `find_clusters`.
"""
function cluster_cells(experiment::SingleCellExperiment; embedding::Union{Nothing,AbstractMatrix}=nothing, method::Symbol=:knn, k::Int=15, shared_threshold::Int=1, n_clusters::Int=8, random_seed::Int=1)
    data = if embedding === nothing
        haskey(experiment.reductions, "pca") ? experiment.reductions["pca"] : run_pca(experiment; n_components=max(n_clusters, 2))
    else
        Matrix{Float64}(embedding)
    end
    method == :knn || method == :kmeans || method == :graph || throw(ArgumentError("method must be :knn, :graph, or :kmeans"))
    labels = if method == :kmeans
        _kmeans(data, n_clusters; random_seed=random_seed)
    elseif method == :graph
        _graph_components(_mutual_knn_graph(data, k))
    else
        _connected_components(_knn_neighbors(data, k); shared_threshold=shared_threshold)
    end
    experiment.clusters[method == :kmeans ? "kmeans" : method == :graph ? "graph" : "knn"] = labels
    update_provenance!(experiment.metadata; source="SingleCell/cluster_cells", notes=["clustered cells from an embedding"], parameters=(method=method, k=Int(k), shared_threshold=Int(shared_threshold), n_clusters=Int(n_clusters), random_seed=Int(random_seed), n_cells=length(labels)))
    _ctx = active_provenance_context()
    _register_singlecell_result!(_ctx, labels, "cluster_cells";
        parameters=(method=method, k=Int(k), n_cells=length(labels), n_clusters_found=length(unique(labels))))

    return labels
end

"""
    find_clusters(experiment; ...)

    Seurat-style clustering wrapper that reuses a cached neighbor graph when
    available and otherwise builds one from the selected reduction.

    `method=:leiden` is the default graph-based community detection path.
    `method=:louvain` uses the same SNN graph but skips the Leiden refinement step.
    `method=:snn` is kept as a compatibility alias for `:leiden`.
"""
function find_clusters(experiment::SingleCellExperiment; neighbors::Union{Nothing,Vector{Vector{Int}}}=nothing, embedding::Union{Nothing,AbstractMatrix}=nothing, reduction::String="pca", graph_name::String="neighbors", method::Symbol=:leiden, k::Int=15, shared_threshold::Int=1, n_clusters::Int=8, random_seed::Int=1)
    method in (:snn, :leiden, :louvain, :graph, :knn, :kmeans) || throw(ArgumentError("method must be :snn, :leiden, :louvain, :graph, :knn, or :kmeans"))
    if method == :kmeans
        return cluster_cells(experiment; embedding=embedding, method=:kmeans, n_clusters=n_clusters, random_seed=random_seed)
    end

    graph = if neighbors !== nothing
        neighbors
    elseif haskey(experiment.neighbors, graph_name)
        experiment.neighbors[graph_name]
    else
        find_neighbors(experiment; embedding=embedding, reduction=reduction, k=k, graph_name=graph_name)
    end
    labels = if method in (:snn, :leiden)
        _cluster_neighbors(graph; method=:leiden, shared_threshold=shared_threshold, random_seed=random_seed)
    elseif method == :louvain
        _cluster_neighbors(graph; method=:louvain, shared_threshold=shared_threshold, random_seed=random_seed)
    else
        _connected_components(graph; shared_threshold=shared_threshold)
    end
    experiment.clusters[graph_name] = labels
    update_provenance!(experiment.metadata; source="SingleCell/find_clusters", notes=["computed cluster labels from a neighbor graph"], parameters=(graph_name=graph_name, method=method, k=Int(k), shared_threshold=Int(shared_threshold), n_clusters=Int(n_clusters), random_seed=Int(random_seed), reused_neighbors=neighbors === nothing && haskey(experiment.neighbors, graph_name)))
    _ctx = active_provenance_context()
    _register_singlecell_result!(_ctx, labels, "find_clusters";
        parameters=(method=method, k=Int(k), n_cells=length(labels), n_clusters_found=length(unique(labels))))

    return labels
end

"""
    run_umap(experiment; embedding=nothing, n_neighbors=15, min_dist=0.3)

Project the data to two dimensions using UMAP when an implementation is
available, otherwise fall back to the first two PCA dimensions.
"""
function run_umap(experiment::SingleCellExperiment; embedding::Union{Nothing,AbstractMatrix}=nothing, n_neighbors::Int=15, min_dist::Real=0.3)
    matrix = if embedding === nothing
        haskey(experiment.reductions, "pca") ? experiment.reductions["pca"] : run_pca(experiment; n_components=10)
    else
        Matrix{Float64}(embedding)
    end
    fallback_notes = String[]
    if isdefined(Main, :UMAP)
        try
            umap_result = getfield(Main, :UMAP).umap(matrix; n_neighbors=n_neighbors, min_dist=min_dist)
            experiment.reductions["umap"] = Matrix{Float64}(umap_result)
            update_provenance!(experiment.metadata; source="SingleCell/run_umap", notes=["computed UMAP embedding"], parameters=(n_neighbors=Int(n_neighbors), min_dist=Float64(min_dist), backend="UMAP.jl"))
            return experiment.reductions["umap"]
        catch err
            @warn "UMAP failed; falling back to the first PCA dimensions" exception=(err, catch_backtrace())
            push!(fallback_notes, "UMAP.jl call failed; used first PCA dimensions instead")
        end
    else
        push!(fallback_notes, "UMAP.jl was unavailable; used first PCA dimensions instead")
    end
    fallback = matrix[:, 1:min(2, size(matrix, 2))]
    experiment.reductions["umap"] = Matrix{Float64}(fallback)
    update_provenance!(experiment.metadata; source="SingleCell/run_umap", status=:warning, notes=vcat(["stored fallback low-dimensional embedding"], fallback_notes), fallbacks=fallback_notes, parameters=(n_neighbors=Int(n_neighbors), min_dist=Float64(min_dist), backend="pca_fallback"))
    _ctx = active_provenance_context()
    _register_singlecell_result!(_ctx, experiment.reductions["umap"], "run_umap";
        parameters=(n_neighbors=Int(n_neighbors), min_dist=Float64(min_dist), backend="pca_fallback", n_cells=size(experiment.reductions["umap"],1)))

    return experiment.reductions["umap"]
end

"""
    summarize_clusters(labels)

Count the number of cells assigned to each cluster label.
"""
function summarize_clusters(labels::AbstractVector{<:Integer})
    counts = Dict{Int,Int}()
    for label in labels
        counts[Int(label)] = get(counts, Int(label), 0) + 1
    end
    return counts
end

"""
    cluster_marker_summary(experiment, labels; min_total=10, shrink=true, top_n=5)

Summarize the top differential-expression markers for each cluster.
"""
function cluster_marker_summary(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; min_total::Integer=10, shrink::Bool=true, top_n::Int=5)
    unique_labels = sort!(unique(Int.(labels)))
    summaries = Dict{Int,Vector{DEResult}}()
    for cluster_id in unique_labels
        markers = find_cluster_markers(experiment, labels; cluster_id=cluster_id, min_total=min_total, shrink=shrink)
        summaries[cluster_id] = first(sort(markers; by = result -> (result.padj, -abs(result.log2_fold_change))), min(top_n, length(markers)))
    end
    update_provenance!(experiment.metadata; source="SingleCell/cluster_marker_summary", notes=["summarized per-cluster marker genes"], parameters=(cluster_count=length(unique_labels), min_total=Int(min_total), shrink=Bool(shrink), top_n=Int(top_n)))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, summaries, "cluster_marker_summary")
end

"""
    find_cluster_markers(experiment, labels; cluster_id, ...)

Compare one cluster against all other cells and return differential-expression
results for the corresponding marker genes.
"""
function find_cluster_markers(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; cluster_id::Integer, min_total::Integer=10, shrink::Bool=true)
    length(labels) == length(experiment.cell_ids) || throw(ArgumentError("labels must match the number of cells"))
    count(==(cluster_id), labels) == 0 && return DEResult[]
    count(!=(cluster_id), labels) == 0 && return DEResult[]
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, find_markers(experiment, labels; ident_1=cluster_id, test=:wilcox, min_total=min_total, shrink=shrink), "find_cluster_markers")
end

"""
    find_markers(experiment, labels; ident_1, ident_2=nothing, ...)

Seurat-style marker wrapper. By default this runs a Wilcoxon rank-sum test on
log-normalized counts, matching Seurat's default marker statistic.

Set `test=:deseq2` or `pseudobulk=true` to route the analysis through the
DifferentialExpression module. Pseudobulk mode aggregates cells by
`sample_labels` before calling `differential_expression`.
"""
function _wilcoxon_rank_sum(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n1 = length(x)
    n2 = length(y)
    if n1 == 0 || n2 == 0
        return (stat=0.0, pvalue=1.0)
    end

    combined = vcat(Float64.(x), Float64.(y))
    n = length(combined)
    order = sortperm(combined)
    ranks = zeros(Float64, n)

    tie_correction = 0.0
    index = 1
    while index <= n
        stop = index
        while stop < n && combined[order[stop + 1]] == combined[order[index]]
            stop += 1
        end
        average_rank = (index + stop) / 2
        for position in index:stop
            ranks[order[position]] = average_rank
        end
        run = stop - index + 1
        run > 1 && (tie_correction += run^3 - run)
        index = stop + 1
    end

    rank_sum = sum(ranks[1:n1])
    u_stat = rank_sum - n1 * (n1 + 1) / 2
    mean_u = n1 * n2 / 2
    variance = n1 * n2 / 12 * ((n + 1) - tie_correction / (n * (n - 1)))
    if variance <= 0
        return (stat=0.0, pvalue=1.0)
    end

    z = (abs(u_stat - mean_u) - 0.5) / sqrt(variance)
    z = max(z, 0.0)
    pvalue = erfc(z / sqrt(2))
    return (stat=(u_stat - mean_u) / sqrt(variance), pvalue=clamp(pvalue, 0.0, 1.0))
end

function _expression_for_cells(matrix::AbstractMatrix{<:Real}, row::Int, columns::AbstractVector{Int})
    values = Vector{Float64}(undef, length(columns))
    @inbounds for (position, column) in enumerate(columns)
        values[position] = Float64(matrix[row, column])
    end
    return values
end

function _expression_for_cells(matrix::CUDA.CuArray{<:Real,2}, row::Int, columns::AbstractVector{Int})
    return _expression_for_cells(Array(matrix), row, columns)
end

function _pseudobulk_count_matrix(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}, ident_1::Integer, ident_2::Integer; sample_labels::AbstractVector{<:String})
    length(labels) == length(sample_labels) || throw(ArgumentError("sample_labels must match labels"))
    selected = [index for (index, label) in enumerate(labels) if label == ident_1 || label == ident_2]
    isempty(selected) && return nothing

    sample_order = String[]
    sample_to_index = Dict{String,Int}()
    sample_to_label = Dict{String,Int}()
    for index in selected
        sample = String(sample_labels[index])
        label = Int(labels[index])
        haskey(sample_to_index, sample) || begin
            push!(sample_order, sample)
            sample_to_index[sample] = length(sample_order)
        end
        if haskey(sample_to_label, sample)
            sample_to_label[sample] == label || throw(ArgumentError("each sample label must map to a single cell group for pseudobulk testing"))
        else
            sample_to_label[sample] = label
        end
    end

    aggregated = zeros(Int, length(experiment.gene_ids), length(sample_order))
    rows, cols, values = _matrix_triplets(experiment.counts)
    selected_set = Set(selected)
    for index in eachindex(values)
        column = cols[index]
        column in selected_set || continue
        sample = String(sample_labels[column])
        aggregated[rows[index], sample_to_index[sample]] += Int(values[index])
    end

    pseudobulk_labels = Symbol[sample_to_label[sample] == ident_1 ? :ident_1 : :ident_2 for sample in sample_order]
    return CountMatrix(aggregated, experiment.gene_ids, sample_order), pseudobulk_labels
end

function _wilcoxon_markers(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}, ident_1::Integer, ident_2::Union{Nothing,Integer}; min_total::Integer=10)
    if ident_2 === nothing
        group1 = [index for (index, label) in enumerate(labels) if label == ident_1]
        group2 = [index for (index, label) in enumerate(labels) if label != ident_1]
    else
        group1 = [index for (index, label) in enumerate(labels) if label == ident_1]
        group2 = [index for (index, label) in enumerate(labels) if label == ident_2]
    end

    isempty(group1) && return DEResult[]
    isempty(group2) && return DEResult[]

    expression = normalize_counts(experiment)
    gene_ids = experiment.gene_ids
    results = Vector{DEResult}(undef, length(gene_ids))
    pvalues = Vector{Float64}(undef, length(gene_ids))

    for gene_index in eachindex(gene_ids)
        x = _expression_for_cells(expression, gene_index, group1)
        y = _expression_for_cells(expression, gene_index, group2)
        base_mean = (sum(x) + sum(y)) / max(length(x) + length(y), 1)
        if base_mean == 0.0 || sum(x) + sum(y) < min_total
            results[gene_index] = DEResult(gene_ids[gene_index], base_mean, 0.0, 1.0, 0.0, 1.0, 1.0)
            pvalues[gene_index] = 1.0
            continue
        end

        stats = _wilcoxon_rank_sum(x, y)
        mean1 = mean(x)
        mean2 = mean(y)
        lfc = log2((mean1 + eps(Float64)) / (mean2 + eps(Float64)))
        lfc_se = sqrt(var(x; corrected=false) / max(length(x), 1) + var(y; corrected=false) / max(length(y), 1)) / log(2)
        results[gene_index] = DEResult(gene_ids[gene_index], base_mean, lfc, lfc_se, stats.stat, stats.pvalue, 1.0)
        pvalues[gene_index] = stats.pvalue
    end

    padj = benjamini_hochberg(pvalues)
    return [DEResult(result.gene_id, result.base_mean, result.log2_fold_change, result.lfc_se, result.stat, result.pvalue, padj[index]) for (index, result) in enumerate(results)]
end

function _deseq2_markers(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}, ident_1::Integer, ident_2::Union{Nothing,Integer}; min_total::Integer=10, shrink::Bool=true, normalization_method::Symbol=:tmm, dispersion_workflow::Symbol=:prior, modelMatrixType::Symbol=:standard, pseudobulk::Bool=false, sample_labels::Union{Nothing,AbstractVector{<:String}}=nothing)
    ident_2 === nothing && return _wilcoxon_markers(experiment, labels, ident_1, ident_2; min_total=min_total)
    if pseudobulk
        sample_labels === nothing && throw(ArgumentError("sample_labels are required for pseudobulk testing"))
        pseudobulk_counts = _pseudobulk_count_matrix(experiment, labels, ident_1, ident_2; sample_labels=sample_labels)
        pseudobulk_counts === nothing && return DEResult[]
        count_matrix, pseudobulk_design = pseudobulk_counts
        return differential_expression(count_matrix, pseudobulk_design; min_total=min_total, shrink=shrink, normalization_method=normalization_method, dispersion_workflow=dispersion_workflow, modelMatrixType=modelMatrixType)
    end

    subset = [index for (index, label) in enumerate(labels) if label == ident_1 || label == ident_2]
    isempty(subset) && return DEResult[]
    subset_matrix = experiment.counts isa CUDA.CuArray ? Array(experiment.counts[:, subset]) : experiment.counts[:, subset]
    subset_counts = CountMatrix(subset_matrix, experiment.gene_ids, experiment.cell_ids[subset])
    design = Symbol[label == ident_1 ? :ident_1 : :ident_2 for label in labels[subset]]
    return differential_expression(subset_counts, design; min_total=min_total, shrink=shrink, normalization_method=normalization_method, dispersion_workflow=dispersion_workflow, modelMatrixType=modelMatrixType)
end

function find_markers(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; ident_1::Integer, ident_2::Union{Nothing,Integer}=nothing, test::Symbol=:wilcox, pseudobulk::Bool=false, sample_labels::Union{Nothing,AbstractVector{<:String}}=nothing, min_total::Integer=10, shrink::Bool=true, normalization_method::Symbol=:tmm, dispersion_workflow::Symbol=:prior, modelMatrixType::Symbol=:standard)
    length(labels) == length(experiment.cell_ids) || throw(ArgumentError("labels must match the number of cells"))
    test in (:wilcox, :deseq2) || throw(ArgumentError("test must be :wilcox or :deseq2"))

    if test == :wilcox && !pseudobulk
        return _wilcoxon_markers(experiment, labels, ident_1, ident_2; min_total=min_total)
    end
    return _deseq2_markers(
        experiment,
        labels,
        ident_1,
        ident_2;
        min_total=min_total,
        shrink=shrink,
        normalization_method=normalization_method,
        dispersion_workflow=dispersion_workflow,
        modelMatrixType=modelMatrixType,
        pseudobulk=pseudobulk || test == :deseq2,
        sample_labels=sample_labels)
end

function _cell_embedding(experiment::SingleCellExperiment; n_components::Int=20)
    embedding = haskey(experiment.reductions, "pca") ? experiment.reductions["pca"] : run_pca(experiment; n_components=n_components)
    return Matrix{Float64}(embedding)
end

function _resolve_root_cell(experiment::SingleCellExperiment, root_cell)
    root_cell === nothing && return 1
    if root_cell isa Integer
        1 <= root_cell <= length(experiment.cell_ids) || throw(ArgumentError("root_cell index is out of bounds"))
        return Int(root_cell)
    elseif root_cell isa String
        index = findfirst(==(String(root_cell)), experiment.cell_ids)
        index === nothing && throw(ArgumentError("root_cell was not found in cell_ids"))
        return index
    end
    throw(ArgumentError("root_cell must be nothing, an integer cell index, or a cell identifier string"))
end

function _trajectory_weighted_graph(embedding::AbstractMatrix{<:Real}, neighbors::Vector{Vector{Int}})
    n_cells = size(embedding, 1)
    graph = SimpleWeightedGraph(n_cells)
    for source in 1:n_cells
        for target in neighbors[source]
            target <= source && continue
            weight = norm(view(embedding, source, :) .- view(embedding, target, :))
            isfinite(weight) || continue
            add_edge!(graph, SimpleWeightedEdge(source, target, max(Float64(weight), eps(Float64))))
        end
    end
    return graph
end

function _mst_pseudotime(graph::SimpleWeightedGraph, root_index::Int)
    mst_edges = kruskal_mst(graph)
    mst = SimpleWeightedGraph(nv(graph))
    for edge in mst_edges
        add_edge!(mst, edge)
    end

    shortest_paths = dijkstra_shortest_paths(mst, root_index, Graphs.weights(mst))
    distances = collect(shortest_paths.dists)
    finite_distances = filter(isfinite, distances)
    max_distance = isempty(finite_distances) ? 0.0 : maximum(finite_distances)
    pseudotime = Vector{Float64}(undef, length(distances))
    for index in eachindex(distances)
        distance = distances[index]
        if !isfinite(distance)
            pseudotime[index] = NaN
        elseif max_distance > 0
            pseudotime[index] = distance / max_distance
        else
            pseudotime[index] = 0.0
        end
    end

    return (pseudotime=pseudotime, mst=mst, distances=distances)
end

"""
    calculate_pseudotime(experiment; root_cell=nothing, embedding=nothing, reduction="pca", graph_name="neighbors", k=15)

Infer a trajectory ordering by building a kNN graph, extracting its minimum spanning tree,
and measuring shortest-path distance from the chosen root cell.

The returned pseudotime is normalized to the range [0, 1] on the connected component
containing the root cell. Cells unreachable from the root are assigned `NaN`.
"""
function calculate_pseudotime(experiment::SingleCellExperiment; root_cell=nothing, embedding::Union{Nothing,AbstractMatrix}=nothing, reduction::String="pca", graph_name::String="neighbors", k::Int=15, use_cuda::Bool=false)
    root_index = _resolve_root_cell(experiment, root_cell)
    matrix = if embedding === nothing
        haskey(experiment.reductions, reduction) ? experiment.reductions[reduction] : run_pca(experiment; n_components=max(2, k), use_cuda=use_cuda)
    else
        Matrix{Float64}(embedding)
    end

    neighbors = haskey(experiment.neighbors, graph_name) ? experiment.neighbors[graph_name] : find_neighbors(experiment; embedding=matrix, k=k, graph_name=graph_name, use_cuda=use_cuda)
    size(matrix, 1) == length(neighbors) || throw(ArgumentError("embedding must have one row per cell"))

    weighted_graph = _trajectory_weighted_graph(matrix, neighbors)
    result = _mst_pseudotime(weighted_graph, root_index)
    experiment.metadata["pseudotime_root"] = experiment.cell_ids[root_index]
    experiment.reductions["pseudotime"] = reshape(result.pseudotime, :, 1)
    provenance = provenance_record(
        "PseudotimeResult",
        "SingleCell/calculate_pseudotime";
        notes=["computed pseudotime from an MST over a neighbor graph"],
        parameters=(root_index=Int(root_index), root_cell=experiment.cell_ids[root_index], graph_name=graph_name, reduction=reduction, k=Int(k), use_cuda=Bool(use_cuda)))
    update_provenance!(experiment.metadata; source="SingleCell/calculate_pseudotime", notes=["computed pseudotime from an MST over a neighbor graph"], parameters=(root_index=Int(root_index), root_cell=experiment.cell_ids[root_index], graph_name=graph_name, reduction=reduction, k=Int(k), use_cuda=Bool(use_cuda)))
    return (
        pseudotime=result.pseudotime,
        mst=result.mst,
        distances=result.distances,
        root_index=root_index,
        root_cell=experiment.cell_ids[root_index],
        graph_name=graph_name,
        provenance=provenance)
end

"""
    detect_doublets(experiment; ...)

Estimate doublet scores by simulating synthetic cells in embedding space and
measuring the density of synthetic neighbors around each observed cell.
"""
function detect_doublets(experiment::SingleCellExperiment; n_simulated::Int=500, k::Int=15, threshold::Union{Nothing,Real}=nothing, n_components::Int=20, random_seed::Int=1)
    n_simulated >= 0 || throw(ArgumentError("n_simulated must be nonnegative"))
    k >= 0 || throw(ArgumentError("k must be nonnegative"))
    embedding = _cell_embedding(experiment; n_components=n_components)
    n_cells = size(embedding, 1)
    if n_cells == 0
        provenance = provenance_record(
            "DoubletDetectionResult",
            "SingleCell/detect_doublets";
            status=:warning,
            notes=["no cells were available for doublet detection"],
            fallbacks=["returned empty doublet scores"],
            parameters=(n_simulated=Int(n_simulated), k=Int(k), n_components=Int(n_components), random_seed=Int(random_seed)))
        update_provenance!(experiment.metadata; source="SingleCell/detect_doublets", status=:warning, notes=["doublet detection received no cells"], fallbacks=["returned empty doublet scores"], parameters=(n_simulated=Int(n_simulated), k=Int(k), n_components=Int(n_components), random_seed=Int(random_seed)))
        return (scores=Float64[], threshold=0.0, doublets=Bool[], simulated_embedding=zeros(Float64, 0, size(embedding, 2)), provenance=provenance)
    end

    rng = MersenneTwister(random_seed)
    simulated = Matrix{Float64}(undef, n_simulated, size(embedding, 2))
    for index in 1:n_simulated
        left = rand(rng, 1:n_cells)
        right = rand(rng, 1:n_cells)
        simulated[index, :] .= (view(embedding, left, :) .+ view(embedding, right, :)) ./ 2
    end

    combined = vcat(embedding, simulated)
    scores = zeros(Float64, n_cells)
    for cell_index in 1:n_cells
        distances = [(other_index, norm(view(combined, cell_index, :) .- view(combined, other_index, :))) for other_index in 1:size(combined, 1) if other_index != cell_index]
        sort!(distances; by = last)
        neighborhood_size = min(k, length(distances))
        nearest = Iterators.take(distances, neighborhood_size)
        synthetic_neighbors = count(pair -> pair[1] > n_cells, nearest)
        scores[cell_index] = neighborhood_size == 0 ? 0.0 : synthetic_neighbors / neighborhood_size
    end

    cutoff = threshold === nothing ? (isempty(scores) ? 0.0 : max(mean(scores) + std(scores), quantile(scores, 0.90))) : Float64(threshold)
    doublets = scores .>= cutoff
    provenance = provenance_record(
        "DoubletDetectionResult",
        "SingleCell/detect_doublets";
        notes=["estimated doublet scores using simulated neighbors"],
        parameters=(n_simulated=Int(n_simulated), k=Int(k), n_components=Int(n_components), random_seed=Int(random_seed), threshold=threshold === nothing ? "auto" : Float64(threshold), n_cells=Int(n_cells)))
    update_provenance!(experiment.metadata; source="SingleCell/detect_doublets", notes=["estimated doublet scores using simulated neighbors"], parameters=(n_simulated=Int(n_simulated), k=Int(k), n_components=Int(n_components), random_seed=Int(random_seed), threshold=threshold === nothing ? "auto" : Float64(threshold), n_cells=Int(n_cells)))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (scores=scores, threshold=cutoff, doublets=doublets, simulated_embedding=simulated, provenance=provenance), "detect_doublets")
end

"""
    integrate_batches(matrix, batch_labels; method=:harmony, n_components=20)

Apply a simple batch-centering integration to a gene-by-cell matrix. The
result includes the corrected matrix, an embedding for downstream analysis, and
per-batch centering offsets.
"""
function integrate_batches(matrix::AbstractMatrix{<:Real}, batch_labels; method::Symbol=:harmony, n_components::Int=20)
    values = Matrix{Float64}(matrix)
    length(batch_labels) == size(values, 2) || throw(ArgumentError("batch_labels must match the number of columns"))

    batches = unique(String.(batch_labels))
    global_mean = mean(values, dims=2)
    corrected = similar(values)
    batch_means = Dict{String,Vector{Float64}}()

    for batch in batches
        indices = findall(label -> String(label) == batch, batch_labels)
        isempty(indices) && continue
        batch_mean = vec(mean(values[:, indices], dims=2))
        batch_means[batch] = batch_mean
        corrected[:, indices] .= values[:, indices] .- batch_mean .+ vec(global_mean)
    end

    embedding = _pca_for_matrix(permutedims(corrected); n_components=n_components)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (corrected_matrix=corrected, embedding=embedding, batch_means=batch_means, method=method), "integrate_batches")
end

"""
    integrate_data(experiments; method=:harmony, n_components=20)

Combine multiple `SingleCellExperiment` objects with matching genes, batch-correct
the concatenated count matrix, and return the corrected matrix, embedding, and
batch labels.
"""
function integrate_data(experiments::AbstractVector{<:SingleCellExperiment}; method::Symbol=:harmony, n_components::Int=20)
    isempty(experiments) && throw(ArgumentError("at least one experiment is required"))
    reference_genes = first(experiments).gene_ids
    for experiment in experiments
        experiment.gene_ids == reference_genes || throw(ArgumentError("all experiments must share the same gene ordering"))
    end

    counts = hcat([Matrix{Float64}(experiment.counts) for experiment in experiments]...)
    batch_labels = String[]
    cell_ids = String[]
    for (index, experiment) in enumerate(experiments)
        append!(batch_labels, fill("batch$(index)", length(experiment.cell_ids)))
        append!(cell_ids, experiment.cell_ids)
    end

    integrated = integrate_batches(counts, batch_labels; method=method, n_components=n_components)
    return (
        corrected_matrix=integrated.corrected_matrix,
        embedding=integrated.embedding,
        batch_means=integrated.batch_means,
        method=integrated.method,
        gene_ids=copy(reference_genes),
        cell_ids=cell_ids,
        batch_labels=batch_labels)
end

function _pca_for_matrix(matrix::AbstractMatrix{<:Real}; n_components::Int=20)
    centered = Matrix{Float64}(matrix)
    centered .-= mean(centered, dims=1)
    svd_result = svd(centered; full=false)
    components = min(n_components, size(svd_result.U, 2))
    return svd_result.U[:, 1:components] * Diagonal(svd_result.S[1:components])
end

function materialize_singlecell_experiment(experiment::SingleCellExperiment)
    counts = experiment.counts isa SparseMatrixCSC ? experiment.counts : sparse(Int.(Array(experiment.counts)))
    metadata = _singlecell_metadata_with_provenance(experiment.metadata; source="materialize_singlecell_experiment", notes=["materialized counts to CPU sparse matrix"], parameters=(n_genes=length(experiment.gene_ids), n_cells=length(experiment.cell_ids)))
    return SingleCellExperiment(counts, experiment.gene_ids, experiment.cell_ids; metadata=metadata, spatial_coords=experiment.spatial_coords)
end

function gpu_singlecell_experiment(experiment::SingleCellExperiment)
    CUDA.functional() || throw(ArgumentError("CUDA is not available"))
    counts = experiment.counts isa CUDA.CuArray ? experiment.counts : CUDA.CuArray{Int}(Array(experiment.counts))
    metadata = _singlecell_metadata_with_provenance(experiment.metadata; source="gpu_singlecell_experiment", notes=["moved counts to CUDA"], parameters=(n_genes=length(experiment.gene_ids), n_cells=length(experiment.cell_ids)))
    return SingleCellExperiment{typeof(counts)}(counts, copy(experiment.gene_ids), copy(experiment.cell_ids), experiment.spatial_coords === nothing ? nothing : copy(experiment.spatial_coords), metadata, Dict{String,Matrix{Float64}}(), Dict{String,Vector{Int}}(), Dict{String,Vector{String}}(), Dict{String,Vector{Vector{Int}}}())
end

function save_singlecell_experiment(path_prefix::String, experiment::SingleCellExperiment)
    counts_path = string(path_prefix, ".counts.bin")
    metadata_path = string(path_prefix, ".metadata.bin")
    dense_counts = Array{Int}(experiment.counts)
    metadata = _singlecell_metadata_with_provenance(experiment.metadata; source="save_singlecell_experiment", notes=["serialized counts and metadata"], parameters=(counts_path=counts_path, metadata_path=metadata_path))
    open(counts_path, "w") do io
        write(io, dense_counts)
    end
    open(metadata_path, "w") do io
        serialize(io, (
            gene_ids=copy(experiment.gene_ids),
            cell_ids=copy(experiment.cell_ids),
            spatial_coords=experiment.spatial_coords === nothing ? nothing : copy(experiment.spatial_coords),
            metadata=metadata,
            size=size(dense_counts)))
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (counts_path=counts_path, metadata_path=metadata_path), "save_singlecell_experiment")
end

function load_singlecell_experiment(path_prefix::String)
    counts_path = string(path_prefix, ".counts.bin")
    metadata_path = string(path_prefix, ".metadata.bin")
    payload = open(metadata_path, "r") do io
        deserialize(io)
    end
    counts = open(counts_path, "r+") do io
        Mmap.mmap(io, Matrix{Int}, payload.size)
    end
    spatial_coords = hasproperty(payload, :spatial_coords) ? payload.spatial_coords : nothing
    metadata = _singlecell_metadata_with_provenance(payload.metadata; source="load_singlecell_experiment", notes=["loaded from $(metadata_path)"], parameters=(counts_path=counts_path, metadata_path=metadata_path))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, SingleCellExperiment{typeof(counts)}(counts, payload.gene_ids, payload.cell_ids, spatial_coords, metadata, Dict{String,Matrix{Float64}}(), Dict{String,Vector{Int}}(), Dict{String,Vector{String}}(), Dict{String,Vector{Vector{Int}}}()), "load_singlecell_experiment")
end

"""
    SingleCellArchive

High-level archive that bundles counts, layers, metadata, embeddings, and cached
analysis state into one serialized artifact.
"""
struct SingleCellArchive
    counts::Any
    gene_ids::Vector{String}
    cell_ids::Vector{String}
    spatial_coords::Union{Nothing,Matrix{Float64}}
    metadata::Dict{String,Any}
    layers::Dict{String,Any}
    embeddings::Dict{String,Matrix{Float64}}
    clusters::Dict{String,Vector{Int}}
    variable_features::Dict{String,Vector{String}}
    neighbors::Dict{String,Vector{Vector{Int}}}
end

function _archive_backend(value)
    value isa CUDA.CuArray && return Array(value)
    value isa CountMatrix && return CountMatrix(_archive_backend(value.counts), copy(value.gene_ids), copy(value.sample_ids))
    value isa AbstractMatrix && return copy(value)
    return value
end

function _archive_embeddings(embeddings::AbstractDict)
    return Dict{String,Matrix{Float64}}(string(name) => Matrix{Float64}(copy(matrix)) for (name, matrix) in embeddings)
end

function _archive_layers(layers::AbstractDict)
    return Dict{String,Any}(string(name) => _archive_backend(layer) for (name, layer) in layers)
end

function save_singlecell_archive(path::String, experiment::SingleCellExperiment; layers::AbstractDict=Dict{String,Any}(), embeddings::AbstractDict=experiment.reductions)
    metadata = _singlecell_metadata_with_provenance(experiment.metadata; source="save_singlecell_archive", notes=["serialized archive contents"], parameters=(path=path, layers=length(layers), embeddings=length(embeddings)))
    archive = SingleCellArchive(
        _archive_backend(experiment.counts),
        copy(experiment.gene_ids),
        copy(experiment.cell_ids),
        experiment.spatial_coords === nothing ? nothing : copy(experiment.spatial_coords),
        metadata,
        _archive_layers(layers),
        _archive_embeddings(embeddings),
        Dict{String,Vector{Int}}(name => copy(labels) for (name, labels) in experiment.clusters),
        Dict{String,Vector{String}}(name => copy(features) for (name, features) in experiment.variable_features),
        Dict{String,Vector{Vector{Int}}}(name => [copy(neighbor_list) for neighbor_list in neighbors] for (name, neighbors) in experiment.neighbors))
    open(path, "w") do io
        serialize(io, archive)
    end
    return archive
end

function load_singlecell_archive(path::String)
    archive = open(path, "r") do io
        deserialize(io)::SingleCellArchive
    end
    update_provenance!(archive.metadata; source="load_singlecell_archive", notes=["loaded from $(path)"], parameters=(path=path, layers=length(archive.layers), embeddings=length(archive.embeddings)))
    return archive
end

function archive_to_singlecell_experiment(archive::SingleCellArchive)
    metadata = _singlecell_metadata_with_provenance(archive.metadata; source="archive_to_singlecell_experiment", notes=["reconstituted from SingleCellArchive"], parameters=(layers=length(archive.layers), embeddings=length(archive.embeddings)))
    experiment = SingleCellExperiment(archive.counts, archive.gene_ids, archive.cell_ids; metadata=metadata, spatial_coords=archive.spatial_coords)
    experiment.reductions = Dict{String,Matrix{Float64}}(name => copy(matrix) for (name, matrix) in archive.embeddings)
    experiment.clusters = Dict{String,Vector{Int}}(name => copy(labels) for (name, labels) in archive.clusters)
    experiment.variable_features = Dict{String,Vector{String}}(name => copy(features) for (name, features) in archive.variable_features)
    experiment.neighbors = Dict{String,Vector{Vector{Int}}}(name => [copy(neighbor_list) for neighbor_list in neighbors] for (name, neighbors) in archive.neighbors)
    if !isempty(archive.layers)
        experiment.metadata["layers"] = archive.layers
    end
    return experiment
end

function _show_singlecell_archive(io::IO, archive::SingleCellArchive)
    print(io, "SingleCellArchive(", length(archive.gene_ids), " genes, ", length(archive.cell_ids), " cells")
    if archive.spatial_coords !== nothing
        print(io, ", spatial=", size(archive.spatial_coords, 1), "x", size(archive.spatial_coords, 2))
    end
    if !isempty(archive.layers)
        print(io, ", layers=", join(sort!(collect(keys(archive.layers))), ", "))
    end
    if !isempty(archive.embeddings)
        print(io, ", embeddings=", join(sort!(collect(keys(archive.embeddings))), ", "))
    end
    provenance = metadata_provenance(archive.metadata)
    provenance !== nothing && print(io, ", provenance=", provenance_summary(provenance))
    print(io, ")")
end

function Base.show(io::IO, archive::SingleCellArchive)
    _show_singlecell_archive(io, archive)
end

function Base.show(io::IO, ::MIME"text/plain", archive::SingleCellArchive)
    _show_singlecell_archive(io, archive)
end

"""
    SingleCellArchiveBrowser

Lightweight archive browser that exposes the available layers and embeddings
while keeping a preview of the selected items ready for display.
"""
struct SingleCellArchiveBrowser
    archive::SingleCellArchive
    selected_layer::Union{Nothing,String}
    selected_embedding::Union{Nothing,String}
    preview_rows::Int
    preview_cols::Int
end

function browse_singlecell_archive(source::Union{SingleCellArchive,String,SingleCellExperiment}; layer::Union{Nothing,String}=nothing, embedding::Union{Nothing,String}=nothing, preview::Int=5)
    archive = if source isa SingleCellArchive
        source
    elseif source isa SingleCellExperiment
        SingleCellArchive(
            _archive_backend(source.counts),
            copy(source.gene_ids),
            copy(source.cell_ids),
            source.spatial_coords === nothing ? nothing : copy(source.spatial_coords),
            Dict{String,Any}(source.metadata),
            Dict{String,Any}(),
            _archive_embeddings(source.reductions),
            Dict{String,Vector{Int}}(name => copy(labels) for (name, labels) in source.clusters),
            Dict{String,Vector{String}}(name => copy(features) for (name, features) in source.variable_features),
            Dict{String,Vector{Vector{Int}}}(name => [copy(neighbor_list) for neighbor_list in neighbors] for (name, neighbors) in source.neighbors))
    else
        load_singlecell_archive(String(source))
    end
    return SingleCellArchiveBrowser(archive, layer === nothing ? nothing : String(layer), embedding === nothing ? nothing : String(embedding), max(preview, 1), max(preview, 1))
end

function _browser_preview(matrix::AbstractMatrix{<:Real}, rows::Int, cols::Int)
    return matrix[1:min(rows, size(matrix, 1)), 1:min(cols, size(matrix, 2))]
end

function archive_layer(browser::SingleCellArchiveBrowser, name::String)
    haskey(browser.archive.layers, String(name)) || throw(ArgumentError("layer '$name' was not found"))
    return browser.archive.layers[String(name)]
end

function archive_embedding(browser::SingleCellArchiveBrowser, name::String)
    haskey(browser.archive.embeddings, String(name)) || throw(ArgumentError("embedding '$name' was not found"))
    return browser.archive.embeddings[String(name)]
end

function _show_singlecell_archive_browser(io::IO, browser::SingleCellArchiveBrowser)
    archive = browser.archive
    println(io, "SingleCellArchiveBrowser(", length(archive.gene_ids), " genes, ", length(archive.cell_ids), " cells)")
    println(io, "  spatial coordinates: ", archive.spatial_coords === nothing ? "none" : string(size(archive.spatial_coords, 1), "x", size(archive.spatial_coords, 2)))
    println(io, "  layers: ", isempty(archive.layers) ? "none" : join(sort!(collect(keys(archive.layers))), ", "))
    println(io, "  embeddings: ", isempty(archive.embeddings) ? "none" : join(sort!(collect(keys(archive.embeddings))), ", "))
    println(io, "  cluster caches: ", isempty(archive.clusters) ? "none" : join(sort!(collect(keys(archive.clusters))), ", "))
    provenance = metadata_provenance(archive.metadata)
    provenance !== nothing && println(io, "  provenance: ", provenance_summary(provenance))
    if browser.selected_layer !== nothing && haskey(archive.layers, browser.selected_layer)
        selected = archive.layers[browser.selected_layer]
        if selected isa AbstractMatrix
            println(io, "  selected layer: ", browser.selected_layer, " preview=", size(_browser_preview(selected, browser.preview_rows, browser.preview_cols)))
        else
            println(io, "  selected layer: ", browser.selected_layer, " type=", typeof(selected))
        end
    end
    if browser.selected_embedding !== nothing && haskey(archive.embeddings, browser.selected_embedding)
        selected = archive.embeddings[browser.selected_embedding]
        println(io, "  selected embedding: ", browser.selected_embedding, " preview=", size(_browser_preview(selected, browser.preview_rows, browser.preview_cols)))
    end
end

function Base.show(io::IO, browser::SingleCellArchiveBrowser)
    _show_singlecell_archive_browser(io, browser)
end

function Base.show(io::IO, ::MIME"text/plain", browser::SingleCellArchiveBrowser)
    _show_singlecell_archive_browser(io, browser)
end

"""
    RNAVelocityResult

Summary of a simple spliced/unspliced RNA velocity fit.
"""
struct RNAVelocityResult <: AbstractAnalysisResult
    gene_ids::Vector{String}
    cell_ids::Vector{String}
    gene_slopes::Vector{Float64}
    velocity::Matrix{Float64}
    cell_scores::Vector{Float64}
    latent_time::Vector{Float64}
    provenance::ResultProvenance
end

RNAVelocityResult(cell_ids, gene_slopes, velocity, cell_scores, latent_time) =
    RNAVelocityResult(cell_ids, gene_slopes, velocity, cell_scores, latent_time, provenance_record("RNAVelocityResult", "singlecell"))

RNAVelocityResult(gene_ids::Vector{String}, cell_ids::Vector{String}, gene_slopes::Vector{Float64}, velocity::Matrix{Float64}, cell_scores::Vector{Float64}, latent_time::Vector{Float64}) =
    RNAVelocityResult(gene_ids, cell_ids, gene_slopes, velocity, cell_scores, latent_time, provenance_record("RNAVelocityResult", "SingleCell/calculate_rna_velocity"))

function Base.show(io::IO, result::RNAVelocityResult)
    print(io, analysis_result_summary(result))
end

function Base.show(io::IO, ::MIME"text/plain", result::RNAVelocityResult)
    print(io, analysis_result_summary(result))
end

function _resolve_velocity_layer(experiment::SingleCellExperiment, layer, key::String)
    if layer !== nothing
        return layer isa CountMatrix ? layer.counts : layer
    end
    haskey(experiment.metadata, String(key)) || throw(ArgumentError("missing velocity layer '$key' in metadata"))
    stored = experiment.metadata[String(key)]
    return stored isa CountMatrix ? stored.counts : stored
end

function attach_velocity_layers!(experiment::SingleCellExperiment; spliced, unspliced, spliced_key::String="spliced_counts", unspliced_key::String="unspliced_counts")
    experiment.metadata[String(spliced_key)] = spliced
    experiment.metadata[String(unspliced_key)] = unspliced
    update_provenance!(experiment.metadata; source="attach_velocity_layers!", notes=["attached velocity layers $(spliced_key)/$(unspliced_key)"], parameters=(spliced_key=spliced_key, unspliced_key=unspliced_key))
    return experiment
end

function _velocity_dense(matrix)
    return Array{Float64}(matrix)
end

function calculate_rna_velocity(experiment::SingleCellExperiment; spliced=nothing, unspliced=nothing, spliced_key::String="spliced_counts", unspliced_key::String="unspliced_counts", normalize::Bool=true, use_cuda::Bool=false)
    spliced_matrix = _resolve_velocity_layer(experiment, spliced, spliced_key)
    unspliced_matrix = _resolve_velocity_layer(experiment, unspliced, unspliced_key)
    size(spliced_matrix) == size(unspliced_matrix) || throw(DimensionMismatch("spliced and unspliced layers must have matching dimensions"))
    size(spliced_matrix, 1) == length(experiment.gene_ids) || throw(DimensionMismatch("velocity layers must match the gene count"))
    size(spliced_matrix, 2) == length(experiment.cell_ids) || throw(DimensionMismatch("velocity layers must match the cell count"))

    if normalize
        spliced_matrix = normalize_counts(SingleCellExperiment(spliced_matrix, experiment.gene_ids, experiment.cell_ids; metadata=experiment.metadata, spatial_coords=experiment.spatial_coords); use_cuda=use_cuda)
        unspliced_matrix = normalize_counts(SingleCellExperiment(unspliced_matrix, experiment.gene_ids, experiment.cell_ids; metadata=experiment.metadata, spatial_coords=experiment.spatial_coords); use_cuda=use_cuda)
    end

    if use_cuda
        CUDA.functional() || throw(ArgumentError("CUDA is not available"))
        spliced_gpu = CUDA.CuArray{Float32}(Array(spliced_matrix))
        unspliced_gpu = CUDA.CuArray{Float32}(Array(unspliced_matrix))
        gene_power = vec(sum(spliced_gpu .* spliced_gpu, dims=2))
        gene_cross = vec(sum(spliced_gpu .* unspliced_gpu, dims=2))
        slopes = Array(gene_cross ./ max.(gene_power, eps(Float32)))
        slope_gpu = reshape(CUDA.CuArray{Float32}(slopes), :, 1)
        velocity = Array(unspliced_gpu .- slope_gpu .* spliced_gpu)
    else
        spliced_dense = _velocity_dense(spliced_matrix)
        unspliced_dense = _velocity_dense(unspliced_matrix)
        slopes = Vector{Float64}(undef, size(spliced_dense, 1))
        velocity = similar(spliced_dense)
        for gene_index in 1:size(spliced_dense, 1)
            x = view(spliced_dense, gene_index, :)
            y = view(unspliced_dense, gene_index, :)
            denominator = sum(abs2, x)
            slopes[gene_index] = denominator <= eps(Float64) ? 0.0 : sum(x .* y) / denominator
            velocity[gene_index, :] .= y .- slopes[gene_index] .* x
        end
    end

    cell_scores = vec(mean(velocity, dims=1))
    if isempty(cell_scores)
        latent_time = Float64[]
    else
        minimum_score = minimum(cell_scores)
        score_range = maximum(cell_scores) - minimum_score
        latent_time = score_range > 0 ? (cell_scores .- minimum_score) ./ score_range : zeros(Float64, length(cell_scores))
    end
    return RNAVelocityResult(
        copy(experiment.gene_ids),
        copy(experiment.cell_ids),
        Float64.(slopes),
        velocity,
        cell_scores,
        latent_time,
        provenance_record(
            "RNAVelocityResult",
            "SingleCell/calculate_rna_velocity";
            notes=["gene-wise slope fit between spliced and unspliced layers"],
            parameters=(gene_count=length(experiment.gene_ids), cell_count=length(experiment.cell_ids), normalize=Bool(normalize), use_cuda=Bool(use_cuda), spliced_source=spliced === nothing ? spliced_key : "explicit_argument", unspliced_source=unspliced === nothing ? unspliced_key : "explicit_argument")))
end

function _scale_to_unit_interval(values::AbstractVector{<:Real})
    isempty(values) && return Float64[]
    dense = Float64.(values)
    minimum_value = minimum(dense)
    value_range = maximum(dense) - minimum_value
    return value_range > 0 ? (dense .- minimum_value) ./ value_range : zeros(Float64, length(dense))
end

function _strictly_increasing_times(times::AbstractVector{<:Real})
    dense = Float64.(times)
    for index in eachindex(dense)
        index == firstindex(dense) && continue
        previous_index = prevind(dense, index)
        previous_time = dense[previous_index]
        dense[index] <= previous_time && (dense[index] = previous_time + 1e-6)
    end
    return dense
end

function _dynamical_velocity_ode!(du, u, p, t)
    alpha, beta, gamma = p
    du[1] = alpha - beta * u[1]
    du[2] = beta * u[1] - gamma * u[2]
    return nothing
end

const _DYNAMICAL_VELOCITY_LOG_PARAM_MIN = -6.0
const _DYNAMICAL_VELOCITY_LOG_PARAM_MAX = 6.0
const _DYNAMICAL_VELOCITY_ABSTOL = 1e-7
const _DYNAMICAL_VELOCITY_RELTOL = 1e-7
const _DYNAMICAL_VELOCITY_MAXITERS = 100_000

@inline function _bounded_dynamical_velocity_parameters(log_parameters)
    bounded = clamp.(Float64.(log_parameters), _DYNAMICAL_VELOCITY_LOG_PARAM_MIN, _DYNAMICAL_VELOCITY_LOG_PARAM_MAX)
    return exp.(bounded)
end

"""
    DynamicalRNAVelocityResult

Per-gene kinetic RNA velocity fit using a simple scVelo-style time-ordered model.
"""
struct DynamicalRNAVelocityResult <: AbstractAnalysisResult
    gene_ids::Vector{String}
    cell_ids::Vector{String}
    alpha::Vector{Float64}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    velocity::Matrix{Float64}
    cell_scores::Vector{Float64}
    latent_time::Vector{Float64}
    provenance::ResultProvenance
end

DynamicalRNAVelocityResult(cell_ids, alpha, beta, gamma, velocity, cell_scores, latent_time) =
    DynamicalRNAVelocityResult(cell_ids, alpha, beta, gamma, velocity, cell_scores, latent_time, provenance_record("DynamicalRNAVelocityResult", "singlecell"))

DynamicalRNAVelocityResult(gene_ids::Vector{String}, cell_ids::Vector{String}, alpha::Vector{Float64}, beta::Vector{Float64}, gamma::Vector{Float64}, velocity::Matrix{Float64}, cell_scores::Vector{Float64}, latent_time::Vector{Float64}) =
    DynamicalRNAVelocityResult(gene_ids, cell_ids, alpha, beta, gamma, velocity, cell_scores, latent_time, provenance_record("DynamicalRNAVelocityResult", "SingleCell/calculate_dynamical_rna_velocity"))

function Base.show(io::IO, result::DynamicalRNAVelocityResult)
    println(io, analysis_result_summary(result))
    alpha_mean = isempty(result.alpha) ? 0.0 : mean(result.alpha)
    beta_mean = isempty(result.beta) ? 0.0 : mean(result.beta)
    gamma_mean = isempty(result.gamma) ? 0.0 : mean(result.gamma)
    println(io, "  mean alpha=", round(alpha_mean; digits=3), ", mean beta=", round(beta_mean; digits=3), ", mean gamma=", round(gamma_mean; digits=3))
end

function _show_dynamical_rna_velocity(io::IO, result::DynamicalRNAVelocityResult)
    println(io, analysis_result_summary(result))
    alpha_mean = isempty(result.alpha) ? 0.0 : mean(result.alpha)
    beta_mean = isempty(result.beta) ? 0.0 : mean(result.beta)
    gamma_mean = isempty(result.gamma) ? 0.0 : mean(result.gamma)
    println(io, "  mean alpha=", round(alpha_mean; digits=3), ", mean beta=", round(beta_mean; digits=3), ", mean gamma=", round(gamma_mean; digits=3))
end

function Base.show(io::IO, ::MIME"text/plain", result::DynamicalRNAVelocityResult)
    _show_dynamical_rna_velocity(io, result)
end

function _kinetic_velocity_fit(spliced::AbstractVector{<:Real}, unspliced::AbstractVector{<:Real}, latent_time::AbstractVector{<:Real})
    n_cells = length(spliced)
    n_cells == length(unspliced) == length(latent_time) || throw(DimensionMismatch("spliced, unspliced, and latent_time must have matching lengths"))
    n_cells < 2 && return 0.0, 0.0, 0.0, zeros(Float64, n_cells)

    order = sortperm(Float64.(latent_time))
    ordered_time = _strictly_increasing_times(latent_time[order])
    ordered_time .-= ordered_time[1]
    ordered_spliced = max.(Float64.(spliced[order]), 0.0)
    ordered_unspliced = max.(Float64.(unspliced[order]), 0.0)
    timespan = ordered_time[end]
    timespan <= eps(Float64) && return 0.0, 0.0, 0.0, zeros(Float64, n_cells)

    initial_state = [ordered_unspliced[1], ordered_spliced[1]]
    observed = vcat(ordered_unspliced, ordered_spliced)
    solver = Rosenbrock23()

    function objective(log_parameters)
        alpha, beta, gamma = _bounded_dynamical_velocity_parameters(log_parameters)
        problem = ODEProblem(_dynamical_velocity_ode!, initial_state, (0.0, timespan), (alpha, beta, gamma))
        solution = try
            solve(problem, solver; saveat=ordered_time, abstol=_DYNAMICAL_VELOCITY_ABSTOL, reltol=_DYNAMICAL_VELOCITY_RELTOL, maxiters=_DYNAMICAL_VELOCITY_MAXITERS)
        catch
            return Inf
        end
        length(solution.u) == n_cells || return Inf
        predicted_unspliced = [state[1] for state in solution.u]
        predicted_spliced = [state[2] for state in solution.u]
        residual = vcat(predicted_unspliced .- ordered_unspliced, predicted_spliced .- ordered_spliced)
        return sum(abs2, residual) / length(residual)
    end

    initial_guess = log.([max(mean(ordered_unspliced), eps(Float64)), max(0.5, 1.0 / max(timespan, 1.0)), max(0.5, 1.0 / max(timespan, 1.0))])
    fit = optimize(objective, initial_guess, NelderMead(), Optim.Options(iterations=200, store_trace=false, show_trace=false))
    alpha, beta, gamma = _bounded_dynamical_velocity_parameters(Optim.minimizer(fit))

    fitted_problem = ODEProblem(_dynamical_velocity_ode!, initial_state, (0.0, timespan), (alpha, beta, gamma))
    fitted_solution = solve(fitted_problem, solver; saveat=ordered_time, abstol=_DYNAMICAL_VELOCITY_ABSTOL, reltol=_DYNAMICAL_VELOCITY_RELTOL, maxiters=_DYNAMICAL_VELOCITY_MAXITERS)
    velocity = zeros(Float64, n_cells)
    for (index, state) in enumerate(fitted_solution.u)
        velocity[index] = beta * state[1] - gamma * state[2]
    end
    return alpha, beta, gamma, velocity[invperm(order)]
end

"""
    calculate_dynamical_rna_velocity(experiment; ...)

Fit a scVelo-style kinetic approximation by ordering cells along latent time,
estimating per-gene alpha/beta/gamma parameters, and returning a velocity
matrix together with the fitted cell-level ordering.
"""
function calculate_dynamical_rna_velocity(experiment::SingleCellExperiment; spliced=nothing, unspliced=nothing, spliced_key::String="spliced_counts", unspliced_key::String="unspliced_counts", latent_time::Union{Nothing,AbstractVector}=nothing, normalize::Bool=true, use_cuda::Bool=false)
    spliced_matrix = _resolve_velocity_layer(experiment, spliced, spliced_key)
    unspliced_matrix = _resolve_velocity_layer(experiment, unspliced, unspliced_key)
    size(spliced_matrix) == size(unspliced_matrix) || throw(DimensionMismatch("spliced and unspliced layers must have matching dimensions"))
    size(spliced_matrix, 1) == length(experiment.gene_ids) || throw(DimensionMismatch("velocity layers must match the gene count"))
    size(spliced_matrix, 2) == length(experiment.cell_ids) || throw(DimensionMismatch("velocity layers must match the cell count"))
    latent_time === nothing || length(latent_time) == length(experiment.cell_ids) || throw(DimensionMismatch("latent_time must match the cell count"))

    if normalize
        spliced_matrix = normalize_counts(SingleCellExperiment(spliced_matrix, experiment.gene_ids, experiment.cell_ids; metadata=experiment.metadata, spatial_coords=experiment.spatial_coords); use_cuda=use_cuda)
        unspliced_matrix = normalize_counts(SingleCellExperiment(unspliced_matrix, experiment.gene_ids, experiment.cell_ids; metadata=experiment.metadata, spatial_coords=experiment.spatial_coords); use_cuda=use_cuda)
    end

    base_result = latent_time === nothing ? calculate_rna_velocity(experiment; spliced=spliced_matrix, unspliced=unspliced_matrix, normalize=false, use_cuda=use_cuda) : nothing
    fitted_latent_time = latent_time === nothing ? base_result.latent_time : _scale_to_unit_interval(latent_time)

    spliced_dense = _velocity_dense(spliced_matrix)
    unspliced_dense = _velocity_dense(unspliced_matrix)
    n_genes, n_cells = size(spliced_dense)
    alpha = zeros(Float64, n_genes)
    beta = zeros(Float64, n_genes)
    gamma = zeros(Float64, n_genes)
    velocity = zeros(Float64, n_genes, n_cells)

    for gene_index in 1:n_genes
        gene_alpha, gene_beta, gene_gamma, gene_velocity = _kinetic_velocity_fit(view(spliced_dense, gene_index, :), view(unspliced_dense, gene_index, :), fitted_latent_time)
        alpha[gene_index] = gene_alpha
        beta[gene_index] = gene_beta
        gamma[gene_index] = gene_gamma
        velocity[gene_index, :] .= gene_velocity
    end

    cell_scores = vec(mean(velocity, dims=1))
    if latent_time === nothing
        fitted_latent_time = base_result.latent_time
    end
    fallback_used = latent_time === nothing
    return DynamicalRNAVelocityResult(
        copy(experiment.gene_ids),
        copy(experiment.cell_ids),
        alpha,
        beta,
        gamma,
        velocity,
        cell_scores,
        fitted_latent_time,
        provenance_record(
            "DynamicalRNAVelocityResult",
            "SingleCell/calculate_dynamical_rna_velocity";
            notes=["per-gene kinetic ODE fit over latent time ordering"],
            fallbacks=fallback_used ? ["latent time inferred from calculate_rna_velocity"] : String[],
            parameters=(gene_count=n_genes, cell_count=n_cells, normalize=Bool(normalize), use_cuda=Bool(use_cuda), used_supplied_latent_time=!fallback_used)))
end

function _wnn_weight_profile(embedding::AbstractMatrix{<:Real}, neighbors::Vector{Vector{Int}})
    weights = zeros(Float64, size(embedding, 1))
    for cell_index in 1:size(embedding, 1)
        neighborhood = neighbors[cell_index]
        isempty(neighborhood) && continue
        distances = [norm(view(embedding, cell_index, :) .- view(embedding, neighbor, :)) for neighbor in neighborhood]
        weights[cell_index] = 1 / (mean(distances) + eps(Float64))
    end
    return weights
end

"""
    WNNResult

Weighted nearest-neighbor integration across multiple modalities.
"""
struct WNNResult <: AbstractAnalysisResult
    cell_ids::Vector{String}
    modalities::Vector{String}
    modality_embeddings::Dict{String,Matrix{Float64}}
    modality_neighbors::Dict{String,Vector{Vector{Int}}}
    modality_weights::Dict{String,Vector{Float64}}
    combined_neighbors::Vector{Vector{Int}}
    provenance::ResultProvenance
end

WNNResult(modalities, modality_embeddings, modality_neighbors, modality_weights, combined_neighbors) =
    WNNResult(modalities, modality_embeddings, modality_neighbors, modality_weights, combined_neighbors, provenance_record("WNNResult", "singlecell"))

WNNResult(cell_ids::Vector{String}, modalities::Vector{String}, modality_embeddings::Dict{String,Matrix{Float64}}, modality_neighbors::Dict{String,Vector{Vector{Int}}}, modality_weights::Dict{String,Vector{Float64}}, combined_neighbors::Vector{Vector{Int}}) =
    WNNResult(cell_ids, modalities, modality_embeddings, modality_neighbors, modality_weights, combined_neighbors, provenance_record("WNNResult", "SingleCell/weighted_nearest_neighbors"))

function Base.show(io::IO, result::WNNResult)
    println(io, analysis_result_summary(result))
    println(io, "  combined neighbor lists: ", length(result.combined_neighbors))
end

function _show_wnn_result(io::IO, result::WNNResult)
    println(io, analysis_result_summary(result))
    println(io, "  combined neighbor lists: ", length(result.combined_neighbors))
end

function Base.show(io::IO, ::MIME"text/plain", result::WNNResult)
    _show_wnn_result(io, result)
end

function weighted_nearest_neighbors(modalities::AbstractDict{<:String,<:SingleCellExperiment}; k::Int=15, n_components::Int=20, use_cuda::Bool=false)
    isempty(modalities) && throw(ArgumentError("at least one modality is required"))
    length(modalities) >= 2 || throw(ArgumentError("at least two modalities are required"))
    k >= 0 || throw(ArgumentError("k must be nonnegative"))

    modality_lookup = Dict{String,SingleCellExperiment}(String(name) => experiment for (name, experiment) in modalities)
    modality_names = sort!(collect(keys(modality_lookup)))
    reference_cells = copy(modality_lookup[modality_names[1]].cell_ids)
    for experiment in values(modality_lookup)
        experiment.cell_ids == reference_cells || throw(ArgumentError("all modalities must use the same cell ordering"))
    end

    modality_embeddings = Dict{String,Matrix{Float64}}()
    modality_neighbors = Dict{String,Vector{Vector{Int}}}()
    modality_weights = Dict{String,Vector{Float64}}()
    n_cells = length(reference_cells)

    for modality in modality_names
        experiment = modality_lookup[modality]
        normalized = normalize_counts(experiment; use_cuda=use_cuda)
        embedding = run_pca(experiment; normalized=normalized, n_components=n_components, use_cuda=use_cuda)
        neighbors = find_neighbors(experiment; embedding=embedding, k=k, graph_name="wnn_$(modality)", use_cuda=use_cuda)
        modality_embeddings[modality] = Matrix{Float64}(embedding)
        modality_neighbors[modality] = neighbors
        modality_weights[modality] = _wnn_weight_profile(embedding, neighbors)
    end

    combined_neighbors = Vector{Vector{Int}}(undef, n_cells)
    for cell_index in 1:n_cells
        candidate_scores = Dict{Int,Float64}()
        weight_total = sum(modality_weights[modality][cell_index] for modality in modality_names)
        for modality in modality_names
            modality_weight = weight_total > 0 ? modality_weights[modality][cell_index] / weight_total : 1 / length(modality_names)
            neighborhood = modality_neighbors[modality][cell_index]
            neighborhood_size = max(length(neighborhood), 1)
            for (rank, neighbor_index) in enumerate(neighborhood)
                neighbor_index == cell_index && continue
                rank_weight = (neighborhood_size - rank + 1) / neighborhood_size
                candidate_scores[neighbor_index] = get(candidate_scores, neighbor_index, 0.0) + modality_weight * rank_weight
            end
        end
        ranked = sort!(collect(candidate_scores); by = pair -> pair[2], rev=true)
        combined_neighbors[cell_index] = [neighbor for (neighbor, _) in Iterators.take(ranked, min(k, length(ranked)))]
    end
    return WNNResult(
        reference_cells,
        modality_names,
        modality_embeddings,
        modality_neighbors,
        modality_weights,
        combined_neighbors,
        provenance_record(
            "WNNResult",
            "SingleCell/weighted_nearest_neighbors";
            notes=["per-modality PCA, neighbor graph construction, and weighted fusion"],
            parameters=(modality_count=length(modality_names), cell_count=n_cells, k=Int(k), n_components=Int(n_components), use_cuda=Bool(use_cuda))))
end

"""
    PerturbationPredictionResult

Ranked perturbation-response summary for a target gene knockout or activation.
"""
struct PerturbationPredictionResult <: AbstractAnalysisResult
    target_gene::String
    direction::Symbol
    gene_ids::Vector{String}
    scores::Vector{Float64}
    predicted_expression::Vector{Float64}
    ranked_genes::Vector{String}
    ranked_scores::Vector{Float64}
    provenance::ResultProvenance
end

PerturbationPredictionResult(direction, gene_ids, scores, predicted_expression, ranked_genes, ranked_scores) =
    PerturbationPredictionResult(direction, gene_ids, scores, predicted_expression, ranked_genes, ranked_scores, provenance_record("PerturbationPredictionResult", "singlecell"))

PerturbationPredictionResult(target_gene::String, direction::Symbol, gene_ids::Vector{String}, scores::Vector{Float64}, predicted_expression::Vector{Float64}, ranked_genes::Vector{String}, ranked_scores::Vector{Float64}) =
    PerturbationPredictionResult(target_gene, direction, gene_ids, scores, predicted_expression, ranked_genes, ranked_scores, provenance_record("PerturbationPredictionResult", "SingleCell/predict_perturbation"; parameters=(target_gene=target_gene, direction=direction)))

function Base.show(io::IO, result::PerturbationPredictionResult)
    println(io, analysis_result_summary(result))
    for (index, gene) in enumerate(result.ranked_genes[1:min(5, length(result.ranked_genes))])
        println(io, "  ", index, ". ", gene, " score=", round(result.ranked_scores[index]; digits=3))
    end
end

function _show_perturbation_prediction(io::IO, result::PerturbationPredictionResult)
    println(io, analysis_result_summary(result))
    for (index, gene) in enumerate(result.ranked_genes[1:min(5, length(result.ranked_genes))])
        println(io, "  ", index, ". ", gene, " score=", round(result.ranked_scores[index]; digits=3))
    end
end

function Base.show(io::IO, ::MIME"text/plain", result::PerturbationPredictionResult)
    _show_perturbation_prediction(io, result)
end

function predict_perturbation(experiment::SingleCellExperiment, target_gene::String; normalize::Bool=true, top_n::Int=20, diffusion_steps::Int=1, direction::Symbol=:knockout)
    direction in (:knockout, :overexpression) || throw(ArgumentError("direction must be :knockout or :overexpression"))
    gene_lookup = Dict(gene => index for (index, gene) in enumerate(experiment.gene_ids))
    haskey(gene_lookup, String(target_gene)) || throw(ArgumentError("target gene '$target_gene' was not found"))

    expression = normalize ? normalize_counts(experiment) : Array{Float64}(experiment.counts)
    matrix = Matrix{Float64}(expression)
    n_genes, n_cells = size(matrix)
    n_cells >= 2 || throw(ArgumentError("at least two cells are required to predict a perturbation response"))

    centered = matrix .- mean(matrix, dims=2)
    covariance = centered * permutedims(centered) / max(n_cells - 1, 1)
    scales = sqrt.(max.(diag(covariance), 0.0) .+ eps(Float64))
    scale_matrix = max.(scales * permutedims(scales), eps(Float64))
    correlation = covariance ./ scale_matrix
    for index in 1:n_genes
        correlation[index, index] = 1.0
    end

    target_index = gene_lookup[String(target_gene)]
    response = direction == :knockout ? -correlation[:, target_index] : correlation[:, target_index]
    response[target_index] = 0.0

    adjacency = abs.(correlation)
    for index in 1:n_genes
        adjacency[index, index] = 0.0
    end
    row_totals = vec(sum(adjacency, dims=2))
    for index in 1:n_genes
        row_totals[index] > 0 && (adjacency[index, :] ./= row_totals[index])
    end

    propagated = copy(response)
    for _ in 1:max(diffusion_steps, 0)
        propagated = adjacency * propagated
        propagated[target_index] = 0.0
    end

    baseline_expression = vec(mean(matrix, dims=2))
    predicted_expression = baseline_expression .+ propagated
    ranked_indices = sort(collect(1:n_genes); by = index -> abs(propagated[index]), rev=true)
    ranked_indices = ranked_indices[1:min(top_n, length(ranked_indices))]
    ranked_genes = [experiment.gene_ids[index] for index in ranked_indices]
    ranked_scores = [propagated[index] for index in ranked_indices]
    return PerturbationPredictionResult(
        String(target_gene),
        direction,
        copy(experiment.gene_ids),
        propagated,
        predicted_expression,
        ranked_genes,
        ranked_scores,
        provenance_record(
            "PerturbationPredictionResult",
            "SingleCell/predict_perturbation";
            notes=["covariance-derived diffusion propagation from target-gene perturbation"],
            parameters=(target_gene=String(target_gene), direction=direction, normalize=Bool(normalize), top_n=Int(top_n), diffusion_steps=Int(diffusion_steps), gene_count=n_genes, cell_count=n_cells)))
end

"""
    LigandReceptorPair

Minimal ligand-receptor pair definition used by the communication scorer.
"""
struct LigandReceptorPair
    ligand::String
    receptor::String
    pathway::String
    evidence::String
end

LigandReceptorPair(ligand::String, receptor::String, evidence::String) = LigandReceptorPair(String(ligand), String(receptor), "general", String(evidence))

"""
    CellCommunicationResult

One sender/receiver interaction scored from cluster-level ligand and receptor expression.
"""
struct CellCommunicationResult <: AbstractAnalysisResult
    sender::String
    receiver::String
    ligand::String
    receptor::String
    pathway::String
    ligand_expression::Float64
    receptor_expression::Float64
    score::Float64
    evidence::String
    provenance::ResultProvenance
end

CellCommunicationResult(sender, receiver, ligand, receptor, pathway, ligand_expression, receptor_expression, score, evidence) =
    CellCommunicationResult(sender, receiver, ligand, receptor, pathway, ligand_expression, receptor_expression, score, evidence, provenance_record("CellCommunicationResult", "singlecell"))

"""
    CellCommunicationNetwork

Cluster-by-cluster aggregate view of the strongest ligand-receptor interactions.
"""
struct CellCommunicationNetwork
    cluster_ids::Vector{String}
    score_matrix::Matrix{Float64}
    interactions::Vector{CellCommunicationResult}
end

function default_ligand_receptor_pairs()
    return LigandReceptorPair[
        LigandReceptorPair("TGFB1", "TGFBR1", "TGF-beta", "canonical"),
        LigandReceptorPair("TGFB2", "TGFBR1", "TGF-beta", "canonical"),
        LigandReceptorPair("VEGFA", "KDR", "VEGF", "canonical"),
        LigandReceptorPair("VEGFB", "FLT1", "VEGF", "canonical"),
        LigandReceptorPair("CXCL12", "CXCR4", "Chemokine", "canonical"),
        LigandReceptorPair("CCL5", "CCR5", "Chemokine", "canonical"),
        LigandReceptorPair("EGF", "EGFR", "EGF", "canonical"),
        LigandReceptorPair("AREG", "EGFR", "EGF", "canonical"),
        LigandReceptorPair("JAG1", "NOTCH1", "Notch", "canonical"),
        LigandReceptorPair("DLL4", "NOTCH1", "Notch", "canonical"),
        LigandReceptorPair("IL6", "IL6R", "IL6/JAK-STAT", "canonical"),
        LigandReceptorPair("IL11", "IL6ST", "IL6/JAK-STAT", "canonical"),
        LigandReceptorPair("TNF", "TNFRSF1A", "TNF", "canonical"),
        LigandReceptorPair("LTA", "TNFRSF1B", "TNF", "canonical"),
        LigandReceptorPair("WNT3A", "FZD1", "WNT", "canonical"),
        LigandReceptorPair("WNT5A", "FZD2", "WNT", "canonical"),
        LigandReceptorPair("FGF2", "FGFR1", "FGF", "canonical"),
        LigandReceptorPair("FGF7", "FGFR2", "FGF", "canonical"),
        LigandReceptorPair("PDGFB", "PDGFRB", "PDGF", "canonical"),
        LigandReceptorPair("HGF", "MET", "HGF", "canonical"),
        LigandReceptorPair("ANGPT1", "TEK", "Angiopoietin", "canonical"),
        LigandReceptorPair("BMP4", "BMPR1A", "BMP", "canonical"),
        LigandReceptorPair("KITLG", "KIT", "Stem cell factor", "canonical"),
    ]
end

function _cluster_memberships(labels::AbstractVector{<:Integer})
    members = Dict{Int,Vector{Int}}()
    for (index, label) in enumerate(labels)
        push!(get!(members, Int(label), Int[]), index)
    end
    return members
end

function find_cell_communication(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; pairs::AbstractVector{<:LigandReceptorPair}=default_ligand_receptor_pairs(), normalize::Bool=true, min_expression::Real=0.0)
    length(labels) == length(experiment.cell_ids) || throw(ArgumentError("labels must match the number of cells"))
    expression = normalize ? normalize_counts(experiment) : Array{Float64}(experiment.counts)
    gene_lookup = Dict(gene => index for (index, gene) in enumerate(experiment.gene_ids))
    cluster_ids = sort!(unique(Int.(labels)))
    cluster_members = _cluster_memberships(labels)
    results = CellCommunicationResult[]

    for sender in cluster_ids
        sender_cells = get(cluster_members, sender, Int[])
        isempty(sender_cells) && continue
        for receiver in cluster_ids
            receiver_cells = get(cluster_members, receiver, Int[])
            isempty(receiver_cells) && continue
            for pair in pairs
                haskey(gene_lookup, pair.ligand) || continue
                haskey(gene_lookup, pair.receptor) || continue
                ligand_expression = mean(_expression_for_cells(expression, gene_lookup[pair.ligand], sender_cells))
                receptor_expression = mean(_expression_for_cells(expression, gene_lookup[pair.receptor], receiver_cells))
                score = sqrt(max(ligand_expression, 0.0) * max(receptor_expression, 0.0))
                score < min_expression && continue
                push!(results, CellCommunicationResult(string(sender), string(receiver), pair.ligand, pair.receptor, pair.pathway, ligand_expression, receptor_expression, score, pair.evidence))
            end
        end
    end
    return results
end

function communication_network(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; pairs::AbstractVector{<:LigandReceptorPair}=default_ligand_receptor_pairs(), normalize::Bool=true, min_expression::Real=0.0)
    interactions = find_cell_communication(experiment, labels; pairs=pairs, normalize=normalize, min_expression=min_expression)
    cluster_ids = sort!(unique(Int.(labels)))
    cluster_names = string.(cluster_ids)
    index_lookup = Dict(name => index for (index, name) in enumerate(cluster_names))
    score_matrix = zeros(Float64, length(cluster_names), length(cluster_names))
    for interaction in interactions
        sender_index = index_lookup[interaction.sender]
        receiver_index = index_lookup[interaction.receiver]
        score_matrix[sender_index, receiver_index] = max(score_matrix[sender_index, receiver_index], interaction.score)
    end
    return CellCommunicationNetwork(cluster_names, score_matrix, interactions)
end

"""
    LigandReceptorReport

Ranked ligand-receptor interaction report with a filtered pathway view.
"""
struct LigandReceptorReport
    pathway::Union{Nothing,String}
    interactions::Vector{CellCommunicationResult}
end

function _filter_ranked_interactions(interactions::AbstractVector{<:CellCommunicationResult}; pathway::Union{Nothing,String}=nothing, sender::Union{Nothing,String}=nothing, receiver::Union{Nothing,String}=nothing)
    ranked = CellCommunicationResult[]
    for interaction in interactions
        pathway !== nothing && interaction.pathway != String(pathway) && continue
        sender !== nothing && interaction.sender != String(sender) && continue
        receiver !== nothing && interaction.receiver != String(receiver) && continue
        push!(ranked, interaction)
    end
    return sort(ranked; by = interaction -> interaction.score, rev=true)
end

function rank_ligand_receptor_report(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; pathway::Union{Nothing,String}=nothing, sender::Union{Nothing,String}=nothing, receiver::Union{Nothing,String}=nothing, top_n::Int=20, pairs::AbstractVector{<:LigandReceptorPair}=default_ligand_receptor_pairs(), normalize::Bool=true, min_expression::Real=0.0)
    interactions = find_cell_communication(experiment, labels; pairs=pairs, normalize=normalize, min_expression=min_expression)
    ranked = _filter_ranked_interactions(interactions; pathway=pathway, sender=sender, receiver=receiver)
    top_interactions = isempty(ranked) ? CellCommunicationResult[] : ranked[1:min(top_n, length(ranked))]
    return LigandReceptorReport(pathway === nothing ? nothing : String(pathway), top_interactions)
end

function rank_ligand_receptor_report(network::CellCommunicationNetwork; pathway::Union{Nothing,String}=nothing, sender::Union{Nothing,String}=nothing, receiver::Union{Nothing,String}=nothing, top_n::Int=20)
    ranked = _filter_ranked_interactions(network.interactions; pathway=pathway, sender=sender, receiver=receiver)
    top_interactions = isempty(ranked) ? CellCommunicationResult[] : ranked[1:min(top_n, length(ranked))]
    return LigandReceptorReport(pathway === nothing ? nothing : String(pathway), top_interactions)
end

function rank_ligand_receptor_report(interactions::AbstractVector{<:CellCommunicationResult}; pathway::Union{Nothing,String}=nothing, sender::Union{Nothing,String}=nothing, receiver::Union{Nothing,String}=nothing, top_n::Int=20)
    ranked = _filter_ranked_interactions(interactions; pathway=pathway, sender=sender, receiver=receiver)
    top_interactions = isempty(ranked) ? CellCommunicationResult[] : ranked[1:min(top_n, length(ranked))]
    return LigandReceptorReport(pathway === nothing ? nothing : String(pathway), top_interactions)
end

function Base.show(io::IO, report::LigandReceptorReport)
    println(io, "LigandReceptorReport(", report.pathway === nothing ? "all pathways" : report.pathway, ", ", length(report.interactions), " interactions)")
    for (index, interaction) in enumerate(report.interactions[1:min(length(report.interactions), 5)])
        println(io, "  ", index, ". ", interaction.pathway, ": ", interaction.ligand, "->", interaction.receptor, " ", interaction.sender, "->", interaction.receiver, " score=", round(interaction.score; digits=3))
    end
end

function _show_ligand_receptor_report(io::IO, report::LigandReceptorReport)
    println(io, "LigandReceptorReport(", report.pathway === nothing ? "all pathways" : report.pathway, ", ", length(report.interactions), " interactions)")
    for (index, interaction) in enumerate(report.interactions[1:min(length(report.interactions), 5)])
        println(io, "  ", index, ". ", interaction.pathway, ": ", interaction.ligand, "->", interaction.receptor, " ", interaction.sender, "->", interaction.receiver, " score=", round(interaction.score; digits=3))
    end
end

function Base.show(io::IO, ::MIME"text/plain", report::LigandReceptorReport)
    _show_ligand_receptor_report(io, report)
end

"""
    CommunicationPathwaySummary

Aggregate view of the communication score per pathway.
"""
struct CommunicationPathwaySummary
    pathway::String
    interaction_count::Int
    total_score::Float64
    mean_score::Float64
    max_score::Float64
    top_pairs::Vector{String}
end

function communication_pathway_summary(interactions::AbstractVector{<:CellCommunicationResult}; top_n::Int=5)
    pathway_groups = Dict{String,Vector{CellCommunicationResult}}()
    for interaction in interactions
        push!(get!(pathway_groups, interaction.pathway, CellCommunicationResult[]), interaction)
    end

    summaries = CommunicationPathwaySummary[]
    for pathway in sort(collect(keys(pathway_groups)))
        group = pathway_groups[pathway]
        scores = [interaction.score for interaction in group]
        ranked = sort(group; by = interaction -> interaction.score, rev=true)
        top_pairs = [string(interaction.ligand, "->", interaction.receptor, " (", interaction.sender, "->", interaction.receiver, ")") for interaction in first(ranked, min(top_n, length(ranked)))]
        push!(summaries, CommunicationPathwaySummary(pathway, length(group), sum(scores), mean(scores), maximum(scores), top_pairs))
    end
    return summaries
end

function communication_pathway_summary(network::CellCommunicationNetwork; top_n::Int=5)
    return communication_pathway_summary(network.interactions; top_n=top_n)
end

function communication_pathway_summary(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; pairs::AbstractVector{<:LigandReceptorPair}=default_ligand_receptor_pairs(), normalize::Bool=true, min_expression::Real=0.0, top_n::Int=5)
    return communication_pathway_summary(find_cell_communication(experiment, labels; pairs=pairs, normalize=normalize, min_expression=min_expression); top_n=top_n)
end

function _embedding_for_plot(experiment::SingleCellExperiment; embedding::Union{Nothing,AbstractMatrix}=nothing, reduction::String="umap", k::Int=15)
    matrix = if embedding === nothing
        if haskey(experiment.reductions, reduction)
            experiment.reductions[reduction]
        elseif reduction == "umap"
            run_umap(experiment; n_neighbors=max(k, 5))
        else
            run_pca(experiment; n_components=max(2, k))
        end
    else
        Matrix{Float64}(embedding)
    end
    matrix = Matrix{Float64}(matrix)
    if size(matrix, 2) >= 2
        return matrix[:, 1:2]
    elseif size(matrix, 2) == 1
        return hcat(matrix[:, 1], zeros(Float64, size(matrix, 1)))
    else
        return zeros(Float64, size(matrix, 1), 2)
    end
end

function _velocity_quiver_vectors(coords::AbstractMatrix{<:Real}, latent_time::AbstractVector{<:Real}; k::Int=5, scale::Real=0.4)
    n_cells = size(coords, 1)
    u = zeros(Float64, n_cells)
    v = zeros(Float64, n_cells)
    for source in 1:n_cells
        source_time = Float64(latent_time[source])
        isfinite(source_time) || continue
        ranked = [(target, norm(view(coords, source, :) .- view(coords, target, :)), Float64(latent_time[target]) - source_time) for target in 1:n_cells if target != source && isfinite(Float64(latent_time[target]))]
        isempty(ranked) && continue
        sort!(ranked; by = item -> item[2])
        direction = zeros(Float64, 2)
        weight_total = 0.0
        for (target, distance, delta_time) in Iterators.take(ranked, min(k, length(ranked)))
            delta_time <= 0 && continue
            weight = delta_time / (distance + eps(Float64))
            direction[1] += (coords[target, 1] - coords[source, 1]) * weight
            direction[2] += (coords[target, 2] - coords[source, 2]) * weight
            weight_total += weight
        end
        if weight_total > 0
            u[source] = scale * direction[1] / weight_total
            v[source] = scale * direction[2] / weight_total
        end
    end
    return u, v
end

function plot_trajectory(experiment::SingleCellExperiment; root_cell=nothing, embedding::Union{Nothing,AbstractMatrix}=nothing, reduction::String="umap", graph_name::String="neighbors", k::Int=15, pseudotime_result=nothing, title::String="Trajectory", show_edges::Bool=true, show_root::Bool=true, kwargs...)
    result = pseudotime_result === nothing ? calculate_pseudotime(experiment; root_cell=root_cell, embedding=embedding, reduction=reduction, graph_name=graph_name, k=k) : pseudotime_result
    coords = _embedding_for_plot(experiment; embedding=embedding, reduction=reduction, k=k)
    colors = result.pseudotime
    plt = scatter(coords[:, 1], coords[:, 2]; marker_z=colors, color=:viridis, legend=false, title=title, xlabel="Dim 1", ylabel="Dim 2", colorbar_title="Pseudotime", kwargs...)
    if show_edges
        for edge in edges(result.mst)
            source = src(edge)
            target = dst(edge)
            plot!(plt, [coords[source, 1], coords[target, 1]], [coords[source, 2], coords[target, 2]]; color=:gray, alpha=0.3, label=false)
        end
    end
    if show_root
        scatter!(plt, [coords[result.root_index, 1]], [coords[result.root_index, 2]]; marker=:star5, color=:red, markersize=8, label="root")
    end
    return plt
end

function plot_rna_velocity(experiment::SingleCellExperiment; velocity_result=nothing, embedding::Union{Nothing,AbstractMatrix}=nothing, reduction::String="umap", title::String="RNA velocity", root_cell=nothing, k::Int=15, show_root::Bool=true, color::Symbol=:latent_time, show_arrows::Bool=true, arrow_scale::Real=0.4, arrow_neighbors::Int=5, kwargs...)
    result = velocity_result === nothing ? calculate_rna_velocity(experiment) : velocity_result
    coords = _embedding_for_plot(experiment; embedding=embedding, reduction=reduction, k=k)
    color_values = color == :cell_scores ? result.cell_scores : result.latent_time
    velocity_magnitude = vec(mean(abs.(result.velocity), dims=1))
    size_range = maximum(velocity_magnitude) - minimum(velocity_magnitude)
    marker_sizes = size_range > 0 ? 4 .+ 8 .* ((velocity_magnitude .- minimum(velocity_magnitude)) ./ size_range) : fill(6.0, length(velocity_magnitude))
    plt = scatter(coords[:, 1], coords[:, 2]; marker_z=color_values, markersize=marker_sizes, color=:plasma, legend=false, title=title, xlabel="Dim 1", ylabel="Dim 2", colorbar_title=String(color), kwargs...)
    if show_arrows
        u, v = _velocity_quiver_vectors(coords, result.latent_time; k=arrow_neighbors, scale=arrow_scale)
        quiver!(plt, coords[:, 1], coords[:, 2], quiver=(u, v); color=:black, alpha=0.45, linewidth=1, label=false)
    end
    if show_root && root_cell !== nothing
        root_index = _resolve_root_cell(experiment, root_cell)
        scatter!(plt, [coords[root_index, 1]], [coords[root_index, 2]]; marker=:star5, color=:white, markerstrokecolor=:black, markersize=8, label="root")
    end
    return plt
end

plot_rna_velocity_quiver(experiment::SingleCellExperiment; kwargs...) = plot_rna_velocity(experiment; show_arrows=true, kwargs...)

function plot_communication_network(network::CellCommunicationNetwork; title::String="Cell-cell communication network", kwargs...)
    return heatmap(network.cluster_ids, network.cluster_ids, network.score_matrix; xlabel="Receiver", ylabel="Sender", title=title, color=:magma, colorbar_title="Score", kwargs...)
end

function plot_communication_pathways(summary::AbstractVector{<:CommunicationPathwaySummary}; top_n::Int=10, title::String="Communication pathways", kwargs...)
    ranked = sort(collect(summary); by = item -> item.total_score, rev=true)
    isempty(ranked) && return bar(String[], Float64[]; title=title, legend=false, xlabel="Pathway", ylabel="Total score", kwargs...)
    top = ranked[1:min(top_n, length(ranked))]
    labels = [item.pathway for item in top]
    scores = [item.total_score for item in top]
    return bar(labels, scores; title=title, legend=false, xlabel="Pathway", ylabel="Total score", color=:steelblue, kwargs...)
end

plot_communication_pathways(network::CellCommunicationNetwork; top_n::Int=10, title::String="Communication pathways", kwargs...) = plot_communication_pathways(communication_pathway_summary(network); top_n=top_n, title=title, kwargs...)

plot_communication_pathways(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; top_n::Int=10, title::String="Communication pathways", kwargs...) = plot_communication_pathways(communication_pathway_summary(experiment, labels); top_n=top_n, title=title, kwargs...)

function plot_ligand_receptor_report(report::LigandReceptorReport; title::String="Ligand-receptor report", top_n::Int=10, kwargs...)
    ranked = report.interactions[1:min(top_n, length(report.interactions))]
    isempty(ranked) && return bar(String[], Float64[]; title=title, legend=false, xlabel="Interaction", ylabel="Score", kwargs...)
    labels = [string(interaction.pathway, ":", interaction.ligand, "->", interaction.receptor, " ", interaction.sender, "->", interaction.receiver) for interaction in ranked]
    scores = [interaction.score for interaction in ranked]
    return bar(labels, scores; title=title, legend=false, xlabel="Interaction", ylabel="Score", color=:darkorange, xrotation=45, kwargs...)
end

plot_ligand_receptor_report(network::CellCommunicationNetwork; kwargs...) = plot_ligand_receptor_report(rank_ligand_receptor_report(network); kwargs...)

plot_ligand_receptor_report(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; kwargs...) = plot_ligand_receptor_report(rank_ligand_receptor_report(experiment, labels); kwargs...)

"""
    score_cell_cycle(experiment; s_genes=..., g2m_genes=...)

Score each cell for S-phase and G2/M-phase gene signatures and assign a coarse
cell-cycle phase label.
"""
function score_cell_cycle(experiment::SingleCellExperiment; s_genes::AbstractVector{<:String}=["MCM5", "PCNA", "TYMS", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7"], g2m_genes::AbstractVector{<:String}=["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "MKI67"])
    expression = normalize_counts(experiment)
    gene_lookup = Dict(gene => index for (index, gene) in pairs(experiment.gene_ids))
    s_indices = [gene_lookup[String(gene)] for gene in s_genes if haskey(gene_lookup, String(gene))]
    g2m_indices = [gene_lookup[String(gene)] for gene in g2m_genes if haskey(gene_lookup, String(gene))]

    n_cells = size(expression, 2)
    s_scores = zeros(Float64, n_cells)
    g2m_scores = zeros(Float64, n_cells)

    for cell_index in 1:n_cells
        s_scores[cell_index] = isempty(s_indices) ? 0.0 : sum(expression[row, cell_index] for row in s_indices) / length(s_indices)
        g2m_scores[cell_index] = isempty(g2m_indices) ? 0.0 : sum(expression[row, cell_index] for row in g2m_indices) / length(g2m_indices)
    end

    phase = Vector{String}(undef, n_cells)
    for cell_index in 1:n_cells
        delta = s_scores[cell_index] - g2m_scores[cell_index]
        phase[cell_index] = delta > 0.1 ? "S" : delta < -0.1 ? "G2M" : "G1"
    end
    return (s_score=s_scores, g2m_score=g2m_scores, phase=phase)
end

@inline _spatial_bh(p::Vector{Float64}) = benjamini_hochberg(clamp.(p, 0.0, 1.0))

function _spatial_zscore_columns(x::AbstractMatrix{<:Real})
    data = Matrix{Float64}(x)
    for j in axes(data, 2)
        col = @view data[:, j]
        μ = mean(col)
        σ = std(col)
        if isfinite(σ) && σ > 0
            col .= (col .- μ) ./ σ
        else
            col .= 0.0
        end
    end
    return data
end

function _spatial_kmeans_lloyd(x::AbstractMatrix{<:Real}, k::Int; max_iter::Int=50, seed::Int=1)
    n = size(x, 1)
    n >= k || throw(ArgumentError("k must be <= number of rows"))
    rng = MersenneTwister(seed)
    centers = Matrix{Float64}(x[rand(rng, 1:n, k), :])
    labels = ones(Int, n)

    for _ in 1:max_iter
        changed = false
        for i in 1:n
            xi = @view x[i, :]
            best = 1
            bestd = Inf
            for c in 1:k
                d = sum(abs2, xi .- @view(centers[c, :]))
                if d < bestd
                    bestd = d
                    best = c
                end
            end
            if labels[i] != best
                labels[i] = best
                changed = true
            end
        end

        for c in 1:k
            idx = findall(==(c), labels)
            if isempty(idx)
                centers[c, :] .= x[rand(rng, 1:n), :]
            else
                centers[c, :] .= vec(mean(x[idx, :], dims=1))
            end
        end
        !changed && break
    end

    return labels, centers
end

"""
    bayesspace_like_domains(experiment; n_domains=7, n_pcs=15, spatial_weight=0.5)

BayesSpace-inspired spatial domain detection using joint expression and spatial
coordinates.
"""
function bayesspace_like_domains(experiment::SingleCellExperiment; n_domains::Int=7, n_pcs::Int=15, spatial_weight::Real=0.5, seed::Int=1)
    experiment.spatial_coords === nothing && throw(ArgumentError("spatial coordinates are required"))
    pca = run_pca(experiment; n_components=max(2, n_pcs), use_variable_features=true)
    pcs = pca[:, 1:min(size(pca, 2), max(2, n_pcs))]
    coords = _spatial_zscore_columns(experiment.spatial_coords)
    joint = hcat(_spatial_zscore_columns(pcs), Float64(spatial_weight) .* coords)
    labels, centers = _spatial_kmeans_lloyd(joint, n_domains; seed=seed)

    table = DataFrame(
        cell_id=copy(experiment.cell_ids),
        domain=labels,
        x=Float64.(experiment.spatial_coords[:, 1]),
        y=Float64.(experiment.spatial_coords[:, 2]))
    return (assignments=table, centers=centers)
end

"""
    spagc_like_domains(experiment; n_domains=7, n_pcs=15)

SpaGC-inspired spatial domain refinement by combining BayesSpace-like
initialization with MRF smoothing over spatial neighbors.
"""
function spagc_like_domains(experiment::SingleCellExperiment; n_domains::Int=7, n_pcs::Int=15, spatial_weight::Real=0.5, beta::Real=1.0, n_iter::Int=6, k::Int=6, seed::Int=1, threaded::Bool=true)
    base = bayesspace_like_domains(experiment; n_domains=n_domains, n_pcs=n_pcs, spatial_weight=spatial_weight, seed=seed)
    initial = Int.(base.assignments.domain)
    refined = spatial_markov_refine_domains(experiment, initial; beta=beta, n_iter=n_iter, k=k, threaded=threaded)
    table = select(base.assignments, :cell_id, :x, :y)
    table[!, :initial_domain] = initial
    table[!, :spagc_domain] = Int.(refined.refined_domain)
    return table
end

"""
    sparkx_spatial_de(experiment; k=6, n_permutations=200)

SPARK-X-style spatially variable gene screening with permutation p-values.
"""
function sparkx_spatial_de(experiment::SingleCellExperiment; k::Int=6, n_permutations::Int=200, seed::Int=1)
    experiment.spatial_coords === nothing && throw(ArgumentError("spatial coordinates are required"))
    expr = Matrix{Float64}(normalize_counts(experiment))
    neighbors = find_spatial_neighbors(experiment; k=k)
    n_genes, n_cells = size(expr)

    rng = MersenneTwister(seed)
    stat = zeros(Float64, n_genes)
    pval = ones(Float64, n_genes)

    function moran_like(x::AbstractVector{<:Real})
        x0 = Float64.(x) .- mean(x)
        denom = sum(abs2, x0)
        denom <= eps(Float64) && return 0.0
        acc = 0.0
        wsum = 0.0
        for i in 1:n_cells
            for j in neighbors[i]
                acc += x0[i] * x0[j]
                wsum += 1.0
            end
        end
        wsum <= eps(Float64) && return 0.0
        return (n_cells / wsum) * (acc / denom)
    end

    for g in 1:n_genes
        observed = moran_like(vec(@view expr[g, :]))
        stat[g] = observed
        if n_permutations > 0
            exceed = 0
            x = vec(@view expr[g, :])
            for _ in 1:n_permutations
                if moran_like(x[randperm(rng, n_cells)]) >= observed
                    exceed += 1
                end
            end
            pval[g] = (exceed + 1) / (n_permutations + 1)
        else
            z = observed * sqrt(n_cells)
            pval[g] = 2 * ccdf(Normal(), abs(z))
        end
    end
    return DataFrame(gene_id=copy(experiment.gene_ids), spatial_stat=stat, pvalue=pval, padj=_spatial_bh(pval))
end

"""
    spatialde_gp_de(experiment; length_scale=1.0, k=6)

SpatialDE-inspired spatially variable gene scoring using an RBF smoother over
spatial coordinates.
"""
function spatialde_gp_de(experiment::SingleCellExperiment; length_scale::Real=1.0, k::Int=6)
    coords = experiment.spatial_coords
    coords === nothing && throw(ArgumentError("spatial coordinates are required"))
    ℓ = max(Float64(length_scale), eps(Float64))
    expr = Matrix{Float64}(normalize_counts(experiment))
    n_genes, n_cells = size(expr)

    d2 = zeros(Float64, n_cells, n_cells)
    for i in 1:n_cells
        for j in i:n_cells
            ci = @view coords[i, :]
            cj = @view coords[j, :]
            dsq = sum(abs2, ci .- cj)
            d2[i, j] = dsq
            d2[j, i] = dsq
        end
    end
    K = exp.(-d2 ./ (2 * ℓ^2))
    K[diagind(K)] .+= 1e-6
    row_sum = vec(sum(K, dims=2))
    for i in 1:n_cells
        s = row_sum[i] > 0 ? row_sum[i] : 1.0
        K[i, :] ./= s
    end

    gp_stat = zeros(Float64, n_genes)
    pval = ones(Float64, n_genes)
    for g in 1:n_genes
        x = vec(@view expr[g, :])
        μ = mean(x)
        xc = x .- μ
        total_var = var(xc)
        if !(isfinite(total_var) && total_var > 0)
            gp_stat[g] = 0.0
            pval[g] = 1.0
            continue
        end
        smooth = K * xc
        explained = var(smooth)
        frac = clamp(explained / total_var, 0.0, 1.0)
        gp_stat[g] = frac
        z = sqrt(n_cells) * frac
        pval[g] = erfc(abs(z) / sqrt(2.0))
    end
    return DataFrame(gene_id=copy(experiment.gene_ids), gp_stat=gp_stat, pvalue=pval, padj=benjamini_hochberg(pval))
end

"""
    spatial_trajectory_graph(experiment; spatial_weight=0.3)

Spatially informed pseudotime by combining PCA embedding and coordinates.
"""
function spatial_trajectory_graph(experiment::SingleCellExperiment; n_components::Int=15, k::Int=15, spatial_weight::Real=0.3, root_cell=nothing)
    coords = experiment.spatial_coords
    coords === nothing && throw(ArgumentError("spatial coordinates are required"))
    embedding = run_pca(experiment; n_components=max(2, n_components))
    joint = hcat(_spatial_zscore_columns(embedding), Float64(spatial_weight) .* _spatial_zscore_columns(coords))
    return calculate_pseudotime(experiment; embedding=joint, k=k, root_cell=root_cell, graph_name="spatial_trajectory")
end

"""
    niche_weighted_communication(experiment, labels; distance_scale=1.0)

Spatially weighted ligand-receptor communication scoring.
"""
function niche_weighted_communication(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; distance_scale::Real=1.0)
    coords = experiment.spatial_coords
    coords === nothing && throw(ArgumentError("spatial coordinates are required"))
    length(labels) == length(experiment.cell_ids) || throw(DimensionMismatch("labels length must match number of cells"))

    interactions = find_cell_communication(experiment, labels)
    clusters = sort!(unique(Int.(labels)))
    centroids = Dict{Int,Tuple{Float64,Float64}}()
    for c in clusters
        idx = findall(==(c), Int.(labels))
        xy = coords[idx, :]
        centroids[c] = (mean(xy[:, 1]), mean(xy[:, 2]))
    end

    table = DataFrame(
        sender=String[],
        receiver=String[],
        ligand=String[],
        receptor=String[],
        pathway=String[],
        score=Float64[],
        distance=Float64[],
        spatial_weight=Float64[],
        weighted_score=Float64[])

    scale = max(Float64(distance_scale), eps(Float64))
    for row in interactions
        s = parse(Int, row.sender)
        r = parse(Int, row.receiver)
        sx, sy = centroids[s]
        rx, ry = centroids[r]
        d = hypot(sx - rx, sy - ry)
        w = exp(-d / scale)
        push!(table, (row.sender, row.receiver, row.ligand, row.receptor, row.pathway, row.score, d, w, row.score * w))
    end

    sort!(table, :weighted_score, rev=true)
    return table
end

"""
    cell2location_like_segmentation(spot_expression, reference_signatures)

Cell2location-inspired nonnegative deconvolution of spot profiles.
"""
function cell2location_like_segmentation(spot_expression::AbstractMatrix{<:Real}, reference_signatures::AbstractMatrix{<:Real}; cell_type_names=nothing)
    S = Matrix{Float64}(spot_expression)
    R = Matrix{Float64}(reference_signatures)
    size(S, 2) == size(R, 2) || throw(DimensionMismatch("spot and reference matrices must share genes in columns"))

    A = permutedims(R)
    n_spots = size(S, 1)
    n_types = size(R, 1)
    abundance = zeros(Float64, n_spots, n_types)

    for i in 1:n_spots
        y = vec(@view S[i, :])
        w = A \ y
        w = max.(w, 0.0)
        total = sum(w)
        abundance[i, :] .= total > 0 ? w ./ total : fill(1 / n_types, n_types)
    end

    names_types = cell_type_names === nothing ? ["celltype_$(i)" for i in 1:n_types] : String.(cell_type_names)
    length(names_types) == n_types || throw(DimensionMismatch("cell_type_names length must match reference rows"))

    table = DataFrame(spot_id=["spot_$(i)" for i in 1:n_spots])
    for j in 1:n_types
        table[!, Symbol(names_types[j])] = abundance[:, j]
    end
    table[!, :dominant_celltype] = [names_types[argmax(@view abundance[i, :])] for i in 1:n_spots]
    return table
end

"""
    mixscape_like_contrast(expression, perturb_labels; control_label="control")

Mixscape-style mean contrast table for perturbation screens.
"""
function mixscape_like_contrast(expression::AbstractMatrix{<:Real}, perturb_labels::AbstractVector{<:AbstractString}; control_label::AbstractString="control")
    X = Matrix{Float64}(expression)
    n_genes, n_cells = size(X)
    length(perturb_labels) == n_cells || throw(DimensionMismatch("perturb_labels must match number of columns"))

    labels = String.(perturb_labels)
    control_idx = findall(==(String(control_label)), labels)
    isempty(control_idx) && throw(ArgumentError("control_label not found"))

    control_mean = vec(mean(X[:, control_idx], dims=2))
    uniq = sort!(unique(labels))

    out = DataFrame(condition=String[], gene_id=String[], delta=Float64[])
    for cond in uniq
        cond == String(control_label) && continue
        idx = findall(==(cond), labels)
        isempty(idx) && continue
        cond_mean = vec(mean(X[:, idx], dims=2))
        delta = cond_mean .- control_mean
        for g in 1:n_genes
            push!(out, (cond, "gene_$(g)", delta[g]))
        end
    end
    return out
end

"""
    mixscape_multibatch_contrast(expression, perturb_labels, batch_labels)

Mixscape-style multi-batch contrast aggregation with batch-level effect pooling.
"""
function mixscape_multibatch_contrast(expression::AbstractMatrix{<:Real}, perturb_labels::AbstractVector{<:AbstractString}, batch_labels::AbstractVector{<:AbstractString}; control_label::AbstractString="control")
    X = Matrix{Float64}(expression)
    n_genes, n_cells = size(X)
    length(perturb_labels) == n_cells || throw(DimensionMismatch("perturb_labels must match expression columns"))
    length(batch_labels) == n_cells || throw(DimensionMismatch("batch_labels must match expression columns"))

    p = String.(perturb_labels)
    b = String.(batch_labels)
    conditions = sort!(setdiff(unique(p), [String(control_label)]))
    batches = sort!(unique(b))

    out = DataFrame(condition=String[], gene_id=String[], effect=Float64[], std_error=Float64[], zscore=Float64[], pvalue=Float64[], n_batches=Int[])
    for cond in conditions
        batch_effects = Vector{Vector{Float64}}()
        for bb in batches
            idx_ctrl = findall(i -> p[i] == String(control_label) && b[i] == bb, 1:n_cells)
            idx_cond = findall(i -> p[i] == cond && b[i] == bb, 1:n_cells)
            isempty(idx_ctrl) && continue
            isempty(idx_cond) && continue
            μc = vec(mean(X[:, idx_ctrl], dims=2))
            μp = vec(mean(X[:, idx_cond], dims=2))
            push!(batch_effects, μp .- μc)
        end
        isempty(batch_effects) && continue

        E = hcat(batch_effects...)
        eff = vec(mean(E, dims=2))
        se = if size(E, 2) > 1
            vec(std(E, dims=2)) ./ sqrt(size(E, 2))
        else
            fill(1.0, n_genes)
        end

        z = eff ./ max.(se, eps(Float64))
        pval = erfc.(abs.(z) ./ sqrt(2.0))
        for g in 1:n_genes
            push!(out, (cond, "gene_$(g)", eff[g], se[g], z[g], pval[g], size(E, 2)))
        end
    end
    out[!, :padj] = benjamini_hochberg(Float64.(out.pvalue))
    sort!(out, [:padj, :pvalue])
    return out
end

"""
    perturbseq_pseudobulk(expression, gene_ids, group_ids)

Aggregate perturb-seq cells into pseudobulk samples per group.
"""
function perturbseq_pseudobulk(expression::AbstractMatrix{<:Real}, gene_ids::AbstractVector{<:AbstractString}, group_ids::AbstractVector{<:AbstractString}; min_cells::Int=5)
    X = Matrix{Float64}(expression)
    n_genes, n_cells = size(X)
    length(gene_ids) == n_genes || throw(DimensionMismatch("gene_ids must match expression rows"))
    length(group_ids) == n_cells || throw(DimensionMismatch("group_ids must match expression columns"))

    groups = String.(group_ids)
    uniq = sort!(unique(groups))
    keep = String[]
    mats = Vector{Vector{Int}}()
    cells_per_group = Int[]

    for g in uniq
        idx = findall(==(g), groups)
        length(idx) < min_cells && continue
        agg = Int.(round.(sum(X[:, idx], dims=2)))
        push!(keep, g)
        push!(mats, vec(agg))
        push!(cells_per_group, length(idx))
    end

    counts = isempty(mats) ? zeros(Int, n_genes, 0) : hcat(mats...)
    meta = DataFrame(sample_id=keep, n_cells=cells_per_group)
    return (counts=counts, gene_ids=String.(gene_ids), sample_ids=keep, metadata=meta)
end

"""
    milo_like_neighborhood_da(experiment, group_labels)

Milo-like neighborhood differential abundance from kNN neighborhoods.
"""
function milo_like_neighborhood_da(experiment::SingleCellExperiment, group_labels::AbstractVector{<:AbstractString}; k::Int=20, reference_group=nothing)
    n_cells = length(experiment.cell_ids)
    length(group_labels) == n_cells || throw(DimensionMismatch("group_labels must match number of cells"))
    grp = String.(group_labels)
    ref = reference_group === nothing ? first(sort!(unique(grp))) : String(reference_group)
    is_case = grp .!= ref

    neigh = find_neighbors(experiment; k=max(1, k), graph_name="milo_like")
    p_global = mean(Float64.(is_case))

    out = DataFrame(neighborhood_id=Int[], center_cell=String[], case_count=Int[], ref_count=Int[], prop_case=Float64[], zscore=Float64[], pvalue=Float64[])
    for i in 1:n_cells
        nb = unique(vcat(i, neigh[i]))
        n = length(nb)
        c = count(j -> is_case[j], nb)
        r = n - c
        p = c / max(n, 1)
        var = max(p_global * (1 - p_global) / max(n, 1), eps(Float64))
        z = (p - p_global) / sqrt(var)
        pv = erfc(abs(z) / sqrt(2.0))
        push!(out, (i, experiment.cell_ids[i], c, r, p, z, pv))
    end
    out[!, :padj] = benjamini_hochberg(Float64.(out.pvalue))
    sort!(out, :padj)
    return out
end

"""
    perturbation_synergy_scores(expression, perturb_a, perturb_b)

Compute gene-wise perturbation synergy: AB - A - B + Control.
"""
function perturbation_synergy_scores(expression::AbstractMatrix{<:Real}, perturb_a::AbstractVector{Bool}, perturb_b::AbstractVector{Bool}; gene_ids=nothing)
    X = Matrix{Float64}(expression)
    n_genes, n_cells = size(X)
    length(perturb_a) == n_cells || throw(DimensionMismatch("perturb_a must match expression columns"))
    length(perturb_b) == n_cells || throw(DimensionMismatch("perturb_b must match expression columns"))

    idx0 = findall(i -> !perturb_a[i] && !perturb_b[i], 1:n_cells)
    idxa = findall(i -> perturb_a[i] && !perturb_b[i], 1:n_cells)
    idxb = findall(i -> !perturb_a[i] && perturb_b[i], 1:n_cells)
    idxab = findall(i -> perturb_a[i] && perturb_b[i], 1:n_cells)
    isempty(idx0) && throw(ArgumentError("control cells are required"))
    isempty(idxa) && throw(ArgumentError("A-only cells are required"))
    isempty(idxb) && throw(ArgumentError("B-only cells are required"))
    isempty(idxab) && throw(ArgumentError("AB cells are required"))

    μ0 = vec(mean(X[:, idx0], dims=2))
    μa = vec(mean(X[:, idxa], dims=2))
    μb = vec(mean(X[:, idxb], dims=2))
    μab = vec(mean(X[:, idxab], dims=2))
    synergy = μab .- μa .- μb .+ μ0

    s0 = vec(std(X[:, idx0], dims=2))
    effect = synergy ./ max.(s0, eps(Float64))
    ids = gene_ids === nothing ? ["gene_$(i)" for i in 1:n_genes] : String.(gene_ids)
    length(ids) == n_genes || throw(DimensionMismatch("gene_ids length must match expression rows"))
    return DataFrame(gene_id=ids, synergy=synergy, effect_size=effect)
end

"""
    spatial_markov_refine_domains(experiment, initial_labels; beta=1.0, n_iter=8, k=6)

Refine initial spatial domains using an MRF-like objective that combines
expression fit and spatial smoothness.
"""
function spatial_markov_refine_domains(experiment::SingleCellExperiment, initial_labels::AbstractVector{<:Integer}; beta::Real=1.0, n_iter::Int=8, k::Int=6, threaded::Bool=true)
    n_cells = length(experiment.cell_ids)
    length(initial_labels) == n_cells || throw(DimensionMismatch("initial_labels must match number of cells"))
    experiment.spatial_coords === nothing && throw(ArgumentError("spatial coordinates required; call attach_spatial_coords! first"))

    X = permutedims(log1p.(Matrix{Float64}(experiment.counts)))
    labels = Int.(initial_labels)
    neigh = find_spatial_neighbors(experiment; k=max(1, k))
    current = copy(labels)

    for _ in 1:max(1, n_iter)
        domains = sort!(unique(current))
        pos = Dict(d => i for (i, d) in enumerate(domains))
        kdom = length(domains)

        centroids = zeros(Float64, kdom, size(X, 2))
        variances = fill(1.0, kdom)
        for d in domains
            idx = findall(==(d), current)
            p = pos[d]
            centroids[p, :] .= vec(mean(X[idx, :], dims=1))
            d2 = [sum(abs2, @view(X[i, :]) .- @view(centroids[p, :])) for i in idx]
            variances[p] = max(isempty(d2) ? 1.0 : mean(d2), eps(Float64))
        end

        proposal = copy(current)
        threaded_foreach(n_cells, i -> begin
            best_label = current[i]
            best_score = -Inf
            for d in domains
                p = pos[d]
                expr = -sum(abs2, @view(X[i, :]) .- @view(centroids[p, :])) / (2 * variances[p])
                same_n = count(j -> current[j] == d, neigh[i])
                score = expr + Float64(beta) * same_n
                if score > best_score
                    best_score = score
                    best_label = d
                end
            end
            proposal[i] = best_label
        end; threaded=threaded)
        changed = any(proposal .!= current)
        current .= proposal
        !changed && break
    end
    return DataFrame(cell_id=experiment.cell_ids, initial_domain=labels, refined_domain=current)
end

"""
    spatial_lr_permutation_test(experiment, labels; n_permutations=200)

Permutation test for spatially weighted ligand-receptor interaction scores.
"""
function spatial_lr_permutation_test(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; n_permutations::Int=200, distance_scale::Real=1.0, seed::Int=1, threaded::Bool=true)
    length(labels) == length(experiment.cell_ids) || throw(DimensionMismatch("labels must match number of cells"))
    n_permutations >= 1 || throw(ArgumentError("n_permutations must be >= 1"))

    observed_tbl = niche_weighted_communication(experiment, Int.(labels); distance_scale=distance_scale)
    observed = combine(groupby(observed_tbl, [:sender, :receiver, :ligand, :receptor, :pathway]), :weighted_score => sum => :observed_score)

    geq = ones(Int, nrow(observed))
    mean_perm = zeros(Float64, nrow(observed))
    label_vec = Int.(labels)
    perm_scores = threaded_map_collect(1:n_permutations; threaded=threaded) do perm_idx
        rng = MersenneTwister(seed + 104729 * perm_idx)
        perm_labels = shuffle(rng, label_vec)
        perm_tbl = niche_weighted_communication(experiment, perm_labels; distance_scale=distance_scale)
        perm = combine(groupby(perm_tbl, [:sender, :receiver, :ligand, :receptor, :pathway]), :weighted_score => sum => :perm_score)
        joined = leftjoin(observed, perm, on=[:sender, :receiver, :ligand, :receptor, :pathway])
        return coalesce.(joined.perm_score, 0.0)
    end

    for scores_any in perm_scores
        scores = Float64.(scores_any)
        geq .+= scores .>= observed.observed_score
        mean_perm .+= scores
    end

    pvals = geq ./ (n_permutations + 1)
    out = copy(observed)
    out[!, :perm_mean] = mean_perm ./ n_permutations
    out[!, :pvalue] = pvals
    out[!, :padj] = benjamini_hochberg(pvals)
    out[!, :enrichment] = out.observed_score ./ max.(out.perm_mean, eps(Float64))
    sort!(out, [:padj, :observed_score])
    return out
end

include("singlecell_interop.jl")

end
