# ==============================================================================
# spatial.jl — Spatial transcriptomics: deconvolution, SVGs, spatial graphs,
#              domain clustering, ligand-receptor scoring, and QC helpers.
#
# References:
#   - Cable et al. (2022) Nat Biotechnol 40:461-471  (RCTD)
#   - Lopez et al. (2022) Nat Methods 19:171-178     (Cell2Location)
#   - Svensson et al. (2018) Science 360:987-990     (SpatialDE)
#   - Moran (1950) Biometrika 37:17-23              (Moran's I)
#   - Geary (1954) Incorp. Statistician 5:115-145   (Geary's C)
#   - Efremova et al. (2020) Nat Methods 17:159-162 (CellPhoneDB / LR scores)
# ==============================================================================

module SpatialDeconvolution

using LinearAlgebra
using Optim
using Random
using Statistics
using DataFrames
using SpecialFunctions: erf

using ..SingleCell: SingleCellExperiment, find_markers, normalize_counts
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!

export SpatialExperiment, DeconvolutionResult
export build_reference_matrix, rctd_deconvolution, cell2location_deconvolution

# Spatial analysis exports
export spatially_variable_genes, spatial_autocorrelation
export spatial_neighborhood_graph, spatial_domain_clustering
export ligand_receptor_spatial_score
export build_spot_deconvolution_qc
export spatial_coexpression_modules
export mark_tissue_boundary_spots
export spatial_pseudotime

# ---------------------------------------------------------------------------
# Core data structures
# ---------------------------------------------------------------------------

struct SpatialExperiment{M<:AbstractMatrix{Int}}
    experiment::SingleCellExperiment{M}
    spatial_coords::Matrix{Float64}
end

struct DeconvolutionResult <: AbstractAnalysisResult
    spot_ids::Vector{String}
    cell_type_ids::Vector{String}
    cell_type_fractions::Matrix{Float64}
    residuals::Vector{Float64}
    method::Symbol
    provenance::ResultProvenance
end

DeconvolutionResult(spot_ids, cell_type_ids, cell_type_fractions, residuals, method) =
    DeconvolutionResult(spot_ids, cell_type_ids, cell_type_fractions, residuals, method, provenance_record("DeconvolutionResult", "spatial"))

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

function SpatialExperiment(experiment::SingleCellExperiment{M}, spatial_coords::AbstractMatrix{<:Real}) where {M<:AbstractMatrix{Int}}
    n_spots = size(experiment.counts, 2)
    size(spatial_coords, 1) == n_spots || throw(ArgumentError("spatial coordinate rows must match number of spots"))
    size(spatial_coords, 2) == 2 || throw(ArgumentError("spatial coordinates must have two columns (x, y)"))
    return SpatialExperiment{M}(experiment, Matrix{Float64}(spatial_coords))
end

function SpatialExperiment(experiment::SingleCellExperiment{M}) where {M<:AbstractMatrix{Int}}
    experiment.spatial_coords === nothing && throw(ArgumentError("SingleCellExperiment does not contain spatial coordinates"))
    return SpatialExperiment(experiment, experiment.spatial_coords)
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

@inline function _softmax(theta::AbstractVector{<:Real})
    shifted = theta .- maximum(theta)
    exps = exp.(shifted)
    return exps ./ sum(exps)
end

function _encode_labels(labels::AbstractVector)
    string_labels = string.(labels)
    unique_labels = sort(unique(string_labels))
    label_to_int = Dict(label => index for (index, label) in enumerate(unique_labels))
    encoded = [label_to_int[label] for label in string_labels]
    return encoded, unique_labels
end

function _safe_dense(matrix)
    return Matrix{Float64}(matrix)
end

@inline function _register_spatial_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

"""
    _spatial_weight_matrix(coords; k=6, radius=nothing)

Build a row-normalised spatial weight matrix W (n_spots × n_spots).
If `radius` is provided, connect spots within that Euclidean distance.
Otherwise connect each spot to its `k` nearest neighbours.
"""
function _spatial_weight_matrix(coords::AbstractMatrix{<:Real}; k::Int=6, radius::Union{Nothing,Real}=nothing)
    n = size(coords, 1)
    C = Matrix{Float64}(coords)
    W = zeros(Float64, n, n)

    if radius !== nothing
        r2 = Float64(radius)^2
        for i in 1:n, j in 1:n
            i == j && continue
            d2 = sum(abs2, C[i, :] .- C[j, :])
            d2 <= r2 && (W[i, j] = 1.0)
        end
    else
        kk = min(k, n - 1)
        for i in 1:n
            dists = [sum(abs2, C[i, :] .- C[j, :]) for j in 1:n]
            dists[i] = Inf
            order = sortperm(dists)
            for j in order[1:kk]
                W[i, j] = 1.0
            end
        end
    end

    # Row-normalise
    row_sums = vec(sum(W, dims=2))
    for i in 1:n
        row_sums[i] > 0 && (W[i, :] ./= row_sums[i])
    end
    return W
end

# ---------------------------------------------------------------------------
# Deconvolution (existing)
# ---------------------------------------------------------------------------

function build_reference_matrix(reference::SingleCellExperiment, reference_labels::AbstractVector; marker_top_n::Integer=30, min_total::Integer=5, test::Symbol=:wilcox, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    size(reference.counts, 2) == length(reference_labels) || throw(ArgumentError("reference labels must match number of cells"))

    encoded_labels, label_names = _encode_labels(reference_labels)
    normalized = _safe_dense(normalize_counts(reference; scale_factor=1e4, log_transform=false))
    n_types = length(label_names)

    selected_markers = Set{String}()
    for type_index in 1:n_types
        markers = find_markers(reference, encoded_labels; ident_1=type_index, ident_2=nothing, test=test, min_total=min_total)
        sort!(markers; by=result -> result.pvalue)
        n_take = min(marker_top_n, length(markers))
        for marker in markers[1:n_take]
            push!(selected_markers, marker.gene_id)
        end
    end

    marker_genes = collect(selected_markers)
    isempty(marker_genes) && (marker_genes = copy(reference.gene_ids))

    marker_index = Dict(gene => idx for (idx, gene) in enumerate(reference.gene_ids))
    marker_indices = [marker_index[gene] for gene in marker_genes if haskey(marker_index, gene)]
    isempty(marker_indices) && (marker_indices = collect(1:length(reference.gene_ids)); marker_genes = copy(reference.gene_ids))

    profiles = zeros(Float64, length(marker_indices), n_types)
    for type_index in 1:n_types
        cells = findall(==(type_index), encoded_labels)
        isempty(cells) && continue
        profile = vec(mean(normalized[marker_indices, cells], dims=2))
        profile .+= eps(Float64)
        profiles[:, type_index] = profile ./ sum(profile)
    end

    result = (
        genes = marker_genes,
        profiles = profiles,
        cell_type_ids = label_names)
    return _register_spatial_result!(_ctx, result, "build_reference_matrix"; parents=provenance_parent_ids(reference), parameters=(marker_top_n=marker_top_n, min_total=min_total, test=String(test), gene_count=length(marker_genes), cell_type_count=length(label_names)))
end

function _shared_profile(reference_tuple, experiment::SingleCellExperiment)
    experiment_index = Dict(gene => idx for (idx, gene) in enumerate(experiment.gene_ids))
    shared_genes = String[]
    reference_rows = Int[]
    experiment_rows = Int[]

    for (row, gene) in enumerate(reference_tuple.genes)
        if haskey(experiment_index, gene)
            push!(shared_genes, gene)
            push!(reference_rows, row)
            push!(experiment_rows, experiment_index[gene])
        end
    end

    isempty(shared_genes) && throw(ArgumentError("no overlapping genes between reference profiles and spatial data"))

    profile = reference_tuple.profiles[reference_rows, :]
    profile = profile .+ eps(Float64)
    profile ./= sum(profile, dims=1)

    return shared_genes, experiment_rows, profile
end

function _optimize_spot(y::Vector{Float64}, profile::Matrix{Float64}; alpha::Union{Nothing,Vector{Float64}}=nothing, maxiter::Integer=250)
    n_types = size(profile, 2)
    depth = max(sum(y), 1.0)

    function objective(theta)
        weights = _softmax(theta)
        mixture = profile * weights
        expected = depth .* max.(mixture, eps(Float64))
        value = sum(expected .- y .* log.(expected))
        if alpha !== nothing
            value -= sum((alpha .- 1.0) .* log.(weights .+ eps(Float64)))
        end
        return value
    end

    result = optimize(objective, zeros(Float64, n_types), LBFGS(), Optim.Options(iterations=maxiter, show_trace=false))
    theta = Optim.minimizer(result)
    weights = _softmax(theta)
    fitted = depth .* (profile * weights)
    residual = sqrt(mean((y .- fitted) .^ 2))
    return weights, residual
end

function _deconvolve(spatial::SpatialExperiment, reference::SingleCellExperiment, reference_labels::AbstractVector; marker_top_n::Integer=30, min_total::Integer=5, test::Symbol=:wilcox, maxiter::Integer=250, alpha::Union{Nothing,Vector{Float64}}=nothing, method::Symbol=:rctd, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    reference_tuple = build_reference_matrix(reference, reference_labels; marker_top_n=marker_top_n, min_total=min_total, test=test, _ctx=_ctx)
    _, spatial_rows, profile = _shared_profile(reference_tuple, spatial.experiment)

    spot_counts = _safe_dense(spatial.experiment.counts[spatial_rows, :])
    n_spots = size(spot_counts, 2)
    n_types = size(profile, 2)

    fractions = zeros(Float64, n_spots, n_types)
    residuals = zeros(Float64, n_spots)

    for spot in 1:n_spots
        y = spot_counts[:, spot]
        weights, residual = _optimize_spot(y, profile; alpha=alpha, maxiter=maxiter)
        fractions[spot, :] = weights
        residuals[spot] = residual
    end

    result = DeconvolutionResult(copy(spatial.experiment.cell_ids), copy(reference_tuple.cell_type_ids), fractions, residuals, method)
    return _register_spatial_result!(_ctx, result, "_deconvolve"; parents=provenance_parent_ids(spatial, reference), parameters=(method=String(method), marker_top_n=marker_top_n, min_total=min_total, test=String(test), maxiter=maxiter, cell_count=length(result.spot_ids), cell_type_count=length(result.cell_type_ids)))
end

function rctd_deconvolution(spatial::SpatialExperiment, reference::SingleCellExperiment, reference_labels::AbstractVector; marker_top_n::Integer=30, min_total::Integer=5, test::Symbol=:wilcox, maxiter::Integer=250, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    return _deconvolve(spatial, reference, reference_labels; marker_top_n=marker_top_n, min_total=min_total, test=test, maxiter=maxiter, alpha=nothing, method=:rctd, _ctx=_ctx)
end

function cell2location_deconvolution(spatial::SpatialExperiment, reference::SingleCellExperiment, reference_labels::AbstractVector; marker_top_n::Integer=30, min_total::Integer=5, test::Symbol=:wilcox, maxiter::Integer=250, dirichlet_alpha::Real=2.0, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    n_types = length(unique(string.(reference_labels)))
    alpha = fill(Float64(dirichlet_alpha), n_types)
    return _deconvolve(spatial, reference, reference_labels; marker_top_n=marker_top_n, min_total=min_total, test=test, maxiter=maxiter, alpha=alpha, method=:cell2location, _ctx=_ctx)
end

# ---------------------------------------------------------------------------
# Spatially Variable Genes — Moran's I
# ---------------------------------------------------------------------------

"""
    spatially_variable_genes(spatial; k=6, radius=nothing, top_n=nothing, permutations=99)

Identify spatially variable genes using Moran's I statistic, analogous to
`spatialDE` and `NNSVG` from Bioconductor.

Returns a `DataFrame` sorted by descending Moran's I with columns:
- `gene_id`, `morans_i`, `expected_i`, `variance_i`, `z_score`, `pvalue`.

Permutation-based p-values are computed when `permutations > 0`.
"""
function spatially_variable_genes(
    spatial::SpatialExperiment;
    k::Int=6,
    radius::Union{Nothing,Real}=nothing,
    top_n::Union{Nothing,Int}=nothing,
    permutations::Int=99,
    min_expr_frac::Real=0.0,
    normalize::Bool=true)
    coords = spatial.spatial_coords
    experiment = spatial.experiment
    n_spots, n_genes = size(experiment.counts, 2), size(experiment.counts, 1)
    W = _spatial_weight_matrix(coords; k=k, radius=radius)

    counts_mat = Matrix{Float64}(experiment.counts)
    if normalize
        col_sums = vec(sum(counts_mat, dims=1))
        col_sums .= max.(col_sums, 1.0)
        counts_mat = counts_mat ./ col_sums' .* 1e4
        counts_mat = log1p.(counts_mat)
    end

    n = n_spots
    E_I = -1.0 / (n - 1)
    S0 = sum(W)
    S1 = 0.5 * sum((W .+ W').^2)
    S2 = sum((vec(sum(W, dims=2)) .+ vec(sum(W, dims=1))).^2)
    b2 = 0.0  # kurtosis term (assume normal)
    var_I = (n^2 * S1 - n * S2 + 3 * S0^2) / ((n^2 - 1) * S0^2) - E_I^2

    min_expr = clamp(Float64(min_expr_frac), 0.0, 1.0)
    rows = NamedTuple{(:gene_id,:morans_i,:expected_i,:variance_i,:z_score,:pvalue),Tuple{String,Float64,Float64,Float64,Float64,Float64}}[]
    for g in 1:n_genes
        x = counts_mat[g, :]
        mean(x .> 0) >= min_expr || continue
        xbar = mean(x)
        z = x .- xbar
        denom = sum(z.^2)
        denom < eps(Float64) && continue
        I = (n / S0) * (z' * (W * z)) / denom
        se = sqrt(max(var_I, 0.0))
        z_score = se > 0 ? (I - E_I) / se : 0.0
        pval = 2.0 * (1.0 - min(0.9999, abs(z_score) > 0 ? 0.5 + 0.5 * erf(abs(z_score) / sqrt(2)) : 0.5))
        if permutations > 0
            extreme = 0
            for _ in 1:permutations
                xp = x[randperm(n)]
                zp = xp .- mean(xp)
                denp = sum(zp.^2)
                Ip = denp < eps(Float64) ? 0.0 : (n / S0) * (zp' * (W * zp)) / denp
                Ip >= I && (extreme += 1)
            end
            pval = (extreme + 1) / (permutations + 1)
        end
        push!(rows, (gene_id=experiment.gene_ids[g], morans_i=I, expected_i=E_I, variance_i=var_I, z_score=z_score, pvalue=pval))
    end

    df = DataFrame(rows)
    sort!(df, :morans_i, rev=true)
    top_n !== nothing && (df = first(df, min(top_n, nrow(df))))

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "spatially_variable_genes";
        parameters=(n_spots=n_spots, n_genes=n_genes, k=k, permutations=permutations, min_expr_frac=min_expr, top_n=top_n))
    end
    return df
end

"""
    spatial_autocorrelation(spatial, gene_id; k=6, radius=nothing, permutations=99)

Compute Moran's I and Geary's C for a single gene.
Returns a NamedTuple with both statistics and p-values.
"""
function spatial_autocorrelation(
    spatial::SpatialExperiment,
    gene_id::String;
    k::Int=6,
    radius::Union{Nothing,Real}=nothing,
    permutations::Int=99,
    min_expr_frac::Real=0.0,
    normalize::Bool=true)
    gene_idx = findfirst(==(gene_id), spatial.experiment.gene_ids)
    gene_idx === nothing && throw(ArgumentError("gene_id $gene_id not found"))
    coords = spatial.spatial_coords
    W = _spatial_weight_matrix(coords; k=k, radius=radius)
    n = size(coords, 1)
    x = Float64.(vec(spatial.experiment.counts[gene_idx, :]))
    if normalize
        cs = max(sum(x), 1.0)
        x = log1p.(x ./ cs .* 1e4)
    end
    xbar = mean(x)
    z = x .- xbar
    S0 = sum(W)
    denom = sum(z.^2)
    morans_i = denom < eps(Float64) ? 0.0 : (n / S0) * (z' * (W * z)) / denom
    E_I = -1.0 / (n - 1)
    # Geary's C
    geary_num = (n - 1) * sum(W[i,j] * (x[i] - x[j])^2 for i in 1:n for j in 1:n)
    gearys_c = denom < eps(Float64) ? 1.0 : geary_num / (2 * S0 * denom)

    pval = permutations > 0 ? begin
        extreme = sum(1 for _ in 1:permutations if begin
            xp = x[randperm(n)]
            zp = xp .- mean(xp)
            denp = sum(zp.^2)
            denp < eps(Float64) ? false : (n / S0) * (zp' * (W * zp)) / denp >= morans_i
        end)
        (extreme + 1) / (permutations + 1)
    end : NaN

    result = (gene_id=gene_id, morans_i=morans_i, gearys_c=gearys_c, expected_i=E_I, pvalue=pval, permutations=permutations)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "spatial_autocorrelation";
        parameters=(gene_id=gene_id, k=k, permutations=permutations))
    end
    return result
end

function spatial_autocorrelation(
    spatial::SpatialExperiment;
    k::Int=6,
    radius::Union{Nothing,Real}=nothing,
    min_expr_frac::Real=0.0,
    normalize::Bool=true)
    coords = spatial.spatial_coords
    experiment = spatial.experiment
    W = _spatial_weight_matrix(coords; k=k, radius=radius)
    counts_mat = Matrix{Float64}(experiment.counts)
    if normalize
        col_sums = vec(sum(counts_mat, dims=1))
        col_sums .= max.(col_sums, 1.0)
        counts_mat = log1p.(counts_mat ./ col_sums' .* 1e4)
    end

    n = size(coords, 1)
    S0 = sum(W)
    expected_i = -1.0 / (n - 1)
    min_expr = clamp(Float64(min_expr_frac), 0.0, 1.0)
    rows = NamedTuple{(:gene_id,:morans_i,:gearys_c,:expected_i),Tuple{String,Float64,Float64,Float64}}[]
    for g in axes(counts_mat, 1)
        x = counts_mat[g, :]
        mean(x .> 0) >= min_expr || continue
        z = x .- mean(x)
        denom = sum(abs2, z)
        denom < eps(Float64) && continue
        morans_i = (n / S0) * (z' * (W * z)) / denom
        geary_num = 0.0
        for i in 1:n, j in 1:n
            geary_num += W[i, j] * (x[i] - x[j])^2
        end
        gearys_c = ((n - 1) / (2 * S0)) * geary_num / denom
        push!(rows, (gene_id=experiment.gene_ids[g], morans_i=morans_i, gearys_c=gearys_c, expected_i=expected_i))
    end
    df = DataFrame(rows)
    _ctx = active_provenance_context()
    _ctx !== nothing && register_provenance!(_ctx, "spatial_autocorrelation"; parameters=(n_genes=size(counts_mat, 1), k=k, min_expr_frac=min_expr))
    return df
end

"""
    spatial_neighborhood_graph(spatial; k=6, radius=nothing)

Build the spatial neighborhood graph. Returns edge and adjacency tables.
"""
function spatial_neighborhood_graph(
    spatial::SpatialExperiment;
    k::Int=6,
    radius::Union{Nothing,Real}=nothing)
    coords = spatial.spatial_coords
    n = size(coords, 1)
    W = _spatial_weight_matrix(coords; k=k, radius=radius)
    adjacency = Int.(W .> 0)
    edges = DataFrame(source=Int[], target=Int[], weight=Float64[])
    for i in 1:n, j in 1:n
        W[i, j] > 0 || continue
        push!(edges, (i, j, W[i, j]))
    end
    result = (edges=edges, adjacency=adjacency)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "spatial_neighborhood_graph";
        parameters=(n_spots=n, k=k, radius=something(radius, NaN), edge_count=nrow(edges)))
    end
    return result
end

"""
    spatial_domain_clustering(spatial; k=6, n_clusters=5, method=:kmeans)

Cluster spatial spots into domains using expression + location features.
Returns a vector of integer cluster labels (1-indexed).
"""
function spatial_domain_clustering(
    spatial::SpatialExperiment;
    k::Int=6,
    n_domains::Union{Nothing,Int}=nothing,
    n_clusters::Int=5,
    method::Symbol=:kmeans,
    normalize::Bool=true,
    seed::Union{Nothing,Int}=nothing,
    random_seed::Int=1)
    n_clusters = n_domains === nothing ? n_clusters : Int(n_domains)
    random_seed = seed === nothing ? random_seed : Int(seed)
    counts_mat = Matrix{Float64}(spatial.experiment.counts)
    if normalize
        cs = vec(sum(counts_mat, dims=1))
        cs .= max.(cs, 1.0)
        counts_mat = log1p.(counts_mat ./ cs' .* 1e4)
    end
    coords = spatial.spatial_coords
    # Normalise coords to same scale as expression PCs
    coord_scale = std(counts_mat) / (std(coords) + eps(Float64))
    scaled_coords = coords .* coord_scale
    features = vcat(counts_mat, scaled_coords')  # (genes+2) × spots

    rng = Random.MersenneTwister(random_seed)
    n_spots = size(features, 2)
    nc = min(n_clusters, n_spots)
    # Simple k-means via random init + assignment
    centroids = features[:, randperm(rng, n_spots)[1:nc]]
    labels = zeros(Int, n_spots)
    for _ in 1:50
        for i in 1:n_spots
            dists = [sum(abs2, features[:, i] .- centroids[:, c]) for c in 1:nc]
            labels[i] = argmin(dists)
        end
        new_centroids = zeros(size(features, 1), nc)
        counts_per_cluster = zeros(Int, nc)
        for i in 1:n_spots
            new_centroids[:, labels[i]] .+= features[:, i]
            counts_per_cluster[labels[i]] += 1
        end
        for c in 1:nc
            counts_per_cluster[c] > 0 && (new_centroids[:, c] ./= counts_per_cluster[c])
        end
        centroids = new_centroids
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "spatial_domain_clustering";
        parameters=(n_spots=n_spots, n_clusters=nc, method=method, random_seed=random_seed))
    end
    return DataFrame(spot_id=spatial.experiment.cell_ids, domain=labels)
end

"""
    ligand_receptor_spatial_score(spatial, lr_pairs; k=6, normalize=true)

Score ligand-receptor pairs using spatial proximity.
`lr_pairs` is a vector of `(ligand_gene, receptor_gene)` tuples.
Returns a DataFrame with columns: ligand, receptor, score, mean_ligand, mean_receptor.
"""
function ligand_receptor_spatial_score(
    spatial::SpatialExperiment,
    lr_pairs::AbstractVector;
    k::Int=6,
    normalize::Bool=true)
    experiment = spatial.experiment
    gene_index = Dict(g => i for (i, g) in enumerate(experiment.gene_ids))
    coords = spatial.spatial_coords
    W = _spatial_weight_matrix(coords; k=k)

    counts_mat = Matrix{Float64}(experiment.counts)
    if normalize
        cs = vec(sum(counts_mat, dims=1))
        cs .= max.(cs, 1.0)
        counts_mat = log1p.(counts_mat ./ cs' .* 1e4)
    end

    rows = NamedTuple[]
    for pair in lr_pairs
        lig, rec = String(pair[1]), String(pair[2])
        lig_idx = get(gene_index, lig, nothing)
        rec_idx = get(gene_index, rec, nothing)
        (lig_idx === nothing || rec_idx === nothing) && continue
        L = vec(counts_mat[lig_idx, :])
        R = vec(counts_mat[rec_idx, :])
        # Score: mean of L_i * sum_j(W_ij * R_j) across spots
        score = mean(L .* (W * R))
        push!(rows, (ligand=lig, receptor=rec, lr_pair="$(lig)|$(rec)", score=score, lr_spatial_score=score, mean_ligand=mean(L), mean_receptor=mean(R)))
    end
    df = isempty(rows) ? DataFrame(ligand=String[], receptor=String[], lr_pair=String[], score=Float64[], lr_spatial_score=Float64[], mean_ligand=Float64[], mean_receptor=Float64[]) : DataFrame(rows)
    sort!(df, :score, rev=true)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "ligand_receptor_spatial_score";
        parameters=(n_pairs=length(lr_pairs), k=k, n_scored=nrow(df)))
    end
    return df
end

"""
    build_spot_deconvolution_qc(result::DeconvolutionResult)

Generate QC summary statistics for a deconvolution result.
Returns a DataFrame with per-spot dominant cell type and residuals.
"""
function build_spot_deconvolution_qc(
    result::DeconvolutionResult)
    n_spots = length(result.spot_ids)
    dominant_types = [result.cell_type_ids[argmax(result.cell_type_fractions[i, :])] for i in 1:n_spots]
    max_fractions = [maximum(result.cell_type_fractions[i, :]) for i in 1:n_spots]
    entropy = [begin
        p = result.cell_type_fractions[i, :]
        -sum(v > 0 ? v * log(v) : 0.0 for v in p)
    end for i in 1:n_spots]
    confidence = [fraction >= 0.7 ? "singlet" : fraction >= 0.45 ? "mixed" : "low_confidence" for fraction in max_fractions]
    df = DataFrame(
        spot_id=result.spot_ids,
        dominant_cell_type=dominant_types,
        dominant_fraction=max_fractions,
        confidence=confidence,
        entropy=entropy,
        residual=result.residuals)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "build_spot_deconvolution_qc";
        parameters=(n_spots=n_spots, method=result.method, n_cell_types=length(result.cell_type_ids)))
    end
    return df
end

"""
    spatial_coexpression_modules(spatial; k=6, n_modules=5, normalize=true)

Identify spatially co-expressed gene modules using a spatially-smoothed
correlation matrix and spectral clustering.
Returns a Vector{Int} of module assignments (1-indexed) for each gene.
"""
function spatial_coexpression_modules(
    spatial::SpatialExperiment;
    k::Int=6,
    n_svgs::Union{Nothing,Int}=nothing,
    n_modules::Int=5,
    min_expr_frac::Real=0.0,
    normalize::Bool=true,
    seed::Union{Nothing,Int}=nothing,
    random_seed::Int=1)
    random_seed = seed === nothing ? random_seed : Int(seed)
    experiment = spatial.experiment
    W = _spatial_weight_matrix(spatial.spatial_coords; k=k)
    counts_mat = Matrix{Float64}(experiment.counts)
    if normalize
        cs = vec(sum(counts_mat, dims=1))
        cs .= max.(cs, 1.0)
        counts_mat = log1p.(counts_mat ./ cs' .* 1e4)
    end
    min_expr = clamp(Float64(min_expr_frac), 0.0, 1.0)
    gene_indices = findall(g -> mean(@view(counts_mat[g, :]) .> 0) >= min_expr, axes(counts_mat, 1))
    isempty(gene_indices) && (gene_indices = collect(axes(counts_mat, 1)))
    if n_svgs !== nothing && length(gene_indices) > n_svgs
        variances = [var(@view counts_mat[g, :]) for g in gene_indices]
        order = sortperm(variances, rev=true)[1:Int(n_svgs)]
        gene_indices = gene_indices[order]
    end
    counts_mat = counts_mat[gene_indices, :]

    # Spatially smooth each gene
    smoothed = counts_mat * W'
    # Correlation matrix across genes
    n_genes = size(smoothed, 1)
    C = cor(smoothed')
    replace!(C, NaN => 0.0)
    # Simple k-means on gene correlation vectors
    rng = Random.MersenneTwister(random_seed)
    nm = min(n_modules, n_genes)
    centroids = C[:, randperm(rng, n_genes)[1:nm]]
    module_labels = zeros(Int, n_genes)
    for _ in 1:50
        for i in 1:n_genes
            dists = [sum(abs2, C[:, i] .- centroids[:, c]) for c in 1:nm]
            module_labels[i] = argmin(dists)
        end
        new_c = zeros(n_genes, nm)
        cnt = zeros(Int, nm)
        for i in 1:n_genes
            new_c[:, module_labels[i]] .+= C[:, i]
            cnt[module_labels[i]] += 1
        end
        for c in 1:nm
            cnt[c] > 0 && (new_c[:, c] ./= cnt[c])
        end
        centroids = new_c
    end
    modules = DataFrame(:gene_id => experiment.gene_ids[gene_indices], :module => module_labels)
    module_scores = zeros(Float64, size(counts_mat, 2), nm)
    for module_id in 1:nm
        rows = findall(==(module_id), module_labels)
        isempty(rows) && continue
        module_scores[:, module_id] .= vec(mean(counts_mat[rows, :], dims=1))
    end
    result = (modules=modules, module_scores=module_scores)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "spatial_coexpression_modules";
        parameters=(n_genes=n_genes, n_modules=nm, k=k, n_svgs=n_svgs === nothing ? n_genes : Int(n_svgs), min_expr_frac=min_expr, random_seed=random_seed))
    end
    return result
end

"""
    mark_tissue_boundary_spots(spatial; k=6, threshold=0.5)

Mark spots that are on the tissue boundary (i.e., have fewer neighbors
than expected).
"""
function mark_tissue_boundary_spots(
    spatial::SpatialExperiment;
    k::Int=6,
    threshold::Real=0.5)
    W = _spatial_weight_matrix(spatial.spatial_coords; k=k)
    n = size(W, 1)
    neighbor_counts = vec(sum(W .> 0, dims=2))
    median_neighbors = median(neighbor_counts)
    boundary = neighbor_counts .< threshold * median_neighbors
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "mark_tissue_boundary_spots";
        parameters=(n_spots=n, k=k, threshold=Float64(threshold), n_boundary=sum(boundary)))
    end
    return DataFrame(spot_id=spatial.experiment.cell_ids, is_boundary=boundary, neighbor_count=neighbor_counts)
end

"""
    spatial_pseudotime(spatial; root_spot=1, k=6)

Compute a spatial pseudotime ordering of spots using a shortest-path
approach on the spatial neighborhood graph (Dijkstra from the root spot).
Returns a Vector{Float64} of pseudotime values per spot (0-indexed from root).
"""
function spatial_pseudotime(
    spatial::SpatialExperiment;
    root_spot::Int=1,
    k::Int=6,
    n_diffusion_steps::Int=0)
    coords = spatial.spatial_coords
    n = size(coords, 1)
    W = _spatial_weight_matrix(coords; k=k)
    # Dijkstra's algorithm on the weight matrix as distances
    dist = fill(Inf, n)
    dist[root_spot] = 0.0
    visited = falses(n)
    for _ in 1:n
        u = argmin([visited[i] ? Inf : dist[i] for i in 1:n])
        visited[u] = true
        for v in 1:n
            W[u, v] > 0 || continue
            edge_dist = 1.0 / (W[u, v] + eps(Float64))
            dist[u] + edge_dist < dist[v] && (dist[v] = dist[u] + edge_dist)
        end
    end
    # Normalize to [0, 1]
    finite_dist = filter(isfinite, dist)
    if !isempty(finite_dist) && maximum(finite_dist) > 0
        dist ./= maximum(finite_dist)
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "spatial_pseudotime";
        parameters=(n_spots=n, root_spot=root_spot, k=k, n_diffusion_steps=Int(n_diffusion_steps)))
    end
    return DataFrame(spot_id=spatial.experiment.cell_ids, pseudotime=dist)
end

end # module SpatialDeconvolution
