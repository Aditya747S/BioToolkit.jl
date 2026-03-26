module SingleCell

using SparseArrays
using Statistics
using LinearAlgebra
using Random

using ..DifferentialExpression: CountMatrix, DEResult, calc_norm_factors, differential_expression, estimate_dispersions, filter_low_counts, vst

export SingleCellExperiment, count_matrix, normalize_counts, sctransform, run_pca, run_umap, cluster_cells, find_cluster_markers, summarize_clusters, cluster_marker_summary
export detect_doublets, integrate_batches, score_cell_cycle

mutable struct SingleCellExperiment
    counts::SparseMatrixCSC{Int,Int}
    gene_ids::Vector{String}
    cell_ids::Vector{String}
    metadata::Dict{String,Any}
    reductions::Dict{String,Matrix{Float64}}
    clusters::Dict{String,Vector{Int}}
end

function SingleCellExperiment(counts::AbstractMatrix{<:Integer}, gene_ids::AbstractVector{<:AbstractString}, cell_ids::AbstractVector{<:AbstractString}; metadata::AbstractDict=Dict{String,Any}())
    matrix = sparse(Int.(counts))
    size(matrix, 1) == length(gene_ids) || throw(ArgumentError("gene_ids must match the number of rows in counts"))
    size(matrix, 2) == length(cell_ids) || throw(ArgumentError("cell_ids must match the number of columns in counts"))
    return SingleCellExperiment(matrix, String.(gene_ids), String.(cell_ids), Dict{String,Any}(string(key) => value for (key, value) in metadata), Dict{String,Matrix{Float64}}(), Dict{String,Vector{Int}}())
end

function count_matrix(experiment::SingleCellExperiment)
    return CountMatrix(experiment.counts, experiment.gene_ids, experiment.cell_ids)
end

function _dense_counts(experiment::SingleCellExperiment)
    return Matrix{Float64}(experiment.counts)
end

function normalize_counts(experiment::SingleCellExperiment; scale_factor::Real=1e4, log_transform::Bool=true)
    counts = _dense_counts(experiment)
    library_sizes = vec(sum(counts, dims=1))
    normalized = counts ./ reshape(max.(library_sizes, eps(Float64)), 1, :)
    normalized .*= Float64(scale_factor)
    log_transform && (normalized = log1p.(normalized))
    return normalized
end

function sctransform(experiment::SingleCellExperiment; min_total::Integer=10)
    cm = count_matrix(experiment)
    filtered = filter_low_counts(cm; min_total=min_total)
    norm_factors = calc_norm_factors(filtered)
    dispersions = estimate_dispersions(filtered, norm_factors)
    transformed = vst(filtered, dispersions; norm_factors=norm_factors)
    return transformed, filtered.gene_ids, filtered.sample_ids
end

function _center_rows(matrix::AbstractMatrix{<:Real})
    dense = Matrix{Float64}(matrix)
    dense .-= mean(dense, dims=1)
    return dense
end

function run_pca(experiment::SingleCellExperiment; normalized::Union{Nothing,AbstractMatrix}=nothing, n_components::Integer=20)
    matrix = normalized === nothing ? normalize_counts(experiment) : Matrix{Float64}(normalized)
    centered = _center_rows(permutedims(matrix))
    svd_result = svd(centered; full=false)
    components = min(n_components, size(svd_result.U, 2))
    embedding = svd_result.U[:, 1:components] * Diagonal(svd_result.S[1:components])
    experiment.reductions["pca"] = embedding
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

function cluster_cells(experiment::SingleCellExperiment; embedding::Union{Nothing,AbstractMatrix}=nothing, method::Symbol=:knn, k::Int=15, shared_threshold::Int=1, n_clusters::Int=8, random_seed::Int=1)
    data = embedding === nothing ? get(experiment.reductions, "pca", run_pca(experiment; n_components=max(n_clusters, 2))) : Matrix{Float64}(embedding)
    method == :knn || method == :kmeans || method == :graph || throw(ArgumentError("method must be :knn, :graph, or :kmeans"))
    labels = if method == :kmeans
        _kmeans(data, n_clusters; random_seed=random_seed)
    elseif method == :graph
        _graph_components(_mutual_knn_graph(data, k))
    else
        _connected_components(_knn_neighbors(data, k); shared_threshold=shared_threshold)
    end
    experiment.clusters[method == :kmeans ? "kmeans" : method == :graph ? "graph" : "knn"] = labels
    return labels
end

function run_umap(experiment::SingleCellExperiment; embedding::Union{Nothing,AbstractMatrix}=nothing, n_neighbors::Int=15, min_dist::Real=0.3)
    matrix = embedding === nothing ? get(experiment.reductions, "pca", run_pca(experiment; n_components=10)) : Matrix{Float64}(embedding)
    if isdefined(Main, :UMAP)
        try
            umap_result = getfield(Main, :UMAP).umap(matrix; n_neighbors=n_neighbors, min_dist=min_dist)
            experiment.reductions["umap"] = Matrix{Float64}(umap_result)
            return experiment.reductions["umap"]
        catch
        end
    end
    fallback = matrix[:, 1:min(2, size(matrix, 2))]
    experiment.reductions["umap"] = Matrix{Float64}(fallback)
    return experiment.reductions["umap"]
end

function summarize_clusters(labels::AbstractVector{<:Integer})
    counts = Dict{Int,Int}()
    for label in labels
        counts[Int(label)] = get(counts, Int(label), 0) + 1
    end
    return counts
end

function cluster_marker_summary(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; min_total::Integer=10, shrink::Bool=true, top_n::Int=5)
    unique_labels = sort!(unique(Int.(labels)))
    summaries = Dict{Int,Vector{DEResult}}()
    for cluster_id in unique_labels
        markers = find_cluster_markers(experiment, labels; cluster_id=cluster_id, min_total=min_total, shrink=shrink)
        summaries[cluster_id] = first(sort(markers; by = result -> (result.padj, -abs(result.log2_fold_change))), min(top_n, length(markers)))
    end
    return summaries
end

function find_cluster_markers(experiment::SingleCellExperiment, labels::AbstractVector{<:Integer}; cluster_id::Integer, min_total::Integer=10, shrink::Bool=true)
    length(labels) == length(experiment.cell_ids) || throw(ArgumentError("labels must match the number of cells"))
    count(==(cluster_id), labels) == 0 && return DEResult[]
    count(!=(cluster_id), labels) == 0 && return DEResult[]
    design = Symbol[label == cluster_id ? :cluster : :background for label in labels]
    results = differential_expression(count_matrix(experiment), design; min_total=min_total, shrink=shrink)
    return [result for result in results if result.padj < 1.0]
end

function _cell_embedding(experiment::SingleCellExperiment; n_components::Int=20)
    embedding = get(experiment.reductions, "pca", nothing)
    embedding === nothing && (embedding = run_pca(experiment; n_components=n_components))
    return Matrix{Float64}(embedding)
end

function detect_doublets(experiment::SingleCellExperiment; n_simulated::Int=500, k::Int=15, threshold::Union{Nothing,Real}=nothing, n_components::Int=20, random_seed::Int=1)
    embedding = _cell_embedding(experiment; n_components=n_components)
    n_cells = size(embedding, 1)
    n_cells == 0 && return (scores=Float64[], threshold=0.0, doublets=Bool[], simulated_embedding=zeros(Float64, 0, size(embedding, 2)))

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
        nearest = Iterators.take(distances, min(k, length(distances)))
        synthetic_neighbors = count(pair -> pair[1] > n_cells, nearest)
        scores[cell_index] = k == 0 ? 0.0 : synthetic_neighbors / min(k, length(distances))
    end

    cutoff = threshold === nothing ? (isempty(scores) ? 0.0 : max(mean(scores) + std(scores), quantile(scores, 0.90))) : Float64(threshold)
    doublets = scores .>= cutoff
    return (scores=scores, threshold=cutoff, doublets=doublets, simulated_embedding=simulated)
end

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
    return (corrected_matrix=corrected, embedding=embedding, batch_means=batch_means, method=method)
end

function _pca_for_matrix(matrix::AbstractMatrix{<:Real}; n_components::Int=20)
    centered = Matrix{Float64}(matrix)
    centered .-= mean(centered, dims=1)
    svd_result = svd(centered; full=false)
    components = min(n_components, size(svd_result.U, 2))
    return svd_result.U[:, 1:components] * Diagonal(svd_result.S[1:components])
end

function score_cell_cycle(experiment::SingleCellExperiment; s_genes::AbstractVector{<:AbstractString}=["MCM5", "PCNA", "TYMS", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7"], g2m_genes::AbstractVector{<:AbstractString}=["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "MKI67"])
    expression = normalize_counts(experiment)
    gene_lookup = Dict(gene => index for (index, gene) in pairs(experiment.gene_ids))
    s_indices = [gene_lookup[String(gene)] for gene in s_genes if haskey(gene_lookup, String(gene))]
    g2m_indices = [gene_lookup[String(gene)] for gene in g2m_genes if haskey(gene_lookup, String(gene))]

    n_cells = size(expression, 2)
    s_scores = zeros(Float64, n_cells)
    g2m_scores = zeros(Float64, n_cells)

    for cell_index in 1:n_cells
        s_scores[cell_index] = isempty(s_indices) ? 0.0 : mean(expression[s_indices, cell_index])
        g2m_scores[cell_index] = isempty(g2m_indices) ? 0.0 : mean(expression[g2m_indices, cell_index])
    end

    phase = Vector{String}(undef, n_cells)
    for cell_index in 1:n_cells
        delta = s_scores[cell_index] - g2m_scores[cell_index]
        phase[cell_index] = delta > 0.1 ? "S" : delta < -0.1 ? "G2M" : "G1"
    end

    return (s_score=s_scores, g2m_score=g2m_scores, phase=phase)
end

end