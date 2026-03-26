module Microbiome

using SparseArrays
using DataFrames
using Statistics
using LinearAlgebra
using Random
using Distributions
using Base.Threads
using Graphs
using Optim
using Plots

using ..DifferentialExpression: CountMatrix, benjamini_hochberg
using ..BioToolkit: PhyloTree, get_terminals

export CommunityProfile, PCoAResult, NMDSResult, ANCOMResult, SongbirdResult, SourceTrackingResult, MicrobiomeNetwork
export clr_transform, ilr_transform, bray_curtis, unifrac, weighted_unifrac, shannon_entropy, simpson_index, faith_pd
export pairwise_bray_curtis, pairwise_unifrac, pcoa, pcoa_plot, nmds, ancom, songbird, cooccurrence_network, network_plot, source_tracking_model, source_tracking, source_tracking_posterior_summary

struct CommunityProfile
    counts::CountMatrix
    taxonomy::DataFrame
    tree::PhyloTree
    metadata::DataFrame

    function CommunityProfile(counts::CountMatrix, taxonomy::DataFrame, tree::PhyloTree, metadata::DataFrame)
        _validate_profile(counts, taxonomy, tree, metadata)
        new(counts, copy(taxonomy), tree, copy(metadata))
    end
end

CommunityProfile(counts::AbstractMatrix{<:Integer}, gene_ids::AbstractVector{<:AbstractString}, sample_ids::AbstractVector{<:AbstractString}, taxonomy::DataFrame, tree::PhyloTree, metadata::DataFrame) =
    CommunityProfile(CountMatrix(counts, gene_ids, sample_ids), taxonomy, tree, metadata)

struct PCoAResult
    coordinates::Matrix{Float64}
    eigenvalues::Vector{Float64}
    variance_explained::Vector{Float64}
end

struct NMDSResult
    coordinates::Matrix{Float64}
    stress::Float64
    iterations::Int
    converged::Bool
end

struct ANCOMResult
    taxon_ids::Vector{String}
    w_stat::Vector{Int}
    min_pvalue::Vector{Float64}
    qvalue::Vector{Float64}
    significant::Vector{Bool}
end

struct SongbirdResult
    coefficients::Matrix{Float64}
    taxon_ids::Vector{String}
    feature_names::Vector{String}
    loglik::Float64
    iterations::Int
    converged::Bool
end

struct SourceTrackingResult
    chain
    mean_proportions::Vector{Float64}
    median_proportions::Vector{Float64}
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
end

struct MicrobiomeNetwork
    graph::SimpleGraph
    weights::Dict{Tuple{Int,Int},Float64}
    taxa::Vector{String}
    coordinates::Matrix{Float64}
end

function Base.show(io::IO, profile::CommunityProfile)
    print(io, "CommunityProfile(", length(profile.counts.gene_ids), " taxa, ", length(profile.counts.sample_ids), " samples)")
end

function _validate_profile(counts::CountMatrix, taxonomy::DataFrame, tree::PhyloTree, metadata::DataFrame)
    hasproperty(taxonomy, :taxon_id) || throw(ArgumentError("taxonomy DataFrame must contain a taxon_id column"))
    hasproperty(metadata, :sample_id) || throw(ArgumentError("metadata DataFrame must contain a sample_id column"))

    taxon_ids = String.(taxonomy.taxon_id)
    sample_ids = String.(metadata.sample_id)
    counts.gene_ids == taxon_ids || throw(ArgumentError("taxonomy taxon_id order must match counts gene_ids"))
    counts.sample_ids == sample_ids || throw(ArgumentError("metadata sample_id order must match counts sample_ids"))

    tree_taxa = Set(node.name for node in get_terminals(tree))
    all(taxon -> taxon in tree_taxa, counts.gene_ids) || throw(ArgumentError("tree must contain every taxon in counts"))
    return nothing
end

function _selection_indices(labels::Vector{String}, selector)
    selector === Colon() && return collect(eachindex(labels))
    selector isa Integer && return [Int(selector)]
    selector isa AbstractVector{Bool} && return findall(selector)
    selector isa AbstractVector{<:Integer} && return Int.(selector)
    selector isa AbstractVector{<:AbstractString} && return [findfirst(==(String(value)), labels) for value in selector]
    selector isa AbstractVector{Symbol} && return [findfirst(==(String(value)), labels) for value in selector]
    selector isa AbstractString && return [findfirst(==(String(selector)), labels)]
    selector isa Symbol && return [findfirst(==(String(selector)), labels)]
    return collect(selector)
end

function _subset_dataframe(df::DataFrame, column::Symbol, ids::Vector{String})
    lookup = Dict(String(value) => index for (index, value) in enumerate(df[!, column]))
    row_indices = [lookup[id] for id in ids if haskey(lookup, id)]
    return df[row_indices, :]
end

function _prune_tree(tree::PhyloTree, keep::Set{String})
    if isempty(tree.children)
        return tree.name in keep ? PhyloTree(tree.name; branch_length=tree.branch_length, support=tree.support, metadata=tree.metadata) : nothing
    end

    children = PhyloTree[]
    for child in tree.children
        pruned = _prune_tree(child, keep)
        pruned === nothing || push!(children, pruned)
    end

    isempty(children) && return nothing
    return PhyloTree(children; name=tree.name, branch_length=tree.branch_length, support=tree.support, metadata=tree.metadata)
end

function _postorder_nodes(tree::PhyloTree)
    stack = Tuple{PhyloTree,Bool}[(tree, false)]
    order = PhyloTree[]
    while !isempty(stack)
        node, visited = pop!(stack)
        if visited
            push!(order, node)
        else
            push!(stack, (node, true))
            for child in Iterators.reverse(node.children)
                push!(stack, (child, false))
            end
        end
    end
    return order
end

function _leaf_lookup(tree::PhyloTree)
    leaves = get_terminals(tree)
    lookup = Dict{String,Int}()
    for (index, leaf) in enumerate(leaves)
        lookup[leaf.name] = index
    end
    return leaves, lookup
end

function _count_matrix_dense(counts::CountMatrix)
    return Matrix{Float64}(counts.counts)
end

function _circle_layout(n::Int)
    n > 0 || return zeros(Float64, 0, 2)
    angles = range(0, 2π; length=n + 1)[1:end-1]
    coordinates = zeros(Float64, n, 2)
    for index in 1:n
        coordinates[index, 1] = cos(angles[index])
        coordinates[index, 2] = sin(angles[index])
    end
    return coordinates
end

function _tree_signal_matrix(tree::PhyloTree, counts::CountMatrix; weighted::Bool=false)
    leaves, _ = _leaf_lookup(tree)
    dense = _count_matrix_dense(counts)
    nsamples = size(dense, 2)
    order = _postorder_nodes(tree)
    node_index = Dict{PhyloTree,Int}(node => index for (index, node) in enumerate(order))
    leaf_index = Dict(leaf.name => node_index[leaf] for leaf in leaves)
    signal = zeros(Float64, length(order), nsamples)

    rows, cols, values = findnz(counts.counts)
    for position in eachindex(values)
        taxon = counts.gene_ids[rows[position]]
        leaf_row = get(leaf_index, taxon, 0)
        leaf_row == 0 && continue
        signal[leaf_row, cols[position]] += weighted ? Float64(values[position]) : 1.0
    end

    if weighted
        totals = vec(sum(signal, dims=1))
        for sample in 1:nsamples
            totals[sample] > 0 || continue
            signal[:, sample] ./= totals[sample]
        end
    end

    for node in order
        node_row = node_index[node]
        isempty(node.children) && continue
        for child in node.children
            signal[node_row, :] .+= signal[node_index[child], :]
        end
    end

    return order, signal, node_index
end

function Base.getindex(profile::CommunityProfile, row_sel, col_sel)
    row_idx = _selection_indices(profile.counts.gene_ids, row_sel)
    col_idx = _selection_indices(profile.counts.sample_ids, col_sel)
    new_counts = CountMatrix(profile.counts.counts[row_idx, col_idx], profile.counts.gene_ids[row_idx], profile.counts.sample_ids[col_idx])
    selected_taxa = profile.counts.gene_ids[row_idx]
    selected_samples = profile.counts.sample_ids[col_idx]
    new_taxonomy = _subset_dataframe(profile.taxonomy, :taxon_id, selected_taxa)
    new_metadata = _subset_dataframe(profile.metadata, :sample_id, selected_samples)
    new_tree = _prune_tree(profile.tree, Set(selected_taxa))
    new_tree === nothing && throw(ArgumentError("requested taxa removed the entire tree"))
    return CommunityProfile(new_counts, new_taxonomy, new_tree, new_metadata)
end

function clr_transform(counts::CountMatrix; pseudocount::Real=0.5)
    pseudocount > 0 || throw(ArgumentError("pseudocount must be positive"))
    dense = _count_matrix_dense(counts)
    transformed = similar(dense)
    @threads for sample in axes(dense, 2)
        values = @view dense[:, sample]
        logged = log.(values .+ pseudocount)
        centered = logged .- mean(logged)
        transformed[:, sample] = centered
    end
    return transformed
end

clr_transform(profile::CommunityProfile; pseudocount::Real=0.5) = clr_transform(profile.counts; pseudocount=pseudocount)

function _helmert_basis(size::Int)
    size >= 2 || throw(ArgumentError("Helmert basis requires at least two taxa"))
    basis = zeros(Float64, size, size - 1)
    for column in 1:(size - 1)
        scale = sqrt(column * (column + 1))
        basis[1:column, column] .= 1 / scale
        basis[column + 1, column] = -column / scale
    end
    return basis
end

function ilr_transform(counts::CountMatrix; pseudocount::Real=0.5)
    clr = clr_transform(counts; pseudocount=pseudocount)
    return _helmert_basis(size(clr, 1))' * clr
end

ilr_transform(profile::CommunityProfile; pseudocount::Real=0.5) = ilr_transform(profile.counts; pseudocount=pseudocount)

function bray_curtis(counts::CountMatrix)
    return pairwise_bray_curtis(counts)
end

function pairwise_bray_curtis(counts::CountMatrix)
    dense = _count_matrix_dense(counts)
    nsamples = size(dense, 2)
    distances = zeros(Float64, nsamples, nsamples)
    @threads for left in 1:nsamples
        for right in left+1:nsamples
            left_counts = @view dense[:, left]
            right_counts = @view dense[:, right]
            denominator = sum(left_counts) + sum(right_counts)
            value = denominator == 0 ? 0.0 : sum(abs.(left_counts .- right_counts)) / denominator
            distances[left, right] = value
            distances[right, left] = value
        end
    end
    return distances
end

pairwise_bray_curtis(profile::CommunityProfile) = pairwise_bray_curtis(profile.counts)

bray_curtis(profile::CommunityProfile) = bray_curtis(profile.counts)

function unifrac(profile::CommunityProfile)
    return pairwise_unifrac(profile)
end

function pairwise_unifrac(profile::CommunityProfile; weighted::Bool=false)
    order, signal, node_index = _tree_signal_matrix(profile.tree, profile.counts; weighted=weighted)
    if weighted
        abundances = signal
        nsamples = size(abundances, 2)
        distances = zeros(Float64, nsamples, nsamples)
        total_length = sum(node.branch_length for node in order if node.branch_length > 0)
        total_length == 0 && return distances
        @threads for left in 1:nsamples
            for right in left+1:nsamples
                distance = 0.0
                for node in order
                    branch_length = node.branch_length
                    branch_length == 0 && continue
                    node_row = node_index[node]
                    distance += branch_length * abs(abundances[node_row, left] - abundances[node_row, right])
                end
                distance /= total_length
                distances[left, right] = distance
                distances[right, left] = distance
            end
        end
        return distances
    end

    order, signal, node_index = _tree_signal_matrix(profile.tree, profile.counts; weighted=false)
    presence = signal .> 0
    nsamples = size(presence, 2)
    distances = zeros(Float64, nsamples, nsamples)
    for node in order
        branch_length = node.branch_length
        branch_length == 0 && continue
        node_presence = presence[node_index[node], :]
        @threads for left in 1:nsamples
            for right in left+1:nsamples
                left_present = node_presence[left]
                right_present = node_presence[right]
                if left_present || right_present
                    distances[left, right] += left_present == right_present ? 0.0 : branch_length
                    distances[right, left] = distances[left, right]
                end
            end
        end
    end
    total_length = sum(node.branch_length for node in order if node.branch_length > 0)
    total_length == 0 && return distances
    return distances ./ total_length
end

weighted_unifrac(profile::CommunityProfile) = pairwise_unifrac(profile; weighted=true)

pairwise_unifrac(counts::CountMatrix, tree; weighted::Bool=false) = pairwise_unifrac(CommunityProfile(counts, String.(counts.gene_ids), String.(counts.sample_ids), nothing, tree, nothing); weighted=weighted)

function shannon_entropy(counts::CountMatrix)
    dense = _count_matrix_dense(counts)
    entropy = zeros(Float64, size(dense, 2))
    for sample in axes(dense, 2)
        values = @view dense[:, sample]
        total = sum(values)
        total == 0 && continue
        probabilities = values ./ total
        entropy[sample] = -sum(probabilities[probabilities .> 0] .* log.(probabilities[probabilities .> 0]))
    end
    return entropy
end

shannon_entropy(profile::CommunityProfile) = shannon_entropy(profile.counts)

function simpson_index(counts::CountMatrix)
    dense = _count_matrix_dense(counts)
    index = zeros(Float64, size(dense, 2))
    for sample in axes(dense, 2)
        values = @view dense[:, sample]
        total = sum(values)
        total == 0 && continue
        probabilities = values ./ total
        index[sample] = sum(probabilities .^ 2)
    end
    return index
end

simpson_index(profile::CommunityProfile) = simpson_index(profile.counts)

function faith_pd(profile::CommunityProfile)
    order, presence, node_index = _tree_signal_matrix(profile.tree, profile.counts; weighted=false)
    pd = zeros(Float64, size(presence, 2))
    for node in order
        branch_length = node.branch_length
        branch_length == 0 && continue
        node_presence = presence[node_index[node], :]
        pd .+= branch_length .* Float64.(node_presence)
    end
    return pd
end

function pcoa(distances::AbstractMatrix{<:Real}; dimensions::Integer=2)
    size(distances, 1) == size(distances, 2) || throw(ArgumentError("distance matrix must be square"))
    matrix = Matrix{Float64}(distances)
    n = size(matrix, 1)
    centering = I - fill(1 / n, n, n)
    gram = -0.5 * centering * (matrix .^ 2) * centering
    decomposition = eigen(Symmetric((gram + gram') / 2))
    order = sortperm(decomposition.values; rev=true)
    values = decomposition.values[order]
    vectors = decomposition.vectors[:, order]
    positive = findall(value -> value > 0, values)
    used = min(Int(dimensions), length(positive))
    used == 0 && throw(ArgumentError("distance matrix produced no positive eigenvalues"))
    coordinates = vectors[:, 1:used] * Diagonal(sqrt.(values[1:used]))
    explained = values[1:used] ./ sum(values[positive])
    return PCoAResult(coordinates, values[1:used], explained)
end

function pcoa_plot(result::PCoAResult; labels=nothing, color=nothing, title::AbstractString="PCoA", kwargs...)
    size(result.coordinates, 2) >= 2 || throw(ArgumentError("PCoAResult must have at least two coordinates"))
    plot_obj = scatter(result.coordinates[:, 1], result.coordinates[:, 2]; legend=false, title=title, xlabel="PCoA 1", ylabel="PCoA 2", color=color, kwargs...)
    if labels !== nothing
        for (index, label) in enumerate(labels)
            annotate!(plot_obj, result.coordinates[index, 1], result.coordinates[index, 2], label)
        end
    end
    return plot_obj
end

function _nmds_stress_and_gradient(flat_coordinates::AbstractVector{<:Real}, distances::Matrix{Float64}, dimensions::Int)
    nsamples = size(distances, 1)
    coordinates = reshape(collect(Float64, flat_coordinates), nsamples, dimensions)
    gradient = zeros(Float64, nsamples, dimensions)
    numerator = 0.0
    denominator = 0.0
    for left in 1:nsamples-1
        for right in left+1:nsamples
            target = distances[left, right]
            denominator += target^2
            difference = coordinates[left, :] .- coordinates[right, :]
            distance = norm(difference)
            if distance > 0
                residual = distance - target
                numerator += residual^2
                factor = 2 * residual / distance
                gradient[left, :] .+= factor .* difference
                gradient[right, :] .-= factor .* difference
            elseif target > 0
                numerator += target^2
            end
        end
    end
    scale = denominator == 0 ? 1.0 : denominator
    return numerator / scale, vec(gradient) ./ scale
end

function nmds(distances::AbstractMatrix{<:Real}; dimensions::Integer=2, n_starts::Integer=50, maxiters::Integer=250, random_seed::Integer=1)
    matrix = Matrix{Float64}(distances)
    nsamples = size(matrix, 1)
    nsamples == size(matrix, 2) || throw(ArgumentError("distance matrix must be square"))
    dimensions > 0 || throw(ArgumentError("dimensions must be positive"))
    n_starts > 0 || throw(ArgumentError("n_starts must be positive"))

    starts = Vector{Vector{Float64}}(undef, n_starts)
    initial = pcoa(matrix; dimensions=min(dimensions, max(1, nsamples - 1))).coordinates
    padded = zeros(Float64, nsamples, dimensions)
    padded[:, 1:min(size(initial, 2), dimensions)] .= initial[:, 1:min(size(initial, 2), dimensions)]
    starts[1] = vec(padded)
    for start in 2:n_starts
        rng = MersenneTwister(random_seed + start)
        starts[start] = vec(randn(rng, nsamples, dimensions))
    end

    best_stress = Inf
    best_coordinates = zeros(Float64, nsamples, dimensions)
    best_iterations = 0
    best_converged = false
    stress_results = fill(Inf, n_starts)
    coordinate_results = Vector{Vector{Float64}}(undef, n_starts)
    iteration_results = fill(0, n_starts)
    converged_results = fill(false, n_starts)

    @threads for start in 1:n_starts
        objective = x -> first(_nmds_stress_and_gradient(x, matrix, Int(dimensions)))
        gradient! = (g, x) -> begin
            _, grad = _nmds_stress_and_gradient(x, matrix, Int(dimensions))
            g .= grad
            return nothing
        end
        result = optimize(objective, gradient!, starts[start], LBFGS())
        stress_results[start] = Optim.minimum(result)
        coordinate_results[start] = Optim.minimizer(result)
        iteration_results[start] = Optim.iterations(result)
        converged_results[start] = Optim.converged(result)
    end

    best_index = argmin(stress_results)
    best_flat = coordinate_results[best_index]
    best_coordinates = reshape(best_flat, nsamples, dimensions)
    best_stress = stress_results[best_index]
    best_iterations = iteration_results[best_index]
    best_converged = converged_results[best_index]
    return NMDSResult(best_coordinates, best_stress, best_iterations, best_converged)
end

function _welch_ttest(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    nx = length(x)
    ny = length(y)
    nx > 1 && ny > 1 || return 1.0
    mean_x = mean(x)
    mean_y = mean(y)
    var_x = var(x; corrected=true)
    var_y = var(y; corrected=true)
    standard_error = sqrt(var_x / nx + var_y / ny)
    standard_error == 0 && return isapprox(mean_x, mean_y) ? 1.0 : 0.0
    t_stat = (mean_x - mean_y) / standard_error
    degrees = (var_x / nx + var_y / ny)^2 / ((var_x^2) / (nx^2 * (nx - 1)) + (var_y^2) / (ny^2 * (ny - 1)))
    isnan(degrees) && return 1.0
    return 2 * ccdf(TDist(degrees), abs(t_stat))
end

function ancom(profile::CommunityProfile, groups::AbstractVector; pseudocount::Real=0.5, alpha::Real=0.05)
    labels = String.(groups)
    unique_labels = unique(labels)
    length(unique_labels) == 2 || throw(ArgumentError("ANCOM implementation currently expects exactly two groups"))
    dense = log.( _count_matrix_dense(profile.counts) .+ pseudocount )
    group_a = findall(==(unique_labels[1]), labels)
    group_b = findall(==(unique_labels[2]), labels)
    ntaxa = size(dense, 1)
    min_pvalues = fill(1.0, ntaxa)
    w_stat = zeros(Int, ntaxa)

    @threads for taxon in 1:ntaxa
        current_min = 1.0
        current_w = 0
        for other in 1:ntaxa
            taxon == other && continue
            ratios = dense[taxon, :] .- dense[other, :]
            pvalue = _welch_ttest(ratios[group_a], ratios[group_b])
            current_min = min(current_min, pvalue)
            current_w += pvalue < alpha ? 1 : 0
        end
        min_pvalues[taxon] = current_min
        w_stat[taxon] = current_w
    end

    qvalues = benjamini_hochberg(min_pvalues)
    significant = qvalues .< alpha
    return ANCOMResult(profile.counts.gene_ids, w_stat, min_pvalues, qvalues, significant)
end

function _logsumexp(values::AbstractVector{<:Real})
    maximum_value = maximum(values)
    return maximum_value + log(sum(exp.(values .- maximum_value)))
end

function songbird(profile::CommunityProfile, design::AbstractMatrix{<:Real}; feature_names=nothing, maxiters::Integer=300)
    dense = _count_matrix_dense(profile.counts)
    nsamples = size(dense, 2)
    ntaxa = size(dense, 1)
    size(design, 1) == nsamples || throw(ArgumentError("design matrix rows must match sample count"))
    nfeatures = size(design, 2)
    nfeatures > 0 || throw(ArgumentError("design matrix must contain at least one feature"))
    reference_taxon = ntaxa
    response = dense[1:ntaxa-1, :]
    totals = vec(sum(dense, dims=1))
    design_matrix = Matrix{Float64}(design)
    initial = zeros(Float64, nfeatures * (ntaxa - 1))

    function objective(flat_parameters)
        coefficients = reshape(flat_parameters, nfeatures, ntaxa - 1)
        linear_predictor = design_matrix * coefficients
        loglikelihood = 0.0
        for sample in 1:nsamples
            logits = vcat(linear_predictor[sample, :], 0.0)
            log_denom = _logsumexp(logits)
            loglikelihood += totals[sample] * log_denom - dot(response[:, sample], linear_predictor[sample, :])
        end
        return loglikelihood
    end

    function gradient!(storage, flat_parameters)
        coefficients = reshape(flat_parameters, nfeatures, ntaxa - 1)
        linear_predictor = design_matrix * coefficients
        errors = zeros(Float64, nsamples, ntaxa - 1)
        for sample in 1:nsamples
            logits = vcat(linear_predictor[sample, :], 0.0)
            log_denom = _logsumexp(logits)
            probabilities = exp.(logits[1:end-1] .- log_denom)
            errors[sample, :] .= totals[sample] .* probabilities .- response[:, sample]
        end
        storage .= vec(design_matrix' * errors)
        return nothing
    end

    result = optimize(objective, gradient!, initial, LBFGS())
    fitted = reshape(Optim.minimizer(result), nfeatures, ntaxa - 1)
    names = feature_names === nothing ? string.(1:nfeatures) : String.(feature_names)
    return SongbirdResult(fitted, profile.counts.gene_ids, names, -Optim.minimum(result), Optim.iterations(result), Optim.converged(result))
end

function source_tracking_posterior_summary(chain)
    samples = Matrix{Float64}(Array(chain))
    mean_proportions = vec(mean(samples, dims=1))
    median_proportions = vec(median(samples, dims=1))
    lower_bounds = [quantile(view(samples, :, column), 0.025) for column in axes(samples, 2)]
    upper_bounds = [quantile(view(samples, :, column), 0.975) for column in axes(samples, 2)]
    return SourceTrackingResult(chain, mean_proportions, median_proportions, Float64.(lower_bounds), Float64.(upper_bounds))
end

const _SOURCE_TRACKING_TURING_LOADED = Ref(false)
const _SOURCE_TRACKING_TURING_MODULE = Ref{Any}(nothing)
const _SOURCE_TRACKING_MODEL_IMPL = Ref{Any}(nothing)

function _source_tracking_turing_module()
    turing = _SOURCE_TRACKING_TURING_MODULE[]
    turing === nothing && throw(ArgumentError("Turing is required for source_tracking"))
    return turing
end

function source_tracking_model(observed::AbstractVector{<:Integer}, source_profiles::AbstractMatrix{<:Real})
    model_impl = _SOURCE_TRACKING_MODEL_IMPL[]
    model_impl === nothing && throw(ArgumentError("Turing is required for source_tracking"))
    return Base.invokelatest(model_impl, observed, source_profiles)
end

function source_tracking(observed::AbstractVector{<:Integer}, source_profiles::AbstractMatrix{<:Real}; draws::Integer=250, rng::AbstractRNG=MersenneTwister(1))
    draws > 0 || throw(ArgumentError("draws must be positive"))
    turing = _source_tracking_turing_module()
    model_impl = _SOURCE_TRACKING_MODEL_IMPL[]
    model_impl === nothing && throw(ArgumentError("Turing is required for source_tracking"))
    model = Base.invokelatest(model_impl, observed, source_profiles)
    chain = Base.invokelatest(turing.sample, rng, model, turing.NUTS(), draws)
    return source_tracking_posterior_summary(chain)
end

function _spring_layout(graph::SimpleGraph; iterations::Int=100, scale::Real=1.0)
    n = nv(graph)
    n == 0 && return zeros(Float64, 0, 2)
    coordinates = _circle_layout(n)
    temperature = 0.1
    for _ in 1:iterations
        displacement = zeros(Float64, n, 2)
        for left in 1:n-1
            for right in left+1:n
                delta = coordinates[left, :] .- coordinates[right, :]
                distance = max(norm(delta), 1e-6)
                force = scale^2 / distance
                direction = delta ./ distance
                displacement[left, :] .+= direction .* force
                displacement[right, :] .-= direction .* force
            end
        end
        for edge in edges(graph)
            left = src(edge)
            right = dst(edge)
            delta = coordinates[left, :] .- coordinates[right, :]
            distance = max(norm(delta), 1e-6)
            force = distance^2 / scale
            direction = delta ./ distance
            displacement[left, :] .-= direction .* force
            displacement[right, :] .+= direction .* force
        end
        for vertex in 1:n
            step = displacement[vertex, :]
            step_norm = norm(step)
            step_norm > 0 || continue
            coordinates[vertex, :] .+= (step ./ step_norm) .* min(step_norm, temperature)
        end
        temperature *= 0.95
    end
    return coordinates
end

function _resolve_layout(graph::SimpleGraph, layout)
    layout === nothing && return _circle_layout(nv(graph))
    layout isa AbstractMatrix{<:Real} && return Matrix{Float64}(layout)
    layout isa Function && return Matrix{Float64}(layout(graph))
    layout === :circle && return _circle_layout(nv(graph))
    layout === :random && return randn(Float64, nv(graph), 2)
    layout === :spring && return _spring_layout(graph)
    throw(ArgumentError("unsupported layout specification"))
end

function cooccurrence_network(profile::CommunityProfile; threshold::Real=0.4, layout=:circle)
    threshold >= 0 || throw(ArgumentError("threshold must be nonnegative"))
    clr = clr_transform(profile)
    correlation = cor(clr')
    ntaxa = size(correlation, 1)
    graph = SimpleGraph(ntaxa)
    weights = Dict{Tuple{Int,Int},Float64}()
    for left in 1:ntaxa-1
        for right in left+1:ntaxa
            value = correlation[left, right]
            isfinite(value) || continue
            abs(value) < threshold && continue
            add_edge!(graph, left, right)
            weights[(left, right)] = value
        end
    end
    return MicrobiomeNetwork(graph, weights, profile.counts.gene_ids, _resolve_layout(graph, layout))
end

function network_plot(network::MicrobiomeNetwork; title::AbstractString="Microbiome co-occurrence", layout=nothing, edge_scale::Real=2.5, show_labels::Bool=true, interactive::Bool=false, kwargs...)
    coordinates = layout === nothing ? network.coordinates : _resolve_layout(network.graph, layout)
    plot_obj = plot(; legend=false, aspect_ratio=:equal, title=title, xaxis=false, yaxis=false, kwargs...)
    for ((left, right), weight) in network.weights
        color = weight >= 0 ? :steelblue : :tomato
        alpha = clamp(abs(weight), 0.15, 1.0)
        linewidth = max(0.5, abs(weight) * edge_scale)
        plot!(plot_obj, [coordinates[left, 1], coordinates[right, 1]], [coordinates[left, 2], coordinates[right, 2]]; color=color, alpha=alpha, linewidth=linewidth)
    end
    scatter!(plot_obj, coordinates[:, 1], coordinates[:, 2]; markerstrokewidth=0, markersize=7, color=:black)
    show_labels || return plot_obj
    for (index, taxon) in enumerate(network.taxa)
        annotate!(plot_obj, coordinates[index, 1], coordinates[index, 2], taxon)
    end
    return plot_obj
end

network_plot(profile::CommunityProfile; kwargs...) = network_plot(cooccurrence_network(profile; kwargs...); kwargs...)

end