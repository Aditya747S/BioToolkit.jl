# ==============================================================================
# microbiome.jl — Microbiome community analysis
#
# Provides alpha/beta diversity metrics, ordination (PCoA, NMDS),
# compositional transforms (CLR, ILR), differential abundance (ANCOM),
# co-occurrence networks, and source tracking.
#
# References:
#   - Bray & Curtis (1957) Pacific Naturalist 28:325-342 (BC distance)
#   - Lozupone & Knight (2005) AEM 71(12):8228-8235 (UniFrac)
#   - Mandal et al. (2015) Microbiome 3:38 (ANCOM)
# ==============================================================================

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
using ..BioToolkit: BioSequence, DNAAlphabet, DNASeq, PhyloTree, get_terminals
using ..BioToolkit: AbstractAnalysisResult, ProvenanceContext, ProvenanceParams, ResultProvenance, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_record, provenance_result!, register_provenance!

@inline function _register_microbiome_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export CommunityProfile, PCoAResult, NMDSResult, ANCOMResult, SongbirdResult, SourceTrackingResult, MicrobiomeNetwork
export clr_transform, ilr_transform, bray_curtis, unifrac, weighted_unifrac, shannon_entropy, simpson_index, faith_pd
export pairwise_bray_curtis, pairwise_unifrac, pcoa, pcoa_plot, nmds, ancom, songbird, cooccurrence_network, network_plot, source_tracking_model, source_tracking, source_tracking_posterior_summary
export mag_bin_contigs, kraken_like_classify, viral_contig_scores, humann_like_pathways, strainge_like_variants
export lca_taxonomy_from_votes, strain_haplotype_profile

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

CommunityProfile(counts::AbstractMatrix{<:Integer}, gene_ids::AbstractVector{<:String}, sample_ids::AbstractVector{<:String}, taxonomy::DataFrame, tree::PhyloTree, metadata::DataFrame) =
    CommunityProfile(CountMatrix(counts, gene_ids, sample_ids), taxonomy, tree, metadata)

struct PCoAResult <: AbstractAnalysisResult
    coordinates::Matrix{Float64}
    eigenvalues::Vector{Float64}
    variance_explained::Vector{Float64}
    provenance::ResultProvenance
end

struct NMDSResult <: AbstractAnalysisResult
    coordinates::Matrix{Float64}
    stress::Float64
    iterations::Int
    converged::Bool
    provenance::ResultProvenance
end

struct ANCOMResult <: AbstractAnalysisResult
    taxon_ids::Vector{String}
    w_stat::Vector{Int}
    min_pvalue::Vector{Float64}
    qvalue::Vector{Float64}
    significant::Vector{Bool}
    provenance::ResultProvenance
end

struct SongbirdResult <: AbstractAnalysisResult
    coefficients::Matrix{Float64}
    taxon_ids::Vector{String}
    feature_names::Vector{String}
    loglik::Float64
    iterations::Int
    converged::Bool
    provenance::ResultProvenance
end

struct SourceTrackingResult <: AbstractAnalysisResult
    chain
    mean_proportions::Vector{Float64}
    median_proportions::Vector{Float64}
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    provenance::ResultProvenance
end

PCoAResult(coordinates::Matrix{Float64}, eigenvalues::Vector{Float64}, variance_explained::Vector{Float64}) =
    PCoAResult(coordinates, eigenvalues, variance_explained, provenance_record("PCoAResult", "Microbiome/pcoa"))

NMDSResult(coordinates::Matrix{Float64}, stress::Float64, iterations::Int, converged::Bool) =
    NMDSResult(coordinates, stress, iterations, converged, provenance_record("NMDSResult", "Microbiome/nmds"; status=converged ? :ok : :warn))

ANCOMResult(taxon_ids::Vector{String}, w_stat::Vector{Int}, min_pvalue::Vector{Float64}, qvalue::Vector{Float64}, significant::Vector{Bool}) =
    ANCOMResult(taxon_ids, w_stat, min_pvalue, qvalue, significant, provenance_record("ANCOMResult", "Microbiome/ancom"))

SongbirdResult(coefficients::Matrix{Float64}, taxon_ids::Vector{String}, feature_names::Vector{String}, loglik::Float64, iterations::Int, converged::Bool) =
    SongbirdResult(coefficients, taxon_ids, feature_names, loglik, iterations, converged, provenance_record("SongbirdResult", "Microbiome/songbird"; status=converged ? :ok : :warn))

SourceTrackingResult(chain, mean_proportions::Vector{Float64}, median_proportions::Vector{Float64}, lower_bounds::Vector{Float64}, upper_bounds::Vector{Float64}) =
    SourceTrackingResult(chain, mean_proportions, median_proportions, lower_bounds, upper_bounds, provenance_record("SourceTrackingResult", "Microbiome/source_tracking"))

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
    selector isa AbstractVector{<:String} && return [findfirst(==(String(value)), labels) for value in selector]
    selector isa AbstractVector{Symbol} && return [findfirst(==(String(value)), labels) for value in selector]
    selector isa String && return [findfirst(==(String(selector)), labels)]
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
    nsamples = size(counts.counts, 2)
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

    all(isfinite, signal) || throw(ArgumentError("tree signal matrix contains non-finite values"))
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

function clr_transform(counts::CountMatrix; pseudocount::Real=0.5, multi_thread::Bool=true)
    pseudocount > 0 || throw(ArgumentError("pseudocount must be positive"))
    dense = _count_matrix_dense(counts)
    transformed = similar(dense)
    if multi_thread && size(dense, 2) > 1 && Threads.nthreads() > 1
        @threads for sample in axes(dense, 2)
            values = @view dense[:, sample]
            logged = log.(values .+ pseudocount)
            centered = logged .- mean(logged)
            transformed[:, sample] = centered
        end
    else
        for sample in axes(dense, 2)
            values = @view dense[:, sample]
            logged = log.(values .+ pseudocount)
            centered = logged .- mean(logged)
            transformed[:, sample] = centered
        end
    end
    _ctx = active_provenance_context()
    return _register_microbiome_result!(_ctx, transformed, "clr_transform"; parents=provenance_parent_ids(counts), parameters=(pseudocount=Float64(pseudocount), multi_thread=multi_thread, n_taxa=size(counts.counts,1)))
end

clr_transform(profile::CommunityProfile; pseudocount::Real=0.5, multi_thread::Bool=true) = clr_transform(profile.counts; pseudocount=pseudocount, multi_thread=multi_thread)

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

function ilr_transform(counts::CountMatrix; pseudocount::Real=0.5, multi_thread::Bool=true)
    clr = clr_transform(counts; pseudocount=pseudocount, multi_thread=multi_thread)
    result = _helmert_basis(size(clr, 1))' * clr
    _ctx = active_provenance_context()
    return _register_microbiome_result!(_ctx, result, "ilr_transform"; parents=provenance_parent_ids(counts), parameters=(pseudocount=Float64(pseudocount)))
end

ilr_transform(profile::CommunityProfile; pseudocount::Real=0.5, multi_thread::Bool=true) = ilr_transform(profile.counts; pseudocount=pseudocount, multi_thread=multi_thread)

function bray_curtis(counts::CountMatrix; multi_thread::Bool=true)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, pairwise_bray_curtis(counts; multi_thread=multi_thread), "bray_curtis")
end

function pairwise_bray_curtis(counts::CountMatrix; multi_thread::Bool=true)
    rows, columns, values = findnz(counts.counts)
    nsamples = size(counts.counts, 2)
    sample_rows = Vector{Vector{Int}}(undef, nsamples)
    sample_values = Vector{Vector{Float64}}(undef, nsamples)

    position = 1
    for sample in 1:nsamples
        sample_start = position
        while position <= length(columns) && columns[position] == sample
            position += 1
        end
        if position == sample_start
            sample_rows[sample] = Int[]
            sample_values[sample] = Float64[]
        else
            sample_rows[sample] = rows[sample_start:position-1]
            sample_values[sample] = Float64.(values[sample_start:position-1])
        end
    end

    distances = zeros(Float64, nsamples, nsamples)
    if multi_thread && nsamples > 1 && Threads.nthreads() > 1
        @threads for left in 1:nsamples
            for right in left+1:nsamples
                left_rows = sample_rows[left]
                right_rows = sample_rows[right]
                left_values = sample_values[left]
                right_values = sample_values[right]
                denominator = sum(left_values) + sum(right_values)
                if denominator == 0
                    distances[left, right] = NaN
                    distances[right, left] = NaN
                    continue
                end

                left_position = 1
                right_position = 1
                numerator = 0.0
                while left_position <= length(left_rows) && right_position <= length(right_rows)
                    left_row = left_rows[left_position]
                    right_row = right_rows[right_position]
                    if left_row == right_row
                        numerator += abs(left_values[left_position] - right_values[right_position])
                        left_position += 1
                        right_position += 1
                    elseif left_row < right_row
                        numerator += left_values[left_position]
                        left_position += 1
                    else
                        numerator += right_values[right_position]
                        right_position += 1
                    end
                end
                while left_position <= length(left_rows)
                    numerator += left_values[left_position]
                    left_position += 1
                end
                while right_position <= length(right_rows)
                    numerator += right_values[right_position]
                    right_position += 1
                end
                value = numerator / denominator
                distances[left, right] = value
                distances[right, left] = value
            end
        end
    else
        for left in 1:nsamples
            for right in left+1:nsamples
                left_rows = sample_rows[left]
                right_rows = sample_rows[right]
                left_values = sample_values[left]
                right_values = sample_values[right]
                denominator = sum(left_values) + sum(right_values)
                if denominator == 0
                    distances[left, right] = NaN
                    distances[right, left] = NaN
                    continue
                end

                left_position = 1
                right_position = 1
                numerator = 0.0
                while left_position <= length(left_rows) && right_position <= length(right_rows)
                    left_row = left_rows[left_position]
                    right_row = right_rows[right_position]
                    if left_row == right_row
                        numerator += abs(left_values[left_position] - right_values[right_position])
                        left_position += 1
                        right_position += 1
                    elseif left_row < right_row
                        numerator += left_values[left_position]
                        left_position += 1
                    else
                        numerator += right_values[right_position]
                        right_position += 1
                    end
                end
                while left_position <= length(left_rows)
                    numerator += left_values[left_position]
                    left_position += 1
                end
                while right_position <= length(right_rows)
                    numerator += right_values[right_position]
                    right_position += 1
                end
                value = numerator / denominator
                distances[left, right] = value
                distances[right, left] = value
            end
        end
    end
    _ctx = active_provenance_context()
    return _register_microbiome_result!(_ctx, distances, "pairwise_bray_curtis"; parents=provenance_parent_ids(counts), parameters=(n_samples=nsamples, multi_thread=multi_thread))
end

pairwise_bray_curtis(profile::CommunityProfile; multi_thread::Bool=true) = pairwise_bray_curtis(profile.counts; multi_thread=multi_thread)

bray_curtis(profile::CommunityProfile; multi_thread::Bool=true) = bray_curtis(profile.counts; multi_thread=multi_thread)

function unifrac(profile::CommunityProfile; multi_thread::Bool=true)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, pairwise_unifrac(profile; multi_thread=multi_thread), "unifrac")
end

function pairwise_unifrac(profile::CommunityProfile; weighted::Bool=false, multi_thread::Bool=true)
    order, signal, node_index = _tree_signal_matrix(profile.tree, profile.counts; weighted=weighted)
    if weighted
        abundances = signal
        nsamples = size(abundances, 2)
        distances = zeros(Float64, nsamples, nsamples)
        total_length = sum(node.branch_length for node in order if node.branch_length > 0)
        total_length == 0 && return distances
        if multi_thread && nsamples > 1 && Threads.nthreads() > 1
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
        else
            for left in 1:nsamples
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
        if multi_thread && nsamples > 1 && Threads.nthreads() > 1
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
        else
            for left in 1:nsamples
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
    end
    total_length = sum(node.branch_length for node in order if node.branch_length > 0)
    total_length == 0 && return distances
    result = distances ./ total_length
    _ctx = active_provenance_context()
    return _register_microbiome_result!(_ctx, result, "pairwise_unifrac"; parents=provenance_parent_ids(profile), parameters=(weighted=weighted, multi_thread=multi_thread, n_samples=size(distances,1)))
end

weighted_unifrac(profile::CommunityProfile; multi_thread::Bool=true) = pairwise_unifrac(profile; weighted=true, multi_thread=multi_thread)

pairwise_unifrac(counts::CountMatrix, tree; weighted::Bool=false, multi_thread::Bool=true) = pairwise_unifrac(CommunityProfile(counts, String.(counts.gene_ids), String.(counts.sample_ids), nothing, tree, nothing); weighted=weighted, multi_thread=multi_thread)

function shannon_entropy(counts::CountMatrix)
    _, columns, values = findnz(counts.counts)
    totals = zeros(Float64, size(counts.counts, 2))
    for position in eachindex(values)
        totals[columns[position]] += Float64(values[position])
    end

    entropy = zeros(Float64, length(totals))
    for sample in eachindex(totals)
        total = totals[sample]
        total == 0 && continue
        entropy_value = 0.0
        for position in eachindex(values)
            columns[position] == sample || continue
            probability = Float64(values[position]) / total
            probability == 0 && continue
            entropy_value -= probability * log(probability)
        end
        entropy[sample] = entropy_value
    end
    _ctx = active_provenance_context()
    return _register_microbiome_result!(_ctx, entropy, "shannon_entropy"; parents=provenance_parent_ids(counts), parameters=(n_taxa=size(counts.counts,1), n_samples=size(counts.counts,2)))
end

shannon_entropy(profile::CommunityProfile) = shannon_entropy(profile.counts)

function simpson_index(counts::CountMatrix)
    _, columns, values = findnz(counts.counts)
    totals = zeros(Float64, size(counts.counts, 2))
    for position in eachindex(values)
        totals[columns[position]] += Float64(values[position])
    end

    index = zeros(Float64, length(totals))
    for sample in eachindex(totals)
        total = totals[sample]
        total == 0 && continue
        index_value = 0.0
        for position in eachindex(values)
            columns[position] == sample || continue
            probability = Float64(values[position]) / total
            index_value += probability^2
        end
        index[sample] = index_value
    end
    _ctx = active_provenance_context()
    return _register_microbiome_result!(_ctx, index, "simpson_index"; parents=provenance_parent_ids(counts), parameters=(n_taxa=size(counts.counts,1), n_samples=size(counts.counts,2)))
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
    _ctx = active_provenance_context()
    return _register_microbiome_result!(_ctx, pd, "faith_pd"; parents=provenance_parent_ids(profile), parameters=(n_taxa=size(profile.counts.counts,1)))
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
    result = PCoAResult(
        coordinates,
        values[1:used],
        explained,
        provenance_record(
            "PCoAResult",
            "Microbiome/pcoa";
            notes=["classical multidimensional scaling from distance matrix"],
            parameters=(dimensions=used, sample_count=n)))
    _ctx = active_provenance_context()
    return _register_microbiome_result!(_ctx, result, "pcoa"; parents=provenance_parent_ids(distances), parameters=(dimensions=used, sample_count=n))
end

function pcoa_plot(result::PCoAResult; labels=nothing, color=nothing, title::String="PCoA", kwargs...)
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
                factor = 2 * residual / max(distance, eps(Float64))
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

function nmds(distances::AbstractMatrix{<:Real}; dimensions::Integer=2, n_starts::Integer=50, maxiters::Integer=250, random_seed::Integer=1, multi_thread::Bool=true)
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

    if multi_thread && n_starts > 1 && Threads.nthreads() > 1
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
    else
        for start in 1:n_starts
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
    end

    best_index = argmin(stress_results)
    best_flat = coordinate_results[best_index]
    best_coordinates = reshape(best_flat, nsamples, dimensions)
    best_stress = stress_results[best_index]
    best_iterations = iteration_results[best_index]
    best_converged = converged_results[best_index]
    warnings = best_converged ? String[] : ["optimizer did not report convergence"]
    fallbacks = best_converged ? String[] : ["best stress solution across random starts was retained"]
    result = NMDSResult(
        best_coordinates,
        best_stress,
        best_iterations,
        best_converged,
        provenance_record(
            "NMDSResult",
            "Microbiome/nmds";
            status=best_converged ? :ok : :warn,
            warnings=warnings,
            fallbacks=fallbacks,
            notes=["best initialization selected from multiple starts"],
            parameters=(dimensions=Int(dimensions), n_starts=Int(n_starts), maxiters=Int(maxiters), random_seed=Int(random_seed), multi_thread=Bool(multi_thread), sample_count=nsamples)))
    _ctx = active_provenance_context()
    return _register_microbiome_result!(_ctx, result, "nmds"; parents=provenance_parent_ids(distances), parameters=(dimensions=Int(dimensions), n_starts=Int(n_starts), converged=best_converged, sample_count=nsamples))
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

function ancom(profile::CommunityProfile, groups::AbstractVector; pseudocount::Real=0.5, alpha::Real=0.05, multi_thread::Bool=true)
    labels = String.(groups)
    unique_labels = unique(labels)
    length(unique_labels) == 2 || throw(ArgumentError("ANCOM implementation currently expects exactly two groups"))
    dense = log.( _count_matrix_dense(profile.counts) .+ pseudocount )
    group_a = findall(==(unique_labels[1]), labels)
    group_b = findall(==(unique_labels[2]), labels)
    ntaxa = size(dense, 1)
    min_pvalues = fill(1.0, ntaxa)
    w_stat = zeros(Int, ntaxa)

    if multi_thread && ntaxa > 1 && Threads.nthreads() > 1
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
    else
        for taxon in 1:ntaxa
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
    end

    qvalues = benjamini_hochberg(min_pvalues)
    significant = qvalues .< alpha
    result = ANCOMResult(
        profile.counts.gene_ids,
        w_stat,
        min_pvalues,
        qvalues,
        significant,
        provenance_record(
            "ANCOMResult",
            "Microbiome/ancom";
            notes=["pairwise log-ratio Welch tests across taxa"],
            parameters=(taxon_count=ntaxa, group_count=length(unique_labels), pseudocount=Float64(pseudocount), alpha=Float64(alpha), multi_thread=Bool(multi_thread))))
    _ctx = active_provenance_context()
    return _register_microbiome_result!(_ctx, result, "ancom"; parents=provenance_parent_ids(profile), parameters=(taxon_count=ntaxa, n_significant=count(identity, significant), alpha=Float64(alpha)))
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
    converged = Optim.converged(result)
    return SongbirdResult(
        fitted,
        profile.counts.gene_ids,
        names,
        -Optim.minimum(result),
        Optim.iterations(result),
        converged,
        provenance_record(
            "SongbirdResult",
            "Microbiome/songbird";
            status=converged ? :ok : :warn,
            warnings=converged ? String[] : ["optimizer stopped before convergence"],
            fallbacks=converged ? String[] : ["best available coefficient estimate was retained"],
            parameters=(sample_count=nsamples, taxon_count=ntaxa, feature_count=nfeatures, maxiters=Int(maxiters))))
end

function source_tracking_posterior_summary(chain; source::AbstractString="Microbiome/source_tracking_posterior_summary", notes::AbstractVector{<:AbstractString}=String[], parameters::NamedTuple=NamedTuple())
    samples = Matrix{Float64}(Array(chain))
    mean_proportions = vec(mean(samples, dims=1))
    median_proportions = vec(median(samples, dims=1))
    lower_bounds = [quantile(view(samples, :, column), 0.025) for column in axes(samples, 2)]
    upper_bounds = [quantile(view(samples, :, column), 0.975) for column in axes(samples, 2)]
    return SourceTrackingResult(
        chain,
        mean_proportions,
        median_proportions,
        Float64.(lower_bounds),
        Float64.(upper_bounds),
        provenance_record("SourceTrackingResult", source; notes=String.(notes), parameters=parameters))
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
    return source_tracking_posterior_summary(
        chain;
        source="Microbiome/source_tracking",
        notes=["Turing NUTS posterior summary"],
        parameters=(draws=Int(draws), source_count=size(source_profiles, 2), feature_count=length(observed)))
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

function network_plot(network::MicrobiomeNetwork; title::String="Microbiome co-occurrence", layout=nothing, edge_scale::Real=2.5, show_labels::Bool=true, interactive::Bool=false, kwargs...)
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

function _meta_zscore_columns(x::AbstractMatrix{<:Real})
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

function _meta_kmeans_lloyd(x::AbstractMatrix{<:Real}, k::Int; max_iter::Int=50, seed::Int=1)
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
    mag_bin_contigs(kmer_features, coverage; n_bins=10)

MetaBAT-style unsupervised MAG binning on composition + coverage features.
"""
function mag_bin_contigs(kmer_features::AbstractMatrix{<:Real}, coverage::AbstractVector{<:Real}; n_bins::Int=10, seed::Int=1)
    size(kmer_features, 1) == length(coverage) || throw(DimensionMismatch("coverage length must match number of contigs"))
    joint = hcat(_meta_zscore_columns(kmer_features), _meta_zscore_columns(reshape(Float64.(coverage), :, 1)))
    k = clamp(n_bins, 1, size(joint, 1))
    labels, _ = _meta_kmeans_lloyd(joint, k; seed=seed)
    return DataFrame(contig_id=["contig_$(i)" for i in 1:length(coverage)], bin=labels, coverage=Float64.(coverage))
end

"""
    kraken_like_classify(reads, kmer_db; k=31)

Kraken-like exact k-mer voting classifier.
"""
function kraken_like_classify(reads::AbstractVector{<:AbstractString}, kmer_db::AbstractDict{<:AbstractString,<:AbstractString}; k::Int=31)
    out = DataFrame(read_id=String[], assigned_taxon=String[], support=Int[])
    for (i, read) in enumerate(reads)
        r = uppercase(String(read))
        votes = Dict{String,Int}()
        if ncodeunits(r) >= k
            for s in 1:(ncodeunits(r) - k + 1)
                km = r[s:(s + k - 1)]
                if haskey(kmer_db, km)
                    tax = String(kmer_db[km])
                    votes[tax] = get(votes, tax, 0) + 1
                end
            end
        end
        if isempty(votes)
            push!(out, ("read_$(i)", "unclassified", 0))
        else
            sup, tax = findmax(votes)
            push!(out, ("read_$(i)", tax, sup))
        end
    end
    return out
end

kraken_like_classify(reads::AbstractVector{<:BioSequence{DNAAlphabet}}, kmer_db::AbstractDict{<:AbstractString,<:AbstractString}; k::Int=31) = kraken_like_classify(String.(reads), kmer_db; k=k)

"""
    viral_contig_scores(contigs)

Heuristic viral-likeness scoring from motif burden and GC content.
"""
function viral_contig_scores(contigs::AbstractVector{<:AbstractString})
    motif = Set(["TATA", "AATAAA", "TTTT"])
    out = DataFrame(contig_id=String[], length=Int[], gc=Float64[], motif_hits=Int[], viral_score=Float64[])
    for (i, contig) in enumerate(contigs)
        seq = uppercase(String(contig))
        len = ncodeunits(seq)
        gc = count(c -> c in ('G', 'C'), collect(seq)) / max(len, 1)
        hits = 0
        for m in motif
            hits += length(collect(eachmatch(Regex(m), seq)))
        end
        score = 0.5 * gc + 0.5 * tanh(hits / 5)
        push!(out, ("contig_$(i)", len, gc, hits, score))
    end
    sort!(out, :viral_score, rev=true)
    return out
end

viral_contig_scores(contigs::AbstractVector{<:BioSequence{DNAAlphabet}}) = viral_contig_scores(String.(contigs))

"""
    humann_like_pathways(uniref_hits, mapping)

Aggregate UniRef abundances to pathway-level totals.
"""
function humann_like_pathways(uniref_hits::DataFrame, mapping::AbstractDict{<:AbstractString,<:AbstractString}; feature_col::Symbol=:feature, abundance_col::Symbol=:abundance)
    hasproperty(uniref_hits, feature_col) || throw(ArgumentError("feature column missing"))
    hasproperty(uniref_hits, abundance_col) || throw(ArgumentError("abundance column missing"))

    pathway = String[]
    abundance = Float64[]
    for row in eachrow(uniref_hits)
        feat = String(row[feature_col])
        haskey(mapping, feat) || continue
        push!(pathway, String(mapping[feat]))
        push!(abundance, Float64(row[abundance_col]))
    end

    df = DataFrame(pathway=pathway, abundance=abundance)
    out = combine(groupby(df, :pathway), :abundance => sum => :pathway_abundance)
    sort!(out, :pathway_abundance, rev=true)
    return out
end

"""
    strainge_like_variants(allele_depth)

StrainGE-style intra-species variant flagging from depth/VAF thresholds.
"""
function strainge_like_variants(allele_depth::DataFrame; ref_col::Symbol=:ref_count, alt_col::Symbol=:alt_count, min_depth::Int=10, min_vaf::Real=0.02)
    hasproperty(allele_depth, ref_col) || throw(ArgumentError("missing ref count column"))
    hasproperty(allele_depth, alt_col) || throw(ArgumentError("missing alt count column"))

    out = copy(allele_depth)
    depth = Int.(out[!, ref_col] .+ out[!, alt_col])
    vaf = Float64.(out[!, alt_col]) ./ max.(depth, 1)
    out[!, :depth] = depth
    out[!, :vaf] = vaf
    out[!, :is_variant] = (depth .>= min_depth) .& (vaf .>= Float64(min_vaf))
    return out
end

function _lca_taxonomy(candidates::AbstractVector{<:AbstractString}; rank_sep::Char=';')
    toks = [split(String(c), rank_sep) for c in candidates if !isempty(strip(String(c))) && lowercase(String(c)) != "unclassified"]
    isempty(toks) && return "unclassified"
    min_depth = minimum(length, toks)
    prefix = String[]
    for i in 1:min_depth
        label = toks[1][i]
        all(t -> t[i] == label, toks) || break
        push!(prefix, label)
    end
    return isempty(prefix) ? "root" : join(prefix, string(rank_sep))
end

"""
    lca_taxonomy_from_votes(votes)

Resolve read-level taxonomy by lowest common ancestor from classifier vote lists.
Votes can be provided as vectors of taxon strings or pipe-delimited strings.
"""
function lca_taxonomy_from_votes(votes::AbstractVector; candidate_delim::Char='|', rank_sep::Char=';')
    out = DataFrame(read_id=String[], lca_taxon=String[], dominant_taxon=String[], support=Int[], n_candidates=Int[])
    for (i, entry) in enumerate(votes)
        candidates = if entry isa AbstractString
            [strip(x) for x in split(String(entry), candidate_delim) if !isempty(strip(x))]
        elseif entry isa AbstractVector
            String.(entry)
        else
            throw(ArgumentError("each vote entry must be a string or vector of strings"))
        end

        if isempty(candidates)
            push!(out, ("read_$(i)", "unclassified", "unclassified", 0, 0))
            continue
        end

        counts = Dict{String,Int}()
        for c in candidates
            counts[c] = get(counts, c, 0) + 1
        end
        support, dominant = findmax(counts)
        lca = _lca_taxonomy(candidates; rank_sep=rank_sep)
        push!(out, ("read_$(i)", lca, dominant, support, length(candidates)))
    end
    return out
end

"""
    strain_haplotype_profile(variant_table)

Construct per-sample strain haplotypes from position/allele calls.
"""
function strain_haplotype_profile(variant_table::DataFrame; sample_col::Symbol=:sample, position_col::Symbol=:position, allele_col::Symbol=:allele)
    hasproperty(variant_table, sample_col) || throw(ArgumentError("missing sample column"))
    hasproperty(variant_table, position_col) || throw(ArgumentError("missing position column"))
    hasproperty(variant_table, allele_col) || throw(ArgumentError("missing allele column"))

    out = DataFrame(sample=String[], n_variants=Int[], haplotype=String[])
    for g in groupby(variant_table, sample_col)
        ord = sortperm(Int.(g[!, position_col]))
        parts = ["$(g[ord[i], position_col]):$(g[ord[i], allele_col])" for i in eachindex(ord)]
        push!(out, (String(g[1, sample_col]), nrow(g), join(parts, "|")))
    end
    sort!(out, :sample)
    return out
end

end
