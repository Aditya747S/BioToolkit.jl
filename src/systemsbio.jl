module SystemsBio

using Graphs
using SimpleWeightedGraphs: SimpleWeightedGraph, SimpleWeightedEdge
using LinearAlgebra
using Statistics
using Random
using Base.Threads
using Distributions

using ..DifferentialExpression: CountMatrix, DEResult, benjamini_hochberg
using ..Enrichment: EnrichmentDatabase, EnrichmentResult
using ..BioToolkit: PhyloTree, get_terminals

export GeneNetwork, SoftThresholdResult, LimmaResult, NetworkInferenceResult, MultiOmicsFactorAnalysisResult
export pick_soft_threshold, find_modules, module_eigengenes, module_dendrogram, limma_fit, limma_deresults, gsea, infer_network, multi_omics_factor_analysis

struct GeneNetwork
    graph::SimpleWeightedGraph{Int,Float64}
    gene_to_node::Dict{String,Int}
    node_to_gene::Vector{String}
    connectivity::Vector{Float64}
    modules::Vector{Int}
end

struct SoftThresholdResult
    powers::Vector{Int}
    scale_free_fit::Vector{Float64}
    mean_connectivity::Vector{Float64}
    best_power::Int
end

struct LimmaResult
    coefficients::Matrix{Float64}
    standard_errors::Matrix{Float64}
    t_statistics::Matrix{Float64}
    pvalues::Matrix{Float64}
    qvalues::Matrix{Float64}
    moderated_variance::Vector{Float64}
    base_mean::Vector{Float64}
    gene_ids::Vector{String}
    coefficient_names::Vector{String}
end

struct NetworkInferenceResult
    graph::SimpleDiGraph
    gene_to_node::Dict{String,Int}
    node_to_gene::Vector{String}
    bic_score::Float64
end

struct MultiOmicsFactorAnalysisResult
    factors::Matrix{Float64}
    loadings::Matrix{Float64}
    assay_loadings::Vector{Matrix{Float64}}
    explained_variance::Vector{Float64}
end

_safe_log1p(x::Real) = log1p(float(x))

function _dense_expression(counts::CountMatrix)
    return _safe_log1p.(Matrix(counts.counts))
end

function _center_rows(matrix::AbstractMatrix{<:Real})
    centered = Matrix{Float64}(matrix)
    centered .-= mean(centered, dims=2)
    return centered
end

function _gene_correlation(counts::CountMatrix)
    centered = _center_rows(_dense_expression(counts))
    correlation = cor(permutedims(centered))
    correlation[.!isfinite.(correlation)] .= 0.0
    for i in axes(correlation, 1)
        correlation[i, i] = 1.0
    end
    return correlation
end

function _scale_free_fit(degrees::AbstractVector{<:Real})
    positive = collect(Float64.(degrees[degrees .> 0]))
    length(positive) < 2 && return 0.0
    rounded = unique(round.(positive; digits=4))
    frequencies = [count(x -> isapprox(x, value; atol=1e-4), positive) / length(positive) for value in rounded]
    valid = frequencies .> 0
    count(valid) < 2 && return 0.0
    x = log.(rounded[valid])
    y = log.(frequencies[valid])
    design = hcat(ones(length(x)), x)
    coefficients = design \ y
    residuals = y .- design * coefficients
    total = y .- mean(y)
    isapprox(sum(abs2, total), 0.0) && return 0.0
    r2 = max(0.0, 1.0 - sum(abs2, residuals) / sum(abs2, total))
    return sign(coefficients[2]) * r2
end

function pick_soft_threshold(counts::CountMatrix; powers::AbstractVector{<:Integer}=2:12)
    correlation = _gene_correlation(counts)
    power_values = collect(Int.(powers))
    fit = zeros(Float64, length(power_values))
    connectivity = zeros(Float64, length(power_values))
    @threads for index in eachindex(power_values)
        power = power_values[index]
        adjacency = abs.(correlation) .^ power
        adjacency[diagind(adjacency)] .= 0.0
        degrees = vec(sum(adjacency, dims=2))
        fit[index] = _scale_free_fit(degrees)
        connectivity[index] = mean(degrees)
    end
    best = power_values[argmax(fit)]
    return SoftThresholdResult(power_values, fit, connectivity, best)
end

function _tom_matrix(adjacency::AbstractMatrix{<:Real})
    weighted = Matrix{Float64}(adjacency)
    shared = weighted * weighted
    degrees = vec(sum(weighted, dims=2))
    numerator = shared .+ weighted
    denominator = min.(degrees, transpose(degrees)) .+ 1 .- weighted
    tom = numerator ./ max.(denominator, eps())
    for i in axes(tom, 1)
        tom[i, i] = 1.0
    end
    return tom
end

function _module_assignments(tom::AbstractMatrix{<:Real}; threshold::Real=0.15, min_module_size::Integer=2)
    n = size(tom, 1)
    graph = SimpleGraph(n)
    for i in 1:n-1
        for j in i+1:n
            tom[i, j] >= threshold && add_edge!(graph, i, j)
        end
    end
    components = connected_components(graph)
    modules = zeros(Int, n)
    next_module = 1
    for component in sort!(components, by=length, rev=true)
        if length(component) < min_module_size
            for node in component
                modules[node] = next_module
                next_module += 1
            end
        else
            for node in component
                modules[node] = next_module
            end
            next_module += 1
        end
    end
    return modules
end

function _module_connectivity(expression::AbstractMatrix{<:Real}, modules::Vector{Int})
    sample_matrix = permutedims(Matrix{Float64}(expression))
    connectivity = zeros(Float64, size(sample_matrix, 2))
    for module_id in unique(modules)
        indices = findall(==(module_id), modules)
        module_matrix = sample_matrix[:, indices]
        eigengene = first(svd(module_matrix; full=false).U) * first(svd(module_matrix; full=false).S)
    end
    return connectivity
end

function module_eigengenes(counts::CountMatrix, modules::AbstractVector{<:Integer})
    expression = permutedims(_dense_expression(counts))
    module_ids = sort(unique(Int.(modules)))
    eigengenes = zeros(Float64, size(expression, 1), length(module_ids))
    for (column, module_id) in enumerate(module_ids)
        indices = findall(==(module_id), modules)
        module_matrix = expression[:, indices]
        module_matrix .-= mean(module_matrix, dims=1)
        decomposition = svd(module_matrix; full=false)
        eigengenes[:, column] .= decomposition.U[:, 1] .* decomposition.S[1]
    end
    return eigengenes
end

function find_modules(counts::CountMatrix; power::Union{Nothing,Integer}=nothing, tom_threshold::Real=0.15, edge_threshold::Real=0.1, min_module_size::Integer=2)
    selected_power = power === nothing ? pick_soft_threshold(counts).best_power : Int(power)
    correlation = _gene_correlation(counts)
    adjacency = abs.(correlation) .^ selected_power
    adjacency[diagind(adjacency)] .= 0.0
    tom = _tom_matrix(adjacency)
    modules = _module_assignments(tom; threshold=tom_threshold, min_module_size=min_module_size)
    eigengenes = module_eigengenes(counts, modules)
    expression = permutedims(_dense_expression(counts))
    connectivity = zeros(Float64, size(adjacency, 1))
    for (gene_index, module_id) in enumerate(modules)
        module_index = findfirst(==(module_id), sort(unique(modules)))
        module_vector = eigengenes[:, module_index]
        gene_vector = expression[:, gene_index]
        if std(gene_vector) == 0 || std(module_vector) == 0
            connectivity[gene_index] = 0.0
        else
            connectivity[gene_index] = cor(gene_vector, module_vector)
        end
    end
    graph = SimpleWeightedGraph(length(counts.gene_ids))
    for i in 1:length(counts.gene_ids)-1
        for j in i+1:length(counts.gene_ids)
            weight = adjacency[i, j]
            if weight >= edge_threshold
                add_edge!(graph, SimpleWeightedEdge(i, j, Float64(weight)))
            end
        end
    end
    return GeneNetwork(graph, Dict(gene => index for (index, gene) in enumerate(counts.gene_ids)), String.(counts.gene_ids), connectivity, modules)
end

function module_dendrogram(network::GeneNetwork)
    children = PhyloTree[]
    for module_id in sort(unique(network.modules))
        module_nodes = findall(==(module_id), network.modules)
        leaf_nodes = [PhyloTree(network.node_to_gene[node]; branch_length=max(0.01, 1 - abs(network.connectivity[node]))) for node in module_nodes]
        push!(children, PhyloTree("module_$(module_id)"; children=leaf_nodes, branch_length=0.1))
    end
    return PhyloTree(children; name="GeneNetwork")
end

function limma_fit(counts::CountMatrix, design::AbstractMatrix{<:Real}; coefficient_names::Union{Nothing,AbstractVector{<:AbstractString}}=nothing)
    expression = permutedims(_dense_expression(counts))
    design_matrix = Matrix{Float64}(design)
    size(design_matrix, 1) == size(expression, 1) || throw(ArgumentError("design must have one row per sample"))
    qr_decomposition = qr(design_matrix)
    coefficients = qr_decomposition \ expression
    fitted = design_matrix * coefficients
    residuals = expression .- fitted
    dof = max(size(design_matrix, 1) - rank(design_matrix), 1)
    residual_variance = vec(sum(abs2, residuals; dims=1)) ./ dof
    positive = residual_variance[residual_variance .> 0]
    prior_mean = isempty(positive) ? 1.0 : mean(positive)
    prior_variance = isempty(positive) ? 1.0 : var(positive)
    prior_df = prior_variance > 0 ? max(2.0, 2.0 + 2.0 * prior_mean^2 / prior_variance) : 4.0
    moderated_variance = (prior_df .* prior_mean .+ dof .* residual_variance) ./ (prior_df + dof)
    xtx_inv = pinv(design_matrix' * design_matrix)
    coef_matrix = permutedims(Matrix{Float64}(coefficients))
    ncoef = size(coef_matrix, 2)
    standard_errors = zeros(Float64, size(coef_matrix, 1), ncoef)
    t_statistics = similar(standard_errors)
    pvalues = similar(standard_errors)
    coefficient_scales = diag(xtx_inv)
    for coef_index in 1:ncoef
        scale = coefficient_scales[coef_index]
        standard_errors[:, coef_index] .= sqrt.(max.(scale .* moderated_variance, eps()))
        t_statistics[:, coef_index] .= coef_matrix[:, coef_index] ./ standard_errors[:, coef_index]
        pvalues[:, coef_index] .= 2 .* ccdf.(TDist(dof + prior_df), abs.(t_statistics[:, coef_index]))
    end
    qvalues = similar(pvalues)
    for coef_index in 1:ncoef
        qvalues[:, coef_index] .= benjamini_hochberg(pvalues[:, coef_index])
    end
    gene_means = vec(mean(expression, dims=1))
    names = coefficient_names === nothing ? ["coef_$(index)" for index in 1:ncoef] : String.(coefficient_names)
    return LimmaResult(coef_matrix, standard_errors, t_statistics, pvalues, qvalues, moderated_variance, gene_means, String.(counts.gene_ids), names)
end

function limma_deresults(result::LimmaResult; coefficient_index::Integer=2)
    1 <= coefficient_index <= size(result.coefficients, 2) || throw(ArgumentError("coefficient_index out of bounds"))
    rows = Vector{DEResult}(undef, length(result.gene_ids))
    for (index, gene_id) in enumerate(result.gene_ids)
        rows[index] = DEResult(
            gene_id,
            result.base_mean[index],
            result.coefficients[index, coefficient_index],
            result.standard_errors[index, coefficient_index],
            result.t_statistics[index, coefficient_index],
            result.pvalues[index, coefficient_index],
            result.qvalues[index, coefficient_index],
        )
    end
    return rows
end

function _term_genes(database::EnrichmentDatabase, term_id::String, universe::Set{String})
    haskey(database.terms, term_id) || return Set{String}()
    return Set(gene for gene in database.terms[term_id].genes if gene in universe)
end

function _gsea_score(ranked_genes::Vector{String}, geneset::Set{String})
    total = length(ranked_genes)
    hits = [gene in geneset for gene in ranked_genes]
    hit_count = count(identity, hits)
    hit_count == 0 && return 0.0
    miss_count = total - hit_count
    hit_step = sqrt(hit_count / total)
    miss_step = miss_count == 0 ? 0.0 : sqrt(miss_count / total)
    running = 0.0
    best = 0.0
    for hit in hits
        running += hit ? hit_step : -miss_step
        if abs(running) > abs(best)
            best = running
        end
    end
    return best
end

function _gsea_pvalue(ranked_genes::Vector{String}, geneset::Set{String}, observed::Float64; permutations::Integer=1000, seed::Integer=1)
    permutation_scores = zeros(Float64, permutations)
    @threads for perm in 1:permutations
        rng = MersenneTwister(seed + perm)
        permutation_scores[perm] = _gsea_score(shuffle(rng, ranked_genes), geneset)
    end
    return (count(score -> abs(score) >= abs(observed), permutation_scores) + 1) / (permutations + 1)
end

function gsea(results::Vector{DEResult}, database::EnrichmentDatabase; permutations::Integer=1000, min_size::Integer=5, max_size::Integer=500, seed::Integer=1)
    sorted = sort(results, by = result -> result.stat, rev=true)
    ranked_genes = [result.gene_id for result in sorted]
    universe = Set(ranked_genes)
    background_size = length(universe)
    output = EnrichmentResult[]
    for term in values(database.terms)
        geneset = _term_genes(database, term.id, universe)
        length(geneset) < min_size && continue
        length(geneset) > max_size && continue
        observed = _gsea_score(ranked_genes, geneset)
        pvalue = _gsea_pvalue(ranked_genes, geneset, observed; permutations=permutations, seed=seed)
        overlap = count(gene -> gene in geneset, ranked_genes)
        odds_ratio = overlap == 0 ? 0.0 : (overlap / length(ranked_genes)) / (length(geneset) / max(background_size, 1))
        push!(output, EnrichmentResult(term.id, term.name, term.namespace, overlap, length(geneset), length(ranked_genes), background_size, pvalue, 1.0, odds_ratio, collect(geneset)))
    end
    isempty(output) && return output
    qvalues = benjamini_hochberg([result.pvalue for result in output])
    for (index, result) in enumerate(output)
        output[index] = EnrichmentResult(result.term_id, result.term_name, result.namespace, result.overlap, result.term_size, result.query_size, result.background_size, result.pvalue, qvalues[index], result.odds_ratio, result.genes)
    end
    return output
end

function gsea(result::LimmaResult, database::EnrichmentDatabase; coefficient_index::Integer=2, kwargs...)
    return gsea(limma_deresults(result; coefficient_index=coefficient_index), database; kwargs...)
end

function _bic_score(response::AbstractVector{<:Real}, predictors::AbstractMatrix{<:Real})
    nobs = length(response)
    if size(predictors, 2) == 0
        centered = response .- mean(response)
        rss = sum(abs2, centered)
        return -0.5 * nobs * log(rss / nobs + eps())
    end
    design = hcat(ones(nobs), predictors)
    coefficients = design \ response
    residuals = response .- design * coefficients
    rss = sum(abs2, residuals)
    nparams = size(design, 2)
    return -0.5 * nobs * log(rss / nobs + eps()) - 0.5 * nparams * log(nobs)
end

function _network_score(data::Matrix{Float64}, graph::SimpleDiGraph)
    score = 0.0
    for node in 1:size(data, 2)
        parents = inneighbors(graph, node)
        predictors = isempty(parents) ? zeros(size(data, 1), 0) : data[:, parents]
        score += _bic_score(data[:, node], predictors)
    end
    return score
end

function infer_network(data::AbstractMatrix{<:Real}; gene_ids::Union{Nothing,AbstractVector{<:AbstractString}}=nothing, max_parents::Integer=3, max_iterations::Integer=100)
    matrix = Matrix{Float64}(data)
    nsamples, nvars = size(matrix)
    labels = gene_ids === nothing ? ["node_$(index)" for index in 1:nvars] : String.(gene_ids)
    graph = SimpleDiGraph(nvars)
    best_score = _network_score(matrix, graph)
    for _ in 1:max_iterations
        improvement = false
        best_candidate = best_score
        best_edge = nothing
        best_graph = graph
        for parent in 1:nvars
            for child in 1:nvars
                parent == child && continue
                has_edge(graph, parent, child) && continue
                outdegree(graph, parent) >= max_parents && continue
                candidate = copy(graph)
                add_edge!(candidate, parent, child)
                is_cyclic(candidate) && continue
                candidate_score = _network_score(matrix, candidate)
                if candidate_score > best_candidate
                    best_candidate = candidate_score
                    best_edge = (parent, child)
                    best_graph = candidate
                end
            end
        end
        if best_edge === nothing
            break
        end
        graph = best_graph
        best_score = best_candidate
        improvement = true
        improvement || break
    end
    return NetworkInferenceResult(graph, Dict(label => index for (index, label) in enumerate(labels)), labels, best_score)
end

function infer_network(counts::CountMatrix; kwargs...)
    return infer_network(permutedims(_dense_expression(counts)); gene_ids=counts.gene_ids, kwargs...)
end

function _zscore(matrix::AbstractMatrix{<:Real})
    data = Matrix{Float64}(matrix)
    for column in axes(data, 2)
        μ = mean(data[:, column])
        σ = std(data[:, column])
        data[:, column] .= σ == 0 ? 0.0 : (data[:, column] .- μ) ./ σ
    end
    return data
end

function multi_omics_factor_analysis(assays::Vector{<:AbstractMatrix{<:Real}}; n_factors::Integer=3)
    length(assays) > 0 || throw(ArgumentError("assays must not be empty"))
    nsamples = size(assays[1], 1)
    all(size(assay, 1) == nsamples for assay in assays) || throw(ArgumentError("all assays must have the same number of samples"))
    standardized = [_zscore(Matrix{Float64}(assay)) for assay in assays]
    concatenated = hcat(standardized...)
    decomposition = svd(concatenated; full=false)
    used = min(Int(n_factors), length(decomposition.S))
    factors = decomposition.U[:, 1:used] * Diagonal(decomposition.S[1:used])
    loadings = decomposition.V[:, 1:used]
    explained_variance = decomposition.S[1:used].^2 ./ sum(decomposition.S .^ 2)
    assay_loadings = Matrix{Float64}[]
    offset = 1
    for assay in standardized
        next_offset = offset + size(assay, 2) - 1
        push!(assay_loadings, loadings[offset:next_offset, 1:used])
        offset = next_offset + 1
    end
    return MultiOmicsFactorAnalysisResult(factors, loadings, assay_loadings, explained_variance)
end

end