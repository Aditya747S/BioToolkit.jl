module SystemsBio

using Graphs
using DataFrames
using SimpleWeightedGraphs: SimpleWeightedGraph, SimpleWeightedEdge
using LinearAlgebra
using Statistics
using Random
using Base.Threads
using Distributions
using SpecialFunctions: digamma, trigamma, polygamma, erfc

using ..DifferentialExpression: CountMatrix, DEResult, benjamini_hochberg, calc_norm_factors
using ..Enrichment: EnrichmentDatabase, EnrichmentResult
using ..BioToolkit: PhyloTree, get_terminals, maybe_to_device, maybe_to_host, resolve_backend, threaded_foreach
using ..BioToolkit: AbstractAnalysisResult, ProvenanceContext, ProvenanceParams, ResultProvenance, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_record, provenance_result!, register_provenance!

@inline function _register_systemsbio_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export GeneNetwork, SoftThresholdResult, LimmaResult, VoomResult, NetworkInferenceResult, MultiOmicsFactorAnalysisResult
export pick_soft_threshold, find_modules, module_eigengenes, module_dendrogram, voom_transform, voom_with_quality_weights, estimate_limma_hyperparameters, limma_moderate_ttest, eBayes, contrasts_fit, duplicateCorrelation, remove_batch_effect_limma, limma_fit, limma_deresults, gsea, infer_network, multi_omics_factor_analysis
export mofa_plus_integration, cca_integration, anchor_based_integration, totalvi_like_integration, causal_grn_inference, gene_program_nmf
export mofa_plus_em, sparse_cca_integration
export causal_ate_regression

struct GeneNetwork
    graph::SimpleWeightedGraph{Int,Float64}
    gene_to_node::Dict{String,Int}
    node_to_gene::Vector{String}
    connectivity::Vector{Float64}
    modules::Vector{Int}
end

struct SoftThresholdResult <: AbstractAnalysisResult
    powers::Vector{Int}
    scale_free_fit::Vector{Float64}
    mean_connectivity::Vector{Float64}
    best_power::Int
    provenance::ResultProvenance
end

SoftThresholdResult(powers, scale_free_fit, mean_connectivity, best_power) =
    SoftThresholdResult(powers, scale_free_fit, mean_connectivity, best_power, provenance_record("SoftThresholdResult", "systemsbio"))

struct LimmaResult <: AbstractAnalysisResult
    coefficients::Matrix{Float64}
    standard_errors::Matrix{Float64}
    t_statistics::Matrix{Float64}
    pvalues::Matrix{Float64}
    qvalues::Matrix{Float64}
    moderated_variance::Vector{Float64}
    base_mean::Vector{Float64}
    gene_ids::Vector{String}
    coefficient_names::Vector{String}
    provenance::ResultProvenance
end

LimmaResult(coefficients, standard_errors, t_statistics, pvalues, qvalues, moderated_variance, base_mean, gene_ids, coefficient_names) =
    LimmaResult(coefficients, standard_errors, t_statistics, pvalues, qvalues, moderated_variance, base_mean, gene_ids, coefficient_names, provenance_record("LimmaResult", "systemsbio"))

struct VoomResult <: AbstractAnalysisResult
    log_cpm::Matrix{Float64}
    weights::Matrix{Float64}
    fitted_logcount::Matrix{Float64}
    effective_lib_sizes::Vector{Float64}
    design::Matrix{Float64}
    provenance::ResultProvenance
end

VoomResult(log_cpm, weights, fitted_logcount, effective_lib_sizes, design) =
    VoomResult(log_cpm, weights, fitted_logcount, effective_lib_sizes, design, provenance_record("VoomResult", "systemsbio"))

struct NetworkInferenceResult <: AbstractAnalysisResult
    graph::SimpleDiGraph
    gene_to_node::Dict{String,Int}
    node_to_gene::Vector{String}
    bic_score::Float64
    provenance::ResultProvenance
end

NetworkInferenceResult(graph, gene_to_node, node_to_gene, bic_score) =
    NetworkInferenceResult(graph, gene_to_node, node_to_gene, bic_score, provenance_record("NetworkInferenceResult", "systemsbio"))

struct MultiOmicsFactorAnalysisResult <: AbstractAnalysisResult
    factors::Matrix{Float64}
    loadings::Matrix{Float64}
    assay_loadings::Vector{Matrix{Float64}}
    explained_variance::Vector{Float64}
    provenance::ResultProvenance
end

MultiOmicsFactorAnalysisResult(factors, loadings, assay_loadings, explained_variance) =
    MultiOmicsFactorAnalysisResult(factors, loadings, assay_loadings, explained_variance, provenance_record("MultiOmicsFactorAnalysisResult", "systemsbio"))

MultiOmicsFactorAnalysisResult(factors::Matrix{Float64}, loadings::Matrix{Float64}, assay_loadings::Vector{Matrix{Float64}}, explained_variance::Vector{Float64}) =
    MultiOmicsFactorAnalysisResult(factors, loadings, assay_loadings, explained_variance, provenance_record("MultiOmicsFactorAnalysisResult", "SystemsBio/multi_omics_factor_analysis"))

_safe_log1p(x::Real) = log1p(float(x))

function _dense_expression(counts::CountMatrix)
    return _safe_log1p.(Matrix(counts.counts))
end

function _center_rows(matrix::AbstractMatrix{<:Real})
    centered = Matrix{Float64}(matrix)
    centered .-= mean(centered, dims=2)
    return centered
end

function _gene_correlation(expression::AbstractMatrix{<:Real})
    centered = _center_rows(expression)
    correlation = cor(permutedims(centered))
    correlation[.!isfinite.(correlation)] .= 0.0
    for i in axes(correlation, 1)
        correlation[i, i] = 1.0
    end
    return correlation
end

_gene_correlation(counts::CountMatrix) = _gene_correlation(_dense_expression(counts))

function _scale_free_fit(degrees::AbstractVector{<:Real})
    positive = collect(Float64.(degrees[degrees .> 0]))
    length(positive) < 2 && return NaN
    rounded = unique(round.(positive; digits=4))
    frequencies = [count(x -> isapprox(x, value; atol=1e-4), positive) / length(positive) for value in rounded]
    valid = frequencies .> 0
    count(valid) < 2 && return NaN
    x = log.(rounded[valid])
    y = log.(frequencies[valid])
    design = hcat(ones(length(x)), x)
    coefficients = design \ y
    residuals = y .- design * coefficients
    total = y .- mean(y)
    isapprox(sum(abs2, total), 0.0) && return NaN
    r2 = max(0.0, 1.0 - sum(abs2, residuals) / sum(abs2, total))
    return sign(coefficients[2]) * r2
end

function pick_soft_threshold(counts::CountMatrix; powers::AbstractVector{<:Integer}=2:12, multi_thread::Bool=true)
    expression = _dense_expression(counts)
    correlation = _gene_correlation(expression)
    power_values = collect(Int.(powers))
    fit = zeros(Float64, length(power_values))
    connectivity = zeros(Float64, length(power_values))
    if multi_thread && length(power_values) > 1 && Threads.nthreads() > 1
        @threads for index in eachindex(power_values)
            power = power_values[index]
            adjacency = abs.(correlation) .^ power
            adjacency[diagind(adjacency)] .= 0.0
            degrees = vec(sum(adjacency, dims=2))
            candidate = _scale_free_fit(degrees)
            fit[index] = isfinite(candidate) ? candidate : -Inf
            connectivity[index] = mean(degrees)
        end
    else
        for index in eachindex(power_values)
            power = power_values[index]
            adjacency = abs.(correlation) .^ power
            adjacency[diagind(adjacency)] .= 0.0
            degrees = vec(sum(adjacency, dims=2))
            candidate = _scale_free_fit(degrees)
            fit[index] = isfinite(candidate) ? candidate : -Inf
            connectivity[index] = mean(degrees)
        end
    end
    best = power_values[argmax(fit)]
    result = SoftThresholdResult(power_values, fit, connectivity, best)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, _register_systemsbio_result!(_ctx, result, "pick_soft_threshold"; parents=provenance_parent_ids(counts), parameters=(powers=power_values, best_power=best)), "pick_soft_threshold")
end

function _tom_matrix(adjacency::AbstractMatrix{<:Real})
    weighted = Matrix{Float64}(adjacency)
    # WGCNA weighted TOM uses matrix multiplication for shared neighborhood overlap:
    # l_ij = sum_u a_iu * a_ju
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

function module_eigengenes(expression::AbstractMatrix{<:Real}, modules::AbstractVector{<:Integer})
    expression = permutedims(Matrix{Float64}(expression))
    module_ids = sort(unique(Int.(modules)))
    eigengenes = zeros(Float64, size(expression, 1), length(module_ids))
    for (column, module_id) in enumerate(module_ids)
        indices = findall(==(module_id), modules)
        module_matrix = expression[:, indices]
        module_matrix .-= mean(module_matrix, dims=1)
        decomposition = svd(module_matrix; full=false)
        eigengenes[:, column] .= decomposition.U[:, 1] .* decomposition.S[1]
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, eigengenes, "module_eigengenes")
end

module_eigengenes(counts::CountMatrix, modules::AbstractVector{<:Integer}) = module_eigengenes(_dense_expression(counts), modules)

function find_modules(counts::CountMatrix; power::Union{Nothing,Integer}=nothing, tom_threshold::Real=0.15, edge_threshold::Real=0.1, min_module_size::Integer=2)
    selected_power = power === nothing ? pick_soft_threshold(counts).best_power : Int(power)
    expression = _dense_expression(counts)
    correlation = _gene_correlation(expression)
    adjacency = abs.(correlation) .^ selected_power
    adjacency[diagind(adjacency)] .= 0.0
    tom = _tom_matrix(adjacency)
    modules = _module_assignments(tom; threshold=tom_threshold, min_module_size=min_module_size)
    eigengenes = module_eigengenes(expression, modules)
    expression = permutedims(Matrix{Float64}(expression))
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
    result = GeneNetwork(graph, Dict(gene => index for (index, gene) in enumerate(counts.gene_ids)), String.(counts.gene_ids), connectivity, modules)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, _register_systemsbio_result!(_ctx, result, "find_modules"; parents=provenance_parent_ids(counts), parameters=(power=selected_power, n_modules=maximum(modules))), "find_modules")
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

"""
    estimate_limma_hyperparameters(S2, V)

Estimate limma-style empirical-Bayes prior parameters `(s0_sq, d0)` from
gene-wise raw residual variances `S2` and expected variance scaling terms `V`.
"""
function _natural_cubic_spline_eval(x_knots::Vector{Float64}, y_knots::Vector{Float64}, x_eval::Vector{Float64})
    knot_count = length(x_knots)
    knot_count == length(y_knots) || throw(ArgumentError("x_knots and y_knots must have the same length"))
    knot_count >= 2 || return fill(y_knots[1], length(x_eval))

    h = diff(x_knots)
    any(h .<= 0) && return fill(median(y_knots), length(x_eval))

    second = zeros(Float64, knot_count)
    if knot_count > 2
        lower = h[2:end-1]
        diagonal = 2.0 .* (h[1:end-1] .+ h[2:end])
        upper = h[2:end-1]
        rhs = 6.0 .* ((y_knots[3:end] .- y_knots[2:end-1]) ./ h[2:end] .- (y_knots[2:end-1] .- y_knots[1:end-2]) ./ h[1:end-1])
        second[2:end-1] .= Tridiagonal(lower, diagonal, upper) \ rhs
    end

    out = similar(x_eval)
    for (index, x) in pairs(x_eval)
        segment = if x <= x_knots[1]
            1
        elseif x >= x_knots[end]
            knot_count - 1
        else
            clamp(searchsortedlast(x_knots, x), 1, knot_count - 1)
        end

        h_segment = x_knots[segment + 1] - x_knots[segment]
        a = (x_knots[segment + 1] - x) / h_segment
        b = (x - x_knots[segment]) / h_segment
        out[index] = a * y_knots[segment] + b * y_knots[segment + 1] + ((a^3 - a) * second[segment] + (b^3 - b) * second[segment + 1]) * (h_segment^2 / 6.0)
    end
    return out
end

function _robust_spline_trend(x::Vector{Float64}, y::Vector{Float64}; min_bins::Int=20, max_bins::Int=200)
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y must have the same length"))
    n == 0 && return Float64[]

    x_span = maximum(x) - minimum(x)
    if n < 6 || x_span <= sqrt(eps(Float64))
        return fill(median(y), n)
    end

    order = sortperm(x)
    x_sorted = x[order]
    y_sorted = y[order]

    nbins = clamp(round(Int, sqrt(n)), min_bins, max_bins)
    nbins = min(nbins, n)
    edges = round.(Int, range(1, n + 1; length=nbins + 1))

    knot_x = Float64[]
    knot_y = Float64[]
    for bin in 1:nbins
        lo = edges[bin]
        hi = edges[bin + 1] - 1
        lo <= hi || continue
        push!(knot_x, median(@view x_sorted[lo:hi]))
        push!(knot_y, median(@view y_sorted[lo:hi]))
    end

    unique_knot_x = Float64[]
    unique_knot_y = Float64[]
    for idx in eachindex(knot_x)
        if isempty(unique_knot_x) || knot_x[idx] > unique_knot_x[end] + sqrt(eps(Float64))
            push!(unique_knot_x, knot_x[idx])
            push!(unique_knot_y, knot_y[idx])
        end
    end

    trend_sorted = if length(unique_knot_x) >= 2
        _natural_cubic_spline_eval(unique_knot_x, unique_knot_y, x_sorted)
    else
        fill(median(y_sorted), n)
    end

    trend = similar(y)
    trend[order] = trend_sorted
    return trend
end

function _robust_spline_predict(x::Vector{Float64}, y::Vector{Float64}, x_eval::Vector{Float64}; min_bins::Int=20, max_bins::Int=200)
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y must have the same length"))
    isempty(x_eval) && return Float64[]
    n == 0 && return fill(0.0, length(x_eval))

    x_span = maximum(x) - minimum(x)
    if n < 6 || x_span <= sqrt(eps(Float64))
        return fill(median(y), length(x_eval))
    end

    order = sortperm(x)
    x_sorted = x[order]
    y_sorted = y[order]

    nbins = clamp(round(Int, sqrt(n)), min_bins, max_bins)
    nbins = min(nbins, n)
    edges = round.(Int, range(1, n + 1; length=nbins + 1))

    knot_x = Float64[]
    knot_y = Float64[]
    for bin in 1:nbins
        lo = edges[bin]
        hi = edges[bin + 1] - 1
        lo <= hi || continue
        push!(knot_x, median(@view x_sorted[lo:hi]))
        push!(knot_y, median(@view y_sorted[lo:hi]))
    end

    unique_knot_x = Float64[]
    unique_knot_y = Float64[]
    for idx in eachindex(knot_x)
        if isempty(unique_knot_x) || knot_x[idx] > unique_knot_x[end] + sqrt(eps(Float64))
            push!(unique_knot_x, knot_x[idx])
            push!(unique_knot_y, knot_y[idx])
        end
    end

    if length(unique_knot_x) < 2
        return fill(median(y), length(x_eval))
    end

    return _natural_cubic_spline_eval(unique_knot_x, unique_knot_y, x_eval)
end

_logmdigamma(x::Float64) = log(x) - digamma(x)

function _trigamma_inverse(x::Float64)
    if !isfinite(x)
        return NaN
    elseif x < 0
        return NaN
    elseif x > 1e7
        return inv(sqrt(x))
    elseif x < 1e-6
        return inv(x)
    end

    y = 0.5 + inv(x)
    for _ in 1:50
        tri = trigamma(y)
        second = polygamma(2, y)
        if !isfinite(tri) || !isfinite(second) || second == 0
            break
        end
        diff = tri * (1 - tri / x) / second
        y += diff
        if abs(diff) <= max(1e-8 * abs(y), 1e-12)
            break
        end
    end
    return y
end

function _fit_fdist(variances::Vector{Float64}, df1::Float64; covariate::Union{Nothing,Vector{Float64}}=nothing)
    n = length(variances)
    n == 0 && return (Float64[], NaN)
    n == 1 && return (copy(variances), 0.0)

    x = copy(variances)
    x .= max.(x, 0.0)
    med = median(x)
    if med == 0.0
        med = 1.0
    end
    x .= max.(x, 1e-5 * med)

    z = log.(x)
    e = z .+ _logmdigamma(df1 / 2)

    emean = if covariate === nothing
        fill(mean(e), n)
    else
        _robust_spline_trend(covariate, e)
    end

    centered = e .- emean
    evar = n > 1 ? sum(abs2, centered) / (n - 1) : 0.0
    evar -= trigamma(df1 / 2)

    if evar > 0
        df2 = 2 * _trigamma_inverse(evar)
        scale = exp.(emean .- _logmdigamma(df2 / 2))
        return (scale, df2)
    else
        df2 = Inf
        if covariate === nothing
            return (fill(mean(x), n), df2)
        else
            return (exp.(emean), df2)
        end
    end
end

function _squeeze_var(variances::Vector{Float64}, df_residual::Float64, var_prior::Vector{Float64}, df_prior::Float64)
    if isfinite(df_prior)
        return (df_residual .* variances .+ df_prior .* var_prior) ./ (df_residual + df_prior)
    else
        return copy(var_prior)
    end
end

function estimate_limma_hyperparameters(S2::AbstractVector{<:Real}, covariate::Union{Nothing,AbstractVector{<:Real}}=nothing; return_trend::Bool=false, df_residual::Real=1.0)
    covariate === nothing || length(S2) == length(covariate) || throw(ArgumentError("S2 and covariate must have the same length"))
    df_residual > 0 || throw(ArgumentError("df_residual must be positive"))

    valid = isfinite.(S2) .& (S2 .> 0)
    if covariate !== nothing
        valid .&= isfinite.(covariate)
    end

    S2_fit = Float64.(S2[valid])
    cov_fit = covariate === nothing ? nothing : Float64.(covariate[valid])

    if length(S2_fit) < 3
        fallback = isempty(S2_fit) ? 1.0 : max(median(S2_fit), eps(Float64))
        if return_trend
            return (fallback, 0.0, fill(fallback, length(S2)))
        end
        return (fallback, 0.0)
    end

    prior_trend, d0 = _fit_fdist(S2_fit, Float64(df_residual); covariate=cov_fit)
    d0 = isfinite(d0) ? max(Float64(d0), 0.0) : Inf
    s0_sq = max(median(prior_trend), eps(Float64))

    if return_trend
        prior_full = fill(s0_sq, length(S2))
        prior_full[valid] = prior_trend
        return (s0_sq, d0, prior_full)
    end
    return (s0_sq, d0)
end

"""
    limma_moderate_ttest(X, y, coef_idx; s0_sq=0.0, d0=0.0)

Compute a limma-style moderated t-statistic for one response vector `y` and
coefficient index `coef_idx` in design matrix `X`.

Returns `(t_moderated, s_posterior, s_raw)`.
"""
function limma_moderate_ttest(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}, coef_idx::Integer; s0_sq::Real=0.0, d0::Real=0.0)
    n, p = size(X)
    length(y) == n || throw(ArgumentError("y must have one value per design row"))
    1 <= coef_idx <= p || throw(ArgumentError("coef_idx must be between 1 and $p"))

    Xf = Matrix{Float64}(X)
    yf = Float64.(y)
    xtx_inv = pinv(Xf' * Xf)
    beta = xtx_inv * (Xf' * yf)

    residuals = yf .- Xf * beta
    df_residual = max(n - rank(Xf), 1)
    s2_raw = sum(abs2, residuals) / df_residual

    prior_var = s0_sq > 0 ? Float64(s0_sq) : s2_raw
    prior_df = max(Float64(d0), 0.0)
    s2_posterior = _squeeze_var([s2_raw], Float64(df_residual), [prior_var], prior_df)[1]

    v_g = max(xtx_inv[coef_idx, coef_idx], eps(Float64))
    t_moderated = beta[coef_idx] / sqrt(max(s2_posterior * v_g, eps(Float64)))
    return (t_moderated, sqrt(max(s2_posterior, 0.0)), sqrt(max(s2_raw, 0.0)))
end

function _weighted_linear_model(X::Matrix{Float64}, y::AbstractVector{<:Real}, w::AbstractVector{<:Real}, df_residual::Int)
    length(y) == size(X, 1) || throw(ArgumentError("y must have one value per design row"))
    length(w) == size(X, 1) || throw(ArgumentError("weights must have one value per design row"))

    w_clean = max.(Float64.(w), eps(Float64))
    sqrt_w = sqrt.(w_clean)
    Xw = X .* sqrt_w
    yw = y .* sqrt_w

    xtx_inv = pinv(Xw' * Xw)
    beta = xtx_inv * (Xw' * yw)
    residuals = y .- X * beta
    s2 = sum(w_clean .* abs2.(residuals)) / max(df_residual, 1)
    stdev_unscaled = sqrt.(max.(diag(xtx_inv), eps(Float64)))
    return beta, s2, stdev_unscaled
end

function voom_transform(counts::CountMatrix, design::AbstractMatrix{<:Real}; normalization_method::Symbol=:tmm, pseudocount::Float64=0.5, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    design_matrix = Matrix{Float64}(design)
    n_genes = size(counts.counts, 1)
    n_samples = size(counts.counts, 2)
    size(design_matrix, 1) == n_samples || throw(ArgumentError("design must have one row per sample"))
    pseudocount > 0 || throw(ArgumentError("pseudocount must be positive"))

    lib_sizes = Float64.(vec(sum(counts.counts, dims=1)))
    all(lib_sizes .> 0) || throw(ArgumentError("all samples must have positive library sizes"))
    norm_factors = calc_norm_factors(counts; method=normalization_method)
    effective_lib_sizes = lib_sizes .* norm_factors

    raw_counts = Matrix{Float64}(counts.counts)
    log_cpm = log2.((raw_counts .+ pseudocount) ./ reshape(effective_lib_sizes .+ 1.0, 1, :) .* 1e6)

    response = permutedims(log_cpm) # samples x genes
    xtx_inv = pinv(design_matrix' * design_matrix)
    beta = xtx_inv * (design_matrix' * response)
    fitted = design_matrix * beta
    residuals = response .- fitted

    df_residual = max(n_samples - rank(design_matrix), 1)
    sigma = sqrt.(vec(sum(abs2, residuals; dims=1)) ./ df_residual)
    amean = vec(mean(log_cpm, dims=2))
    sx = amean .+ mean(log2.(effective_lib_sizes .+ 1.0)) .- log2(1e6)

    fitted_values = permutedims(fitted) # genes x samples
    fitted_cpm = 2 .^ fitted_values
    fitted_count = 1e-6 .* fitted_cpm .* reshape(effective_lib_sizes .+ 1.0, 1, :)
    fitted_logcount = log2.(max.(fitted_count, eps(Float64)))

    allzero = vec(sum(raw_counts, dims=2)) .== 0.0
    trend_keep = .!allzero .& isfinite.(sx) .& isfinite.(sigma)
    if count(trend_keep) < 3
        finite_sigma = sigma[isfinite.(sigma)]
        base_sd = isempty(finite_sigma) ? 1.0 : max(median(finite_sigma), sqrt(eps(Float64)))
        weights = fill(1.0 / (base_sd^4), size(log_cpm))
        return VoomResult(log_cpm, weights, fitted_logcount, effective_lib_sizes, design_matrix)
    end

    trend_sd = _robust_spline_predict(Float64.(sx[trend_keep]), Float64.(sigma[trend_keep]), vec(fitted_logcount))
    trend_sd = max.(trend_sd, sqrt(eps(Float64)))
    weights = reshape(1.0 ./ (trend_sd .^ 4), size(log_cpm))
    result = VoomResult(log_cpm, weights, fitted_logcount, effective_lib_sizes, design_matrix)


    return _register_systemsbio_result!(_ctx, result, "voom_transform"; parents=provenance_parent_ids(counts), parameters=(normalization_method=normalization_method, pseudocount=pseudocount, n_genes=n_genes, n_samples=n_samples))
end

function voom_transform(counts::AbstractMatrix{<:Integer}, design::AbstractMatrix{<:Real}; gene_ids::Union{Nothing,AbstractVector{<:String}}=nothing, sample_ids::Union{Nothing,AbstractVector{<:String}}=nothing, kwargs...)
    genes = gene_ids === nothing ? ["gene$(index)" for index in 1:size(counts, 1)] : String.(gene_ids)
    samples = sample_ids === nothing ? ["sample$(index)" for index in 1:size(counts, 2)] : String.(sample_ids)
    cm = CountMatrix(counts, genes, samples)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, voom_transform(cm, design; _ctx=_ctx, kwargs...), "voom_transform")
end

"""
    voom_with_quality_weights(counts, design; kwargs...)

Two-pass voom that estimates sample-level quality weights `q_i` and returns
combined observation weights `w_gi * q_i`.
"""
function voom_with_quality_weights(counts::CountMatrix, design::AbstractMatrix{<:Real}; normalization_method::Symbol=:tmm, pseudocount::Float64=0.5, weight_floor::Float64=1e-6, weight_cap::Float64=1e10, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    base_voom = voom_transform(counts, design; normalization_method=normalization_method, pseudocount=pseudocount)
    X = base_voom.design
    n_genes, n_samples = size(base_voom.log_cpm)
    df_residual = max(n_samples - rank(X), 1)

    sample_mse = zeros(Float64, n_samples)
    for gene in 1:n_genes
        y = vec(@view base_voom.log_cpm[gene, :])
        w = vec(@view base_voom.weights[gene, :])
        beta, _, _ = _weighted_linear_model(X, y, w, df_residual)
        residuals = y .- X * beta
        sample_mse .+= (sqrt.(max.(w, eps(Float64))) .* residuals) .^ 2
    end
    sample_mse ./= max(n_genes, 1)

    sample_weights = 1.0 ./ max.(sample_mse, eps(Float64))
    sample_weights ./= exp(mean(log.(sample_weights)))
    sample_weights = clamp.(sample_weights, weight_floor, weight_cap)

    combined_weights = clamp.(base_voom.weights .* reshape(sample_weights, 1, :), weight_floor, weight_cap)
    adjusted_voom = VoomResult(base_voom.log_cpm, combined_weights, base_voom.fitted_logcount, base_voom.effective_lib_sizes, base_voom.design)
    result = (voom=adjusted_voom, sample_weights=sample_weights)


    return _register_systemsbio_result!(_ctx, result, "voom_with_quality_weights"; parents=provenance_parent_ids(counts), parameters=(normalization_method=normalization_method, pseudocount=pseudocount, n_samples=length(sample_weights)))
end

function voom_with_quality_weights(counts::AbstractMatrix{<:Integer}, design::AbstractMatrix{<:Real}; gene_ids::Union{Nothing,AbstractVector{<:String}}=nothing, sample_ids::Union{Nothing,AbstractVector{<:String}}=nothing, kwargs...)
    genes = gene_ids === nothing ? ["gene$(index)" for index in 1:size(counts, 1)] : String.(gene_ids)
    samples = sample_ids === nothing ? ["sample$(index)" for index in 1:size(counts, 2)] : String.(sample_ids)
    cm = CountMatrix(counts, genes, samples)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, voom_with_quality_weights(cm, design; _ctx=_ctx, kwargs...), "voom_with_quality_weights")
end

function _ebayes_raw_fit(coefficients::Matrix{Float64}, stdev_unscaled::Matrix{Float64}, sigma2::Vector{Float64}, df_residual::Float64, gene_ids::Vector{String}, coefficient_names::Vector{String}, base_mean::Vector{Float64})
    return (
        coefficients=coefficients,
        stdev_unscaled=stdev_unscaled,
        sigma2=sigma2,
        df_residual=df_residual,
        gene_ids=gene_ids,
        coefficient_names=coefficient_names,
        base_mean=base_mean)
end

function _winsorize(values::AbstractVector{<:Real}, tails::Tuple{<:Real,<:Real})
    lower_tail, upper_tail = Float64(tails[1]), Float64(tails[2])
    (0.0 <= lower_tail < 0.5 && 0.0 <= upper_tail < 0.5 && lower_tail + upper_tail < 1.0) ||
        throw(ArgumentError("winsor_tail_p must satisfy 0 <= tails < 0.5 and sum < 1"))

    data = Float64.(values)
    isempty(data) && return data
    lo = quantile(data, lower_tail)
    hi = quantile(data, 1.0 - upper_tail)
    return clamp.(data, lo, hi)
end

"""
    eBayes(fit; robust=false, winsor_tail_p=(0.05, 0.1), trend=true)

Apply limma-style empirical-Bayes moderation to a linear-model fit payload.
"""
function eBayes(fit;
                robust::Bool=false,
                winsor_tail_p::Tuple{<:Real,<:Real}=(0.05, 0.1),
                trend::Bool=true,
                proportion::Real=0.01,
                return_hyperparameters::Bool=false)
    hasproperty(fit, :coefficients) || throw(ArgumentError("fit must provide coefficients"))
    hasproperty(fit, :stdev_unscaled) || throw(ArgumentError("fit must provide stdev_unscaled"))
    hasproperty(fit, :sigma2) || throw(ArgumentError("fit must provide sigma2"))
    hasproperty(fit, :df_residual) || throw(ArgumentError("fit must provide df_residual"))
    hasproperty(fit, :gene_ids) || throw(ArgumentError("fit must provide gene_ids"))
    hasproperty(fit, :coefficient_names) || throw(ArgumentError("fit must provide coefficient_names"))
    hasproperty(fit, :base_mean) || throw(ArgumentError("fit must provide base_mean"))

    coefficients = Matrix{Float64}(fit.coefficients)
    stdev_unscaled = Matrix{Float64}(fit.stdev_unscaled)
    sigma2 = max.(Float64.(fit.sigma2), eps(Float64))
    df_residual = max(Float64(fit.df_residual), 1.0)
    gene_ids = String.(fit.gene_ids)
    coefficient_names = String.(fit.coefficient_names)
    base_mean = Float64.(fit.base_mean)

    prior_source = copy(sigma2)
    use_robust = robust && length(prior_source) >= 10
    if use_robust
        log_variance = log.(prior_source)
        log_variance = _winsorize(log_variance, winsor_tail_p)
        prior_source .= exp.(log_variance)
    end

    covariate = trend ? base_mean : nothing
    s0_sq, d0, prior_trend = estimate_limma_hyperparameters(prior_source, covariate; return_trend=true, df_residual=df_residual)
    moderated_variance = _squeeze_var(sigma2, df_residual, prior_trend, Float64(d0))
    moderated_variance = max.(moderated_variance, eps(Float64))

    n_genes, n_coef = size(coefficients)
    standard_errors = sqrt.(max.(moderated_variance, eps(Float64))) .* stdev_unscaled
    t_statistics = coefficients ./ max.(standard_errors, eps(Float64))
    df_total = min(df_residual + d0, n_genes * df_residual)
    t_distribution = TDist(max(df_total, 1.0))
    pvalues = 2.0 .* ccdf.(Ref(t_distribution), abs.(t_statistics))

    qvalues = similar(pvalues)
    for coef_index in 1:n_coef
        qvalues[:, coef_index] .= benjamini_hochberg(pvalues[:, coef_index])
    end

    result = LimmaResult(coefficients, standard_errors, t_statistics, pvalues, qvalues, moderated_variance, base_mean, gene_ids, coefficient_names)
    if return_hyperparameters
        return (fit=result, s0_sq=s0_sq, d0=d0, prior_trend=prior_trend, robust_used=use_robust, proportion=Float64(proportion))
    end
    return result
end

"""
    contrasts_fit(fit, contrast_matrix; coefficient_names=nothing)

Project limma coefficients onto a contrast matrix using pseudoinverse-safe
variance propagation. Redundant (non-estimable) contrast columns are returned
as `NaN` statistics.
"""
function contrasts_fit(fit::LimmaResult, contrast_matrix::AbstractMatrix{<:Real}; coefficient_names::Union{Nothing,AbstractVector{<:String}}=nothing)
    C = Matrix{Float64}(contrast_matrix)
    ncoef = size(fit.coefficients, 2)
    size(C, 1) == ncoef || throw(ArgumentError("contrast_matrix must have $ncoef rows"))

    ncontrasts = size(C, 2)
    raw_names = coefficient_names === nothing ? ["contrast_$(index)" for index in 1:ncontrasts] : String.(coefficient_names)
    length(raw_names) == ncontrasts || throw(ArgumentError("coefficient_names must match number of contrast columns"))

    estimable = trues(ncontrasts)
    if rank(C) < ncontrasts
        decomposition = qr(C, ColumnNorm())
        independent = Set(Int.(decomposition.p[1:rank(C)]))
        for index in 1:ncontrasts
            estimable[index] = index in independent
        end
    end

    n_genes = size(fit.coefficients, 1)
    coefficients = fill(NaN, n_genes, ncontrasts)
    standard_errors = fill(NaN, n_genes, ncontrasts)
    t_statistics = fill(NaN, n_genes, ncontrasts)
    pvalues = fill(NaN, n_genes, ncontrasts)

    coef_var = fit.standard_errors .^ 2
    for contrast_index in 1:ncontrasts
        estimable[contrast_index] || continue
        c = C[:, contrast_index]
        for gene in 1:n_genes
            beta_row = vec(@view fit.coefficients[gene, :])
            any(!isfinite, beta_row) && continue
            coefficients[gene, contrast_index] = dot(beta_row, c)

            variance = dot(c .^ 2, vec(@view coef_var[gene, :]))
            variance = max(variance, eps(Float64))
            se = sqrt(variance)
            standard_errors[gene, contrast_index] = se
            t_statistics[gene, contrast_index] = coefficients[gene, contrast_index] / se
            pvalues[gene, contrast_index] = 2.0 * ccdf(Normal(), abs(t_statistics[gene, contrast_index]))
        end
    end

    qvalues = fill(NaN, n_genes, ncontrasts)
    for contrast_index in 1:ncontrasts
        finite_mask = isfinite.(pvalues[:, contrast_index])
        if any(finite_mask)
            adjusted = benjamini_hochberg(pvalues[finite_mask, contrast_index])
            qvalues[finite_mask, contrast_index] .= adjusted
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, LimmaResult(coefficients, standard_errors, t_statistics, pvalues, qvalues, copy(fit.moderated_variance), copy(fit.base_mean), copy(fit.gene_ids), raw_names), "contrasts_fit")
end

function _duplicatecorrelation_objective(residuals::Matrix{Float64}, block_indices::Vector{Vector{Int}}, rho::Float64)
    total = 0.0
    for indices in block_indices
        m = length(indices)
        inv_one_minus = 1.0 / max(1.0 - rho, eps(Float64))
        denom = max(1.0 + (m - 1) * rho, eps(Float64))
        logdet = (m - 1) * log(max(1.0 - rho, eps(Float64))) + log(denom)

        view_block = @view residuals[:, indices]
        block_sum = vec(sum(view_block; dims=2))
        block_sq = vec(sum(abs2, view_block; dims=2))
        quad = inv_one_minus .* (block_sq .- (rho / denom) .* (block_sum .^ 2))
        total += 0.5 * sum(logdet .+ quad)
    end
    return total
end

function _bounded_golden_section_minimize(objective::Function, lower::Float64, upper::Float64; tol::Float64=1e-6, maxit::Int=200)
    phi = (sqrt(5.0) - 1.0) / 2.0
    a, b = lower, upper
    c = b - phi * (b - a)
    d = a + phi * (b - a)
    fc = objective(c)
    fd = objective(d)
    for _ in 1:maxit
        if abs(b - a) <= tol
            break
        end
        if fc < fd
            b, d, fd = d, c, fc
            c = b - phi * (b - a)
            fc = objective(c)
        else
            a, c, fc = c, d, fd
            d = a + phi * (b - a)
            fd = objective(d)
        end
    end
    return (a + b) / 2.0
end

"""
    duplicateCorrelation(expression, design, block)

Estimate consensus intra-block correlation for repeated measures using bounded
optimization of a compound-symmetry likelihood.
"""
function duplicateCorrelation(expression::AbstractMatrix{<:Real}, design::AbstractMatrix{<:Real}, block::AbstractVector)
    matrix = Matrix{Float64}(expression)
    n_genes, n_samples = size(matrix)
    size(design, 1) == n_samples || throw(ArgumentError("design must have one row per sample"))
    length(block) == n_samples || throw(ArgumentError("block must have one value per sample"))

    block_values = collect(block)
    counts = Dict{Any,Int}()
    for value in block_values
        counts[value] = get(counts, value, 0) + 1
    end

    keep_samples = [counts[value] > 1 for value in block_values]
    if count(keep_samples) < 2
        return (consensus_correlation=0.0, correlation=0.0, at_boundary=false)
    end

    sample_indices = findall(keep_samples)
    block_filtered = block_values[sample_indices]
    design_filtered = Matrix{Float64}(design[sample_indices, :])
    matrix_filtered = matrix[:, sample_indices]

    block_map = Dict{Any,Vector{Int}}()
    for (position, value) in enumerate(block_filtered)
        push!(get!(block_map, value, Int[]), position)
    end
    block_indices = [indices for indices in values(block_map) if length(indices) > 1]
    isempty(block_indices) && return (consensus_correlation=0.0, correlation=0.0, at_boundary=false)

    X = Matrix{Float64}(design_filtered)
    projector = pinv(X)
    fitted = X * (projector * permutedims(matrix_filtered))
    residuals = permutedims(matrix_filtered) .- fitted
    residuals = permutedims(residuals)

    max_block = maximum(length.(block_indices))
    lower = -1.0 / (max_block - 1) + 1e-6
    upper = 1.0 / (max_block - 1) - 1e-6
    if lower >= upper
        return (consensus_correlation=0.0, correlation=0.0, at_boundary=false)
    end

    objective(rho) = _duplicatecorrelation_objective(residuals, block_indices, rho)
    rho_hat = _bounded_golden_section_minimize(objective, lower, upper; tol=1e-6, maxit=300)
    rho_hat = clamp(rho_hat, lower, upper)
    at_boundary = abs(rho_hat - lower) < 1e-4 || abs(rho_hat - upper) < 1e-4
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (consensus_correlation=rho_hat, correlation=rho_hat, at_boundary=at_boundary), "duplicateCorrelation")
end

function _batch_design_matrix(batch::AbstractVector)
    levels = unique(batch)
    matrix = zeros(Float64, length(batch), length(levels))
    for (index, level) in enumerate(levels)
        for sample in eachindex(batch)
            matrix[sample, index] = batch[sample] == level ? 1.0 : 0.0
        end
    end
    return matrix
end

"""
    remove_batch_effect_limma(expression, batch; covariates=nothing)

Remove additive batch effects using pseudoinverse projection, while optionally
preserving biological covariates.
"""
function remove_batch_effect_limma(expression::AbstractMatrix{<:Real}, batch::AbstractVector; covariates::Union{Nothing,AbstractVector,AbstractMatrix}=nothing)
    Y = permutedims(Matrix{Float64}(expression)) # samples x features
    n_samples = size(Y, 1)
    length(batch) == n_samples || throw(ArgumentError("batch must have one value per sample"))

    batch_design = _batch_design_matrix(batch)
    if covariates === nothing
        residual = Y
        preserved = zeros(Float64, size(Y))
    else
        cov_matrix = covariates isa AbstractVector ? reshape(Float64.(covariates), :, 1) : Matrix{Float64}(covariates)
        size(cov_matrix, 1) == n_samples || throw(ArgumentError("covariates must have one row per sample"))
        centered_covariates = cov_matrix .- mean(cov_matrix, dims=1)
        bio_design = hcat(ones(Float64, n_samples), centered_covariates)
        preserved = bio_design * (pinv(bio_design) * Y)
        residual = Y .- preserved
    end

    batch_effect = batch_design * (pinv(batch_design) * residual)
    corrected = residual .- batch_effect .+ preserved
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, permutedims(corrected), "remove_batch_effect_limma")
end

function limma_fit(counts::CountMatrix, design::AbstractMatrix{<:Real}; coefficient_names::Union{Nothing,AbstractVector{<:String}}=nothing, normalization_method::Symbol=:tmm, pseudocount::Float64=0.5, use_voom::Bool=true, trend::Bool=true, robust::Bool=false, winsor_tail_p::Tuple{<:Real,<:Real}=(0.05, 0.1), proportion::Real=0.01, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    design_matrix = Matrix{Float64}(design)
    n_samples, n_coef = size(design_matrix)
    size(counts.counts, 2) == n_samples || throw(ArgumentError("design must have one row per sample"))
    pseudocount > 0 || throw(ArgumentError("pseudocount must be positive"))

    log_cpm = Matrix{Float64}(undef, size(counts.counts, 1), size(counts.counts, 2))
    observation_weights = Matrix{Float64}(undef, size(counts.counts, 1), size(counts.counts, 2))
    if use_voom
        voom = voom_transform(counts, design_matrix; normalization_method=normalization_method, pseudocount=pseudocount, _ctx=_ctx)
        log_cpm .= voom.log_cpm
        observation_weights .= voom.weights
    else
        lib_sizes = Float64.(vec(sum(counts.counts, dims=1)))
        all(lib_sizes .> 0) || throw(ArgumentError("all samples must have positive library sizes"))
        norm_factors = calc_norm_factors(counts; method=normalization_method)
        effective_lib_sizes = lib_sizes .* norm_factors
        raw_counts = Matrix{Float64}(counts.counts)
        log_cpm .= log2.(raw_counts .+ pseudocount) .- log2.(reshape(effective_lib_sizes ./ 1e6 .+ pseudocount, 1, :))
        fill!(observation_weights, 1.0)
    end

    gene_means = vec(mean(log_cpm, dims=2))

    n_genes = size(log_cpm, 1)
    coefficient_matrix = zeros(Float64, n_genes, n_coef)
    stdev_unscaled = zeros(Float64, n_genes, n_coef)
    s2_all = zeros(Float64, n_genes)

    df_residual = max(n_samples - rank(design_matrix), 1)
    for gene in 1:n_genes
        y = vec(@view log_cpm[gene, :])
        w = vec(@view observation_weights[gene, :])
        beta_gene, s2_gene, stdev_gene = _weighted_linear_model(design_matrix, y, w, df_residual)
        coefficient_matrix[gene, :] .= beta_gene
        stdev_unscaled[gene, :] .= stdev_gene
        s2_all[gene] = s2_gene
    end

    names = coefficient_names === nothing ? ["coef_$(index)" for index in 1:n_coef] : String.(coefficient_names)
    raw_fit = _ebayes_raw_fit(coefficient_matrix, stdev_unscaled, s2_all, Float64(df_residual), String.(counts.gene_ids), names, gene_means)
    moderated = eBayes(raw_fit; robust=robust, winsor_tail_p=winsor_tail_p, trend=trend, proportion=proportion)
    _register_systemsbio_result!(_ctx, moderated, "limma_fit"; parents=provenance_parent_ids(counts), parameters=(use_voom=use_voom, trend=trend, robust=robust, n_genes=n_genes))
    return moderated
end

function limma_fit(counts::AbstractMatrix{<:Integer}, design::AbstractMatrix{<:Real}; coefficient_names::Union{Nothing,AbstractVector{<:String}}=nothing, gene_ids::Union{Nothing,AbstractVector{<:String}}=nothing, sample_ids::Union{Nothing,AbstractVector{<:String}}=nothing, kwargs...)
    genes = gene_ids === nothing ? ["gene$(index)" for index in 1:size(counts, 1)] : String.(gene_ids)
    samples = sample_ids === nothing ? ["sample$(index)" for index in 1:size(counts, 2)] : String.(sample_ids)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, limma_fit(CountMatrix(counts, genes, samples), design; coefficient_names=coefficient_names, _ctx=_ctx, kwargs...), "limma_fit")
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
            result.qvalues[index, coefficient_index])
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, rows, "limma_deresults")
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

function _gsea_pvalue(ranked_genes::Vector{String}, geneset::Set{String}, observed::Float64; permutations::Integer=1000, seed::Integer=1, multi_thread::Bool=true)
    permutation_scores = zeros(Float64, permutations)
    if multi_thread && permutations > 1 && Threads.nthreads() > 1
        @threads for perm in 1:permutations
            rng = MersenneTwister(seed + perm)
            permutation_scores[perm] = _gsea_score(shuffle(rng, ranked_genes), geneset)
        end
    else
        for perm in 1:permutations
            rng = MersenneTwister(seed + perm)
            permutation_scores[perm] = _gsea_score(shuffle(rng, ranked_genes), geneset)
        end
    end
    return (count(score -> abs(score) >= abs(observed), permutation_scores) + 1) / (permutations + 1)
end

function gsea(results::Vector{DEResult}, database::EnrichmentDatabase; permutations::Integer=1000, min_size::Integer=5, max_size::Integer=500, seed::Integer=1, multi_thread::Bool=true)
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
        pvalue = _gsea_pvalue(ranked_genes, geneset, observed; permutations=permutations, seed=seed, multi_thread=multi_thread)
        overlap = count(gene -> gene in geneset, ranked_genes)
        odds_ratio = overlap == 0 ? 0.0 : (overlap / length(ranked_genes)) / (length(geneset) / max(background_size, 1))
        push!(output, EnrichmentResult(term.id, term.name, term.namespace, overlap, length(geneset), length(ranked_genes), background_size, pvalue, 1.0, odds_ratio, collect(geneset)))
    end
    isempty(output) && return output
    qvalues = benjamini_hochberg([result.pvalue for result in output])
    for (index, result) in enumerate(output)
        output[index] = EnrichmentResult(result.term_id, result.term_name, result.namespace, result.overlap, result.term_size, result.query_size, result.background_size, result.pvalue, qvalues[index], result.odds_ratio, result.genes)
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, output, "gsea")
end

function gsea(result::LimmaResult, database::EnrichmentDatabase; coefficient_index::Integer=2, kwargs...)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, gsea(limma_deresults(result; coefficient_index=coefficient_index), database; kwargs...), "gsea")
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

function infer_network(data::AbstractMatrix{<:Real}; gene_ids::Union{Nothing,AbstractVector{<:String}}=nothing, max_parents::Integer=3, max_iterations::Integer=100, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
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
                # max_parents constrains incoming regulators per child node.
                indegree(graph, child) >= max_parents && continue
                # Adding parent -> child creates a cycle iff parent is already reachable from child.
                has_path(graph, child, parent) && continue
                candidate = copy(graph)
                add_edge!(candidate, parent, child)
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
    result = NetworkInferenceResult(graph, Dict(label => index for (index, label) in enumerate(labels)), labels, best_score)


    return _register_systemsbio_result!(_ctx, result, "infer_network"; parents=provenance_parent_ids(data), parameters=(n_nodes=nvars, max_parents=Int(max_parents), max_iterations=Int(max_iterations)))
end

function infer_network(counts::CountMatrix; kwargs...)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, infer_network(permutedims(_dense_expression(counts)); gene_ids=counts.gene_ids, _ctx=_ctx, kwargs...), "infer_network")
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
    result = MultiOmicsFactorAnalysisResult(
        factors,
        loadings,
        assay_loadings,
        explained_variance,
        provenance_record(
            "MultiOmicsFactorAnalysisResult",
            "SystemsBio/multi_omics_factor_analysis";
            notes=["shared SVD across z-scored assays"],
            parameters=(assay_count=length(assays), sample_count=nsamples, factor_count=used)))
    _ctx = active_provenance_context()
    _register_systemsbio_result!(_ctx, result, "multi_omics_factor_analysis"; parameters=(n_modalities=length(assays), n_factors=used, n_samples=nsamples))


    return result
end

function _impute_missing_columns(x::AbstractMatrix{<:Real})
    y = Matrix{Float64}(x)
    for j in axes(y, 2)
        col = @view y[:, j]
        finite = col[isfinite.(col)]
        fillv = isempty(finite) ? 0.0 : mean(finite)
        for i in eachindex(col)
            isfinite(col[i]) || (col[i] = fillv)
        end
    end
    return y
end

"""
    mofa_plus_integration(assays; n_factors=5)

MOFA+-style integration with per-view missing-value imputation.
"""
function mofa_plus_integration(assays::Vector{<:AbstractMatrix{<:Real}}; n_factors::Int=5)
    cleaned = [_impute_missing_columns(assay) for assay in assays]
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, multi_omics_factor_analysis(cleaned; n_factors=n_factors), "mofa_plus_integration")
end

"""
    cca_integration(x, y; n_components=10, ridge=1e-5)

Canonical correlation analysis between two omics matrices (samples x features).
"""
function cca_integration(x::AbstractMatrix{<:Real}, y::AbstractMatrix{<:Real}; n_components::Int=10, ridge::Real=1e-5)
    size(x, 1) == size(y, 1) || throw(DimensionMismatch("x and y must have same sample count (rows)"))
    X = Matrix{Float64}(x)
    Y = Matrix{Float64}(y)
    X .-= mean(X, dims=1)
    Y .-= mean(Y, dims=1)

    n = size(X, 1)
    Sxx = (X' * X) / max(n - 1, 1) + Float64(ridge) * I
    Syy = (Y' * Y) / max(n - 1, 1) + Float64(ridge) * I
    Sxy = (X' * Y) / max(n - 1, 1)

    Ex = eigen(Symmetric(Sxx))
    Ey = eigen(Symmetric(Syy))
    Wx = Ex.vectors * Diagonal(1.0 ./ sqrt.(max.(Ex.values, eps(Float64)))) * Ex.vectors'
    Wy = Ey.vectors * Diagonal(1.0 ./ sqrt.(max.(Ey.values, eps(Float64)))) * Ey.vectors'

    M = Wx * Sxy * Wy
    U, S, V = svd(M)
    used = min(n_components, length(S))
    x_scores = X * (Wx * U[:, 1:used])
    y_scores = Y * (Wy * V[:, 1:used])
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (x_scores=x_scores, y_scores=y_scores, canonical_correlations=S[1:used]), "cca_integration")
end

"""
    anchor_based_integration(reference, query; n_components=20, k=5)

Seurat-style anchor candidate discovery in a shared latent embedding.
"""
function anchor_based_integration(reference::AbstractMatrix{<:Real}, query::AbstractMatrix{<:Real}; n_components::Int=20, k::Int=5)
    size(reference, 2) == size(query, 2) || throw(DimensionMismatch("reference and query must have same feature count"))
    concat = vcat(Matrix{Float64}(reference), Matrix{Float64}(query))
    concat .-= mean(concat, dims=1)
    fac = svd(concat; full=false)
    used = min(n_components, size(fac.U, 2))
    latent = fac.U[:, 1:used] * Diagonal(fac.S[1:used])

    n_ref = size(reference, 1)
    ref_latent = latent[1:n_ref, :]
    qry_latent = latent[(n_ref + 1):end, :]

    anchor_q = Int[]
    anchor_r = Int[]
    dist = Float64[]

    for q in 1:size(qry_latent, 1)
        d = [(r, sum(abs2, @view(qry_latent[q, :]) .- @view(ref_latent[r, :]))) for r in 1:size(ref_latent, 1)]
        sort!(d; by=last)
        for (r, v) in Iterators.take(d, min(k, length(d)))
            push!(anchor_q, q)
            push!(anchor_r, r)
            push!(dist, sqrt(v))
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, DataFrame(query_index=anchor_q, reference_index=anchor_r, distance=dist), "anchor_based_integration")
end

"""
    totalvi_like_integration(rna_counts, protein_counts; n_latent=16)

totalVI-like latent embedding from joint RNA/protein counts with background
correction on protein channels.
"""
function totalvi_like_integration(rna_counts::AbstractMatrix{<:Real}, protein_counts::AbstractMatrix{<:Real}; n_latent::Int=16, background_quantile::Real=0.1)
    size(rna_counts, 2) == size(protein_counts, 2) || throw(DimensionMismatch("RNA and protein matrices must share cell count in columns"))
    rna = log1p.(Float64.(rna_counts))
    protein = log1p.(Float64.(protein_counts))

    bg = mapslices(col -> quantile(vec(col), clamp(Float64(background_quantile), 0.0, 1.0)), protein; dims=2)
    corrected = max.(protein .- bg, 0.0)

    joint = vcat(rna, corrected)
    cells_by_feat = permutedims(joint)
    cells_by_feat .-= mean(cells_by_feat, dims=1)
    fac = svd(cells_by_feat; full=false)
    used = min(n_latent, size(fac.U, 2))
    latent = fac.U[:, 1:used] * Diagonal(fac.S[1:used])
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (latent=latent, protein_corrected=corrected), "totalvi_like_integration")
end

"""
    causal_grn_inference(expression; top_k=200)

Precision-matrix based causal graph approximation for gene networks.
"""
function causal_grn_inference(expression::AbstractMatrix{<:Real}; top_k::Int=200)
    X = Matrix{Float64}(expression)
    X .-= mean(X, dims=2)
    Σ = cov(permutedims(X))
    Ω = pinv(Σ)
    p = size(Ω, 1)

    edges = DataFrame(source=String[], target=String[], weight=Float64[])
    for i in 1:p
        for j in (i + 1):p
            denom = sqrt(max(Ω[i, i] * Ω[j, j], eps(Float64)))
            ρ = -Ω[i, j] / denom
            isfinite(ρ) || continue
            push!(edges, ("gene_$(i)", "gene_$(j)", abs(ρ)))
        end
    end
    sort!(edges, :weight, rev=true)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, edges[1:min(top_k, nrow(edges)), :], "causal_grn_inference")
end

"""
    gene_program_nmf(expression; n_programs=5, n_iter=200)

Nonnegative matrix factorization for gene program discovery.
"""
function gene_program_nmf(expression::AbstractMatrix{<:Real}; n_programs::Int=5, n_iter::Int=200, seed::Int=1, backend::Symbol=:auto)
    X = max.(Float64.(expression), 0.0)
    n_genes, n_cells = size(X)
    k = clamp(n_programs, 1, min(n_genes, n_cells))
    selected = resolve_backend(; backend=backend)
    Xdev = maybe_to_device(X; backend=selected)

    rng = MersenneTwister(seed)
    W = maybe_to_device(rand(rng, n_genes, k) .+ 1e-3; backend=selected)
    H = maybe_to_device(rand(rng, k, n_cells) .+ 1e-3; backend=selected)

    for _ in 1:n_iter
        H .*= (W' * Xdev) ./ max.(W' * W * H, 1e-12)
        W .*= (Xdev * H') ./ max.(W * (H * H'), 1e-12)
    end

    W_host = Matrix{Float64}(maybe_to_host(W))
    H_host = Matrix{Float64}(maybe_to_host(H))
    reconstruction = W_host * H_host
    err = norm(X - reconstruction) / max(norm(X), eps(Float64))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (program_loadings=W_host, program_activity=H_host, relative_error=err, backend=selected), "gene_program_nmf")
end

"""
    mofa_plus_em(assays; n_factors=8, n_iter=100)

EM-style low-rank integration with explicit missing-value reconstruction across
multiple omics assays (features x samples).
"""
function mofa_plus_em(assays::Vector{<:AbstractMatrix{<:Real}}; n_factors::Int=8, n_iter::Int=100, backend::Symbol=:auto, threaded::Bool=true)
    isempty(assays) && throw(ArgumentError("assays must contain at least one matrix"))
    n_samples = size(first(assays), 2)
    all(size(a, 2) == n_samples for a in assays) || throw(DimensionMismatch("all assays must share sample count in columns"))
    selected = resolve_backend(; backend=backend)

    blocks = [Matrix{Float64}(a) for a in assays]
    for b in blocks
        b[.!isfinite.(b)] .= NaN
    end
    X = vcat(blocks...)
    obs = .!isnan.(X)

    Ximp = copy(X)
    threaded_foreach(size(Ximp, 1), i -> begin
        row = @view Ximp[i, :]
        mask = @view obs[i, :]
        μ = any(mask) ? mean(row[mask]) : 0.0
        for j in axes(row, 1)
            mask[j] || (row[j] = μ)
        end
    end; threaded=threaded)

    used = clamp(n_factors, 1, min(size(Ximp, 1), size(Ximp, 2)))
    for _ in 1:max(1, n_iter)
        row_mean = mean(Ximp, dims=2)
        centered = Ximp .- row_mean
        fac = svd(maybe_to_device(permutedims(centered); backend=selected); full=false)
        U = Matrix{Float64}(maybe_to_host(fac.U))
        S = Vector{Float64}(maybe_to_host(fac.S))
        V = Matrix{Float64}(maybe_to_host(fac.V))
        r = min(used, length(S))
        factors = U[:, 1:r] * Diagonal(S[1:r])
        loadings = V[:, 1:r]
        recon = loadings * permutedims(factors) .+ row_mean
        Ximp[.!obs] .= recon[.!obs]
    end

    row_mean = mean(Ximp, dims=2)
    fac = svd(maybe_to_device(permutedims(Ximp .- row_mean); backend=selected); full=false)
    U = Matrix{Float64}(maybe_to_host(fac.U))
    S = Vector{Float64}(maybe_to_host(fac.S))
    V = Matrix{Float64}(maybe_to_host(fac.V))
    r = min(used, length(S))
    factors = U[:, 1:r] * Diagonal(S[1:r])
    loadings = V[:, 1:r]
    total_var = sum(abs2, S)
    explained = total_var > 0 ? (S[1:r] .^ 2) ./ total_var : fill(0.0, r)

    assay_loadings = Matrix{Float64}[]
    offset = 1
    for b in blocks
        n_rows = size(b, 1)
        push!(assay_loadings, loadings[offset:(offset + n_rows - 1), :])
        offset += n_rows
    end
    return MultiOmicsFactorAnalysisResult(
        factors,
        loadings,
        assay_loadings,
        explained,
        provenance_record(
            "MultiOmicsFactorAnalysisResult",
            "SystemsBio/mofa_plus_em";
            notes=["EM-style low-rank integration with per-assay block loadings"],
            parameters=(assay_count=length(blocks), sample_count=n_samples, factor_count=r, n_iter=Int(n_iter), backend=backend, threaded=Bool(threaded))))
end

"""
    sparse_cca_integration(x, y; n_components=10, sparsity=0.2)

CCA with hard-thresholded feature loadings to produce sparse projections.
"""
function sparse_cca_integration(x::AbstractMatrix{<:Real}, y::AbstractMatrix{<:Real}; n_components::Int=10, ridge::Real=1e-5, sparsity::Real=0.2, threaded::Bool=true)
    0 < sparsity <= 1 || throw(ArgumentError("sparsity must be in (0, 1]"))
    base = cca_integration(x, y; n_components=n_components, ridge=ridge)

    X = Matrix{Float64}(x)
    Y = Matrix{Float64}(y)
    X .-= mean(X, dims=1)
    Y .-= mean(Y, dims=1)
    n = max(size(X, 1) - 1, 1)

    x_load = (X' * base.x_scores) ./ n
    y_load = (Y' * base.y_scores) ./ n

    keep_x = falses(size(x_load))
    keep_y = falses(size(y_load))
    threaded_foreach(size(x_load, 2), c -> begin
        tx = quantile(abs.(x_load[:, c]), 1 - Float64(sparsity))
        ty = quantile(abs.(y_load[:, c]), 1 - Float64(sparsity))
        keep_x[:, c] .= abs.(x_load[:, c]) .>= tx
        keep_y[:, c] .= abs.(y_load[:, c]) .>= ty
    end; threaded=threaded)

    x_sparse_load = x_load .* keep_x
    y_sparse_load = y_load .* keep_y
    x_sparse_scores = X * x_sparse_load
    y_sparse_scores = Y * y_sparse_load
    return (
        x_scores=x_sparse_scores,
        y_scores=y_sparse_scores,
        x_loadings=x_sparse_load,
        y_loadings=y_sparse_load,
        selected_x=[findall(@view keep_x[:, c]) for c in axes(keep_x, 2)],
        selected_y=[findall(@view keep_y[:, c]) for c in axes(keep_y, 2)],
        canonical_correlations=base.canonical_correlations)
end

"""
    causal_ate_regression(df, treatment, outcome; covariates=Symbol[])

Estimate an adjusted average treatment effect with linear regression,
DoWhy-style backdoor adjustment.
"""
function causal_ate_regression(df::DataFrame, treatment::Symbol, outcome::Symbol; covariates::Vector{Symbol}=Symbol[])
    hasproperty(df, treatment) || throw(ArgumentError("missing treatment column"))
    hasproperty(df, outcome) || throw(ArgumentError("missing outcome column"))
    for c in covariates
        hasproperty(df, c) || throw(ArgumentError("missing covariate column $(c)"))
    end

    y = Float64.(df[!, outcome])
    t = Float64.(df[!, treatment])
    X = hcat(ones(Float64, nrow(df)), t)
    if !isempty(covariates)
        Z = hcat([Float64.(df[!, c]) for c in covariates]...)
        X = hcat(X, Z)
    end

    β = X \ y
    resid = y .- X * β
    dof = max(length(y) - size(X, 2), 1)
    σ2 = sum(abs2, resid) / dof
    vcov = σ2 * inv(X' * X)
    se = sqrt(max(vcov[2, 2], eps(Float64)))
    ate = β[2]
    z = ate / se
    p = erfc(abs(z) / sqrt(2.0))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (ate=ate, std_error=se, zscore=z, pvalue=p, coefficients=β, covariates=covariates), "causal_ate_regression")
end

end
