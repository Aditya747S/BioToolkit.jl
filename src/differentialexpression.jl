module DifferentialExpression

using SparseArrays
using Statistics
using LinearAlgebra
using Distributions
using Random
using DataFrames
using ProgressMeter

export CountMatrix, DEResult, GLMSolver, calc_norm_factors, estimate_dispersions, benjamini_hochberg, differential_expression, filter_low_counts, shrink_lfc, vst, fit_gene_fast!, remove_batch_effect, combat_correction, estimate_surrogates

"""
    CountMatrix(counts, gene_ids, sample_ids)

Sparse count-matrix container with gene and sample labels. The matrix stores
integer counts in gene-by-sample orientation.
"""
struct CountMatrix
    counts::SparseMatrixCSC{Int,Int}
    gene_ids::Vector{String}
    sample_ids::Vector{String}

    function CountMatrix(counts::SparseMatrixCSC{Int,Int}, gene_ids::Vector{String}, sample_ids::Vector{String})
        size(counts, 1) == length(gene_ids) || throw(ArgumentError("gene_ids must match the number of rows in counts"))
        size(counts, 2) == length(sample_ids) || throw(ArgumentError("sample_ids must match the number of columns in counts"))
        new(counts, gene_ids, sample_ids)
    end
end

function CountMatrix(counts::SparseMatrixCSC{Int,Int}, gene_ids::AbstractVector{<:AbstractString}, sample_ids::AbstractVector{<:AbstractString})
    return CountMatrix(counts, String.(gene_ids), String.(sample_ids))
end

function CountMatrix(counts::AbstractMatrix{<:Integer}, gene_ids::AbstractVector{<:AbstractString}, sample_ids::AbstractVector{<:AbstractString})
    return CountMatrix(sparse(Int.(counts)), gene_ids, sample_ids)
end

"""
    filter_low_counts(counts::CountMatrix; min_total=10)

Drop genes whose total count does not exceed `min_total`.
"""
function filter_low_counts(counts::CountMatrix; min_total::Integer=10)
    min_total >= 0 || throw(ArgumentError("min_total must be nonnegative"))
    row_sums = vec(sum(counts.counts, dims=2))
    keep = row_sums .> min_total
    return CountMatrix(counts.counts[keep, :], counts.gene_ids[keep], counts.sample_ids)
end

"""
    DEResult

Single-gene differential-expression result record.
"""
struct DEResult
    gene_id::String
    base_mean::Float64
    log2_fold_change::Float64
    lfc_se::Float64
    stat::Float64
    pvalue::Float64
    padj::Float64
end

"""
    DataFrames.DataFrame(results::AbstractVector{<:DEResult})

Convert a vector of differential-expression results into a tabular data frame.
"""
function DataFrames.DataFrame(results::AbstractVector{<:DEResult})
    rows = collect(results)
    return DataFrames.DataFrame(
        gene_id = [result.gene_id for result in rows],
        base_mean = [result.base_mean for result in rows],
        log2_fold_change = [result.log2_fold_change for result in rows],
        lfc_se = [result.lfc_se for result in rows],
        stat = [result.stat for result in rows],
        pvalue = [result.pvalue for result in rows],
        padj = [result.padj for result in rows],
    )
end

"""
    Base.show(io::IO, result::DEResult)

Print a compact one-line summary for a differential-expression result.
"""
function Base.show(io::IO, result::DEResult)
    print(io, "DEResult(", result.gene_id, ", log2FC=", result.log2_fold_change, ", pvalue=", result.pvalue, ", padj=", result.padj, ")")
end

"""
    GLMSolver(X)

Internal working state for iteratively fitting the gene-wise GLM.
"""
struct GLMSolver
    X::Matrix{Float64}
    xtwx::Matrix{Float64}
    xtwz::Vector{Float64}
    beta::Vector{Float64}
    eta::Vector{Float64}
    mu::Vector{Float64}
    weights::Vector{Float64}
    z::Vector{Float64}

    function GLMSolver(X::AbstractMatrix{<:Real})
        design = Matrix{Float64}(X)
        nobs, ncoef = size(design)
        new(design, zeros(Float64, ncoef, ncoef), zeros(Float64, ncoef), zeros(Float64, ncoef), zeros(Float64, nobs), zeros(Float64, nobs), zeros(Float64, nobs), zeros(Float64, nobs))
    end
end

function _as_dense_row(counts::AbstractMatrix{<:Integer}, gene_index::Integer)
    return Float64[counts[gene_index, sample] for sample in axes(counts, 2)]
end

function _as_dense_row(counts::AbstractMatrix{<:Real}, gene_index::Integer)
    return Float64[counts[gene_index, sample] for sample in axes(counts, 2)]
end

function _library_sizes(counts::AbstractMatrix{<:Integer})
    return vec(sum(counts, dims=1))
end

function _row_totals(counts::CountMatrix)
    return vec(sum(counts.counts, dims=2))
end

function _trimmed_weighted_mean(m_values::Vector{Float64}, a_values::Vector{Float64}, weights::Vector{Float64}; logratio_trim::Float64=0.3, sum_trim::Float64=0.05)
    isempty(m_values) && return 0.0
    low_m, high_m = quantile(m_values, [logratio_trim / 2, 1 - logratio_trim / 2])
    low_a, high_a = quantile(a_values, [sum_trim / 2, 1 - sum_trim / 2])
    keep = (m_values .>= low_m) .& (m_values .<= high_m) .& (a_values .>= low_a) .& (a_values .<= high_a)
    any(keep) || return 0.0
    kept_weights = weights[keep]
    kept_values = m_values[keep]
    total_weight = sum(kept_weights)
    total_weight > 0 || return 0.0
    return sum(kept_weights .* kept_values) / total_weight
end

"""
    calc_norm_factors(counts::AbstractMatrix{<:Integer}; threaded=true)

Estimate sample-wise normalization factors using a trimmed mean of M-values.
"""
function calc_norm_factors(counts::AbstractMatrix{<:Integer}; threaded::Bool=true)
    ngenes, nsamples = size(counts)
    nsamples > 0 || throw(ArgumentError("counts must contain at least one sample"))

    lib_sizes = _library_sizes(counts)
    reference = argmin(abs.(lib_sizes .- median(lib_sizes)))
    ref_counts = Float64[counts[gene, reference] for gene in 1:ngenes]
    ref_lib = Float64(lib_sizes[reference])

    factors = ones(Float64, nsamples)

    if threaded && nsamples > 1 && Threads.nthreads() > 1
        Threads.@threads for sample in 1:nsamples
            sample == reference && continue

            sample_counts = Float64[counts[gene, sample] for gene in 1:ngenes]
            sample_lib = Float64(lib_sizes[sample])
            valid = (sample_counts .> 0) .& (ref_counts .> 0)
            any(valid) || continue

            y = sample_counts[valid]
            r = ref_counts[valid]
            normalized_sample = y ./ sample_lib
            normalized_ref = r ./ ref_lib
            m_values = log2.(normalized_sample) .- log2.(normalized_ref)
            a_values = 0.5 .* (log2.(normalized_sample) .+ log2.(normalized_ref))
            weights = 1.0 ./ (1.0 ./ y .+ 1.0 ./ r)
            factors[sample] = 2.0 ^ _trimmed_weighted_mean(m_values, a_values, weights)
        end
    else
        for sample in 1:nsamples
            sample == reference && continue

            sample_counts = Float64[counts[gene, sample] for gene in 1:ngenes]
            sample_lib = Float64(lib_sizes[sample])
            valid = (sample_counts .> 0) .& (ref_counts .> 0)
            any(valid) || continue

            y = sample_counts[valid]
            r = ref_counts[valid]
            normalized_sample = y ./ sample_lib
            normalized_ref = r ./ ref_lib
            m_values = log2.(normalized_sample) .- log2.(normalized_ref)
            a_values = 0.5 .* (log2.(normalized_sample) .+ log2.(normalized_ref))
            weights = 1.0 ./ (1.0 ./ y .+ 1.0 ./ r)
            factors[sample] = 2.0 ^ _trimmed_weighted_mean(m_values, a_values, weights)
        end
    end

    geometric_mean = exp(mean(log.(factors)))
    factors ./= geometric_mean
    return factors
end

"""
    calc_norm_factors(counts::CountMatrix)

Estimate normalization factors for a sparse count matrix.
"""
calc_norm_factors(counts::CountMatrix) = calc_norm_factors(counts.counts)

"""
    estimate_dispersions(counts::AbstractMatrix{<:Integer}, norm_factors; threaded=true)

Estimate gene-wise dispersions from normalized count profiles.
"""
function estimate_dispersions(counts::AbstractMatrix{<:Integer}, norm_factors::AbstractVector{<:Real}; threaded::Bool=true)
    ngenes, nsamples = size(counts)
    nsamples == length(norm_factors) || throw(ArgumentError("norm_factors must match the number of samples"))

    means = zeros(Float64, ngenes)
    raw = zeros(Float64, ngenes)

    if threaded && ngenes > 1 && Threads.nthreads() > 1
        Threads.@threads for gene in 1:ngenes
            normalized = [Float64(counts[gene, sample]) / Float64(norm_factors[sample]) for sample in 1:nsamples]
            means[gene] = mean(normalized)
            if nsamples > 1 && means[gene] > 0
                variance = var(normalized; corrected=true)
                raw[gene] = max((variance - means[gene]) / (means[gene]^2), 0.0)
            end
        end
    else
        for gene in 1:ngenes
            normalized = [Float64(counts[gene, sample]) / Float64(norm_factors[sample]) for sample in 1:nsamples]
            means[gene] = mean(normalized)
            if nsamples > 1 && means[gene] > 0
                variance = var(normalized; corrected=true)
                raw[gene] = max((variance - means[gene]) / (means[gene]^2), 0.0)
            end
        end
    end

    positive = (means .> 0) .& (raw .> 0)
    trend = similar(raw)
    if count(positive) >= 2
        x = log.(means[positive] .+ eps())
        y = log.(raw[positive] .+ eps())
        design = hcat(ones(length(x)), x)
        coefficients = design \ y
        trend .= exp.(coefficients[1] .+ coefficients[2] .* log.(means .+ eps()))
    else
        fallback = count(raw .> 0) > 0 ? median(raw[raw .> 0]) : 0.1
        trend .= fallback
    end

    shrink = clamp.(means ./ (means .+ 10.0), 0.05, 0.95)
    dispersions = shrink .* raw .+ (1 .- shrink) .* trend
    return max.(dispersions, eps())
end

"""
    estimate_dispersions(counts::CountMatrix, norm_factors)

Estimate dispersions for a sparse count matrix.
"""
estimate_dispersions(counts::CountMatrix, norm_factors::AbstractVector{<:Real}) = estimate_dispersions(counts.counts, norm_factors)

function _solve_irls!(solver::GLMSolver)
    if size(solver.X, 2) == 2
        a = solver.xtwx[1, 1]
        b = solver.xtwx[1, 2]
        c = solver.xtwx[2, 2]
        det = a * c - b * b
        abs(det) <= eps() && return solver.beta
        x0 = solver.xtwz[1]
        x1 = solver.xtwz[2]
        solver.beta[1] = (c * x0 - b * x1) / det
        solver.beta[2] = (-b * x0 + a * x1) / det
        return solver.beta
    end

    solver.beta .= solver.xtwx \ solver.xtwz
    return solver.beta
end

"""
    fit_gene_fast!(solver::GLMSolver, counts, offset, phi)

Fit a negative-binomial GLM for one gene using iterative reweighted least squares.
"""
function fit_gene_fast!(solver::GLMSolver, counts::AbstractVector{<:Real}, offset::AbstractVector{<:Real}, phi::Real)
    nobs = length(counts)
    length(offset) == nobs || throw(ArgumentError("offset must match counts"))

    y = Float64.(counts)
    offset_values = Float64.(offset)
    solver.beta .= 0.0
    solver.beta[1] = log(mean(y .+ 0.1))

    for _ in 1:30
        mul!(solver.eta, solver.X, solver.beta)
        @inbounds for i in 1:nobs
            solver.eta[i] += offset_values[i]
        end
        @. solver.mu = exp(clamp(solver.eta, -30.0, 30.0))
        @. solver.weights = max(solver.mu / (1 + Float64(phi) * solver.mu), 1e-8)
        @. solver.z = solver.eta + (y - solver.mu) / max(solver.mu, 1e-8) - offset_values

        fill!(solver.xtwx, 0.0)
        fill!(solver.xtwz, 0.0)
        ncoef = size(solver.X, 2)
        @inbounds for row in 1:nobs
            w = solver.weights[row]
            z = solver.z[row]
            for col in 1:ncoef
                xcol = solver.X[row, col]
                solver.xtwz[col] += xcol * w * z
                for inner in col:ncoef
                    solver.xtwx[col, inner] += xcol * solver.X[row, inner] * w
                end
            end
        end
        for col in 1:ncoef, inner in 1:col-1
            solver.xtwx[col, inner] = solver.xtwx[inner, col]
        end

        previous_beta = copy(solver.beta)
        beta_new = _solve_irls!(solver)
        if norm(beta_new - previous_beta) <= 1e-8 * (1 + norm(previous_beta))
            break
        end
    end

    mul!(solver.eta, solver.X, solver.beta)
    @inbounds for i in 1:nobs
        solver.eta[i] += offset_values[i]
    end
    @. solver.mu = exp(clamp(solver.eta, -30.0, 30.0))
    @. solver.weights = max(solver.mu / (1 + Float64(phi) * solver.mu), 1e-8)

    fill!(solver.xtwx, 0.0)
    ncoef = size(solver.X, 2)
    @inbounds for row in 1:nobs
        w = solver.weights[row]
        for col in 1:ncoef
            xcol = solver.X[row, col]
            for inner in col:ncoef
                solver.xtwx[col, inner] += xcol * solver.X[row, inner] * w
            end
        end
    end
    for col in 1:ncoef, inner in 1:col-1
        solver.xtwx[col, inner] = solver.xtwx[inner, col]
    end

    covariance = inv(Symmetric(solver.xtwx))
    lfc_se = size(solver.X, 2) > 1 ? sqrt(max(covariance[2, 2], 0.0)) : 0.0
    stat = lfc_se > 0 ? (solver.beta[2] / lfc_se)^2 : 0.0
    return solver.beta, lfc_se, stat
end

"""
    shrink_lfc(results; prior_sd=nothing)

Apply empirical shrinkage to log2 fold-changes in a vector of results.
"""
function shrink_lfc(results::AbstractVector{<:DEResult}; prior_sd::Union{Nothing,Real}=nothing)
    isempty(results) && return DEResult[]
    lfc_values = [result.log2_fold_change for result in results]
    estimated_prior_sd = prior_sd === nothing ? max(median(abs.(lfc_values)) / 0.6745, eps()) : Float64(prior_sd)

    shrunk = DEResult[]
    for result in results
        weight = estimated_prior_sd^2 / (estimated_prior_sd^2 + max(result.lfc_se, eps())^2)
        shrunk_lfc = result.log2_fold_change * weight
        shrunk_se = sqrt(max((result.lfc_se^2 * estimated_prior_sd^2) / (result.lfc_se^2 + estimated_prior_sd^2), eps()))
        push!(shrunk, DEResult(result.gene_id, result.base_mean, shrunk_lfc, shrunk_se, result.stat, result.pvalue, result.padj))
    end
    return shrunk
end

"""
    vst(counts::CountMatrix, dispersions; norm_factors=calc_norm_factors(counts))

Apply a variance-stabilizing transform to sparse count data.
"""
function vst(counts::CountMatrix, dispersions::AbstractVector{<:Real}; norm_factors::AbstractVector{<:Real}=calc_norm_factors(counts))
    ngenes, nsamples = size(counts.counts)
    length(dispersions) == ngenes || throw(ArgumentError("dispersions must match the number of genes"))
    length(norm_factors) == nsamples || throw(ArgumentError("norm_factors must match the number of samples"))

    transformed = Matrix{Float64}(undef, ngenes, nsamples)
    for gene in 1:ngenes
        disp = max(Float64(dispersions[gene]), eps())
        for sample in 1:nsamples
            normalized = Float64(counts.counts[gene, sample]) / Float64(norm_factors[sample])
            transformed[gene, sample] = asinh(sqrt(normalized / (1 + disp)))
        end
    end
    return transformed
end

function _fit_gene_model(y::AbstractVector{<:Real}, X::AbstractMatrix{<:Real}, offset::AbstractVector{<:Real}, phi::Real)
    y_values = Float64.(y)
    offset_values = Float64.(offset)
    beta = zeros(Float64, size(X, 2))
    beta[1] = log(mean(y_values .+ 0.1))

    for _ in 1:50
        eta = offset_values .+ X * beta
        mu = exp.(clamp.(eta, -30.0, 30.0))
        weights = max.(mu ./ (1 .+ Float64(phi) .* mu), 1e-8)
        target = eta .+ (y_values .- mu) ./ max.(mu, 1e-8) .- offset_values
        weighted_X = X .* reshape(sqrt.(weights), :, 1)
        weighted_target = target .* sqrt.(weights)
        beta_new = weighted_X \ weighted_target
        if norm(beta_new - beta) <= 1e-8 * (1 + norm(beta))
            beta = beta_new
            break
        end
        beta = beta_new
    end

    eta = offset_values .+ X * beta
    mu = exp.(clamp.(eta, -30.0, 30.0))
    weights = max.(mu ./ (1 .+ Float64(phi) .* mu), 1e-8)
    xtwx = X' * (X .* reshape(weights, :, 1))
    covariance = inv(Symmetric(xtwx))
    lfc_se = sqrt(max(covariance[2, 2], 0.0))
    stat = lfc_se > 0 ? (beta[2] / lfc_se)^2 : 0.0
    return beta, lfc_se, stat
end

"""
    benjamini_hochberg(pvalues)

Adjust p-values using the Benjamini-Hochberg false-discovery-rate procedure.
"""
function benjamini_hochberg(pvalues::AbstractVector{<:Real})
    n = length(pvalues)
    n == 0 && return Float64[]

    sanitized = [clamp(isfinite(Float64(pvalue)) ? Float64(pvalue) : 1.0, 0.0, 1.0) for pvalue in pvalues]
    order = sortperm(sanitized)
    sorted = sanitized[order]
    adjusted_sorted = similar(sorted)

    running_min = 1.0
    for index in n:-1:1
        running_min = min(running_min, sorted[index] * n / index)
        adjusted_sorted[index] = min(running_min, 1.0)
    end

    adjusted = similar(sanitized)
    for (rank, original_index) in enumerate(order)
        adjusted[original_index] = adjusted_sorted[rank]
    end

    return adjusted
end

function _fit_differential_expression(count_matrix::CountMatrix, design::AbstractVector{<:Symbol}, norm_factors::AbstractVector{<:Real}, dispersions::AbstractVector{<:Real}, treatment_level::Symbol; show_progress::Bool=false)
    ngenes, nsamples = size(count_matrix.counts)
    condition = Float64[design[index] == treatment_level ? 1.0 : 0.0 for index in 1:nsamples]
    X = hcat(ones(Float64, nsamples), condition)
    chunk_size = max(1, cld(ngenes, max(Base.Threads.nthreads(), 1)))

    results = Vector{DEResult}(undef, ngenes)
    if show_progress
        solver = GLMSolver(X)
        @showprogress for gene in 1:ngenes
            y = _as_dense_row(count_matrix.counts, gene)
            base_mean = mean(y ./ norm_factors)

            if all(iszero, y)
                results[gene] = DEResult(count_matrix.gene_ids[gene], base_mean, 0.0, Inf, 0.0, 1.0, 1.0)
                continue
            end

            beta, lfc_se_natural, stat = fit_gene_fast!(solver, y, log.(norm_factors), dispersions[gene])
            log2_fold_change = beta[2] / log(2)
            lfc_se = lfc_se_natural / log(2)
            pvalue = ccdf(Chisq(1), stat)
            results[gene] = DEResult(count_matrix.gene_ids[gene], base_mean, log2_fold_change, lfc_se, stat, pvalue, 1.0)
        end
    else
        Base.Threads.@threads for chunk_start in 1:chunk_size:ngenes
            solver = GLMSolver(X)
            chunk_stop = min(chunk_start + chunk_size - 1, ngenes)
            for gene in chunk_start:chunk_stop
                y = _as_dense_row(count_matrix.counts, gene)
                base_mean = mean(y ./ norm_factors)

                if all(iszero, y)
                    results[gene] = DEResult(count_matrix.gene_ids[gene], base_mean, 0.0, Inf, 0.0, 1.0, 1.0)
                    continue
                end

                beta, lfc_se_natural, stat = fit_gene_fast!(solver, y, log.(norm_factors), dispersions[gene])
                log2_fold_change = beta[2] / log(2)
                lfc_se = lfc_se_natural / log(2)
                pvalue = ccdf(Chisq(1), stat)
                results[gene] = DEResult(count_matrix.gene_ids[gene], base_mean, log2_fold_change, lfc_se, stat, pvalue, 1.0)
            end
        end
    end

    return results
end

"""
    differential_expression(count_matrix::CountMatrix, design; min_total=10, shrink=false, show_progress=false)

Run the full two-group differential-expression workflow for a count matrix.
"""
function differential_expression(count_matrix::CountMatrix, design::AbstractVector{<:Symbol}; min_total::Integer=10, shrink::Bool=false, show_progress::Bool=false)
    ngenes, nsamples = size(count_matrix.counts)
    length(design) == nsamples || throw(ArgumentError("design must match the number of samples"))

    levels = unique(design)
    length(levels) == 2 || throw(ArgumentError("design must contain exactly two conditions"))

    row_totals = _row_totals(count_matrix)
    keep = row_totals .> min_total
    norm_factors = calc_norm_factors(count_matrix)
    dispersions = estimate_dispersions(count_matrix, norm_factors)

    results = Vector{DEResult}(undef, ngenes)
    if any(keep)
        filtered = CountMatrix(count_matrix.counts[keep, :], count_matrix.gene_ids[keep], count_matrix.sample_ids)
        filtered_results = _fit_differential_expression(filtered, design, norm_factors, dispersions[keep], levels[2]; show_progress=show_progress)
        filtered_padj = benjamini_hochberg([result.pvalue for result in filtered_results])
        filtered_results = [DEResult(result.gene_id, result.base_mean, result.log2_fold_change, result.lfc_se, result.stat, result.pvalue, filtered_padj[index]) for (index, result) in enumerate(filtered_results)]
        iterator = 1
        for gene in 1:ngenes
            if keep[gene]
                results[gene] = filtered_results[iterator]
                iterator += 1
            else
                results[gene] = DEResult(count_matrix.gene_ids[gene], 0.0, 0.0, Inf, 0.0, 1.0, 1.0)
            end
        end
    else
        for gene in 1:ngenes
            results[gene] = DEResult(count_matrix.gene_ids[gene], 0.0, 0.0, Inf, 0.0, 1.0, 1.0)
        end
    end

    if shrink
        results = shrink_lfc(results)
    end

    return results
end

function _design_basis(design::AbstractMatrix{<:Real})
    matrix = Matrix{Float64}(design)
    size(matrix, 1) > 0 || return zeros(Float64, size(matrix, 1), 0)
    basis = qr(matrix, ColumnNorm())
    diagonal = abs.(diag(Matrix(basis.R)))
    scale = isempty(diagonal) ? 1.0 : max(maximum(abs, matrix), 1.0)
    tolerance = max(size(matrix)...)*eps(scale)
    rank_value = count(value -> value > tolerance, diagonal)
    rank_value == 0 && return zeros(Float64, size(matrix, 1), 0)
    return Matrix(basis.Q)[:, 1:rank_value]
end

function _batch_design_matrix(batch_labels::AbstractVector)
    samples = length(batch_labels)
    samples > 0 || throw(ArgumentError("batch_labels must contain at least one sample"))
    levels = unique(batch_labels)
    design = zeros(Float64, samples, length(levels))
    lookup = Dict(level => index for (index, level) in enumerate(levels))
    for (sample, label) in enumerate(batch_labels)
        design[sample, lookup[label]] = 1.0
    end
    if !isempty(levels)
        for column in axes(design, 2)
            design[:, column] .-= mean(design[:, column])
        end
    end
    return design
end

function _project_onto_design(data::AbstractMatrix{<:Real}, design::AbstractMatrix{<:Real})
    matrix = Matrix{Float64}(design)
    isempty(matrix) && return zeros(Float64, size(data))
    projector = matrix * pinv(matrix)
    return Matrix{Float64}(data) * projector
end

function _residualize_design(design::AbstractMatrix{<:Real}, bio_design::Union{Nothing,AbstractMatrix{<:Real}})
    matrix = Matrix{Float64}(design)
    bio_design === nothing && return matrix
    bio = hcat(ones(Float64, size(matrix, 1)), Matrix{Float64}(bio_design))
    basis = _design_basis(bio)
    isempty(basis) && return matrix
    return matrix - basis * (basis' * matrix)
end

"""
    remove_batch_effect(data, batch_labels; bio_design=nothing)

Remove batch effects from an expression matrix using batch labels.
"""
function remove_batch_effect(data::AbstractMatrix{<:Real}, batch_labels::AbstractVector; bio_design::Union{Nothing,AbstractMatrix{<:Real}}=nothing)
    design = _batch_design_matrix(batch_labels)
    return remove_batch_effect(data, design; bio_design=bio_design)
end

"""
    remove_batch_effect(data, batch_design; bio_design=nothing)

Remove batch effects from an expression matrix using an explicit batch design.
"""
function remove_batch_effect(data::AbstractMatrix{<:Real}, batch_design::AbstractMatrix{<:Real}; bio_design::Union{Nothing,AbstractMatrix{<:Real}}=nothing)
    matrix = Matrix{Float64}(data)
    design = _residualize_design(batch_design, bio_design)
    isempty(design) && return matrix
    projector = design * pinv(design)
    return matrix - matrix * projector
end

"""
    remove_batch_effect(counts::CountMatrix, batch_labels; bio_design=nothing)

Remove batch effects after converting a sparse count matrix to dense form.
"""
remove_batch_effect(counts::CountMatrix, batch_labels::AbstractVector; bio_design::Union{Nothing,AbstractMatrix{<:Real}}=nothing) = remove_batch_effect(Matrix{Float64}(counts.counts), batch_labels; bio_design=bio_design)

"""
    remove_batch_effect(counts::CountMatrix, batch_design; bio_design=nothing)

Remove batch effects from sparse count data using an explicit batch design.
"""
remove_batch_effect(counts::CountMatrix, batch_design::AbstractMatrix{<:Real}; bio_design::Union{Nothing,AbstractMatrix{<:Real}}=nothing) = remove_batch_effect(Matrix{Float64}(counts.counts), batch_design; bio_design=bio_design)

"""
    combat_correction(data, batch_labels; bio_design=nothing)

Apply ComBat-style batch correction using batch labels.
"""
function combat_correction(data::AbstractMatrix{<:Real}, batch_labels::AbstractVector; bio_design::Union{Nothing,AbstractMatrix{<:Real}}=nothing)
    design = _batch_design_matrix(batch_labels)
    return combat_correction(data, design; bio_design=bio_design)
end

"""
    combat_correction(data, batch_design; bio_design=nothing, threaded=true)

Apply ComBat-style batch correction using an explicit batch design.
"""
function combat_correction(data::AbstractMatrix{<:Real}, batch_design::AbstractMatrix{<:Real}; bio_design::Union{Nothing,AbstractMatrix{<:Real}}=nothing, threaded::Bool=true)
    matrix = Matrix{Float64}(data)
    size(matrix, 2) == size(batch_design, 1) || throw(ArgumentError("batch design must have one row per sample"))

    bio_fit = bio_design === nothing ? zeros(Float64, size(matrix)) : _project_onto_design(matrix, hcat(ones(Float64, size(matrix, 2)), Matrix{Float64}(bio_design)))
    residual = matrix - bio_fit

    gene_means = vec(mean(residual, dims=2))
    gene_sds = vec(std(residual, dims=2; corrected=true))
    gene_sds = max.(gene_sds, eps(Float64))
    standardized = (residual .- gene_means) ./ gene_sds

    design = _residualize_design(batch_design, bio_design)
    basis = _design_basis(design)
    isempty(basis) && return matrix

    batch_assignments = [argmax(batch_design[sample, :]) for sample in axes(batch_design, 1)]
    n_batches = size(design, 2)
    corrected = copy(standardized)

    if threaded && n_batches > 1 && Threads.nthreads() > 1
        Threads.@threads for batch_index in 1:n_batches
            sample_indices = findall(==(batch_index), batch_assignments)
            isempty(sample_indices) && continue
            batch_slice = standardized[:, sample_indices]
            gamma_hat = vec(mean(batch_slice, dims=2))
            sigma2_hat = vec(var(batch_slice, dims=2; corrected=true))
            sigma2_hat = max.(sigma2_hat, eps(Float64))
            gamma_bar = mean(gamma_hat)
            tau2 = max(var(gamma_hat; corrected=true), eps(Float64))
            sample_count = length(sample_indices)
            post_mean = (tau2 .* gamma_hat .+ (sigma2_hat ./ sample_count) .* gamma_bar) ./ (tau2 .+ sigma2_hat ./ sample_count)
            corrected[:, sample_indices] .-= post_mean
        end
    else
        for batch_index in 1:n_batches
            sample_indices = findall(==(batch_index), batch_assignments)
            isempty(sample_indices) && continue
            batch_slice = standardized[:, sample_indices]
            gamma_hat = vec(mean(batch_slice, dims=2))
            sigma2_hat = vec(var(batch_slice, dims=2; corrected=true))
            sigma2_hat = max.(sigma2_hat, eps(Float64))
            gamma_bar = mean(gamma_hat)
            tau2 = max(var(gamma_hat; corrected=true), eps(Float64))
            sample_count = length(sample_indices)
            post_mean = (tau2 .* gamma_hat .+ (sigma2_hat ./ sample_count) .* gamma_bar) ./ (tau2 .+ sigma2_hat ./ sample_count)
            corrected[:, sample_indices] .-= post_mean
        end
    end

    return corrected .* gene_sds .+ gene_means .+ bio_fit
end

"""
    combat_correction(counts::CountMatrix, batch_labels; bio_design=nothing)

Apply ComBat-style batch correction to sparse count data using batch labels.
"""
combat_correction(counts::CountMatrix, batch_labels::AbstractVector; bio_design::Union{Nothing,AbstractMatrix{<:Real}}=nothing) = combat_correction(Matrix{Float64}(counts.counts), batch_labels; bio_design=bio_design)

"""
    combat_correction(counts::CountMatrix, batch_design; bio_design=nothing)

Apply ComBat-style batch correction to sparse count data using a batch design.
"""
combat_correction(counts::CountMatrix, batch_design::AbstractMatrix{<:Real}; bio_design::Union{Nothing,AbstractMatrix{<:Real}}=nothing) = combat_correction(Matrix{Float64}(counts.counts), batch_design; bio_design=bio_design)

"""
    estimate_surrogates(data; bio_design=nothing, nperm=20, alpha=0.1, max_components=10, rng=Random.default_rng())

Estimate surrogate variables from the residual structure in a matrix of data.
"""
function estimate_surrogates(data::AbstractMatrix{<:Real}; bio_design::Union{Nothing,AbstractMatrix{<:Real}}=nothing, nperm::Integer=20, alpha::Real=0.1, max_components::Integer=10, rng::Random.AbstractRNG=Random.default_rng())
    matrix = Matrix{Float64}(data)
    nperm >= 0 || throw(ArgumentError("nperm must be nonnegative"))
    0 < alpha < 1 || throw(ArgumentError("alpha must be between 0 and 1"))

    null_fit = bio_design === nothing ? repeat(mean(matrix, dims=2), 1, size(matrix, 2)) : _project_onto_design(matrix, hcat(ones(Float64, size(matrix, 2)), Matrix{Float64}(bio_design)))
    residual = matrix - null_fit
    centered = residual .- mean(residual, dims=2)

    full_svd = svd(centered; full=false)
    singular_values = full_svd.S
    component_limit = min(max_components, length(singular_values), max(0, size(matrix, 2) - 1))
    component_limit == 0 && return zeros(Float64, size(matrix, 2), 0)

    permuted_singular_values = zeros(Float64, component_limit, nperm)
    if nperm > 0
        for perm_index in 1:nperm
            permuted = similar(centered)
            for gene in axes(centered, 1)
                permuted[gene, :] = centered[gene, randperm(rng, size(centered, 2))]
            end
            permuted_singular_values[:, perm_index] .= svd(permuted; full=false).S[1:component_limit]
        end
    end

    significant = Int[]
    for component in 1:component_limit
        if nperm == 0
            push!(significant, component)
            continue
        end
        empirical_p = (1 + count(permuted_singular_values[component, :] .>= singular_values[component])) / (nperm + 1)
        empirical_p <= alpha && push!(significant, component)
    end

    isempty(significant) && return zeros(Float64, size(matrix, 2), 0)
    return full_svd.V[:, significant]
end

"""
    estimate_surrogates(counts::CountMatrix; bio_design=nothing, nperm=20, alpha=0.1, max_components=10, rng=Random.default_rng())

Estimate surrogate variables from sparse count data.
"""
estimate_surrogates(counts::CountMatrix; bio_design::Union{Nothing,AbstractMatrix{<:Real}}=nothing, nperm::Integer=20, alpha::Real=0.1, max_components::Integer=10, rng::Random.AbstractRNG=Random.default_rng()) = estimate_surrogates(Matrix{Float64}(counts.counts); bio_design=bio_design, nperm=nperm, alpha=alpha, max_components=max_components, rng=rng)

end # module DifferentialExpression