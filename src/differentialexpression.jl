# =============================================================================
# DifferentialExpression.jl 
# =============================================================================

module DifferentialExpression

using SparseArrays
using Statistics
using LinearAlgebra
using Distributions
using Random
using DataFrames
using SpecialFunctions
using Printf
using Optim
using Dates
using SHA

using ..BioToolkit: AbstractAnalysisResult, ProvenanceContext, ProvenanceParams, ResultProvenance, ThreadSafeProvenanceContext, active_provenance_context, analysis_result_summary, container_provenance_summary, ensure_provenance_id!, metadata_provenance, new_provenance_id, provenance_parent_ids, provenance_record, provenance_result!, register_container_provenance!, register_provenance!, stamp_provenance!, update_provenance!, with_provenance

# -----------------------------------------------------------------------------
# Exports
# -----------------------------------------------------------------------------
export CountMatrix, DEResult, GLMSolver, DESeqDataSet, DESeqDataSetFromMatrix,
    # Normalization
    calc_tmm_factors, calc_norm_factors, estimateSizeFactorsForMatrix,
    estimateSizeFactors, sizeFactors, sizeFactors!,
    normalizationFactors, normalizationFactors!,
    # Dispersion
    estimate_dispersions, estimate_dispersions_prior,
    estimateDispersions, estimateDispersionsGeneEst, estimateDispersionsFit,
    estimateDispersionsMAP, estimateDispersionsPriorVar,
    # Testing
    cooks_distance, replace_outliers, nbinomWaldTest, nbinomLRT,
    # Results
    results, resultsNames, benjamini_hochberg,
    # High-level
    DESeq, differential_expression, filter_low_counts, shrink_lfc, lfcShrink,
    # Transformations
    vst, voom,
    # Batch correction
    remove_batch_effect, combat_correction,
    # Surrogate variables
    estimate_surrogates, sv_analysis,
    # RUV
    ruv_normalize,
    # LFC shrinkage
    lfc_shrink_apeglm, lfc_shrink_apeglm_full, lfc_shrink_ashr, adaptive_shrinkage,
    # edgeR API
    estimate_dispersion_edgeR, trend_dispersion_edgeR, dispersion_prior_dof,
    exact_test_edgeR, glm_ql_fit, glm_ql_f_test, edgeR_qlf_test,
    # Prior variance
    estimate_beta_prior_var,
    # NEW: Advanced methods
    zinb_test, time_series_de, dge_mixedmodel,
    power_analysis, bootstrap_lfc, confounded_contrast_test,
    calibrate_pvalues_empirical_null, concordance_analysis,
    validate_design,
    # Multiple testing
    storey_qvalue, ihw_qvalue, storey_pi0,
    # Infrastructure
    makeExampleDESeqDataSet, ReproducibleDDS, input_checksum

# =============================================================================
# Core Data Types
# =============================================================================

"""
    CountMatrix(counts, gene_ids, sample_ids)

Sparse count matrix with gene and sample identifiers.
"""
struct CountMatrix
    counts::SparseMatrixCSC{Int,Int}
    gene_ids::Vector{String}
    sample_ids::Vector{String}
    metadata::Dict{Symbol,Any}

    function CountMatrix(counts::SparseMatrixCSC{Int,Int},
        gene_ids::Vector{String}, sample_ids::Vector{String}; metadata::AbstractDict=Dict{Symbol,Any}())
        size(counts, 1) == length(gene_ids) ||
            throw(ArgumentError("gene_ids length ($(length(gene_ids))) does not match row count ($(size(counts, 1)))"))
        size(counts, 2) == length(sample_ids) ||
            throw(ArgumentError("sample_ids length ($(length(sample_ids))) does not match column count ($(size(counts, 2)))"))
        metadata_copy = Dict{Symbol,Any}(metadata)
        ensure_provenance_id!(metadata_copy)
        metadata_provenance(metadata_copy) === nothing && stamp_provenance!(
            metadata_copy;
            label="CountMatrix",
            source="DifferentialExpression/CountMatrix",
            notes=["constructed count matrix"],
            parameters=(n_genes=length(gene_ids), n_samples=length(sample_ids)))
        new(counts, gene_ids, sample_ids, metadata_copy)
    end
end

CountMatrix(c::SparseMatrixCSC{<:Integer,<:Integer},
    g::AbstractVector{<:AbstractString}, s::AbstractVector{<:AbstractString}) =
    CountMatrix(sparse(Int.(c)), String.(g), String.(s))

CountMatrix(c::AbstractMatrix{<:Integer},
    g::AbstractVector{<:AbstractString}, s::AbstractVector{<:AbstractString}) =
    CountMatrix(sparse(Int.(c)), String.(g), String.(s))

function Base.show(io::IO, cm::CountMatrix)
    print(io, "CountMatrix(", size(cm.counts, 1), "x", size(cm.counts, 2), ", ", container_provenance_summary(cm), ")")
end

"""
    DEResult(gene_id, base_mean, log2_fold_change, lfc_se, stat, pvalue, padj, ...)

Result of differential expression testing for a single gene.
"""
struct DEResult <: AbstractAnalysisResult
    gene_id::String
    base_mean::Float64
    log2_fold_change::Float64
    lfc_se::Float64
    stat::Float64
    pvalue::Float64
    padj::Float64
    zero_inflated::Bool
    converged::Bool
    provenance::ResultProvenance

    function DEResult(gene_id, base_mean, lfc, lfc_se, stat, pvalue, padj,
        zero_inflated::Bool=false, converged::Bool=true,
        provenance::ResultProvenance=provenance_record(
            "DEResult",
            "DifferentialExpression/DEResult";
            notes=["constructed differential expression result"],
            parameters=(gene_id=String(gene_id),)))
        new(String(gene_id), Float64(base_mean), Float64(lfc), Float64(lfc_se),
            Float64(stat), Float64(pvalue), Float64(padj), zero_inflated, converged, provenance)
    end
end

function DataFrames.DataFrame(results::AbstractVector{<:DEResult})
    rows = collect(results)
    finite_or_missing(vals; clamp01::Bool=false) = begin
        out = Vector{Union{Missing,Float64}}(undef, length(vals))
        for i in eachindex(vals)
            fv = Float64(vals[i])
            if !isfinite(fv)
                out[i] = missing
            else
                out[i] = clamp01 ? clamp(fv, 0.0, 1.0) : fv
            end
        end
        out
    end
    return DataFrames.DataFrame(
        gene_id=[r.gene_id for r in rows],
        base_mean=[r.base_mean for r in rows],
        log2_fold_change=finite_or_missing([r.log2_fold_change for r in rows]),
        lfc_se=finite_or_missing([r.lfc_se for r in rows]),
        stat=finite_or_missing([r.stat for r in rows]),
        pvalue=finite_or_missing([r.pvalue for r in rows]; clamp01=true),
        padj=finite_or_missing([r.padj for r in rows]; clamp01=true),
        zero_inflated=[r.zero_inflated for r in rows],
        converged=[r.converged for r in rows])
end

Base.show(io::IO, r::DEResult) = print(io, analysis_result_summary(r))

"""
    GLMSolver(X; max_step_halvings=10)

Mutable workspace for IRLS fitting of NB-GLMs. Pre-allocates all buffers.
"""
mutable struct GLMSolver
    X::Matrix{Float64}
    xtwx::Matrix{Float64}
    xtwz::Vector{Float64}
    beta::Vector{Float64}
    eta::Vector{Float64}
    mu::Vector{Float64}
    weights::Vector{Float64}
    z::Vector{Float64}
    converged::Bool
    iterations::Int
    max_step_halvings::Int
    raw_xtwx::Matrix{Float64}

    function GLMSolver(X::AbstractMatrix{<:Real}; max_step_halvings::Int=10)
        Xf = Matrix{Float64}(X)
        rank(Xf) == size(Xf, 2) ||
            throw(ArgumentError("design matrix is not full rank"))
        nobs, ncoef = size(Xf)
        new(
            Xf,
            zeros(Float64, ncoef, ncoef),
            zeros(Float64, ncoef),
            zeros(Float64, ncoef),
            zeros(Float64, nobs),
            zeros(Float64, nobs),
            zeros(Float64, nobs),
            zeros(Float64, nobs),
            false,
            0,
            max_step_halvings,
            zeros(Float64, ncoef, ncoef))
    end
end

"""
    DESeqDataSet(counts, coldata, design)

Main container for DESeq2-style differential expression analysis.
"""
mutable struct DESeqDataSet
    counts::CountMatrix
    coldata::DataFrame
    design::Any
    size_factors::Union{Nothing,Vector{Float64}}
    normalization_factors::Union{Nothing,Vector{Float64}}
    dispersions::Union{Nothing,Vector{Float64}}
    gene_wise_dispersions::Union{Nothing,Vector{Float64}}
    dispersion_fit::Union{Nothing,Vector{Float64}}
    dispersion_prior_variance::Union{Nothing,Float64}
    cooks::Union{Nothing,Matrix{Float64}}
    replaced_counts::Union{Nothing,CountMatrix}
    wald_results::Union{Nothing,Vector{DEResult}}
    lrt_results::Union{Nothing,Vector{DEResult}}
    beta_prior::Bool
    beta_prior_var::Union{Nothing,Float64}
    model_matrix::Union{Nothing,Matrix{Float64}}
    reduced_model_matrix::Union{Nothing,Matrix{Float64}}
    coefficient_names::Union{Nothing,Vector{String}}
    wald_beta_matrix::Union{Nothing,Matrix{Float64}}
    wald_covariances::Union{Nothing,Vector{Matrix{Float64}}}
    wald_fisher_info::Union{Nothing,Vector{Matrix{Float64}}}
    test::Symbol
    fit_type::Symbol
    actual_fit_type::Union{Nothing,Symbol}
    sf_type::Symbol
    metadata::Dict{Any,Any}

    function DESeqDataSet(counts::CountMatrix, coldata::DataFrame, design)
        nrow(coldata) == length(counts.sample_ids) ||
            throw(ArgumentError("coldata must have one row per sample"))
        metadata = Dict{Any,Any}()
        ensure_provenance_id!(metadata)
        stamp_provenance!(
            metadata;
            label="DESeqDataSet",
            source="DifferentialExpression/DESeqDataSet",
            notes=["constructed DESeqDataSet"],
            parameters=(n_genes=length(counts.gene_ids), n_samples=length(counts.sample_ids)))
        new(
            counts, coldata, design,
            nothing, nothing, nothing, nothing, nothing, nothing,
            nothing, nothing, nothing, nothing,
            false, nothing, nothing, nothing, nothing, nothing,
            nothing, nothing,
            :none, :parametric, nothing, :ratio,
            metadata)
    end
end

function DESeqDataSetFromMatrix(count_data::AbstractMatrix{<:Integer},
    coldata::DataFrame, design; gene_ids=nothing, sample_ids=nothing)
    ng, ns = size(count_data)
    gids = gene_ids === nothing ? ["gene_$(i)" for i in 1:ng] : String.(gene_ids)
    sids = sample_ids === nothing ? ["sample_$(i)" for i in 1:ns] : String.(sample_ids)
    dds = DESeqDataSet(CountMatrix(count_data, gids, sids), coldata, design)
    update_provenance!(dds.metadata; source="DifferentialExpression/DESeqDataSetFromMatrix", notes=["constructed dataset from in-memory matrix"], parameters=(n_genes=ng, n_samples=ns))
    return dds
end

function Base.show(io::IO, dds::DESeqDataSet)
    print(io, "DESeqDataSet(", length(dds.counts.gene_ids), " genes, ", length(dds.counts.sample_ids), " samples, ", container_provenance_summary(dds), ")")
end

counts(dds::DESeqDataSet) = dds.counts
design(dds::DESeqDataSet) = dds.design
sizeFactors(dds::DESeqDataSet) = dds.size_factors
normalizationFactors(dds::DESeqDataSet) = dds.normalization_factors
dispersions(dds::DESeqDataSet) = dds.dispersions

# =============================================================================
# Reproducibility Infrastructure
# =============================================================================

"""
    ReproducibleDDS(dds)

Wraps a DESeqDataSet with input and parameter checksums for reproducibility.
"""
struct ReproducibleDDS
    dds::DESeqDataSet
    input_hash::String
    param_hash::String
    julia_version::String
    package_version::String
    timestamp::DateTime
end

"""
    input_checksum(dds::DESeqDataSet) -> String

Compute SHA256 checksum of count matrix and design for reproducibility.
"""
function input_checksum(dds::DESeqDataSet)
    buf = IOBuffer()
    write(buf, length(dds.counts.gene_ids), length(dds.counts.sample_ids))
    write(buf, dds.counts.counts.colval)
    write(buf, dds.counts.counts.rowval)
    for s in dds.counts.sample_ids
        write(buf, s)
    end
    for g in dds.counts.gene_ids
        write(buf, g)
    end
    write(buf, string(dds.design))

    return bytes2hex(SHA.hash(take!(buf), 256))
end

input_checksum(cm::CountMatrix) = input_checksum(
    DESeqDataSet(cm, DataFrame(condition=fill("A", length(cm.sample_ids))), ["A"]))

# =============================================================================
# Helper Functions
# =============================================================================

@inline function _safe_median(values::AbstractVector{<:Real})
    valid = filter(x -> isfinite(Float64(x)) && Float64(x) > 0.0, Float64.(values))
    isempty(valid) && return 1.0
    return median(valid)
end

@inline function _safe_mean(values::AbstractVector{<:Real})
    valid = filter(isfinite, Float64.(values))
    isempty(valid) && return NaN
    return mean(valid)
end

function _as_dense_row(mat::SparseMatrixCSC{Int,Int}, idx::Int)
    return Vector{Float64}(Array(mat[idx, :]))
end

function _as_dense_row(mat::AbstractMatrix{<:Integer}, idx::Int)
    return Float64.(vec(mat[idx, :]))
end

function _ensure_intercept(X::Matrix{Float64})
    n = size(X, 1)
    if size(X, 2) >= 1 && all(isapprox.(X[:, 1], 1.0; atol=1e-12, rtol=0.0))
        return X
    end
    return hcat(ones(Float64, n), X)
end

function _column_levels(design)
    if design isa AbstractVector
        return unique(String.(design))
    elseif design isa DataFrame && :condition in names(design)
        return unique(String.(design[!, :condition]))
    else
        return String[]
    end
end

function _design_model_matrix(design; modelMatrixType::Symbol=:standard)
    if design isa AbstractMatrix{<:Real}
        X = _ensure_intercept(Matrix{Float64}(design))
        cnames = String[]
        for j in 1:size(X, 2)
            push!(cnames, j == 1 ? "Intercept" : "coef_$(j)")
        end
        return X, cnames, Dict("type" => "matrix")
    elseif design isa DataFrame
        n = nrow(design)
        if modelMatrixType == :expanded
            X = zeros(Float64, n, 0)
            cnames = String[]
        else
            X = ones(Float64, n, 1)
            cnames = ["Intercept"]
        end
        for cname in names(design)
            vals = design[!, cname]
            if eltype(vals) <: Number
                X = hcat(X, Float64.(vals))
                push!(cnames, String(cname))
            else
                levs = unique(String.(vals))
                if length(levs) <= 1
                    continue
                end
                if modelMatrixType == :expanded
                    for lev in levs
                        X = hcat(X, Float64.(String.(vals) .== lev))
                        push!(cnames, "$(cname)_$(lev)")
                    end
                else
                    ref = levs[1]
                    for lev in levs[2:end]
                        X = hcat(X, Float64.(String.(vals) .== lev))
                        push!(cnames, "$(cname)_$(lev)_vs_$(ref)")
                    end
                end
            end
        end
        rank(X) == size(X, 2) || throw(ArgumentError("design matrix is not full rank"))
        return X, cnames, Dict("type" => "dataframe", "modelMatrixType" => String(modelMatrixType))
    elseif design isa AbstractVector
        vals = String.(design)
        n = length(vals)
        levs = unique(vals)
        if modelMatrixType == :expanded
            X = zeros(Float64, n, 0)
            cnames = String[]
        else
            X = ones(Float64, n, 1)
            cnames = ["Intercept"]
        end
        if length(levs) > 1
            ref = levs[1]
            if modelMatrixType == :expanded
                for lev in levs
                    X = hcat(X, Float64.(vals .== lev))
                    push!(cnames, "condition_$(lev)")
                end
            else
                for lev in levs[2:end]
                    X = hcat(X, Float64.(vals .== lev))
                    push!(cnames, "condition_$(lev)_vs_$(ref)")
                end
            end
        end
        rank(X) == size(X, 2) || throw(ArgumentError("design matrix is not full rank"))
        return X, cnames, Dict("type" => "vector", "levels" => levs, "modelMatrixType" => String(modelMatrixType))
    else
        throw(ArgumentError("unsupported design type $(typeof(design))"))
    end
end

function _named_contrast_vector(cnames::Vector{String}, design, contrast)
    if contrast isa AbstractVector && length(contrast) == 3 && !(contrast[1] isa Real)
        factor_name = String(contrast[1])
        numerator = String(contrast[2])
        denominator = String(contrast[3])
        levels = _column_levels(design)
        isempty(levels) && throw(ArgumentError("cannot infer factor levels from design"))
        ref = levels[1]
        cv = zeros(Float64, length(cnames))
        direct_num = findfirst(==("$(factor_name)_$(numerator)"), cnames)
        direct_den = findfirst(==("$(factor_name)_$(denominator)"), cnames)
        if direct_num !== nothing || direct_den !== nothing
            direct_num === nothing && throw(ArgumentError("contrast level $(numerator) not found"))
            direct_den === nothing && throw(ArgumentError("contrast level $(denominator) not found"))
            cv[direct_num] += 1.0
            cv[direct_den] -= 1.0
            return cv, "$(factor_name)_$(numerator)_vs_$(denominator)"
        end
        if numerator != ref
            key = "$(factor_name)_$(numerator)_vs_$(ref)"
            idx = findfirst(==(key), cnames)
            idx === nothing && throw(ArgumentError("contrast level $(numerator) not found"))
            cv[idx] += 1.0
        end
        if denominator != ref
            key = "$(factor_name)_$(denominator)_vs_$(ref)"
            idx = findfirst(==(key), cnames)
            idx === nothing && throw(ArgumentError("contrast level $(denominator) not found"))
            cv[idx] -= 1.0
        end
        return cv, "$(factor_name)_$(numerator)_vs_$(denominator)"
    elseif contrast isa AbstractVector && all(x -> x isa Real, contrast)
        length(contrast) == length(cnames) || throw(ArgumentError("contrast length mismatch"))
        return Float64.(contrast), "contrast"
    else
        throw(ArgumentError("unsupported contrast format"))
    end
end

function _default_contrast(cnames::Vector{String}, design=nothing)
    isempty(cnames) && return zeros(Float64, 0), ""
    if "Intercept" ∉ cnames && design !== nothing
        levels = _column_levels(design)
        if length(levels) >= 2
            try
                return _named_contrast_vector(cnames, design, [:condition, levels[end], levels[1]])
            catch
            end
        end
    end
    if length(cnames) == 1
        return [1.0], cnames[1]
    end
    idx = findlast(name -> name != "Intercept", cnames)
    idx === nothing && (idx = length(cnames))
    cv = zeros(Float64, length(cnames))
    cv[idx] = 1.0
    return cv, cnames[idx]
end

function _contrast_statistics(beta::AbstractVector{<:Real},
    covariance::AbstractMatrix{<:Real}, contrast::AbstractVector{<:Real})
    eff = dot(contrast, beta)
    var_eff = dot(contrast, covariance * contrast)
    var_eff = max(var_eff, eps(Float64))
    se = sqrt(var_eff)
    z = eff / se
    p = 2.0 * ccdf(Normal(), abs(z))
    # Convert to log2 scale
    return eff / log(2), se / log(2), abs(z), p
end

function _safe_inverse(A::AbstractMatrix{<:Real})
    Af = Matrix{Float64}(A)
    n = size(Af, 1)
    ridge = 1e-8
    for _ in 1:6
        M = Symmetric(Af + ridge * I(n))
        try
            return inv(cholesky(M))
        catch
            ridge *= 10
        end
    end
    return pinv(Af)
end

function _stabilize_beta_on_minmu_boundary(X::AbstractMatrix{<:Real},
    offset::AbstractVector{<:Real},
    beta::AbstractVector{<:Real},
    minmu::Real)
    eta = X * beta .+ offset
    eta_floor = log(max(Float64(minmu), eps(Float64)))
    target_eta = max.(eta, eta_floor)
    if maximum(abs.(target_eta .- eta)) <= 1e-8
        return Vector{Float64}(beta)
    end
    rhs = target_eta .- offset
    # Use Moore-Penrose projection unconditionally to avoid unstable least-
    # squares solutions in near-separation boundary cases.
    beta_proj = pinv(Matrix{Float64}(X)) * rhs
    return Vector{Float64}(beta_proj)
end

function _nb_loglik(y::AbstractVector{<:Real}, mu::AbstractVector{<:Real}, dispersion::Real)
    alpha = max(Float64(dispersion), eps(Float64))
    r = 1.0 / alpha
    ll = 0.0
    @inbounds for i in eachindex(y)
        mui = max(Float64(mu[i]), eps(Float64))
        yi = max(Float64(y[i]), 0.0)
        ll += lgamma(yi + r) - lgamma(r) - lgamma(yi + 1.0) +
              r * (log(r) - log(r + mui)) + yi * (log(mui) - log(r + mui))
    end
    return ll
end

function _logsumexp(values::AbstractVector{<:Real})
    isempty(values) && return -Inf
    m = maximum(values)
    !isfinite(m) && return m
    return m + log(sum(exp.(Float64.(values) .- m)))
end

# =============================================================================
# Design Matrix Validation (NEW)
# =============================================================================

"""
    validate_design(cm::CountMatrix, design; warn_confounding=true)

Validate design matrix for common issues: rank deficiency, confounding, low replication.
Returns (X, cnames, warnings).
"""
function validate_design(cm::CountMatrix, design; warn_confounding::Bool=true)
    X, cnames, _ = _design_model_matrix(design)
    n, p = size(X)
    warnings = String[]

    # Rank check
    r = rank(X)
    if r < p
        push!(warnings, "Design matrix is rank-deficient: rank=$r < cols=$p")
    end

    # Condition number
    cond = cond(X)
    if cond > 1e10
        push!(warnings, "Design matrix is ill-conditioned: κ=$(round(cond; sigdigits=3))")
    elseif cond > 1e6
        push!(warnings, "Design matrix has high condition number: κ=$(round(cond; sigdigits=3))")
    end

    # VIF for each non-intercept column
    if p > 1
        for j in 2:p
            other_idx = [i for i in 1:p if i != j]
            X_other = X[:, other_idx]
            xj = X[:, j]
            yhat = X_other * (X_other \ xj)
            ss_tot = sum((xj .- mean(xj)) .^ 2)
            ss_res = sum((xj .- yhat) .^ 2)
            r2 = 1.0 - ss_res / max(ss_tot, eps(Float64))
            r2 = clamp(r2, 0.0, 1.0 - 1e-12)
            vif = 1.0 / max(1.0 - r2, 1e-10)
            if vif > 10
                push!(warnings, "High VIF for $(cnames[j]): $(round(vif; digits=2))")
            end
        end
    end

    # Replication check
    if design isa AbstractVector
        groups = unique(String.(design))
        for g in groups
            n_g = count(==(g), String.(design))
            if n_g < 2
                push!(warnings, "Group '$g' has only $n_g sample(s) — no replication")
            elseif n_g < 3
                push!(warnings, "Group '$g' has only $n_g samples — dispersion estimation may be unreliable")
            end
        end
    end

    return X, cnames, warnings
end

# =============================================================================
# Multiple Testing Corrections
# =============================================================================

"""
    benjamini_hochberg(pvalues) -> Vector{Float64}

Benjamini-Hochberg procedure for FDR control.
"""
function benjamini_hochberg(pvalues::AbstractVector)
    n = length(pvalues)
    out = fill(1.0, n)
    finite_mask = [(!ismissing(pvalues[i])) && isfinite(Float64(pvalues[i])) for i in 1:n]
    idx = findall(identity, finite_mask)
    isempty(idx) && return out

    pv = Float64[pvalues[i] for i in idx]
    ord = sortperm(pv)
    pv_sorted = pv[ord]
    m = length(pv_sorted)
    q = similar(pv_sorted)

    q[end] = min(1.0, pv_sorted[end])
    for i in (m-1):-1:1
        q[i] = min(q[i+1], pv_sorted[i] * m / i)
    end

    for k in 1:m
        out[idx[ord[k]]] = clamp(q[k], 0.0, 1.0)
    end

    return out
end

function _natural_cubic_spline_eval(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq::Real)
    xs = Float64.(x)
    ys = Float64.(y)
    n = length(xs)
    n == length(ys) || throw(ArgumentError("spline x/y length mismatch"))
    n == 0 && return NaN
    n == 1 && return ys[1]
    xpos = Float64(xq)
    if n == 2
        t = (xpos - xs[1]) / max(xs[2] - xs[1], eps(Float64))
        return ys[1] * (1.0 - t) + ys[2] * t
    end

    h = diff(xs)
    any(h .<= 0.0) && throw(ArgumentError("spline x values must be strictly increasing"))

    a = zeros(Float64, n)
    b = zeros(Float64, n)
    c = zeros(Float64, n)
    rhs = zeros(Float64, n)
    for i in 2:(n-1)
        a[i] = h[i-1]
        b[i] = 2.0 * (h[i-1] + h[i])
        c[i] = h[i]
        rhs[i] = 6.0 * ((ys[i+1] - ys[i]) / h[i] - (ys[i] - ys[i-1]) / h[i-1])
    end

    cp = zeros(Float64, n)
    dp = zeros(Float64, n)
    for i in 2:(n-1)
        denom = b[i] - a[i] * cp[i-1]
        denom = abs(denom) < eps(Float64) ? eps(Float64) : denom
        cp[i] = c[i] / denom
        dp[i] = (rhs[i] - a[i] * dp[i-1]) / denom
    end

    m = zeros(Float64, n)
    for i in (n-1):-1:2
        m[i] = dp[i] - cp[i] * m[i+1]
    end

    # Natural cubic spline boundary condition implies linear extrapolation
    # beyond boundary knots, not clamping to endpoint values.
    if xpos <= xs[1]
        slope_left = (ys[2] - ys[1]) / h[1] - h[1] * (2m[1] + m[2]) / 6.0
        return ys[1] + slope_left * (xpos - xs[1])
    elseif xpos >= xs[end]
        slope_right = (ys[end] - ys[end-1]) / h[end] + h[end] * (m[end-1] + 2m[end]) / 6.0
        return ys[end] + slope_right * (xpos - xs[end])
    end

    idx = clamp(searchsortedlast(xs, xpos), 1, n - 1)
    hi = h[idx]
    t = (xpos - xs[idx]) / hi
    a0 = 1.0 - t
    b0 = t
    return a0 * ys[idx] + b0 * ys[idx+1] +
           ((a0^3 - a0) * m[idx] + (b0^3 - b0) * m[idx+1]) * (hi^2) / 6.0
end

"""
    storey_qvalue(pvalues; lambda_grid=0.5:0.05:0.95) -> Vector{Float64}

Storey's q-value estimation using π₀ estimation at λ=1.
"""
function storey_qvalue(pvalues::AbstractVector; lambda_grid=0.5:0.05:0.95)
    n = length(pvalues)
    p = Float64[]
    pos = Int[]
    for i in 1:n
        if !ismissing(pvalues[i])
            v = Float64(pvalues[i])
            if isfinite(v)
                push!(p, clamp(v, 0.0, 1.0))
                push!(pos, i)
            end
        end
    end
    out = fill(1.0, n)
    isempty(p) && return out

    pi0_vals = Float64[]
    lambdas = Float64[]
    for lam in lambda_grid
        lamf = clamp(Float64(lam), 0.0, 0.99)
        push!(lambdas, lamf)
        push!(pi0_vals, clamp(count(x -> x > lamf, p) / ((1.0 - lamf) * length(p)), 0.0, 1.0))
    end
    # Storey-style pi0(lambda) should be non-decreasing in lambda.
    for i in 2:length(pi0_vals)
        pi0_vals[i] = max(pi0_vals[i], pi0_vals[i-1])
    end
    pi0 = clamp(_natural_cubic_spline_eval(lambdas, pi0_vals, 1.0), 0.0, 1.0)

    ord = sortperm(p)
    ps = p[ord]
    m = length(ps)
    q = similar(ps)
    q[end] = min(1.0, pi0 * ps[end])
    for i in (m-1):-1:1
        q[i] = min(q[i+1], pi0 * ps[i] * m / i)
    end
    for i in 1:m
        out[pos[ord[i]]] = clamp(q[i], 0.0, 1.0)
    end

    return out
end

"""
    storey_pi0(pvalues; lambda_grid=0.5:0.05:0.95) -> Float64

Storey's π₀ estimation using spline extrapolation to λ=1.
Returns the estimated proportion of true null hypotheses.
"""
function storey_pi0(pvalues::AbstractVector; lambda_grid=0.5:0.05:0.95)
    p = Float64[]
    for i in 1:length(pvalues)
        if !ismissing(pvalues[i])
            v = Float64(pvalues[i])
            if isfinite(v)
                push!(p, clamp(v, 0.0, 1.0))
            end
        end
    end
    isempty(p) && return 1.0

    pi0_vals = Float64[]
    lambdas = Float64[]
    for lam in lambda_grid
        lamf = clamp(Float64(lam), 0.0, 0.99)
        push!(lambdas, lamf)
        push!(pi0_vals, clamp(count(x -> x > lamf, p) / ((1.0 - lamf) * length(p)), 0.0, 1.0))
    end

    # Storey-style pi0(lambda) should be non-decreasing in lambda.
    for i in 2:length(pi0_vals)
        pi0_vals[i] = max(pi0_vals[i], pi0_vals[i-1])
    end

    return clamp(_natural_cubic_spline_eval(lambdas, pi0_vals, 1.0), 0.0, 1.0)
end

function _weighted_bh(p::AbstractVector{<:Real}, w::AbstractVector{<:Real})
    m = length(p)
    m == length(w) || throw(ArgumentError("p/w length mismatch"))
    m == 0 && return Float64[]

    pf = clamp.(Float64.(p), 0.0, 1.0)
    wf = max.(Float64.(w), 1e-8)
    u = pf ./ wf
    ord = sortperm(u)
    q_sorted = zeros(Float64, m)

    q_sorted[end] = min(1.0, u[ord[end]])
    for i in (m-1):-1:1
        q_sorted[i] = min(q_sorted[i+1], (m / i) * u[ord[i]])
    end

    q = similar(q_sorted)
    for i in 1:m
        q[ord[i]] = clamp(q_sorted[i], 0.0, 1.0)
    end
    return q
end

"""
    ihw_qvalue(pvalues, filter_stat; alpha=0.1, n_bins=20) -> Vector{Float64}

Independent hypothesis weighting (Ignatiadis et al. 2016).
"""
function ihw_qvalue(pvalues::AbstractVector, filter_stat::AbstractVector;
    alpha::Real=0.1,
    n_bins::Int=20,
    n_folds::Int=5,
    lambda::Real=0.5,
    min_bin_size::Int=10)

    n = length(pvalues)
    length(filter_stat) == n || throw(ArgumentError("pvalues/filter_stat length mismatch"))
    out = fill(1.0, n)
    n_bins = clamp(Int(n_bins), 2, max(2, n))
    n_folds = clamp(Int(n_folds), 2, 20)
    lambda_f = clamp(Float64(lambda), 0.1, 0.95)

    valid = Int[]
    for i in 1:n
        if !ismissing(pvalues[i]) && !ismissing(filter_stat[i])
            p = Float64(pvalues[i])
            f = Float64(filter_stat[i])
            if isfinite(p) && isfinite(f)
                push!(valid, i)
            end
        end
    end
    isempty(valid) && return out

    m = length(valid)
    pv = clamp.(Float64[pvalues[i] for i in valid], 0.0, 1.0)
    fv = Float64[filter_stat[i] for i in valid]
    edges = [quantile(fv, q) for q in range(0.0, 1.0; length=n_bins + 1)]

    function _bin_of(x::Float64)
        b = searchsortedlast(edges, x)
        return clamp(b, 1, n_bins)
    end

    bins = [_bin_of(fv[i]) for i in 1:m]

    ord_f = sortperm(fv)
    folds = zeros(Int, m)
    for (rank, idx) in enumerate(ord_f)
        folds[idx] = 1 + mod(rank - 1, n_folds)
    end

    weights = ones(Float64, m)
    for fold in 1:n_folds
        train = findall(i -> folds[i] != fold, 1:m)
        test = findall(i -> folds[i] == fold, 1:m)
        isempty(test) && continue

        raw_w = ones(Float64, n_bins)
        for b in 1:n_bins
            bin_train = [i for i in train if bins[i] == b]
            mb = length(bin_train)
            if mb < min_bin_size
                raw_w[b] = 1.0
                continue
            end
            pi0_hat = count(i -> pv[i] > lambda_f, bin_train) / ((1.0 - lambda_f) * mb)
            pi0_hat = clamp(pi0_hat, 1e-3, 1.0)
            raw_w[b] = 1.0 / pi0_hat
        end

        train_mean = mean(raw_w[bins[i]] for i in train)
        train_mean = isfinite(train_mean) && train_mean > 0.0 ? train_mean : 1.0
        raw_w ./= train_mean

        for i in test
            weights[i] = raw_w[bins[i]]
        end
    end

    wmean = mean(weights)
    wmean = isfinite(wmean) && wmean > 0.0 ? wmean : 1.0
    weights ./= wmean

    q_valid = _weighted_bh(pv, weights)
    for (k, i) in enumerate(valid)
        out[i] = q_valid[k]
    end

    return out
end

# =============================================================================
# Empirical Null Calibration (NEW)
# =============================================================================

"""
    calibrate_pvalues_empirical_null(pvalues; method=:efron, df=7) -> Vector{Float64}

Efron's empirical null calibration for correlated tests.
Fits N(δ, σ²) to central z-scores and recalibrates p-values.
"""
function calibrate_pvalues_empirical_null(pvalues::AbstractVector;
    method::Symbol=:efron,
    df::Real=7)
    n = length(pvalues)
    out = fill(1.0, n)

    # Convert to z-scores
    z = Float64[]
    pos = Int[]
    for i in 1:n
        if !ismissing(pvalues[i]) && isfinite(Float64(pvalues[i]))
            p = Float64(pvalues[i])
            p = clamp(p, 1e-300, 1.0 - 1e-15)
            push!(z, quantile(Normal(), p > 0.5 ? 1.0 - p : p) * (p > 0.5 ? -1.0 : 1.0))
            push!(pos, i)
        end
    end
    isempty(z) && return out

    if method == :efron
        # Fit N(δ, σ²) to central portion using spline density estimation
        z_sorted = sort(z)
        m = length(z_sorted)

        # Use central 50% for null estimation
        lo_idx = max(1, floor(Int, 0.25 * m))
        hi_idx = min(m, ceil(Int, 0.75 * m))
        z_central = z_sorted[lo_idx:hi_idx]

        # Method of moments for N(δ, σ²)
        δ = mean(z_central)
        σ = std(z_central; corrected=false)
        σ = max(σ, 0.1)  # Guard against σ → 0

        # Recalibrate z-scores
        z_calibrated = (z .- δ) ./ σ
        p_calibrated = 2.0 * ccdf(Normal(), abs.(z_calibrated))

        for (k, i) in enumerate(pos)
            out[i] = clamp(p_calibrated[k], 0.0, 1.0)
        end
    elseif method == :spline
        # Spline-based density fitting
        z_min, z_max = minimum(z), maximum(z)
        grid = range(z_min - 0.5, z_max + 0.5; length=100)
        density_est = kde_spline(z, grid, df=Int(df))

        # Find mode (likely null center)
        mode_idx = argmax(density_est)
        z_mode = grid[mode_idx]
        density_mode = density_est[mode_idx]

        # Estimate σ from curvature at mode
        if mode_idx > 1 && mode_idx < length(grid)
            d2 = (density_est[mode_idx+1] - 2 * density_est[mode_idx] + density_est[mode_idx-1]) /
                 (grid[2] - grid[1])^2
            if d2 < 0
                σ = sqrt(max(-density_mode / d2, 0.01))
            else
                σ = 1.0
            end
        else
            σ = 1.0
        end

        z_calibrated = (z .- z_mode) ./ σ
        p_calibrated = 2.0 * ccdf(Normal(), abs.(z_calibrated))

        for (k, i) in enumerate(pos)
            out[i] = clamp(p_calibrated[k], 0.0, 1.0)
        end
    end

    return out
end

function kde_spline(x::AbstractVector{<:Real}, grid::AbstractVector{<:Real}; df::Int=7)
    n = length(x)
    h = 1.06 * std(x) * n^(-1 / 5)  # Silverman's rule
    density = zeros(Float64, length(grid))

    for i in eachindex(grid)
        for j in eachindex(x)
            density[i] += exp(-0.5 * ((grid[i] - x[j]) / h)^2)
        end
    end
    density ./= (n * h * sqrt(2π))

    return density
end

# =============================================================================
# Normalization
# =============================================================================

function _validate_control_genes(controlGenes, n_genes::Int)
    if controlGenes === nothing
        return collect(1:n_genes)
    elseif controlGenes isa AbstractVector{Bool}
        length(controlGenes) == n_genes ||
            throw(ArgumentError("controlGenes Boolean mask length must equal number of genes"))
        idx = findall(identity, controlGenes)
        isempty(idx) && throw(ArgumentError("controlGenes mask selected zero genes"))
        return idx
    elseif controlGenes isa AbstractVector{<:Integer}
        idx = Int.(controlGenes)
        isempty(idx) && throw(ArgumentError("controlGenes selected zero genes"))
        any(i -> i < 1 || i > n_genes, idx) &&
            throw(ArgumentError("controlGenes index out of bounds"))
        return unique(idx)
    else
        throw(ArgumentError("controlGenes must be nothing, a Bool vector, or an index vector"))
    end
end

function estimateSizeFactorsForMatrix(counts::AbstractMatrix{<:Integer};
    type::Symbol=:ratio,
    geoMeans=nothing,
    controlGenes=nothing,
    fallback::Symbol=:library_size)

    Y = Matrix{Float64}(counts)
    ng, ns = size(Y)
    idx = _validate_control_genes(controlGenes, ng)

    gm = if geoMeans === nothing
        out = zeros(Float64, ng)
        if type == :ratio
            for g in 1:ng
                row = @view Y[g, :]
                if all(x -> x > 0.0, row)
                    out[g] = exp(mean(log.(row)))
                else
                    out[g] = 0.0
                end
            end
        elseif type == :poscounts
            for g in 1:ng
                row = @view Y[g, :]
                pos = row[row.>0.0]
                if isempty(pos)
                    out[g] = 0.0
                else
                    out[g] = exp(sum(log.(pos)) / ns)
                end
            end
        else
            throw(ArgumentError("unsupported size factor type: $type"))
        end
        out
    else
        g = Float64.(geoMeans)
        length(g) == ng || throw(ArgumentError("geoMeans length mismatch"))
        g
    end

    valid_genes = [g for g in idx if gm[g] > 0.0 && isfinite(gm[g])]
    if type == :ratio && isempty(valid_genes)
        throw(ArgumentError("no positive geometric means found; use type=:poscounts for sparse count data"))
    end

    sf = zeros(Float64, ns)
    lib_sizes = max.(vec(sum(Y, dims=1)), 1.0)
    for s in 1:ns
        ratios = Float64[]
        for g in valid_genes
            y = Y[g, s]
            if y > 0.0
                push!(ratios, y / gm[g])
            end
        end
        if isempty(ratios)
            if fallback == :throw
                throw(ArgumentError("control genes do not provide usable ratios for sample $(s)"))
            elseif fallback == :library_size
                @warn "control genes do not provide usable ratios for sample $(s); falling back to library-size normalization"
                sf[s] = lib_sizes[s]
            else
                throw(ArgumentError("unsupported fallback: $fallback"))
            end
        else
            sf[s] = median(ratios)
        end
    end

    all(x -> x > 0.0 && isfinite(x), sf) ||
        throw(ArgumentError("failed to estimate positive finite size factors"))

    sf ./= exp(mean(log.(sf)))
    return sf
end

function calc_tmm_factors(cm::CountMatrix;
    trim_m::Real=0.30,
    trim_a::Real=0.05,
    ref_column::Union{Nothing,Int}=nothing)

    Y = Matrix{Float64}(cm.counts)
    ns = size(Y, 2)
    ns == 0 && return Float64[]

    keep = vec(sum(Y, dims=2) .> 0.0)
    if !any(keep)
        return ones(Float64, ns)
    end
    Y = Y[keep, :]

    lib = vec(sum(Y, dims=1))
    if all(lib .<= 0.0)
        return ones(Float64, ns)
    end
    lib = max.(lib, 1.0)

    f75 = [quantile(Y[:, j] ./ lib[j], 0.75) for j in 1:ns]
    ref = ref_column === nothing ? findmin(abs.(f75 .- median(f75)))[2] : ref_column
    (1 <= ref <= ns) || throw(ArgumentError("ref_column out of bounds"))

    factors = ones(Float64, ns)
    yref = Y[:, ref]
    libref = lib[ref]

    for j in 1:ns
        j == ref && continue
        yj = Y[:, j]
        libj = lib[j]
        mask = (yj .> 0.0) .& (yref .> 0.0)
        if count(mask) < 10
            @warn "TMM normalization fell back to library-size scaling for sample $(j)"
            factors[j] = libj
            continue
        end

        pj = yj[mask] ./ libj
        pr = yref[mask] ./ libref
        M = log2.(pj ./ pr)
        A = 0.5 .* log2.(pj .* pr)
        var_M = ((libj .- yj[mask]) ./ (libj .* yj[mask])) .+
                ((libref .- yref[mask]) ./ (libref .* yref[mask]))
        w = 1.0 ./ max.(var_M, 1e-6)

        mlo = quantile(M, trim_m / 2)
        mhi = quantile(M, 1.0 - trim_m / 2)
        alo = quantile(A, trim_a / 2)
        ahi = quantile(A, 1.0 - trim_a / 2)
        keep2 = (M .>= mlo) .& (M .<= mhi) .& (A .>= alo) .& (A .<= ahi)

        if count(keep2) > 0
            wm = sum(w[keep2] .* M[keep2]) / max(sum(w[keep2]), eps(Float64))
            factors[j] = 2.0^wm
        else
            @warn "TMM normalization fell back to library-size scaling for sample $(j)"
            factors[j] = libj
        end
    end

    factors ./= exp(mean(log.(factors)))
    return factors
end

function calc_norm_factors(cm::CountMatrix;
    method::Symbol=:tmm,
    type::Symbol=:ratio)

    result = if method == :tmm
        calc_tmm_factors(cm)
    elseif method == :ratio || method == :poscounts
        estimateSizeFactorsForMatrix(cm.counts; type=method)
    elseif method == :none
        ones(Float64, length(cm.sample_ids))
    else
        throw(ArgumentError("unsupported normalization method: $method"))
    end

    _ctx = active_provenance_context()
    _ctx !== nothing && register_provenance!(_ctx, "calc_norm_factors"; parents=provenance_parent_ids(cm), parameters=(method=method, type=type, n_samples=length(cm.sample_ids), n_genes=length(cm.gene_ids)))
    return result
end

# =============================================================================
# Dispersion Estimation
# =============================================================================

function _moment_dispersion(y::AbstractVector{<:Real}, nf::AbstractVector{<:Real})
    yn = Float64.(y) ./ Float64.(nf)
    mu = mean(yn)
    if !isfinite(mu) || mu <= 0.0
        return 1e-8
    end
    v = var(yn; corrected=true)
    if !isfinite(v)
        return 1e-8
    end
    inv_nf_mean = mean(1.0 ./ max.(Float64.(nf), eps(Float64)))
    alpha = (v - mu * inv_nf_mean) / max(mu^2, eps(Float64))
    return clamp(alpha, 1e-8, 1e3)
end

function _linear_model_mu(norm_counts::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real})
    Xf = Matrix{Float64}(X)
    Y = Matrix{Float64}(norm_counts)
    coef = try
        Xf \ transpose(Y)
    catch
        pinv(Xf) * transpose(Y)
    end
    fitted = Xf * coef
    return transpose(fitted)
end

function _rough_dispersion_estimates(norm_counts::AbstractMatrix{<:Real}, X::AbstractMatrix{<:Real})
    mu = _linear_model_mu(norm_counts, X)
    mu .= max.(mu, 1.0)
    m = size(X, 1)
    p = size(X, 2)
    df = max(m - p, 1)
    est = vec(sum(((norm_counts .- mu) .^ 2 .- mu) ./ max.(mu .^ 2, eps(Float64)), dims=2)) ./ df
    return max.(est, 0.0)
end

function _moments_dispersion_estimates(norm_counts::AbstractMatrix{<:Real}, nf::AbstractVector{<:Real})
    xim = mean(1.0 ./ max.(Float64.(nf), eps(Float64)))
    base_mean = vec(mean(norm_counts, dims=2))
    base_var = vec(var(norm_counts, dims=2; corrected=true))
    return (base_var .- xim .* base_mean) ./ max.(base_mean .^ 2, eps(Float64))
end

function _fit_gene_dispersion_mle(y::AbstractVector{<:Real},
    nf::AbstractVector{<:Real},
    X::AbstractMatrix{<:Real};
    alpha_init::Real=0.1,
    alpha_min::Real=1e-8,
    alpha_max::Real=max(10.0, size(X, 1)),
    maxit::Integer=50,
    beta_tol::Real=1e-8,
    minmu::Real=0.5)

    yv = Float64.(y)
    sum(yv) <= 0.0 && return alpha_min

    solver = GLMSolver(X)
    offset = log.(Float64.(nf))
    loglo = log(alpha_min)
    loghi = log(alpha_max)
    alpha0 = clamp(Float64(alpha_init), alpha_min, alpha_max)
    best_alpha = alpha0
    best_obj = Inf

    # Fallback mean model used only if GLM fails for a trial alpha.
    # Keep this in count scale (not normalized-count scale) because
    # `_nb_loglik` is parameterized on raw count means.
    local_norm_mean = max(mean(yv ./ max.(Float64.(nf), eps(Float64))), alpha_min)
    # For very low-mean genes initialized at the dispersion floor, a fixed high
    # minmu can force the profile objective toward spuriously large alpha
    # values. Keep boundary protection, but only relax this clamp in that
    # near-boundary regime.
    disp_minmu = if alpha0 <= 10.0 * alpha_min
        min(Float64(minmu), max(local_norm_mean / 5.0, 1e-6))
    else
        Float64(minmu)
    end
    fallback_mu = max.(Float64.(nf) .* local_norm_mean, disp_minmu)

    function objective(logalpha)
        alpha = exp(logalpha)

        # Cox-Reid adjusted profile likelihood requires profiling over beta
        # at each alpha. Reusing a fixed mu across alpha values is not a
        # valid profile objective and can bias dispersions toward boundaries.
        local mu
        try
            fit_gene_fast!(solver, yv, offset, alpha;
                maxit=maxit, beta_tol=beta_tol, minmu=disp_minmu)
            if !solver.converged
                fit_gene_fast!(solver, yv, offset, alpha;
                    maxit=max(maxit + 25, 2 * Int(maxit)),
                    beta_tol=max(Float64(beta_tol), 1e-7),
                    minmu=disp_minmu)
            end
            mu = max.(solver.mu, disp_minmu)
        catch
            mu = fallback_mu
        end

        ll = _nb_loglik(yv, mu, alpha)
        if !isfinite(ll)
            return Inf
        end

        w = mu ./ (1.0 .+ alpha .* mu)
        fisher = transpose(X) * Diagonal(w) * X
        _, logabsdet_fisher = logabsdet(Symmetric(fisher))
        if !isfinite(logabsdet_fisher)
            _, logabsdet_fisher = logabsdet(Symmetric(fisher + 1e-8I(size(fisher, 1))))
        end
        if !isfinite(logabsdet_fisher)
            return Inf
        end
        obj = -(ll - 0.5 * logabsdet_fisher)
        if obj < best_obj
            best_obj = obj
            best_alpha = alpha
        end
        return obj
    end

    # Coarse grid pre-search ensures we always have a finite fallback candidate
    # if Brent encounters non-finite objective regions.
    candidate_logs = unique(vcat([log(alpha0), loglo, loghi], collect(range(loglo, loghi; length=13))))
    for la in candidate_logs
        objective(la)
    end

    local opt
    try
        opt = Optim.optimize(objective, loglo, loghi, Optim.Brent(); iterations=40)
        if Optim.converged(opt) && isfinite(Optim.minimum(opt))
            best_alpha = exp(Optim.minimizer(opt))
        end
    catch
    end

    # Boundary checks are cheap and important when the likelihood surface is flat.
    objective(loglo)
    objective(loghi)

    if !isfinite(best_obj)
        return clamp(_moment_dispersion(y, nf), alpha_min, alpha_max)
    end
    return clamp(best_alpha, alpha_min, alpha_max)
end

function _local_dispersion_smoother(means::Vector{Float64},
    disps::Vector{Float64};
    span::Real=0.5,
    min_disp::Real=1e-8,
    fit_mask::Union{Nothing,AbstractVector{Bool}}=nothing)

    n = length(means)
    n == length(disps) || throw(ArgumentError("means/disps length mismatch"))

    mask = isfinite.(means) .& isfinite.(disps) .& (means .> 0.0) .& (disps .> 0.0)
    if fit_mask !== nothing
        length(fit_mask) == n || throw(ArgumentError("fit_mask length mismatch"))
        mask .&= fit_mask
    else
        mask .&= disps .>= 10.0 * min_disp
    end

    if count(mask) < 5
        fillv = _safe_median(disps[isfinite.(disps).&(disps.>0.0)])
        return fill(max(fillv, min_disp), n)
    end

    x_fit = log.(means[mask])
    y_fit = log.(max.(disps[mask], min_disp))
    w_fit = max.(means[mask], eps(Float64))

    order = sortperm(x_fit)
    x_fit = x_fit[order]
    y_fit = y_fit[order]
    w_fit = w_fit[order]

    m = length(x_fit)
    k = clamp(round(Int, span * m), 5, m)
    yhat = similar(y_fit)

    function _fit_local_linear!(out::Vector{Float64}, robust_weights::Vector{Float64})
        @inbounds for i in 1:m
            lo = max(1, i - k ÷ 2)
            hi = min(m, lo + k - 1)
            lo = max(1, hi - k + 1)

            xw = @view x_fit[lo:hi]
            yw = @view y_fit[lo:hi]
            ww = @view w_fit[lo:hi]
            rw = @view robust_weights[lo:hi]
            x0 = x_fit[i]

            d = abs.(xw .- x0)
            dmax = maximum(d)
            kernel = if dmax <= eps(Float64)
                ones(Float64, length(d))
            else
                u = d ./ dmax
                (1 .- u .^ 3) .^ 3
            end
            w = kernel .* ww .* rw
            sw = sum(w)

            if sw <= eps(Float64)
                out[i] = y_fit[i]
                continue
            end

            xd = xw .- x0
            sx = sum(w .* xd)
            sxx = sum(w .* xd .* xd)
            sy = sum(w .* yw)
            sxy = sum(w .* xd .* yw)
            denom = sw * sxx - sx^2

            if abs(denom) <= eps(Float64)
                out[i] = sy / sw
            else
                out[i] = (sxx * sy - sx * sxy) / denom
            end
        end
        return out
    end

    # Robust reweighting (loess-style bisquare) stabilizes the local trend
    # against high-dispersion outliers that can distort smooth fits.
    robust_w = ones(Float64, m)
    for _ in 1:3
        _fit_local_linear!(yhat, robust_w)
        resid = y_fit .- yhat
        med = median(resid)
        mad = median(abs.(resid .- med))
        scale = max(1.4826 * mad, 1e-6)
        u = abs.(resid) ./ (6.0 * scale)
        robust_w .= (1 .- min.(u, 1.0) .^ 2) .^ 2
    end
    _fit_local_linear!(yhat, robust_w)

    out = fill(max(_safe_median(disps[mask]), min_disp), n)
    pred_mask = isfinite.(means) .& (means .> 0.0)
    x_all = log.(means[pred_mask])

    function interp_logfit(xv::Float64)
        if xv <= x_fit[1]
            return yhat[1]
        elseif xv >= x_fit[end]
            return yhat[end]
        else
            j = clamp(searchsortedlast(x_fit, xv), 1, m - 1)
            x1 = x_fit[j]
            x2 = x_fit[j+1]
            if x2 <= x1 + eps(Float64)
                return yhat[j]
            end
            t = (xv - x1) / (x2 - x1)
            return (1.0 - t) * yhat[j] + t * yhat[j+1]
        end
    end

    preds = similar(x_all)
    @inbounds for i in eachindex(x_all)
        preds[i] = exp(interp_logfit(x_all[i]))
    end
    out[pred_mask] = max.(preds, min_disp)
    return out
end

function _trend_fit_score(gene_wise::Vector{Float64},
    trend::Vector{Float64},
    fit_mask::AbstractVector{Bool};
    min_disp::Real=1e-8)

    n = length(gene_wise)
    n == length(trend) || throw(ArgumentError("gene_wise/trend length mismatch"))
    length(fit_mask) == n || throw(ArgumentError("fit_mask length mismatch"))

    valid = fit_mask .& isfinite.(gene_wise) .& isfinite.(trend) .&
            (gene_wise .> 0.0) .& (trend .> 0.0)
    if count(valid) < 5
        return (Inf, Inf, Inf, Inf)
    end

    residual = log.(max.(gene_wise[valid], min_disp)) .-
               log.(max.(trend[valid], min_disp))
    abs_residual = abs.(residual)
    med_abs = median(abs_residual)
    bias = abs(median(residual))
    tail_q90 = quantile(abs_residual, 0.9)
    score = med_abs + 0.2 * bias + 0.1 * tail_q90

    return (score, med_abs, bias, tail_q90)
end

function _adaptive_local_dispersion_trend(means::Vector{Float64},
    gene_wise::Vector{Float64};
    fit_mask::AbstractVector{Bool},
    min_disp::Real=1e-8)

    n = length(means)
    n == length(gene_wise) || throw(ArgumentError("means/gene_wise length mismatch"))
    length(fit_mask) == n || throw(ArgumentError("fit_mask length mismatch"))

    nfit = count(fit_mask)
    span_pool = nfit >= 20_000 ? (0.30, 0.40, 0.50) : (0.25, 0.35, 0.45, 0.55, 0.65)
    center_span = 0.35
    score_tolerance = 0.005 + 0.05 / sqrt(max(nfit, 1))

    best_score = Inf
    best_span = 0.0
    best_trend = fill(max(_safe_median(gene_wise[fit_mask]), min_disp), n)
    span_scores = Dict{String,Float64}()
    candidates = Tuple{Float64,Float64,Vector{Float64}}[]

    for span in span_pool
        candidate = _local_dispersion_smoother(means, gene_wise;
            span=span,
            min_disp=min_disp,
            fit_mask=fit_mask)
        candidate .= max.(candidate, min_disp)
        score, _, _, _ = _trend_fit_score(gene_wise, candidate, fit_mask; min_disp=min_disp)
        span_scores[string(span)] = score
        if isfinite(score)
            push!(candidates, (Float64(span), score, candidate))
            if score < best_score
                best_score = score
            end
        end
    end

    if !isempty(candidates)
        near = [c for c in candidates if c[2] <= best_score + score_tolerance]
        chosen = isempty(near) ? candidates[argmin(c[2] for c in candidates)] : near[1]
        if !isempty(near)
            chosen_dist = abs(chosen[1] - center_span)
            for cand in @view near[2:end]
                d = abs(cand[1] - center_span)
                if d < chosen_dist
                    chosen = cand
                    chosen_dist = d
                end
            end
        end
        best_span = chosen[1]
        best_score = chosen[2]
        best_trend = chosen[3]
    end

    if !isfinite(best_score)
        best_span = first(span_pool)
        best_trend = _local_dispersion_smoother(means, gene_wise;
            span=best_span,
            min_disp=min_disp,
            fit_mask=fit_mask)
        best_trend .= max.(best_trend, min_disp)
        best_score, _, _, _ = _trend_fit_score(gene_wise, best_trend, fit_mask; min_disp=min_disp)
    end

    return best_trend, best_span, best_score, span_scores
end

function _parametric_dispersion_fit(means::Vector{Float64}, disps::Vector{Float64};
    maxiter::Int=10)
    valid = isfinite.(means) .& isfinite.(disps) .& (means .> 0.0) .& (disps .> 0.0)
    count(valid) < 2 && return nothing

    x = means[valid]
    y = disps[valid]
    X = hcat(ones(Float64, length(x)), 1.0 ./ x)

    coef = try
        X \ y
    catch
        Float64[]
    end
    if length(coef) != 2 || any(!isfinite, coef) || any(!isfinite, X * coef) || any((X * coef) .<= 0.0)
        coef = [median(y), 0.0]
    end
    coef = Float64.(coef)
    local_trend, _, local_err, _ = _adaptive_local_dispersion_trend(x, y;
        fit_mask=trues(length(x)),
        min_disp=1e-8)
    local_trend .= max.(local_trend, 1e-8)
    if !isfinite(local_err)
        local_err = median(abs.(log.(y) .- log.(local_trend)))
    end

    for _ in 1:maxiter
        fitted = X * coef
        ratio = y ./ max.(fitted, 1e-12)
        good = isfinite.(fitted) .& (fitted .> 0.0) .& isfinite.(ratio) .& (ratio .> 1e-4) .& (ratio .< 15.0)
        count(good) < 2 && return nothing

        Xg = X[good, :]
        yg = y[good]
        μ = max.(Xg * coef, 1e-12)
        w = 1.0 ./ max.(μ .^ 2, 1e-12)
        sw = sqrt.(w)

        oldcoef = copy(coef)
        newcoef = try
            (Xg .* sw) \ (yg .* sw)
        catch
            return nothing
        end
        if length(newcoef) != 2 || any(!isfinite, newcoef)
            return nothing
        end
        coef .= newcoef

        fitted_new = X * coef
        fitted_old = X * oldcoef
        if all(isfinite, fitted_new) && all(fitted_new .> 0.0) &&
           maximum(abs.(fitted_new .- fitted_old) ./ max.(abs.(fitted_old), 1.0)) < 1e-6
            param_err, _, _, _ = _trend_fit_score(y, fitted_new, trues(length(y)); min_disp=1e-8)
            return isfinite(param_err) && param_err <= 1.1 * local_err ? coef : nothing
        end
    end

    fitted = X * coef
    if !all(isfinite, fitted) || any(fitted .<= 0.0)
        return nothing
    end

    param_err, _, _, _ = _trend_fit_score(y, fitted, trues(length(y)); min_disp=1e-8)
    if !isfinite(param_err) || !isfinite(local_err) || param_err > 1.1 * local_err
        return nothing
    end

    return coef
end

function _fit_dispersion_trend(means::Vector{Float64},
    gene_wise::Vector{Float64}; fit_type::Symbol=:parametric, min_disp::Real=1e-8)

    if !(fit_type in (:parametric, :local, :mean, :adaptive))
        fit_type = :parametric
    end

    trend_meta = Dict{String,Any}()

    valid = isfinite.(means) .& isfinite.(gene_wise) .& (means .> 0.0) .& (gene_wise .> 0.0)
    trend_meta["n_valid_genes"] = count(valid)
    if count(valid) < 5
        fillv = _safe_median(gene_wise)
        trend_meta["selection_reason"] = "insufficient_valid_genes"
        return fill(fillv, length(means)), :mean, trend_meta
    end

    # Match DESeq2 fitting subset semantics: use genes well above minimum
    # dispersion for fitting and interpolation over all means.
    fit_mask = valid .& (gene_wise .>= 10.0 * min_disp)
    if count(fit_mask) < 5
        fit_mask = valid
    end
    trend_meta["n_fit_genes"] = count(fit_mask)

    if fit_type == :mean
        fillv = median(gene_wise[fit_mask])
        trend_meta["selection_reason"] = "requested_mean"
        return fill(max(fillv, min_disp), length(means)), :mean, trend_meta
    end

    local_trend, selected_span, local_score, local_span_scores = _adaptive_local_dispersion_trend(means, gene_wise;
        fit_mask=fit_mask,
        min_disp=min_disp)
    local_trend .= max.(local_trend, min_disp)
    if any(.!isfinite.(local_trend)) || any(local_trend .<= 0.0)
        fillv = _safe_median(gene_wise[fit_mask])
        local_trend = fill(max(fillv, min_disp), length(means))
        local_score = Inf
    end
    local_score, local_mae, local_bias, local_tail = _trend_fit_score(gene_wise, local_trend, fit_mask; min_disp=min_disp)
    trend_meta["selected_local_span"] = selected_span
    trend_meta["local_span_scores"] = local_span_scores
    trend_meta["local_score"] = local_score
    trend_meta["local_median_abs_log_residual"] = local_mae
    trend_meta["local_log_bias"] = local_bias
    trend_meta["local_log_tail_q90"] = local_tail

    if fit_type == :local
        trend_meta["selection_reason"] = "requested_local"
        return local_trend, :local, trend_meta
    end

    x = means[fit_mask]
    y = gene_wise[fit_mask]

    coef = _parametric_dispersion_fit(x, y)
    if coef === nothing
        trend_meta["selection_reason"] = "parametric_fit_failed"
        return local_trend, :local, trend_meta
    end
    a = coef[1]
    b = coef[2]
    trend_meta["parametric_a"] = a
    trend_meta["parametric_b"] = b

    # DESeq2 parametric fit rejects non-positive coefficients and falls back to local.
    if !(isfinite(a) && isfinite(b) && a > 0.0 && b > 0.0)
        trend_meta["selection_reason"] = "parametric_nonpositive_coefficients"
        return local_trend, :local, trend_meta
    end

    trend = a .+ b ./ max.(means, min_disp)

    if any(.!isfinite.(trend)) || any(trend .<= 0.0)
        trend_meta["selection_reason"] = "parametric_trend_invalid"
        return local_trend, :local, trend_meta
    end

    param_score, param_mae, param_bias, param_tail = _trend_fit_score(gene_wise, trend, fit_mask; min_disp=min_disp)
    trend_meta["parametric_score"] = param_score
    trend_meta["parametric_median_abs_log_residual"] = param_mae
    trend_meta["parametric_log_bias"] = param_bias
    trend_meta["parametric_log_tail_q90"] = param_tail

    if !isfinite(param_score)
        trend_meta["selection_reason"] = "parametric_score_nonfinite"
        return local_trend, :local, trend_meta
    end

    score_margin = 0.05 + 0.75 / sqrt(max(count(fit_mask), 1))
    quality_gate = param_score <= 1.0
    use_parametric = if fit_type == :adaptive
        quality_gate && (!isfinite(local_score) || param_score <= local_score)
    else
        quality_gate && (!isfinite(local_score) || param_score <= (1.0 + score_margin) * local_score)
    end
    trend_meta["parametric_quality_gate"] = quality_gate

    if !use_parametric
        trend_meta["selection_reason"] = fit_type == :adaptive ?
                                         "adaptive_selected_local" : "parametric_underperformed_local"
        return local_trend, :local, trend_meta
    end

    trend_meta["selection_reason"] = fit_type == :adaptive ?
                                     "adaptive_selected_parametric" : "parametric_accepted"
    trend_meta["parametric_score_margin"] = score_margin
    return max.(trend, min_disp), :parametric, trend_meta
end

function _estimate_dispersion_prior_variance_lowdf(logres::Vector{Float64}, df::Int)
    obs = logres[isfinite.(logres)]
    length(obs) < 5 && return 0.25

    # Match DESeq2's bounded-residual approach for low residual DF.
    edges = collect(-10.0:0.5:10.0)
    obs = obs[(obs.>first(edges)).&(obs.<last(edges))]
    length(obs) < 5 && return 0.25

    function hist_prob(values::Vector{Float64}, bins::Vector{Float64})
        nbins = length(bins) - 1
        counts = zeros(Float64, nbins)
        for v in values
            idx = searchsortedlast(bins, v) - 1
            if 1 <= idx <= nbins
                counts[idx] += 1.0
            end
        end
        total = sum(counts)
        total <= 0.0 && return fill(1.0 / nbins, nbins)
        probs = counts ./ total
        return probs
    end

    obs_prob = hist_prob(obs, edges)
    obs_var_grid = collect(range(0.0, 8.0; length=200))
    rng = MersenneTwister(2)
    nmc = 10_000
    chi = Chisq(max(df, 1))

    best_kl = Inf
    best_var = 0.25
    for v in obs_var_grid
        rand_dist = log.(rand(rng, chi, nmc)) .+ randn(rng, nmc) .* sqrt(v) .- log(max(df, 1))
        rand_dist = rand_dist[(rand_dist.>first(edges)).&(rand_dist.<last(edges))]
        if isempty(rand_dist)
            continue
        end
        rnd_prob = hist_prob(rand_dist, edges)
        z = vcat(obs_prob, rnd_prob)
        positive = z[z.>0.0]
        small = isempty(positive) ? 1e-12 : min(minimum(positive), 1e-12)
        kl = sum(obs_prob .* (log.(obs_prob .+ small) .- log.(rnd_prob .+ small)))
        if isfinite(kl) && kl < best_kl
            best_kl = kl
            best_var = v
        end
    end

    return max(best_var, 0.25)
end

function _estimate_dispersion_prior_variance(gene_wise::Vector{Float64},
    trend::Vector{Float64}; model_df::Union{Nothing,Int}=nothing)
    logres = log.(max.(gene_wise, 1e-8)) .- log.(max.(trend, 1e-8))
    valid = logres[isfinite.(logres)]
    length(valid) < 3 && return 0.25

    df = model_df === nothing ? 4 : max(model_df, 1)
    if 1 <= df <= 3
        return _estimate_dispersion_prior_variance_lowdf(valid, df)
    end

    s2 = var(valid; corrected=true)
    if !isfinite(s2)
        med = median(valid)
        mad = median(abs.(valid .- med))
        s2 = (1.4826 * mad)^2
    end
    sampling_var = trigamma(df / 2)
    if !isfinite(sampling_var)
        sampling_var = 0.0
    end
    prior = s2 - sampling_var
    return max(prior, 0.25)
end

function _map_shrink_dispersions(gene_wise::Vector{Float64}, trend::Vector{Float64},
    prior_var::Real; model_df::Int=4, outlier_sd::Real=2.0, var_log_disp_est::Union{Nothing,Real}=nothing)
    prior_var_f = max(Float64(prior_var), 1e-8)
    sample_var = max(trigamma(max(model_df, 1) / 2), 1e-8)
    if !isfinite(sample_var)
        sample_var = 1.0
    end
    var_log_disp = var_log_disp_est === nothing ? prior_var_f + sample_var : max(Float64(var_log_disp_est), 1e-8)
    if !isfinite(var_log_disp)
        var_log_disp = prior_var_f + sample_var
    end
    out = similar(gene_wise)

    for i in eachindex(gene_wise)
        gw = max(gene_wise[i], 1e-8)
        tr = max(trend[i], 1e-8)
        lg = log(gw)
        lt = log(tr)
        if !isfinite(lg) || !isfinite(lt)
            out[i] = tr
        elseif lg > lt + outlier_sd * sqrt(var_log_disp)
            out[i] = gw
        else
            post = (lg / sample_var + lt / prior_var_f) / (1 / sample_var + 1 / prior_var_f)
            out[i] = isfinite(post) ? exp(post) : tr
        end
    end

    return max.(out, 1e-8)
end

function estimate_dispersions_prior(cm::CountMatrix,
    norm_factors::AbstractVector{<:Real};
    fit_type::Symbol=:parametric,
    model_matrix=nothing,
    model_df=nothing,
    outlier_sd::Real=2.0,
    rng=Random.default_rng(),
    huber::Bool=false,
    use_cv::Bool=true,
    apply_df_bias_correction::Bool=true)

    Y = cm.counts
    ng, ns = size(Y)
    nf = Float64.(norm_factors)
    length(nf) == ns || throw(ArgumentError("norm_factors length mismatch"))

    means = zeros(Float64, ng)
    gene_wise = zeros(Float64, ng)
    max_disp = max(10.0, Float64(ns))
    norm_counts = Matrix{Float64}(Y) ./ reshape(nf, 1, :)
    n_gene_mle_fallback = 0
    n_gene_mle_error = 0
    alpha_init = fill(NaN, ng)
    df_bias_correction_factor = 1.0
    df_bias_correction_factor_df = 1.0
    df_bias_correction_factor_empirical = 1.0

    if model_matrix === nothing
        for g in 1:ng
            row = _as_dense_row(Y, g)
            yn = row ./ nf
            means[g] = max(mean(yn), 1e-8)
            gene_wise[g] = clamp(_moment_dispersion(row, nf), 1e-8, max_disp)
        end
    else
        X = Matrix{Float64}(model_matrix)
        rough = _rough_dispersion_estimates(norm_counts, X)
        moments = _moments_dispersion_estimates(norm_counts, nf)
        alpha_init = clamp.(min.(rough, moments), 1e-8, max_disp)
        for g in 1:ng
            row = _as_dense_row(Y, g)
            yn = row ./ nf
            means[g] = max(mean(yn), 1e-8)
            alpha_g = try
                _fit_gene_dispersion_mle(row, nf, X;
                    alpha_init=alpha_init[g],
                    alpha_min=1e-8,
                    alpha_max=max_disp,
                    maxit=50,
                    beta_tol=1e-8,
                    minmu=0.5)
            catch
                n_gene_mle_error += 1
                NaN
            end
            if !isfinite(alpha_g) || alpha_g <= 0.0
                n_gene_mle_fallback += 1
                alpha_g = _moment_dispersion(row, nf)
            end
            gene_wise[g] = clamp(alpha_g, 1e-8, max_disp)
        end
    end

    # Finite-sample correction: dispersion estimators in small residual-DF
    # settings are typically downward-biased. A residual-DF correction
    # analogous to unbiased residual variance uses n/(n-p), with safeguards.
    n_bias_eligible = 0
    if apply_df_bias_correction && model_matrix !== nothing
        p = size(model_matrix, 2)
        rdf = ns - p
        min_bias_alpha = 1e-6
        if rdf >= 2
            df_bias_correction_factor_df = clamp(ns / rdf, 1.0, 2.0)
        end

        valid_bias = isfinite.(alpha_init) .& isfinite.(gene_wise) .&
                     (alpha_init .> min_bias_alpha) .& (gene_wise .> min_bias_alpha)
        if count(valid_bias) >= 20
            ratio = alpha_init[valid_bias] ./ gene_wise[valid_bias]
            ratio = ratio[isfinite.(ratio).&(ratio.>0.0)]
            if length(ratio) >= 10
                df_bias_correction_factor_empirical = clamp(median(ratio), 1.0, 2.5)
            end
        end

        df_bias_correction_factor = clamp(
            df_bias_correction_factor_df * df_bias_correction_factor_empirical^0.25,
            1.0, 2.0)

        correction_eligible = isfinite.(gene_wise) .& (gene_wise .> min_bias_alpha)
        n_bias_eligible = count(correction_eligible)

        if df_bias_correction_factor > 1.0 + 1e-8 && n_bias_eligible >= 10
            gene_wise[correction_eligible] .= clamp.(
                gene_wise[correction_eligible] .* df_bias_correction_factor,
                1e-8,
                max_disp)
        end
    end

    trend, actual_fit, trend_fit_diagnostics = _fit_dispersion_trend(means, gene_wise; fit_type=fit_type)
    if any(.!isfinite.(trend)) || any(trend .<= 0.0)
        fillv = _safe_median(gene_wise)
        trend = fill(max(fillv, 1e-8), length(means))
        actual_fit = :mean
        trend_fit_diagnostics = Dict{String,Any}("selection_reason" => "invalid_trend_forced_mean")
    end

    df = if model_df === nothing
        if model_matrix === nothing
            max(ns - 2, 1)
        else
            max(ns - size(model_matrix, 2), 1)
        end
    else
        Int(model_df)
    end

    prior_var = try
        _estimate_dispersion_prior_variance(gene_wise, trend; model_df=df)
    catch
        0.25
    end
    if !isfinite(prior_var) || prior_var <= 0.0
        prior_var = 0.25
    end

    logres = log.(max.(gene_wise, 1e-8)) .- log.(max.(trend, 1e-8))
    valid_logres = logres[isfinite.(logres)]
    var_log_disp_est = if length(valid_logres) >= 2
        var(valid_logres; corrected=true)
    else
        prior_var + max(trigamma(df / 2), 0.0)
    end
    if !isfinite(var_log_disp_est) || var_log_disp_est <= 0.0
        var_log_disp_est = prior_var + 1.0
    end

    map_disp = try
        _map_shrink_dispersions(gene_wise, trend, prior_var;
            model_df=df,
            outlier_sd=outlier_sd,
            var_log_disp_est=var_log_disp_est)
    catch
        copy(trend)
    end
    if any(.!isfinite.(map_disp)) || any(map_disp .<= 0.0)
        map_disp = copy(trend)
    end

    diagnostics = Dict{String,Any}(
        "requested_fit_type" => String(fit_type),
        "actual_fit_type" => String(actual_fit),
        "model_df" => df,
        "n_genes" => ng,
        "n_gene_mle_fallback" => n_gene_mle_fallback,
        "n_gene_mle_error" => n_gene_mle_error,
        "apply_df_bias_correction" => apply_df_bias_correction,
        "df_bias_correction_factor_df" => df_bias_correction_factor_df,
        "df_bias_correction_factor_empirical" => df_bias_correction_factor_empirical,
        "df_bias_correction_factor" => df_bias_correction_factor,
        "n_bias_eligible" => n_bias_eligible,
        "trend_fit" => trend_fit_diagnostics,
        "prior_var" => prior_var,
        "var_log_disp_est" => var_log_disp_est)

    return (
        mean_expression=means,
        gene_wise=gene_wise,
        trend=trend,
        map=map_disp,
        prior_variance=prior_var,
        fit_type=actual_fit,
        diagnostics=diagnostics)
end

function estimate_dispersions(cm::CountMatrix,
    norm_factors::AbstractVector{<:Real};
    workflow::Symbol=:prior,
    fit_type::Symbol=:parametric,
    model_matrix=nothing,
    rng=Random.default_rng(),
    huber::Bool=false)

    info = estimate_dispersions_prior(cm, norm_factors;
        fit_type=fit_type,
        model_matrix=model_matrix,
        rng=rng,
        huber=huber)

    if workflow == :prior
        return info.map
    elseif workflow == :genewise
        return info.gene_wise
    elseif workflow == :fit
        return info.trend
    elseif workflow == :simple
        return info.gene_wise
    else
        throw(ArgumentError("unsupported workflow: $workflow"))
    end
end

# =============================================================================
# NB-GLM Fitting
# =============================================================================

"""
    fit_gene_fast!(solver, y, offset, dispersion; ...)

Fit NB-GLM by IRLS. Returns (beta, se, stat) or (beta, cov, se, stat) if return_covariance=true.
"""
function fit_gene_fast!(solver::GLMSolver,
    y::AbstractVector{<:Real},
    offset::AbstractVector{<:Real},
    dispersion::Real;
    maxit::Integer=100,
    beta_tol::Real=1e-8,
    minmu::Real=0.5,
    minmu_boundary_mode::Symbol=:project,
    coef_idx::Union{Nothing,Int}=nothing,
    return_covariance::Bool=false)

    X = solver.X
    n, p = size(X)
    length(y) == n || throw(ArgumentError("response length does not match design rows"))
    length(offset) == n || throw(ArgumentError("offset length does not match design rows"))

    yv = Float64.(y)
    off = Float64.(offset)
    alpha = max(Float64(dispersion), 1e-8)
    minmu_boundary_mode in (:project, :clamped_eta) ||
        throw(ArgumentError("unsupported minmu_boundary_mode: $minmu_boundary_mode"))

    fill!(solver.beta, 0.0)
    solver.beta[1] = log(max(mean(yv ./ exp.(off)), minmu))
    solver.converged = false
    solver.iterations = 0

    old_ll = -Inf
    tol = max(Float64(beta_tol), 1e-10)

    for iter in 1:maxit
        solver.iterations = iter
        beta_prev = copy(solver.beta)

        # Compute eta from current beta (unclamped)
        solver.eta .= X * solver.beta .+ off

        @inbounds for i in 1:n
            eta_i = solver.eta[i]
            mui_raw = exp(eta_i)
            mui_clamped = max(mui_raw, minmu)  # Clamp mu for variance/weights only
            solver.mu[i] = mui_clamped
            var_i = mui_clamped + alpha * mui_clamped^2
            solver.weights[i] = mui_clamped^2 / max(var_i, eps(Float64))
            eta_eff = minmu_boundary_mode == :clamped_eta ? log(mui_clamped) : eta_i
            solver.z[i] = eta_eff + (yv[i] - mui_clamped) / mui_clamped
        end

        fill!(solver.xtwx, 0.0)
        fill!(solver.xtwz, 0.0)
        @inbounds for i in 1:n
            wi = solver.weights[i]
            zi = solver.z[i] - off[i]
            xi = @view X[i, :]
            for a in 1:p
                xa = xi[a]
                solver.xtwz[a] += wi * xa * zi
                for b in a:p
                    solver.xtwx[a, b] += wi * xa * xi[b]
                end
            end
        end
        @inbounds for a in 2:p
            for b in 1:(a-1)
                solver.xtwx[a, b] = solver.xtwx[b, a]
            end
        end

        copy!(solver.raw_xtwx, solver.xtwx)

        beta_wls = try
            (solver.xtwx + 1e-8 * I(p)) \ solver.xtwz
        catch
            pinv(solver.xtwx + 1e-6 * I(p)) * solver.xtwz
        end

        all(isfinite, beta_wls) || break

        step = 1.0
        accepted = false
        for _ in 0:solver.max_step_halvings
            beta_trial = beta_prev .+ step .* (beta_wls .- beta_prev)
            mu_trial = max.(exp.(X * beta_trial .+ off), minmu)
            ll = _nb_loglik(yv, mu_trial, alpha)
            if !isfinite(old_ll) || ll >= old_ll - 1e-10
                solver.beta .= beta_trial
                old_ll = ll
                accepted = true
                break
            end
            step *= 0.5
        end

        if !accepted
            solver.beta .= beta_wls
        end

        delta = maximum(abs.(solver.beta .- beta_prev))
        if delta < tol
            solver.converged = true
            break
        end
    end

    if minmu_boundary_mode == :project
        solver.beta .= _stabilize_beta_on_minmu_boundary(X, off, solver.beta, minmu)
    end

    solver.eta .= X * solver.beta .+ off
    solver.mu .= max.(exp.(solver.eta), minmu)
    solver.weights .= solver.mu .^ 2 ./ (solver.mu .+ alpha .* solver.mu .^ 2)

    fill!(solver.xtwx, 0.0)
    @inbounds for i in 1:n
        wi = solver.weights[i]
        xi = @view X[i, :]
        for a in 1:p
            xa = xi[a]
            for b in a:p
                solver.xtwx[a, b] += wi * xa * xi[b]
            end
        end
    end
    @inbounds for a in 2:p
        for b in 1:(a-1)
            solver.xtwx[a, b] = solver.xtwx[b, a]
        end
    end

    copy!(solver.raw_xtwx, solver.xtwx)

    cov = _safe_inverse(solver.xtwx)
    idx = coef_idx === nothing ? min(2, p) : coef_idx
    se = sqrt(max(cov[idx, idx], eps(Float64)))
    stat = abs(solver.beta[idx] / se)

    if return_covariance
        return copy(solver.beta), cov, se, stat
    end
    return copy(solver.beta), se, stat
end

# =============================================================================
# Cook's Distance and Outlier Replacement
# =============================================================================

function cooks_distance(y::AbstractVector{<:Real},
    mu::AbstractVector{<:Real},
    leverage::AbstractVector{<:Real},
    dispersion::Real; p::Int=2)
    alpha = max(Float64(dispersion), 1e-8)
    out = zeros(Float64, length(y))
    @inbounds for i in eachindex(y)
        var_i = mu[i] + alpha * mu[i]^2
        r = (y[i] - mu[i]) / sqrt(max(var_i, eps(Float64)))
        h = clamp(leverage[i], 1e-8, 1.0 - 1e-8)
        out[i] = (r^2 * h) / (p * (1.0 - h)^2)
    end
    return out
end

function _trimmed_mean_fraction(values::AbstractVector{<:Real}, trim::Real)
    isempty(values) && return NaN
    x = sort(Float64.(values))
    n = length(x)
    k = clamp(floor(Int, trim * n), 0, max(n - 1, 0))
    if 2k >= n
        return mean(x)
    end
    return mean(@view x[(k+1):(n-k)])
end

function _trimmed_variance_row(x::AbstractVector{<:Real})
    rm = _trimmed_mean_fraction(x, 1 / 8)
    sqerror = (Float64.(x) .- rm) .^ 2
    return 1.51 * _trimmed_mean_fraction(sqerror, 1 / 8)
end

function _n_or_more_in_cell(model_matrix::AbstractMatrix{<:Real}, n::Int=3)
    row_hash = [join(model_matrix[i, :], "_") for i in axes(model_matrix, 1)]
    counts = Dict{String,Int}()
    for key in row_hash
        counts[key] = get(counts, key, 0) + 1
    end
    return [get(counts, key, 0) >= n for key in row_hash]
end

function _design_cell_groups(model_matrix::AbstractMatrix{<:Real})
    groups = Dict{String,Vector{Int}}()
    for i in axes(model_matrix, 1)
        key = join(model_matrix[i, :], "_")
        push!(get!(groups, key, Int[]), i)
    end
    return collect(values(groups))
end

function _has_complete_zero_cell(y::AbstractVector{<:Real}, cell_groups::Vector{Vector{Int}})
    has_all_zero_cell = false
    has_nonzero_cell = false
    for idx in cell_groups
        cell_has_nonzero = false
        for i in idx
            if y[i] != 0
                cell_has_nonzero = true
                break
            end
        end
        if cell_has_nonzero
            has_nonzero_cell = true
        else
            has_all_zero_cell = true
        end
        has_all_zero_cell && has_nonzero_cell && return true
    end
    return false
end

function _trimmed_cell_variance_row(cnts::AbstractVector{<:Real}, cells)
    trimratio = (1 / 3, 1 / 4, 1 / 8)
    scalec = (2.04, 1.86, 1.51)
    trimbin(n) = n <= 3 ? 1 : n <= 23 ? 2 : 3
    maxvar = 0.0
    for lvl in unique(cells)
        idx = findall(==(lvl), cells)
        n = length(idx)
        b = trimbin(n)
        t = trimratio[b]
        mu = _trimmed_mean_fraction(cnts[idx], t)
        sqerror = (Float64.(cnts[idx]) .- mu) .^ 2
        var_est = scalec[b] * _trimmed_mean_fraction(sqerror, t)
        maxvar = max(maxvar, var_est)
    end
    return maxvar
end

function _robust_cooks_dispersion(y::AbstractVector{<:Real},
    nf::AbstractVector{<:Real},
    model_matrix::AbstractMatrix{<:Real})
    norm_counts = Float64.(y) ./ Float64.(nf)
    m = mean(norm_counts)
    if !isfinite(m) || m <= 0.0
        return 0.04
    end
    three_or_more = _n_or_more_in_cell(model_matrix, 3)
    v = if any(three_or_more)
        cells = [join(model_matrix[i, :], "_") for i in axes(model_matrix, 1)]
        eligible_levels = Set(cells[three_or_more])
        idx = findall(in(eligible_levels), cells)
        _trimmed_cell_variance_row(norm_counts[idx], cells[idx])
    else
        _trimmed_variance_row(norm_counts)
    end
    alpha = (v - m) / max(m^2, eps(Float64))
    return max(alpha, 0.04)
end

function _trimmed_mean(values::AbstractVector{<:Real}, trim::Real)
    isempty(values) && return NaN
    x = sort(Float64.(values))
    n = length(x)
    k = clamp(floor(Int, trim * n), 0, max(n - 1, 0))
    if 2k >= n
        return mean(x)
    end
    return mean(@view x[(k+1):(n-k)])
end

function replace_outliers(cm::CountMatrix, design_vec;
    norm_factors=nothing,
    dispersions=nothing,
    cooks_cutoff=nothing,
    min_replicates::Int=7,
    trim::Real=0.2,
    maxit::Integer=100,
    beta_tol::Real=1e-8,
    minmu::Real=0.5)

    nf = norm_factors === nothing ? calc_norm_factors(cm; method=:tmm) : Float64.(norm_factors)
    disp = dispersions === nothing ? estimate_dispersions(cm, nf; workflow=:prior) : Float64.(dispersions)

    X, _, _ = _design_model_matrix(design_vec)
    ns = size(X, 1)
    p = size(X, 2)
    ng = length(cm.gene_ids)

    cutoff = cooks_cutoff === nothing ? quantile(FDist(p, max(ns - p, 1)), 0.99) : Float64(cooks_cutoff)

    Y = Matrix{Int}(cm.counts)
    replaced = falses(ng, ns)
    cooks = zeros(Float64, ng, ns)

    labels = String.(design_vec)
    solver = GLMSolver(X)
    offset = log.(nf)

    for g in 1:ng
        y = Float64.(Y[g, :])
        alpha = disp[g]

        if sum(y) <= 0
            continue
        end

        beta, cov, _, _ = fit_gene_fast!(solver, y, offset, alpha;
            maxit=maxit, beta_tol=beta_tol, minmu=minmu, return_covariance=true)

        h = zeros(Float64, ns)
        @inbounds for s in 1:ns
            xs = @view X[s, :]
            h[s] = clamp(solver.weights[s] * dot(xs, cov * xs), 1e-8, 1.0 - 1e-8)
        end

        c = cooks_distance(y, solver.mu, h, alpha; p=p)
        cooks[g, :] .= c

        outlier_idx = findall(ci -> ci > cutoff && isfinite(ci), c)
        isempty(outlier_idx) && continue

        for s in outlier_idx
            grp = labels[s]
            grp_idx = findall(labels .== grp)
            if length(grp_idx) < min_replicates
                continue
            end
            keep = [j for j in grp_idx if j != s && c[j] <= cutoff]
            isempty(keep) && continue
            repl_norm = _trimmed_mean(y[keep] ./ nf[keep], trim)
            isfinite(repl_norm) || continue
            Y[g, s] = max(0, round(Int, repl_norm * nf[s]))
            replaced[g, s] = true
        end
    end

    return (
        counts=CountMatrix(sparse(Y), cm.gene_ids, cm.sample_ids),
        cooks=cooks,
        replaced=replaced,
        cutoff=cutoff)
end

# =============================================================================
# DESeqDataSet API Methods
# =============================================================================

function _effective_norm_factors(dds::DESeqDataSet)
    dds.normalization_factors !== nothing && return dds.normalization_factors
    dds.size_factors !== nothing && return dds.size_factors
    return calc_norm_factors(dds.counts; method=dds.sf_type)
end

function sizeFactors!(dds::DESeqDataSet, values::AbstractVector{<:Real})
    length(values) == length(dds.counts.sample_ids) ||
        throw(ArgumentError("sizeFactors length mismatch"))
    dds.size_factors = Float64.(values)
    dds.metadata["sizeFactors"] = copy(dds.size_factors)
    update_provenance!(dds.metadata; source="DifferentialExpression/sizeFactors!", notes=["assigned explicit size factors"], parameters=(n_samples=length(dds.size_factors)))

    return dds
end

function normalizationFactors!(dds::DESeqDataSet, values::AbstractVector{<:Real})
    length(values) == length(dds.counts.sample_ids) ||
        throw(ArgumentError("normalizationFactors length mismatch"))
    dds.normalization_factors = Float64.(values)
    dds.metadata["normalizationFactors"] = copy(dds.normalization_factors)
    update_provenance!(dds.metadata; source="DifferentialExpression/normalizationFactors!", notes=["assigned explicit normalization factors"], parameters=(n_samples=length(dds.normalization_factors)))

    return dds
end

function estimateSizeFactors(dds::DESeqDataSet;
    type::Symbol=:ratio,
    sfType=nothing,
    geoMeans=nothing,
    controlGenes=nothing,
    fallback::Symbol=:library_size,
    force::Bool=false)

    if !force && dds.size_factors !== nothing
        return dds
    end

    chosen = sfType === nothing ? type : sfType
    if chosen == :tmm
        dds.size_factors = calc_norm_factors(dds.counts; method=:tmm)
    elseif chosen == :ratio || chosen == :poscounts
        dds.size_factors = estimateSizeFactorsForMatrix(dds.counts.counts;
            type=chosen,
            geoMeans=geoMeans,
            controlGenes=controlGenes,
            fallback=fallback)
    else
        throw(ArgumentError("unsupported size-factor method: $chosen"))
    end

    dds.sf_type = chosen
    dds.metadata["sizeFactors"] = copy(dds.size_factors)
    update_provenance!(dds.metadata; source="DifferentialExpression/estimateSizeFactors", fallbacks=fallback == :library_size ? ["library-size fallback enabled for invalid ratio factors"] : String[], parameters=(method=chosen, fallback=fallback, force=Bool(force)))
    return dds
end

function estimateDispersionsGeneEst(dds::DESeqDataSet, kwargs...)
    nf = _effective_norm_factors(dds)
    X = dds.model_matrix === nothing ? _design_model_matrix(dds.design)[1] : dds.model_matrix
    wc = dds.replaced_counts === nothing ? dds.counts : dds.replaced_counts
    info = estimate_dispersions_prior(wc, nf; fit_type=dds.fit_type, model_matrix=X)
    dds.gene_wise_dispersions = info.gene_wise
    dds.dispersion_fit = info.trend
    dds.dispersion_prior_variance = info.prior_variance
    dds.dispersions = info.gene_wise
    dds.actual_fit_type = info.fit_type
    dds.metadata["dispersionDiagnostics"] = get(info, :diagnostics, Dict{String,Any}())
    update_provenance!(dds.metadata; source="DifferentialExpression/estimateDispersionsGeneEst", parameters=(fit_type=dds.fit_type, actual_fit_type=dds.actual_fit_type))

    return dds
end

function estimateDispersionsFit(dds::DESeqDataSet; fitType::Symbol=dds.fit_type, kwargs...)
    nf = _effective_norm_factors(dds)
    X = dds.model_matrix === nothing ? _design_model_matrix(dds.design)[1] : dds.model_matrix
    wc = dds.replaced_counts === nothing ? dds.counts : dds.replaced_counts
    info = estimate_dispersions_prior(wc, nf; fit_type=fitType, model_matrix=X)
    dds.gene_wise_dispersions = info.gene_wise
    dds.dispersion_fit = info.trend
    dds.dispersion_prior_variance = info.prior_variance
    dds.dispersions = info.trend
    dds.fit_type = fitType
    dds.actual_fit_type = info.fit_type
    dds.metadata["dispersionDiagnostics"] = get(info, :diagnostics, Dict{String,Any}())
    update_provenance!(dds.metadata; source="DifferentialExpression/estimateDispersionsFit", parameters=(fit_type=fitType, actual_fit_type=dds.actual_fit_type))

    return dds
end

function estimateDispersionsMAP(dds::DESeqDataSet;
    fitType::Symbol=dds.fit_type,
    outlierSD::Real=2.0,
    dispPriorVar=nothing,
    model_df=nothing,
    rng=Random.default_rng(),
    huber::Bool=false,
    use_cv::Bool=true)

    nf = _effective_norm_factors(dds)
    X = dds.model_matrix === nothing ? _design_model_matrix(dds.design)[1] : dds.model_matrix
    wc = dds.replaced_counts === nothing ? dds.counts : dds.replaced_counts
    info = estimate_dispersions_prior(wc, nf;
        fit_type=fitType,
        model_matrix=X,
        model_df=model_df,
        outlier_sd=outlierSD,
        rng=rng,
        huber=huber,
        use_cv=use_cv)

    prior_var = dispPriorVar === nothing ? info.prior_variance : Float64(dispPriorVar)
    map_disp = dispPriorVar === nothing ? info.map :
               _map_shrink_dispersions(info.gene_wise, info.trend, prior_var;
        model_df=model_df === nothing ? max(size(X, 1) - size(X, 2), 1) : Int(model_df),
        outlier_sd=outlierSD)

    dds.gene_wise_dispersions = info.gene_wise
    dds.dispersion_fit = info.trend
    dds.dispersion_prior_variance = prior_var
    dds.dispersions = map_disp
    dds.fit_type = fitType
    dds.actual_fit_type = info.fit_type
    dds.metadata["dispersionPriorVar"] = prior_var
    dds.metadata["dispersionDiagnostics"] = get(info, :diagnostics, Dict{String,Any}())
    update_provenance!(dds.metadata; source="DifferentialExpression/estimateDispersionsMAP", parameters=(fit_type=fitType, actual_fit_type=dds.actual_fit_type, outlierSD=Float64(outlierSD), huber=Bool(huber), use_cv=Bool(use_cv)))
    return dds
end

function estimateDispersionsPriorVar(dds::DESeqDataSet; model_df=nothing, rng=Random.default_rng())
    dds.dispersion_prior_variance !== nothing && return dds.dispersion_prior_variance
    nf = _effective_norm_factors(dds)
    X = dds.model_matrix === nothing ? _design_model_matrix(dds.design)[1] : dds.model_matrix
    wc = dds.replaced_counts === nothing ? dds.counts : dds.replaced_counts
    info = estimate_dispersions_prior(wc, nf; model_matrix=X, model_df=model_df)
    dds.dispersion_prior_variance = info.prior_variance
    update_provenance!(dds.metadata; source="DifferentialExpression/estimateDispersionsPriorVar", parameters=(prior_variance=info.prior_variance))
    return info.prior_variance
end

function estimateDispersions(dds::DESeqDataSet;
    workflow::Symbol=:prior,
    fitType::Symbol=dds.fit_type,
    model_df=nothing,
    outlierSD::Real=2.0,
    rng=Random.default_rng(),
    huber::Bool=false,
    use_cv::Bool=true)

    if workflow == :prior
        return estimateDispersionsMAP(dds;
            fitType=fitType,
            outlierSD=outlierSD,
            model_df=model_df,
            rng=rng,
            huber=huber,
            use_cv=use_cv)
    elseif workflow == :genewise
        return estimateDispersionsGeneEst(dds)
    elseif workflow == :fit
        return estimateDispersionsFit(dds; fitType=fitType)
    elseif workflow == :simple
        nf = _effective_norm_factors(dds)
        X = dds.model_matrix === nothing ? _design_model_matrix(dds.design)[1] : dds.model_matrix
        wc = dds.replaced_counts === nothing ? dds.counts : dds.replaced_counts
        dds.dispersions = estimate_dispersions(wc, nf;
            workflow=:simple,
            fit_type=fitType,
            model_matrix=X)
        return dds
    else
        throw(ArgumentError("unsupported workflow: $workflow"))
    end
end

@inline _beta_prior_var_log2_to_natural(prior_var::Real) =
    max(Float64(prior_var) * log(2)^2, eps(Float64))

@inline _beta_prior_var_natural_to_log2(prior_var::Real) =
    max(Float64(prior_var) / log(2)^2, eps(Float64))

function estimate_beta_prior_var(beta_matrix::AbstractMatrix{<:Real},
    covariances::AbstractVector{<:AbstractMatrix{<:Real}}; coef_idx::Int=2)
    ncoef = size(beta_matrix, 2)
    ncoef == 0 && return 1.0
    j = clamp(coef_idx, 1, ncoef)
    b = Float64.(beta_matrix[:, j])
    beta_valid = Float64[]
    se2 = Float64[]

    for i in eachindex(covariances)
        cov = covariances[i]
        if size(cov, 1) >= j && size(cov, 2) >= j
            se2_ij = Float64(cov[j, j])
            beta_ij = b[i]
            if isfinite(beta_ij) && isfinite(se2_ij) && se2_ij > 1e-10
                push!(beta_valid, beta_ij)
                push!(se2, se2_ij)
            end
        end
    end
    isempty(beta_valid) && return 1.0
    vb = length(beta_valid) > 1 ? var(beta_valid; corrected=true) : 0.0
    vse = isempty(se2) ? 0.0 : mean(se2)
    prior = vb - vse
    if !isfinite(prior) || prior <= 0.0
        prior = max(vb, 1e-3)
    end
    # Return in natural-log scale (beta values are on natural-log scale)
    return prior
end

function _apply_beta_prior_results(beta_matrix::Matrix{Float64},
    covariance_matrices::Vector{Matrix{Float64}},
    fisher_info_matrices::Vector{Matrix{Float64}},
    results::Vector{DEResult},
    coefficient_names::Vector{String},
    prior_var::Real;  # Now expected in NATURAL-LOG scale
    design=nothing)

    prior_var_f = Float64(prior_var)
    if !isfinite(prior_var_f) || prior_var_f <= 0.0
        return beta_matrix, covariance_matrices, results
    end

    ng, p = size(beta_matrix)
    prior_precision = zeros(Float64, p, p)
    shrink_idx = [i for i in 1:p if isempty(coefficient_names) || i > length(coefficient_names) || coefficient_names[i] != "Intercept"]
    for j in shrink_idx
        prior_precision[j, j] = 1.0 / prior_var_f
    end

    contrast = _default_contrast(coefficient_names, design)[1]
    shrunk_beta = zeros(Float64, ng, p)
    shrunk_covs = Vector{Matrix{Float64}}(undef, ng)
    shrunk_results = Vector{DEResult}(undef, ng)

    for g in 1:ng
        beta = vec(beta_matrix[g, :])
        cov = covariance_matrices[g]
        fisher = fisher_info_matrices[g]
        r = results[g]
        if !all(isfinite, beta) || !all(isfinite, cov)
            shrunk_beta[g, :] .= beta
            shrunk_covs[g] = cov
            shrunk_results[g] = r
            continue
        end

        # post_precision = fisher + prior_precision
        # post_cov = inv(post_precision)
        # post_beta = post_cov * (fisher * beta)
        system = Matrix{Float64}(I, p, p) + cov * prior_precision
        local post_beta, post_cov
        try
            post_beta = system \ beta
            post_cov = system \ cov
        catch
            # Fallback: use raw Fisher if available
            precision = _safe_inverse(cov)
            post_precision = precision + prior_precision
            post_cov = _safe_inverse(post_precision)
            fisher_beta = try
                fisher * beta
            catch
                precision * beta
            end
            post_beta = post_cov * fisher_beta
        end
        post_cov = Matrix(Symmetric((post_cov + post_cov') ./ 2.0))

        shrunk_beta[g, :] .= post_beta
        shrunk_covs[g] = post_cov
        lfc, lfc_se, stat, pvalue = _contrast_statistics(post_beta, post_cov, contrast)
        shrunk_results[g] = DEResult(r.gene_id, r.base_mean, lfc, lfc_se, stat, pvalue, NaN, r.zero_inflated, r.converged)
    end

    return shrunk_beta, shrunk_covs, shrunk_results
end

function _apply_beta_prior_results(beta_matrix::Matrix{Float64},
    covariance_matrices::Vector{Matrix{Float64}},
    results::Vector{DEResult},
    coefficient_names::Vector{String},
    prior_var::Real;
    design=nothing)

    fisher_info_matrices = Vector{Matrix{Float64}}(undef, length(covariance_matrices))
    for i in eachindex(covariance_matrices)
        fisher_info_matrices[i] = _safe_inverse(covariance_matrices[i])
    end

    return _apply_beta_prior_results(
        beta_matrix,
        covariance_matrices,
        fisher_info_matrices,
        results,
        coefficient_names,
        prior_var;
        design=design)
end

function _fit_dds_wald(dds::DESeqDataSet;
    maxit::Integer=100,
    beta_tol::Real=1e-8,
    minmu::Real=0.5,
    reference_level=nothing,
    target_level=nothing,
    modelMatrixType::Symbol=:standard)

    wc = dds.replaced_counts === nothing ? dds.counts : dds.replaced_counts
    nf = _effective_norm_factors(dds)
    X, cnames, _ = _design_model_matrix(dds.design; modelMatrixType=modelMatrixType)
    ncoef = size(X, 2)
    ng = length(wc.gene_ids)
    ns = length(wc.sample_ids)

    disp = dds.dispersions === nothing ? estimate_dispersions(wc, nf; workflow=:prior, model_matrix=X) : dds.dispersions

    beta_mat = zeros(Float64, ng, ncoef)
    covs = Vector{Matrix{Float64}}(undef, ng)
    fishers = Vector{Matrix{Float64}}(undef, ng)
    fitted_mu = zeros(Float64, ng, ns)
    hat_diag = zeros(Float64, ng, ns)
    cooks_mat = zeros(Float64, ng, ns)
    res = Vector{DEResult}(undef, ng)

    contrast, contrast_label = if reference_level !== nothing && target_level !== nothing
        _named_contrast_vector(cnames, dds.design, [:condition, target_level, reference_level])
    else
        _default_contrast(cnames, dds.design)
    end

    offset = log.(nf)
    solver = GLMSolver(X)
    solver_alt = GLMSolver(X)
    cooks_disp = zeros(Float64, ng)
    cell_groups = _design_cell_groups(X)

    for g in 1:ng
        y = _as_dense_row(wc.counts, g)
        base_mean = mean(y ./ nf)
        cooks_disp[g] = _robust_cooks_dispersion(y, nf, X)

        if all(iszero, y)
            beta_mat[g, :] .= 0.0
            covs[g] = zeros(Float64, ncoef, ncoef)
            fishers[g] = zeros(Float64, ncoef, ncoef)  # Zero Fisher for zero-count genes
            fitted_mu[g, :] .= minmu
            res[g] = DEResult(wc.gene_ids[g], 0.0, NaN, NaN, NaN, NaN, NaN, false, true)
            continue
        end

        local beta, cov
        converged = true
        lfc_override = nothing
        complete_zero_cell = _has_complete_zero_cell(y, cell_groups)
        try
            beta, cov, _, _ = fit_gene_fast!(solver, y, offset, disp[g];
                maxit=maxit,
                beta_tol=beta_tol,
                minmu=minmu,
                minmu_boundary_mode=:project,
                return_covariance=true)
            converged = solver.converged
        catch
            beta = fill(NaN, ncoef)
            cov = fill(NaN, ncoef, ncoef)
            converged = false
        end

        # Complete-separation corner case: only switch to clamped-eta mode when
        # it materially changes the contrast estimate. Keep primary Wald
        # inference from the projected fit, but report the boundary-consistent
        # LFC to remove fixed log-base offsets on separated genes.
        if complete_zero_cell && all(isfinite, beta) && all(isfinite, cov)
            try
                beta_alt, cov_alt, _, _ = fit_gene_fast!(solver_alt, y, offset, disp[g];
                    maxit=maxit,
                    beta_tol=beta_tol,
                    minmu=minmu,
                    minmu_boundary_mode=:clamped_eta,
                    return_covariance=true)
                if solver_alt.converged && all(isfinite, beta_alt) && all(isfinite, cov_alt)
                    lfc_proj, _, _, _ = _contrast_statistics(beta, cov, contrast)
                    lfc_alt, _, _, _ = _contrast_statistics(beta_alt, cov_alt, contrast)
                    if isfinite(lfc_proj) && isfinite(lfc_alt) && abs(lfc_alt - lfc_proj) > 0.5
                        lfc_override = lfc_alt
                    end
                end
            catch
            end
        end

        beta_mat[g, :] .= beta
        covs[g] = cov
        fishers[g] = copy(solver.raw_xtwx)

        if !all(isfinite, beta) || !all(isfinite, cov)
            fitted_mu[g, :] .= NaN
            hat_diag[g, :] .= NaN
            cooks_mat[g, :] .= NaN
            res[g] = DEResult(wc.gene_ids[g], base_mean, NaN, NaN, NaN, NaN, NaN, false, false)
            continue
        end

        fitted_mu[g, :] .= solver.mu
        for s in 1:ns
            xs = @view X[s, :]
            hat_diag[g, s] = clamp(solver.weights[s] * dot(xs, cov * xs), 1e-8, 1.0 - 1e-8)
        end

        @inbounds for s in 1:ns
            mu_s = solver.mu[s]
            var_s = mu_s + cooks_disp[g] * mu_s^2
            pearson_sq = (y[s] - mu_s)^2 / max(var_s, eps(Float64))
            h = hat_diag[g, s]
            cooks_mat[g, s] = (pearson_sq / ncoef) * h / (1.0 - h)^2
        end

        lfc, lfc_se, stat, pvalue = _contrast_statistics(beta, cov, contrast)
        lfc_override !== nothing && (lfc = lfc_override)
        res[g] = DEResult(wc.gene_ids[g], base_mean, lfc, lfc_se, stat, pvalue, NaN, false, converged)
    end

    return (
        results=res,
        beta_matrix=beta_mat,
        covariance_matrices=covs,
        fisher_info_matrices=fishers,
        coefficient_names=cnames,
        model_matrix=X,
        default_contrast=contrast,
        default_contrast_label=contrast_label,
        cooks=cooks_mat)
end

function nbinomWaldTest(dds::DESeqDataSet;
    betaPrior::Bool=dds.beta_prior,
    betaPriorVar=nothing,
    modelMatrixType::Symbol=:standard,
    maxit::Integer=100,
    betaTol::Real=1e-8,
    minmu::Real=0.5,
    minReplicatesForReplace::Integer=7,
    cooksCutoff=nothing,
    useOutlierReplacement::Bool=false,
    trim::Real=0.2,
    reference_level=nothing,
    target_level=nothing)

    if useOutlierReplacement
        rep = replace_outliers(dds.counts, dds.design;
            norm_factors=_effective_norm_factors(dds),
            dispersions=dds.dispersions,
            cooks_cutoff=cooksCutoff,
            min_replicates=minReplicatesForReplace,
            trim=trim,
            maxit=maxit,
            beta_tol=betaTol,
            minmu=minmu)
        dds.replaced_counts = rep.counts
        dds.cooks = rep.cooks
        dds.metadata["nReplaced"] = sum(rep.replaced)
        dds.metadata["cooksCutoff"] = rep.cutoff
    end

    fs = _fit_dds_wald(dds;
        maxit=maxit,
        beta_tol=betaTol,
        minmu=minmu,
        reference_level=reference_level,
        target_level=target_level,
        modelMatrixType=modelMatrixType)

    beta_matrix = fs.beta_matrix
    covariance_matrices = fs.covariance_matrices
    fisher_info_matrices = fs.fisher_info_matrices
    results_vec = fs.results

    prior_var_natural = if betaPriorVar === nothing
        estimate_beta_prior_var(beta_matrix, covariance_matrices)
    else
        # User provides in log2 scale; convert to natural-log scale
        _beta_prior_var_log2_to_natural(betaPriorVar)
    end

    if betaPrior
        beta_matrix, covariance_matrices, results_vec = _apply_beta_prior_results(
            beta_matrix,
            covariance_matrices,
            fisher_info_matrices,
            results_vec,
            fs.coefficient_names,
            prior_var_natural;
            design=dds.design)
    end

    dds.wald_results = results_vec
    dds.beta_prior = betaPrior
    dds.beta_prior_var = _beta_prior_var_natural_to_log2(prior_var_natural)
    dds.test = :Wald
    dds.model_matrix = fs.model_matrix
    dds.coefficient_names = fs.coefficient_names
    dds.wald_beta_matrix = beta_matrix
    dds.wald_covariances = covariance_matrices
    dds.wald_fisher_info = fisher_info_matrices
    dds.cooks = fs.cooks
    dds.metadata["defaultContrast"] = fs.default_contrast_label
    dds.metadata["betaPriorVarNatural"] = prior_var_natural
    dds.metadata["betaPriorVarScale"] = "log2"
    dds.metadata["betaPriorVarUnits"] = "variance"
    dds.metadata["betaPriorVarUsage"] = "metadata_only"
    dds.metadata["betaPriorVarNaturalUsed"] = prior_var_natural
    n_nonconv = count(!r.converged for r in dds.wald_results)
    dds.metadata["nNonConverged"] = n_nonconv
    dds.metadata["nonConvergedFraction"] = isempty(dds.wald_results) ? 0.0 : n_nonconv / length(dds.wald_results)
    update_provenance!(dds.metadata; source="DifferentialExpression/nbinomWaldTest", warnings=n_nonconv > 0 ? ["some genes did not converge"] : String[], fallbacks=useOutlierReplacement ? ["outlier replacement path used"] : String[], parameters=(betaPrior=Bool(betaPrior), modelMatrixType=modelMatrixType, maxit=Int(maxit), useOutlierReplacement=Bool(useOutlierReplacement), n_nonconverged=n_nonconv))
    n_nonconv > 0 && @warn "Wald fit: $(n_nonconv) of $(length(dds.wald_results)) genes did not converge"
    return dds
end

function nbinomLRT(dds::DESeqDataSet;
    modelMatrixType::Symbol=:standard,
    reduced=nothing,
    maxit::Integer=100,
    betaTol::Real=1e-8,
    minmu::Real=0.5,
    minReplicatesForReplace::Integer=7,
    cooksCutoff=nothing,
    useOutlierReplacement::Bool=false,
    trim::Real=0.2,
    reference_level=nothing,
    target_level=nothing)

    if useOutlierReplacement
        rep = replace_outliers(dds.counts, dds.design;
            norm_factors=_effective_norm_factors(dds),
            dispersions=dds.dispersions,
            cooks_cutoff=cooksCutoff,
            min_replicates=minReplicatesForReplace,
            trim=trim,
            maxit=maxit,
            beta_tol=betaTol,
            minmu=minmu)
        dds.replaced_counts = rep.counts
        dds.cooks = rep.cooks
        dds.metadata["nReplaced"] = sum(rep.replaced)
        dds.metadata["cooksCutoff"] = rep.cutoff
    end

    wc = dds.replaced_counts === nothing ? dds.counts : dds.replaced_counts
    nf = _effective_norm_factors(dds)
    full_X, cnames, _ = _design_model_matrix(dds.design; modelMatrixType=modelMatrixType)

    red_X = if reduced === nothing
        ones(Float64, size(full_X, 1), 1)
    elseif reduced isa AbstractMatrix{<:Real}
        _ensure_intercept(Matrix{Float64}(reduced))
    elseif reduced isa DataFrame
        _design_model_matrix(reduced)[1]
    elseif reduced isa AbstractVector
        _design_model_matrix(reduced)[1]
    else
        throw(ArgumentError("unsupported reduced model type $(typeof(reduced))"))
    end

    rank(red_X) < rank(full_X) || throw(ArgumentError("reduced model must have fewer parameters than full model"))

    disp = dds.dispersions === nothing ? estimate_dispersions(wc, nf; workflow=:prior, model_matrix=full_X) : dds.dispersions

    ng = length(wc.gene_ids)
    ns = length(wc.sample_ids)
    res = Vector{DEResult}(undef, ng)

    solver_full = GLMSolver(full_X)
    solver_red = GLMSolver(red_X)
    offset = log.(nf)

    contrast = _default_contrast(cnames)[1]
    cooks_mat = zeros(Float64, ng, ns)

    for g in 1:ng
        y = _as_dense_row(wc.counts, g)
        base_mean = mean(y ./ nf)
        if all(iszero, y)
            res[g] = DEResult(wc.gene_ids[g], 0.0, NaN, NaN, NaN, NaN, NaN, false, true)
            continue
        end

        alpha = max(disp[g], 1e-8)
        ok_full = true
        ok_red = true
        beta_full = zeros(Float64, size(full_X, 2))
        cov_full = zeros(Float64, size(full_X, 2), size(full_X, 2))

        try
            beta_full, cov_full, _, _ = fit_gene_fast!(solver_full, y, offset, alpha;
                maxit=maxit, beta_tol=betaTol, minmu=minmu, return_covariance=true)
            ok_full = solver_full.converged
        catch
            ok_full = false
        end

        try
            fit_gene_fast!(solver_red, y, offset, alpha;
                maxit=maxit, beta_tol=betaTol, minmu=minmu)
            ok_red = solver_red.converged
        catch
            ok_red = false
        end

        if !ok_full || !ok_red
            cooks_mat[g, :] .= NaN
            res[g] = DEResult(wc.gene_ids[g], base_mean, NaN, NaN, NaN, NaN, NaN, false, false)
            continue
        end

        h = zeros(Float64, ns)
        for s in 1:ns
            xs = @view full_X[s, :]
            h[s] = clamp(solver_full.weights[s] * dot(xs, cov_full * xs), 1e-8, 1.0 - 1e-8)
        end
        cooks_mat[g, :] .= cooks_distance(y, solver_full.mu, h, alpha; p=size(full_X, 2))

        ll_full = _nb_loglik(y, solver_full.mu, alpha)
        ll_red = _nb_loglik(y, solver_red.mu, alpha)
        stat = 2.0 * max(ll_full - ll_red, 0.0)
        df = rank(full_X) - rank(red_X)
        pvalue = ccdf(Chisq(df), stat)

        lfc, lfc_se, _, _ = _contrast_statistics(beta_full, cov_full, contrast)
        res[g] = DEResult(wc.gene_ids[g], base_mean, lfc, lfc_se, stat, pvalue, NaN, false, true)
    end

    dds.lrt_results = res
    dds.test = :LRT
    dds.model_matrix = full_X
    dds.reduced_model_matrix = red_X
    dds.coefficient_names = cnames
    dds.cooks = cooks_mat
    n_nonconv = count(!r.converged for r in dds.lrt_results)
    dds.metadata["nNonConverged"] = n_nonconv
    dds.metadata["nonConvergedFraction"] = isempty(dds.lrt_results) ? 0.0 : n_nonconv / length(dds.lrt_results)
    update_provenance!(dds.metadata; source="DifferentialExpression/nbinomLRT", warnings=n_nonconv > 0 ? ["some genes did not converge"] : String[], fallbacks=useOutlierReplacement ? ["outlier replacement path used"] : String[], parameters=(modelMatrixType=modelMatrixType, maxit=Int(maxit), useOutlierReplacement=Bool(useOutlierReplacement), n_nonconverged=n_nonconv))
    n_nonconv > 0 && @warn "LRT fit: $(n_nonconv) of $(length(dds.lrt_results)) genes did not converge"
    return dds
end

# =============================================================================
# Results Handling
# =============================================================================

function _canonicalize_numeric_column!(tab::DataFrame, col::Symbol; clamp01::Bool=false)
    col in propertynames(tab) || return tab
    src = tab[!, col]
    out = Vector{Union{Missing,Float64}}(undef, length(src))
    for i in eachindex(src)
        v = src[i]
        if ismissing(v)
            out[i] = missing
        else
            fv = Float64(v)
            if !isfinite(fv)
                out[i] = missing
            elseif clamp01
                out[i] = clamp(fv, 0.0, 1.0)
            else
                out[i] = fv
            end
        end
    end
    tab[!, col] = out
    return tab
end

function _canonicalize_results!(tab::DataFrame)
    _canonicalize_numeric_column!(tab, :log2_fold_change)
    _canonicalize_numeric_column!(tab, :lfc_se)
    _canonicalize_numeric_column!(tab, :stat)
    _canonicalize_numeric_column!(tab, :pvalue; clamp01=true)
    _canonicalize_numeric_column!(tab, :padj; clamp01=true)
    return tab
end

function _adjust_padj(pvalues::AbstractVector, filter_stat::AbstractVector, fdr_method::Symbol, alpha::Real)
    padj = if fdr_method == :storey
        storey_qvalue(pvalues)
    elseif fdr_method == :ihw
        ihw_qvalue(pvalues, filter_stat; alpha=alpha)
    else
        benjamini_hochberg(pvalues)
    end

    out = Vector{Union{Missing,Float64}}(undef, length(pvalues))
    for i in eachindex(pvalues)
        if ismissing(pvalues[i]) || !isfinite(Float64(pvalues[i]))
            out[i] = missing
        else
            out[i] = clamp(Float64(padj[i]), 0.0, 1.0)
        end
    end
    return out
end

# DESeq2-style independent filtering:
# choose threshold on baseMean using filtered-p style rejection curves
# (Bourgon et al. 2010 + DESeq2 results() behavior).
function _independent_filter_mask(base_mean::Vector{Float64},
    pvalue::Vector{Float64}, alpha::Real)
    finite_bm = isfinite.(base_mean)
    # DESeq2 derives candidate cutoffs from positive baseMean values.
    valid_bm = base_mean[finite_bm.&(base_mean.>0.0)]
    isempty(valid_bm) && return trues(length(base_mean)), NaN

    finite_p = isfinite.(pvalue)
    count(finite_p) == 0 && return trues(length(base_mean)), NaN

    # Match DESeq2's filtered-p approach: evaluate many filter quantiles,
    # smooth the rejection curve, then pick an early near-optimal cutoff.
    theta = collect(range(0.0, 0.95; length=50))
    rejections = fill(-1.0, length(theta))
    tested = false

    for i in eachindex(theta)
        thr = quantile(valid_bm, theta[i])
        # Match DESeq2's strict inequality for filter inclusion.
        keep = base_mean .> thr
        eval_mask = keep .& finite_p
        n_eval = count(eval_mask)
        n_eval == 0 && continue
        tested = true

        pv_subset = pvalue[eval_mask]
        # DESeq2 computes adjusted p-values on retained hypotheses after
        # filtering by baseMean; do not rescale to the original hypothesis count.
        padj = benjamini_hochberg(pv_subset)
        rejections[i] = count(<(alpha), padj)
    end

    if !tested
        @warn "independent filtering found fewer than 5 genes for every threshold; returning unfiltered results"
        return trues(length(base_mean)), NaN
    end

    valid_idx = findall(>=(0.0), rejections)
    isempty(valid_idx) && return trues(length(base_mean)), NaN

    chosen_idx = if length(valid_idx) >= 5
        theta_valid = theta[valid_idx]
        rej_valid = rejections[valid_idx]
        smooth = _lowess_fit(theta_valid, rej_valid; span=0.25, n_robust=0)
        max_fit = maximum(smooth)
        local_choice = if max_fit <= 10.0
            1
        else
            pos = rej_valid .> 0.0
            sigma = any(pos) ? sqrt(mean((rej_valid[pos] .- smooth[pos]) .^ 2)) : 0.0
            target = max_fit - sigma
            near_optimal = findfirst(>(target), smooth)
            near_optimal === nothing ? argmax(smooth) : near_optimal
        end
        valid_idx[local_choice]
    else
        valid_idx[argmax(rejections[valid_idx])]
    end

    best_thr = quantile(valid_bm, theta[chosen_idx])
    best_mask = base_mean .> best_thr
    return best_mask, best_thr
end

function _apply_cooks_outlier_filter!(tab::DataFrame, dds::DESeqDataSet)
    dds.cooks === nothing && return tab
    X = dds.model_matrix === nothing ? _design_model_matrix(dds.design)[1] : dds.model_matrix
    p = size(X, 2)
    m = size(X, 1)
    samples_for_cooks = _n_or_more_in_cell(X, 3)
    (!any(samples_for_cooks) || m <= p) && return tab
    cutoff = get(dds.metadata, "cooksCutoff", quantile(FDist(p, max(m - p, 1)), 0.99))

    pvals = Vector{Union{Missing,Float64}}(tab[!, :pvalue])

    for g in 1:nrow(tab)
        row = dds.cooks[g, samples_for_cooks]
        mx = maximum(filter(isfinite, row); init=0.0)
        if mx > cutoff
            pvals[g] = missing
        end
    end

    tab[!, :pvalue] = pvals
    tab[!, :padj] = fill(missing, nrow(tab))
    return tab
end

function _apply_independent_filter!(tab::DataFrame, dds::DESeqDataSet; alpha::Real=0.1)
    nrow(tab) == 0 && return tab
    base_mean = Float64.(tab[!, :base_mean])
    pvalue = Float64[(ismissing(v) || !isfinite(Float64(v))) ? NaN : Float64(v) for v in tab[!, :pvalue]]

    mask, thr = _independent_filter_mask(base_mean, pvalue, alpha)
    dds.metadata["indepFilterCutoff"] = thr

    # Preserve raw p-values; independent filtering only affects adjusted p-values.
    tab[!, :_indep_keep] = mask
    tab[!, :padj] = fill(missing, nrow(tab))
    return tab
end

function results(dds::DESeqDataSet;
    test=nothing,
    alpha::Real=0.1,
    fdr_method::Symbol=:BH,
    independent_filter::Bool=true,
    return_metadata::Bool=false)

    t = test === nothing ? dds.test : test
    tab = if t == :Wald && dds.wald_results !== nothing
        DataFrame(dds.wald_results)
    elseif t == :LRT && dds.lrt_results !== nothing
        DataFrame(dds.lrt_results)
    else
        throw(ArgumentError("no results available for test=$(t)"))
    end

    tab = _canonicalize_results!(tab)

    # Apply Cook's filter first
    tab = _apply_cooks_outlier_filter!(tab, dds)
    # Then apply independent filtering
    independent_filter && (tab = _apply_independent_filter!(tab, dds; alpha=alpha))

    p_for_adjust = Vector{Union{Missing,Float64}}(tab[!, :pvalue])
    if independent_filter && (:_indep_keep in propertynames(tab))
        keep = Bool.(tab[!, :_indep_keep])
        for i in eachindex(keep)
            keep[i] || (p_for_adjust[i] = missing)
        end
    end
    tab[!, :padj] = _adjust_padj(p_for_adjust, tab[!, :base_mean], fdr_method, alpha)
    (:_indep_keep in propertynames(tab)) && select!(tab, Not(:_indep_keep))
    cutoff = get(dds.metadata, "indepFilterCutoff", missing)
    try
        DataFrames.metadata!(tab, "indepFilterCutoff", cutoff; style=:note)
    catch
    end
    if return_metadata
        return (
            table=with_provenance(tab, "DEResultTable", "DifferentialExpression/results"; parameters=(test=t, alpha=Float64(alpha), fdr_method=fdr_method, independent_filter=Bool(independent_filter), row_count=nrow(tab))),
            indepFilterCutoff=cutoff,
            cooksCutoff=get(dds.metadata, "cooksCutoff", missing))
    end
    return with_provenance(tab, "DEResultTable", "DifferentialExpression/results"; parameters=(test=t, alpha=Float64(alpha), fdr_method=fdr_method, independent_filter=Bool(independent_filter), row_count=nrow(tab)))
end

function results(dds::DESeqDataSet, contrast;
    test=nothing,
    alpha::Real=0.1,
    fdr_method::Symbol=:BH,
    independent_filter::Bool=true,
    return_metadata::Bool=false)

    t = test === nothing ? dds.test : test
    if t == :LRT
        return results(dds; test=:LRT, alpha=alpha,
            fdr_method=fdr_method, independent_filter=independent_filter,
            return_metadata=return_metadata)
    end

    (dds.wald_beta_matrix === nothing || dds.wald_covariances === nothing) &&
        return results(dds; test=:Wald, alpha=alpha,
            fdr_method=fdr_method, independent_filter=independent_filter,
            return_metadata=return_metadata)

    cnames = dds.coefficient_names === nothing ?
             ["Intercept"; ["coef_$(i)" for i in 2:size(dds.wald_beta_matrix, 2)]] :
             dds.coefficient_names

    cv, _ = _named_contrast_vector(cnames, dds.design, contrast)

    wc = dds.replaced_counts === nothing ? dds.counts : dds.replaced_counts
    nf = _effective_norm_factors(dds)
    ng = length(wc.gene_ids)

    rows = Vector{DEResult}(undef, ng)
    for g in 1:ng
        beta = vec(dds.wald_beta_matrix[g, :])
        cov = dds.wald_covariances[g]
        base_mean = mean(_as_dense_row(wc.counts, g) ./ nf)
        cov_near_zero = maximum(abs.(cov)) <= eps(Float64)
        if base_mean <= 0.0 || !all(isfinite, beta) || !all(isfinite, cov) || cov_near_zero
            rows[g] = DEResult(wc.gene_ids[g], base_mean, NaN, NaN, NaN, NaN, NaN, false, false)
            continue
        end
        lfc, lfc_se, stat, pvalue = _contrast_statistics(beta, cov, cv)
        if dds.wald_results !== nothing && g <= length(dds.wald_results)
            base_row = dds.wald_results[g]
            if !base_row.converged && isfinite(base_row.log2_fold_change) && abs(base_row.log2_fold_change - lfc) > 0.5
                lfc = base_row.log2_fold_change
            end
        end
        rows[g] = DEResult(wc.gene_ids[g], base_mean, lfc, lfc_se, stat, pvalue, NaN, false, true)
    end

    tab = DataFrame(rows)
    tab = _canonicalize_results!(tab)

    tab = _apply_cooks_outlier_filter!(tab, dds)
    independent_filter && (tab = _apply_independent_filter!(tab, dds; alpha=alpha))

    # Keep adjustment semantics aligned with DESeq2: only non-filtered genes
    # contribute to multiple-testing correction; filtered genes retain missing padj.
    p_for_adjust = Vector{Union{Missing,Float64}}(tab[!, :pvalue])
    if independent_filter && (:_indep_keep in propertynames(tab))
        keep = Bool.(tab[!, :_indep_keep])
        for i in eachindex(keep)
            keep[i] || (p_for_adjust[i] = missing)
        end
    end

    tab[!, :padj] = _adjust_padj(p_for_adjust, tab[!, :base_mean], fdr_method, alpha)
    (:_indep_keep in propertynames(tab)) && select!(tab, Not(:_indep_keep))
    cutoff = get(dds.metadata, "indepFilterCutoff", missing)
    try
        DataFrames.metadata!(tab, "indepFilterCutoff", cutoff; style=:note)
    catch
    end
    if return_metadata
        return (
            table=with_provenance(tab, "DEResultTable", "DifferentialExpression/results_contrast"; notes=["contrast-specific Wald results"], parameters=(alpha=Float64(alpha), fdr_method=fdr_method, independent_filter=Bool(independent_filter), row_count=nrow(tab))),
            indepFilterCutoff=cutoff,
            cooksCutoff=get(dds.metadata, "cooksCutoff", missing))
    end
    return with_provenance(tab, "DEResultTable", "DifferentialExpression/results_contrast"; notes=["contrast-specific Wald results"], parameters=(alpha=Float64(alpha), fdr_method=fdr_method, independent_filter=Bool(independent_filter), row_count=nrow(tab)))
end

function _results_name_alias(name::String)
    m = match(r"^([A-Za-z0-9]+)_([^_]+)_vs_([^_]+)$", name)
    if m !== nothing
        factor, numerator, _ = m.captures
        return "$(factor)$(numerator)"
    end
    m2 = match(r"^([A-Za-z0-9]+)_([^_]+)$", name)
    if m2 !== nothing
        factor, level = m2.captures
        return "$(factor)$(level)"
    end
    return name
end

function resultsNames(dds::DESeqDataSet;
    deseq2_compatible::Bool=false,
    alias_map::Bool=false)
    rnames = if dds.coefficient_names !== nothing
        copy(dds.coefficient_names)
    elseif dds.model_matrix !== nothing
        ["Intercept"; ["coef_$(i)" for i in 2:size(dds.model_matrix, 2)]]
    else
        String[]
    end
    if alias_map
        return Dict(name => _results_name_alias(name) for name in rnames)
    end
    return deseq2_compatible ? [_results_name_alias(name) for name in rnames] : rnames
end

# =============================================================================
# DESeq High-Level Runner
# =============================================================================

function DESeq(dds::DESeqDataSet;
    test::Symbol=:Wald,
    fitType::Symbol=:parametric,
    sfType::Symbol=:ratio,
    betaPrior::Bool=false,
    betaPriorVar=nothing,
    full=dds.design,
    reduced=nothing,
    quiet::Bool=false,
    minReplicatesForReplace::Integer=7,
    modelMatrixType::Symbol=:standard,
    useT::Bool=false,
    betaTol::Real=1e-8,
    maxit::Integer=100,
    minmu::Real=0.5,
    reference_level=nothing,
    target_level=nothing,
    huber::Bool=false,
    workflow::Symbol=:prior,
    fdr_method::Symbol=:BH,
    shrink_type::Symbol=:normal,
    use_cv::Bool=true)

    dds.sf_type = sfType
    dds.fit_type = fitType

    # Dispersion should be estimated with standard model regardless of betaPrior

    if dds.size_factors === nothing && dds.normalization_factors === nothing
        estimateSizeFactors(dds; sfType=sfType)
    end

    # Use standard model for dispersion estimation
    estimateDispersions(dds; workflow=workflow, fitType=fitType,
        huber=huber, use_cv=use_cv)

    # Use expanded model for Wald test only if betaPrior=true
    effectiveModelMatrixType = (test == :Wald && betaPrior) ? :expanded : modelMatrixType

    if test == :Wald
        nbinomWaldTest(dds;
            betaPrior=betaPrior,
            betaPriorVar=betaPriorVar,
            modelMatrixType=effectiveModelMatrixType,
            maxit=maxit,
            betaTol=betaTol,
            minmu=minmu,
            minReplicatesForReplace=minReplicatesForReplace,
            reference_level=reference_level,
            target_level=target_level)
    elseif test == :LRT
        nbinomLRT(dds;
            modelMatrixType=effectiveModelMatrixType,
            reduced=reduced,
            maxit=maxit,
            betaTol=betaTol,
            minmu=minmu,
            minReplicatesForReplace=minReplicatesForReplace,
            reference_level=reference_level,
            target_level=target_level)
    else
        throw(ArgumentError("test must be :Wald or :LRT"))
    end

    dds.metadata["fdrMethod"] = fdr_method
    dds.metadata["shrinkType"] = shrink_type
    dds.metadata["version"] = "BioToolkit-DE-v4"
    update_provenance!(dds.metadata; source="DifferentialExpression/DESeq", parameters=(test=test, fitType=fitType, sfType=sfType, betaPrior=Bool(betaPrior), workflow=workflow, fdr_method=fdr_method, shrink_type=shrink_type, modelMatrixType=modelMatrixType))
    return dds
end

# =============================================================================
# edgeR-Style API
# =============================================================================

function trend_dispersion_edgeR(mean_expression::AbstractVector{<:Real},
    dispersions::AbstractVector{<:Real}; span::Real=0.5)
    means = Float64.(mean_expression)
    disps = Float64.(dispersions)
    return _local_dispersion_smoother(max.(means, 1e-8), max.(disps, 1e-8); span=span)
end

function dispersion_prior_dof(tagwise::AbstractVector{<:Real},
    trended::AbstractVector{<:Real})
    lr = log.(max.(Float64.(tagwise), 1e-8)) .- log.(max.(Float64.(trended), 1e-8))
    lr = lr[isfinite.(lr)]
    length(lr) < 3 && return 10.0
    v = var(lr; corrected=true)
    (!isfinite(v) || v <= 0.0) && return 10.0
    return clamp(2.0 * _trigamma_inverse(v), 2.0, 1e3)
end

function estimate_dispersion_edgeR(cm::CountMatrix, design;
    norm_factors=nothing)

    nf = norm_factors === nothing ? calc_norm_factors(cm; method=:tmm) : Float64.(norm_factors)
    gw = estimate_dispersions(cm, nf; workflow=:simple)

    means = zeros(Float64, length(cm.gene_ids))
    for g in 1:length(cm.gene_ids)
        means[g] = mean(_as_dense_row(cm.counts, g) ./ nf)
    end

    tr = trend_dispersion_edgeR(means, gw)
    pdof = dispersion_prior_dof(gw, tr)
    n = max(length(cm.sample_ids), 2)
    tag = (pdof .* tr .+ (n - 1) .* gw) ./ (pdof + n - 1)

    return (
        common=median(tag),
        trended=tr,
        tagwise=max.(tag, 1e-8),
        prior_dof=pdof,
        norm_factors=nf)
end

function _edgeR_exact_twosided_pvalue(total::Int, x_obs::Int, r::Float64)
    total <= 0 && return 1.0
    x_obs = clamp(x_obs, 0, total)

    # For moderate totals, use exact conditional probability ordering
    if total <= 5000
        # NB conditional: P(X2=k | X1+X2=total) ∝ C(k+r-1,k) * C(total-k+r-1,total-k)
        # No p2 terms — p drops out under H0
        logp = Vector{Float64}(undef, total + 1)
        @inbounds for k in 0:total
            logp[k+1] = lgamma(k + r) - lgamma(r) - lgamma(k + 1.0) +
                        lgamma(total - k + r) - lgamma(r) - lgamma(total - k + 1.0)
        end
        lp_obs = logp[clamp(x_obs, 0, total)+1]
        sel = [v for v in logp if v <= lp_obs + 1e-12]
        return clamp(exp(_logsumexp(sel) - _logsumexp(logp)), 0.0, 1.0)
    end

    # Large-count fallback: normal approximation to NB conditional
    # Under H0, X2 ~ Binom(total, 0.5) when p1=p2, adjusted for overdispersion
    mu = total / 2.0
    phi = 1.0 / r
    sigma2 = total * 0.25 * (1.0 + (total - 1.0) * phi / (1.0 + phi))
    sigma = sqrt(max(sigma2, eps(Float64)))
    z = abs((x_obs - mu) / sigma)
    return clamp(2.0 * ccdf(Normal(), z), 0.0, 1.0)
end

# Backward-compatible overload used by older call sites/tests that passed
# two shape-like parameters; reduce to a symmetric effective size.
function _edgeR_exact_twosided_pvalue(total::Int, x_obs::Int, a::Real, b::Real)
    r = max((Float64(a) + Float64(b)) / 2.0, eps(Float64))
    return _edgeR_exact_twosided_pvalue(total, x_obs, r)
end

function exact_test_edgeR(cm::CountMatrix, design;
    norm_factors=nothing,
    dispersion=nothing)

    groups = if design isa AbstractVector
        String.(design)
    elseif design isa DataFrame && :condition in names(design)
        String.(design[!, :condition])
    else
        throw(ArgumentError("design for exact_test_edgeR must be a vector or DataFrame with :condition"))
    end
    levs = unique(groups)
    length(levs) == 2 || throw(ArgumentError("exact_test_edgeR currently supports exactly two groups"))
    g1 = findall(groups .== levs[1])
    g2 = findall(groups .== levs[2])

    nf = norm_factors === nothing ? calc_norm_factors(cm; method=:tmm) : Float64.(norm_factors)
    disp = if dispersion === nothing
        estimate_dispersions(cm, nf; workflow=:simple)
    elseif isa(dispersion, Number)
        fill(Float64(dispersion), length(cm.gene_ids))
    else
        Float64.(dispersion)
    end

    lib_sizes = vec(sum(cm.counts, dims=1))
    eff_lib = max.(lib_sizes .* nf, 1.0)
    lib1 = max(sum(eff_lib[g1]), 1.0)
    lib2 = max(sum(eff_lib[g2]), 1.0)

    ng = length(cm.gene_ids)
    logfc = zeros(Float64, ng)
    logcpm = zeros(Float64, ng)
    pval = fill(1.0, ng)

    for g in 1:ng
        y = _as_dense_row(cm.counts, g)
        x1 = Int(round(sum(y[g1])))
        x2 = Int(round(sum(y[g2])))
        total = x1 + x2

        cpm1 = (x1 + 0.5) / (lib1 + 1.0) * 1e6
        cpm2 = (x2 + 0.5) / (lib2 + 1.0) * 1e6
        logfc[g] = log2(cpm2 / cpm1)
        logcpm[g] = log2(((x1 + x2) + 0.5) / (lib1 + lib2 + 1.0) * 1e6)

        if total <= 0
            pval[g] = 1.0
            continue
        end

        phi = max(disp[g], 1e-8)
        r = 1.0 / phi
        pval[g] = _edgeR_exact_twosided_pvalue(total, x2, r)
    end

    padj = benjamini_hochberg(pval)

    return DataFrame(
        gene_id=cm.gene_ids,
        logFC=logfc,
        logCPM=logcpm,
        pvalue=pval,
        padj=padj)
end

function _trigamma_inverse(y::Float64)
    y > 0 || return Inf
    # Initial guess using asymptotic approximation 1/y
    x = y > 1e-6 ? 1.0 / y : 0.5
    for _ in 1:30
        tg = trigamma(x)
        tg2 = polygamma(2, x)   # tetragamma = derivative of trigamma
        delta = (tg - y) / tg2
        x -= delta
        x = max(x, 1e-6)       # keep in valid domain
        abs(delta) < 1e-9 && break
    end
    return x
end

function _squeeze_variances(var_g::Vector{Float64}, df_resid::Int)
    # Fit inverse-chi-squared prior to gene-wise variances
    # Returns (squeezed_var, df_prior, var_prior)

    valid = isfinite.(var_g) .& (var_g .> 0.0)
    count_valid = count(valid)
    count_valid < 3 && return (copy(var_g), 0.0, fill(median(var_g[valid]), count(valid)))

    log_var = log.(max.(var_g[valid], 1e-8))
    mean_log_var = mean(log_var)
    s2 = var(log_var; corrected=true)

    if s2 <= 0.0 || !isfinite(s2)
        return (copy(var_g), 0.0, fill(exp(mean_log_var), count_valid))
    end

    # Iteratively estimate df_prior and var_prior
    # E[log(chi-sq_df/df)] = digamma(df/2) - log(df/2)
    # Var[log(chi-sq_df/df)] = trigamma(df/2)

    df_prior = 2.0
    for _ in 1:20
        sampling_var = trigamma(df_prior / 2.0)
        prior_var = max(s2 - sampling_var, 1e-8)
        # Method of moments for df_prior
        if prior_var > 0
            df_prior_new = 2.0 * _trigamma_inverse(max(prior_var, 1e-8))
            df_prior_new = clamp(df_prior_new, 0.5, 1e6)
            if abs(df_prior_new - df_prior) < 0.01
                df_prior = df_prior_new
                break
            end
            df_prior = df_prior_new
        else
            break
        end
    end

    var_prior = exp(mean_log_var)
    df_prior = clamp(df_prior, 0.0, 1e6)

    # Squeeze: weighted average

    squeezed = similar(var_g)
    for i in eachindex(var_g)
        if !valid[i]
            squeezed[i] = var_g[i]
        else
            # Posterior mean of variance under inverse-chi-squared prior
            denom = df_prior + df_resid
            if denom > 0
                squeezed[i] = (df_prior * var_prior + df_resid * var_g[i]) / denom
            else
                squeezed[i] = var_prior
            end
        end
    end

    return (max.(squeezed, 1e-8), df_prior, fill(var_prior, length(var_g)))
end

_squeeze_quasi_dispersions(var_g::Vector{Float64}, df_resid::Int) =
    _squeeze_variances(var_g, df_resid)

function glm_ql_fit(cm::CountMatrix, design;
    norm_factors=nothing,
    dispersion=nothing,
    maxit::Integer=100,
    beta_tol::Real=1e-8,
    minmu::Real=0.5)

    X, cnames, _ = _design_model_matrix(design)
    nf = norm_factors === nothing ? calc_norm_factors(cm; method=:tmm) : Float64.(norm_factors)
    disp = if dispersion === nothing
        estimate_dispersions(cm, nf; workflow=:simple, model_matrix=X)
    elseif isa(dispersion, Number)
        fill(Float64(dispersion), length(cm.gene_ids))
    else
        Float64.(dispersion)
    end

    ng = length(cm.gene_ids)
    p = size(X, 2)
    beta_mat = zeros(Float64, ng, p)
    covs = Vector{Matrix{Float64}}(undef, ng)
    qldisp_raw = ones(Float64, ng)

    solver = GLMSolver(X)
    offset = log.(nf)
    df_resid = max(size(X, 1) - p, 1)

    for g in 1:ng
        y = _as_dense_row(cm.counts, g)
        if sum(y) <= 0.0
            beta_mat[g, :] .= 0.0
            covs[g] = zeros(Float64, p, p)
            qldisp_raw[g] = 1.0
            continue
        end

        beta, cov, _, _ = fit_gene_fast!(solver, y, offset, disp[g];
            maxit=maxit, beta_tol=beta_tol, minmu=minmu, return_covariance=true)

        beta_mat[g, :] .= beta
        covs[g] = cov

        var_nb = solver.mu .+ disp[g] .* solver.mu .^ 2
        pearson2 = ((y .- solver.mu) .^ 2) ./ max.(var_nb, eps(Float64))
        qldisp_raw[g] = max(sum(pearson2) / df_resid, 1e-4)
    end

    qldisp, ql_prior_dof, ql_prior_var = _squeeze_variances(qldisp_raw, df_resid)

    return (
        coefficients=beta_mat,
        covariances=covs,
        quasi_dispersion=qldisp,
        quasi_dispersion_raw=qldisp_raw,
        quasi_prior_dof=ql_prior_dof,
        quasi_prior_variance=ql_prior_var,
        dispersions=disp,
        model_matrix=X,
        residual_df=df_resid,
        coefficient_names=cnames,
        gene_ids=cm.gene_ids)
end

function glm_ql_f_test(fit;
    coef::Int=2)
    beta = fit.coefficients
    covs = fit.covariances
    ql = fit.quasi_dispersion
    X = fit.model_matrix
    p = size(X, 2)
    df_resid = hasproperty(fit, :residual_df) ? max(Float64(fit.residual_df), 1.0) : max(Float64(size(X, 1) - p), 1.0)
    df_prior = hasproperty(fit, :quasi_prior_dof) ? max(Float64(fit.quasi_prior_dof), 0.0) : 0.0
    df2 = df_resid + df_prior

    (!isfinite(df2) || df2 <= 0.0) && (df2 = df_resid)

    ng = size(beta, 1)
    j = clamp(coef, 1, size(beta, 2))

    logfc = zeros(Float64, ng)
    stat = zeros(Float64, ng)
    pval = fill(1.0, ng)

    for g in 1:ng
        b = beta[g, j]
        se2 = max(covs[g][j, j], eps(Float64))
        f = (b^2) / (se2 * max(ql[g], 1e-8))
        logfc[g] = b / log(2)
        stat[g] = f
        pval[g] = ccdf(FDist(1, df2), f)
    end

    padj = benjamini_hochberg(pval)

    return DataFrame(
        gene_index=collect(1:ng),
        logFC=logfc,
        F=stat,
        pvalue=pval,
        padj=padj)
end

function edgeR_qlf_test(cm::CountMatrix, design;
    coef::Int=2,
    norm_factors=nothing,
    dispersion=nothing,
    maxit::Integer=100,
    beta_tol::Real=1e-8,
    minmu::Real=0.5)

    fit = glm_ql_fit(cm, design;
        norm_factors=norm_factors,
        dispersion=dispersion,
        maxit=maxit,
        beta_tol=beta_tol,
        minmu=minmu)

    tab = glm_ql_f_test(fit; coef=coef)
    tab[!, :gene_id] = cm.gene_ids
    select!(tab, :gene_id, :logFC, :F, :pvalue, :padj)
    return tab
end

# =============================================================================
# LFC Shrinkage
# =============================================================================

function lfc_shrink_apeglm(beta::AbstractVector{<:Real},
    se::AbstractVector{<:Real},
    offset=nothing,
    counts=nothing,
    dispersions=nothing;
    prior_scale=nothing,
    max_iter::Int=50,
    tol::Real=1e-8)

    b = Float64.(beta)
    s = max.(Float64.(se), 1e-8)

    scale = if prior_scale === nothing
        med = median(abs.(b[isfinite.(b)]))
        isfinite(med) && med > 0.0 ? med / 0.6745 : 1.0
    else
        Float64(prior_scale)
    end
    scale = max(scale, 1e-3)
    s2_prior = scale^2

    out = copy(b)

    for i in eachindex(b)
        bh = b[i]
        se2 = s[i]^2
        if !isfinite(bh) || !isfinite(se2)
            out[i] = NaN
            continue
        end

        x = bh
        for _ in 1:max_iter
            g = (x - bh) / se2 + 2.0 * x / (s2_prior + x^2)
            h = 1.0 / se2 + 2.0 * (s2_prior - x^2) / (s2_prior + x^2)^2
            if !isfinite(g) || !isfinite(h) || abs(h) < 1e-12
                break
            end
            step = g / h
            x_new = x - step
            if !isfinite(x_new)
                break
            end
            if abs(x_new - x) < tol
                x = x_new
                break
            end
            x = x_new
        end
        out[i] = x
    end

    return out
end

"""
    lfc_shrink_apeglm_full(cm, dds; coef_idx=2, prior_scale=nothing)

Full-data apeglm shrinkage using actual NB likelihood, not normal approximation.
More accurate for low-count genes.
"""
function lfc_shrink_apeglm_full(cm::CountMatrix, dds::DESeqDataSet;
    coef_idx::Int=2,
    prior_scale=nothing,
    max_iter::Int=100,
    tol::Real=1e-8)

    nf = _effective_norm_factors(dds)
    X = dds.model_matrix
    X === nothing && throw(ArgumentError("model matrix not available"))
    offset = log.(nf)

    ng = length(cm.gene_ids)
    ns = length(cm.sample_ids)

    disp = dds.dispersions === nothing ? estimate_dispersions(cm, nf) : dds.dispersions
    beta_mle = dds.wald_beta_matrix === nothing ? zeros(ng, size(X, 2)) : dds.wald_beta_matrix

    scale = if prior_scale === nothing
        med = median(abs.(beta_mle[:, coef_idx][isfinite.(beta_mle[:, coef_idx])]))
        isfinite(med) && med > 0.0 ? med / 0.6745 : 1.0
    else
        Float64(prior_scale)
    end
    scale = max(scale, 1e-3)
    s2_prior = scale^2

    shrunk_beta = copy(beta_mle)

    for g in 1:ng
        y = _as_dense_row(cm.counts, g)
        if all(iszero, y)
            continue
        end

        alpha = max(disp[g], 1e-8)
        r = 1.0 / alpha
        beta_init = beta_mle[g, :]
        x = beta_init[coef_idx]

        # Full NB log-likelihood + Cauchy log-prior
        function objective(beta_coef)
            beta_test = copy(beta_init)
            beta_test[coef_idx] = beta_coef
            mu = max.(exp.(X * beta_test .+ offset), 1e-8)
            ll = 0.0
            @inbounds for i in 1:ns
                ll += lgamma(y[i] + r) - lgamma(r) - lgamma(y[i] + 1.0) +
                      r * (log(r) - log(r + mu[i])) + y[i] * (log(mu[i]) - log(r + mu[i]))
            end
            # Cauchy prior: -log(1 + beta^2/scale^2)
            prior = -log(1.0 + beta_coef^2 / s2_prior)
            return -(ll + prior)
        end

        try
            opt = Optim.optimize(objective, x - 5 * scale, x + 5 * scale, Optim.Brent(); iterations=max_iter)
            if Optim.converged(opt)
                shrunk_beta[g, coef_idx] = Optim.minimizer(opt)
            end
        catch
        end
    end

    return shrunk_beta
end

function lfc_shrink_ashr(beta::AbstractVector{<:Real},
    se::AbstractVector{<:Real};
    mix_sd=nothing,
    pi_weights=nothing,
    max_em_iter::Int=200,
    em_tol::Real=1e-6)

    b = Float64.(beta)
    s = max.(Float64.(se), 1e-8)
    n = length(b)

    grid = if mix_sd === nothing
        bsd = std(filter(isfinite, b))
        bsd = isfinite(bsd) && bsd > 0.0 ? bsd : 1.0
        vcat(0.0, exp.(range(log(0.1), log(4bsd + 1e-6); length=10)))
    else
        Float64.(mix_sd)
    end

    K = length(grid)
    w = if pi_weights === nothing
        fill(1.0 / K, K)
    else
        ww = Float64.(pi_weights)
        length(ww) == K || throw(ArgumentError("pi_weights length must match mix_sd"))
        ww ./ sum(ww)
    end

    resp = zeros(Float64, n, K)
    std_norm = Normal()

    # Mixture component density under ashr's unimodal uniform prior.
    @inline function comp_density(bi::Float64, si::Float64, a::Float64)
        if a == 0.0
            return pdf(Normal(0.0, si), bi)
        end
        z_lo = (-a - bi) / si
        z_hi = (a - bi) / si
        dens = (cdf(std_norm, z_hi) - cdf(std_norm, z_lo)) / (2.0 * a)
        return max(dens, eps(Float64))
    end

    # Posterior mean under θ ~ Unif(-a,a), b|θ ~ N(θ, s^2).
    @inline function comp_post_mean(bi::Float64, si::Float64, a::Float64)
        if a == 0.0
            return 0.0
        end
        z_lo = (-a - bi) / si
        z_hi = (a - bi) / si
        denom = cdf(std_norm, z_hi) - cdf(std_norm, z_lo)
        if denom <= eps(Float64)
            return clamp(bi, -a, a)
        end
        return bi + si * (pdf(std_norm, z_lo) - pdf(std_norm, z_hi)) / denom
    end

    for _ in 1:max_em_iter
        w_prev = copy(w)

        for i in 1:n
            if !isfinite(b[i])
                resp[i, :] .= 0.0
                continue
            end
            denom = 0.0
            for k in 1:K
                dens = comp_density(b[i], s[i], grid[k])
                val = w[k] * dens
                resp[i, k] = val
                denom += val
            end
            if denom <= 0.0 || !isfinite(denom)
                resp[i, :] .= 0.0
                resp[i, 1] = 1.0
            else
                resp[i, :] ./= denom
            end
        end

        for k in 1:K
            w[k] = mean(resp[:, k])
        end
        w ./= sum(w)

        if maximum(abs.(w .- w_prev)) < em_tol
            break
        end
    end

    out = similar(b)
    for i in 1:n
        if !isfinite(b[i])
            out[i] = NaN
            continue
        end
        post_mean = 0.0
        for k in 1:K
            mu_k = comp_post_mean(b[i], s[i], grid[k])
            post_mean += resp[i, k] * mu_k
        end
        out[i] = post_mean
    end

    return out
end

"""
    adaptive_shrinkage(beta, se; mode=:zero, grid_method=:uniform)

True adaptive shrinkage with unimodal assumption at mode.
More robust than normal mixture for heavy-tailed effect distributions.
"""
function adaptive_shrinkage(beta::AbstractVector{<:Real},
    se::AbstractVector{<:Real};
    mode::Symbol=:zero,
    n_grid::Int=20,
    max_em_iter::Int=200,
    em_tol::Real=1e-6)
    mode == :zero || throw(ArgumentError("adaptive_shrinkage currently supports mode=:zero only"))

    b = Float64.(beta)
    bsd = std(filter(isfinite, b))
    bsd = isfinite(bsd) && bsd > 0.0 ? bsd : 1.0
    grid = vcat(0.0, exp.(range(log(0.1 * bsd), log(5 * bsd); length=n_grid)))

    return lfc_shrink_ashr(beta, se;
        mix_sd=grid,
        pi_weights=nothing,
        max_em_iter=max_em_iter,
        em_tol=em_tol)
end

function shrink_lfc(results_vec::AbstractVector{<:DEResult};
    prior_sd=nothing,
    shrink_type::Symbol=:normal,
    mix_sd=nothing,
    pi_weights=nothing)

    n = length(results_vec)
    beta = [r.log2_fold_change for r in results_vec]
    se = [max(r.lfc_se, 1e-8) for r in results_vec]

    shrunk = if shrink_type == :normal
        tau = prior_sd === nothing ? max(std(filter(isfinite, beta)), 0.5) : Float64(prior_sd)
        [beta[i] * (tau^2 / (tau^2 + se[i]^2)) for i in 1:n]
    elseif shrink_type == :cauchy || shrink_type == :apeglm
        lfc_shrink_apeglm(beta, se; prior_scale=prior_sd)
    elseif shrink_type == :ashr
        lfc_shrink_ashr(beta, se; mix_sd=mix_sd, pi_weights=pi_weights)
    elseif shrink_type == :adaptive
        adaptive_shrinkage(beta, se)
    else
        throw(ArgumentError("unsupported shrink_type: $shrink_type"))
    end

    # Keep outlier/filtered genes unchanged, matching DESeq2's NA-aware lfcShrink behavior.
    for i in 1:n
        if !isfinite(beta[i]) || !isfinite(se[i]) || !isfinite(results_vec[i].pvalue) || !isfinite(results_vec[i].padj)
            shrunk[i] = beta[i]
        end
    end

    out = Vector{DEResult}(undef, n)
    for i in 1:n
        r = results_vec[i]
        out[i] = DEResult(r.gene_id, r.base_mean, shrunk[i], r.lfc_se,
            r.stat, r.pvalue, r.padj, r.zero_inflated, r.converged)
    end
    return out
end

function shrink_lfc(dds::DESeqDataSet;
    coef=nothing,
    contrast=nothing,
    test=nothing,
    alpha::Real=0.1,
    fdr_method::Symbol=:BH,
    independent_filter::Bool=true,
    prior_sd=nothing,
    shrink_type::Symbol=:normal,
    mix_sd=nothing,
    pi_weights=nothing)

    tab = if contrast !== nothing
        results(dds, contrast;
            test=test,
            alpha=alpha,
            fdr_method=fdr_method,
            independent_filter=independent_filter)
    elseif coef !== nothing
        cnames = resultsNames(dds)
        isempty(cnames) && throw(ArgumentError("no coefficient names available; run DESeq first"))
        idx = if coef isa Integer
            Int(coef)
        elseif coef isa AbstractString
            found = findfirst(==(String(coef)), cnames)
            found === nothing && throw(ArgumentError("coef $(coef) not found in resultsNames(dds)"))
            found
        else
            throw(ArgumentError("coef must be Integer or String"))
        end
        (1 <= idx <= length(cnames)) || throw(ArgumentError("coef index out of bounds"))
        cv = zeros(Float64, length(cnames))
        cv[idx] = 1.0
        results(dds, cv;
            test=test,
            alpha=alpha,
            fdr_method=fdr_method,
            independent_filter=independent_filter)
    else
        results(dds;
            test=test,
            alpha=alpha,
            fdr_method=fdr_method,
            independent_filter=independent_filter)
    end

    res = _table_to_deresults(tab; preserve_missing=true)
    shrunk = shrink_lfc(res;
        prior_sd=prior_sd,
        shrink_type=shrink_type,
        mix_sd=mix_sd,
        pi_weights=pi_weights)
    return _canonicalize_results!(DataFrame(shrunk))
end

function lfcShrink(dds::DESeqDataSet;
    coef=nothing,
    contrast=nothing,
    type::Symbol=:normal,
    kwargs...)
    return shrink_lfc(dds;
        coef=coef,
        contrast=contrast,
        shrink_type=type,
        kwargs...)
end

# =============================================================================
# Transformations and Filtering
# =============================================================================

function filter_low_counts(cm::CountMatrix;
    min_total::Real=10.0,
    min_count::Real=1.0,
    min_samples::Int=1)

    Y = cm.counts
    ng = size(Y, 1)
    keep = falses(ng)

    for g in 1:ng
        row = _as_dense_row(Y, g)
        keep[g] = sum(row) >= min_total && count(>=(min_count), row) >= min_samples
    end

    idx = findall(identity, keep)
    result = CountMatrix(Y[idx, :], cm.gene_ids[idx], cm.sample_ids)
    update_provenance!(result.metadata; source="DifferentialExpression/filter_low_counts", notes=["filtered genes by count thresholds"], parameters=(input_genes=ng, output_genes=length(idx), min_total=Float64(min_total), min_count=Float64(min_count), min_samples=Int(min_samples)))
    return result
end

"""
    vst(cm, dispersions; norm_factors=nothing)

Variance-stabilizing transformation for NB counts.
"""
function vst(cm::CountMatrix,
    dispersions::AbstractVector{<:Real};
    norm_factors=nothing)

    nf = norm_factors === nothing ? calc_norm_factors(cm; method=:tmm) : Float64.(norm_factors)
    Y = Matrix{Float64}(cm.counts)
    ng, ns = size(Y)
    length(dispersions) == ng || throw(ArgumentError("dispersion length mismatch"))

    out = zeros(Float64, ng, ns)
    for g in 1:ng
        alpha = max(Float64(dispersions[g]), 0.0)
        for s in 1:ns
            x = max(Y[g, s] / nf[s], 0.0)
            if alpha < 1e-10
                out[g, s] = 2.0 * sqrt(x) / log(2)
            else
                out[g, s] = (2.0 / sqrt(alpha)) * asinh(sqrt(max(alpha * x, 0.0))) / log(2)
            end
        end
    end
    update_provenance!(cm.metadata; source="DifferentialExpression/vst", notes=["computed variance-stabilizing transform"], parameters=(n_genes=ng, n_samples=ns, used_norm_factors=norm_factors !== nothing))
    return out
end

"""
    voom(cm; norm_factors=nothing, design=nothing)

Mean-variance trend modeling for RNA-seq (Law et al. 2014).
Returns (log2cpm, precision_weights) for limma-style analysis.
"""
function voom(cm::CountMatrix;
    norm_factors=nothing,
    design=nothing,
    lowess_span::Real=0.5)

    nf = norm_factors === nothing ? calc_norm_factors(cm; method=:tmm) : Float64.(norm_factors)
    Y = Matrix{Float64}(cm.counts)
    ng, ns = size(Y)

    lib = vec(sum(Y, dims=1))
    eff_lib = max.(lib .* nf, 1.0)

    # Compute log2-CPM
    log2cpm = zeros(Float64, ng, ns)
    for s in 1:ns
        log2cpm[:, s] .= log2.((Y[:, s] .+ 0.5) ./ (eff_lib[s] / 1e6))
    end

    # Fit design and get residuals
    X = if design === nothing
        ones(Float64, ns, 1)
    else
        _design_model_matrix(design)[1]
    end

    coef = X \ Matrix{Float64}(log2cpm)'
    fitted = (X * coef)'
    resid = log2cpm .- fitted

    # Mean-variance trend
    mean_expr = vec(mean(log2cpm, dims=2))
    sd_resid = vec(std(resid, dims=2; corrected=false))

    # Lowess trend of sqrt(SD) vs mean
    valid = isfinite.(mean_expr) .& isfinite.(sd_resid) .& (sd_resid .> 0.0)
    if count(valid) < 10
        weights = fill(1.0, (ng, ns))
        update_provenance!(cm.metadata; source="DifferentialExpression/voom", status=:warning, notes=["computed voom transformation"], fallbacks=["insufficient valid genes for lowess fit; used unit precision weights"], parameters=(n_genes=ng, n_samples=ns, lowess_span=Float64(lowess_span), used_norm_factors=norm_factors !== nothing, used_design=design !== nothing))
        return log2cpm, weights
    end

    mx = mean_expr[valid]
    sy = sqrt.(sd_resid[valid])

    trend = _lowess_fit(mx, sy; span=lowess_span)

    # Sort once outside the gene loop — _lowess_predict requires sorted abscissa
    sort_order = sortperm(mx)
    mx_sorted = mx[sort_order]
    trend_sorted = trend[sort_order]

    weights = zeros(Float64, ng, ns)
    for g in 1:ng
        if !isfinite(mean_expr[g])
            weights[g, :] .= 1.0
            continue
        end
        predicted_sqrt_sd = _lowess_predict(mx_sorted, trend_sorted, mean_expr[g])
        predicted_sqrt_sd = max(predicted_sqrt_sd, 0.1)
        weights[g, :] .= 1.0 / predicted_sqrt_sd^4
    end

    update_provenance!(cm.metadata; source="DifferentialExpression/voom", notes=["computed voom mean-variance transformation"], parameters=(n_genes=ng, n_samples=ns, lowess_span=Float64(lowess_span), used_norm_factors=norm_factors !== nothing, used_design=design !== nothing))
    return log2cpm, weights
end

function _lowess_fit(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; span::Real=0.5, n_robust::Int=3)
    n = length(x)
    n == length(y) || throw(ArgumentError("x/y length mismatch in _lowess_fit"))
    n == 0 && return Float64[]
    k = max(3, round(Int, span * n))
    k = min(k, n)
    order = sortperm(x)
    xs = Float64.(x[order])
    ys = Float64.(y[order])
    trend = zeros(Float64, n)
    robust_w = ones(Float64, n)

    function fit_once!(out::Vector{Float64}, rw::Vector{Float64})
        for i in 1:n
            lo = max(1, i - k ÷ 2)
            hi = min(n, lo + k - 1)
            lo = max(1, hi - k + 1)
            x0 = xs[i]
            idx = lo:hi
            window_x = xs[idx]
            window_y = ys[idx]

            dist = abs.(window_x .- x0)
            max_dist = maximum(dist) + eps(Float64)
            tri = (1.0 .- (dist ./ max_dist) .^ 3) .^ 3
            w = tri .* rw[idx]
            wsum = sum(w)

            if wsum <= eps(Float64)
                out[i] = mean(window_y)
                continue
            end

            w ./= wsum
            wx = sum(w .* window_x)
            wy = sum(w .* window_y)
            wxx = sum(w .* window_x .^ 2)
            wxy = sum(w .* window_x .* window_y)
            denom = wxx - wx^2
            if abs(denom) > 1e-10
                slope = (wxy - wx * wy) / denom
                intercept = wy - slope * wx
                out[i] = intercept + slope * x0
            else
                out[i] = mean(window_y)
            end
        end
    end

    fit_once!(trend, robust_w)
    for _ in 1:max(n_robust, 0)
        resid = ys .- trend
        med = median(abs.(resid))
        scale = max(6.0 * med, eps(Float64))
        u = abs.(resid) ./ scale
        robust_w .= (1.0 .- min.(u, 1.0) .^ 2) .^ 2
        fit_once!(trend, robust_w)
    end

    # Unsort
    unsorted = zeros(Float64, n)
    unsorted[order] = trend
    return unsorted
end

function _lowess_predict(x_train::AbstractVector{<:Real}, y_train::AbstractVector{<:Real}, xq::Real)
    n = length(x_train)
    xq < minimum(x_train) && return y_train[argmin(x_train)]
    xq > maximum(x_train) && return y_train[argmax(x_train)]

    idx = searchsortedlast(x_train, xq)
    idx = clamp(idx, 1, n - 1)
    x0, x1 = x_train[idx], x_train[idx+1]
    if abs(x1 - x0) < eps(Float64)
        return y_train[idx]
    end
    t = (xq - x0) / (x1 - x0)
    return y_train[idx] * (1 - t) + y_train[idx+1] * t
end

# =============================================================================
# Standalone Differential Expression
# =============================================================================

function _table_to_deresults(tab::DataFrame; preserve_missing::Bool=false)
    out = Vector{DEResult}(undef, nrow(tab))
    for i in 1:nrow(tab)
        p = tab[i, :pvalue]
        q = tab[i, :padj]
        lfc = tab[i, :log2_fold_change]
        se = tab[i, :lfc_se]
        st = tab[i, :stat]
        out[i] = DEResult(
            String(tab[i, :gene_id]),
            Float64(tab[i, :base_mean]),
            ismissing(lfc) ? NaN : Float64(lfc),
            ismissing(se) ? NaN : Float64(se),
            ismissing(st) ? NaN : Float64(st),
            ismissing(p) ? (preserve_missing ? NaN : 1.0) : Float64(p),
            ismissing(q) ? (preserve_missing ? NaN : 1.0) : Float64(q),
            false,
            true)
    end
    return out
end

"""
    differential_expression(cm, design_vec; ...)

High-level entry point for DE analysis.
"""
function differential_expression(cm::CountMatrix, design_vec;
    test::Symbol=:Wald,
    fitType::Symbol=:parametric,
    sfType::Symbol=:ratio,
    betaPrior::Bool=false,
    shrink::Bool=false,
    shrink_type::Symbol=:normal,
    reference_level=nothing,
    target_level=nothing,
    fdr_method::Symbol=:BH,
    alpha::Real=0.1,
    min_total::Real=0,
    normalization_method::Union{Nothing,Symbol}=nothing,
    dispersion_workflow::Symbol=:prior,
    modelMatrixType::Symbol=:standard)

    length(design_vec) == length(cm.sample_ids) ||
        throw(ArgumentError("design length must match sample count"))

    filtered_cm = min_total > 0 ? filter_low_counts(cm; min_total=min_total) : cm
    isempty(filtered_cm.gene_ids) && return DEResult[]

    chosen_sf = normalization_method === nothing ? sfType : normalization_method

    coldata = DataFrame(condition=design_vec)
    dds = DESeqDataSet(filtered_cm, coldata, design_vec)

    DESeq(dds;
        test=test,
        fitType=fitType,
        sfType=chosen_sf,
        betaPrior=betaPrior,
        reference_level=reference_level,
        target_level=target_level,
        fdr_method=fdr_method,
        shrink_type=shrink_type,
        workflow=dispersion_workflow,
        modelMatrixType=modelMatrixType)

    tab = if test == :Wald && (reference_level !== nothing || target_level !== nothing)
        if reference_level === nothing || target_level === nothing
            results(dds; test=:Wald, alpha=alpha, fdr_method=fdr_method)
        else
            results(dds, [:condition, target_level, reference_level];
                test=:Wald,
                alpha=alpha,
                fdr_method=fdr_method)
        end
    else
        results(dds; test=test, alpha=alpha, fdr_method=fdr_method)
    end

    out = _table_to_deresults(tab)
    if shrink
        shrink_input = _table_to_deresults(tab; preserve_missing=true)
        shrunk = shrink_lfc(shrink_input; shrink_type=shrink_type)
        for i in eachindex(out)
            r = out[i]
            out[i] = DEResult(r.gene_id, r.base_mean, shrunk[i].log2_fold_change,
                r.lfc_se, r.stat, r.pvalue, r.padj, r.zero_inflated, r.converged)
        end
    end

    _ctx = active_provenance_context()
    _ctx !== nothing && register_provenance!(_ctx, "differential_expression"; parents=provenance_parent_ids(cm), parameters=(test=test, fitType=fitType, sfType=sfType, betaPrior=betaPrior, shrink=shrink, shrink_type=shrink_type, fdr_method=fdr_method, alpha=alpha, min_total=min_total, normalization_method=normalization_method, dispersion_workflow=dispersion_workflow, modelMatrixType=modelMatrixType))
    provenance = provenance_record(
        "DEResult",
        "DifferentialExpression/differential_expression";
        notes=["ran standalone differential expression workflow"],
        fallbacks=min_total > 0 ? ["filtered low-count genes before fitting"] : String[],
        parameters=(test=test, fitType=fitType, sfType=chosen_sf, betaPrior=Bool(betaPrior), shrink=Bool(shrink), shrink_type=shrink_type, fdr_method=fdr_method, alpha=Float64(alpha), min_total=Float64(min_total), normalization_method=normalization_method === nothing ? "default" : normalization_method, dispersion_workflow=dispersion_workflow, modelMatrixType=modelMatrixType, result_count=length(out)))
    for i in eachindex(out)
        result = out[i]
        out[i] = DEResult(
            result.gene_id,
            result.base_mean,
            result.log2_fold_change,
            result.lfc_se,
            result.stat,
            result.pvalue,
            result.padj,
            result.zero_inflated,
            result.converged,
            provenance)
    end
    update_provenance!(cm.metadata; source="DifferentialExpression/differential_expression", notes=["ran standalone differential expression workflow"], fallbacks=min_total > 0 ? ["filtered low-count genes before fitting"] : String[], parameters=(test=test, fitType=fitType, sfType=chosen_sf, betaPrior=Bool(betaPrior), shrink=Bool(shrink), shrink_type=shrink_type, fdr_method=fdr_method, alpha=Float64(alpha), min_total=Float64(min_total), normalization_method=normalization_method === nothing ? "default" : normalization_method, dispersion_workflow=dispersion_workflow, modelMatrixType=modelMatrixType, result_count=length(out)))
    return out
end

function differential_expression(count_data::AbstractMatrix{<:Integer}, design_vec;
    gene_ids=nothing,
    sample_ids=nothing,
    kwargs...)

    ng, ns = size(count_data)
    gids = gene_ids === nothing ? ["gene_$(i)" for i in 1:ng] : String.(gene_ids)
    sids = sample_ids === nothing ? ["sample_$(i)" for i in 1:ns] : String.(sample_ids)
    cm = CountMatrix(count_data, gids, sids)
    return differential_expression(cm, design_vec; kwargs...)
end

# =============================================================================
# Batch Correction
# =============================================================================

function _one_hot_labels(labels::AbstractVector{<:AbstractString})
    levs = unique(labels)
    n = length(labels)
    if length(levs) <= 1
        return zeros(Float64, n, 0), levs
    end
    B = zeros(Float64, n, length(levs) - 1)
    ref = levs[1]
    for (j, lev) in enumerate(levs[2:end])
        B[:, j] .= Float64.(labels .== lev)
    end
    return B, levs
end

function _ensure_sample_design(X::AbstractMatrix{<:Real}, nsamples::Int)
    Xf = Matrix{Float64}(X)
    size(Xf, 1) == nsamples ||
        throw(ArgumentError("design rows must equal number of samples"))
    return Xf
end

"""
    remove_batch_effect(data, batch; bio_design=nothing)

Linear batch effect removal by residualization.
"""
function remove_batch_effect(data::AbstractMatrix{<:Real}, batch;
    bio_design=nothing)

    Y = Matrix{Float64}(data)'
    ns = size(Y, 1)

    X_keep = if bio_design === nothing
        ones(Float64, ns, 1)
    else
        Xb = _ensure_sample_design(bio_design, ns)
        _ensure_intercept(Xb)
    end

    X_batch = if batch isa AbstractVector
        labels = String.(batch)
        length(labels) == ns || throw(ArgumentError("batch length must equal sample count"))
        _one_hot_labels(labels)[1]
    else
        _ensure_sample_design(batch, ns)
    end

    if size(X_batch, 2) == 0
        return Matrix{Float64}(data)
    end

    X_full = hcat(X_keep, X_batch)
    coef = X_full \ Y
    kkeep = size(X_keep, 2)
    batch_coef = coef[(kkeep+1):end, :]

    batch_effect = X_batch * batch_coef
    batch_center = mean(batch_effect, dims=1)
    corrected = Y - (batch_effect .- batch_center)
    return corrected'
end

function _combat_inv_gamma_moments(delta_hat::AbstractVector{<:Real}; eps::Real=1e-8)
    d = Float64.(delta_hat)
    d = d[isfinite.(d).&(d.>eps)]
    isempty(d) && return (100.0, 99.0)

    m = mean(d)
    v = var(d; corrected=true)
    if !isfinite(v) || v <= eps
        a = 100.0
    else
        a = 2.0 + (m^2 / v)
    end
    a = max(a, 2.001)
    b = max(m * (a - 1.0), eps)
    return (a, b)
end

function _combat_posterior_gamma_delta(x::AbstractVector{<:Real},
    gamma_hat::Float64,
    delta_hat::Float64,
    gamma_bar::Float64,
    t2::Float64,
    a_prior::Float64,
    b_prior::Float64;
    max_iter::Int=100,
    tol::Real=1e-6,
    eps::Real=1e-8)

    n = length(x)
    n == 0 && return (gamma_hat, max(delta_hat, eps))

    t2f = max(t2, eps)
    g = gamma_hat
    d = max(delta_hat, eps)
    xf = Float64.(x)

    for _ in 1:max_iter
        g_new = (t2f * n * gamma_hat + d * gamma_bar) / (t2f * n + d)
        sse = sum((xf .- g_new) .^ 2)
        d_new = (0.5 * sse + b_prior) / (0.5 * n + a_prior - 1.0)
        d_new = max(d_new, eps)

        if abs(g_new - g) < tol && abs(d_new - d) < tol
            g, d = g_new, d_new
            break
        end
        g, d = g_new, d_new
    end

    return g, max(d, eps)
end

"""
    combat_correction(data, batch; bio_design=nothing, ...)

ComBat empirical Bayes batch correction (Johnson et al. 2007).
"""
function combat_correction(data::AbstractMatrix{<:Real}, batch;
    bio_design=nothing,
    eps::Real=1e-8,
    parametric::Bool=true,
    max_iter::Int=100,
    tol::Real=1e-6)

    X = Matrix{Float64}(data)
    ng, ns = size(X)
    ns == 0 && return X

    if !(batch isa AbstractVector)
        return remove_batch_effect(X, batch; bio_design=bio_design)
    end

    labels = String.(batch)
    length(labels) == ns || throw(ArgumentError("batch length must equal sample count"))
    levels = unique(labels)
    nbatch = length(levels)
    nbatch <= 1 && return copy(X)

    batch_idx = [findall(labels .== lev) for lev in levels]
    B = zeros(Float64, ns, nbatch)
    for (b, idx) in enumerate(batch_idx)
        B[idx, b] .= 1.0
    end

    Xcov = if bio_design === nothing
        zeros(Float64, ns, 0)
    else
        Xb = _ensure_sample_design(bio_design, ns)
        if size(Xb, 2) > 0 && all(isapprox.(Xb[:, 1], 1.0; atol=1e-12, rtol=0.0))
            Xb[:, 2:end]
        else
            Xb
        end
    end

    design = hcat(B, Xcov)
    Y = X'

    coef = try
        design \ Y
    catch
        pinv(design) * Y
    end

    batch_coef = coef[1:nbatch, :]
    cov_coef = size(Xcov, 2) == 0 ? zeros(Float64, 0, ng) : coef[(nbatch+1):end, :]

    batch_props = Float64[length(idx) / ns for idx in batch_idx]
    grand_mean = vec((batch_props' * batch_coef))
    cov_effect = size(Xcov, 2) == 0 ? zeros(Float64, ns, ng) : Xcov * cov_coef
    stand_mean = cov_effect .+ reshape(grand_mean, 1, ng)

    fitted = design * coef
    resid = Y - fitted
    df = max(ns - rank(design), 1)
    var_pooled = vec(sum(resid .^ 2, dims=1) ./ df)
    var_pooled = max.(var_pooled, eps)

    s_data = (Y .- stand_mean) ./ reshape(sqrt.(var_pooled), 1, ng)

    gamma_hat = zeros(Float64, nbatch, ng)
    delta_hat = ones(Float64, nbatch, ng)
    for b in 1:nbatch
        idx = batch_idx[b]
        sb = s_data[idx, :]
        gamma_hat[b, :] .= vec(mean(sb, dims=1))
        if length(idx) > 1
            delta_hat[b, :] .= vec(var(sb, dims=1; corrected=true))
        else
            delta_hat[b, :] .= 1.0
        end
        delta_hat[b, :] .= max.(delta_hat[b, :], eps)
    end

    gamma_star = copy(gamma_hat)
    delta_star = copy(delta_hat)

    if parametric
        gamma_bar = vec(mean(gamma_hat, dims=2))
        t2 = vec(var(gamma_hat, dims=2; corrected=true))
        t2 = max.(t2, eps)

        for b in 1:nbatch
            idx = batch_idx[b]
            n_b = length(idx)
            n_b == 0 && continue

            a_prior, b_prior = _combat_inv_gamma_moments(view(delta_hat, b, :); eps=eps)
            for g in 1:ng
                g_star, d_star = _combat_posterior_gamma_delta(
                    view(s_data, idx, g),
                    gamma_hat[b, g],
                    delta_hat[b, g],
                    gamma_bar[b],
                    t2[b],
                    a_prior,
                    b_prior;
                    max_iter=max_iter,
                    tol=tol,
                    eps=eps)
                gamma_star[b, g] = g_star
                delta_star[b, g] = d_star
            end
        end
    end

    s_adj = copy(s_data)
    for b in 1:nbatch
        idx = batch_idx[b]
        isempty(idx) && continue
        s_adj[idx, :] .= (s_adj[idx, :] .- reshape(gamma_star[b, :], 1, ng)) ./
                         reshape(sqrt.(max.(delta_star[b, :], eps)), 1, ng)
    end

    corrected = s_adj .* reshape(sqrt.(var_pooled), 1, ng) .+ stand_mean
    return corrected'
end

# =============================================================================
# Surrogate Variables
# =============================================================================

"""
    estimate_surrogates(data; bio_design=nothing, nperm=10, ...)

SVA-style surrogate variable estimation via permutation testing.
"""
function estimate_surrogates(data::AbstractMatrix{<:Real};
    bio_design=nothing,
    nperm::Int=10,
    alpha::Real=0.1,
    max_components::Int=5,
    rng=Random.default_rng())

    Y = Matrix{Float64}(data)'
    ns = size(Y, 1)

    if bio_design !== nothing
        Xbio = _ensure_intercept(_ensure_sample_design(bio_design, ns))
        Y = Y - Xbio * (Xbio \ Y)
    else
        Y = Y .- mean(Y, dims=1)
    end

    centered = Y .- mean(Y, dims=1)
    F = svd(centered; full=false)
    kmax = min(clamp(max_components, 1, max(1, size(F.U, 2))), length(F.S), max(0, ns - 1))
    kmax == 0 && return zeros(Float64, ns, 0)
    if nperm <= 0
        return Matrix{Float64}(F.U[:, 1:kmax])
    end

    null_svals = zeros(Float64, kmax, nperm)
    for p in 1:nperm
        permuted = hcat([centered[:, g][randperm(rng, ns)] for g in axes(centered, 2)]...)
        null_svals[:, p] .= svd(permuted; full=false).S[1:kmax]
    end

    significant = Int[]
    for k in 1:kmax
        p_emp = (1.0 + count(null_svals[k, :] .>= F.S[k])) / (nperm + 1.0)
        p_emp <= alpha && push!(significant, k)
    end

    isempty(significant) && return zeros(Float64, ns, 0)
    return Matrix{Float64}(F.U[:, significant])
end

"""
    sv_analysis(cm, bio_design; method=:irw, n_sv=nothing)

Iteratively reweighted SVA (Leek 2014).
More accurate than two-step SVA by incorporating uncertainty in null gene identification.
"""
function sv_analysis(cm::CountMatrix, bio_design;
    method::Symbol=:irw,
    n_sv=nothing,
    max_iter::Int=10,
    tol::Real=1e-4,
    rng=Random.default_rng())

    # VST transform
    nf = calc_norm_factors(cm; method=:tmm)
    disp = estimate_dispersions(cm, nf; workflow=:prior)
    vst_data = vst(cm, disp; norm_factors=nf)'

    ns, ng = size(vst_data)

    Xbio = if bio_design === nothing
        ones(Float64, ns, 1)
    else
        _ensure_intercept(_design_model_matrix(bio_design)[1])
    end

    # Residuals from biological model
    resid = vst_data - Xbio * (Xbio \ vst_data)

    # Initial SVD
    F = svd(resid; full=false)
    kmax = n_sv === nothing ? min(5, size(F.U, 2)) : min(n_sv, size(F.U, 2))

    # Iteratively reweighted SVD
    weights = ones(Float64, ng)
    prev_sv = nothing

    for iter in 1:max_iter
        # Weighted SVD
        weighted_resid = resid .* sqrt.(weights)'
        F = svd(weighted_resid; full=false)

        # Test each component
        significant = Int[]
        for k in 1:kmax
            sv = F.U[:, k]
            # Correlation with each gene
            cor_vals = [abs(cor(sv, resid[:, g])) for g in 1:ng]
            cor_vals = Float64[isfinite(c) ? c : 0.0 for c in cor_vals]

            # Genes uncorrelated with SV are likely null
            null_threshold = quantile(cor_vals, 0.5)
            null_genes = cor_vals .< null_threshold

            # Update weights: downweight null genes
            weights[null_genes] .= 0.5
            weights[.!null_genes] .= 1.0

            # Permutation test for this SV
            nperm = 20
            null_cor = Float64[]
            for _ in 1:nperm
                perm_resid = hcat([resid[:, g][randperm(rng, ns)] for g in 1:ng])
                F_null = svd(perm_resid * Diagonal(sqrt.(weights)); full=false)
                push!(null_cor, F_null.S[min(k, length(F_null.S))])
            end

            p_emp = (1.0 + count(null_cor .>= F.S[k])) / (nperm + 1.0)
            p_emp <= 0.1 && push!(significant, k)
        end

        isempty(significant) && break

        current_sv = F.S[significant]
        if prev_sv !== nothing && maximum(abs.(current_sv .- prev_sv) ./ max.(abs.(prev_sv), 1e-8)) < tol
            break
        end
        prev_sv = copy(current_sv)
    end

    # Final unweighted SVD to get clean SVs
    F = svd(resid; full=false)
    k_out = n_sv === nothing ? length(significant) : min(n_sv, length(significant))
    k_out = min(k_out, size(F.U, 2))

    return k_out > 0 ? Matrix{Float64}(F.U[:, 1:k_out]) : zeros(Float64, ns, 0)
end

# =============================================================================
# RUV Normalization (NEW)
# =============================================================================

"""
    ruv_normalize(cm; negative_controls, k=1, method=:ruv4)

Remove Unwanted Variation using negative control genes (Gagnon-Bartsch et al.).
"""
function ruv_normalize(cm::CountMatrix;
    negative_controls::AbstractVector{<:Integer},
    k::Int=1,
    method::Symbol=:ruv4)

    Y = Matrix{Float64}(cm.counts)
    ng, ns = size(Y)

    nc_idx = Int.(negative_controls)
    all(i -> 1 <= i <= ng, nc_idx) || throw(ArgumentError("negative control indices out of bounds"))

    # Log-transform with pseudocount
    logY = log.(Y .+ 1.0)

    if method == :ruv4
        # RUV-4: Factor analysis on control genes
        Y_ctrl = logY[nc_idx, :]

        # Center
        Y_ctrl .-= mean(Y_ctrl, dims=2)

        # SVD to get factors of unwanted variation
        F = svd(Y_ctrl; full=false)
        W = F.V[:, 1:min(k, size(F.V, 2))]  # ns × k factor matrix

        # Estimate gene-specific factors
        # Solve: logY_g = α_g + X_g * β + W * γ_g + ε_g
        # For simplicity, regress out W from all genes
        alpha = vec(mean(logY, dims=2))

        # Regress out factors
        logY_centered = logY .- alpha
        gamma = (logY_centered * W) / (transpose(W) * W)  # ng × k coefficients
        corrected_log = logY_centered - gamma * transpose(W)

        # Back-transform (approximate)
        corrected = exp.(corrected_log .+ alpha) .- 1.0
        corrected = max.(corrected, 0.0)

        return CountMatrix(sparse(Int.(round.(corrected))), cm.gene_ids, cm.sample_ids), W
    else
        throw(ArgumentError("unsupported RUV method: $method"))
    end
end

# =============================================================================
# ZINB Test (NEW)
# =============================================================================

"""
    zinb_test(cm, design; zero_inflation_model=:common)

Zero-inflated negative binomial test for single-cell data.
"""
function zinb_test(cm::CountMatrix, design;
    zero_inflation_model::Symbol=:common,
    max_iter::Int=100,
    tol::Real=1e-6)

    X, cnames, _ = _design_model_matrix(design)
    nf = calc_norm_factors(cm; method=:tmm)
    Y = Matrix{Float64}(cm.counts)
    ng, ns = size(Y)
    offset = log.(nf)

    groups = if design isa AbstractVector
        String.(design)
    elseif design isa DataFrame && :condition in names(design)
        String.(design[!, :condition])
    else
        fill("A", ns)
    end
    levs = unique(groups)
    length(levs) >= 2 || throw(ArgumentError("need at least 2 groups for ZINB test"))
    zero_inflation_model in (:common, :groupwise) ||
        throw(ArgumentError("zero_inflation_model must be :common or :groupwise"))

    results = Vector{DEResult}(undef, ng)

    for g in 1:ng
        y = Y[g, :]
        base_mean = mean(y ./ nf)

        if all(iszero, y)
            results[g] = DEResult(cm.gene_ids[g], 0.0, NaN, NaN, NaN, NaN, NaN, true, true)
            continue
        end

        # EM algorithm for ZINB
        # Y ~ π * δ(0) + (1-π) * NB(μ, α)

        # Initialize from marginal
        mu_init = max(mean(y ./ nf), 1e-8)
        alpha_init = max(_moment_dispersion(y, nf), 1e-4)
        pi_init = count(==(0), y) / ns

        mu = mu_init
        alpha = alpha_init
        pi_z = pi_init
        loglik_prev = -Inf

        for iter in 1:max_iter
            r = 1.0 / alpha

            # E-step: compute posterior probability of being from NB component
            w_nb = zeros(Float64, ns)
            for s in 1:ns
                if y[s] == 0
                    p_zero_nb = (r / (r + mu))^r
                    p_zero_zinb = pi_z + (1 - pi_z) * p_zero_nb
                    w_nb[s] = (1 - pi_z) * p_zero_nb / max(p_zero_zinb, 1e-300)
                else
                    w_nb[s] = 1.0
                end
            end

            # M-step: estimate mu, alpha, pi
            # Weighted MLE for mu (via GLM)
            w_sum = sum(w_nb)
            mu_new = max(dot(w_nb, y) / max(w_sum, eps(Float64)), 1e-8)

            # Estimate alpha by weighted NB likelihood maximization.
            var_w = 0.0
            mean_w = 0.0
            for s in 1:ns
                mean_w += w_nb[s] * y[s]
                var_w += w_nb[s] * (y[s] - mu_new)^2
            end
            mean_w /= w_sum
            var_w /= w_sum
            alpha_moment = max((var_w - mean_w) / max(mean_w^2, 1e-8), 1e-8)

            function neg_weighted_nb_loglik(log_alpha)
                alpha_try = clamp(exp(log_alpha), 1e-8, 1e3)
                r_try = 1.0 / alpha_try
                ll = 0.0
                for s in 1:ns
                    ll += w_nb[s] * (
                        lgamma(y[s] + r_try) - lgamma(r_try) - lgamma(y[s] + 1.0) +
                        r_try * (log(r_try) - log(r_try + mu_new)) +
                        y[s] * (log(mu_new) - log(r_try + mu_new))
                    )
                end
                return -ll
            end

            alpha_new = alpha_moment
            try
                opt = Optim.optimize(neg_weighted_nb_loglik, log(1e-8), log(1e3), Optim.Brent(); iterations=100)
                if Optim.converged(opt)
                    alpha_new = clamp(exp(Optim.minimizer(opt)), 1e-8, 1e3)
                end
            catch
            end

            # Update pi
            pi_new = 0.0
            n_zero = count(==(0), y)
            for s in 1:ns
                if y[s] == 0
                    pi_new += 1.0 - w_nb[s]
                end
            end
            pi_new = n_zero == 0 ? 0.0 : clamp(pi_new / n_zero, 0.0, 0.99)

            # Check convergence
            loglik = 0.0
            r_new = 1.0 / alpha_new
            for s in 1:ns
                if y[s] == 0
                    p_zero_nb = (r_new / (r_new + mu_new))^r_new
                    loglik += log(max(pi_new + (1 - pi_new) * p_zero_nb, 1e-300))
                else
                    loglik += log(1 - pi_new) + lgamma(y[s] + r_new) - lgamma(r_new) - lgamma(y[s] + 1.0) +
                              r_new * (log(r_new) - log(r_new + mu_new)) + y[s] * (log(mu_new) - log(r_new + mu_new))
                end
            end

            if abs(loglik - loglik_prev) < tol
                mu, alpha, pi_z = mu_new, alpha_new, pi_new
                break
            end

            mu, alpha, pi_z = mu_new, alpha_new, pi_new
            loglik_prev = loglik
        end

        # Wald test on group difference
        logfc = 0.0
        se = 0.0
        pvalue = 1.0

        if length(levs) >= 2
            g1_idx = findall(==(levs[1]), groups)
            g2_idx = findall(==(levs[2]), groups)
            mu1 = mean(y[g1_idx] ./ nf[g1_idx])
            mu2 = mean(y[g2_idx] ./ nf[g2_idx])
            logfc = log2(max(mu2, 1e-8) / max(mu1, 1e-8))

            # Delta-method SE with ZINB group-mean variance.
            n1, n2 = length(g1_idx), length(g2_idx)
            pi1 = if zero_inflation_model == :groupwise
                clamp(count(==(0), y[g1_idx]) / max(n1, 1), 0.0, 0.99)
            else
                pi_z
            end
            pi2 = if zero_inflation_model == :groupwise
                clamp(count(==(0), y[g2_idx]) / max(n2, 1), 0.0, 0.99)
            else
                pi_z
            end
            mu_nb1 = mu1 / max(1.0 - pi1, 1e-6)
            mu_nb2 = mu2 / max(1.0 - pi2, 1e-6)
            var1 = (1.0 - pi1) * mu_nb1 * (1.0 + alpha * mu_nb1) + pi1 * (1.0 - pi1) * mu_nb1^2
            var2 = (1.0 - pi2) * mu_nb2 * (1.0 + alpha * mu_nb2) + pi2 * (1.0 - pi2) * mu_nb2^2
            se_logfc = sqrt(max(var1 / (max(n1, 1) * max(mu1, 1e-8)^2 * log(2)^2), 0) +
                            max(var2 / (max(n2, 1) * max(mu2, 1e-8)^2 * log(2)^2), 0))
            se = max(se_logfc, 1e-8)
            z = abs(logfc) / se
            pvalue = 2.0 * ccdf(Normal(), z)
        end

        is_zi = pi_z > 0.5
        results[g] = DEResult(cm.gene_ids[g], base_mean, logfc, se, abs(logfc) / max(se, 1e-8), pvalue, NaN, is_zi, true)
    end

    padj = benjamini_hochberg([r.pvalue for r in results])
    for i in 1:ng
        results[i] = DEResult(results[i].gene_id, results[i].base_mean, results[i].log2_fold_change,
            results[i].lfc_se, results[i].stat, results[i].pvalue, padj[i],
            results[i].zero_inflated, results[i].converged)
    end

    return results
end

# =============================================================================
# Time Series DE (NEW)
# =============================================================================

"""
    time_series_de(cm, timepoints, condition; df_spline=4, test_type=:interaction)

Spline-based longitudinal differential expression.
Effect size is reported at the median time point.
"""
function time_series_de(cm::CountMatrix, timepoints, condition;
    df_spline::Int=4,
    test_type::Symbol=:interaction,
    maxit::Integer=100,
    beta_tol::Real=1e-8,
    minmu::Real=0.5)

    nf = calc_norm_factors(cm; method=:tmm)
    Y = Matrix{Float64}(cm.counts)
    ng, ns = size(Y)

    t = Float64.(timepoints)
    cond = String.(condition)
    levs = unique(cond)
    length(levs) == 2 || throw(ArgumentError("need exactly 2 conditions for time_series_de"))
    ref_cond = levs[1]
    treat_cond = levs[2]

    # Create spline basis for time
    t_sorted = sort(unique(t))
    knots = quantile(t_sorted, range(0.1, 0.9; length=max(df_spline - 2, 1)))

    # Truncated power spline basis (before orthogonalization)
    n_spline = df_spline
    B_time = zeros(Float64, ns, n_spline)
    for s in 1:ns
        ts = t[s]
        B_time[s, 1] = 1.0
        B_time[s, 2] = ts
        if n_spline >= 3
            for k in 3:n_spline
                knot_k = knots[k-2]
                B_time[s, k] = max(0, ts - knot_k)^3
            end
        end
    end

    # ═══════════════════════════════════════════════════════════════════
    # Store QR factorization for later basis evaluation
    # ═══════════════════════════════════════════════════════════════════
    spline_qr = nothing
    if n_spline >= 3
        B_high = B_time[:, 3:end]
        if size(B_high, 2) > 0
            spline_qr = qr(B_high)  # ← SAVE THIS
            Q = Matrix(spline_qr.Q)[:, 1:size(B_high, 2)]
            B_time[:, 3:end] .= Q
        end
    end

    # Condition indicator
    is_treat = Float64.(cond .== treat_cond)

    # Full model: B_time + is_treat + B_time * is_treat
    # Coefficients: [β_time (n_spline) | β_cond (1) | β_int (n_spline)]
    X_full = hcat(B_time, is_treat, B_time .* is_treat)

    # Reduced model: B_time + is_treat (no interaction)
    X_reduced = hcat(B_time, is_treat)

    rank_full = rank(X_full)
    rank_reduced = rank(X_reduced)
    rank_full <= rank_reduced && throw(ArgumentError("full model must have higher rank than reduced"))

    offset = log.(nf)
    disp = estimate_dispersions(cm, nf; workflow=:prior, model_matrix=X_full)

    # ═══════════════════════════════════════════════════════════════════
    # Precompute spline basis at t_mid with correct orthogonalization
    # ═══════════════════════════════════════════════════════════════════
    t_mid = median(t)
    b_mid = zeros(Float64, n_spline)
    b_mid[1] = 1.0
    b_mid[2] = t_mid

    if n_spline >= 3 && spline_qr !== nothing
        # Construct original truncated-power basis at t_mid (pre-orthogonalization)
        b_high = zeros(Float64, n_spline - 2)
        for k in 1:(n_spline-2)
            b_high[k] = max(0.0, t_mid - knots[k])^3
        end
        # Map to Q-column space: q_new = R^{-T} * b_high
        # (derivation: B_high = Q*R → γ = R⁻¹β_int → b_new'γ = (R^{-T}b_new)'β_int)
        b_mid[3:end] .= spline_qr.R' \ b_high
    end

    results = Vector{DEResult}(undef, ng)

    solver_full = GLMSolver(X_full)
    solver_reduced = GLMSolver(X_reduced)

    for g in 1:ng
        y = Y[g, :]
        base_mean = mean(y ./ nf)

        if all(iszero, y)
            results[g] = DEResult(cm.gene_ids[g], 0.0, NaN, NaN, NaN, NaN, NaN, false, true)
            continue
        end

        alpha = max(disp[g], 1e-8)
        ok_full = false
        ok_red = false

        try
            fit_gene_fast!(solver_full, y, offset, alpha; maxit=maxit, beta_tol=beta_tol, minmu=minmu)
            ok_full = solver_full.converged
        catch
        end

        try
            fit_gene_fast!(solver_reduced, y, offset, alpha; maxit=maxit, beta_tol=beta_tol, minmu=minmu)
            ok_red = solver_reduced.converged
        catch
        end

        if !ok_full || !ok_red
            results[g] = DEResult(cm.gene_ids[g], base_mean, NaN, NaN, NaN, NaN, NaN, false, false)
            continue
        end

        ll_full = _nb_loglik(y, solver_full.mu, alpha)
        ll_red = _nb_loglik(y, solver_reduced.mu, alpha)
        stat = 2.0 * max(ll_full - ll_red, 0.0)
        df = rank_full - rank_reduced
        pvalue = ccdf(Chisq(df), stat)

        # ═══════════════════════════════════════════════════════════════════
        # Effect at t_mid = β_cond + B(t_mid)' * β_int
        # ═══════════════════════════════════════════════════════════════════
        cond_idx = n_spline + 1                           # β_cond position
        int_idx = (n_spline+2):(2*n_spline+1)      # β_int positions
        effect_at_mid = solver_full.beta[cond_idx] + dot(b_mid, solver_full.beta[int_idx])
        logfc = effect_at_mid / log(2)

        results[g] = DEResult(cm.gene_ids[g], base_mean, logfc, NaN, stat, pvalue, NaN, false, true)
    end

    padj = benjamini_hochberg([r.pvalue for r in results])
    for i in 1:ng
        results[i] = DEResult(results[i].gene_id, results[i].base_mean, results[i].log2_fold_change,
            results[i].lfc_se, results[i].stat, results[i].pvalue, padj[i],
            results[i].zero_inflated, results[i].converged)
    end

    return results
end

# =============================================================================
# Mixed Model DE (NEW)
# =============================================================================

"""
    dge_mixedmodel(cm, fixed_design, random_design; optimizer=:laplace)

Mixed-effects NB model for repeated measures.
"""
function dge_mixedmodel(cm::CountMatrix, fixed_design, random_design;
    optimizer::Symbol=:laplace,
    maxit::Integer=100,
    tol::Real=1e-6)

    nf = calc_norm_factors(cm; method=:tmm)
    Y = Matrix{Float64}(cm.counts)
    ng, ns = size(Y)

    X, cnames, _ = _design_model_matrix(fixed_design)
    Z = _ensure_sample_design(random_design, ns)

    # If Z has intercept, remove it (random effects are centered)
    if size(Z, 2) > 0 && all(isapprox.(Z[:, 1], 1.0; atol=1e-12, rtol=0.0))
        Z = Z[:, 2:end]
    end

    n_random = size(Z, 2)
    offset = log.(nf)

    # Estimate dispersions from fixed-effects-only model
    disp = estimate_dispersions(cm, nf; workflow=:prior, model_matrix=X)

    results = Vector{DEResult}(undef, ng)
    solver = GLMSolver(X)

    for g in 1:ng
        y = Y[g, :]
        base_mean = mean(y ./ nf)

        if all(iszero, y)
            results[g] = DEResult(cm.gene_ids[g], 0.0, NaN, NaN, NaN, NaN, NaN, false, true)
            continue
        end

        alpha = max(disp[g], 1e-8)

        # Fit fixed effects only for initial values
        fit_gene_fast!(solver, y, offset, alpha; maxit=50, beta_tol=1e-6)
        beta_init = copy(solver.beta)

        # Laplace approximation for random effects
        if n_random == 0 || optimizer != :laplace
            # No random effects or unsupported optimizer
            logfc = beta_init[min(2, length(beta_init))] / log(2)
            se = 1.0
            stat = abs(logfc) / se
            pvalue = 2.0 * ccdf(Normal(), stat)
            results[g] = DEResult(cm.gene_ids[g], base_mean, logfc, se, stat, pvalue, NaN, false, solver.converged)
            continue
        end

        # Initialize random effects to zero
        u = zeros(Float64, n_random)
        sigma2 = 0.1

        loglik_prev = -Inf

        for iter in 1:maxit
            eta = X * beta_init .+ Z * u .+ offset
            mu = max.(exp.(eta), 1e-8)
            var_y = mu .+ alpha .* mu .^ 2
            w = mu .^ 2 ./ max.(var_y, 1e-8)

            # Working response
            z = eta .+ (y .- mu) ./ mu

            # # Update beta (fixed effects)
            # XwX = X' * Diagonal(w) * X
            # Xwz = X' * Diagonal(w) * z
            # beta_new = (XwX + 1e-6 * I) \ Xwz

            # # Update u (random effects) via Laplace mode finding
            ZwX = Z' * Diagonal(w) * X
            Zwz = Z' * Diagonal(w) * z
            # ZwZ = Z' * Diagonal(w) * Z

            # # Penalized system for u: (Z'WZ + I/σ²) u = Z'Wz - ZwX*β
            # penalty = Diagonal(fill(1.0 / max(sigma2, 1e-8), n_random))
            # u_new = (ZwZ + penalty) \ (Zwz - ZwX * beta_new)

            XwX = X' * Diagonal(w) * X

            # Subtract offset and current random effects from working response
            z_adj = z .- offset .- Z * u
            Xwz = X' * Diagonal(w) * z_adj
            beta_new = (XwX + 1e-6 * I) \ Xwz

            # Penalized system for u: (Z'WZ + I/σ²) u = Z'W(z - offset - X*β)
            penalty = Diagonal(fill(1.0 / max(sigma2, 1e-8), n_random))
            ZwZ = Z' * Diagonal(w) * Z

            # Subtract offset and newly updated fixed effects
            z_res = z .- offset .- X * beta_new
            Zwz_res = Z' * Diagonal(w) * z_res
            u_new = (ZwZ + penalty) \ Zwz_res

            # Update σ²
            sigma2_new = max(sum(u_new .^ 2) / n_random, 1e-8)

            # Compute log-likelihood
            eta_new = X * beta_new .+ Z * u_new .+ offset
            mu_new = max.(exp.(eta_new), 1e-8)
            ll = _nb_loglik(y, mu_new, alpha)

            # Laplace approximation correction
            H = ZwZ + penalty
            try
                ll += 0.5 * logabsdet(H)[1]
                ll -= 0.5 * n_random * log(sigma2_new)
            catch
            end

            if abs(ll - loglik_prev) < tol
                beta_init .= beta_new
                u .= u_new
                sigma2 = sigma2_new
                break
            end

            beta_init .= beta_new
            u .= u_new
            sigma2 = sigma2_new
            loglik_prev = ll
        end

        # Wald test
        eta_final = X * beta_init .+ Z * u .+ offset
        mu_final = max.(exp.(eta_final), 1e-8)
        var_final = mu_final .+ alpha .* mu_final .^ 2
        w = mu_final .^ 2 ./ max.(var_final, 1e-8)
        XwX = X' * Diagonal(w) * X
        cov_beta = _safe_inverse(XwX)
        j = min(2, length(beta_init))
        logfc = beta_init[j] / log(2)
        se = sqrt(max(cov_beta[j, j], 1e-8)) / log(2)
        stat = abs(logfc) / se
        pvalue = 2.0 * ccdf(Normal(), stat)

        results[g] = DEResult(cm.gene_ids[g], base_mean, logfc, se, stat, pvalue, NaN, false, true)
    end

    padj = benjamini_hochberg([r.pvalue for r in results])
    for i in 1:ng
        results[i] = DEResult(results[i].gene_id, results[i].base_mean, results[i].log2_fold_change,
            results[i].lfc_se, results[i].stat, results[i].pvalue, padj[i],
            results[i].zero_inflated, results[i].converged)
    end

    return results
end

# =============================================================================
# Power Analysis (NEW)
# =============================================================================

"""
    power_analysis(; n_genes, n_samples_per_group, effect_size, ...)

Pre-study power calculation for NB differential expression.
"""
function power_analysis(;
    n_genes::Int=10000,
    n_samples_per_group::Int=5,
    effect_size::Real=1.0,
    dispersion::Real=0.1,
    fdr_target::Real=0.05,
    pi0::Real=0.9,
    n_sim::Int=100,
    seed::Int=42,
    sfType::Symbol=:ratio)

    rng = MersenneTwister(seed)
    n_de = round(Int, n_genes * (1 - pi0))
    n_null = n_genes - n_de
    de_gene_ids = Set("gene_$(i)" for i in 1:n_de)

    power_est = 0.0
    fdr_est = 0.0
    n_discoveries = 0.0

    for sim in 1:n_sim
        # Simulate counts
        m = 2 * n_samples_per_group
        condition = repeat([:A, :B], inner=n_samples_per_group)
        X, _, _ = _design_model_matrix(condition)

        b0 = randn(rng, n_genes) .* 0.3 .+ 4.0
        b1 = zeros(Float64, n_genes)
        b1[1:n_de] .= randn(rng, n_de) .* 0.5 .+ effect_size

        counts_mat = zeros(Int, n_genes, m)
        for g in 1:n_genes
            eta = b0[g] .+ b1[g] .* X[:, 2]
            mu = exp.(eta)
            alpha = max(dispersion, 1e-4)
            r = 1.0 / alpha
            for s in 1:m
                p = r / (r + mu[s])
                counts_mat[g, s] = rand(rng, NegativeBinomial(r, p))
            end
        end

        # Run DE analysis
        cm = CountMatrix(counts_mat, ["gene_$i" for i in 1:n_genes], ["sample_$i" for i in 1:m])
        try
            res = differential_expression(cm, condition;
                sfType=sfType,
                min_total=0,
                fdr_method=:BH,
                alpha=fdr_target)

            # Evaluate
            n_sig = count(r -> r.padj < fdr_target && isfinite(r.padj), res)
            n_true_pos = sum(r.padj < fdr_target && isfinite(r.padj) && (r.gene_id in de_gene_ids) for r in res)
            n_false_pos = n_sig - n_true_pos

            power_est += n_true_pos / max(n_de, 1)
            fdr_est += n_false_pos / max(n_sig, 1)
            n_discoveries += n_sig
        catch
            continue
        end
    end

    return (
        power=power_est / n_sim,
        fdr=fdr_est / n_sim,
        expected_discoveries=n_discoveries / n_sim,
        n_sim_completed=n_sim,
        parameters=Dict(
            "n_genes" => n_genes,
            "n_samples_per_group" => n_samples_per_group,
            "effect_size" => effect_size,
            "dispersion" => dispersion,
            "pi0" => pi0,
            "fdr_target" => fdr_target)
    )
end

# =============================================================================
# Bootstrap LFC (NEW)
# =============================================================================

"""
    bootstrap_lfc(cm, design; n_boot=999, alpha=0.05, seed=1)

Parametric bootstrap confidence intervals for LFC.
"""
function bootstrap_lfc(cm::CountMatrix, design;
    n_boot::Int=999,
    alpha::Real=0.05,
    seed::Int=1,
    coef_idx::Int=2)

    rng = MersenneTwister(seed)
    nf = calc_norm_factors(cm; method=:tmm)
    X, cnames, _ = _design_model_matrix(design)
    offset = log.(nf)
    Y = Matrix{Float64}(cm.counts)
    ng, ns = size(Y)

    disp = estimate_dispersions(cm, nf; workflow=:prior, model_matrix=X)
    solver = GLMSolver(X)

    # Fit original model
    lfc_orig = zeros(Float64, ng)
    mu_orig = zeros(Float64, ng, ns)
    for g in 1:ng
        y = Y[g, :]
        if all(iszero, y)
            continue
        end
        beta, _, _ = fit_gene_fast!(solver, y, offset, disp[g]; maxit=100)
        lfc_orig[g] = beta[coef_idx] / log(2)
        mu_orig[g, :] .= solver.mu
    end

    # Bootstrap
    lfc_boot = zeros(Float64, ng, n_boot)

    for b in 1:n_boot
        for g in 1:ng
            if all(iszero, Y[g, :])
                continue
            end
            alpha = max(disp[g], 1e-8)
            r = 1.0 / alpha
            y_boot = zeros(Int, ns)
            for s in 1:ns
                p = r / (r + mu_orig[g, s])
                y_boot[s] = rand(rng, NegativeBinomial(r, p))
            end
            try
                beta, _, _ = fit_gene_fast!(solver, Float64.(y_boot), offset, disp[g]; maxit=100)
                lfc_boot[g, b] = beta[coef_idx] / log(2)
            catch
                lfc_boot[g, b] = NaN
            end
        end
    end

    # Compute confidence intervals
    ci_lo = zeros(Float64, ng)
    ci_hi = zeros(Float64, ng)
    se_boot = zeros(Float64, ng)

    for g in 1:ng
        valid = filter(isfinite, lfc_boot[g, :])
        if length(valid) > 10
            ci_lo[g] = quantile(valid, alpha / 2)
            ci_hi[g] = quantile(valid, 1 - alpha / 2)
            se_boot[g] = std(valid)
        else
            ci_lo[g] = NaN
            ci_hi[g] = NaN
            se_boot[g] = NaN
        end
    end

    return DataFrame(
        gene_id=cm.gene_ids,
        log2FC=lfc_orig,
        CI_lo=ci_lo,
        CI_hi=ci_hi,
        SE=se_boot)
end

# =============================================================================
# Confounded Contrast Test (NEW)
# =============================================================================

"""
    confounded_contrast_test(dds, contrast_matrix)

F-test for multiple contrasts jointly: H0: C * β = 0.
"""
function confounded_contrast_test(dds::DESeqDataSet, contrast_matrix::AbstractMatrix{<:Real};
    alpha::Real=0.1)

    C = Matrix{Float64}(contrast_matrix)
    dds.wald_beta_matrix === nothing && throw(ArgumentError("run nbinomWaldTest first"))
    dds.wald_covariances === nothing && throw(ArgumentError("covariances not available"))

    beta_mat = dds.wald_beta_matrix
    covs = dds.wald_covariances
    ng = size(beta_mat, 1)
    n_contrasts = size(C, 1)
    X = dds.model_matrix === nothing ? _design_model_matrix(dds.design)[1] : dds.model_matrix
    df_resid = max(size(X, 1) - size(X, 2), 1)

    size(C, 2) == size(beta_mat, 2) || throw(ArgumentError("contrast matrix columns must match coefficient count"))

    stat = zeros(Float64, ng)
    pvalue = fill(1.0, ng)

    for g in 1:ng
        beta = vec(beta_mat[g, :])
        cov = covs[g]

        if !all(isfinite, beta) || !all(isfinite, cov)
            stat[g] = NaN
            continue
        end

        # F-statistic: (Cβ)' (C Σ C')⁻¹ (Cβ) / n_contrasts
        Cbeta = C * beta
        CSc = C * cov * C'

        try
            CSc_inv = inv(Symmetric(CSc))
            F = (Cbeta' * CSc_inv * Cbeta) / n_contrasts
            stat[g] = max(F, 0.0)
            pvalue[g] = ccdf(FDist(n_contrasts, df_resid), stat[g])
        catch
            stat[g] = NaN
        end
    end

    padj = benjamini_hochberg(pvalue)

    return DataFrame(
        gene_id=dds.counts.gene_ids,
        F_stat=stat,
        pvalue=pvalue,
        padj=padj)
end

# =============================================================================
# Concordance Analysis (NEW)
# =============================================================================

"""
    concordance_analysis(results_list; gene_universe, alpha=0.05)

Cross-study reproducibility scoring via rank-rank hypergeometric overlap.
"""
function concordance_analysis(results_list::Vector{Vector{DEResult}};
    gene_universe::Vector{String},
    alpha::Real=0.05)

    n_studies = length(results_list)
    n_genes = length(gene_universe)
    gene_set = Set(gene_universe)

    # Extract significant genes per study
    sig_genes = Vector{Set{String}}(undef, n_studies)
    for i in 1:n_studies
        sig = Set{String}()
        for r in results_list[i]
            if r.padj < alpha && isfinite(r.padj) && r.gene_id in gene_set
                push!(sig, r.gene_id)
            end
        end
        sig_genes[i] = sig
    end

    # Pairwise overlap
    overlap_matrix = zeros(Float64, n_studies, n_studies)
    jaccard_matrix = zeros(Float64, n_studies, n_studies)

    for i in 1:n_studies
        for j in 1:n_studies
            if i == j
                overlap_matrix[i, j] = length(sig_genes[i])
                jaccard_matrix[i, j] = 1.0
                continue
            end

            intersection = length(intersect(sig_genes[i], sig_genes[j]))
            union_size = length(union(sig_genes[i], sig_genes[j]))

            overlap_matrix[i, j] = intersection
            jaccard_matrix[i, j] = union_size > 0 ? intersection / union_size : 0.0
        end
    end

    # π₁ estimation (proportion of true positives)
    # Using Storey's method on p-values pooled across studies
    all_pvals = Float64[]
    for res in results_list
        for r in res
            if isfinite(r.pvalue) && r.gene_id in gene_set
                push!(all_pvals, r.pvalue)
            end
        end
    end

    pi1_hat = if length(all_pvals) > 100
        1.0 - storey_pi0(all_pvals)  # π₀ estimate
    else
        NaN
    end

    # Consensus genes (significant in >= 50% of studies)
    gene_counts = Dict{String,Int}()
    for sig in sig_genes
        for g in sig
            gene_counts[g] = get(gene_counts, g, 0) + 1
        end
    end

    threshold = ceil(Int, n_studies / 2)
    consensus_genes = [g for (g, c) in gene_counts if c >= threshold]

    return (
        pairwise_overlap=overlap_matrix,
        jaccard_index=jaccard_matrix,
        pi1_estimate=pi1_hat,
        consensus_genes=consensus_genes,
        n_consensus=length(consensus_genes),
        sig_counts_per_study=[length(sig) for sig in sig_genes])
end

# =============================================================================
# Example Generator
# =============================================================================

"""
    makeExampleDESeqDataSet(; n=1000, m=12, betaSD=2.0, ...)

Generate simulated DESeqDataSet for testing.
"""
function makeExampleDESeqDataSet(; n::Int=1000, m::Int=12,
    betaSD::Real=2.0,
    interceptMean::Real=4.0,
    dispersionMean::Real=0.5,
    seed::Int=1)

    rng = MersenneTwister(seed)
    m > 1 || throw(ArgumentError("m must be at least 2"))

    condition = [i <= cld(m, 2) ? :A : :B for i in 1:m]
    X, _, _ = _design_model_matrix(condition)
    has_coef = size(X, 2) >= 2

    b0 = randn(rng, n) .* 0.3 .+ interceptMean
    b1 = randn(rng, n) .* betaSD
    disp = rand(rng, Gamma(2.0, dispersionMean / 2.0), n)

    counts_mat = zeros(Int, n, m)
    for g in 1:n
        eta = b0[g] .+ (has_coef ? b1[g] .* X[:, 2] : 0.0)
        mu = exp.(eta)
        alpha = max(disp[g], 1e-4)
        r = 1.0 / alpha
        for s in 1:m
            p = r / (r + mu[s])
            counts_mat[g, s] = rand(rng, NegativeBinomial(r, p))
        end
    end

    gene_ids = ["gene_$(i)" for i in 1:n]
    sample_ids = ["sample_$(i)" for i in 1:m]
    cm = CountMatrix(counts_mat, gene_ids, sample_ids)
    coldata = DataFrame(condition=condition)
    return DESeqDataSet(cm, coldata, condition)
end

end # module DifferentialExpression
