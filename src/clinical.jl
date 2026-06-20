# ==============================================================================
# clinical.jl — Clinical genomics and survival analysis
#
# Provides Kaplan-Meier estimation, Cox proportional hazards modeling,
# log-rank testing, MAF file parsing, TCGA data access, survival ROC,
# competing risks, and oncoprint visualization.
#
# References:
#   - Kaplan & Meier (1958) JASA 53(282):457-481 (KM estimator)
#   - Cox (1972) JRSS B 34(2):187-220 (proportional hazards)
# ==============================================================================

module Clinical

using SparseArrays
using DataFrames
using Statistics
using LinearAlgebra
using Random
using Distributions
using JSON
using Downloads
using Printf
using Plots

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_container_provenance!, register_provenance!
using ..DifferentialExpression: CountMatrix, benjamini_hochberg

export PatientCohort, KaplanMeierResult, CoxResult, CoxTermResult, MAFRecord, MAFSummary
export read_maf, summarize_maf, tcga_query, tcga_download_files, merge_tcga_count_files, tcga_ingest, kaplan_meier, logrank_test, cox_ph
export forest_plot, survival_roc, cif_curve, neural_cox, dose_response_curve, oncoprint
export pharmacogenomics_star_alleles, omop_visit_summary, synthpop_like_cohort, trial_suitability_scores
export cpic_metabolizer_phenotype, pharmgkb_like_recommendations, propensity_score_match

"""
    PatientCohort

Clinical cohort with matched genomics and patient identifiers.
"""
struct PatientCohort
    clinical::DataFrame
    genomics::Union{CountMatrix,SparseMatrixCSC{Int,Int},Matrix{Float32}}
    patient_ids::Vector{String}

    function PatientCohort(clinical::DataFrame, genomics::Union{CountMatrix,SparseMatrixCSC{Int,Int},Matrix{Float32}}, patient_ids::AbstractVector{<:String})
        ids = String.(patient_ids)
        :patient_id in Symbol.(names(clinical)) || throw(ArgumentError("clinical DataFrame must contain a patient_id column"))
        clinical_ids = String.(clinical.patient_id)
        clinical_ids == ids || throw(ArgumentError("patient_ids must match clinical patient_id order"))
        if genomics isa CountMatrix
            genomics.sample_ids == ids || throw(ArgumentError("CountMatrix sample IDs must match patient_ids"))
        elseif genomics isa SparseMatrixCSC{Int,Int} || genomics isa Matrix{Float32}
            size(genomics, 2) == length(ids) || throw(ArgumentError("genomics columns must match patient_ids"))
        end
        new(clinical, genomics, ids)
    end
end

"""
    KaplanMeierResult

Kaplan-Meier survival curve summary with event and censoring counts.
"""
struct KaplanMeierResult <: AbstractAnalysisResult
    time::Vector{Float64}
    survival::Vector{Float64}
    at_risk::Vector{Int}
    events::Vector{Int}
    censored::Vector{Int}
    censor_times::Vector{Float64}
    censor_survival::Vector{Float64}
    provenance::ResultProvenance
end

KaplanMeierResult(time, survival, at_risk, events, censored, censor_times, censor_survival) =
    KaplanMeierResult(time, survival, at_risk, events, censored, censor_times, censor_survival, provenance_record("KaplanMeierResult", "clinical"))

"""
    CoxTermResult

Single coefficient summary from a Cox proportional hazards model.
"""
struct CoxTermResult <: AbstractAnalysisResult
    term::String
    beta::Float64
    hazard_ratio::Float64
    standard_error::Float64
    z_score::Float64
    pvalue::Float64
    ci_lower::Float64
    ci_upper::Float64
    provenance::ResultProvenance
end

CoxTermResult(term, beta, hazard_ratio, se, z_score, pvalue, ci_lower, ci_upper) =
    CoxTermResult(term, beta, hazard_ratio, se, z_score, pvalue, ci_lower, ci_upper, provenance_record("CoxTermResult", "clinical"))

"""
    CoxResult

Cox proportional hazards fit result with baseline hazard and convergence state.
"""
struct CoxResult <: AbstractAnalysisResult
    terms::Vector{CoxTermResult}
    baseline_times::Vector{Float64}
    baseline_hazard::Vector{Float64}
    loglik::Float64
    iterations::Int
    converged::Bool
    provenance::ResultProvenance
end

CoxResult(terms, baseline_times, baseline_hazard, loglik, iterations, converged) =
    CoxResult(terms, baseline_times, baseline_hazard, loglik, iterations, converged, provenance_record("CoxResult", "clinical"))

"""
    MAFRecord

Single mutation annotation record from a MAF file.
"""
struct MAFRecord
    gene::String
    sample::String
    variant_classification::String
    variant_type::String
    chromosome::String
    start_position::Int
    end_position::Int
    reference_allele::String
    tumor_seq_allele2::String
end

"""
    MAFSummary

Aggregate mutation summary derived from a MAF table.
"""
struct MAFSummary
    total_mutations::Int
    per_sample::Dict{String,Int}
    per_gene::Dict{String,Int}
    variant_classes::Dict{String,Int}
end

struct ROCResult <: AbstractAnalysisResult
    time::Float64
    thresholds::Vector{Float64}
    tpr::Vector{Float64}
    fpr::Vector{Float64}
    auc::Float64
    provenance::ResultProvenance
end

ROCResult(time, thresholds, tpr, fpr, auc) =
    ROCResult(time, thresholds, tpr, fpr, auc, provenance_record("ROCResult", "clinical"))

struct CIFResult <: AbstractAnalysisResult
    time::Vector{Float64}
    cumulative_incidence::Dict{Int,Vector{Float64}}
    censoring::Vector{Float64}
    provenance::ResultProvenance
end

CIFResult(time, cumulative_incidence, censoring) =
    CIFResult(time, cumulative_incidence, censoring, provenance_record("CIFResult", "clinical"))

struct NeuralCoxResult <: AbstractAnalysisResult
    weights::Vector{Matrix{Float64}}
    biases::Vector{Vector{Float64}}
    risk_scores::Vector{Float64}
    loglik::Float64
    provenance::ResultProvenance
end

NeuralCoxResult(weights, biases, risk_scores, loglik) =
    NeuralCoxResult(weights, biases, risk_scores, loglik, provenance_record("NeuralCoxResult", "clinical"))

struct DoseResponseResult <: AbstractAnalysisResult
    concentrations::Vector{Float64}
    responses::Vector{Float64}
    emax::Float64
    ec50::Float64
    hill::Float64
    fitted::Vector{Float64}
    ic50::Float64
    provenance::ResultProvenance
end

DoseResponseResult(concentrations, responses, emax, ec50, hill, fitted, ic50) =
    DoseResponseResult(concentrations, responses, emax, ec50, hill, fitted, ic50, provenance_record("DoseResponseResult", "clinical"))

struct OncoprintResult <: AbstractAnalysisResult
    genes::Vector{String}
    samples::Vector{String}
    matrix::Matrix{Int}
    mutation_labels::Matrix{String}
    provenance::ResultProvenance
end

OncoprintResult(genes, samples, matrix, mutation_labels) =
    OncoprintResult(genes, samples, matrix, mutation_labels, provenance_record("OncoprintResult", "clinical"))

function _clamp_probability(value::Real)
    return clamp(isfinite(Float64(value)) ? Float64(value) : 1.0, eps(Float64), 1.0)
end

function _clinical_index(cohort::PatientCohort, patient_id::String)
    index = findfirst(==(String(patient_id)), cohort.patient_ids)
    index === nothing && throw(KeyError(String(patient_id)))
    return index
end

function _mask_to_bool(mask::AbstractVector)
    boolmask = Vector{Bool}(undef, length(mask))
    for (index, value) in enumerate(mask)
        if value isa Bool
            boolmask[index] = value
        elseif value isa Integer || value isa AbstractFloat
            boolmask[index] = !iszero(value)
        else
            throw(ArgumentError("mask must contain Bool or numeric values"))
        end
    end
    return boolmask
end

function _subset_genomics(genomics::CountMatrix, columns::Vector{Int})
    return CountMatrix(genomics.counts[:, columns], genomics.gene_ids, genomics.sample_ids[columns])
end

function _subset_genomics(genomics::SparseMatrixCSC{Int,Int}, columns::Vector{Int})
    return genomics[:, columns]
end

function _subset_genomics(genomics::Matrix{Float32}, columns::Vector{Int})
    return genomics[:, columns]
end

function _cohort_view(cohort::PatientCohort, columns::Vector{Int})
    clinical = cohort.clinical[columns, :]
    genomics = _subset_genomics(cohort.genomics, columns)
    ids = cohort.patient_ids[columns]
    return PatientCohort(clinical, genomics, ids)
end

function Base.getindex(cohort::PatientCohort, patient_id::String)
    index = _clinical_index(cohort, patient_id)
    clinical_row = cohort.clinical[index, :]
    genomics_column = cohort.genomics isa CountMatrix ? cohort.genomics.counts[:, index] : cohort.genomics[:, index]
    return clinical_row, genomics_column
end

function Base.getindex(cohort::PatientCohort, mask::AbstractVector{Bool})
    length(mask) == length(cohort.patient_ids) || throw(ArgumentError("mask must match patient count"))
    columns = findall(mask)
    return _cohort_view(cohort, columns)
end

function Base.getindex(cohort::PatientCohort, mask::AbstractVector)
    return cohort[_mask_to_bool(mask)]
end

function _sorted_event_data(time::AbstractVector{<:Real}, status::AbstractVector{<:Integer})
    length(time) == length(status) || throw(ArgumentError("time and status must have the same length"))
    order = sortperm(Float64.(time), rev=false)
    return Float64.(time)[order], Int.(status)[order]
end

function _validated_group_levels(groups)
    levels = unique(groups)
    length(levels) == 2 || throw(ArgumentError("groups must contain exactly two levels"))
    return levels
end

function _count_tie_events(time::Vector{Float64}, status::Vector{Int}, event_time::Float64)
    event_mask = (time .== event_time) .& (status .> 0)
    censor_mask = (time .== event_time) .& (status .== 0)
    return count(event_mask), count(censor_mask)
end

function _km_survival_at(result::KaplanMeierResult, query_time::Real)
    survival = 1.0
    for (index, event_time) in enumerate(result.time)
        event_time <= query_time || break
        survival = result.survival[index]
    end
    return survival
end

function _km_plot_data(result::KaplanMeierResult)
    x = [0.0]
    y = [1.0]
    current = 1.0
    for (event_time, survival) in zip(result.time, result.survival)
        push!(x, event_time)
        push!(y, current)
        push!(x, event_time)
        push!(y, survival)
        current = survival
    end
    if !isempty(result.time)
        push!(x, last(result.time))
        push!(y, last(result.survival))
    end
    return x, y
end

"""
    kaplan_meier(time, status)

Compute a Kaplan-Meier survival estimate from event times and censoring status.
"""
function kaplan_meier(time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    sorted_time, sorted_status = _sorted_event_data(time, status)
    unique_times = sort(unique(sorted_time[sorted_status .> 0]))
    survival = Float64[]
    at_risk = Int[]
    events = Int[]
    censored = Int[]
    censor_times = Float64[]
    censor_survival = Float64[]
    current_survival = 1.0

    for event_time in unique_times
        risk = count(t -> t >= event_time, sorted_time)
        event_count, censor_count = _count_tie_events(sorted_time, sorted_status, event_time)
        risk > 0 || continue
        current_survival *= (1 - event_count / risk)
        push!(survival, current_survival)
        push!(at_risk, risk)
        push!(events, event_count)
        push!(censored, censor_count)
        for _ in 1:censor_count
            push!(censor_times, event_time)
            push!(censor_survival, current_survival)
        end
    end

    result = KaplanMeierResult(unique_times, survival, at_risk, events, censored, censor_times, censor_survival)
    return provenance_result!(_ctx, result, "kaplan_meier"; parents=provenance_parent_ids(time, status), parameters=(n=length(time), event_times=length(unique_times)))
end

"""
    kaplan_meier_plot(result; title="Kaplan-Meier", xlabel="Time", ylabel="Survival probability", show_censors=true)

Plot a Kaplan-Meier survival curve.
"""
function kaplan_meier_plot(result::KaplanMeierResult; title::String="Kaplan-Meier", xlabel::String="Time", ylabel::String="Survival probability", show_censors::Bool=true, kwargs...)
    _ctx = active_provenance_context()

    x, y = _km_plot_data(result)
    plt = plot(x, y; seriestype=:steppost, linewidth=2.5, color=:black, title=title, xlabel=xlabel, ylabel=ylabel, ylim=(0, 1.05), legend=false, kwargs...)
    if show_censors && !isempty(result.censor_times)
        censor_y = [_km_survival_at(result, t) for t in result.censor_times]
        scatter!(plt, result.censor_times, censor_y; markershape=_KM_CENSOR_MARKER, markercolor=:black, markersize=6, label=nothing)
    end
    return provenance_result!(_ctx, plt, "kaplan_meier_plot"; parents=provenance_parent_ids(result), parameters=(title=title, show_censors=show_censors))
end

function _logrank_components(time::Vector{Float64}, status::Vector{Int}, groups::Vector{Int})
    unique_times = sort(unique(time[status .> 0]))
    observed = 0.0
    expected = 0.0
    variance = 0.0
    for event_time in unique_times
        at_risk = time .>= event_time
        events_now = (time .== event_time) .& (status .> 0)
        n = count(at_risk)
        d = count(events_now)
        n1 = count(at_risk .& (groups .== 1))
        d1 = count(events_now .& (groups .== 1))
        n > 0 || continue
        e1 = d * n1 / n
        v1 = n > 1 ? d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1 + eps()) : 0.0
        observed += d1
        expected += e1
        variance += v1
    end
    return observed, expected, variance
end

"""
    logrank_test(time, status, groups)

Compare two survival groups with a log-rank test.
"""
function logrank_test(time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}, groups::AbstractVector; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    sorted_time, sorted_status = _sorted_event_data(time, status)
    length(groups) == length(time) || throw(ArgumentError("groups must match time and status"))
    levels = _validated_group_levels(groups)
    group_map = Dict(levels[1] => 1, levels[2] => 2)
    reorder = sortperm(Float64.(time))
    reordered_groups = Int[group_map[groups[index]] for index in reorder]
    observed, expected, variance = _logrank_components(sorted_time, sorted_status, reordered_groups)
    statistic = variance > 0 ? (observed - expected)^2 / variance : 0.0
    pvalue = ccdf(Chisq(1), statistic)
    result = (statistic=statistic, pvalue=pvalue, observed=observed, expected=expected, variance=variance)
    return provenance_result!(_ctx, result, "logrank_test"; parents=provenance_parent_ids(time, status, groups), parameters=(n=length(time), group_count=length(levels), pvalue=pvalue))
end

function _encode_covariate(values::AbstractVector, term_name::String)
    any(ismissing, values) && throw(ArgumentError("covariate $(term_name) contains missing values"))
    if all(v -> v isa Real || v isa Bool, values)
        return reshape(Float64.(values), :, 1), [String(term_name)]
    end
    labels = String.(values)
    levels = sort(unique(labels))
    if length(levels) <= 1
        return zeros(Float64, length(labels), 0), String[]
    end
    matrix = zeros(Float64, length(labels), length(levels) - 1)
    names = String[]
    for (column_index, level) in enumerate(levels[2:end])
        matrix[:, column_index] = Float64.(labels .== level)
        push!(names, string(term_name, "=", level))
    end
    return matrix, names
end

function _cox_design_matrix(cohort::PatientCohort, covariates)
    if covariates isa AbstractVector{<:String}
        design = ones(Float64, nrow(cohort.clinical), 1)
        names = String[]
        for covariate in covariates
            encoded, encoded_names = _encode_covariate(cohort.clinical[!, Symbol(covariate)], covariate)
            design = hcat(design, encoded)
            append!(names, encoded_names)
        end
        return design, names
    elseif covariates isa AbstractMatrix{<:Real}
        return hcat(ones(Float64, size(covariates, 1)), Matrix{Float64}(covariates)), [string("x", index) for index in 1:size(covariates, 2)]
    else
        throw(ArgumentError("covariates must be a matrix or a vector of clinical column names"))
    end
end

function _parse_surv_formula(formula)
    formula isa Expr && formula.head == :call && formula.args[1] == :~ || throw(ArgumentError("formula must look like Surv(time, status) ~ covariates"))
    lhs = formula.args[2]
    rhs = formula.args[3]
    lhs isa Expr && lhs.head == :call && lhs.args[1] == :Surv && length(lhs.args) == 3 || throw(ArgumentError("formula must use Surv(time, status) on the left hand side"))
    time_name = Symbol(lhs.args[2])
    status_name = Symbol(lhs.args[3])
    covariate_terms = String[]
    function _collect_terms(expr)
        if expr isa Symbol
            expr == :1 && return
            push!(covariate_terms, String(expr))
        elseif expr isa Number
            expr == 1 && return
            throw(ArgumentError("unsupported numeric term in formula"))
        elseif expr isa Expr && expr.head == :call && expr.args[1] == :+
            for arg in expr.args[2:end]
                _collect_terms(arg)
            end
        else
            push!(covariate_terms, string(expr))
        end
    end
    _collect_terms(rhs)
    return time_name, status_name, covariate_terms
end

function _cox_order(time::Vector{Float64}, status::Vector{Int})
    order = sortperm(time; rev=true)
    return time[order], status[order], order
end

function _cox_loglik_gradient_hessian(beta::Vector{Float64}, X::Matrix{Float64}, time::Vector{Float64}, status::Vector{Int})
    eta = X * beta
    risk = exp.(clamp.(eta, -40.0, 40.0))
    p = size(X, 2)
    loglik = 0.0
    gradient = zeros(Float64, p)
    hessian = zeros(Float64, p, p)
    for i in eachindex(time)
        status[i] > 0 || continue
        risk_set = time .>= time[i]
        risk_indices = findall(risk_set)
        denom = sum(risk[risk_indices])
        denom <= 0 && continue
        weighted_x = zeros(Float64, p)
        for j in risk_indices
            weighted_x .+= risk[j] .* X[j, :]
        end
        weighted_x ./= denom
        loglik += dot(X[i, :], beta) - log(denom)
        gradient .+= X[i, :] .- weighted_x
        for j in risk_indices
            centered = X[j, :] .- weighted_x
            hessian .-= (risk[j] / denom) .* (centered * centered')
        end
    end
    return loglik, gradient, hessian
end

"""
    cox_ph(formula, cohort; max_iter=50, tol=1e-7)

Fit a Cox proportional hazards model from a formula and cohort.
"""
function cox_ph(formula, cohort::PatientCohort; max_iter::Int=50, tol::Real=1e-7)
    _ctx = active_provenance_context()

    time_name, status_name, covariate_names = _parse_surv_formula(formula)
    X, term_names = _cox_design_matrix(cohort, covariate_names)
    time = Float64.(cohort.clinical[!, time_name])
    status = Int.(cohort.clinical[!, status_name])
    time, status, order = _cox_order(time, status)
    X = X[order, :]
    beta = zeros(Float64, size(X, 2))
    converged = false
    loglik = -Inf
    iterations = 0

    for iteration in 1:max_iter
        current_loglik, gradient, hessian = _cox_loglik_gradient_hessian(beta, X, time, status)
        regularized_hessian = hessian - 1e-6 * I
        step = -(regularized_hessian \ gradient)
        beta_new = beta .+ step
        iterations = iteration
        if norm(beta_new - beta) <= tol * (1 + norm(beta))
            beta = beta_new
            loglik = current_loglik
            converged = true
            break
        end
        beta = beta_new
        loglik = current_loglik
    end

    _, _, hessian = _cox_loglik_gradient_hessian(beta, X, time, status)
    covariance = try
        inv(Matrix(-Symmetric(hessian) + 1e-6 * I))
    catch
        pinv(Matrix(-hessian + 1e-6 * I))
    end
    terms = CoxTermResult[]
    for index in 2:size(X, 2)
        se = sqrt(max(covariance[index, index], eps(Float64)))
        z = beta[index] / se
        pvalue = 2 * ccdf(Normal(), abs(z))
        hr = exp(beta[index])
        push!(terms, CoxTermResult(term_names[index - 1], beta[index], hr, se, z, pvalue, exp(beta[index] - 1.96 * se), exp(beta[index] + 1.96 * se)))
    end

    base_times = sort(unique(time[status .> 0]))
    baseline_hazard = Float64[]
    cumulative = 0.0
    for event_time in base_times
        risk_set = time .>= event_time
        events_now = count((time .== event_time) .& (status .> 0))
        risk_score = sum(exp.(clamp.(X[risk_set, :] * beta, -40.0, 40.0)))
        increment = risk_score > 0 ? events_now / risk_score : 0.0
        cumulative += increment
        push!(baseline_hazard, cumulative)
    end

    result = CoxResult(terms, base_times, baseline_hazard, loglik, iterations, converged)
    return provenance_result!(_ctx, result, "cox_ph"; parents=provenance_parent_ids(cohort), parameters=(max_iter=max_iter, tol=Float64(tol), term_count=length(terms), iterations=iterations, converged=converged))
end

"""
    forest_plot(cox_result)

Plot coefficient estimates and confidence intervals from a Cox model fit.
"""
function forest_plot(cox_result::CoxResult)
    _ctx = active_provenance_context()
    n = length(cox_result.terms)
    if n == 0
        plt = plot(title="Cox forest plot", legend=false)
        return provenance_result!(_ctx, plt, "forest_plot"; parents=provenance_parent_ids(cox_result), parameters=(term_count=0))
    end
    terms = reverse(cox_result.terms)
    labels = [term.term for term in terms]
    hazard_ratios = [term.hazard_ratio for term in terms]
    lower_errors = [max(term.hazard_ratio - term.ci_lower, eps(Float64)) for term in terms]
    upper_errors = [max(term.ci_upper - term.hazard_ratio, eps(Float64)) for term in terms]
    pvalues = [term.pvalue for term in terms]
    p1 = plot(hazard_ratios, 1:n; xerror=(lower_errors, upper_errors), seriestype=:scatter, marker=:circle, markersize=7, color=:black, legend=false, xscale=:log10, xlabel="Hazard ratio", yticks=(1:n, labels), title="Cox forest plot", ylim=(0.5, n + 0.5), framestyle=:box)
    vline!(p1, [1.0]; linestyle=:dash, color=:gray)
    p2 = plot(; xlim=(0, 1), ylim=(0.5, n + 0.5), framestyle=:none, grid=false, legend=false, xticks=false, yticks=false, title="Summary")
    for (index, term) in enumerate(terms)
        row = n - index + 1
        annotate!(p2, 0.02, row, Plots.text(term.term, 8, :black, :left))
        annotate!(p2, 0.40, row, Plots.text(@sprintf("HR %.2f", term.hazard_ratio), 8, :black, :left))
        annotate!(p2, 0.68, row, Plots.text(@sprintf("95%% CI %.2f-%.2f", term.ci_lower, term.ci_upper), 8, :black, :left))
        annotate!(p2, 0.98, row, Plots.text(@sprintf("p=%.3g", term.pvalue), 8, :black, :right))
    end
    plt = plot(p1, p2; layout=(1, 2), size=(1100, 420))
    return provenance_result!(_ctx, plt, "forest_plot"; parents=provenance_parent_ids(cox_result), parameters=(term_count=n))
end

"""
    survival_roc(time, status, marker, predict_time)

Compute a time-dependent ROC curve for survival prediction.
"""
function survival_roc(time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}, marker::AbstractVector{<:Real}, predict_time::Real)
    _ctx = active_provenance_context()

    length(time) == length(status) == length(marker) || throw(ArgumentError("time, status, and marker must have the same length"))
    timef = Float64.(time)
    statusi = Int.(status)
    markerf = Float64.(marker)
    event = (timef .<= Float64(predict_time)) .& (statusi .> 0)
    control = timef .> Float64(predict_time)
    usable = event .| control
    thresholds = sort(unique(markerf[usable]))
    isempty(thresholds) && (thresholds = sort(unique(markerf)))
    tpr = Float64[]
    fpr = Float64[]
    for threshold in thresholds
        predicted = markerf .>= threshold
        tp = count(predicted .& event)
        fp = count(predicted .& control)
        fn = count((.!predicted) .& event)
        tn = count((.!predicted) .& control)
        push!(tpr, tp + fn > 0 ? tp / (tp + fn) : 0.0)
        push!(fpr, fp + tn > 0 ? fp / (fp + tn) : 0.0)
    end
    order = sortperm(fpr)
    fpr_curve = vcat(0.0, fpr[order], 1.0)
    tpr_curve = vcat(0.0, tpr[order], 1.0)
    auc = sum((fpr_curve[2:end] .- fpr_curve[1:end-1]) .* (tpr_curve[2:end] .+ tpr_curve[1:end-1]) ./ 2)
    result = ROCResult(Float64(predict_time), thresholds, tpr, fpr, auc)
    return provenance_result!(_ctx, result, "survival_roc"; parents=provenance_parent_ids(time, status, marker), parameters=(n=length(time), predict_time=Float64(predict_time), auc=auc))
end

"""
    cif_curve(time, status, cause)

Compute cumulative incidence curves for competing risks data.
"""
function cif_curve(time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}, cause::AbstractVector{<:Integer})
    _ctx = active_provenance_context()

    length(time) == length(status) == length(cause) || throw(ArgumentError("time, status, and cause must have the same length"))
    timef = Float64.(time)
    statusi = Int.(status)
    causei = Int.(cause)
    event_times = sort(unique(timef[statusi .> 0]))
    causes = sort(unique(causei[(statusi .> 0) .& (causei .> 0)]))
    cumulative = Dict{Int,Vector{Float64}}(cause_id => Float64[] for cause_id in causes)
    current = Dict{Int,Float64}(cause_id => 0.0 for cause_id in causes)
    censoring = Float64[]
    survival = 1.0
    for event_time in event_times
        at_risk = timef .>= event_time
        n = count(at_risk)
        n == 0 && continue
        total_events = count((timef .== event_time) .& (statusi .> 0))
        event_by_cause = Dict(cause_id => count((timef .== event_time) .& (statusi .> 0) .& (causei .== cause_id)) for cause_id in causes)
        for cause_id in causes
            current[cause_id] += survival * (event_by_cause[cause_id] / n)
            push!(cumulative[cause_id], current[cause_id])
        end
        survival *= max(1 - total_events / n, 0.0)
        push!(censoring, survival)
    end
    result = CIFResult(event_times, cumulative, censoring)
    return provenance_result!(_ctx, result, "cif_curve"; parents=provenance_parent_ids(time, status, cause), parameters=(n=length(time), cause_count=length(causes), event_times=length(event_times)))
end

"""
    read_maf(path::String)

Read a MAF file into a vector of mutation records.
"""
function read_maf(path::String)
    _ctx = active_provenance_context()
    open(path, "r") do io
        header = String[]
        fieldmap = Dict{String,Int}()
        records = MAFRecord[]
        for line in eachline(io)
            stripped = strip(line)
            isempty(stripped) && continue
            startswith(line, "#") && continue
            if isempty(header)
                header = split(stripped, '\t')
                fieldmap = Dict(name => index for (index, name) in enumerate(header))
                required = ["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "Variant_Type", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
                all(haskey(fieldmap, key) for key in required) || throw(ArgumentError("MAF header is missing required fields"))
                continue
            end
            fields = split(stripped, '\t')
            length(fields) >= length(header) || throw(ArgumentError("MAF row has fewer fields than the header"))
            push!(records, MAFRecord(
                fields[fieldmap["Hugo_Symbol"]],
                fields[fieldmap["Tumor_Sample_Barcode"]],
                fields[fieldmap["Variant_Classification"]],
                fields[fieldmap["Variant_Type"]],
                fields[fieldmap["Chromosome"]],
                parse(Int, fields[fieldmap["Start_Position"]]),
                parse(Int, fields[fieldmap["End_Position"]]),
                fields[fieldmap["Reference_Allele"]],
                fields[fieldmap["Tumor_Seq_Allele2"]]))
        end
        return provenance_result!(_ctx, records, "read_maf"; parents=String[], parameters=(path=path, record_count=length(records)))
    end
end

"""
    summarize_maf(maf::AbstractVector{<:MAFRecord})

Summarize mutation counts per sample, gene, and variant class.
"""
function summarize_maf(maf::AbstractVector{<:MAFRecord})
    _ctx = active_provenance_context()

    per_sample = Dict{String,Int}()
    per_gene = Dict{String,Int}()
    variant_classes = Dict{String,Int}()
    for record in maf
        per_sample[record.sample] = get(per_sample, record.sample, 0) + 1
        per_gene[record.gene] = get(per_gene, record.gene, 0) + 1
        variant_classes[record.variant_classification] = get(variant_classes, record.variant_classification, 0) + 1
    end
    result = MAFSummary(length(maf), per_sample, per_gene, variant_classes)
    return provenance_result!(_ctx, result, "summarize_maf"; parents=provenance_parent_ids(maf), parameters=(mutation_count=length(maf), sample_count=length(per_sample), gene_count=length(per_gene)))
end

function _maf_to_matrix(maf::AbstractVector{<:MAFRecord})
    genes = sort(unique(record.gene for record in maf))
    samples = sort(unique(record.sample for record in maf))
    matrix = zeros(Int, length(genes), length(samples))
    gene_index = Dict(gene => index for (index, gene) in enumerate(genes))
    sample_index = Dict(sample => index for (index, sample) in enumerate(samples))
    labels = fill("", length(genes), length(samples))
    for record in maf
        row = gene_index[record.gene]
        col = sample_index[record.sample]
        matrix[row, col] += 1
        labels[row, col] = isempty(labels[row, col]) ? record.variant_classification : labels[row, col] * ";" * record.variant_classification
    end
    return genes, samples, matrix, labels
end

"""
    oncoprint(maf::AbstractVector{<:MAFRecord})

Convert MAF records into an oncoprint matrix representation.
"""
function oncoprint(maf::AbstractVector{<:MAFRecord})
    _ctx = active_provenance_context()
    genes, samples, matrix, labels = _maf_to_matrix(maf)
    result = OncoprintResult(genes, samples, matrix, labels)
    return provenance_result!(_ctx, result, "oncoprint"; parents=provenance_parent_ids(maf), parameters=(gene_count=length(genes), sample_count=length(samples), mutation_count=length(maf)))
end

function _mutation_level(labels::String)
    isempty(labels) && return 0
    first_label = split(labels, ';')[1]
    return first_label == "Missense_Mutation" ? 1 : first_label == "Nonsense_Mutation" ? 2 : first_label == "Frame_Shift_Del" ? 3 : first_label == "Frame_Shift_Ins" ? 4 : first_label == "Splice_Site" ? 5 : 6
end

function _oncoprint_plot(result::OncoprintResult; title::String="Oncoprint", kwargs...)
    if isempty(result.genes) || isempty(result.samples)
        return plot(title=title, legend=false)
    end
    codes = zeros(Int, size(result.matrix))
    for row in axes(result.matrix, 1), col in axes(result.matrix, 2)
        codes[row, col] = result.matrix[row, col] > 0 ? _mutation_level(result.mutation_labels[row, col]) : 0
    end
    heat = plot(codes; seriestype=:heatmap, colorbar=true, c=:magma, xlabel="Patients", ylabel="Genes", title=title, yticks=(1:length(result.genes), reverse(result.genes)), xticks=(1:length(result.samples), result.samples), yflip=true, kwargs...)
    for row in axes(result.mutation_labels, 1), col in axes(result.mutation_labels, 2)
        label = result.mutation_labels[row, col]
        isempty(label) && continue
        annotate!(heat, col, row, text(label, 6, :white, :center))
    end
    return heat
end

"""
    oncoprint_plot(result; kwargs...)

Plot an oncoprint from a prepared result object.
"""
function oncoprint_plot(result::OncoprintResult, kwargs...)
    _ctx = active_provenance_context()
    plt = _oncoprint_plot(result; kwargs...)
    return provenance_result!(_ctx, plt, "oncoprint_plot"; parents=provenance_parent_ids(result), parameters=(gene_count=length(result.genes), sample_count=length(result.samples)))
end

"""
    tcga_query(; project, data_type, base_url="https://api.gdc.cancer.gov", limit=50)

Query the GDC/TCGA API for matching files.
"""
function tcga_query(; project::String, data_type::String, base_url::String="https://api.gdc.cancer.gov", limit::Int=50)
    _ctx = active_provenance_context()
    filters = Dict("op" => "and", "content" => [Dict("op" => "in", "content" => Dict("field" => "cases.project.project_id", "value" => [project])), Dict("op" => "in", "content" => Dict("field" => "data_type", "value" => [data_type]))])
    endpoint = string(base_url, "/files")
    payload = Dict("filters" => filters, "format" => "JSON", "size" => limit)
    io = IOBuffer()
    JSON.print(io, payload)
    body = String(take!(io))
    response_path = tempname()
    try
        Downloads.request(endpoint; method="POST", headers=Dict("Content-Type" => "application/json"), input=body, output=response_path)
        parsed = JSON.parse(read(response_path, String))
        results = get(parsed, "data", Dict())
        hits = get(results, "hits", Any[])
        result = (project=project, data_type=data_type, hits=hits)
        return provenance_result!(_ctx, result, "tcga_query"; parents=String[], parameters=(project=project, data_type=data_type, base_url=base_url, limit=limit, hit_count=length(hits)))
    catch err
        result = (project=project, data_type=data_type, hits=Any[], error=string(err))
        return provenance_result!(_ctx, result, "tcga_query"; parents=String[], parameters=(project=project, data_type=data_type, base_url=base_url, limit=limit, hit_count=0, error=string(err)))
    finally
        isfile(response_path) && rm(response_path; force=true)
    end
end

function _tcga_hit_string(hit::AbstractDict, key::String, fallback::String="")
    value = get(hit, key, nothing)
    value === nothing && return fallback
    return String(value)
end

function _tcga_case_field(hit::AbstractDict, field::String, fallback::String="")
    cases = get(hit, "cases", Any[])
    isempty(cases) && return fallback
    first_case = first(cases)
    first_case isa AbstractDict || return fallback
    value = get(first_case, field, nothing)
    value === nothing && return fallback
    return String(value)
end

function _tcga_file_id(hit::AbstractDict)
    return _tcga_hit_string(hit, "file_id", "")
end

function _tcga_sample_id(hit::AbstractDict)
    sample_id = _tcga_case_field(hit, "submitter_id", "")
    isempty(sample_id) && (sample_id = _tcga_case_field(hit, "case_id", ""))
    isempty(sample_id) && (sample_id = _tcga_hit_string(hit, "file_name", ""))
    isempty(sample_id) && (sample_id = _tcga_file_id(hit))
    return sample_id
end

function _tcga_file_name(hit::AbstractDict)
    name = _tcga_hit_string(hit, "file_name", "")
    isempty(name) && (name = _tcga_file_id(hit))
    return name
end

"""
    tcga_download_files(hits; base_url="https://api.gdc.cancer.gov", download_dir=tempdir(), fetcher=Downloads.download, _ctx=nothing)

Download TCGA/GDC files for a set of query hits.
"""
function tcga_download_files(hits::AbstractVector; base_url::String="https://api.gdc.cancer.gov", download_dir::String=tempdir(), fetcher=Downloads.download, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    mkpath(download_dir)
    file_paths = String[]
    sample_ids = String[]
    for hit_any in hits
        hit = hit_any isa AbstractDict ? hit_any : continue
        file_id = _tcga_file_id(hit)
        isempty(file_id) && continue
        sample_id = _tcga_sample_id(hit)
        file_name = _tcga_file_name(hit)
        output_name = isempty(sample_id) ? file_name : string(sample_id, "_", file_name)
        output_path = joinpath(download_dir, output_name)
        fetcher(string(base_url, "/data/", file_id), output_path)
        push!(file_paths, output_path)
        push!(sample_ids, sample_id)
    end
    _ctx === nothing || register_provenance!(_ctx, "tcga_download_files"; parameters=(base_url=String(base_url), download_dir=String(download_dir), file_count=length(file_paths), sample_ids=sample_ids))
    return (file_paths=file_paths, sample_ids=sample_ids)
end

function _parse_tcga_count_file(path::String; has_header::Bool=true, gene_column::Int=1, count_column::Int=2)
    genes = String[]
    counts = Int[]
    open(path, "r") do io
        first_data_line = true
        for line in eachline(io)
            stripped = strip(line)
            isempty(stripped) && continue
            startswith(stripped, "#") && continue
            if has_header && first_data_line
                first_data_line = false
                continue
            end
            first_data_line = false
            fields = split(stripped, '\t')
            length(fields) >= max(gene_column, count_column) || continue
            gene = fields[gene_column]
            count_value = try
                parse(Int, fields[count_column])
            catch
                round(Int, parse(Float64, fields[count_column]))
            end
            push!(genes, gene)
            push!(counts, count_value)
        end
    end
    return genes, counts
end

"""
    merge_tcga_count_files(file_paths, sample_ids; has_header=true, gene_column=1, count_column=2)

Merge multiple TCGA count files into a sparse count matrix.
"""
function merge_tcga_count_files(file_paths::AbstractVector{<:String}, sample_ids::AbstractVector{<:String}; has_header::Bool=true, gene_column::Int=1, count_column::Int=2)
    _ctx = active_provenance_context()

    length(file_paths) == length(sample_ids) || throw(ArgumentError("file_paths and sample_ids must have the same length"))
    gene_order = String[]
    gene_index = Dict{String,Int}()
    parsed_files = Vector{Tuple{Vector{String},Vector{Int}}}(undef, length(file_paths))

    for (index, path) in enumerate(file_paths)
        genes, counts = _parse_tcga_count_file(path; has_header=has_header, gene_column=gene_column, count_column=count_column)
        parsed_files[index] = (genes, counts)
        for gene in genes
            haskey(gene_index, gene) && continue
            gene_index[gene] = length(gene_order) + 1
            push!(gene_order, gene)
        end
    end

    matrix = zeros(Int, length(gene_order), length(sample_ids))
    for (sample_index, (genes, counts)) in enumerate(parsed_files)
        for (gene, count_value) in zip(genes, counts)
            row = gene_index[gene]
            matrix[row, sample_index] += count_value
        end
    end

    result = CountMatrix(sparse(matrix), gene_order, String.(sample_ids))
    return provenance_result!(_ctx, result, "merge_tcga_count_files"; parents=provenance_parent_ids(file_paths), parameters=(file_count=length(file_paths), sample_count=length(sample_ids), gene_count=length(gene_order), has_header=has_header, gene_column=gene_column, count_column=count_column))
end

"""
    tcga_ingest(file_paths, sample_ids; has_header=true, gene_column=1, count_column=2)

Ingest TCGA count files directly from local paths.
"""
function tcga_ingest(file_paths::AbstractVector{<:String}, sample_ids::AbstractVector{<:String}; has_header::Bool=true, gene_column::Int=1, count_column::Int=2)
    _ctx = active_provenance_context()
    matrix = merge_tcga_count_files(file_paths, sample_ids; has_header=has_header, gene_column=gene_column, count_column=count_column)
    return provenance_result!(_ctx, matrix, "tcga_ingest"; parents=provenance_parent_ids(file_paths), parameters=(file_count=length(file_paths), sample_count=length(sample_ids), has_header=has_header, gene_column=gene_column, count_column=count_column))
end

"""
    tcga_ingest(; project, data_type, base_url="https://api.gdc.cancer.gov", limit=50, download_dir=tempdir(), has_header=true, gene_column=1, count_column=2, _ctx=nothing)

Query, download, and merge TCGA count files in one step.
"""
function tcga_ingest(; project::String, data_type::String, base_url::String="https://api.gdc.cancer.gov", limit::Int=50, download_dir::String=tempdir(), has_header::Bool=true, gene_column::Int=1, count_column::Int=2)
    _ctx = active_provenance_context()
    query = tcga_query(project=project, data_type=data_type, base_url=base_url, limit=limit)
    isempty(query.hits) && throw(ArgumentError("TCGA query returned no hits for project=$(project), data_type=$(data_type)"))
    downloads = tcga_download_files(query.hits; base_url=base_url, download_dir=download_dir)
    isempty(downloads.file_paths) && throw(ArgumentError("TCGA query returned downloadable hits, but no files were downloaded"))
    matrix = merge_tcga_count_files(downloads.file_paths, downloads.sample_ids; has_header=has_header, gene_column=gene_column, count_column=count_column)
    _ctx === nothing || register_container_provenance!(_ctx, matrix, "tcga_ingest"; parameters=(project=project, data_type=data_type, base_url=base_url, file_count=length(downloads.file_paths), sample_count=length(downloads.sample_ids)))
    return matrix
end

function _gradient_step!(weights::Vector{Matrix{Float64}}, biases::Vector{Vector{Float64}}, activations, deltas, learning_rate::Float64)
    for layer in eachindex(weights)
        weights[layer] .-= learning_rate .* (deltas[layer] * activations[layer]')
        biases[layer] .-= learning_rate .* vec(sum(deltas[layer], dims=2))
    end
end

function _forward_pass(X::Matrix{Float64}, weights::Vector{Matrix{Float64}}, biases::Vector{Vector{Float64}})
    activations = Vector{Matrix{Float64}}(undef, length(weights) + 1)
    activations[1] = X'
    for layer in eachindex(weights)
        z = weights[layer] * activations[layer] .+ biases[layer]
        activations[layer + 1] = layer == length(weights) ? z : tanh.(z)
    end
    return activations
end

function _cox_partial_loglik(risk::Vector{Float64}, time::Vector{Float64}, status::Vector{Int})
    loglik = 0.0
    for index in eachindex(time)
        status[index] > 0 || continue
        risk_set = time .>= time[index]
        loglik += risk[index] - log(sum(exp.(clamp.(risk[risk_set], -40.0, 40.0))))
    end
    return loglik
end

"""
    neural_cox(X, time, status; hidden_units=8, learning_rate=0.01, epochs=200, seed=42)

Fit a small neural-network Cox model for survival prediction.
"""
function neural_cox(X::AbstractMatrix{<:Real}, time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}; hidden_units::Int=8, learning_rate::Real=0.01, epochs::Int=200, seed::Int=42)
    _ctx = active_provenance_context()

    Xf = Matrix{Float64}(X)
    timef = Float64.(time)
    statusi = Int.(status)
    rng = MersenneTwister(seed)
    weights = [randn(rng, hidden_units, size(Xf, 2)) / sqrt(size(Xf, 2)), randn(rng, 1, hidden_units) / sqrt(hidden_units)]
    biases = [zeros(Float64, hidden_units), zeros(Float64, 1)]
    best_loglik = -Inf
    best_weights = deepcopy(weights)
    best_biases = deepcopy(biases)

    for _ in 1:epochs
        activations = _forward_pass(Xf, weights, biases)
        risk = vec(activations[end])
        loglik = _cox_partial_loglik(risk, timef, statusi)
        if loglik > best_loglik
            best_loglik = loglik
            best_weights = deepcopy(weights)
            best_biases = deepcopy(biases)
        end
        output_delta = zeros(Float64, 1, length(timef))
        for index in eachindex(timef)
            statusi[index] > 0 || continue
            risk_set = timef .>= timef[index]
            risk_indices = findall(risk_set)
            scores = exp.(risk[risk_indices])
            denom = sum(scores)
            denom <= 0 && continue
            probabilities = scores / denom
            output_delta[1, index] -= 1.0
            output_delta[1, risk_indices] .+= probabilities
        end
        hidden = activations[2]
        hidden_delta = (weights[end]' * output_delta) .* (1 .- hidden .^ 2)
        _gradient_step!(weights, biases, activations, [hidden_delta, output_delta], Float64(learning_rate))
    end

    final_activations = _forward_pass(Xf, best_weights, best_biases)
    risk_scores = vec(final_activations[end])
    result = NeuralCoxResult(best_weights, best_biases, risk_scores, best_loglik)
    return provenance_result!(_ctx, result, "neural_cox"; parents=provenance_parent_ids(X, time, status), parameters=(hidden_units=hidden_units, learning_rate=Float64(learning_rate), epochs=epochs, seed=seed, loglik=best_loglik))
end

"""
    dose_response_curve(drug, cell_lines, concentrations, responses)

Fit a Hill-style dose-response curve from concentration/response measurements.
"""
function dose_response_curve(drug::String, cell_lines::AbstractVector{<:String}, concentrations::AbstractVector{<:Real}, responses::AbstractVector{<:Real})
    _ctx = active_provenance_context()
    length(concentrations) == length(responses) || throw(ArgumentError("concentrations and responses must have the same length"))
    length(cell_lines) == length(concentrations) || throw(ArgumentError("cell_lines must match concentrations and responses"))
    x = Float64.(concentrations)
    y = Float64.(responses)
    all(isfinite, x) || throw(ArgumentError("concentrations must be finite"))
    any(x .<= 0) && throw(ArgumentError("concentrations must be positive"))
    if all(y .== first(y))
        fitted = fill(first(y), length(y))
        result = DoseResponseResult(x, y, first(y), median(x), 1.0, fitted, median(x))
        return provenance_result!(_ctx, result, "dose_response_curve"; parents=provenance_parent_ids(cell_lines, concentrations, responses), parameters=(drug=drug, n=length(x), flat=true))
    end
    emin = minimum(y)
    emax = maximum(y)
    log_candidates = exp.(range(log(minimum(x)), log(maximum(x)), length=25))
    ec50_candidates = unique(sort(vcat([median(x)], x, log_candidates)))
    hill_candidates = collect(range(0.5, 4.0, length=15))
    best_sse = Inf
    best_ec50 = median(x)
    best_hill = 1.0
    best_fitted = fill(mean(y), length(y))
    for ec50 in ec50_candidates
        for hill in hill_candidates
            preds = emin .+ (emax - emin) ./ (1 .+ (ec50 ./ x).^hill)
            sse = sum((y .- preds).^2)
            if sse < best_sse
                best_sse = sse
                best_ec50 = ec50
                best_hill = hill
                best_fitted = preds
            end
        end
    end
    result = DoseResponseResult(x, y, emax, best_ec50, best_hill, best_fitted, best_ec50)
    return provenance_result!(_ctx, result, "dose_response_curve"; parents=provenance_parent_ids(cell_lines, concentrations, responses), parameters=(drug=drug, n=length(x), ec50=best_ec50, hill=best_hill))
end

"""
    pharmacogenomics_star_alleles(variant_calls)

Map PGx variants to a compact star-allele assignment table.
"""
function pharmacogenomics_star_alleles(variant_calls::DataFrame; gene_col::Symbol=:gene, variant_col::Symbol=:variant, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    hasproperty(variant_calls, gene_col) || throw(ArgumentError("missing gene column"))
    hasproperty(variant_calls, variant_col) || throw(ArgumentError("missing variant column"))

    star_map = Dict(
        "CYP2D6" => Dict("rs1065852" => "*10", "rs3892097" => "*4", "rs16947" => "*2"),
        "CYP2C19" => Dict("rs4244285" => "*2", "rs12248560" => "*17"))

    out = DataFrame(gene=String[], variant=String[], star_allele=String[])
    for row in eachrow(variant_calls)
        gene = String(row[gene_col])
        var = String(row[variant_col])
        allele = haskey(star_map, gene) ? get(star_map[gene], var, "unknown") : "unknown"
        push!(out, (gene, var, allele))
    end
    return provenance_result!(_ctx, out, "pharmacogenomics_star_alleles"; parents=provenance_parent_ids(variant_calls), parameters=(row_count=nrow(out), gene_col=gene_col, variant_col=variant_col))
end

"""
    omop_visit_summary(person, visit_occurrence, condition_occurrence)

Summarize OMOP person-level visit and condition burden.
"""
function omop_visit_summary(person::DataFrame, visit_occurrence::DataFrame, condition_occurrence::DataFrame)
    hasproperty(person, :person_id) || throw(ArgumentError("person table must contain person_id"))
    hasproperty(visit_occurrence, :person_id) || throw(ArgumentError("visit table must contain person_id"))
    hasproperty(condition_occurrence, :person_id) || throw(ArgumentError("condition table must contain person_id"))
    _ctx = active_provenance_context()

    visits = combine(groupby(visit_occurrence, :person_id), nrow => :n_visits)
    conds = combine(groupby(condition_occurrence, :person_id), nrow => :n_conditions)
    out = leftjoin(person, visits, on=:person_id)
    out = leftjoin(out, conds, on=:person_id)
    out[!, :n_visits] = coalesce.(out.n_visits, 0)
    out[!, :n_conditions] = coalesce.(out.n_conditions, 0)
    return provenance_result!(_ctx, out, "omop_visit_summary"; parents=provenance_parent_ids(person, visit_occurrence, condition_occurrence), parameters=(person_count=nrow(person), visit_count=nrow(visit_occurrence), condition_count=nrow(condition_occurrence)))
end

"""
    synthpop_like_cohort(df; n=nrow(df))

Generate a synthetic cohort via bootstrap with mild Gaussian jitter for
continuous variables.
"""
function synthpop_like_cohort(df::DataFrame; n::Int=nrow(df), seed::Int=1)
    _ctx = active_provenance_context()

    n >= 1 || throw(ArgumentError("n must be positive"))
    rng = MersenneTwister(seed)
    idx = rand(rng, 1:nrow(df), n)
    out = copy(df[idx, :])

    for name in names(out)
        col = out[!, name]
        if eltype(col) <: Real
            σ = std(skipmissing(Float64.(col)))
            noise = randn(rng, length(col)) .* (isfinite(σ) ? 0.05 * σ : 0.0)
            out[!, name] = Float64.(col) .+ noise
        end
    end
    return provenance_result!(_ctx, out, "synthpop_like_cohort"; parents=provenance_parent_ids(df), parameters=(input_rows=nrow(df), output_rows=n, seed=seed))
end

"""
    trial_suitability_scores(clinical)

Compute simple trial suitability scores from ECOG, age, and biomarker columns.
"""
function trial_suitability_scores(clinical::DataFrame; age_col::Symbol=:age, ecog_col::Symbol=:ecog, biomarker_cols::Vector{Symbol}=Symbol[], prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    hasproperty(clinical, age_col) || throw(ArgumentError("clinical table missing age column"))
    hasproperty(clinical, ecog_col) || throw(ArgumentError("clinical table missing ecog column"))

    n = nrow(clinical)
    score = zeros(Float64, n)
    for i in 1:n
        age = Float64(clinical[i, age_col])
        ecog = Float64(clinical[i, ecog_col])
        s = 0.0
        s += clamp((75.0 - age) / 40.0, 0.0, 1.0)
        s += clamp((2.0 - ecog) / 2.0, 0.0, 1.0)
        for b in biomarker_cols
            hasproperty(clinical, b) || continue
            v = clinical[i, b]
            if v isa Bool
                s += v ? 0.6 : 0.0
            else
                s += clamp(Float64(v), 0.0, 1.0)
            end
        end
        score[i] = s / max(2 + length(biomarker_cols), 1)
    end

    out = copy(clinical)
    out[!, :trial_score] = score
    out[!, :eligible] = score .>= 0.5
    return provenance_result!(_ctx, out, "trial_suitability_scores"; parents=provenance_parent_ids(clinical), parameters=(row_count=nrow(out), age_col=age_col, ecog_col=ecog_col, biomarker_count=length(biomarker_cols)))
end

"""
    cpic_metabolizer_phenotype(diplotypes; gene="CYP2D6")

Assign CPIC-like metabolizer categories from star-allele diplotypes.
"""
function cpic_metabolizer_phenotype(diplotypes::AbstractVector{<:AbstractString}; gene::AbstractString="CYP2D6", prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    activity = Dict(
        "*1" => 1.0,
        "*2" => 1.0,
        "*4" => 0.0,
        "*5" => 0.0,
        "*10" => 0.25,
        "*17" => 0.5,
        "*41" => 0.5)

    out = DataFrame(gene=String[], diplotype=String[], activity_score=Float64[], phenotype=String[])
    for d in diplotypes
        tokens = split(replace(String(d), '+' => '/'), '/')
        score = sum(get(activity, strip(t), 0.5) for t in tokens)
        pheno = score > 2.0 ? "ultrarapid_metabolizer" : score >= 1.25 ? "normal_metabolizer" : score >= 0.25 ? "intermediate_metabolizer" : "poor_metabolizer"
        push!(out, (String(gene), String(d), score, pheno))
    end
    return provenance_result!(_ctx, out, "cpic_metabolizer_phenotype"; parents=provenance_parent_ids(diplotypes), parameters=(gene=String(gene), row_count=nrow(out)))
end

"""
    pharmgkb_like_recommendations(genotypes)

Generate coarse PharmGKB/CPIC-like treatment recommendations from PGx diplotypes.
"""
function pharmgkb_like_recommendations(genotypes::DataFrame; gene_col::Symbol=:gene, diplotype_col::Symbol=:diplotype, sample_col::Union{Nothing,Symbol}=nothing)
    _ctx = active_provenance_context()
    hasproperty(genotypes, gene_col) || throw(ArgumentError("missing gene column"))
    hasproperty(genotypes, diplotype_col) || throw(ArgumentError("missing diplotype column"))

    rec_map = Dict(
        "CYP2D6" => Dict(
            "poor_metabolizer" => "consider alternative therapy or lower dose",
            "intermediate_metabolizer" => "consider reduced dose",
            "normal_metabolizer" => "standard dosing",
            "ultrarapid_metabolizer" => "consider higher dose or alternate drug"),
        "CYP2C19" => Dict(
            "poor_metabolizer" => "avoid prodrug activation-dependent therapies",
            "intermediate_metabolizer" => "consider alternative antiplatelet",
            "normal_metabolizer" => "standard dosing",
            "ultrarapid_metabolizer" => "monitor efficacy and adjust dose"))
    drug_map = Dict(
        "CYP2D6" => "codeine/tamoxifen",
        "CYP2C19" => "clopidogrel/PPIs")

    out = DataFrame(sample_id=String[], gene=String[], diplotype=String[], phenotype=String[], drug_context=String[], recommendation=String[])
    for row in eachrow(genotypes)
        gene = String(row[gene_col])
        diplotype = String(row[diplotype_col])
        ph = cpic_metabolizer_phenotype([diplotype]; gene=gene, prov_ctx=nothing).phenotype[1]
        recommendation = haskey(rec_map, gene) ? get(rec_map[gene], ph, "manual curation required") : "manual curation required"
        sample_id = sample_col === nothing || !hasproperty(genotypes, sample_col) ? "sample_unknown" : String(row[sample_col])
        push!(out, (sample_id, gene, diplotype, ph, get(drug_map, gene, "general"), recommendation))
    end
    return provenance_result!(_ctx, out, "pharmgkb_like_recommendations"; parents=provenance_parent_ids(genotypes), parameters=(row_count=nrow(out), gene_col=gene_col, diplotype_col=diplotype_col, sample_col=sample_col === nothing ? "none" : String(sample_col)))
end

function _fit_logistic_irls(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}; max_iter::Int=100, tol::Real=1e-7)
    β = zeros(Float64, size(X, 2))
    for _ in 1:max_iter
        η = X * β
        η = clamp.(η, -25.0, 25.0)
        μ = 1.0 ./ (1.0 .+ exp.(-η))
        w = max.(μ .* (1 .- μ), 1e-6)
        H = X' * (X .* w) + 1e-6 * I
        g = X' * (y .- μ)
        Δ = H \ g
        β_new = β .+ Δ
        if norm(β_new - β) <= Float64(tol)
            β = β_new
            break
        end
        β = β_new
    end
    return β
end

"""
    propensity_score_match(clinical; treatment_col=:treatment, covariates=Symbol[], ratio=1)

Nearest-neighbor propensity score matching without replacement.
"""
function propensity_score_match(clinical::DataFrame; treatment_col::Symbol=:treatment, covariates::Vector{Symbol}=Symbol[], ratio::Int=1)
    _ctx = active_provenance_context()

    hasproperty(clinical, treatment_col) || throw(ArgumentError("clinical table missing treatment column"))
    ratio >= 1 || throw(ArgumentError("ratio must be >= 1"))

    if isempty(covariates)
        covariates = [Symbol(name) for name in names(clinical) if Symbol(name) != treatment_col && eltype(clinical[!, Symbol(name)]) <: Real]
    end
    isempty(covariates) && throw(ArgumentError("no numeric covariates available for propensity model"))

    Xraw = hcat([Float64.(clinical[!, c]) for c in covariates]...)
    for j in axes(Xraw, 2)
        μ = mean(@view Xraw[:, j])
        σ = std(@view Xraw[:, j])
        if isfinite(σ) && σ > 0
            Xraw[:, j] .= (Xraw[:, j] .- μ) ./ σ
        else
            Xraw[:, j] .= 0.0
        end
    end

    X = hcat(ones(Float64, nrow(clinical)), Xraw)
    y = Float64.(clinical[!, treatment_col] .!= 0)
    β = _fit_logistic_irls(X, y)
    ps = 1.0 ./ (1.0 .+ exp.(-clamp.(X * β, -25.0, 25.0)))

    treated = findall(==(1.0), y)
    controls = Set(findall(==(0.0), y))
    matches = DataFrame(treated_index=Int[], control_index=Int[], treated_ps=Float64[], control_ps=Float64[], abs_distance=Float64[])

    for ti in sort(treated; by=i -> ps[i], rev=true)
        isempty(controls) && break
        for _ in 1:ratio
            isempty(controls) && break
            best = 0
            best_dist = Inf
            for ci in controls
                d = abs(ps[ti] - ps[ci])
                if d < best_dist
                    best_dist = d
                    best = ci
                end
            end
            best == 0 && break
            delete!(controls, best)
            push!(matches, (ti, best, ps[ti], ps[best], best_dist))
        end
    end

    result = (matches=matches, propensity_score=ps, coefficients=β, covariates=vcat(:intercept, covariates))
    return provenance_result!(_ctx, result, "propensity_score_match"; parents=provenance_parent_ids(clinical), parameters=(row_count=nrow(clinical), treatment_col=treatment_col, covariate_count=length(covariates), ratio=ratio, match_count=nrow(matches)))
end

end # module Clinical
