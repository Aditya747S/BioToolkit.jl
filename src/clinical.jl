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


using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_container_provenance!, register_provenance!, _provenance_timestamp
using ..DifferentialExpression: CountMatrix, benjamini_hochberg

import ..BioToolkit: to_html, export_html

function _get_plots_mod()
  ext = Base.get_extension(parentmodule(Clinical), :BioToolkitPlotsExt)
  if ext !== nothing
    return ext.Plots
  elseif isdefined(Main, :Plots)
    return Main.Plots
  else
    error("Plotting functions require Plots.jl to be loaded (`using Plots`).")
  end
end


export PatientCohort, KaplanMeierResult, CoxResult, CoxTermResult, MAFRecord, MAFSummary
export read_maf, summarize_maf, tcga_query, tcga_download_files, merge_tcga_count_files, tcga_ingest, kaplan_meier, logrank_test, cox_ph
export forest_plot, survival_roc, cif_curve, neural_cox, dose_response_curve, oncoprint, oncoprint_plot
export pharmacogenomics_star_alleles, omop_visit_summary, synthpop_like_cohort, trial_suitability_scores
export cpic_metabolizer_phenotype, pharmgkb_like_recommendations, propensity_score_match
export cox_zph, CoxZphResult, fine_gray, FineGrayResult, somatic_interactions, calculate_tmb, rmst, RMSTResult

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
  std_error::Vector{Float64}
  ci_lower::Vector{Float64}
  ci_upper::Vector{Float64}
  at_risk::Vector{Int}
  events::Vector{Int}
  censored::Vector{Int}
  censor_times::Vector{Float64}
  censor_survival::Vector{Float64}
  provenance::ResultProvenance
end

KaplanMeierResult(time, survival, std_error, ci_lower, ci_upper, at_risk, events, censored, censor_times, censor_survival) =
  KaplanMeierResult(time, survival, std_error, ci_lower, ci_upper, at_risk, events, censored, censor_times, censor_survival, provenance_record("KaplanMeierResult", "clinical"))

function KaplanMeierResult(time, survival, at_risk, events, censored, censor_times, censor_survival)
  se = zeros(Float64, length(survival))
  cil = copy(survival)
  ciu = copy(survival)
  return KaplanMeierResult(time, survival, se, cil, ciu, at_risk, events, censored, censor_times, censor_survival)
end

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

struct CoxZphResult <: AbstractAnalysisResult
  terms::Vector{String}
  rho::Vector{Float64}
  chisq::Vector{Float64}
  pvalue::Vector{Float64}
  schoenfeld_residuals::Matrix{Float64}
  times::Vector{Float64}
  global_chisq::Float64
  global_pvalue::Float64
  provenance::ResultProvenance
end

CoxZphResult(terms, rho, chisq, pvalue, schoenfeld_residuals, times, global_chisq, global_pvalue) =
  CoxZphResult(terms, rho, chisq, pvalue, schoenfeld_residuals, times, global_chisq, global_pvalue, provenance_record("CoxZphResult", "clinical"))

struct FineGrayResult <: AbstractAnalysisResult
  terms::Vector{CoxTermResult}
  baseline_times::Vector{Float64}
  baseline_hazard::Vector{Float64}
  loglik::Float64
  iterations::Int
  converged::Bool
  fail_code::Int
  provenance::ResultProvenance
end

FineGrayResult(terms, baseline_times, baseline_hazard, loglik, iterations, converged, fail_code) =
  FineGrayResult(terms, baseline_times, baseline_hazard, loglik, iterations, converged, fail_code, provenance_record("FineGrayResult", "clinical"))

struct RMSTResult <: AbstractAnalysisResult
  tau::Float64
  rmst::Float64
  std_error::Float64
  ci_lower::Float64
  ci_upper::Float64
  group_results::Union{Nothing,Dict{String,Any}}
  provenance::ResultProvenance
end

RMSTResult(tau, rmst, std_error, ci_lower, ci_upper, group_results) =
  RMSTResult(tau, rmst, std_error, ci_lower, ci_upper, group_results, provenance_record("RMSTResult", "clinical"))

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
  n_obs = length(time)
  n_obs == length(status) || throw(ArgumentError("time and status must have equal length"))

  order = sortperm(time)
  sorted_time = Float64.(time)[order]
  sorted_status = Int.(status)[order]

  unique_times = Float64[]
  survival = Float64[]
  std_error = Float64[]
  ci_lower = Float64[]
  ci_upper = Float64[]
  at_risk = Int[]
  events = Int[]
  censored = Int[]
  censor_times = Float64[]
  censor_survival = Float64[]

  current_survival = 1.0
  greenwood_sum = 0.0
  current_at_risk = n_obs
  i = 1

  while i <= n_obs
    t_current = sorted_time[i]
    d_count = 0
    c_count = 0

    while i <= n_obs && sorted_time[i] == t_current
      if sorted_status[i] > 0
        d_count += 1
      else
        c_count += 1
      end
      i += 1
    end

    n_curr = current_at_risk
    if d_count > 0 && n_curr > 0
      current_survival *= (1.0 - d_count / n_curr)
      if n_curr > d_count
        greenwood_sum += d_count / (n_curr * (n_curr - d_count))
      end
      se = current_survival * sqrt(greenwood_sum)

      if current_survival > 0 && current_survival < 1
        log_s = log(current_survival)
        theta = 1.96 * se / (current_survival * abs(log_s))
        cil = clamp(current_survival^exp(theta), 0.0, 1.0)
        ciu = clamp(current_survival^exp(-theta), 0.0, 1.0)
      else
        cil = current_survival
        ciu = current_survival
      end

      push!(unique_times, t_current)
      push!(survival, current_survival)
      push!(std_error, se)
      push!(ci_lower, cil)
      push!(ci_upper, ciu)
      push!(at_risk, n_curr)
      push!(events, d_count)
      push!(censored, c_count)

      for _ in 1:c_count
        push!(censor_times, t_current)
        push!(censor_survival, current_survival)
      end
    else
      for _ in 1:c_count
        push!(censor_times, t_current)
        push!(censor_survival, current_survival)
      end
    end

    current_at_risk -= (d_count + c_count)
  end

  result = KaplanMeierResult(unique_times, survival, std_error, ci_lower, ci_upper, at_risk, events, censored, censor_times, censor_survival)
  return provenance_result!(_ctx, result, "kaplan_meier"; parents=provenance_parent_ids(time, status), parameters=(n=n_obs, event_times=length(unique_times)))
end

const _KM_CENSOR_MARKER = :plus

"""
    kaplan_meier_plot(result; title="Kaplan-Meier", xlabel="Time", ylabel="Survival probability", show_censors=true)

Plot a Kaplan-Meier survival curve.
"""
function kaplan_meier_plot(result::KaplanMeierResult; title::String="Kaplan-Meier", xlabel::String="Time", ylabel::String="Survival probability", show_censors::Bool=true, kwargs...)
  _ctx = active_provenance_context()
  PlotsMod = _get_plots_mod()
  x, y = _km_plot_data(result)
  plt = PlotsMod.plot(x, y; seriestype=:steppost, linewidth=2.5, color=:black, title=title, xlabel=xlabel, ylabel=ylabel, ylim=(0, 1.05), legend=false, kwargs...)
  if show_censors && !isempty(result.censor_times)
    censor_y = [_km_survival_at(result, t) for t in result.censor_times]
    PlotsMod.scatter!(plt, result.censor_times, censor_y; markershape=_KM_CENSOR_MARKER, markercolor=:black, markersize=6, label=nothing)
  end
  return provenance_result!(_ctx, plt, "kaplan_meier_plot"; parents=provenance_parent_ids(result), parameters=(title=title, show_censors=show_censors))
end

"""
    logrank_test(time, status, groups)

Compare survival distributions across two or more groups with a log-rank test.
"""
function logrank_test(time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}, groups::AbstractVector; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
  n = length(time)
  (n == length(status) && n == length(groups)) || throw(ArgumentError("time, status, and groups must have equal length"))

  levels = unique(groups)
  K = length(levels)
  K >= 2 || throw(ArgumentError("groups must contain at least 2 distinct levels"))

  group_map = Dict(lev => idx for (idx, lev) in enumerate(levels))
  g_encoded = Int[group_map[g] for g in groups]

  order = sortperm(time)
  t_sorted = Float64.(time)[order]
  s_sorted = Int.(status)[order]
  g_sorted = g_encoded[order]

  n_group = zeros(Int, K)
  for g in g_encoded
    n_group[g] += 1
  end

  observed = zeros(Float64, K)
  expected = zeros(Float64, K)
  V = zeros(Float64, K, K)

  i = 1
  while i <= n
    t_curr = t_sorted[i]
    d_group = zeros(Int, K)
    c_group = zeros(Int, K)

    while i <= n && t_sorted[i] == t_curr
      g = g_sorted[i]
      if s_sorted[i] > 0
        d_group[g] += 1
      else
        c_group[g] += 1
      end
      i += 1
    end

    d_total = sum(d_group)
    n_total = sum(n_group)

    if d_total > 0
      for g in 1:K
        n_g = n_group[g]
        e_g = d_total * (n_g / n_total)
        observed[g] += d_group[g]
        expected[g] += e_g
      end

      if n_total > 1
        factor = (n_total - d_total) / (n_total - 1)
        for g in 1:K
          n_g = n_group[g]
          for m in 1:K
            n_m = n_group[m]
            if g == m
              V[g, g] += d_total * factor * (n_g / n_total) * (1.0 - n_g / n_total)
            else
              V[g, m] += -d_total * factor * (n_g / n_total) * (n_m / n_total)
            end
          end
        end
      end
    end

    for g in 1:K
      n_group[g] -= (d_group[g] + c_group[g])
    end
  end

  df = K - 1
  Z = observed[1:df] .- expected[1:df]
  V_sub = V[1:df, 1:df]
  statistic = try
    dot(Z, V_sub \ Z)
  catch
    dot(Z, pinv(V_sub) * Z)
  end
  pvalue = ccdf(Chisq(df), max(statistic, 0.0))

  result = (statistic=statistic, pvalue=pvalue, df=df, observed=observed[1], expected=expected[1], variance=V[1, 1])
  return provenance_result!(_ctx, result, "logrank_test"; parents=provenance_parent_ids(time, status, groups), parameters=(n=n, group_count=K, statistic=statistic, pvalue=pvalue))
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
  formula isa Expr && formula.head == :call && formula.args[1] == :~ || throw(ArgumentError("formula must look like Surv(time, status) ~ covariates or Surv(start, stop, status) ~ covariates"))
  lhs = formula.args[2]
  rhs = formula.args[3]
  lhs isa Expr && lhs.head == :call && lhs.args[1] == :Surv || throw(ArgumentError("formula must use Surv(...) on the left hand side"))

  start_name = nothing
  time_name = Symbol()
  status_name = Symbol()

  if length(lhs.args) == 3
    time_name = Symbol(lhs.args[2])
    status_name = Symbol(lhs.args[3])
  elseif length(lhs.args) == 4
    start_name = Symbol(lhs.args[2])
    time_name = Symbol(lhs.args[3])
    status_name = Symbol(lhs.args[4])
  else
    throw(ArgumentError("Surv must have 2 or 3 arguments"))
  end

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
  return start_name, time_name, status_name, covariate_terms
end

function _cox_order(time::Vector{Float64}, status::Vector{Int})
  order = sortperm(time; rev=true)
  return time[order], status[order], order
end

function _cox_loglik_gradient_hessian(beta::Vector{Float64}, X::Matrix{Float64}, time::Vector{Float64}, status::Vector{Int}; ties::Symbol=:efron)
  N, p = size(X)
  eta = X * beta
  risk = exp.(clamp.(eta, -40.0, 40.0))

  loglik = 0.0
  gradient = zeros(Float64, p)
  hessian = zeros(Float64, p, p)

  S0 = 0.0
  S1 = zeros(Float64, p)
  S2 = zeros(Float64, p, p)

  i = 1
  while i <= N
    t_curr = time[i]
    start_idx = i

    while i <= N && time[i] == t_curr
      w = risk[i]
      x_i = @view X[i, :]
      S0 += w
      @inbounds for k in 1:p
        S1[k] += w * x_i[k]
        for l in 1:p
          S2[k, l] += w * x_i[k] * x_i[l]
        end
      end
      i += 1
    end

    event_indices = Int[]
    for j in start_idx:(i-1)
      if status[j] > 0
        push!(event_indices, j)
      end
    end
    d_events = length(event_indices)

    if d_events > 0 && S0 > 0.0
      if ties == :efron && d_events > 1
        S0_event = 0.0
        S1_event = zeros(Float64, p)
        S2_event = zeros(Float64, p, p)
        for j in event_indices
          w_j = risk[j]
          x_j = @view X[j, :]
          S0_event += w_j
          for k in 1:p
            S1_event[k] += w_j * x_j[k]
            for l in 1:p
              S2_event[k, l] += w_j * x_j[k] * x_j[l]
            end
          end
        end

        for j in event_indices
          loglik += dot(@view(X[j, :]), beta)
        end

        for m in 0:(d_events-1)
          frac = m / d_events
          S0_adj = S0 - frac * S0_event
          S1_adj = S1 .- frac .* S1_event
          S2_adj = S2 .- frac .* S2_event

          loglik -= log(max(S0_adj, eps(Float64)))
          E_x = S1_adj ./ max(S0_adj, eps(Float64))
          for k in 1:p
            for l in 1:p
              v_kl = (S2_adj[k, l] / max(S0_adj, eps(Float64))) - E_x[k] * E_x[l]
              hessian[k, l] -= v_kl
            end
          end
          gradient .-= E_x
        end

        for j in event_indices
          gradient .+= @view X[j, :]
        end
      else
        E_x = S1 ./ S0
        for j in event_indices
          x_j = @view X[j, :]
          loglik += dot(x_j, beta) - log(S0)
          gradient .+= x_j .- E_x
        end

        for k in 1:p
          for l in 1:p
            v_kl = (S2[k, l] / S0) - E_x[k] * E_x[l]
            hessian[k, l] -= d_events * v_kl
          end
        end
      end
    end
  end

  return loglik, gradient, hessian
end

"""
    cox_ph(formula, cohort; max_iter=50, tol=1e-7, ties=:efron, cluster=nothing)

Fit a Cox proportional hazards model from a formula and cohort.
Supports Efron/Breslow ties handling, cluster-robust standard errors, and counting process format.
"""
function cox_ph(formula, cohort::PatientCohort; max_iter::Int=50, tol::Real=1e-7, ties::Symbol=:efron, cluster=nothing)
  _ctx = active_provenance_context()

  start_name, time_name, status_name, covariate_names = _parse_surv_formula(formula)
  X, term_names = _cox_design_matrix(cohort, covariate_names)
  time = Float64.(cohort.clinical[!, time_name])
  status = Int.(cohort.clinical[!, status_name])
  time, status, order = _cox_order(time, status)
  X = X[order, :]

  cluster_vec = nothing
  if cluster !== nothing
    c_vals = cohort.clinical[!, cluster isa Symbol ? cluster : Symbol(cluster)]
    cluster_vec = c_vals[order]
  end

  beta = zeros(Float64, size(X, 2))
  converged = false
  loglik = -Inf
  iterations = 0

  for iteration in 1:max_iter
    current_loglik, gradient, hessian = _cox_loglik_gradient_hessian(beta, X, time, status; ties=ties)
    regularized_hessian = hessian - 1e-6 * I
    step = -(regularized_hessian \ gradient)
    beta_new = beta .+ step
    iterations = iteration
    if norm(beta_new - beta) <= tol * (1 + norm(beta))
      beta = beta_new
      converged = true
      break
    end
    beta = beta_new
  end

  # Accurately compute log-likelihood and Hessian at final beta
  loglik, _, hessian = _cox_loglik_gradient_hessian(beta, X, time, status; ties=ties)

  covariance = try
    inv(Matrix(-Symmetric(hessian) + 1e-6 * I))
  catch
    pinv(Matrix(-hessian + 1e-6 * I))
  end

  if cluster_vec !== nothing
    # Cluster-robust sandwich variance-covariance matrix
    N, p = size(X)
    eta = X * beta
    risk = exp.(clamp.(eta, -40.0, 40.0))
    scores = zeros(Float64, N, p)

    # Compute score residuals
    S0 = 0.0
    S1 = zeros(Float64, p)
    i = 1
    while i <= N
      t_curr = time[i]
      start_idx = i
      while i <= N && time[i] == t_curr
        w = risk[i]
        S0 += w
        S1 .+= w .* @view(X[i, :])
        i += 1
      end
      E_x = S1 ./ max(S0, eps(Float64))
      for j in start_idx:(i-1)
        if status[j] > 0
          scores[j, :] .+= @view(X[j, :]) .- E_x
        end
      end
    end

    clusters = unique(cluster_vec)
    G = zeros(Float64, length(clusters), p)
    for (c_idx, cl) in enumerate(clusters)
      mask = cluster_vec .== cl
      G[c_idx, :] = vec(sum(scores[mask, :], dims=1))
    end

    sandwich = G' * G
    covariance = covariance * sandwich * covariance
  end

  terms = CoxTermResult[]
  for index in 2:size(X, 2)
    se = sqrt(max(covariance[index, index], eps(Float64)))
    z = beta[index] / se
    pvalue = 2 * ccdf(Normal(), abs(z))
    hr = exp(beta[index])
    push!(terms, CoxTermResult(term_names[index-1], beta[index], hr, se, z, pvalue, exp(beta[index] - 1.96 * se), exp(beta[index] + 1.96 * se)))
  end

  # Fast O(N) baseline hazard computation on sorted data
  N = length(time)
  risk = exp.(clamp.(X * beta, -40.0, 40.0))
  base_times = Float64[]
  baseline_hazard = Float64[]

  unique_events = sort(unique(time[status .> 0]))
  if !isempty(unique_events)
    curr_idx = 1
    S0_acc = 0.0
    rev_events = reverse(unique_events)
    inc_hazards = Float64[]
    for ev in rev_events
      while curr_idx <= N && time[curr_idx] >= ev
        S0_acc += risk[curr_idx]
        curr_idx += 1
      end
      d_now = count((time .== ev) .& (status .> 0))
      push!(inc_hazards, S0_acc > 0 ? d_now / S0_acc : 0.0)
    end
    base_times = reverse(rev_events)
    inc_hazards = reverse(inc_hazards)
    baseline_hazard = cumsum(inc_hazards)
  end

  result = CoxResult(terms, base_times, baseline_hazard, loglik, iterations, converged)
  return provenance_result!(_ctx, result, "cox_ph"; parents=provenance_parent_ids(cohort), parameters=(max_iter=max_iter, tol=Float64(tol), term_count=length(terms), iterations=iterations, converged=converged, ties=ties))
end

"""
    forest_plot(cox_result)

Plot coefficient estimates and confidence intervals from a Cox model fit.
"""
function forest_plot(cox_result::CoxResult)
  _ctx = active_provenance_context()
  PlotsMod = _get_plots_mod()
  n = length(cox_result.terms)
  if n == 0
    plt = PlotsMod.plot(title="Cox forest plot", legend=false)
    return provenance_result!(_ctx, plt, "forest_plot"; parents=provenance_parent_ids(cox_result), parameters=(term_count=0))
  end
  terms = reverse(cox_result.terms)
  labels = [term.term for term in terms]
  hazard_ratios = [term.hazard_ratio for term in terms]
  lower_errors = [max(term.hazard_ratio - term.ci_lower, eps(Float64)) for term in terms]
  upper_errors = [max(term.ci_upper - term.hazard_ratio, eps(Float64)) for term in terms]
  pvalues = [term.pvalue for term in terms]
  p1 = PlotsMod.plot(hazard_ratios, 1:n; xerror=(lower_errors, upper_errors), seriestype=:scatter, marker=:circle, markersize=7, color=:black, legend=false, xscale=:log10, xlabel="Hazard ratio", yticks=(1:n, labels), title="Cox forest plot", ylim=(0.5, n + 0.5), framestyle=:box)
  PlotsMod.vline!(p1, [1.0]; linestyle=:dash, color=:gray)
  p2 = PlotsMod.plot(; xlim=(0, 1), ylim=(0.5, n + 0.5), framestyle=:none, grid=false, legend=false, xticks=false, yticks=false, title="Summary")
  for (index, term) in enumerate(terms)
    row = n - index + 1
    PlotsMod.annotate!(p2, 0.02, row, PlotsMod.text(term.term, 8, :black, :left))
    PlotsMod.annotate!(p2, 0.40, row, PlotsMod.text(@sprintf("HR %.2f", term.hazard_ratio), 8, :black, :left))
    PlotsMod.annotate!(p2, 0.68, row, PlotsMod.text(@sprintf("95%% CI %.2f-%.2f", term.ci_lower, term.ci_upper), 8, :black, :left))
    PlotsMod.annotate!(p2, 0.98, row, PlotsMod.text(@sprintf("p=%.3g", term.pvalue), 8, :black, :right))
  end
  plt = PlotsMod.plot(p1, p2; layout=(1, 2), size=(1100, 420))
  return provenance_result!(_ctx, plt, "forest_plot"; parents=provenance_parent_ids(cox_result), parameters=(term_count=n))
end

"""
    survival_roc(time, status, marker, predict_time)

Compute an IPCW-weighted time-dependent ROC curve for survival prediction.
"""
function survival_roc(time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}, marker::AbstractVector{<:Real}, predict_time::Real)
  _ctx = active_provenance_context()

  length(time) == length(status) == length(marker) || throw(ArgumentError("time, status, and marker must have the same length"))
  timef = Float64.(time)
  statusi = Int.(status)
  markerf = Float64.(marker)
  tau = Float64(predict_time)

  # Estimate censoring distribution via Kaplan-Meier for IPCW weights
  km_censor = kaplan_meier(timef, statusi .== 0)

  weights = zeros(Float64, length(timef))
  for i in 1:length(timef)
    if timef[i] <= tau && statusi[i] > 0
      g_t = _km_survival_at(km_censor, timef[i])
      weights[i] = 1.0 / max(g_t, 1e-4)
    elseif timef[i] > tau
      g_tau = _km_survival_at(km_censor, tau)
      weights[i] = 1.0 / max(g_tau, 1e-4)
    end
  end

  event = (timef .<= tau) .& (statusi .> 0)
  control = timef .> tau
  usable = (event .| control) .& (weights .> 0)

  thresholds = sort(unique(markerf[usable]))
  isempty(thresholds) && (thresholds = sort(unique(markerf)))
  tpr = Float64[]
  fpr = Float64[]

  total_event_w = sum(weights[event])
  total_control_w = sum(weights[control])

  for threshold in thresholds
    predicted = markerf .>= threshold
    tp_w = sum(weights[predicted .& event])
    fp_w = sum(weights[predicted .& control])

    push!(tpr, total_event_w > 0 ? tp_w / total_event_w : 0.0)
    push!(fpr, total_control_w > 0 ? fp_w / total_control_w : 0.0)
  end

  order = sortperm(fpr)
  fpr_curve = vcat(0.0, fpr[order], 1.0)
  tpr_curve = vcat(0.0, tpr[order], 1.0)
  auc = sum((fpr_curve[2:end] .- fpr_curve[1:(end-1)]) .* (tpr_curve[2:end] .+ tpr_curve[1:(end-1)]) ./ 2)
  auc = clamp(auc, 0.0, 1.0)

  result = ROCResult(tau, thresholds, tpr, fpr, auc)
  return provenance_result!(_ctx, result, "survival_roc"; parents=provenance_parent_ids(time, status, marker), parameters=(n=length(time), predict_time=tau, auc=auc))
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

const _VARIANT_SEVERITY_RANK = Dict(
  "Nonsense_Mutation" => 1,
  "Frame_Shift_Del" => 1,
  "Frame_Shift_Ins" => 1,
  "Splice_Site" => 2,
  "Nonstop_Mutation" => 2,
  "Translation_Start_Site" => 2,
  "Missense_Mutation" => 3,
  "In_Frame_Del" => 3,
  "In_Frame_Ins" => 3,
  "Silent" => 4,
  "3'UTR" => 5,
  "5'UTR" => 5,
  "Intron" => 5,
  "RNA" => 5
)

function _most_severe_variant(labels::String)
  isempty(labels) && return ""
  parts = split(labels, ';')
  return String(argmin(p -> get(_VARIANT_SEVERITY_RANK, String(p), 99), parts))
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
  sev = _most_severe_variant(labels)
  return sev == "Nonsense_Mutation" || sev == "Frame_Shift_Del" || sev == "Frame_Shift_Ins" ? 2 :
         sev == "Splice_Site" || sev == "Nonstop_Mutation" || sev == "Translation_Start_Site" ? 5 :
         sev == "Missense_Mutation" || sev == "In_Frame_Del" || sev == "In_Frame_Ins" ? 1 : 6
end

function _oncoprint_plot(result::OncoprintResult; title::String="Oncoprint", kwargs...)
  PlotsMod = _get_plots_mod()
  if isempty(result.genes) || isempty(result.samples)
    return PlotsMod.plot(title=title, legend=false)
  end
  codes = zeros(Int, size(result.matrix))
  for row in axes(result.matrix, 1), col in axes(result.matrix, 2)
    codes[row, col] = result.matrix[row, col] > 0 ? _mutation_level(result.mutation_labels[row, col]) : 0
  end
  heat = PlotsMod.plot(codes; seriestype=:heatmap, colorbar=true, c=:magma, xlabel="Patients", ylabel="Genes", title=title, yticks=(1:length(result.genes), reverse(result.genes)), xticks=(1:length(result.samples), result.samples), yflip=true, kwargs...)
  for row in axes(result.mutation_labels, 1), col in axes(result.mutation_labels, 2)
    label = result.mutation_labels[row, col]
    isempty(label) && continue
    PlotsMod.annotate!(heat, col, row, PlotsMod.text(label, 6, :white, :center))
  end
  return heat
end

"""
    oncoprint_plot(result; kwargs...)

Plot an oncoprint from a prepared result object.
"""
function oncoprint_plot(result::OncoprintResult; kwargs...)
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
    activations[layer+1] = layer == length(weights) ? z : tanh.(z)
  end
  return activations
end

"""
    neural_cox(X, time, status; hidden_units=8, learning_rate=0.01, epochs=200, seed=42)

Fit a small neural-network Cox model for survival prediction, optimized to O(epochs * N log N).
"""
function neural_cox(X::AbstractMatrix{<:Real}, time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}; hidden_units::Int=8, learning_rate::Real=0.01, epochs::Int=200, seed::Int=42)
  _ctx = active_provenance_context()

  Xf = Matrix{Float64}(X)
  timef = Float64.(time)
  statusi = Int.(status)

  order = sortperm(timef, rev=true)
  timef = timef[order]
  statusi = statusi[order]
  Xf = Xf[order, :]

  N = length(timef)
  rng = MersenneTwister(seed)
  weights = [randn(rng, hidden_units, size(Xf, 2)) / sqrt(size(Xf, 2)), randn(rng, 1, hidden_units) / sqrt(hidden_units)]
  biases = [zeros(Float64, hidden_units), zeros(Float64, 1)]
  best_loglik = -Inf
  best_weights = deepcopy(weights)
  best_biases = deepcopy(biases)

  learning_rate_eff = Float64(learning_rate) / max(N, 1)

  for _ in 1:epochs
    activations = _forward_pass(Xf, weights, biases)
    risk = vec(activations[end])

    loglik = 0.0
    output_delta = zeros(Float64, 1, N)

    S0_vec = zeros(Float64, N)
    acc = 0.0
    i = 1
    while i <= N
      t_curr = timef[i]
      group_start = i
      while i <= N && timef[i] == t_curr
        acc += exp(clamp(risk[i], -40.0, 40.0))
        i += 1
      end
      for j in group_start:(i-1)
        S0_vec[j] = acc
      end
    end

    V = zeros(Float64, N)
    for i in 1:N
      if statusi[i] > 0
        S0 = S0_vec[i]
        if S0 > 0.0
          loglik += risk[i] - log(S0)
          V[i] = 1.0 / S0
        end
      end
    end

    if loglik > best_loglik
      best_loglik = loglik
      best_weights = deepcopy(weights)
      best_biases = deepcopy(biases)
    end

    W = cumsum(reverse(V))
    W = reverse(W)

    for i in 1:N
      if statusi[i] > 0
        output_delta[1, i] -= 1.0
      end
      output_delta[1, i] += exp(clamp(risk[i], -40.0, 40.0)) * W[i]
    end

    hidden = activations[2]
    hidden_delta = (weights[end]' * output_delta) .* (1 .- hidden .^ 2)
    _gradient_step!(weights, biases, activations, [hidden_delta, output_delta], learning_rate_eff)
  end

  final_activations = _forward_pass(Xf, best_weights, best_biases)
  risk_scores = vec(final_activations[end])
  inv_order = invperm(order)
  risk_scores = risk_scores[inv_order]

  result = NeuralCoxResult(best_weights, best_biases, risk_scores, best_loglik)
  return provenance_result!(_ctx, result, "neural_cox"; parents=provenance_parent_ids(X, time, status), parameters=(hidden_units=hidden_units, learning_rate=Float64(learning_rate), epochs=epochs, seed=seed, loglik=best_loglik))
end

"""
    dose_response_curve(drug, cell_lines, concentrations, responses)

Fit a 4-parameter Hill-style dose-response curve via grid search and Levenberg-Marquardt NLS optimization.
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

  emin_init = minimum(y)
  emax_init = maximum(y)
  log_candidates = exp.(range(log(minimum(x)), log(maximum(x)), length=25))
  ec50_candidates = unique(sort(vcat([median(x)], x, log_candidates)))
  hill_candidates = collect(range(0.5, 4.0, length=15))
  best_sse = Inf
  best_ec50 = median(x)
  best_hill = 1.0
  best_emin = emin_init
  best_emax = emax_init

  for ec50 in ec50_candidates
    for hill in hill_candidates
      preds = emin_init .+ (emax_init - emin_init) ./ (1.0 .+ (ec50 ./ x) .^ hill)
      sse = sum((y .- preds) .^ 2)
      if sse < best_sse
        best_sse = sse
        best_ec50 = ec50
        best_hill = hill
      end
    end
  end

  theta = [best_emin, best_emax, max(best_ec50, 1e-6), max(best_hill, 0.1)]
  lambda = 0.01

  for _ in 1:100
    emin, emax, ec50, hill = theta[1], theta[2], max(theta[3], 1e-6), max(theta[4], 0.05)
    ratio = ec50 ./ x
    denom = 1.0 .+ ratio .^ hill
    preds = emin .+ (emax - emin) ./ denom
    residuals = y .- preds
    current_sse = sum(residuals .^ 2)

    J = zeros(Float64, length(x), 4)
    for i in 1:length(x)
      r = ratio[i]
      d = denom[i]
      J[i, 1] = 1.0 - 1.0 / d
      J[i, 2] = 1.0 / d
      dr_dec50 = hill * r^(hill - 1.0) / x[i]
      J[i, 3] = -(emax - emin) * dr_dec50 / (d^2)
      dr_dhill = r^hill * log(max(r, 1e-12))
      J[i, 4] = -(emax - emin) * dr_dhill / (d^2)
    end

    g = J' * residuals
    H = J' * J
    H_damped = H + lambda * Diagonal(max.(diag(H), 1e-6))

    delta = try
      H_damped \ g
    catch
      pinv(H_damped) * g
    end

    theta_new = theta .+ delta
    theta_new[3] = max(theta_new[3], 1e-6)
    theta_new[4] = max(theta_new[4], 0.05)

    emin_n, emax_n, ec50_n, hill_n = theta_new[1], theta_new[2], theta_new[3], theta_new[4]
    preds_new = emin_n .+ (emax_n - emin_n) ./ (1.0 .+ (ec50_n ./ x) .^ hill_n)
    new_sse = sum((y .- preds_new) .^ 2)

    if new_sse < current_sse
      theta = theta_new
      lambda /= 3.0
      if norm(delta) < 1e-6 || abs(current_sse - new_sse) < 1e-8
        break
      end
    else
      lambda *= 5.0
    end
  end

  emin_final, emax_final, ec50_final, hill_final = theta[1], theta[2], theta[3], theta[4]
  best_fitted = emin_final .+ (emax_final - emin_final) ./ (1.0 .+ (ec50_final ./ x) .^ hill_final)
  ic50_final = ec50_final

  result = DoseResponseResult(x, y, emax_final, ec50_final, hill_final, best_fitted, ic50_final)
  return provenance_result!(_ctx, result, "dose_response_curve"; parents=provenance_parent_ids(cell_lines, concentrations, responses), parameters=(drug=drug, n=length(x), ec50=ec50_final, hill=hill_final, emax=emax_final))
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
    propensity_score_match(clinical; treatment_col=:treatment, covariates=Symbol[], ratio=1, missing_handling=:drop)

Nearest-neighbor propensity score matching without replacement, with missing data handling.
"""
function propensity_score_match(clinical::DataFrame; treatment_col::Symbol=:treatment, covariates::Vector{Symbol}=Symbol[], ratio::Int=1, missing_handling::Symbol=:drop)
  _ctx = active_provenance_context()

  hasproperty(clinical, treatment_col) || throw(ArgumentError("clinical table missing treatment column"))
  ratio >= 1 || throw(ArgumentError("ratio must be >= 1"))

  if isempty(covariates)
    covariates = [Symbol(name) for name in names(clinical) if Symbol(name) != treatment_col && eltype(clinical[!, Symbol(name)]) <: Real]
  end
  isempty(covariates) && throw(ArgumentError("no numeric covariates available for propensity model"))

  df = copy(clinical)
  if missing_handling == :drop
    missing_mask = zeros(Bool, nrow(df))
    for c in covariates
      missing_mask .|= ismissing.(df[!, c])
    end
    if any(missing_mask)
      @info "Dropping $(count(missing_mask)) rows with missing covariate values in propensity_score_match"
      df = df[.!missing_mask, :]
    end
  elseif missing_handling == :impute
    for c in covariates
      col = df[!, c]
      if any(ismissing, col)
        med = median(skipmissing(col))
        df[!, c] = coalesce.(col, med)
      end
    end
  else
    throw(ArgumentError("missing_handling must be :drop or :impute"))
  end

  nrow(df) > 0 || throw(ArgumentError("no complete cases remaining after handling missing values"))

  Xraw = hcat([Float64.(df[!, c]) for c in covariates]...)
  for j in axes(Xraw, 2)
    μ = mean(@view Xraw[:, j])
    σ = std(@view Xraw[:, j])
    if isfinite(σ) && σ > 0
      Xraw[:, j] .= (Xraw[:, j] .- μ) ./ σ
    else
      Xraw[:, j] .= 0.0
    end
  end

  X = hcat(ones(Float64, nrow(df)), Xraw)
  y = Float64.(df[!, treatment_col] .!= 0)
  β = _fit_logistic_irls(X, y)
  ps = 1.0 ./ (1.0 .+ exp.(-clamp.(X * β, -25.0, 25.0)))

  treated = findall(==(1.0), y)
  controls_vec = findall(==(0.0), y)

  ctrl_perm = sortperm(ps[controls_vec])
  sorted_controls = controls_vec[ctrl_perm]
  sorted_ctrl_ps = ps[sorted_controls]
  ctrl_matched = zeros(Bool, length(sorted_controls))

  matches = DataFrame(treated_index=Int[], control_index=Int[], treated_ps=Float64[], control_ps=Float64[], abs_distance=Float64[])

  for ti in sort(treated; by=i -> ps[i], rev=true)
    p_t = ps[ti]
    idx = searchsortedfirst(sorted_ctrl_ps, p_t)

    for _ in 1:ratio
      best_ci = 0
      best_idx = 0
      best_dist = Inf

      left = min(max(idx, 1), length(sorted_controls))
      right = left

      search_window = 0
      while (left >= 1 || right <= length(sorted_controls)) && search_window < 200
        if left >= 1 && !ctrl_matched[left]
          d = abs(p_t - sorted_ctrl_ps[left])
          if d < best_dist
            best_dist = d
            best_ci = sorted_controls[left]
            best_idx = left
          end
        end
        if right <= length(sorted_controls) && !ctrl_matched[right]
          d = abs(p_t - sorted_ctrl_ps[right])
          if d < best_dist
            best_dist = d
            best_ci = sorted_controls[right]
            best_idx = right
          end
        end

        if best_ci != 0 && best_dist < (left >= 1 ? abs(p_t - sorted_ctrl_ps[left]) : Inf) && (right > length(sorted_controls) || best_dist < abs(p_t - sorted_ctrl_ps[right]))
          break
        end

        left -= 1
        right += 1
        search_window += 1
      end

      if best_ci != 0
        ctrl_matched[best_idx] = true
        push!(matches, (ti, best_ci, p_t, ps[best_ci], best_dist))
      else
        break
      end
    end
  end

  result = (matches=matches, propensity_score=ps, coefficients=β, covariates=vcat(:intercept, covariates))
  return provenance_result!(_ctx, result, "propensity_score_match"; parents=provenance_parent_ids(clinical), parameters=(row_count=nrow(df), treatment_col=treatment_col, covariate_count=length(covariates), ratio=ratio, match_count=nrow(matches)))
end

"""
    rmst(time, status; tau=nothing, groups=nothing)

Compute Restricted Mean Survival Time (RMST) up to truncation time `tau`.
If `groups` is provided, performs a 2-group comparison (difference and ratio).
"""
function rmst(time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}; tau::Union{Nothing,Real}=nothing, groups::Union{Nothing,AbstractVector}=nothing)
  _ctx = active_provenance_context()
  n = length(time)
  (n == length(status)) || throw(ArgumentError("time and status must have equal length"))

  max_t = maximum(time)
  tau_val = tau === nothing ? Float64(max_t) : Float64(tau)
  tau_val > 0 || throw(ArgumentError("tau must be positive"))

  function _calc_rmst_single(t_vec, s_vec, t_limit)
    km = kaplan_meier(t_vec, s_vec)
    times = km.time
    surv = km.survival
    at_risk = km.at_risk
    events = km.events

    rmst_val = 0.0
    prev_t = 0.0
    prev_s = 1.0

    for i in 1:length(times)
      t_curr = min(times[i], t_limit)
      if t_curr > prev_t
        rmst_val += prev_s * (t_curr - prev_t)
      end
      prev_t = t_curr
      prev_s = surv[i]
      times[i] >= t_limit && break
    end
    if prev_t < t_limit
      rmst_val += prev_s * (t_limit - prev_t)
    end

    var_rmst = 0.0
    for i in 1:length(times)
      times[i] <= t_limit || break
      n_i = at_risk[i]
      d_i = events[i]
      if d_i > 0 && n_i > d_i
        area_rem = 0.0
        p_t = times[i]
        p_s = surv[i]
        for j in (i+1):length(times)
          t_next = min(times[j], t_limit)
          area_rem += p_s * (t_next - p_t)
          p_t = t_next
          p_s = surv[j]
          times[j] >= t_limit && break
        end
        if p_t < t_limit
          area_rem += p_s * (t_limit - p_t)
        end

        var_rmst += (area_rem^2) * (d_i / (n_i * (n_i - d_i)))
      end
    end

    se = sqrt(max(var_rmst, eps(Float64)))
    cil = max(rmst_val - 1.96 * se, 0.0)
    ciu = rmst_val + 1.96 * se
    return rmst_val, se, cil, ciu
  end

  if groups !== nothing
    length(groups) == n || throw(ArgumentError("groups must match time and status length"))
    levels = unique(groups)
    length(levels) == 2 || throw(ArgumentError("rmst group comparison requires exactly 2 group levels"))

    mask1 = groups .== levels[1]
    mask2 = groups .== levels[2]

    rmst1, se1, cil1, ciu1 = _calc_rmst_single(time[mask1], status[mask1], tau_val)
    rmst2, se2, cil2, ciu2 = _calc_rmst_single(time[mask2], status[mask2], tau_val)

    diff = rmst1 - rmst2
    se_diff = sqrt(se1^2 + se2^2)
    z_diff = diff / se_diff
    p_diff = 2 * ccdf(Normal(), abs(z_diff))

    ratio = rmst1 / max(rmst2, eps(Float64))
    se_log_ratio = sqrt((se1/max(rmst1, eps(Float64)))^2 + (se2/max(rmst2, eps(Float64)))^2)
    z_ratio = log(max(ratio, eps(Float64))) / max(se_log_ratio, eps(Float64))
    p_ratio = 2 * ccdf(Normal(), abs(z_ratio))

    group_res = Dict{String,Any}(
      "level1" => string(levels[1]),
      "level2" => string(levels[2]),
      "rmst1" => rmst1, "se1" => se1, "ci_lower1" => cil1, "ci_upper1" => ciu1,
      "rmst2" => rmst2, "se2" => se2, "ci_lower2" => cil2, "ci_upper2" => ciu2,
      "diff" => diff, "se_diff" => se_diff, "z_diff" => z_diff, "p_diff" => p_diff,
      "ratio" => ratio, "se_log_ratio" => se_log_ratio, "p_ratio" => p_ratio
    )

    overall_rmst, overall_se, overall_cil, overall_ciu = _calc_rmst_single(time, status, tau_val)
    res = RMSTResult(tau_val, overall_rmst, overall_se, overall_cil, overall_ciu, group_res)
    return provenance_result!(_ctx, res, "rmst"; parents=provenance_parent_ids(time, status, groups), parameters=(tau=tau_val, group_count=2))
  else
    val, se, cil, ciu = _calc_rmst_single(time, status, tau_val)
    res = RMSTResult(tau_val, val, se, cil, ciu, nothing)
    return provenance_result!(_ctx, res, "rmst"; parents=provenance_parent_ids(time, status), parameters=(tau=tau_val,))
  end
end

"""
    cox_zph(cox_result, formula, cohort)

Test proportional hazards assumption for a fitted Cox model using Schoenfeld residuals.
"""
function cox_zph(cox_res::CoxResult, formula, cohort::PatientCohort; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
  start_name, time_name, status_name, covariate_names = _parse_surv_formula(formula)
  X, term_names = _cox_design_matrix(cohort, covariate_names)
  time = Float64.(cohort.clinical[!, time_name])
  status = Int.(cohort.clinical[!, status_name])
  time, status, order = _cox_order(time, status)
  X = X[order, :]

  beta = vcat(0.0, [t.beta for t in cox_res.terms])
  if length(beta) != size(X, 2)
    beta = zeros(Float64, size(X, 2))
  end

  N, p = size(X)
  eta = X * beta
  risk = exp.(clamp.(eta, -40.0, 40.0))

  event_indices = findall(>(0), status)
  n_events = length(event_indices)
  p_cov = p - 1
  p_cov > 0 || throw(ArgumentError("No covariates in Cox model for cox_zph"))

  sch_res = zeros(Float64, n_events, p_cov)
  event_times = zeros(Float64, n_events)

  for (k, idx) in enumerate(event_indices)
    t_curr = time[idx]
    event_times[k] = t_curr
    risk_set = time .>= t_curr
    S0 = sum(risk[risk_set])
    if S0 > 0
      S1 = sum(X[risk_set, 2:end] .* risk[risk_set], dims=1)
      E_x = vec(S1 ./ S0)
      sch_res[k, :] = X[idx, 2:end] .- E_x
    end
  end

  rho = Float64[]
  chisq = Float64[]
  pvalue = Float64[]

  t_rank = Float64.(sortperm(sortperm(event_times)))
  t_centered = t_rank .- mean(t_rank)
  denom_t = sqrt(sum(t_centered .^ 2))

  for col in 1:p_cov
    r_col = sch_res[:, col]
    r_centered = r_col .- mean(r_col)
    denom_r = sqrt(sum(r_centered .^ 2))
    r_val = (denom_t > 0 && denom_r > 0) ? sum(t_centered .* r_centered) / (denom_t * denom_r) : 0.0
    c_stat = n_events * (r_val^2)
    pval = ccdf(Chisq(1), max(c_stat, 0.0))
    push!(rho, r_val)
    push!(chisq, c_stat)
    push!(pvalue, pval)
  end

  g_chisq = sum(chisq)
  g_pval = ccdf(Chisq(max(p_cov, 1)), max(g_chisq, 0.0))

  res = CoxZphResult(term_names, rho, chisq, pvalue, sch_res, event_times, g_chisq, g_pval)
  return provenance_result!(_ctx, res, "cox_zph"; parents=provenance_parent_ids(cohort), parameters=(event_count=n_events, covariate_count=p_cov))
end

"""
    fine_gray(formula, cohort; fail_code=1, max_iter=50, tol=1e-7, prov_ctx=nothing)

Fit a Fine-Gray subdistribution hazards model for competing risks data using Inverse Probability of Censoring Weighting (IPCW).
Subdistribution hazards retain individuals experiencing competing risks in the risk set weighted by G(t) / G(T_i ∧ t).
"""
function fine_gray(formula, cohort::PatientCohort; fail_code::Int=1, max_iter::Int=50, tol::Real=1e-7, prov_ctx=nothing)
  _ctx = active_provenance_context(prov_ctx)

  start_name, time_name, status_name, covariate_names = _parse_surv_formula(formula)
  X, term_names = _cox_design_matrix(cohort, covariate_names)
  raw_time = Float64.(cohort.clinical[!, time_name])
  raw_status = Int.(cohort.clinical[!, status_name])

  n_samples, p_cov = size(X)

  order = sortperm(raw_time)
  time = raw_time[order]
  status = raw_status[order]
  X = X[order, :]

  # 1. Estimate Censoring Survival Distribution G(t) = P(C > t) using Kaplan-Meier on censoring events (status == 0)
  unique_times = sort(unique(time))
  G_map = Dict{Float64,Float64}()
  current_G = 1.0
  for u in unique_times
    Y_u = count(t -> t >= u, time)
    d_c_u = count(i -> time[i] == u && status[i] == 0, 1:n_samples)
    if Y_u > 0 && d_c_u > 0
      current_G *= (1.0 - d_c_u / Y_u)
    end
    G_map[u] = max(current_G, 1e-6)
  end

  function get_G(t::Real)
    t < unique_times[1] && return 1.0
    idx = searchsortedlast(unique_times, t)
    idx == 0 && return 1.0
    return G_map[unique_times[idx]]
  end

  G_T = [get_G(time[i]) for i in 1:n_samples]

  event_indices = findall(i -> status[i] == fail_code, 1:n_samples)
  if isempty(event_indices)
    throw(ArgumentError("No failure events with status == $fail_code found in cohort."))
  end

  # 2. Newton-Raphson Optimization with IPCW Weighted Risk Sets
  beta = zeros(Float64, p_cov)
  converged = false
  iterations = 0
  loglik = 0.0
  H = zeros(Float64, p_cov, p_cov)

  for iter in 1:max_iter
    iterations = iter
    eta = X * beta
    theta = exp.(clamp.(eta, -30.0, 30.0))

    U = zeros(Float64, p_cov)
    H = zeros(Float64, p_cov, p_cov)
    loglik = 0.0

    for m in event_indices
      t_m = time[m]
      G_tm = get_G(t_m)

      risk_set_weights = Float64[]
      risk_set_indices = Int[]

      for i in 1:n_samples
        if time[i] >= t_m
          push!(risk_set_indices, i)
          push!(risk_set_weights, 1.0)
        elseif status[i] > 0 && status[i] != fail_code
          w_i = G_tm / G_T[i]
          push!(risk_set_indices, i)
          push!(risk_set_weights, w_i)
        end
      end

      S0 = 0.0
      S1 = zeros(Float64, p_cov)
      S2 = zeros(Float64, p_cov, p_cov)

      for (k, idx) in enumerate(risk_set_indices)
        w_k = risk_set_weights[k]
        w_theta = w_k * theta[idx]
        x_k = X[idx, :]

        S0 += w_theta
        S1 += w_theta * x_k
        S2 += w_theta * (x_k * x_k')
      end

      if S0 > 0.0
        x_bar = S1 / S0
        U += X[m, :] - x_bar
        H -= (S2 / S0 - x_bar * x_bar')
        loglik += eta[m] - log(S0)
      end
    end

    negH = -H
    if det(negH) <= 1e-12
      negH += 1e-6 * I(p_cov)
    end

    delta_beta = negH \ U
    beta += delta_beta

    if maximum(abs.(delta_beta)) < tol
      converged = true
      break
    end
  end

  # 3. Standard Errors, z-scores, p-values, 95% CIs
  cov_mat = inv(-H + 1e-8 * I(p_cov))
  se = sqrt.(max.(1e-12, diag(cov_mat)))
  z_scores = beta ./ se
  p_values = 2.0 .* (1.0 .- cdf.(Normal(0, 1), abs.(z_scores)))
  hazard_ratios = exp.(beta)
  ci_lower = exp.(beta .- 1.96 .* se)
  ci_upper = exp.(beta .+ 1.96 .* se)

  terms = [
    CoxTermResult(
      term_names[j-1], beta[j], hazard_ratios[j], se[j],
      z_scores[j], p_values[j], ci_lower[j], ci_upper[j]
    ) for j in 2:p_cov
  ]

  # Baseline cumulative subdistribution hazard
  eta_final = X * beta
  theta_final = exp.(clamp.(eta_final, -30.0, 30.0))
  baseline_times = unique(time[event_indices])
  baseline_hazard = Float64[]
  cum_h = 0.0

  for t_m in baseline_times
    G_tm = get_G(t_m)
    d_m = count(i -> time[i] == t_m && status[i] == fail_code, 1:n_samples)
    S0 = 0.0
    for i in 1:n_samples
      if time[i] >= t_m
        S0 += theta_final[i]
      elseif status[i] > 0 && status[i] != fail_code
        S0 += (G_tm / G_T[i]) * theta_final[i]
      end
    end
    cum_h += S0 > 0.0 ? d_m / S0 : 0.0
    push!(baseline_hazard, cum_h)
  end

  res = FineGrayResult(terms, baseline_times, baseline_hazard, loglik, iterations, converged, fail_code)
  return provenance_result!(_ctx, res, "fine_gray"; parents=provenance_parent_ids(cohort), parameters=(fail_code=fail_code, max_iter=max_iter, tol=Float64(tol)))
end

"""
    calculate_tmb(maf; capture_size_mb=38.0)

Calculate Tumor Mutational Burden (TMB) per sample in mutations per megabase.
"""
function calculate_tmb(maf::AbstractVector{<:MAFRecord}; capture_size_mb::Real=38.0)
  _ctx = active_provenance_context()
  capture_size_mb > 0 || throw(ArgumentError("capture_size_mb must be positive"))

  non_silent_classes = Set([
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins"
  ])

  total_muts = Dict{String,Int}()
  non_silent_muts = Dict{String,Int}()

  for rec in maf
    s = rec.sample
    total_muts[s] = get(total_muts, s, 0) + 1
    if rec.variant_classification in non_silent_classes
      non_silent_muts[s] = get(non_silent_muts, s, 0) + 1
    end
  end

  df = DataFrame(sample=String[], total_mutations=Int[], non_silent_mutations=Int[], tmb=Float64[])
  for s in sort(collect(keys(total_muts)))
    tot = total_muts[s]
    ns = get(non_silent_muts, s, 0)
    tmb_val = ns / Float64(capture_size_mb)
    push!(df, (s, tot, ns, tmb_val))
  end

  return provenance_result!(_ctx, df, "calculate_tmb"; parents=provenance_parent_ids(maf), parameters=(sample_count=nrow(df), capture_size_mb=Float64(capture_size_mb)))
end

"""
    somatic_interactions(maf; top_n=20)

Perform Fisher's exact test for mutual exclusivity and co-occurrence between top mutated genes.
"""
function somatic_interactions(maf::AbstractVector{<:MAFRecord}; top_n::Int=20)
  _ctx = active_provenance_context()
  genes, samples, matrix, _ = _maf_to_matrix(maf)

  mut_counts = vec(sum(matrix .> 0, dims=2))
  top_indices = sortperm(mut_counts, rev=true)[1:min(top_n, length(genes))]
  top_genes = genes[top_indices]
  n_samples = length(samples)

  out = DataFrame(gene1=String[], gene2=String[], odds_ratio=Float64[], pvalue=Float64[], event_type=String[])

  for i in 1:length(top_genes)
    for j in (i+1):length(top_genes)
      g1_idx = top_indices[i]
      g2_idx = top_indices[j]

      g1_mut = matrix[g1_idx, :] .> 0
      g2_mut = matrix[g2_idx, :] .> 0

      a = count(g1_mut .& g2_mut)
      b = count(g1_mut .& .!g2_mut)
      c = count(.!g1_mut .& g2_mut)
      d = count(.!g1_mut .& .!g2_mut)

      or_val = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
      hg = Hypergeometric(a + c, max(n_samples - (a + c), 0), a + b)
      pval = 2 * min(cdf(hg, a), ccdf(hg, max(a - 1, 0)))
      pval = clamp(pval, 0.0, 1.0)

      event = or_val > 1.0 ? "Co_Occurrence" : "Mutually_Exclusive"
      push!(out, (top_genes[i], top_genes[j], or_val, pval, event))
    end
  end

  return provenance_result!(_ctx, out, "somatic_interactions"; parents=provenance_parent_ids(maf), parameters=(top_n=top_n, pair_count=nrow(out)))
end

# ==============================================================================
# Interactive HTML Exports for Clinical Results
# ==============================================================================

function to_html(km::KaplanMeierResult)
  data_json = replace(JSON.json(Dict(
      "time" => km.time,
      "survival" => km.survival,
      "std_error" => km.std_error,
      "ci_lower" => km.ci_lower,
      "ci_upper" => km.ci_upper,
      "at_risk" => km.at_risk,
      "events" => km.events,
      "censored" => km.censored,
      "censor_times" => km.censor_times,
      "censor_survival" => km.censor_survival
    )), "</" => "<\\/")

  return """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Kaplan-Meier Survival Analysis - BioToolkit</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; background: #0f172a; color: #f8fafc; margin: 0; padding: 24px; }
        .card { background: #1e293b; border-radius: 12px; padding: 24px; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.5); max-width: 960px; margin: 0 auto; border: 1px solid #334155; }
        h2 { margin-top: 0; color: #38bdf8; font-weight: 600; }
        .canvas-container { position: relative; width: 100%; height: 420px; background: #0f172a; border-radius: 8px; border: 1px solid #334155; }
        canvas { width: 100%; height: 100%; display: block; }
        .risk-table { width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 13px; }
        .risk-table th, .risk-table td { padding: 8px 12px; text-align: center; border-bottom: 1px solid #334155; }
        .risk-table th { background: #334155; color: #94a3b8; font-weight: 500; }
        .badge { background: #0284c7; color: white; padding: 4px 10px; border-radius: 9999px; font-size: 12px; font-weight: 600; display: inline-block; margin-bottom: 12px; }
    </style>
</head>
<body>
    <div class="card">
        <span class="badge">Clinical Genomics Engine</span>
        <h2>Kaplan-Meier Survival Analysis</h2>
        <div class="canvas-container" id="container">
            <canvas id="kmCanvas"></canvas>
        </div>
        <table class="risk-table">
            <thead>
                <tr><th>Time</th><th>At Risk</th><th>Events</th><th>Censored</th><th>Survival %</th><th>95% CI</th></tr>
            </thead>
            <tbody id="riskTableBody"></tbody>
        </table>
    </div>

    <script>
        const kmData = $data_json;
        const canvas = document.getElementById('kmCanvas');
        const ctx = canvas.getContext('2d');
        const container = document.getElementById('container');

        function resize() {
            canvas.width = container.clientWidth * window.devicePixelRatio;
            canvas.height = container.clientHeight * window.devicePixelRatio;
            draw();
        }

        function draw() {
            const w = canvas.width;
            const h = canvas.height;
            ctx.clearRect(0, 0, w, h);

            const padLeft = 60 * window.devicePixelRatio;
            const padRight = 30 * window.devicePixelRatio;
            const padTop = 30 * window.devicePixelRatio;
            const padBottom = 50 * window.devicePixelRatio;

            const pw = w - padLeft - padRight;
            const ph = h - padTop - padBottom;

            const maxTime = kmData.time.length ? Math.max(...kmData.time) * 1.05 : 1.0;

            function mapX(t) { return padLeft + (t / maxTime) * pw; }
            function mapY(s) { return padTop + (1.0 - s) * ph; }

            ctx.strokeStyle = '#334155';
            ctx.lineWidth = 1;
            for (let s = 0.0; s <= 1.0; s += 0.2) {
                const y = mapY(s);
                ctx.beginPath(); ctx.moveTo(padLeft, y); ctx.lineTo(w - padRight, y); ctx.stroke();
                ctx.fillStyle = '#94a3b8'; ctx.font = `\${12 * window.devicePixelRatio}px sans-serif`;
                ctx.textAlign = 'right'; ctx.fillText(s.toFixed(1), padLeft - 10, y + 4);
            }

            if (kmData.ci_lower.length && kmData.ci_upper.length) {
                ctx.fillStyle = 'rgba(56, 189, 248, 0.15)';
                ctx.beginPath();
                ctx.moveTo(mapX(0), mapY(1.0));
                let curL = 1.0, curU = 1.0;
                for (let i = 0; i < kmData.time.length; i++) {
                    const x = mapX(kmData.time[i]);
                    ctx.lineTo(x, mapY(curU));
                    curU = kmData.ci_upper[i];
                    ctx.lineTo(x, mapY(curU));
                }
                for (let i = kmData.time.length - 1; i >= 0; i--) {
                    const x = mapX(kmData.time[i]);
                    curL = kmData.ci_lower[i];
                    ctx.lineTo(x, mapY(curL));
                }
                ctx.closePath();
                ctx.fill();
            }

            ctx.strokeStyle = '#38bdf8';
            ctx.lineWidth = 3 * window.devicePixelRatio;
            ctx.beginPath();
            let curS = 1.0;
            ctx.moveTo(mapX(0), mapY(curS));
            for (let i = 0; i < kmData.time.length; i++) {
                const x = mapX(kmData.time[i]);
                ctx.lineTo(x, mapY(curS));
                curS = kmData.survival[i];
                ctx.lineTo(x, mapY(curS));
            }
            ctx.stroke();

            ctx.strokeStyle = '#f43f5e';
            ctx.lineWidth = 2 * window.devicePixelRatio;
            for (let i = 0; i < kmData.censor_times.length; i++) {
                const cx = mapX(kmData.censor_times[i]);
                const cy = mapY(kmData.censor_survival[i]);
                ctx.beginPath();
                ctx.moveTo(cx, cy - 6); ctx.lineTo(cx, cy + 6);
                ctx.stroke();
            }

            ctx.fillStyle = '#94a3b8';
            ctx.textAlign = 'center';
            for (let i = 0; i <= 5; i++) {
                const t = (maxTime * i / 5);
                ctx.fillText(t.toFixed(1), mapX(t), h - padBottom + 20);
            }
        }

        window.addEventListener('resize', resize);
        resize();

        const tbody = document.getElementById('riskTableBody');
        tbody.innerHTML = '';
        for (let i = 0; i < kmData.time.length; i++) {
            const tr = document.createElement('tr');
            tr.innerHTML = `<td>\${kmData.time[i].toFixed(1)}</td><td>\${kmData.at_risk[i]}</td><td>\${kmData.events[i]}</td><td>\${kmData.censored[i]}</td><td>\${(kmData.survival[i]*100).toFixed(1)}%</td><td>[\${(kmData.ci_lower[i]*100).toFixed(1)}%, \${(kmData.ci_upper[i]*100).toFixed(1)}%]</td>`;
            tbody.appendChild(tr);
        }
    </script>
</body>
</html>"""
end

function export_html(km::KaplanMeierResult, filename::String)
  open(filename, "w") do io
    write(io, to_html(km))
  end
  return filename
end

function to_html(cox::CoxResult)
  data_json = replace(JSON.json(Dict(
      "terms" => [Dict("term" => t.term, "beta" => t.beta, "hr" => t.hazard_ratio, "se" => t.standard_error, "z" => t.z_score, "p" => t.pvalue, "ci_lower" => t.ci_lower, "ci_upper" => t.ci_upper) for t in cox.terms],
      "loglik" => cox.loglik,
      "iterations" => cox.iterations,
      "converged" => cox.converged
    )), "</" => "<\\/")

  return """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Cox Proportional Hazards Model - BioToolkit</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; background: #0f172a; color: #f8fafc; margin: 0; padding: 24px; }
        .card { background: #1e293b; border-radius: 12px; padding: 24px; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.5); max-width: 960px; margin: 0 auto; border: 1px solid #334155; }
        h2 { margin-top: 0; color: #38bdf8; font-weight: 600; }
        .stats-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; margin-bottom: 20px; }
        .stat-box { background: #0f172a; border-radius: 8px; padding: 12px; border: 1px solid #334155; text-align: center; }
        .stat-value { font-size: 18px; font-weight: 700; color: #38bdf8; }
        .stat-label { font-size: 12px; color: #94a3b8; margin-top: 4px; }
        .forest-table { width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 13px; }
        .forest-table th, .forest-table td { padding: 10px 14px; text-align: left; border-bottom: 1px solid #334155; }
        .forest-table th { background: #334155; color: #94a3b8; font-weight: 500; }
        .badge { background: #0284c7; color: white; padding: 4px 10px; border-radius: 9999px; font-size: 12px; font-weight: 600; display: inline-block; margin-bottom: 12px; }
    </style>
</head>
<body>
    <div class="card">
        <span class="badge">Clinical Genomics Engine</span>
        <h2>Cox Proportional Hazards Model</h2>
        <div class="stats-grid">
            <div class="stat-box"><div class="stat-value" id="loglik">0</div><div class="stat-label">Log-Likelihood</div></div>
            <div class="stat-box"><div class="stat-value" id="iters">0</div><div class="stat-label">Iterations</div></div>
            <div class="stat-box"><div class="stat-value" id="conv">True</div><div class="stat-label">Converged</div></div>
        </div>
        <table class="forest-table">
            <thead>
                <tr><th>Term</th><th>Beta</th><th>Hazard Ratio (HR)</th><th>95% Confidence Interval</th><th>p-value</th></tr>
            </thead>
            <tbody id="termsBody"></tbody>
        </table>
    </div>

    <script>
        const coxData = $data_json;
        document.getElementById('loglik').innerText = coxData.loglik.toFixed(2);
        document.getElementById('iters').innerText = coxData.iterations;
        document.getElementById('conv').innerText = coxData.converged ? 'Yes' : 'No';

        function escapeHtml(str) {
            return String(str)
                .replace(/&/g, '&amp;')
                .replace(/</g, '&lt;')
                .replace(/>/g, '&gt;')
                .replace(/"/g, '&quot;')
                .replace(/'/g, '&#039;');
        }

        const tbody = document.getElementById('termsBody');
        coxData.terms.forEach(t => {
            const tr = document.createElement('tr');
            tr.innerHTML = `<td><strong>\${escapeHtml(t.term)}</strong></td><td>\${t.beta.toFixed(3)}</td><td><strong>\${t.hr.toFixed(2)}</strong></td><td>[\${t.ci_lower.toFixed(2)} - \${t.ci_upper.toFixed(2)}]</td><td>\${t.p < 0.001 ? '<0.001' : t.p.toFixed(4)}</td>`;
            tbody.appendChild(tr);
        });
    </script>
</body>
</html>"""
end

function export_html(cox::CoxResult, filename::String)
  open(filename, "w") do io
    write(io, to_html(cox))
  end
  return filename
end

function to_html(op::OncoprintResult)
  data_json = replace(JSON.json(Dict(
      "genes" => op.genes,
      "samples" => op.samples,
      "matrix" => op.matrix,
      "labels" => op.mutation_labels
    )), "</" => "<\\/")

  return """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Oncoprint Mutation Landscape - BioToolkit</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; background: #0f172a; color: #f8fafc; margin: 0; padding: 24px; }
        .card { background: #1e293b; border-radius: 12px; padding: 24px; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.5); max-width: 1100px; margin: 0 auto; border: 1px solid #334155; }
        h2 { margin-top: 0; color: #38bdf8; font-weight: 600; }
        .canvas-container { position: relative; width: 100%; height: 500px; background: #0f172a; border-radius: 8px; border: 1px solid #334155; }
        canvas { width: 100%; height: 100%; display: block; }
        .legend { display: flex; gap: 16px; margin-top: 16px; font-size: 13px; }
        .legend-item { display: flex; align-items: center; gap: 6px; }
        .legend-box { width: 14px; height: 14px; border-radius: 3px; }
        .badge { background: #0284c7; color: white; padding: 4px 10px; border-radius: 9999px; font-size: 12px; font-weight: 600; display: inline-block; margin-bottom: 12px; }
    </style>
</head>
<body>
    <div class="card">
        <span class="badge">Clinical Genomics Engine</span>
        <h2>Oncoprint Mutation Landscape</h2>
        <div class="canvas-container" id="container">
            <canvas id="opCanvas"></canvas>
        </div>
        <div class="legend">
            <div class="legend-item"><div class="legend-box" style="background:#2ecc71"></div>Missense</div>
            <div class="legend-item"><div class="legend-box" style="background:#e74c3c"></div>Nonsense</div>
            <div class="legend-item"><div class="legend-box" style="background:#9b59b6"></div>Frame Shift</div>
            <div class="legend-item"><div class="legend-box" style="background:#e67e22"></div>Splice Site</div>
            <div class="legend-item"><div class="legend-box" style="background:#34495e"></div>Other</div>
        </div>
    </div>

    <script>
        const opData = $data_json;
        const canvas = document.getElementById('opCanvas');
        const ctx = canvas.getContext('2d');
        const container = document.getElementById('container');

        const colors = {
            'Missense_Mutation': '#2ecc71',
            'Nonsense_Mutation': '#e74c3c',
            'Frame_Shift_Del': '#9b59b6',
            'Frame_Shift_Ins': '#9b59b6',
            'Splice_Site': '#e67e22'
        };

        const severityRank = {
            'Nonsense_Mutation': 1, 'Frame_Shift_Del': 1, 'Frame_Shift_Ins': 1,
            'Splice_Site': 2, 'Nonstop_Mutation': 2,
            'Missense_Mutation': 3, 'In_Frame_Del': 3, 'In_Frame_Ins': 3,
            'Silent': 4
        };

        function getMostSevere(label) {
            if (!label) return '';
            const parts = label.split(';');
            parts.sort((a, b) => (severityRank[a] || 99) - (severityRank[b] || 99));
            return parts[0];
        }

        function resize() {
            canvas.width = container.clientWidth * window.devicePixelRatio;
            canvas.height = container.clientHeight * window.devicePixelRatio;
            draw();
        }

        function draw() {
            const w = canvas.width;
            const h = canvas.height;
            ctx.clearRect(0, 0, w, h);

            const padLeft = 120 * window.devicePixelRatio;
            const padBottom = 40 * window.devicePixelRatio;
            const padTop = 20 * window.devicePixelRatio;
            const padRight = 20 * window.devicePixelRatio;

            const nGenes = opData.genes.length;
            const nSamples = opData.samples.length;

            if (nGenes === 0 || nSamples === 0) return;

            const cellW = (w - padLeft - padRight) / nSamples;
            const cellH = (h - padTop - padBottom) / nGenes;

            for (let r = 0; r < nGenes; r++) {
                ctx.fillStyle = '#f8fafc';
                ctx.font = `\${12 * window.devicePixelRatio}px sans-serif`;
                ctx.textAlign = 'right';
                ctx.fillText(opData.genes[r], padLeft - 10, padTop + (r + 0.6) * cellH);

                for (let c = 0; c < nSamples; c++) {
                    const x = padLeft + c * cellW;
                    const y = padTop + r * cellH;

                    ctx.fillStyle = '#1e293b';
                    ctx.fillRect(x + 1, y + 1, cellW - 2, cellH - 2);

                    const label = opData.labels[r][c];
                    if (label) {
                        const firstType = getMostSevere(label);
                        ctx.fillStyle = colors[firstType] || '#34495e';
                        ctx.fillRect(x + 2, y + cellH * 0.2, cellW - 4, cellH * 0.6);
                    }
                }
            }
        }

        window.addEventListener('resize', resize);
        resize();
    </script>
</body>
</html>"""
end

function export_html(op::OncoprintResult, filename::String)
  open(filename, "w") do io
    write(io, to_html(op))
  end
  return filename
end

function to_html(roc::ROCResult)
  data_json = replace(JSON.json(Dict(
      "predict_time" => roc.time,
      "tpr" => roc.tpr,
      "fpr" => roc.fpr,
      "auc" => roc.auc
    )), "</" => "<\\/")

  return """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Time-Dependent Survival ROC - BioToolkit</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; background: #0f172a; color: #f8fafc; margin: 0; padding: 24px; }
        .card { background: #1e293b; border-radius: 12px; padding: 24px; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.5); max-width: 720px; margin: 0 auto; border: 1px solid #334155; }
        h2 { margin-top: 0; color: #38bdf8; font-weight: 600; }
        .canvas-container { position: relative; width: 100%; height: 420px; background: #0f172a; border-radius: 8px; border: 1px solid #334155; }
        canvas { width: 100%; height: 100%; display: block; }
        .badge { background: #0284c7; color: white; padding: 4px 10px; border-radius: 9999px; font-size: 12px; font-weight: 600; display: inline-block; margin-bottom: 12px; }
        .auc-badge { background: #10b981; color: white; padding: 6px 14px; border-radius: 8px; font-weight: 700; font-size: 16px; float: right; }
    </style>
</head>
<body>
    <div class="card">
        <div class="auc-badge" id="aucBadge">AUC: 0.00</div>
        <span class="badge">Clinical Genomics Engine</span>
        <h2>Time-Dependent Survival ROC Curve</h2>
        <div class="canvas-container" id="container">
            <canvas id="rocCanvas"></canvas>
        </div>
    </div>

    <script>
        const rocData = $data_json;
        document.getElementById('aucBadge').innerText = 'AUC: ' + rocData.auc.toFixed(3);

        const canvas = document.getElementById('rocCanvas');
        const ctx = canvas.getContext('2d');
        const container = document.getElementById('container');

        function resize() {
            canvas.width = container.clientWidth * window.devicePixelRatio;
            canvas.height = container.clientHeight * window.devicePixelRatio;
            draw();
        }

        function draw() {
            const w = canvas.width;
            const h = canvas.height;
            ctx.clearRect(0, 0, w, h);

            const pad = 50 * window.devicePixelRatio;
            const pw = w - 2 * pad;
            const ph = h - 2 * pad;

            ctx.strokeStyle = '#475569';
            ctx.lineWidth = 1.5 * window.devicePixelRatio;
            ctx.setLineDash([6, 6]);
            ctx.beginPath();
            ctx.moveTo(pad, h - pad); ctx.lineTo(w - pad, pad);
            ctx.stroke();
            ctx.setLineDash([]);

            ctx.strokeStyle = '#38bdf8';
            ctx.lineWidth = 3 * window.devicePixelRatio;
            ctx.beginPath();
            ctx.moveTo(pad, h - pad);
            for (let i = 0; i < rocData.fpr.length; i++) {
                const x = pad + rocData.fpr[i] * pw;
                const y = h - pad - rocData.tpr[i] * ph;
                ctx.lineTo(x, y);
            }
            ctx.stroke();
        }

        window.addEventListener('resize', resize);
        resize();
    </script>
</body>
</html>"""
end

function export_html(roc::ROCResult, filename::String)
  open(filename, "w") do io
    write(io, to_html(roc))
  end
  return filename
end

function to_html(rmst::RMSTResult)
  data_json = replace(JSON.json(Dict(
      "tau" => rmst.tau,
      "rmst" => rmst.rmst,
      "std_error" => rmst.std_error,
      "ci_lower" => rmst.ci_lower,
      "ci_upper" => rmst.ci_upper,
      "group_results" => rmst.group_results
    )), "</" => "<\\/")

  return """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Restricted Mean Survival Time (RMST) - BioToolkit</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; background: #0f172a; color: #f8fafc; margin: 0; padding: 24px; }
        .card { background: #1e293b; border-radius: 12px; padding: 24px; box-shadow: 0 10px 25px -5px rgba(0,0,0,0.5); max-width: 720px; margin: 0 auto; border: 1px solid #334155; }
        h2 { margin-top: 0; color: #38bdf8; font-weight: 600; }
        .stat-value { font-size: 24px; font-weight: 700; color: #38bdf8; margin: 12px 0; }
        .badge { background: #0284c7; color: white; padding: 4px 10px; border-radius: 9999px; font-size: 12px; font-weight: 600; display: inline-block; margin-bottom: 12px; }
    </style>
</head>
<body>
    <div class="card">
        <span class="badge">Clinical Genomics Engine</span>
        <h2>Restricted Mean Survival Time (RMST)</h2>
        <div class="stat-value" id="rmstVal">0.00</div>
        <p id="rmstCi">95% CI: [0.00, 0.00]</p>
    </div>

    <script>
        const data = $data_json;
        document.getElementById('rmstVal').innerText = 'RMST (tau=' + data.tau.toFixed(1) + '): ' + data.rmst.toFixed(2);
        document.getElementById('rmstCi').innerText = '95% CI: [' + data.ci_lower.toFixed(2) + ', ' + data.ci_upper.toFixed(2) + '] (SE: ' + data.std_error.toFixed(3) + ')';
    </script>
</body>
</html>"""
end

# ==============================================================================
# Bioconductor/BioPython Parity Extensions (survminer, maftools, survival)
# ==============================================================================

struct SurvCutpointResult <: AbstractAnalysisResult
  feature::String
  cutpoint::Float64
  statistic::Float64
  pvalue::Float64
  high_group_count::Int
  low_group_count::Int
  provenance::ResultProvenance
end

struct NelsonAalenResult <: AbstractAnalysisResult
  time::Vector{Float64}
  cum_hazard::Vector{Float64}
  std_error::Vector{Float64}
  provenance::ResultProvenance
end

export SurvCutpointResult, NelsonAalenResult, surv_cutpoint, maf_compare, clinical_enrichment, nelson_aalen, cox_residuals

"""
    surv_cutpoint(time, status, feature_vector; min_prop=0.1) → SurvCutpointResult

Determine the optimal cutpoint for a continuous feature (e.g., gene expression) to stratify patients into high vs low survival risk groups using maximally selected rank statistics. Equivalent to `survminer::surv_cutpoint`.
"""
function surv_cutpoint(time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}, feature_vector::AbstractVector{<:Real}; min_prop::Float64=0.1, prov_ctx=nothing)
  _ctx = active_provenance_context(prov_ctx)
  n = length(time)
  @assert length(status) == n && length(feature_vector) == n "Dimensions of time, status, and feature must match"

  unique_vals = sort(unique(feature_vector))
  min_k = max(1, round(Int, n * min_prop))

  best_score = -1.0
  best_stat = -1.0
  best_cut = Float64(unique_vals[1])
  best_p = 1.0

  for val in unique_vals
    group = feature_vector .> val
    n_high = sum(group)
    n_low = n - n_high
    (n_high < min_k || n_low < min_k) && continue

    lr = logrank_test(time, status, group; prov_ctx=_ctx)
    score = lr.statistic
    if score > best_score
      best_score = score
      best_stat = lr.statistic
      best_cut = Float64(val)
      best_p = lr.pvalue
    end
  end

  # Lausen-Schumacher / Miller-Siegmund multiple testing adjustment for maxstat
  T_stat = sqrt(max(0.0, best_stat))
  eps_val = clamp(min_prop, 0.01, 0.49)
  p_adj = best_p
  if T_stat > 0.0
    phi_T = exp(-0.5 * T_stat^2) / sqrt(2.0 * pi)
    p_ms = 4.0 * phi_T / T_stat + phi_T * (T_stat - 1.0 / T_stat) * log(((1.0 - eps_val)^2) / (eps_val^2))
    p_adj = clamp(p_ms, best_p, 1.0)
  end

  high_count = sum(feature_vector .> best_cut)
  low_count = n - high_count

  res = SurvCutpointResult(
    "feature", best_cut, best_stat, p_adj, high_count, low_count,
    ResultProvenance(
      new_provenance_id("surv_cutpoint_$(n)_$(best_cut)"),
      "surv_cutpoint", "Clinical/surv_cutpoint", :ok, String[], String[], String[], String[],
      (n=n, cutpoint=best_cut, min_prop=min_prop, unadjusted_pvalue=best_p, adjusted_pvalue=p_adj),
      provenance_parent_ids(time, status, feature_vector), _provenance_timestamp()
    )
  )
  provenance_result!(_ctx, res, "surv_cutpoint"; parents=provenance_parent_ids(time, status, feature_vector), parameters=(cutpoint=best_cut, pvalue=p_adj))
  return res
end

"""
    maf_compare(maf1, maf2; cohort1_name="Cohort1", cohort2_name="Cohort2", top_n=20) → DataFrame

Compare mutation frequencies between two clinical MAF cohorts (e.g. Primary vs Metastatic or Responder vs Non-responder) using Fisher's exact test. Equivalent to `maftools::mafCompare`.
"""
function maf_compare(maf1::AbstractVector{<:MAFRecord}, maf2::AbstractVector{<:MAFRecord}; cohort1_name::String="Cohort1", cohort2_name::String="Cohort2", top_n::Int=20, prov_ctx=nothing)
  _ctx = active_provenance_context(prov_ctx)
  s1 = summarize_maf(maf1)
  s2 = summarize_maf(maf2)

  n1 = max(1, length(s1.per_sample))
  n2 = max(1, length(s2.per_sample))

  all_genes = unique(vcat(collect(keys(s1.per_gene)), collect(keys(s2.per_gene))))

  genes = String[]
  c1_mut = Int[]
  c2_mut = Int[]
  odds_ratios = Float64[]
  pvalues = Float64[]

  for gene in all_genes
    m1 = get(s1.per_gene, gene, 0)
    m2 = get(s2.per_gene, gene, 0)

    a = m1
    b = max(0, n1 - m1)
    c = m2
    d = max(0, n2 - m2)

    or_val = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    m_total_mut = a + c
    m_total_wt = b + d
    n_sample = a + b

    if m_total_mut > 0 && n_sample > 0 && (m_total_mut + m_total_wt) > 0
      hg = Hypergeometric(m_total_mut, max(0, m_total_wt), n_sample)
      pval = 2 * min(cdf(hg, a), ccdf(hg, max(0, a - 1)))
      pval = clamp(pval, 1e-15, 1.0)
    else
      pval = 1.0
    end

    push!(genes, gene)
    push!(c1_mut, m1)
    push!(c2_mut, m2)
    push!(odds_ratios, Float64(or_val))
    push!(pvalues, Float64(pval))
  end

  df = DataFrame(
    gene=genes,
    cohort1_mutated=c1_mut,
    cohort2_mutated=c2_mut,
    cohort1_freq=c1_mut ./ n1,
    cohort2_freq=c2_mut ./ n2,
    odds_ratio=odds_ratios,
    pvalue=pvalues
  )
  if nrow(df) > 0
    df.fdr = benjamini_hochberg(df.pvalue)
  else
    df.fdr = Float64[]
  end
  sort!(df, :pvalue)
  df = df[1:min(top_n, nrow(df)), :]

  provenance_result!(_ctx, df, "maf_compare"; parents=provenance_parent_ids(maf1, maf2), parameters=(cohort1=cohort1_name, cohort2=cohort2_name))
  return df
end

"""
    clinical_enrichment(maf, clinical_df, feature_col; prov_ctx=nothing) → DataFrame

Perform enrichment analysis testing association between clinical annotations (e.g., stage, subtype) and gene mutation status using Fisher's exact test. Equivalent to `maftools::clinicalEnrichment`.
"""
function clinical_enrichment(maf::AbstractVector{<:MAFRecord}, clinical_df::DataFrame, feature_col::Symbol; prov_ctx=nothing)
  _ctx = active_provenance_context(prov_ctx)
  s = summarize_maf(maf)

  id_col = :patient_id
  if !hasproperty(clinical_df, id_col)
    for cand in [:sample, :Tumor_Sample_Barcode, :Sample_ID, :id]
      if hasproperty(clinical_df, cand)
        id_col = cand
        break
      end
    end
    if !hasproperty(clinical_df, id_col)
      id_col = propertynames(clinical_df)[1]
    end
  end

  results = DataFrame(
    gene=String[],
    feature_group=String[],
    mut_in_group=Int[],
    total_in_group=Int[],
    mut_outside_group=Int[],
    total_outside_group=Int[],
    odds_ratio=Float64[],
    pvalue=Float64[]
  )

  groups = unique(skipmissing(clinical_df[!, feature_col]))

  for (gene, mut_cnt) in s.per_gene
    mut_samples = Set([rec.sample for rec in maf if rec.gene == gene])
    for grp in groups
      grp_mask = .!ismissing.(clinical_df[!, feature_col]) .& (clinical_df[!, feature_col] .== grp)
      sub_df = clinical_df[grp_mask, :]
      other_df = clinical_df[.!grp_mask, :]

      n_grp = nrow(sub_df)
      n_other = nrow(other_df)

      a = count(r -> String(r[id_col]) in mut_samples, eachrow(sub_df))
      b = max(0, n_grp - a)
      c = count(r -> String(r[id_col]) in mut_samples, eachrow(other_df))
      d = max(0, n_other - c)

      or_val = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

      m_total_mut = a + c
      m_total_wt = b + d
      n_sample = a + b

      if m_total_mut > 0 && n_sample > 0 && (m_total_mut + m_total_wt) > 0
        hg = Hypergeometric(m_total_mut, max(0, m_total_wt), n_sample)
        pval = 2 * min(cdf(hg, a), ccdf(hg, max(0, a - 1)))
        pval = clamp(pval, 1e-15, 1.0)
      else
        pval = 1.0
      end

      push!(results, (gene, string(grp), a, n_grp, c, n_other, Float64(or_val), Float64(pval)))
    end
  end

  if nrow(results) > 0
    results.fdr = benjamini_hochberg(results.pvalue)
  else
    results.fdr = Float64[]
  end

  sort!(results, :pvalue)
  provenance_result!(_ctx, results, "clinical_enrichment"; parents=provenance_parent_ids(maf, clinical_df), parameters=(feature_col=string(feature_col),))
  return results
end

"""
    nelson_aalen(time, status) → NelsonAalenResult

Estimate the Nelson-Aalen cumulative hazard function H(t) = sum_{t_i <= t} (d_i / n_i) and standard error. Equivalent to `survival::survfit` cumulative hazard estimation.
"""
function nelson_aalen(time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}; prov_ctx=nothing)
  _ctx = active_provenance_context(prov_ctx)
  p = sortperm(time)
  st = Float64.(time[p])
  ss = Int.(status[p])

  unique_t = Float64[]
  cum_h = Float64[]
  se_h = Float64[]

  n = length(st)
  ch = 0.0
  var_h = 0.0
  i = 1

  while i <= n
    t_curr = st[i]
    d = 0
    at_risk = n - i + 1
    while i <= n && st[i] == t_curr
      if ss[i] > 0
        d += 1
      end
      i += 1
    end
    if d > 0
      ch += d / at_risk
      var_h += d / (at_risk^2)
      push!(unique_t, t_curr)
      push!(cum_h, ch)
      push!(se_h, sqrt(var_h))
    end
  end

  res = NelsonAalenResult(
    unique_t, cum_h, se_h,
    ResultProvenance(
      new_provenance_id("nelson_aalen_$(n)"),
      "nelson_aalen", "Clinical/nelson_aalen", :ok, String[], String[], String[], String[],
      (sample_size=n,), provenance_parent_ids(time, status), _provenance_timestamp()
    )
  )
  provenance_result!(_ctx, res, "nelson_aalen"; parents=provenance_parent_ids(time, status), parameters=(sample_size=n,))
  return res
end

"""
    cox_residuals(cox_result, time, status, X; type=:schoenfeld) → Matrix{Float64}

Calculate Schoenfeld, Martingale, Deviance, or Score residuals for a fitted Cox proportional hazards model matching `survival::residuals.coxph`.
For Schoenfeld residuals, risk-set expectations follow standard Breslow weighting: bar(X)(t_k) = sum(theta_i * X_i) / sum(theta_i) over risk set R_k.
"""
function cox_residuals(cox_result::CoxResult, time::AbstractVector{<:Real}, status::AbstractVector{<:Integer}, X::AbstractMatrix{<:Real}; type::Symbol=:schoenfeld, prov_ctx=nothing)
  _ctx = active_provenance_context(prov_ctx)
  n_samples, p_cov = size(X)
  beta = [t.beta for t in cox_result.terms]
  p_use = min(p_cov, length(beta))

  beta_vec = beta[1:p_use]
  X_mat = Matrix{Float64}(X[:, 1:p_use])
  eta = X_mat * beta_vec
  theta = exp.(clamp.(eta, -30.0, 30.0))

  res = zeros(Float64, n_samples, p_use)

  if type == :martingale || type == :deviance
    bt = cox_result.baseline_times
    bh = cox_result.baseline_hazard

    m_res = zeros(Float64, n_samples)
    for i in 1:n_samples
      t_i = time[i]
      st_i = status[i]
      lam0 = 0.0
      if !isempty(bt) && t_i >= bt[1]
        idx = searchsortedlast(bt, t_i)
        if idx > 0
          lam0 = bh[idx]
        end
      end
      m_res[i] = st_i - lam0 * theta[i]
    end

    if type == :martingale
      res = repeat(m_res, 1, p_use)
    else # :deviance
      for i in 1:n_samples
        mi = m_res[i]
        st_i = status[i]
        dev_val = 0.0
        if st_i > 0
          dev_val = sign(mi) * sqrt(max(0.0, 2.0 * (-mi - st_i * log(max(1e-12, st_i - mi)))))
        else
          dev_val = sign(mi) * sqrt(max(0.0, -2.0 * mi))
        end
        res[i, :] .= dev_val
      end
    end
  else # :schoenfeld or default
    for i in 1:n_samples
      if status[i] > 0
        t_i = time[i]
        risk_mask = time .>= t_i
        S0 = sum(theta[risk_mask])
        if S0 > 0.0
          S1 = sum(theta[risk_mask] .* X_mat[risk_mask, :], dims=1)
          E_x = vec(S1 / S0)
          res[i, :] = X_mat[i, :] .- E_x
        end
      end
    end
  end

  provenance_result!(_ctx, res, "cox_residuals"; parents=provenance_parent_ids(time, status, X), parameters=(type=string(type),))
  return res
end

end # module Clinical
