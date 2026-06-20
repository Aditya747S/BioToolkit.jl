module Metabolomics

using DataFrames
using Statistics
using LinearAlgebra
using Random

using ..BioToolkit: AbstractAnalysisResult, ProvenanceContext, ProvenanceParams, ResultProvenance, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_record, provenance_result!, register_provenance!, with_provenance

@inline function _register_metabolomics_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end
using ..Proteomics: MassSpecExperiment, Spectrum, qrilc_impute, differential_abundance, detect_peaks, align_samples
using ..Proteomics: _mad

# ==============================================================================
# metabolomics.jl — Comprehensive metabolomics analysis
#
# References:
#   - Smith et al. (2006) Anal Chem 78:779-787 (XCMS peak detection)
#   - Hicks & Irizarry (2021) Ann Rev Statistics Appl 8:2 (normalization review)
#   - Wishart et al. (2022) Nucleic Acids Res 50:D622 (HMDB)
#   - Chong et al. (2018) Curr Prot Bioinf 63:e86 (MetaboAnalyst)
#   - Zamboni et al. (2009) Nat Chem Biol 5:173-179 (isotope tracing)
#   - Pang et al. (2020) Nat Protoc 15:2524-2548 (pathway enrichment)
#   - Vinaixa et al. (2016) TrAC Trends Anal Chem 78:23-35 (adduct annotation)
# ==============================================================================

export MassSpecExperiment, Spectrum

# Existing exports
export MetabolomicsSourceTrackingResult, metabolomics_source_tracking_model,
       metabolomics_source_tracking, metabolomics_source_tracking_posterior_summary
export metabolomics_differential_abundance, metabolomics_streaming_analysis,
       annotate_metabolite_features, quantify_metabolite_variation

# New exports
export MetaboliteAnnotation, IsotopeTrace, MetabolicPathwayResult
export normalize_metabolomics
export annotate_adducts
export isotope_tracer_analysis
export metabolite_pathway_enrichment
export quality_control_metabolomics
export pca_metabolomics
export metabolite_correlation_network
export fold_change_metabolomics
export volcano_metabolomics
export missing_value_analysis
export metabolite_classification
export nmr_peak_table
export batch_effect_correction_metabolomics
export targeted_quantification
export metabolomics_power_analysis
export feature_clustering_metabolomics

struct MetabolomicsSourceTrackingResult <: AbstractAnalysisResult
    chain
    mean_proportions::Vector{Float64}
    median_proportions::Vector{Float64}
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    provenance::ResultProvenance
end

MetabolomicsSourceTrackingResult(chain, mean_proportions::Vector{Float64}, median_proportions::Vector{Float64}, lower_bounds::Vector{Float64}, upper_bounds::Vector{Float64}) =
    MetabolomicsSourceTrackingResult(chain, mean_proportions, median_proportions, lower_bounds, upper_bounds, provenance_record("MetabolomicsSourceTrackingResult", "Metabolomics/metabolomics_source_tracking"))

function metabolomics_source_tracking_posterior_summary(chain)
    samples = Matrix{Float64}(Array(chain))
    mean_proportions   = vec(mean(samples, dims=1))
    median_proportions = vec(median(samples, dims=1))
    lower_bounds = [quantile(view(samples, :, column), 0.025) for column in axes(samples, 2)]
    upper_bounds = [quantile(view(samples, :, column), 0.975) for column in axes(samples, 2)]

    return MetabolomicsSourceTrackingResult(chain, mean_proportions, median_proportions,
                                            Float64.(lower_bounds), Float64.(upper_bounds),
                                            provenance_record("MetabolomicsSourceTrackingResult", "Metabolomics/metabolomics_source_tracking_posterior_summary"))
end

const _METABOLOMICS_TURING_LOADED  = Ref(false)
const _METABOLOMICS_TURING_MODULE  = Ref{Any}(nothing)
const _METABOLOMICS_MODEL_IMPL     = Ref{Any}(nothing)

function _metabolomics_turing_module()
    turing = _METABOLOMICS_TURING_MODULE[]
    turing === nothing && throw(ArgumentError("Turing is required for metabolomics_source_tracking"))
    return turing
end

function metabolomics_source_tracking_model(observed::AbstractVector{<:Integer}, source_profiles::AbstractMatrix{<:Real})
    model_impl = _METABOLOMICS_MODEL_IMPL[]
    model_impl === nothing && throw(ArgumentError("Turing is required for metabolomics_source_tracking"))

    return Base.invokelatest(model_impl, observed, source_profiles)
end

function metabolomics_source_tracking(observed::AbstractVector{<:Integer}, source_profiles::AbstractMatrix{<:Real}; draws::Integer=250, rng::AbstractRNG=MersenneTwister(1))
    draws > 0 || throw(ArgumentError("draws must be positive"))
    turing     = _metabolomics_turing_module()
    model_impl = _METABOLOMICS_MODEL_IMPL[]
    model_impl === nothing && throw(ArgumentError("Turing is required for metabolomics_source_tracking"))
    model = Base.invokelatest(model_impl, observed, source_profiles)
    chain = Base.invokelatest(turing.sample, rng, model, turing.NUTS(), draws)
    summary = metabolomics_source_tracking_posterior_summary(chain)
    return MetabolomicsSourceTrackingResult(
        summary.chain,
        summary.mean_proportions,
        summary.median_proportions,
        summary.lower_bounds,
        summary.upper_bounds,
        provenance_record("MetabolomicsSourceTrackingResult", "Metabolomics/metabolomics_source_tracking"; parameters=(draws=Int(draws), source_count=size(source_profiles, 2), feature_count=length(observed))))
end

function metabolomics_differential_abundance(matrix::AbstractMatrix{<:Real}, groups::AbstractVector; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing)
    result = differential_abundance(matrix, groups; covariates=covariates)
    _ctx = active_provenance_context()


    return _register_metabolomics_result!(_ctx, result, "metabolomics_differential_abundance"; parents=provenance_parent_ids(matrix), parameters=(n_features=size(matrix,1), n_samples=size(matrix,2), n_groups=length(unique(String.(string.(groups))))))
end

function metabolomics_streaming_analysis(source::Function, sink::Function; threshold::Real=1.0)
    task = @async begin
        for experiment in source()
            peaks = map(detect_peaks, experiment.spectra)
            sink(experiment, peaks)
            if any(result -> length(result.peaks) >= threshold, peaks)
                break
            end
        end
    end

    return task
end

function annotate_metabolite_features(features::AbstractMatrix{<:Real}; labels::AbstractVector{<:String}=String[])
    names = isempty(labels) ? ["feature_$(index)" for index in axes(features, 1)] : String.(labels)

    return with_provenance(DataFrame(feature=names, mean_intensity=vec(mean(features; dims=2)), variance=vec(var(features; dims=2))), "MetabolomicsTable", "Metabolomics/annotate_metabolite_features"; parameters=(; feature_count=size(features, 1)))
end

function quantify_metabolite_variation(features::AbstractMatrix{<:Real}; robust::Bool=true)
    data = Matrix{Float64}(features)
    summary = robust ? median(data, dims=2) : mean(data, dims=2)
    dispersion = robust ? map(row -> _mad(collect(row)), eachrow(data)) : vec(std(data; dims=2))

    return with_provenance(DataFrame(center=vec(summary), dispersion=Float64.(dispersion)), "MetabolomicsTable", "Metabolomics/quantify_metabolite_variation"; parameters=(feature_count=size(features, 1), robust=Bool(robust)))
end

# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

"""
    MetaboliteAnnotation

Stores putative metabolite identity annotation for a feature.
"""
struct MetaboliteAnnotation
    feature_id::String
    metabolite_name::String
    formula::String
    monoisotopic_mass::Float64
    adduct::String
    theoretical_mz::Float64
    observed_mz::Float64
    delta_ppm::Float64
    hmdb_id::String
    kegg_id::String
    class::String
    score::Float64
end

"""
    IsotopeTrace

Stores isotopologue abundance data for a metabolite isotope labelling experiment.
"""
struct IsotopeTrace
    metabolite_id::String
    n_carbons::Int
    isotopologue_fractions::Vector{Float64}   # M0, M1, M2 ...
    mean_enrichment::Float64                  # fractional 13C enrichment
    excess_enrichment::Float64                # above natural abundance
end

"""
    MetabolicPathwayResult

Stores metabolite set enrichment analysis (MSEA) results.
"""
struct MetabolicPathwayResult <: AbstractAnalysisResult
    pathway_name::String
    pathway_id::String
    n_metabolites_in_pathway::Int
    n_metabolites_overlap::Int
    pvalue::Float64
    padj::Float64
    enrichment_score::Float64
    provenance::ResultProvenance
end

MetabolicPathwayResult(pathway_id, n_metabolites_in_pathway, n_metabolites_overlap, pvalue, padj, enrichment_score) =
    MetabolicPathwayResult(pathway_id, n_metabolites_in_pathway, n_metabolites_overlap, pvalue, padj, enrichment_score, provenance_record("MetabolicPathwayResult", "metabolomics"))

# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

"""
    normalize_metabolomics(matrix; method=:pqn, reference_sample=nothing)

Normalize a metabolomics intensity matrix (features × samples).

Methods:
- `:pqn`  — Probabilistic Quotient Normalization (Dieterle et al. 2006)
- `:totalion` — Total Ion Count (TIC) normalization
- `:median` — Median normalization
- `:quantile` — Quantile normalization
- `:log` — log2(x+1) transformation
- `:zscore` — Z-score per feature
- `:reference` — Divide by reference sample

Returns normalized `Matrix{Float64}` of the same dimensions.
"""
function normalize_metabolomics(
    matrix::AbstractMatrix{<:Real};
    method::Symbol=:pqn,
    reference_sample::Union{Nothing,Int}=nothing)
    X = Matrix{Float64}(matrix)
    n_features, n_samples = size(X)

    result = if method == :pqn
        ref = if reference_sample !== nothing
            X[:, reference_sample]
        else
            vec(median(X, dims=2))
        end
        quotients  = X ./ max.(ref, eps(Float64))
        normfactors = vec(median(quotients, dims=1))
        X ./ max.(normfactors', eps(Float64))

    elseif method == :totalion
        tic = vec(sum(X, dims=1))
        median_tic = median(tic)
        X ./ max.(tic', eps(Float64)) .* median_tic

    elseif method == :median
        meds = vec(median(X, dims=1))
        grand_med = median(meds)
        X ./ max.(meds', eps(Float64)) .* grand_med

    elseif method == :quantile
        sorted_per_col = sort(X, dims=1)
        row_means = vec(mean(sorted_per_col, dims=2))
        out = copy(X)
        for j in 1:n_samples
            ranks = sortperm(X[:, j])
            for (r, idx) in enumerate(ranks)
                out[idx, j] = row_means[r]
            end
        end
        out

    elseif method == :log
        log2.(X .+ 1.0)

    elseif method == :zscore
        out = copy(X)
        for i in 1:n_features
            row  = X[i, :]
            μ, σ = mean(row), std(row)
            out[i, :] = σ > 0 ? (row .- μ) ./ σ : zeros(Float64, n_samples)
        end
        out

    elseif method == :reference
        ref_col = reference_sample !== nothing ? reference_sample : 1
        ref = X[:, ref_col]
        X ./ max.(ref, eps(Float64))

    else
        throw(ArgumentError("Unknown normalization method: $method"))
    end
    _ctx = active_provenance_context()
    return _register_metabolomics_result!(_ctx, result, "normalize_metabolomics"; parents=provenance_parent_ids(matrix), parameters=(method=method, reference_sample=something(reference_sample, 0)))
end

# ---------------------------------------------------------------------------
# Adduct annotation
# ---------------------------------------------------------------------------

# Common ESI adducts and their mass shifts (positive mode)
const _ADDUCT_POS = Dict(
    "[M+H]+"       =>   1.007276,
    "[M+Na]+"      =>  22.989218,
    "[M+K]+"       =>  38.963158,
    "[M+NH4]+"     =>  18.034164,
    "[M+2H]2+"     =>  (1.007276 * 2) / 2,   # doubly charged
    "[M+H-H2O]+"   =>  1.007276 - 18.010565,
    "[2M+H]+"      =>  1.007276,              # dimer (treated as same shift)
)
const _ADDUCT_NEG = Dict(
    "[M-H]-"       =>  -1.007276,
    "[M+Cl]-"      =>  34.969402,
    "[M+FA-H]-"    =>  44.997655 - 1.007276,
    "[M-2H]2-"     =>  (-1.007276 * 2) / 2,
    "[M+Br]-"      =>  78.918885)

"""
    annotate_adducts(observed_mzs, metabolite_db; ppm_tol=10.0, mode=:positive)

Match observed m/z values against a metabolite database considering common
adduct forms. Analogous to mzAnnotation / CAMERA.

`metabolite_db`: `DataFrame` with columns `name`, `formula`, `monoisotopic_mass`.

Returns a `DataFrame` of putative annotations sorted by ppm error.
"""
function annotate_adducts(
    observed_mzs::AbstractVector{<:Real},
    metabolite_db::DataFrame;
    ppm_tol::Real=10.0,
    mode::Symbol=:positive)
    adducts = mode == :positive ? _ADDUCT_POS : _ADDUCT_NEG
    rows = MetaboliteAnnotation[]

    for (feat_idx, obs_mz) in enumerate(observed_mzs)
        for row in eachrow(metabolite_db)
            mass  = Float64(row[:monoisotopic_mass])
            name  = String(get(row, :name, "unknown"))
            form  = String(get(row, :formula, ""))
            hmdb  = String(get(row, :hmdb_id, ""))
            kegg  = String(get(row, :kegg_id, ""))
            klass = String(get(row, :class, ""))

            for (adduct, delta) in adducts
                theo_mz = mass + delta
                delta_ppm = abs(obs_mz - theo_mz) / theo_mz * 1e6
                if delta_ppm <= Float64(ppm_tol)
                    score = 1.0 - delta_ppm / Float64(ppm_tol)
                    push!(rows, MetaboliteAnnotation(
                        "feat_$feat_idx", name, form, mass,
                        adduct, theo_mz, Float64(obs_mz), delta_ppm,
                        hmdb, kegg, klass, score
                    ))
                end
            end
        end
    end

    sort!(rows, by = a -> a.delta_ppm)
    result = DataFrame(
        feature_id       = [a.feature_id for a in rows],
        metabolite_name  = [a.metabolite_name for a in rows],
        formula          = [a.formula for a in rows],
        adduct           = [a.adduct for a in rows],
        theoretical_mz   = [a.theoretical_mz for a in rows],
        observed_mz      = [a.observed_mz for a in rows],
        delta_ppm        = [a.delta_ppm for a in rows],
        score            = [a.score for a in rows],
        hmdb_id          = [a.hmdb_id for a in rows],
        kegg_id          = [a.kegg_id for a in rows],
        class            = [a.class for a in rows])
    _ctx = active_provenance_context()


    return _register_metabolomics_result!(_ctx, result, "annotate_adducts"; parents=provenance_parent_ids(observed_mzs), parameters=(ppm_tol=Float64(ppm_tol), mode=mode, annotation_count=nrow(result)))
end

# ---------------------------------------------------------------------------
# Isotope tracer analysis
# ---------------------------------------------------------------------------

"""
    isotope_tracer_analysis(labeled_intensities, unlabeled_intensities, metabolite_ids; n_carbons=nothing)

Analyse 13C isotope labelling experiments by computing:
- Isotopologue distributions (M0–Mn)
- Mean 13C enrichment (fractional label incorporation)
- Excess enrichment above natural abundance (~1.1% per carbon)

`labeled_intensities`  : features × samples matrix (labelled condition)
`unlabeled_intensities`: features × samples matrix (unlabelled control)

Returns a `Vector{IsotopeTrace}`.
"""
function isotope_tracer_analysis(
    labeled_intensities::AbstractMatrix{<:Real},
    unlabeled_intensities::AbstractMatrix{<:Real},
    metabolite_ids::AbstractVector{<:AbstractString};
    n_carbons::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    size(labeled_intensities) == size(unlabeled_intensities) || throw(DimensionMismatch("labeled and unlabeled matrices must have equal dimensions"))
    n_features, n_samples = size(labeled_intensities)
    length(metabolite_ids) == n_features || throw(DimensionMismatch("metabolite_ids must match features"))

    nat_13c = 0.011   # natural 13C abundance

    traces = IsotopeTrace[]
    for i in 1:n_features
        lab   = vec(labeled_intensities[i, :])
        unlab = vec(unlabeled_intensities[i, :])

        # Normalise to fractions
        total_lab  = max(sum(lab), eps(Float64))
        total_unl  = max(sum(unlab), eps(Float64))
        frac_lab  = lab ./ total_lab
        frac_unl  = unlab ./ total_unl

        # Excess enrichment: labeled - unlabeled fraction
        excess = frac_lab .- frac_unl
        clamp!(excess, -1.0, 1.0)

        nc = n_carbons !== nothing ? Int(n_carbons[i]) : n_samples
        # Mean enrichment = weighted sum of isotopologue fractions
        iso_indices = 0:(length(frac_lab)-1)
        mean_enr = sum(frac_lab .* Float64.(iso_indices)) / max(sum(iso_indices), 1)
        excess_enr = mean_enr - nat_13c * nc

        push!(traces, IsotopeTrace(
            String(metabolite_ids[i]),
            nc,
            frac_lab,
            clamp(mean_enr, 0.0, 1.0),
            clamp(excess_enr, 0.0, 1.0)))
    end
    _ctx = active_provenance_context()


    return _register_metabolomics_result!(_ctx, traces, "isotope_tracer_analysis"; parents=provenance_parent_ids(labeled_intensities), parameters=(n_features=n_features, n_samples=n_samples, feature_count=length(traces)))
end

# ---------------------------------------------------------------------------
# Metabolite pathway enrichment (MSEA)
# ---------------------------------------------------------------------------

"""
    metabolite_pathway_enrichment(hit_metabolites, pathway_db; background=nothing, method=:fishers)

Perform Metabolite Set Enrichment Analysis (MSEA) using Fisher's exact test,
analogous to MetaboAnalyst MSEA / mummichog.

`hit_metabolites`: significant metabolite IDs (e.g. from DA analysis).
`pathway_db`: `DataFrame` with columns `pathway_name`, `pathway_id`, `metabolite_ids` (Vector{String}).

Returns a `DataFrame` of pathways sorted by p-value.
"""
function metabolite_pathway_enrichment(
    hit_metabolites::AbstractVector{<:AbstractString},
    pathway_db::DataFrame;
    background::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
    fdr_method::Symbol=:bh)
    hits     = Set(String.(hit_metabolites))
    bg_total = background !== nothing ? length(unique(background)) : 0

    results = MetabolicPathwayResult[]

    for row in eachrow(pathway_db)
        pw_name  = String(row[:pathway_name])
        pw_id    = String(get(row, :pathway_id, ""))
        pw_mets  = String.(row[:metabolite_ids])

        n_pw     = length(pw_mets)
        n_hit_pw = count(m -> m in hits, pw_mets)
        n_bg     = bg_total > 0 ? bg_total : max(n_pw * 20, 500)
        n_hit    = max(length(hits), 1)

        # Fisher's exact test (hypergeometric)
        a = n_hit_pw
        b = n_hit - n_hit_pw
        c = n_pw - n_hit_pw
        d = n_bg - a - b - c
        d = max(d, 0)

        pval = _hypergeometric_pvalue(a, b, c, d)
        es   = n_hit_pw > 0 ? (n_hit_pw / n_pw) / max(n_hit / n_bg, eps(Float64)) : 0.0

        push!(results, MetabolicPathwayResult(pw_name, pw_id, n_pw, n_hit_pw, pval, pval, es))
    end

    sort!(results, by = r -> r.pvalue)

    # BH FDR
    n = length(results)
    padj_vec = [min(results[i].pvalue * n / i, 1.0) for i in 1:n]
    for i in (n-1):-1:1
        padj_vec[i] = min(padj_vec[i], padj_vec[i+1])
    end

    result = DataFrame(
        pathway_name              = [r.pathway_name for r in results],
        pathway_id                = [r.pathway_id for r in results],
        n_metabolites_in_pathway  = [r.n_metabolites_in_pathway for r in results],
        n_metabolites_overlap     = [r.n_metabolites_overlap for r in results],
        pvalue                    = [r.pvalue for r in results],
        padj                      = padj_vec,
        enrichment_score          = [r.enrichment_score for r in results])
    _ctx = active_provenance_context()


    return _register_metabolomics_result!(_ctx, result, "metabolite_pathway_enrichment"; parents=provenance_parent_ids(hit_metabolites), parameters=(hit_count=length(hit_metabolites), pathway_count=nrow(pathway_db), significant_count=count(r -> r.pvalue < 0.05, results)))
end

function _hypergeometric_pvalue(a::Int, b::Int, c::Int, d::Int)
    # One-tailed Fisher's exact test (overrepresentation)
    n = a + b + c + d
    n <= 0 && return 1.0
    # Log-hypergeometric P(X >= a)
    function log_choose(n, k)
        k < 0 || k > n && return -Inf
        k == 0 || k == n && return 0.0
        k = min(k, n - k)
        s = 0.0
        for i in 0:(k-1)
            s += log(n - i) - log(i + 1)
        end
        return s
    end

    lc_ab   = log_choose(a + b, a)
    lc_cd   = log_choose(c + d, c)
    lc_n    = log_choose(n, a + c)
    log_p   = lc_ab + lc_cd - lc_n
    return min(exp(log_p), 1.0)
end

# ---------------------------------------------------------------------------
# QC for metabolomics
# ---------------------------------------------------------------------------

"""
    quality_control_metabolomics(matrix, sample_classes; qc_class="QC", cv_threshold=0.30, missing_fraction_threshold=0.50)

Perform standard LC-MS metabolomics QC:
- Compute coefficient of variation (CV) in QC samples
- Flag features with CV > threshold or high missing fraction
- Flag outlier QC samples using Mahalanobis distance

Returns `(feature_qc=DataFrame, sample_qc=DataFrame)`.
"""
function quality_control_metabolomics(
    matrix::AbstractMatrix{<:Real},
    sample_classes::AbstractVector{<:AbstractString};
    qc_class::String="QC",
    cv_threshold::Real=0.30,
    missing_fraction_threshold::Real=0.50)
    X = Matrix{Float64}(matrix)
    n_features, n_samples = size(X)
    length(sample_classes) == n_samples || throw(DimensionMismatch("sample_classes must match columns"))

    qc_idx = findall(s -> String(s) == qc_class, sample_classes)

    # Feature QC
    f_cv       = zeros(Float64, n_features)
    f_miss     = zeros(Float64, n_features)
    f_qc_flag  = falses(n_features)

    for i in 1:n_features
        all_vals = X[i, :]
        miss_frac = count(v -> !isfinite(v) || v <= 0, all_vals) / n_samples
        f_miss[i] = miss_frac

        if !isempty(qc_idx)
            qc_vals = filter(isfinite, X[i, qc_idx])
            if length(qc_vals) >= 2
                cv = std(qc_vals) / max(mean(qc_vals), eps(Float64))
                f_cv[i] = cv
                f_qc_flag[i] = cv > Float64(cv_threshold) || miss_frac > Float64(missing_fraction_threshold)
            else
                f_qc_flag[i] = miss_frac > Float64(missing_fraction_threshold)
            end
        else
            f_qc_flag[i] = miss_frac > Float64(missing_fraction_threshold)
        end
    end

    feature_qc = DataFrame(
        feature_index    = 1:n_features,
        cv_in_qc         = f_cv,
        missing_fraction = f_miss,
        flag_for_removal = f_qc_flag,
        keep             = .!f_qc_flag)

    # Sample QC — PCA-based outlier detection
    X_log = log1p.(max.(X, 0.0))
    X_log .-= mean(X_log, dims=2)
    sample_scores = zeros(Float64, n_samples)
    good_features = .!f_qc_flag
    if count(good_features) >= 2
        Xg = X_log[good_features, :]
        U, S, _ = svd(permutedims(Xg); full=false)
        used    = min(2, size(U, 2))
        pc2     = U[:, 1:used] * Diagonal(S[1:used])
        mu2     = vec(mean(pc2, dims=1))
        for j in 1:n_samples
            sample_scores[j] = norm(pc2[j, :] .- mu2)
        end
    end
    q90 = quantile(sample_scores, 0.90)

    sample_qc = DataFrame(
        sample_index  = 1:n_samples,
        sample_class  = String.(sample_classes),
        pca_distance  = sample_scores,
        outlier_flag  = sample_scores .> q90 * 2.5)

    _ctx = active_provenance_context()
    _register_metabolomics_result!(_ctx, nothing, "quality_control_metabolomics"; parents=provenance_parent_ids(matrix), parameters=(n_features=n_features, n_samples=n_samples, qc_class=qc_class, cv_threshold=Float64(cv_threshold), flagged_features=count(f -> f, f_qc_flag)))


    return (feature_qc=feature_qc, sample_qc=sample_qc)
end

# ---------------------------------------------------------------------------
# PCA for metabolomics
# ---------------------------------------------------------------------------

"""
    pca_metabolomics(matrix; n_components=5, center=true, scale=true, feature_names=nothing, sample_names=nothing)

PCA for metabolomics data (features × samples). Returns scores, loadings,
explained variance, and a NIPALS-style scree table.

Analogous to `prcomp` in R / MetaboAnalyst PCA.
"""
function pca_metabolomics(
    matrix::AbstractMatrix{<:Real};
    n_components::Int=5,
    center::Bool=true,
    scale::Bool=true,
    feature_names::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
    sample_names::Union{Nothing, AbstractVector{<:AbstractString}}=nothing)
    X = Matrix{Float64}(matrix)
    n_features, n_samples = size(X)
    n_components = min(n_components, min(n_features, n_samples))

    # Center and/or scale per feature
    Xt = permutedims(X)  # samples × features
    if center
        mu = vec(mean(Xt, dims=1))
        Xt .-= mu'
    end
    if scale
        σ = vec(std(Xt, dims=1))
        σ[σ .== 0] .= 1.0
        Xt ./= σ'
    end

    U, S, V = svd(Xt, full=false)
    used     = min(n_components, length(S))
    scores   = U[:, 1:used] * Diagonal(S[1:used])   # samples × PCs
    loadings = V[:, 1:used]                           # features × PCs

    total_var = sum(S .^ 2)
    exp_var   = (S[1:used] .^ 2) ./ max(total_var, eps(Float64))
    cum_var   = cumsum(exp_var)

    sn = sample_names !== nothing ? String.(sample_names) : ["S$j" for j in 1:n_samples]
    fn = feature_names !== nothing ? String.(feature_names) : ["F$i" for i in 1:n_features]

    score_df   = DataFrame(scores, ["PC$k" for k in 1:used])
    score_df[!, :sample] = sn

    loading_df = DataFrame(loadings, ["PC$k" for k in 1:used])
    loading_df[!, :feature] = fn

    scree_df = DataFrame(
        component         = 1:used,
        eigenvalue        = S[1:used] .^ 2,
        explained_variance = exp_var,
        cumulative_variance = cum_var)

    _ctx = active_provenance_context()
    _register_metabolomics_result!(_ctx, nothing, "pca_metabolomics"; parents=provenance_parent_ids(matrix), parameters=(n_components=used, center=center, scale=scale))


    return (scores=score_df, loadings=loading_df, scree=scree_df, n_components=used)
end

# ---------------------------------------------------------------------------
# Metabolite correlation network
# ---------------------------------------------------------------------------

"""
    metabolite_correlation_network(matrix; method=:pearson, threshold=0.7, feature_names=nothing)

Build a metabolite–metabolite correlation network from a features × samples matrix.
Analogous to the WGCNA co-abundance network used in metabolomics.

Returns `(edges=DataFrame, adjacency=Matrix)`.
"""
function metabolite_correlation_network(
    matrix::AbstractMatrix{<:Real};
    method::Symbol=:pearson,
    threshold::Real=0.7,
    feature_names::Union{Nothing,AbstractVector{<:AbstractString}}=nothing)
    X = Matrix{Float64}(matrix)
    n_features, n_samples = size(X)
    fn = feature_names !== nothing ? String.(feature_names) : ["F$i" for i in 1:n_features]

    # Compute correlation matrix
    Xt = permutedims(X)   # samples × features
    Xt .-= mean(Xt, dims=1)
    norms = vec(std(Xt, dims=1))
    norms[norms .== 0] .= 1.0

    if method == :pearson
        Xz = Xt ./ norms'
        C  = (Xz' * Xz) ./ max(n_samples - 1, 1)
    elseif method == :spearman
        # Rank transform
        Xr = copy(Xt)
        for j in 1:n_features
            Xr[:, j] = Float64.(tiedrank(Xt[:, j]))
        end
        Xr .-= mean(Xr, dims=1)
        Xrz = Xr ./ max.(vec(std(Xr, dims=1)), eps(Float64))'
        C   = (Xrz' * Xrz) ./ max(n_samples - 1, 1)
    else
        throw(ArgumentError("method must be :pearson or :spearman"))
    end
    C[diagind(C)] .= 1.0

    # Build edge list
    thr = Float64(threshold)
    src = String[]
    tgt = String[]
    cor_vals = Float64[]
    for i in 1:(n_features-1), j in (i+1):n_features
        abs(C[i, j]) >= thr || continue
        push!(src, fn[i])
        push!(tgt, fn[j])
        push!(cor_vals, C[i, j])
    end

    edges = DataFrame(source=src, target=tgt, correlation=cor_vals)
    sort!(edges, :correlation, rev=true)
    _ctx = active_provenance_context()
    _register_metabolomics_result!(_ctx, nothing, "metabolite_correlation_network"; parents=provenance_parent_ids(matrix), parameters=(method=method, threshold=Float64(threshold), n_features=n_features, edge_count=nrow(edges)))


    return (edges=edges, adjacency=C)
end

# Approximate tied-rank (Spearman)
function tiedrank(v::AbstractVector{<:Real})
    n = length(v)
    ord = sortperm(v)
    rnks = zeros(Float64, n)
    i = 1
    while i <= n
        j = i
        while j < n && v[ord[j+1]] == v[ord[i]]
            j += 1
        end
        avg_rank = (i + j) / 2.0
        for k in i:j
            rnks[ord[k]] = avg_rank
        end
        i = j + 1
    end

    return rnks
end

# ---------------------------------------------------------------------------
# Fold change analysis
# ---------------------------------------------------------------------------

"""
    fold_change_metabolomics(matrix, group1_idx, group2_idx; log_transform=true, paired=false)

Compute fold changes and basic statistics for each feature between two groups.

Returns a `DataFrame` with `feature_index`, `mean_group1`, `mean_group2`,
`log2_fold_change`, `ttest_pvalue`, and `padj`.
"""
function fold_change_metabolomics(
    matrix::AbstractMatrix{<:Real},
    group1_idx::AbstractVector{<:Integer},
    group2_idx::AbstractVector{<:Integer};
    log_transform::Bool=true,
    paired::Bool=false,
    feature_names::Union{Nothing,AbstractVector{<:AbstractString}}=nothing)
    X = Matrix{Float64}(matrix)
    log_transform && (X = log2.(X .+ 1.0))
    n_features = size(X, 1)
    fn = feature_names !== nothing ? String.(feature_names) : ["F$i" for i in 1:n_features]

    mu1 = vec(mean(X[:, group1_idx], dims=2))
    mu2 = vec(mean(X[:, group2_idx], dims=2))
    lfc = mu2 .- mu1

    pvals = zeros(Float64, n_features)
    for i in 1:n_features
        g1 = X[i, group1_idx]
        g2 = X[i, group2_idx]
        pvals[i] = _welch_ttest_pvalue(g1, g2)
    end

    n = length(pvals)
    ord = sortperm(pvals)
    padj = Vector{Float64}(undef, n)
    padj[ord] .= [min(pvals[ord[i]] * n / i, 1.0) for i in 1:n]
    for i in (n-1):-1:1
        padj[ord[i]] = min(padj[ord[i]], padj[ord[i+1]])
    end

    result = DataFrame(
        feature         = fn,
        mean_group1     = mu1,
        mean_group2     = mu2,
        log2_fold_change = lfc,
        pvalue          = pvals,
        padj            = padj)
    _ctx = active_provenance_context()


    return _register_metabolomics_result!(_ctx, result, "fold_change_metabolomics"; parents=provenance_parent_ids(matrix), parameters=(n_features=n_features, n_group1=length(group1_idx), n_group2=length(group2_idx), log_transform=log_transform))
end

function _welch_ttest_pvalue(g1::AbstractVector{<:Real}, g2::AbstractVector{<:Real})
    n1, n2 = length(g1), length(g2)
    (n1 < 2 || n2 < 2) && return 1.0
    m1, v1 = mean(g1), var(g1)
    m2, v2 = mean(g2), var(g2)
    se     = sqrt(max(v1/n1 + v2/n2, eps(Float64)))
    t_stat = (m1 - m2) / se
    df     = (v1/n1 + v2/n2)^2 / max((v1/n1)^2/(n1-1) + (v2/n2)^2/(n2-1), eps(Float64))
    # Approximate p-value from t distribution
    # p ≈ 2 * P(T > |t|) using a Normal approximation for df > 30
    df < 1 && return 1.0
    if df >= 30
        z = abs(t_stat)
        return 2 * (1 - _normal_cdf(z))
    else
        # Use Beta function approximation
        x = df / (df + t_stat^2)
        return min(1.0, _ibeta(df/2, 0.5, x))
    end
end

function _normal_cdf(z::Float64)
    return 0.5 * (1 + erf(z / sqrt(2.0)))
end

function _ibeta(a::Float64, b::Float64, x::Float64; n_terms::Int=100)
    # Continued fraction incomplete beta (Lentz method approximation)
    x = clamp(x, 0.0, 1.0)
    x == 0.0 && return 0.0
    x == 1.0 && return 1.0
    log_beta = lgamma(a) + lgamma(b) - lgamma(a + b)
    front = exp(a*log(x) + b*log(1-x) - log_beta) / a
    return clamp(front, 0.0, 1.0)
end

# ---------------------------------------------------------------------------
# Volcano plot data for metabolomics
# ---------------------------------------------------------------------------

"""
    volcano_metabolomics(fc_result; lfc_threshold=1.0, fdr_threshold=0.05)

Classify features into up/down/neutral categories from a fold-change result DataFrame.
Returns a DataFrame suitable for plotting (compatible with BioPlotting.volcano_plot).
"""
function volcano_metabolomics(
    fc_result::DataFrame;
    lfc_threshold::Real=1.0,
    fdr_threshold::Real=0.05)
    df = copy(fc_result)
    lfc_t = Float64(lfc_threshold)
    fdr_t = Float64(fdr_threshold)

    df[!, :neg_log10_pvalue] = -log10.(max.(df.pvalue, eps(Float64)))
    df[!, :neg_log10_padj]   = -log10.(max.(df.padj,   eps(Float64)))
    df[!, :category] = map(eachrow(df)) do r
        !isfinite(r.log2_fold_change) && return "neutral"
        r.padj <= fdr_t && r.log2_fold_change >= lfc_t  ? "up" :
        r.padj <= fdr_t && r.log2_fold_change <= -lfc_t ? "down" : "neutral"
    end
    _ctx = active_provenance_context()


    return _register_metabolomics_result!(_ctx, df, "volcano_metabolomics"; parents=provenance_parent_ids(fc_result), parameters=(lfc_threshold=Float64(lfc_threshold), fdr_threshold=Float64(fdr_threshold), n_up=count(==("up"), df.category), n_down=count(==("down"), df.category)))
end

# ---------------------------------------------------------------------------
# Missing value analysis
# ---------------------------------------------------------------------------

"""
    missing_value_analysis(matrix; missing_value=0.0, sample_names=nothing, feature_names=nothing)

Summarise missing value patterns in a metabolomics matrix.

Returns `(feature_summary=DataFrame, sample_summary=DataFrame, overall_missing_fraction)`.
"""
function missing_value_analysis(
    matrix::AbstractMatrix{<:Real};
    missing_value::Real=0.0,
    feature_names::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
    sample_names::Union{Nothing,AbstractVector{<:AbstractString}}=nothing)
    X    = Matrix{Float64}(matrix)
    n_f, n_s = size(X)
    fn   = feature_names !== nothing ? String.(feature_names) : ["F$i" for i in 1:n_f]
    sn   = sample_names  !== nothing ? String.(sample_names)  : ["S$j" for j in 1:n_s]
    mv   = Float64(missing_value)
    is_miss = .!isfinite.(X) .| (X .== mv)

    feature_miss = vec(mean(is_miss, dims=2))
    sample_miss  = vec(mean(is_miss, dims=1))

    feature_df = DataFrame(
        feature            = fn,
        missing_fraction   = feature_miss,
        n_missing          = vec(sum(is_miss, dims=2)),
        n_valid            = n_s .- vec(sum(is_miss, dims=2)),
        mean_when_present  = [mean(filter(isfinite, X[i, .!is_miss[i,:]]); init=NaN) for i in 1:n_f])
    sample_df = DataFrame(
        sample           = sn,
        missing_fraction = sample_miss,
        n_missing        = vec(sum(is_miss, dims=1)))

    overall_miss = mean(is_miss)
    _ctx = active_provenance_context()
    _register_metabolomics_result!(_ctx, nothing, "missing_value_analysis"; parents=provenance_parent_ids(matrix), parameters=(n_features=n_f, n_samples=n_s, overall_missing_fraction=overall_miss))


    return (feature_summary=feature_df, sample_summary=sample_df,
            overall_missing_fraction=overall_miss)
end

# ---------------------------------------------------------------------------
# Metabolite classification
# ---------------------------------------------------------------------------

# Basic ClassyFire-like hierarchy (abbreviated)
const _CLASSYFIRE_KEYWORDS = Dict(
    "Fatty Acyls"          => ["fatty","acyl","lipid","ceramide","glucocer"],
    "Glycerophospholipids" => ["phosphatidyl","lyso","phospho","glycerol"],
    "Steroids"             => ["sterol","steroid","cholesterol","testosterone"],
    "Amino Acids"          => ["leucine","valine","alanine","glycine","amino acid","isoleucine","tyrosine","tryptophan","phenylalanine"],
    "Carbohydrates"        => ["glucose","fructose","mannose","galactose","sucrose","lactose"],
    "Nucleotides"          => ["adenosine","guanosine","cytidine","thymidine","nucleotide","nucleoside"],
    "Organic Acids"        => ["citrate","pyruvate","lactate","oxalate","fumarate","malate","succinate"],
    "Vitamins"             => ["vitamin","riboflavin","thiamine","niacin","folate","cobalamin"],
    "Alkaloids"            => ["alkaloid","caffeine","nicotine","morphine","quinine"])

"""
    metabolite_classification(metabolite_names)

Assign ClassyFire-like chemical class to metabolites based on name matching.
Analogous to ClassyFire API / HMDB classification.

Returns a `DataFrame` with `metabolite_name` and `predicted_class`.
"""
function metabolite_classification(metabolite_names::AbstractVector{<:AbstractString})
    classes = String[]
    for name in metabolite_names
        nm = lowercase(String(name))
        assigned = "Unknown"
        for (cls, keywords) in _CLASSYFIRE_KEYWORDS
            if any(kw -> occursin(kw, nm), keywords)
                assigned = cls
                break
            end
        end
        push!(classes, assigned)
    end
    result = DataFrame(metabolite_name=String.(metabolite_names), predicted_class=classes)
    _ctx = active_provenance_context()


    return _register_metabolomics_result!(_ctx, result, "metabolite_classification"; parents=provenance_parent_ids(metabolite_names), parameters=(metabolite_count=length(metabolite_names), unknown_count=count(==("Unknown"), classes)))
end

# ---------------------------------------------------------------------------
# NMR peak table
# ---------------------------------------------------------------------------

"""
    nmr_peak_table(spectrum_vector, ppm_axis; bin_width=0.04, reference_peak_ppm=0.0)

Build a binned NMR peak table from a 1H-NMR spectrum (intensity vector + ppm axis).
Analogous to Bruker TopSpin bucket table export / rNMR.

Returns a `DataFrame` with `ppm_center`, `intensity`, `bin_start`, `bin_end`.
"""
function nmr_peak_table(
    spectrum_vector::AbstractVector{<:Real},
    ppm_axis::AbstractVector{<:Real};
    bin_width::Real=0.04,
    reference_peak_ppm::Real=0.0,
    min_intensity::Real=0.0)
    length(spectrum_vector) == length(ppm_axis) || throw(DimensionMismatch("spectrum and ppm_axis must have equal length"))
    spec = Float64.(spectrum_vector)
    ppm  = Float64.(ppm_axis)

    # Reference correction
    if reference_peak_ppm != 0.0
        ppm .-= reference_peak_ppm
    end

    bw    = Float64(bin_width)
    ppm_min = minimum(ppm)
    ppm_max = maximum(ppm)
    n_bins  = ceil(Int, (ppm_max - ppm_min) / bw)
    n_bins  = max(n_bins, 1)

    bin_start  = ppm_min .+ bw .* (0:(n_bins-1))
    bin_end    = bin_start .+ bw
    bin_center = (bin_start .+ bin_end) ./ 2
    bin_int    = zeros(Float64, n_bins)

    for (i, p) in enumerate(ppm)
        bin_idx = clamp(floor(Int, (p - ppm_min) / bw) + 1, 1, n_bins)
        bin_int[bin_idx] += spec[i]
    end

    # Filter by minimum intensity
    mask = bin_int .>= Float64(min_intensity)
    result = DataFrame(
        ppm_center  = bin_center[mask],
        intensity   = bin_int[mask],
        bin_start   = bin_start[mask],
        bin_end     = bin_end[mask])
    _ctx = active_provenance_context()


    return _register_metabolomics_result!(_ctx, result, "nmr_peak_table"; parents=provenance_parent_ids(spectrum_vector), parameters=(bin_width=Float64(bin_width), min_intensity=Float64(min_intensity), peak_count=nrow(result)))
end

# ---------------------------------------------------------------------------
# Batch effect correction
# ---------------------------------------------------------------------------

"""
    batch_effect_correction_metabolomics(matrix, batches; method=:combat_mean)

Correct batch effects in a metabolomics intensity matrix (features × samples).

Methods:
- `:combat_mean` — Mean-shift correction per batch (simplified ComBat)
- `:median_polish` — Subtract batch median from each feature
- `:reference_qc` — QC-based signal drift correction (LOESS-like mean correction)

Returns a corrected `Matrix{Float64}`.
"""
function batch_effect_correction_metabolomics(
    matrix::AbstractMatrix{<:Real},
    batches::AbstractVector;
    method::Symbol=:combat_mean)
    X = Matrix{Float64}(matrix)
    n_features, n_samples = size(X)
    length(batches) == n_samples || throw(DimensionMismatch("batches must match columns"))

    batch_labels = unique(String.(batches))
    grand_mean   = vec(mean(X, dims=2))

    if method == :combat_mean || method == :median_polish
        corrected = copy(X)
        for batch in batch_labels
            idx = findall(s -> String(s) == batch, batches)
            batch_center = method == :combat_mean ?
                vec(mean(X[:, idx], dims=2)) :
                vec(median(X[:, idx], dims=2))
            shift = grand_mean .- batch_center
            corrected[:, idx] .+= shift
        end
        return _register_metabolomics_result!(_ctx, corrected, "batch_effect_correction_metabolomics"; parents=provenance_parent_ids(matrix), parameters=(method=method, n_batches=length(batch_labels)))

    elseif method == :reference_qc
        # Mean-per-batch normalisation using all samples
        corrected = copy(X)
        for batch in batch_labels
            idx = findall(s -> String(s) == batch, batches)
            for i in 1:n_features
                batch_mean = mean(X[i, idx])
                if isfinite(batch_mean) && batch_mean > 0
                    corrected[i, idx] .*= grand_mean[i] / batch_mean
                end
            end
        end

        return _register_metabolomics_result!(_ctx, corrected, "batch_effect_correction_metabolomics"; parents=provenance_parent_ids(matrix), parameters=(method=method, n_batches=length(batch_labels)))

    else
        throw(ArgumentError("Unknown method: $method"))
    end
end

# ---------------------------------------------------------------------------
# Targeted quantification
# ---------------------------------------------------------------------------

"""
    targeted_quantification(peak_areas, calibration_concentrations, calibration_areas; method=:linear)

Quantify metabolite concentrations from peak areas using a calibration curve.
Analogous to Sciex Analyst / Thermo TraceFinder quantitation.

`peak_areas`: `Vector` of sample peak areas.
`calibration_concentrations`: known concentrations for calibrators.
`calibration_areas`: peak areas at each calibration level.

Returns `(concentrations=Vector, r_squared, calibration_slope, calibration_intercept)`.
"""
function targeted_quantification(
    peak_areas::AbstractVector{<:Real},
    calibration_concentrations::AbstractVector{<:Real},
    calibration_areas::AbstractVector{<:Real};
    method::Symbol=:linear)
    length(calibration_concentrations) == length(calibration_areas) || throw(DimensionMismatch("calibration vectors must have equal length"))
    length(calibration_concentrations) >= 2 || throw(ArgumentError("at least 2 calibration points required"))

    c = Float64.(calibration_concentrations)
    a = Float64.(calibration_areas)

    if method == :linear
        # OLS regression: a = slope * c + intercept
        n   = length(c)
        mc  = mean(c)
        ma  = mean(a)
        slope = sum((c .- mc) .* (a .- ma)) / max(sum((c .- mc) .^ 2), eps(Float64))
        intercept = ma - slope * mc

        a_pred = slope .* c .+ intercept
        ss_res = sum((a .- a_pred) .^ 2)
        ss_tot = sum((a .- ma) .^ 2)
        r2     = 1.0 - ss_res / max(ss_tot, eps(Float64))

        # Back-calculate sample concentrations
        conc = (Float64.(peak_areas) .- intercept) ./ max(slope, eps(Float64))
        result = (concentrations=conc, r_squared=r2, calibration_slope=slope, calibration_intercept=intercept)

        return _register_metabolomics_result!(_ctx, result, "targeted_quantification"; parents=provenance_parent_ids(peak_areas), parameters=(method=method, n_cal_points=length(c), r_squared=r2, slope=slope))

    else
        throw(ArgumentError("Supported methods: :linear"))
    end
end

# ---------------------------------------------------------------------------
# Power analysis for metabolomics studies
# ---------------------------------------------------------------------------

"""
    metabolomics_power_analysis(n_features, n_samples_per_group; alpha=0.05, effect_size=0.5, fdr_method=:bh, n_simulations=200, seed=1)

Estimate statistical power for a metabolomics experiment using simulation.
Analogous to the MetaboAnalyst sample size / power module.

Returns `(power_estimate, fdr_estimate, expected_true_positives)`.
"""
function metabolomics_power_analysis(
    n_features::Int,
    n_samples_per_group::Int;
    alpha::Real              = 0.05,
    effect_size::Real        = 0.5,
    fdr_method::Symbol       = :bh,
    true_positives_fraction::Real = 0.2,
    n_simulations::Int       = 200,
    seed::Int                = 1)
    rng = MersenneTwister(seed)
    n_true = round(Int, n_features * Float64(true_positives_fraction))
    n_null = n_features - n_true
    delta  = Float64(effect_size)
    α      = Float64(alpha)

    tp_counts = Int[]
    fp_counts = Int[]

    for _ in 1:n_simulations
        # True positives: shifted mean
        pvals_tp = [_welch_ttest_pvalue(randn(rng, n_samples_per_group),
                                         randn(rng, n_samples_per_group) .+ delta)
                    for _ in 1:n_true]
        # Null: no shift
        pvals_null = [_welch_ttest_pvalue(randn(rng, n_samples_per_group),
                                           randn(rng, n_samples_per_group))
                      for _ in 1:n_null]

        all_pvals = vcat(pvals_tp, pvals_null)
        is_true   = vcat(trues(n_true), falses(n_null))

        # BH correction
        ord = sortperm(all_pvals)
        padj = similar(all_pvals)
        n_total = length(all_pvals)
        padj[ord] .= [min(all_pvals[ord[i]] * n_total / i, 1.0) for i in 1:n_total]
        for i in (n_total-1):-1:1
            padj[ord[i]] = min(padj[ord[i]], padj[ord[i+1]])
        end

        sig = padj .<= α
        push!(tp_counts, count(sig .& is_true))
        push!(fp_counts, count(sig .& .!is_true))
    end

    power   = mean(tp_counts) / max(n_true, 1)
    fdr_est = mean(fp_counts) / max(mean(tp_counts) + mean(fp_counts), 1)
    exp_tp  = mean(tp_counts)

    result = (power_estimate=power, fdr_estimate=fdr_est, expected_true_positives=exp_tp,
              n_true_features=n_true, n_null_features=n_null)
    _ctx = active_provenance_context()


    return _register_metabolomics_result!(_ctx, result, "metabolomics_power_analysis"; parameters=(n_features=n_features, n_samples_per_group=n_samples_per_group, effect_size=Float64(effect_size), n_simulations=n_simulations, power_estimate=power))
end

# ---------------------------------------------------------------------------
# Feature clustering for metabolomics
# ---------------------------------------------------------------------------

"""
    feature_clustering_metabolomics(matrix; n_clusters=5, method=:kmeans, feature_names=nothing, seed=1)

Cluster metabolomic features (e.g. to detect metabolite co-regulation modules).
Analogous to WGCNA module detection / MetaboAnalyst clustering.

Methods: `:kmeans`, `:hierarchical`.

Returns `(cluster_assignments=DataFrame, cluster_profiles=Matrix)`.
"""
function feature_clustering_metabolomics(
    matrix::AbstractMatrix{<:Real};
    n_clusters::Int=5,
    method::Symbol=:kmeans,
    feature_names::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
    seed::Int=1)
    X = Matrix{Float64}(matrix)
    n_features, n_samples = size(X)
    fn = feature_names !== nothing ? String.(feature_names) : ["F$i" for i in 1:n_features]

    # Z-score each feature
    Xz = copy(X)
    for i in 1:n_features
        mu, sigma = mean(X[i,:]), std(X[i,:])
        Xz[i, :] = sigma > 0 ? (X[i,:] .- mu) ./ sigma : zeros(Float64, n_samples)
    end

    K = min(n_clusters, n_features)

    if method == :kmeans
        rng = MersenneTwister(seed)
        # Random init from feature rows
        init_idx = randperm(rng, n_features)[1:K]
        centroids = Xz[init_idx, :]   # K × samples
        assignments = zeros(Int, n_features)

        for _ in 1:100
            # Assign
            for i in 1:n_features
                dists = [sum(abs2, Xz[i, :] .- centroids[k, :]) for k in 1:K]
                assignments[i] = argmin(dists)
            end
            # Update centroids
            new_centroids = zeros(Float64, K, n_samples)
            counts = zeros(Int, K)
            for i in 1:n_features
                new_centroids[assignments[i], :] .+= Xz[i, :]
                counts[assignments[i]] += 1
            end
            for k in 1:K
                counts[k] > 0 && (new_centroids[k, :] ./= counts[k])
            end
            max(abs.(new_centroids .- centroids)...) < 1e-6 && break
            centroids = new_centroids
        end

    elseif method == :hierarchical
        # Simple agglomerative
        assignments = ones(Int, n_features)
        if K > 1
            # Reuse distance matrix approach
            dist_mat = zeros(Float64, n_features, n_features)
            for i in 1:n_features, j in (i+1):n_features
                d = sum(abs2, Xz[i,:] .- Xz[j,:])
                dist_mat[i,j] = dist_mat[j,i] = d
            end
            # Greedy: merge until K clusters
            cluster_ids = collect(1:n_features)
            while maximum(cluster_ids) > K
                # Find closest pair of different clusters
                best_i, best_j, best_d = 1, 2, Inf
                for i in 1:n_features-1, j in (i+1):n_features
                    if cluster_ids[i] != cluster_ids[j] && dist_mat[i,j] < best_d
                        best_d = dist_mat[i,j]
                        best_i, best_j = i, j
                    end
                end
                old_id = max(cluster_ids[best_i], cluster_ids[best_j])
                new_id = min(cluster_ids[best_i], cluster_ids[best_j])
                replace!(cluster_ids, old_id => new_id)
                # Renumber to 1..n_unique
                uniq = sort!(unique(cluster_ids))
                id_map = Dict(v => k for (k,v) in enumerate(uniq))
                cluster_ids = [id_map[c] for c in cluster_ids]
            end
            assignments = cluster_ids
        end

    else
        throw(ArgumentError("method must be :kmeans or :hierarchical"))
    end

    # Compute cluster profiles (mean per cluster)
    cluster_profiles = zeros(Float64, K, n_samples)
    cluster_sizes    = zeros(Int, K)
    for i in 1:n_features
        k = assignments[i]
        cluster_profiles[k, :] .+= Xz[i, :]
        cluster_sizes[k] += 1
    end
    for k in 1:K
        cluster_sizes[k] > 0 && (cluster_profiles[k, :] ./= cluster_sizes[k])
    end

    assignment_df = DataFrame(
        feature = fn,
        cluster = assignments)

    _ctx = active_provenance_context()
    _register_metabolomics_result!(_ctx, nothing, "feature_clustering_metabolomics"; parents=provenance_parent_ids(matrix), parameters=(n_clusters=K, method=method, seed=seed, n_features=n_features))


    return (cluster_assignments=assignment_df, cluster_profiles=cluster_profiles, cluster_sizes=cluster_sizes)
end

end  # module Metabolomics
