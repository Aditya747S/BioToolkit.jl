# ==============================================================================
# somatic.jl — Germline/somatic calling helpers and SV graph summaries
#
# References:
#   - Cibulskis et al. (2013) Nat Biotechnol 31:213-219 (Mutect concepts)
#   - Obenchain et al. (2014) Bioinformatics 30:2076-2078 (variant annotation)
#   - Alexandrov et al. (2020) Nature 578:94-101 (COSMIC signatures)
#   - Roth et al. (2014) Nat Methods 11:396-398 (PyClone CCF)
#   - Van Loo et al. (2010) PNAS 107:16910-16915 (ASCAT allelic imbalance)
#   - Olshen et al. (2004) Biostatistics 5:557-572 (CBS)
# ==============================================================================

module Somatic

using DataFrames
using Distributions
using LinearAlgebra
using Random
using Statistics
using ..BioToolkit: maybe_to_device, maybe_to_host, resolve_backend
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!

export germline_bayesian_call, somatic_tumor_normal_call, haplotypecaller_call, mutect2_call, strelka2_call
export sv_breakpoint_graph_table, vep_like_annotation
export tumor_mutational_burden, mutational_signature_nmf
export cbs_like_segments

# New exports
export cosmic_signature_attribution
export clonal_deconvolution_pyclone
export clonal_evolution_tree
export allelic_imbalance_test
export copy_number_from_coverage
export sv_fusion_candidates
export driver_enrichment
export variant_tier_classification
export somatic_hotspot_scan
export mutational_spectrum_96

_tool_cmd(program::String, args::Vector{String}) = `$(program) $(args...)`

@inline function _register_result_provenance!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

@inline function _expected_alt_fraction(
    mutated_copies::Real,
    total_copy_number::Real,
    purity::Real;
    normal_copy_number::Real=2.0,
    background_mutated_copies::Real=0.0)
    ρ = clamp(Float64(purity), 0.0, 1.0)
    tumor_cn = max(Float64(total_copy_number), eps(Float64))
    normal_cn = max(Float64(normal_copy_number), eps(Float64))
    denominator = ρ * tumor_cn + (1.0 - ρ) * normal_cn
    denominator > 0 || return eps(Float64)
    numerator = ρ * Float64(mutated_copies) + (1.0 - ρ) * Float64(background_mutated_copies)
    return clamp(numerator / denominator, eps(Float64), 1.0 - eps(Float64))
end

@inline function _binomial_loglikelihood(depth::Integer, alt::Integer, p::Real)
    d = Int(depth)
    a = Int(alt)
    0 <= a <= d || return -Inf
    q = clamp(Float64(p), eps(Float64), 1.0 - eps(Float64))
    return logpdf(Binomial(d, q), a)
end

function germline_bayesian_call(ref_count::Integer, alt_count::Integer; error_rate::Real=1e-3, het_prior::Real=1e-3)
    r = Int(ref_count)
    a = Int(alt_count)
    total = r + a
    total > 0 || return (genotype="./.", posterior_het=0.0)

    l00 = pdf(Binomial(total, Float64(error_rate)), a)
    l01 = pdf(Binomial(total, 0.5), a)
    l11 = pdf(Binomial(total, 1.0 - Float64(error_rate)), a)

    p00 = (1 - Float64(het_prior)) * 0.5
    p11 = (1 - Float64(het_prior)) * 0.5
    p01 = Float64(het_prior)

    z = l00 * p00 + l01 * p01 + l11 * p11
    post00 = (l00 * p00) / max(z, eps(Float64))
    post01 = (l01 * p01) / max(z, eps(Float64))
    post11 = (l11 * p11) / max(z, eps(Float64))

    gt = post01 >= max(post00, post11) ? "0/1" : post11 > post00 ? "1/1" : "0/0"
    result = (genotype=gt, posterior_het=post01, posterior_homref=post00, posterior_homalt=post11)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, result, "germline_bayesian_call"; parents=String[], parameters=(ref_count=r, alt_count=a, error_rate=Float64(error_rate), het_prior=Float64(het_prior), total=total))
end

function somatic_tumor_normal_call(
    tumor_ref::Integer,
    tumor_alt::Integer,
    normal_ref::Integer,
    normal_alt::Integer;
    min_tumor_vaf::Real=0.05,
    max_normal_vaf::Real=0.02,
    tumor_copy_number::Real=2.0,
    normal_copy_number::Real=2.0,
    tumor_purity::Real=1.0,
    somatic_multiplicity::Real=1.0,
    germline_multiplicity::Real=1.0,
    error_rate::Real=1e-3)
    tdepth = Int(tumor_ref + tumor_alt)
    ndepth = Int(normal_ref + normal_alt)
    tvaf = tdepth == 0 ? 0.0 : Float64(tumor_alt) / tdepth
    nvaf = ndepth == 0 ? 0.0 : Float64(normal_alt) / ndepth

    tumor_cn = max(Float64(tumor_copy_number), eps(Float64))
    normal_cn = max(Float64(normal_copy_number), eps(Float64))
    ρ = clamp(Float64(tumor_purity), 0.0, 1.0)
    somatic_mult = max(Float64(somatic_multiplicity), eps(Float64))
    germline_mult = max(Float64(germline_multiplicity), eps(Float64))

    tumor_expected_vaf_somatic = _expected_alt_fraction(
        somatic_mult,
        tumor_cn,
        ρ;
        normal_copy_number=normal_cn,
        background_mutated_copies=0.0)
    normal_expected_vaf_somatic = clamp(Float64(error_rate), eps(Float64), 1.0 - eps(Float64))

    tumor_expected_vaf_germline = _expected_alt_fraction(
        germline_mult,
        tumor_cn,
        ρ;
        normal_copy_number=normal_cn,
        background_mutated_copies=germline_mult)
    normal_expected_vaf_germline = clamp(germline_mult / normal_cn, eps(Float64), 1.0 - eps(Float64))

    loglik_somatic = _binomial_loglikelihood(tdepth, tumor_alt, tumor_expected_vaf_somatic) +
                     _binomial_loglikelihood(ndepth, normal_alt, normal_expected_vaf_somatic)
    loglik_germline = _binomial_loglikelihood(tdepth, tumor_alt, tumor_expected_vaf_germline) +
                      _binomial_loglikelihood(ndepth, normal_alt, normal_expected_vaf_germline)
    log10_lod = (loglik_somatic - loglik_germline) / log(10.0)
    delta = loglik_somatic - loglik_germline
    posterior_somatic = delta >= 0 ? 1.0 / (1.0 + exp(-delta)) : exp(delta) / (1.0 + exp(delta))
    somatic = posterior_somatic >= 0.5
    passes_thresholds = (tvaf >= Float64(min_tumor_vaf)) && (nvaf <= Float64(max_normal_vaf))

    result = (
        somatic=somatic,
        passes_thresholds=passes_thresholds,
        tumor_vaf=tvaf,
        normal_vaf=nvaf,
        tumor_expected_vaf_somatic=tumor_expected_vaf_somatic,
        normal_expected_vaf_somatic=normal_expected_vaf_somatic,
        tumor_expected_vaf_germline=tumor_expected_vaf_germline,
        normal_expected_vaf_germline=normal_expected_vaf_germline,
        posterior_somatic=posterior_somatic,
        log_lod=log10_lod)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, result, "somatic_tumor_normal_call"; parents=String[], parameters=(tumor_ref=Int(tumor_ref), tumor_alt=Int(tumor_alt), normal_ref=Int(normal_ref), normal_alt=Int(normal_alt), min_tumor_vaf=Float64(min_tumor_vaf), max_normal_vaf=Float64(max_normal_vaf), tumor_copy_number=tumor_cn, normal_copy_number=normal_cn, tumor_purity=ρ, somatic_multiplicity=somatic_mult, germline_multiplicity=germline_mult, error_rate=Float64(error_rate)))
end

function haplotypecaller_call(reference_fasta::AbstractString, bam_path::AbstractString; output_vcf::AbstractString="haplotypecaller.vcf.gz", intervals::Union{Nothing,AbstractString}=nothing, dry_run::Bool=true)
    _ctx = active_provenance_context()
    args = String["HaplotypeCaller", "-R", String(reference_fasta), "-I", String(bam_path), "-O", String(output_vcf)]
    intervals === nothing || append!(args, ["-L", String(intervals)])
    cmd = _tool_cmd("gatk", args)
    if dry_run
        result = (cmd=string(cmd), output=String(output_vcf), status=:planned)
        return _register_result_provenance!(_ctx, result, "haplotypecaller_call"; parents=String[], parameters=(reference=String(reference_fasta), input=String(bam_path), output=String(output_vcf), intervals=intervals === nothing ? "none" : String(intervals), dry_run=dry_run))
    end
    run(cmd)
    result = (cmd=string(cmd), output=String(output_vcf), status=:ok)

    return _register_result_provenance!(_ctx, result, "haplotypecaller_call"; parents=String[], parameters=(reference=String(reference_fasta), input=String(bam_path), output=String(output_vcf), intervals=intervals === nothing ? "none" : String(intervals), dry_run=dry_run))
end

function mutect2_call(reference_fasta::AbstractString, tumor_bam::AbstractString; normal_bam::Union{Nothing,AbstractString}=nothing, output_vcf::AbstractString="mutect2.vcf.gz", tumor_sample::AbstractString="TUMOR", normal_sample::AbstractString="NORMAL", dry_run::Bool=true)
    _ctx = active_provenance_context()

    args = String[
        "Mutect2", "-R", String(reference_fasta),
        "-I", String(tumor_bam), "-tumor", String(tumor_sample),
        "-O", String(output_vcf),
    ]
    if normal_bam !== nothing
        append!(args, ["-I", String(normal_bam), "-normal", String(normal_sample)])
    end
    cmd = _tool_cmd("gatk", args)
    if dry_run
        result = (cmd=string(cmd), output=String(output_vcf), status=:planned)
        return _register_result_provenance!(_ctx, result, "mutect2_call"; parents=String[], parameters=(reference=String(reference_fasta), tumor_bam=String(tumor_bam), normal_bam=normal_bam === nothing ? "none" : String(normal_bam), output=String(output_vcf), tumor_sample=String(tumor_sample), normal_sample=String(normal_sample), dry_run=dry_run))
    end
    run(cmd)
    result = (cmd=string(cmd), output=String(output_vcf), status=:ok)

    return _register_result_provenance!(_ctx, result, "mutect2_call"; parents=String[], parameters=(reference=String(reference_fasta), tumor_bam=String(tumor_bam), normal_bam=normal_bam === nothing ? "none" : String(normal_bam), output=String(output_vcf), tumor_sample=String(tumor_sample), normal_sample=String(normal_sample), dry_run=dry_run))
end

function strelka2_call(reference_fasta::AbstractString, tumor_bam::AbstractString, normal_bam::AbstractString; run_dir::AbstractString="strelka_run", output_vcf::AbstractString="somatic.snvs.vcf.gz", threads::Int=4, dry_run::Bool=true)
    _ctx = active_provenance_context()
    configure_cmd = _tool_cmd("configureStrelkaSomaticWorkflow.py", ["--referenceFasta", String(reference_fasta), "--tumorBam", String(tumor_bam), "--normalBam", String(normal_bam), "--runDir", String(run_dir)])
    run_cmd = _tool_cmd("python", [joinpath(String(run_dir), "runWorkflow.py"), "-m", "local", "-j", string(threads)])
    if dry_run
        result = (configure_cmd=string(configure_cmd), run_cmd=string(run_cmd), output=joinpath(String(run_dir), "results", "variants", String(output_vcf)), status=:planned)
        return _register_result_provenance!(_ctx, result, "strelka2_call"; parents=String[], parameters=(reference=String(reference_fasta), tumor_bam=String(tumor_bam), normal_bam=String(normal_bam), run_dir=String(run_dir), output=String(output_vcf), threads=threads, dry_run=dry_run))
    end
    run(configure_cmd)
    run(run_cmd)
    result = (configure_cmd=string(configure_cmd), run_cmd=string(run_cmd), output=joinpath(String(run_dir), "results", "variants", String(output_vcf)), status=:ok)

    return _register_result_provenance!(_ctx, result, "strelka2_call"; parents=String[], parameters=(reference=String(reference_fasta), tumor_bam=String(tumor_bam), normal_bam=String(normal_bam), run_dir=String(run_dir), output=String(output_vcf), threads=threads, dry_run=dry_run))
end

function sv_breakpoint_graph_table(sv_calls::DataFrame; chrom1_col::Symbol=:chrom1, pos1_col::Symbol=:pos1, chrom2_col::Symbol=:chrom2, pos2_col::Symbol=:pos2)
    hasproperty(sv_calls, chrom1_col) || throw(ArgumentError("missing chrom1 column"))
    hasproperty(sv_calls, pos1_col) || throw(ArgumentError("missing pos1 column"))
    hasproperty(sv_calls, chrom2_col) || throw(ArgumentError("missing chrom2 column"))
    hasproperty(sv_calls, pos2_col) || throw(ArgumentError("missing pos2 column"))

    node_a = String[]
    node_b = String[]
    weight = Int[]

    for row in eachrow(sv_calls)
        a = "$(row[chrom1_col]):$(row[pos1_col])"
        b = "$(row[chrom2_col]):$(row[pos2_col])"
        push!(node_a, a)
        push!(node_b, b)
        push!(weight, 1)
    end

    df = DataFrame(node_a=node_a, node_b=node_b, weight=weight)
    result = combine(groupby(df, [:node_a, :node_b]), :weight => sum => :support)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, result, "sv_breakpoint_graph_table"; parents=provenance_parent_ids(sv_calls), parameters=(chrom1_col=String(chrom1_col), pos1_col=String(pos1_col), chrom2_col=String(chrom2_col), pos2_col=String(pos2_col), edge_count=nrow(result)))
end

function vep_like_annotation(variants::DataFrame; consequence_col::Symbol=:consequence)
    hasproperty(variants, consequence_col) || throw(ArgumentError("variants table must include consequence column"))
    severity_rank = Dict(
        "stop_gained"         => 5,
        "frameshift_variant"  => 5,
        "missense_variant"    => 4,
        "splice_region_variant" => 3,
        "synonymous_variant"  => 2,
        "intron_variant"      => 1)
    out = copy(variants)
    out[!, :impact_score] = [get(severity_rank, String(v), 0) for v in out[!, consequence_col]]
    out[!, :impact_label] = ifelse.(out.impact_score .>= 5, "HIGH", ifelse.(out.impact_score .>= 3, "MODERATE", "LOW"))
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, out, "vep_like_annotation"; parents=provenance_parent_ids(variants), parameters=(consequence_col=String(consequence_col), row_count=nrow(out)))
end

"""
    tumor_mutational_burden(variants; callable_mb=30.0)

Compute per-sample TMB as non-synonymous variant count per callable megabase.
"""
function tumor_mutational_burden(variants::DataFrame; sample_col::Symbol=:sample, consequence_col::Symbol=:consequence, callable_mb::Real=30.0)
    hasproperty(variants, sample_col) || throw(ArgumentError("variants table missing sample column"))
    hasproperty(variants, consequence_col) || throw(ArgumentError("variants table missing consequence column"))
    Float64(callable_mb) > 0 || throw(ArgumentError("callable_mb must be positive"))

    nonsyn = Set(["missense_variant", "frameshift_variant", "stop_gained", "start_lost", "splice_region_variant"])
    df = copy(variants)
    df[!, :is_nonsynonymous] = [String(v) in nonsyn for v in df[!, consequence_col]]
    agg = combine(groupby(df, sample_col), :is_nonsynonymous => sum => :n_nonsynonymous)
    agg[!, :tmb_per_mb] = Float64.(agg.n_nonsynonymous) ./ Float64(callable_mb)
    sort!(agg, :tmb_per_mb, rev=true)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, agg, "tumor_mutational_burden"; parents=provenance_parent_ids(variants), parameters=(sample_col=String(sample_col), consequence_col=String(consequence_col), callable_mb=Float64(callable_mb), row_count=nrow(agg)))
end

"""
    mutational_signature_nmf(catalog; n_signatures=5, n_iter=300)

NMF decomposition of mutation catalogs into signatures and exposures.
Expected catalog shape is contexts x samples.
"""
function mutational_signature_nmf(catalog::AbstractMatrix{<:Real}; n_signatures::Int=5, n_iter::Int=300, seed::Int=1, backend::Symbol=:auto)
    X = max.(Float64.(catalog), 0.0)
    n_ctx, n_samples = size(X)
    k = clamp(n_signatures, 1, min(n_ctx, n_samples))
    selected = resolve_backend(; backend=backend)
    Xdev = maybe_to_device(X; backend=selected)

    rng = MersenneTwister(seed)
    W = maybe_to_device(rand(rng, n_ctx, k) .+ 1e-6; backend=selected)
    H = maybe_to_device(rand(rng, k, n_samples) .+ 1e-6; backend=selected)

    for _ in 1:max(1, n_iter)
        H .*= (W' * Xdev) ./ max.(W' * W * H, 1e-12)
        W .*= (Xdev * H') ./ max.(W * (H * H'), 1e-12)
    end

    signatures = maybe_to_host(W ./ max.(sum(W, dims=1), eps(Float64)))
    exposures = maybe_to_host(H)
    signatures_host = Matrix{Float64}(signatures)
    exposures_host = Matrix{Float64}(exposures)
    rel_err = norm(X - signatures_host * exposures_host) / max(norm(X), eps(Float64))
    result = (signatures=signatures_host, exposures=exposures_host, relative_error=rel_err, backend=selected)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, result, "mutational_signature_nmf"; parents=provenance_parent_ids(catalog), parameters=(n_signatures=k, n_iter=n_iter, seed=seed, backend=String(selected), relative_error=rel_err))
end

"""
    cbs_like_segments(log2_ratio; min_bins=10)

Simple binary segmentation of log2 copy-ratio profiles.
"""
function cbs_like_segments(log2_ratio::AbstractVector{<:Real}; min_bins::Int=10, threshold::Real=0.25)
    _ctx = active_provenance_context()
    x = Float64.(log2_ratio)
    n = length(x)
    if n == 0
        result = DataFrame(segment_start=Int[], segment_end=Int[], n_bins=Int[], mean_log2_ratio=Float64[])
        return _register_result_provenance!(_ctx, result, "cbs_like_segments"; parents=provenance_parent_ids(log2_ratio), parameters=(min_bins=min_bins, threshold=Float64(threshold), segment_count=0))
    end
    min_bins >= 2 || throw(ArgumentError("min_bins must be >= 2"))

    breaks = [1]
    i = 1
    while i <= n - min_bins
        remaining = n - i + 1
        if remaining < 2 * min_bins
            break
        end
        left = i:(i + min_bins - 1)
        right = (i + min_bins):(i + 2 * min_bins - 1)
        Δ = abs(mean(x[left]) - mean(x[right]))
        if Δ >= Float64(threshold)
            push!(breaks, i + min_bins)
            i += min_bins
        else
            i += 1
        end
    end
    push!(breaks, n + 1)

    out = DataFrame(segment_start=Int[], segment_end=Int[], n_bins=Int[], mean_log2_ratio=Float64[])
    for b in 1:(length(breaks) - 1)
        s = breaks[b]
        e = breaks[b + 1] - 1
        push!(out, (s, e, e - s + 1, mean(@view x[s:e])))
    end
    return _register_result_provenance!(_ctx, out, "cbs_like_segments"; parents=provenance_parent_ids(log2_ratio), parameters=(min_bins=min_bins, threshold=Float64(threshold), segment_count=nrow(out)))
end

# ---------------------------------------------------------------------------
# COSMIC Signature Attribution (cosine-similarity)
# ---------------------------------------------------------------------------

"""
    cosmic_signature_attribution(sample_catalog, cosmic_signatures; signature_names=nothing)

Assign COSMIC SBS signatures to tumour samples by cosine similarity, then
refine exposures via constrained least-squares (NNLS).

`sample_catalog`   : 96-context mutation counts vector (length 96) or matrix
                     (96 × n_samples).
`cosmic_signatures`: 96 × n_sigs reference matrix (columns sum to 1).

Returns a `DataFrame` with one row per sample and columns for every signature
exposure plus `reconstruction_cosine` quality metric.
"""
function cosmic_signature_attribution(
    sample_catalog::AbstractMatrix{<:Real},
    cosmic_signatures::AbstractMatrix{<:Real};
    signature_names::Union{Nothing,AbstractVector}=nothing,
    min_exposure::Real=0.01)
    C = Matrix{Float64}(sample_catalog)   # 96 × n_samples
    S = Matrix{Float64}(cosmic_signatures) # 96 × n_sigs
    size(C, 1) == size(S, 1) || throw(DimensionMismatch("catalog and signature matrix must have the same number of trinucleotide contexts (rows)"))

    n_sigs    = size(S, 2)
    n_samples = size(C, 2)

    sig_names = signature_names !== nothing ? String.(signature_names) :
                ["SBS$(i)" for i in 1:n_sigs]
    length(sig_names) == n_sigs || throw(DimensionMismatch("signature_names length must match number of signatures"))

    results = DataFrame()
    results[!, :sample_index] = 1:n_samples

    exposures_mat = zeros(Float64, n_samples, n_sigs)
    reconstruction_cosine = zeros(Float64, n_samples)

    for s in 1:n_samples
        y = max.(C[:, s], 0.0)
        total = sum(y)
        y_norm = total > 0 ? y ./ total : y

        # NNLS via projected gradient descent
        x = ones(Float64, n_sigs) ./ n_sigs
        lr = 1e-2
        for _ in 1:500
            grad = S' * (S * x .- y_norm)
            x = max.(x .- lr .* grad, 0.0)
            s_sum = sum(x)
            s_sum > 0 && (x ./= s_sum)
        end

        # Zero out negligible exposures
        x[x .< Float64(min_exposure)] .= 0.0
        s_sum = sum(x)
        s_sum > 0 && (x ./= s_sum)

        exposures_mat[s, :] = x
        recon = S * x
        cos_sim = dot(y_norm, recon) / max(norm(y_norm) * norm(recon), eps(Float64))
        reconstruction_cosine[s] = cos_sim
    end

    for (k, name) in enumerate(sig_names)
        results[!, Symbol(name)] = exposures_mat[:, k]
    end
    results[!, :reconstruction_cosine] = reconstruction_cosine
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, results, "cosmic_signature_attribution"; parents=provenance_parent_ids(sample_catalog, cosmic_signatures), parameters=(signature_count=n_sigs, sample_count=n_samples, min_exposure=Float64(min_exposure), reconstruction_cosine_mean=mean(reconstruction_cosine)))
end

# Convenience: single-sample vector overload
function cosmic_signature_attribution(
    sample_catalog::AbstractVector{<:Real},
    cosmic_signatures::AbstractMatrix{<:Real},
    kw...)
    return cosmic_signature_attribution(reshape(Float64.(sample_catalog), :, 1), cosmic_signatures; kw...)
end

# ---------------------------------------------------------------------------
# Clonal Deconvolution — PyClone-like EM for CCF estimation
# ---------------------------------------------------------------------------

"""
    clonal_deconvolution_pyclone(tumor_vafs, normal_vafs; n_clones=3, purity=nothing, n_iter=200)

EM-based cancer cell fraction (CCF) decomposition, inspired by PyClone.

Fits a mixture of `n_clones` Beta-distributed clonal CCF components to the
observed VAF distribution under a given tumour purity.

Returns `(ccf_means, mixing_weights, responsibilities, bic)`.
"""
function clonal_deconvolution_pyclone(
    tumor_vafs::AbstractVector{<:Real};
    normal_vafs::Union{Nothing,AbstractVector{<:Real}}=nothing,
    n_clones::Int=3,
    purity::Real=0.7,
    n_iter::Int=200,
    seed::Int=1,
    copy_numbers::Union{Nothing,AbstractVector{<:Real}}=nothing,
    mutation_multiplicity::Union{Nothing,AbstractVector{<:Real}}=nothing,
    normal_copy_number::Real=2.0)
    n_clones >= 1 || throw(ArgumentError("n_clones must be >= 1"))
    0 < purity <= 1 || throw(ArgumentError("purity must be in (0, 1]"))

    vafs = clamp.(Float64.(tumor_vafs), 1e-6, 1.0 - 1e-6)
    n = length(vafs)
    ρ = Float64(purity)
    tumor_copy_numbers = copy_numbers === nothing ? fill(2.0, n) : Float64.(copy_numbers)
    mutation_multiplicities = mutation_multiplicity === nothing ? fill(1.0, n) : Float64.(mutation_multiplicity)

    length(tumor_copy_numbers) == n || throw(DimensionMismatch("copy_numbers must match tumor_vafs length"))
    length(mutation_multiplicities) == n || throw(DimensionMismatch("mutation_multiplicity must match tumor_vafs length"))
    all(tumor_copy_numbers .> 0) || throw(ArgumentError("copy_numbers must be positive"))
    all(mutation_multiplicities .> 0) || throw(ArgumentError("mutation_multiplicity must be positive"))

    # Convert VAF → CCF using a copy-number-aware mixture model.
    ccf = similar(vafs)
    normal_cn = max(Float64(normal_copy_number), eps(Float64))
    for i in eachindex(vafs)
        denom = ρ * tumor_copy_numbers[i] + (1.0 - ρ) * normal_cn
        ccf[i] = clamp(vafs[i] * denom / (ρ * mutation_multiplicities[i]), 1e-4, 1.0 - 1e-4)
    end

    K = min(n_clones, n)

    # Initialise cluster means from the empirical CCF quantiles for stability.
    sorted_ccf = sort(ccf)
    μ = K == 1 ? [mean(sorted_ccf)] : [quantile(sorted_ccf, (k - 0.5) / K) for k in 1:K]
    μ = clamp.(μ .+ 1e-6 .* rand(MersenneTwister(seed), K), 1e-4, 1.0 - 1e-4)
    π_k = fill(1.0 / K, K)
    conc = 5.0   # Beta concentration (shared)

    r = zeros(Float64, n, K)  # responsibilities

    log_lik_prev = -Inf
    for iter in 1:n_iter
        # E-step
        for i in 1:n
            for k in 1:K
                α = max(μ[k] * conc, 1e-4)
                β = max((1.0 - μ[k]) * conc, 1e-4)
                r[i, k] = log(π_k[k]) + logpdf(Beta(α, β), ccf[i])
            end
            # Log-sum-exp normalise
            row_max = maximum(r[i, :])
            r[i, :] = exp.(r[i, :] .- row_max)
            r[i, :] ./= max(sum(r[i, :]), eps(Float64))
        end

        # M-step
        Nk = vec(sum(r, dims=1))
        for k in 1:K
            π_k[k] = Nk[k] / n
            μ[k] = dot(r[:, k], ccf) / max(Nk[k], eps(Float64))
            μ[k] = clamp(μ[k], 1e-4, 1.0 - 1e-4)
        end

        # Log-likelihood
        log_lik = sum(
            log(max(sum(π_k[k] * pdf(Beta(max(μ[k]*conc, 1e-4), max((1-μ[k])*conc, 1e-4)), ccf[i]) for k in 1:K), eps(Float64)))
            for i in 1:n
        )
        abs(log_lik - log_lik_prev) < 1e-6 && break
        log_lik_prev = log_lik
    end

    n_params = 2 * K - 1
    bic = -2 * log_lik_prev + n_params * log(max(n, 2))
    sort_order = sortperm(μ)

    result = (
        ccf_means         = μ[sort_order],
        mixing_weights    = π_k[sort_order],
        responsibilities  = r[:, sort_order],
        bic               = bic,
        n_clones          = K)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, result, "clonal_deconvolution_pyclone"; parents=provenance_parent_ids(tumor_vafs, normal_vafs === nothing ? Float64[] : normal_vafs), parameters=(n_clones=n_clones, purity=Float64(purity), n_iter=n_iter, seed=seed, normal_copy_number=Float64(normal_copy_number)))
end

# ---------------------------------------------------------------------------
# Clonal Evolution Tree
# ---------------------------------------------------------------------------

"""
    clonal_evolution_tree(ccf_table; sample_col=:sample, clone_col=:clone, ccf_col=:ccf)

Build a simple clonal phylogeny from CCF values across samples using a
parent-assignment heuristic (CITUP / ClonEvol style):
A clone C is the parent of clone D in a sample if CCF(C) > CCF(D) and
there is no other clone E with CCF(C) > CCF(E) > CCF(D).

Returns a `DataFrame` edge list with `parent_clone`, `child_clone`,
and per-sample CCF columns.
"""
function clonal_evolution_tree(
    ccf_table::DataFrame;
    sample_col::Symbol=:sample,
    clone_col::Symbol=:clone,
    ccf_col::Symbol=:ccf)
    hasproperty(ccf_table, sample_col) || throw(ArgumentError("missing sample column"))
    hasproperty(ccf_table, clone_col) || throw(ArgumentError("missing clone column"))
    hasproperty(ccf_table, ccf_col) || throw(ArgumentError("missing ccf column"))

    clones  = sort!(unique(String.(ccf_table[!, clone_col])))
    samples = sort!(unique(String.(ccf_table[!, sample_col])))

    # Build CCF matrix: rows=clones, cols=samples
    ccf_mat = zeros(Float64, length(clones), length(samples))
    clone_idx  = Dict(c => i for (i, c) in enumerate(clones))
    sample_idx = Dict(s => j for (j, s) in enumerate(samples))

    for row in eachrow(ccf_table)
        c = String(row[clone_col])
        s = String(row[sample_col])
        haskey(clone_idx, c) && haskey(sample_idx, s) &&
            (ccf_mat[clone_idx[c], sample_idx[s]] = Float64(row[ccf_col]))
    end

    # Mean CCF across samples for ordering
    mean_ccf = vec(mean(ccf_mat, dims=2))
    order = sortperm(mean_ccf, rev=true)

    parent = Dict{String,String}()
    for i in 2:length(order)
        child_idx  = order[i]
        child_ccf  = mean_ccf[child_idx]
        best_parent = order[1]   # root as fallback
        best_diff  = Inf
        for j in 1:(i-1)
            p_idx = order[j]
            diff  = mean_ccf[p_idx] - child_ccf
            if diff > 0 && diff < best_diff
                best_diff   = diff
                best_parent = p_idx
            end
        end
        parent[clones[child_idx]] = clones[best_parent]
    end

    edges = DataFrame(parent_clone=String[], child_clone=String[])
    for (child, par) in parent
        push!(edges, (par, child))
    end

    # Attach CCF columns
    for (j, s) in enumerate(samples)
        edges[!, Symbol("ccf_$(s)_parent")] = [ccf_mat[clone_idx[r.parent_clone], j] for r in eachrow(edges)]
        edges[!, Symbol("ccf_$(s)_child")]  = [ccf_mat[clone_idx[r.child_clone], j] for r in eachrow(edges)]
    end

    result = (edges=edges, clones=clones, ccf_matrix=ccf_mat, samples=samples)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, result, "clonal_evolution_tree"; parents=provenance_parent_ids(ccf_table), parameters=(sample_col=String(sample_col), clone_col=String(clone_col), ccf_col=String(ccf_col), clone_count=length(clones), sample_count=length(samples)))
end

# ---------------------------------------------------------------------------
# Allelic Imbalance Test (ASCAT-style BAF)
# ---------------------------------------------------------------------------

"""
    allelic_imbalance_test(ref_counts, alt_counts; segments=nothing, min_depth=10)

Test genomic positions for allelic imbalance (loss of heterozygosity / CN gain)
using a Binomial test against the null B-allele frequency of 0.5.

`ref_counts` / `alt_counts`: Integer vectors of per-locus read counts.
`segments`                 : Optional `DataFrame` with `segment_start` and
                             `segment_end` columns (1-indexed) used to
                             compute segment-level BAF summaries.

Returns per-locus `DataFrame` with `ref`, `alt`, `depth`, `baf`, `pvalue`,
`imbalanced` and (if segments given) a segment summary.
"""
function allelic_imbalance_test(
    ref_counts::AbstractVector{<:Integer},
    alt_counts::AbstractVector{<:Integer};
    segments::Union{Nothing,DataFrame}=nothing,
    min_depth::Int=10,
    fdr_threshold::Real=0.05)
    n = length(ref_counts)
    n == length(alt_counts) || throw(DimensionMismatch("ref_counts and alt_counts must have equal length"))

    depth = Int.(ref_counts) .+ Int.(alt_counts)
    baf   = [depth[i] > 0 ? Float64(alt_counts[i]) / depth[i] : 0.5 for i in 1:n]
    pvals = [
        depth[i] >= min_depth ?
            2.0 * min(cdf(Binomial(depth[i], 0.5), alt_counts[i]),
                      ccdf(Binomial(depth[i], 0.5), alt_counts[i] - 1)) : 1.0
        for i in 1:n
    ]

    # BH FDR correction
    m = n
    order = sortperm(pvals)
    padj = similar(pvals)
    running = 1.0
    for rank in m:-1:1
        idx = order[rank]
        running = min(running, pvals[idx] * m / rank)
        padj[idx] = running
    end
    padj = clamp.(padj, 0.0, 1.0)

    locus_df = DataFrame(
        locus      = 1:n,
        ref        = Int.(ref_counts),
        alt        = Int.(alt_counts),
        depth      = depth,
        baf        = baf,
        pvalue     = pvals,
        padj       = padj,
        imbalanced = padj .< Float64(fdr_threshold))

    _ctx = active_provenance_context()
    if segments === nothing
        return _register_result_provenance!(_ctx, locus_df, "allelic_imbalance_test"; parents=provenance_parent_ids(ref_counts, alt_counts), parameters=(min_depth=min_depth, fdr_threshold=Float64(fdr_threshold), locus_count=nrow(locus_df)))
    end

    # Segment-level summary
    hasproperty(segments, :segment_start) || throw(ArgumentError("segments must have segment_start column"))
    hasproperty(segments, :segment_end)   || throw(ArgumentError("segments must have segment_end column"))

    seg_df = copy(segments)
    seg_baf        = Float64[]
    seg_n_imbal    = Int[]
    seg_mean_depth = Float64[]
    seg_lengths    = Float64[]

    for row in eachrow(segments)
        s = max(Int(row.segment_start), 1)
        e = min(Int(row.segment_end), n)
        if s > e
            push!(seg_baf, 0.5); push!(seg_n_imbal, 0); push!(seg_mean_depth, 0.0); push!(seg_lengths, 1.0)
            continue
        end
        idx = s:e
        push!(seg_baf, mean(baf[idx]))
        push!(seg_n_imbal, count(locus_df.imbalanced[idx]))
        push!(seg_mean_depth, mean(depth[idx]))
        push!(seg_lengths, Float64(e - s + 1))
    end

    seg_df[!, :mean_baf]        = seg_baf
    seg_df[!, :n_imbalanced]    = seg_n_imbal
    seg_df[!, :mean_depth]      = seg_mean_depth
    seg_df[!, :loh_fraction]    = seg_n_imbal ./ seg_lengths

    result = (loci=locus_df, segments=seg_df)

    return _register_result_provenance!(_ctx, result, "allelic_imbalance_test"; parents=provenance_parent_ids(ref_counts, alt_counts, segments), parameters=(min_depth=min_depth, fdr_threshold=Float64(fdr_threshold), locus_count=nrow(locus_df), segment_count=nrow(seg_df)))
end

# ---------------------------------------------------------------------------
# Copy Number from Coverage
# ---------------------------------------------------------------------------

"""
    copy_number_from_coverage(coverage; window=50, reference_coverage=nothing, ploidy=2)

Infer relative copy-number from a per-base or per-bin coverage vector.
Applies rolling-window median smoothing then converts to log2 copy-ratio
(as in CNVkit / Control-FREEC).

`reference_coverage`: If provided, normalise tumour by matched normal.

Returns a `DataFrame` with `bin`, `raw_ratio`, `log2_ratio`, `cn_call`.
"""
function copy_number_from_coverage(
    coverage::AbstractVector{<:Real};
    window::Int=50,
    reference_coverage::Union{Nothing,AbstractVector{<:Real}}=nothing,
    ploidy::Real=2.0,
    threshold_gain::Real=0.3,
    threshold_loss::Real=-0.3)
    window >= 1 || throw(ArgumentError("window must be >= 1"))
    x = Float64.(coverage)
    n = length(x)

    ref = reference_coverage !== nothing ? Float64.(reference_coverage) : ones(Float64, n)
    length(ref) == n || throw(DimensionMismatch("reference_coverage length must match coverage"))

    half = div(window, 2)
    smooth = zeros(Float64, n)
    for i in 1:n
        lo = max(1, i - half)
        hi = min(n, i + half)
        smooth[i] = median(x[lo:hi])
    end

    ref_mean = mean(ref)
    ref_mean = ref_mean > 0 ? ref_mean : 1.0
    tumor_mean = mean(smooth)
    tumor_mean = tumor_mean > 0 ? tumor_mean : 1.0

    # Normalise
    ratio = (smooth ./ tumor_mean) ./ (ref ./ ref_mean)
    ratio = max.(ratio, 1e-6)
    log2_ratio = log2.(ratio)

    cn_call = [
        lr >= Float64(threshold_gain) ? "gain" :
        lr <= Float64(threshold_loss) ? "loss" : "neutral"
        for lr in log2_ratio
    ]

    result = DataFrame(
        bin        = 1:n,
        raw_ratio  = ratio,
        log2_ratio = log2_ratio,
        cn_call    = cn_call)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, result, "copy_number_from_coverage"; parents=provenance_parent_ids(coverage, reference_coverage === nothing ? Float64[] : reference_coverage), parameters=(window=window, ploidy=Float64(ploidy), threshold_gain=Float64(threshold_gain), threshold_loss=Float64(threshold_loss), bin_count=nrow(result)))
end

# ---------------------------------------------------------------------------
# SV Fusion Candidates
# ---------------------------------------------------------------------------

"""
    sv_fusion_candidates(sv_calls, gene_annotations; max_distance=100_000)

Identify putative gene fusions from structural variant calls by intersecting
breakpoints with a gene annotation `DataFrame`.

`sv_calls`       : `DataFrame` with at minimum `chrom1`, `pos1`, `chrom2`, `pos2`.
`gene_annotations`: `DataFrame` with `chrom`, `gene_start`, `gene_end`, `gene_name`.

Returns a `DataFrame` of candidate fusions with `gene5p`, `gene3p`,
`sv_type` (if present), `distance_5p`, `distance_3p`.
"""
function sv_fusion_candidates(
    sv_calls::DataFrame,
    gene_annotations::DataFrame;
    max_distance::Int=100_000)
    required_sv = [:chrom1, :pos1, :chrom2, :pos2]
    for col in required_sv
        hasproperty(sv_calls, col) || throw(ArgumentError("sv_calls missing column $col"))
    end
    required_ann = [:chrom, :gene_start, :gene_end, :gene_name]
    for col in required_ann
        hasproperty(gene_annotations, col) || throw(ArgumentError("gene_annotations missing column $col"))
    end

    function nearest_gene(chrom, pos)
        same = filter(r -> String(r.chrom) == String(chrom), gene_annotations)
        isempty(same) && return ("intergenic", typemax(Int))
        best_dist = typemax(Int)
        best_gene = "intergenic"
        for row in eachrow(same)
            gs = Int(row.gene_start)
            ge = Int(row.gene_end)
            dist = max(0, max(gs - Int(pos), Int(pos) - ge))
            if dist < best_dist
                best_dist = dist
                best_gene = String(row.gene_name)
            end
        end
        return (best_gene, best_dist)
    end

    gene5p = String[]
    gene3p = String[]
    dist5p = Int[]
    dist3p = Int[]
    sv_type_out = String[]

    for row in eachrow(sv_calls)
        g1, d1 = nearest_gene(row.chrom1, Int(row.pos1))
        g2, d2 = nearest_gene(row.chrom2, Int(row.pos2))
        (d1 > max_distance && d2 > max_distance) && continue
        push!(gene5p, g1)
        push!(gene3p, g2)
        push!(dist5p, d1)
        push!(dist3p, d2)
        t = hasproperty(sv_calls, :sv_type) ? String(row.sv_type) : "unknown"
        push!(sv_type_out, t)
    end

    result = DataFrame(
        gene5p      = gene5p,
        gene3p      = gene3p,
        sv_type     = sv_type_out,
        distance_5p = dist5p,
        distance_3p = dist3p)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, result, "sv_fusion_candidates"; parents=provenance_parent_ids(sv_calls, gene_annotations), parameters=(max_distance=max_distance, fusion_count=nrow(result)))
end

# ---------------------------------------------------------------------------
# Driver Gene Enrichment
# ---------------------------------------------------------------------------

"""
    driver_enrichment(variants, driver_genes; gene_col=:gene, consequence_col=:consequence)

Test for enrichment of non-synonymous variants in known cancer driver genes
using a Fisher's exact test (one-sided), analogous to dNdS / OncodriveFML
output interpretation.

Returns a `DataFrame` sorted by `pvalue` with columns:
`gene`, `n_driver_nonsynonymous`, `n_background_nonsynonymous`,
`odds_ratio`, `pvalue`, `padj`, `is_driver`.
"""
function driver_enrichment(
    variants::DataFrame,
    driver_genes::AbstractVector{<:AbstractString};
    gene_col::Symbol=:gene,
    consequence_col::Symbol=:consequence,
    fdr_threshold::Real=0.05)
    hasproperty(variants, gene_col) || throw(ArgumentError("variants missing gene column"))
    hasproperty(variants, consequence_col) || throw(ArgumentError("variants missing consequence column"))

    nonsyn = Set(["missense_variant", "frameshift_variant", "stop_gained", "start_lost", "splice_region_variant", "stop_lost"])
    drivers = Set(String.(driver_genes))

    df = copy(variants)
    df[!, :is_driver]   = [String(g) in drivers for g in df[!, gene_col]]
    df[!, :is_nonsyn]   = [String(c) in nonsyn  for c in df[!, consequence_col]]

    total_driver_nonsyn = count(r -> r.is_driver && r.is_nonsyn, eachrow(df))
    total_bg_nonsyn     = count(r -> !r.is_driver && r.is_nonsyn, eachrow(df))

    genes = sort!(unique(String.(df[!, gene_col])))
    out = DataFrame(gene=String[], n_nonsyn=Int[], is_driver_gene=Bool[], odds_ratio=Float64[], pvalue=Float64[])

    for g in genes
        sub = filter(r -> String(r[gene_col]) == g, df)
        a = count(r -> r.is_nonsyn, eachrow(sub))
        b = total_driver_nonsyn - (String(sub[!, gene_col][1]) in drivers ? a : 0)
        c = total_bg_nonsyn - (String(sub[!, gene_col][1]) in drivers ? 0 : a)
        d_val = nrow(df) - a - b - c
        or = ((a + 0.5) * (max(d_val, 0) + 0.5)) / ((max(b, 0) + 0.5) * (max(c, 0) + 0.5))
        # Fisher right tail
        N = a + max(b, 0) + max(c, 0) + max(d_val, 0)
        K = a + max(c, 0)
        nn = a + max(b, 0)
        pval = (N > 0 && K > 0 && nn > 0) ? ccdf(Hypergeometric(N, K, nn), a - 1) : 1.0
        push!(out, (g, a, g in drivers, or, pval))
    end

    sort!(out, :pvalue)
    m = nrow(out)
    padj = similar(out.pvalue)
    running = 1.0
    for rank in m:-1:1
        running = min(running, out.pvalue[rank] * m / rank)
        padj[rank] = running
    end
    out[!, :padj] = clamp.(padj, 0.0, 1.0)
    out[!, :significant] = out.padj .< Float64(fdr_threshold)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, out, "driver_enrichment"; parents=provenance_parent_ids(variants), parameters=(gene_col=String(gene_col), consequence_col=String(consequence_col), driver_count=length(driver_genes), row_count=nrow(out), fdr_threshold=Float64(fdr_threshold)))
end

# ---------------------------------------------------------------------------
# Variant Tier Classification (ACMG-inspired)
# ---------------------------------------------------------------------------

"""
    variant_tier_classification(variants; consequence_col=:consequence, vaf_col=nothing, depth_col=nothing)

Multi-tier variant classification analogous to CGI / ClinVar tier labelling.

Tiers:
- Tier I: High-impact coding variants in known drivers (stop_gained, frameshift)
- Tier II: Potentially damaging (missense, splice)
- Tier III: Low-impact or non-coding
- Tier IV: Synonymous / intergenic

Returns the input `DataFrame` with added columns `tier`, `tier_label`.
"""
function variant_tier_classification(
    variants::DataFrame;
    consequence_col::Symbol=:consequence,
    vaf_col::Union{Nothing,Symbol}=nothing,
    depth_col::Union{Nothing,Symbol}=nothing,
    min_vaf::Real=0.0,
    min_depth::Int=0)
    hasproperty(variants, consequence_col) || throw(ArgumentError("missing consequence column"))

    tier_map = Dict(
        "stop_gained"           => 1,
        "frameshift_variant"    => 1,
        "start_lost"            => 1,
        "stop_lost"             => 1,
        "missense_variant"      => 2,
        "splice_region_variant" => 2,
        "splice_acceptor_variant" => 2,
        "splice_donor_variant"  => 2,
        "inframe_insertion"     => 3,
        "inframe_deletion"      => 3,
        "synonymous_variant"    => 4,
        "intron_variant"        => 4,
        "intergenic_variant"    => 4)
    tier_label = Dict(1 => "Tier_I", 2 => "Tier_II", 3 => "Tier_III", 4 => "Tier_IV")

    out = copy(variants)
    out[!, :tier] = [get(tier_map, String(c), 3) for c in out[!, consequence_col]]
    out[!, :tier_label] = [tier_label[t] for t in out.tier]

    # Filter on VAF/depth if requested
    if vaf_col !== nothing && hasproperty(out, vaf_col)
        out = filter(r -> Float64(r[vaf_col]) >= Float64(min_vaf), out)
    end
    if depth_col !== nothing && hasproperty(out, depth_col)
        out = filter(r -> Int(r[depth_col]) >= min_depth, out)
    end

    sort!(out, :tier)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, out, "variant_tier_classification"; parents=provenance_parent_ids(variants), parameters=(consequence_col=String(consequence_col), min_vaf=Float64(min_vaf), min_depth=min_depth, row_count=nrow(out)))
end

# ---------------------------------------------------------------------------
# Somatic Hotspot Scan
# ---------------------------------------------------------------------------

"""
    somatic_hotspot_scan(variants; gene_col=:gene, position_col=:position, min_recurrence=3)

Identify recurrently mutated hotspot positions across samples, analogous to
MutSig2CV hotspot output.

Returns a `DataFrame` with `gene`, `position`, `n_samples`, `recurrence_rank`.
"""
function somatic_hotspot_scan(
    variants::DataFrame;
    gene_col::Symbol=:gene,
    position_col::Symbol=:position,
    sample_col::Union{Nothing,Symbol}=:sample,
    min_recurrence::Int=3)
    hasproperty(variants, gene_col) || throw(ArgumentError("missing gene column"))
    hasproperty(variants, position_col) || throw(ArgumentError("missing position column"))

    group_cols = [gene_col, position_col]
    if sample_col !== nothing && hasproperty(variants, sample_col)
        # Count unique samples per position
        df = combine(
            groupby(variants, group_cols),
            sample_col => (x -> length(unique(String.(x)))) => :n_samples
        )
    else
        df = combine(groupby(variants, group_cols), nrow => :n_samples)
    end

    filter!(r -> r.n_samples >= min_recurrence, df)
    sort!(df, :n_samples, rev=true)
    df[!, :recurrence_rank] = 1:nrow(df)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, df, "somatic_hotspot_scan"; parents=provenance_parent_ids(variants), parameters=(gene_col=String(gene_col), position_col=String(position_col), sample_col=sample_col === nothing ? "none" : String(sample_col), min_recurrence=min_recurrence, hotspot_count=nrow(df)))
end

# ---------------------------------------------------------------------------
# Mutational Spectrum (96 trinucleotide contexts summary)
# ---------------------------------------------------------------------------

"""
    mutational_spectrum_96(variants; context_col=:trinucleotide_context, ref_col=:ref, alt_col=:alt)

Summarise SNV mutations into SBS-96 trinucleotide context counts per sample.
The `context_col` should contain 3-character trinucleotide context strings (e.g. \"ACA\").

Returns a `DataFrame` with `context`, `substitution`, `count`, `fraction`.
"""
function mutational_spectrum_96(
    variants::DataFrame;
    context_col::Symbol=:trinucleotide_context,
    ref_col::Symbol=:ref,
    alt_col::Symbol=:alt,
    sample_col::Union{Nothing,Symbol}=nothing)
    hasproperty(variants, context_col) || throw(ArgumentError("missing trinucleotide_context column"))
    hasproperty(variants, ref_col) || throw(ArgumentError("missing ref column"))
    hasproperty(variants, alt_col) || throw(ArgumentError("missing alt column"))

    df = copy(variants)
    df[!, :mutation_type] = [
        "$(uppercase(String(r[ref_col])))>$(uppercase(String(r[alt_col])))"
        for r in eachrow(df)
    ]
    df[!, :sbs_key] = [
        "$(uppercase(string(r[context_col][1])))[$(r.mutation_type)]$(uppercase(string(r[context_col][3])))"
        for r in eachrow(df) if ncodeunits(String(r[context_col])) >= 3
    ]

    group_cols = sample_col !== nothing && hasproperty(df, sample_col) ? [sample_col, :sbs_key] : [:sbs_key]
    counts = combine(groupby(df, group_cols), nrow => :count)

    if sample_col === nothing || !hasproperty(df, sample_col)
        total = sum(counts.count)
        counts[!, :fraction] = counts.count ./ max(total, 1)
    else
        counts = transform(
            groupby(counts, sample_col),
            :count => (c -> c ./ max(sum(c), 1)) => :fraction
        )
    end
    sort!(counts, :count, rev=true)
    _ctx = active_provenance_context()

    return _register_result_provenance!(_ctx, counts, "mutational_spectrum_96"; parents=provenance_parent_ids(variants), parameters=(context_col=String(context_col), ref_col=String(ref_col), alt_col=String(alt_col), sample_col=sample_col === nothing ? "none" : String(sample_col), row_count=nrow(counts)))
end

end
