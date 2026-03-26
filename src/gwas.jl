module GWAS

using Mmap
using DataFrames
using Statistics
using LinearAlgebra
using Random
using Distributions
using Base.Threads
using PooledArrays

using ..DifferentialExpression: benjamini_hochberg
using ..GenomicRanges: GenomicInterval, build_collection, find_overlaps
using ..Epigenetics: PeakSet, Peak
using ..SystemsBio: GeneNetwork
using ..Enrichment: EnrichmentDatabase, EnrichmentTerm, build_annotation_database, enrichment_test

export GenotypeMatrix, GWASResult, MetaAnalysisResult
export read_plink, write_plink, gwas_linear_scan, gwas_lmm_scan
export calculate_prs, prs_ldpred, prs_cross_validation, ld_clumping, meta_analyze
export overlap_gwas_peaks, gene_based_test

"""
    GenotypeMatrix

Decoded PLINK genotype matrix together with BIM/FAM metadata and source prefix.
"""
struct GenotypeMatrix <: AbstractMatrix{Float64}
    raw::Vector{UInt8}
    decoded::Matrix{Float64}
    bim::DataFrame
    fam::DataFrame
    prefix::String
end

GenotypeMatrix(decoded::AbstractMatrix{<:Real}, bim::DataFrame, fam::DataFrame; prefix::AbstractString="") =
    GenotypeMatrix(UInt8[], Matrix{Float64}(decoded), bim, fam, String(prefix))

"""
    GWASResult

Association scan result table for a single GWAS analysis.
"""
struct GWASResult
    snp_ids::PooledVector{String,UInt32,Vector{UInt32}}
    chromosomes::PooledVector{String,UInt32,Vector{UInt32}}
    positions::Vector{Int}
    alleles::Vector{Tuple{String,String}}
    gene_ids::PooledVector{String,UInt32,Vector{UInt32}}
    beta::Vector{Float64}
    standard_error::Vector{Float64}
    zscore::Vector{Float64}
    pvalue::Vector{Float64}
    sample_size::Int
    covariate_names::PooledVector{String,UInt32,Vector{UInt32}}
    phenotype_name::String
    method::String
end

"""
    MetaAnalysisResult

Meta-analysis summary across multiple GWAS result tables.
"""
struct MetaAnalysisResult
    snp_ids::PooledVector{String,UInt32,Vector{UInt32}}
    chromosomes::PooledVector{String,UInt32,Vector{UInt32}}
    positions::Vector{Int}
    alleles::Vector{Tuple{String,String}}
    gene_ids::PooledVector{String,UInt32,Vector{UInt32}}
    beta::Vector{Float64}
    standard_error::Vector{Float64}
    zscore::Vector{Float64}
    pvalue::Vector{Float64}
    qvalue::Vector{Float64}
    tau2::Vector{Float64}
    i2::Vector{Float64}
    study_count::Int
    method::String
end

_pooled_strings(values) = PooledArray(String.(values))

function GWASResult(
    snp_ids::AbstractVector{<:AbstractString},
    chromosomes::AbstractVector{<:AbstractString},
    positions::AbstractVector{<:Integer},
    alleles::AbstractVector{<:Tuple},
    gene_ids::AbstractVector{<:AbstractString},
    beta::AbstractVector{<:Real},
    standard_error::AbstractVector{<:Real},
    zscore::AbstractVector{<:Real},
    pvalue::AbstractVector{<:Real},
    sample_size::Integer,
    covariate_names::AbstractVector{<:AbstractString},
    phenotype_name::AbstractString,
    method::AbstractString,
)
    return GWASResult(
        _pooled_strings(snp_ids),
        _pooled_strings(chromosomes),
        Int.(positions),
        [(String(left), String(right)) for (left, right) in alleles],
        _pooled_strings(gene_ids),
        Float64.(beta),
        Float64.(standard_error),
        Float64.(zscore),
        Float64.(pvalue),
        Int(sample_size),
        _pooled_strings(covariate_names),
        String(phenotype_name),
        String(method),
    )
end

function MetaAnalysisResult(
    snp_ids::AbstractVector{<:AbstractString},
    chromosomes::AbstractVector{<:AbstractString},
    positions::AbstractVector{<:Integer},
    alleles::AbstractVector{<:Tuple},
    gene_ids::AbstractVector{<:AbstractString},
    beta::AbstractVector{<:Real},
    standard_error::AbstractVector{<:Real},
    zscore::AbstractVector{<:Real},
    pvalue::AbstractVector{<:Real},
    qvalue::AbstractVector{<:Real},
    tau2::AbstractVector{<:Real},
    i2::AbstractVector{<:Real},
    study_count::Integer,
    method::AbstractString,
)
    return MetaAnalysisResult(
        _pooled_strings(snp_ids),
        _pooled_strings(chromosomes),
        Int.(positions),
        [(String(left), String(right)) for (left, right) in alleles],
        _pooled_strings(gene_ids),
        Float64.(beta),
        Float64.(standard_error),
        Float64.(zscore),
        Float64.(pvalue),
        Float64.(qvalue),
        Float64.(tau2),
        Float64.(i2),
        Int(study_count),
        String(method),
    )
end

Base.size(genotypes::GenotypeMatrix) = size(genotypes.decoded)
Base.getindex(genotypes::GenotypeMatrix, row::Int, column::Int) = genotypes.decoded[row, column]
Base.IndexStyle(::Type{GenotypeMatrix}) = IndexCartesian()
Base.eltype(::Type{GenotypeMatrix}) = Float64
Base.Matrix(genotypes::GenotypeMatrix) = copy(genotypes.decoded)

_matrix(genotypes::GenotypeMatrix) = genotypes.decoded
_matrix(genotypes::AbstractMatrix{<:Real}) = Matrix{Float64}(genotypes)

function _table_rows(path::AbstractString)
    rows = Vector{Vector{String}}()
    open(path, "r") do io
        for line in eachline(io)
            stripped = strip(line)
            isempty(stripped) && continue
            push!(rows, split(stripped))
        end
    end
    return rows
end

function _read_bim(path::AbstractString)
    rows = _table_rows(path)
    return DataFrame(
        chromosome = [row[1] for row in rows],
        snp_id = [row[2] for row in rows],
        genetic_distance = parse.(Float64, [row[3] for row in rows]),
        position = parse.(Int, [row[4] for row in rows]),
        allele1 = [row[5] for row in rows],
        allele2 = [row[6] for row in rows],
    )
end

function _read_fam(path::AbstractString)
    rows = _table_rows(path)
    return DataFrame(
        family_id = [row[1] for row in rows],
        sample_id = [row[2] for row in rows],
        paternal_id = [row[3] for row in rows],
        maternal_id = [row[4] for row in rows],
        sex = parse.(Int, [row[5] for row in rows]),
        phenotype = parse.(Float64, [row[6] for row in rows]),
    )
end

_plink_code(value) = ismissing(value) ? 0x01 : UInt8(clamp(Int(round(Float64(value))), 0, 2)) == 0x00 ? 0x00 : UInt8(clamp(Int(round(Float64(value))), 0, 2)) == 0x01 ? 0x02 : 0x03

function _decode_plink_code(byte::UInt8)
    code = byte & 0x03
    code == 0x00 && return 0.0
    code == 0x01 && return NaN
    code == 0x02 && return 1.0
    return 2.0
end

function _decode_bed(raw::Vector{UInt8}, nsamples::Int, nsnps::Int)
    length(raw) >= 3 || throw(ArgumentError("invalid .bed file"))
    raw[1:3] == UInt8[0x6c, 0x1b, 0x01] || throw(ArgumentError("invalid PLINK .bed header"))
    bytes = raw[4:end]
    bytes_per_snp = cld(nsamples, 4)
    length(bytes) >= nsnps * bytes_per_snp || throw(ArgumentError("bed file too short for declared dimensions"))
    decoded = fill(NaN, nsamples, nsnps)
    for snp in 1:nsnps
        base = (snp - 1) * bytes_per_snp
        for sample in 1:nsamples
            byte_index = base + cld(sample, 4)
            shift = 2 * ((sample - 1) % 4)
            decoded[sample, snp] = _decode_plink_code((bytes[byte_index] >> shift) & 0x03)
        end
    end
    return decoded
end

function _encode_bed(matrix::AbstractMatrix{<:Real})
    nsamples, nsnps = size(matrix)
    bytes_per_snp = cld(nsamples, 4)
    payload = fill(UInt8(0), nsnps * bytes_per_snp)
    for snp in 1:nsnps
        base = (snp - 1) * bytes_per_snp
        for sample in 1:nsamples
            byte_index = base + cld(sample, 4)
            shift = 2 * ((sample - 1) % 4)
            payload[byte_index] |= UInt8(_plink_code(matrix[sample, snp])) << shift
        end
    end
    return vcat(UInt8[0x6c, 0x1b, 0x01], payload)
end

"""
    read_plink(prefix::AbstractString)

Read a PLINK dataset given its file prefix.
"""
function read_plink(prefix::AbstractString)
    bim = _read_bim(String(prefix) * ".bim")
    fam = _read_fam(String(prefix) * ".fam")
    raw = open(String(prefix) * ".bed", "r") do io
        Mmap.mmap(io)
    end
    decoded = _decode_bed(raw, nrow(fam), nrow(bim))
    return GenotypeMatrix(raw, decoded, bim, fam, String(prefix))
end

"""
    write_plink(prefix, genotypes, bim, fam)

Write a PLINK dataset to disk.
"""
function write_plink(prefix::AbstractString, genotypes::AbstractMatrix{<:Real}, bim::DataFrame, fam::DataFrame)
    open(String(prefix) * ".bim", "w") do io
        for row in eachrow(bim)
            println(io, join((row.chromosome, row.snp_id, row.genetic_distance, row.position, row.allele1, row.allele2), '\t'))
        end
    end
    open(String(prefix) * ".fam", "w") do io
        for row in eachrow(fam)
            println(io, join((row.family_id, row.sample_id, row.paternal_id, row.maternal_id, row.sex, row.phenotype), ' '))
        end
    end
    open(String(prefix) * ".bed", "w") do io
        write(io, _encode_bed(genotypes))
    end
    return prefix
end

function _covariate_matrix(covariates::Union{Nothing,AbstractMatrix{<:Real}}, nobs::Int)
    if covariates === nothing
        return ones(Float64, nobs, 1), _pooled_strings(["intercept"])
    end
    matrix = Matrix{Float64}(covariates)
    size(matrix, 1) == nobs || throw(ArgumentError("covariates must have one row per sample"))
    names = ["intercept"; ["cov_$(index)" for index in 1:size(matrix, 2)]...]
    return hcat(ones(Float64, nobs, 1), matrix), _pooled_strings(names)
end

function _project_out(design::AbstractMatrix{<:Real}, values::AbstractMatrix{<:Real})
    basis = Matrix(qr(Matrix{Float64}(design)).Q)[:, 1:rank(design)]
    return Matrix{Float64}(values) .- basis * (basis' * Matrix{Float64}(values)), basis
end

function _project_out(design::AbstractMatrix{<:Real}, values::AbstractVector{<:Real})
    basis = Matrix(qr(Matrix{Float64}(design)).Q)[:, 1:rank(design)]
    residual = Vector{Float64}(values) .- basis * (basis' * Vector{Float64}(values))
    return residual, basis
end

function _result_from_matrix(genotypes::GenotypeMatrix, beta::Vector{Float64}, se::Vector{Float64}, z::Vector{Float64}, p::Vector{Float64}; method::AbstractString, phenotype_name::AbstractString, covariate_names)
    gene_ids = hasproperty(genotypes.bim, :gene_id) ? String.(genotypes.bim.gene_id) : fill("", length(beta))
    alleles = collect(zip(String.(genotypes.bim.allele1), String.(genotypes.bim.allele2)))
    return GWASResult(_pooled_strings(genotypes.bim.snp_id), _pooled_strings(genotypes.bim.chromosome), Int.(genotypes.bim.position), alleles, _pooled_strings(gene_ids), beta, se, z, p, size(genotypes, 1), _pooled_strings(covariate_names), String(phenotype_name), String(method))
end

function _result_from_matrix(genotypes::AbstractMatrix{<:Real}, beta::Vector{Float64}, se::Vector{Float64}, z::Vector{Float64}, p::Vector{Float64}; method::AbstractString, phenotype_name::AbstractString, covariate_names)
    nsnps = size(genotypes, 2)
    snp_ids = ["snp_$(index)" for index in 1:nsnps]
    chromosomes = fill("", nsnps)
    positions = collect(1:nsnps)
    alleles = fill(("", ""), nsnps)
    gene_ids = fill("", nsnps)
    return GWASResult(_pooled_strings(snp_ids), _pooled_strings(chromosomes), positions, alleles, _pooled_strings(gene_ids), beta, se, z, p, size(genotypes, 1), _pooled_strings(covariate_names), String(phenotype_name), String(method))
end

function DataFrames.DataFrame(result::GWASResult)
    return DataFrames.DataFrame(
        snp_id = result.snp_ids,
        chromosome = result.chromosomes,
        position = result.positions,
        allele1 = [allele[1] for allele in result.alleles],
        allele2 = [allele[2] for allele in result.alleles],
        gene_id = result.gene_ids,
        beta = result.beta,
        standard_error = result.standard_error,
        zscore = result.zscore,
        pvalue = result.pvalue,
        sample_size = fill(result.sample_size, length(result.snp_ids)),
        phenotype_name = PooledArray(fill(result.phenotype_name, length(result.snp_ids))),
        method = PooledArray(fill(result.method, length(result.snp_ids))),
    )
end

function DataFrames.DataFrame(result::MetaAnalysisResult)
    return DataFrames.DataFrame(
        snp_id = result.snp_ids,
        chromosome = result.chromosomes,
        position = result.positions,
        allele1 = [allele[1] for allele in result.alleles],
        allele2 = [allele[2] for allele in result.alleles],
        gene_id = result.gene_ids,
        beta = result.beta,
        standard_error = result.standard_error,
        zscore = result.zscore,
        pvalue = result.pvalue,
        qvalue = result.qvalue,
        tau2 = result.tau2,
        i2 = result.i2,
        study_count = fill(result.study_count, length(result.snp_ids)),
        method = PooledArray(fill(result.method, length(result.snp_ids))),
    )
end

"""
    gwas_linear_scan(genotypes, phenotype; covariates=nothing, phenotype_name="phenotype")

Run a linear-regression GWAS scan across all variants.
"""
function gwas_linear_scan(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, phenotype::AbstractVector{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, phenotype_name::AbstractString="phenotype")
    G = _matrix(genotypes)
    nobs, nsnps = size(G)
    length(phenotype) == nobs || throw(ArgumentError("phenotype length must match sample count"))
    design, covariate_names = _covariate_matrix(covariates, nobs)
    y_res, basis = _project_out(design, phenotype)
    G_res = G .- basis * (basis' * G)
    dof = max(nobs - size(design, 2) - 1, 1)
    beta = zeros(Float64, nsnps)
    se = zeros(Float64, nsnps)
    z = zeros(Float64, nsnps)
    p = ones(Float64, nsnps)
    @threads for snp in 1:nsnps
        x = G_res[:, snp]
        xxt = dot(x, x)
        if !isfinite(xxt) || xxt <= eps()
            continue
        end
        b = dot(x, y_res) / xxt
        resid = y_res .- x .* b
        sigma2 = max(dot(resid, resid) / dof, eps())
        stderr = sqrt(sigma2 / xxt)
        beta[snp] = b
        se[snp] = stderr
        z[snp] = b / stderr
        p[snp] = 2 * ccdf(TDist(dof), abs(z[snp]))
    end
    return genotypes isa GenotypeMatrix ? _result_from_matrix(genotypes, beta, se, z, p; method="linear_scan", phenotype_name=phenotype_name, covariate_names=covariate_names) : _result_from_matrix(G, beta, se, z, p; method="linear_scan", phenotype_name=phenotype_name, covariate_names=covariate_names)
end

function _kinship_from_genotypes(genotypes::AbstractMatrix{<:Real})
    X = Matrix{Float64}(genotypes)
    for column in 1:size(X, 2)
        centered = X[:, column] .- mean(X[:, column])
        spread = std(centered)
        X[:, column] = spread > 0 ? centered ./ spread : centered
    end
    return (X * X') / max(size(X, 2), 1)
end

"""
    gwas_lmm_scan(genotypes, phenotype; covariates=nothing, kinship=nothing, phenotype_name="phenotype")

Run a mixed-model GWAS scan using an optional kinship matrix.
"""
function gwas_lmm_scan(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, phenotype::AbstractVector{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, kinship::Union{Nothing,AbstractMatrix{<:Real}}=nothing, phenotype_name::AbstractString="phenotype")
    G = _matrix(genotypes)
    nobs = size(G, 1)
    length(phenotype) == nobs || throw(ArgumentError("phenotype length must match sample count"))
    K = kinship === nothing ? _kinship_from_genotypes(G) : Matrix{Float64}(kinship)
    symK = Symmetric((K + K') / 2)
    decomp = eigen(symK)
    weights = 1.0 ./ sqrt.(max.(decomp.values .+ 1e-6, 1e-6))
    transform = Diagonal(weights) * decomp.vectors'
    return gwas_linear_scan(transform * G, transform * Vector{Float64}(phenotype); covariates=covariates === nothing ? nothing : transform * Matrix{Float64}(covariates), phenotype_name=phenotype_name)
end

function _pairwise_ld(genotypes::AbstractMatrix{<:Real}, locus_1::Int, locus_2::Int)
    x = Vector{Float64}(genotypes[:, locus_1])
    y = Vector{Float64}(genotypes[:, locus_2])
    valid = isfinite.(x) .& isfinite.(y)
    sum(valid) < 3 && return (0.0, 0.0, 0.0)
    x = x[valid]
    y = y[valid]
    sx = std(x)
    sy = std(y)
    (sx <= eps() || sy <= eps()) && return (0.0, 0.0, 0.0)
    r = cor(x, y)
    d = r * sx * sy
    return (d, r, r^2)
end

"""
    ld_clumping(result, genotypes; p_threshold=5e-8, r2_threshold=0.2, window=250_000)

Prune correlated GWAS hits by linkage disequilibrium.
"""
function ld_clumping(result::GWASResult, genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; p_threshold::Real=5e-8, r2_threshold::Real=0.2, window::Int=250_000)
    G = _matrix(genotypes)
    order = sortperm(result.pvalue)
    kept = Int[]
    for index in order
        result.pvalue[index] <= p_threshold || continue
        blocked = false
        for existing in kept
            result.chromosomes[existing] == result.chromosomes[index] || continue
            abs(result.positions[existing] - result.positions[index]) <= window || continue
            _, _, r2 = _pairwise_ld(G, existing, index)
            if r2 >= r2_threshold
                blocked = true
                break
            end
        end
        blocked || push!(kept, index)
    end
    return GWASResult(result.snp_ids[kept], result.chromosomes[kept], result.positions[kept], result.alleles[kept], result.gene_ids[kept], result.beta[kept], result.standard_error[kept], result.zscore[kept], result.pvalue[kept], result.sample_size, result.covariate_names, result.phenotype_name, result.method)
end

function _summary_lookup(summary_stats)
    hasproperty(summary_stats, :snp_ids) && hasproperty(summary_stats, :beta) || throw(ArgumentError("summary stats must expose snp_ids and beta"))
    return Dict(String(summary_stats.snp_ids[index]) => index for index in eachindex(summary_stats.snp_ids))
end

"""
    calculate_prs(genotypes, summary_stats; standardize=false)

Compute a polygenic risk score from summary statistics.
"""
function calculate_prs(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, summary_stats; standardize::Bool=false)
    G = _matrix(genotypes)
    lookup = _summary_lookup(summary_stats)
    weights = zeros(Float64, size(G, 2))
    if genotypes isa GenotypeMatrix
        for (index, snp_id) in enumerate(String.(genotypes.bim.snp_id))
            if haskey(lookup, snp_id)
                weights[index] = Float64(summary_stats.beta[lookup[snp_id]])
            end
        end
    else
        for index in 1:size(G, 2)
            key = "snp_$(index)"
            if haskey(lookup, key)
                weights[index] = Float64(summary_stats.beta[lookup[key]])
            end
        end
    end
    matrix = standardize ? (G .- mean(G, dims=1)) ./ max.(std(G, dims=1), eps()) : G
    return matrix * weights
end

function _solve_ld_system(ld_matrix::AbstractMatrix{<:Real}, weights::AbstractVector{<:Real}; ridge::Real=0.1)
    A = Matrix{Float64}(ld_matrix) + Float64(ridge) * I
    return A \ Vector{Float64}(weights)
end

"""
    prs_ldpred(genotypes, ld_matrix, summary_stats; ridge=0.1)

Compute an LD-adjusted PRS using ridge regularization.
"""
function prs_ldpred(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, ld_matrix::AbstractMatrix{<:Real}, summary_stats; ridge::Real=0.1)
    weights = hasproperty(summary_stats, :beta) ? Float64.(summary_stats.beta) : throw(ArgumentError("summary stats must expose beta"))
    adjusted = _solve_ld_system(ld_matrix, weights; ridge=ridge)
    return _matrix(genotypes) * adjusted
end

"""
    prs_cross_validation(genotypes, phenotype, summary_stats; p_thresholds=[...], ld_windows=[...])

Choose PRS tuning parameters by cross-validation.
"""
function prs_cross_validation(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, phenotype::AbstractVector{<:Real}, summary_stats; p_thresholds::AbstractVector{<:Real}=[1e-5, 1e-4, 1e-3], ld_windows::AbstractVector{<:Integer}=[50_000, 100_000, 250_000])
    G = _matrix(genotypes)
    best_score = -Inf
    best_threshold = first(p_thresholds)
    best_window = first(ld_windows)
    scores = Dict{Tuple{Float64,Int},Float64}()
    for threshold in p_thresholds, window in ld_windows
        clumped = ld_clumping(summary_stats, G; p_threshold=threshold, window=window)
        prs = calculate_prs(G, clumped)
        score = iszero(std(prs)) || iszero(std(phenotype)) ? 0.0 : abs(cor(prs, phenotype))
        scores[(Float64(threshold), Int(window))] = score
        if score > best_score
            best_score = score
            best_threshold = threshold
            best_window = window
        end
    end
    return (; best_threshold=Float64(best_threshold), best_window=Int(best_window), best_score=best_score, scores=scores)
end

function _common_snps(results::Vector{GWASResult})
    isempty(results) && return String[]
    common = Set(results[1].snp_ids)
    for result in results[2:end]
        intersect!(common, Set(result.snp_ids))
    end
    return [snp for snp in results[1].snp_ids if snp in common]
end

"""
    meta_analyze(results; random_effects=false)

Combine multiple GWAS result tables using fixed- or random-effects meta-analysis.
"""
function meta_analyze(results::Vector{GWASResult}; random_effects::Bool=false)
    isempty(results) && return MetaAnalysisResult(_pooled_strings(String[]), _pooled_strings(String[]), Int[], Tuple{String,String}[], _pooled_strings(String[]), Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], 0, random_effects ? "random_effects" : "fixed_effects")
    snp_ids = _common_snps(results)
    chromosomes = String[]
    positions = Int[]
    alleles = Tuple{String,String}[]
    gene_ids = String[]
    beta = Float64[]
    se = Float64[]
    z = Float64[]
    p = Float64[]
    tau2 = Float64[]
    i2 = Float64[]
    for snp_id in snp_ids
        per_study = [(result, findfirst(==(snp_id), result.snp_ids)) for result in results]
        per_study = [(result, index) for (result, index) in per_study if index !== nothing]
        betas = [Float64(result.beta[index]) for (result, index) in per_study]
        ses = [Float64(result.standard_error[index]) for (result, index) in per_study]
        weights = 1.0 ./ (ses .^ 2)
        fixed_beta = sum(weights .* betas) / sum(weights)
        q_stat = sum(weights .* (betas .- fixed_beta).^2)
        c_stat = sum(weights) - sum(weights .^ 2) / sum(weights)
        tau = random_effects ? max((q_stat - (length(betas) - 1)) / max(c_stat, eps()), 0.0) : 0.0
        weights_star = 1.0 ./ (ses .^ 2 .+ tau)
        meta_beta = sum(weights_star .* betas) / sum(weights_star)
        meta_se = sqrt(1 / sum(weights_star))
        meta_z = meta_beta / meta_se
        meta_p = 2 * ccdf(Normal(), abs(meta_z))
        push!(beta, meta_beta)
        push!(se, meta_se)
        push!(z, meta_z)
        push!(p, meta_p)
        push!(tau2, tau)
        push!(i2, q_stat <= 0 ? 0.0 : max((q_stat - (length(betas) - 1)) / q_stat, 0.0))
        first_result = results[1]
        first_index = findfirst(==(snp_id), first_result.snp_ids)
        push!(chromosomes, first_result.chromosomes[first_index])
        push!(positions, first_result.positions[first_index])
        push!(alleles, first_result.alleles[first_index])
        push!(gene_ids, first_result.gene_ids[first_index])
    end
    qvalue = benjamini_hochberg(p)
    return MetaAnalysisResult(_pooled_strings(snp_ids), _pooled_strings(chromosomes), positions, alleles, _pooled_strings(gene_ids), beta, se, z, p, qvalue, tau2, i2, length(results), random_effects ? "random_effects" : "fixed_effects")
end

function _peak_intervals(peakset::PeakSet)
    intervals = Vector{GenomicInterval}(undef, length(peakset.peaks))
    for (index, peak) in enumerate(peakset.peaks)
        intervals[index] = GenomicInterval(peak.chrom, peak.left, peak.right, '+', Dict{String,Any}("peak_index" => index, "score" => peak.score, "pvalue" => peak.pvalue, "qvalue" => peak.qvalue))
    end
    return intervals
end

"""
    overlap_gwas_peaks(result, peakset; flank=50_000, pvalue_threshold=5e-8)

Find GWAS loci that overlap peak intervals.
"""
function overlap_gwas_peaks(result::GWASResult, peakset::PeakSet; flank::Int=50_000, pvalue_threshold::Real=5e-8)
    intervals = build_collection(_peak_intervals(peakset))
    rows = NamedTuple[]
    for index in eachindex(result.snp_ids)
        result.pvalue[index] <= pvalue_threshold || continue
        query = GenomicInterval(result.chromosomes[index], max(1, result.positions[index] - flank), result.positions[index] + flank, '+', Dict{String,Any}("snp_id" => result.snp_ids[index], "beta" => result.beta[index], "pvalue" => result.pvalue[index]))
        for hit in find_overlaps(query, intervals)
            push!(rows, (
                snp_id = result.snp_ids[index],
                chromosome = result.chromosomes[index],
                position = result.positions[index],
                beta = result.beta[index],
                pvalue = result.pvalue[index],
                peak_chromosome = hit.chrom,
                peak_left = hit.left,
                peak_right = hit.right,
                peak_index = get(hit.metadata, "peak_index", 0),
                peak_score = get(hit.metadata, "score", NaN),
                peak_pvalue = get(hit.metadata, "pvalue", NaN),
                peak_qvalue = get(hit.metadata, "qvalue", NaN),
            ))
        end
    end
    return DataFrame(rows)
end

"""
    gene_based_test(result, network; snp_to_gene=Dict(), pvalue_threshold=5e-8)

Run a gene-based association test using a gene network and SNP-to-gene mapping.
"""
function gene_based_test(result::GWASResult, network::GeneNetwork; snp_to_gene::Dict{String,String}=Dict{String,String}(), pvalue_threshold::Real=5e-8)
    genes = String[]
    for index in eachindex(result.snp_ids)
        result.pvalue[index] <= pvalue_threshold || continue
        gene = result.gene_ids[index]
        isempty(gene) && (gene = get(snp_to_gene, result.snp_ids[index], ""))
        isempty(gene) || push!(genes, gene)
    end
    query_genes = unique(genes)
    terms = EnrichmentTerm[]
    for module_id in sort(unique(network.modules))
        gene_indices = findall(==(module_id), network.modules)
        module_genes = [network.node_to_gene[index] for index in gene_indices]
        push!(terms, EnrichmentTerm("MODULE:$(module_id)", "Module $(module_id)", "MODULE", module_genes, String[]))
    end
    database = build_annotation_database(terms)
    return enrichment_test(query_genes, database; namespace="MODULE")
end

end