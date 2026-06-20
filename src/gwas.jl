# ==============================================================================
# gwas.jl — Genome-wide association studies
#
# Provides PLINK format I/O, linear/LMM GWAS scans, polygenic risk
# score computation (PRS, LDpred), LD clumping, meta-analysis, and
# gene-based testing.
#
# References:
#   - Purcell et al. (2007) AJHG 81(3):559-575 (PLINK)
#   - Vilhjálmsson et al. (2015) AJHG 97(4):576-592 (LDpred)
# ==============================================================================

module GWAS

using Mmap
using DataFrames
using Statistics
using LinearAlgebra
using Random
using Distributions
using Base.Threads
using PooledArrays
using SparseArrays
using Arrow
using Downloads
using JSON

using ..BioToolkit: AbstractAnalysisResult, ProvenanceContext, ProvenanceParams, ResultProvenance, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_record, provenance_result!, register_provenance!

@inline function _register_gwas_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end
using ..DifferentialExpression: benjamini_hochberg
using ..GenomicRanges: GenomicInterval, build_collection, find_overlaps
using ..Epigenetics: PeakSet, Peak
using ..SystemsBio: GeneNetwork
using ..Enrichment: EnrichmentDatabase, EnrichmentTerm, build_annotation_database, enrichment_test

export GenotypeMatrix, GWASResult, MetaAnalysisResult
export read_plink, write_plink, gwas_linear_scan, gwas_lmm_scan
export calculate_prs, prs_ldpred, prs_cross_validation, ld_clumping, meta_analyze
export overlap_gwas_peaks, gene_based_test
export BedReader, BgenReader, BgenVariant
export read_bed_genotypes, read_bgen, write_bgen
export calculate_maf, hwe_exact, calculate_hwe_pvalues, hwe_filter
export calculate_missingness, missingness_filter, info_score_filter, gwas_qc_report
export genomic_control_lambda, apply_genomic_control
export calculate_ld_matrix, calculate_grm, loco_projection
export gwas_logistic_scan, gwas_gxe_interaction
export to_plink_dataframe, write_plink_sumstats
export gwas_migration_guide
export LDSCResult, ldsc_heritability, estimate_heritability_greml, partitioned_heritability, ldsc_genetic_correlation
export SuSiEResult, fine_map_susie, calculate_credible_set, posterior_inclusion_probability
export MRResult, mr_two_sample, mr_egger, mr_pleiotropy_test
export conditional_analysis, cojo_stepwise, joint_analysis
export gwas_pca, project_pca, calculate_ibd, calculate_king_kinship, detect_related_pairs
export coloc_result, coloc_abf
export ihs_score, xp_ehh_score, fst_outlier_test
export rank_inverse_normal, rank_inverse_normal!
export calculate_sex_check, calculate_heterozygosity_outliers, sample_qc_report
export normalise_chromosome, filter_variants, filter_samples, merge_genotype_matrices
export read_vcf, write_vcf, read_oxford_gen, compute_ld_scores, prune_ld, flip_alleles, harmonise_alleles
export sample_prune_relatedness, ancestry_inference, manhattan_data, manhattan_plot, qq_data, qq_plot, locus_zoom
export score_test_linear, likelihood_ratio_test, gwas_survival_scan, gwas_ordinal_scan, gwas_multivariate_scan
export burden_test, skat_test, conditional_fdr, storey_pi0_estimate, gwas_power_calculation, permutation_pvalue, genomic_sem_fit
export twas_scan, smr_test, liftover, functional_annotation, ld_expand_credible_set, estimate_effective_n, he_regression_variance_components, reml_variance_components
export simulate_genotypes, simulate_phenotype, dosage_to_hardcall, info_score_from_dosage, popcorn_genetic_correlation, ebi_lookup
export save_gwas_result, load_gwas_result

const _GWAS_MIGRATION_GUIDE = """
### Migrating from SnpArrays.jl and JWAS.jl

BioToolkit provides a pure-Julia GWAS stack with explicit data structures.

| Task | SnpArrays/JWAS | BioToolkit |
|------|----------------|------------|
| Read PLINK prefix | `SnpArray(prefix)` | `read_plink(prefix)` |
| Read BGEN stream | package-specific wrappers | `BgenReader` and `read_bgen` |
| MAF filtering | custom code | `calculate_maf` |
| HWE exact test | custom code | `hwe_exact`, `calculate_hwe_pvalues` |
| Missingness / call rate | PLINK `--geno/--mind` | `calculate_missingness` |
| Full LD matrix | package-specific | `calculate_ld_matrix` |
| GRM construction | package-specific | `calculate_grm` |
| Binary trait GWAS | GLM wrappers | `gwas_logistic_scan` |
"""

"""
    gwas_migration_guide()

Return a markdown migration guide comparing GWAS workflows with SnpArrays/JWAS.
"""
gwas_migration_guide() = _GWAS_MIGRATION_GUIDE

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

GenotypeMatrix(decoded::AbstractMatrix{<:Real}, bim::DataFrame, fam::DataFrame; prefix::String="") =
    GenotypeMatrix(UInt8[], Matrix{Float64}(decoded), bim, fam, String(prefix))

"""
    BedReader

Streaming parser for tab-delimited BED-like dosage matrices with columns:
`CHR POS ID REF ALT SAMPLE1 SAMPLE2 ...`.

The reader yields one marker at a time and avoids loading the full file in memory.
"""
mutable struct BedReader
    io::IOStream
    header::Vector{String}
    marker_ids::Vector{String}
    sample_ids::Vector{String}
    n_samples::Int
    n_snps::Int
    data_start::Int64
    closed::Bool
end

"""
    BgenVariant

Single BGEN marker decoded as allele dosages.
"""
struct BgenVariant
    chromosome::String
    position::Int
    snp_id::String
    rsid::String
    alleles::Vector{String}
    dosage::Vector{Float64}
end

"""
    BgenReader

Streaming reader for BGEN files (layout 2, uncompressed blocks).
"""
mutable struct BgenReader
    io::IOStream
    path::String
    header::NamedTuple
    n_samples::Int
    n_snps::Int
    sample_ids::Vector{String}
    compression::Int
    layout::Int
    variant_index::Int
    closed::Bool
end

"""
    GWASResult

Association scan result table for a single GWAS analysis.
"""
struct GWASResult <: AbstractAnalysisResult
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
    provenance::ResultProvenance
end

GWASResult(snp_ids, chromosomes, positions, alleles, gene_ids, beta, standard_error, zscore, pvalue, sample_size, covariate_names, phenotype_name, method) =
    GWASResult(snp_ids, chromosomes, positions, alleles, gene_ids, beta, standard_error, zscore, pvalue, sample_size, covariate_names, phenotype_name, method, provenance_record("GWASResult", "gwas"))

"""
    MetaAnalysisResult

Meta-analysis summary across multiple GWAS result tables.
"""
struct MetaAnalysisResult <: AbstractAnalysisResult
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
    provenance::ResultProvenance
end

MetaAnalysisResult(snp_ids, chromosomes, positions, alleles, gene_ids, beta, standard_error, zscore, pvalue, qvalue, tau2, i2, study_count, method) =
    MetaAnalysisResult(snp_ids, chromosomes, positions, alleles, gene_ids, beta, standard_error, zscore, pvalue, qvalue, tau2, i2, study_count, method, provenance_record("MetaAnalysisResult", "gwas"))

_pooled_strings(values) = PooledArray(String.(values))

function GWASResult(
    snp_ids::AbstractVector{<:String},
    chromosomes::AbstractVector{<:String},
    positions::AbstractVector{<:Integer},
    alleles::AbstractVector{<:Tuple},
    gene_ids::AbstractVector{<:String},
    beta::AbstractVector{<:Real},
    standard_error::AbstractVector{<:Real},
    zscore::AbstractVector{<:Real},
    pvalue::AbstractVector{<:Real},
    sample_size::Integer,
    covariate_names::AbstractVector{<:String},
    phenotype_name::String,
    method::String)
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
        provenance_record("GWASResult", "gwas"))
end

function MetaAnalysisResult(
    snp_ids::AbstractVector{<:String},
    chromosomes::AbstractVector{<:String},
    positions::AbstractVector{<:Integer},
    alleles::AbstractVector{<:Tuple},
    gene_ids::AbstractVector{<:String},
    beta::AbstractVector{<:Real},
    standard_error::AbstractVector{<:Real},
    zscore::AbstractVector{<:Real},
    pvalue::AbstractVector{<:Real},
    qvalue::AbstractVector{<:Real},
    tau2::AbstractVector{<:Real},
    i2::AbstractVector{<:Real},
    study_count::Integer,
    method::String)
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
        provenance_record("MetaAnalysisResult", "gwas"))
end

Base.size(genotypes::GenotypeMatrix) = size(genotypes.decoded)
Base.getindex(genotypes::GenotypeMatrix, row::Int, column::Int) = genotypes.decoded[row, column]
Base.IndexStyle(::Type{GenotypeMatrix}) = IndexCartesian()
Base.eltype(::Type{GenotypeMatrix}) = Float64
Base.Matrix(genotypes::GenotypeMatrix) = copy(genotypes.decoded)

_matrix(genotypes::GenotypeMatrix) = genotypes.decoded
_matrix(genotypes::AbstractMatrix{<:Real}) = Matrix{Float64}(genotypes)

function normalise_chromosome(chr)
    value = uppercase(strip(String(chr)))
    isempty(value) && return ""
    startswith(value, "CHR") && (value = value[4:end])

    if value in ("23", "X")
        return "X"
    elseif value in ("24", "Y")
        return "Y"
    elseif value in ("25", "XY")
        return "XY"
    elseif value in ("26", "M", "MT")
        return "MT"
    end

    parsed = tryparse(Int, value)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, parsed === nothing ? value : string(parsed), "normalise_chromosome")
end

function _indices_from_mask(mask, n::Int)
    if mask isa AbstractVector{Bool}
        length(mask) == n || throw(DimensionMismatch("mask length must match axis length"))
        return findall(mask)
    elseif mask isa AbstractVector{<:Integer}
        idx = Int.(mask)
        all(i -> 1 <= i <= n, idx) || throw(BoundsError("index mask out of bounds"))
        return idx
    end
    throw(ArgumentError("mask must be a Bool vector or integer index vector"))
end

@inline _thread_enabled(nwork::Int, multi_thread::Bool; min_grain::Int=512) =
    multi_thread && Threads.nthreads() > 1 && nwork >= min_grain

@inline function _progress_tick!(verbose::Bool, label::String, index::Int, total::Int; step::Int=50_000)
    verbose || return nothing
    if index == 1 || index == total || (step > 0 && index % step == 0)
        @info("$label progress", done=index, total=total)
    end
    return nothing
end

function _validate_gwas_inputs(
    G::AbstractMatrix{<:Real};
    phenotype::Union{Nothing,AbstractVector{<:Real}}=nothing,
    covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
    kinship::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
    ld_matrix::Union{Nothing,AbstractMatrix{<:Real}}=nothing)
    nobs, nsnps = size(G)
    nobs > 0 || throw(ArgumentError("genotype matrix must contain at least one sample"))
    nsnps > 0 || throw(ArgumentError("genotype matrix must contain at least one variant"))

    if phenotype !== nothing
        length(phenotype) == nobs || throw(DimensionMismatch("phenotype length must match sample count"))
    end
    if covariates !== nothing
        size(covariates, 1) == nobs || throw(DimensionMismatch("covariates must have one row per sample"))
    end
    if kinship !== nothing
        size(kinship, 1) == nobs == size(kinship, 2) || throw(DimensionMismatch("kinship matrix must be square with side length equal to sample count"))
    end
    if ld_matrix !== nothing
        size(ld_matrix, 1) == nsnps == size(ld_matrix, 2) || throw(DimensionMismatch("LD matrix must be square with side length equal to variant count"))
    end
    return nothing
end

function _table_rows(path::String)
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

function _read_bim(path::String)
    rows = _table_rows(path)
    return DataFrame(
        chromosome = [normalise_chromosome(row[1]) for row in rows],
        snp_id = [row[2] for row in rows],
        genetic_distance = parse.(Float64, [row[3] for row in rows]),
        position = parse.(Int, [row[4] for row in rows]),
        allele1 = [row[5] for row in rows],
        allele2 = [row[6] for row in rows])
end

function _read_fam(path::String)
    rows = _table_rows(path)
    return DataFrame(
        family_id = [row[1] for row in rows],
        sample_id = [row[2] for row in rows],
        paternal_id = [row[3] for row in rows],
        maternal_id = [row[4] for row in rows],
        sex = parse.(Int, [row[5] for row in rows]),
        phenotype = parse.(Float64, [row[6] for row in rows]))
end

"""
    filter_variants(gm, mask)

Subset variants while keeping BIM metadata aligned.
"""
function filter_variants(gm::GenotypeMatrix, mask)
    idx = _indices_from_mask(mask, size(gm, 2))
    bim = gm.bim[idx, :]
    decoded = gm.decoded[:, idx]
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GenotypeMatrix(UInt8[], Matrix{Float64}(decoded), DataFrame(bim), copy(gm.fam), gm.prefix), "filter_variants")
end

"""
    filter_samples(gm, mask)

Subset samples while keeping FAM metadata aligned.
"""
function filter_samples(gm::GenotypeMatrix, mask)
    idx = _indices_from_mask(mask, size(gm, 1))
    fam = gm.fam[idx, :]
    decoded = gm.decoded[idx, :]
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GenotypeMatrix(UInt8[], Matrix{Float64}(decoded), copy(gm.bim), DataFrame(fam), gm.prefix), "filter_samples")
end

"""
    merge_genotype_matrices(gm_list; by=:variants)

Merge multiple `GenotypeMatrix` objects either by concatenating variants
or by concatenating samples.
"""
function merge_genotype_matrices(gm_list::AbstractVector{<:GenotypeMatrix}; by::Symbol=:variants)
    isempty(gm_list) && throw(ArgumentError("gm_list cannot be empty"))
    if by == :variants
        sample_ids = String.(gm_list[1].fam.sample_id)
        for gm in gm_list[2:end]
            String.(gm.fam.sample_id) == sample_ids || throw(DimensionMismatch("all matrices must share identical sample IDs for by=:variants"))
        end
        decoded = hcat([gm.decoded for gm in gm_list]...)
        bim = vcat([gm.bim for gm in gm_list]...)
        fam = copy(gm_list[1].fam)
        return GenotypeMatrix(UInt8[], Matrix{Float64}(decoded), DataFrame(bim), fam, gm_list[1].prefix)
    elseif by == :samples
        snp_ids = String.(gm_list[1].bim.snp_id)
        for gm in gm_list[2:end]
            String.(gm.bim.snp_id) == snp_ids || throw(DimensionMismatch("all matrices must share identical SNP IDs and order for by=:samples"))
        end
        decoded = vcat([gm.decoded for gm in gm_list]...)
        fam = vcat([gm.fam for gm in gm_list]...)
        bim = copy(gm_list[1].bim)

        return GenotypeMatrix(UInt8[], Matrix{Float64}(decoded), bim, DataFrame(fam), gm_list[1].prefix)
    end
    throw(ArgumentError("by must be :variants or :samples"))
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

@inline function _split_bed_line(line::AbstractString)
    stripped = strip(line)
    isempty(stripped) && return String[]
    return occursin('\t', stripped) ? split(stripped, '\t'; keepempty=false) : split(stripped)
end

@inline function _parse_dosage_or_nan(token::AbstractString)
    normalized = lowercase(strip(String(token)))
    if normalized in ("na", "nan", ".", "missing", "")
        return NaN
    end
    return parse(Float64, normalized)
end

function _read_bed_index(path::String)
    isfile(path) || return nothing
    ids = String[]
    open(path, "r") do io
        for line in eachline(io)
            stripped = strip(line)
            isempty(stripped) && continue
            push!(ids, stripped)
        end
    end
    return isempty(ids) ? nothing : ids
end

"""
    BedReader(path; index_path=nothing)

Open a BED-like tabular dosage file for streaming marker-by-marker iteration.

If an index file is provided (or `<path>.snpidx` exists), marker IDs are read from
the index and the constructor avoids a full scan to count rows.
"""
function BedReader(path::String; index_path::Union{Nothing,String}=nothing)
    io = open(path, "r")
    header_line = ""
    while !eof(io)
        candidate = strip(readline(io))
        isempty(candidate) && continue
        startswith(candidate, "#") && (candidate = strip(candidate[2:end]))
        isempty(candidate) && continue
        header_line = candidate
        break
    end
    isempty(header_line) && throw(ArgumentError("empty BED-like dosage file"))
    header = _split_bed_line(header_line)
    length(header) >= 6 || throw(ArgumentError("BED-like dosage files must provide CHR POS ID REF ALT and at least one sample column"))
    sample_ids = String.(header[6:end])
    data_start = position(io)

    effective_index_path = if index_path === nothing
        auto = String(path) * ".snpidx"
        isfile(auto) ? auto : nothing
    else
        index_path
    end

    marker_ids = effective_index_path === nothing ? nothing : _read_bed_index(effective_index_path)
    n_snps = 0

    if marker_ids === nothing
        marker_ids = String[]
        for line in eachline(io)
            parts = _split_bed_line(line)
            isempty(parts) && continue
            n_snps += 1
            push!(marker_ids, length(parts) >= 3 ? String(parts[3]) : "snp_$(n_snps)")
        end
    else
        n_snps = length(marker_ids)
    end

    seek(io, data_start)
    return BedReader(io, String.(header[1:5]), marker_ids, sample_ids, length(sample_ids), n_snps, data_start, false)
end

function Base.close(reader::BedReader)
    if !reader.closed
        close(reader.io)
        reader.closed = true
    end
    return nothing
end

Base.length(reader::BedReader) = reader.n_snps

function Base.iterate(reader::BedReader, state::Int=1)
    reader.closed && return nothing
    if state > reader.n_snps
        close(reader)
        return nothing
    end

    while !eof(reader.io)
        parts = _split_bed_line(readline(reader.io))
        isempty(parts) && continue
        length(parts) >= 5 + reader.n_samples || throw(ArgumentError("malformed BED-like dosage row with insufficient sample columns"))
        dosage = Vector{Float64}(undef, reader.n_samples)
        for sample in 1:reader.n_samples
            dosage[sample] = _parse_dosage_or_nan(parts[5 + sample])
        end
        variant = (
            chromosome = normalise_chromosome(parts[1]),
            position = parse(Int, parts[2]),
            snp_id = String(parts[3]),
            ref = String(parts[4]),
            alt = String(parts[5]),
            dosage = dosage)
        return variant, state + 1
    end

    close(reader)
    return nothing
end

"""
    read_bed_genotypes(path)

Read a BED-like dosage table into a `GenotypeMatrix`.
"""
function read_bed_genotypes(path::String)
    reader = BedReader(path)
    matrix = fill(NaN, reader.n_samples, reader.n_snps)
    chromosomes = String[]
    positions = Int[]
    snp_ids = String[]
    allele1 = String[]
    allele2 = String[]

    col = 0
    try
        for variant in reader
            col += 1
            matrix[:, col] = variant.dosage
            push!(chromosomes, variant.chromosome)
            push!(positions, variant.position)
            push!(snp_ids, variant.snp_id)
            push!(allele1, variant.ref)
            push!(allele2, variant.alt)
        end
    finally
        close(reader)
    end

    if col < size(matrix, 2)
        matrix = matrix[:, 1:col]
    end

    bim = DataFrame(
        chromosome = chromosomes,
        snp_id = snp_ids,
        genetic_distance = zeros(Float64, length(snp_ids)),
        position = positions,
        allele1 = allele1,
        allele2 = allele2)
    fam = DataFrame(
        family_id = fill("0", reader.n_samples),
        sample_id = copy(reader.sample_ids),
        paternal_id = fill("0", reader.n_samples),
        maternal_id = fill("0", reader.n_samples),
        sex = fill(0, reader.n_samples),
        phenotype = fill(NaN, reader.n_samples))
    prefix = endswith(lowercase(path), ".bed") ? path[1:end-4] : path
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GenotypeMatrix(matrix, bim, fam; prefix=prefix), "read_bed_genotypes")
end

@inline function _read_u16_le(io::IO)
    bytes = read(io, 2)
    length(bytes) == 2 || throw(EOFError())
    return UInt16(bytes[1]) | (UInt16(bytes[2]) << 8)
end

@inline function _read_u32_le(io::IO)
    bytes = read(io, 4)
    length(bytes) == 4 || throw(EOFError())
    return UInt32(bytes[1]) | (UInt32(bytes[2]) << 8) | (UInt32(bytes[3]) << 16) | (UInt32(bytes[4]) << 24)
end

@inline function _write_u16_le(io::IO, value::Integer)
    raw = UInt16(value)
    write(io, UInt8(raw & 0xff))
    write(io, UInt8((raw >> 8) & 0xff))
    return nothing
end

@inline function _write_u32_le(io::IO, value::Integer)
    raw = UInt32(value)
    write(io, UInt8(raw & 0xff))
    write(io, UInt8((raw >> 8) & 0xff))
    write(io, UInt8((raw >> 16) & 0xff))
    write(io, UInt8((raw >> 24) & 0xff))
    return nothing
end

@inline function _bgen_probability_count(ploidy::Int, n_alleles::Int, phased::Bool)
    if phased
        return ploidy * max(n_alleles - 1, 0)
    end
    return Int(binomial(ploidy + n_alleles - 1, n_alleles - 1) - 1)
end

function _read_packed_bits!(payload::Vector{UInt8}, bit_offset::Base.RefValue{Int}, nbits::Int)
    value = UInt32(0)
    for bit in 0:(nbits - 1)
        offset = bit_offset[] + bit
        byte_index = (offset >>> 3) + 1
        bit_index = offset & 0x07
        byte_index <= length(payload) || throw(EOFError())
        value |= UInt32((payload[byte_index] >>> bit_index) & 0x01) << bit
    end
    bit_offset[] += nbits
    return Int(value)
end

function _decode_bgen_layout2_dosage(block::Vector{UInt8}, n_samples::Int, n_alleles::Int)
    io = IOBuffer(block)
    block_samples = Int(_read_u32_le(io))
    block_samples == n_samples || throw(ArgumentError("BGEN block sample count does not match header"))
    block_alleles = Int(_read_u16_le(io))
    block_alleles == n_alleles || throw(ArgumentError("BGEN block allele count does not match marker metadata"))

    _min_ploidy = Int(read(io, UInt8))
    _max_ploidy = Int(read(io, UInt8))
    ploidy_bytes = read(io, n_samples)
    phased_flag = Int(read(io, UInt8))
    bits_per_probability = Int(read(io, UInt8))
    bits_per_probability > 0 || throw(ArgumentError("invalid BGEN bits-per-probability field"))
    bits_per_probability <= 16 || throw(ArgumentError("BGEN bits-per-probability > 16 is unsupported"))

    probability_payload = read(io)
    max_probability = (1 << bits_per_probability) - 1
    dosage = fill(NaN, n_samples)
    bit_offset = Ref(0)

    for sample in 1:n_samples
        ploidy = Int(ploidy_bytes[sample] & 0x3f)
        missing = (ploidy_bytes[sample] & 0x80) != 0
        phased = phased_flag == 1
        n_probabilities = _bgen_probability_count(ploidy, n_alleles, phased)

        if missing
            bit_offset[] += n_probabilities * bits_per_probability
            continue
        end

        if n_alleles == 2 && ploidy == 2
            if phased
                p_hap1_alt = _read_packed_bits!(probability_payload, bit_offset, bits_per_probability) / max_probability
                p_hap2_alt = _read_packed_bits!(probability_payload, bit_offset, bits_per_probability) / max_probability
                dosage[sample] = p_hap1_alt + p_hap2_alt
            else
                p_hom_ref = _read_packed_bits!(probability_payload, bit_offset, bits_per_probability) / max_probability
                p_het = _read_packed_bits!(probability_payload, bit_offset, bits_per_probability) / max_probability
                p_hom_alt = max(1.0 - p_hom_ref - p_het, 0.0)
                dosage[sample] = p_het + 2.0 * p_hom_alt
            end
        else
            for _ in 1:n_probabilities
                _read_packed_bits!(probability_payload, bit_offset, bits_per_probability)
            end
            dosage[sample] = NaN
        end
    end

    return dosage
end

"""
    BgenReader(path)

Create a streaming BGEN reader.

Pure-Julia decoding targets layout-2 blocks. zlib/zstd decompression is
available when optional codecs are installed.
"""
function BgenReader(path::String)
    io = open(path, "r")
    try
        offset = Int(_read_u32_le(io))
        header_length = Int(_read_u32_le(io))
        n_snps = Int(_read_u32_le(io))
        n_samples = Int(_read_u32_le(io))
        magic = String(read(io, 4))
        lowercase(magic) == "bgen" || throw(ArgumentError("not a BGEN file"))

        free_data_length = max(header_length - 20, 0)
        free_data_length > 0 && read(io, free_data_length)
        flags = _read_u32_le(io)
        compression = Int(flags & UInt32(0x03))
        layout = Int((flags >> 2) & UInt32(0x0f))
        sample_ids_present = ((flags >> 31) & UInt32(0x01)) == UInt32(0x01)

        sample_ids = ["sample_$(index)" for index in 1:n_samples]
        if sample_ids_present
            sample_block_length = Int(_read_u32_le(io))
            sample_block_start = position(io)
            block_samples = Int(_read_u32_le(io))
            block_samples == n_samples || throw(ArgumentError("BGEN sample identifier block does not match header sample count"))
            sample_ids = Vector{String}(undef, n_samples)
            for index in 1:n_samples
                length_id = Int(_read_u16_le(io))
                sample_ids[index] = String(read(io, length_id))
            end
            consumed = position(io) - sample_block_start
            remaining = sample_block_length - consumed
            remaining > 0 && read(io, remaining)
        end

        header = (
            offset = offset,
            header_length = header_length,
            magic = magic,
            flags = flags,
            sample_ids_present = sample_ids_present)
        return BgenReader(io, path, header, n_samples, n_snps, sample_ids, compression, layout, 1, false)
    catch
        close(io)
        rethrow()
    end
end

function Base.close(reader::BgenReader)
    if !reader.closed
        close(reader.io)
        reader.closed = true
    end
    return nothing
end

Base.length(reader::BgenReader) = reader.n_snps

function _load_optional_module(module_name::Symbol)
    isdefined(@__MODULE__, module_name) && return getfield(@__MODULE__, module_name)
    try
        Base.eval(@__MODULE__, Expr(:import, module_name))
    catch
        return nothing
    end
    return getfield(@__MODULE__, module_name)
end

function _decode_bgen_payload(payload::Vector{UInt8}, compression::Int, expected_uncompressed::Int)
    if compression == 0
        decoded = payload
    elseif compression == 1
        codec = _load_optional_module(:CodecZlib)
        codec === nothing && throw(ArgumentError("compressed BGEN (zlib) requires CodecZlib.jl"))
        stream = codec.ZlibDecompressorStream(IOBuffer(payload))
        try
            decoded = read(stream)
        finally
            close(stream)
        end
    elseif compression == 2
        codec = _load_optional_module(:CodecZstd)
        codec === nothing && throw(ArgumentError("compressed BGEN (zstd) requires CodecZstd.jl"))
        stream = codec.ZstdDecompressorStream(IOBuffer(payload))
        try
            decoded = read(stream)
        finally
            close(stream)
        end
    else
        throw(ArgumentError("unsupported BGEN compression mode $(compression)"))
    end

    if expected_uncompressed > 0 && length(decoded) != expected_uncompressed
        throw(ArgumentError("BGEN decoded payload size does not match block header"))
    end
    return decoded
end

function Base.iterate(reader::BgenReader, state::Int=1)
    reader.closed && return nothing
    if state > reader.n_snps
        close(reader)
        return nothing
    end
    eof(reader.io) && (close(reader); return nothing)

    snp_id_length = Int(_read_u16_le(reader.io))
    snp_id = String(read(reader.io, snp_id_length))
    rsid_length = Int(_read_u16_le(reader.io))
    rsid = String(read(reader.io, rsid_length))
    chrom_length = Int(_read_u16_le(reader.io))
    chromosome = normalise_chromosome(String(read(reader.io, chrom_length)))
    position = Int(_read_u32_le(reader.io))
    n_alleles = Int(_read_u16_le(reader.io))
    n_alleles >= 2 || throw(ArgumentError("BGEN marker must contain at least two alleles"))

    alleles = Vector{String}(undef, n_alleles)
    for allele in 1:n_alleles
        allele_length = Int(_read_u32_le(reader.io))
        alleles[allele] = String(read(reader.io, allele_length))
    end

    reader.layout == 2 || throw(ArgumentError("only BGEN layout 2 is supported"))

    compressed_size = Int(_read_u32_le(reader.io))
    compressed_size >= 0 || throw(ArgumentError("invalid BGEN compressed block size"))
    uncompressed_size = Int(_read_u32_le(reader.io))
    payload_compressed = read(reader.io, compressed_size)
    length(payload_compressed) == compressed_size || throw(EOFError())
    payload = _decode_bgen_payload(payload_compressed, reader.compression, uncompressed_size)

    dosage = _decode_bgen_layout2_dosage(payload, reader.n_samples, n_alleles)
    reader.variant_index = state + 1
    return BgenVariant(chromosome, position, snp_id, rsid, alleles, dosage), state + 1
end

"""
    read_bgen(path)

Read a BGEN dataset into a `GenotypeMatrix`.
"""
function read_bgen(path::String)
    reader = BgenReader(path)
    matrix = fill(NaN, reader.n_samples, reader.n_snps)
    chromosomes = String[]
    positions = Int[]
    snp_ids = String[]
    allele1 = String[]
    allele2 = String[]

    col = 0
    try
        for variant in reader
            col += 1
            matrix[:, col] = variant.dosage
            push!(chromosomes, variant.chromosome)
            push!(positions, variant.position)
            push!(snp_ids, isempty(variant.snp_id) ? variant.rsid : variant.snp_id)
            push!(allele1, length(variant.alleles) >= 1 ? variant.alleles[1] : "")
            push!(allele2, length(variant.alleles) >= 2 ? variant.alleles[2] : "")
        end
    finally
        close(reader)
    end

    if col < size(matrix, 2)
        matrix = matrix[:, 1:col]
    end

    bim = DataFrame(
        chromosome = chromosomes,
        snp_id = snp_ids,
        genetic_distance = zeros(Float64, length(snp_ids)),
        position = positions,
        allele1 = allele1,
        allele2 = allele2)
    fam = DataFrame(
        family_id = fill("0", reader.n_samples),
        sample_id = copy(reader.sample_ids),
        paternal_id = fill("0", reader.n_samples),
        maternal_id = fill("0", reader.n_samples),
        sex = fill(0, reader.n_samples),
        phenotype = fill(NaN, reader.n_samples))
    prefix = endswith(lowercase(path), ".bgen") ? path[1:end-5] : path
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GenotypeMatrix(matrix, bim, fam; prefix=prefix), "read_bgen")
end

function _encode_bgen_layout2_payload(dosage::AbstractVector{<:Real}; bits_per_probability::Int=8)
    bits_per_probability == 8 || throw(ArgumentError("write_bgen currently supports bits_per_probability=8 only"))
    n_samples = length(dosage)

    payload_io = IOBuffer()
    _write_u32_le(payload_io, n_samples)
    _write_u16_le(payload_io, 2)  # bi-allelic
    write(payload_io, UInt8(2))   # min ploidy
    write(payload_io, UInt8(2))   # max ploidy

    ploidy_bytes = fill(UInt8(0x02), n_samples)  # diploid, observed
    probability_bytes = UInt8[]
    max_q = 255.0

    @inbounds for i in 1:n_samples
        value = Float64(dosage[i])
        if !isfinite(value)
            ploidy_bytes[i] = UInt8(0x82)  # missing + diploid
            push!(probability_bytes, UInt8(0x00), UInt8(0x00))
            continue
        end

        d = clamp(value, 0.0, 2.0)
        p0 = d <= 1.0 ? 1.0 - d : 0.0
        p1 = d <= 1.0 ? d : 2.0 - d
        p0 = clamp(p0, 0.0, 1.0)
        p1 = clamp(p1, 0.0, 1.0 - p0)

        q0_int = clamp(round(Int, p0 * max_q), 0, Int(max_q))
        q1_int = clamp(round(Int, p1 * max_q), 0, Int(max_q))
        if q0_int + q1_int > Int(max_q)
            q1_int = Int(max_q) - q0_int
        end
        q0 = UInt8(q0_int)
        q1 = UInt8(q1_int)
        push!(probability_bytes, q0, q1)
    end

    write(payload_io, ploidy_bytes)
    write(payload_io, UInt8(0))                    # unphased
    write(payload_io, UInt8(bits_per_probability))
    write(payload_io, probability_bytes)
    return take!(payload_io)
end

"""
    write_bgen(path, gm; bits_per_probability=8)

Write a minimal layout-2 BGEN file (uncompressed) from a `GenotypeMatrix`.
"""
function write_bgen(path::String, gm::GenotypeMatrix; bits_per_probability::Int=8)
    G = gm.decoded
    n_samples, n_snps = size(G)
    sample_ids = String.(gm.fam.sample_id)
    length(sample_ids) == n_samples || throw(DimensionMismatch("sample metadata does not match genotype row count"))

    sample_block_length = 4 + sum(2 + ncodeunits(id) for id in sample_ids)
    # BGEN offset points to the first variant block from file start.
    # Current layout: 24-byte fixed header + 4-byte sample block length field + sample block payload.
    first_variant_offset = 24 + 4 + sample_block_length

    open(path, "w") do io
        _write_u32_le(io, first_variant_offset)
        _write_u32_le(io, 20)
        _write_u32_le(io, n_snps)
        _write_u32_le(io, n_samples)
        write(io, "bgen")
        _write_u32_le(io, 0x80000008)  # sample IDs present, layout 2, uncompressed

        _write_u32_le(io, sample_block_length)
        _write_u32_le(io, n_samples)
        for id in sample_ids
            _write_u16_le(io, ncodeunits(id))
            write(io, id)
        end

        for snp in 1:n_snps
            snp_id = String(gm.bim.snp_id[snp])
            rsid = isempty(snp_id) ? "snp_$(snp)" : snp_id
            chr = normalise_chromosome(gm.bim.chromosome[snp])
            pos = Int(gm.bim.position[snp])
            ref = String(gm.bim.allele1[snp])
            alt = String(gm.bim.allele2[snp])

            _write_u16_le(io, ncodeunits(snp_id)); write(io, snp_id)
            _write_u16_le(io, ncodeunits(rsid)); write(io, rsid)
            _write_u16_le(io, ncodeunits(chr)); write(io, chr)
            _write_u32_le(io, pos)
            _write_u16_le(io, 2)
            _write_u32_le(io, ncodeunits(ref)); write(io, ref)
            _write_u32_le(io, ncodeunits(alt)); write(io, alt)

            payload = _encode_bgen_layout2_payload(@view(G[:, snp]); bits_per_probability=bits_per_probability)
            _write_u32_le(io, length(payload))
            _write_u32_le(io, length(payload))
            write(io, payload)
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, path, "write_bgen")
end

function write_bgen(path::String, genotypes::AbstractMatrix{<:Real}, bim::DataFrame, fam::DataFrame; bits_per_probability::Int=8)
    gm = GenotypeMatrix(UInt8[], Matrix{Float64}(genotypes), copy(bim), copy(fam), replace(path, r"\.bgen$"i => ""))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, write_bgen(path, gm; bits_per_probability=bits_per_probability), "write_bgen")
end

@inline function _vcf_gt_to_dosage(gt::AbstractString)
    gt in (".", "./.", ".|.") && return NaN
    sep = occursin('|', gt) ? '|' : '/'
    alleles = split(gt, sep)
    length(alleles) == 2 || return NaN
    any(a -> a == "." || isempty(a), alleles) && return NaN
    a1 = tryparse(Int, alleles[1])
    a2 = tryparse(Int, alleles[2])
    (a1 === nothing || a2 === nothing) && return NaN
    # Multi-allelic records are represented as NaN in this lightweight parser.
    ((a1 > 1) || (a2 > 1)) && return NaN
    return Float64(a1 + a2)
end

@inline function _vcf_gt_to_dosage(format_field::AbstractString, sample_field::AbstractString)
    format_tokens = split(String(format_field), ':')
    gt_index = findfirst(==("GT"), format_tokens)
    gt_index === nothing && return _vcf_gt_to_dosage(sample_field)

    sample_tokens = split(String(sample_field), ':')
    gt_index > length(sample_tokens) && return NaN
    return _vcf_gt_to_dosage(sample_tokens[gt_index])
end

"""
    read_vcf(path)

Read a text VCF file into a `GenotypeMatrix` using GT dosages (0/1/2, NaN for missing).
"""
function read_vcf(path::String)
    chromosomes = String[]
    positions = Int[]
    snp_ids = String[]
    allele1 = String[]
    allele2 = String[]
    dosage_cols = Vector{Vector{Float64}}()
    sample_ids = String[]

    function _parse_stream!(io::IO)
        for line in eachline(io)
            startswith(line, "##") && continue
            if startswith(line, "#CHROM")
                header = split(chomp(line), '\t')
                length(header) >= 10 || throw(ArgumentError("VCF header must contain FORMAT and at least one sample"))
                sample_ids = String.(header[10:end])
                continue
            end
            isempty(strip(line)) && continue

            parts = split(chomp(line), '\t')
            length(parts) >= 10 || continue
            chr = normalise_chromosome(parts[1])
            pos = tryparse(Int, parts[2])
            pos === nothing && continue
            snp_id = parts[3] == "." ? "$(chr):$(pos)" : parts[3]
            ref = parts[4]
            alt = split(parts[5], ',')[1]
            format = length(parts) >= 9 ? parts[9] : ""

            dosage = Vector{Float64}(undef, length(parts) - 9)
            @inbounds for i in eachindex(dosage)
                dosage[i] = _vcf_gt_to_dosage(format, parts[9 + i])
            end

            if isempty(sample_ids)
                sample_ids = ["sample_$(i)" for i in 1:length(dosage)]
            elseif length(dosage) != length(sample_ids)
                throw(DimensionMismatch("VCF record sample count mismatch"))
            end

            push!(chromosomes, chr)
            push!(positions, pos)
            push!(snp_ids, snp_id)
            push!(allele1, ref)
            push!(allele2, alt)
            push!(dosage_cols, dosage)
        end
        return nothing
    end

    if endswith(lowercase(path), ".gz")
        codec = _load_optional_module(:CodecZlib)
        codec === nothing && throw(ArgumentError("reading .vcf.gz requires CodecZlib.jl"))
        open(path, "r") do raw
            stream = codec.GzipDecompressorStream(raw)
            try
                _parse_stream!(stream)
            finally
                close(stream)
            end
        end
    else
        open(path, "r") do io
            _parse_stream!(io)
        end
    end

    matrix = isempty(dosage_cols) ? zeros(Float64, length(sample_ids), 0) : hcat(dosage_cols...)
    bim = DataFrame(
        chromosome = chromosomes,
        snp_id = snp_ids,
        genetic_distance = zeros(Float64, length(snp_ids)),
        position = positions,
        allele1 = allele1,
        allele2 = allele2)
    fam = DataFrame(
        family_id = fill("0", length(sample_ids)),
        sample_id = sample_ids,
        paternal_id = fill("0", length(sample_ids)),
        maternal_id = fill("0", length(sample_ids)),
        sex = fill(0, length(sample_ids)),
        phenotype = fill(NaN, length(sample_ids)))
    prefix = replace(path, r"\.vcf(\.gz)?$"i => "")
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GenotypeMatrix(matrix, bim, fam; prefix=prefix), "read_vcf")
end

@inline function _dosage_to_gt(value::Float64; threshold::Real=0.1)
    isfinite(value) || return "./."
    hard = round(Int, clamp(value, 0.0, 2.0))
    abs(value - hard) <= threshold || return "./."
    hard == 0 && return "0/0"
    hard == 1 && return "0/1"
    return "1/1"
end

"""
    write_vcf(path, gm; threshold=0.1)

Write a `GenotypeMatrix` to a minimal VCF representation using hard-called GT fields.
"""
function write_vcf(path::String, gm::GenotypeMatrix; threshold::Real=0.1)
    G = gm.decoded
    sample_ids = String.(gm.fam.sample_id)
    open(path, "w") do io
        println(io, "##fileformat=VCFv4.2")
        println(io, join(vcat(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"], sample_ids), '\t'))
        for j in 1:size(G, 2)
            chr = normalise_chromosome(gm.bim.chromosome[j])
            pos = Int(gm.bim.position[j])
            snp_id = String(gm.bim.snp_id[j])
            ref = String(gm.bim.allele1[j])
            alt = String(gm.bim.allele2[j])
            gts = Vector{String}(undef, size(G, 1))
            @inbounds for i in 1:size(G, 1)
                gts[i] = _dosage_to_gt(Float64(G[i, j]); threshold=threshold)
            end
            println(io, join(vcat([chr, string(pos), snp_id, ref, alt, ".", "PASS", ".", "GT"], gts), '\t'))
        end
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, path, "write_vcf")
end

"""
    read_oxford_gen(path)

Read an Oxford GEN file into a `GenotypeMatrix` by converting genotype probabilities to dosages.
"""
function read_oxford_gen(path::String)
    chromosomes = String[]
    positions = Int[]
    snp_ids = String[]
    allele1 = String[]
    allele2 = String[]
    dosage_cols = Vector{Vector{Float64}}()
    n_samples = 0

    open(path, "r") do io
        for (line_no, line) in enumerate(eachline(io))
            stripped = strip(line)
            isempty(stripped) && continue
            parts = split(stripped)
            length(parts) >= 8 || throw(ArgumentError("invalid GEN line $(line_no): expected at least 8 fields"))

            chr = normalise_chromosome(parts[1])
            snp_id = isempty(parts[2]) ? parts[3] : parts[2]
            pos = parse(Int, parts[4])
            a1 = parts[5]
            a2 = parts[6]

            prob_tokens = parts[7:end]
            (length(prob_tokens) % 3 == 0) || throw(ArgumentError("invalid GEN line $(line_no): probability fields must be in triplets"))
            n_line_samples = div(length(prob_tokens), 3)
            if n_samples == 0
                n_samples = n_line_samples
            elseif n_line_samples != n_samples
                throw(DimensionMismatch("inconsistent sample count across GEN records"))
            end

            dosage = Vector{Float64}(undef, n_samples)
            idx = 1
            @inbounds for s in 1:n_samples
                p_aa = parse(Float64, prob_tokens[idx]); idx += 1
                p_ab = parse(Float64, prob_tokens[idx]); idx += 1
                p_bb = parse(Float64, prob_tokens[idx]); idx += 1
                total = p_aa + p_ab + p_bb
                if total <= eps(Float64)
                    dosage[s] = NaN
                else
                    dosage[s] = (p_ab + 2.0 * p_bb) / total
                end
            end

            push!(chromosomes, chr)
            push!(positions, pos)
            push!(snp_ids, snp_id)
            push!(allele1, a1)
            push!(allele2, a2)
            push!(dosage_cols, dosage)
        end
    end

    sample_ids = ["sample_$(i)" for i in 1:n_samples]
    matrix = isempty(dosage_cols) ? zeros(Float64, n_samples, 0) : hcat(dosage_cols...)
    bim = DataFrame(
        chromosome = chromosomes,
        snp_id = snp_ids,
        genetic_distance = zeros(Float64, length(snp_ids)),
        position = positions,
        allele1 = allele1,
        allele2 = allele2)
    fam = DataFrame(
        family_id = fill("0", n_samples),
        sample_id = sample_ids,
        paternal_id = fill("0", n_samples),
        maternal_id = fill("0", n_samples),
        sex = fill(0, n_samples),
        phenotype = fill(NaN, n_samples))
    prefix = replace(path, r"\.gen$"i => "")
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GenotypeMatrix(matrix, bim, fam; prefix=prefix), "read_oxford_gen")
end

"""
    read_plink(prefix::String)

Read a PLINK dataset given its file prefix.
"""
function read_plink(prefix::String)
    bim = _read_bim(String(prefix) * ".bim")
    fam = _read_fam(String(prefix) * ".fam")
    raw = open(String(prefix) * ".bed", "r") do io
        Mmap.mmap(io)
    end
    decoded = _decode_bed(raw, nrow(fam), nrow(bim))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GenotypeMatrix(raw, decoded, bim, fam, String(prefix)), "read_plink")
end

"""
    write_plink(prefix, genotypes, bim, fam)

Write a PLINK dataset to disk.
"""
function write_plink(prefix::String, genotypes::AbstractMatrix{<:Real}, bim::DataFrame, fam::DataFrame)
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
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, prefix, "write_plink")
end

"""
    calculate_maf(genotypes; min_maf=0.01, max_maf=1.0, return_frequency=false)

Compute per-variant minor allele frequencies and return a boolean pass mask.
"""
function calculate_maf(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; min_maf::Real=0.01, max_maf::Real=1.0, return_frequency::Bool=false)
    G = _matrix(genotypes)
    nobs, nsnps = size(G)
    maf = fill(NaN, nsnps)
    mask = trues(nsnps)

    if nsnps > 1 && Threads.nthreads() > 1
        @threads for snp in 1:nsnps
            alt_sum = 0.0
            n_valid = 0
            @inbounds for sample in 1:nobs
                value = G[sample, snp]
                if isfinite(value)
                    alt_sum += value
                    n_valid += 1
                end
            end
            if n_valid == 0
                mask[snp] = false
                continue
            end
            alt_freq = alt_sum / (2.0 * n_valid)
            alt_freq = clamp(alt_freq, 0.0, 1.0)
            maf[snp] = min(alt_freq, 1.0 - alt_freq)
            (maf[snp] >= min_maf && maf[snp] <= max_maf) || (mask[snp] = false)
        end
    else
        for snp in 1:nsnps
            alt_sum = 0.0
            n_valid = 0
            @inbounds for sample in 1:nobs
                value = G[sample, snp]
                if isfinite(value)
                    alt_sum += value
                    n_valid += 1
                end
            end
            if n_valid == 0
                mask[snp] = false
                continue
            end
            alt_freq = alt_sum / (2.0 * n_valid)
            alt_freq = clamp(alt_freq, 0.0, 1.0)
            maf[snp] = min(alt_freq, 1.0 - alt_freq)
            (maf[snp] >= min_maf && maf[snp] <= max_maf) || (mask[snp] = false)
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, return_frequency ? (mask, maf) : mask, "calculate_maf")
end

@inline function _qc_axis_dim(by::Symbol)
    by == :variant && return 1
    by == :sample && return 2
    throw(ArgumentError("`by` must be either :variant or :sample"))
end

"""
    calculate_missingness(genotypes; by=:variant, return_call_rate=false)

Compute missingness rates akin PLINK/ Bioconductor QC workflows.

- `by=:variant` returns one rate per SNP.
- `by=:sample` returns one rate per sample.
"""
function calculate_missingness(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; by::Symbol=:variant, return_call_rate::Bool=false)
    G = _matrix(genotypes)
    nobs, nsnps = size(G)
    missing = if by == :variant
        rates = zeros(Float64, nsnps)
        if nsnps > 1 && Threads.nthreads() > 1
            @threads for snp in 1:nsnps
                missing_count = 0
                @inbounds for sample in 1:nobs
                    missing_count += isfinite(G[sample, snp]) ? 0 : 1
                end
                rates[snp] = missing_count / max(nobs, 1)
            end
        else
            for snp in 1:nsnps
                missing_count = 0
                @inbounds for sample in 1:nobs
                    missing_count += isfinite(G[sample, snp]) ? 0 : 1
                end
                rates[snp] = missing_count / max(nobs, 1)
            end
        end
        rates
    elseif by == :sample
        rates = zeros(Float64, nobs)
        if nobs > 1 && Threads.nthreads() > 1
            @threads for sample in 1:nobs
                missing_count = 0
                @inbounds for snp in 1:nsnps
                    missing_count += isfinite(G[sample, snp]) ? 0 : 1
                end
                rates[sample] = missing_count / max(nsnps, 1)
            end
        else
            for sample in 1:nobs
                missing_count = 0
                @inbounds for snp in 1:nsnps
                    missing_count += isfinite(G[sample, snp]) ? 0 : 1
                end
                rates[sample] = missing_count / max(nsnps, 1)
            end
        end
        rates
    else
        _qc_axis_dim(by)
        throw(ArgumentError("`by` must be either :variant or :sample"))
    end
    call_rate = 1.0 .- missing
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, return_call_rate ? (missing, call_rate) : missing, "calculate_missingness")
end

"""
    missingness_filter(genotypes; by=:variant, max_missing=0.05, return_rates=false)

Return a pass mask based on maximum missingness threshold.
"""
function missingness_filter(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; by::Symbol=:variant, max_missing::Real=0.05, return_rates::Bool=false)
    rates = calculate_missingness(genotypes; by=by)
    mask = map(rate -> isfinite(rate) && rate <= max_missing, rates)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, return_rates ? (mask, rates) : mask, "missingness_filter")
end

"""
    info_score_filter(info_scores; min_info=0.8, max_info=1.0, return_scores=false)

Filter imputation INFO scores (common post-imputation QC step).
"""
function info_score_filter(info_scores::AbstractVector{<:Real}; min_info::Real=0.8, max_info::Real=1.0, return_scores::Bool=false)
    scores = Float64.(info_scores)
    mask = map(score -> isfinite(score) && score >= min_info && score <= max_info, scores)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, return_scores ? (mask, scores) : mask, "info_score_filter")
end

"""
    gwas_qc_report(genotypes; min_maf=0.01, max_missing=0.05, hwe_threshold=1e-6, info_scores=nothing, min_info=0.8)

Build a variant-level QC report inspired by PLINK/Bioconductor workflows.
"""
function gwas_qc_report(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; min_maf::Real=0.01, max_missing::Real=0.05, hwe_threshold::Real=1e-6, info_scores::Union{Nothing,AbstractVector{<:Real}}=nothing, min_info::Real=0.8)
    G = _matrix(genotypes)
    nsnps = size(G, 2)

    _, maf = calculate_maf(G; min_maf=min_maf, max_maf=1.0, return_frequency=true)
    hwe_p = calculate_hwe_pvalues(G)
    missing_rate, call_rate = calculate_missingness(G; by=:variant, return_call_rate=true)

    pass_maf = map(value -> isfinite(value) && value >= min_maf, maf)
    pass_missing = map(value -> isfinite(value) && value <= max_missing, missing_rate)
    pass_hwe = map(value -> isfinite(value) && value >= hwe_threshold, hwe_p)

    info_vec = fill(NaN, nsnps)
    pass_info = trues(nsnps)
    if info_scores !== nothing
        length(info_scores) == nsnps || throw(DimensionMismatch("info_scores length must match number of variants"))
        pass_info, info_vec = info_score_filter(info_scores; min_info=min_info, return_scores=true)
    end

    pass_all = pass_maf .& pass_missing .& pass_hwe .& pass_info

    if genotypes isa GenotypeMatrix
        chr = String.(genotypes.bim.chromosome)
        pos = Int.(genotypes.bim.position)
        ids = String.(genotypes.bim.snp_id)
    else
        chr = fill("", nsnps)
        pos = collect(1:nsnps)
        ids = ["snp_$(index)" for index in 1:nsnps]
    end

    df = DataFrame(
        CHR = chr,
        POS = pos,
        ID = ids,
        CALL_RATE = call_rate,
        MISSING_RATE = missing_rate,
        MAF = maf,
        HWE_P = hwe_p,
        INFO = info_vec,
        PASS_MAF = pass_maf,
        PASS_MISSING = pass_missing,
        PASS_HWE = pass_hwe,
        PASS_INFO = pass_info,
        PASS = pass_all)
    _ctx = active_provenance_context()
    return _register_gwas_result!(_ctx, df, "gwas_qc_report"; parents=genotypes isa GenotypeMatrix ? provenance_parent_ids(genotypes) : String[], parameters=(min_maf=Float64(min_maf), max_missing=Float64(max_missing), hwe_threshold=Float64(hwe_threshold), n_variants=nrow(df), n_pass=count(identity, pass_all)))
end

"""
    hwe_exact(obs_hom_ref, obs_het, obs_hom_alt, obs_unknown=0)

Exact Hardy-Weinberg equilibrium test for genotype counts.
"""
function hwe_exact(obs_hom_ref::Int, obs_het::Int, obs_hom_alt::Int, obs_unknown::Int=0)
    _ = obs_unknown
    n_called = obs_hom_ref + obs_het + obs_hom_alt
    n_called <= 0 && return NaN
    obs_hom_ref < 0 && throw(ArgumentError("HWE counts must be non-negative"))
    obs_het < 0 && throw(ArgumentError("HWE counts must be non-negative"))
    obs_hom_alt < 0 && throw(ArgumentError("HWE counts must be non-negative"))

    obs_homr = min(obs_hom_ref, obs_hom_alt)
    obs_homc = max(obs_hom_ref, obs_hom_alt)
    rare_copies = 2 * obs_homr + obs_het
    genotypes = n_called

    probs = zeros(Float64, rare_copies + 1)
    midpoint = Int(floor(rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes)))
    if isodd(rare_copies - midpoint)
        midpoint += 1
    end
    midpoint = clamp(midpoint, 0, rare_copies)

    probs[midpoint + 1] = 1.0
    normalizer = 1.0

    het = midpoint
    homr = Int((rare_copies - het) / 2)
    homc = genotypes - homr - het
    while het > 1
        probability = probs[het + 1] * het * (het - 1) / (4.0 * (homr + 1) * (homc + 1))
        probs[het - 1 + 1] = probability
        normalizer += probability
        het -= 2
        homr += 1
        homc += 1
    end

    het = midpoint
    homr = Int((rare_copies - het) / 2)
    homc = genotypes - homr - het
    while het <= rare_copies - 2
        probability = probs[het + 1] * 4.0 * homr * homc / ((het + 2) * (het + 1))
        probs[het + 2 + 1] = probability
        normalizer += probability
        het += 2
        homr -= 1
        homc -= 1
    end

    probs ./= max(normalizer, eps(Float64))
    obs_probability = obs_het <= rare_copies ? probs[obs_het + 1] : 0.0
    pvalue = sum(probability for probability in probs if probability <= obs_probability + 1e-12)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, clamp(pvalue, 0.0, 1.0), "hwe_exact")
end

"""
    calculate_hwe_pvalues(genotypes)

Compute exact HWE p-values for each variant.
"""
function calculate_hwe_pvalues(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}})
    G = _matrix(genotypes)
    nsnps = size(G, 2)
    pvalues = fill(NaN, nsnps)

    for snp in 1:nsnps
        col = @view G[:, snp]
        valid = isfinite.(col)
        n_valid = sum(valid)
        if n_valid == 0
            continue
        end
        rounded = round.(Int, clamp.(col[valid], 0.0, 2.0))
        obs_hom_ref = count(==(0), rounded)
        obs_het = count(==(1), rounded)
        obs_hom_alt = count(==(2), rounded)
        obs_unknown = size(G, 1) - n_valid
        pvalues[snp] = hwe_exact(obs_hom_ref, obs_het, obs_hom_alt, obs_unknown)
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, pvalues, "calculate_hwe_pvalues")
end

"""
    hwe_filter(genotypes; p_threshold=1e-6, return_pvalues=false)

Return a pass mask of variants satisfying exact HWE p-value threshold.
"""
function hwe_filter(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; p_threshold::Real=1e-6, return_pvalues::Bool=false)
    pvalues = calculate_hwe_pvalues(genotypes)
    mask = map(value -> isfinite(value) && value >= p_threshold, pvalues)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, return_pvalues ? (mask, pvalues) : mask, "hwe_filter")
end

"""
    genomic_control_lambda(result)

Compute genomic-control lambda from GWAS p-values.
"""
function genomic_control_lambda(result::GWASResult)
    valid = [clamp(Float64(value), eps(Float64), 1.0) for value in result.pvalue if isfinite(value)]
    isempty(valid) && return 1.0
    chi2 = quantile.(Ref(Chisq(1)), 1.0 .- valid)
    median_chi2 = median(chi2)
    expected = quantile(Chisq(1), 0.5)
    expected <= 0 && return 1.0
    lambda = median_chi2 / expected
    λ = isfinite(lambda) && lambda > 0 ? lambda : 1.0
    _ctx = active_provenance_context()
    return _register_gwas_result!(_ctx, λ, "genomic_control_lambda"; parents=provenance_parent_ids(result), parameters=(lambda=λ, n_variants=length(result.pvalue)))
end

"""
    apply_genomic_control(result; lambda=genomic_control_lambda(result))

Return a lambda-adjusted GWAS result table.
"""
function apply_genomic_control(result::GWASResult; lambda::Real=genomic_control_lambda(result))
    scale = sqrt(max(Float64(lambda), 1.0))
    adjusted_z = result.zscore ./ scale
    adjusted_p = 2 .* ccdf.(Ref(Normal()), abs.(adjusted_z))
    adjusted_se = result.standard_error .* scale
    return GWASResult(
        result.snp_ids,
        result.chromosomes,
        result.positions,
        result.alleles,
        result.gene_ids,
        result.beta,
        adjusted_se,
        adjusted_z,
        adjusted_p,
        result.sample_size,
        result.covariate_names,
        result.phenotype_name,
        string(result.method, "+gc"))
end

"""
    calculate_ld_matrix(genotypes; min_maf=0.01, max_maf=1.0, multi_thread=true, return_snp_indices=false, max_snps=50_000, window_kb=nothing, sparse=false, verbose=false)

Compute SNP-by-SNP LD correlation after MAF filtering.

For very large variant counts, set `window_kb` for band-limited LD or `sparse=true`
to avoid dense matrix materialization.
"""
function calculate_ld_matrix(
    genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}};
    min_maf::Real=0.01,
    max_maf::Real=1.0,
    multi_thread::Bool=true,
    return_snp_indices::Bool=false,
    max_snps::Int=50_000,
    window_kb::Union{Nothing,Int}=nothing,
    sparse::Bool=false,
    verbose::Bool=false)
    G = _matrix(genotypes)
    mask, maf = calculate_maf(G; min_maf=min_maf, max_maf=max_maf, return_frequency=true)
    snp_indices = findall(mask)
    isempty(snp_indices) && return return_snp_indices ? (zeros(Float64, 0, 0), Int[], Float64[]) : zeros(Float64, 0, 0)

    filtered = Matrix{Float64}(G[:, snp_indices])
    nobs, nsnps = size(filtered)

    if window_kb === nothing && nsnps > max_snps
        throw(ArgumentError("requested dense LD matrix for $(nsnps) SNPs exceeds max_snps=$(max_snps); set window_kb or increase max_snps explicitly"))
    end

    @inline function _normalize_column!(col::AbstractVector{<:AbstractFloat}, nrows::Int)
        finite_sum = 0.0
        n_valid = 0
        @inbounds for row in eachindex(col)
            value = col[row]
            if isfinite(value)
                finite_sum += value
                n_valid += 1
            end
        end
        if n_valid == 0
            fill!(col, 0.0)
            return nothing
        end

        mean_value = finite_sum / n_valid
        ss = 0.0
        @inbounds for row in eachindex(col)
            value = col[row]
            centered = isfinite(value) ? (value - mean_value) : 0.0
            col[row] = centered
            ss += centered * centered
        end

        std_value = sqrt(ss / max(nrows - 1, 1))
        if !(isfinite(std_value) && std_value > eps(Float64))
            fill!(col, 0.0)
        else
            col ./= std_value
        end
        return nothing
    end

    if _thread_enabled(nsnps, multi_thread)
        @threads for snp in 1:nsnps
            col = @view filtered[:, snp]
            _normalize_column!(col, nobs)
        end
    else
        for snp in 1:nsnps
            col = @view filtered[:, snp]
            _normalize_column!(col, nobs)
            _progress_tick!(verbose, "calculate_ld_matrix:normalise", snp, nsnps; step=25_000)
        end
    end

    denom = max(size(filtered, 1) - 1, 1)
    if window_kb === nothing
        ld = Symmetric((filtered' * filtered) / denom)
        ld_matrix = sparse ? SparseArrays.sparse(Matrix(ld)) : Matrix(ld)
        for index in 1:nsnps
            ld_matrix[index, index] = 1.0
        end
        return return_snp_indices ? (ld_matrix, snp_indices, maf[snp_indices]) : ld_matrix
    end

    window_bp = max(Int(window_kb), 0) * 1000
    chrom = genotypes isa GenotypeMatrix ? String.(genotypes.bim.chromosome[snp_indices]) : fill("", nsnps)
    pos = genotypes isa GenotypeMatrix ? Int.(genotypes.bim.position[snp_indices]) : collect(1:nsnps)

    if sparse
        rows = Int[]
        cols = Int[]
        vals = Float64[]
        sizehint!(rows, 4 * nsnps)
        sizehint!(cols, 4 * nsnps)
        sizehint!(vals, 4 * nsnps)

        for i in 1:nsnps
            push!(rows, i); push!(cols, i); push!(vals, 1.0)
            xi = @view filtered[:, i]
            for j in (i + 1):nsnps
                chrom[i] == chrom[j] || continue
                abs(pos[j] - pos[i]) <= window_bp || continue
                xj = @view filtered[:, j]
                r = dot(xi, xj) / denom
                push!(rows, i); push!(cols, j); push!(vals, r)
                push!(rows, j); push!(cols, i); push!(vals, r)
            end
            _progress_tick!(verbose, "calculate_ld_matrix:windowed", i, nsnps; step=5_000)
        end

        ld_sparse = SparseArrays.sparse(rows, cols, vals, nsnps, nsnps)
        return return_snp_indices ? (ld_sparse, snp_indices, maf[snp_indices]) : ld_sparse
    end

    ld_dense = Matrix{Float64}(I, nsnps, nsnps)
    for i in 1:nsnps
        xi = @view filtered[:, i]
        for j in (i + 1):nsnps
            chrom[i] == chrom[j] || continue
            abs(pos[j] - pos[i]) <= window_bp || continue
            xj = @view filtered[:, j]
            r = dot(xi, xj) / denom
            ld_dense[i, j] = r
            ld_dense[j, i] = r
        end
        _progress_tick!(verbose, "calculate_ld_matrix:windowed", i, nsnps; step=5_000)
    end

    return return_snp_indices ? (ld_dense, snp_indices, maf[snp_indices]) : ld_dense
end

"""
    calculate_grm(genotypes; block_size=2048, center=false, verbose=false)

Construct a PLINK-style uncentered GRM normalized by its mean diagonal.
"""
function calculate_grm(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; block_size::Int=2048, center::Bool=false, verbose::Bool=false)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, _kinship_from_genotypes(_matrix(genotypes); block_size=block_size, center=center, verbose=verbose), "calculate_grm")
end

"""
    loco_projection(covariates, values; ridge=1e-6, return_projector=false)

Project values onto the orthogonal complement of covariates (LoCO-style projection).
"""
function loco_projection(covariates::AbstractMatrix{<:Real}, values::AbstractMatrix{<:Real}; ridge::Real=1e-6, return_projector::Bool=false)
    X = Matrix{Float64}(covariates)
    V = Matrix{Float64}(values)
    size(X, 1) == size(V, 1) || throw(DimensionMismatch("covariates and values must share sample dimension"))
    XtX = Symmetric(X' * X + Float64(ridge) * I)
    projected = V - X * (XtX \ (X' * V))
    if return_projector
        projector = Matrix{Float64}(I, size(X, 1), size(X, 1)) - X * (XtX \ X')
        return projected, projector
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, projected, "loco_projection")
end

function loco_projection(covariates::AbstractMatrix{<:Real}, values::AbstractVector{<:Real}; ridge::Real=1e-6, return_projector::Bool=false)
    projected = loco_projection(covariates, reshape(Vector{Float64}(values), :, 1); ridge=ridge, return_projector=return_projector)
    if return_projector
        projected_values, projector = projected
        return vec(projected_values), projector
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, vec(projected), "loco_projection")
end

function _covariate_matrix(covariates::Union{Nothing,AbstractMatrix{<:Real}}, nobs::Int)
    if covariates === nothing
        return ones(Float64, nobs, 1), _pooled_strings(["intercept"])
    end
    matrix = Matrix{Float64}(covariates)
    size(matrix, 1) == nobs || throw(ArgumentError("covariates must have one row per sample"))

    # Avoid duplicating the intercept if callers already provided an all-ones column.
    intercept_col = findfirst(col -> all(abs.(col .- 1.0) .<= 1e-8), eachcol(matrix))
    if intercept_col === nothing
        names = ["intercept"; ["cov_$(index)" for index in 1:size(matrix, 2)]...]
        return hcat(ones(Float64, nobs, 1), matrix), _pooled_strings(names)
    end

    ordered = vcat([intercept_col], [idx for idx in 1:size(matrix, 2) if idx != intercept_col])
    matrix_ordered = matrix[:, ordered]
    names = ["intercept"; ["cov_$(index)" for index in 1:(size(matrix_ordered, 2) - 1)]...]
    return matrix_ordered, _pooled_strings(names)
end

function _project_out(design::AbstractMatrix{<:Real}, values::AbstractMatrix{<:Real})
    X = Matrix{Float64}(design)
    V = Matrix{Float64}(values)
    design_qr = qr(X)
    coeff = design_qr \ V
    residual = copy(V)
    mul!(residual, X, coeff, -1.0, 1.0)
    return residual, design_qr
end

function _project_out(design::AbstractMatrix{<:Real}, values::AbstractVector{<:Real})
    X = Matrix{Float64}(design)
    v = Vector{Float64}(values)
    design_qr = qr(X)
    coeff = design_qr \ v
    residual = copy(v)
    mul!(residual, X, coeff, -1.0, 1.0)
    return residual, design_qr
end

function _result_from_matrix(genotypes::GenotypeMatrix, beta::Vector{Float64}, se::Vector{Float64}, z::Vector{Float64}, p::Vector{Float64}; method::String, phenotype_name::String, covariate_names)
    gene_ids = hasproperty(genotypes.bim, :gene_id) ? String.(genotypes.bim.gene_id) : fill("", length(beta))
    alleles = collect(zip(String.(genotypes.bim.allele1), String.(genotypes.bim.allele2)))
    return GWASResult(_pooled_strings(genotypes.bim.snp_id), _pooled_strings(genotypes.bim.chromosome), Int.(genotypes.bim.position), alleles, _pooled_strings(gene_ids), beta, se, z, p, size(genotypes, 1), _pooled_strings(covariate_names), String(phenotype_name), String(method))
end

function _result_from_matrix(genotypes::AbstractMatrix{<:Real}, beta::Vector{Float64}, se::Vector{Float64}, z::Vector{Float64}, p::Vector{Float64}; method::String, phenotype_name::String, covariate_names)
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
        method = PooledArray(fill(result.method, length(result.snp_ids))))
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
        method = PooledArray(fill(result.method, length(result.snp_ids))))
end

"""
    to_plink_dataframe(result, genotypes=nothing)

Format GWAS summary statistics using PLINK-style schema:
`CHR POS ID REF ALT BETA SE P`.
"""
function to_plink_dataframe(result::GWASResult, genotypes::Union{Nothing,GenotypeMatrix}=nothing)
    allele_lookup = Dict{String,Tuple{String,String}}()
    if genotypes !== nothing
        for index in eachindex(genotypes.bim.snp_id)
            allele_lookup[String(genotypes.bim.snp_id[index])] = (String(genotypes.bim.allele1[index]), String(genotypes.bim.allele2[index]))
        end
    end

    refs = Vector{String}(undef, length(result.snp_ids))
    alts = Vector{String}(undef, length(result.snp_ids))
    for index in eachindex(result.snp_ids)
        snp_id = String(result.snp_ids[index])
        allele_pair = haskey(allele_lookup, snp_id) ? allele_lookup[snp_id] : result.alleles[index]
        refs[index] = allele_pair[1]
        alts[index] = allele_pair[2]
    end
    return DataFrame(
        CHR = String.(result.chromosomes),
        POS = Int.(result.positions),
        ID = String.(result.snp_ids),
        REF = refs,
        ALT = alts,
        BETA = Float64.(result.beta),
        SE = Float64.(result.standard_error),
        P = Float64.(result.pvalue))
end

"""
    write_plink_sumstats(path, result; genotypes=nothing)

Write PLINK-style summary statistics file with columns:
`CHR POS ID REF ALT BETA SE P`.
"""
function write_plink_sumstats(path::String, result::GWASResult; genotypes::Union{Nothing,GenotypeMatrix}=nothing)
    table = to_plink_dataframe(result, genotypes)
    open(path, "w") do io
        buffer = IOBuffer()
        println(buffer, "CHR\tPOS\tID\tREF\tALT\tBETA\tSE\tP")
        for row in eachrow(table)
            println(buffer, join((row.CHR, row.POS, row.ID, row.REF, row.ALT, row.BETA, row.SE, row.P), '\t'))
            if position(buffer) >= 1_000_000
                write(io, take!(buffer))
            end
        end
        write(io, take!(buffer))
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, path, "write_plink_sumstats")
end

function _gwas_dot_products_cuda(G::Matrix{Float64}, y::Vector{Float64})
    if !isdefined(@__MODULE__, :CUDA)
        try
            @eval import CUDA
        catch
            return nothing
        end
    end
    CUDA.functional() || return nothing
    G_cu = CUDA.CuArray(G)
    y_cu = CUDA.CuArray(y)
    dotxy = vec(Array(transpose(G_cu) * y_cu))
    xxt = vec(Array(sum(abs2, G_cu; dims=1)))
    return dotxy, xxt
end

"""
    gwas_linear_scan(genotypes, phenotype; covariates=nothing, phenotype_name="phenotype")

Run a linear-regression GWAS scan across all variants.
"""
function gwas_linear_scan(
    genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}},
    phenotype::AbstractVector{<:Real};
    covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
    phenotype_name::String="phenotype",
    multi_thread::Bool=true,
    use_cuda::Bool=false,
    dof_override::Union{Nothing,Integer}=nothing,
    precision::DataType=Float64,
    verbose::Bool=false)
    base_G = _matrix(genotypes)
    _validate_gwas_inputs(base_G; phenotype=phenotype, covariates=covariates)
    T = precision == Float32 ? Float32 : Float64
    G = Matrix{T}(base_G)
    nobs, nsnps = size(G)

    design, covariate_names = _covariate_matrix(covariates, nobs)
    X = Matrix{Float64}(design)
    design_qr = qr(X)

    y_res = Vector{Float64}(phenotype)
    coeff_y = design_qr \ y_res
    mul!(y_res, X, coeff_y, -1.0, 1.0)

    dof = dof_override === nothing ? max(nobs - size(design, 2) - 1, 1) : max(Int(dof_override), 1)
    beta = zeros(Float64, nsnps)
    se = fill(Inf, nsnps)
    z = zeros(Float64, nsnps)
    p = ones(Float64, nsnps)

    use_cuda && @warn("CUDA path is disabled for streaming gwas_linear_scan; running on CPU")

    y_norm2 = dot(y_res, y_res)
    t_dist = TDist(dof)
    thread_buffers = [Vector{Float64}(undef, nobs) for _ in 1:Threads.nthreads()]

    @inline function _scan_one_snp!(snp::Int)
        g = thread_buffers[Threads.threadid()]
        @inbounds for i in 1:nobs
            g[i] = Float64(G[i, snp])
        end

        finite_sum = 0.0
        n_valid = 0
        @inbounds for i in 1:nobs
            value = g[i]
            if isfinite(value)
                finite_sum += value
                n_valid += 1
            end
        end
        n_valid == 0 && return nothing

        mean_value = finite_sum / n_valid
        @inbounds for i in 1:nobs
            isfinite(g[i]) || (g[i] = mean_value)
        end

        coeff_g = design_qr \ g
        mul!(g, X, coeff_g, -1.0, 1.0)
        denom = dot(g, g)
        (!isfinite(denom) || denom <= eps(Float64)) && return nothing

        dotxy = dot(g, y_res)
        b = dotxy / denom
        rss = y_norm2 - 2 * b * dotxy + b * b * denom
        sigma2 = max(rss / dof, eps(Float64))
        stderr = sqrt(sigma2 / denom)
        beta[snp] = b
        se[snp] = stderr
        z[snp] = b / stderr
        p[snp] = 2 * ccdf(t_dist, abs(z[snp]))
        return nothing
    end

    if _thread_enabled(nsnps, multi_thread)
        @threads for snp in 1:nsnps
            _scan_one_snp!(snp)
        end
    else
        for snp in 1:nsnps
            _scan_one_snp!(snp)
            _progress_tick!(verbose, "gwas_linear_scan", snp, nsnps; step=25_000)
        end
    end

    return genotypes isa GenotypeMatrix ?
           _result_from_matrix(genotypes, beta, se, z, p; method="linear_scan", phenotype_name=phenotype_name, covariate_names=covariate_names) :
           _result_from_matrix(base_G, beta, se, z, p; method="linear_scan", phenotype_name=phenotype_name, covariate_names=covariate_names)
end

function _binary_response(phenotype::AbstractVector{<:Real})
    y = Float64.(phenotype)
    valid = [value for value in y if isfinite(value)]
    isempty(valid) && throw(ArgumentError("binary phenotype has no finite values"))
    values = sort!(unique(valid))
    if all(value -> value in (0.0, 1.0), values)
        return y
    elseif length(values) == 2
        low, high = values[1], values[2]
        return [!isfinite(value) ? NaN : (value == high ? 1.0 : 0.0) for value in y]
    end
    throw(ArgumentError("binary trait GWAS requires a phenotype with two unique finite values"))
end

function _logistic_fit_last_coefficient(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}; firth::Bool=true, max_iter::Int=50, tol::Real=1e-8, ridge::Real=1e-6)
    Xf = Matrix{Float64}(X)
    yf = Vector{Float64}(y)
    nobs, p = size(Xf)
    nobs > p || return 0.0, Inf, 0.0, 1.0

    beta = zeros(Float64, p)
    converged = false

    for _ in 1:max_iter
        eta = clamp.(Xf * beta, -35.0, 35.0)
        mu = 1.0 ./ (1.0 .+ exp.(-eta))
        W = clamp.(mu .* (1.0 .- mu), eps(Float64), Inf)
        XtWX = Xf' * (Xf .* W)
        hessian = Symmetric(XtWX + Float64(ridge) * I)
        factor = try
            cholesky(hessian)
        catch
            return 0.0, Inf, 0.0, 1.0
        end
        score = Xf' * (yf .- mu)

        if firth
            hinv_xt = factor \ transpose(Xf)
            leverage = W .* vec(sum(Xf .* transpose(hinv_xt), dims=2))
            score .+= Xf' * ((0.5 .- mu) .* leverage)
        end

        step = factor \ score
        beta_new = beta + step
        if maximum(abs.(step)) <= tol * (1.0 + maximum(abs.(beta_new)))
            beta = beta_new
            converged = true
            break
        end
        beta = beta_new
    end

    if !converged && any(!isfinite, beta)
        return 0.0, Inf, 0.0, 1.0
    end

    eta = clamp.(Xf * beta, -35.0, 35.0)
    mu = 1.0 ./ (1.0 .+ exp.(-eta))
    W = clamp.(mu .* (1.0 .- mu), eps(Float64), Inf)
    XtWX = Xf' * (Xf .* W)
    hessian = Symmetric(XtWX + Float64(ridge) * I)
    factor = try
        cholesky(hessian)
    catch
        return 0.0, Inf, 0.0, 1.0
    end
    e_last = zeros(Float64, p)
    e_last[end] = 1.0
    v_last = factor \ e_last
    stderr = sqrt(max(v_last[end], eps(Float64)))
    zscore = beta[end] / stderr
    pvalue = 2.0 * ccdf(Normal(), abs(zscore))
    return beta[end], stderr, zscore, pvalue
end

function _linear_fit_last_coefficient(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}; ridge::Real=1e-8)
    Xf = Matrix{Float64}(X)
    yf = Vector{Float64}(y)
    nobs, p = size(Xf)
    nobs > p || return 0.0, Inf, 0.0, 1.0
    xtx = Symmetric(Xf' * Xf + Float64(ridge) * I)
    beta = try
        xtx \ (Xf' * yf)
    catch
        return 0.0, Inf, 0.0, 1.0
    end
    residual = yf - Xf * beta
    dof = max(nobs - p, 1)
    sigma2 = max(sum(abs2, residual) / dof, eps(Float64))
    xtx_inv = try
        inv(xtx)
    catch
        return 0.0, Inf, 0.0, 1.0
    end
    stderr = sqrt(max(sigma2 * xtx_inv[end, end], eps(Float64)))
    zscore = beta[end] / stderr
    pvalue = 2.0 * ccdf(TDist(dof), abs(zscore))
    return beta[end], stderr, zscore, pvalue
end

"""
    gwas_logistic_scan(genotypes, phenotype; covariates=nothing, phenotype_name="phenotype", firth=true)

Run a binary-trait GWAS logistic scan with optional Firth bias correction.
"""
function gwas_logistic_scan(
    genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}},
    phenotype::AbstractVector{<:Real};
    covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
    phenotype_name::String="phenotype",
    firth::Bool=true,
    multi_thread::Bool=true,
    max_iter::Int=50,
    tol::Real=1e-8,
    ridge::Real=1e-6,
    precision::DataType=Float64,
    verbose::Bool=false)
    base_G = _matrix(genotypes)
    _validate_gwas_inputs(base_G; phenotype=phenotype, covariates=covariates)
    T = precision == Float32 ? Float32 : Float64
    G = Matrix{T}(base_G)
    nobs, nsnps = size(G)

    y = _binary_response(phenotype)
    y_finite = isfinite.(y)
    design, covariate_names = _covariate_matrix(covariates, nobs)
    covariate_names_vec = String.(covariate_names)
    n_cov = size(design, 2)
    design_valid = vec(all(isfinite.(design), dims=2))

    beta = zeros(Float64, nsnps)
    se = fill(Inf, nsnps)
    z = zeros(Float64, nsnps)
    p = ones(Float64, nsnps)

    @inline function _scan_one_snp!(snp::Int)
        g = @view G[:, snp]
        valid_buf = falses(nobs)
        n_valid = 0
        @inbounds for i in 1:nobs
            keep = design_valid[i] && y_finite[i] && isfinite(g[i])
            valid_buf[i] = keep
            n_valid += keep
        end
        n_valid > n_cov || return nothing

        X = Matrix{Float64}(undef, n_valid, n_cov + 1)
        y_local = Vector{Float64}(undef, n_valid)
        row = 0
        @inbounds for i in 1:nobs
            valid_buf[i] || continue
            row += 1
            X[row, 1:n_cov] = design[i, :]
            X[row, end] = Float64(g[i])
            y_local[row] = y[i]
        end

        coeff, stderr, zscore, pvalue = _logistic_fit_last_coefficient(X, y_local; firth=firth, max_iter=max_iter, tol=tol, ridge=ridge)
        beta[snp] = coeff
        se[snp] = stderr
        z[snp] = zscore
        p[snp] = pvalue
        return nothing
    end

    if _thread_enabled(nsnps, multi_thread)
        @threads for snp in 1:nsnps
            _scan_one_snp!(snp)
        end
    else
        for snp in 1:nsnps
            _scan_one_snp!(snp)
            _progress_tick!(verbose, "gwas_logistic_scan", snp, nsnps; step=10_000)
        end
    end

    method = firth ? "logistic_scan_firth" : "logistic_scan"
    return genotypes isa GenotypeMatrix ?
           _result_from_matrix(genotypes, beta, se, z, p; method=method, phenotype_name=phenotype_name, covariate_names=covariate_names_vec) :
           _result_from_matrix(base_G, beta, se, z, p; method=method, phenotype_name=phenotype_name, covariate_names=covariate_names_vec)
end

"""
    gwas_gxe_interaction(genotypes, phenotype, snp_of_interest, snps_of_effect_modifier; ...)

Test SNP-by-SNP interaction effects using either linear or logistic models.
"""
function gwas_gxe_interaction(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, phenotype::AbstractVector{<:Real}, snp_of_interest::Integer, snps_of_effect_modifier::Union{Integer,AbstractVector{<:Integer}}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, phenotype_name::String="phenotype", logistic::Bool=false, firth::Bool=true, multi_thread::Bool=true, max_iter::Int=50, tol::Real=1e-8, ridge::Real=1e-6)
    G = _matrix(genotypes)
    nobs, nsnps = size(G)
    length(phenotype) == nobs || throw(ArgumentError("phenotype length must match sample count"))
    1 <= snp_of_interest <= nsnps || throw(BoundsError("snp_of_interest is out of bounds"))

    modifiers = snps_of_effect_modifier isa Integer ? [Int(snps_of_effect_modifier)] : unique(Int.(collect(snps_of_effect_modifier)))
    all(index -> 1 <= index <= nsnps, modifiers) || throw(BoundsError("modifier SNP index out of bounds"))

    design, covariate_names = _covariate_matrix(covariates, nobs)
    covariate_names_vec = String.(covariate_names)
    design_valid = vec(all(isfinite.(design), dims=2))
    y = logistic ? _binary_response(phenotype) : Float64.(phenotype)

    m = length(modifiers)
    beta = zeros(Float64, m)
    se = fill(Inf, m)
    z = zeros(Float64, m)
    p = ones(Float64, m)

    focal = Vector{Float64}(G[:, snp_of_interest])

    if multi_thread && m > 1 && Threads.nthreads() > 1
        @threads for idx in 1:m
            modifier_index = modifiers[idx]
            modifier = Vector{Float64}(G[:, modifier_index])
            interaction = focal .* modifier
            valid = design_valid .& isfinite.(focal) .& isfinite.(modifier) .& isfinite.(interaction) .& isfinite.(y)
            sum(valid) > size(design, 2) + 3 || continue
            X = hcat(design[valid, :], focal[valid], modifier[valid], interaction[valid])
            coeff, stderr, zscore, pvalue = if logistic
                _logistic_fit_last_coefficient(X, y[valid]; firth=firth, max_iter=max_iter, tol=tol, ridge=ridge)
            else
                _linear_fit_last_coefficient(X, y[valid]; ridge=ridge)
            end
            beta[idx] = coeff
            se[idx] = stderr
            z[idx] = zscore
            p[idx] = pvalue
        end
    else
        for idx in 1:m
            modifier_index = modifiers[idx]
            modifier = Vector{Float64}(G[:, modifier_index])
            interaction = focal .* modifier
            valid = design_valid .& isfinite.(focal) .& isfinite.(modifier) .& isfinite.(interaction) .& isfinite.(y)
            sum(valid) > size(design, 2) + 3 || continue
            X = hcat(design[valid, :], focal[valid], modifier[valid], interaction[valid])
            coeff, stderr, zscore, pvalue = if logistic
                _logistic_fit_last_coefficient(X, y[valid]; firth=firth, max_iter=max_iter, tol=tol, ridge=ridge)
            else
                _linear_fit_last_coefficient(X, y[valid]; ridge=ridge)
            end
            beta[idx] = coeff
            se[idx] = stderr
            z[idx] = zscore
            p[idx] = pvalue
        end
    end

    if genotypes isa GenotypeMatrix
        snp_ids = [string(genotypes.bim.snp_id[snp_of_interest], "_x_", genotypes.bim.snp_id[index]) for index in modifiers]
        chromosomes = String.(genotypes.bim.chromosome[modifiers])
        positions = Int.(genotypes.bim.position[modifiers])
        alleles = [(String(genotypes.bim.allele1[index]), String(genotypes.bim.allele2[index])) for index in modifiers]
        gene_ids = hasproperty(genotypes.bim, :gene_id) ? String.(genotypes.bim.gene_id[modifiers]) : fill("", m)
    else
        snp_ids = ["snp_$(snp_of_interest)_x_snp_$(index)" for index in modifiers]
        chromosomes = fill("", m)
        positions = modifiers
        alleles = fill(("", ""), m)
        gene_ids = fill("", m)
    end

    model_name = logistic ? (firth ? "gxe_logistic_firth" : "gxe_logistic") : "gxe_linear"
    return GWASResult(
        snp_ids,
        chromosomes,
        positions,
        alleles,
        gene_ids,
        beta,
        se,
        z,
        p,
        nobs,
        vcat(covariate_names_vec, ["focal_snp", "modifier_snp", "interaction_term"]),
        phenotype_name,
        model_name)
end

function _impute_missing_by_mean!(matrix::AbstractMatrix{<:Real})
    nobs, ncols = size(matrix)
    for column in 1:ncols
        col = @view matrix[:, column]
        finite_sum = 0.0
        n_valid = 0
        @inbounds for row in 1:nobs
            value = col[row]
            if isfinite(value)
                finite_sum += value
                n_valid += 1
            end
        end
        if n_valid == 0
            col .= 0.0
            continue
        end
        mean_value = finite_sum / n_valid
        @inbounds for row in 1:nobs
            isfinite(col[row]) || (col[row] = mean_value)
        end
    end
    return matrix
end

function _kinship_from_genotypes(genotypes::AbstractMatrix{<:Real}; block_size::Int=2048, center::Bool=false, verbose::Bool=false)
    X = Matrix{Float64}(genotypes)
    _impute_missing_by_mean!(X)
    n, p = size(X)
    p > 0 || return Matrix{Float64}(I, n, n)

    K = zeros(Float64, n, n)
    b = max(1, Int(block_size))
    for start_col in 1:b:p
        end_col = min(start_col + b - 1, p)
        block = @view X[:, start_col:end_col]
        if center
            centered = block .- mean(block, dims=1)
            K .+= centered * centered'
        else
            K .+= block * block'
        end
        _progress_tick!(verbose, "calculate_grm", end_col, p; step=max(10_000, 10 * b))
    end

    K ./= p
    if any(value -> !isfinite(value), K)
        return Matrix{Float64}(I, n, n)
    end

    denominator = mean(diag(K))
    if !isfinite(denominator) || denominator <= eps(Float64)
        return Matrix{Float64}(I, n, n)
    end
    return K ./ denominator
end

# he_regression_variance_components is defined below with full _ctx support.

function _estimate_variance_components(y::AbstractVector{<:Real}, K::AbstractMatrix{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing)
    yv = Vector{Float64}(y)
    n = length(yv)
    n > 1 || return (sigma_g=0.0, sigma_e=1.0, h2=0.0, dof=1)

    X, _ = _covariate_matrix(covariates, n)
    size(X, 1) == n || throw(DimensionMismatch("covariates must match phenotype length"))

    reml_fit = try
        reml_variance_components(yv, [Matrix{Float64}(K)]; covariates=X, add_intercept=false, max_iter=60, tol=1e-6)
    catch
        he_regression_variance_components(yv, [Matrix{Float64}(K)])
    end

    sigma_g = max(Float64(reml_fit.sigma2_random[1]), 0.0)
    sigma_e = max(Float64(reml_fit.sigma2_residual), eps(Float64))
    total = sigma_g + sigma_e
    h2 = total <= eps(Float64) ? 0.0 : sigma_g / total
    dof = max(n - LinearAlgebra.rank(X) - 1, 1)
    return (sigma_g=sigma_g, sigma_e=sigma_e, h2=h2, dof=dof)
end

"""
    gwas_lmm_scan(genotypes, phenotype; covariates=nothing, kinship=nothing, phenotype_name="phenotype")

Run a mixed-model GWAS scan using an optional kinship matrix.
"""
function gwas_lmm_scan(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, phenotype::AbstractVector{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, kinship::Union{Nothing,AbstractMatrix{<:Real}}=nothing, phenotype_name::String="phenotype", multi_thread::Bool=true, use_cuda::Bool=false, precision::DataType=Float64, verbose::Bool=false)
    base_G = _matrix(genotypes)
    _validate_gwas_inputs(base_G; phenotype=phenotype, covariates=covariates, kinship=kinship)
    nobs = size(base_G, 1)

    T = precision == Float32 ? Float32 : Float64
    G = Matrix{T}(base_G)
    K = kinship === nothing ? _kinship_from_genotypes(base_G; verbose=verbose) : Matrix{Float64}(kinship)
    if any(value -> !isfinite(value), K)
        K = Matrix{Float64}(I, nobs, nobs)
    end

    design = covariates === nothing ? ones(Float64, nobs, 1) : hcat(ones(Float64, nobs, 1), Matrix{Float64}(covariates))
    vc = _estimate_variance_components(phenotype, K; covariates=design)

    symK = Symmetric((K + K') / 2)
    decomp = eigen(symK)
    lambdas = vc.sigma_g .* decomp.values .+ vc.sigma_e
    weights = 1.0 ./ sqrt.(max.(lambdas, eps(Float64)))

    Ut = decomp.vectors'
    G_transformed = Ut * Matrix{Float64}(G)
    G_transformed .*= reshape(weights, :, 1)

    y_transformed = Ut * Vector{Float64}(phenotype)
    y_transformed .*= weights

    covariates_transformed = if covariates === nothing
        nothing
    else
        C = Ut * Matrix{Float64}(covariates)
        C .*= reshape(weights, :, 1)
        C
    end

    scan = gwas_linear_scan(
        G_transformed,
        y_transformed;
        covariates=covariates_transformed,
        phenotype_name=phenotype_name,
        multi_thread=multi_thread,
        use_cuda=use_cuda,
        dof_override=vc.dof,
        precision=precision,
        verbose=verbose)
    return GWASResult(
        scan.snp_ids,
        scan.chromosomes,
        scan.positions,
        scan.alleles,
        scan.gene_ids,
        scan.beta,
        scan.standard_error,
        scan.zscore,
        scan.pvalue,
        scan.sample_size,
        scan.covariate_names,
        scan.phenotype_name,
        "lmm_scan")
end

function _pairwise_ld(genotypes::AbstractMatrix{<:Real}, locus_1::Int, locus_2::Int)
    n = 0
    sum_x = 0.0
    sum_y = 0.0
    sum_xx = 0.0
    sum_yy = 0.0
    sum_xy = 0.0

    @inbounds for row in axes(genotypes, 1)
        x = genotypes[row, locus_1]
        y = genotypes[row, locus_2]
        if isfinite(x) && isfinite(y)
            n += 1
            sum_x += x
            sum_y += y
            sum_xx += x * x
            sum_yy += y * y
            sum_xy += x * y
        end
    end

    n < 3 && return (0.0, 0.0, 0.0)
    n_minus_one = max(n - 1, 1)
    var_x = (sum_xx - (sum_x * sum_x) / n) / n_minus_one
    var_y = (sum_yy - (sum_y * sum_y) / n) / n_minus_one
    (var_x <= eps(Float64) || var_y <= eps(Float64)) && return (0.0, 0.0, 0.0)

    cov_xy = (sum_xy - (sum_x * sum_y) / n) / n_minus_one
    r = clamp(cov_xy / sqrt(var_x * var_y), -1.0, 1.0)
    d = cov_xy
    return (d, r, r^2)
end

function _standardize_genotype_columns!(X::AbstractMatrix{<:AbstractFloat})
    nobs, ncols = size(X)
    for column in 1:ncols
        col = @view X[:, column]
        finite_sum = 0.0
        n_valid = 0
        @inbounds for row in 1:nobs
            value = col[row]
            if isfinite(value)
                finite_sum += value
                n_valid += 1
            end
        end

        if n_valid == 0
            fill!(col, 0.0)
            continue
        end

        mean_value = finite_sum / n_valid
        ss = 0.0
        @inbounds for row in 1:nobs
            value = col[row]
            centered = isfinite(value) ? (value - mean_value) : 0.0
            col[row] = centered
            ss += centered * centered
        end

        std_value = sqrt(ss / max(nobs - 1, 1))
        if !(isfinite(std_value) && std_value > eps(Float64))
            fill!(col, 0.0)
        else
            col ./= std_value
        end
    end
    return X
end

function _conditional_analysis_from_standardized(result::GWASResult, X::AbstractMatrix{<:Real}, cond::AbstractVector{Int}; method::String=string(result.method, "+conditional"))
    p = size(X, 2)
    isempty(cond) && return _subset_gwas_result(result, collect(1:p), method)

    Xc = X[:, cond]
    nobs = max(size(X, 1) - 1, 1)
    Rcc = Symmetric((Xc' * Xc) / nobs + 1e-6 * I)
    Rjc = (X' * Xc) / nobs
    adj = Rcc \ Float64.(result.beta[cond])

    beta_adj = Float64.(result.beta) .- Rjc * adj
    se_adj = Float64.(result.standard_error)
    z_adj = beta_adj ./ max.(se_adj, eps(Float64))
    p_adj = 2.0 .* ccdf.(Ref(Normal()), abs.(z_adj))

    return GWASResult(
        result.snp_ids,
        result.chromosomes,
        result.positions,
        result.alleles,
        result.gene_ids,
        beta_adj,
        se_adj,
        z_adj,
        p_adj,
        result.sample_size,
        result.covariate_names,
        result.phenotype_name,
        method)
end

function _joint_analysis_from_standardized(result::GWASResult, X::AbstractMatrix{<:Real}, idx::AbstractVector{Int}; ridge::Real=1e-4, method::String="COJO_joint")
    isempty(idx) && return _subset_gwas_result(result, Int[], method)

    Xjoint = X[:, idx]
    nobs = max(size(X, 1) - 1, 1)
    R = Symmetric((Xjoint' * Xjoint) / nobs + Float64(ridge) * I)
    beta_marg = Float64.(result.beta[idx])
    beta_joint = R \ beta_marg

    invR = try
        inv(R)
    catch
        Matrix{Float64}(I, length(idx), length(idx))
    end
    se_scale = median(Float64.(result.standard_error[idx]))
    se_joint = sqrt.(abs.(diag(invR))) .* max(se_scale, eps(Float64))
    z_joint = beta_joint ./ max.(se_joint, eps(Float64))
    p_joint = 2.0 .* ccdf.(Ref(Normal()), abs.(z_joint))

    joint = _subset_gwas_result(result, idx, method)
    return GWASResult(
        joint.snp_ids,
        joint.chromosomes,
        joint.positions,
        joint.alleles,
        joint.gene_ids,
        beta_joint,
        se_joint,
        z_joint,
        p_joint,
        joint.sample_size,
        joint.covariate_names,
        joint.phenotype_name,
        method)
end

"""
    ld_clumping(result, genotypes; p_threshold=5e-8, r2_threshold=0.2, window=250_000)

Prune correlated GWAS hits by linkage disequilibrium.
"""
function ld_clumping(result::GWASResult, genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; p_threshold::Real=5e-8, r2_threshold::Real=0.2, window::Int=250_000)
    G = _matrix(genotypes)
    snp_to_col = Dict{String,Int}()
    if genotypes isa GenotypeMatrix
        for (col, snp_id) in enumerate(String.(genotypes.bim.snp_id))
            snp_to_col[snp_id] = col
        end
    elseif size(G, 2) == length(result.snp_ids)
        for col in 1:size(G, 2)
            snp_to_col[String(result.snp_ids[col])] = col
        end
    else
        for col in 1:size(G, 2)
            snp_to_col["snp_$(col)"] = col
        end
    end

    result_to_col = [get(snp_to_col, String(result.snp_ids[i]), 0) for i in eachindex(result.snp_ids)]
    order = sortperm(result.pvalue)
    kept = Int[]
    for index in order
        result.pvalue[index] <= p_threshold || continue
        blocked = false
        idx_col = result_to_col[index]
        for existing in kept
            normalise_chromosome(result.chromosomes[existing]) == normalise_chromosome(result.chromosomes[index]) || continue
            abs(result.positions[existing] - result.positions[index]) <= window || continue
            ex_col = result_to_col[existing]
            if ex_col == 0 || idx_col == 0
                continue
            end
            _, _, r2 = _pairwise_ld(G, ex_col, idx_col)
            if r2 >= r2_threshold
                blocked = true
                break
            end
        end
        blocked || push!(kept, index)
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GWASResult(result.snp_ids[kept], result.chromosomes[kept], result.positions[kept], result.alleles[kept], result.gene_ids[kept], result.beta[kept], result.standard_error[kept], result.zscore[kept], result.pvalue[kept], result.sample_size, result.covariate_names, result.phenotype_name, result.method), "ld_clumping")
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
    _validate_gwas_inputs(G)
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
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, matrix * weights, "calculate_prs")
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
    G = _matrix(genotypes)
    _validate_gwas_inputs(G; ld_matrix=ld_matrix)
    weights = hasproperty(summary_stats, :beta) ? Float64.(summary_stats.beta) : throw(ArgumentError("summary stats must expose beta"))
    length(weights) == size(ld_matrix, 1) || throw(DimensionMismatch("summary-stat beta length must match LD matrix dimensions"))
    adjusted = _solve_ld_system(ld_matrix, weights; ridge=ridge)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, G * adjusted, "prs_ldpred")
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
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (; best_threshold=Float64(best_threshold), best_window=Int(best_window), best_score=best_score, scores=scores), "prs_cross_validation")
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
    index_maps = [Dict{String,Int}(String(result.snp_ids[i]) => i for i in eachindex(result.snp_ids)) for result in results]
    first_map = index_maps[1]
    meta_snp_ids = String[]
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
        snp_key = String(snp_id)
        per_study = [(results[i], get(index_maps[i], snp_key, 0)) for i in eachindex(results)]
        per_study = [(result, idx) for (result, idx) in per_study if idx != 0]
        first_result = results[1]
        first_index = get(first_map, snp_key, 0)
        first_index == 0 && continue
        betas = [Float64(result.beta[index]) for (result, index) in per_study]
        ses = [Float64(result.standard_error[index]) for (result, index) in per_study]
        isempty(betas) && continue
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
        push!(meta_snp_ids, snp_key)
        push!(chromosomes, first_result.chromosomes[first_index])
        push!(positions, first_result.positions[first_index])
        push!(alleles, first_result.alleles[first_index])
        push!(gene_ids, first_result.gene_ids[first_index])
    end
    qvalue = benjamini_hochberg(p)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, MetaAnalysisResult(_pooled_strings(meta_snp_ids), _pooled_strings(chromosomes), positions, alleles, _pooled_strings(gene_ids), beta, se, z, p, qvalue, tau2, i2, length(results), random_effects ? "random_effects" : "fixed_effects"), "meta_analyze")
end

function _peak_intervals(peakset::PeakSet)
    intervals = Vector{GenomicInterval}(undef, length(peakset.peaks))
    for (index, peak) in enumerate(peakset.peaks)
        intervals[index] = GenomicInterval(normalise_chromosome(peak.chrom), peak.left, peak.right, '+', Dict{String,Any}("peak_index" => index, "score" => peak.score, "pvalue" => peak.pvalue, "qvalue" => peak.qvalue))
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
        query = GenomicInterval(normalise_chromosome(result.chromosomes[index]), max(1, result.positions[index] - flank), result.positions[index] + flank, '+', Dict{String,Any}("snp_id" => result.snp_ids[index], "beta" => result.beta[index], "pvalue" => result.pvalue[index]))
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
                peak_qvalue = get(hit.metadata, "qvalue", NaN)))
        end
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, DataFrame(rows), "overlap_gwas_peaks")
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

"""
    LDSCResult

Result from LD-score regression style heritability estimation.
"""
struct LDSCResult <: AbstractAnalysisResult
    h2::Float64
    h2_se::Float64
    intercept::Float64
    intercept_se::Float64
    mean_chi2::Float64
    lambda_gc::Float64
    n_snps_used::Int
    n_h2_blocks::Int
    provenance::ResultProvenance
end

LDSCResult(h2, h2_se, intercept, intercept_se, mean_chi2, lambda_gc, n_snps_used, n_h2_blocks) =
    LDSCResult(h2, h2_se, intercept, intercept_se, mean_chi2, lambda_gc, n_snps_used, n_h2_blocks, provenance_record("LDSCResult", "gwas"))

LDSCResult(h2::Float64, h2_se::Float64, intercept::Float64, intercept_se::Float64, mean_chi2::Float64, lambda_gc::Float64, n_snps_used::Int, n_h2_blocks::Int) =
    LDSCResult(h2, h2_se, intercept, intercept_se, mean_chi2, lambda_gc, n_snps_used, n_h2_blocks, provenance_record("LDSCResult", "GWAS/ldsc_heritability"))

@inline function _ldsc_empty_result(n_blocks::Int=0; reason::AbstractString="insufficient valid SNPs", parameters::NamedTuple=NamedTuple())
    return LDSCResult(
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        0,
        n_blocks,
        provenance_record(
            "LDSCResult",
            "GWAS/ldsc_heritability";
            status=:warn,
            warnings=[String(reason)],
            fallbacks=["returned NaN placeholder result"],
            parameters=parameters))
end

"""
    ldsc_heritability(sumstats, ld_scores; n_blocks=200, intercept_only=false)

Estimate SNP heritability using an LDSC-style weighted regression.
"""
function ldsc_heritability(sumstats::GWASResult, ld_scores::AbstractVector{<:Real}; n_blocks::Int=200, intercept_only::Bool=false)
    length(sumstats.pvalue) == length(ld_scores) || throw(DimensionMismatch("ld_scores length must match summary statistics"))

    pvalues = clamp.(Float64.(sumstats.pvalue), eps(Float64), 1.0)
    ld = Float64.(ld_scores)
    valid = isfinite.(pvalues) .& isfinite.(ld) .& (ld .>= 0.0)
    valid_count = sum(valid)
    valid_count >= 10 || return _ldsc_empty_result(n_blocks; reason="insufficient valid SNPs after filtering", parameters=(n_blocks=Int(n_blocks), intercept_only=Bool(intercept_only), valid_snp_count=Int(valid_count)))

    p = pvalues[valid]
    l = ld[valid]
    chi2 = quantile.(Ref(Chisq(1)), 1.0 .- p)
    mean_chi2 = mean(chi2)
    lambda_gc = median(chi2) / quantile(Chisq(1), 0.5)

    weights = 1.0 ./ max.(l .+ 1.0, 1e-8)
    weights ./= mean(weights)

    if intercept_only
        intercept = sum(weights .* chi2) / sum(weights)
        intercept_se = sqrt(max(var(chi2; corrected=true), eps(Float64)) / length(chi2))
        return LDSCResult(
            NaN,
            NaN,
            intercept,
            intercept_se,
            mean_chi2,
            lambda_gc,
            length(chi2),
            1,
            provenance_record(
                "LDSCResult",
                "GWAS/ldsc_heritability";
                notes=["intercept-only LDSC path used"],
                parameters=(n_blocks=Int(n_blocks), intercept_only=true, valid_snp_count=length(chi2))))
    end

    X = hcat(ones(Float64, length(l)), l)
    Xw = X .* sqrt.(weights)
    yw = chi2 .* sqrt.(weights)
    beta = Xw \ yw

    cov_beta = try
        inv(Symmetric(Xw' * Xw))
    catch
        return _ldsc_empty_result(n_blocks; reason="weighted regression design was singular", parameters=(n_blocks=Int(n_blocks), intercept_only=false, valid_snp_count=length(chi2)))
    end

    N = max(sumstats.sample_size, 1)
    M = length(chi2)
    slope = beta[2]
    h2 = max(slope * M / N, 0.0)
    h2_se = sqrt(max(cov_beta[2, 2], 0.0)) * M / N
    intercept = beta[1]
    intercept_se = sqrt(max(cov_beta[1, 1], 0.0))
    return LDSCResult(
        h2,
        h2_se,
        intercept,
        intercept_se,
        mean_chi2,
        lambda_gc,
        M,
        n_blocks,
        provenance_record(
            "LDSCResult",
            "GWAS/ldsc_heritability";
            notes=["weighted regression on LD scores"],
            parameters=(n_blocks=Int(n_blocks), intercept_only=false, valid_snp_count=M, sample_size=Int(N))))
end

"""
    estimate_heritability_greml(genotypes, phenotype; kinship=nothing, constrained=true)

Estimate SNP heritability using a Haseman-Elston / GREML-style moment estimator.
"""
function estimate_heritability_greml(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, phenotype::AbstractVector{<:Real}; kinship::Union{Nothing,AbstractMatrix{<:Real}}=nothing, constrained::Bool=true)
    G = _matrix(genotypes)
    n = size(G, 1)
    length(phenotype) == n || throw(DimensionMismatch("phenotype length must match sample count"))

    y = Float64.(phenotype)
    valid_samples = isfinite.(y)
    sum(valid_samples) >= 5 || return NaN

    Gv = G[valid_samples, :]
    yv = y[valid_samples]
    yv .-= mean(yv)

    K = kinship === nothing ? calculate_grm(Gv) : Matrix{Float64}(kinship)
    size(K, 1) == length(yv) == size(K, 2) || throw(DimensionMismatch("kinship matrix must match phenotype/sample dimensions"))

    kvals = Float64[]
    yprod = Float64[]
    for i in 1:(length(yv) - 1), j in (i + 1):length(yv)
        push!(kvals, K[i, j])
        push!(yprod, yv[i] * yv[j])
    end

    var(kvals; corrected=true) > eps(Float64) || return 0.0
    slope = cov(kvals, yprod; corrected=true) / var(kvals; corrected=true)
    total_var = max(var(yv; corrected=false), eps(Float64))
    h2 = slope / total_var
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, constrained ? clamp(h2, 0.0, 1.0) : h2, "estimate_heritability_greml")
end

"""
    partitioned_heritability(sumstats, ld_scores; annotations=nothing, min_snps=20, n_blocks=200)

Estimate LDSC-style heritability partitioned by category labels.
"""
function partitioned_heritability(sumstats::GWASResult, ld_scores::AbstractVector{<:Real}; annotations::Union{Nothing,AbstractVector{<:AbstractString}}=nothing, min_snps::Int=20, n_blocks::Int=200)
    length(sumstats.pvalue) == length(ld_scores) || throw(DimensionMismatch("ld_scores length must match summary statistics"))
    labels = annotations === nothing ? String.(sumstats.gene_ids) : String.(annotations)
    length(labels) == length(sumstats.pvalue) || throw(DimensionMismatch("annotation length must match summary statistics"))

    output = Dict{String,LDSCResult}()
    for label in unique(labels)
        mask = labels .== label
        if sum(mask) < min_snps
            output[label] = _ldsc_empty_result(n_blocks)
            continue
        end
        subset = GWASResult(
            sumstats.snp_ids[mask],
            sumstats.chromosomes[mask],
            sumstats.positions[mask],
            sumstats.alleles[mask],
            sumstats.gene_ids[mask],
            sumstats.beta[mask],
            sumstats.standard_error[mask],
            sumstats.zscore[mask],
            sumstats.pvalue[mask],
            sumstats.sample_size,
            sumstats.covariate_names,
            sumstats.phenotype_name,
            sumstats.method)
        output[label] = ldsc_heritability(subset, ld_scores[mask]; n_blocks=n_blocks)
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, output, "partitioned_heritability")
end

"""
    ldsc_genetic_correlation(traits)

Estimate genetic correlation between two GWAS result sets from matched z-scores.
"""
function ldsc_genetic_correlation(traits::Vector{GWASResult}; min_snps::Int=100)
    length(traits) == 2 || throw(ArgumentError("genetic correlation requires exactly two traits"))
    t1, t2 = traits

    idx2 = Dict{String,Int}(String(snp) => index for (index, snp) in enumerate(t2.snp_ids))
    z1 = Float64[]
    z2 = Float64[]
    for index in eachindex(t1.snp_ids)
        snp = String(t1.snp_ids[index])
        haskey(idx2, snp) || continue
        j = idx2[snp]
        z1_i = Float64(t1.zscore[index])
        z2_i = Float64(t2.zscore[j])
        (isfinite(z1_i) && isfinite(z2_i)) || continue
        push!(z1, z1_i)
        push!(z2, z2_i)
    end

    n = length(z1)
    n >= min_snps || return (NaN, NaN)
    rg = cor(z1, z2)
    isfinite(rg) || return (NaN, NaN)
    se = sqrt(max((1.0 - rg^2) / max(n - 2, 1), 0.0))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (rg, se), "ldsc_genetic_correlation")
end

"""
    SuSiEResult

Lightweight fine-mapping result object with PIPs and credible sets.
"""
struct SuSiEResult <: AbstractAnalysisResult
    pip::Vector{Float64}
    posterior_mean::Vector{Float64}
    posterior_sd::Vector{Float64}
    credible_sets::Vector{Vector{Int}}
    n_effects::Int
    residual_variance::Float64
    heritability::Float64
    ldp_logbf::Vector{Float64}
    provenance::ResultProvenance
end

SuSiEResult(pip, posterior_mean, posterior_sd, credible_sets, n_effects, residual_variance, heritability, ldp_logbf) =
    SuSiEResult(pip, posterior_mean, posterior_sd, credible_sets, n_effects, residual_variance, heritability, ldp_logbf, provenance_record("SuSiEResult", "gwas"))

"""
    posterior_inclusion_probability(beta, se, ld_matrix; prior_variance=0.01)

Compute approximate PIPs with Wakefield ABFs and a simple LD shrinkage factor.
"""
function posterior_inclusion_probability(beta::AbstractVector{<:Real}, se::AbstractVector{<:Real}, ld_matrix::AbstractMatrix{<:Real}; prior_variance::Real=0.01)
    length(beta) == length(se) || throw(DimensionMismatch("beta and se lengths must match"))
    size(ld_matrix, 1) == size(ld_matrix, 2) == length(beta) || throw(DimensionMismatch("ld_matrix must be square with side length equal to number of variants"))

    p = length(beta)
    pip = zeros(Float64, p)
    valid = isfinite.(beta) .& isfinite.(se) .& (Float64.(se) .> 0.0)
    sum(valid) > 0 || return pip

    b = Float64.(beta[valid])
    s = Float64.(se[valid])
    v = s .^ 2
    W = Float64(prior_variance)
    z = b ./ s
    abf = sqrt.(v ./ (v .+ W)) .* exp.((z .^ 2 .* W) ./ (2.0 .* (v .+ W)))

    ld_diag = abs.(diag(Matrix{Float64}(ld_matrix))[valid])
    abf ./= max.(1.0 .+ ld_diag, 1.0)
    pip_valid = abf ./ (sum(abf) + 1.0)
    pip[valid] = clamp.(pip_valid, 0.0, 1.0)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, pip, "posterior_inclusion_probability")
end

"""
    calculate_credible_set(pip; coverage=0.95)

Build an index set that reaches target posterior coverage.
"""
function calculate_credible_set(pip::AbstractVector{<:Real}; coverage::Real=0.95)
    order = sortperm(Float64.(pip); rev=true)
    cs = Int[]
    cumulative = 0.0
    for index in order
        value = clamp(Float64(pip[index]), 0.0, 1.0)
        value <= 0 && continue
        push!(cs, index)
        cumulative += value
        cumulative >= coverage && break
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, cs, "calculate_credible_set")
end

"""
    fine_map_susie(beta, se, ld_matrix; n_effects=10, coverage=0.95, max_iter=100, tol=1e-6)

Run a lightweight SuSiE-style approximation from summary statistics.
"""
function fine_map_susie(beta::AbstractVector{<:Real}, se::AbstractVector{<:Real}, ld_matrix::AbstractMatrix{<:Real}; n_effects::Int=10, coverage::Real=0.95, max_iter::Int=100, tol::Real=1e-6)
    _ = max_iter
    _ = tol
    p = length(beta)
    p > 0 || throw(ArgumentError("empty summary statistics"))

    pip = posterior_inclusion_probability(beta, se, ld_matrix)
    post_mean = Float64.(beta) .* pip
    post_sd = sqrt.(max.(Float64.(se) .^ 2 .* pip .* (1.0 .- pip), eps(Float64)))
    credible = [calculate_credible_set(pip; coverage=coverage)]
    resid_var = var(Float64.(beta) .- post_mean; corrected=false)
    h2 = sum(post_mean .^ 2) / max(sum(Float64.(beta) .^ 2), eps(Float64))
    logbf = log.(pip .+ eps(Float64)) .- log.(1.0 .- pip .+ eps(Float64))

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, SuSiEResult(pip, post_mean, post_sd, credible, min(n_effects, p), resid_var, h2, logbf), "fine_map_susie")
end

"""
    MRResult

Result for Mendelian randomization analyses.
"""
struct MRResult <: AbstractAnalysisResult
    estimate::Float64
    se::Float64
    pvalue::Float64
    method::String
    egger_intercept::Union{Float64,Nothing}
    egger_intercept_p::Union{Float64,Nothing}
    provenance::ResultProvenance
end

MRResult(estimate, se, pvalue, method, egger_intercept, egger_intercept_p) =
    MRResult(estimate, se, pvalue, method, egger_intercept, egger_intercept_p, provenance_record("MRResult", "gwas"))

@inline function _mr_valid(beta_exp, se_exp, beta_out, se_out)
    return isfinite(beta_exp) && isfinite(se_exp) && isfinite(beta_out) && isfinite(se_out) && abs(beta_exp) > eps(Float64) && se_out > 0
end

"""
    mr_two_sample(beta_exp, se_exp, beta_out, se_out)

Two-sample IVW Mendelian randomization.
"""
function mr_two_sample(beta_exp::AbstractVector{<:Real}, se_exp::AbstractVector{<:Real}, beta_out::AbstractVector{<:Real}, se_out::AbstractVector{<:Real}; exposure_name::String="exposure", outcome_name::String="outcome")
    _ = exposure_name
    _ = outcome_name
    length(beta_exp) == length(se_exp) == length(beta_out) == length(se_out) || throw(DimensionMismatch("all MR vectors must have equal length"))

    valid = [_mr_valid(Float64(beta_exp[i]), Float64(se_exp[i]), Float64(beta_out[i]), Float64(se_out[i])) for i in eachindex(beta_exp)]
    sum(valid) >= 3 || return MRResult(NaN, NaN, NaN, "IVW", nothing, nothing)

    b_exp = Float64.(beta_exp[valid])
    s_exp = Float64.(se_exp[valid])
    b_out = Float64.(beta_out[valid])
    s_out = Float64.(se_out[valid])

    ratio = b_out ./ b_exp
    var_ratio = (s_out .^ 2 ./ (b_exp .^ 2)) .+ ((b_out .^ 2 .* s_exp .^ 2) ./ (b_exp .^ 4))
    w = 1.0 ./ max.(var_ratio, eps(Float64))

    estimate = sum(w .* ratio) / sum(w)
    se = sqrt(1.0 / sum(w))
    z = estimate / se
    p = 2.0 * ccdf(Normal(), abs(z))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, MRResult(estimate, se, p, "IVW", nothing, nothing), "mr_two_sample")
end

"""
    mr_egger(beta_exp, se_exp, beta_out, se_out)

MR-Egger weighted regression with intercept.
"""
function mr_egger(beta_exp::AbstractVector{<:Real}, se_exp::AbstractVector{<:Real}, beta_out::AbstractVector{<:Real}, se_out::AbstractVector{<:Real})
    length(beta_exp) == length(se_exp) == length(beta_out) == length(se_out) || throw(DimensionMismatch("all MR vectors must have equal length"))
    valid = [_mr_valid(Float64(beta_exp[i]), Float64(se_exp[i]), Float64(beta_out[i]), Float64(se_out[i])) for i in eachindex(beta_exp)]
    sum(valid) >= 3 || return MRResult(NaN, NaN, NaN, "MR-Egger", NaN, NaN)

    x = Float64.(beta_exp[valid])
    y = Float64.(beta_out[valid])
    w = 1.0 ./ max.(Float64.(se_out[valid]) .^ 2, eps(Float64))

    X = hcat(ones(Float64, length(x)), x)
    Xw = X .* sqrt.(w)
    yw = y .* sqrt.(w)
    beta_hat = Xw \ yw

    cov_beta = try
        inv(Symmetric(Xw' * Xw))
    catch
        return MRResult(NaN, NaN, NaN, "MR-Egger", NaN, NaN)
    end

    intercept = beta_hat[1]
    slope = beta_hat[2]
    intercept_se = sqrt(max(cov_beta[1, 1], 0.0))
    slope_se = sqrt(max(cov_beta[2, 2], 0.0))

    p_slope = 2.0 * ccdf(Normal(), abs(slope / max(slope_se, eps(Float64))))
    p_intercept = 2.0 * ccdf(Normal(), abs(intercept / max(intercept_se, eps(Float64))))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, MRResult(slope, slope_se, p_slope, "MR-Egger", intercept, p_intercept), "mr_egger")
end

"""
    mr_pleiotropy_test(beta_exp, se_exp, beta_out, se_out)

Global pleiotropy test based on Cochran's Q statistic around IVW estimate.
"""
function mr_pleiotropy_test(beta_exp::AbstractVector{<:Real}, se_exp::AbstractVector{<:Real}, beta_out::AbstractVector{<:Real}, se_out::AbstractVector{<:Real})
    length(beta_exp) == length(se_exp) == length(beta_out) == length(se_out) || throw(DimensionMismatch("all MR vectors must have equal length"))
    valid = [_mr_valid(Float64(beta_exp[i]), Float64(se_exp[i]), Float64(beta_out[i]), Float64(se_out[i])) for i in eachindex(beta_exp)]
    sum(valid) >= 3 || return (NaN, NaN)

    b_exp = Float64.(beta_exp[valid])
    s_exp = Float64.(se_exp[valid])
    b_out = Float64.(beta_out[valid])
    s_out = Float64.(se_out[valid])
    ivw = mr_two_sample(b_exp, s_exp, b_out, s_out)

    ratio = b_out ./ b_exp
    var_ratio = (s_out .^ 2 ./ (b_exp .^ 2)) .+ ((b_out .^ 2 .* s_exp .^ 2) ./ (b_exp .^ 4))
    q = sum((ratio .- ivw.estimate) .^ 2 ./ max.(var_ratio, eps(Float64)))
    df = max(length(ratio) - 1, 1)
    p = ccdf(Chisq(df), q)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (q, p), "mr_pleiotropy_test")
end

@inline function _subset_gwas_result(result::GWASResult, indices::AbstractVector{<:Integer}, method::String)
    idx = Int.(indices)
    return GWASResult(
        result.snp_ids[idx],
        result.chromosomes[idx],
        result.positions[idx],
        result.alleles[idx],
        result.gene_ids[idx],
        result.beta[idx],
        result.standard_error[idx],
        result.zscore[idx],
        result.pvalue[idx],
        result.sample_size,
        result.covariate_names,
        result.phenotype_name,
        method)
end

"""
    conditional_analysis(result, genotypes, conditional_indices; window=1_000_000)

Approximate COJO-style conditional effect adjustment using LD with conditioning SNPs.
"""
function conditional_analysis(result::GWASResult, genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, conditional_indices::AbstractVector{Int}; window::Int=1_000_000)
    _ = window
    G = _matrix(genotypes)
    p = size(G, 2)
    cond = unique(Int.(conditional_indices))
    all(index -> 1 <= index <= p, cond) || throw(BoundsError("conditional index out of bounds"))
    X = Matrix{Float64}(G)
    _standardize_genotype_columns!(X)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, _conditional_analysis_from_standardized(result, X, cond; method=string(result.method, "+conditional")), "conditional_analysis")
end

"""
    joint_analysis(result, genotypes, joint_indices; ridge=1e-4)

Joint effect estimate for a selected SNP set using LD inversion.
"""
function joint_analysis(result::GWASResult, genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, joint_indices::AbstractVector{Int}; ridge::Real=1e-4)
    G = _matrix(genotypes)
    p = size(G, 2)
    idx = unique(Int.(joint_indices))
    all(index -> 1 <= index <= p, idx) || throw(BoundsError("joint index out of bounds"))
    X = Matrix{Float64}(G)
    _standardize_genotype_columns!(X)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, _joint_analysis_from_standardized(result, X, idx; ridge=ridge, method="COJO_joint"), "joint_analysis")
end

"""
    cojo_stepwise(result, genotypes; p_threshold=5e-8, r2_threshold=0.01, window=10_000_000, max_snps=100)

Stepwise conditional selection akin COJO.
"""
function cojo_stepwise(result::GWASResult, genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; p_threshold::Real=5e-8, r2_threshold::Real=0.01, window::Int=10_000_000, max_snps::Int=100)
    G = _matrix(genotypes)
    X = Matrix{Float64}(G)
    _standardize_genotype_columns!(X)

    result_index = Dict{String,Int}(String(result.snp_ids[i]) => i for i in eachindex(result.snp_ids))
    genotype_index = Dict{String,Int}()
    if genotypes isa GenotypeMatrix
        for (col, snp_id) in enumerate(String.(genotypes.bim.snp_id))
            genotype_index[snp_id] = col
        end
    else
        for i in 1:size(G, 2)
            genotype_index["snp_$(i)"] = i
            i <= length(result.snp_ids) && (genotype_index[String(result.snp_ids[i])] = i)
        end
    end

    selected_ids = String[]
    selected_set = Set{String}()
    current = result
    nobs = max(size(X, 1) - 1, 1)

    for _ in 1:max_snps
        best = 0
        best_p = Inf
        for i in eachindex(current.pvalue)
            snp_id = String(current.snp_ids[i])
            snp_id in selected_set && continue
            pval = current.pvalue[i]
            if isfinite(pval) && pval < best_p
                best_p = pval
                best = i
            end
        end

        (best == 0 || best_p > p_threshold) && break
        best_id = String(current.snp_ids[best])
        best_result_idx = get(result_index, best_id, 0)
        best_geno_idx = get(genotype_index, best_id, 0)
        (best_result_idx == 0 || best_geno_idx == 0) && break

        blocked = false
        for existing_id in selected_ids
            ex_result_idx = get(result_index, existing_id, 0)
            ex_geno_idx = get(genotype_index, existing_id, 0)
            (ex_result_idx == 0 || ex_geno_idx == 0) && continue

            if normalise_chromosome(result.chromosomes[ex_result_idx]) == normalise_chromosome(result.chromosomes[best_result_idx]) &&
               abs(result.positions[ex_result_idx] - result.positions[best_result_idx]) <= window
                r = dot(@view(X[:, ex_geno_idx]), @view(X[:, best_geno_idx])) / nobs
                r2 = r * r
                if r2 >= r2_threshold
                    blocked = true
                    break
                end
            end
        end

        if !blocked
            push!(selected_ids, best_id)
            push!(selected_set, best_id)
        end

        if !isempty(selected_ids)
            cond_idx = [result_index[snp_id] for snp_id in selected_ids if haskey(result_index, snp_id)]
            current = _conditional_analysis_from_standardized(result, X, cond_idx; method=string(result.method, "+conditional"))
        end
    end

    selected_idx = [result_index[snp_id] for snp_id in selected_ids if haskey(result_index, snp_id)]
    isempty(selected_idx) && return _subset_gwas_result(result, Int[], "COJO_stepwise")
    joint = _joint_analysis_from_standardized(result, X, selected_idx; ridge=1e-4, method="COJO_joint")
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GWASResult(joint.snp_ids, joint.chromosomes, joint.positions, joint.alleles, joint.gene_ids, joint.beta, joint.standard_error, joint.zscore, joint.pvalue, joint.sample_size, joint.covariate_names, joint.phenotype_name, "COJO_stepwise"), "cojo_stepwise")
end

"""
    gwas_pca(genotypes; n_components=20, maf_threshold=0.01, return_loadings=false)

Perform PCA on genotype data for population structure covariates.

Set `return_loadings=true` to also return PCA loadings and scaling parameters
for out-of-sample projection via `project_pca`.
"""
function gwas_pca(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; n_components::Int=20, maf_threshold::Real=0.01, return_loadings::Bool=false)
    G = _matrix(genotypes)
    mask, _ = calculate_maf(G; min_maf=maf_threshold, max_maf=1.0, return_frequency=true)
    idx = findall(mask)
    isempty(idx) && throw(ArgumentError("no SNPs passed MAF threshold"))

    X = Matrix{Float64}(G[:, idx])
    _impute_missing_by_mean!(X)
    center = vec(mean(X, dims=1))
    scale = vec(std(X, dims=1))
    scale = max.(scale, eps(Float64))
    X .-= reshape(center, 1, :)
    X ./= reshape(scale, 1, :)

    decomp = svd(X; full=false)
    k = min(n_components, length(decomp.S))
    pcs = decomp.U[:, 1:k] * Diagonal(decomp.S[1:k])
    loadings = decomp.V[:, 1:k]
    total_var = max(sum(decomp.S .^ 2), eps(Float64))
    var_exp = (decomp.S[1:k] .^ 2) ./ total_var
    if return_loadings
        return (pcs, var_exp, loadings, center, scale, idx)
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (pcs, var_exp), "gwas_pca")
end

"""
    project_pca(new_genotypes, loadings; center=nothing, scale=nothing, snp_indices=nothing)

Project out-of-sample genotypes onto previously fitted PCA loadings.
"""
function project_pca(
    new_genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}},
    loadings::AbstractMatrix{<:Real};
    center::Union{Nothing,AbstractVector{<:Real}}=nothing,
    scale::Union{Nothing,AbstractVector{<:Real}}=nothing,
    snp_indices::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    G = _matrix(new_genotypes)
    X = snp_indices === nothing ? Matrix{Float64}(G) : Matrix{Float64}(G[:, Int.(snp_indices)])
    _impute_missing_by_mean!(X)

    size(X, 2) == size(loadings, 1) || throw(DimensionMismatch("loadings row count must match number of projected SNPs"))
    if center !== nothing
        length(center) == size(X, 2) || throw(DimensionMismatch("center length must match number of projected SNPs"))
        X .-= reshape(Float64.(center), 1, :)
    end
    if scale !== nothing
        length(scale) == size(X, 2) || throw(DimensionMismatch("scale length must match number of projected SNPs"))
        denom = max.(Float64.(scale), eps(Float64))
        X ./= reshape(denom, 1, :)
    end
    return X * Matrix{Float64}(loadings)
end

"""
    calculate_ibd(genotypes)

Compute pairwise IBS similarity matrix.
"""
function calculate_ibd(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}})
    G = _matrix(genotypes)
    n = size(G, 1)
    ibd = Matrix{Float64}(I, n, n)

    if _thread_enabled(n, true; min_grain=64)
        @threads for i in 1:(n - 1)
            xi = @view G[i, :]
            for j in (i + 1):n
                xj = @view G[j, :]
                n_valid = 0
                n_equal = 0
                @inbounds for k in eachindex(xi)
                    left = xi[k]
                    right = xj[k]
                    if isfinite(left) && isfinite(right)
                        n_valid += 1
                        n_equal += (left == right) ? 1 : 0
                    end
                end
                ibd[i, j] = n_valid == 0 ? NaN : n_equal / n_valid
            end
        end
    else
        for i in 1:(n - 1)
            xi = @view G[i, :]
            for j in (i + 1):n
                xj = @view G[j, :]
                n_valid = 0
                n_equal = 0
                @inbounds for k in eachindex(xi)
                    left = xi[k]
                    right = xj[k]
                    if isfinite(left) && isfinite(right)
                        n_valid += 1
                        n_equal += (left == right) ? 1 : 0
                    end
                end
                ibd[i, j] = n_valid == 0 ? NaN : n_equal / n_valid
            end
        end
    end

    for i in 1:(n - 1)
        for j in (i + 1):n
            ibd[j, i] = ibd[i, j]
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, Symmetric(ibd), "calculate_ibd")
end

"""
    calculate_king_kinship(genotypes; ignore_hwe=true)

Approximate KING-style kinship matrix from centered genotypes.
"""
function calculate_king_kinship(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; ignore_hwe::Bool=true)
    _ = ignore_hwe
    G = _matrix(genotypes)
    n, p = size(G)

    allele_freq = zeros(Float64, p)
    for j in 1:p
        valid_count = 0
        allele_sum = 0.0
        @inbounds for i in 1:n
            value = G[i, j]
            if isfinite(value)
                valid_count += 1
                allele_sum += value
            end
        end
        allele_freq[j] = valid_count == 0 ? 0.0 : clamp(allele_sum / (2.0 * valid_count), 0.0, 1.0)
    end
    var_term = 2.0 .* allele_freq .* (1.0 .- allele_freq)
    informative_locus = var_term .> eps(Float64)

    kinship = zeros(Float64, n, n)

    # Diagonal terms are deterministic and race-free.
    for i in 1:n
        xi = @view G[i, :]
        num = 0.0
        den_sum = 0.0
        @inbounds for k in eachindex(xi)
            informative_locus[k] || continue
            value = xi[k]
            if isfinite(value)
                centered = value - 2.0 * allele_freq[k]
                num += centered * centered
                den_sum += var_term[k]
            end
        end
        den = 2.0 * den_sum
        kinship[i, i] = den <= eps(Float64) ? 0.0 : num / den
    end

    if _thread_enabled(n, true; min_grain=64)
        @threads for i in 1:(n - 1)
            xi = @view G[i, :]
            for j in (i + 1):n
                xj = @view G[j, :]
                num = 0.0
                den_sum = 0.0
                @inbounds for k in eachindex(xi)
                    informative_locus[k] || continue
                    left = xi[k]
                    right = xj[k]
                    if isfinite(left) && isfinite(right)
                        centered_left = left - 2.0 * allele_freq[k]
                        centered_right = right - 2.0 * allele_freq[k]
                        num += centered_left * centered_right
                        den_sum += var_term[k]
                    end
                end
                den = 2.0 * den_sum
                kinship[i, j] = den <= eps(Float64) ? 0.0 : num / den
            end
        end
    else
        for i in 1:(n - 1)
            xi = @view G[i, :]
            for j in (i + 1):n
                xj = @view G[j, :]
                num = 0.0
                den_sum = 0.0
                @inbounds for k in eachindex(xi)
                    informative_locus[k] || continue
                    left = xi[k]
                    right = xj[k]
                    if isfinite(left) && isfinite(right)
                        centered_left = left - 2.0 * allele_freq[k]
                        centered_right = right - 2.0 * allele_freq[k]
                        num += centered_left * centered_right
                        den_sum += var_term[k]
                    end
                end
                den = 2.0 * den_sum
                kinship[i, j] = den <= eps(Float64) ? 0.0 : num / den
            end
        end
    end

    for i in 1:(n - 1)
        for j in (i + 1):n
            kinship[j, i] = kinship[i, j]
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, Symmetric((kinship + kinship') ./ 2), "calculate_king_kinship")
end

"""
    detect_related_pairs(kinship; threshold=0.125, max_pairs=100000)

Return sample-index pairs above kinship threshold.
"""
function detect_related_pairs(kinship::AbstractMatrix{<:Real}; threshold::Real=0.125, max_pairs::Int=100000)
    n = size(kinship, 1)
    pairs = Tuple{Int,Int,Float64}[]
    for i in 1:n
        for j in (i + 1):n
            value = Float64(kinship[i, j])
            (isfinite(value) && value >= threshold) || continue
            push!(pairs, (i, j, value))
            length(pairs) >= max_pairs && break
        end
        length(pairs) >= max_pairs && break
    end
    sort!(pairs; by=pair -> pair[3], rev=true)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, pairs, "detect_related_pairs")
end

"""
    coloc_result

Result from approximate Bayesian colocalization.
"""
struct coloc_result
    posterior_prob_ab::Float64
    posterior_prob_a::Float64
    posterior_prob_b::Float64
    n_snps::Int
    hypothesis::String
end

"""
    coloc_abf(beta_a, se_a, beta_b, se_b, n_samples; prior_std=0.2)

Approximate colocalization posterior probabilities from ABFs.
"""
function coloc_abf(beta_a::AbstractVector{<:Real}, se_a::AbstractVector{<:Real}, beta_b::AbstractVector{<:Real}, se_b::AbstractVector{<:Real}, n_samples::Int; prior_std::Real=0.2)
    _ = n_samples
    length(beta_a) == length(se_a) == length(beta_b) == length(se_b) || throw(DimensionMismatch("all coloc vectors must have equal length"))

    valid = isfinite.(beta_a) .& isfinite.(se_a) .& isfinite.(beta_b) .& isfinite.(se_b) .& (Float64.(se_a) .> 0.0) .& (Float64.(se_b) .> 0.0)
    n = sum(valid)
    n > 0 || return coloc_result(NaN, NaN, NaN, 0, "H0")

    ba = Float64.(beta_a[valid])
    sa = Float64.(se_a[valid])
    bb = Float64.(beta_b[valid])
    sb = Float64.(se_b[valid])

    W = prior_std^2
    va = sa .^ 2
    vb = sb .^ 2
    za = ba ./ sa
    zb = bb ./ sb
    abf_a = sqrt.(va ./ (va .+ W)) .* exp.((za .^ 2 .* W) ./ (2.0 .* (va .+ W)))
    abf_b = sqrt.(vb ./ (vb .+ W)) .* exp.((zb .^ 2 .* W) ./ (2.0 .* (vb .+ W)))

    p1 = 1e-4
    p2 = 1e-4
    p12 = 1e-5

    h0 = 1.0
    h1 = p1 * sum(abf_a)
    h2 = p2 * sum(abf_b)
    h4 = p12 * sum(abf_a .* abf_b)
    h3 = max(p1 * p2 * sum(abf_a) * sum(abf_b) - h4, eps(Float64))

    total = h0 + h1 + h2 + h3 + h4
    p_h0, p_h1, p_h2, p_h3, p_h4 = h0 / total, h1 / total, h2 / total, h3 / total, h4 / total
    labels = ["H0", "H1", "H2", "H3", "H4"]
    probs = [p_h0, p_h1, p_h2, p_h3, p_h4]
    hypothesis = labels[argmax(probs)]

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, coloc_result(p_h4, p_h1, p_h2, n, hypothesis), "coloc_abf")
end

"""
    ihs_score(genotypes, positions, haplotype_a, haplotype_b; window=100000)

Compute a simple iHS-like score from local haplotype homozygosity asymmetry.
"""
function ihs_score(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, positions::AbstractVector{<:Integer}, haplotype_a::AbstractVector{<:Integer}, haplotype_b::AbstractVector{<:Integer}; window::Int=100_000)
    _ = genotypes
    length(positions) == length(haplotype_a) == length(haplotype_b) || throw(DimensionMismatch("positions and haplotypes must have equal length"))
    p = length(positions)
    raw = fill(NaN, p)

    for i in 1:p
        idx = findall(abs.(Int.(positions) .- Int(positions[i])) .<= window)
        isempty(idx) && continue
        ehh_a = mean(haplotype_a[idx] .== haplotype_a[i])
        ehh_b = mean(haplotype_b[idx] .== haplotype_b[i])
        raw[i] = log((ehh_a + eps(Float64)) / (ehh_b + eps(Float64)))
    end

    valid = isfinite.(raw)
    sum(valid) < 3 && return raw
    mu = mean(raw[valid])
    sd = max(std(raw[valid]), eps(Float64))
    raw[valid] = (raw[valid] .- mu) ./ sd
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, raw, "ihs_score")
end

"""
    xp_ehh_score(genotypes, positions, pops_a, pops_b; window=200000, min_maf=0.01)

Compute a cross-population EHH-like score based on local homozygosity contrast.
"""
function xp_ehh_score(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, positions::AbstractVector{<:Integer}, pops_a::AbstractVector{Int}, pops_b::AbstractVector{Int}; window::Int=200_000, min_maf::Real=0.01)
    G = _matrix(genotypes)
    n, p = size(G)
    length(positions) == p || throw(DimensionMismatch("positions length must match number of variants"))
    all(index -> 1 <= index <= n, pops_a) || throw(BoundsError("population A index out of bounds"))
    all(index -> 1 <= index <= n, pops_b) || throw(BoundsError("population B index out of bounds"))

    mask = calculate_maf(G; min_maf=min_maf, max_maf=1.0)
    scores = fill(NaN, p)

    for i in 1:p
        mask[i] || continue
        idx = findall(abs.(Int.(positions) .- Int(positions[i])) .<= window)
        isempty(idx) && continue

        GA = Matrix{Float64}(G[pops_a, idx])
        GB = Matrix{Float64}(G[pops_b, idx])
        _impute_missing_by_mean!(GA)
        _impute_missing_by_mean!(GB)

        pA = vec(mean(GA, dims=1) ./ 2.0)
        pB = vec(mean(GB, dims=1) ./ 2.0)
        ehh_a = mean(1.0 .- 2.0 .* pA .* (1.0 .- pA))
        ehh_b = mean(1.0 .- 2.0 .* pB .* (1.0 .- pB))
        scores[i] = log((ehh_a + eps(Float64)) / (ehh_b + eps(Float64)))
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, scores, "xp_ehh_score")
end

"""
    fst_outlier_test(genotypes, positions; window=500000, top_frac=0.01, groups=nothing)

Estimate per-variant Fst and report top outlier loci.
"""
function fst_outlier_test(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, positions::AbstractVector{<:Integer}; window::Int=500_000, top_frac::Real=0.01, groups::Union{Nothing,AbstractVector{<:Integer}}=nothing, min_maf::Real=0.01)
    _ = window
    G = _matrix(genotypes)
    n, p = size(G)
    length(positions) == p || throw(DimensionMismatch("positions length must match number of variants"))

    grp = if groups === nothing
        vcat(fill(1, n ÷ 2), fill(2, n - n ÷ 2))
    else
        length(groups) == n || throw(DimensionMismatch("groups length must match number of samples"))
        Int.(groups)
    end
    unique_groups = unique(grp)
    length(unique_groups) >= 2 || throw(ArgumentError("at least two groups are required for Fst"))
    g1 = findall(==(unique_groups[1]), grp)
    g2 = findall(==(unique_groups[2]), grp)

    fst = fill(NaN, p)
    maf_mask = calculate_maf(G; min_maf=min_maf, max_maf=1.0)
    for j in 1:p
        maf_mask[j] || continue
        x1 = @view G[g1, j]
        x2 = @view G[g2, j]
        v1 = x1[isfinite.(x1)]
        v2 = x2[isfinite.(x2)]
        (isempty(v1) || isempty(v2)) && continue

        p1 = clamp(mean(v1) / 2.0, 0.0, 1.0)
        p2 = clamp(mean(v2) / 2.0, 0.0, 1.0)
        pbar = (p1 + p2) / 2.0
        ht = 2.0 * pbar * (1.0 - pbar)
        hs = (2.0 * p1 * (1.0 - p1) + 2.0 * p2 * (1.0 - p2)) / 2.0
        fst[j] = ht <= eps(Float64) ? 0.0 : max((ht - hs) / ht, 0.0)
    end

    valid = fst[isfinite.(fst)]
    isempty(valid) && return (outlier_positions=Int[], fst_values=fst, threshold=NaN, outlier_indices=Int[])
    threshold = quantile(valid, clamp(1.0 - top_frac, 0.0, 1.0))
    outlier_idx = [index for index in eachindex(fst) if isfinite(fst[index]) && fst[index] >= threshold]
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (outlier_positions=Int.(positions[outlier_idx]), fst_values=fst, threshold=threshold, outlier_indices=outlier_idx), "fst_outlier_test")
end

function _rank_with_ties(values::AbstractVector{<:Real}; ties::String="average")
    vals = Float64.(values)
    order = sortperm(vals)
    ranks = zeros(Float64, length(vals))
    i = 1
    while i <= length(order)
        j = i
        while j < length(order) && vals[order[j + 1]] == vals[order[i]]
            j += 1
        end
        rank_value = if ties == "average"
            (i + j) / 2.0
        elseif ties == "min"
            Float64(i)
        elseif ties == "max"
            Float64(j)
        else
            throw(ArgumentError("unsupported ties mode '$ties'"))
        end
        for k in i:j
            ranks[order[k]] = rank_value
        end
        i = j + 1
    end
    return ranks
end

"""
    rank_inverse_normal(phenotype; ties="average")

Rank-based inverse normal transform for quantitative traits.
"""
function rank_inverse_normal(phenotype::AbstractVector{<:Real}; ties::String="average")
    y = Float64.(phenotype)
    valid = isfinite.(y)
    n_valid = sum(valid)
    n_valid > 0 || return fill(NaN, length(y))

    ranks = _rank_with_ties(y[valid]; ties=ties)
    q = (ranks .- 0.375) ./ (n_valid + 0.25)
    q = clamp.(q, eps(Float64), 1.0 - eps(Float64))
    transformed = quantile.(Ref(Normal()), q)

    out = fill(NaN, length(y))
    out[valid] = transformed
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, out, "rank_inverse_normal")
end

"""
    rank_inverse_normal!(phenotype; ties="average")

In-place rank-based inverse normal transform for Float64 vectors.
"""
function rank_inverse_normal!(phenotype::AbstractVector{Float64}; ties::String="average")
    phenotype .= rank_inverse_normal(phenotype; ties=ties)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, phenotype, "rank_inverse_normal!")
end

function rank_inverse_normal!(phenotype::AbstractVector{<:Real}; ties::String="average")
    transformed = rank_inverse_normal(phenotype; ties=ties)
    if eltype(phenotype) <: AbstractFloat
        try
            phenotype .= transformed
            return phenotype
        catch
        end
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, transformed, "rank_inverse_normal!")
end

"""
    calculate_sex_check(genotypes)

Compute per-sample X/Y heterozygosity summary when chromosome labels are available.
"""
function calculate_sex_check(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}})
    if !(genotypes isa GenotypeMatrix)
        n = size(_matrix(genotypes), 1)
        return DataFrame(SAMPLE=["sample_$(index)" for index in 1:n], X_HET=fill(NaN, n), Y_HET=fill(NaN, n), INFERRED_SEX=fill("UNKNOWN", n), STATUS=fill("UNKNOWN", n))
    end

    chr = lowercase.(String.(genotypes.bim.chromosome))
    xmask = map(value -> value in ("x", "chrx", "23"), chr)
    ymask = map(value -> value in ("y", "chry", "24"), chr)

    G = _matrix(genotypes)
    n = size(G, 1)
    sample_ids = String.(genotypes.fam.sample_id)
    x_het = fill(NaN, n)
    y_het = fill(NaN, n)
    inferred = fill("UNKNOWN", n)
    status = fill("UNKNOWN", n)

    for i in 1:n
        if any(xmask)
            x = @view G[i, xmask]
            valid_x = isfinite.(x)
            if sum(valid_x) > 0
                x_het[i] = sum(x[valid_x] .== 1.0) / sum(valid_x)
            end
        end
        if any(ymask)
            y = @view G[i, ymask]
            valid_y = isfinite.(y)
            if sum(valid_y) > 0
                y_het[i] = sum(y[valid_y] .== 1.0) / sum(valid_y)
            end
        end

        inferred[i] = isfinite(x_het[i]) ? (x_het[i] < 0.10 ? "XY" : "XX") : "UNKNOWN"
        status[i] = if isfinite(y_het[i]) && inferred[i] == "XX" && y_het[i] > 0.20
            "FAIL"
        else
            "PASS"
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, DataFrame(SAMPLE=sample_ids, X_HET=x_het, Y_HET=y_het, INFERRED_SEX=inferred, STATUS=status), "calculate_sex_check")
end

"""
    calculate_heterozygosity_outliers(genotypes; threshold=3.0)

Identify samples with heterozygosity rates far from the cohort mean.
"""
function calculate_heterozygosity_outliers(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; threshold::Real=3.0)
    G = _matrix(genotypes)
    n = size(G, 1)
    het = fill(NaN, n)

    for i in 1:n
        row = @view G[i, :]
        valid = isfinite.(row)
        if sum(valid) == 0
            continue
        end
        het[i] = sum(row[valid] .== 1.0) / sum(valid)
    end

    finite_mask = isfinite.(het)
    if sum(finite_mask) == 0
        return (outlier=fill(false, n), heterozygosity=het, zscore=fill(NaN, n))
    end

    mu = mean(het[finite_mask])
    sigma = max(std(het[finite_mask]), eps(Float64))
    z = (het .- mu) ./ sigma
    outlier = map(value -> isfinite(value) && abs(value) > threshold, z)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (outlier=outlier, heterozygosity=het, zscore=z), "calculate_heterozygosity_outliers")
end

"""
    sample_qc_report(genotypes; missing_threshold=0.05, het_threshold=3.0, kinship=nothing, related_threshold=0.125)

Build a per-sample QC report combining missingness, heterozygosity, and optional relatedness flags.
"""
function sample_qc_report(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; missing_threshold::Real=0.05, het_threshold::Real=3.0, kinship::Union{Nothing,AbstractMatrix{<:Real}}=nothing, related_threshold::Real=0.125)
    G = _matrix(genotypes)
    n = size(G, 1)
    sample_ids = genotypes isa GenotypeMatrix ? String.(genotypes.fam.sample_id) : ["sample_$(index)" for index in 1:n]

    missing_rate, call_rate = calculate_missingness(G; by=:sample, return_call_rate=true)
    het = calculate_heterozygosity_outliers(G; threshold=het_threshold)
    sex_df = calculate_sex_check(genotypes)

    related_count = zeros(Int, n)
    if kinship !== nothing
        pairs = detect_related_pairs(kinship; threshold=related_threshold)
        for (i, j, _) in pairs
            1 <= i <= n && (related_count[i] += 1)
            1 <= j <= n && (related_count[j] += 1)
        end
    end

    status = fill("PASS", n)
    for i in 1:n
        if !isfinite(missing_rate[i]) || missing_rate[i] > missing_threshold || het.outlier[i] || (i <= nrow(sex_df) && sex_df.STATUS[i] == "FAIL")
            status[i] = "FAIL"
        end
    end
    return DataFrame(
        SAMPLE = sample_ids,
        CALL_RATE = call_rate,
        MISSING_RATE = missing_rate,
        HETEROZYGOSITY = het.heterozygosity,
        HET_Z = het.zscore,
        HET_OUTLIER = het.outlier,
        X_HET = nrow(sex_df) == n ? sex_df.X_HET : fill(NaN, n),
        Y_HET = nrow(sex_df) == n ? sex_df.Y_HET : fill(NaN, n),
        INFERRED_SEX = nrow(sex_df) == n ? sex_df.INFERRED_SEX : fill("UNKNOWN", n),
        RELATED_COUNT = related_count,
        STATUS = status)
end

@inline function _allele_complement(base::AbstractString)
    b = uppercase(String(base))
    b == "A" && return "T"
    b == "T" && return "A"
    b == "C" && return "G"
    b == "G" && return "C"
    return b
end

@inline function _is_palindromic_pair(a1::AbstractString, a2::AbstractString)
    pair = Set((uppercase(String(a1)), uppercase(String(a2))))
    return pair == Set(("A", "T")) || pair == Set(("C", "G"))
end

"""
    harmonise_alleles(result_a, result_b; drop_ambiguous_palindromic=true)

Align two GWAS result sets on effect alleles, flipping signs in `result_b` when
required. Returns harmonized result pairs in a named tuple.
"""
function harmonise_alleles(result_a::GWASResult, result_b::GWASResult; drop_ambiguous_palindromic::Bool=true)
    idx_b = Dict{String,Int}(String(result_b.snp_ids[i]) => i for i in eachindex(result_b.snp_ids))

    keep_a = Int[]
    keep_b = Int[]
    flip_b = Bool[]

    for i in eachindex(result_a.snp_ids)
        snp_id = String(result_a.snp_ids[i])
        j = get(idx_b, snp_id, 0)
        j == 0 && continue

        a1 = uppercase(String(result_a.alleles[i][1]))
        a2 = uppercase(String(result_a.alleles[i][2]))
        b1 = uppercase(String(result_b.alleles[j][1]))
        b2 = uppercase(String(result_b.alleles[j][2]))

        drop_ambiguous_palindromic && (_is_palindromic_pair(a1, a2) || _is_palindromic_pair(b1, b2)) && continue

        b1c = _allele_complement(b1)
        b2c = _allele_complement(b2)

        aligned = false
        flip = false
        if a1 == b1 && a2 == b2
            aligned = true
        elseif a1 == b2 && a2 == b1
            aligned = true
            flip = true
        elseif a1 == b1c && a2 == b2c
            aligned = true
        elseif a1 == b2c && a2 == b1c
            aligned = true
            flip = true
        end
        aligned || continue

        push!(keep_a, i)
        push!(keep_b, j)
        push!(flip_b, flip)
    end

    beta_b = Float64.(result_b.beta[keep_b])
    z_b = Float64.(result_b.zscore[keep_b])
    for k in eachindex(keep_b)
        if flip_b[k]
            beta_b[k] = -beta_b[k]
            z_b[k] = -z_b[k]
        end
    end

    aligned_a = GWASResult(
        result_a.snp_ids[keep_a],
        result_a.chromosomes[keep_a],
        result_a.positions[keep_a],
        result_a.alleles[keep_a],
        result_a.gene_ids[keep_a],
        result_a.beta[keep_a],
        result_a.standard_error[keep_a],
        result_a.zscore[keep_a],
        result_a.pvalue[keep_a],
        result_a.sample_size,
        result_a.covariate_names,
        result_a.phenotype_name,
        string(result_a.method, "+harmonised"))

    aligned_b = GWASResult(
        result_b.snp_ids[keep_b],
        result_b.chromosomes[keep_b],
        result_b.positions[keep_b],
        result_a.alleles[keep_a],
        result_b.gene_ids[keep_b],
        beta_b,
        result_b.standard_error[keep_b],
        z_b,
        result_b.pvalue[keep_b],
        result_b.sample_size,
        result_b.covariate_names,
        result_b.phenotype_name,
        string(result_b.method, "+harmonised"))

    flipped_ids = [String(result_b.snp_ids[keep_b[k]]) for k in eachindex(keep_b) if flip_b[k]]
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (; result_a=aligned_a, result_b=aligned_b, flipped=flipped_ids), "harmonise_alleles")
end

"""
    flip_alleles(gm, reference_alleles)

Flip dosage encoding for variants where allele2 matches the target reference allele.
"""
function flip_alleles(gm::GenotypeMatrix, reference_alleles::Dict{String,String})
    decoded = copy(gm.decoded)
    bim = copy(gm.bim)
    flipped = falses(size(decoded, 2))

    for j in 1:size(decoded, 2)
        snp_id = String(bim.snp_id[j])
        haskey(reference_alleles, snp_id) || continue
        target_ref = uppercase(strip(reference_alleles[snp_id]))
        a1 = uppercase(String(bim.allele1[j]))
        a2 = uppercase(String(bim.allele2[j]))
        (target_ref == a2 && target_ref != a1) || continue

        @inbounds for i in 1:size(decoded, 1)
            value = decoded[i, j]
            isfinite(value) && (decoded[i, j] = 2.0 - value)
        end
        bim.allele1[j], bim.allele2[j] = bim.allele2[j], bim.allele1[j]
        flipped[j] = true
    end

    out = GenotypeMatrix(UInt8[], decoded, bim, copy(gm.fam), gm.prefix)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (; genotypes=out, flipped=flipped), "flip_alleles")
end

"""
    compute_ld_scores(genotypes; window_kb=1000, min_maf=0.01)

Compute per-SNP LD scores (sum of r^2 within a genomic window).
"""
function compute_ld_scores(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; window_kb::Int=1000, min_maf::Real=0.01)
    G = _matrix(genotypes)
    n_snps = size(G, 2)
    scores = zeros(Float64, n_snps)
    n_snps == 0 && return scores

    maf_mask = calculate_maf(G; min_maf=min_maf, max_maf=1.0)
    idx = findall(maf_mask)
    isempty(idx) && return scores

    chr = genotypes isa GenotypeMatrix ? String.(genotypes.bim.chromosome) : fill("", n_snps)
    pos = genotypes isa GenotypeMatrix ? Int.(genotypes.bim.position) : collect(1:n_snps)
    window_bp = max(window_kb, 0) * 1000

    ordered = sort(idx; by=i -> (normalise_chromosome(chr[i]), pos[i]))

    for ii in eachindex(ordered)
        i = ordered[ii]
        chr_i = normalise_chromosome(chr[i])
        pos_i = pos[i]
        for jj in ii:length(ordered)
            j = ordered[jj]
            chr_j = normalise_chromosome(chr[j])
            if chr_i != chr_j
                break
            end
            delta = pos[j] - pos_i
            delta > window_bp && break
            delta < -window_bp && continue
            _, _, r2 = _pairwise_ld(G, i, j)
            isfinite(r2) || continue
            scores[i] += r2
            i == j || (scores[j] += r2)
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, scores, "compute_ld_scores")
end

"""
    prune_ld(genotypes; r2=0.2, window_kb=250, step_kb=50, return_indices=false)

Greedy LD pruning akin PLINK's indep-pairwise workflow.
"""
function prune_ld(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; r2::Real=0.2, window_kb::Int=250, step_kb::Int=50, return_indices::Bool=false)
    _ = step_kb
    G = _matrix(genotypes)
    n_snps = size(G, 2)
    n_snps == 0 && return return_indices ? Int[] : genotypes

    chr = genotypes isa GenotypeMatrix ? String.(genotypes.bim.chromosome) : fill("", n_snps)
    pos = genotypes isa GenotypeMatrix ? Int.(genotypes.bim.position) : collect(1:n_snps)
    order = sortperm(1:n_snps; by=i -> (normalise_chromosome(chr[i]), pos[i]))
    keep = trues(n_snps)
    window_bp = max(window_kb, 0) * 1000

    for ii in eachindex(order)
        i = order[ii]
        keep[i] || continue
        for jj in (ii + 1):length(order)
            j = order[jj]
            normalise_chromosome(chr[i]) == normalise_chromosome(chr[j]) || continue
            if abs(pos[j] - pos[i]) > window_bp
                pos[j] > pos[i] && break
                continue
            end
            _, _, pair_r2 = _pairwise_ld(G, i, j)
            if isfinite(pair_r2) && pair_r2 >= r2
                keep[j] = false
            end
        end
    end

    keep_idx = findall(keep)
    if return_indices
        return keep_idx
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, genotypes isa GenotypeMatrix ? filter_variants(genotypes, keep_idx) : G[:, keep_idx], "prune_ld")
end

"""
    sample_prune_relatedness(kinship; threshold=0.125)

Greedy sample pruning so that no retained pair exceeds the kinship threshold.
"""
function sample_prune_relatedness(kinship::AbstractMatrix{<:Real}; threshold::Real=0.125)
    n = size(kinship, 1)
    size(kinship, 2) == n || throw(DimensionMismatch("kinship must be square"))

    adj = falses(n, n)
    degree = zeros(Int, n)
    for i in 1:n
        for j in (i + 1):n
            value = Float64(kinship[i, j])
            if isfinite(value) && value >= threshold
                adj[i, j] = true
                adj[j, i] = true
                degree[i] += 1
                degree[j] += 1
            end
        end
    end

    removed = falses(n)
    while true
        max_degree, idx = findmax(degree)
        max_degree == 0 && break
        removed[idx] = true
        for j in 1:n
            if adj[idx, j]
                adj[idx, j] = false
                adj[j, idx] = false
                degree[j] = max(degree[j] - 1, 0)
            end
        end
        degree[idx] = 0
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (; keep_indices=findall(.!removed), removed_indices=findall(removed)), "sample_prune_relatedness")
end

"""
    ancestry_inference(pcs, reference_pcs, labels; k=5)

Assign ancestry labels by k-NN projection in PCA space.
"""
function ancestry_inference(pcs::AbstractMatrix{<:Real}, reference_pcs::AbstractMatrix{<:Real}, labels::AbstractVector{<:AbstractString}; k::Int=5)
    size(pcs, 2) == size(reference_pcs, 2) || throw(DimensionMismatch("pcs and reference_pcs must have the same component count"))
    size(reference_pcs, 1) == length(labels) || throw(DimensionMismatch("labels must match reference sample count"))
    n = size(pcs, 1)
    nref = size(reference_pcs, 1)
    nref > 0 || throw(ArgumentError("reference_pcs cannot be empty"))

    k_eff = clamp(k, 1, nref)
    ref = Matrix{Float64}(reference_pcs)
    qry = Matrix{Float64}(pcs)
    lbl = String.(labels)

    assigned = Vector{String}(undef, n)
    confidence = zeros(Float64, n)

    for i in 1:n
        d2 = [sum((qry[i, :] .- ref[j, :]) .^ 2) for j in 1:nref]
        nn = partialsortperm(d2, 1:k_eff)
        votes = Dict{String,Int}()
        for j in nn
            votes[lbl[j]] = get(votes, lbl[j], 0) + 1
        end
        top_label, top_count = first(sort!(collect(votes); by=x -> x[2], rev=true))
        assigned[i] = top_label
        confidence[i] = top_count / k_eff
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, DataFrame(sample=collect(1:n), ancestry=assigned, confidence=confidence), "ancestry_inference")
end

"""
    manhattan_data(result; threshold=5e-8, annotate_top=10)

Return Manhattan-plot-ready coordinates and top hit annotations.
"""
function manhattan_data(result::GWASResult; threshold::Real=5e-8, annotate_top::Int=10)
    n = length(result.snp_ids)
    n == 0 && return (; data=DataFrame(), threshold_line=-log10(clamp(Float64(threshold), eps(Float64), 1.0)), top_hits=DataFrame())

    chrom = normalise_chromosome.(String.(result.chromosomes))
    pos = Int.(result.positions)
    pvals = clamp.(Float64.(result.pvalue), eps(Float64), 1.0)
    neglog10 = -log10.(pvals)

    chrom_order = sort(unique(chrom))
    chrom_max = Dict(ch => maximum(pos[findall(==(ch), chrom)]) for ch in chrom_order)
    offsets = Dict{String,Int}()
    running = 0
    for ch in chrom_order
        offsets[ch] = running
        running += chrom_max[ch]
    end

    cumulative = [offsets[chrom[i]] + pos[i] for i in 1:n]
    data = DataFrame(
        snp_id = String.(result.snp_ids),
        chromosome = chrom,
        position = pos,
        cumulative_position = cumulative,
        pvalue = pvals,
        logp = neglog10,
        significant = pvals .<= threshold)
    top = data[sortperm(data.pvalue)[1:min(annotate_top, n)], :]
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (; data=data, threshold_line=-log10(clamp(Float64(threshold), eps(Float64), 1.0)), top_hits=top), "manhattan_data")
end

"""
    manhattan_plot(result; threshold=5e-8, annotate_top=10)

Back-compat alias for `manhattan_data`. For rendered figures, prefer
`BioToolkit.BioPlotting.manhattan_plot`.
"""
function manhattan_plot(result::GWASResult; threshold::Real=5e-8, annotate_top::Int=10)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, manhattan_data(result; threshold=threshold, annotate_top=annotate_top), "manhattan_plot")
end

"""
    qq_data(pvalues; lambda=NaN)

Return QQ-plot-ready observed/expected coordinates and genomic inflation lambda.
"""
function qq_data(pvalues::AbstractVector{<:Real}; lambda::Real=NaN)
    p = sort(clamp.(Float64.(pvalues), eps(Float64), 1.0))
    n = length(p)
    expected = n == 0 ? Float64[] : [-log10((i - 0.5) / n) for i in 1:n]
    observed = n == 0 ? Float64[] : -log10.(p)

    lambda_gc = if isfinite(lambda)
        Float64(lambda)
    elseif n == 0
        NaN
    else
        chisq = quantile.(Ref(Chisq(1)), 1 .- p)
        median(chisq) / 0.4549364
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (; data=DataFrame(expected=expected, observed=observed, pvalue=p), lambda_gc=lambda_gc), "qq_data")
end

"""
    qq_plot(pvalues; lambda=NaN)

Back-compat alias for `qq_data`. For rendered figures, prefer
`BioToolkit.BioPlotting.qq_plot`.
"""
function qq_plot(pvalues::AbstractVector{<:Real}; lambda::Real=NaN)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, qq_data(pvalues; lambda=lambda), "qq_plot")
end

"""
    locus_zoom(result, chr, start, stop)

Return regional association table for a genomic locus.
"""
function locus_zoom(result::GWASResult, chr, start::Integer, stop::Integer)
    chr_norm = normalise_chromosome(chr)
    left = Int(min(start, stop))
    right = Int(max(start, stop))

    idx = [i for i in eachindex(result.snp_ids) if normalise_chromosome(result.chromosomes[i]) == chr_norm && left <= result.positions[i] <= right]
    data = DataFrame(
        snp_id = String.(result.snp_ids[idx]),
        chromosome = String.(result.chromosomes[idx]),
        position = Int.(result.positions[idx]),
        beta = Float64.(result.beta[idx]),
        se = Float64.(result.standard_error[idx]),
        pvalue = Float64.(result.pvalue[idx]),
        logp = -log10.(clamp.(Float64.(result.pvalue[idx]), eps(Float64), 1.0)))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, data, "locus_zoom")
end

"""
    save_gwas_result(path, result)

Persist a GWAS result to Arrow for fast reload.
"""
function save_gwas_result(path::String, result::GWASResult)
    table = DataFrame(result)
    table.covariate_names = fill(join(String.(result.covariate_names), ';'), nrow(table))
    Arrow.write(path, table)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, path, "save_gwas_result")
end

"""
    load_gwas_result(path)

Load a GWAS result saved by `save_gwas_result`.
"""
function load_gwas_result(path::String)
    table = DataFrame(Arrow.Table(path))
    isempty(table) && return GWASResult(String[], String[], Int[], Tuple{String,String}[], String[], Float64[], Float64[], Float64[], Float64[], 0, String[], "", "")

    covariates = hasproperty(table, :covariate_names) ? String.(split(String(table.covariate_names[1]), ';')) : String[]
    alleles = [(String(table.allele1[i]), String(table.allele2[i])) for i in 1:nrow(table)]
    sample_size = hasproperty(table, :sample_size) ? Int(table.sample_size[1]) : 0
    phenotype_name = hasproperty(table, :phenotype_name) ? String(table.phenotype_name[1]) : "phenotype"
    method = hasproperty(table, :method) ? String(table.method[1]) : "loaded"
    return GWASResult(
        String.(table.snp_id),
        String.(table.chromosome),
        Int.(table.position),
        alleles,
        hasproperty(table, :gene_id) ? String.(table.gene_id) : fill("", nrow(table)),
        Float64.(table.beta),
        Float64.(table.standard_error),
        Float64.(table.zscore),
        Float64.(table.pvalue),
        sample_size,
        covariates,
        phenotype_name,
        method)
end

function _coerce_gene_sets(gene_sets)
    if gene_sets isa Dict
        return [(String(key), Int.(value)) for (key, value) in pairs(gene_sets)]
    elseif gene_sets isa AbstractVector
        return [(String(name), Int.(idx)) for (name, idx) in gene_sets]
    end
    throw(ArgumentError("gene_sets must be a Dict or a vector of (name, indices) tuples"))
end

"""
    score_test_linear(genotypes, phenotype; covariates=nothing, phenotype_name="phenotype")

Compute per-SNP Rao score tests for quantitative traits.
"""
function score_test_linear(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, phenotype::AbstractVector{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, phenotype_name::String="phenotype", multi_thread::Bool=true)
    G = _matrix(genotypes)
    _validate_gwas_inputs(G; phenotype=phenotype, covariates=covariates)
    nobs, nsnps = size(G)

    y_all = Float64.(phenotype)
    design, covariate_names = _covariate_matrix(covariates, nobs)
    X_all = Matrix{Float64}(design)
    row_valid = vec(all(isfinite.(X_all), dims=2)) .& isfinite.(y_all)
    n_valid = sum(row_valid)
    n_cov = size(X_all, 2)
    n_valid > n_cov || throw(ArgumentError("not enough finite samples for score_test_linear"))

    X = X_all[row_valid, :]
    y = y_all[row_valid]
    qrX = qr(X)

    y_res = copy(y)
    coeff_y = qrX \ y_res
    mul!(y_res, X, coeff_y, -1.0, 1.0)
    dof = max(n_valid - n_cov, 1)
    sigma2 = max(dot(y_res, y_res) / dof, eps(Float64))
    t_dist = TDist(dof)

    beta = zeros(Float64, nsnps)
    se = fill(Inf, nsnps)
    z = zeros(Float64, nsnps)
    p = ones(Float64, nsnps)

    @inline function _scan_snp!(snp::Int)
        g = Float64.(G[row_valid, snp])
        finite_mask = isfinite.(g)
        sum(finite_mask) > 0 || return nothing
        mean_g = mean(g[finite_mask])
        g[.!finite_mask] .= mean_g

        coeff_g = qrX \ g
        mul!(g, X, coeff_g, -1.0, 1.0)

        v = dot(g, g)
        v <= eps(Float64) && return nothing
        u = dot(g, y_res)
        b = u / v
        stderr = sqrt(sigma2 / v)
        zscore = u / sqrt(max(sigma2 * v, eps(Float64)))
        pvalue = 2.0 * ccdf(t_dist, abs(zscore))

        beta[snp] = b
        se[snp] = stderr
        z[snp] = zscore
        p[snp] = pvalue
        return nothing
    end

    if _thread_enabled(nsnps, multi_thread)
        @threads for snp in 1:nsnps
            _scan_snp!(snp)
        end
    else
        for snp in 1:nsnps
            _scan_snp!(snp)
        end
    end

    method = "score_test_linear"
    return genotypes isa GenotypeMatrix ?
           _result_from_matrix(genotypes, beta, se, z, p; method=method, phenotype_name=phenotype_name, covariate_names=covariate_names) :
           _result_from_matrix(G, beta, se, z, p; method=method, phenotype_name=phenotype_name, covariate_names=covariate_names)
end

"""
    likelihood_ratio_test(genotypes, phenotype; covariates=nothing, phenotype_name="phenotype")

Per-SNP linear-model LRT against a covariate-only null model.
"""
function likelihood_ratio_test(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, phenotype::AbstractVector{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, phenotype_name::String="phenotype", ridge::Real=1e-8)
    G = _matrix(genotypes)
    _validate_gwas_inputs(G; phenotype=phenotype, covariates=covariates)
    nobs, nsnps = size(G)

    y_all = Float64.(phenotype)
    design, covariate_names = _covariate_matrix(covariates, nobs)
    X0_all = Matrix{Float64}(design)
    row_valid = vec(all(isfinite.(X0_all), dims=2)) .& isfinite.(y_all)
    n_valid = sum(row_valid)
    n_cov = size(X0_all, 2)
    n_valid > n_cov || throw(ArgumentError("not enough finite samples for likelihood_ratio_test"))

    X0 = X0_all[row_valid, :]
    y = y_all[row_valid]
    Gv = G[row_valid, :]

    beta0 = (X0' * X0 + Float64(ridge) * I) \ (X0' * y)
    rss0 = max(sum(abs2, y - X0 * beta0), eps(Float64))

    beta = zeros(Float64, nsnps)
    se = fill(Inf, nsnps)
    z = zeros(Float64, nsnps)
    p = ones(Float64, nsnps)

    for snp in 1:nsnps
        g = Float64.(Gv[:, snp])
        finite_mask = isfinite.(g)
        sum(finite_mask) > 0 || continue
        mean_g = mean(g[finite_mask])
        g[.!finite_mask] .= mean_g

        X1 = hcat(X0, g)
        XtX = Symmetric(X1' * X1 + Float64(ridge) * I)
        beta1 = try
            XtX \ (X1' * y)
        catch
            continue
        end
        resid1 = y - X1 * beta1
        rss1 = max(sum(abs2, resid1), eps(Float64))
        lrt = max(n_valid * log(rss0 / rss1), 0.0)
        pvalue = ccdf(Chisq(1), lrt)

        inv_xtx = try
            inv(XtX)
        catch
            Matrix{Float64}(I, size(XtX, 1), size(XtX, 2))
        end
        sigma2 = rss1 / max(n_valid - size(X1, 2), 1)
        stderr = sqrt(max(sigma2 * inv_xtx[end, end], eps(Float64)))

        beta[snp] = beta1[end]
        se[snp] = stderr
        z[snp] = sign(beta1[end]) * sqrt(lrt)
        p[snp] = pvalue
    end

    method = "likelihood_ratio_test"
    return genotypes isa GenotypeMatrix ?
           _result_from_matrix(genotypes, beta, se, z, p; method=method, phenotype_name=phenotype_name, covariate_names=covariate_names) :
           _result_from_matrix(G, beta, se, z, p; method=method, phenotype_name=phenotype_name, covariate_names=covariate_names)
end

function _cox_group_boundaries(times_sorted::AbstractVector{<:Real}, events_sorted::AbstractVector{Bool})
    starts = Int[]
    ends = Int[]
    event_counts = Int[]
    n = length(times_sorted)
    i = 1
    while i <= n
        j = i
        t = times_sorted[i]
        while j < n && times_sorted[j + 1] == t
            j += 1
        end
        push!(starts, i)
        push!(ends, j)
        push!(event_counts, count(events_sorted[i:j]))
        i = j + 1
    end
    return starts, ends, event_counts
end

function _cox_fit_null(covariates::AbstractMatrix{<:Real}, times::AbstractVector{<:Real}, events::AbstractVector{Bool}; max_iter::Int=60, tol::Real=1e-8, ridge::Real=1e-6)
    n = length(times)
    size(covariates, 1) == n || throw(DimensionMismatch("covariates and times must have matching lengths"))

    order = sortperm(times)
    t_sorted = Float64.(times[order])
    e_sorted = events[order]
    X_sorted = Matrix{Float64}(covariates[order, :])
    p = size(X_sorted, 2)
    starts, ends, event_counts = _cox_group_boundaries(t_sorted, e_sorted)

    if p == 0
        return (beta=zeros(Float64, 0), order=order, times=t_sorted, events=e_sorted, weights=ones(Float64, n), starts=starts, ends=ends, event_counts=event_counts)
    end

    beta = zeros(Float64, p)
    converged = false

    for _ in 1:max_iter
        eta = clamp.(X_sorted * beta, -30.0, 30.0)
        w = exp.(eta)

        score = zeros(Float64, p)
        info = zeros(Float64, p, p)

        S0_at = zeros(Float64, length(starts))
        S1_at = zeros(Float64, length(starts), p)
        S2_at = [zeros(Float64, p, p) for _ in eachindex(starts)]

        S0_cur = 0.0
        S1_cur = zeros(Float64, p)
        S2_cur = zeros(Float64, p, p)
        g = length(starts)

        # Ascending-time order uses suffix risk sets: R(t) = { i : t_i >= t }.
        for idx in n:-1:1
            wi = w[idx]
            xi = @view X_sorted[idx, :]
            S0_cur += wi
            @inbounds for a in 1:p
                xa = xi[a]
                S1_cur[a] += wi * xa
                for b in 1:p
                    S2_cur[a, b] += wi * xa * xi[b]
                end
            end

            if g >= 1 && idx == starts[g]
                S0_at[g] = S0_cur
                @inbounds for a in 1:p
                    S1_at[g, a] = S1_cur[a]
                end
                S2_at[g] = copy(S2_cur)
                g -= 1
            end
        end

        x_event = zeros(Float64, p)
        for g in eachindex(starts)
            left = starts[g]
            right = ends[g]
            d = event_counts[g]
            d == 0 && continue

            fill!(x_event, 0.0)
            for idx in left:right
                e_sorted[idx] || continue
                @inbounds for a in 1:p
                    x_event[a] += X_sorted[idx, a]
                end
            end

            denom = max(S0_at[g], eps(Float64))
            S1 = @view S1_at[g, :]
            S2 = S2_at[g]
            mu = S1 ./ denom
            score .+= x_event .- d .* mu
            info .+= d .* (S2 ./ denom .- (S1 * S1') ./ max(denom^2, eps(Float64)))
        end

        for a in 1:p
            info[a, a] += ridge
        end

        step = try
            info \ score
        catch
            zeros(Float64, p)
        end
        beta_new = beta + step
        if maximum(abs.(step)) <= tol * (1.0 + maximum(abs.(beta_new)))
            beta = beta_new
            converged = true
            break
        end
        beta = beta_new
    end

    eta = p == 0 ? zeros(Float64, n) : clamp.(X_sorted * beta, -30.0, 30.0)
    w = exp.(eta)
    converged || @debug("cox null fit reached max iterations")

    return (beta=beta, order=order, times=t_sorted, events=e_sorted, weights=w, starts=starts, ends=ends, event_counts=event_counts)
end

"""
    gwas_survival_scan(genotypes, time, event; covariates=nothing, phenotype_name="survival")

Cox proportional hazards score scan (Breslow ties) with optional covariate adjustment.
"""
function gwas_survival_scan(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, time::AbstractVector{<:Real}, event::AbstractVector{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, phenotype_name::String="survival", multi_thread::Bool=true)
    length(time) == length(event) || throw(DimensionMismatch("time and event vectors must have equal length"))
    G = _matrix(genotypes)
    _validate_gwas_inputs(G)
    nobs, nsnps = size(G)
    length(time) == nobs || throw(DimensionMismatch("time length must match sample count"))

    times_all = Float64.(time)
    event_all = Float64.(event)
    cov_matrix = covariates === nothing ? zeros(Float64, nobs, 0) : Matrix{Float64}(covariates)
    size(cov_matrix, 1) == nobs || throw(DimensionMismatch("covariates must have one row per sample"))

    cov_valid = size(cov_matrix, 2) == 0 ? trues(nobs) : vec(all(isfinite.(cov_matrix), dims=2))
    row_valid = isfinite.(times_all) .& (times_all .> 0.0) .& isfinite.(event_all) .& cov_valid
    n_valid = sum(row_valid)
    n_valid > size(cov_matrix, 2) || throw(ArgumentError("not enough valid samples for survival scan"))

    times = times_all[row_valid]
    events = Bool.(event_all[row_valid] .> 0.0)
    count(events) > 0 || throw(ArgumentError("survival scan requires at least one event"))
    cov_valid_matrix = cov_matrix[row_valid, :]

    null_fit = _cox_fit_null(cov_valid_matrix, times, events)
    order = null_fit.order
    w_sorted = null_fit.weights
    e_sorted = null_fit.events
    starts = null_fit.starts
    ends = null_fit.ends
    event_counts = null_fit.event_counts
    cw = reverse(cumsum(reverse(w_sorted)))
    e_sorted_float = Float64.(e_sorted)

    beta = zeros(Float64, nsnps)
    se = fill(Inf, nsnps)
    z = zeros(Float64, nsnps)
    p = ones(Float64, nsnps)

    @inline function _scan_one_snp!(snp::Int)
        g_local = Float64.(G[row_valid, snp])
        finite_mask = isfinite.(g_local)
        sum(finite_mask) > 0 || return nothing
        mean_g = mean(g_local[finite_mask])
        g_local[.!finite_mask] .= mean_g

        gs = g_local[order]
        cwg = reverse(cumsum(reverse(w_sorted .* gs)))
        cwg2 = reverse(cumsum(reverse(w_sorted .* gs .* gs)))
        cevent_g = cumsum(gs .* e_sorted_float)

        U = 0.0
        V = 0.0
        for idx in eachindex(starts)
            d = event_counts[idx]
            d == 0 && continue
            left = starts[idx]
            right = ends[idx]
            s0 = max(cw[left], eps(Float64))
            s1 = cwg[left] / s0
            s2 = cwg2[left] / s0
            event_sum = cevent_g[right] - (left > 1 ? cevent_g[left - 1] : 0.0)
            U += event_sum - d * s1
            V += d * max(s2 - s1^2, 0.0)
        end

        V > eps(Float64) || return nothing
        stderr = 1.0 / sqrt(V)
        zscore = U * stderr
        beta[snp] = U / V
        se[snp] = stderr
        z[snp] = zscore
        p[snp] = 2.0 * ccdf(Normal(), abs(zscore))
        return nothing
    end

    if _thread_enabled(nsnps, multi_thread)
        @threads for snp in 1:nsnps
            _scan_one_snp!(snp)
        end
    else
        for snp in 1:nsnps
            _scan_one_snp!(snp)
        end
    end

    covariate_names = size(cov_valid_matrix, 2) == 0 ? String[] : ["cov_$(i)" for i in 1:size(cov_valid_matrix, 2)]
    base = genotypes isa GenotypeMatrix ?
           _result_from_matrix(genotypes, beta, se, z, p; method="cox_score_scan", phenotype_name=phenotype_name, covariate_names=covariate_names) :
           _result_from_matrix(G, beta, se, z, p; method="cox_score_scan", phenotype_name=phenotype_name, covariate_names=covariate_names)

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GWASResult(base.snp_ids, base.chromosomes, base.positions, base.alleles, base.gene_ids, base.beta, base.standard_error, base.zscore, base.pvalue, n_valid, base.covariate_names, base.phenotype_name, "cox_score_scan"), "gwas_survival_scan")
end

"""
    gwas_ordinal_scan(genotypes, phenotype; covariates=nothing, phenotype_name="ordinal")

Approximate proportional-odds GWAS by rank-transforming the ordinal phenotype.
"""
function gwas_ordinal_scan(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, phenotype::AbstractVector{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, phenotype_name::String="ordinal", multi_thread::Bool=true)
    y = Float64.(phenotype)
    valid = isfinite.(y)
    y_rank = fill(NaN, length(y))
    y_rank[valid] = _rank_with_ties(collect(y[valid]); ties="average")

    scan = gwas_linear_scan(genotypes, y_rank; covariates=covariates, phenotype_name=phenotype_name, multi_thread=multi_thread)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GWASResult(scan.snp_ids, scan.chromosomes, scan.positions, scan.alleles, scan.gene_ids, scan.beta, scan.standard_error, scan.zscore, scan.pvalue, scan.sample_size, scan.covariate_names, scan.phenotype_name, "ordinal_scan_approx"), "gwas_ordinal_scan")
end

"""
    gwas_multivariate_scan(genotypes, phenotypes; covariates=nothing, phenotype_name="multivariate")

Combine univariate linear-scan z-scores across multiple phenotypes using Stouffer's method.
"""
function gwas_multivariate_scan(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, phenotypes::AbstractMatrix{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, phenotype_name::String="multivariate", multi_thread::Bool=true)
    size(phenotypes, 2) >= 1 || throw(ArgumentError("phenotypes must contain at least one column"))
    scans = [gwas_linear_scan(genotypes, phenotypes[:, j]; covariates=covariates, phenotype_name="$(phenotype_name)_$(j)", multi_thread=multi_thread) for j in 1:size(phenotypes, 2)]
    base = scans[1]
    nsnps = length(base.snp_ids)

    beta = zeros(Float64, nsnps)
    se = fill(Inf, nsnps)
    z = zeros(Float64, nsnps)
    p = ones(Float64, nsnps)

    for i in 1:nsnps
        zvals = [scan.zscore[i] for scan in scans if isfinite(scan.zscore[i])]
        bvals = [scan.beta[i] for scan in scans if isfinite(scan.beta[i])]
        isempty(zvals) && continue
        z_comb = sum(zvals) / sqrt(length(zvals))
        b_comb = mean(bvals)
        beta[i] = b_comb
        z[i] = z_comb
        se[i] = abs(b_comb) / max(abs(z_comb), eps(Float64))
        p[i] = 2.0 * ccdf(Normal(), abs(z_comb))
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GWASResult(base.snp_ids, base.chromosomes, base.positions, base.alleles, base.gene_ids, beta, se, z, p, base.sample_size, base.covariate_names, phenotype_name, "multivariate_scan"), "gwas_multivariate_scan")
end

"""
    burden_test(genotypes, gene_sets, phenotype; covariates=nothing, weights=nothing)

Gene/set-based burden tests from aggregated variant dosages.
"""
function burden_test(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, gene_sets, phenotype::AbstractVector{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, weights=nothing)
    G = _matrix(genotypes)
    _validate_gwas_inputs(G; phenotype=phenotype, covariates=covariates)
    nobs, nsnps = size(G)
    sets = _coerce_gene_sets(gene_sets)
    design, _ = _covariate_matrix(covariates, nobs)
    X_all = Matrix{Float64}(design)
    y_all = Float64.(phenotype)
    design_valid = vec(all(isfinite.(X_all), dims=2))

    rows = NamedTuple[]
    for (set_name, idx_raw) in sets
        idx = unique([i for i in idx_raw if 1 <= i <= nsnps])
        isempty(idx) && continue

        w = if weights === nothing
            ones(Float64, length(idx))
        elseif weights isa AbstractVector
            length(weights) == length(idx) ? Float64.(weights) : ones(Float64, length(idx))
        elseif weights isa Dict
            [Float64(get(weights, i, 1.0)) for i in idx]
        else
            ones(Float64, length(idx))
        end

        burden = zeros(Float64, nobs)
        for (k, snp) in enumerate(idx)
            col = Float64.(G[:, snp])
            finite_mask = isfinite.(col)
            sum(finite_mask) > 0 || continue
            mean_col = mean(col[finite_mask])
            col[.!finite_mask] .= mean_col
            burden .+= w[k] .* col
        end

        valid = design_valid .& isfinite.(y_all) .& isfinite.(burden)
        sum(valid) > size(X_all, 2) || continue
        X = hcat(X_all[valid, :], burden[valid])
        y = y_all[valid]
        coeff, stderr, zscore, pvalue = _linear_fit_last_coefficient(X, y)
        push!(rows, (set=set_name, n_variants=length(idx), beta=coeff, se=stderr, zscore=zscore, pvalue=pvalue))
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, DataFrame(rows), "burden_test")
end

"""
    skat_test(genotypes, gene_sets, phenotype; covariates=nothing, kernel=:linear, weights=nothing)

Approximate SKAT-style set tests using residualized phenotypes.
"""
function skat_test(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, gene_sets, phenotype::AbstractVector{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, kernel::Symbol=:linear, weights=nothing)
    _ = kernel
    G = _matrix(genotypes)
    _validate_gwas_inputs(G; phenotype=phenotype, covariates=covariates)
    nobs, nsnps = size(G)
    sets = _coerce_gene_sets(gene_sets)

    design, _ = _covariate_matrix(covariates, nobs)
    X_all = Matrix{Float64}(design)
    y_all = Float64.(phenotype)
    valid = vec(all(isfinite.(X_all), dims=2)) .& isfinite.(y_all)
    sum(valid) > size(X_all, 2) || return DataFrame(set=String[], n_variants=Int[], q_stat=Float64[], dof=Float64[], pvalue=Float64[])

    X = X_all[valid, :]
    y = y_all[valid]
    beta0 = (X' * X + 1e-8I) \ (X' * y)
    residual = y - X * beta0
    G_valid = G[valid, :]

    rows = NamedTuple[]
    for (set_name, idx_raw) in sets
        idx = unique([i for i in idx_raw if 1 <= i <= nsnps])
        isempty(idx) && continue

        Gs = Matrix{Float64}(G_valid[:, idx])
        _impute_missing_by_mean!(Gs)
        Gs .-= mean(Gs, dims=1)
        if weights !== nothing && weights isa AbstractVector && length(weights) == size(Gs, 2)
            Gs .*= reshape(Float64.(weights), 1, :)
        end

        eigvals = Float64[]
        q_stat = 0.0
        if size(Gs, 1) > size(Gs, 2)
            # Kernel trick: non-zero eigenvalues of G*G' equal those of G'*G.
            small_kernel = Symmetric(Gs' * Gs)
            eigvals = max.(Float64.(real.(LinearAlgebra.eigvals(small_kernel))), 0.0)
            proj = Gs' * residual
            q_stat = dot(proj, proj)
        else
            full_kernel = Symmetric(Gs * Gs')
            eigvals = max.(Float64.(real.(LinearAlgebra.eigvals(full_kernel))), 0.0)
            q_stat = dot(residual, full_kernel * residual)
        end

        m1 = sum(eigvals)
        m2 = 2.0 * sum(eigvals .^ 2)
        if m1 <= eps(Float64) || m2 <= eps(Float64)
            push!(rows, (set=set_name, n_variants=length(idx), q_stat=q_stat, dof=1.0, pvalue=1.0))
            continue
        end
        scale = m2 / max(2.0 * m1, eps(Float64))
        dof = max(2.0 * m1^2 / max(m2, eps(Float64)), 1.0)
        pvalue = ccdf(Chisq(dof), q_stat / max(scale, eps(Float64)))

        push!(rows, (set=set_name, n_variants=length(idx), q_stat=q_stat, dof=dof, pvalue=pvalue))
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, DataFrame(rows), "skat_test")
end

"""
    conditional_fdr(pvalues, prior_pvalues; pi0=1.0)

Compute a simple conditional FDR estimate using primary and prior p-values.
"""
function conditional_fdr(pvalues::AbstractVector{<:Real}, prior_pvalues::AbstractVector{<:Real}; pi0::Real=1.0)
    length(pvalues) == length(prior_pvalues) || throw(DimensionMismatch("pvalues and prior_pvalues must have equal length"))
    primary = clamp.(Float64.(pvalues), eps(Float64), 1.0)
    prior = clamp.(Float64.(prior_pvalues), eps(Float64), 1.0)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, clamp.(primary .* prior ./ max(Float64(pi0), eps(Float64)), 0.0, 1.0), "conditional_fdr")
end

"""
    storey_pi0_estimate(pvalues; lambda=0.5)

Estimate the null fraction pi0 with Storey's estimator.
"""
function storey_pi0_estimate(pvalues::AbstractVector{<:Real}; lambda=0.5)
    p = clamp.(Float64.(pvalues), 0.0, 1.0)
    isempty(p) && return NaN

    if lambda isa AbstractVector
        lambdas = clamp.(Float64.(lambda), 0.0, 0.95)
        isempty(lambdas) && return NaN
        estimates = [count(>(l), p) / (max(1.0 - l, eps(Float64)) * length(p)) for l in lambdas]
        return clamp(mean(estimates), 0.0, 1.0)
    end

    l = clamp(Float64(lambda), 0.0, 0.95)
    pi0 = count(>(l), p) / (max(1.0 - l, eps(Float64)) * length(p))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, clamp(pi0, 0.0, 1.0), "storey_pi0_estimate")
end

"""
    gwas_power_calculation(n, maf, beta, alpha)

Approximate two-sided power for additive linear GWAS models.
"""
function gwas_power_calculation(n::Real, maf::Real, beta::Real, alpha::Real)
    n > 0 || throw(ArgumentError("n must be positive"))
    p = clamp(Float64(maf), eps(Float64), 1.0 - eps(Float64))
    b = Float64(beta)
    a = clamp(Float64(alpha), eps(Float64), 1.0)

    ncp = Float64(n) * (2.0 * p * (1.0 - p)) * b^2
    zcrit = quantile(Normal(), 1.0 - a / 2.0)
    delta = sqrt(max(ncp, 0.0))
    power = ccdf(Normal(), zcrit - delta) + cdf(Normal(), -zcrit - delta)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, clamp(power, 0.0, 1.0), "gwas_power_calculation")
end

"""
    permutation_pvalue(stat, null_dist; n_perm=length(null_dist))

Compute empirical two-sided permutation p-value.
"""
function permutation_pvalue(stat::Real, null_dist::AbstractVector{<:Real}; n_perm::Int=length(null_dist))
    n = min(max(n_perm, 0), length(null_dist))
    n == 0 && return NaN
    obs = abs(Float64(stat))
    tail = count(i -> abs(Float64(null_dist[i])) >= obs, 1:n)
    return (tail + 1) / (n + 1)
end

"""
    genomic_sem_fit(ldsc_results; model=:common_factor)

Lightweight genomic SEM summary fit from a vector of LDSC-like results.
"""
function genomic_sem_fit(ldsc_results; model=:common_factor)
    h2 = Float64[]
    for item in ldsc_results
        if item isa LDSCResult
            isfinite(item.h2) && push!(h2, item.h2)
        elseif hasproperty(item, :h2)
            value = Float64(getproperty(item, :h2))
            isfinite(value) && push!(h2, value)
        elseif item isa Real
            value = Float64(item)
            isfinite(value) && push!(h2, value)
        end
    end

    isempty(h2) && return (; model=model, h2=Float64[], covariance=zeros(Float64, 0, 0), loadings=Float64[])
    covariance = Diagonal(h2)
    loadings = fill(1.0 / sqrt(length(h2)), length(h2))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (; model=model, h2=h2, covariance=Matrix(covariance), loadings=loadings), "genomic_sem_fit")
end

"""
    dosage_to_hardcall(dosage; threshold=0.1)

Convert dosage values to hard calls (0/1/2), marking uncertain values as NaN.
"""
function dosage_to_hardcall(dosage::AbstractArray{<:Real}; threshold::Real=0.1)
    out = fill(NaN, size(dosage))
    thr = max(Float64(threshold), 0.0)
    for idx in eachindex(dosage)
        value = Float64(dosage[idx])
        isfinite(value) || continue
        hard = round(value)
        abs(value - hard) <= thr || continue
        out[idx] = clamp(hard, 0.0, 2.0)
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, out, "dosage_to_hardcall")
end

"""
    info_score_from_dosage(dosage, gp)

Compute an imputation INFO-like score from dosage and genotype probabilities.
"""
function info_score_from_dosage(dosage::AbstractVector{<:Real}, gp::AbstractMatrix{<:Real})
    length(dosage) == size(gp, 1) || throw(DimensionMismatch("dosage and gp must have matching sample counts"))
    size(gp, 2) >= 3 || throw(DimensionMismatch("gp must have at least 3 probability columns"))

    d = Float64.(dosage)
    expected = Float64.(gp[:, 2]) .+ 2.0 .* Float64.(gp[:, 3])
    valid = isfinite.(d) .& isfinite.(expected)
    sum(valid) > 1 || return NaN

    d_valid = d[valid]
    expected_valid = expected[valid]
    p = clamp(mean(expected_valid) / 2.0, 0.0, 1.0)
    expected_var = 2.0 * p * (1.0 - p)
    expected_var <= eps(Float64) && return NaN

    mse = mean((d_valid .- expected_valid) .^ 2)
    info = 1.0 - mse / expected_var
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, clamp(info, 0.0, 1.0), "info_score_from_dosage")
end

"""
    he_regression_variance_components(y, K_list)

Haseman-Elston (HE) regression estimator for multiple random-effect components.
"""
function he_regression_variance_components(y::AbstractVector{<:Real}, K_list::AbstractVector{<:AbstractMatrix{<:Real}})
    n = length(y)
    n > 1 || throw(ArgumentError("y must contain at least two samples"))
    isempty(K_list) && throw(ArgumentError("K_list cannot be empty"))
    all(size(K, 1) == n == size(K, 2) for K in K_list) || throw(DimensionMismatch("all kinship matrices must be n x n"))

    yv = Float64.(y)
    yv .-= mean(yv)

    components = [Matrix{Float64}(K) for K in K_list]
    push!(components, Matrix{Float64}(I, n, n))
    m = length(components)

    A = zeros(Float64, m, m)
    b = zeros(Float64, m)
    for i in 1:m
        Ki = components[i]
        b[i] = dot(yv, Ki * yv)
        for j in i:m
            value = tr(Ki * components[j])
            A[i, j] = value
            A[j, i] = value
        end
    end

    vc = try
        A \ b
    catch
        fill(0.0, m)
    end
    vc = max.(vc, 0.0)
    total = sum(vc)
    h2 = total <= eps(Float64) ? fill(0.0, m - 1) : vc[1:(m - 1)] ./ total
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (; sigma2_random=vc[1:(m - 1)], sigma2_residual=vc[end], h2=h2), "he_regression_variance_components")
end

"""
    reml_variance_components(y, K_list; covariates=nothing, add_intercept=true, max_iter=100, tol=1e-6, ridge=1e-8, verbose=false)

Estimate variance components using restricted maximum likelihood (REML)
under a linear mixed model with covariance:

`V = sum(σ²ᵢ Kᵢ) + σ²ₑ I`
"""
function reml_variance_components(
    y::AbstractVector{<:Real},
    K_list::AbstractVector{<:AbstractMatrix{<:Real}};
    covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
    add_intercept::Bool=true,
    max_iter::Int=100,
    tol::Real=1e-6,
    ridge::Real=1e-8,
    verbose::Bool=false)
    n = length(y)
    n > 1 || throw(ArgumentError("y must contain at least two samples"))
    isempty(K_list) && throw(ArgumentError("K_list cannot be empty"))
    all(size(K, 1) == n == size(K, 2) for K in K_list) || throw(DimensionMismatch("all kinship matrices must be n x n"))

    yv = Float64.(y)
    isfinite.(yv) |> all || throw(ArgumentError("y must contain only finite values"))

    X = if covariates === nothing
        add_intercept ? ones(Float64, n, 1) : zeros(Float64, n, 0)
    else
        C = Matrix{Float64}(covariates)
        size(C, 1) == n || throw(DimensionMismatch("covariates must have one row per sample"))
        all(isfinite, C) || throw(ArgumentError("covariates must contain only finite values"))
        add_intercept ? hcat(ones(Float64, n, 1), C) : C
    end

    kinship = [Matrix{Float64}(Symmetric((Matrix{Float64}(K) + Matrix{Float64}(K)') / 2)) for K in K_list]
    n_random = length(kinship)
    deriv_mats = vcat(kinship, [Matrix{Float64}(I, n, n)])
    n_components = n_random + 1

    he = he_regression_variance_components(yv, kinship)
    state = vcat(max.(Float64.(he.sigma2_random), eps(Float64)), max(Float64(he.sigma2_residual), eps(Float64)))

    function _eval_state(sigma2::Vector{Float64})
        V = sigma2[end] .* Matrix{Float64}(I, n, n)
        for i in 1:n_random
            V .+= sigma2[i] .* kinship[i]
        end
        V = (V + V') ./ 2

        jitter = max(Float64(ridge), eps(Float64))
        cholV = nothing
        for _ in 1:8
            try
                cholV = cholesky(Symmetric(V + jitter * I))
                break
            catch
                jitter *= 10.0
            end
        end
        cholV === nothing && return nothing

        Vinv = Matrix{Float64}(cholV \ Matrix{Float64}(I, n, n))
        logdetV = 2.0 * sum(log, diag(cholV.U))

        if size(X, 2) == 0
            P = Vinv
            logdetX = 0.0
        else
            VinvX = Vinv * X
            XtVinvX = Symmetric((X' * VinvX + (X' * VinvX)') / 2)

            jitter_x = max(Float64(ridge), eps(Float64))
            cholX = nothing
            for _ in 1:8
                try
                    cholX = cholesky(Symmetric(Matrix{Float64}(XtVinvX) + jitter_x * I))
                    break
                catch
                    jitter_x *= 10.0
                end
            end
            cholX === nothing && return nothing

            Xt_inv = Matrix{Float64}(cholX \ Matrix{Float64}(I, size(X, 2), size(X, 2)))
            P = Vinv - VinvX * Xt_inv * transpose(VinvX)
            logdetX = 2.0 * sum(log, diag(cholX.U))
        end

        P = (P + P') ./ 2
        Py = P * yv
        dof = max(n - size(X, 2), 1)
        loglik = -0.5 * (dof * log(2π) + logdetV + logdetX + dot(yv, Py))

        score = zeros(Float64, n_components)
        PD = Vector{Matrix{Float64}}(undef, n_components)
        for i in 1:n_components
            Di = deriv_mats[i]
            PDi = P * Di
            PD[i] = PDi
            score[i] = -0.5 * (tr(PDi) - dot(Py, Di * Py))
        end

        info = zeros(Float64, n_components, n_components)
        for i in 1:n_components
            for j in i:n_components
                value = 0.5 * tr(PD[i] * PD[j])
                info[i, j] = value
                info[j, i] = value
            end
        end

        return (; loglik=loglik, score=score, info=info)
    end

    evaluation = _eval_state(state)
    if evaluation === nothing
        return (; sigma2_random=Float64.(he.sigma2_random), sigma2_residual=Float64(he.sigma2_residual), h2=Float64.(he.h2), converged=false, iterations=0, loglik=NaN)
    end

    converged = false
    loglik = evaluation.loglik
    iterations = 0

    for iter in 1:max_iter
        iterations = iter
        info_reg = evaluation.info + Float64(ridge) * Matrix{Float64}(I, n_components, n_components)
        delta = try
            info_reg \ evaluation.score
        catch
            zeros(Float64, n_components)
        end
        all(isfinite, delta) || break

        step = 1.0
        best_state = nothing
        best_eval = nothing
        while step > 1e-6
            candidate = state .+ step .* delta
            if any(candidate .<= eps(Float64))
                step *= 0.5
                continue
            end
            candidate_eval = _eval_state(candidate)
            if candidate_eval !== nothing && candidate_eval.loglik >= loglik - 1e-8
                best_state = candidate
                best_eval = candidate_eval
                break
            end
            step *= 0.5
        end

        best_state === nothing && break
        rel_change = maximum(abs.(best_state .- state) ./ (abs.(state) .+ 1e-8))
        state = best_state
        evaluation = best_eval
        loglik = evaluation.loglik

        if verbose && (iter == 1 || iter % 10 == 0)
            @info("reml_variance_components", iteration=iter, loglik=loglik, sigma2=state)
        end

        if rel_change < tol
            converged = true
            break
        end
    end

    total = sum(state)
    h2 = total <= eps(Float64) ? fill(0.0, n_random) : state[1:n_random] ./ total
    return (; sigma2_random=state[1:n_random], sigma2_residual=state[end], h2=h2, converged=converged, iterations=iterations, loglik=loglik)
end

"""
    simulate_genotypes(n, p; maf_dist=Beta(0.8, 0.8), ld_blocks=nothing)

Simulate genotype dosages with optional coarse LD block structure.
"""
function simulate_genotypes(n::Int, p::Int; maf_dist=Beta(0.8, 0.8), ld_blocks::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    (n > 0 && p > 0) || throw(ArgumentError("n and p must be positive"))
    mafs = clamp.(rand(maf_dist, p), 0.01, 0.5)
    G = zeros(Float64, n, p)

    if ld_blocks === nothing
        for j in 1:p
            G[:, j] .= rand(Binomial(2, mafs[j]), n)
        end
        return G
    end

    block_sizes = Int.(ld_blocks)
    sum(block_sizes) == p || throw(DimensionMismatch("sum(ld_blocks) must equal p"))
    offset = 0
    for block_size in block_sizes
        block_idx = (offset + 1):(offset + block_size)
        block_maf = mean(mafs[block_idx])
        latent = rand(Binomial(2, block_maf), n)
        for j in block_idx
            noise = rand(Binomial(2, mafs[j]), n)
            G[:, j] .= clamp.(0.75 .* latent .+ 0.25 .* noise, 0.0, 2.0)
        end
        offset += block_size
    end
    return G
end

"""
    simulate_phenotype(genotypes, causal_idx, beta, h2; binary=false)

Simulate phenotypes from a genotype matrix and a causal architecture.
"""
function simulate_phenotype(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, causal_idx::AbstractVector{<:Integer}, beta, h2::Real; binary::Bool=false)
    G = _matrix(genotypes)
    n, p = size(G)
    idx = [Int(i) for i in causal_idx if 1 <= i <= p]
    isempty(idx) && throw(ArgumentError("causal_idx must contain valid variant indices"))

    b = beta isa Real ? fill(Float64(beta), length(idx)) : Float64.(beta)
    length(b) == length(idx) || throw(DimensionMismatch("beta length must match number of causal variants"))

    X = Matrix{Float64}(G[:, idx])
    _impute_missing_by_mean!(X)
    genetic = X * b
    genetic .-= mean(genetic)

    target_h2 = clamp(Float64(h2), 0.0, 1.0)
    vg = max(var(genetic; corrected=false), eps(Float64))
    ve = target_h2 <= eps(Float64) ? vg : vg * (1.0 - target_h2) / target_h2
    noise = randn(n) .* sqrt(max(ve, eps(Float64)))
    y = genetic + noise

    if binary
        prob = 1.0 ./ (1.0 .+ exp.(-y))
        return Float64.(rand.(Bernoulli.(clamp.(prob, eps(Float64), 1.0 - eps(Float64)))))
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, y, "simulate_phenotype")
end

"""
    ld_expand_credible_set(cs_indices, genotypes; r2=0.8)

Expand a credible set by adding variants in LD above threshold with any seed SNP.
"""
function ld_expand_credible_set(cs_indices::AbstractVector{<:Integer}, genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}; r2::Real=0.8)
    G = _matrix(genotypes)
    p = size(G, 2)
    seeds = unique([Int(i) for i in cs_indices if 1 <= i <= p])
    isempty(seeds) && return Int[]

    expanded = Set(seeds)
    for i in seeds
        for j in 1:p
            _, _, pair_r2 = _pairwise_ld(G, i, j)
            isfinite(pair_r2) && pair_r2 >= r2 && push!(expanded, j)
        end
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, sort!(collect(expanded)), "ld_expand_credible_set")
end

"""
    estimate_effective_n(result)

Estimate effective sample size from summary statistics.
"""
function estimate_effective_n(result::GWASResult)
    if result.sample_size > 0
        return Float64(result.sample_size)
    end
    valid = isfinite.(result.beta) .& isfinite.(result.standard_error) .& (result.standard_error .> 0.0)
    sum(valid) > 0 || return NaN
    neff = median((result.beta[valid] ./ result.standard_error[valid]) .^ 2)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, max(neff, 0.0), "estimate_effective_n")
end

"""
    popcorn_genetic_correlation(z1, z2, ld_scores_1, ld_scores_2)

Cross-population genetic correlation with LD-aware weighting.
"""
function popcorn_genetic_correlation(z1::AbstractVector{<:Real}, z2::AbstractVector{<:Real}, ld_scores_1::AbstractVector{<:Real}, ld_scores_2::AbstractVector{<:Real})
    length(z1) == length(z2) == length(ld_scores_1) == length(ld_scores_2) || throw(DimensionMismatch("all vectors must have equal length"))
    valid = isfinite.(z1) .& isfinite.(z2) .& isfinite.(ld_scores_1) .& isfinite.(ld_scores_2)
    n = sum(valid)
    n >= 3 || return (NaN, NaN)

    x = Float64.(z1[valid])
    y = Float64.(z2[valid])
    w = 1.0 ./ max.(Float64.(ld_scores_1[valid]) .+ Float64.(ld_scores_2[valid]), eps(Float64))
    w ./= sum(w)

    mx = sum(w .* x)
    my = sum(w .* y)
    cov_xy = sum(w .* (x .- mx) .* (y .- my))
    var_x = sum(w .* (x .- mx) .^ 2)
    var_y = sum(w .* (y .- my) .^ 2)
    rg = cov_xy / sqrt(max(var_x * var_y, eps(Float64)))
    se = sqrt(max((1.0 - rg^2) / max(n - 2, 1), 0.0))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (rg, se), "popcorn_genetic_correlation")
end

"""
    twas_scan(genotypes, expression_weights, phenotype; covariates=nothing)

Perform a simple TWAS by predicting expression and testing against phenotype.
"""
function twas_scan(genotypes::Union{GenotypeMatrix,AbstractMatrix{<:Real}}, expression_weights::Dict{<:AbstractString,<:AbstractVector{<:Real}}, phenotype::AbstractVector{<:Real}; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing)
    G = _matrix(genotypes)
    _validate_gwas_inputs(G; phenotype=phenotype, covariates=covariates)
    n, p = size(G)
    design, _ = _covariate_matrix(covariates, n)
    y_all = Float64.(phenotype)
    design_valid = vec(all(isfinite.(design), dims=2))
    Xg = Matrix{Float64}(G)
    _impute_missing_by_mean!(Xg)

    rows = NamedTuple[]
    for (gene, w_raw) in expression_weights
        w = Float64.(w_raw)
        length(w) == p || continue
        expr = Xg * w
        valid = isfinite.(y_all) .& isfinite.(expr) .& design_valid
        sum(valid) > size(design, 2) || continue
        X = hcat(design[valid, :], expr[valid])
        coeff, stderr, zscore, pvalue = _linear_fit_last_coefficient(X, y_all[valid])
        push!(rows, (gene=String(gene), beta=coeff, se=stderr, zscore=zscore, pvalue=pvalue))
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, DataFrame(rows), "twas_scan")
end

"""
    smr_test(gwas_result, eqtl_result)

Summary-data MR across matched SNPs between GWAS and eQTL summary statistics.
"""
function smr_test(gwas_result::GWASResult, eqtl_result::GWASResult)
    aligned = harmonise_alleles(gwas_result, eqtl_result)
    g = aligned.result_a
    e = aligned.result_b

    rows = NamedTuple[]
    for i in eachindex(g.snp_ids)
        be = e.beta[i]
        se_be = e.standard_error[i]
        bg = g.beta[i]
        se_bg = g.standard_error[i]
        (isfinite(be) && isfinite(bg) && isfinite(se_be) && isfinite(se_bg) && abs(be) > eps(Float64)) || continue

        beta_smr = bg / be
        var_smr = (se_bg^2 / be^2) + (bg^2 * se_be^2 / be^4)
        se_smr = sqrt(max(var_smr, eps(Float64)))
        z_smr = beta_smr / se_smr
        p_smr = 2.0 * ccdf(Normal(), abs(z_smr))
        push!(rows, (snp_id=String(g.snp_ids[i]), beta=beta_smr, se=se_smr, zscore=z_smr, pvalue=p_smr))
    end

    table = DataFrame(rows)
    if nrow(table) == 0
        return (; summary=DataFrame(), top_hit=nothing)
    end
    top = table[argmin(table.pvalue), :]
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (; summary=table, top_hit=top), "smr_test")
end

"""
    liftover(result, from_build, to_build; position_map=nothing, allow_passthrough=false)

Coordinate liftover using a custom `(chromosome, position) => new_position` map.
By default this throws when no map is supplied; set `allow_passthrough=true`
to return unchanged positions explicitly.
"""
function liftover(result::GWASResult, from_build, to_build; position_map::Union{Nothing,Dict}=nothing, allow_passthrough::Bool=false)
    from_str = String(from_build)
    to_str = String(to_build)
    if from_str == to_str
        return GWASResult(result.snp_ids, result.chromosomes, result.positions, result.alleles, result.gene_ids, result.beta, result.standard_error, result.zscore, result.pvalue, result.sample_size, result.covariate_names, result.phenotype_name, string(result.method, "+liftover:" , from_str, "->", to_str))
    end

    lifted_positions = similar(result.positions)
    if position_map === nothing
        allow_passthrough || throw(ArgumentError("no liftover map supplied for $(from_str) -> $(to_str); pass position_map or set allow_passthrough=true"))
        @warn "No liftover map supplied; returning unchanged coordinates" from_build=from_str to_build=to_str
        lifted_positions .= result.positions
    else
        for i in eachindex(result.positions)
            key = (String(result.chromosomes[i]), Int(result.positions[i]))
            lifted_positions[i] = Int(get(position_map, key, result.positions[i]))
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, GWASResult(result.snp_ids, result.chromosomes, lifted_positions, result.alleles, result.gene_ids, result.beta, result.standard_error, result.zscore, result.pvalue, result.sample_size, result.covariate_names, result.phenotype_name, string(result.method, "+liftover:" , from_str, "->", to_str)), "liftover")
end

"""
    functional_annotation(result, vep_cache; pvalue_threshold=5e-8)

Annotate GWAS hits using a precomputed annotation dictionary keyed by SNP ID.
"""
function functional_annotation(result::GWASResult, vep_cache::Dict{<:AbstractString,<:Any}; pvalue_threshold::Real=5e-8)
    threshold = clamp(Float64(pvalue_threshold), 0.0, 1.0)
    rows = DataFrame(
        snp_id=String[],
        chromosome=String[],
        position=Int[],
        beta=Float64[],
        pvalue=Float64[],
        gene=String[],
        consequence=String[])
    for i in eachindex(result.snp_ids)
        pval = Float64(result.pvalue[i])
        (isfinite(pval) && pval <= threshold) || continue
        snp_id = String(result.snp_ids[i])
        ann = get(vep_cache, snp_id, nothing)
        consequence = ann === nothing ? "unknown" : (ann isa AbstractDict && haskey(ann, "consequence") ? String(ann["consequence"]) : String(ann))
        gene = ann isa AbstractDict && haskey(ann, "gene") ? String(ann["gene"]) : String(result.gene_ids[i])
        push!(rows, (
            snp_id = snp_id,
            chromosome = String(result.chromosomes[i]),
            position = Int(result.positions[i]),
            beta = Float64(result.beta[i]),
            pvalue = pval,
            gene = gene,
            consequence = consequence))
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, rows, "functional_annotation")
end

"""
    ebi_lookup(snp_ids; max_requests=50, max_concurrency=16)

Fetch SNP annotations from the EBI GWAS Catalog REST endpoint.
"""
function ebi_lookup(snp_ids::AbstractVector{<:AbstractString}; max_requests::Int=50, max_concurrency::Int=16)
    n = min(length(snp_ids), max_requests)
    n <= 0 && return DataFrame(snp_id=String[], found=Bool[], rsid=String[], functional_class=String[], error=String[])

    requests = String.(snp_ids[1:n])
    ntasks = max(1, min(max_concurrency, n))
    rows = asyncmap(requests; ntasks=ntasks) do snp
        url = "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/$(snp)"
        response_path = tempname()
        try
            Downloads.download(url, response_path)
            payload = JSON.parse(read(response_path, String))
            rsid = haskey(payload, "rsId") ? String(payload["rsId"]) : snp
            func_class = haskey(payload, "functionalClass") ? String(payload["functionalClass"]) : ""
            return (snp_id=snp, found=true, rsid=rsid, functional_class=func_class, error="")
        catch err
        return (snp_id=snp, found=false, rsid="", functional_class="", error=sprint(showerror, err))
        finally
        isfile(response_path) && rm(response_path; force=true)
        end
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, DataFrame(rows), "ebi_lookup")
end

end
