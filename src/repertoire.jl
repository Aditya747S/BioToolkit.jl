# ==============================================================================
# repertoire.jl — Immune Repertoire Contig Parsing & Clonotype Analytics
#
# Provides:
#   - Contig parsing for 10x Genomics, MiXCR, TRUST4, AIRR formats
#   - Clonotype combination across chains (scRepertoire-style)
#   - Overlap/diversity metrics (Jaccard, Morisita-Horn, Chao1, ACE, Shannon, Gini-Simpson)
#   - CDR3 clustering with Levenshtein distance (GLIPH/CoNGA-style)
#   - Convergent clonotype detection
#
# References:
#   - Borcherding et al. (2020) Front Immunol 11:569710 (scRepertoire)
#   - Glanville et al. (2017) Nature 547:94-98 (GLIPH)
#   - Ramesh et al. (2023) Nat Methods 20:1533-1541 (CoNGA)
# ==============================================================================

module Repertoire

using DataFrames
using Statistics
using LinearAlgebra
using SparseArrays
using Random
using Distributions
using SpecialFunctions
using Dates

using ..BioToolkit: BioSequence, AminoAcidAlphabet, DNAAlphabet
using ..BioToolkit: ProvenanceContext, ThreadSafeProvenanceContext, active_provenance_context
using ..BioToolkit: provenance_result!, provenance_parent_ids, provenance_record, register_provenance!
import ..BioToolkit: shannon_entropy

@inline function _register_repertoire_result!(_ctx, result, operation; parents=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export RepertoireContig, RepertoireSample, ContigData
export ContigClonotype, ClonotypeAssignment
export AIRRReader, airr_read, airr_parse_row
export tenx_contig_reader, mixcr_reader, trust4_reader
export combine_clonotypes, clonotype_network
export clonotype_overlap_matrix, overlap_jaccard, overlap_morisita_horn
export clonality_index, shannon_entropy, gini_simpson, chao1_richness, ace_richness
export diversity_curve, hill_diversity, rarefaction_curve
export cdr3_clustering, cdr3_levenshtein_distance
export convergent_clonotypes, public_clonotype_analysis
export repertoire_statistics, clonal_distribution_summary
export clone_size_distribution, expanded_clonotype_test
export vj_usage_matrix, vj_pairing_analysis
export repertoire_distance, repertoire_divergence
export repertoire_heatmap, clonotype_venn

# ---- Core Types ---------------------------------------------------------------

"""
    RepertoireContig

Single reconstructed immune receptor contig from V(D)J assembly.
"""
struct RepertoireContig
    contig_id::String
    barcode::String
    chain::Symbol           # :TRA, :TRB, :TRG, :TRD, :IGH, :IGL, :IGK
    v_gene::String
    d_gene::String
    j_gene::String
    c_gene::String
    cdr1_nt::String
    cdr2_nt::String
    cdr3_nt::String
    cdr1_aa::String
    cdr2_aa::String
    cdr3_aa::String
    sequence::String
    length::Int
    umi_count::Int
    read_count::Int
    quality::Float64
    productive::Bool
    full_length::Bool
    annotations::Dict{String,Any}
end

"""
    ContigData

Container for a collection of contigs from a single sample/library.
"""
struct ContigData
    contigs::Vector{RepertoireContig}
    sample_id::String
    metadata::Dict{String,Any}
    provenance::Vector{Dict{String,Any}}
end

Base.length(c::ContigData) = length(c.contigs)
Base.isempty(c::ContigData) = isempty(c.contigs)

"""
    ClonotypeAssignment

Minimal fingerprint of a cell's immune receptor for clonotyping.
"""
struct ClonotypeAssignment
    barcode::String
    sample_id::String
    clonotype_id::String
    v_gene::String
    j_gene::String
    cdr3_aa::String
    cdr3_nt::String
    cdr3_length::Int
    chain::Symbol
    umi_count::Int
    clone_size::Int
end

"""
    RepertoireSample

Full repertoire for one sample: per-cell clonotype assignments + per-clonotype summary.
"""
struct RepertoireSample
    assignments::Vector{ClonotypeAssignment}
    clonotypes::DataFrame
    sample_id::String
    metadata::Dict{String,Any}
end

Base.length(r::RepertoireSample) = length(r.assignments)

function ContigData(contigs::AbstractVector{<:RepertoireContig}, sample_id::String; metadata::AbstractDict=Dict{String,Any}(), provenance=Dict{String,Any}[])
    return ContigData(collect(contigs), String(sample_id), Dict{String,Any}(metadata), collect(Dict{String,Any}, provenance))
end

# ---- AIRR Format Parser -------------------------------------------------------

function airr_parse_row(row::AbstractDict)
    chain_raw = get(row, "locus", get(row, "chain", ""))
    chain = if occursin(r"TRA|TR.alpha"i, String(chain_raw))
        :TRA
    elseif occursin(r"TRB|TR.beta"i, String(chain_raw))
        :TRB
    elseif occursin(r"TRG|TR.gamma"i, String(chain_raw))
        :TRG
    elseif occursin(r"TRD|TR.delta"i, String(chain_raw))
        :TRD
    elseif occursin(r"IGH|IgH"i, String(chain_raw))
        :IGH
    elseif occursin(r"IGL|IgL"i, String(chain_raw))
        :IGL
    elseif occursin(r"IGK|IgK"i, String(chain_raw))
        :IGK
    else
        :UNKNOWN
    end

    seq = String(get(row, "sequence", ""))
    return RepertoireContig(
        String(get(row, "sequence_id", "")),
        String(get(row, "cell_id", get(row, "barcode", ""))),
        chain,
        String(get(row, "v_call", "")),
        String(get(row, "d_call", "")),
        String(get(row, "j_call", "")),
        String(get(row, "c_call", "")),
        String(get(row, "cdr1", "")),
        String(get(row, "cdr2", "")),
        String(get(row, "junction_aa", "")),
        String(get(row, "cdr1_aa", "")),
        String(get(row, "cdr2_aa", "")),
        String(get(row, "junction_aa", "")),
        uppercase(seq),
        length(seq),
        get(row, "duplicate_count", get(row, "umi_count", 1)),
        get(row, "read_count", get(row, "read_qual", 1)),
        get(row, "qual", get(row, "quality", 0.0)),
        get(row, "productive", true),
        get(row, "full_length", false),
        Dict{String,Any}(
            "airr_fields" => row,
            "formatted" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
        )
    )
end

function airr_read(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    isfile(path) || throw(ArgumentError("AIRR file not found: $path"))
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]

    rows = if endswith(lowercase(path), ".tsv") || endswith(lowercase(path), ".airr")
        df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
        eachrow(df)
    elseif endswith(lowercase(path), ".csv")
        df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
        eachrow(df)
    elseif endswith(lowercase(path), ".json") || endswith(lowercase(path), ".jsonl")
        lines = readlines(path)
        [JSON.parse(line) for line in lines if !isempty(strip(line))]
    else
        df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
        eachrow(df)
    end

    contigs = RepertoireContig[airr_parse_row(row) for row in rows]
    metadata = Dict{String,Any}(
        "source" => "airr",
        "path" => String(path),
        "contig_count" => length(contigs)
    )
    result = ContigData(contigs, sample; metadata=metadata)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "airr_read";
            parameters=(path=String(path), sample=sample, contig_count=length(contigs)))
    end
    return result
end

# ---- 10x Genomics V(D)J Parsers ----------------------------------------------

function tenx_contig_reader(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]
    isfile(path) || throw(ArgumentError("10x contig file not found: $path"))

    df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
    contigs = RepertoireContig[]

    for row in eachrow(df)
        chain = if occursin(r"TRA|TR A|alpha", String(get(row, :chain, "")))
            :TRA
        elseif occursin(r"TRB|TR B|beta", String(get(row, :chain, "")))
            :TRB
        elseif occursin(r"IGH|Ig H|heavy", String(get(row, :chain, "")))
            :IGH
        elseif occursin(r"IGL|Ig L|lambda", String(get(row, :chain, "")))
            :IGL
        elseif occursin(r"IGK|Ig K|kappa", String(get(row, :chain, "")))
            :IGK
        else
            :UNKNOWN
        end

        seq = String(get(row, :sequence, ""))
        push!(contigs, RepertoireContig(
            String(get(row, :contig_id, "")),
            String(get(row, :barcode, "")),
            chain,
            String(get(row, :v_gene, "")),
            String(get(row, :d_gene, "")),
            String(get(row, :j_gene, "")),
            String(get(row, :c_gene, "")),
            String(get(row, :cdr1_seq, "")),
            String(get(row, :cdr2_seq, "")),
            String(get(row, :cdr3_seq, "")),
            String(get(row, :cdr1_aa, "")),
            String(get(row, :cdr2_aa, "")),
            String(get(row, :cdr3_aa, "")),
            uppercase(seq),
            length(seq),
            get(row, :umis, get(row, :umi_count, 1)),
            get(row, :read_count, 1),
            get(row, :high_confidence, get(row, :quality, 1.0)),
            get(row, :productive, true),
            get(row, :full_length, false),
            Dict{String,Any}("source" => "10x", "format" => "contig")
        ))
    end

    result = ContigData(contigs, sample; metadata=Dict{String,Any}("source" => "10x", "path" => String(path), "contig_count" => length(contigs)))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "tenx_contig_reader";
            parameters=(path=String(path), sample=sample, contig_count=length(contigs)))
    end
    return result
end

# ---- MiXCR Reader ------------------------------------------------------------

function mixcr_reader(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]
    isfile(path) || throw(ArgumentError("MiXCR file not found: $path"))

    df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
    contigs = RepertoireContig[]

    for row in eachrow(df)
        chain = if occursin(r"TRA", String(get(row, :allele, get(row, :chain, ""))))
            :TRA
        elseif occursin(r"TRB", String(get(row, :allele, get(row, :chain, ""))))
            :TRB
        elseif occursin(r"IGH", String(get(row, :allele, get(row, :chain, ""))))
            :IGH
        else
            :UNKNOWN
        end

        seq = String(get(row, :sequence, get(row, :clonesequence, "")))
        push!(contigs, RepertoireContig(
            String(get(row, :id, "")),
            String(get(row, :count, get(row, :barcode, ""))),
            chain,
            String(get(row, :v_gene, "")),
            String(get(row, :d_gene, "")),
            String(get(row, :j_gene, "")),
            String(get(row, :c_gene, "")),
            String(get(row, :cdr1, "")),
            String(get(row, :cdr2, "")),
            String(get(row, :aaSeqCDR3, get(row, :cdr3_aa, ""))),
            String(get(row, :cdr1_aa, "")),
            String(get(row, :cdr2_aa, "")),
            String(get(row, :aaSeqCDR3, get(row, :cdr3_aa, ""))),
            uppercase(seq),
            length(seq),
            get(row, :count, 1),
            get(row, :read_count, 1),
            1.0,
            true,
            false,
            Dict{String,Any}("source" => "mixcr", "path" => String(path))
        ))
    end

    result = ContigData(contigs, sample; metadata=Dict{String,Any}("source" => "mixcr", "path" => String(path), "contig_count" => length(contigs)))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "mixcr_reader";
            parameters=(path=String(path), sample=sample, contig_count=length(contigs)))
    end
    return result
end

# ---- TRUST4 Reader -----------------------------------------------------------

function trust4_reader(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]
    isfile(path) || throw(ArgumentError("TRUST4 file not found: $path"))

    contigs = RepertoireContig[]
    open(String(path)) do io
        for line in eachline(io)
            isempty(strip(line)) && continue
            fields = split(strip(line), '\t')
            length(fields) < 8 && continue

            contig_id = fields[1]
            barcode = length(fields) > 8 ? fields[9] : contig_id
            chain_str = fields[2]
            v_gene = fields[3]
            d_gene = length(fields) > 4 ? fields[5] : ""
            j_gene = length(fields) > 6 ? fields[7] : ""
            cdr3_nt = length(fields) > 8 ? fields[9] : ""
            seq = length(fields) > 9 ? fields[10] : ""

            chain = if occursin(r"TRA", chain_str)
                :TRA
            elseif occursin(r"TRB", chain_str)
                :TRB
            elseif occursin(r"IGH", chain_str)
                :IGH
            elseif occursin(r"IGL|IGK", chain_str)
                :IGL
            else
                :UNKNOWN
            end

            push!(contigs, RepertoireContig(
                contig_id, barcode, chain, v_gene, d_gene, j_gene, "",
                "", "", cdr3_nt, "", "", cdr3_nt,
                uppercase(seq), length(seq), 1, 1, 1.0, true, false,
                Dict{String,Any}("source" => "trust4")
            ))
        end
    end

    result = ContigData(contigs, sample; metadata=Dict{String,Any}("source" => "trust4", "path" => String(path), "contig_count" => length(contigs)))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "trust4_reader";
            parameters=(path=String(path), sample=sample, contig_count=length(contigs)))
    end
    return result
end

# ---- Clonotype Combination ---------------------------------------------------

function combine_clonotypes(data::ContigData; chain_pair::Symbol=:TRB, fallback_chain::Symbol=:TRA)
    filtered = filter(c -> c.chain in (chain_pair, fallback_chain), data.contigs)
    isempty(filtered) && return ContigData(RepertoireContig[], data.sample_id; metadata=data.metadata)

    by_barcode = Dict{String,Vector{RepertoireContig}}()
    for c in filtered
        push!(get!(by_barcode, c.barcode, RepertoireContig[]), c)
    end

    clonotypes = Dict{String,Vector{RepertoireContig}}()
    for (barcode, contigs) in by_barcode
        primary = filter(c -> c.chain == chain_pair, contigs)
        isempty(primary) && (primary = filter(c -> c.chain == fallback_chain, contigs))
        isempty(primary) && continue

        key = join(sort!([c.cdr3_aa for c in primary if !isempty(c.cdr3_aa)]), "|")
        isempty(key) && continue
        push!(get!(clonotypes, key, RepertoireContig[]), primary[1])
    end

    combined = RepertoireContig[]
    for contigs in values(clonotypes)
        if !isempty(contigs)
            c = contigs[1]
            push!(combined, RepertoireContig(
                c.contig_id, c.barcode, c.chain, c.v_gene, c.d_gene, c.j_gene, c.c_gene,
                c.cdr1_nt, c.cdr2_nt, c.cdr3_nt, c.cdr1_aa, c.cdr2_aa, c.cdr3_aa,
                c.sequence, c.length, c.umi_count, c.read_count, c.quality, c.productive, c.full_length,
                Dict{String,Any}("clone_size" => length(contigs), "combined_from" => [x.contig_id for x in contigs])
            ))
        end
    end

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "combine_clonotypes";
            parameters=(sample_id=data.sample_id, input_contigs=length(filtered), output_contigs=length(combined)))
    end
    return ContigData(combined, data.sample_id; metadata=merge(data.metadata, Dict{String,Any}("combined" => true, "chain_pair" => String(chain_pair))))
end

# ---- Clonotype Assignment -----------------------------------------------------

function _make_clonotype_key(cdr3_aa::String, v_gene::String, j_gene::String)
    isempty(cdr3_aa) && return "unresolved"
    v_family = first(split(v_gene, '*'))
    j_family = first(split(j_gene, '*'))
    return "$(v_family)_$(j_family)_$(cdr3_aa)"
end

function RepertoireSample(data::ContigData; clone_threshold::Int=1)
    by_barcode = Dict{String,Vector{RepertoireContig}}()
    for c in data.contigs
        push!(get!(by_barcode, c.barcode, RepertoireContig[]), c)
    end

    assignments = ClonotypeAssignment[]
    clone_counts = Dict{String,Int}()

    for (barcode, contigs) in by_barcode
        best = contigs[1]
        for c in contigs[2:end]
            if !isempty(c.cdr3_aa) && (isempty(best.cdr3_aa) || c.umi_count > best.umi_count)
                best = c
            end
        end

        key = _make_clonotype_key(best.cdr3_aa, best.v_gene, best.j_gene)
        clone_counts[key] = get(clone_counts, key, 0) + 1

        push!(assignments, ClonotypeAssignment(
            barcode, data.sample_id, key, best.v_gene, best.j_gene,
            best.cdr3_aa, best.cdr3_nt, length(best.cdr3_aa), best.chain,
            best.umi_count, clone_counts[key]
        ))
    end

    df_rows = [(clonotype_id=k, v_gene=first(split(assignments[findfirst(a -> a.clonotype_id == k, assignments)].v_gene, '*')),
                j_gene=first(split(assignments[findfirst(a -> a.clonotype_id == k, assignments)].j_gene, '*')),
                cdr3_aa=assignments[findfirst(a -> a.clonotype_id == k, assignments)].cdr3_aa,
                cdr3_length=assignments[findfirst(a -> a.clonotype_id == k, assignments)].cdr3_length,
                clone_size=clone_counts[k],
                frequency=clone_counts[k] / max(length(assignments), 1))
               for k in keys(clone_counts)]
    clonotypes_df = DataFrames.DataFrame(df_rows)
    sort!(clonotypes_df, :clone_size, rev=true)

    result = RepertoireSample(assignments, clonotypes_df, data.sample_id, merge(data.metadata, Dict{String,Any}("n_cells" => length(assignments), "n_clonotypes" => length(clone_counts))))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "repertoire_sample";
            parameters=(sample_id=data.sample_id, n_cells=length(assignments), n_clonotypes=length(clone_counts)))
    end
    return result
end

# ---- Diversity Metrics -------------------------------------------------------

"""
    shannon_entropy(clone_sizes)

Shannon entropy H = -Σ pᵢ log(pᵢ). Uses natural log.
"""
function shannon_entropy(clone_sizes::AbstractVector{<:Integer})
    n = sum(clone_sizes)
    n == 0 && return 0.0
    p = Float64.(clone_sizes) ./ n
    return -sum(pi > 0 ? pi * log(pi) : 0.0 for pi in p)
end

"""
    gini_simpson(clone_sizes)

Gini-Simpson index = 1 - Σ pᵢ². Range [0, 1].
"""
function gini_simpson(clone_sizes::AbstractVector{<:Integer})
    n = sum(clone_sizes)
    n == 0 && return 0.0
    p = Float64.(clone_sizes) ./ n
    return 1.0 - sum(abs2, p)
end

"""
    clonality_index(clone_sizes)

Clonality = 1 - normalized Shannon entropy. Range [0, 1].
"""
function clonality_index(clone_sizes::AbstractVector{<:Integer})
    h = shannon_entropy(clone_sizes)
    h_max = log(max(length(clone_sizes), 1))
    return h_max > 0 ? 1.0 - (h / h_max) : 0.0
end

"""
    chao1_richness(clone_sizes)

Chao1 richness estimator for unseen species.
"""
function chao1_richness(clone_sizes::AbstractVector{<:Integer})
    n = sum(clone_sizes)
    f1 = count(==(1), clone_sizes)
    f2 = count(==(2), clone_sizes)
    s_obs = count(>(0), clone_sizes)

    if f2 > 0
        return s_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))
    elseif f1 > 0
        return s_obs + (f1 * (f1 - 1)) / 2
    else
        return s_obs
    end
end

"""
    ace_richness(clone_sizes; cutoff=10)

ACE (Abundance-based Coverage Estimator) richness.
"""
function ace_richness(clone_sizes::AbstractVector{<:Integer}; cutoff::Int=10)
    n = sum(clone_sizes)
    counts = [count(==(i), clone_sizes) for i in 1:cutoff]
    f1, fc = counts[1], sum(counts[i] * (i - 1) for i in 1:cutoff)
    s_obs = count(>(0), clone_sizes)

    if fc == 0 || n == 0
        return Float64(s_obs)
    end

    gamma_sq = max(0.0, (cutoff / fc) * sum(i * (i - 1) * counts[i] for i in 1:cutoff) - 1.0)
    C_ace = 1.0 - f1 / n
    return C_ace > 0 ? s_obs + (f1 / C_ace) * gamma_sq : Float64(s_obs)
end

"""
    hill_diversity(clone_sizes, q)

Hill numbers of order q. q=0 → richness, q=1 → exp(Shannon), q=2 → 1/Simpson.
"""
function hill_diversity(clone_sizes::AbstractVector{<:Integer}, q::Real)
    n = sum(clone_sizes)
    n == 0 && return 0.0
    p = Float64.(clone_sizes) ./ n

    if q == 1
        exp_shannon = exp(-sum(pi > 0 ? pi * log(pi) : 0.0 for pi in p))
        return isfinite(exp_shannon) ? exp_shannon : 0.0
    elseif q == 0
        return Float64(count(>(0), clone_sizes))
    elseif q == 2
        return 1.0 / max(sum(abs2, p), eps(Float64))
    else
        sum_pq = sum(pi > 0 ? pi^q : 0.0 for pi in p)
        return sum_pq > 0 ? sum_pq^(1.0 / (1.0 - q)) : 0.0
    end
end

"""
    diversity_curve(clone_sizes; q_values=[0.0, 1.0, 2.0, 3.0])

Compute Hill diversity for a range of q orders.
"""
function diversity_curve(clone_sizes::AbstractVector{<:Integer}; q_values::AbstractVector{<:Real}=[0.0, 1.0, 2.0, 3.0])
    hill_vals = Float64[hill_diversity(clone_sizes, q) for q in q_values]
    return DataFrames.DataFrame(q=q_values, hill_number=hill_vals)
end

"""
    rarefaction_curve(clone_sizes; n_points=50)

Compute sample-size-corrected rarefaction and extrapolation curves.
"""
function rarefaction_curve(clone_sizes::AbstractVector{<:Integer}; n_points::Int=50)
    n_total = sum(clone_sizes)
    n_clones = count(>(0), clone_sizes)
    n_points = clamp(n_points, 2, n_total)

    rarefaction = Float64[]
    extrapolation = Float64[]
    sample_sizes = Int[]

    for m in range(1, n_total, length=n_points)
        if m <= n_total
            push!(sample_sizes, round(Int, m))
            push!(rarefaction, _rarefaction_observed(clone_sizes, round(Int, m)))
            push!(extrapolation, _extrapolation_chao(clone_sizes, round(Int, m)))
        end
    end

    return DataFrames.DataFrame(sample_size=sample_sizes, rarefaction=rarefaction, extrapolation=extrapolation)
end

function _rarefaction_observed(clone_sizes::AbstractVector{<:Integer}, m::Int)
    n = sum(clone_sizes)
    m >= n && return Float64(count(>(0), clone_sizes))

    u = 0.0
    for ni in clone_sizes
        ni > 0 && (u += 1.0 - exp(log1p(-ni / n) * m / n))
    end
    return u
end

function _extrapolation_chao(clone_sizes::AbstractVector{<:Integer}, m::Int)
    f1 = count(==(1), clone_sizes)
    f2 = count(==(2), clone_sizes)
    s_obs = count(>(0), clone_sizes)

    if f2 > 0
        f0_hat = (f1^2) / (2 * f2)
        return s_obs + f0_hat * (1 - exp(-f1 / (f0_hat + f1) * (m / sum(clone_sizes))))
    else
        return s_obs
    end
end

"""
    repertoire_statistics(sample)

Compute comprehensive repertoire statistics.
"""
function repertoire_statistics(sample::RepertoireSample)
    cs = Int.(sample.clonotypes.clone_size)
    n = sum(cs)
    n == 0 && return DataFrames.DataFrame()

    results = Dict{Symbol,Any}(
        :sample_id => sample.sample_id,
        :n_cells => length(sample.assignments),
        :n_clonotypes => length(cs),
        :n_singletons => count(==(1), cs),
        :n_expanded => count(>=(2), cs),
        :clone_fraction_expanded => sum(cs[cs .>= 2]) / max(n, 1),
        :shannon_entropy => shannon_entropy(cs),
        :gini_simpson => gini_simpson(cs),
        :clonality => clonality_index(cs),
        :chao1 => chao1_richness(cs),
        :ace => ace_richness(cs),
        :d50 => _d50(cs),
        :top10_clone_fraction => sum(sort(cs, rev=true)[1:min(10, length(cs))]) / max(n, 1),
    )

    if length(cs) >= 2
        results[:cv] = std(cs) / mean(cs)
    else
        results[:cv] = 0.0
    end

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "repertoire_statistics";
            parameters=(sample_id=sample.sample_id, n_cells=results[:n_cells], n_clonotypes=results[:n_clonotypes], clonality=results[:clonality]))
    end

    return DataFrames.DataFrame([(; results...)])
end

function _d50(clone_sizes::AbstractVector{<:Integer})
    n_total = sum(clone_sizes)
    n_total == 0 && return 0.0
    target = n_total / 2
    cumsum = 0
    sorted = sort(clone_sizes, rev=true)
    for (i, cs) in enumerate(sorted)
        cumsum += cs
        if cumsum >= target
            return i / max(length(sorted), 1)
        end
    end
    return 1.0
end

# ---- Overlap Metrics ---------------------------------------------------------

"""
    clonotype_overlap_matrix(samples; metric=:jaccard)

Build pairwise overlap matrix between repertoire samples.
"""
function clonotype_overlap_matrix(samples::AbstractVector{<:RepertoireSample}; metric::Symbol=:jaccard)
    n = length(samples)
    result = zeros(Float64, n, n)

    clonotypes_per_sample = [Set{String}(sample.clonotypes.clonotype_id) for sample in samples]

    for i in 1:n
        result[i, i] = 1.0
        for j in (i+1):n
            overlap = if metric == :jaccard
                overlap_jaccard(clonotypes_per_sample[i], clonotypes_per_sample[j])
            elseif metric == :morisita_horn
                overlap_morisita_horn(samples[i], samples[j])
            else
                0.0
            end
            result[i, j] = overlap
            result[j, i] = overlap
        end
    end

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "clonotype_overlap_matrix";
            parameters=(n_samples=n, metric=String(metric)))
    end
    return result
end

function overlap_jaccard(set_a::Set{String}, set_b::Set{String})
    union_size = length(union(set_a, set_b))
    union_size == 0 && return 0.0
    return length(intersect(set_a, set_b)) / union_size
end

function overlap_morisita_horn(sample_a::RepertoireSample, sample_b::RepertoireSample)
    df_a = sample_a.clonotypes
    df_b = sample_b.clonotypes

    all_ids = sort!(union(df_a.clonotype_id, df_b.clonotype_id))
    n_a = length(sample_a.assignments)
    n_b = length(sample_b.assignments)

    dict_a = Dict(zip(df_a.clonotype_id, df_a.clone_size))
    dict_b = Dict(zip(df_b.clonotype_id, df_b.clone_size))
    counts_a = [get(dict_a, id, 0) for id in all_ids]
    counts_b = [get(dict_b, id, 0) for id in all_ids]

    sum_a = sum(counts_a)
    sum_b = sum(counts_b)

    if sum_a == 0 || sum_b == 0
        return 0.0
    end

    p_a = Float64.(counts_a) ./ sum_a
    p_b = Float64.(counts_b) ./ sum_b

    numerator = 2.0 * sum(p_a .* p_b)
    denominator = sum(p_a .^ 2) + sum(p_b .^ 2)

    return denominator > 0 ? numerator / denominator : 0.0
end

# ---- CDR3 Clustering ---------------------------------------------------------

function _levenshtein_distance(a::AbstractString, b::AbstractString)
    la = ncodeunits(a)
    lb = ncodeunits(b)

    la == 0 && return lb
    lb == 0 && return la

    prev = zeros(Int, lb + 1)
    curr = zeros(Int, lb + 1)

    for j in 0:lb
        prev[j + 1] = j
    end

    for i in 1:la
        curr[1] = i
        for j in 1:lb
            cost = a[i] == b[j] ? 0 : 1
            curr[j + 1] = min(prev[j + 1] + 1, curr[j] + 1, prev[j] + cost)
        end
        prev, curr = curr, prev
    end

    return prev[lb + 1]
end

"""
    cdr3_levenshtein_distance(cdr3_list)

Compute pairwise Levenshtein distance matrix for CDR3 sequences.
"""
function cdr3_levenshtein_distance(cdr3_list::AbstractVector{<:AbstractString}; threaded::Bool=true)
    n = length(cdr3_list)
    D = zeros(Float64, n, n)

    @inbounds for i in 1:n
        D[i, i] = 0.0
        for j in (i+1):n
            d = Float64(_levenshtein_distance(cdr3_list[i], cdr3_list[j]))
            D[i, j] = d
            D[j, i] = d
        end
    end

    return D
end

"""
    cdr3_clustering(cdr3_sequences; threshold::Real=2, min_size::Int=2)

Single-linkage hierarchical clustering on Levenshtein distances.
Returns cluster assignments and cluster summary.
"""
function cdr3_clustering(
    cdr3_sequences::AbstractVector{<:AbstractString};
    threshold::Real=2,
    min_size::Int=2,
    v_genes::Union{Nothing,AbstractVector{<:AbstractString}}=nothing)
    seqs = uppercase.(String.(cdr3_sequences))
    n = length(seqs)
    n == 0 && return (assignments=Int[], clusters=DataFrames.DataFrame())

    D = cdr3_levenshtein_distance(seqs)
    parent = collect(1:n)

    function find!(p, x)
        while p[x] != x
            p[x] = p[p[x]]
            x = p[x]
        end
        return x
    end
    function union!(p, x, y)
        rx = find!(p, x)
        ry = find!(p, y)
        rx != ry && (p[ry] = rx)
    end

    for i in 1:n
        for j in (i+1):n
            if D[i, j] <= Float64(threshold)
                union!(parent, i, j)
            end
        end
    end

    labels = [find!(parent, i) for i in 1:n]
    unique_roots = sort!(unique(labels))
    root_map = Dict(r => k for (k, r) in enumerate(unique_roots))
    cluster_ids = [root_map[l] for l in labels]

    df = DataFrames.DataFrame(sequence=seqs, cluster_id=cluster_ids)
    if v_genes !== nothing && length(v_genes) == n
        df[!, :v_gene] = String.(v_genes)
    end

    clusters = combine(groupby(df, :cluster_id), nrow => :size)
    filter!(r -> r.size >= min_size, clusters)
    sort!(clusters, :size, rev=true)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "cdr3_clustering";
            parameters=(n_sequences=n, threshold=Float64(threshold), min_size=min_size, n_clusters=nrow(clusters)))
    end
    return (assignments=cluster_ids, clusters=clusters)
end

# ---- Convergent Clonotypes ---------------------------------------------------

"""
    convergent_clonotypes(samples; max_hamming::Int=1, min_shared::Int=2)

Find public clonotypes shared across multiple samples.
"""
function convergent_clonotypes(
    samples::AbstractVector{<:RepertoireSample};
    max_hamming::Int=1,
    min_shared::Int=2)
    n_samples = length(samples)
    sample_ids = [s.sample_id for s in samples]

    clonotype_sets = [Set{String}(s.clonotypes.clonotype_id) for s in samples]

    direct_shared = Set{String}()
    for i in 1:n_samples
        for j in (i+1):n_samples
            for ct in clonotype_sets[i]
                ct in clonotype_sets[j] && push!(direct_shared, ct)
            end
        end
    end

    convergent = DataFrames.DataFrame(
        clonotype_id=String[],
        n_samples=Int[],
        shared_in=String[],
        min_clone_size=Int[],
        max_clone_size=Int[]
    )

    for ct in direct_shared
        shared_in = String[]
        clone_sizes = Int[]
        for (i, sample) in enumerate(samples)
            if ct in clonotype_sets[i]
                push!(shared_in, sample_ids[i])
                row_idx = findfirst(==(ct), sample.clonotypes.clonotype_id)
                if row_idx !== nothing
                    push!(clone_sizes, Int(sample.clonotypes.clone_size[row_idx]))
                end
            end
        end

        if length(shared_in) >= min_shared
            push!(convergent, (ct, length(shared_in), join(shared_in, ", "),
                minimum(clone_sizes), maximum(clone_sizes)))
        end
    end

    sort!(convergent, :n_samples, rev=true)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "convergent_clonotypes";
            parameters=(n_samples=n_samples, min_shared=min_shared, n_convergent=nrow(convergent)))
    end
    return convergent
end

# ---- V(J) Usage Analysis -----------------------------------------------------

"""
    vj_usage_matrix(sample)

Compute V and J gene usage frequencies for a sample.
"""
function vj_usage_matrix(sample::RepertoireSample)
    v_usage = combine(groupby(DataFrame(v_gene=first.(split.(String.([a.v_gene for a in sample.assignments]), '*'))), :v_gene), nrow => :count)
    v_usage[!, :frequency] = v_usage.count ./ max(sum(v_usage.count), 1)
    sort!(v_usage, :count, rev=true)

    j_usage = combine(groupby(DataFrame(j_gene=first.(split.(String.([a.j_gene for a in sample.assignments]), '*'))), :j_gene), nrow => :count)
    j_usage[!, :frequency] = j_usage.count ./ max(sum(j_usage.count), 1)
    sort!(j_usage, :count, rev=true)

    return (v_usage=v_usage, j_usage=j_usage)
end

"""
    vj_pairing_analysis(sample)

Compute V-J pairing frequencies for paired-chain repertoires.
"""
function vj_pairing_analysis(sample::RepertoireSample)
    pairs = [(v_gene=first(split(a.v_gene, '*')), j_gene=first(split(a.j_gene, '*')), cdr3=a.cdr3_aa) for a in sample.assignments]
    df = DataFrames.DataFrame(pairs)
    result = combine(groupby(df, [:v_gene, :j_gene]), nrow => :count)
    result[!, :frequency] = result.count ./ max(sum(result.count), 1)
    sort!(result, :count, rev=true)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "vj_pairing_analysis";
            parameters=(sample_id=sample.sample_id, n_pairs=nrow(result)))
    end
    return result
end

# ---- Clone Size Distribution -------------------------------------------------

"""
    clonal_distribution_summary(sample)

Summarize the clone size distribution.
"""
function clonal_distribution_summary(sample::RepertoireSample)
    cs = Int.(sample.clonotypes.clone_size)

    return DataFrames.DataFrame(
        metric=["n_cells", "n_clones", "singletons", "doubletons", "expanded_(2-9)", "large_(10-99)", "hyperexpanded_(>=100)",
                "shannon", "gini_simpson", "clonality", "chao1", "d50", "top1_freq", "top10_freq", "gini_coefficient"],
        value=[
            length(sample.assignments), length(cs), count(==(1), cs), count(==(2), cs),
            count(c -> 2 <= c <= 9, cs), count(c -> 10 <= c <= 99, cs), count(>=(100), cs),
            shannon_entropy(cs), gini_simpson(cs), clonality_index(cs), chao1_richness(cs), _d50(cs),
            maximum(cs) / max(sum(cs), 1),
            sum(sort(cs, rev=true)[1:min(10, length(cs))]) / max(sum(cs), 1),
            _gini_coefficient(cs)
        ]
    )
end

function _gini_coefficient(clone_sizes::AbstractVector{<:Integer})
    n = length(clone_sizes)
    n <= 1 && return 0.0
    sorted = sort(clone_sizes)
    return (2.0 / (n * sum(sorted))) * sum((2*i - n - 1) * sorted[i] for i in 1:n)
end

"""
    expanded_clonotype_test(sample; threshold::Int=5)

Identify significantly expanded clonotypes using a multinomial test.
"""
function expanded_clonotype_test(sample::RepertoireSample; threshold::Int=5)
    cs = Int.(sample.clonotypes.clone_size)
    total = sum(cs)
    n_clones = length(cs)

    n_expected = total / n_clones
    results = DataFrames.DataFrame(
        clonotype_id=String[],
        clone_size=Int[],
        expected=Float64[],
        fold_change=Float64[],
        pvalue=Float64[]
    )

    for (i, ct_row) in enumerate(eachrow(sample.clonotypes))
        obs = ct_row.clone_size
        expected = n_expected
        fc = expected > 0 ? obs / expected : Inf
        pval = if obs >= threshold
            1.0 - cdf(Poisson(expected), obs - 1)
        else
            1.0
        end
        push!(results, (ct_row.clonotype_id, obs, expected, fc, pval))
    end

    m = nrow(results)
    padj = similar(results.pvalue)
    perm = sortperm(results.pvalue)
    running = 1.0
    for i in m:-1:1
        j = perm[i]
        running = min(running, results.pvalue[j] * m / i)
        padj[j] = min(running, 1.0)
    end
    results[!, :padj] = padj

    sort!(results, :pvalue)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "expanded_clonotype_test";
            parameters=(sample_id=sample.sample_id, n_clones=n_clones, threshold=threshold))
    end
    return results
end

# ---- Repertoire Distance ------------------------------------------------------

"""
    repertoire_distance(sample_a, sample_b; metric=:jaccard)

Compute distance between two repertoires.
"""
function repertoire_distance(sample_a::RepertoireSample, sample_b::RepertoireSample; metric::Symbol=:jaccard)
    set_a = Set{String}(sample_a.clonotypes.clonotype_id)
    set_b = Set{String}(sample_b.clonotypes.clonotype_id)

    d = if metric == :jaccard
        1.0 - overlap_jaccard(set_a, set_b)
    elseif metric == :morisita_horn
        1.0 - overlap_morisita_horn(sample_a, sample_b)
    else
        1.0
    end
    return clamp(d, 0.0, 1.0)
end

"""
    repertoire_divergence(samples)

Compute pairwise divergence matrix.
"""
function repertoire_divergence(samples::AbstractVector{<:RepertoireSample}; metric::Symbol=:jaccard)
    n = length(samples)
    D = zeros(Float64, n, n)

    for i in 1:n
        D[i, i] = 0.0
        for j in (i+1):n
            d = repertoire_distance(samples[i], samples[j]; metric=metric)
            D[i, j] = d
            D[j, i] = d
        end
    end

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "repertoire_divergence";
            parameters=(n_samples=n, metric=String(metric)))
    end
    return D
end

end
