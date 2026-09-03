# ==============================================================================
# repertoire.jl — Comprehensive Immune Repertoire & Single-Cell AIRR Analytics
#
# Features:
#   - Full AIRR, 10x Genomics, MiXCR, TRUST4, Adaptive readers/writers
#   - Advanced single-cell chain pairing, multi-chain/orphan chain analytics
#   - BK-Tree metric spatial indexing & zero-allocation thread-local Levenshtein DP
#   - Comprehensive diversity, richness & evenness indices + Bootstrap 95% CIs
#   - Sample-size corrected analytical rarefaction & extrapolation curves (Chao et al. 2014)
#   - Pairwise repertoire overlap & Monte Carlo permutation significance tests
#   - High-performance multi-threaded CDR3 Levenshtein, Hamming & TCRdist matrices
#   - GLIPH2-style motif discovery & enrichment + Positional CDR3 sequence logos (PWM)
#   - Sequence clustering & Clonotype similarity network analysis
#   - Cohort public & convergent clonotype identification
#   - V(D)J usage, pairing, JSD/KLD divergence & differential usage hypothesis tests
#   - Somatic Hypermutation (SHM) rate & BCR mutation hotspot/coldspot profiling
#   - Atchley (1-5) & Kidera (1-10) biophysical factor profiling of CDR3s
#   - Clonotype lineage Minimum Spanning Trees (MST) & spectratyping analytics
#   - Visualization helpers for heatmaps and Venn/Upset diagrams
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
using CSV
using JSON

using ..BioToolkit: BioSequence, AminoAcidAlphabet, DNAAlphabet
using ..BioToolkit: ProvenanceContext, ThreadSafeProvenanceContext, active_provenance_context
using ..BioToolkit: provenance_result!, provenance_parent_ids, provenance_record, register_provenance!
import ..BioToolkit: shannon_entropy, differential_vj_usage

@inline function _register_repertoire_result!(_ctx, result, operation; parents=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

# Export core types
export RepertoireContig, RepertoireSample, ContigData, ContigClonotype, ClonotypeAssignment, AIRRReader, BKTree

# Export parsers and writers
export airr_read, airr_write, airr_parse_row, tenx_contig_reader, tenx_clonotype_reader, mixcr_reader, trust4_reader, adaptive_reader

# Export clonotype grouping & multi-chain pairing
export group_clonotypes, combine_clonotypes, clone_size_distribution, clonotype_network, multi_chain_pairing_analysis

# Export diversity, richness & bootstrap CIs
export shannon_entropy, gini_simpson, inverse_simpson, clonality_index, chao1_richness, ace_richness, d50_index, d75_index, pielou_evenness, hill_diversity, diversity_curve, rarefaction_curve, repertoire_statistics, clonal_distribution_summary, bootstrap_diversity_ci

# Export overlap, distance & permutation tests
export overlap_jaccard, overlap_morisita_horn, overlap_sorensen_dice, overlap_overlap_coefficient, overlap_cosine, clonotype_overlap_matrix, repertoire_distance, repertoire_divergence, repertoire_overlap_permutation_test

# Export sequence clustering, trees & fast metrics
export cdr3_levenshtein_distance, cdr3_hamming_distance, tcrdist_matrix, cdr3_clustering, bktree_range_query

# Export motif discovery & positional logos
export gliph2_motif_enrichment, positional_motif_logo

# Export biophysical CDR3 properties & Atchley/Kidera factors
export repertoire_cdr3_physicochemical_properties

# Export spectratype & lineage tree analytics
export spectratype_analysis, clonotype_lineage_tree

# Export convergent & public clonotypes
export convergent_clonotypes, public_clonotype_analysis

# Export V(D)J gene usage & tests
export vj_usage_matrix, vj_pairing_analysis, vj_usage_divergence, differential_vj_usage, expanded_clonotype_test

# Export BCR-specific analytics & hotspot profiling
export bcr_shm_rate, isotype_usage_summary, bcr_mutation_hotspot_profile

# Export visualization helpers
export repertoire_heatmap, clonotype_venn

# ---- Core Types ---------------------------------------------------------------

"""
    RepertoireContig

Single reconstructed immune receptor contig from V(D)J assembly.
"""
struct RepertoireContig
    contig_id::String
    barcode::String
    chain::Symbol           # :TRA, :TRB, :TRG, :TRD, :IGH, :IGL, :IGK, :UNKNOWN
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
Base.iterate(c::ContigData, state=1) = state > length(c.contigs) ? nothing : (c.contigs[state], state + 1)
Base.getindex(c::ContigData, i::Int) = c.contigs[i]

function ContigData(contigs::AbstractVector{<:RepertoireContig}, sample_id::String; metadata::AbstractDict=Dict{String,Any}(), provenance=Dict{String,Any}[])
    return ContigData(collect(contigs), String(sample_id), Dict{String,Any}(metadata), collect(Dict{String,Any}, provenance))
end

"""
    ContigClonotype

Summary struct for a single unique clonotype across cells.
"""
struct ContigClonotype
    clonotype_id::String
    cdr3_aa::String
    cdr3_nt::String
    v_gene::String
    d_gene::String
    j_gene::String
    c_gene::String
    chain::Symbol
    clone_size::Int
    frequency::Float64
    barcodes::Vector{String}
    contig_ids::Vector{String}
end

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
Base.isempty(r::RepertoireSample) = isempty(r.assignments)

"""
    AIRRReader

Streaming iterator for AIRR standard files.
"""
struct AIRRReader
    path::String
    sample_id::String
    rows::Vector{Dict{String,Any}}
end

function AIRRReader(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    isfile(path) || throw(ArgumentError("AIRR file not found: $path"))
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]
    df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
    rows = [Dict{String,Any}(string(k) => (ismissing(v) ? "" : v) for (k, v) in pairs(r)) for r in eachrow(df)]
    return AIRRReader(String(path), sample, rows)
end

Base.length(r::AIRRReader) = length(r.rows)
Base.iterate(r::AIRRReader, state=1) = state > length(r.rows) ? nothing : (airr_parse_row(r.rows[state]), state + 1)

# ---- Fast Spatial Metric Indexing (BK-Tree) -----------------------------------

struct BKTreeNode
    sequence::String
    units::Vector{UInt8}
    index::Int
    children::Dict{Int, BKTreeNode}
end

"""
    BKTree

Burkhard-Keller tree for sub-quadratic spatial sequence metric search.
"""
struct BKTree
    root::Union{Nothing, BKTreeNode}
end

@inline function _fast_levenshtein(a::AbstractVector{UInt8}, b::AbstractVector{UInt8}, buf1::Vector{Int}, buf2::Vector{Int})
    la, lb = length(a), length(b)
    la == 0 && return lb
    lb == 0 && return la
    if la > lb
        return _fast_levenshtein(b, a, buf1, buf2)
    end

    @inbounds for i in 1:(la + 1)
        buf1[i] = i - 1
    end

    @inbounds for j in 1:lb
        buf2[1] = j
        bj = b[j]
        @simd for i in 1:la
            cost = a[i] == bj ? 0 : 1
            buf2[i + 1] = min(buf1[i + 1] + 1, buf2[i] + 1, buf1[i] + cost)
        end
        copyto!(buf1, 1, buf2, 1, la + 1)
    end

    return buf1[la + 1]
end

function BKTree(sequences::AbstractVector{<:AbstractString})
    isempty(sequences) && return BKTree(nothing)
    units_list = [codeunits(String(s)) for s in sequences]
    root = BKTreeNode(String(sequences[1]), units_list[1], 1, Dict{Int, BKTreeNode}())

    buf1 = zeros(Int, 128)
    buf2 = zeros(Int, 128)

    for i in 2:length(sequences)
        seq = String(sequences[i])
        u = units_list[i]
        curr = root
        while true
            d = _fast_levenshtein(curr.units, u, buf1, buf2)
            if haskey(curr.children, d)
                curr = curr.children[d]
            else
                curr.children[d] = BKTreeNode(seq, u, i, Dict{Int, BKTreeNode}())
                break
            end
        end
    end
    return BKTree(root)
end

"""
    bktree_range_query(tree, query, max_dist)

Find indices of sequences in BKTree within `max_dist` Levenshtein distance.
"""
function bktree_range_query(tree::BKTree, query::AbstractString, max_dist::Int)
    results = Int[]
    tree.root === nothing && return results
    q_units = codeunits(String(query))

    buf1 = zeros(Int, 128)
    buf2 = zeros(Int, 128)

    stack = [tree.root]
    while !isempty(stack)
        node = pop!(stack)
        d = _fast_levenshtein(node.units, q_units, buf1, buf2)
        if d <= max_dist
            push!(results, node.index)
        end
        low = d - max_dist
        high = d + max_dist
        for (dist, child) in node.children
            if low <= dist <= high
                push!(stack, child)
            end
        end
    end
    return results
end

# ---- AIRR & External Format Parsers / Writers ---------------------------------

"""
    airr_parse_row(row)

Parse a single AIRR-compliant dictionary or row into a `RepertoireContig`.
"""
function airr_parse_row(row::Union{AbstractDict, DataFrames.DataFrameRow})
    get_str(k, default="") = begin
        v = get(row, k, get(row, Symbol(k), default))
        ismissing(v) ? default : String(string(v))
    end
    get_int(k, default=1) = begin
        v = get(row, k, get(row, Symbol(k), default))
        ismissing(v) ? default : (v isa Integer ? Int(v) : parse(Int, string(v)))
    end
    get_float(k, default=0.0) = begin
        v = get(row, k, get(row, Symbol(k), default))
        ismissing(v) ? default : (v isa Real ? Float64(v) : parse(Float64, string(v)))
    end
    get_bool(k, default=true) = begin
        v = get(row, k, get(row, Symbol(k), default))
        if ismissing(v)
            return default
        elseif v isa Bool
            return v
        else
            s = lowercase(strip(string(v)))
            return s in ("true", "t", "1", "yes")
        end
    end

    chain_raw = get_str("locus", get_str("chain", ""))
    chain = if occursin(r"TRA|TR.alpha"i, chain_raw)
        :TRA
    elseif occursin(r"TRB|TR.beta"i, chain_raw)
        :TRB
    elseif occursin(r"TRG|TR.gamma"i, chain_raw)
        :TRG
    elseif occursin(r"TRD|TR.delta"i, chain_raw)
        :TRD
    elseif occursin(r"IGH|IgH"i, chain_raw)
        :IGH
    elseif occursin(r"IGL|IgL"i, chain_raw)
        :IGL
    elseif occursin(r"IGK|IgK"i, chain_raw)
        :IGK
    else
        :UNKNOWN
    end

    seq = get_str("sequence", "")
    cdr3_aa = get_str("junction_aa", get_str("cdr3_aa", ""))
    cdr3_nt = get_str("junction", get_str("cdr3", get_str("cdr3_nt", "")))

    return RepertoireContig(
        get_str("sequence_id", get_str("contig_id", "")),
        get_str("cell_id", get_str("barcode", "")),
        chain,
        get_str("v_call", get_str("v_gene", "")),
        get_str("d_call", get_str("d_gene", "")),
        get_str("j_call", get_str("j_gene", "")),
        get_str("c_call", get_str("c_gene", "")),
        get_str("cdr1", get_str("cdr1_nt", "")),
        get_str("cdr2", get_str("cdr2_nt", "")),
        cdr3_nt,
        get_str("cdr1_aa", ""),
        get_str("cdr2_aa", ""),
        cdr3_aa,
        uppercase(seq),
        length(seq),
        get_int("duplicate_count", get_int("umi_count", 1)),
        get_int("read_count", 1),
        get_float("quality", get_float("qual", 1.0)),
        get_bool("productive", true),
        get_bool("full_length", false),
        Dict{String,Any}(
            "source" => "airr",
            "formatted" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
        )
    )
end

"""
    airr_read(path; sample_id=nothing)

Read AIRR TSV/CSV/JSON file into `ContigData`.
"""
function airr_read(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    isfile(path) || throw(ArgumentError("AIRR file not found: $path"))
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]

    path_lower = lowercase(path)
    contigs = RepertoireContig[]

    if endswith(path_lower, ".json") || endswith(path_lower, ".jsonl")
        lines = readlines(path)
        for line in lines
            sline = strip(line)
            isempty(sline) && continue
            obj = JSON.parse(sline)
            if obj isa Dict
                push!(contigs, airr_parse_row(obj))
            end
        end
    else
        df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
        for r in eachrow(df)
            push!(contigs, airr_parse_row(r))
        end
    end

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

"""
    airr_write(path, data)

Export `ContigData` or `RepertoireSample` to an AIRR standard TSV file.
"""
function airr_write(path::AbstractString, data::ContigData)
    df = DataFrames.DataFrame(
        sequence_id = [c.contig_id for c in data.contigs],
        cell_id = [c.barcode for c in data.contigs],
        locus = [string(c.chain) for c in data.contigs],
        v_call = [c.v_gene for c in data.contigs],
        d_call = [c.d_gene for c in data.contigs],
        j_call = [c.j_gene for c in data.contigs],
        c_call = [c.c_gene for c in data.contigs],
        junction = [c.cdr3_nt for c in data.contigs],
        junction_aa = [c.cdr3_aa for c in data.contigs],
        sequence = [c.sequence for c in data.contigs],
        duplicate_count = [c.umi_count for c in data.contigs],
        read_count = [c.read_count for c in data.contigs],
        productive = [c.productive for c in data.contigs]
    )
    CSV.write(path, df; delim='\t')
    return path
end

function airr_write(path::AbstractString, sample::RepertoireSample)
    df = DataFrames.DataFrame(
        sequence_id = [a.clonotype_id for a in sample.assignments],
        cell_id = [a.barcode for a in sample.assignments],
        locus = [string(a.chain) for a in sample.assignments],
        v_call = [a.v_gene for a in sample.assignments],
        j_call = [a.j_gene for a in sample.assignments],
        junction = [a.cdr3_nt for a in sample.assignments],
        junction_aa = [a.cdr3_aa for a in sample.assignments],
        duplicate_count = [a.umi_count for a in sample.assignments],
        productive = fill(true, length(sample.assignments))
    )
    CSV.write(path, df; delim='\t')
    return path
end

"""
    tenx_contig_reader(path; sample_id=nothing)

Read 10x Genomics `filtered_contig_annotations.csv` / `.json`.
"""
function tenx_contig_reader(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]
    isfile(path) || throw(ArgumentError("10x contig file not found: $path"))

    df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
    contigs = RepertoireContig[]

    for row in eachrow(df)
        chain_str = String(get(row, :chain, ""))
        chain = if occursin(r"TRA|TR A|alpha"i, chain_str)
            :TRA
        elseif occursin(r"TRB|TR B|beta"i, chain_str)
            :TRB
        elseif occursin(r"TRG|TR G|gamma"i, chain_str)
            :TRG
        elseif occursin(r"TRD|TR D|delta"i, chain_str)
            :TRD
        elseif occursin(r"IGH|Ig H|heavy"i, chain_str)
            :IGH
        elseif occursin(r"IGL|Ig L|lambda"i, chain_str)
            :IGL
        elseif occursin(r"IGK|Ig K|kappa"i, chain_str)
            :IGK
        else
            :UNKNOWN
        end

        seq = String(get(row, :sequence, ""))
        umis = get(row, :umis, get(row, :umi_count, 1))
        reads = get(row, :reads, get(row, :read_count, 1))
        high_conf = get(row, :high_confidence, true)
        high_conf_val = high_conf isa Bool ? (high_conf ? 1.0 : 0.0) : Float64(high_conf)
        prod = get(row, :productive, true)
        prod_val = prod isa Bool ? prod : lowercase(string(prod)) in ("true", "t", "1")
        full_l = get(row, :full_length, false)
        full_l_val = full_l isa Bool ? full_l : lowercase(string(full_l)) in ("true", "t", "1")

        push!(contigs, RepertoireContig(
            String(get(row, :contig_id, "")),
            String(get(row, :barcode, "")),
            chain,
            String(get(row, :v_gene, "")),
            String(get(row, :d_gene, "")),
            String(get(row, :j_gene, "")),
            String(get(row, :c_gene, "")),
            String(get(row, :cdr1, get(row, :cdr1_seq, ""))),
            String(get(row, :cdr2, get(row, :cdr2_seq, ""))),
            String(get(row, :cdr3, get(row, :cdr3_seq, ""))),
            String(get(row, :cdr1_aa, "")),
            String(get(row, :cdr2_aa, "")),
            String(get(row, :cdr3_aa, "")),
            uppercase(seq),
            length(seq),
            umis isa Integer ? Int(umis) : parse(Int, string(umis)),
            reads isa Integer ? Int(reads) : parse(Int, string(reads)),
            high_conf_val,
            prod_val,
            full_l_val,
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

"""
    tenx_clonotype_reader(path; sample_id=nothing)

Read 10x Genomics `clonotypes.csv`.
"""
function tenx_clonotype_reader(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]
    isfile(path) || throw(ArgumentError("10x clonotypes file not found: $path"))

    df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
    return (sample_id=sample, clonotypes=df)
end

"""
    mixcr_reader(path; sample_id=nothing)

Read MiXCR export file into `ContigData`.
"""
function mixcr_reader(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]
    isfile(path) || throw(ArgumentError("MiXCR file not found: $path"))

    df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
    contigs = RepertoireContig[]

    for row in eachrow(df)
        allele_str = String(get(row, :bestVGene, get(row, :v_gene, get(row, :chain, ""))))
        chain = if occursin(r"TRA"i, allele_str)
            :TRA
        elseif occursin(r"TRB"i, allele_str)
            :TRB
        elseif occursin(r"TRG"i, allele_str)
            :TRG
        elseif occursin(r"TRD"i, allele_str)
            :TRD
        elseif occursin(r"IGH"i, allele_str)
            :IGH
        elseif occursin(r"IGL"i, allele_str)
            :IGL
        elseif occursin(r"IGK"i, allele_str)
            :IGK
        else
            :UNKNOWN
        end

        seq = String(get(row, :targetSequences, get(row, :sequence, get(row, :clonesequence, ""))))
        cdr3_aa = String(get(row, :aaSeqCDR3, get(row, :cdr3_aa, "")))
        cdr3_nt = String(get(row, :nSeqCDR3, get(row, :cdr3_nt, "")))
        clone_count = get(row, :cloneCount, get(row, :count, 1))

        push!(contigs, RepertoireContig(
            String(get(row, :cloneId, get(row, :id, ""))),
            String(get(row, :cloneId, get(row, :barcode, ""))),
            chain,
            String(get(row, :bestVGene, get(row, :v_gene, ""))),
            String(get(row, :bestDGene, get(row, :d_gene, ""))),
            String(get(row, :bestJGene, get(row, :j_gene, ""))),
            String(get(row, :bestCGene, get(row, :c_gene, ""))),
            String(get(row, :nSeqCDR1, get(row, :cdr1, ""))),
            String(get(row, :nSeqCDR2, get(row, :cdr2, ""))),
            cdr3_nt,
            String(get(row, :aaSeqCDR1, get(row, :cdr1_aa, ""))),
            String(get(row, :aaSeqCDR2, get(row, :cdr2_aa, ""))),
            cdr3_aa,
            uppercase(seq),
            length(seq),
            clone_count isa Integer ? Int(clone_count) : parse(Int, string(clone_count)),
            1, 1.0, true, false,
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

"""
    trust4_reader(path; sample_id=nothing)

Read TRUST4 TSV report into `ContigData`.
"""
function trust4_reader(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]
    isfile(path) || throw(ArgumentError("TRUST4 file not found: $path"))

    contigs = RepertoireContig[]
    open(String(path)) do io
        for line in eachline(io)
            sline = strip(line)
            (isempty(sline) || startswith(sline, "#")) && continue
            fields = split(sline, '\t')
            length(fields) < 4 && continue

            contig_id = fields[1]
            barcode = length(fields) >= 9 ? fields[9] : contig_id
            chain_str = fields[2]
            v_gene = length(fields) >= 3 ? fields[3] : ""
            d_gene = length(fields) >= 5 ? fields[5] : ""
            j_gene = length(fields) >= 7 ? fields[7] : ""
            cdr3_aa = length(fields) >= 9 ? fields[9] : ""
            cdr3_nt = length(fields) >= 10 ? fields[10] : ""
            seq = length(fields) >= 11 ? fields[11] : ""

            chain = if occursin(r"TRA"i, chain_str)
                :TRA
            elseif occursin(r"TRB"i, chain_str)
                :TRB
            elseif occursin(r"TRG"i, chain_str)
                :TRG
            elseif occursin(r"TRD"i, chain_str)
                :TRD
            elseif occursin(r"IGH"i, chain_str)
                :IGH
            elseif occursin(r"IGL"i, chain_str)
                :IGL
            elseif occursin(r"IGK"i, chain_str)
                :IGK
            else
                :UNKNOWN
            end

            push!(contigs, RepertoireContig(
                contig_id, barcode, chain, v_gene, d_gene, j_gene, "",
                "", "", cdr3_nt, "", "", cdr3_aa,
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

"""
    adaptive_reader(path; sample_id=nothing)

Read ImmunoSEQ / Adaptive Biotechnologies TSV format into `ContigData`.
"""
function adaptive_reader(path::AbstractString; sample_id::Union{Nothing,String}=nothing)
    sample = sample_id !== nothing ? String(sample_id) : splitext(basename(path))[1]
    isfile(path) || throw(ArgumentError("Adaptive Biotechnologies file not found: $path"))

    df = DataFrames.DataFrame(CSV.read(path, DataFrames.DataFrame; comment="#"))
    contigs = RepertoireContig[]

    for row in eachrow(df)
        v_gene = String(get(row, :vMaxResolved, get(row, :vGeneName, get(row, :v_gene, ""))))
        j_gene = String(get(row, :jMaxResolved, get(row, :jGeneName, get(row, :j_gene, ""))))
        d_gene = String(get(row, :dMaxResolved, get(row, :dGeneName, get(row, :d_gene, ""))))
        cdr3_aa = String(get(row, :aminoAcid, get(row, :cdr3_aa, "")))
        cdr3_nt = String(get(row, :nucleotide, get(row, :cdr3_nt, "")))
        count = get(row, :templates, get(row, :sequenceCount, get(row, :count, 1)))

        chain = if occursin(r"TCRB|TRB"i, v_gene)
            :TRB
        elseif occursin(r"TCRA|TRA"i, v_gene)
            :TRA
        elseif occursin(r"IGH"i, v_gene)
            :IGH
        else
            :UNKNOWN
        end

        push!(contigs, RepertoireContig(
            String(get(row, :readId, get(row, :cloneId, ""))),
            String(get(row, :readId, "")),
            chain, v_gene, d_gene, j_gene, "",
            "", "", cdr3_nt, "", "", cdr3_aa,
            uppercase(cdr3_nt), length(cdr3_nt),
            count isa Integer ? Int(count) : parse(Int, string(count)),
            1, 1.0, true, false,
            Dict{String,Any}("source" => "adaptive")
        ))
    end

    result = ContigData(contigs, sample; metadata=Dict{String,Any}("source" => "adaptive", "path" => String(path), "contig_count" => length(contigs)))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "adaptive_reader";
            parameters=(path=String(path), sample=sample, contig_count=length(contigs)))
    end
    return result
end

# ---- Clonotype Grouping & Multi-Chain Analytics -----------------------------

"""
    group_clonotypes(data::ContigData; by::Symbol=:v_j_cdr3_aa)

Group contigs into unique clonotype definitions.
"""
function group_clonotypes(data::ContigData; by::Symbol=:v_j_cdr3_aa)
    by_key = Dict{String, Vector{RepertoireContig}}()

    for c in data.contigs
        key = if by == :v_j_cdr3_aa
            _make_clonotype_key(c.cdr3_aa, c.v_gene, c.j_gene)
        elseif by == :cdr3_aa
            isempty(c.cdr3_aa) ? "unresolved" : c.cdr3_aa
        elseif by == :v_j_cdr3_nt
            "$(c.v_gene)_$(c.j_gene)_$(c.cdr3_nt)"
        elseif by == :cdr3_nt
            isempty(c.cdr3_nt) ? "unresolved" : c.cdr3_nt
        else
            c.barcode
        end
        push!(get!(by_key, key, RepertoireContig[]), c)
    end

    total_cells = length(data.contigs)
    clonotypes = ContigClonotype[]

    for (key, group) in by_key
        c1 = group[1]
        sz = length(group)
        freq = total_cells > 0 ? sz / total_cells : 0.0
        barcodes = [c.barcode for c in group]
        contig_ids = [c.contig_id for c in group]
        push!(clonotypes, ContigClonotype(key, c1.cdr3_aa, c1.cdr3_nt, c1.v_gene, c1.d_gene, c1.j_gene, c1.c_gene, c1.chain, sz, freq, barcodes, contig_ids))
    end

    sort!(clonotypes, by=c->c.clone_size, rev=true)
    return clonotypes
end

"""
    combine_clonotypes(data::ContigData; chain_pair::Symbol=:TRB, fallback_chain::Symbol=:TRA)

Combine single-cell contigs by barcode across paired chains (scRepertoire style).
"""
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

"""
    multi_chain_pairing_analysis(data::ContigData)

Profile single-cell receptor pairing, dual-alpha/beta chains, and orphan chains.
"""
function multi_chain_pairing_analysis(data::ContigData)
    by_barcode = Dict{String,Vector{RepertoireContig}}()
    for c in data.contigs
        push!(get!(by_barcode, c.barcode, RepertoireContig[]), c)
    end

    n_cells = length(by_barcode)
    n_dual_alpha = 0
    n_dual_beta = 0
    n_canonical_pair = 0
    n_orphan_alpha = 0
    n_orphan_beta = 0

    for (barcode, contigs) in by_barcode
        tra = count(c -> c.chain == :TRA, contigs)
        trb = count(c -> c.chain == :TRB, contigs)

        if tra == 1 && trb == 1
            n_canonical_pair += 1
        elseif tra >= 2 && trb == 1
            n_dual_alpha += 1
        elseif tra == 1 && trb >= 2
            n_dual_beta += 1
        elseif tra >= 1 && trb == 0
            n_orphan_alpha += 1
        elseif tra == 0 && trb >= 1
            n_orphan_beta += 1
        end
    end

    df = DataFrames.DataFrame(
        pairing_type=["Canonical (1 TRA + 1 TRB)", "Dual Alpha (>=2 TRA + 1 TRB)", "Dual Beta (1 TRA + >=2 TRB)", "Orphan Alpha (TRA only)", "Orphan Beta (TRB only)"],
        cell_count=[n_canonical_pair, n_dual_alpha, n_dual_beta, n_orphan_alpha, n_orphan_beta],
        percentage=[
            (n_canonical_pair / max(n_cells, 1)) * 100,
            (n_dual_alpha / max(n_cells, 1)) * 100,
            (n_dual_beta / max(n_cells, 1)) * 100,
            (n_orphan_alpha / max(n_cells, 1)) * 100,
            (n_orphan_beta / max(n_cells, 1)) * 100
        ]
    )
    return df
end

function _make_clonotype_key(cdr3_aa::String, v_gene::String, j_gene::String)
    isempty(cdr3_aa) && return "unresolved"
    v_family = first(split(v_gene, '*'))
    j_family = first(split(j_gene, '*'))
    return "$(v_family)_$(j_family)_$(cdr3_aa)"
end

"""
    RepertoireSample(data::ContigData; clone_threshold=1)

Construct full `RepertoireSample` representation with clonotype statistics.
"""
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

    df_rows = [(clonotype_id=k,
                v_gene=first(split(assignments[findfirst(a -> a.clonotype_id == k, assignments)].v_gene, '*')),
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

"""
    clone_size_distribution(sample::RepertoireSample)

Summarize clonotypes into size category bins.
"""
function clone_size_distribution(sample::RepertoireSample)
    cs = Int.(sample.clonotypes.clone_size)
    total_cells = sum(cs)

    bins = [
        ("Singleton (n=1)", count(==(1), cs)),
        ("Doubleton (n=2)", count(==(2), cs)),
        ("Small (3<=n<=5)", count(c -> 3 <= c <= 5, cs)),
        ("Medium (6<=n<=20)", count(c -> 6 <= c <= 20, cs)),
        ("Large (21<=n<=100)", count(c -> 21 <= c <= 100, cs)),
        ("Hyperexpanded (n>100)", count(>(100), cs))
    ]

    bin_names = String[b[1] for b in bins]
    n_clones = Int[b[2] for b in bins]
    proportions = total_cells > 0 ? n_clones ./ sum(n_clones) : zeros(Float64, length(bins))

    return DataFrames.DataFrame(bin=bin_names, clone_count=n_clones, proportion=proportions)
end

# ---- Diversity & Richness Metrics + Bootstrap --------------------------------

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
    inverse_simpson(clone_sizes)

Inverse Simpson index = 1 / Σ pᵢ².
"""
function inverse_simpson(clone_sizes::AbstractVector{<:Integer})
    n = sum(clone_sizes)
    n == 0 && return 0.0
    p = Float64.(clone_sizes) ./ n
    sum_p2 = sum(abs2, p)
    return sum_p2 > 0 ? 1.0 / sum_p2 : 0.0
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

Chao1 richness estimator for unseen clonotypes.
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
    d50_index(clone_sizes) / _d50

Percentage of top clonotypes required to account for 50% of the total cells.
"""
function d50_index(clone_sizes::AbstractVector{<:Integer})
    return _d50(clone_sizes)
end

"""
    d75_index(clone_sizes)

Percentage of top clonotypes required to account for 75% of the total cells.
"""
function d75_index(clone_sizes::AbstractVector{<:Integer})
    n_total = sum(clone_sizes)
    n_total == 0 && return 0.0
    target = 0.75 * n_total
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

"""
    pielou_evenness(clone_sizes)

Pielou's evenness J' = H / ln(S).
"""
function pielou_evenness(clone_sizes::AbstractVector{<:Integer})
    s = count(>(0), clone_sizes)
    s <= 1 && return 0.0
    return shannon_entropy(clone_sizes) / log(s)
end

"""
    hill_diversity(clone_sizes, q)

Hill numbers of order q.
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
    diversity_curve(clone_sizes; q_values=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
"""
function diversity_curve(clone_sizes::AbstractVector{<:Integer}; q_values::AbstractVector{<:Real}=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    hill_vals = Float64[hill_diversity(clone_sizes, q) for q in q_values]
    return DataFrames.DataFrame(q=q_values, hill_number=hill_vals)
end

"""
    rarefaction_curve(clone_sizes; n_points=50)

Compute analytical sample-size-corrected rarefaction and extrapolation curves (Chao et al. 2014).
"""
function rarefaction_curve(clone_sizes::AbstractVector{<:Integer}; n_points::Int=50)
    n_total = sum(clone_sizes)
    n_clones = count(>(0), clone_sizes)
    n_points = clamp(n_points, 2, max(n_total, 2))

    rarefaction = Float64[]
    extrapolation = Float64[]
    sample_sizes = Int[]

    for m in range(1, n_total, length=n_points)
        sm = round(Int, m)
        push!(sample_sizes, sm)
        push!(rarefaction, _rarefaction_observed(clone_sizes, sm))
        push!(extrapolation, _extrapolation_chao(clone_sizes, sm))
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
        return Float64(s_obs)
    end
end

"""
    bootstrap_diversity_ci(sample; n_bootstraps=500, ci=0.95)

Multinomial bootstrap estimation of 95% confidence intervals for diversity metrics.
"""
function bootstrap_diversity_ci(sample::RepertoireSample; n_bootstraps::Int=500, ci::Float64=0.95)
    cs = Int.(sample.clonotypes.clone_size)
    n_total = sum(cs)
    n_total == 0 && return DataFrames.DataFrame()

    probs = Float64.(cs) ./ n_total
    dist = Multinomial(n_total, probs)

    shannon_boots = zeros(Float64, n_bootstraps)
    simpson_boots = zeros(Float64, n_bootstraps)
    chao1_boots = zeros(Float64, n_bootstraps)
    clonality_boots = zeros(Float64, n_bootstraps)

    for b in 1:n_bootstraps
        resample = rand(dist)
        shannon_boots[b] = shannon_entropy(resample)
        simpson_boots[b] = gini_simpson(resample)
        chao1_boots[b] = chao1_richness(resample)
        clonality_boots[b] = clonality_index(resample)
    end

    alpha = 1.0 - ci
    lower_p = alpha / 2.0
    upper_p = 1.0 - lower_p

    metrics = ["Shannon Entropy", "Gini-Simpson", "Chao1 Richness", "Clonality"]
    means = [mean(shannon_boots), mean(simpson_boots), mean(chao1_boots), mean(clonality_boots)]
    lowers = [quantile(shannon_boots, lower_p), quantile(simpson_boots, lower_p), quantile(chao1_boots, lower_p), quantile(clonality_boots, lower_p)]
    uppers = [quantile(shannon_boots, upper_p), quantile(simpson_boots, upper_p), quantile(chao1_boots, upper_p), quantile(clonality_boots, upper_p)]

    return DataFrames.DataFrame(metric=metrics, mean=means, ci_lower=lowers, ci_upper=uppers)
end

"""
    repertoire_statistics(sample)

Compute comprehensive repertoire diversity & expansion statistics.
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
        :inverse_simpson => inverse_simpson(cs),
        :clonality => clonality_index(cs),
        :pielou_evenness => pielou_evenness(cs),
        :chao1 => chao1_richness(cs),
        :ace => ace_richness(cs),
        :d50 => _d50(cs),
        :d75 => d75_index(cs),
        :top10_clone_fraction => sum(sort(cs, rev=true)[1:min(10, length(cs))]) / max(n, 1),
    )

    results[:cv] = length(cs) >= 2 ? std(cs) / mean(cs) : 0.0

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "repertoire_statistics";
            parameters=(sample_id=sample.sample_id, n_cells=results[:n_cells], n_clonotypes=results[:n_clonotypes], clonality=results[:clonality]))
    end

    return DataFrames.DataFrame([(; results...)])
end

"""
    clonal_distribution_summary(sample)

Summarize repertoire metrics into key-value summary table.
"""
function clonal_distribution_summary(sample::RepertoireSample)
    cs = Int.(sample.clonotypes.clone_size)

    return DataFrames.DataFrame(
        metric=["n_cells", "n_clones", "singletons", "doubletons", "expanded_(2-9)", "large_(10-99)", "hyperexpanded_(>=100)",
                "shannon", "gini_simpson", "inverse_simpson", "clonality", "pielou_evenness", "chao1", "d50", "d75", "top1_freq", "top10_freq", "gini_coefficient"],
        value=[
            length(sample.assignments), length(cs), count(==(1), cs), count(==(2), cs),
            count(c -> 2 <= c <= 9, cs), count(c -> 10 <= c <= 99, cs), count(>=(100), cs),
            shannon_entropy(cs), gini_simpson(cs), inverse_simpson(cs), clonality_index(cs), pielou_evenness(cs), chao1_richness(cs), _d50(cs), d75_index(cs),
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
    s = sum(sorted)
    s == 0 && return 0.0
    return (2.0 / (n * s)) * sum((2*i - n - 1) * sorted[i] for i in 1:n)
end

# ---- Overlap, Distance & Permutation Tests -----------------------------------

"""
    overlap_jaccard(set_a, set_b)

Jaccard overlap coefficient between two sets of clonotype IDs.
"""
function overlap_jaccard(set_a::Set{String}, set_b::Set{String})
    union_size = length(union(set_a, set_b))
    union_size == 0 && return 0.0
    return length(intersect(set_a, set_b)) / union_size
end

"""
    overlap_sorensen_dice(set_a, set_b)

Sørensen-Dice overlap coefficient.
"""
function overlap_sorensen_dice(set_a::Set{String}, set_b::Set{String})
    total = length(set_a) + length(set_b)
    total == 0 && return 0.0
    return 2.0 * length(intersect(set_a, set_b)) / total
end

"""
    overlap_overlap_coefficient(set_a, set_b)

Szymkiewicz-Simpson overlap coefficient.
"""
function overlap_overlap_coefficient(set_a::Set{String}, set_b::Set{String})
    min_sz = min(length(set_a), length(set_b))
    min_sz == 0 && return 0.0
    return length(intersect(set_a, set_b)) / min_sz
end

"""
    overlap_morisita_horn(sample_a, sample_b)

Morisita-Horn quantitative overlap index.
"""
function overlap_morisita_horn(sample_a::RepertoireSample, sample_b::RepertoireSample)
    df_a = sample_a.clonotypes
    df_b = sample_b.clonotypes

    all_ids = sort!(union(df_a.clonotype_id, df_b.clonotype_id))

    dict_a = Dict(zip(df_a.clonotype_id, df_a.clone_size))
    dict_b = Dict(zip(df_b.clonotype_id, df_b.clone_size))
    counts_a = [get(dict_a, id, 0) for id in all_ids]
    counts_b = [get(dict_b, id, 0) for id in all_ids]

    sum_a = sum(counts_a)
    sum_b = sum(counts_b)

    (sum_a == 0 || sum_b == 0) && return 0.0

    p_a = Float64.(counts_a) ./ sum_a
    p_b = Float64.(counts_b) ./ sum_b

    numerator = 2.0 * sum(p_a .* p_b)
    denominator = sum(p_a .^ 2) + sum(p_b .^ 2)

    return denominator > 0 ? numerator / denominator : 0.0
end

"""
    overlap_cosine(sample_a, sample_b)

Cosine similarity between clonotype frequency vectors.
"""
function overlap_cosine(sample_a::RepertoireSample, sample_b::RepertoireSample)
    df_a = sample_a.clonotypes
    df_b = sample_b.clonotypes

    all_ids = sort!(union(df_a.clonotype_id, df_b.clonotype_id))

    dict_a = Dict(zip(df_a.clonotype_id, df_a.frequency))
    dict_b = Dict(zip(df_b.clonotype_id, df_b.frequency))
    v_a = [get(dict_a, id, 0.0) for id in all_ids]
    v_b = [get(dict_b, id, 0.0) for id in all_ids]

    norm_a = norm(v_a)
    norm_b = norm(v_b)
    (norm_a == 0 || norm_b == 0) && return 0.0
    return dot(v_a, v_b) / (norm_a * norm_b)
end

"""
    clonotype_overlap_matrix(samples; metric=:jaccard)

Build pairwise overlap matrix between repertoire samples using multi-threading.
"""
function clonotype_overlap_matrix(samples::AbstractVector{<:RepertoireSample}; metric::Symbol=:jaccard)
    n = length(samples)
    result = zeros(Float64, n, n)

    clonotypes_per_sample = [Set{String}(sample.clonotypes.clonotype_id) for sample in samples]

    if Threads.nthreads() > 1
        Threads.@threads for i in 1:n
            result[i, i] = 1.0
            for j in (i+1):n
                overlap = if metric == :jaccard
                    overlap_jaccard(clonotypes_per_sample[i], clonotypes_per_sample[j])
                elseif metric == :sorensen_dice
                    overlap_sorensen_dice(clonotypes_per_sample[i], clonotypes_per_sample[j])
                elseif metric == :overlap_coeff
                    overlap_overlap_coefficient(clonotypes_per_sample[i], clonotypes_per_sample[j])
                elseif metric == :morisita_horn
                    overlap_morisita_horn(samples[i], samples[j])
                elseif metric == :cosine
                    overlap_cosine(samples[i], samples[j])
                else
                    overlap_jaccard(clonotypes_per_sample[i], clonotypes_per_sample[j])
                end
                result[i, j] = overlap
                result[j, i] = overlap
            end
        end
    else
        for i in 1:n
            result[i, i] = 1.0
            for j in (i+1):n
                overlap = if metric == :jaccard
                    overlap_jaccard(clonotypes_per_sample[i], clonotypes_per_sample[j])
                elseif metric == :sorensen_dice
                    overlap_sorensen_dice(clonotypes_per_sample[i], clonotypes_per_sample[j])
                elseif metric == :overlap_coeff
                    overlap_overlap_coefficient(clonotypes_per_sample[i], clonotypes_per_sample[j])
                elseif metric == :morisita_horn
                    overlap_morisita_horn(samples[i], samples[j])
                elseif metric == :cosine
                    overlap_cosine(samples[i], samples[j])
                else
                    overlap_jaccard(clonotypes_per_sample[i], clonotypes_per_sample[j])
                end
                result[i, j] = overlap
                result[j, i] = overlap
            end
        end
    end

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "clonotype_overlap_matrix";
            parameters=(n_samples=n, metric=String(metric)))
    end
    return result
end

"""
    repertoire_distance(sample_a, sample_b; metric=:jaccard)

Distance = 1.0 - similarity between two repertoires.
"""
function repertoire_distance(sample_a::RepertoireSample, sample_b::RepertoireSample; metric::Symbol=:jaccard)
    set_a = Set{String}(sample_a.clonotypes.clonotype_id)
    set_b = Set{String}(sample_b.clonotypes.clonotype_id)

    d = if metric == :jaccard
        1.0 - overlap_jaccard(set_a, set_b)
    elseif metric == :sorensen_dice
        1.0 - overlap_sorensen_dice(set_a, set_b)
    elseif metric == :morisita_horn
        1.0 - overlap_morisita_horn(sample_a, sample_b)
    elseif metric == :cosine
        1.0 - overlap_cosine(sample_a, sample_b)
    else
        1.0 - overlap_jaccard(set_a, set_b)
    end
    return clamp(d, 0.0, 1.0)
end

"""
    repertoire_divergence(samples; metric=:jaccard)

Compute pairwise divergence/distance matrix across samples.
"""
function repertoire_divergence(samples::AbstractVector{<:RepertoireSample}; metric::Symbol=:jaccard)
    n = length(samples)
    D = zeros(Float64, n, n)

    if Threads.nthreads() > 1
        Threads.@threads for i in 1:n
            D[i, i] = 0.0
            for j in (i+1):n
                d = repertoire_distance(samples[i], samples[j]; metric=metric)
                D[i, j] = d
                D[j, i] = d
            end
        end
    else
        for i in 1:n
            D[i, i] = 0.0
            for j in (i+1):n
                d = repertoire_distance(samples[i], samples[j]; metric=metric)
                D[i, j] = d
                D[j, i] = d
            end
        end
    end

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "repertoire_divergence";
            parameters=(n_samples=n, metric=String(metric)))
    end
    return D
end

"""
    repertoire_overlap_permutation_test(sample_a, sample_b; n_permutations=1000, metric=:jaccard)

Monte Carlo permutation significance test for clonotype sharing between repertoires.
"""
function repertoire_overlap_permutation_test(sample_a::RepertoireSample, sample_b::RepertoireSample; n_permutations::Int=1000, metric::Symbol=:jaccard)
    set_a = Set{String}(sample_a.clonotypes.clonotype_id)
    set_b = Set{String}(sample_b.clonotypes.clonotype_id)

    obs_overlap = if metric == :sorensen_dice
        overlap_sorensen_dice(set_a, set_b)
    elseif metric == :morisita_horn
        overlap_morisita_horn(sample_a, sample_b)
    else
        overlap_jaccard(set_a, set_b)
    end

    pooled_clonotypes = collect(union(set_a, set_b))
    n_pooled = length(pooled_clonotypes)
    sz_a = length(set_a)
    sz_b = length(set_b)

    perm_overlaps = zeros(Float64, n_permutations)
    for p in 1:n_permutations
        shuffled = shuffle(pooled_clonotypes)
        perm_a = Set(shuffled[1:sz_a])
        perm_b = Set(shuffled[(sz_a + 1):min(sz_a + sz_b, n_pooled)])
        perm_overlaps[p] = overlap_jaccard(perm_a, perm_b)
    end

    mean_perm = mean(perm_overlaps)
    std_perm = std(perm_overlaps)
    z_score = std_perm > 0 ? (obs_overlap - mean_perm) / std_perm : 0.0
    p_value = count(>=(obs_overlap), perm_overlaps) / max(n_permutations, 1)

    return (observed_overlap=obs_overlap, expected_perm_mean=mean_perm, z_score=z_score, pvalue=p_value)
end

# ---- High-Performance Distance & CDR3 Spatial Clustering ---------------------

"""
    cdr3_levenshtein_distance(cdr3_list; threaded=true)

Compute pairwise Levenshtein distance matrix for CDR3 sequences using pre-allocated DP buffers.
"""
function cdr3_levenshtein_distance(cdr3_list::AbstractVector{<:AbstractString}; threaded::Bool=true)
    n = length(cdr3_list)
    D = zeros(Float64, n, n)
    units = [codeunits(String(s)) for s in cdr3_list]

    if threaded && Threads.nthreads() > 1
        Threads.@threads for i in 1:n
            buf1 = zeros(Int, 128)
            buf2 = zeros(Int, 128)
            D[i, i] = 0.0
            for j in (i+1):n
                d = Float64(_fast_levenshtein(units[i], units[j], buf1, buf2))
                D[i, j] = d
                D[j, i] = d
            end
        end
    else
        buf1 = zeros(Int, 128)
        buf2 = zeros(Int, 128)
        @inbounds for i in 1:n
            D[i, i] = 0.0
            for j in (i+1):n
                d = Float64(_fast_levenshtein(units[i], units[j], buf1, buf2))
                D[i, j] = d
                D[j, i] = d
            end
        end
    end

    return D
end

"""
    cdr3_hamming_distance(cdr3_list; threaded=true)

Compute pairwise Hamming distance matrix for CDR3 sequences of equal length.
"""
function cdr3_hamming_distance(cdr3_list::AbstractVector{<:AbstractString}; threaded::Bool=true)
    n = length(cdr3_list)
    D = zeros(Float64, n, n)
    units = [codeunits(String(s)) for s in cdr3_list]

    function _hamming(u1, u2)
        l1, l2 = length(u1), length(u2)
        if l1 != l2
            return Float64(abs(l1 - l2) + sum(u1[k] != u2[k] for k in 1:min(l1, l2)))
        end
        diff = 0
        @inbounds for k in 1:l1
            diff += (u1[k] != u2[k])
        end
        return Float64(diff)
    end

    if threaded && Threads.nthreads() > 1
        Threads.@threads for i in 1:n
            for j in (i+1):n
                d = _hamming(units[i], units[j])
                D[i, j] = d
                D[j, i] = d
            end
        end
    else
        for i in 1:n
            for j in (i+1):n
                d = _hamming(units[i], units[j])
                D[i, j] = d
                D[j, i] = d
            end
        end
    end

    return D
end

# Default BLOSUM62 matrix lookup for TCRdist
const _BLOSUM62_SCORES = Dict{Tuple{Char,Char},Int}(
    ('A','A')=>4, ('A','R')=>-1, ('A','N')=>-2, ('A','D')=>-2, ('A','C')=>0, ('A','Q')=>-1, ('A','E')=>-1, ('A','G')=>0, ('A','H')=>-2, ('A','I')=>-1, ('A','L')=>-1, ('A','K')=>-1, ('A','M')=>-1, ('A','F')=>-2, ('A','P')=>-1, ('A','S')=>1, ('A','T')=>0, ('A','W')=>-3, ('A','Y')=>-2, ('A','V')=>0,
    ('R','R')=>5, ('R','N')=>0, ('R','D')=>-2, ('R','C')=>-3, ('R','Q')=>1, ('R','E')=>0, ('R','G')=>-2, ('R','H')=>0, ('R','I')=>-3, ('R','L')=>-2, ('R','K')=>2, ('R','M')=>-1, ('R','F')=>-3, ('R','P')=>-2, ('R','S')=>-1, ('R','T')=>-1, ('R','W')=>-3, ('R','Y')=>-2, ('R','V')=>-3,
    ('N','N')=>6, ('N','D')=>1, ('N','C')=>-3, ('N','Q')=>0, ('N','E')=>0, ('N','G')=>0, ('N','H')=>1, ('N','I')=>-3, ('N','L')=>-3, ('N','K')=>0, ('N','M')=>-2, ('N','F')=>-3, ('N','P')=>-2, ('N','S')=>1, ('N','T')=>0, ('N','W')=>-4, ('N','Y')=>-2, ('N','V')=>-3,
    ('D','D')=>6, ('D','C')=>-3, ('D','Q')=>0, ('D','E')=>2, ('D','G')=>-1, ('D','H')=>-1, ('D','I')=>-3, ('D','L')=>-4, ('D','K')=>-1, ('D','M')=>-3, ('D','F')=>-3, ('D','P')=>-1, ('D','S')=>0, ('D','T')=>-1, ('D','W')=>-4, ('D','Y')=>-3, ('D','V')=>-3,
    ('C','C')=>9, ('C','Q')=>-3, ('C','E')=>-4, ('C','G')=>-3, ('C','H')=>-3, ('C','I')=>-1, ('C','L')=>-1, ('C','K')=>-3, ('C','M')=>-1, ('C','F')=>-2, ('C','P')=>-3, ('C','S')=>-1, ('C','T')=>-1, ('C','W')=>-2, ('C','Y')=>-2, ('C','V')=>-1,
    ('Q','Q')=>5, ('Q','E')=>2, ('Q','G')=>-2, ('Q','H')=>0, ('Q','I')=>-3, ('Q','L')=>-2, ('Q','K')=>1, ('Q','M')=>0, ('Q','F')=>-3, ('Q','P')=>-1, ('Q','S')=>0, ('Q','T')=>-1, ('Q','W')=>-2, ('Q','Y')=>-1, ('Q','V')=>-2,
    ('E','E')=>5, ('E','G')=>-2, ('E','H')=>0, ('E','I')=>-3, ('E','L')=>-3, ('E','K')=>1, ('E','M')=>-2, ('E','F')=>-3, ('E','P')=>-1, ('E','S')=>0, ('E','T')=>-1, ('E','W')=>-3, ('E','Y')=>-2, ('E','V')=>-2,
    ('G','G')=>6, ('G','H')=>-2, ('G','I')=>-4, ('G','L')=>-4, ('G','K')=>-2, ('G','M')=>-3, ('G','F')=>-3, ('G','P')=>-2, ('G','S')=>0, ('G','T')=>-2, ('G','W')=>-2, ('G','Y')=>-3, ('G','V')=>-3,
    ('H','H')=>8, ('H','I')=>-3, ('H','L')=>-3, ('H','K')=>-1, ('H','M')=>-2, ('H','F')=>-1, ('H','P')=>-2, ('H','S')=>-1, ('H','T')=>-2, ('H','W')=>-2, ('H','Y')=>2, ('H','V')=>-3,
    ('I','I')=>4, ('I','L')=>2, ('I','K')=>-3, ('I','M')=>1, ('I','F')=>0, ('I','P')=>-3, ('I','S')=>-2, ('I','T')=>-1, ('I','W')=>-1, ('I','Y')=>-1, ('I','V')=>3,
    ('L','L')=>4, ('L','K')=>-2, ('L','M')=>2, ('L','F')=>0, ('L','P')=>-3, ('L','S')=>-2, ('L','T')=>-1, ('L','W')=>-2, ('L','Y')=>-1, ('L','V')=>1,
    ('K','K')=>5, ('K','M')=>-1, ('K','F')=>-3, ('K','P')=>-1, ('K','S')=>0, ('K','T')=>-1, ('K','W')=>-3, ('K','Y')=>-2, ('K','V')=>-2,
    ('M','M')=>5, ('M','F')=>0, ('M','P')=>-2, ('M','S')=>-1, ('M','T')=>-1, ('M','W')=>-1, ('M','Y')=>-1, ('M','V')=>1,
    ('F','F')=>6, ('F','P')=>-4, ('F','S')=>-2, ('F','T')=>-2, ('F','W')=>1, ('F','Y')=>3, ('F','V')=>0,
    ('P','P')=>7, ('P','S')=>-1, ('P','T')=>-1, ('P','W')=>-4, ('P','Y')=>-3, ('P','V')=>-2,
    ('S','S')=>4, ('S','T')=>1, ('S','W')=>-3, ('S','Y')=>-2, ('S','V')=>-2,
    ('T','T')=>5, ('T','W')=>-2, ('T','Y')=>-2, ('T','V')=>0,
    ('W','W')=>11, ('W','Y')=>2, ('W','V')=>-3,
    ('Y','Y')=>7, ('Y','V')=>-1,
    ('V','V')=>4
)

function _blosum62_dist(c1::Char, c2::Char)
    c1_u, c2_u = uppercase(c1), uppercase(c2)
    c1_u == c2_u && return 0.0
    key = (c1_u <= c2_u) ? (c1_u, c2_u) : (c2_u, c1_u)
    score = get(_BLOSUM62_SCORES, key, -1)
    return max(0.0, 4.0 - Float64(score))
end

"""
    tcrdist_matrix(cdr3_list; v_genes=nothing, threaded=true)

TCRdist distance matrix (Dash et al. 2017 / Steele et al. 2022).
"""
function tcrdist_matrix(
    cdr3_list::AbstractVector{<:AbstractString};
    v_genes::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
    threaded::Bool=true)

    n = length(cdr3_list)
    D = zeros(Float64, n, n)
    v_vec = v_genes !== nothing ? String.(v_genes) : String[]

    function _pair_tcrdist(s1::String, s2::String, v1::String, v2::String)
        dist = 0.0
        if !isempty(v1) && !isempty(v2) && v1 != v2
            dist += 50.0
        end

        l1, l2 = length(s1), length(s2)
        length_penalty = 12.0 * abs(l1 - l2)
        min_len = min(l1, l2)

        aa_dist = 0.0
        for k in 1:min_len
            aa_dist += _blosum62_dist(s1[k], s2[k])
        end
        return 3.0 * aa_dist + length_penalty + dist
    end

    if threaded && Threads.nthreads() > 1
        Threads.@threads for i in 1:n
            v1 = (length(v_vec) == n) ? v_vec[i] : ""
            for j in (i+1):n
                v2 = (length(v_vec) == n) ? v_vec[j] : ""
                d = _pair_tcrdist(String(cdr3_list[i]), String(cdr3_list[j]), v1, v2)
                D[i, j] = d
                D[j, i] = d
            end
        end
    else
        for i in 1:n
            v1 = (length(v_vec) == n) ? v_vec[i] : ""
            for j in (i+1):n
                v2 = (length(v_vec) == n) ? v_vec[j] : ""
                d = _pair_tcrdist(String(cdr3_list[i]), String(cdr3_list[j]), v1, v2)
                D[i, j] = d
                D[j, i] = d
            end
        end
    end

    return D
end

"""
    cdr3_clustering(cdr3_sequences; threshold=2, min_size=2, v_genes=nothing, use_bktree=true)

Fast sequence clustering using BKTree spatial metric indexing for large cohorts (GLIPH / CoNGA style).
"""
function cdr3_clustering(
    cdr3_sequences::AbstractVector{<:AbstractString};
    threshold::Real=2,
    min_size::Int=2,
    v_genes::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
    use_bktree::Bool=true)

    seqs = uppercase.(String.(cdr3_sequences))
    n = length(seqs)
    n == 0 && return (assignments=Int[], clusters=DataFrames.DataFrame())

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

    thresh = Int(floor(threshold))

    if use_bktree && n >= 20
        tree = BKTree(seqs)
        for i in 1:n
            neighbors = bktree_range_query(tree, seqs[i], thresh)
            for j in neighbors
                if j > i
                    if v_genes === nothing || (length(v_genes) == n && v_genes[i] == v_genes[j])
                        union!(parent, i, j)
                    end
                end
            end
        end
    else
        D = cdr3_levenshtein_distance(seqs)
        for i in 1:n
            for j in (i+1):n
                if D[i, j] <= Float64(threshold)
                    if v_genes === nothing || (length(v_genes) == n && v_genes[i] == v_genes[j])
                        union!(parent, i, j)
                    end
                end
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

"""
    clonotype_network(sample; threshold=2, metric=:levenshtein)

Construct network graph of clonotype sequences within a threshold distance.
"""
function clonotype_network(sample::RepertoireSample; threshold::Real=2, metric::Symbol=:levenshtein)
    df = sample.clonotypes
    n = nrow(df)
    n == 0 && return (adjacency=spzeros(0,0), nodes=DataFrames.DataFrame(), edges=DataFrames.DataFrame(), components=DataFrames.DataFrame(), summary=Dict{String,Any}())

    seqs = String.(df.cdr3_aa)
    D = if metric == :hamming
        cdr3_hamming_distance(seqs)
    elseif metric == :tcrdist
        tcrdist_matrix(seqs; v_genes=String.(df.v_gene))
    else
        cdr3_levenshtein_distance(seqs)
    end

    I = Int[]
    J = Int[]
    V = Float64[]

    edge_sources = String[]
    edge_targets = String[]
    edge_dists = Float64[]

    for i in 1:n
        for j in (i+1):n
            if D[i, j] <= Float64(threshold)
                push!(I, i); push!(J, j); push!(V, D[i, j])
                push!(I, j); push!(J, i); push!(V, D[i, j])
                push!(edge_sources, df.clonotype_id[i])
                push!(edge_targets, df.clonotype_id[j])
                push!(edge_dists, D[i, j])
            end
        end
    end

    adj = sparse(I, J, V, n, n)
    degrees = [count(>(0), adj[i, :]) for i in 1:n]

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

    for k in 1:length(edge_sources)
        i = findfirst(==(edge_sources[k]), df.clonotype_id)
        j = findfirst(==(edge_targets[k]), df.clonotype_id)
        (i !== nothing && j !== nothing) && union!(parent, i, j)
    end

    comp_roots = [find!(parent, i) for i in 1:n]
    uniq_comp = sort!(unique(comp_roots))
    comp_map = Dict(r => k for (k, r) in enumerate(uniq_comp))
    comp_ids = [comp_map[r] for r in comp_roots]

    nodes_df = DataFrames.DataFrame(
        clonotype_id=df.clonotype_id,
        cdr3_aa=df.cdr3_aa,
        v_gene=df.v_gene,
        j_gene=df.j_gene,
        clone_size=df.clone_size,
        degree=degrees,
        component_id=comp_ids
    )

    edges_df = DataFrames.DataFrame(source=edge_sources, target=edge_targets, distance=edge_dists)
    comp_df = combine(groupby(nodes_df, :component_id), nrow => :n_nodes, :clone_size => sum => :total_cells)
    sort!(comp_df, :n_nodes, rev=true)

    summary = Dict{String,Any}(
        "sample_id" => sample.sample_id,
        "n_nodes" => n,
        "n_edges" => length(edge_sources),
        "n_components" => nrow(comp_df),
        "max_component_size" => isempty(comp_df) ? 0 : maximum(comp_df.n_nodes),
        "hub_clonotype" => isempty(nodes_df) ? "" : nodes_df.clonotype_id[argmax(nodes_df.degree)]
    )

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "clonotype_network";
            parameters=(sample_id=sample.sample_id, n_nodes=n, n_edges=length(edge_sources)))
    end

    return (adjacency=adj, nodes=nodes_df, edges=edges_df, components=comp_df, summary=summary)
end

# ---- Motif Discovery & Positional Sequence Logos -----------------------------

"""
    gliph2_motif_enrichment(cdr3_list; k_range=3:4, control_list=nothing, p_cutoff=0.05)

Discover conserved amino acid k-mer motifs enriched in CDR3 sequences (GLIPH2 style).
"""
function gliph2_motif_enrichment(
    cdr3_list::AbstractVector{<:AbstractString};
    k_range::AbstractVector{Int}=3:4,
    control_list::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
    p_cutoff::Float64=0.05)

    sample_seqs = uppercase.(String.(cdr3_list))
    n_sample = length(sample_seqs)

    motif_sample_counts = Dict{String,Int}()
    for seq in sample_seqs
        L = length(seq)
        seen_in_seq = Set{String}()
        for k in k_range
            for i in 1:(L - k + 1)
                push!(seen_in_seq, seq[i:(i + k - 1)])
            end
        end
        for m in seen_in_seq
            motif_sample_counts[m] = get(motif_sample_counts, m, 0) + 1
        end
    end

    motifs = String[]
    sample_counts = Int[]
    sample_freqs = Float64[]
    pvalues = Float64[]

    for (m, count_s) in motif_sample_counts
        freq_s = count_s / max(n_sample, 1)

        # Right-tailed Binomial / Fisher test vs expected background (default 0.01 per kmer)
        expected_p = 0.01 ^ length(m)
        pval = 1.0 - cdf(Binomial(n_sample, clamp(expected_p, 1e-6, 0.5)), count_s - 1)

        if pval <= p_cutoff
            push!(motifs, m)
            push!(sample_counts, count_s)
            push!(sample_freqs, freq_s)
            push!(pvalues, pval)
        end
    end

    res = DataFrames.DataFrame(motif=motifs, count=sample_counts, frequency=sample_freqs, pvalue=pvalues)
    sort!(res, :pvalue)
    return res
end

"""
    positional_motif_logo(cdr3_list; method=:information_content)

Compute Position Weight Matrix (PWM) and positional information content (bits) for CDR3s.
"""
function positional_motif_logo(cdr3_list::AbstractVector{<:AbstractString}; method::Symbol=:information_content)
    seqs = uppercase.(String.(cdr3_list))
    filter!(!isempty, seqs)
    isempty(seqs) && return DataFrames.DataFrame()

    mode_len = mode(length.(seqs))
    fixed_seqs = filter(s -> length(s) == mode_len, seqs)
    N = length(fixed_seqs)
    N == 0 && return DataFrames.DataFrame()

    amino_acids = collect("ACDEFGHIKLMNPQRSTVWY")
    n_aa = length(amino_acids)

    pwm = zeros(Float64, mode_len, n_aa)
    for s in fixed_seqs
        for (i, char) in enumerate(s)
            idx = findfirst(==(char), amino_acids)
            if idx !== nothing
                pwm[i, idx] += 1.0
            end
        end
    end
    pwm ./= N

    info_content = Float64[]
    for i in 1:mode_len
        pi = pwm[i, :]
        h = -sum(p > 0 ? p * log2(p) : 0.0 for p in pi)
        R_i = log2(20) - h
        push!(info_content, max(R_i, 0.0))
    end

    df = DataFrames.DataFrame(position=1:mode_len, information_bits=info_content)
    for (j, aa) in enumerate(amino_acids)
        df[!, Symbol("freq_$(aa)")] = pwm[:, j]
    end
    return df
end

# ---- Biophysical CDR3 Properties & Atchley/Kidera Factors --------------------

const _KYTE_DOOLITTLE = Dict{Char,Float64}(
    'A'=>1.8, 'R'=>-4.5, 'N'=>-3.5, 'D'=>-3.5, 'C'=>2.5, 'Q'=>-3.5, 'E'=>-3.5, 'G'=>-0.4, 'H'=>-3.2, 'I'=>4.5,
    'L'=>3.8, 'K'=>-3.9, 'M'=>1.9, 'F'=>2.8, 'P'=>-1.6, 'S'=>-0.8, 'T'=>-0.7, 'W'=>-0.9, 'Y'=>-1.3, 'V'=>4.2
)

# Atchley Factors (Atchley et al. 2005, PNAS): Polarity, SecStruct, Size, CodonComposition, Charge
const _ATCHLEY_FACTORS = Dict{Char, Vector{Float64}}(
    'A' => [-0.591, -1.302, -0.733, 1.570, -0.146],
    'R' => [1.538, -0.055, 1.502, 0.440, 2.897],
    'N' => [0.945, 0.828, 1.299, -0.169, 0.933],
    'D' => [1.050, 0.302, -0.365, -1.238, -0.259],
    'C' => [-0.108, 1.522, 1.719, -0.551, -0.718],
    'Q' => [0.931, -0.179, -0.300, -0.380, -0.134],
    'E' => [1.357, -1.453, 1.477, 0.113, -0.837],
    'G' => [-0.384, 1.652, 1.330, 1.045, 2.064],
    'H' => [0.336, -0.417, -1.673, -1.474, -0.078],
    'I' => [-1.239, -0.547, 2.131, 0.393, 0.816],
    'L' => [-1.019, -0.987, -1.505, 1.266, -0.912],
    'K' => [1.831, -0.561, 0.533, -0.277, 1.648],
    'M' => [-0.663, -1.524, 2.219, -1.005, 1.212],
    'F' => [-1.006, -0.590, 1.891, -0.397, 0.412],
    'P' => [0.189, 2.081, -1.628, 0.421, -1.392],
    'S' => [0.016, 0.870, -0.814, 0.238, -0.413],
    'T' => [-0.188, 0.170, -0.468, -0.267, -0.474],
    'W' => [-0.595, 0.009, 0.672, -2.128, -0.184],
    'Y' => [-0.260, 0.246, 0.620, -1.099, -0.329],
    'V' => [-1.337, -0.279, -0.544, 1.242, -0.126]
)

"""
    repertoire_cdr3_physicochemical_properties(cdr3_list)

Compute biophysical profiles of CDR3 sequences including Kyte-Doolittle GRAVY, Charge, pI, and Atchley Factors (1-5).
"""
function repertoire_cdr3_physicochemical_properties(cdr3_list::AbstractVector{<:AbstractString})
    seqs = uppercase.(String.(cdr3_list))

    len_vec = Int[]
    gravy_vec = Float64[]
    charge_vec = Float64[]
    aliphatic_vec = Float64[]
    aromatic_vec = Float64[]

    a1_vec = Float64[]
    a2_vec = Float64[]
    a3_vec = Float64[]
    a4_vec = Float64[]
    a5_vec = Float64[]

    for s in seqs
        L = length(s)
        push!(len_vec, L)
        if L == 0
            push!(gravy_vec, 0.0); push!(charge_vec, 0.0)
            push!(aliphatic_vec, 0.0); push!(aromatic_vec, 0.0)
            push!(a1_vec, 0.0); push!(a2_vec, 0.0); push!(a3_vec, 0.0); push!(a4_vec, 0.0); push!(a5_vec, 0.0)
            continue
        end

        hydros = [get(_KYTE_DOOLITTLE, c, 0.0) for c in s]
        push!(gravy_vec, mean(hydros))

        # Net charge at pH 7.4
        pos = count(c -> c in ('K', 'R', 'H'), s)
        neg = count(c -> c in ('D', 'E'), s)
        push!(charge_vec, Float64(pos - neg))

        # Aliphatic index & Aromaticity
        ala = count(==('A'), s); val = count(==('V'), s); ile = count(==('I'), s); leu = count(==('L'), s)
        push!(aliphatic_vec, (ala + 2.9*val + 3.9*(ile + leu)) / L)

        aro = count(c -> c in ('F', 'W', 'Y'), s)
        push!(aromatic_vec, aro / L)

        # Atchley Factors
        af_vals = [get(_ATCHLEY_FACTORS, c, zeros(5)) for c in s]
        mean_af = mean(af_vals)
        push!(a1_vec, mean_af[1]); push!(a2_vec, mean_af[2]); push!(a3_vec, mean_af[3]); push!(a4_vec, mean_af[4]); push!(a5_vec, mean_af[5])
    end

    return DataFrames.DataFrame(
        cdr3_aa=seqs,
        length=len_vec,
        gravy_hydrophobicity=gravy_vec,
        net_charge=charge_vec,
        aliphatic_index=aliphatic_vec,
        aromaticity=aromatic_vec,
        atchley_polarity=a1_vec,
        atchley_sec_struct=a2_vec,
        atchley_size=a3_vec,
        atchley_codon=a4_vec,
        atchley_charge=a5_vec
    )
end

# ---- Spectratype & Lineage Tree Analytics ------------------------------------

"""
    spectratype_analysis(sample; by=:length)

Generate CDR3 length spectratype distribution and compute spectratype perturbation entropy.
"""
function spectratype_analysis(sample::RepertoireSample; by::Symbol=:length)
    df = sample.clonotypes
    lengths = Int.(df.cdr3_length)
    sizes = Int.(df.clone_size)

    length_counts = Dict{Int,Int}()
    for (l, sz) in zip(lengths, sizes)
        length_counts[l] = get(length_counts, l, 0) + sz
    end

    uniq_lens = sort!(collect(keys(length_counts)))
    counts = [length_counts[l] for l in uniq_lens]
    total = sum(counts)
    freqs = total > 0 ? counts ./ total : zeros(Float64, length(counts))

    shannon_spectratype = -sum(p > 0 ? p * log(p) : 0.0 for p in freqs)

    res = DataFrames.DataFrame(cdr3_length=uniq_lens, cell_count=counts, frequency=freqs)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "spectratype_analysis";
            parameters=(sample_id=sample.sample_id, spectratype_entropy=shannon_spectratype))
    end
    return (spectratype=res, spectratype_entropy=shannon_spectratype)
end

"""
    clonotype_lineage_tree(sample, clonotype_id)

Construct Minimum Spanning Tree (MST) for intraclonal sequence variants within a clonotype.
"""
function clonotype_lineage_tree(sample::RepertoireSample, clonotype_id::AbstractString)
    assigns = filter(a -> a.clonotype_id == clonotype_id, sample.assignments)
    isempty(assigns) && return (tree_nodes=DataFrames.DataFrame(), tree_edges=DataFrames.DataFrame())

    unique_seqs = sort!(unique([a.cdr3_aa for a in assigns]))
    n = length(unique_seqs)
    n == 0 && return (tree_nodes=DataFrames.DataFrame(), tree_edges=DataFrames.DataFrame())

    D = cdr3_levenshtein_distance(unique_seqs)

    # Prim's Algorithm for Minimum Spanning Tree
    in_tree = fill(false, n)
    in_tree[1] = true

    edge_sources = String[]
    edge_targets = String[]
    edge_weights = Float64[]

    for step in 1:(n - 1)
        min_d = Inf
        min_i = 1
        min_j = 1
        for i in 1:n
            if in_tree[i]
                for j in 1:n
                    if !in_tree[j] && D[i, j] < min_d
                        min_d = D[i, j]
                        min_i = i
                        min_j = j
                    end
                end
            end
        end
        if isfinite(min_d)
            in_tree[min_j] = true
            push!(edge_sources, unique_seqs[min_i])
            push!(edge_targets, unique_seqs[min_j])
            push!(edge_weights, min_d)
        end
    end

    nodes_df = DataFrames.DataFrame(sequence=unique_seqs, is_root=[i == 1 for i in 1:n])
    edges_df = DataFrames.DataFrame(source=edge_sources, target=edge_targets, distance=edge_weights)
    return (tree_nodes=nodes_df, tree_edges=edges_df)
end

# ---- Public & Convergent Clonotypes -----------------------------------------

"""
    convergent_clonotypes(samples; max_hamming=1, min_shared=2)

Identify public clonotypes shared across multiple samples.
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

"""
    public_clonotype_analysis(samples; min_samples=2)

Generate cohort public clonotype sharing matrix and prevalence summary.
"""
function public_clonotype_analysis(samples::AbstractVector{<:RepertoireSample}; min_samples::Int=2)
    n_samples = length(samples)
    sample_ids = [s.sample_id for s in samples]

    all_clonotypes = Set{String}()
    for s in samples
        for id in s.clonotypes.clonotype_id
            push!(all_clonotypes, id)
        end
    end

    clonotype_freqs = Dict{String, Vector{Float64}}()
    for ct in all_clonotypes
        clonotype_freqs[ct] = zeros(Float64, n_samples)
    end

    for (s_idx, s) in enumerate(samples)
        dict_freq = Dict(zip(s.clonotypes.clonotype_id, s.clonotypes.frequency))
        for ct in keys(clonotype_freqs)
            clonotype_freqs[ct][s_idx] = get(dict_freq, ct, 0.0)
        end
    end

    res = DataFrames.DataFrame(
        clonotype_id=String[],
        prevalence=Int[],
        prevalence_pct=Float64[],
        mean_frequency=Float64[],
        max_frequency=Float64[]
    )

    for (ct, freqs) in clonotype_freqs
        prev = count(>(0.0), freqs)
        if prev >= min_samples
            push!(res, (ct, prev, (prev / n_samples) * 100.0, mean(freqs), maximum(freqs)))
        end
    end

    sort!(res, :prevalence, rev=true)
    return res
end

# ---- V(D)J Usage & Statistical Tests -----------------------------------------

"""
    vj_usage_matrix(sample)

Compute V and J gene usage tables and frequencies.
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

Compute paired V-J gene usage matrix and relative frequencies.
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

"""
    vj_usage_divergence(sample_a, sample_b; metric=:jsd)

Compute Jensen-Shannon Divergence (JSD) or Kullback-Leibler (KL) divergence between V/J distributions.
"""
function vj_usage_divergence(sample_a::RepertoireSample, sample_b::RepertoireSample; metric::Symbol=:jsd)
    v_a = vj_usage_matrix(sample_a).v_usage
    v_b = vj_usage_matrix(sample_b).v_usage

    all_v = sort!(union(v_a.v_gene, v_b.v_gene))
    dict_a = Dict(zip(v_a.v_gene, v_a.frequency))
    dict_b = Dict(zip(v_b.v_gene, v_b.frequency))

    p = Float64[get(dict_a, g, 1e-10) for g in all_v]
    q = Float64[get(dict_b, g, 1e-10) for g in all_v]
    p ./= sum(p)
    q ./= sum(q)

    if metric == :kld
        return sum(p[i] * log(p[i] / q[i]) for i in 1:length(p))
    else
        m = 0.5 .* (p .+ q)
        kl_pm = sum(p[i] * log(p[i] / m[i]) for i in 1:length(p))
        kl_qm = sum(q[i] * log(q[i] / m[i]) for i in 1:length(q))
        return 0.5 * kl_pm + 0.5 * kl_qm
    end
end

"""
    differential_vj_usage(samples_group1, samples_group2)

Statistical test comparing V/J gene counts between two condition groups using Welch's t-test and FDR control.
"""
function differential_vj_usage(group1::AbstractVector{<:RepertoireSample}, group2::AbstractVector{<:RepertoireSample})
    n1 = length(group1)
    n2 = length(group2)

    all_v = Set{String}()
    for s in vcat(group1, group2)
        for g in vj_usage_matrix(s).v_usage.v_gene
            push!(all_v, g)
        end
    end

    genes = sort!(collect(all_v))
    res = DataFrames.DataFrame(
        v_gene=String[],
        mean_freq_g1=Float64[],
        mean_freq_g2=Float64[],
        fold_change=Float64[],
        pvalue=Float64[]
    )

    for g in genes
        freqs1 = Float64[]
        freqs2 = Float64[]
        for s in group1
            v_u = vj_usage_matrix(s).v_usage
            idx = findfirst(==(g), v_u.v_gene)
            push!(freqs1, idx !== nothing ? v_u.frequency[idx] : 0.0)
        end
        for s in group2
            v_u = vj_usage_matrix(s).v_usage
            idx = findfirst(==(g), v_u.v_gene)
            push!(freqs2, idx !== nothing ? v_u.frequency[idx] : 0.0)
        end

        m1, m2 = mean(freqs1), mean(freqs2)
        fc = (m1 + 1e-6) / (m2 + 1e-6)

        s1, s2 = var(freqs1), var(freqs2)
        se = sqrt(s1/max(n1,1) + s2/max(n2,1))
        t_stat = se > 0 ? (m1 - m2) / se : 0.0
        df_welch = se > 0 ? (s1/n1 + s2/n2)^2 / ((s1/n1)^2/(n1-1 + 1e-5) + (s2/n2)^2/(n2-1 + 1e-5)) : 1.0
        pval = se > 0 ? 2.0 * (1.0 - cdf(TDist(max(df_welch, 1.0)), abs(t_stat))) : 1.0

        push!(res, (g, m1, m2, fc, pval))
    end

    m = nrow(res)
    padj = similar(res.pvalue)
    perm = sortperm(res.pvalue)
    running = 1.0
    for i in m:-1:1
        j = perm[i]
        running = min(running, res.pvalue[j] * m / i)
        padj[j] = min(running, 1.0)
    end
    res[!, :padj] = padj

    sort!(res, :pvalue)
    return res
end

"""
    expanded_clonotype_test(sample; threshold=5)

Poisson statistical test identifying significantly expanded clonotypes.
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

# ---- BCR-Specific Analytics & Hotspot Profiling -----------------------------

"""
    bcr_shm_rate(contigs::Vector{RepertoireContig})

Compute Somatic Hypermutation (SHM) rate density across BCR V regions.
"""
function bcr_shm_rate(contigs::AbstractVector{<:RepertoireContig})
    bcr_contigs = filter(c -> c.chain in (:IGH, :IGL, :IGK), contigs)
    n = length(bcr_contigs)
    n == 0 && return (mean_shm=0.0, shm_df=DataFrames.DataFrame())

    shm_rates = Float64[]
    for c in bcr_contigs
        shm_val = get(c.annotations, "shm_rate", get(c.annotations, "v_identity", 1.0))
        rate = shm_val isa Real ? (shm_val <= 1.0 ? (1.0 - Float64(shm_val)) * 100.0 : Float64(shm_val)) : 0.0
        push!(shm_rates, rate)
    end

    df = DataFrames.DataFrame(
        contig_id=[c.contig_id for c in bcr_contigs],
        barcode=[c.barcode for c in bcr_contigs],
        chain=[c.chain for c in bcr_contigs],
        v_gene=[c.v_gene for c in bcr_contigs],
        shm_rate_pct=shm_rates
    )

    return (mean_shm=mean(shm_rates), shm_df=df)
end

"""
    bcr_mutation_hotspot_profile(sequences)

Analyze AID mutation hotspots (WRC/GYW, WA/TW) vs coldspots (SYC/GRS) in BCR sequences (SHazaM style).
"""
function bcr_mutation_hotspot_profile(sequences::AbstractVector{<:AbstractString})
    seqs = uppercase.(String.(sequences))
    n = length(seqs)

    aid_hotspots = [r"WRC", r"GYW", r"WA", r"TW"]
    coldspots = [r"SYC", r"GRS"]

    hot_counts = Int[]
    cold_counts = Int[]
    ratios = Float64[]

    for s in seqs
        L = length(s)
        n_hot = 0
        for pat in (r"[AT][AG]C", r"C[CT][AT]", r"[AT]A", r"T[AT]")
            pos = 1
            while pos <= L
                m = match(pat, s, pos)
                m === nothing && break
                n_hot += 1
                pos = m.offset + 1
            end
        end
        n_cold = 0
        for pat in (r"[CG][CT]C", r"G[AG][CG]")
            pos = 1
            while pos <= L
                m = match(pat, s, pos)
                m === nothing && break
                n_cold += 1
                pos = m.offset + 1
            end
        end

        push!(hot_counts, n_hot)
        push!(cold_counts, n_cold)
        push!(ratios, (n_hot + 0.1) / (n_cold + 0.1))
    end

    return DataFrames.DataFrame(
        sequence=seqs,
        hotspot_count=hot_counts,
        coldspot_count=cold_counts,
        hotspot_to_coldspot_ratio=ratios
    )
end

"""
    isotype_usage_summary(contigs::Vector{RepertoireContig})

Summarize BCR Heavy chain constant region isotype frequencies (IgM, IgD, IgG, IgA, IgE).
"""
function isotype_usage_summary(contigs::AbstractVector{<:RepertoireContig})
    igh_contigs = filter(c -> c.chain == :IGH, contigs)
    n = length(igh_contigs)
    n == 0 && return DataFrames.DataFrame(isotype=String[], count=Int[], frequency=Float64[])

    isotype_counts = Dict{String,Int}()
    for c in igh_contigs
        iso = if occursin(r"IGHM|IgM"i, c.c_gene)
            "IgM"
        elseif occursin(r"IGHD|IgD"i, c.c_gene)
            "IgD"
        elseif occursin(r"IGHG1|IgG1"i, c.c_gene)
            "IgG1"
        elseif occursin(r"IGHG2|IgG2"i, c.c_gene)
            "IgG2"
        elseif occursin(r"IGHG3|IgG3"i, c.c_gene)
            "IgG3"
        elseif occursin(r"IGHG4|IgG4"i, c.c_gene)
            "IgG4"
        elseif occursin(r"IGHG|IgG"i, c.c_gene)
            "IgG"
        elseif occursin(r"IGHA|IgA"i, c.c_gene)
            "IgA"
        elseif occursin(r"IGHE|IgE"i, c.c_gene)
            "IgE"
        else
            "Other/Unassigned"
        end
        isotype_counts[iso] = get(isotype_counts, iso, 0) + 1
    end

    df = DataFrames.DataFrame(
        isotype=collect(keys(isotype_counts)),
        count=collect(values(isotype_counts))
    )
    df[!, :frequency] = df.count ./ n
    sort!(df, :count, rev=true)
    return df
end

# ---- Visualization Helpers ---------------------------------------------------

"""
    repertoire_heatmap(samples; metric=:jaccard)

Generate long-format DataFrame (sample1, sample2, similarity) formatted for heatmaps.
"""
function repertoire_heatmap(samples::AbstractVector{<:RepertoireSample}; metric::Symbol=:jaccard)
    n = length(samples)
    names = [s.sample_id for s in samples]
    mat = clonotype_overlap_matrix(samples; metric=metric)

    df = DataFrames.DataFrame(sample1=String[], sample2=String[], value=Float64[])
    for i in 1:n
        for j in 1:n
            push!(df, (names[i], names[j], mat[i, j]))
        end
    end

    return df
end

"""
    clonotype_venn(samples)

Compute set intersections for 2-way, 3-way, or N-way Venn/Upset diagrams.
"""
function clonotype_venn(samples::AbstractVector{<:RepertoireSample})
    n = length(samples)
    ids = [s.sample_id for s in samples]
    sets = [Set{String}(s.clonotypes.clonotype_id) for s in samples]

    df = DataFrames.DataFrame(set_combination=String[], count=Int[], clonotypes=String[])

    for k in 1:((1 << n) - 1)
        active_indices = [i for i in 1:n if (k & (1 << (i - 1))) != 0]
        comb_name = join([ids[i] for i in active_indices], " & ")

        inter = copy(sets[active_indices[1]])
        for i in active_indices[2:end]
            intersect!(inter, sets[i])
        end

        inactive_indices = [i for i in 1:n if (k & (1 << (i - 1))) == 0]
        for i in inactive_indices
            setdiff!(inter, sets[i])
        end

        cts = collect(inter)
        push!(df, (comb_name, length(cts), join(cts[1:min(10, length(cts))], ", ")))
    end

    sort!(df, :count, rev=true)
    return df
end

end
