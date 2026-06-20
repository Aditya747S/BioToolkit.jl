module LongRead

using DataFrames
using LinearAlgebra
using Statistics
using Random

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, BamCigarOp, BamFile, BamRecord, BioAlphabet, BioSequence, DNAAlphabet, DNASeq, FastqRecord, IntervalTree, query_overlaps, reverse_complement, threaded_foreach
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!

@inline function _register_longread_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export NanoporeRead, PacBioRead, StructuralVariant
export MinimizerSeed, CandidateRegion, BandedAlignmentResult
export canonical_kmer_hash, minimizer_sketch, build_minimizer_index, find_candidate_regions
export banded_semiglobal_alignment
export split_read_evidence, discordant_pair_evidence, cluster_structural_variants, call_structural_variants
export isoseq_consensus, minimap2_align, graphaligner_align, sniffles_call, svim_call, whatshap_phase, phase_reads_by_alleles
export overlap_layout_consensus, read_n50
export overlap_graph_table
export porec_contact_table, omnic_contact_matrix

struct NanoporeRead{A<:BioAlphabet}
    sequence::BioSequence{A}
    identifier::String
    description::String
    quality::String
    signal::Vector{Float32}
    channel_info::Dict{String,Any}
end

struct PacBioRead{A<:BioAlphabet}
    sequence::BioSequence{A}
    identifier::String
    description::String
    quality::String
    signal::Vector{Float32}
    channel_info::Dict{String,Any}
end

struct StructuralVariant
    chrom1::String
    pos1::Int
    chrom2::String
    pos2::Int
    sv_type::Symbol
    supporting_reads::Vector{String}
    genotype::String
end

struct MinimizerSeed
    hash::UInt64
    position::Int
end

struct CandidateRegion
    chrom::String
    start::Int
    stop::Int
    support::Int
end

struct BandedAlignmentResult <: AbstractAnalysisResult
    score::Int
    query_start::Int
    query_end::Int
    ref_start::Int
    ref_end::Int
    cigar::Vector{BamCigarOp}
    aligned_query::String
    aligned_reference::String
    provenance::ResultProvenance
end

BandedAlignmentResult(score, query_start, query_end, ref_start, ref_end, cigar, aligned_query, aligned_reference) =
    BandedAlignmentResult(score, query_start, query_end, ref_start, ref_end, cigar, aligned_query, aligned_reference, provenance_record("BandedAlignmentResult", "longread"))

struct _SVEvidence
    chrom1::String
    pos1::Int
    chrom2::String
    pos2::Int
    sv_type::Symbol
    read_name::String
end

mutable struct _SVCluster
    chrom1::String
    chrom2::String
    sv_type::Symbol
    pos1_values::Vector{Int}
    pos2_values::Vector{Int}
    reads::Vector{String}
end

@inline function _record_sequence(record)
    if record isa NanoporeRead || record isa PacBioRead
        return String(record.sequence)
    elseif record isa FastqRecord
        return String(record.sequence)
    elseif record isa BioSequence
        return String(record)
    else
        return String(record)
    end
end

@inline function _base_code(byte::UInt8)
    if byte == 0x41 || byte == 0x61
        return UInt8(0)
    elseif byte == 0x43 || byte == 0x63
        return UInt8(1)
    elseif byte == 0x47 || byte == 0x67
        return UInt8(2)
    elseif byte == 0x54 || byte == 0x74
        return UInt8(3)
    end
    return UInt8(4)
end

function NanoporeRead(record::FastqRecord{A}; signal::AbstractVector{<:Real}=Float32[], channel_info::AbstractDict=Dict{String,Any}()) where {A<:BioAlphabet}
    return NanoporeRead{A}(record.sequence, record.identifier, record.description, record.quality, Float32.(signal), Dict{String,Any}(string(key) => value for (key, value) in channel_info))
end

function PacBioRead(record::FastqRecord{A}; signal::AbstractVector{<:Real}=Float32[], channel_info::AbstractDict=Dict{String,Any}()) where {A<:BioAlphabet}
    return PacBioRead{A}(record.sequence, record.identifier, record.description, record.quality, Float32.(signal), Dict{String,Any}(string(key) => value for (key, value) in channel_info))
end

function canonical_kmer_hash(kmer::BioSequence{DNAAlphabet})
    text = uppercase(String(kmer))
    rc = uppercase(reverse_complement(text))

    function _hash(seq::String)
        value = UInt64(0xcbf29ce484222325)
        for byte in codeunits(seq)
            value = (value ⊻ UInt64(_base_code(byte) + 1)) * UInt64(0x100000001b3)
        end
        return value
    end

    forward = _hash(text)
    reverse = _hash(rc)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, min(forward, reverse), "canonical_kmer_hash")
end

function canonical_kmer_hash(kmer::String)
    text = uppercase(String(kmer))
    rc = uppercase(reverse_complement(text))

    function _hash(seq::String)
        value = UInt64(0xcbf29ce484222325)
        for byte in codeunits(seq)
            value = (value ⊻ UInt64(_base_code(byte) + 1)) * UInt64(0x100000001b3)
        end
        return value
    end

    forward = _hash(text)
    reverse = _hash(rc)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, min(forward, reverse), "canonical_kmer_hash")
end

function minimizer_sketch(sequence; k::Integer=15, w::Integer=10)
    k > 0 || throw(ArgumentError("k must be positive"))
    w > 0 || throw(ArgumentError("w must be positive"))

    text = uppercase(_record_sequence(sequence))
    n = ncodeunits(text)
    n < k && return MinimizerSeed[]

    n_kmers = n - k + 1
    hashes = Vector{UInt64}(undef, n_kmers)
    @inbounds for index in eachindex(hashes)
        hashes[index] = canonical_kmer_hash(String(SubString(text, index, index + k - 1)))
    end

    window_width = min(Int(w), n_kmers)
    minimizers = MinimizerSeed[]
    last_position = 0

    for start in 1:(n_kmers - window_width + 1)
        stop = start + window_width - 1
        window = @view hashes[start:stop]
        local_offset = argmin(window)
        position = start + local_offset - 1
        if position != last_position
            push!(minimizers, MinimizerSeed(hashes[position], position))
            last_position = position
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, minimizers, "minimizer_sketch")
end

function build_minimizer_index(reference::AbstractDict{<:AbstractString,<:Any}; k::Integer=15, w::Integer=10)
    index = Dict{UInt64,Vector{Tuple{String,Int}}}()
    for (chromosome, sequence) in reference
        for seed in minimizer_sketch(sequence; k=k, w=w)
            push!(get!(index, seed.hash, Tuple{String,Int}[]), (String(chromosome), seed.position))
        end
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, index, "build_minimizer_index")
end

function build_minimizer_index(reference_sequence; chrom::AbstractString="ref", k::Integer=15, w::Integer=10)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, build_minimizer_index(Dict(String(chrom) => _record_sequence(reference_sequence)); k=k, w=w), "build_minimizer_index")
end

function find_candidate_regions(read, index::Dict{UInt64,Vector{Tuple{String,Int}}}; k::Integer=15, w::Integer=10, max_hits_per_minimizer::Integer=100)
    read_text = _record_sequence(read)
    votes = Dict{Tuple{String,Int},Int}()
    read_length = ncodeunits(read_text)

    for seed in minimizer_sketch(read_text; k=k, w=w)
        hits = get(index, seed.hash, Tuple{String,Int}[])
        if length(hits) > max_hits_per_minimizer
            hits = @view hits[1:max_hits_per_minimizer]
        end
        for hit in hits
            chromosome, hit_position = hit
            offset = hit_position - seed.position
            key = (chromosome, offset)
            votes[key] = get(votes, key, 0) + 1
        end
    end

    candidates = CandidateRegion[]
    for ((chromosome, offset), support) in votes
        start = max(1, offset + 1)
        stop = start + max(read_length - 1, 0)
        push!(candidates, CandidateRegion(chromosome, start, stop, support))
    end

    sort!(candidates; by=region -> (-region.support, region.chrom, region.start))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, candidates, "find_candidate_regions")
end

function _alignment_to_cigar(aligned_query::String, aligned_reference::String)
    ncodeunits(aligned_query) == ncodeunits(aligned_reference) || throw(ArgumentError("aligned strings must have equal length"))
    isempty(aligned_query) && return BamCigarOp[]

    cigar = BamCigarOp[]
    current_op = '\0'
    current_len = 0

    for index in eachindex(aligned_query)
        q = aligned_query[index]
        r = aligned_reference[index]
        op = if q == '-'
            'D'
        elseif r == '-'
            'I'
        else
            'M'
        end

        if op == current_op
            current_len += 1
        else
            if current_len > 0
                push!(cigar, BamCigarOp(current_len, current_op))
            end
            current_op = op
            current_len = 1
        end
    end

    current_len > 0 && push!(cigar, BamCigarOp(current_len, current_op))
    return cigar
end

function banded_semiglobal_alignment(query, reference; bandwidth::Union{Nothing,Integer}=nothing, match_score::Integer=2, mismatch_score::Integer=-3, gap_score::Integer=-4)
    q = uppercase(_record_sequence(query))
    r = uppercase(_record_sequence(reference))
    n = ncodeunits(q)
    m = ncodeunits(r)
    bw = bandwidth === nothing ? max(n, m) : Int(bandwidth)
    bw >= 1 || throw(ArgumentError("bandwidth must be at least one"))

    neg_inf = -10^9
    dp = fill(neg_inf, n + 1, m + 1)
    trace = fill(UInt8(0), n + 1, m + 1)

    for j in 0:m
        dp[1, j + 1] = 0
        trace[1, j + 1] = UInt8(3)
    end
    for i in 1:n
        dp[i + 1, 1] = i * gap_score
        trace[i + 1, 1] = UInt8(2)
    end

    for i in 1:n
        j_lo = max(1, i - bw)
        j_hi = min(m, i + bw)
        for j in j_lo:j_hi
            score_diag = dp[i, j] + (q[i] == r[j] ? match_score : mismatch_score)
            score_up = dp[i, j + 1] + gap_score
            score_left = dp[i + 1, j] + gap_score

            best = score_diag
            direction = UInt8(1)
            if score_up > best
                best = score_up
                direction = UInt8(2)
            end
            if score_left > best
                best = score_left
                direction = UInt8(3)
            end

            dp[i + 1, j + 1] = best
            trace[i + 1, j + 1] = direction
        end
    end

    best_col = 0
    best_score = neg_inf
    for j in 0:m
        score = dp[n + 1, j + 1]
        if score > best_score
            best_score = score
            best_col = j
        end
    end

    aligned_query = Vector{UInt8}()
    aligned_reference = Vector{UInt8}()
    i = n
    j = best_col

    while i > 0
        if j == 0
            push!(aligned_query, UInt8(q[i]))
            push!(aligned_reference, UInt8('-'))
            i -= 1
            continue
        end

        direction = trace[i + 1, j + 1]
        if direction == UInt8(1)
            push!(aligned_query, UInt8(q[i]))
            push!(aligned_reference, UInt8(r[j]))
            i -= 1
            j -= 1
        elseif direction == UInt8(2)
            push!(aligned_query, UInt8(q[i]))
            push!(aligned_reference, UInt8('-'))
            i -= 1
        else
            push!(aligned_query, UInt8('-'))
            push!(aligned_reference, UInt8(r[j]))
            j -= 1
        end
    end

    reverse!(aligned_query)
    reverse!(aligned_reference)
    query_start = 1
    query_end = n
    ref_start = j + 1
    ref_end = best_col

    aligned_query_text = String(aligned_query)
    aligned_reference_text = String(aligned_reference)
    cigar = _alignment_to_cigar(aligned_query_text, aligned_reference_text)

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, BandedAlignmentResult(best_score, query_start, query_end, ref_start, ref_end, cigar, aligned_query_text, aligned_reference_text), "banded_semiglobal_alignment")
end

function _make_evidence(record::BamRecord, sv_type::Symbol, pos1::Int, pos2::Int; chrom2::Union{Nothing,String}=nothing)
    chrom1 = record.refname === nothing ? "*" : record.refname
    partner = chrom2 === nothing ? chrom1 : String(chrom2)
    return _SVEvidence(chrom1, max(1, pos1), partner, max(1, pos2), sv_type, record.qname)
end

function split_read_evidence(record::BamRecord; min_sv_size::Integer=30, min_softclip::Integer=20)
    events = _SVEvidence[]
    record.refname === nothing && return events

    ref_pos = Int(record.pos) + 1
    query_pos = 1

    for op in record.cigar
        if op.op == 'M' || op.op == '=' || op.op == 'X'
            ref_pos += op.length
            query_pos += op.length
        elseif op.op == 'I'
            if op.length >= min_sv_size
                push!(events, _make_evidence(record, :DUP, ref_pos, ref_pos + op.length))
            end
            query_pos += op.length
        elseif op.op == 'D' || op.op == 'N'
            if op.length >= min_sv_size
                push!(events, _make_evidence(record, :DEL, ref_pos, ref_pos + op.length))
            end
            ref_pos += op.length
        elseif op.op == 'S'
            if op.length >= min_softclip
                partner_pos = record.mate_pos >= 0 ? Int(record.mate_pos) + 1 : ref_pos
                partner_chrom = record.mate_refname === nothing ? record.refname : record.mate_refname
                push!(events, _make_evidence(record, :TRA, ref_pos, partner_pos; chrom2=partner_chrom))
            end
            query_pos += op.length
        elseif op.op == 'H' || op.op == 'P'
            continue
        else
            ref_pos += op.length
            query_pos += op.length
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, events, "split_read_evidence")
end

@inline _flag_is_set(flag::UInt16, bit::UInt16) = (flag & bit) != 0

function discordant_pair_evidence(record::BamRecord; insert_size_mean::Real=500.0, insert_size_std::Real=100.0, z_threshold::Real=3.0, min_sv_size::Integer=50)
    events = _SVEvidence[]
    record.refname === nothing && return events
    record.mate_refname === nothing && return events
    record.mate_pos < 0 && return events

    read_pos = Int(record.pos) + 1
    mate_pos = Int(record.mate_pos) + 1

    if record.refname != record.mate_refname
        push!(events, _make_evidence(record, :TRA, read_pos, mate_pos; chrom2=record.mate_refname))
        return events
    end

    insert_size = abs(Int(record.template_length))
    discordant_cutoff = insert_size_mean + z_threshold * insert_size_std

    if insert_size >= max(Int(round(discordant_cutoff)), min_sv_size)
        left = min(read_pos, mate_pos)
        right = max(read_pos, mate_pos)
        sv_type = Int(record.template_length) >= 0 ? :DEL : :DUP
        push!(events, _make_evidence(record, sv_type, left, right))
    end

    read_reverse = _flag_is_set(record.flag, UInt16(0x10))
    mate_reverse = _flag_is_set(record.flag, UInt16(0x20))
    if read_reverse == mate_reverse
        left = min(read_pos, mate_pos)
        right = max(read_pos, mate_pos)
        push!(events, _make_evidence(record, :INV, left, right))
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, events, "discordant_pair_evidence")
end

function _cluster_key(event::_SVEvidence)
    return (event.sv_type, event.chrom1, event.chrom2)
end

function _sv_genotype(support::Int, min_support::Int)
    if support >= max(8, 2 * min_support)
        return "1/1"
    elseif support >= min_support
        return "0/1"
    end
    return "./."
end

function cluster_structural_variants(events::AbstractVector{<:_SVEvidence}; window::Integer=100, min_support::Integer=2)
    grouped = Dict{Tuple{Symbol,String,String},Vector{_SVEvidence}}()
    for event in events
        push!(get!(grouped, _cluster_key(event), _SVEvidence[]), event)
    end

    calls = StructuralVariant[]
    for ((sv_type, chrom1, chrom2), bucket) in grouped
        sort!(bucket; by=event -> event.pos1)
        tree = IntervalTree{Int}()
        clusters = Dict{Int,_SVCluster}()
        next_cluster = 1

        for event in bucket
            candidate_clusters = unique(query_overlaps(tree, event.pos1 - window, event.pos1 + window))
            selected = 0
            for cluster_id in candidate_clusters
                cluster = clusters[cluster_id]
                center2 = Int(round(mean(cluster.pos2_values)))
                if abs(event.pos2 - center2) <= window
                    selected = cluster_id
                    break
                end
            end

            if selected == 0
                selected = next_cluster
                next_cluster += 1
                clusters[selected] = _SVCluster(chrom1, chrom2, sv_type, Int[event.pos1], Int[event.pos2], String[event.read_name])
            else
                cluster = clusters[selected]
                push!(cluster.pos1_values, event.pos1)
                push!(cluster.pos2_values, event.pos2)
                push!(cluster.reads, event.read_name)
            end

            Base.insert!(tree, event.pos1 - window, event.pos1 + window, selected)
        end

        cluster_values = collect(values(clusters))
        for cluster_idx in eachindex(cluster_values)
            cluster = cluster_values[cluster_idx]
            support_reads = unique(cluster.reads)
            support = length(support_reads)
            support < min_support && continue
            pos1 = Int(round(mean(cluster.pos1_values)))
            pos2 = Int(round(mean(cluster.pos2_values)))
            genotype = _sv_genotype(support, min_support)
            push!(calls, StructuralVariant(cluster.chrom1, pos1, cluster.chrom2, pos2, cluster.sv_type, support_reads, genotype))
        end
    end

    sort!(calls; by=call -> (call.chrom1, call.pos1, call.sv_type, call.chrom2, call.pos2))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, calls, "cluster_structural_variants")
end

function call_structural_variants(records::AbstractVector{<:BamRecord}; min_sv_size::Integer=30, min_softclip::Integer=20, insert_size_mean::Real=500.0, insert_size_std::Real=100.0, z_threshold::Real=3.0, cluster_window::Integer=100, min_support::Integer=2)
    events = _SVEvidence[]
    for record in records
        append!(events, split_read_evidence(record; min_sv_size=min_sv_size, min_softclip=min_softclip))
        append!(events, discordant_pair_evidence(record; insert_size_mean=insert_size_mean, insert_size_std=insert_size_std, z_threshold=z_threshold, min_sv_size=min_sv_size))
    end
    result = cluster_structural_variants(events; window=cluster_window, min_support=min_support)
    _ctx = active_provenance_context()
    return _register_longread_result!(_ctx, result, "call_structural_variants"; parameters=(n_records=length(records), min_sv_size=Int(min_sv_size), min_support=Int(min_support), n_svs=length(result)))
end

function call_structural_variants(file::BamFile, kwargs...)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, call_structural_variants(file.records; kwargs...), "call_structural_variants")
end

function _longread_kmeans_lloyd(x::AbstractMatrix{<:Real}, k::Int; max_iter::Int=50, seed::Int=1)
    n = size(x, 1)
    n >= k || throw(ArgumentError("k must be <= number of rows"))
    rng = MersenneTwister(seed)
    centers = Matrix{Float64}(x[rand(rng, 1:n, k), :])
    labels = ones(Int, n)

    for _ in 1:max_iter
        changed = false
        for i in 1:n
            xi = @view x[i, :]
            best = 1
            bestd = Inf
            for c in 1:k
                d = sum(abs2, xi .- @view(centers[c, :]))
                if d < bestd
                    bestd = d
                    best = c
                end
            end
            if labels[i] != best
                labels[i] = best
                changed = true
            end
        end

        for c in 1:k
            idx = findall(==(c), labels)
            if isempty(idx)
                centers[c, :] .= x[rand(rng, 1:n), :]
            else
                centers[c, :] .= vec(mean(x[idx, :], dims=1))
            end
        end
        !changed && break
    end
    return labels, centers
end

function _consensus_from_cluster(reads::Vector{String})
    maxlen = maximum(ncodeunits.(reads))
    consensus = Vector{Char}(undef, maxlen)
    for pos in 1:maxlen
        counts = Dict{Char,Int}()
        for r in reads
            if pos <= ncodeunits(r)
                c = r[pos]
                counts[c] = get(counts, c, 0) + 1
            end
        end
        consensus[pos] = isempty(counts) ? 'N' : findmax(counts)[2]
    end
    return String(consensus)
end

"""
    isoseq_consensus(reads; n_clusters=3)

IsoSeq-style consensus generation using clustering in length/GC space.
"""
function isoseq_consensus(reads::AbstractVector{<:BioSequence{DNAAlphabet}}; n_clusters::Int=3, seed::Int=1)
    seqs = String.(reads)
    lengths = Float64.(ncodeunits.(seqs))
    gc = [count(c -> c in ('G', 'C'), collect(uppercase(s))) / max(ncodeunits(s), 1) for s in seqs]
    features = hcat(lengths, gc)

    k = clamp(n_clusters, 1, length(seqs))
    labels, _ = _longread_kmeans_lloyd(features, k; seed=seed)

    out = DataFrame(cluster=Int[], n_reads=Int[], consensus=String[])
    for c in 1:k
        idx = findall(==(c), labels)
        isempty(idx) && continue
        cons = _consensus_from_cluster(seqs[idx])
        push!(out, (c, length(idx), cons))
    end
    _ctx = active_provenance_context()
    return _register_longread_result!(_ctx, out, "isoseq_consensus"; parameters=(n_reads=length(reads), n_clusters=k, n_consensus=nrow(out)))
end

function isoseq_consensus(reads::Vector{String}; n_clusters::Int=3, seed::Int=1)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, isoseq_consensus(DNASeq.(reads); n_clusters=n_clusters, seed=seed), "isoseq_consensus")
end

_tool_cmd(program::String, args::Vector{String}) = `$(program) $(args...)`

function minimap2_align(reference_fasta::AbstractString, reads_fastq::AbstractString; output_sam::AbstractString="alignment.sam", preset::AbstractString="map-ont", threads::Int=4, dry_run::Bool=true)
    cmd = _tool_cmd("minimap2", ["-ax", String(preset), "-t", string(threads), String(reference_fasta), String(reads_fastq)])
    if dry_run
        return (cmd=string(cmd), output=String(output_sam), status=:planned)
    end
    open(String(output_sam), "w") do io
        run(pipeline(cmd, stdout=io))
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (cmd=string(cmd), output=String(output_sam), status=:ok), "minimap2_align")
end

function graphaligner_align(graph_gfa::AbstractString, reads_fastq::AbstractString; output_gaf::AbstractString="graph_alignment.gaf", threads::Int=4, dry_run::Bool=true)
    cmd = _tool_cmd("GraphAligner", ["-g", String(graph_gfa), "-f", String(reads_fastq), "-t", string(threads), "-a", String(output_gaf)])
    if dry_run
        return (cmd=string(cmd), output=String(output_gaf), status=:planned)
    end
    run(cmd)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (cmd=string(cmd), output=String(output_gaf), status=:ok), "graphaligner_align")
end

function sniffles_call(bam_path::AbstractString; output_vcf::AbstractString="sniffles.vcf", min_support::Int=5, dry_run::Bool=true)
    cmd = _tool_cmd("sniffles", ["-i", String(bam_path), "-v", String(output_vcf), "--minsupport", string(min_support)])
    if dry_run
        return (cmd=string(cmd), output=String(output_vcf), status=:planned)
    end
    run(cmd)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (cmd=string(cmd), output=String(output_vcf), status=:ok), "sniffles_call")
end

function svim_call(bam_path::AbstractString, reference_fasta::AbstractString; output_dir::AbstractString="svim_out", min_sv_size::Int=50, dry_run::Bool=true)
    cmd = _tool_cmd("svim", ["alignment", String(output_dir), String(bam_path), String(reference_fasta), "--min_sv_size", string(min_sv_size)])
    if dry_run
        return (cmd=string(cmd), output=String(output_dir), status=:planned)
    end
    run(cmd)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (cmd=string(cmd), output=String(output_dir), status=:ok), "svim_call")
end

function whatshap_phase(vcf_path::AbstractString, bam_path::AbstractString; output_vcf::AbstractString="phased.vcf", dry_run::Bool=true)
    cmd = _tool_cmd("whatshap", ["phase", "-o", String(output_vcf), String(vcf_path), String(bam_path)])
    if dry_run
        return (cmd=string(cmd), output=String(output_vcf), status=:planned)
    end
    run(cmd)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (cmd=string(cmd), output=String(output_vcf), status=:ok), "whatshap_phase")
end

"""
    phase_reads_by_alleles(read_by_variant)

Simple EM-like read-backed phasing from a read x variant allele matrix.
"""
function phase_reads_by_alleles(read_by_variant::AbstractMatrix{<:Integer}; max_iter::Int=25, seed::Int=1)
    X = Int.(read_by_variant)
    n_reads, n_vars = size(X)
    rng = MersenneTwister(seed)
    z = rand(rng, Bool, n_reads)

    hap1 = zeros(Float64, n_vars)
    hap2 = zeros(Float64, n_vars)

    for _ in 1:max_iter
        for v in 1:n_vars
            a = Int[]
            b = Int[]
            for r in 1:n_reads
                val = X[r, v]
                val == -1 && continue
                if z[r]
                    push!(a, val)
                else
                    push!(b, val)
                end
            end
            hap1[v] = isempty(a) ? 0.5 : mean(a)
            hap2[v] = isempty(b) ? 0.5 : mean(b)
        end

        changed = false
        for r in 1:n_reads
            row = @view X[r, :]
            idx = findall(!=(-1), row)
            isempty(idx) && continue
            d1 = sum(abs2, Float64.(row[idx]) .- hap1[idx])
            d2 = sum(abs2, Float64.(row[idx]) .- hap2[idx])
            newz = d1 <= d2
            if newz != z[r]
                z[r] = newz
                changed = true
            end
        end
        !changed && break
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (read_haplotype=Int.(z) .+ 1, haplotype1=hap1, haplotype2=hap2), "phase_reads_by_alleles")
end

function _suffix_prefix_overlap(a::BioSequence, b::BioSequence, min_overlap::Int)
    max_k = min(ncodeunits(a), ncodeunits(b))
    best = 0
    for k in max(min_overlap, 1):max_k
        if a[(end - k + 1):end] == b[1:k]
            best = k
        end
    end
    return best
end

function _suffix_prefix_overlap(a::String, b::String, min_overlap::Int)
    max_k = min(ncodeunits(a), ncodeunits(b))
    best = 0
    for k in max(min_overlap, 1):max_k
        if a[(end - k + 1):end] == b[1:k]
            best = k
        end
    end
    return best
end

"""
    overlap_layout_consensus(reads; min_overlap=30)

Greedy overlap-layout-consensus assembly from long reads.
"""
function overlap_layout_consensus(reads::AbstractVector; min_overlap::Int=30)
    seqs = [uppercase(_record_sequence(r)) for r in reads]
    if isempty(seqs)
        _register_longread_result!(_ctx, nothing, "overlap_layout_consensus"; parameters=(n_reads=0, min_overlap=min_overlap))
        return (consensus="", used_reads=Int[], unused_reads=Int[])
    end

    lengths = ncodeunits.(seqs)
    start = argmax(lengths)
    used = Set([start])
    consensus = seqs[start]

    while true
        best_j = 0
        best_ov = 0
        for j in eachindex(seqs)
            j in used && continue
            ov = _suffix_prefix_overlap(consensus, seqs[j], min_overlap)
            if ov > best_ov
                best_ov = ov
                best_j = j
            end
        end

        if best_j == 0
            break
        end

        push!(used, best_j)
        consensus *= seqs[best_j][(best_ov + 1):end]
    end

    used_idx   = sort!(collect(used))
    unused_idx = sort!(setdiff(collect(eachindex(seqs)), used_idx))
    result = (consensus=consensus, used_reads=used_idx, unused_reads=unused_idx)
    _ctx = active_provenance_context()
    return _register_longread_result!(_ctx, result, "overlap_layout_consensus"; parameters=(n_reads=length(reads), min_overlap=min_overlap, consensus_length=ncodeunits(consensus), n_used=length(used_idx)))
end

"""
    read_n50(reads)

Compute N50 for read lengths.
"""
function read_n50(reads::AbstractVector)
    lengths = sort(ncodeunits.(uppercase.(_record_sequence.(reads))); rev=true)
    if isempty(lengths)
        _register_longread_result!(_ctx, 0, "read_n50"; parameters=(n_reads=0, n50=0))
        return 0
    end
    target = sum(lengths) / 2
    running = 0.0
    n50 = 0
    for len in lengths
        running += len
        if running >= target
            n50 = len
            break
        end
    end
    n50 = n50 == 0 ? last(lengths) : n50
    _ctx = active_provenance_context()
    return _register_longread_result!(_ctx, n50, "read_n50"; parameters=(n_reads=length(reads), n50=n50, total_bases=sum(lengths)))
end

"""
    overlap_graph_table(reads; min_overlap=30)

Construct a directed overlap graph table between reads.
"""
function overlap_graph_table(reads::AbstractVector; min_overlap::Int=30, threaded::Bool=true)
    seqs = [uppercase(_record_sequence(r)) for r in reads]
    n = length(seqs)
    rows = Vector{Vector{NTuple{4,Any}}}(undef, n)

    threaded_foreach(n, i -> begin
        local_edges = NTuple{4,Any}[]
        for j in 1:n
            i == j && continue
            ov = _suffix_prefix_overlap(seqs[i], seqs[j], min_overlap)
            ov >= min_overlap || continue
            push!(local_edges, ("read_$(i)", "read_$(j)", ov, ov / max(min(ncodeunits(seqs[i]), ncodeunits(seqs[j])), 1)))
        end
        rows[i] = local_edges
    end; threaded=threaded)

    src = String[]
    dst = String[]
    overlap = Int[]
    frac = Float64[]
    for edges in rows
        for e in edges
            push!(src, e[1])
            push!(dst, e[2])
            push!(overlap, Int(e[3]))
            push!(frac, Float64(e[4]))
        end
    end
    result = DataFrame(source=src, target=dst, overlap_bp=overlap, overlap_fraction=frac)
    _ctx = active_provenance_context()
    return _register_longread_result!(_ctx, result, "overlap_graph_table"; parameters=(n_reads=n, min_overlap=min_overlap, n_edges=nrow(result)))
end

"""
    porec_contact_table(segments)

Build long-range contact pairs from per-read multi-segment alignments
(Pore-C/Omni-C style).
"""
function porec_contact_table(segments::DataFrame; read_col::Symbol=:read_id, chrom_col::Symbol=:chrom, start_col::Symbol=:start, stop_col::Symbol=:stop)
    hasproperty(segments, read_col) || throw(ArgumentError("missing read column"))
    hasproperty(segments, chrom_col) || throw(ArgumentError("missing chrom column"))
    hasproperty(segments, start_col) || throw(ArgumentError("missing start column"))
    hasproperty(segments, stop_col) || throw(ArgumentError("missing stop column"))

    out = DataFrame(read_id=String[], chrom1=String[], pos1=Int[], chrom2=String[], pos2=Int[], distance=Int[])
    for g in groupby(segments, read_col)
        n = nrow(g)
        n >= 2 || continue
        mids = [Int(round((Int(g[i, start_col]) + Int(g[i, stop_col])) / 2)) for i in 1:n]
        for i in 1:(n - 1)
            for j in (i + 1):n
                c1 = String(g[i, chrom_col])
                c2 = String(g[j, chrom_col])
                p1 = mids[i]
                p2 = mids[j]
                d = c1 == c2 ? abs(p2 - p1) : -1
                push!(out, (String(g[1, read_col]), c1, p1, c2, p2, d))
            end
        end
    end
    _ctx = active_provenance_context()
    return _register_longread_result!(_ctx, out, "porec_contact_table"; parents=provenance_parent_ids(segments), parameters=(n_segments=nrow(segments), n_contacts=nrow(out)))
end

"""
    omnic_contact_matrix(contacts; chrom="chr1", bin_size=100_000)

Generate an intra-chromosomal contact matrix from long-range contacts.
"""
function omnic_contact_matrix(contacts::DataFrame; chrom::AbstractString="chr1", bin_size::Int=100_000)
    hasproperty(contacts, :chrom1) || throw(ArgumentError("contacts must contain chrom1"))
    hasproperty(contacts, :pos1) || throw(ArgumentError("contacts must contain pos1"))
    hasproperty(contacts, :chrom2) || throw(ArgumentError("contacts must contain chrom2"))
    hasproperty(contacts, :pos2) || throw(ArgumentError("contacts must contain pos2"))

    sub = filter(row -> String(row.chrom1) == String(chrom) && String(row.chrom2) == String(chrom), contacts)
    if nrow(sub) == 0
        result = (matrix=zeros(Int, 0, 0), bins=Int[], chrom=String(chrom), bin_size=bin_size)
        return _register_longread_result!(_ctx, result, "omnic_contact_matrix"; parameters=(n_contacts=0, chrom=String(chrom), bin_size=bin_size))
    end

    b1 = Int.(fld.(sub.pos1, bin_size)) .+ 1
    b2 = Int.(fld.(sub.pos2, bin_size)) .+ 1
    nb = max(maximum(b1), maximum(b2))
    M = zeros(Int, nb, nb)

    for i in eachindex(b1)
        x = b1[i]
        y = b2[i]
        M[x, y] += 1
        x == y || (M[y, x] += 1)
    end

    result = (matrix=M, bins=collect(1:nb), chrom=String(chrom), bin_size=Int(bin_size))
    _ctx = active_provenance_context()
    return _register_longread_result!(_ctx, result, "omnic_contact_matrix"; parents=provenance_parent_ids(contacts), parameters=(n_contacts=nrow(sub), chrom=String(chrom), bin_size=bin_size, n_bins=nb))
end

end
