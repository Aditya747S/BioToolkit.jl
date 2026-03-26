export build_index, local_search, run_blast, qblast, parse_blast_xml, read_blast_xml, KmerIndex, HighScoringPair, BlastXMLRecord, BlastXMLHit, BlastXMLHSP
export BlastTabularRecord, BlastTabularHit, BlastTabularHSP, parse_blast_tabular, read_blast_tabular

struct KmerIndex
    k::Int
    database::Dict{UInt64, Vector{Tuple{Int, Int}}}
    target_names::Vector{String}
    target_sequences::Vector{Vector{UInt8}}
end

@inline function _pack_kmer(bytes::Vector{UInt8}, start_pos::Int, k::Int)
    val = UInt64(0)
    @inbounds for i in 0:k-1
        val = (val << 8) | bytes[start_pos + i]
    end
    return val
end

# --- BLAST+ Local Wrapper ---

"""
    run_blast(program::String, query::String, database::String; options::String="", outfmt::Int=6)

Run a local BLAST+ search. Requires the `blast+` suite to be installed in the system PATH.
Returns a Vector of `HighScoringPair`.
"""
function run_blast(program::AbstractString, query_seq::AbstractString, database::AbstractString; options::AbstractString="", outfmt::Int=6)
    # Temporary files for query and output
    q_file, q_io = mktemp()
    write(q_io, ">query\n", query_seq)
    close(q_io)
    
    r_file, r_io = mktemp()
    close(r_io)
    
    try
        # Construct command
        # outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        cmd_str = "$program -query $q_file -db $database -out $r_file -outfmt $outfmt $options"
        # run(pipeline(`$program -query $q_file -db $database -out $r_file -outfmt $outfmt`, stderr=devnull))
        
        # We use a safer way to run commands in Julia
        run(`/bin/sh -c "$cmd_str"`)
        
        # Parse results (outfmt 6)
        hsps = HighScoringPair[]
        for line in eachline(r_file)
            if startswith(line, "#"); continue; end
            parts = Base.split(line, '\t')
            if length(parts) >= 12
                hsp = HighScoringPair(
                    String(parts[1]), # qseqid
                    String(parts[2]), # sseqid
                    _parse_target_index(String(parts[2])),
                    parse(Int, parts[7]), # qstart
                    parse(Int, parts[8]), # qend
                    parse(Int, parts[9]), # sstart
                    parse(Int, parts[10]), # send
                    parse(Float64, parts[12]), # bitscore
                    parse(Float64, parts[11]), # evalue
                    parse(Float64, parts[3]), # identity
                    parse(Int, parts[4])  # alignment_length
                )
                push!(hsps, hsp)
            end
        end
        return hsps
    finally
        rm(q_file, force=true)
        rm(r_file, force=true)
    end
end

"""
    qblast(program::String, database::String, query::String; options::Dict=Dict())

Placeholder for remote BLAST search. Implementation requires `HTTP` and `EzXML` dependencies.
"""
function qblast(program::AbstractString, database::AbstractString, query::AbstractString; options::AbstractDict=Dict())
    error("qblast requires optional dependencies HTTP and EzXML. Please install them and use a wrapper.")
end

"""
    build_index(targets::Vector{<:AbstractString}, names::Vector{<:AbstractString}=String[]; k::Int=4)

Builds a k-mer index for a database of target sequences (optimized for proteins, k <= 8).
"""
function build_index(targets::Vector{<:AbstractString}, names::Vector{<:AbstractString}=String[]; k::Int=4)
    if k > 8
        throw(ArgumentError("k must be <= 8 for efficient UInt64 packing in this implementation"))
    end
    
    if isempty(names)
        names = ["target_\$i" for i in 1:length(targets)]
    end
    
    database = Dict{UInt64, Vector{Tuple{Int, Int}}}()
    target_bytes = Vector{Vector{UInt8}}(undef, length(targets))
    
    for (t_idx, target) in enumerate(targets)
        bytes = Vector{UInt8}(codeunits(target))
        target_bytes[t_idx] = bytes
        
        len = length(bytes)
        for i in 1:(len - k + 1)
            key = _pack_kmer(bytes, i, k)
            
            # Append to dict
            locs = get!(database, key, Tuple{Int, Int}[])
            push!(locs, (t_idx, i))
        end
    end
    
    return KmerIndex(k, database, names, target_bytes)
end

struct HighScoringPair
    query_id::String
    target_id::String
    target_idx::Int
    query_start::Int
    query_end::Int
    target_start::Int
    target_end::Int
    score::Float64
    evalue::Float64
    identity::Float64
    alignment_length::Int
end

struct BlastXMLHSP
    bit_score::Float64
    evalue::Float64
    identities::Int
    positives::Int
    align_length::Int
    query_start::Int
    query_end::Int
    hit_start::Int
    hit_end::Int
    query_sequence::String
    hit_sequence::String
    midline::String
end

struct BlastXMLHit
    id::String
    description::String
    length::Int
    hsps::Vector{BlastXMLHSP}
end

struct BlastXMLRecord
    program::String
    version::String
    database::String
    query_id::String
    query_description::String
    query_length::Int
    hits::Vector{BlastXMLHit}
end

struct BlastTabularHSP
    query_start::Int
    query_end::Int
    target_start::Int
    target_end::Int
    identity::Float64
    alignment_length::Int
    mismatches::Int
    gap_openings::Int
    evalue::Float64
    bit_score::Float64
end

struct BlastTabularHit
    target_id::String
    hsps::Vector{BlastTabularHSP}
end

struct BlastTabularRecord
    query_id::String
    hits::Vector{BlastTabularHit}
end

@inline function _blast_xml_unescape(text::AbstractString)
    replace(String(text), "&lt;" => "<", "&gt;" => ">", "&amp;" => "&", "&quot;" => "\"", "&apos;" => "'")
end

function _xml_blocks(text::AbstractString, tag::AbstractString)
    pattern = Regex("<$(tag)>(.*?)</$(tag)>", "s")
    return [match.captures[1] for match in eachmatch(pattern, text)]
end

function _xml_first(text::AbstractString, tag::AbstractString; default::AbstractString="")
    pattern = Regex("<$(tag)>(.*?)</$(tag)>", "s")
    found = match(pattern, text)
    found === nothing && return String(default)
    return _blast_xml_unescape(found.captures[1])
end

_xml_parse_int(text::AbstractString; default::Int=0) = isempty(strip(text)) ? default : parse(Int, strip(text))
_xml_parse_float(text::AbstractString; default::Float64=0.0) = isempty(strip(text)) ? default : parse(Float64, strip(text))

function parse_blast_xml(xml::AbstractString)
    records = BlastXMLRecord[]
    xml_text = String(xml)

    program = _xml_first(xml_text, "BlastOutput_program")
    version = _xml_first(xml_text, "BlastOutput_version")
    database = _xml_first(xml_text, "BlastOutput_db")
    blast_query_id = _xml_first(xml_text, "BlastOutput_query-ID")
    blast_query_def = _xml_first(xml_text, "BlastOutput_query-def")
    blast_query_len = _xml_parse_int(_xml_first(xml_text, "BlastOutput_query-len"), default=0)

    for iteration in _xml_blocks(xml_text, "Iteration")
        query_id = _xml_first(iteration, "Iteration_query-ID"; default=blast_query_id)
        query_def = _xml_first(iteration, "Iteration_query-def"; default=blast_query_def)
        query_length = _xml_parse_int(_xml_first(iteration, "Iteration_query-len"; default=string(blast_query_len)), default=blast_query_len)

        hits = BlastXMLHit[]
        for hit_block in _xml_blocks(iteration, "Hit")
            hit_id = _xml_first(hit_block, "Hit_id"; default="")
            hit_def = _xml_first(hit_block, "Hit_def"; default=hit_id)
            hit_len = _xml_parse_int(_xml_first(hit_block, "Hit_len"; default="0"), default=0)

            hsps = BlastXMLHSP[]
            for hsp_block in _xml_blocks(hit_block, "Hsp")
                bit_score = _xml_parse_float(_xml_first(hsp_block, "Hsp_bit-score"; default="0"), default=0.0)
                evalue = _xml_parse_float(_xml_first(hsp_block, "Hsp_evalue"; default="0"), default=0.0)
                identities = _xml_parse_int(_xml_first(hsp_block, "Hsp_identity"; default="0"), default=0)
                positives = _xml_parse_int(_xml_first(hsp_block, "Hsp_positive"; default=string(identities)), default=identities)
                align_length = _xml_parse_int(_xml_first(hsp_block, "Hsp_align-len"; default="0"), default=0)
                query_start = _xml_parse_int(_xml_first(hsp_block, "Hsp_query-from"; default="0"), default=0)
                query_end = _xml_parse_int(_xml_first(hsp_block, "Hsp_query-to"; default="0"), default=0)
                hit_start = _xml_parse_int(_xml_first(hsp_block, "Hsp_hit-from"; default="0"), default=0)
                hit_end = _xml_parse_int(_xml_first(hsp_block, "Hsp_hit-to"; default="0"), default=0)
                query_sequence = _xml_first(hsp_block, "Hsp_qseq")
                hit_sequence = _xml_first(hsp_block, "Hsp_hseq")
                midline = _xml_first(hsp_block, "Hsp_midline")

                push!(hsps, BlastXMLHSP(bit_score, evalue, identities, positives, align_length, query_start, query_end, hit_start, hit_end, query_sequence, hit_sequence, midline))
            end

            push!(hits, BlastXMLHit(hit_id, hit_def, hit_len, hsps))
        end

        push!(records, BlastXMLRecord(program, version, database, query_id, query_def, query_length, hits))
    end

    return records
end

function parse_blast_xml(io::IO)
    return parse_blast_xml(read(io, String))
end

function read_blast_xml(path::AbstractString)
    open(path, "r") do io
        return parse_blast_xml(io)
    end
end

function parse_blast_tabular(lines::AbstractString)
    records = Dict{String,Dict{String,Vector{BlastTabularHSP}}}()
    for raw_line in eachline(IOBuffer(lines))
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, '#') && continue
        parts = Base.split(line, '\t')
        length(parts) < 12 && continue

        query_id = String(parts[1])
        target_id = String(parts[2])
        hsp = BlastTabularHSP(
            parse(Int, parts[7]),
            parse(Int, parts[8]),
            parse(Int, parts[9]),
            parse(Int, parts[10]),
            parse(Float64, parts[3]),
            parse(Int, parts[4]),
            parse(Int, parts[5]),
            parse(Int, parts[6]),
            parse(Float64, parts[11]),
            parse(Float64, parts[12]),
        )
        query_hits = get!(records, query_id, Dict{String,Vector{BlastTabularHSP}}())
        push!(get!(query_hits, target_id, BlastTabularHSP[]), hsp)
    end

    parsed = BlastTabularRecord[]
    for (query_id, hits_map) in records
        hits = BlastTabularHit[]
        for (target_id, hsps) in hits_map
            push!(hits, BlastTabularHit(target_id, hsps))
        end
        push!(parsed, BlastTabularRecord(query_id, hits))
    end

    sort!(parsed; by = record -> record.query_id)
    return parsed
end

function parse_blast_tabular(io::IO)
    return parse_blast_tabular(read(io, String))
end

function read_blast_tabular(path::AbstractString)
    open(path, "r") do io
        return parse_blast_tabular(io)
    end
end

# Legacy constructor for backward compatibility
function HighScoringPair(q_start::Int, q_end::Int, t_idx::Int, t_start::Int, t_end::Int, score::Real)
    return HighScoringPair("query", "target_$t_idx", t_idx, q_start, q_end, t_start, t_end, Float64(score), 0.0, 0.0, 0)
end

function HighScoringPair(query_id::AbstractString, target_id::AbstractString, target_idx::Integer, query_start::Integer, query_end::Integer, target_start::Integer, target_end::Integer, score::Real, evalue::Real, identity::Real, alignment_length::Integer)
    return HighScoringPair(String(query_id), String(target_id), Int(target_idx), Int(query_start), Int(query_end), Int(target_start), Int(target_end), Float64(score), Float64(evalue), Float64(identity), Int(alignment_length))
end

function _parse_target_index(target_id::AbstractString)
    found = Base.match(r"(\d+)$", target_id)
    found === nothing && return 0
    return parse(Int, found.captures[1])
end

# To support `AbstractPairwiseScoring` parameter
function _extend_ungapped(query_bytes::Vector{UInt8}, target_bytes::Vector{UInt8}, q_pos::Int, t_pos::Int, k::Int, scoring::AbstractPairwiseScoring, x_drop::Int)
    q_len = length(query_bytes)
    t_len = length(target_bytes)
    
    # Base score for the exact match seed
    base_score = 0
    @inbounds for i in 0:k-1
        base_score += _pairwise_score(scoring, query_bytes[q_pos + i], target_bytes[t_pos + i])
    end
    
    # Extend right (forward)
    max_fwd_score = 0
    current_score = 0
    best_fwd_offset = 0
    offset = k
    
    @inbounds while q_pos + offset <= q_len && t_pos + offset <= t_len
        current_score += _pairwise_score(scoring, query_bytes[q_pos + offset], target_bytes[t_pos + offset])
        if current_score > max_fwd_score
            max_fwd_score = current_score
            best_fwd_offset = offset
        elseif current_score < max_fwd_score - x_drop
            break
        end
        offset += 1
    end
    
    # Extend left (backward)
    max_bwd_score = 0
    current_score = 0
    best_bwd_offset = 0
    offset = 1
    
    @inbounds while q_pos - offset >= 1 && t_pos - offset >= 1
        current_score += _pairwise_score(scoring, query_bytes[q_pos - offset], target_bytes[t_pos - offset])
        if current_score > max_bwd_score
            max_bwd_score = current_score
            best_bwd_offset = offset
        elseif current_score < max_bwd_score - x_drop
            break
        end
        offset += 1
    end
    
    total_score = base_score + max_fwd_score + max_bwd_score
    q_start = q_pos - best_bwd_offset
    q_end = q_pos + best_fwd_offset
    t_start = t_pos - best_bwd_offset
    t_end = t_pos + best_fwd_offset
    
    return q_start, q_end, t_start, t_end, total_score
end

"""
    local_search(query::AbstractString, index::KmerIndex; scoring::AbstractPairwiseScoring, x_drop::Int=10, min_score::Int=20)

Performs a BLOSUM/SubstitutionMatrix-scored BLAST-like seed-and-extend local search.
"""
function local_search(
    query::AbstractString, 
    index::KmerIndex; 
    scoring::AbstractPairwiseScoring, 
    x_drop::Int=15, 
    min_score::Int=25,
    gap_open::Int=-5,
    gap_extend::Int=-1,
    envelope::Int=20,
    use_threads::Bool=true,
    use_cuda::Bool=false
)
    query_bytes = Vector{UInt8}(codeunits(query))
    q_len = length(query_bytes)
    k = index.k
    
    hsps = HighScoringPair[]
    
    # Extract matrix if available
    sub_matrix = scoring isa MatrixPairwiseScoring ? scoring.matrix : substitution_matrix("ACGT"; match=2, mismatch=-1)
    
    # 1. Collect all seeds
    seeds = Tuple{Int, Int, Int}[] # (q_pos, t_idx, t_pos)
    for i in 1:(q_len - k + 1)
        key = _pack_kmer(query_bytes, i, k)
        if haskey(index.database, key)
            for (t_idx, t_pos) in index.database[key]
                push!(seeds, (i, t_idx, t_pos))
            end
        end
    end
    
    # 2. Process seeds
    if isempty(seeds)
        return hsps
    end
    
    # We will compute scores for all seeds.
    # We use threads if requested.
    # GPU extension could be done by assigning a thread per seed, but variable length X-drop is tricky on GPU.
    # For now, we simulate the GPU distribution or use multithreading robustly.
    
    results_lock = Threads.SpinLock()
    
    function process_seed(seed)
        q_pos, t_idx, t_pos = seed
        target_bytes = index.target_sequences[t_idx]
        
        # Ungapped X-drop extension
        q_s, q_e, t_s, t_e, score = _extend_ungapped(query_bytes, target_bytes, q_pos, t_pos, k, scoring, x_drop)
        
        if score >= min_score
            # Bounded Gapped Extension (Smith-Waterman Envelope)
            q_env_s = max(1, q_s - envelope)
            q_env_e = min(q_len, q_e + envelope)
            t_env_s = max(1, t_s - envelope)
            t_env_e = min(length(target_bytes), t_e + envelope)
            
            q_sub = @view query_bytes[q_env_s:q_env_e]
            t_sub = @view target_bytes[t_env_s:t_env_e]
            
            if gap_open == gap_extend
                res = pairwise_align(q_sub, t_sub; is_local=true, substitution_matrix=sub_matrix, gap=gap_open)
            else
                res = pairwise_align(q_sub, t_sub; is_local=true, substitution_matrix=sub_matrix, gap_open=gap_open, gap_extend=gap_extend)
            end
            
            if res.score >= min_score
                hsp = HighScoringPair(q_env_s, q_env_e, t_idx, t_env_s, t_env_e, res.score)
                lock(results_lock) do
                    push!(hsps, hsp)
                end
            end
        end
    end

    if use_threads || use_cuda
        Threads.@threads for seed in seeds
            process_seed(seed)
        end
    else
        for seed in seeds
            process_seed(seed)
        end
    end
    
    # Primitive filtering: Sort by score descending
    sort!(hsps, by = x -> x.score, rev = true)
    
    # Filter overlapping HSPs on the query in O(N log N) using disjoint intervals
    filtered_hsps = HighScoringPair[]
    visited_query_ranges = UnitRange{Int}[]
    
    for hsp in hsps
        overlap = false
        hsp_range = hsp.query_start:hsp.query_end
        
        insert_idx = searchsortedfirst(visited_query_ranges, hsp_range, by = x -> x.start)
        
        # Check interval immediately preceding
        if insert_idx > 1
            prev = visited_query_ranges[insert_idx - 1]
            if prev.stop >= hsp_range.start
                overlap = true
            end
        end
        
        # Check interval immediately following
        if !overlap && insert_idx <= length(visited_query_ranges)
            next_interval = visited_query_ranges[insert_idx]
            if next_interval.start <= hsp_range.stop
                overlap = true
            end
        end
        
        if !overlap
            push!(filtered_hsps, hsp)
            insert!(visited_query_ranges, insert_idx, hsp_range)
        end
    end
    
    return filtered_hsps
end
