# ==============================================================================
# io.jl — Biological Format I/O
#
# Provides high-performance, scientific-grade parsers and writers for common
# bioinformatics formats (FASTA, FASTQ, VCF, BED, GFF3, GenBank).
#
# Design decisions:
#   - Specialized BioAlphabet dispatch for sequence formats.
#   - Automatic alphabet inference for untyped reads.
#   - Byte-level scanning for performance.
#   - Samtools-compatible .fai indexing for random access.
# ==============================================================================

# ─── FASTA I/O ───────────────────────────

const _FASTA_IM = UInt32(139968)
const _FASTA_IA = UInt32(3877)
const _FASTA_IC = UInt32(29573)
const _FASTA_LINE_LENGTH = 60

# Direct 139,968-byte lookup table for O(1) L2-cached random nucleotide selection
function _make_fasta_lut(mapping::Vector{Tuple{UInt8, Float64}})
    lut = Vector{UInt8}(undef, _FASTA_IM)
    cumprob = Vector{UInt32}(undef, length(mapping))
    acc = 0.0
    for (i, (char, prob)) in enumerate(mapping)
        acc += prob
        cumprob[i] = floor(UInt32, acc * _FASTA_IM)
    end
    
    for s in 0:(_FASTA_IM-1)
        cnt = 1
        for cp in cumprob
            if cp <= s
                cnt += 1
            else
                break
            end
        end
        lut[s + 1] = mapping[cnt][1]
    end
    return lut
end

# Lock-free LCG Jump-ahead in O(log k) using binary matrix exponentiation
@inline function _lcg_jump(seed::UInt32, k::Int)
    ra, rc = UInt64(1), UInt64(0)
    ca, cc = UInt64(_FASTA_IA), UInt64(_FASTA_IC)
    m = UInt64(_FASTA_IM)
    n = k
    while n > 0
        if (n & 1) == 1
            ra = (ra * ca) % m
            rc = (rc * ca + cc) % m
        end
        cc = (cc * ca + cc) % m
        ca = (ca * ca) % m
        n >>= 1
    end
    return UInt32((ra * UInt64(seed) + rc) % m)
end

function _fasta_repeat(io::IO, seq::Vector{UInt8}, n::Int)
    len = length(seq)
    buf = Vector{UInt8}(undef, _FASTA_LINE_LENGTH + 1)
    buf[end] = UInt8('\n')
    
    pos = 1
    rem = n
    while rem > 0
        to_write = min(rem, _FASTA_LINE_LENGTH)
        for i in 1:to_write
            buf[i] = seq[pos]
            pos = pos == len ? 1 : pos + 1
        end
        if to_write < _FASTA_LINE_LENGTH
            write(io, view(buf, 1:to_write), UInt8('\n'))
        else
            write(io, buf)
        end
        rem -= to_write
    end
end

function _fasta_random_par(io::IO, initial_seed::UInt32, lut::Vector{UInt8}, n::Int; block_lines=1024)
    block_chars = block_lines * _FASTA_LINE_LENGTH
    num_blocks = cld(n, block_chars)
    
    results = Vector{Vector{UInt8}}(undef, num_blocks)
    
    Threads.@threads for b in 0:(num_blocks-1)
        chars_in_block = min(n - b * block_chars, block_chars)
        lines_in_block = cld(chars_in_block, _FASTA_LINE_LENGTH)
        
        block_seed = _lcg_jump(initial_seed, b * block_chars)
        buf = Vector{UInt8}(undef, lines_in_block * (_FASTA_LINE_LENGTH + 1))
        
        s = block_seed
        out_idx = 1
        chars_left = chars_in_block
        
        for l in 1:lines_in_block
            line_len = min(chars_left, _FASTA_LINE_LENGTH)
            for j in 1:line_len
                s = (s * _FASTA_IA + _FASTA_IC) % _FASTA_IM
                @inbounds buf[out_idx] = lut[s + 1]
                out_idx += 1
            end
            @inbounds buf[out_idx] = UInt8('\n')
            out_idx += 1
            chars_left -= line_len
        end
        results[b + 1] = buf
    end
    
    for buf in results
        write(io, buf)
    end
    
    return _lcg_jump(initial_seed, n)
end

"""
    fasta_benchmark(io::IO, n::Integer)
    fasta_benchmark(path::AbstractString, n::Integer)

High-performance FASTA generator benchmark matching the Computer Language Benchmarks Game.
Uses lock-free parallel LCG jump-ahead generation and O(1) L2-cached lookup tables.
"""
function fasta_benchmark(io::IO, n::Integer)
    n_int = Int(n)
    alu = Vector{UInt8}("GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATACAAAAATTAGCCGGGCGTGGTGGCGCGCGCCTGTAATCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCGCTTGAACCCGGGAGGCGGAGGTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAA")

    iub = [
        (UInt8('a'), 0.27), (UInt8('c'), 0.12), (UInt8('g'), 0.12), (UInt8('t'), 0.27),
        (UInt8('B'), 0.02), (UInt8('D'), 0.02), (UInt8('H'), 0.02), (UInt8('K'), 0.02),
        (UInt8('M'), 0.02), (UInt8('N'), 0.02), (UInt8('R'), 0.02), (UInt8('S'), 0.02),
        (UInt8('V'), 0.02), (UInt8('W'), 0.02), (UInt8('Y'), 0.02)
    ]
    iub_lut = _make_fasta_lut(iub)

    homosapiens = [
        (UInt8('a'), 0.3029549426680), (UInt8('c'), 0.1979883004921),
        (UInt8('g'), 0.1975473066391), (UInt8('t'), 0.3015094502008)
    ]
    hs_lut = _make_fasta_lut(homosapiens)

    write(io, ">ONE Homo sapiens alu\n")
    _fasta_repeat(io, alu, n_int * 2)

    write(io, ">TWO IUB ambiguity codes\n")
    seed = _fasta_random_par(io, UInt32(42), iub_lut, n_int * 3)

    write(io, ">THREE Homo sapiens frequency\n")
    _fasta_random_par(io, seed, hs_lut, n_int * 5)
    return io
end

function fasta_benchmark(path::AbstractString, n::Integer)
    open(path, "w") do io
        fasta_benchmark(io, n)
    end
    return path
end

function fasta_benchmark_to_bytes(n::Integer)
    io = IOBuffer()
    fasta_benchmark(io, n)
    return take!(io)
end

function _fast_parse_fasta_bytes(bytes::Vector{UInt8}, ::Type{A}; preserve_case::Bool=false) where {A <: BioAlphabet}
    len = length(bytes)
    records = SeqRecord{A}[]
    i = 1
    seen_header = false
    
    header = ""
    seq_buf = Vector{UInt8}(undef, len)
    seq_len = 0
    
    while i <= len
        b = @inbounds bytes[i]
        if b == UInt8('>')
            if seen_header
                seq_data = @inbounds seq_buf[1:seq_len]
                push!(records, SeqRecord(BioSequence{A}(seq_data; validate=false), identifier=header))
                seq_len = 0
            end
            i += 1
            h_start = i
            while i <= len && @inbounds(bytes[i]) != UInt8('\n') && @inbounds(bytes[i]) != UInt8('\r')
                i += 1
            end
            h_end = i - 1
            while h_end >= h_start && (@inbounds(bytes[h_end]) == UInt8(' ') || @inbounds(bytes[h_end]) == UInt8('\t'))
                h_end -= 1
            end
            header = String(view(bytes, h_start:h_end))
            seen_header = true
            if i <= len && @inbounds(bytes[i]) == UInt8('\r')
                i += 1
            end
            if i <= len && @inbounds(bytes[i]) == UInt8('\n')
                i += 1
            end
        elseif b == UInt8('\n') || b == UInt8('\r') || b == UInt8(' ') || b == UInt8('\t')
            i += 1
        else
            if !seen_header
                throw(ArgumentError("FASTA sequence data before header"))
            end
            u = (!preserve_case && b >= 0x61 && b <= 0x7a) ? (b - 0x20) : b
            seq_len += 1
            @inbounds seq_buf[seq_len] = u
            i += 1
        end
    end
    if seen_header
        seq_data = @inbounds seq_buf[1:seq_len]
        push!(records, SeqRecord(BioSequence{A}(seq_data; validate=false), identifier=header))
    end
    return records
end

"""
    read_fasta(path; alphabet=DNAAlphabet, preserve_case=false) -> Vector{SeqRecord{A}}
    read_fasta(io::IO; alphabet=DNAAlphabet, preserve_case=false) -> Vector{SeqRecord{A}}

Read all FASTA records using high-performance zero-allocation byte scanning.
Setting `preserve_case=true` preserves soft-masked lowercase characters.
"""
function read_fasta(path::String; alphabet::Type{<:BioAlphabet}=DNAAlphabet, preserve_case::Bool=false, prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    if _ctx === nothing
        open(path, "r") do io
            return read_fasta(io; alphabet=alphabet, preserve_case=preserve_case, _ctx=nothing)
        end
    end
    raw_bytes = read(path)
    provenance_hash = bytes2hex(sha256(raw_bytes))
    return read_fasta(IOBuffer(raw_bytes); alphabet=alphabet, preserve_case=preserve_case, _ctx=_ctx, provenance_hash=provenance_hash, provenance_source=path)
end

function read_fasta(io::IO; alphabet::Type{A}=DNAAlphabet, preserve_case::Bool=false, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx), provenance_hash::Union{Nothing,AbstractString}=nothing, provenance_source::AbstractString="read_fasta") where {A <: BioAlphabet}
    bytes = read(io)
    records = _fast_parse_fasta_bytes(bytes, A; preserve_case=preserve_case)

    if provenance_hash !== nothing
        for record in records
            record.metadata[PROVENANCE_HASH_KEY] = String(provenance_hash)
        end
    end

    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        root = register_provenance!(_ctx, "read_fasta"; parents=String[], parameters=(source=provenance_source, alphabet=string(alphabet), record_count=length(records), hash=provenance_hash))
        for (index, record) in enumerate(records)
            register_container_provenance!(_ctx, record, "read_fasta_record"; parents=[root.id], parameters=(source=provenance_source, record_index=index, identifier=record.identifier, alphabet=string(alphabet)), provenance_hash=provenance_hash)
        end
    end

    return records
end

"""
    write_fasta(path, records; width=60)
    write_fasta(io::IO, records; width=60)

Write FASTA records using chunk-buffered vectorized output streams.
"""
@inline function _materialize_records_for_provenance(records, _ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext})
    _ctx === nothing && return records
    return records isa AbstractArray ? records : collect(records)
end

function _register_path_write_provenance!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, operation::AbstractString, output::AbstractString, parents::AbstractVector{<:AbstractString}; record_count::Union{Nothing,Int}=nothing, provenance_hash::Union{Nothing,AbstractString}=nothing)
    _ctx === nothing && return nothing
    parameters = Dict{Symbol,Any}(:output => String(output))
    record_count === nothing || (parameters[:record_count] = record_count)
    provenance_hash === nothing || (parameters[:hash] = String(provenance_hash))
    register_provenance!(_ctx, operation; parents=parents, parameters=parameters)
    return nothing
end

function write_fasta(io::IO, records; width::Integer=60)
    width > 0 || throw(ArgumentError("width must be positive"))
    w = Int(width)
    buf = Vector{UInt8}(undef, 65536)
    buf_pos = 0

    @inline function flush_buf()
        if buf_pos > 0
            write(io, view(buf, 1:buf_pos))
            buf_pos = 0
        end
    end

    @inline function write_byte(b::UInt8)
        buf_pos += 1
        @inbounds buf[buf_pos] = b
        if buf_pos == length(buf)
            flush_buf()
        end
    end

    @inline function write_bytes(data::AbstractVector{UInt8})
        n = length(data)
        src_pos = 1
        while src_pos <= n
            avail = length(buf) - buf_pos
            to_copy = min(avail, n - src_pos + 1)
            copyto!(buf, buf_pos + 1, data, src_pos, to_copy)
            buf_pos += to_copy
            src_pos += to_copy
            if buf_pos == length(buf)
                flush_buf()
            end
        end
    end

    for record in records
        write_byte(UInt8('>'))
        write_bytes(codeunits(record.identifier))
        write_byte(UInt8('\n'))

        data = record.sequence.data
        len = length(data)
        i = 1
        while i <= len
            chunk_len = min(w, len - i + 1)
            write_bytes(view(data, i:(i + chunk_len - 1)))
            write_byte(UInt8('\n'))
            i += chunk_len
        end
    end
    flush_buf()
    return io
end

function write_fasta(path::String, records; width::Integer=60, prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    materialized = _materialize_records_for_provenance(records, _ctx)
    open(path, "w") do io
        write_fasta(io, materialized; width=width)
    end
    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(path)))
        _register_path_write_provenance!(_ctx, "write_fasta", path, provenance_parent_ids(materialized); record_count=length(materialized), provenance_hash=provenance_hash)
    end
    
    return path
end

# ─── FASTQ I/O ────────────────────────────────────────────────────────────────

"""
    read_fastq(path; alphabet=DNAAlphabet) -> Vector{FastqRecord{A}}
"""
function _fast_parse_fastq_bytes(bytes::Vector{UInt8}, ::Type{A}; preserve_case::Bool=false) where {A <: BioAlphabet}
    len = length(bytes)
    records = FastqRecord{A}[]
    i = 1
    
    while i <= len
        while i <= len && (@inbounds(bytes[i]) == UInt8('\n') || @inbounds(bytes[i]) == UInt8('\r') || @inbounds(bytes[i]) == UInt8(' '))
            i += 1
        end
        i > len && break
        
        @inbounds(bytes[i]) == UInt8('@') || throw(ArgumentError("Malformed FASTQ: expected '@'"))
        i += 1
        h_start = i
        while i <= len && @inbounds(bytes[i]) != UInt8('\n') && @inbounds(bytes[i]) != UInt8('\r')
            i += 1
        end
        header = String(view(bytes, h_start:(i-1)))
        
        if i <= len && @inbounds(bytes[i]) == UInt8('\r') i += 1 end
        if i <= len && @inbounds(bytes[i]) == UInt8('\n') i += 1 end
        
        seq_buf = UInt8[]
        sizehint!(seq_buf, 300)
        seq_start = i
        while i <= len
            b = @inbounds bytes[i]
            if b == UInt8('+')
                if i > seq_start
                    prev_b = @inbounds bytes[i-1]
                    if prev_b == UInt8('\n') || prev_b == UInt8('\r')
                        break
                    end
                end
            end
            if !preserve_case && b >= 0x61 && b <= 0x7a
                push!(seq_buf, b - 0x20)
            elseif b != UInt8('\n') && b != UInt8('\r') && b != UInt8(' ') && b != UInt8('\t')
                push!(seq_buf, b)
            end
            i += 1
        end
        
        i <= len && @inbounds(bytes[i]) == UInt8('+') || throw(ArgumentError("Malformed FASTQ: expected '+' line"))
        while i <= len && @inbounds(bytes[i]) != UInt8('\n') && @inbounds(bytes[i]) != UInt8('\r')
            i += 1
        end
        
        if i <= len && @inbounds(bytes[i]) == UInt8('\r') i += 1 end
        if i <= len && @inbounds(bytes[i]) == UInt8('\n') i += 1 end
        
        seq_len = length(seq_buf)
        qual_buf = UInt8[]
        sizehint!(qual_buf, seq_len)
        while i <= len && length(qual_buf) < seq_len
            b = @inbounds bytes[i]
            if b != UInt8('\n') && b != UInt8('\r')
                push!(qual_buf, b)
            end
            i += 1
        end
        
        if i <= len && @inbounds(bytes[i]) == UInt8('\r') i += 1 end
        if i <= len && @inbounds(bytes[i]) == UInt8('\n') i += 1 end
        
        identifier = _fastq_identifier(header)
        sequence = BioSequence{A}(seq_buf; validate=false)
        push!(records, FastqRecord(sequence, String(qual_buf); identifier=identifier, description=header))
    end
    return records
end

function read_fastq(path::String; alphabet::Type{A}=DNAAlphabet, preserve_case::Bool=false, prov_ctx=nothing) where {A <: BioAlphabet}
    _ctx = active_provenance_context(prov_ctx)
    if _ctx === nothing
        open(path, "r") do io
            return read_fastq(io; alphabet=alphabet, preserve_case=preserve_case, _ctx=nothing)
        end
    end
    raw_bytes = read(path)
    provenance_hash = bytes2hex(sha256(raw_bytes))
    return read_fastq(IOBuffer(raw_bytes); alphabet=alphabet, preserve_case=preserve_case, _ctx=_ctx, provenance_hash=provenance_hash, provenance_source=path)
end

function read_fastq(io::IO; alphabet::Type{A}=DNAAlphabet, preserve_case::Bool=false, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx), provenance_hash::Union{Nothing,AbstractString}=nothing, provenance_source::AbstractString="read_fastq") where {A <: BioAlphabet}
    bytes = read(io)
    records = _fast_parse_fastq_bytes(bytes, A; preserve_case=preserve_case)

    if provenance_hash !== nothing
        for record in records
            record.metadata[PROVENANCE_HASH_KEY] = String(provenance_hash)
        end
    end

    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        root = register_provenance!(_ctx, "read_fastq"; parents=String[], parameters=(source=provenance_source, alphabet=string(alphabet), record_count=length(records), hash=provenance_hash))
        for (index, record) in enumerate(records)
            register_container_provenance!(_ctx, record, "read_fastq_record"; parents=[root.id], parameters=(source=provenance_source, record_index=index, identifier=record.identifier, alphabet=string(alphabet)), provenance_hash=provenance_hash)
        end
    end

    return records
end

"""
    each_fasta_record(io::IO; alphabet=DNAAlphabet, preserve_case=false)
    each_fasta_record(path::AbstractString; alphabet=DNAAlphabet, preserve_case=false)

Stream FASTA records one at a time without loading the full file into memory.
"""
function each_fasta_record(io::IO; alphabet::Type{A}=DNAAlphabet, preserve_case::Bool=false) where {A <: BioAlphabet}
    return Channel{SeqRecord{A}}(32) do ch
        header = ""
        seq_buf = UInt8[]
        seen_header = false
        
        for line in eachline(io)
            stripped = strip(line)
            isempty(stripped) && continue
            if startswith(stripped, '>')
                if seen_header
                    sequence = BioSequence{A}(copy(seq_buf); validate=false)
                    put!(ch, SeqRecord(sequence, identifier=String(header)))
                    empty!(seq_buf)
                end
                header = String(strip(line[2:end]))
                seen_header = true
            elseif seen_header
                for b in codeunits(stripped)
                    if (!preserve_case && b >= 0x61 && b <= 0x7a)
                        push!(seq_buf, b - 0x20)
                    elseif b != UInt8(' ') && b != UInt8('\t')
                        push!(seq_buf, b)
                    end
                end
            end
        end
        if seen_header
            sequence = BioSequence{A}(copy(seq_buf); validate=false)
            put!(ch, SeqRecord(sequence, identifier=String(header)))
        end
    end
end

function each_fasta_record(path::AbstractString; alphabet::Type{A}=DNAAlphabet, preserve_case::Bool=false) where {A <: BioAlphabet}
    return Channel{SeqRecord{A}}(32) do ch
        open(path, "r") do io
            for record in each_fasta_record(io; alphabet=alphabet, preserve_case=preserve_case)
                put!(ch, record)
            end
        end
    end
end

@inline function _fastq_identifier(header::String)
    space_index = findfirst(isspace, header)
    space_index === nothing && return String(header)
    return String(header[firstindex(header):prevind(header, space_index)])
end

@inline function _fastq_quality_string(quality)
    quality isa String && return String(quality)
    return String(UInt8[UInt8(character) for character in quality])
end

@inline function _fastq_components(record)
    throw(ArgumentError("unsupported FASTQ record type: $(typeof(record))"))
end

"""
    write_fastq(path, records)
    write_fastq(io::IO, records)
"""
function write_fastq(io::IO, records)
    buf = Vector{UInt8}(undef, 65536)
    buf_pos = 0

    @inline function flush_buf()
        if buf_pos > 0
            write(io, view(buf, 1:buf_pos))
            buf_pos = 0
        end
    end

    @inline function write_byte(b::UInt8)
        buf_pos += 1
        @inbounds buf[buf_pos] = b
        if buf_pos == length(buf)
            flush_buf()
        end
    end

    @inline function write_bytes(data::AbstractVector{UInt8})
        n = length(data)
        src_pos = 1
        while src_pos <= n
            avail = length(buf) - buf_pos
            to_copy = min(avail, n - src_pos + 1)
            copyto!(buf, buf_pos + 1, data, src_pos, to_copy)
            buf_pos += to_copy
            src_pos += to_copy
            if buf_pos == length(buf)
                flush_buf()
            end
        end
    end

    for record in records
        identifier, description, sequence, quality = _fastq_components(record)
        header = isempty(description) ? identifier : description
        write_byte(UInt8('@'))
        write_bytes(codeunits(header))
        write_byte(UInt8('\n'))

        seq_bytes = sequence isa BioSequence ? sequence.data : codeunits(String(sequence))
        write_bytes(seq_bytes)

        write_byte(UInt8('\n'))
        write_byte(UInt8('+'))
        write_byte(UInt8('\n'))

        qual_bytes = codeunits(quality)
        write_bytes(qual_bytes)
        write_byte(UInt8('\n'))
    end
    flush_buf()
    return io
end

function write_fastq(path::String, records; prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    materialized = _materialize_records_for_provenance(records, _ctx)
    open(path, "w") do io
        write_fastq(io, materialized)
    end
    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(path)))
        _register_path_write_provenance!(_ctx, "write_fastq", path, provenance_parent_ids(materialized); record_count=length(materialized), provenance_hash=provenance_hash)
    end
    
    return path
end

# ─── Variant/Interval I/O (VCF, BED, GFF3) ────────────────────────────────────

function _load_optional_module(module_name::Symbol)
    isdefined(@__MODULE__, module_name) && return getfield(@__MODULE__, module_name)
    try
        Base.eval(@__MODULE__, Expr(:import, module_name))
    catch
        return nothing
    end
    return getfield(@__MODULE__, module_name)
end

function _open_vcf_input(path::String, f)
    if endswith(lowercase(path), ".gz")
        codec = _load_optional_module(:CodecZlib)
        codec === nothing && throw(ArgumentError("reading .vcf.gz requires CodecZlib.jl"))
        open(path, "r") do raw
            stream = codec.GzipDecompressorStream(raw)
            try
                return f(stream)
            finally
                close(stream)
            end
        end
    end

    open(path, "r") do io
        return f(io)
    end
end

function _open_vcf_output(path::String, f)
    if endswith(lowercase(path), ".gz")
        codec = _load_optional_module(:CodecZlib)
        codec === nothing && throw(ArgumentError("writing .vcf.gz requires CodecZlib.jl"))
        open(path, "w") do raw
            stream = codec.GzipCompressorStream(raw)
            try
                return f(stream)
            finally
                close(stream)
            end
        end
    end

    open(path, "w") do io
        return f(io)
    end
end



@inline function _parse_vcf_qual(field::AbstractString)
    stripped = strip(String(field))
    isempty(stripped) && return missing
    stripped == "." && return missing
    parsed = tryparse(Float64, stripped)
    parsed === nothing && throw(ArgumentError("invalid VCF QUAL value: $(field)"))
    return Float32(parsed)
end

function _vcf_header_columns(header::Union{Nothing,VcfHeader}, sample_names::Vector{String})
    if header === nothing || isempty(header.columns)
        columns = String["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    else
        columns = copy(header.columns)
    end

    if isempty(sample_names)
        return columns
    end

    if length(columns) == 8
        return vcat(columns, ["FORMAT"], sample_names)
    elseif length(columns) == 9
        return vcat(columns, sample_names)
    elseif length(columns) >= 10
        return columns
    end

    return vcat(String["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"], ["FORMAT"], sample_names)
end

function _vcf_sample_names(header::Union{Nothing,VcfHeader}, records::AbstractVector{<:VariantTextRecord})
    if header !== nothing && !isempty(header.sample_names)
        return copy(header.sample_names)
    end

    sample_count = isempty(records) ? 0 : maximum(length(record.samples) for record in records)
    return ["sample_$(index)" for index in 1:sample_count]
end

function _write_vcf_header(io::IO, header::Union{Nothing,VcfHeader}, sample_names::Vector{String})
    if header === nothing || isempty(header.meta_lines)
        println(io, "##fileformat=VCFv4.2")
    else
        has_fileformat = any(startswith(line, "##fileformat=") for line in header.meta_lines)
        for line in header.meta_lines
            println(io, line)
        end
        has_fileformat || println(io, "##fileformat=VCFv4.2")
    end

    println(io, join(_vcf_header_columns(header, sample_names), '\t'))
    return nothing
end

function _vcf_record_fields(record::VariantTextRecord, sample_count::Int)
    fields = String[
        record.chrom,
        string(record.pos),
        record.id,
        record.ref,
        record.alt,
        record.qual === missing ? "." : string(record.qual),
        isempty(record.filter) ? "PASS" : record.filter,
        isempty(record.info) ? "." : record.info,
    ]

    if sample_count > 0
        push!(fields, isempty(record.format) ? "GT" : record.format)
        samples = copy(record.samples)
        if length(samples) > sample_count
            throw(DimensionMismatch("VCF record sample count mismatch"))
        elseif length(samples) < sample_count
            append!(samples, fill(".", sample_count - length(samples)))
        end
        append!(fields, samples)
    end

    return fields
end

function _write_vcf_records(io::IO, records::AbstractVector{<:VariantTextRecord}; header::Union{Nothing,VcfHeader}=nothing)
    sample_names = _vcf_sample_names(header, records)
    _write_vcf_header(io, header, sample_names)
    sample_count = length(sample_names)

    buf = Vector{UInt8}(undef, 65536)
    buf_pos = 0

    @inline function flush_buf()
        if buf_pos > 0
            write(io, view(buf, 1:buf_pos))
            buf_pos = 0
        end
    end

    @inline function write_byte(b::UInt8)
        buf_pos += 1
        @inbounds buf[buf_pos] = b
        if buf_pos == length(buf)
            flush_buf()
        end
    end

    @inline function write_bytes(data::AbstractVector{UInt8})
        n = length(data)
        src_pos = 1
        while src_pos <= n
            avail = length(buf) - buf_pos
            to_copy = min(avail, n - src_pos + 1)
            copyto!(buf, buf_pos + 1, data, src_pos, to_copy)
            buf_pos += to_copy
            src_pos += to_copy
            if buf_pos == length(buf)
                flush_buf()
            end
        end
    end

    for record in records
        write_bytes(codeunits(record.chrom))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(string(record.pos)))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(record.id))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(record.ref))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(record.alt))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(record.qual === missing ? "." : string(record.qual)))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(isempty(record.filter) ? "PASS" : record.filter))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(isempty(record.info) ? "." : record.info))

        if sample_count > 0
            write_byte(UInt8('\t'))
            write_bytes(codeunits(isempty(record.format) ? "GT" : record.format))
            samples = record.samples
            for (idx, s) in enumerate(samples)
                write_byte(UInt8('\t'))
                write_bytes(codeunits(s))
            end
            if length(samples) < sample_count
                for _ in (length(samples)+1):sample_count
                    write_byte(UInt8('\t'))
                    write_byte(UInt8('.'))
                end
            end
        end
        write_byte(UInt8('\n'))
    end
    flush_buf()
    return nothing
end

"""
    parse_vcf_record(line::AbstractString) -> Union{Nothing, VariantTextRecord}

Parse a single VCF variant record line into a structured record. By design, accepts 
truncated VCF lines (with a minimum of CHROM, POS, ID, REF, ALT, QUAL) and populates 
unspecified trailing fields with standard VCF defaults (`QUAL=missing`, `FILTER="PASS"`, `INFO="."`).
"""
function parse_vcf_record(line::AbstractString; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    s = strip(line)
    isempty(s) && return nothing
    
    idx1 = findnext('\t', s, 1)
    idx1 === nothing && return nothing
    chrom = String(SubString(s, 1, idx1 - 1))
    
    idx2 = findnext('\t', s, idx1 + 1)
    idx2 === nothing && return nothing
    pos_str = SubString(s, idx1 + 1, idx2 - 1)
    pos_val = tryparse(Int, pos_str)
    (pos_val === nothing || pos_val <= 0) && throw(ArgumentError("VCF position must be positive"))
    pos = Int32(pos_val)
    
    idx3 = findnext('\t', s, idx2 + 1)
    idx3 === nothing && return nothing
    id = String(SubString(s, idx2 + 1, idx3 - 1))
    
    idx4 = findnext('\t', s, idx3 + 1)
    idx4 === nothing && return nothing
    ref = String(SubString(s, idx3 + 1, idx4 - 1))
    
    idx5 = findnext('\t', s, idx4 + 1)
    idx5 === nothing && return nothing
    alt = String(SubString(s, idx4 + 1, idx5 - 1))
    
    idx6 = findnext('\t', s, idx5 + 1)
    if idx6 === nothing
        qual_str = SubString(s, idx5 + 1)
        qual = _parse_vcf_qual(qual_str)
        return VariantTextRecord(chrom, pos, id, ref, alt, qual; filter="PASS", info=".", format="", samples=String[])
    end
    qual_str = SubString(s, idx5 + 1, idx6 - 1)
    qual = _parse_vcf_qual(qual_str)
    
    idx7 = findnext('\t', s, idx6 + 1)
    if idx7 === nothing
        filter_str = SubString(s, idx6 + 1)
        filter = isempty(filter_str) ? "PASS" : String(filter_str)
        return VariantTextRecord(chrom, pos, id, ref, alt, qual; filter=filter, info=".", format="", samples=String[])
    end
    filter_str = SubString(s, idx6 + 1, idx7 - 1)
    filter = isempty(filter_str) ? "PASS" : String(filter_str)
    
    idx8 = findnext('\t', s, idx7 + 1)
    if idx8 === nothing
        info_str = SubString(s, idx7 + 1)
        info = isempty(info_str) ? "." : String(info_str)
        return VariantTextRecord(chrom, pos, id, ref, alt, qual; filter=filter, info=info, format="", samples=String[])
    end
    info_str = SubString(s, idx7 + 1, idx8 - 1)
    info = isempty(info_str) ? "." : String(info_str)
    
    idx9 = findnext('\t', s, idx8 + 1)
    if idx9 === nothing
        format_str = SubString(s, idx8 + 1)
        return VariantTextRecord(chrom, pos, id, ref, alt, qual; filter=filter, info=info, format=String(format_str), samples=String[])
    end
    format_str = SubString(s, idx8 + 1, idx9 - 1)
    
    samples_str = SubString(s, idx9 + 1)
    samples = isempty(samples_str) ? String[] : String[String(sub) for sub in Base.eachsplit(samples_str, '\t'; keepempty=true)]
    return VariantTextRecord(chrom, pos, id, ref, alt, qual; filter=filter, info=info, format=String(format_str), samples=samples)
end

"""
    read_vcf_document(input)

Read a full VCF document, including header metadata and parsed records.
"""
function read_vcf_document(input::String; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    if _ctx === nothing
        return _open_vcf_input(input, io -> begin
            return read_vcf_document(io, nothing, input; _ctx=nothing)
        end)
    end
    raw_bytes = read(input)
    provenance_hash = bytes2hex(sha256(raw_bytes))
    return _open_vcf_input(input, io -> begin
        return read_vcf_document(io, provenance_hash, input; _ctx=_ctx)
    end)
end

"""
    read_vcf_document(io)

Read a full VCF document, including header metadata and parsed records.
"""
function read_vcf_document(io::IO, provenance_hash::Union{Nothing,AbstractString}=nothing, provenance_source::AbstractString="read_vcf_document"; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    meta_lines = String[]
    columns = String[]
    records = VariantTextRecord[]

    for raw_line in eachline(io)
        line = strip(raw_line)
        isempty(line) && continue
        if startswith(line, "##")
            push!(meta_lines, line)
            continue
        elseif startswith(line, "#CHROM")
            columns = String[String(sub) for sub in Base.eachsplit(line, '\t'; keepempty=true)]
            continue
        end

        record = parse_vcf_record(line)
        record === nothing && continue
        push!(records, record)
    end

    sample_names = length(columns) >= 10 ? String.(columns[10:end]) : String[]
    inferred_samples = isempty(records) ? 0 : maximum(length(record.samples) for record in records)
    if isempty(columns)
        if inferred_samples > 0
            sample_names = ["sample_$(index)" for index in 1:inferred_samples]
            columns = vcat(String["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"], sample_names)
        else
            columns = String["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        end
    elseif length(columns) == 8 && inferred_samples > 0
        sample_names = ["sample_$(index)" for index in 1:inferred_samples]
        columns = vcat(columns, ["FORMAT"], sample_names)
    elseif length(columns) == 9 && inferred_samples > 0 && isempty(sample_names)
        sample_names = ["sample_$(index)" for index in 1:inferred_samples]
        columns = vcat(columns, sample_names)
    end

    header = VcfHeader(meta_lines, columns, sample_names)
    document = VcfDocument(header, records)
    if provenance_hash !== nothing
        document.metadata[PROVENANCE_HASH_KEY] = String(provenance_hash)
        for record in document.records
            record.metadata[PROVENANCE_HASH_KEY] = String(provenance_hash)
        end
    end
    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        root = register_provenance!(_ctx, "read_vcf_document"; parents=String[], parameters=(source=provenance_source, record_count=length(records), sample_count=length(sample_names), hash=provenance_hash))
        register_container_provenance!(_ctx, document, "read_vcf_document"; parents=[root.id], parameters=(source=provenance_source, record_count=length(records), sample_count=length(sample_names)), provenance_hash=provenance_hash)
        for (index, record) in enumerate(document.records)
        register_container_provenance!(_ctx, record, "read_vcf_record"; parents=[root.id], parameters=(source=provenance_source, record_index=index, chrom=record.chrom, pos=record.pos), provenance_hash=provenance_hash)
        end
    end
    return document
end

"""
    read_vcf(input)

Read VCF records from a file path or IO stream.
"""
function read_vcf(input::String; prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)


    return read_vcf_document(input; _ctx=_ctx).records
end

function read_vcf(io::IO; prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)


    return read_vcf_document(io; _ctx=_ctx).records
end

"""
    write_vcf_document(output, doc)

Write a full VCF document to a file path.
"""
function write_vcf_document(output::String, doc::VcfDocument; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    result = _open_vcf_output(output, io -> begin
        return write_vcf_document(io, doc; _ctx=nothing)
    end)
    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(output)))
        _register_path_write_provenance!(_ctx, "write_vcf_document", output, provenance_parent_ids(doc, doc.records); record_count=length(doc.records), provenance_hash=provenance_hash)
    end
    
    return result
end

"""
    write_vcf_document(io, doc)

Write a full VCF document to an IO stream.
"""
function write_vcf_document(io::IO, doc::VcfDocument; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    _write_vcf_records(io, doc.records; header=doc.header)
    _ctx = active_provenance_context(_ctx)
    return nothing
end

"""
    write_vcf(output, records)

Write VCF records to a file path.
"""
function write_vcf(output::String, records::AbstractVector{<:VariantTextRecord}; header::Union{Nothing,VcfHeader}=nothing, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    materialized = _materialize_records_for_provenance(records, _ctx)
    result = _open_vcf_output(output, io -> begin
        return write_vcf(io, materialized; header=header, _ctx=_ctx)
    end)
    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(output)))
        _register_path_write_provenance!(_ctx, "write_vcf", output, provenance_parent_ids(materialized); record_count=length(materialized), provenance_hash=provenance_hash)
    end
    
    return result
end

function write_vcf(io::IO, records::AbstractVector{<:VariantTextRecord}; header::Union{Nothing,VcfHeader}=nothing, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    _write_vcf_records(io, records; header=header)
    _ctx = active_provenance_context(_ctx)
    return nothing
end

function write_vcf(output::String, doc::VcfDocument; prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)


    return write_vcf_document(output, doc; _ctx=_ctx)
end

function write_vcf(io::IO, doc::VcfDocument; prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    return write_vcf_document(io, doc; _ctx=nothing)
end

"""
    BedRecord

Compact BED interval record.
"""
struct BedRecord
    chrom::String
    start::Int32
    stop::Int32
end

"""
    GffRecord

Structured GFF3 record with parsed attributes.
"""
struct GffRecord
    chrom::String
    source::String
    feature::String
    start::Int32
    stop::Int32
    score::Union{Missing,Float32}
    strand::String
    phase::Union{Missing,Int8}
    attributes::String
    attribute_map::Dict{String,Vector{String}}
end

"""
    Base.:(==)(left, right)

Test two GFF records for value equality.
"""
function Base.:(==)(left::GffRecord, right::GffRecord)
    return isequal(left.chrom, right.chrom) && isequal(left.source, right.source) && isequal(left.feature, right.feature) && isequal(left.start, right.start) && isequal(left.stop, right.stop) && isequal(left.score, right.score) && isequal(left.strand, right.strand) && isequal(left.phase, right.phase) && isequal(left.attributes, right.attributes) && isequal(left.attribute_map, right.attribute_map)
end

"""
    _parse_gff_attributes(attributes)

Parse the attribute column of a GFF record into a dictionary.
"""
function _parse_gff_attributes(attributes::AbstractString)
    parsed = Dict{String,Vector{String}}()
    s = strip(attributes)
    isempty(s) && return parsed

    start_idx = 1
    len = ncodeunits(s)
    while start_idx <= len
        semi_idx = findnext(';', s, start_idx)
        pair_end = semi_idx === nothing ? len : semi_idx - 1
        pair = strip(SubString(s, start_idx, pair_end))
        if !isempty(pair)
            eq_idx = findnext('=', pair, 1)
            if eq_idx !== nothing
                key = String(strip(SubString(pair, 1, eq_idx - 1)))
                val_str = strip(SubString(pair, eq_idx + 1))
                if isempty(val_str)
                    parsed[key] = String[""]
                else
                    vals = String[]
                    vstart = 1
                    vlen = ncodeunits(val_str)
                    while vstart <= vlen
                        comma_idx = findnext(',', val_str, vstart)
                        vend = comma_idx === nothing ? vlen : comma_idx - 1
                        push!(vals, String(strip(SubString(val_str, vstart, vend))))
                        vstart = comma_idx === nothing ? vlen + 1 : comma_idx + 1
                    end
                    parsed[key] = vals
                end
            else
                parsed[String(pair)] = String[""]
            end
        end
        start_idx = semi_idx === nothing ? len + 1 : semi_idx + 1
    end

    return parsed
end

"""
    _render_gff_attributes(attributes, attribute_map)

Render GFF attributes back into a text column.
"""
function _render_gff_attributes(attributes::String, attribute_map::Dict{String,Vector{String}})
    isempty(strip(attributes)) || return attributes
    isempty(attribute_map) && return "."

    parts = String[]
    for (key, values) in attribute_map
        isempty(values) ? push!(parts, key) : push!(parts, string(key, "=", join(values, ',')))
    end

    return join(parts, ';')
end

"""
    GffRecord(chrom, source, feature, start, stop, score, strand, phase, attributes, attribute_map=...)

Construct a normalized GFF record from text or parsed components.
"""
function GffRecord(
    chrom::String,
    source::String,
    feature::String,
    start::Integer,
    stop::Integer,
    score::Union{Missing,Float32},
    strand::String,
    phase::Union{Missing,Int8},
    attributes::String,
    attribute_map::AbstractDict=Dict{String,Vector{String}}())
    return GffRecord(
        String(chrom),
        String(source),
        String(feature),
        Int32(start),
        Int32(stop),
        score,
        String(strand),
        phase,
        String(attributes),
        Dict{String,Vector{String}}(attribute_map))
end

"""
    GenBankFeature

Parsed GenBank feature with raw and parsed location information.
"""
struct GenBankFeature
    key::String
    location::String
    qualifiers::Dict{String,Vector{String}}
    parsed_location::Any
end

"""
    Base.:(==)(left, right)

Test two GenBank features for value equality.
"""
function Base.:(==)(left::GenBankFeature, right::GenBankFeature)
    return isequal(left.key, right.key) && isequal(left.location, right.location) && isequal(left.qualifiers, right.qualifiers) && isequal(left.parsed_location, right.parsed_location)
end

GenBankFeature(key::String, location::String, qualifiers::AbstractDict=Dict{String,Vector{String}}()) = GenBankFeature(
    String(key),
    String(location),
    Dict{String,Vector{String}}(qualifiers),
    nothing)

GenBankFeature(key::String, location::String, qualifiers::AbstractDict, parsed_location) = GenBankFeature(
    String(key),
    String(location),
    Dict{String,Vector{String}}(qualifiers),
    parsed_location)

"""
    GenBankRecord

Structured GenBank record containing metadata, sequence, and features.
"""
struct GenBankRecord
    locus::String
    locus_line::String
    definition::String
    accession::String
    version::String
    keywords::String
    source::String
    organism::String
    comment::String
    sequence::BioSequence
    features::Vector{GenBankFeature}
    metadata::Dict{Symbol,Any}
end

"""
    Base.:(==)(left, right)

Test two GenBank records for value equality.
"""
function Base.:(==)(left::GenBankRecord, right::GenBankRecord)
    return isequal(left.locus, right.locus) && isequal(left.locus_line, right.locus_line) && isequal(left.definition, right.definition) && isequal(left.accession, right.accession) && isequal(left.version, right.version) && isequal(left.keywords, right.keywords) && isequal(left.source, right.source) && isequal(left.organism, right.organism) && isequal(left.comment, right.comment) && isequal(left.sequence, right.sequence) && isequal(left.features, right.features) && isequal(left.metadata, right.metadata)
end

function Base.show(io::IO, record::GenBankRecord)
    print(io, "GenBankRecord(", record.locus, ", ", length(record.sequence), " bp, ", container_provenance_summary(record), ")")
end

"""
    GenBankRecord(...)

Construct a normalized GenBank record from parsed components.
"""
function GenBankRecord(
    locus::AbstractString,
    definition::AbstractString,
    accession::AbstractString,
    version::AbstractString,
    source::AbstractString,
    organism::AbstractString,
    sequence,
    features::AbstractVector{GenBankFeature},
    locus_line::AbstractString="",
    keywords::AbstractString="",
    comment::AbstractString="",
    metadata::AbstractDict=Dict{Symbol,Any}())
    sequence_typed = if sequence isa BioSequence
        sequence
    else
        sequence_text = String(sequence)
        seq_bytes = Vector{UInt8}(sequence_text)
        if isempty(sequence_text) || validate_sequence(DNAAlphabet, sequence_text)
            BioSequence{DNAAlphabet}(seq_bytes; validate=false)
        elseif validate_sequence(RNAAlphabet, sequence_text)
            BioSequence{RNAAlphabet}(seq_bytes; validate=false)
        elseif validate_sequence(AminoAcidAlphabet, sequence_text)
            BioSequence{AminoAcidAlphabet}(seq_bytes; validate=false)
        else
            throw(ArgumentError("cannot infer alphabet for GenBank sequence"))
        end
    end

    return GenBankRecord(
        String(locus),
        String(locus_line),
        String(definition),
        String(accession),
        String(version),
        String(keywords),
        String(source),
        String(organism),
        String(comment),
        sequence_typed,
        features isa Vector{GenBankFeature} ? features : GenBankFeature[feature for feature in features],
        begin
            metadata_copy = Dict{Symbol,Any}(metadata)
            ensure_provenance_id!(metadata_copy)
            metadata_copy
        end)
end

"""
    GenBankArrowRecord

Flattened GenBank record used for Arrow ingestion.
"""
struct GenBankArrowRecord
    locus::String
    locus_line::String
    accession::String
    version::String
    definition::String
    keywords::String
    source::String
    organism::String
    comment::String
    sequence_text::String
    feature_count::Int32
    feature_keys::String
    feature_locations::String
end

const _GENBANK_FEATURE_KEY_PREFIX = "     "
const _GENBANK_FEATURE_QUALIFIER_PREFIX = "                     "

"""
    parse_bed_record(line)

Parse a single BED text line into a BED record.
"""
function parse_bed_record(line::AbstractString; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    s = strip(line)
    (isempty(s) || startswith(s, '#') || startswith(s, "track") || startswith(s, "browser")) && return nothing

    idx1 = findnext('\t', s, 1)
    idx1 === nothing && return nothing
    chrom = String(SubString(s, 1, idx1 - 1))

    idx2 = findnext('\t', s, idx1 + 1)
    idx2 === nothing && return nothing
    start_str = SubString(s, idx1 + 1, idx2 - 1)
    start_val = tryparse(Int, start_str)
    start_val === nothing && throw(ArgumentError("malformed BED record: $(line)"))

    idx3 = findnext('\t', s, idx2 + 1)
    stop_str = idx3 === nothing ? SubString(s, idx2 + 1) : SubString(s, idx2 + 1, idx3 - 1)
    stop_val = tryparse(Int, stop_str)
    stop_val === nothing && throw(ArgumentError("malformed BED record: $(line)"))

    start_val >= 0 || throw(ArgumentError("BED start must be nonnegative"))
    stop_val > start_val || throw(ArgumentError("BED stop must be greater than start"))
    return BedRecord(chrom, Int32(start_val), Int32(stop_val))
end

"""
    read_bed(input)

Read BED records from a file path.
"""
function read_bed(input::String; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    if _ctx === nothing
        open(input, "r") do io
            return read_bed(io)
        end
    end
    raw_bytes = read(input)
    provenance_hash = bytes2hex(sha256(raw_bytes))
    return read_bed(IOBuffer(raw_bytes), provenance_hash, input; _ctx=_ctx)
end

"""
    read_bed(io)

Read BED records from an IO stream.
"""
function read_bed(io::IO, provenance_hash::Union{Nothing,AbstractString}=nothing, provenance_source::AbstractString="read_bed"; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    records = BedRecord[]

    for (line_number, raw_line) in enumerate(eachline(io))
        startswith(raw_line, '#') && continue
        record = try
            parse_bed_record(raw_line)
        catch err
            throw(ArgumentError("malformed BED record on line $(line_number): $(err.msg)"))
        end
        record === nothing && continue
        push!(records, record)
    end

    _ctx = active_provenance_context(_ctx)
    _ctx !== nothing && register_provenance!(_ctx, "read_bed"; parents=String[], parameters=(source=provenance_source, record_count=length(records), hash=provenance_hash))
    return records
end

"""
    write_bed(output, records)

Write BED records to a file path.
"""
function write_bed(output::String, records; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    materialized = _materialize_records_for_provenance(records, _ctx)
    open(output, "w") do io
        write_bed(io, materialized; _ctx=nothing)
    end

    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(output)))
        _register_path_write_provenance!(_ctx, "write_bed", output, provenance_parent_ids(materialized); record_count=length(materialized), provenance_hash=provenance_hash)
    end
    
    return output
end

"""
    write_bed(io, records)

Write BED records to an IO stream.
"""
function write_bed(io::IO, records; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    buf = Vector{UInt8}(undef, 65536)
    buf_pos = 0

    @inline function flush_buf()
        if buf_pos > 0
            write(io, view(buf, 1:buf_pos))
            buf_pos = 0
        end
    end

    @inline function write_byte(b::UInt8)
        buf_pos += 1
        @inbounds buf[buf_pos] = b
        if buf_pos == length(buf)
            flush_buf()
        end
    end

    @inline function write_bytes(data::AbstractVector{UInt8})
        n = length(data)
        src_pos = 1
        while src_pos <= n
            avail = length(buf) - buf_pos
            to_copy = min(avail, n - src_pos + 1)
            copyto!(buf, buf_pos + 1, data, src_pos, to_copy)
            buf_pos += to_copy
            src_pos += to_copy
            if buf_pos == length(buf)
                flush_buf()
            end
        end
    end

    for record in records
        write_bytes(codeunits(record.chrom))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(string(record.start)))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(string(record.stop)))
        write_byte(UInt8('\n'))
    end
    flush_buf()
    _ctx = active_provenance_context(_ctx)
    return nothing
end

"""
    parse_gff_record(line)

Parse a single GFF3 text line into a structured record.
"""
function parse_gff_record(line::String; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    first_index = firstindex(line)
    tab1 = findnext('\t', line, first_index)
    tab1 === nothing && return nothing
    tab2 = findnext('\t', line, nextind(line, tab1))
    tab2 === nothing && return nothing
    tab3 = findnext('\t', line, nextind(line, tab2))
    tab3 === nothing && return nothing
    tab4 = findnext('\t', line, nextind(line, tab3))
    tab4 === nothing && return nothing
    tab5 = findnext('\t', line, nextind(line, tab4))
    tab5 === nothing && return nothing
    tab6 = findnext('\t', line, nextind(line, tab5))
    tab6 === nothing && return nothing
    tab7 = findnext('\t', line, nextind(line, tab6))
    tab7 === nothing && return nothing
    tab8 = findnext('\t', line, nextind(line, tab7))
    tab8 === nothing && return nothing

    try
        chrom = SubString(line, first_index, prevind(line, tab1))
        source = SubString(line, nextind(line, tab1), prevind(line, tab2))
        feature = SubString(line, nextind(line, tab2), prevind(line, tab3))
        start = parse(Int, SubString(line, nextind(line, tab3), prevind(line, tab4)))
        stop = parse(Int, SubString(line, nextind(line, tab4), prevind(line, tab5)))
        start > 0 || throw(ArgumentError("GFF start must be positive"))
        stop >= start || throw(ArgumentError("GFF stop must be >= start"))

        score_field = SubString(line, nextind(line, tab5), prevind(line, tab6))
        score = score_field == "." ? missing : Float32(parse(Float64, score_field))

        strand = SubString(line, nextind(line, tab6), prevind(line, tab7))
        strand in ("+", "-", ".") || throw(ArgumentError("GFF strand must be +, -, or ."))

        phase_field = SubString(line, nextind(line, tab7), prevind(line, tab8))
        if phase_field == "."
            phase = missing
        else
            phase_value = parse(Int, phase_field)
            0 <= phase_value <= 2 || throw(ArgumentError("GFF phase must be 0, 1, or 2"))
            phase = Int8(phase_value)
        end

        attributes = SubString(line, nextind(line, tab8), lastindex(line))
        attribute_map = _parse_gff_attributes(attributes)
        return GffRecord(chrom, source, feature, Int32(start), Int32(stop), score, strand, phase, attributes, attribute_map)
    catch err
        err isa ArgumentError && rethrow()
        throw(ArgumentError("malformed GFF record: $(line)"))
    end
end

"""
    read_gff(input)

Read GFF records from a file path.
"""
function read_gff(input::String; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    if _ctx === nothing
        open(input, "r") do io
            return read_gff(io)
        end
    end
    raw_bytes = read(input)
    provenance_hash = bytes2hex(sha256(raw_bytes))
    return read_gff(IOBuffer(raw_bytes), provenance_hash, input; _ctx=_ctx)
end

"""
    read_gff(io)

Read GFF records from an IO stream.
"""
function read_gff(io::IO, provenance_hash::Union{Nothing,AbstractString}=nothing, provenance_source::AbstractString="read_gff"; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    records = GffRecord[]

    for (line_number, raw_line) in enumerate(eachline(io))
        startswith(raw_line, '#') && continue
        record = try
            parse_gff_record(raw_line)
        catch err
            throw(ArgumentError("malformed GFF record on line $(line_number): $(err.msg)"))
        end
        record === nothing && continue
        push!(records, record)
    end

    _ctx = active_provenance_context(_ctx)
    _ctx !== nothing && register_provenance!(_ctx, "read_gff"; parents=String[], parameters=(source=provenance_source, record_count=length(records), hash=provenance_hash))
    return records
end

"""
    write_gff(output, records)

Write GFF records to a file path.
"""
function write_gff(output::String, records; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    materialized = _materialize_records_for_provenance(records, _ctx)
    open(output, "w") do io
        write_gff(io, materialized; _ctx=nothing)
    end

    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(output)))
        _register_path_write_provenance!(_ctx, "write_gff", output, provenance_parent_ids(materialized); record_count=length(materialized), provenance_hash=provenance_hash)
    end
    
    return output
end

"""
    write_gff(io, records)

Write GFF records to an IO stream.
"""
function write_gff(io::IO, records; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    buf = Vector{UInt8}(undef, 65536)
    buf_pos = 0

    @inline function flush_buf()
        if buf_pos > 0
            write(io, view(buf, 1:buf_pos))
            buf_pos = 0
        end
    end

    @inline function write_byte(b::UInt8)
        buf_pos += 1
        @inbounds buf[buf_pos] = b
        if buf_pos == length(buf)
            flush_buf()
        end
    end

    @inline function write_bytes(data::AbstractVector{UInt8})
        n = length(data)
        src_pos = 1
        while src_pos <= n
            avail = length(buf) - buf_pos
            to_copy = min(avail, n - src_pos + 1)
            copyto!(buf, buf_pos + 1, data, src_pos, to_copy)
            buf_pos += to_copy
            src_pos += to_copy
            if buf_pos == length(buf)
                flush_buf()
            end
        end
    end

    for record in records
        score = record.score === missing ? "." : string(record.score)
        phase = record.phase === missing ? "." : string(record.phase)
        attributes = _render_gff_attributes(record.attributes, record.attribute_map)
        
        write_bytes(codeunits(record.chrom))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(record.source))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(record.feature))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(string(record.start)))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(string(record.stop)))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(score))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(record.strand))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(phase))
        write_byte(UInt8('\t'))
        write_bytes(codeunits(attributes))
        write_byte(UInt8('\n'))
    end
    flush_buf()
    _ctx = active_provenance_context(_ctx)
    return nothing
end

"""
    _genbank_push_qualifier!(qualifiers, key, value)

Append a qualifier value to a GenBank qualifier dictionary.
"""
function _genbank_push_qualifier!(qualifiers::Dict{String,Vector{String}}, key::String, value::String)
    push!(get!(qualifiers, key, String[]), value)
    return nothing
end

"""
    _genbank_parse_qualifier(text)

Parse a single GenBank qualifier line into a key-value pair.
"""
function _genbank_parse_qualifier(text::AbstractString)
    stripped = strip(text)
    startswith(stripped, "/") || return nothing

    payload = SubString(stripped, 2)
    eq_idx = findnext('=', payload, 1)
    if eq_idx !== nothing
        key = strip(SubString(payload, 1, eq_idx - 1))
        val_str = strip(SubString(payload, eq_idx + 1))
        if startswith(val_str, "\"") && endswith(val_str, "\"") && length(val_str) >= 2
            val_str = SubString(val_str, 2, length(val_str) - 1)
        end
        return key, val_str
    end

    return strip(payload), SubString("")
end

"""
    parse_genbank_record(lines)

Parse a GenBank flat-file record from its raw text lines.
"""
function parse_genbank_record(lines::AbstractVector{<:String}; preserve_case::Bool=false)
    isempty(lines) && return nothing

    locus = ""
    locus_line = ""
    definition = IOBuffer()
    accession = ""
    version = ""
    keywords = IOBuffer()
    source = IOBuffer()
    organism = IOBuffer()
    comment = IOBuffer()
    sequence = IOBuffer()
    features = GenBankFeature[]

    current_field = :none
    in_features = false
    in_origin = false
    current_feature_key = ""
    current_feature_location = IOBuffer()
    current_feature_qualifiers = Dict{String,Vector{String}}()
    current_qualifier_key = nothing

    function flush_feature!()
        if !isempty(current_feature_key)
            qualifiers = Dict{String,Vector{String}}()
            for (key, values) in current_feature_qualifiers
                clean_values = String[]
                for val in values
                    v = strip(val)
                    if startswith(v, '"') && endswith(v, '"') && length(v) >= 2
                        v = v[2:end-1]
                    elseif startswith(v, '"')
                        v = v[2:end]
                    elseif endswith(v, '"')
                        v = v[1:end-1]
                    end
                    push!(clean_values, v)
                end
                qualifiers[key] = clean_values
            end
            location_text = strip(String(take!(current_feature_location)))
            parsed_location = isempty(location_text) ? nothing : parse_feature_location(location_text)
            push!(features, GenBankFeature(current_feature_key, location_text, qualifiers, parsed_location))
        end
        empty!(current_feature_qualifiers)
        current_feature_location = IOBuffer()
        current_feature_key = ""
        current_qualifier_key = nothing
        return nothing
    end

    for raw_line in lines
        line = rstrip(raw_line)
        line == "//" && break

        if in_origin
            for byte in codeunits(line)
                if (!preserve_case && UInt8('a') <= byte <= UInt8('z'))
                    write(sequence, UInt8(byte - 0x20))
                elseif (UInt8('a') <= byte <= UInt8('z')) || (UInt8('A') <= byte <= UInt8('Z'))
                    write(sequence, UInt8(byte))
                end
            end
            continue
        end

        if in_features
            if startswith(line, "ORIGIN")
                flush_feature!()
                in_features = false
                in_origin = true
                current_field = :origin
                continue
            end

            if startswith(line, _GENBANK_FEATURE_QUALIFIER_PREFIX)
                qualifier_line = line[length(_GENBANK_FEATURE_QUALIFIER_PREFIX) + 1:end]
                parsed = _genbank_parse_qualifier(qualifier_line)
                if parsed !== nothing
                    key, value = parsed
                    _genbank_push_qualifier!(current_feature_qualifiers, String(key), String(value))
                    current_qualifier_key = String(key)
                elseif current_qualifier_key !== nothing
                    values = current_feature_qualifiers[current_qualifier_key]
                    stripped = strip(qualifier_line)
                    if !isempty(stripped)
                        clean_str = startswith(stripped, "\"") && endswith(stripped, "\"") && length(stripped) >= 2 ? stripped[2:end-1] : (endswith(stripped, "\"") ? rstrip(stripped, '"') : stripped)
                        values[end] = string(values[end], clean_str)
                    end
                end
                continue
            end

            if startswith(line, _GENBANK_FEATURE_KEY_PREFIX)
                s = strip(line)
                sp_idx = findnext(isspace, s, 1)
                if sp_idx !== nothing
                    flush_feature!()
                    current_feature_key = String(SubString(s, 1, sp_idx - 1))
                    print(current_feature_location, strip(SubString(s, sp_idx + 1)))
                    current_qualifier_key = nothing
                    continue
                end
            end

            if current_qualifier_key !== nothing
                values = current_feature_qualifiers[current_qualifier_key]
                stripped = strip(line)
                if !isempty(stripped)
                    # For qualifiers like /translation, remove quotes or space insertion
                    clean_str = startswith(stripped, "\"") && endswith(stripped, "\"") && length(stripped) >= 2 ? stripped[2:end-1] : (endswith(stripped, "\"") ? rstrip(stripped, '"') : stripped)
                    values[end] = string(values[end], clean_str)
                end
                continue
            end

            stripped = strip(line)
            if !isempty(stripped)
                write(current_feature_location, ' ')
                write(current_feature_location, stripped)
            end
            continue
        end

        if startswith(line, "LOCUS")
            s = strip(line)
            sp1 = findnext(isspace, s, 1)
            if sp1 !== nothing
                p2 = findnext(!isspace, s, sp1 + 1)
                if p2 !== nothing
                    sp2 = findnext(isspace, s, p2)
                    locus = String(sp2 === nothing ? SubString(s, p2) : SubString(s, p2, sp2 - 1))
                end
            end
            locus_line = String(s)
            current_field = :none
            continue
        elseif startswith(line, "DEFINITION")
            print(definition, strip(line[11:end]))
            current_field = :definition
            continue
        elseif startswith(line, "ACCESSION")
            accession = strip(line[10:end])
            current_field = :accession
            continue
        elseif startswith(line, "VERSION")
            version = strip(line[8:end])
            current_field = :version
            continue
        elseif startswith(line, "KEYWORDS")
            print(keywords, strip(line[9:end]))
            current_field = :keywords
            continue
        elseif startswith(line, "SOURCE")
            print(source, strip(line[7:end]))
            current_field = :source
            continue
        elseif startswith(line, "  ORGANISM")
            print(organism, strip(line[11:end]))
            current_field = :organism
            continue
        elseif startswith(line, "COMMENT")
            print(comment, strip(line[8:end]))
            current_field = :comment
            continue
        elseif startswith(line, "FEATURES")
            in_features = true
            in_origin = false
            current_field = :features
            continue
        elseif startswith(line, "ORIGIN")
            flush_feature!()
            in_features = false
            in_origin = true
            current_field = :origin
            continue
        end

        if current_field == :definition
            print(definition, " ", strip(line))
        elseif current_field == :keywords
            print(keywords, " ", strip(line))
        elseif current_field == :source
            print(source, " ", strip(line))
        elseif current_field == :organism
            print(organism, " ", strip(line))
        elseif current_field == :comment
            print(comment, " ", strip(line))
        end
    end

    flush_feature!()

    keywords_text = strip(String(take!(keywords)))
    if endswith(keywords_text, ".")
        keywords_text = strip(keywords_text[1:end-1])
    end
    return GenBankRecord(
        locus,
        String(take!(definition)),
        accession,
        version,
        String(take!(source)),
        String(take!(organism)),
        String(take!(sequence)),
        features,
        locus_line,
        keywords_text,
        String(take!(comment)))
end

"""
    _wrap_genbank_text(prefix, text; continuation_prefix=..., width=79, empty_text="")

Wrap a GenBank text field to the expected line width.
"""
function _wrap_genbank_text(prefix::String, text::String; continuation_prefix::String=repeat(" ", ncodeunits(prefix)), width::Int=79, empty_text::String="")
    payload = strip(text)
    isempty(payload) && return isempty(empty_text) ? String[] : [string(prefix, empty_text)]

    lines = String[]
    current_prefix = String(prefix)
    current = String(prefix)
    current_length = ncodeunits(current_prefix)

    for word in Base.eachsplit(payload)
        needs_space = current_length > ncodeunits(current_prefix)
        projected = current_length + (needs_space ? 1 : 0) + ncodeunits(word)
        if projected > width && current != current_prefix
            push!(lines, current)
            current_prefix = String(continuation_prefix)
            current = string(current_prefix, word)
            current_length = ncodeunits(current)
        else
            if needs_space
                current = string(current, " ")
                current_length += 1
            end
            current = string(current, word)
            current_length += ncodeunits(word)
        end
    end

    push!(lines, current)
    return lines
end

"""
    _format_genbank_locus_line(record)

Render the LOCUS line for a GenBank record.
"""
function _format_genbank_locus_line(record::GenBankRecord)
    isempty(strip(record.locus_line)) && return string("LOCUS       ", rpad(record.locus, 16), lpad(string(length(record.sequence)), 11), " bp    DNA     linear   UNK")
    return record.locus_line
end

"""
    _write_genbank_wrapped_section(io, prefix, text; continuation_prefix=..., width=79, empty_text="")

Write a wrapped GenBank text section to an IO stream.
"""
function _write_genbank_wrapped_section(io::IO, prefix::String, text::String; continuation_prefix::String=repeat(" ", ncodeunits(prefix)), width::Int=79, empty_text::String="")
    for line in _wrap_genbank_text(prefix, text; continuation_prefix=continuation_prefix, width=width, empty_text=empty_text)
        isempty(line) && continue
        println(io, line)
    end
end

"""
    _render_genbank_qualifier(key, values)

Render a GenBank qualifier and its values as text.
"""
function _render_genbank_qualifier(key::String, values::Vector{String})
    isempty(values) && return "                     /$(key)"
    if length(values) == 1
        value = values[1]
        isempty(value) && return "                     /$(key)"
        return "                     /$(key)=\"$(value)\""
    end
    return join(("                     /$(key)=\"$(value)\"" for value in values), '\n')
end

"""
    _render_genbank_sequence(sequence)

Render a GenBank ORIGIN sequence block.
"""
function _render_genbank_sequence(sequence::Union{BioSequence, AbstractString})
    seq_bytes = sequence isa BioSequence ? sequence.data : codeunits(String(sequence))
    n = length(seq_bytes)
    buffer = IOBuffer()
    pos = 1
    while pos <= n
        print(buffer, lpad(string(pos), 9), " ")
        end_pos = min(pos + 59, n)
        p = pos
        while p <= end_pos
            block_end = min(p + 9, end_pos)
            for i in p:block_end
                b = @inbounds seq_bytes[i]
                ch = (b >= 0x41 && b <= 0x5a) ? UInt8(b + 0x20) : b
                write(buffer, ch)
            end
            write(buffer, UInt8(' '))
            p += 10
        end
        write(buffer, UInt8('\n'))
        pos += 60
    end
    return String(take!(buffer))
end

"""
    write_genbank(output_path, records)

Write GenBank records to a file path.
"""
function write_genbank(output_path::String, records; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    materialized = _materialize_records_for_provenance(records, _ctx)
    open(output_path, "w") do io
        write_genbank(io, materialized; _ctx=nothing)
    end

    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(output_path)))
        _register_path_write_provenance!(_ctx, "write_genbank", output_path, provenance_parent_ids(materialized); record_count=length(materialized), provenance_hash=provenance_hash)
    end
    
    return output_path
end

"""
    write_genbank(io, records)

Write GenBank records to an IO stream.
"""
function write_genbank(io::IO, records; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    for record in records
        println(io, _format_genbank_locus_line(record))
        _write_genbank_wrapped_section(io, "DEFINITION  ", record.definition; continuation_prefix="            ")
        println(io, "ACCESSION   ", record.accession)
        println(io, "VERSION     ", record.version)
        keyword_text = strip(record.keywords)
        isempty(keyword_text) ? println(io, "KEYWORDS    .") : _write_genbank_wrapped_section(io, "KEYWORDS    ", endswith(keyword_text, ".") ? keyword_text : string(keyword_text, "."); continuation_prefix="            ")
        _write_genbank_wrapped_section(io, "SOURCE      ", record.source; continuation_prefix="            ")
        _write_genbank_wrapped_section(io, "  ORGANISM  ", record.organism; continuation_prefix="            ")
        isempty(strip(record.comment)) || _write_genbank_wrapped_section(io, "COMMENT     ", record.comment; continuation_prefix="            ")
        println(io, "FEATURES             Location/Qualifiers")
        for feature in record.features
            println(io, "     ", rpad(feature.key, 15), feature.location)
            for (key, values) in feature.qualifiers
                for rendered in Base.eachsplit(_render_genbank_qualifier(key, values), '\n')
                    println(io, rendered)
                end
            end
        end
        println(io, "ORIGIN")
        print(io, _render_genbank_sequence(record.sequence))
        println(io, "//")
    end
    _ctx = active_provenance_context(_ctx)
    return nothing
end

"""
    _genbank_flatten_record(record)

Flatten a GenBank record into an Arrow-ready storage layout.
"""
function _genbank_flatten_record(record::GenBankRecord)
    feature_keys = join((feature.key for feature in record.features), ";")
    feature_locations = join((feature.location for feature in record.features), ";")
    return GenBankArrowRecord(
        record.locus,
        record.locus_line,
        record.accession,
        record.version,
        record.definition,
        record.keywords,
        record.source,
        record.organism,
        record.comment,
        String(record.sequence),
        Int32(length(record.features)),
        feature_keys,
        feature_locations)
end

"""
    read_genbank(input_path)

Read GenBank records from a file path.
"""
function read_genbank(input_path::String; preserve_case::Bool=false, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    if _ctx === nothing
        open(input_path, "r") do io
            return read_genbank(io, nothing, input_path; preserve_case=preserve_case, _ctx=nothing)
        end
    end
    raw_bytes = read(input_path)
    provenance_hash = bytes2hex(sha256(raw_bytes))
    return read_genbank(IOBuffer(raw_bytes), provenance_hash, input_path; preserve_case=preserve_case, _ctx=_ctx)
end

function read_genbank(io::IO, provenance_hash::Union{Nothing,AbstractString}=nothing, provenance_source::AbstractString="read_genbank"; preserve_case::Bool=false, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    records = GenBankRecord[]
    current_lines = String[]

    for raw_line in eachline(io)
        push!(current_lines, raw_line)
        if strip(raw_line) == "//"
            record = parse_genbank_record(current_lines; preserve_case=preserve_case)
            record === nothing || push!(records, record)
            empty!(current_lines)
        end
    end

    if !isempty(current_lines)
        rec = parse_genbank_record(current_lines; preserve_case=preserve_case)
        rec === nothing || push!(records, rec)
    end
    if provenance_hash !== nothing
        for record in records
            record.metadata[PROVENANCE_HASH_KEY] = String(provenance_hash)
        end
    end
    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        root = register_provenance!(_ctx, "read_genbank"; parents=String[], parameters=(source=provenance_source, record_count=length(records), hash=provenance_hash))
        for (index, record) in enumerate(records)
        register_container_provenance!(_ctx, record, "read_genbank_record"; parents=[root.id], parameters=(source=provenance_source, record_index=index, locus=record.locus, accession=record.accession), provenance_hash=provenance_hash)
        end
    end
    return records
end

"""
    _write_genbank_chunk!(writer, loci, locus_lines, accessions, versions, definitions, keywords, sources, organisms, comments, sequences, feature_counts, feature_keys, feature_locations)

Write a chunk of flattened GenBank records to an Arrow writer.
"""
function _write_genbank_chunk!(writer, loci, locus_lines, accessions, versions, definitions, keywords, sources, organisms, comments, sequences, feature_counts, feature_keys, feature_locations)
    isempty(loci) && return nothing

    Arrow.write(
        writer,
        (
            locus = copy(loci),
            locus_line = copy(locus_lines),
            accession = copy(accessions),
            version = copy(versions),
            definition = copy(definitions),
            keywords = copy(keywords),
            source = copy(sources),
            organism = copy(organisms),
            comment = copy(comments),
            sequence = copy(sequences),
            feature_count = copy(feature_counts),
            feature_keys = copy(feature_keys),
            feature_locations = copy(feature_locations)))

    empty!(loci)
    empty!(locus_lines)
    empty!(accessions)
    empty!(versions)
    empty!(definitions)
    empty!(keywords)
    empty!(sources)
    empty!(organisms)
    empty!(comments)
    empty!(sequences)
    empty!(feature_counts)
    empty!(feature_keys)
    empty!(feature_locations)

    return nothing
end

"""
    ingest_genbank(input_path, output_path; chunk_size=100)

Convert GenBank records into a chunked Arrow table.
"""
function ingest_genbank(
    input_path::String,
    output_path::String;
    chunk_size::Integer=100,
    preserve_case::Bool=false)
    loci = String[]
    locus_lines = String[]
    accessions = String[]
    versions = String[]
    definitions = String[]
    keywords = String[]
    sources = String[]
    organisms = String[]
    comments = String[]
    sequences = String[]
    feature_counts = Int32[]
    feature_keys = String[]
    feature_locations = String[]

    current_lines = String[]

    open(Arrow.Writer, output_path; file=true) do writer
        open(input_path, "r") do io
            for raw_line in eachline(io)
                push!(current_lines, raw_line)
                if strip(raw_line) == "//"
                    record = parse_genbank_record(current_lines; preserve_case=preserve_case)
                    if record !== nothing
                        flattened = _genbank_flatten_record(record)
                        push!(loci, flattened.locus)
                        push!(locus_lines, flattened.locus_line)
                        push!(accessions, flattened.accession)
                        push!(versions, flattened.version)
                        push!(definitions, flattened.definition)
                        push!(keywords, flattened.keywords)
                        push!(sources, flattened.source)
                        push!(organisms, flattened.organism)
                        push!(comments, flattened.comment)
                        push!(sequences, flattened.sequence_text)
                        push!(feature_counts, flattened.feature_count)
                        push!(feature_keys, flattened.feature_keys)
                        push!(feature_locations, flattened.feature_locations)

                        if length(loci) >= chunk_size
                            _write_genbank_chunk!(writer, loci, locus_lines, accessions, versions, definitions, keywords, sources, organisms, comments, sequences, feature_counts, feature_keys, feature_locations)
                        end
                    end
                    empty!(current_lines)
                end
            end
        end

        if !isempty(current_lines)
            record = parse_genbank_record(current_lines; preserve_case=preserve_case)
            if record !== nothing
                flattened = _genbank_flatten_record(record)
                push!(loci, flattened.locus)
                push!(locus_lines, flattened.locus_line)
                push!(accessions, flattened.accession)
                push!(versions, flattened.version)
                push!(definitions, flattened.definition)
                push!(keywords, flattened.keywords)
                push!(sources, flattened.source)
                push!(organisms, flattened.organism)
                push!(comments, flattened.comment)
                push!(sequences, flattened.sequence_text)
                push!(feature_counts, flattened.feature_count)
                push!(feature_keys, flattened.feature_keys)
                push!(feature_locations, flattened.feature_locations)
            end
        end

        _write_genbank_chunk!(writer, loci, locus_lines, accessions, versions, definitions, keywords, sources, organisms, comments, sequences, feature_counts, feature_keys, feature_locations)
    end

    return output_path
end

"""
    _write_vcf_chunk!(writer, chroms, positions, ids, refs, alts, quals, filters, infos, formats, sample_counts, samples)

Write a chunk of flattened VCF records to an Arrow writer.
"""
function _write_vcf_chunk!(writer, chroms, positions, ids, refs, alts, quals, filters, infos, formats, sample_counts, samples)
    isempty(chroms) && return nothing

    chunk = (
        chrom = copy(chroms),
        pos = copy(positions),
        id = copy(ids),
        ref = copy(refs),
        alt = copy(alts),
        qual = copy(quals),
        filter = copy(filters),
        info = copy(infos),
        format = copy(formats),
        sample_count = copy(sample_counts),
        samples = copy(samples))

    Arrow.write(
        writer,
        chunk)

    empty!(chroms)
    empty!(positions)
    empty!(ids)
    empty!(refs)
    empty!(alts)
    empty!(quals)
    empty!(filters)
    empty!(infos)
    empty!(formats)
    empty!(sample_counts)
    empty!(samples)

    return nothing
end

"""
    _write_bed_chunk!(writer, chroms, starts, stops)

Write a chunk of flattened BED records to an Arrow writer.
"""
function _write_bed_chunk!(writer, chroms, starts, stops)
    isempty(chroms) && return nothing

    Arrow.write(
        writer,
        (
            chrom = copy(chroms),
            start = copy(starts),
            stop = copy(stops)))

    empty!(chroms)
    empty!(starts)
    empty!(stops)

    return nothing
end

"""
    ingest_vcf(input_path, output_path; chunk_size=10_000)

Convert VCF records into a chunked Arrow table.
"""
function ingest_vcf(
    input_path::String,
    output_path::String;
    chunk_size::Integer=10_000)
    chroms = String[]
    positions = Int32[]
    ids = String[]
    refs = String[]
    alts = String[]
    quals = Union{Missing,Float32}[]
    filters = String[]
    infos = String[]
    formats = String[]
    sample_counts = Int32[]
    samples = String[]

    _open_vcf_input(input_path, io -> begin
        open(Arrow.Writer, output_path; file=true) do writer
            for raw_line in eachline(io)
                startswith(raw_line, '#') && continue
                record = parse_vcf_record(raw_line)
                record === nothing && continue

                push!(chroms, record.chrom)
                push!(positions, record.pos)
                push!(ids, record.id)
                push!(refs, record.ref)
                push!(alts, record.alt)
                push!(quals, record.qual)
                push!(filters, record.filter)
                push!(infos, record.info)
                push!(formats, record.format)
                push!(sample_counts, Int32(length(record.samples)))
                push!(samples, join(record.samples, '\t'))

                if length(chroms) >= chunk_size
                    _write_vcf_chunk!(writer, chroms, positions, ids, refs, alts, quals, filters, infos, formats, sample_counts, samples)
                end
            end

            _write_vcf_chunk!(writer, chroms, positions, ids, refs, alts, quals, filters, infos, formats, sample_counts, samples)
        end
    end)

    return output_path
end

"""
    ingest_bed(input_path, output_path; chunk_size=10_000)

Convert BED records into a chunked Arrow table.
"""
function ingest_bed(
    input_path::String,
    output_path::String;
    chunk_size::Integer=10_000)
    chroms = String[]
    starts = Int32[]
    stops = Int32[]

    open(Arrow.Writer, output_path; file=true) do writer
        open(input_path, "r") do io
            for raw_line in eachline(io)
                startswith(raw_line, '#') && continue
                record = parse_bed_record(raw_line)
                record === nothing && continue

                push!(chroms, record.chrom)
                push!(starts, record.start)
                push!(stops, record.stop)

                if length(chroms) >= chunk_size
                    _write_bed_chunk!(writer, chroms, starts, stops)
                end
            end
        end

        _write_bed_chunk!(writer, chroms, starts, stops)
    end

    return output_path
end

"""
    _write_gff_chunk!(writer, chroms, sources, features, starts, stops, scores, strands, phases, attributes)

Write a chunk of flattened GFF records to an Arrow writer.
"""
function _write_gff_chunk!(writer, chroms, sources, features, starts, stops, scores, strands, phases, attributes)
    isempty(chroms) && return nothing

    Arrow.write(
        writer,
        (
            chrom = copy(chroms),
            source = copy(sources),
            feature = copy(features),
            start = copy(starts),
            stop = copy(stops),
            score = copy(scores),
            strand = copy(strands),
            phase = copy(phases),
            attributes = copy(attributes)))

    empty!(chroms)
    empty!(sources)
    empty!(features)
    empty!(starts)
    empty!(stops)
    empty!(scores)
    empty!(strands)
    empty!(phases)
    empty!(attributes)

    return nothing
end

"""
    ingest_gff(input_path, output_path; chunk_size=10_000)

Convert GFF records into a chunked Arrow table.
"""
function ingest_gff(
    input_path::String,
    output_path::String;
    chunk_size::Integer=10_000)
    chroms = String[]
    sources = String[]
    features = String[]
    starts = Int32[]
    stops = Int32[]
    scores = Union{Missing,Float32}[]
    strands = String[]
    phases = Union{Missing,Int8}[]
    attributes = String[]

    open(Arrow.Writer, output_path; file=true) do writer
        open(input_path, "r") do io
            for raw_line in eachline(io)
                startswith(raw_line, '#') && continue
                record = parse_gff_record(raw_line)
                record === nothing && continue

                push!(chroms, record.chrom)
                push!(sources, record.source)
                push!(features, record.feature)
                push!(starts, record.start)
                push!(stops, record.stop)
                push!(scores, record.score)
                push!(strands, record.strand)
                push!(phases, record.phase)
                push!(attributes, record.attributes)

                if length(chroms) >= chunk_size
                    _write_gff_chunk!(writer, chroms, sources, features, starts, stops, scores, strands, phases, attributes)
                end
            end
        end

        _write_gff_chunk!(writer, chroms, sources, features, starts, stops, scores, strands, phases, attributes)
    end

    return output_path
end

"""
    load_arrow_table(path)

Load an Arrow table from disk.
"""
function load_arrow_table(path::String)
    return Arrow.Table(path)
end

"""
    write_arrow_table(output_path, table)

Write a table to an Arrow file.
"""
function write_arrow_table(output_path::String, table)
    Arrow.write(output_path, table)
    return output_path
end

# 10. EMBL & SwissProt Support
# ------------------------------------------------------------------------------

"""
    read_embl(filepath::String)

Parses an EMBL flat file into a vector of AnnotatedSeqRecord objects.
"""
function read_embl(filepath::String; preserve_case::Bool=false)
    open(filepath, "r") do io
        return read_embl(io; preserve_case=preserve_case)
    end
end

"""
    read_embl(io; preserve_case=false)

Read EMBL records from an IO stream.
"""
function read_embl(io::IO; preserve_case::Bool=false)
    records = AnnotatedSeqRecord[]
    current_lines = String[]

    for raw_line in eachline(io)
        push!(current_lines, String(raw_line))
        if startswith(raw_line, "//")
            record = _parse_embl_record(current_lines; preserve_case=preserve_case)
            record === nothing || push!(records, record)
            empty!(current_lines)
        end
    end

    if !isempty(current_lines)
        record = _parse_embl_record(current_lines; preserve_case=preserve_case)
        record === nothing || push!(records, record)
    end
    return records
end

"""
    _embl_first_token(text)

Return the first whitespace-delimited token from an EMBL field.
"""
@inline function _embl_first_token(text::AbstractString)
    token_end = findfirst(isspace, text)
    token_end === nothing && return String(text)
    return String(text[firstindex(text):prevind(text, token_end)])
end

"""
    _embl_compact_tail(text)

Strip spaces from an EMBL continuation field.
"""
@inline function _embl_compact_tail(text::AbstractString)
    isempty(text) && return ""
    buffer = IOBuffer()
    for character in text
        character == ' ' && continue
        print(buffer, character)
    end
    return String(take!(buffer))
end

"""
    _parse_embl_record(lines)

Parse a single EMBL record from raw text lines.
"""
function _parse_embl_record(lines::AbstractVector{<:String}; preserve_case::Bool=false)
    id = ""
    accession = ""
    description = IOBuffer()
    sequence = IOBuffer()
    features = SeqFeatureLite[]
    annotations = Dict{Symbol, Any}()
    
    current_feature_type = ""
    current_feature_loc = ""
    current_feature_qualifiers = Dict{String, Vector{String}}()
    
    in_sequence = false
    
    for line in lines
        if length(line) < 2; continue; end
        code = line[1:2]
        content = length(line) > 5 ? strip(line[6:end]) : ""
        
        if code == "ID"
            id = _embl_first_token(content)
        elseif code == "AC"
            accession = isempty(accession) ? content : accession * " " * content
        elseif code == "DE"
            print(description, content, " ")
        elseif code == "FT"
            # Feature Table
            # EMBL format: 
            # FT   feature_name    location
            # FT                   /qualifier="value"
            # FT                   location_continuation
            
            # Check if this is a qualifier line
            line_content = line[6:end]
            trimmed_content = strip(line_content)
            
            if startswith(trimmed_content, "/")
                # Qualifier
                qual_part = trimmed_content[2:end] # remove /
                eq_index = findfirst(==( '=' ), qual_part)
                if eq_index !== nothing
                    q_key = String(qual_part[firstindex(qual_part):prevind(qual_part, eq_index)])
                    q_val = strip(String(qual_part[nextind(qual_part, eq_index):end]), '"')
                    push!(get!(current_feature_qualifiers, q_key, String[]), q_val)
                else
                    push!(get!(current_feature_qualifiers, String(qual_part), String[]), "")
                end
            elseif !isempty(trimmed_content) && isspace(line[6]) && isspace(line[7]) && isspace(line[8])
                current_feature_loc *= trimmed_content
            elseif !isempty(trimmed_content)
                space_index = findfirst(isspace, trimmed_content)
                if space_index !== nothing
                    if !isempty(current_feature_type)
                        loc_obj = parse_feature_location(current_feature_loc)
                        push!(features, SeqFeatureLite(current_feature_type, loc_obj, qualifiers=current_feature_qualifiers))
                    end
                    current_feature_type = String(trimmed_content[firstindex(trimmed_content):prevind(trimmed_content, space_index)])
                    current_feature_loc = _embl_compact_tail(trimmed_content[nextind(trimmed_content, space_index):end])
                    current_feature_qualifiers = Dict{String, Vector{String}}()
                end
            end
        elseif code == "SQ"
            in_sequence = true
        elseif in_sequence && code != "//"
            for b in codeunits(content)
                if (!preserve_case && UInt8('a') <= b <= UInt8('z'))
                    write(sequence, UInt8(b - 0x20))
                elseif (UInt8('a') <= b <= UInt8('z')) || (UInt8('A') <= b <= UInt8('Z'))
                    write(sequence, b)
                end
            end
        end
    end
    
    # Flush last feature
    if !isempty(current_feature_type)
        loc_obj = parse_feature_location(current_feature_loc)
        push!(features, SeqFeatureLite(current_feature_type, loc_obj, qualifiers=current_feature_qualifiers))
    end
    
    if !isempty(accession)
        annotations[:accession] = accession
    end

    sequence_bytes = take!(sequence)
    sequence_text = String(copy(sequence_bytes))
    sequence_typed = if isempty(sequence_text) || validate_sequence(DNAAlphabet, sequence_text)
        BioSequence{DNAAlphabet}(sequence_bytes; validate=false)
    elseif validate_sequence(RNAAlphabet, sequence_text)
        BioSequence{RNAAlphabet}(sequence_bytes; validate=false)
    elseif validate_sequence(AminoAcidAlphabet, sequence_text)
        BioSequence{AminoAcidAlphabet}(sequence_bytes; validate=false)
    else
        throw(ArgumentError("cannot infer alphabet for EMBL sequence"))
    end

    return AnnotatedSeqRecord(
        sequence_typed;
        identifier=id,
        name=id,
        description=String(strip(String(take!(description)))),
        annotations=annotations,
        features=features)
end

"""
    read_swissprot(filepath::String)

Parses a SwissProt/UniProt .dat file into a vector of AnnotatedSeqRecord objects.
"""
function read_swissprot(filepath::String)
    return read_embl(filepath) # Formats are extremely similar in structure
end

# --- ABIF (.ab1) Sanger Trace Parser ---

"""
    read_abif(filepath)

Read an ABIF chromatogram file from a path.
"""
function read_abif(filepath::String)
    open(filepath, "r") do io
        return read_abif(io)
    end
end

"""
    read_abif(io)

Read an ABIF chromatogram file from an IO stream.
"""
function read_abif(io::IO)
    magic = String(read(io, 4))
    if magic != "ABIF"
        throw(ArgumentError("Not a valid ABIF file (magic: \$magic)"))
    end
    
    version = ntoh(read(io, UInt16))
    
    # Directory entry (28 bytes) begins at offset 6
    dir_name = String(read(io, 4))
    dir_number = ntoh(read(io, UInt32))
    dir_type = ntoh(read(io, UInt16))
    dir_elem_size = ntoh(read(io, UInt16))
    dir_elements = ntoh(read(io, UInt32))
    dir_data_size = ntoh(read(io, UInt32))
    dir_data_offset = ntoh(read(io, UInt32))
    
    # Header is 128 bytes total, seek past it
    seek(io, 128)
    
    # Then seek to the actual directory location
    seek(io, dir_data_offset)
    sequence_bytes = UInt8[]
    qualities = UInt8[]
    fwo = "GATC"
    trace_data = Dict{Int, Vector{UInt16}}()

    @inline function read_tag_bytes(size::UInt32, data_offset::UInt32, data_pos::Int)
        if size <= 4
            seek(io, data_pos - 4)
            data = read(io, 4)[1:size]
        else
            seek(io, data_offset)
            data = read(io, size)
        end
        seek(io, data_pos)
        return data
    end
    
    for i in 1:dir_elements
        tag_name = String(read(io, 4))
        tag_num = ntoh(read(io, UInt32))
        tag_type = ntoh(read(io, UInt16))
        tag_elem_size = ntoh(read(io, UInt16))
        tag_num_elems = ntoh(read(io, UInt32))
        tag_data_size = ntoh(read(io, UInt32))
        tag_data_offset = ntoh(read(io, UInt32))
        data_handle = ntoh(read(io, UInt32))
        
        pos = position(io)

        if tag_name == "PBAS" && (tag_num == 1 || tag_num == 2)
            sequence_bytes = read_tag_bytes(tag_data_size, tag_data_offset, pos)
        elseif (tag_name == "PQC" || tag_name == "PCON") && (tag_num == 1 || tag_num == 2)
            qualities = read_tag_bytes(tag_data_size, tag_data_offset, pos)
        elseif tag_name == "FWO_" && tag_num == 1
            fwo = String(read_tag_bytes(tag_data_size, tag_data_offset, pos))
        elseif tag_name == "DATA"
            trace_data[Int(tag_num)] = [ntoh(val) for val in reinterpret(UInt16, read_tag_bytes(tag_data_size, tag_data_offset, pos))]
        else
            seek(io, pos)
        end
    end
    
    sequence = isempty(sequence_bytes) ? BioSequence{DNAAlphabet}(UInt8[]; validate=false) : BioSequence{DNAAlphabet}(String(sequence_bytes); validate=false)
    qualities = qualities
    
    trace_map = Dict{Char, Vector{UInt16}}()
    for (i, char) in enumerate(fwo)
        trace_map[char] = get(trace_data, i, UInt16[])
    end
    
    return SangerTrace(
        sequence,
        qualities,
        get(trace_map, 'A', UInt16[]),
        get(trace_map, 'C', UInt16[]),
        get(trace_map, 'G', UInt16[]),
        get(trace_map, 'T', UInt16[]),
        Dict{String, Any}("fwo" => fwo, "version" => version)
    )
end
