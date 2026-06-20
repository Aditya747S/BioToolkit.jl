# ==============================================================================
# bam.jl — BAM/SAM binary alignment format I/O
#
# Provides BAM reading/writing with BGZF compression, CIGAR encoding,
# sequence nibble-packing, quality encoding, auxiliary tag parsing, and
# BAI index support for random-access region queries.
#
# References:
#   - SAM/BAM specification v1.6 (hts-specs, samtools)
#   - Li et al. (2009) Bioinformatics 25(16):2078-2079 (SAM format)
#   - BGZF: blocked gzip with virtual file offsets
# ==============================================================================

const _BAM_MAGIC = UInt8[0x42, 0x41, 0x4d, 0x01]
const _BAM_BAI_MAGIC = UInt8[0x42, 0x41, 0x49, 0x01]
const _BAM_LINEAR_INDEX_WINDOW = 14

const _BAM_CIGAR_CODES = Dict(
    'M' => UInt32(0),
    'I' => UInt32(1),
    'D' => UInt32(2),
    'N' => UInt32(3),
    'S' => UInt32(4),
    'H' => UInt32(5),
    'P' => UInt32(6),
    '=' => UInt32(7),
    'X' => UInt32(8),
)

const _BAM_CIGAR_SYMBOLS = Dict(value => key for (key, value) in _BAM_CIGAR_CODES)

const _BAM_SEQUENCE_CODES = Dict(
    0x1 => 'A',
    0x2 => 'C',
    0x3 => 'M',
    0x4 => 'G',
    0x5 => 'R',
    0x6 => 'S',
    0x7 => 'V',
    0x8 => 'T',
    0x9 => 'W',
    0xA => 'Y',
    0xB => 'H',
    0xC => 'K',
    0xD => 'D',
    0xE => 'B',
    0xF => 'N',
)

const _BAM_SEQUENCE_ENCODE = Dict(
    'A' => 0x1,
    'C' => 0x2,
    'G' => 0x4,
    'T' => 0x8,
    'N' => 0xF,
)

"""
    BamReference

Reference sequence entry from a BAM header.
"""
struct BamReference
    name::String
    length::Int
end

"""
    BamHeader

Textual BAM header plus its reference list.
"""
struct BamHeader
    text::String
    references::Vector{BamReference}
end

"""
    BamCigarOp

Single CIGAR operation and its run length.
"""
struct BamCigarOp
    op::Char
    length::Int
end

"""
    BamRecord

In-memory BAM alignment record with sequence, qualities, and auxiliary tags.
"""
struct BamRecord
    qname::String
    flag::UInt16
    refname::Union{Nothing,String}
    pos::Int32
    mapq::UInt8
    cigar::Vector{BamCigarOp}
    mate_refname::Union{Nothing,String}
    mate_pos::Int32
    template_length::Int32
    sequence::BioSequence{DNAAlphabet}
    quality::Union{Missing,String}
    tags::Dict{String,Any}
end

"""
    BamFile

In-memory BAM file container holding a header and records.
"""
struct BamFile
    header::BamHeader
    records::Vector{BamRecord}
    metadata::Dict{Symbol,Any}
end

"""
    BamChunk

Compressed file chunk used for BAM random-access indexing.
"""
struct BamChunk
    start::UInt64
    stop::UInt64
end

"""
    BamIndex

Bin and linear index tables used for BAM region queries.
"""
struct BamIndex
    bins::Vector{Dict{UInt32,Vector{BamChunk}}}
    linear::Vector{Dict{Int,UInt64}}
end

"""
    AbstractBamRecordReader

Abstract supertype for lazy BAM record iterators.
"""
abstract type AbstractBamRecordReader end

"""
    BamReader

Sequential lazy BAM record reader.
"""
mutable struct BamReader <: AbstractBamRecordReader
    io::BGZFStreams.BGZFStream
    header::BamHeader
    done::Bool
end

"""
    BamRegionScanReader

Sequential lazy reader that yields only records overlapping a region.
"""
mutable struct BamRegionScanReader <: AbstractBamRecordReader
    io::BGZFStreams.BGZFStream
    header::BamHeader
    region::GenomicRanges.GenomicInterval
    done::Bool
end

"""
    BamIndexedRegionReader

Indexed lazy reader that jumps to candidate chunks for a region.
"""
mutable struct BamIndexedRegionReader <: AbstractBamRecordReader
    io::BGZFStreams.BGZFStream
    header::BamHeader
    region::GenomicRanges.GenomicInterval
    chunks::Vector{BamChunk}
    chunk_index::Int
    chunk_stop::UInt64
    seen::Set{UInt64}
    done::Bool
end

"""
    BamEmptyReader

Empty lazy BAM record reader with a known header.
"""
struct BamEmptyReader <: AbstractBamRecordReader
    header::BamHeader
end

Base.IteratorEltype(::Type{<:AbstractBamRecordReader}) = Base.HasEltype()
Base.IteratorSize(::Type{<:AbstractBamRecordReader}) = Base.SizeUnknown()
Base.eltype(::Type{<:AbstractBamRecordReader}) = BamRecord

_close_bam_io!(io::BGZFStreams.BGZFStream) = (
    try
        close(io)
    catch
    end
)

@inline function _register_bam_result!(_explicit_ctx, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    _ctx = active_provenance_context(_explicit_ctx)
    _ctx === nothing || register_provenance!(_ctx, operation; parents=parents, parameters=parameters)
    return nothing
end

"""
    Base.:(==)(left, right)

Test two BAM references for value equality.
"""
function Base.:(==)(left::BamReference, right::BamReference)
    return isequal(left.name, right.name) && isequal(left.length, right.length)
end

"""
    Base.:(==)(left, right)

Test two BAM CIGAR operations for value equality.
"""
function Base.:(==)(left::BamCigarOp, right::BamCigarOp)
    return isequal(left.op, right.op) && isequal(left.length, right.length)
end

"""
    Base.:(==)(left, right)

Test two BAM records for value equality.
"""
function Base.:(==)(left::BamRecord, right::BamRecord)
    return isequal(left.qname, right.qname) && isequal(left.flag, right.flag) && isequal(left.refname, right.refname) && isequal(left.pos, right.pos) && isequal(left.mapq, right.mapq) && isequal(left.cigar, right.cigar) && isequal(left.mate_refname, right.mate_refname) && isequal(left.mate_pos, right.mate_pos) && isequal(left.template_length, right.template_length) && isequal(left.sequence, right.sequence) && isequal(left.quality, right.quality) && isequal(left.tags, right.tags)
end

"""
    Base.:(==)(left, right)

Test two BAM headers for value equality.
"""
function Base.:(==)(left::BamHeader, right::BamHeader)
    return isequal(left.text, right.text) && isequal(left.references, right.references)
end

Base.iterate(file::BamFile, state...) = iterate(file.records, state...)
Base.length(file::BamFile) = length(file.records)
Base.getindex(file::BamFile, index::Integer) = file.records[index]

"""
    Base.show(io, reference)

Render a compact BAM reference summary.
"""
function Base.show(io::IO, reference::BamReference)
    print(io, "BamReference(", reference.name, ", length=", reference.length, ")")
end

"""
    Base.show(io, op)

Render a compact BAM CIGAR summary.
"""
function Base.show(io::IO, op::BamCigarOp)
    print(io, "BamCigarOp(", op.length, op.op, ")")
end

"""
    Base.show(io, record)

Render a compact BAM record summary.
"""
function Base.show(io::IO, record::BamRecord)
    chrom = record.refname === nothing ? "*" : record.refname
    print(io, "BamRecord(", record.qname, ", ", chrom, ":", record.pos + 1, ", ", length(record.sequence), " bp)")
end

"""
    Base.show(io, file)

Render a compact BAM file summary.
"""
function Base.show(io::IO, file::BamFile)
    print(io, "BamFile(", length(file.records), " records, ", length(file.header.references), " references, ", container_provenance_summary(file), ")")
end

Base.:(==)(left::BamFile, right::BamFile) = isequal(left.header, right.header) && isequal(left.records, right.records) && isequal(left.metadata, right.metadata)

"""
    _bam_header_text(references)

Generate a minimal SAM-style header text block from BAM references.
"""
function _bam_header_text(references::AbstractVector{<:BamReference})
    buffer = IOBuffer()
    println(buffer, "@HD\tVN:1.6\tSO:unknown")
    for reference in references
        println(buffer, "@SQ\tSN:", reference.name, "\tLN:", reference.length)
    end
    return String(take!(buffer))
end

"""
    BamReference(name, length)

Construct a normalized BAM reference entry.
"""
function BamReference(name::String, length::Integer)
    return BamReference(String(name), Int(length))
end

"""
    BamHeader(references; text="")

Construct a BAM header from a reference list and optional header text.
"""
function BamHeader(references::AbstractVector{<:BamReference}; text::String="")
    reference_list = BamReference[BamReference(reference.name, reference.length) for reference in references]
    header_text = isempty(strip(text)) ? _bam_header_text(reference_list) : String(text)
    return BamHeader(header_text, reference_list)
end

"""
    BamHeader(text, references)

Construct a BAM header from explicit text and references.
"""
function BamHeader(text::String, references::AbstractVector{<:BamReference})
    return BamHeader(String(text), BamReference[BamReference(reference.name, reference.length) for reference in references])
end

"""
    BamCigarOp(length, op)

Construct a BAM CIGAR operation from its length and symbol.
"""
function BamCigarOp(length::Integer, op::Char)
    return BamCigarOp(op, Int(length))
end

"""
    BamRecord(...)

Construct a fully typed BAM alignment record.
"""
function BamRecord(
    qname::String,
    refname::Union{Nothing,String},
    pos::Integer,
    cigar::AbstractVector{<:BamCigarOp},
    sequence::BioSequence{DNAAlphabet};
    flag::Integer=0,
    mapq::Integer=60,
    mate_refname::Union{Nothing,String}=nothing,
    mate_pos::Integer=-1,
    template_length::Integer=0,
    quality::Union{Missing,String}=missing,
    tags::AbstractDict=Dict{String,Any}(),
)
    ref_value = refname === nothing ? nothing : String(refname)
    mate_ref_value = mate_refname === nothing ? nothing : String(mate_refname)
    cigar_ops = BamCigarOp[BamCigarOp(op.length, op.op) for op in cigar]
    quality_value = quality === missing ? missing : String(quality)
    return BamRecord(
        String(qname),
        UInt16(flag),
        ref_value,
        Int32(pos),
        UInt8(mapq),
        cigar_ops,
        mate_ref_value,
        Int32(mate_pos),
        Int32(template_length),
        sequence,
        quality_value,
        Dict{String,Any}(tags),
    )
end

function BamRecord(
    qname::String,
    refname::Union{Nothing,String},
    pos::Integer,
    cigar::AbstractVector{<:BamCigarOp},
    sequence::BioSequence;
    kwargs...
)
    return BamRecord(qname, refname, pos, cigar, BioSequence{DNAAlphabet}(sequence.data; validate=false); kwargs...)
end

@inline _bam_coerce_dna_sequence(sequence::BioSequence{DNAAlphabet}) = sequence

function _bam_coerce_dna_sequence(sequence::BioSequence)
    sequence_text = String(sequence)
    validate_sequence(DNAAlphabet, sequence_text) || throw(ArgumentError("BamRecord sequence must be DNA-compatible"))
    return BioSequence{DNAAlphabet}(sequence_text)
end

function _bam_coerce_dna_sequence(sequence)
    sequence_text = String(sequence)
    validate_sequence(DNAAlphabet, sequence_text) || throw(ArgumentError("BamRecord sequence must be DNA-compatible"))
    return BioSequence{DNAAlphabet}(sequence_text)
end

function BamRecord(
    qname::String,
    refname::Union{Nothing,String},
    pos::Integer,
    cigar::AbstractVector{<:BamCigarOp},
    sequence;
    kwargs...
)
    return BamRecord(qname, refname, pos, cigar, _bam_coerce_dna_sequence(sequence); kwargs...)
end

"""
    BamFile(records; header=nothing)

Wrap a list of BAM records in an in-memory file container.
"""
function BamFile(records::AbstractVector{<:BamRecord}; header::Union{Nothing,BamHeader}=nothing, metadata::AbstractDict=Dict{Symbol,Any}())
    bam_records = BamRecord[record for record in records]
    bam_header = header === nothing ? _infer_bam_header(bam_records) : header
    metadata_copy = Dict{Symbol,Any}(metadata)
    ensure_provenance_id!(metadata_copy)
    return BamFile(bam_header, bam_records, metadata_copy)
end

function BamFile(header::BamHeader, records::AbstractVector{<:BamRecord}; metadata::AbstractDict=Dict{Symbol,Any}())
    bam_records = BamRecord[record for record in records]
    metadata_copy = Dict{Symbol,Any}(metadata)
    ensure_provenance_id!(metadata_copy)
    return BamFile(header, bam_records, metadata_copy)
end

"""
    _read_cstring(io, nbytes)

Read a NUL-terminated string from a BAM binary stream.
"""
function _read_cstring(io::IO, nbytes::Integer)
    nbytes > 0 || return ""
    bytes = read(io, Int(nbytes))
    isempty(bytes) && return ""
    if bytes[end] == 0x00
        pop!(bytes)
    end
    return String(bytes)
end

"""
    _read_bam_name(io, nbytes)

Read a BAM name field from a binary stream.
"""
function _read_bam_name(io::IO, nbytes::Integer)
    bytes = read(io, Int(nbytes))
    isempty(bytes) && return ""
    if bytes[end] == 0x00
        pop!(bytes)
    end
    return String(bytes)
end

"""
    _decode_bam_sequence(bytes, length)

Decode BAM-packed nucleotide bytes into a typed DNA sequence.
"""
function _decode_bam_sequence(bytes::Vector{UInt8}, length::Integer)
    if length == 0
        return BioSequence{DNAAlphabet}(UInt8[]; validate=false)
    end

    buffer = Vector{UInt8}(undef, Int(length))
    expected = Int(length)
    written = 0
    for byte in bytes
        high = (byte >> 4) & 0x0f
        low = byte & 0x0f
        if written < expected
            written += 1
            buffer[written] = UInt8(get(_BAM_SEQUENCE_CODES, high, 'N'))
        end
        if written < expected
            written += 1
            buffer[written] = UInt8(get(_BAM_SEQUENCE_CODES, low, 'N'))
        end
    end

    return BioSequence{DNAAlphabet}(buffer; validate=false)
end

"""
    _encode_bam_sequence(sequence)

Encode a nucleotide sequence into BAM nibble-packed bytes.
"""
function _encode_bam_sequence(sequence::BioSequence{DNAAlphabet})
    encoded = UInt8[]
    buffer = UInt8(0)
    use_high = true

    for character in String(sequence)
        code = get(_BAM_SEQUENCE_ENCODE, character, 0xF)
        if use_high
            buffer = UInt8(code << 4)
            use_high = false
        else
            push!(encoded, buffer | UInt8(code))
            use_high = true
        end
    end

    if !use_high
        push!(encoded, buffer)
    end

    return encoded
end

"""
    _decode_bam_quality(bytes)

Decode BAM quality bytes into an ASCII quality string.
"""
function _decode_bam_quality(bytes::Vector{UInt8})
    isempty(bytes) && return missing
    all(byte -> byte == 0xff, bytes) && return missing
    buffer = IOBuffer()
    for byte in bytes
        write(buffer, Char(byte + 33))
    end
    return String(take!(buffer))
end

"""
    _encode_bam_quality(quality, length)

Encode an ASCII quality string for BAM storage.
"""
function _encode_bam_quality(quality::Union{Missing,String}, length::Integer)
    quality === missing && return fill(UInt8(0xff), Int(length))
    ncodeunits(quality) == length || throw(ArgumentError("quality length must match sequence length"))
    return UInt8.(codeunits(String(quality)) .- UInt8(33))
end

"""
    _decode_bam_cigar(value)

Decode a packed BAM CIGAR value into a `BamCigarOp`.
"""
function _decode_bam_cigar(value::UInt32)
    op = get(_BAM_CIGAR_SYMBOLS, value & 0x0f, 'M')
    length = Int(value >> 4)
    return BamCigarOp(op, length)
end

"""
    _encode_bam_cigar(op)

Encode a `BamCigarOp` into a packed BAM CIGAR value.
"""
function _encode_bam_cigar(op::BamCigarOp)
    code = get(_BAM_CIGAR_CODES, op.op, nothing)
    code === nothing && throw(ArgumentError("unsupported BAM CIGAR op $(op.op)"))
    return UInt32(op.length) << 4 | code
end

"""
    _read_bam_tag_value(io, type)

Read a single BAM auxiliary tag value from a binary stream.
"""
function _read_bam_array_values(io::IO, subtype::Char, count::Int)
    count >= 0 || throw(ArgumentError("BAM B-array count must be non-negative"))

    if subtype == 'c'
        values = Vector{Int8}(undef, count)
        read!(io, values)
        return values
    elseif subtype == 'C'
        values = Vector{UInt8}(undef, count)
        read!(io, values)
        return values
    elseif subtype == 's'
        values = Vector{Int16}(undef, count)
        read!(io, values)
        return values
    elseif subtype == 'S'
        values = Vector{UInt16}(undef, count)
        read!(io, values)
        return values
    elseif subtype == 'i'
        values = Vector{Int32}(undef, count)
        read!(io, values)
        return values
    elseif subtype == 'I'
        values = Vector{UInt32}(undef, count)
        read!(io, values)
        return values
    elseif subtype == 'f'
        values = Vector{Float32}(undef, count)
        read!(io, values)
        return values
    elseif subtype == 'd'
        values = Vector{Float64}(undef, count)
        read!(io, values)
        return values
    end

    throw(ArgumentError("unsupported BAM B-array subtype '$subtype'"))
end

function _read_bam_tag_value(io::IO, type::Char)
    if type == 'A'
        return Char(read(io, UInt8))
    elseif type == 'c'
        return read(io, Int8)
    elseif type == 'C'
        return read(io, UInt8)
    elseif type == 's'
        return read(io, Int16)
    elseif type == 'S'
        return read(io, UInt16)
    elseif type == 'i'
        return read(io, Int32)
    elseif type == 'I'
        return read(io, UInt32)
    elseif type == 'f'
        return read(io, Float32)
    elseif type == 'd'
        return read(io, Float64)
    elseif type == 'Z' || type == 'H'
        buffer = IOBuffer()
        while !eof(io)
            byte = read(io, UInt8)
            byte == 0x00 && break
            write(buffer, byte)
        end
        return String(take!(buffer))
    elseif type == 'B'
        subtype = Char(read(io, UInt8))
        count = Int(read(io, Int32))
        return _read_bam_array_values(io, subtype, count)
    end

    throw(ArgumentError("unsupported BAM tag type '$type'"))
end

"""
    _parse_bam_tags(bytes)

Parse packed BAM auxiliary tags into a dictionary.
"""
function _parse_bam_tags(bytes::Vector{UInt8})
    tags = Dict{String,Any}()
    io = IOBuffer(bytes)

    while !eof(io)
        remaining = bytesavailable(io)
        remaining == 0 && break
        remaining < 3 && throw(ArgumentError("truncated BAM auxiliary tag payload"))

        key = String(read(io, 2))
        type_char = Char(read(io, UInt8))
        tags[key] = _read_bam_tag_value(io, type_char)
    end

    return tags
end

"""
    _encode_bam_tag_value(value)

Encode a BAM auxiliary tag value into binary form.
"""
function _encode_bam_tag_value(value)
    if value isa Char
        return 'A', UInt8(value)
    elseif value isa Int8
        return 'c', value
    elseif value isa UInt8
        return 'C', value
    elseif value isa Int16
        return 's', value
    elseif value isa UInt16
        return 'S', value
    elseif value isa Int32 || value isa Int || value isa Int64
        return 'i', Int32(value)
    elseif value isa UInt32 || value isa UInt || value isa UInt64
        return 'I', UInt32(value)
    elseif value isa Float64
        return 'd', value
    elseif value isa AbstractFloat
        return 'f', Float32(value)
    elseif value isa String
        return 'Z', String(value)
    elseif value isa AbstractVector
        if isempty(value)
            return 'B', ('i', Int32[])
        elseif eltype(value) <: Integer
            if eltype(value) <: Int8
                return 'B', ('c', Int8.(value))
            elseif eltype(value) <: UInt8
                return 'B', ('C', UInt8.(value))
            elseif eltype(value) <: Int16
                return 'B', ('s', Int16.(value))
            elseif eltype(value) <: UInt16
                return 'B', ('S', UInt16.(value))
            elseif eltype(value) <: UInt32 || eltype(value) <: UInt || eltype(value) <: UInt64
                return 'B', ('I', UInt32.(value))
            else
                return 'B', ('i', Int32.(value))
            end
        elseif eltype(value) <: AbstractFloat
            if eltype(value) <: Float64
                return 'B', ('d', Float64.(value))
            end
            return 'B', ('f', Float32.(value))
        end
    end

    throw(ArgumentError("unsupported BAM tag value $(typeof(value))"))
end

"""
    _write_bam_tags(io, tags)

Write BAM auxiliary tags to a binary stream.
"""
function _write_bam_tags(io::IO, tags::Dict{String,Any})
    for (tag, value) in tags
        length(tag) == 2 || throw(ArgumentError("BAM tag keys must be two characters"))
        type, encoded = _encode_bam_tag_value(value)
        write(io, codeunits(tag))
        write(io, UInt8(type))
        if type == 'A'
            write(io, encoded)
        elseif type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I' || type == 'f' || type == 'd'
            write(io, encoded)
        elseif type == 'Z' || type == 'H'
            write(io, codeunits(encoded))
            write(io, UInt8(0x00))
        elseif type == 'B'
            subtype, values = encoded
            write(io, UInt8(subtype))
            write(io, Int32(length(values)))
            for entry in values
                write(io, entry)
            end
        end
    end
    return nothing
end

"""
    _bam_reference_name_map(header)

Build a mapping from BAM reference names to reference indices.
"""
function _bam_reference_name_map(header::BamHeader)
    return Dict(reference.name => index for (index, reference) in pairs(header.references))
end

"""
    _reference_index(header, refname)

Resolve a reference name to its BAM index.
"""
function _reference_index(header::BamHeader, refname::Union{Nothing,String})
    refname === nothing && return 0
    for (index, reference) in pairs(header.references)
        reference.name == refname && return index
    end
    return 0
end

"""
    _reference_name(header, refindex)

Resolve a BAM reference index back to its name.
"""
function _reference_name(header::BamHeader, refindex::Integer)
    refindex < 0 && return nothing
    index = refindex + 1
    index > length(header.references) && return nothing
    return header.references[index].name
end

"""
    _bam_reference_span(cigar)

Compute the reference span covered by a BAM CIGAR string.
"""
function _bam_reference_span(cigar::AbstractVector{<:BamCigarOp})
    span = 0
    for op in cigar
        op.op in ('M', 'D', 'N', '=', 'X') && (span += op.length)
    end
    return span
end

"""
    _bam_query_span(record)

Compute the query span covered by a BAM record.
"""
function _bam_query_span(record::BamRecord)
    record.pos < 0 && return 0, -1
    span = max(1, _bam_reference_span(record.cigar))
    return Int(record.pos), Int(record.pos) + span
end

"""
    _bam_overlaps(record, region)

Test whether a BAM record overlaps a genomic interval.
"""
function _bam_overlaps(record::BamRecord, region::GenomicRanges.GenomicInterval)
    record.refname === nothing && return false
    record.refname == region.chrom || return false
    beg, end_ = _bam_query_span(record)
    end_ <= beg && return false
    region_left = region.left - 1
    region_right = region.right
    return beg < region_right && end_ > region_left
end

"""
    _bam_record_from_stream(io, header)

Decode one BAM record from a BGZF stream.
"""
function _bam_record_from_stream(io::BGZFStreams.BGZFStream, header::BamHeader)
    record_start = convert(UInt64, BGZFStreams.virtualoffset(io))
    try
        block_size = Int(read(io, Int32))
        block_size >= 32 || throw(ArgumentError("invalid BAM record block size: $block_size"))
        refid = Int(read(io, Int32))
        pos = read(io, Int32)
        l_read_name = Int(read(io, UInt8))
        l_read_name >= 0 || throw(ArgumentError("invalid BAM read-name length: $l_read_name"))
        mapq = read(io, UInt8)
        _bin = read(io, UInt16)
        n_cigar = Int(read(io, UInt16))
        flag = read(io, UInt16)
        l_seq = Int(read(io, Int32))
        l_seq >= 0 || throw(ArgumentError("invalid BAM sequence length: $l_seq"))
        next_refid = Int(read(io, Int32))
        next_pos = read(io, Int32)
        tlen = read(io, Int32)

        qname = l_read_name == 0 ? "" : _read_bam_name(io, l_read_name)
        cigar = BamCigarOp[_decode_bam_cigar(read(io, UInt32)) for _ in 1:n_cigar]

        seq_bytes = read(io, cld(l_seq, 2))
        sequence = _decode_bam_sequence(seq_bytes, l_seq)

        quality_bytes = read(io, l_seq)
        quality = _decode_bam_quality(quality_bytes)

        consumed = 32 + l_read_name + 4 * n_cigar + cld(l_seq, 2) + l_seq
        aux_bytes = block_size - consumed
        aux_bytes >= 0 || throw(ArgumentError("malformed BAM record: block size ($block_size) is smaller than consumed bytes ($consumed)"))
        tags = aux_bytes == 0 ? Dict{String,Any}() : _parse_bam_tags(read(io, aux_bytes))

        refname = _reference_name(header, refid)
        mate_refname = _reference_name(header, next_refid)
        return record_start, BamRecord(qname, flag, refname, pos, mapq, cigar, mate_refname, next_pos, tlen, sequence, quality, tags)
    catch err
        if err isa EOFError
            return nothing, nothing
        end
        rethrow(err)
    end
end

"""
    _bam_reg2bins(beg, end_)

Compute the bin list covering a BAM genomic region.
"""
function _bam_reg2bins(beg::Integer, end_::Integer)
    beg_int = Int(beg)
    end_int = Int(end_) - 1
    end_int < beg_int && return UInt32[0]

    bins = UInt32[0]
    for (shift, offset) in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681))
        start_bin = offset + (beg_int >> shift)
        end_bin = offset + (end_int >> shift)
        for bin in start_bin:end_bin
            push!(bins, UInt32(bin))
        end
    end

    return bins
end

"""
    _bam_region_chunks(index, refindex, left, right)

Return candidate BAM chunks for a genomic region.
"""
function _bam_linear_floor_offset(index::BamIndex, refindex::Integer, left::Integer)
    refindex <= 0 && return UInt64(0)
    refindex > length(index.linear) && return UInt64(0)

    ref_linear = index.linear[refindex]
    isempty(ref_linear) && return UInt64(0)
    window = max((left - 1) >> _BAM_LINEAR_INDEX_WINDOW, 0)

    while window >= 0
        offset = get(ref_linear, window, UInt64(0))
        offset != 0 && return offset
        window -= 1
    end

    return UInt64(0)
end

function _bam_merge_chunks(chunks::Vector{BamChunk})
    isempty(chunks) && return BamChunk[]

    merged = BamChunk[]
    current_start = chunks[1].start
    current_stop = chunks[1].stop

    for chunk in @view chunks[2:end]
        if chunk.start <= current_stop
            current_stop = max(current_stop, chunk.stop)
        else
            push!(merged, BamChunk(current_start, current_stop))
            current_start = chunk.start
            current_stop = chunk.stop
        end
    end
    push!(merged, BamChunk(current_start, current_stop))

    return merged
end

function _bam_region_chunks(index::BamIndex, refindex::Integer, left::Integer, right::Integer)
    refindex <= 0 && return BamChunk[]
    refindex > length(index.bins) && return BamChunk[]
    region_beg = max(left - 1, 0)
    region_end = max(right, region_beg + 1)
    chunks = BamChunk[]
    seen = Set{Tuple{UInt64,UInt64}}()
    ref_bins = index.bins[refindex]

    for bin in _bam_reg2bins(region_beg, region_end)
        haskey(ref_bins, bin) || continue
        for chunk in ref_bins[bin]
            key = (chunk.start, chunk.stop)
            key in seen && continue
            push!(seen, key)
            push!(chunks, chunk)
        end
    end

    sort!(chunks, by=chunk -> (chunk.start, chunk.stop))
    min_offset = _bam_linear_floor_offset(index, refindex, left)
    if min_offset != 0
        filtered = BamChunk[]
        for chunk in chunks
            chunk.stop <= min_offset && continue
            push!(filtered, BamChunk(max(chunk.start, min_offset), chunk.stop))
        end
        chunks = filtered
    end

    return _bam_merge_chunks(chunks)
end

"""
    _bam_infer_reference_lengths(records)

Infer reference lengths from a set of BAM records.
"""
function _bam_infer_reference_lengths(records::AbstractVector{<:BamRecord})
    names = String[]
    lengths = Dict{String,Int}()

    for record in records
        record.refname === nothing && continue
        if !(record.refname in names)
            push!(names, record.refname)
        end
        beg, end_ = _bam_query_span(record)
        end_ > beg || continue
        current = get(lengths, record.refname, 0)
        lengths[record.refname] = max(current, end_)
    end

    return BamReference[BamReference(name, max(get(lengths, name, 0), 1)) for name in names]
end

"""
    _infer_bam_header(records)

Infer a BAM header from the records when no header is provided.
"""
function _infer_bam_header(records::AbstractVector{<:BamRecord})
    references = _bam_infer_reference_lengths(records)
    return BamHeader(references)
end

"""
    _write_bam_record(io, header, record, index=nothing)

Write a single BAM record to a BGZF stream.
"""
function _write_bam_record(io::BGZFStreams.BGZFStream, header::BamHeader, record::BamRecord, index::Union{Nothing,BamIndex}=nothing)
    record_start = convert(UInt64, BGZFStreams.virtualoffset(io))
    refindex = _reference_index(header, record.refname)
    mate_refindex = _reference_index(header, record.mate_refname)
    cigar_bytes = UInt32[_encode_bam_cigar(op) for op in record.cigar]
    sequence_bytes = _encode_bam_sequence(record.sequence)
    quality_bytes = _encode_bam_quality(record.quality, length(record.sequence))
    tag_io = IOBuffer()
    _write_bam_tags(tag_io, record.tags)
    tags_bytes = take!(tag_io)

    block_size = Int32(32 + ncodeunits(record.qname) + 1 + 4 * length(cigar_bytes) + length(sequence_bytes) + length(quality_bytes) + length(tags_bytes))
    write(io, block_size)
    write(io, Int32(refindex == 0 ? -1 : refindex - 1))
    write(io, record.pos)
    write(io, UInt8(ncodeunits(record.qname) + 1))
    write(io, record.mapq)
    write(io, UInt16(0))
    write(io, UInt16(length(cigar_bytes)))
    write(io, record.flag)
    write(io, Int32(length(record.sequence)))
    write(io, Int32(mate_refindex == 0 ? -1 : mate_refindex - 1))
    write(io, record.mate_pos)
    write(io, record.template_length)
    write(io, codeunits(record.qname))
    write(io, UInt8(0x00))
    for value in cigar_bytes
        write(io, value)
    end
    write(io, sequence_bytes)
    write(io, quality_bytes)
    write(io, tags_bytes)

    record_end = convert(UInt64, BGZFStreams.virtualoffset(io))

    if index !== nothing && refindex > 0
        beg, end_ = _bam_query_span(record)
        if end_ > beg
            ref_bins = index.bins[refindex]
            ref_linear = index.linear[refindex]
            for bin in _bam_reg2bins(beg, end_)
                push!(get!(ref_bins, bin, BamChunk[]), BamChunk(record_start, record_end))
            end
            window_start = beg >> _BAM_LINEAR_INDEX_WINDOW
            window_stop = (end_ - 1) >> _BAM_LINEAR_INDEX_WINDOW
            for window in window_start:window_stop
                current = get(ref_linear, window, typemax(UInt64))
                current > record_start && (ref_linear[window] = record_start)
            end
        end
    end

    return nothing
end

"""
    _write_bam_index(path, header, index)

Write a BAM index file to disk.
"""
function _write_bam_index(path::String, header::BamHeader, index::BamIndex)
    open(path, "w") do io
        write(io, _BAM_BAI_MAGIC)
        write(io, Int32(length(header.references)))
        for refindex in eachindex(header.references)
            ref_bins = index.bins[refindex]
            write(io, Int32(length(ref_bins)))
            for (bin_id, chunks) in sort(collect(ref_bins); by=first)
                write(io, UInt32(bin_id))
                write(io, Int32(length(chunks)))
                for chunk in chunks
                    write(io, chunk.start)
                    write(io, chunk.stop)
                end
            end

            ref_linear = index.linear[refindex]
            if isempty(ref_linear)
                write(io, Int32(0))
            else
                max_window = maximum(keys(ref_linear))
                write(io, Int32(max_window + 1))
                for window in 0:max_window
                    write(io, get(ref_linear, window, UInt64(0)))
                end
            end
        end
    end
    return path
end

"""
    _read_bam_index(path)

Read a BAM index file from disk.
"""
function _read_bam_index(path::String)
    open(path, "r") do io
        magic = read(io, 4)
        magic == _BAM_BAI_MAGIC || throw(ArgumentError("not a BAI index file"))
        n_ref = Int(read(io, Int32))
        bins = Vector{Dict{UInt32,Vector{BamChunk}}}(undef, n_ref)
        linear = Vector{Dict{Int,UInt64}}(undef, n_ref)

        for refindex in 1:n_ref
            n_bin = Int(read(io, Int32))
            ref_bins = Dict{UInt32,Vector{BamChunk}}()
            for _ in 1:n_bin
                bin_id = read(io, UInt32)
                n_chunk = Int(read(io, Int32))
                chunks = BamChunk[]
                for _ in 1:n_chunk
                    start = read(io, UInt64)
                    stop = read(io, UInt64)
                    push!(chunks, BamChunk(start, stop))
                end
                ref_bins[bin_id] = chunks
            end
            n_intv = Int(read(io, Int32))
            ref_linear = Dict{Int,UInt64}()
            for window in 0:max(n_intv - 1, -1)
                offset = read(io, UInt64)
                offset != 0 && (ref_linear[window] = offset)
            end
            bins[refindex] = ref_bins
            linear[refindex] = ref_linear
        end

        return BamIndex(bins, linear)
    end
end

"""
    _read_bam_header!(io)

Read a BAM header from an already-open BGZF stream.
"""
function _read_bam_header!(io::BGZFStreams.BGZFStream)
    magic = read(io, 4)
    magic == _BAM_MAGIC || throw(ArgumentError("not a BAM file"))
    text_length = Int(read(io, Int32))
    text = text_length > 0 ? String(read(io, text_length)) : ""
    n_ref = Int(read(io, Int32))
    references = BamReference[]
    for _ in 1:n_ref
        name_length = Int(read(io, Int32))
        name = _read_bam_name(io, name_length)
        length = Int(read(io, Int32))
        push!(references, BamReference(name, length))
    end
    return BamHeader(text, references)
end

function _assert_pure_julia_alignment_path(path::String)
    endswith(lowercase(path), ".cram") && throw(ArgumentError("CRAM input is not supported yet in pure-Julia mode"))
    return nothing
end

function _open_bam_stream(path::String)
    _assert_pure_julia_alignment_path(path)
    io = open(BGZFStreams.BGZFStream, path, "r")
    try
        header = _read_bam_header!(io)
        return io, header
    catch
        _close_bam_io!(io)
        rethrow()
    end
end

"""
    _read_bam_header(path)

Read a BAM header from disk.
"""
function _read_bam_header(path::String)
    _assert_pure_julia_alignment_path(path)
    open(BGZFStreams.BGZFStream, path, "r") do io
        return _read_bam_header!(io)
    end
end

function BamReader(path::String)
    io, header = _open_bam_stream(path)
    return BamReader(io, header, false)
end

function BamRegionScanReader(path::String, region::GenomicRanges.GenomicInterval)
    io, header = _open_bam_stream(path)
    return BamRegionScanReader(io, header, region, false)
end

function Base.close(reader::BamReader)
    reader.done = true
    _close_bam_io!(reader.io)
    return nothing
end

function Base.close(reader::BamRegionScanReader)
    reader.done = true
    _close_bam_io!(reader.io)
    return nothing
end

function Base.close(reader::BamIndexedRegionReader)
    reader.done = true
    _close_bam_io!(reader.io)
    return nothing
end

Base.close(::BamEmptyReader) = nothing

function Base.iterate(reader::BamReader, _state=nothing)
    reader.done && return nothing

    _record_start, record = _bam_record_from_stream(reader.io, reader.header)
    if record === nothing
        close(reader)
        return nothing
    end

    return record, nothing
end

function Base.iterate(reader::BamRegionScanReader, _state=nothing)
    reader.done && return nothing

    while true
        _record_start, record = _bam_record_from_stream(reader.io, reader.header)
        if record === nothing
            close(reader)
            return nothing
        end
        _bam_overlaps(record, reader.region) && return record, nothing
    end
end

function _bam_indexed_seek_next_chunk!(reader::BamIndexedRegionReader)
    while reader.chunk_index <= length(reader.chunks)
        chunk = reader.chunks[reader.chunk_index]
        reader.chunk_index += 1
        chunk.stop > chunk.start || continue
        seek(reader.io, convert(BGZFStreams.VirtualOffset, chunk.start))
        reader.chunk_stop = chunk.stop
        return true
    end

    reader.chunk_stop = UInt64(0)
    return false
end

function Base.iterate(reader::BamIndexedRegionReader, _state=nothing)
    reader.done && return nothing

    while true
        if reader.chunk_stop == 0 || convert(UInt64, BGZFStreams.virtualoffset(reader.io)) >= reader.chunk_stop
            _bam_indexed_seek_next_chunk!(reader) || begin
                close(reader)
                return nothing
            end
        end

        record_start, record = _bam_record_from_stream(reader.io, reader.header)
        if record === nothing
            close(reader)
            return nothing
        end
        record_start >= reader.chunk_stop && continue
        record_start in reader.seen && continue
        push!(reader.seen, record_start)
        _bam_overlaps(record, reader.region) || continue
        return record, nothing
    end
end

Base.iterate(::BamEmptyReader, _state=nothing) = nothing

function _bam_collect(reader::AbstractBamRecordReader)
    records = BamRecord[]
    try
        for record in reader
            push!(records, record)
        end
    finally
        close(reader)
    end
    return BamFile(reader.header, records)
end

"""
    stream_bam(path)

Open a BAM file as a lazy pure-Julia record iterator.
"""
function stream_bam(path::String)

    return BamReader(path)
end

"""
    stream_bam(path, region)

Open a BAM file as a lazy iterator of records overlapping a region.
"""
function stream_bam(path::String, region::GenomicRanges.GenomicInterval)

    io, header = _open_bam_stream(path)
    try
        refindex = _reference_index(header, region.chrom === nothing ? nothing : String(region.chrom))

        if refindex == 0
            _close_bam_io!(io)
            return BamEmptyReader(header)
        end

        index_path = string(path, ".bai")
        if !isfile(index_path)
            return BamRegionScanReader(io, header, region, false)
        end

        index = _read_bam_index(index_path)
        chunks = _bam_region_chunks(index, refindex, region.left, region.right)
        isempty(chunks) && return BamRegionScanReader(io, header, region, false)

        return BamIndexedRegionReader(io, header, region, chunks, 1, UInt64(0), Set{UInt64}(), false)
    catch
        _close_bam_io!(io)
        rethrow()
    end
end

"""
    _bam_scan_region(path, header, region)

Scan a BAM file sequentially for alignments overlapping a region.
"""
function _bam_scan_region(path::String, header::BamHeader, region::GenomicRanges.GenomicInterval)
    bam = _bam_collect(BamRegionScanReader(path, region))
    return header == bam.header ? bam : BamFile(header, bam.records)
end

"""
    _bam_region_from_index(path, header, region, index)

Use a BAM index to read only alignments overlapping a region.
"""
function _bam_region_from_index(path::String, header::BamHeader, region::GenomicRanges.GenomicInterval, index::BamIndex)
    refindex = _reference_index(header, region.chrom === nothing ? nothing : String(region.chrom))
    refindex == 0 && return BamFile(header, BamRecord[])

    chunks = _bam_region_chunks(index, refindex, region.left, region.right)
    isempty(chunks) && return _bam_scan_region(path, header, region)

    io, stream_header = _open_bam_stream(path)
    reader = BamIndexedRegionReader(io, stream_header, region, chunks, 1, UInt64(0), Set{UInt64}(), false)
    bam = _bam_collect(reader)
    return header == bam.header ? bam : BamFile(header, bam.records)
end

"""
    read_bam(path)

Read a BAM file into an in-memory `BamFile`.
"""
function read_bam(path::String; materialize::Bool=true, prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    reader = stream_bam(path)
    if materialize
        provenance_hash = _ctx === nothing ? nothing : file_provenance_hash(path)
        bam = _bam_collect(reader)
        provenance_hash === nothing || (bam.metadata[PROVENANCE_HASH_KEY] = provenance_hash)
        _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
            root = register_provenance!(_ctx, "read_bam"; parents=String[], parameters=(source=path, materialized=materialize, hash=provenance_hash, record_count=length(bam.records)))
            register_container_provenance!(_ctx, bam, "read_bam"; parents=[root.id], parameters=(source=path, materialized=materialize, record_count=length(bam.records)), provenance_hash=provenance_hash)
        end
        return bam
    end
    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(path)))
        register_provenance!(_ctx, "read_bam"; parents=String[], parameters=(source=path, materialized=materialize, hash=provenance_hash))
    end
    return reader
end

"""
    read_bam(path, region)

Read only alignments overlapping a genomic interval.
"""
function read_bam(path::String, region::GenomicRanges.GenomicInterval; materialize::Bool=true, prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    reader = stream_bam(path, region)
    if materialize
        provenance_hash = _ctx === nothing ? nothing : file_provenance_hash(path)
        bam = _bam_collect(reader)
        provenance_hash === nothing || (bam.metadata[PROVENANCE_HASH_KEY] = provenance_hash)
        _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
            root = register_provenance!(_ctx, "read_bam"; parents=String[], parameters=(source=path, region=string(region), materialized=materialize, hash=provenance_hash, record_count=length(bam.records)))
            register_container_provenance!(_ctx, bam, "read_bam"; parents=[root.id], parameters=(source=path, region=string(region), materialized=materialize, record_count=length(bam.records)), provenance_hash=provenance_hash)
        end
        return bam
    end
    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(path)))
        register_provenance!(_ctx, "read_bam"; parents=String[], parameters=(source=path, region=string(region), materialized=materialize, hash=provenance_hash))
    end
    return reader
end

"""
    write_bam(path, records; header=nothing, write_index=true)

Write BAM records to disk and optionally create an index.
"""
function write_bam(path::String, records; header::Union{Nothing,BamHeader}=nothing, write_index::Bool=true, prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    _assert_pure_julia_alignment_path(path)
    bam = records isa BamFile ? records : BamFile(records; header=header)
    index = BamIndex([Dict{UInt32,Vector{BamChunk}}() for _ in bam.header.references], [Dict{Int,UInt64}() for _ in bam.header.references])

    open(BGZFStreams.BGZFStream, path, "w") do io
        write(io, _BAM_MAGIC)
        text = String(bam.header.text)
        write(io, Int32(ncodeunits(text)))
        write(io, codeunits(text))
        write(io, Int32(length(bam.header.references)))
        for reference in bam.header.references
            write(io, Int32(ncodeunits(reference.name) + 1))
            write(io, codeunits(reference.name))
            write(io, UInt8(0x00))
            write(io, Int32(reference.length))
        end

        for record in bam.records
            _write_bam_record(io, bam.header, record, write_index ? index : nothing)
        end
    end

    if write_index
        _write_bam_index(string(path, ".bai"), bam.header, index)
    end

    _ctx = active_provenance_context(_ctx)
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(path)))
        _register_bam_result!(_ctx, "write_bam"; parents=provenance_parent_ids(bam), parameters=(output=path, header_reference_count=length(bam.header.references), record_count=length(bam.records), write_index=write_index, hash=provenance_hash))
    end

    _ctx = active_provenance_context(_ctx)

    return path
end
