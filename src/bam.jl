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

const _BAM_CIGAR_OP_SYMBOLS = ('M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X')

const _BAM_CIGAR_OP_ENCODE_LUT = let lut = fill(UInt32(0xFFFFFFFF), 256)
  lut[Int('M')] = UInt32(0)
  lut[Int('I')] = UInt32(1)
  lut[Int('D')] = UInt32(2)
  lut[Int('N')] = UInt32(3)
  lut[Int('S')] = UInt32(4)
  lut[Int('H')] = UInt32(5)
  lut[Int('P')] = UInt32(6)
  lut[Int('=')] = UInt32(7)
  lut[Int('X')] = UInt32(8)
  lut
end

const _BAM_SEQUENCE_DECODE_LUT = ('=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N')

const _BAM_SEQUENCE_ENCODE_LUT = let lut = fill(0x0f, 256)
  lut[Int('A')] = 0x1;
  lut[Int('a')] = 0x1
  lut[Int('C')] = 0x2;
  lut[Int('c')] = 0x2
  lut[Int('M')] = 0x3;
  lut[Int('m')] = 0x3
  lut[Int('G')] = 0x4;
  lut[Int('g')] = 0x4
  lut[Int('R')] = 0x5;
  lut[Int('r')] = 0x5
  lut[Int('S')] = 0x6;
  lut[Int('s')] = 0x6
  lut[Int('V')] = 0x7;
  lut[Int('v')] = 0x7
  lut[Int('T')] = 0x8;
  lut[Int('t')] = 0x8
  lut[Int('W')] = 0x9;
  lut[Int('w')] = 0x9
  lut[Int('Y')] = 0xA;
  lut[Int('y')] = 0xA
  lut[Int('H')] = 0xB;
  lut[Int('h')] = 0xB
  lut[Int('K')] = 0xC;
  lut[Int('k')] = 0xC
  lut[Int('D')] = 0xD;
  lut[Int('d')] = 0xD
  lut[Int('B')] = 0xE;
  lut[Int('b')] = 0xE
  lut[Int('N')] = 0xF;
  lut[Int('n')] = 0xF
  lut
end

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

Textual BAM header plus its reference list and fast cached reference index map.
"""
struct BamHeader
  text::String
  references::Vector{BamReference}
  ref_map::Dict{String,Int}

  function BamHeader(text::String, references::Vector{BamReference})
    ref_map = Dict{String,Int}(ref.name => i for (i, ref) in enumerate(references))
    return new(text, references, ref_map)
  end
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

Base.hash(x::BamReference, h::UInt) = hash(x.name, hash(x.length, h))

"""
    Base.:(==)(left, right)

Test two BAM CIGAR operations for value equality.
"""
function Base.:(==)(left::BamCigarOp, right::BamCigarOp)
  return isequal(left.op, right.op) && isequal(left.length, right.length)
end

Base.hash(x::BamCigarOp, h::UInt) = hash(x.op, hash(x.length, h))

"""
    Base.:(==)(left, right)

Test two BAM records for value equality.
"""
function Base.:(==)(left::BamRecord, right::BamRecord)
  return isequal(left.qname, right.qname) && isequal(left.flag, right.flag) && isequal(left.refname, right.refname) && isequal(left.pos, right.pos) && isequal(left.mapq, right.mapq) && isequal(left.cigar, right.cigar) && isequal(left.mate_refname, right.mate_refname) && isequal(left.mate_pos, right.mate_pos) && isequal(left.template_length, right.template_length) && isequal(left.sequence, right.sequence) && isequal(left.quality, right.quality) && isequal(left.tags, right.tags)
end

Base.hash(x::BamRecord, h::UInt) = hash((x.qname, x.flag, x.refname, x.pos, x.mapq, x.cigar,
    x.mate_refname, x.mate_pos, x.template_length, x.sequence, x.quality, x.tags), h)

"""
    Base.:(==)(left, right)

Test two BAM headers for value equality.
"""
function Base.:(==)(left::BamHeader, right::BamHeader)
  return isequal(left.text, right.text) && isequal(left.references, right.references)
end

Base.hash(x::BamHeader, h::UInt) = hash(x.text, hash(x.references, h))

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
  len = Int(length)
  if len == 0
    return BioSequence{DNAAlphabet}(UInt8[]; validate=false)
  end

  buffer = Vector{UInt8}(undef, len)
  written = 0
  @inbounds for byte in bytes
    high = (byte >> 4) & 0x0f
    low = byte & 0x0f
    if written < len
      written += 1
      buffer[written] = UInt8(_BAM_SEQUENCE_DECODE_LUT[high+1])
    end
    if written < len
      written += 1
      buffer[written] = UInt8(_BAM_SEQUENCE_DECODE_LUT[low+1])
    end
  end

  return BioSequence{DNAAlphabet}(buffer; validate=false)
end

"""
    _encode_bam_sequence(sequence)

Encode a nucleotide sequence into BAM nibble-packed bytes.
"""
function _encode_bam_sequence(sequence::BioSequence{DNAAlphabet})
  len = length(sequence.data)
  encoded = Vector{UInt8}(undef, cld(len, 2))
  buffer = UInt8(0)
  byte_idx = 0
  use_high = true

  @inbounds for i in 1:len
    char_val = sequence.data[i]
    code = _BAM_SEQUENCE_ENCODE_LUT[char_val]
    if use_high
      buffer = UInt8(code << 4)
      use_high = false
    else
      byte_idx += 1
      encoded[byte_idx] = buffer | UInt8(code)
      use_high = true
    end
  end

  if !use_high
    byte_idx += 1
    encoded[byte_idx] = buffer
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
  return String(bytes .+ UInt8(33))
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
@inline function _decode_bam_cigar(value::UInt32)
  op_code = value & 0x0f
  op = op_code < 9 ? _BAM_CIGAR_OP_SYMBOLS[op_code+1] : 'M'
  len = Int(value >> 4)
  return BamCigarOp(op, len)
end

"""
    _encode_bam_cigar(op)

Encode a `BamCigarOp` into a packed BAM CIGAR value.
"""
@inline function _encode_bam_cigar(op::BamCigarOp)
  char_val = UInt8(op.op)
  code = _BAM_CIGAR_OP_ENCODE_LUT[char_val]
  code == 0xFFFFFFFF && throw(ArgumentError("unsupported BAM CIGAR op $(op.op)"))
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
    bytes = readuntil(io, 0x00)
    return String(bytes)
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

    key_bytes = read(io, 2)
    key = String(key_bytes)
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
  return header.ref_map
end

"""
    _reference_index(header, refname)

Resolve a reference name to its BAM index in O(1) time.
"""
@inline function _reference_index(header::BamHeader, refname::Union{Nothing,String})
  refname === nothing && return 0
  return get(header.ref_map, refname, 0)
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
  span = _bam_reference_span(record.cigar)
  span <= 0 && return 0, -1
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
  local block_size
  try
    block_size = Int(read(io, Int32))
  catch err
    err isa EOFError && return nothing, nothing
    rethrow()
  end
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

  if n_cigar == 2 && haskey(tags, "CG")
    cg_val = tags["CG"]
    if cg_val isa AbstractVector
      cigar = BamCigarOp[_decode_bam_cigar(UInt32(op)) for op in cg_val]
    end
  end

  refname = _reference_name(header, refid)
  mate_refname = _reference_name(header, next_refid)
  return record_start, BamRecord(qname, flag, refname, pos, mapq, cigar, mate_refname, next_pos, tlen, sequence, quality, tags)
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
  refindex = record.refname === nothing ? 0 : get(header.ref_map, record.refname) do
    throw(ArgumentError("reference '$(record.refname)' not found in BAM header"))
  end
  mate_refindex = record.mate_refname === nothing ? 0 : get(header.ref_map, record.mate_refname) do
    throw(ArgumentError("reference '$(record.mate_refname)' not found in BAM header"))
  end
  record_tags = record.tags
  cigar_to_encode = record.cigar
  if length(record.cigar) > 65535
    record_tags = copy(record.tags)
    record_tags["CG"] = UInt32[_encode_bam_cigar(op) for op in record.cigar]
    cigar_to_encode = BamCigarOp[BamCigarOp(length(record.sequence), 'S'), BamCigarOp(0, 'N')]
  end

  cigar_bytes = UInt32[_encode_bam_cigar(op) for op in cigar_to_encode]
  sequence_bytes = _encode_bam_sequence(record.sequence)
  quality_bytes = _encode_bam_quality(record.quality, length(record.sequence))
  tag_io = IOBuffer()
  _write_bam_tags(tag_io, record_tags)
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
      for window in 0:max(n_intv-1, -1)
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
    if isempty(chunks)
      _close_bam_io!(io)
      return BamEmptyReader(header)
    end

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
  isempty(chunks) && return BamFile(header, BamRecord[])

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
    register_provenance!(_ctx, "read_bam"; parents=String[], parameters=(source=path, materialized=materialize, hash=nothing))
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
    register_provenance!(_ctx, "read_bam"; parents=String[], parameters=(source=path, region=string(region), materialized=materialize, hash=nothing))
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

# ==============================================================================
# SAM Flag Predicates & Utilities (Bioconductor / Rsamtools Parity)
# ==============================================================================

const SAM_FLAG_PAIRED = UInt16(0x0001)
const SAM_FLAG_PROPER_PAIR = UInt16(0x0002)
const SAM_FLAG_UNMAPPED = UInt16(0x0004)
const SAM_FLAG_MATE_UNMAPPED = UInt16(0x0008)
const SAM_FLAG_REVERSE = UInt16(0x0010)
const SAM_FLAG_MATE_REVERSE = UInt16(0x0020)
const SAM_FLAG_READ1 = UInt16(0x0040)
const SAM_FLAG_READ2 = UInt16(0x0080)
const SAM_FLAG_SECONDARY = UInt16(0x0100)
const SAM_FLAG_QC_FAIL = UInt16(0x0200)
const SAM_FLAG_DUPLICATE = UInt16(0x0400)
const SAM_FLAG_SUPPLEMENTARY = UInt16(0x0800)

is_paired(record::BamRecord) = (record.flag & SAM_FLAG_PAIRED) != 0
is_proper_pair(record::BamRecord) = (record.flag & SAM_FLAG_PROPER_PAIR) != 0
is_unmapped(record::BamRecord) = (record.flag & SAM_FLAG_UNMAPPED) != 0
is_mate_unmapped(record::BamRecord) = (record.flag & SAM_FLAG_MATE_UNMAPPED) != 0
is_reverse_strand(record::BamRecord) = (record.flag & SAM_FLAG_REVERSE) != 0
is_mate_reverse_strand(record::BamRecord) = (record.flag & SAM_FLAG_MATE_REVERSE) != 0
is_read1(record::BamRecord) = (record.flag & SAM_FLAG_READ1) != 0
is_read2(record::BamRecord) = (record.flag & SAM_FLAG_READ2) != 0
is_secondary(record::BamRecord) = (record.flag & SAM_FLAG_SECONDARY) != 0
is_qc_failed(record::BamRecord) = (record.flag & SAM_FLAG_QC_FAIL) != 0
is_duplicate(record::BamRecord) = (record.flag & SAM_FLAG_DUPLICATE) != 0
is_supplementary(record::BamRecord) = (record.flag & SAM_FLAG_SUPPLEMENTARY) != 0
ismapped(record::BamRecord) = !is_unmapped(record)

alignlength(record::BamRecord) = _bam_reference_span(record.cigar)
leftposition(record::BamRecord) = record.pos < 0 ? 0 : Int(record.pos) + 1
rightposition(record::BamRecord) = record.pos < 0 ? 0 : Int(record.pos) + alignlength(record)
readname(record::BamRecord) = record.qname
cigar_rle(record::BamRecord) = ([op.op for op in record.cigar], [op.length for op in record.cigar])

"""
    filter_bam(reader; flags_req=0, flags_no=0, min_mapq=0)

Stream alignments from a BAM reader filtering by required flags (`flags_req`),
excluded flags (`flags_no`), and minimum mapping quality (`min_mapq`).
"""
function filter_bam(reader::AbstractBamRecordReader; flags_req::Integer=0, flags_no::Integer=0, min_mapq::Integer=0)
  records = BamRecord[]
  for record in reader
    (record.flag & UInt16(flags_req)) == UInt16(flags_req) || continue
    (record.flag & UInt16(flags_no)) == 0 || continue
    record.mapq >= UInt8(min_mapq) || continue
    push!(records, record)
  end
  return records
end

function filter_bam(path::String; kwargs...)
  reader = stream_bam(path)
  try
    return filter_bam(reader; kwargs...)
  finally
    close(reader)
  end
end

"""
    count_bam(path; region=nothing)

Fast alignment record count for a BAM file or region.
"""
function count_bam(path::String; region::Union{Nothing,GenomicRanges.GenomicInterval}=nothing)
  reader = region === nothing ? stream_bam(path) : stream_bam(path, region)
  count = 0
  try
    for _ in reader
      count += 1
    end
  finally
    close(reader)
  end
  return count
end

"""
    bam_coverage(records, references; exclude_flags=SAM_FLAG_SECONDARY | SAM_FLAG_QC_FAIL | SAM_FLAG_DUPLICATE | SAM_FLAG_SUPPLEMENTARY) -> Dict{String, Vector{Int}}

Compute per-base reference coverage across all mapped alignments, respecting CIGAR ops ('M', 'D', '=', 'X').
"""
function bam_coverage(records::AbstractVector{<:BamRecord}, references::AbstractVector{<:BamReference}; exclude_flags::Integer=SAM_FLAG_SECONDARY | SAM_FLAG_QC_FAIL | SAM_FLAG_DUPLICATE | SAM_FLAG_SUPPLEMENTARY)
  cov = Dict{String,Vector{Int}}(ref.name => zeros(Int, ref.length) for ref in references)

  for rec in records
    (rec.refname === nothing || rec.pos < 0 || is_unmapped(rec)) && continue
    (rec.flag & UInt16(exclude_flags)) == 0 || continue
    ref_cov = get(cov, rec.refname, nothing)
    ref_cov === nothing && continue
    curr_pos = Int(rec.pos) + 1
    for op in rec.cigar
      if op.op in ('M', '=', 'X')
        stop_pos = min(curr_pos + op.length - 1, length(ref_cov))
        if curr_pos <= length(ref_cov) && stop_pos >= 1
          start_k = max(curr_pos, 1)
          @inbounds for k in start_k:stop_pos
            ref_cov[k] += 1
          end
        end
        curr_pos += op.length
      elseif op.op in ('D', 'N')
        curr_pos += op.length
      end
    end
  end
  return cov
end

bam_coverage(file::BamFile; kwargs...) = bam_coverage(file.records, file.header.references; kwargs...)

"""
    alignments_to_interval_collection(records) -> IntervalCollection

Convert BAM alignment records to an `IntervalCollection` of `GenomicInterval` objects.
"""
function alignments_to_interval_collection(records::AbstractVector{<:BamRecord})
  intervals = GenomicRanges.GenomicInterval[]
  for rec in records
    (rec.refname === nothing || rec.pos < 0) && continue
    span = _bam_reference_span(rec.cigar)
    span <= 0 && continue
    strand = is_reverse_strand(rec) ? '-' : '+'
    left = Int(rec.pos) + 1
    right = left + span - 1
    metadata = Dict{String,Any}("qname" => rec.qname, "mapq" => rec.mapq, "flag" => rec.flag)
    push!(intervals, GenomicRanges.GenomicInterval(rec.refname, left, right, strand, metadata))
  end
  return GenomicRanges.build_collection(intervals)
end

alignments_to_interval_collection(file::BamFile) = alignments_to_interval_collection(file.records)

# Alias matching Bioconductor's granges() convention
const granges = alignments_to_interval_collection

# ==============================================================================
# SAM Text Format I/O (Bioconductor / SAMtools Parity)
# ==============================================================================

function _parse_sam_cigar(cigar_str::AbstractString)
  (cigar_str == "*" || isempty(cigar_str)) && return BamCigarOp[]
  ops = BamCigarOp[]
  len_val = 0
  for c in cigar_str
    if isdigit(c)
      len_val = len_val * 10 + Int(c - '0')
    else
      push!(ops, BamCigarOp(c, len_val))
      len_val = 0
    end
  end
  return ops
end

function _parse_sam_tag_value(type_char::AbstractString, value_str::AbstractString)
  if type_char == "i"
    return parse(Int32, value_str)
  elseif type_char == "c"
    return parse(Int8, value_str)
  elseif type_char == "C"
    return parse(UInt8, value_str)
  elseif type_char == "s"
    return parse(Int16, value_str)
  elseif type_char == "S"
    return parse(UInt16, value_str)
  elseif type_char == "I"
    return parse(UInt32, value_str)
  elseif type_char == "f"
    return parse(Float32, value_str)
  elseif type_char == "d"
    return parse(Float64, value_str)
  elseif type_char == "A"
    return isempty(value_str) ? ' ' : value_str[1]
  elseif type_char == "B"
    length(value_str) < 2 && return Int32[]
    subtype = value_str[1]
    raw_nums = value_str[3:end]
    isempty(raw_nums) && return subtype in ('f', 'd') ? Float32[] : Int32[]
    nums = Base.split(raw_nums, ',')
    if subtype == 'c'
      return parse.(Int8, nums)
    elseif subtype == 'C'
      return parse.(UInt8, nums)
    elseif subtype == 's'
      return parse.(Int16, nums)
    elseif subtype == 'S'
      return parse.(UInt16, nums)
    elseif subtype == 'I'
      return parse.(UInt32, nums)
    elseif subtype == 'f'
      return parse.(Float32, nums)
    elseif subtype == 'd'
      return parse.(Float64, nums)
    else
      return parse.(Int32, nums)
    end
  else  # "Z" or "H"
    return String(value_str)
  end
end

"""
    read_sam(io_or_path) -> BamFile

Read SAM format textual alignment records into an in-memory `BamFile`.
"""
function read_sam(io::IO)
  header_lines = String[]
  references = BamReference[]
  records = BamRecord[]

  for line in eachline(io)
    isempty(line) && continue
    if startswith(line, '@')
      push!(header_lines, line)
      if startswith(line, "@SQ")
        parts = Base.split(line, '\t')
        sn = ""
        ln = 0
        for part in parts[2:end]
          if startswith(part, "SN:")
            sn = String(part[4:end])
          elseif startswith(part, "LN:")
            ln = parse(Int, part[4:end])
          end
        end
        !isempty(sn) && push!(references, BamReference(sn, ln))
      end
      continue
    end

    fields = Base.split(line, '\t')
    length(fields) >= 11 || continue
    qname = String(fields[1])
    flag = parse(UInt16, fields[2])
    rname = fields[3] == "*" ? nothing : String(fields[3])
    pos = parse(Int32, fields[4]) - Int32(1)
    mapq = parse(UInt8, fields[5])
    cigar_str = fields[6]
    cigar_ops = _parse_sam_cigar(cigar_str)
    rnext = fields[7] == "*" ? nothing : (fields[7] == "=" ? rname : String(fields[7]))
    pnext = parse(Int32, fields[8]) - Int32(1)
    tlen = parse(Int32, fields[9])
    seq_str = fields[10]
    qual_str = fields[11]

    seq = seq_str == "*" ? BioSequence{DNAAlphabet}(UInt8[]; validate=false) : BioSequence{DNAAlphabet}(String(seq_str))
    qual = qual_str == "*" ? missing : String(qual_str)

    tags = Dict{String,Any}()
    for tag_str in fields[12:end]
      tparts = Base.split(tag_str, ':', limit=3)
      length(tparts) == 3 || continue
      tags[String(tparts[1])] = _parse_sam_tag_value(tparts[2], tparts[3])
    end

    push!(records, BamRecord(qname, flag, rname, pos, mapq, cigar_ops, rnext, pnext, tlen, seq, qual, tags))
  end

  hdr_text = join(header_lines, "\n")
  if isempty(references)
    references = _bam_infer_reference_lengths(records)
  end
  hdr = BamHeader(hdr_text, references)
  return BamFile(hdr, records)
end

function read_sam(path::String)
  open(path, "r") do io
    return read_sam(io)
  end
end

"""
    write_sam(path_or_io, file)

Write BAM records to SAM text format.
"""
function write_sam(io::IO, file::BamFile)
  println(io, file.header.text)
  for rec in file.records
    rname = rec.refname === nothing ? "*" : rec.refname
    pos_1based = rec.pos < 0 ? 0 : rec.pos + 1
    cigar_str = isempty(rec.cigar) ? "*" : join(["$(op.length)$(op.op)" for op in rec.cigar])
    rnext = rec.mate_refname === nothing ? "*" : (rec.mate_refname == rec.refname ? "=" : rec.mate_refname)
    pnext_1based = rec.mate_pos < 0 ? 0 : rec.mate_pos + 1
    seq_str = isempty(rec.sequence) ? "*" : String(rec.sequence)
    qual_str = rec.quality === missing ? "*" : rec.quality

    tag_strs = String[]
    for (k, v) in rec.tags
      if v isa AbstractVector
        subtype_char = if eltype(v) <: Integer
          eltype(v) <: Int8 ? "c" : eltype(v) <: UInt8 ? "C" :
                                    eltype(v) <: Int16 ? "s" : eltype(v) <: UInt16 ? "S" :
                                                               eltype(v) <: Unsigned ? "I" : "i"
        else
          "f"
        end
        push!(tag_strs, "$k:B:$subtype_char," * join(v, ","))
      elseif v isa Integer
        push!(tag_strs, "$k:i:$v")
      elseif v isa AbstractFloat
        push!(tag_strs, "$k:f:$v")
      elseif v isa Char
        push!(tag_strs, "$k:A:$v")
      else
        push!(tag_strs, "$k:Z:$v")
      end
    end
    tags_line = isempty(tag_strs) ? "" : "\t" * join(tag_strs, "\t")

    println(io, "$(rec.qname)\t$(rec.flag)\t$rname\t$pos_1based\t$(rec.mapq)\t$cigar_str\t$rnext\t$pnext_1based\t$(rec.template_length)\t$seq_str\t$qual_str$tags_line")
  end
end

function write_sam(path::String, file::BamFile)
  open(path, "w") do io
    write_sam(io, file)
  end
end
