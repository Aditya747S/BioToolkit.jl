# ==============================================================================
# biotypes.jl — Native Biological Type System
#
# Provides parametric type safety for biological alphabets and sequences
# without any external dependencies. Sequences are stored as compact
# byte vectors with typed wrappers that prevent accidental cross-alphabet
# operations (e.g., translating a protein sequence as if it were DNA).
#
# Design decisions:
#   - Alphabets are singleton types used purely for dispatch, not storage.
#   - Sequences store raw bytes (one byte per symbol) for simplicity and
#     interoperability with existing String-based code. A 2-bit encoding
#     would halve DNA memory but complicate slicing and printing; given
#     that BioToolkit targets interactive analysis workloads (not assembly),
#     the byte-per-symbol layout is the right trade-off.
#   - All constructors validate input against the alphabet.
#   - String conversion is zero-copy where possible.
# ==============================================================================

# ---- Alphabets ---------------------------------------------------------------

"""
    BioAlphabet

Abstract supertype for all biological sequence alphabets.

Concrete subtypes exist solely for dispatch — they carry no data.
Use `symbols(A)` to retrieve the set of valid byte values for alphabet `A`.
"""
abstract type BioAlphabet end

"""
    DNAAlphabet <: BioAlphabet

Standard IUPAC DNA alphabet (A, C, G, T, plus ambiguity codes N, R, Y, S, W, K, M, B, D, H, V).
Case-insensitive: lowercase inputs are accepted and stored as uppercase.
"""
struct DNAAlphabet <: BioAlphabet end

"""
    RNAAlphabet <: BioAlphabet

Standard IUPAC RNA alphabet (A, C, G, U, plus ambiguity codes).
Case-insensitive: lowercase inputs are accepted and stored as uppercase.
"""
struct RNAAlphabet <: BioAlphabet end

"""
    AminoAcidAlphabet <: BioAlphabet

Standard IUPAC amino acid alphabet (20 canonical + B, Z, X, *, U, O).
Case-insensitive.
"""
struct AminoAcidAlphabet <: BioAlphabet end

# Valid symbol sets — used for input validation.  The sets contain uppercase
# bytes; validation normalises input to uppercase before checking membership.

const _DNA_VALID_BYTES = Set{UInt8}(UInt8.(collect("ACGTNRYSWKMBDHV-")))
const _RNA_VALID_BYTES = Set{UInt8}(UInt8.(collect("ACGUNRYSWKMBDHV-")))
const _AA_VALID_BYTES = Set{UInt8}(UInt8.(collect("ACDEFGHIKLMNPQRSTVWYXBZUO*-")))

"""
    symbols(::Type{A}) -> Set{UInt8}

Return the set of valid uppercase byte symbols for alphabet `A`.
"""
symbols(::Type{DNAAlphabet}) = _DNA_VALID_BYTES
symbols(::Type{RNAAlphabet}) = _RNA_VALID_BYTES
symbols(::Type{AminoAcidAlphabet}) = _AA_VALID_BYTES

"""
    gap(::Type{A}) -> UInt8

Return the gap character byte (fixed as '-') for alphabet `A`.
"""
gap(::Type{<:BioAlphabet}) = UInt8('-')

"""
    isvalid_symbol(::Type{A}, byte) -> Bool

Check whether `byte` (after uppercasing) is a member of alphabet `A`.
"""
@inline function isvalid_symbol(::Type{A}, byte::UInt8) where {A<:BioAlphabet}
    upper = byte < 0x61 ? byte : (byte <= 0x7a ? byte - 0x20 : byte)
    return upper in symbols(A)
end

# ---- Biological Symbol System & 1-Hot Bitwise Compatibility ------------------

export
    BioSymbol, DNA, RNA, AminoAcid,
    DNA_Gap, DNA_A, DNA_C, DNA_G, DNA_T, DNA_M, DNA_R, DNA_W, DNA_S, DNA_Y, DNA_K, DNA_V, DNA_H, DNA_D, DNA_B, DNA_N,
    RNA_Gap, RNA_A, RNA_C, RNA_G, RNA_U, RNA_M, RNA_R, RNA_W, RNA_S, RNA_Y, RNA_K, RNA_V, RNA_H, RNA_D, RNA_B, RNA_N,
    ACGT, ACGU,
    AA_A, AA_R, AA_N, AA_D, AA_C, AA_Q, AA_E, AA_G, AA_H, AA_I, AA_L, AA_K, AA_M, AA_F, AA_P, AA_S, AA_T, AA_W, AA_Y, AA_V, AA_O, AA_U, AA_B, AA_J, AA_Z, AA_X, AA_Term, AA_Gap,
    parse_amino_acid, iscompatible, compatbits, isambiguous, iscertain, isgap, isGC, ispurine, ispyrimidine, complement

"""
    BioSymbol

Abstract supertype representing single biological characters (nucleotides or amino acids).
Concrete primitive types `DNA`, `RNA`, and `AminoAcid` provide 1-hot bitwise ambiguity
matching and single-cycle operations.
"""
abstract type BioSymbol end

primitive type DNA <: BioSymbol 8 end
primitive type RNA <: BioSymbol 8 end
primitive type AminoAcid <: BioSymbol 8 end

encoded_data(x::BioSymbol) = reinterpret(UInt8, x)
encode(::Type{T}, x::UInt8) where {T<:BioSymbol} = reinterpret(T, x)

Base.broadcastable(x::BioSymbol) = (x,)

const _DNA_TO_CHAR = Char['-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N']
const _RNA_TO_CHAR = Char['-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'U', 'W', 'Y', 'H', 'K', 'D', 'B', 'N']

const _AA_TO_CHAR = Char[
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
    'O', 'U', 'B', 'J', 'Z', 'X', '*', '-'
]

const _AA_COMPATBITS = UInt32[
    1<<0, 1<<1, 1<<2, 1<<3, 1<<4, 1<<5, 1<<6, 1<<7, 1<<8, 1<<9,
    1<<10, 1<<11, 1<<12, 1<<13, 1<<14, 1<<15, 1<<16, 1<<17, 1<<18, 1<<19,
    1<<20, 1<<21, (1<<2)|(1<<3), (1<<9)|(1<<10), (1<<5)|(1<<6), (1<<22)-1, 1<<26, 0
]

Base.convert(::Type{Char}, x::DNA) = @inbounds _DNA_TO_CHAR[encoded_data(x) + 1]
Base.convert(::Type{Char}, x::RNA) = @inbounds _RNA_TO_CHAR[encoded_data(x) + 1]
Base.convert(::Type{Char}, x::AminoAcid) = @inbounds _AA_TO_CHAR[min(Int(encoded_data(x)) + 1, 28)]

Base.Char(x::BioSymbol) = convert(Char, x)
Base.show(io::IO, x::BioSymbol) = print(io, Char(x))
Base.print(io::IO, x::BioSymbol) = print(io, Char(x))

@inline compatbits(x::DNA) = encoded_data(x)
@inline compatbits(x::RNA) = encoded_data(x)
@inline compatbits(x::AminoAcid) = @inbounds _AA_COMPATBITS[min(Int(encoded_data(x)) + 1, 28)]

"""
    iscompatible(x::BioSymbol, y::BioSymbol) -> Bool

Test if two biological symbols are compatible using bitwise 1-hot bitmasks.
"""
@inline iscompatible(x::S, y::S) where {S<:BioSymbol} = (compatbits(x) & compatbits(y)) != 0

@inline isambiguous(x::DNA) = count_ones(encoded_data(x)) > 1
@inline isambiguous(x::RNA) = count_ones(encoded_data(x)) > 1
@inline isambiguous(x::AminoAcid) = 0x16 <= encoded_data(x) <= 0x19

@inline iscertain(x::BioSymbol) = !isambiguous(x) && !isgap(x)
@inline isgap(x::BioSymbol) = encoded_data(x) == 0x00 || Char(x) == '-'

@inline isGC(x::Union{DNA, RNA}) = (encoded_data(x) != 0) && ((encoded_data(x) & 0b1001) == 0)
@inline ispurine(x::Union{DNA, RNA}) = (encoded_data(x) != 0) && ((encoded_data(x) & 0b1010) == 0)
@inline ispyrimidine(x::Union{DNA, RNA}) = (encoded_data(x) != 0) && ((encoded_data(x) & 0b0101) == 0)

"""
    complement(nt::DNA) / complement(nt::RNA)

Single-cycle bitwise complement of a nucleotide symbol.
"""
function complement(nt::DNA)
    u64 = 0xf7b3d591e6a2c480 >>> ((4 * encoded_data(nt)) & 63)
    return encode(DNA, UInt8(u64 & 0x0f))
end

function complement(nt::RNA)
    u64 = 0xf7b3d591e6a2c480 >>> ((4 * encoded_data(nt)) & 63)
    return encode(RNA, UInt8(u64 & 0x0f))
end

function _register_symbol_complement_interop!()
    if isdefined(Main, :complement)
        try
            f = getfield(Main, :complement)
            if f !== complement
                f(nt::DNA) = complement(nt)
                f(nt::RNA) = complement(nt)
            end
        catch
        end
    end
end

# Symbol constants
const DNA_Gap = encode(DNA, 0b0000)
const DNA_A   = encode(DNA, 0b0001)
const DNA_C   = encode(DNA, 0b0010)
const DNA_G   = encode(DNA, 0b0100)
const DNA_T   = encode(DNA, 0b1000)
const DNA_M   = encode(DNA, 0b0011)
const DNA_R   = encode(DNA, 0b0101)
const DNA_W   = encode(DNA, 0b1001)
const DNA_S   = encode(DNA, 0b0110)
const DNA_Y   = encode(DNA, 0b1010)
const DNA_K   = encode(DNA, 0b1100)
const DNA_V   = encode(DNA, 0b0111)
const DNA_H   = encode(DNA, 0b1011)
const DNA_D   = encode(DNA, 0b1101)
const DNA_B   = encode(DNA, 0b1110)
const DNA_N   = encode(DNA, 0b1111)

const RNA_Gap = encode(RNA, 0b0000)
const RNA_A   = encode(RNA, 0b0001)
const RNA_C   = encode(RNA, 0b0010)
const RNA_G   = encode(RNA, 0b0100)
const RNA_U   = encode(RNA, 0b1000)
const RNA_M   = encode(RNA, 0b0011)
const RNA_R   = encode(RNA, 0b0101)
const RNA_W   = encode(RNA, 0b1001)
const RNA_S   = encode(RNA, 0b0110)
const RNA_Y   = encode(RNA, 0b1010)
const RNA_K   = encode(RNA, 0b1100)
const RNA_V   = encode(RNA, 0b0111)
const RNA_H   = encode(RNA, 0b1011)
const RNA_D   = encode(RNA, 0b1101)
const RNA_B   = encode(RNA, 0b1110)
const RNA_N   = encode(RNA, 0b1111)

const ACGT = (DNA_A, DNA_C, DNA_G, DNA_T)
const ACGU = (RNA_A, RNA_C, RNA_G, RNA_U)

const AA_A = encode(AminoAcid, 0x00)
const AA_R = encode(AminoAcid, 0x01)
const AA_N = encode(AminoAcid, 0x02)
const AA_D = encode(AminoAcid, 0x03)
const AA_C = encode(AminoAcid, 0x04)
const AA_Q = encode(AminoAcid, 0x05)
const AA_E = encode(AminoAcid, 0x06)
const AA_G = encode(AminoAcid, 0x07)
const AA_H = encode(AminoAcid, 0x08)
const AA_I = encode(AminoAcid, 0x09)
const AA_L = encode(AminoAcid, 0x0a)
const AA_K = encode(AminoAcid, 0x0b)
const AA_M = encode(AminoAcid, 0x0c)
const AA_F = encode(AminoAcid, 0x0d)
const AA_P = encode(AminoAcid, 0x0e)
const AA_S = encode(AminoAcid, 0x0f)
const AA_T = encode(AminoAcid, 0x10)
const AA_W = encode(AminoAcid, 0x11)
const AA_Y = encode(AminoAcid, 0x12)
const AA_V = encode(AminoAcid, 0x13)
const AA_O = encode(AminoAcid, 0x14)
const AA_U = encode(AminoAcid, 0x15)
const AA_B = encode(AminoAcid, 0x16)
const AA_J = encode(AminoAcid, 0x17)
const AA_Z = encode(AminoAcid, 0x18)
const AA_X = encode(AminoAcid, 0x19)
const AA_Term = encode(AminoAcid, 0x1a)
const AA_Gap  = encode(AminoAcid, 0x1b)

# ---- 3-Letter Amino Acid Parser Engine ----------------------------------------

const _THREE_LETTER_TO_AA = Dict{String, Char}(
    "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D', "CYS" => 'C',
    "GLN" => 'Q', "GLU" => 'E', "GLY" => 'G', "HIS" => 'H', "ILE" => 'I',
    "LEU" => 'L', "LYS" => 'K', "MET" => 'M', "PHE" => 'F', "PRO" => 'P',
    "SER" => 'S', "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
    "SEC" => 'U', "PYL" => 'O', "ASX" => 'B', "GLX" => 'Z', "XLE" => 'J', "XAA" => 'X'
)

"""
    parse_amino_acid(s::Union{AbstractString, Char}) -> Char

Parse a 1-letter or 3-letter amino acid code into a standard single uppercase amino acid character.
Supports PDB codes ("ALA", "CYS", "MET", "SEC", "PYL", "ASX", etc.).
"""
function parse_amino_acid(s::AbstractString)
    clean = strip(s)
    if length(clean) == 1
        c = uppercase(clean[1])
        c in _AA_VALID_BYTES && return c
        throw(ArgumentError("invalid amino acid symbol '$(s)'"))
    elseif length(clean) == 3
        key = uppercase(clean)
        haskey(_THREE_LETTER_TO_AA, key) && return _THREE_LETTER_TO_AA[key]
        throw(ArgumentError("invalid 3-letter amino acid code '$(s)'"))
    else
        throw(ArgumentError("amino acid code must be 1 or 3 characters long"))
    end
end

parse_amino_acid(c::Char) = parse_amino_acid(string(c))

# ---- 2-Bit / Bit-Packed Vector Primitive --------------------------------------

const _BYTE_TO_2BIT = fill(0x03, 256)
for (byte, val) in [(UInt8('A'), 0x00), (UInt8('a'), 0x00), (UInt8('C'), 0x01), (UInt8('c'), 0x01), (UInt8('G'), 0x02), (UInt8('g'), 0x02), (UInt8('T'), 0x03), (UInt8('t'), 0x03), (UInt8('U'), 0x03), (UInt8('u'), 0x03)]
    _BYTE_TO_2BIT[Int(byte)+1] = val
end
const _2BIT_TO_BYTE = UInt8[UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T')]

struct BitPackedVector{Bits} <: AbstractVector{UInt8}
    chunks::Vector{UInt64}
    len::Int

    function BitPackedVector{2}(bytes::AbstractVector{UInt8})
        n = length(bytes)
        num_chunks = (n + 31) ÷ 32
        chunks = zeros(UInt64, num_chunks)
        @inbounds for i in 1:n
            c = UInt64(_BYTE_TO_2BIT[Int(bytes[i]) + 1])
            chunk_idx = ((i - 1) >>> 5) + 1
            shift = ((i - 1) & 0x1f) << 1
            chunks[chunk_idx] |= (c << shift)
        end
        return new{2}(chunks, n)
    end

    function BitPackedVector{2}(chunks::Vector{UInt64}, len::Int)
        return new{2}(chunks, len)
    end
end

Base.size(v::BitPackedVector) = (v.len,)
Base.length(v::BitPackedVector) = v.len
Base.IndexStyle(::Type{<:BitPackedVector}) = IndexLinear()

@inline function Base.getindex(v::BitPackedVector{2}, i::Int)
    @boundscheck checkbounds(v, i)
    chunk_idx = ((i - 1) >>> 5) + 1
    shift = ((i - 1) & 0x1f) << 1
    @inbounds code = (v.chunks[chunk_idx] >>> shift) & 0x03
    return _2BIT_TO_BYTE[code + 1]
end

@inline function Base.setindex!(v::BitPackedVector{2}, val::UInt8, i::Int)
    @boundscheck checkbounds(v, i)
    chunk_idx = ((i - 1) >>> 5) + 1
    shift = ((i - 1) & 0x1f) << 1
    code = UInt64(_BYTE_TO_2BIT[Int(val) + 1])
    @inbounds v.chunks[chunk_idx] = (v.chunks[chunk_idx] & ~(UInt64(0x03) << shift)) | (code << shift)
    return val
end

Base.copy(v::BitPackedVector{2}) = BitPackedVector{2}(copy(v.chunks), v.len)

# ---- Sequence type -----------------------------------------------------------

"""
    BioSequence{A <: BioAlphabet, S <: AbstractVector{UInt8}}

A typed biological sequence over alphabet `A`, stored as a contiguous byte
vector (`Vector{UInt8}`) or 2-bit bit-packed vector (`BitPackedVector{2}`).
"""
struct BioSequence{A<:BioAlphabet, S<:AbstractVector{UInt8}}
    data::S

    function BioSequence{A}(data::S; validate::Bool=true) where {A<:BioAlphabet, S<:AbstractVector{UInt8}}
        if validate
            @inbounds for i in 1:length(data)
                isvalid_symbol(A, data[i]) || throw(ArgumentError(
                    "invalid symbol '$(Char(data[i]))' (byte $(data[i])) for $(A); " *
                    "expected one of: $(join(sort!([Char(b) for b in symbols(A)]), ", "))"))
            end
        end
        return new{A, S}(data)
    end
end

# Convenience aliases
const DNASeq = BioSequence{DNAAlphabet, Vector{UInt8}}
const RNASeq = BioSequence{RNAAlphabet, Vector{UInt8}}
const AASeq = BioSequence{AminoAcidAlphabet, Vector{UInt8}}
const PackedDNASeq = BioSequence{DNAAlphabet, BitPackedVector{2}}

function PackedDNASeq(s::AbstractString; validate::Bool=true)
    bytes = Vector{UInt8}(undef, ncodeunits(s))
    @inbounds for (i, byte) in enumerate(codeunits(s))
        bytes[i] = byte < 0x61 ? byte : (byte <= 0x7a ? byte - 0x20 : byte)
    end
    if validate
        for b in bytes
            b in (UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T'), UInt8('U')) ||
                throw(ArgumentError("PackedDNASeq only supports A/C/G/T/U; got '$(Char(b))' — use BioSequence{DNAAlphabet} for ambiguity codes or gaps"))
        end
    end
    packed = BitPackedVector{2}(bytes)
    return BioSequence{DNAAlphabet}(packed; validate=false)
end

# ---- Constructors ------------------------------------------------------------

"""
    BioSequence{A}(s::String; validate=true)

Construct a `BioSequence{A}` from a string, normalising to uppercase.
"""
function BioSequence{A, Vector{UInt8}}(data::Vector{UInt8}; validate::Bool=true) where {A<:BioAlphabet}
    return BioSequence{A}(data; validate=validate)
end

function BioSequence{A, Vector{UInt8}}(s::AbstractString; validate::Bool=true) where {A<:BioAlphabet}
    bytes = Vector{UInt8}(undef, ncodeunits(s))
    @inbounds for (i, byte) in enumerate(codeunits(s))
        bytes[i] = byte < 0x61 ? byte : (byte <= 0x7a ? byte - 0x20 : byte)
    end
    return BioSequence{A}(bytes; validate=validate)
end

function BioSequence{AminoAcidAlphabet, Vector{UInt8}}(s::AbstractString; validate::Bool=true)
    if occursin('-', s) || occursin(' ', s) || occursin(',', s)
        tokens = Base.split(s, Regex("[- ,]+"))
        if !isempty(tokens) && all(t -> length(strip(t)) == 3, tokens)
            parsed_chars = [parse_amino_acid(t) for t in tokens]
            bytes = UInt8.(parsed_chars)
            return BioSequence{AminoAcidAlphabet}(bytes; validate=validate)
        end
    end

    bytes = Vector{UInt8}(undef, ncodeunits(s))
    @inbounds for (i, byte) in enumerate(codeunits(s))
        bytes[i] = byte < 0x61 ? byte : (byte <= 0x7a ? byte - 0x20 : byte)
    end
    return BioSequence{AminoAcidAlphabet}(bytes; validate=validate)
end

function BioSequence{A, BitPackedVector{2}}(s::AbstractString; validate::Bool=true) where {A<:BioAlphabet}
    return PackedDNASeq(s; validate=validate)
end

function BioSequence{A}(s::AbstractString; validate::Bool=true) where {A<:BioAlphabet}
    return BioSequence{A, Vector{UInt8}}(s; validate=validate)
end

# Note: Named constructors are NOT needed here. Since `DNASeq`, `RNASeq`, and
# `AASeq` are const aliases for `BioSequence{DNAAlphabet}` etc., Julia already
# dispatches `DNASeq("ACGT")` to `BioSequence{DNAAlphabet}(::String)`.
# Explicitly defining `DNASeq(s) = BioSequence{DNAAlphabet}(s)` would create
# infinite recursion because Julia resolves both sides to the same method.

# ---- String-like interface -------------------------------------------

Base.length(seq::BioSequence) = length(seq.data)
Base.sizeof(seq::BioSequence) = length(seq.data)
Base.isempty(seq::BioSequence) = isempty(seq.data)
Base.ncodeunits(seq::BioSequence) = length(seq.data)
Base.codeunits(seq::BioSequence) = seq.data
Base.firstindex(seq::BioSequence) = 1
Base.lastindex(seq::BioSequence) = length(seq.data)
Base.getindex(seq::BioSequence, i::Integer) = Char(seq.data[i])
Base.iterate(seq::BioSequence) = isempty(seq.data) ? nothing : (Char(seq.data[1]), 2)
Base.iterate(seq::BioSequence, i::Int) = i > length(seq.data) ? nothing : (Char(seq.data[i]), i + 1)
Base.eltype(::Type{<:BioSequence}) = Char

function Base.getindex(seq::BioSequence{A}, r::UnitRange{<:Integer}) where {A}
    return BioSequence{A}(seq.data[r]; validate=false)
end

function Base.:(==)(a::BioSequence{A}, b::BioSequence{A}) where {A}
    return a.data == b.data
end

# Compatibility with legacy String-based code and tests.
Base.:(==)(a::BioSequence, b::String) = String(a) == String(b)
Base.:(==)(a::String, b::BioSequence) = String(a) == String(b)
Base.startswith(sequence::BioSequence, prefix::AbstractString) = startswith(String(sequence), String(prefix))
Base.startswith(sequence::BioSequence, prefix::BioSequence) = startswith(String(sequence), String(prefix))

function Base.hash(seq::BioSequence, h::UInt)
    return hash(seq.data, hash(:BioSequence, h))
end

function Base.show(io::IO, seq::BioSequence{A, S}) where {A, S}
    name = S <: BitPackedVector ? "Packed$(A === DNAAlphabet ? "DNASeq" : "BioSequence")" :
           A === DNAAlphabet ? "DNASeq" :
           A === RNAAlphabet ? "RNASeq" :
           A === AminoAcidAlphabet ? "AASeq" : "BioSequence{$(A)}"
    n = length(seq.data)
    if n <= 60
        print(io, name, "(\"", String(collect(seq.data)), "\")")
    else
        bytes = collect(seq.data[1:30])
        bytes_end = collect(seq.data[end-29:end])
        print(io, name, "(\"", String(bytes), "…", String(bytes_end), "\") [", n, " nt]")
    end
end

"""
    String(seq::BioSequence) -> String

Convert a typed biological sequence back to a plain string.
"""
Base.String(seq::BioSequence{A, Vector{UInt8}}) where {A} = String(copy(seq.data))
Base.String(seq::BioSequence{A, BitPackedVector{2}}) where {A} = String(collect(seq.data))

"""
    convert(::Type{String}, seq::BioSequence)

Enable implicit conversion to `String` for interoperability with existing code.
"""
Base.convert(::Type{String}, seq::BioSequence) = String(seq)

# ---- Alphabet queries --------------------------------------------------------

"""
    alphabet(seq::BioSequence{A}) -> Type{A}

Return the alphabet type of a sequence.
"""
alphabet(::BioSequence{A}) where {A} = A

"""
    isdna(seq) -> Bool

Check whether a sequence is a DNA sequence.
"""
isdna(::BioSequence{DNAAlphabet}) = true
isdna(::BioSequence) = false
isdna(::String) = false

"""
    isrna(seq) -> Bool

Check whether a sequence is an RNA sequence.
"""
isrna(::BioSequence{RNAAlphabet}) = true
isrna(::BioSequence) = false
isrna(::String) = false

"""
    isaminoacid(seq) -> Bool

Check whether a sequence is an amino acid sequence.
"""
isaminoacid(::BioSequence{AminoAcidAlphabet}) = true
isaminoacid(::BioSequence) = false
isaminoacid(::String) = false

# ---- Validation helpers ------------------------------------------------------

"""
    validate_sequence(::Type{A}, s::String) -> Bool

Check whether all characters in `s` are valid for alphabet `A`.
"""
function validate_sequence(::Type{A}, s::String) where {A<:BioAlphabet}
    @inbounds for byte in codeunits(s)
        isvalid_symbol(A, byte) || return false
    end
    return true
end

validate_aa_sequence(s::String) = validate_sequence(AminoAcidAlphabet, s)

# ---- Run-Length Encoding (Rle) -----------------------------------------------

"""
    Rle{T}

Run-length encoded vector. Stores repeated values efficiently as (value, length)
pairs. Equivalent to R/Bioconductor's `IRanges::Rle`.

Useful for representing genomic coverage vectors, GC content tracks, and
other signals with long runs of identical values over a large genomic space.

# Fields
- `values::Vector{T}`: The unique values in order.
- `lengths::Vector{Int}`: The number of times each corresponding value is repeated.
"""
struct Rle{T}
    values::Vector{T}
    lengths::Vector{Int}

    function Rle{T}(values::Vector{T}, lengths::Vector{Int}) where {T}
        length(values) == length(lengths) || throw(ArgumentError("values and lengths must have same length"))
        all(lengths .> 0) || throw(ArgumentError("all lengths must be positive"))
        return new{T}(values, lengths)
    end
end

"""
    Rle(data::AbstractVector)

Construct an Rle from a raw vector of values by compressing adjacent identical values.
"""
function Rle(data::AbstractVector{T}) where {T}
    isempty(data) && return Rle{T}(T[], Int[])

    values = T[]
    lengths = Int[]

    current = data[1]
    count = 1

    @inbounds for i in 2:length(data)
        if data[i] == current
            count += 1
        else
            push!(values, current)
            push!(lengths, count)
            current = data[i]
            count = 1
        end
    end
    push!(values, current)
    push!(lengths, count)

    return Rle{T}(values, lengths)
end

Base.length(rle::Rle) = isempty(rle.lengths) ? 0 : sum(rle.lengths)
Base.isempty(rle::Rle) = isempty(rle.values)
Base.eltype(::Type{Rle{T}}) where {T} = T

function Base.getindex(rle::Rle{T}, i::Integer) where {T}
    cumulative = 0
    @inbounds for j in eachindex(rle.values)
        cumulative += rle.lengths[j]
        cumulative >= i && return rle.values[j]
    end
    throw(BoundsError(rle, i))
end

"""
    decode(rle::Rle{T}) -> Vector{T}

Expand an Rle back to its full vector representation.
"""
function decode(rle::Rle{T}) where {T}
    result = Vector{T}(undef, length(rle))
    idx = 1
    @inbounds for j in eachindex(rle.values)
        len = rle.lengths[j]
        val = rle.values[j]
        for _ in 1:len
            result[idx] = val
            idx += 1
        end
    end
    return result
end

function Base.show(io::IO, rle::Rle{T}) where {T}
    n = length(rle.values)
    l = length(rle)
    print(io, "Rle{$T}(", n, " runs, total length ", l, ")")
end

# ---- SummarizedExperiment ----------------------------------------------------

"""
    SummarizedExperiment

A matrix-like container that coordinates assay data (e.g., counts) with row (feature/gene)
and column (sample) metadata. Equivalent to R/Bioconductor's `SummarizedExperiment`.

Ensures that subsetting samples or features automatically maintains metadata alignment.

# Fields
- `assays::Dict{String, Matrix{<:Real}}`: Named matrices where rows = features, cols = samples.
- `rowData::Dict{Symbol, Vector}`: Metadata for each row (feature).
- `colData::Dict{Symbol, Vector}`: Metadata for each column (sample).
- `metadata::Dict{Symbol, Any}`: General experiment-level metadata.
"""
struct SummarizedExperiment{T<:Real}
    assays::Dict{String,Matrix{T}}
    rowData::Dict{Symbol,Vector}
    colData::Dict{Symbol,Vector}
    metadata::Dict{Symbol,Any}

    function SummarizedExperiment(
        assays::AbstractDict{String, <:AbstractMatrix{T}},
        rowData::AbstractDict{Symbol, <:AbstractVector},
        colData::AbstractDict{Symbol, <:AbstractVector},
        metadata::AbstractDict=Dict{Symbol,Any}()
    ) where {T<:Real}
        assays_conv = Dict{String, Matrix{T}}(k => Matrix{T}(v) for (k, v) in assays)
        rowData_conv = Dict{Symbol, Vector}(k => Vector(v) for (k, v) in rowData)
        colData_conv = Dict{Symbol, Vector}(k => Vector(v) for (k, v) in colData)
        metadata_conv = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in metadata)

        isempty(assays_conv) && throw(ArgumentError("at least one assay must be provided"))
        # Validate dimensions
        n_features, n_samples = size(first(values(assays_conv)))
        for (name, matrix) in assays_conv
            size(matrix) == (n_features, n_samples) || throw(DimensionMismatch("assay '$name' has dimensions $(size(matrix)), expected ($n_features, $n_samples)"))
        end
        for (key, vec) in rowData_conv
            length(vec) == n_features || throw(DimensionMismatch("rowData column '$key' has length $(length(vec)), expected $n_features"))
        end
        for (key, vec) in colData_conv
            length(vec) == n_samples || throw(DimensionMismatch("colData column '$key' has length $(length(vec)), expected $n_samples"))
        end
        return new{T}(assays_conv, rowData_conv, colData_conv, metadata_conv)
    end
end


"""
    SummarizedExperiment(counts::Matrix; rowData=..., colData=..., metadata=...)

Main constructor for `SummarizedExperiment`.
"""
function SummarizedExperiment(
    counts::Matrix{T};
    assay_name::String="counts",
    rowData::Dict{Symbol,Vector}=Dict{Symbol,Vector}(),
    colData::Dict{Symbol,Vector}=Dict{Symbol,Vector}(),
    metadata::Dict{Symbol,Any}=Dict{Symbol,Any}()
) where {T<:Real}
    return SummarizedExperiment(
        Dict(assay_name => counts),
        rowData,
        colData,
        metadata
    )
end

# Accessors
assay(se::SummarizedExperiment, name::String) = se.assays[name]
assay(se::SummarizedExperiment) = first(values(se.assays))
rowData(se::SummarizedExperiment) = se.rowData
colData(se::SummarizedExperiment) = se.colData
metadata(se::SummarizedExperiment) = se.metadata

# Basic subsetting
function subset_features(se::SummarizedExperiment, feature_indices)
    new_assays = Dict(name => matrix[feature_indices, :] for (name, matrix) in se.assays)
    new_rowData = Dict(key => vec[feature_indices] for (key, vec) in se.rowData)
    return SummarizedExperiment(new_assays, new_rowData, copy(se.colData), copy(se.metadata))
end

function subset_samples(se::SummarizedExperiment, sample_indices)
    new_assays = Dict(name => matrix[:, sample_indices] for (name, matrix) in se.assays)
    new_colData = Dict(key => vec[sample_indices] for (key, vec) in se.colData)
    return SummarizedExperiment(new_assays, copy(se.rowData), new_colData, copy(se.metadata))
end

function Base.show(io::IO, se::SummarizedExperiment)
    nf, ns = size(assay(se))
    na = length(se.assays)
    print(io, "SummarizedExperiment($nf features, $ns samples, $na assay(s))")
end

# ---- Interval Tree for Genomic Ranges ----------------------------------------

"""
    IntervalTreeNode{T}

Node in an augmented interval tree (centered interval tree / red-black tree
variant).  Each node stores a single interval and caches `max_end` — the
maximum right endpoint in its subtree — enabling O(log n + k) overlap queries.

Reference: Cormen, Leiserson, Rivest, Stein. "Introduction to Algorithms",
Chapter 14.3 — Augmenting Data Structures (Interval Trees).
"""
mutable struct IntervalTreeNode{T}
    left_endpoint::Int
    right_endpoint::Int
    max_end::Int         # max right_endpoint in this subtree
    payload::T
    left::Union{Nothing,IntervalTreeNode{T}}
    right::Union{Nothing,IntervalTreeNode{T}}
    height::Int          # AVL balance factor
end

"""
    IntervalTree{T}

Self-balancing (AVL) augmented interval tree supporting O(log n) insertion
and O(log n + k) overlap queries, where k is the number of results.

This replaces the sorted-array + binary-search approach that degraded to
O(n) for overlap queries on large collections.
"""
struct IntervalTree{T}
    root::Base.RefValue{Union{Nothing,IntervalTreeNode{T}}}
end

IntervalTree{T}() where {T} = IntervalTree{T}(Ref{Union{Nothing,IntervalTreeNode{T}}}(nothing))

function _itn_height(node::Union{Nothing,IntervalTreeNode})
    return node === nothing ? 0 : node.height
end

function _itn_max_end(node::Union{Nothing,IntervalTreeNode})
    return node === nothing ? typemin(Int) : node.max_end
end

function _itn_update!(node::IntervalTreeNode)
    node.height = 1 + max(_itn_height(node.left), _itn_height(node.right))
    node.max_end = max(node.right_endpoint,
        _itn_max_end(node.left),
        _itn_max_end(node.right))
    return node
end

function _itn_balance(node::IntervalTreeNode)
    return _itn_height(node.left) - _itn_height(node.right)
end

function _itn_rotate_right(y::IntervalTreeNode{T}) where {T}
    x = y.left::IntervalTreeNode{T}
    t2 = x.right
    x.right = y
    y.left = t2
    _itn_update!(y)
    _itn_update!(x)
    return x
end

function _itn_rotate_left(x::IntervalTreeNode{T}) where {T}
    y = x.right::IntervalTreeNode{T}
    t2 = y.left
    y.left = x
    x.right = t2
    _itn_update!(x)
    _itn_update!(y)
    return y
end

"""
    insert!(tree::IntervalTree, left, right, payload)

Insert an interval `[left, right]` with associated `payload` into the tree.
Maintains AVL balance and augmented `max_end` invariant.
"""
function Base.insert!(tree::IntervalTree{T}, left::Int, right::Int, payload::T) where {T}
    tree.root[] = _itn_insert(tree.root[], left, right, payload)
    return tree
end

function _itn_insert(node::Nothing, left::Int, right::Int, payload::T) where {T}
    return IntervalTreeNode{T}(left, right, right, payload, nothing, nothing, 1)
end

function _itn_insert(node::IntervalTreeNode{T}, left::Int, right::Int, payload::T) where {T}
    if left < node.left_endpoint || (left == node.left_endpoint && right < node.right_endpoint)
        node.left = _itn_insert(node.left, left, right, payload)
    else
        node.right = _itn_insert(node.right, left, right, payload)
    end

    _itn_update!(node)

    balance = _itn_balance(node)

    # Left-heavy
    if balance > 1
        if node.left !== nothing && _itn_balance(node.left) >= 0
            return _itn_rotate_right(node)
        elseif node.left !== nothing
            node.left = _itn_rotate_left(node.left)
            return _itn_rotate_right(node)
        end
    end

    # Right-heavy
    if balance < -1
        if node.right !== nothing && _itn_balance(node.right) <= 0
            return _itn_rotate_left(node)
        elseif node.right !== nothing
            node.right = _itn_rotate_right(node.right)
            return _itn_rotate_left(node)
        end
    end

    return node
end

"""
    query_overlaps(tree::IntervalTree, query_left, query_right) -> Vector

Return all payloads whose intervals overlap `[query_left, query_right]`.
Runs in O(log n + k) time where k is the number of overlapping intervals,
thanks to the augmented `max_end` field that prunes entire subtrees.
"""
function query_overlaps(tree::IntervalTree{T}, query_left::Int, query_right::Int) where {T}
    results = T[]
    _itn_query(tree.root[], query_left, query_right, results)
    return results
end

function _itn_query(node::Nothing, ::Int, ::Int, ::Vector)
    return nothing
end

function _itn_query(node::IntervalTreeNode{T}, ql::Int, qr::Int, results::Vector{T}) where {T}
    # Prune: if the maximum endpoint in this subtree is below our query start,
    # no interval in the subtree can overlap.
    node.max_end < ql && return nothing

    # Search left subtree first (may contain overlapping intervals)
    _itn_query(node.left, ql, qr, results)

    # Check current node
    if node.left_endpoint <= qr && node.right_endpoint >= ql
        push!(results, node.payload)
    end

    # Prune right subtree: if the current node starts after query end,
    # all nodes in the right subtree also start after query end.
    node.left_endpoint > qr && return nothing

    _itn_query(node.right, ql, qr, results)
    return nothing
end

"""
    Base.length(tree::IntervalTree) -> Int

Count the number of intervals stored in the tree.
"""
function Base.length(tree::IntervalTree)
    return _itn_count(tree.root[])
end

function _itn_count(node::Nothing)
    return 0
end

function _itn_count(node::IntervalTreeNode)
    return 1 + _itn_count(node.left) + _itn_count(node.right)
end

Base.isempty(tree::IntervalTree) = tree.root[] === nothing

# ---- Bioconductor Parity Types ----------------------------------------------

abstract type GeneIdType end
struct EnsemblID <: GeneIdType end
struct EntrezID <: GeneIdType end
struct SymbolID <: GeneIdType end
struct RefSeqID <: GeneIdType end
struct UniProtID <: GeneIdType end

const ensembl = EnsemblID()
const entrez = EntrezID()
const symbol = SymbolID()
const refseq = RefSeqID()
const uniprot = UniProtID()

function GeneIdType(s::Symbol)
    if s === :ensembl return EnsemblID()
    elseif s === :entrez return EntrezID()
    elseif s === :symbol return SymbolID()
    elseif s === :refseq return RefSeqID()
    elseif s === :uniprot return UniProtID()
    else throw(ArgumentError("Unknown GeneIdType symbol: $s"))
    end
end

Base.Symbol(::EnsemblID) = :ensembl
Base.Symbol(::EntrezID) = :entrez
Base.Symbol(::SymbolID) = :symbol
Base.Symbol(::RefSeqID) = :refseq
Base.Symbol(::UniProtID) = :uniprot

struct OrganismDb{Organism}
    metadata::Dict{Symbol,Any}
    mapper::Any
    term_indices::Dict{Symbol,Any}
end

function OrganismDb(organism::Symbol, metadata::Dict{Symbol,Any}=Dict{Symbol,Any}(), mapper=nothing, term_indices::Dict{Symbol,Any}=Dict{Symbol,Any}())
    return OrganismDb{organism}(metadata, mapper, term_indices)
end

struct TxQuantRecord
    tx_ids::Vector{String}
    counts::Vector{Float64}
    tpm::Vector{Float64}
    efflength::Vector{Float64}
    tool::String
    sample_id::String
end

struct Tx2GeneMap <: AbstractAnalysisResult
    map::Dict{String,String}
    gene_id_type::GeneIdType
    bias_flags::Dict{Symbol,Any}
    provenance::ResultProvenance
end

Tx2GeneMap(map, gene_id_type, bias_flags) = Tx2GeneMap(map, gene_id_type, bias_flags, provenance_record("Tx2GeneMap", "biotypes/Tx2GeneMap"))

struct AnnotatedHeatmapSpec
    row_annotations::Dict{Symbol,Vector}
    col_annotations::Dict{Symbol,Vector}
    row_splits::Vector
    col_splits::Vector
end

# Accept DataFrame by converting each column to a Dict entry
_df_to_dict(d::Dict) = d
function _df_to_dict(df)
    out = Dict{Symbol,Vector}()
    for n in propertynames(df)
        out[n] = Vector(df[!, n])
    end
    out
end

function AnnotatedHeatmapSpec(; row_annotations=Dict{Symbol,Vector}(), col_annotations=Dict{Symbol,Vector}(), row_splits=String[], col_splits=String[])
    return AnnotatedHeatmapSpec(_df_to_dict(row_annotations), _df_to_dict(col_annotations), row_splits, col_splits)
end
