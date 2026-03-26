"""
    PairwiseAlignmentResult

Result object returned by the pairwise alignment routines. It stores the aligned
left and right sequences, the final score, the number of exact matches, and the
identity fraction.
"""
struct PairwiseAlignmentResult
    left::String
    right::String
    score::Int
    matches::Int
    identity::Float64
end

"""
    AbstractPairwiseScoring

Abstract parent type for all pairwise alignment scoring models.
"""
abstract type AbstractPairwiseScoring end

"""
    LinearPairwiseScoring

Simple match/mismatch scoring model for pairwise alignment.
"""
struct LinearPairwiseScoring <: AbstractPairwiseScoring
    match::Int
    mismatch::Int
end

"""
    SubstitutionMatrix

General substitution-matrix container with alphabet lookup tables, score
storage, and a fallback score for unknown symbols.
"""
struct SubstitutionMatrix
    alphabet::Vector{UInt8}
    scores::Matrix{Int}
    lookup::Dict{UInt8,Int}
    lookup_table::Vector{Int}
    default::Int
end

"""
    MatrixPairwiseScoring

Pairwise scoring wrapper around a `SubstitutionMatrix`.
"""
struct MatrixPairwiseScoring <: AbstractPairwiseScoring
    matrix::SubstitutionMatrix
end

"""
    CodonSubstitutionMatrix

Substitution matrix specialized for codon tokens encoded as packed byte values.
"""
struct CodonSubstitutionMatrix
    alphabet::Vector{UInt8}
    scores::Matrix{Int}
    lookup::Dict{UInt8,Int}
    lookup_table::Vector{Int}
    default::Int
end

"""
    CodonMatrixPairwiseScoring

Pairwise scoring wrapper around a `CodonSubstitutionMatrix`.
"""
struct CodonMatrixPairwiseScoring <: AbstractPairwiseScoring
    matrix::CodonSubstitutionMatrix
end

const _STANDARD_SUBSTITUTION_MATRIX_TEXT = Dict{String,String}(
    "BLOSUM62" => raw"""
        A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
    A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
    R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
    N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
    D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
    C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
    Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
    E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
    G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
    H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
    I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
    L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
    K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
    M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
    F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
    P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
    S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
    T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
    W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
    Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
    V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
    B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
    Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
    X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
    * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
    """,
    "BLOSUM80" => raw"""
        A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
    A  7 -3 -3 -3 -1 -2 -2  0 -3 -3 -3 -1 -2 -4 -1  2  0 -5 -4 -1 -3 -2 -1 -8
    R -3  9 -1 -3 -6  1 -1 -4  0 -5 -4  3 -3 -5 -3 -2 -2 -5 -4 -4 -2  0 -2 -8
    N -3 -1  9  2 -5  0 -1 -1  1 -6 -6  0 -4 -6 -4  1  0 -7 -4 -5  5 -1 -2 -8
    D -3 -3  2 10 -7 -1  2 -3 -2 -7 -7 -2 -6 -6 -3 -1 -2 -8 -6 -6  6  1 -3 -8
    C -1 -6 -5 -7 13 -5 -7 -6 -7 -2 -3 -6 -3 -4 -6 -2 -2 -5 -5 -2 -6 -7 -4 -8
    Q -2  1  0 -1 -5  9  3 -4  1 -5 -4  2 -1 -5 -3 -1 -1 -4 -3 -4 -1  5 -2 -8
    E -2 -1 -1  2 -7  3  8 -4  0 -6 -6  1 -4 -6 -2 -1 -2 -6 -5 -4  1  6 -2 -8
    G  0 -4 -1 -3 -6 -4 -4  9 -4 -7 -7 -3 -5 -6 -5 -1 -3 -6 -6 -6 -2 -4 -3 -8
    H -3  0  1 -2 -7  1  0 -4 12 -6 -5 -1 -4 -2 -4 -2 -3 -4  3 -5 -1  0 -2 -8
    I -3 -5 -6 -7 -2 -5 -6 -7 -6  7  2 -5  2 -1 -5 -4 -2 -5 -3  4 -6 -6 -2 -8
    L -3 -4 -6 -7 -3 -4 -6 -7 -5  2  6 -4  3  0 -5 -4 -3 -4 -2  1 -7 -5 -2 -8
    K -1  3  0 -2 -6  2  1 -3 -1 -5 -4  8 -3 -5 -2 -1 -1 -6 -4 -4 -1  1 -2 -8
    M -2 -3 -4 -6 -3 -1 -4 -5 -4  2  3 -3  9  0 -4 -3 -1 -3 -3  1 -5 -3 -2 -8
    F -4 -5 -6 -6 -4 -5 -6 -6 -2 -1  0 -5  0 10 -6 -4 -4  0  4 -2 -6 -6 -3 -8
    P -1 -3 -4 -3 -6 -3 -2 -5 -4 -5 -5 -2 -4 -6 12 -2 -3 -7 -6 -4 -4 -2 -3 -8
    S  2 -2  1 -1 -2 -1 -1 -1 -2 -4 -4 -1 -3 -4 -2  7  2 -6 -3 -3  0 -1 -1 -8
    T  0 -2  0 -2 -2 -1 -2 -3 -3 -2 -3 -1 -1 -4 -3  2  8 -5 -3  0 -1 -2 -1 -8
    W -5 -5 -7 -8 -5 -4 -6 -6 -4 -5 -4 -6 -3  0 -7 -6 -5 16  3 -5 -8 -5 -5 -8
    Y -4 -4 -4 -6 -5 -3 -5 -6  3 -3 -2 -4 -3  4 -6 -3 -3  3 11 -3 -5 -4 -3 -8
    V -1 -4 -5 -6 -2 -4 -4 -6 -5  4  1 -4  1 -2 -4 -3  0 -5 -3  7 -6 -4 -2 -8
    B -3 -2  5  6 -6 -1  1 -2 -1 -6 -7 -1 -5 -6 -4  0 -1 -8 -5 -6  6  0 -3 -8
    Z -2  0 -1  1 -7  5  6 -4  0 -6 -5  1 -3 -6 -2 -1 -2 -5 -4 -4  0  6 -1 -8
    X -1 -2 -2 -3 -4 -2 -2 -3 -2 -2 -2 -2 -2 -3 -3 -1 -1 -5 -3 -2 -3 -1 -2 -8
    * -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8  1
    """,
    "PAM250" => raw"""
        A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
    A  2 -2  0  0 -2  0  0  1 -1 -1 -2 -1 -1 -3  1  1  1 -6 -3  0  0  0  0 -8
    R -2  6  0 -1 -4  1 -1 -3  2 -2 -3  3  0 -4  0  0 -1  2 -4 -2 -1  0 -1 -8
    N  0  0  2  2 -4  1  1  0  2 -2 -3  1 -2 -3  0  1  0 -4 -2 -2  2  1  0 -8
    D  0 -1  2  4 -5  2  3  1  1 -2 -4  0 -3 -6 -1  0  0 -7 -4 -2  3  3 -1 -8
    C -2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3  0 -2 -8  0 -2 -4 -5 -3 -8
    Q  0  1  1  2 -5  4  2 -1  3 -2 -2  1 -1 -5  0 -1 -1 -5 -4 -2  1  3 -1 -8
    E  0 -1  1  3 -5  2  4  0  1 -2 -3  0 -2 -5 -1  0  0 -7 -4 -2  3  3 -1 -8
    G  1 -3  0  1 -3 -1  0  5 -2 -3 -4 -2 -3 -5  0  1  0 -7 -5 -1  0  0 -1 -8
    H -1  2  2  1 -3  3  1 -2  6 -2 -2  0 -2 -2  0 -1 -1 -3  0 -2  1  2 -1 -8
    I -1 -2 -2 -2 -2 -2 -2 -3 -2  5  2 -2  2  1 -2 -1  0 -5 -1  4 -2 -2 -1 -8
    L -2 -3 -3 -4 -6 -2 -3 -4 -2  2  6 -3  4  2 -3 -3 -2 -2 -1  2 -3 -3 -1 -8
    K -1  3  1  0 -5  1  0 -2  0 -2 -3  5  0 -5 -1  0  0 -3 -4 -2  1  0 -1 -8
    M -1  0 -2 -3 -5 -1 -2 -3 -2  2  4  0  6  0 -2 -2 -1 -4 -2  2 -2 -2 -1 -8
    F -3 -4 -3 -6 -4 -5 -5 -5 -2  1  2 -5  0  9 -5 -3 -3  0  7 -1 -4 -5 -2 -8
    P  1  0  0 -1 -3  0 -1  0  0 -2 -3 -1 -2 -5  6  1  0 -6 -5 -1 -1  0 -1 -8
    S  1  0  1  0  0 -1  0  1 -1 -1 -3  0 -2 -3  1  2  1 -2 -3 -1  0  0  0 -8
    T  1 -1  0  0 -2 -1  0  0 -1  0 -2  0 -1 -3  0  1  3 -5 -3  0  0 -1  0 -8
    W -6  2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4  0 -6 -2 -5 17  0 -6 -5 -6 -4 -8
    Y -3 -4 -2 -4  0 -4 -4 -5  0 -1 -1 -4 -2  7 -5 -3 -3  0 10 -2 -3 -4 -2 -8
    V  0 -2 -2 -2 -2 -2 -2 -1 -2  4  2 -2  2 -1 -1 -1  0 -6 -2  4 -2 -2 -1 -8
    B  0 -1  2  3 -4  1  3  0  1 -2 -3  1 -2 -4 -1  0  0 -5 -3 -2  3  2 -1 -8
    Z  0  0  1  3 -5  3  3  0  2 -2 -3  0 -2 -5  0  0 -1 -6 -4 -2  2  3 -1 -8
    X  0 -1  0 -1 -3 -1 -1 -1 -1 -1 -1 -1 -1 -2 -1  0  0 -4 -2 -1 -1 -1 -1 -8
    * -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8  1
    """
)

include("substitution_matrices_data.jl")
merge!(_STANDARD_SUBSTITUTION_MATRIX_TEXT, _BIOPYTHON_NAMED_SUBSTITUTION_MATRIX_TEXT)

const _STANDARD_SUBSTITUTION_MATRIX_CACHE = Dict{String,SubstitutionMatrix}()

"""
    _normalize_substitution_matrix_name(name)

Normalize a matrix name by uppercasing it and removing spacing and separator
characters.
"""
@inline function _normalize_substitution_matrix_name(name::AbstractString)
    return replace(uppercase(strip(name)), r"[\s._-]+" => "")
end

"""
    _parse_standard_substitution_matrix(spec; threaded=true)

Parse a standard amino-acid substitution matrix from its textual representation.
"""
function _parse_standard_substitution_matrix(spec::AbstractString; threaded::Bool=true)
    lines = filter(!isempty, Base.split(chomp(spec), '\n'))
    length(lines) >= 2 || throw(ArgumentError("substitution matrix spec must include a header and at least one row"))

    header = Base.split(strip(lines[1]))
    alphabet = join(header)
    row_count = length(header)
    scores = Matrix{Int}(undef, row_count, row_count)

    if threaded && row_count > 1 && Threads.nthreads() > 1
        Threads.@threads for row_index in 1:row_count
            fields = Base.split(strip(lines[row_index + 1]))
            length(fields) == row_count + 1 || throw(ArgumentError("substitution matrix row has the wrong width"))
            fields[1] == header[row_index] || throw(ArgumentError("substitution matrix row order does not match its header"))
            @inbounds for col_index in 1:row_count
                scores[row_index, col_index] = parse(Int, fields[col_index + 1])
            end
        end
    else
        for (row_index, line) in enumerate(@view lines[2:end])
            fields = Base.split(strip(line))
            length(fields) == row_count + 1 || throw(ArgumentError("substitution matrix row has the wrong width"))
            fields[1] == header[row_index] || throw(ArgumentError("substitution matrix row order does not match its header"))
            @inbounds for col_index in 1:row_count
                scores[row_index, col_index] = parse(Int, fields[col_index + 1])
            end
        end
    end

    return SubstitutionMatrix(alphabet, scores; default=minimum(scores))
end

const _CODON_DECODE_TABLE = let table = Vector{String}(undef, 64)
    for first in 0:3
        for second in 0:3
            for third in 0:3
                code = (first << 4) | (second << 2) | third
                table[code + 1] = String(UInt8[_kmer_base_from_code(first), _kmer_base_from_code(second), _kmer_base_from_code(third)])
            end
        end
    end
    table
end

"""
    _codon_code(byte)

Map a nucleotide byte to its compact DNA code.
"""
@inline function _codon_code(byte::UInt8)
    return _DNA_CODE[Int(byte) + 1]
end

"""
    _pack_codon_token(first, second, third)

Pack three nucleotide bytes into a single codon token.
"""
@inline function _pack_codon_token(first::UInt8, second::UInt8, third::UInt8)
    code1 = _codon_code(first)
    code2 = _codon_code(second)
    code3 = _codon_code(third)
    return code1 > 3 || code2 > 3 || code3 > 3 ? UInt8(255) : UInt8((Int(code1) << 4) | (Int(code2) << 2) | Int(code3))
end

"""
    _encode_codon_sequence(bytes)

Convert a nucleotide byte sequence into packed codon tokens.
"""
function _encode_codon_sequence(bytes::AbstractVector{UInt8})
    length(bytes) % 3 == 0 || throw(ArgumentError("codon sequence length must be divisible by 3"))
    token_count = length(bytes) ÷ 3
    tokens = Vector{UInt8}(undef, token_count)

    @inbounds for token_index in 1:token_count
        byte_index = (token_index - 1) * 3 + 1
        tokens[token_index] = _pack_codon_token(bytes[byte_index], bytes[byte_index + 1], bytes[byte_index + 2])
    end

    return tokens
end

"""
    _decode_codon_token(token)

Convert a packed codon token back into a three-character codon string.
"""
function _decode_codon_token(token::UInt8)
    return token > 63 ? "NNN" : _CODON_DECODE_TABLE[Int(token) + 1]
end

"""
    _parse_codon_substitution_matrix(spec; threaded=true)

Parse a codon substitution matrix from its textual table form.
"""
function _parse_codon_substitution_matrix(spec::AbstractString; threaded::Bool=true)
    lines = filter(!isempty, Base.split(chomp(spec), '\n'))
    length(lines) >= 2 || throw(ArgumentError("substitution matrix spec must include a header and at least one row"))

    header = Base.split(strip(lines[1]))
    length(header) % 3 == 0 || throw(ArgumentError("codon substitution matrix header must contain codon triplets"))
    row_count = length(header) ÷ 3
    row_count == 64 || throw(ArgumentError("codon substitution matrix must contain exactly 64 codons"))

    alphabet = Vector{UInt8}(undef, row_count)
    lookup_table = fill(0, 64)
    lookup = Dict{UInt8,Int}()

    for index in 1:row_count
        token_start = (index - 1) * 3 + 1
        token1 = header[token_start]
        token2 = header[token_start + 1]
        token3 = header[token_start + 2]
        token = join((token1, token2, token3), " ")
        code = _pack_codon_token(UInt8(token1[1]), UInt8(token2[1]), UInt8(token3[1]))
        code == 255 && throw(ArgumentError("invalid codon token in matrix header: $token"))
        alphabet[index] = code
        lookup[code] = index
        lookup_table[Int(code) + 1] = index
    end

    scores = Matrix{Int}(undef, row_count, row_count)

    if threaded && row_count > 1 && Threads.nthreads() > 1
        Threads.@threads for row_index in 1:row_count
            fields = Base.split(strip(lines[row_index + 1]))
            length(fields) == row_count + 1 || throw(ArgumentError("substitution matrix row has the wrong width"))
            @inbounds for col_index in 1:row_count
                scores[row_index, col_index] = parse(Int, fields[col_index + 1])
            end
        end
    else
        for (row_index, line) in enumerate(@view lines[2:end])
            fields = Base.split(strip(line))
            length(fields) == row_count + 1 || throw(ArgumentError("substitution matrix row has the wrong width"))
            @inbounds for col_index in 1:row_count
                scores[row_index, col_index] = parse(Int, fields[col_index + 1])
            end
        end
    end

    return CodonSubstitutionMatrix(alphabet, scores, lookup, lookup_table, minimum(scores))
end

const _PAIRWISE_STATE_MATCH = UInt8(0x01)
const _PAIRWISE_STATE_GAP_LEFT = UInt8(0x02)
const _PAIRWISE_STATE_GAP_RIGHT = UInt8(0x03)
const _PAIRWISE_TRACE_NONE = UInt8(0x00)

"""
    Base.show(io, result)

Render a compact summary of a `PairwiseAlignmentResult`.
"""
function Base.show(io::IO, result::PairwiseAlignmentResult)
    print(
        io,
        "PairwiseAlignmentResult(score=",
        result.score,
        ", matches=",
        result.matches,
        ", identity=",
        round(result.identity, digits=4),
        ")",
    )
end

"""
    SubstitutionMatrix(alphabet; match=1, mismatch=-1, default=mismatch, threaded=true)

Build a simple square substitution matrix from an alphabet string.
"""
function SubstitutionMatrix(alphabet::AbstractString; match::Int=1, mismatch::Int=-1, default::Int=mismatch, threaded::Bool=true)
    alphabet_bytes = collect(codeunits(alphabet))
    scores = Matrix{Int}(undef, length(alphabet_bytes), length(alphabet_bytes))
    lookup_table = fill(0, 256)

    if threaded && length(alphabet_bytes) > 1 && Threads.nthreads() > 1
        Threads.@threads for i in eachindex(alphabet_bytes)
            lookup_table[Int(alphabet_bytes[i]) + 1] = i
            @inbounds for j in eachindex(alphabet_bytes)
                scores[i, j] = alphabet_bytes[i] == alphabet_bytes[j] ? match : mismatch
            end
        end
    else
        @inbounds for i in eachindex(alphabet_bytes)
            lookup_table[Int(alphabet_bytes[i]) + 1] = i
            for j in eachindex(alphabet_bytes)
                scores[i, j] = alphabet_bytes[i] == alphabet_bytes[j] ? match : mismatch
            end
        end
    end

    lookup = Dict{UInt8,Int}(byte => index for (index, byte) in pairs(alphabet_bytes))
    return SubstitutionMatrix(alphabet_bytes, scores, lookup, lookup_table, default)
end

"""
    SubstitutionMatrix(alphabet, scores; default=0, threaded=true)

Wrap an explicit square score matrix for the given alphabet.
"""
function SubstitutionMatrix(alphabet::AbstractString, scores::AbstractMatrix{<:Integer}; default::Int=0, threaded::Bool=true)
    alphabet_bytes = collect(codeunits(alphabet))
    size(scores, 1) == length(alphabet_bytes) || throw(ArgumentError("score matrix row count must match alphabet length"))
    size(scores, 2) == length(alphabet_bytes) || throw(ArgumentError("score matrix column count must match alphabet length"))
    lookup = Dict{UInt8,Int}(byte => index for (index, byte) in pairs(alphabet_bytes))
    lookup_table = fill(0, 256)
    if threaded && length(alphabet_bytes) > 1 && Threads.nthreads() > 1
        Threads.@threads for index in eachindex(alphabet_bytes)
            byte = alphabet_bytes[index]
            lookup_table[Int(byte) + 1] = index
        end
    else
        @inbounds for (index, byte) in pairs(alphabet_bytes)
            lookup_table[Int(byte) + 1] = index
        end
    end
    return SubstitutionMatrix(alphabet_bytes, Matrix{Int}(scores), lookup, lookup_table, default)
end

"""
    substitution_matrix(alphabet; kwargs...)

Convenience alias for `SubstitutionMatrix(alphabet; kwargs...)`.
"""
substitution_matrix(alphabet::AbstractString; kwargs...) = SubstitutionMatrix(alphabet; kwargs...)

"""
    substitution_matrix(alphabet, scores; kwargs...)

Convenience alias for `SubstitutionMatrix(alphabet, scores; kwargs...)`.
"""
substitution_matrix(alphabet::AbstractString, scores::AbstractMatrix{<:Integer}; kwargs...) = SubstitutionMatrix(alphabet, scores; kwargs...)

"""
    named_substitution_matrix(name)

Load one of the built-in named substitution matrices or a DNA identity matrix.
"""
function named_substitution_matrix(name::AbstractString)
    normalized = _normalize_substitution_matrix_name(name)
    if normalized == "DNA" || normalized == "NUC" || normalized == "NUCLEOTIDE" || normalized == "IDENTITY"
        return SubstitutionMatrix("ACGT"; match=1, mismatch=-1, default=-1)
    end

    matrix = get(_STANDARD_SUBSTITUTION_MATRIX_CACHE, normalized, nothing)
    if matrix !== nothing
        return matrix
    end

    spec = get(_STANDARD_SUBSTITUTION_MATRIX_TEXT, String(name), nothing)
    spec === nothing && (spec = get(_STANDARD_SUBSTITUTION_MATRIX_TEXT, normalized, nothing))
    spec === nothing && throw(ArgumentError("unknown named substitution matrix: $name"))

    parsed = _parse_standard_substitution_matrix(spec)
    _STANDARD_SUBSTITUTION_MATRIX_CACHE[normalized] = parsed
    return parsed
end

"""
    named_substitution_matrix(name::Symbol)

Symbol-based wrapper for `named_substitution_matrix`.
"""
named_substitution_matrix(name::Symbol) = named_substitution_matrix(String(name))

"""
    substitution_matrix(name::Symbol)

Symbol-based alias for `named_substitution_matrix`.
"""
substitution_matrix(name::Symbol) = named_substitution_matrix(name)

"""
    available_named_substitution_matrices()

List the built-in named substitution matrices available in this module.
"""
available_named_substitution_matrices() = sort!(collect(keys(_STANDARD_SUBSTITUTION_MATRIX_TEXT)))

"""
    named_codon_substitution_matrix(name)

Load a built-in codon substitution matrix by name.
"""
function named_codon_substitution_matrix(name::AbstractString)
    normalized = _normalize_substitution_matrix_name(name)
    normalized == "SCHNEIDER" || throw(ArgumentError("unknown named codon substitution matrix: $name"))

    spec = get(_STANDARD_SUBSTITUTION_MATRIX_TEXT, normalized, nothing)
    spec === nothing && throw(ArgumentError("unknown named codon substitution matrix: $name"))

    return _parse_codon_substitution_matrix(spec)
end

"""
    named_codon_substitution_matrix(name::Symbol)

Symbol-based wrapper for `named_codon_substitution_matrix`.
"""
named_codon_substitution_matrix(name::Symbol) = named_codon_substitution_matrix(String(name))

"""
    codon_substitution_matrix(alphabet, scores; default=0, threaded=true)

Construct a codon-aware substitution matrix from an explicit codon alphabet and
score matrix.
"""
function codon_substitution_matrix(alphabet::AbstractString, scores::AbstractMatrix{<:Integer}; default::Int=0, threaded::Bool=true)
    codon_tokens = Base.split(strip(alphabet))
    length(codon_tokens) == 64 || throw(ArgumentError("codon alphabet must contain exactly 64 codons"))

    alphabet_codes = Vector{UInt8}(undef, 64)
    lookup_table = fill(0, 64)

    if threaded && length(codon_tokens) > 1 && Threads.nthreads() > 1
        Threads.@threads for index in eachindex(codon_tokens)
            token = codon_tokens[index]
            length(token) == 3 || throw(ArgumentError("codon alphabet entries must be triplets"))
            code = _pack_codon_token(UInt8(token[1]), UInt8(token[2]), UInt8(token[3]))
            code == 255 && throw(ArgumentError("invalid codon token in alphabet: $token"))
            alphabet_codes[index] = code
            lookup_table[Int(code) + 1] = index
        end
    else
        for (index, token) in enumerate(codon_tokens)
            length(token) == 3 || throw(ArgumentError("codon alphabet entries must be triplets"))
            code = _pack_codon_token(UInt8(token[1]), UInt8(token[2]), UInt8(token[3]))
            code == 255 && throw(ArgumentError("invalid codon token in alphabet: $token"))
            alphabet_codes[index] = code
            lookup_table[Int(code) + 1] = index
        end
    end

    lookup = Dict{UInt8,Int}(code => index for (index, code) in pairs(alphabet_codes))

    size(scores, 1) == 64 || throw(ArgumentError("codon score matrix row count must be 64"))
    size(scores, 2) == 64 || throw(ArgumentError("codon score matrix column count must be 64"))
    return CodonSubstitutionMatrix(alphabet_codes, Matrix{Int}(scores), lookup, lookup_table, default)
end

"""
    _pairwise_score(scoring, left_byte, right_byte)

Score a nucleotide pair under linear match/mismatch scoring.
"""
@inline function _pairwise_score(scoring::LinearPairwiseScoring, left_byte::UInt8, right_byte::UInt8)
    return left_byte == right_byte ? scoring.match : scoring.mismatch
end

"""
    _pairwise_score(scoring, left_byte, right_byte)

Score a nucleotide pair using a substitution matrix lookup.
"""
@inline function _pairwise_score(scoring::MatrixPairwiseScoring, left_byte::UInt8, right_byte::UInt8)
    left_index = scoring.matrix.lookup_table[Int(left_byte) + 1]
    right_index = scoring.matrix.lookup_table[Int(right_byte) + 1]
    left_index == 0 || right_index == 0 ? scoring.matrix.default : scoring.matrix.scores[left_index, right_index]
end

"""
    _pairwise_score(scoring, left_token, right_token)

Score a codon pair using a codon substitution matrix lookup.
"""
@inline function _pairwise_score(scoring::CodonMatrixPairwiseScoring, left_token::UInt8, right_token::UInt8)
    left_index = left_token > 63 ? 0 : scoring.matrix.lookup_table[Int(left_token) + 1]
    right_index = right_token > 63 ? 0 : scoring.matrix.lookup_table[Int(right_token) + 1]
    left_index == 0 || right_index == 0 ? scoring.matrix.default : scoring.matrix.scores[left_index, right_index]
end

"""
    _pairwise_gap_parameters(gap, gap_open, gap_extend)

Resolve linear and affine gap penalties into explicit open and extend values.
"""
@inline function _pairwise_gap_parameters(gap::Int, gap_open::Union{Nothing,Int}, gap_extend::Union{Nothing,Int})
    if gap_open === nothing && gap_extend === nothing
        return gap, gap
    end

    open_score = gap_open === nothing ? something(gap_extend, gap) : gap_open
    extend_score = gap_extend === nothing ? something(gap_open, gap) : gap_extend
    return open_score, extend_score
end

"""
    pairwise_align(left_bytes, right_bytes; kwargs...)

Align two byte sequences using linear or affine gap scoring.
"""
function pairwise_align(
    left_bytes::AbstractVector{UInt8},
    right_bytes::AbstractVector{UInt8};
    match::Int=1,
    mismatch::Int=-1,
    substitution_matrix::Union{Nothing,SubstitutionMatrix}=nothing,
    gap::Int=-1,
    gap_open::Union{Nothing,Int}=nothing,
    gap_extend::Union{Nothing,Int}=nothing,
    is_local::Bool=false,
)
    scoring = substitution_matrix === nothing ? LinearPairwiseScoring(match, mismatch) : MatrixPairwiseScoring(substitution_matrix)
    if gap_open === nothing && gap_extend === nothing
        return is_local ? _pairwise_align_local(left_bytes, right_bytes, scoring, gap) : _pairwise_align_global(left_bytes, right_bytes, scoring, gap)
    end

    affine_gap_open, affine_gap_extend = _pairwise_gap_parameters(gap, gap_open, gap_extend)
    return is_local ? _pairwise_align_affine_local(left_bytes, right_bytes, scoring, affine_gap_open, affine_gap_extend) : _pairwise_align_affine_global(left_bytes, right_bytes, scoring, affine_gap_open, affine_gap_extend)
end

"""
    pairwise_align(left, right; kwargs...)

String-based convenience wrapper for `pairwise_align`.
"""
function pairwise_align(
    left::AbstractString,
    right::AbstractString;
    kwargs...
)
    # Using simple Arrays to avoid SubArray performance hits entirely inside DP
    # For small sequences like local_search envelopes, direct views might be fine
    # but the inner kernels accept AbstractVector{UInt8}. 
    return pairwise_align(codeunits(left), codeunits(right); kwargs...)
end

"""
    pairwise_align(left::SeqRecordLite, right::SeqRecordLite; kwargs...)

Align the sequence payloads stored in two `SeqRecordLite` values.
"""
function pairwise_align(left::SeqRecordLite, right::SeqRecordLite; kwargs...)
    return pairwise_align(left.sequence, right.sequence; kwargs...)
end

"""
    pairwise_align(left::FastqRecord, right::FastqRecord; kwargs...)

Align the sequence payloads stored in two FASTQ records.
"""
function pairwise_align(left::FastqRecord, right::FastqRecord; kwargs...)
    return pairwise_align(left.sequence, right.sequence; kwargs...)
end

"""
    _write_codon_token!(buffer, write_index, codon)

Write a three-character codon into a preallocated byte buffer from the back.
"""
@inline function _write_codon_token!(buffer::Vector{UInt8}, write_index::Int, codon::AbstractString)
    codon_bytes = codeunits(codon)
    buffer[write_index - 2] = codon_bytes[1]
    buffer[write_index - 1] = codon_bytes[2]
    buffer[write_index] = codon_bytes[3]
    return write_index - 3
end

"""
    _pairwise_traceback_codon(...)

Reconstruct a codon alignment from a standard traceback matrix.
"""
function _pairwise_traceback_codon(
    left_tokens::AbstractVector{UInt8},
    right_tokens::AbstractVector{UInt8},
    scores::Union{Nothing,Matrix{Int}},
    trace::Matrix{UInt8},
    i::Int,
    j::Int,
    best_score::Int,
)
    aligned_left = Vector{UInt8}(undef, 3 * (length(left_tokens) + length(right_tokens)))
    aligned_right = Vector{UInt8}(undef, 3 * (length(left_tokens) + length(right_tokens)))
    left_write_index = length(aligned_left)
    right_write_index = length(aligned_right)
    aligned_length = 0
    matches = 0

    current_i = i
    current_j = j

    while current_i > 1 || current_j > 1
        if scores !== nothing && scores[current_i, current_j] == 0 && trace[current_i, current_j] == 0x00
            break
        end

        direction = trace[current_i, current_j]
        direction == 0x00 && break

        if direction == 0x01
            left_token = left_tokens[current_i - 1]
            right_token = right_tokens[current_j - 1]
            left_write_index = _write_codon_token!(aligned_left, left_write_index, _decode_codon_token(left_token))
            right_write_index = _write_codon_token!(aligned_right, right_write_index, _decode_codon_token(right_token))
            matches += left_token == right_token ? 1 : 0
            current_i -= 1
            current_j -= 1
        elseif direction == 0x02
            left_write_index = _write_codon_token!(aligned_left, left_write_index, _decode_codon_token(left_tokens[current_i - 1]))
            right_write_index = _write_codon_token!(aligned_right, right_write_index, "---")
            current_i -= 1
        else
            left_write_index = _write_codon_token!(aligned_left, left_write_index, "---")
            right_write_index = _write_codon_token!(aligned_right, right_write_index, _decode_codon_token(right_tokens[current_j - 1]))
            current_j -= 1
        end

        aligned_length += 1
    end

    left_start_index = left_write_index + 1
    right_start_index = right_write_index + 1
    left_string = String(aligned_left[left_start_index:end])
    right_string = String(aligned_right[right_start_index:end])
    identity = aligned_length == 0 ? 0.0 : matches / aligned_length

    return PairwiseAlignmentResult(left_string, right_string, best_score, matches, identity)
end

"""
    _pairwise_traceback_codon_affine(...)

Reconstruct a codon alignment from affine-gap traceback matrices.
"""
function _pairwise_traceback_codon_affine(
    left_tokens::AbstractVector{UInt8},
    right_tokens::AbstractVector{UInt8},
    match_trace::Matrix{UInt8},
    gap_left_trace::Matrix{UInt8},
    gap_right_trace::Matrix{UInt8},
    i::Int,
    j::Int,
    start_state::UInt8,
    best_score::Int,
)
    aligned_left = Vector{UInt8}(undef, 3 * (length(left_tokens) + length(right_tokens)))
    aligned_right = Vector{UInt8}(undef, 3 * (length(left_tokens) + length(right_tokens)))
    left_write_index = length(aligned_left)
    right_write_index = length(aligned_right)
    aligned_length = 0
    matches = 0

    current_i = i
    current_j = j
    current_state = start_state

    while current_state != _PAIRWISE_TRACE_NONE && (current_i > 1 || current_j > 1)
        if current_state == _PAIRWISE_STATE_MATCH
            left_token = left_tokens[current_i - 1]
            right_token = right_tokens[current_j - 1]
            left_write_index = _write_codon_token!(aligned_left, left_write_index, _decode_codon_token(left_token))
            right_write_index = _write_codon_token!(aligned_right, right_write_index, _decode_codon_token(right_token))
            matches += left_token == right_token ? 1 : 0
            current_state = match_trace[current_i, current_j]
            current_i -= 1
            current_j -= 1
        elseif current_state == _PAIRWISE_STATE_GAP_LEFT
            left_write_index = _write_codon_token!(aligned_left, left_write_index, _decode_codon_token(left_tokens[current_i - 1]))
            right_write_index = _write_codon_token!(aligned_right, right_write_index, "---")
            current_state = gap_left_trace[current_i, current_j]
            current_i -= 1
        else
            left_write_index = _write_codon_token!(aligned_left, left_write_index, "---")
            right_write_index = _write_codon_token!(aligned_right, right_write_index, _decode_codon_token(right_tokens[current_j - 1]))
            current_state = gap_right_trace[current_i, current_j]
            current_j -= 1
        end

        aligned_length += 1
    end

    left_start_index = left_write_index + 1
    right_start_index = right_write_index + 1
    left_string = String(aligned_left[left_start_index:end])
    right_string = String(aligned_right[right_start_index:end])
    identity = aligned_length == 0 ? 0.0 : matches / aligned_length

    return PairwiseAlignmentResult(left_string, right_string, best_score, matches, identity)
end

"""
    _pairwise_align_codon_global(left_tokens, right_tokens, scoring, gap)

Run global codon alignment with a linear gap model.
"""
function _pairwise_align_codon_global(left_tokens::AbstractVector{UInt8}, right_tokens::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap::Int)
    left_length = length(left_tokens)
    right_length = length(right_tokens)

    previous_scores = Vector{Int}(undef, left_length + 1)
    current_scores = similar(previous_scores)
    trace = Matrix{UInt8}(undef, left_length + 1, right_length + 1)

    @inbounds begin
        previous_scores[1] = 0
        trace[1, 1] = 0x00

        for i in 2:left_length + 1
            previous_scores[i] = previous_scores[i - 1] + gap
            trace[i, 1] = 0x02
        end

        for j in 2:right_length + 1
            current_scores[1] = previous_scores[1] + gap
            trace[1, j] = 0x03
            right_token = right_tokens[j - 1]
            for i in 2:left_length + 1
                left_token = left_tokens[i - 1]

                diag_score = previous_scores[i - 1] + _pairwise_score(scoring, left_token, right_token)
                up_score = previous_scores[i] + gap
                left_score = current_scores[i - 1] + gap

                cell_score = diag_score
                cell_trace = 0x01

                if up_score > cell_score
                    cell_score = up_score
                    cell_trace = 0x02
                end
                if left_score > cell_score
                    cell_score = left_score
                    cell_trace = 0x03
                end

                current_scores[i] = cell_score
                trace[i, j] = cell_trace
            end

            previous_scores, current_scores = current_scores, previous_scores
        end
    end

    return _pairwise_traceback_codon(left_tokens, right_tokens, nothing, trace, left_length + 1, right_length + 1, previous_scores[left_length + 1])
end

"""
    _pairwise_align_codon_local(left_tokens, right_tokens, scoring, gap)

Run local codon alignment with a linear gap model.
"""
function _pairwise_align_codon_local(left_tokens::AbstractVector{UInt8}, right_tokens::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap::Int)
    left_length = length(left_tokens)
    right_length = length(right_tokens)

    scores = Matrix{Int}(undef, left_length + 1, right_length + 1)
    trace = Matrix{UInt8}(undef, left_length + 1, right_length + 1)

    best_score = 0
    best_i = 1
    best_j = 1

    @inbounds begin
        scores[1, 1] = 0
        trace[1, 1] = 0x00

        for i in 2:left_length + 1
            scores[i, 1] = 0
            trace[i, 1] = 0x00
        end

        for j in 2:right_length + 1
            scores[1, j] = 0
            trace[1, j] = 0x00
        end

        for j in 2:right_length + 1
            right_token = right_tokens[j - 1]
            for i in 2:left_length + 1
                left_token = left_tokens[i - 1]

                diag_score = scores[i - 1, j - 1] + _pairwise_score(scoring, left_token, right_token)
                up_score = scores[i - 1, j] + gap
                left_score = scores[i, j - 1] + gap

                cell_score = 0
                cell_trace = 0x00

                if diag_score > cell_score
                    cell_score = diag_score
                    cell_trace = 0x01
                end
                if up_score > cell_score
                    cell_score = up_score
                    cell_trace = 0x02
                end
                if left_score > cell_score
                    cell_score = left_score
                    cell_trace = 0x03
                end

                scores[i, j] = cell_score
                trace[i, j] = cell_trace

                if cell_score > best_score
                    best_score = cell_score
                    best_i = i
                    best_j = j
                end
            end
        end
    end

    return _pairwise_traceback_codon(left_tokens, right_tokens, scores, trace, best_i, best_j, best_score)
end

"""
    _pairwise_align_codon_affine_global(left_tokens, right_tokens, scoring, gap_open, gap_extend)

Run global codon alignment with an affine gap model.
"""
function _pairwise_align_codon_affine_global(left_tokens::AbstractVector{UInt8}, right_tokens::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap_open::Int, gap_extend::Int)
    left_length = length(left_tokens)
    right_length = length(right_tokens)
    negative_infinity = typemin(Int) ÷ 4

    previous_match_scores = fill(negative_infinity, left_length + 1)
    previous_gap_left_scores = fill(negative_infinity, left_length + 1)
    previous_gap_right_scores = fill(negative_infinity, left_length + 1)
    current_match_scores = similar(previous_match_scores)
    current_gap_left_scores = similar(previous_gap_left_scores)
    current_gap_right_scores = similar(previous_gap_right_scores)

    match_trace = zeros(UInt8, left_length + 1, right_length + 1)
    gap_left_trace = zeros(UInt8, left_length + 1, right_length + 1)
    gap_right_trace = zeros(UInt8, left_length + 1, right_length + 1)

    previous_match_scores[1] = 0
    previous_gap_left_scores[1] = negative_infinity
    previous_gap_right_scores[1] = negative_infinity

    @inbounds begin
        for i in 2:left_length + 1
            from_match = previous_match_scores[i - 1] + gap_open
            from_gap = previous_gap_left_scores[i - 1] + gap_extend
            if from_match >= from_gap
                previous_gap_left_scores[i] = from_match
                gap_left_trace[i, 1] = _PAIRWISE_STATE_MATCH
            else
                previous_gap_left_scores[i] = from_gap
                gap_left_trace[i, 1] = _PAIRWISE_STATE_GAP_LEFT
            end
        end

        for j in 2:right_length + 1
            current_match_scores[1] = negative_infinity
            current_gap_left_scores[1] = negative_infinity

            from_match = previous_match_scores[1] + gap_open
            from_gap = previous_gap_right_scores[1] + gap_extend
            if from_match >= from_gap
                current_gap_right_scores[1] = from_match
                gap_right_trace[1, j] = _PAIRWISE_STATE_MATCH
            else
                current_gap_right_scores[1] = from_gap
                gap_right_trace[1, j] = _PAIRWISE_STATE_GAP_RIGHT
            end

            right_token = right_tokens[j - 1]
            for i in 2:left_length + 1
                left_token = left_tokens[i - 1]

                best_prev = previous_match_scores[i - 1]
                best_state = _PAIRWISE_STATE_MATCH
                if previous_gap_left_scores[i - 1] > best_prev
                    best_prev = previous_gap_left_scores[i - 1]
                    best_state = _PAIRWISE_STATE_GAP_LEFT
                end
                if previous_gap_right_scores[i - 1] > best_prev
                    best_prev = previous_gap_right_scores[i - 1]
                    best_state = _PAIRWISE_STATE_GAP_RIGHT
                end
                current_match_scores[i] = best_prev + _pairwise_score(scoring, left_token, right_token)
                match_trace[i, j] = best_state

                from_match = current_match_scores[i - 1] + gap_open
                from_gap = current_gap_left_scores[i - 1] + gap_extend
                if from_match >= from_gap
                    current_gap_left_scores[i] = from_match
                    gap_left_trace[i, j] = _PAIRWISE_STATE_MATCH
                else
                    current_gap_left_scores[i] = from_gap
                    gap_left_trace[i, j] = _PAIRWISE_STATE_GAP_LEFT
                end

                from_match = previous_match_scores[i] + gap_open
                from_gap = previous_gap_right_scores[i] + gap_extend
                if from_match >= from_gap
                    current_gap_right_scores[i] = from_match
                    gap_right_trace[i, j] = _PAIRWISE_STATE_MATCH
                else
                    current_gap_right_scores[i] = from_gap
                    gap_right_trace[i, j] = _PAIRWISE_STATE_GAP_RIGHT
                end
            end

            previous_match_scores, current_match_scores = current_match_scores, previous_match_scores
            previous_gap_left_scores, current_gap_left_scores = current_gap_left_scores, previous_gap_left_scores
            previous_gap_right_scores, current_gap_right_scores = current_gap_right_scores, previous_gap_right_scores
        end
    end

    best_score = previous_match_scores[left_length + 1]
    best_state = _PAIRWISE_STATE_MATCH
    if previous_gap_left_scores[left_length + 1] > best_score
        best_score = previous_gap_left_scores[left_length + 1]
        best_state = _PAIRWISE_STATE_GAP_LEFT
    end
    if previous_gap_right_scores[left_length + 1] > best_score
        best_score = previous_gap_right_scores[left_length + 1]
        best_state = _PAIRWISE_STATE_GAP_RIGHT
    end

    return _pairwise_traceback_codon_affine(
        left_tokens,
        right_tokens,
        match_trace,
        gap_left_trace,
        gap_right_trace,
        left_length + 1,
        right_length + 1,
        best_state,
        best_score,
    )
end

"""
    _pairwise_align_codon_affine_local(left_tokens, right_tokens, scoring, gap_open, gap_extend)

Run local codon alignment with an affine gap model.
"""
function _pairwise_align_codon_affine_local(left_tokens::AbstractVector{UInt8}, right_tokens::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap_open::Int, gap_extend::Int)
    left_length = length(left_tokens)
    right_length = length(right_tokens)

    previous_match_scores = zeros(Int, left_length + 1)
    previous_gap_left_scores = zeros(Int, left_length + 1)
    previous_gap_right_scores = zeros(Int, left_length + 1)
    current_match_scores = zeros(Int, left_length + 1)
    current_gap_left_scores = zeros(Int, left_length + 1)
    current_gap_right_scores = zeros(Int, left_length + 1)

    match_trace = zeros(UInt8, left_length + 1, right_length + 1)
    gap_left_trace = zeros(UInt8, left_length + 1, right_length + 1)
    gap_right_trace = zeros(UInt8, left_length + 1, right_length + 1)

    best_score = 0
    best_i = 1
    best_j = 1
    best_state = _PAIRWISE_TRACE_NONE

    @inbounds for j in 2:right_length + 1
        right_token = right_tokens[j - 1]
        current_match_scores[1] = 0
        current_gap_left_scores[1] = 0
        current_gap_right_scores[1] = 0

        for i in 2:left_length + 1
            left_token = left_tokens[i - 1]

            best_prev = 0
            prev_state = _PAIRWISE_TRACE_NONE

            candidate = previous_match_scores[i - 1]
            if candidate > best_prev
                best_prev = candidate
                prev_state = _PAIRWISE_STATE_MATCH
            end
            candidate = previous_gap_left_scores[i - 1]
            if candidate > best_prev
                best_prev = candidate
                prev_state = _PAIRWISE_STATE_GAP_LEFT
            end
            candidate = previous_gap_right_scores[i - 1]
            if candidate > best_prev
                best_prev = candidate
                prev_state = _PAIRWISE_STATE_GAP_RIGHT
            end

            match_score = best_prev + _pairwise_score(scoring, left_token, right_token)
            if match_score > 0
                current_match_scores[i] = match_score
                match_trace[i, j] = prev_state
                if match_score > best_score
                    best_score = match_score
                    best_i = i
                    best_j = j
                    best_state = _PAIRWISE_STATE_MATCH
                end
            else
                current_match_scores[i] = 0
            end

            from_match = current_match_scores[i - 1] + gap_open
            from_gap = current_gap_left_scores[i - 1] + gap_extend
            gap_score = 0
            gap_state = _PAIRWISE_TRACE_NONE
            if from_match > gap_score
                gap_score = from_match
                gap_state = _PAIRWISE_STATE_MATCH
            end
            if from_gap > gap_score
                gap_score = from_gap
                gap_state = _PAIRWISE_STATE_GAP_LEFT
            end
            if gap_score > 0
                current_gap_left_scores[i] = gap_score
                gap_left_trace[i, j] = gap_state
                if gap_score > best_score
                    best_score = gap_score
                    best_i = i
                    best_j = j
                    best_state = _PAIRWISE_STATE_GAP_LEFT
                end
            else
                current_gap_left_scores[i] = 0
            end

            from_match = previous_match_scores[i] + gap_open
            from_gap = previous_gap_right_scores[i] + gap_extend
            gap_score = 0
            gap_state = _PAIRWISE_TRACE_NONE
            if from_match > gap_score
                gap_score = from_match
                gap_state = _PAIRWISE_STATE_MATCH
            end
            if from_gap > gap_score
                gap_score = from_gap
                gap_state = _PAIRWISE_STATE_GAP_RIGHT
            end
            if gap_score > 0
                current_gap_right_scores[i] = gap_score
                gap_right_trace[i, j] = gap_state
                if gap_score > best_score
                    best_score = gap_score
                    best_i = i
                    best_j = j
                    best_state = _PAIRWISE_STATE_GAP_RIGHT
                end
            else
                current_gap_right_scores[i] = 0
            end
        end

        previous_match_scores, current_match_scores = current_match_scores, previous_match_scores
        previous_gap_left_scores, current_gap_left_scores = current_gap_left_scores, previous_gap_left_scores
        previous_gap_right_scores, current_gap_right_scores = current_gap_right_scores, previous_gap_right_scores
    end

    return _pairwise_traceback_codon_affine(
        left_tokens,
        right_tokens,
        match_trace,
        gap_left_trace,
        gap_right_trace,
        best_i,
        best_j,
        best_state,
        best_score,
    )
end

"""
    pairwise_align_codons(left, right; kwargs...)

Align two nucleotide strings in codon space, preserving triplet boundaries.
"""
function pairwise_align_codons(
    left::AbstractString,
    right::AbstractString;
    match::Int=1,
    mismatch::Int=-1,
    substitution_matrix::Union{Nothing,CodonSubstitutionMatrix}=nothing,
    gap::Int=-1,
    gap_open::Union{Nothing,Int}=nothing,
    gap_extend::Union{Nothing,Int}=nothing,
    is_local::Bool=false,
)
    left_tokens = _encode_codon_sequence(codeunits(left))
    right_tokens = _encode_codon_sequence(codeunits(right))
    scoring = substitution_matrix === nothing ? LinearPairwiseScoring(match, mismatch) : CodonMatrixPairwiseScoring(substitution_matrix)
    if gap_open === nothing && gap_extend === nothing
        return is_local ? _pairwise_align_codon_local(left_tokens, right_tokens, scoring, gap) : _pairwise_align_codon_global(left_tokens, right_tokens, scoring, gap)
    end

    affine_gap_open, affine_gap_extend = _pairwise_gap_parameters(gap, gap_open, gap_extend)
    return is_local ? _pairwise_align_codon_affine_local(left_tokens, right_tokens, scoring, affine_gap_open, affine_gap_extend) : _pairwise_align_codon_affine_global(left_tokens, right_tokens, scoring, affine_gap_open, affine_gap_extend)
end

"""
    needleman_wunsch_codons(left, right; kwargs...)

Wrapper for global codon alignment.
"""
needleman_wunsch_codons(left, right; kwargs...) = pairwise_align_codons(left, right; is_local=false, kwargs...)

"""
    smith_waterman_codons(left, right; kwargs...)

Wrapper for local codon alignment.
"""
smith_waterman_codons(left, right; kwargs...) = pairwise_align_codons(left, right; is_local=true, kwargs...)

"""
    local_align_codons(left, right; kwargs...)

Backward-compatible alias for `smith_waterman_codons`.
"""
local_align_codons(left, right; kwargs...) = smith_waterman_codons(left, right; kwargs...)

"""
    needleman_wunsch(left, right; kwargs...)

Explicit wrapper for the global pairwise alignment path.
"""
needleman_wunsch(left, right; kwargs...) = pairwise_align(left, right; is_local=false, kwargs...)

"""
    smith_waterman(left, right; kwargs...)

Explicit wrapper for the local pairwise alignment path.
"""
smith_waterman(left, right; kwargs...) = pairwise_align(left, right; is_local=true, kwargs...)

"""
    local_align(left, right; kwargs...)

Backward-compatible alias for `smith_waterman`.
"""
local_align(left, right; kwargs...) = smith_waterman(left, right; kwargs...)

"""
    _pairwise_align_global(left_bytes, right_bytes, scoring, gap)

Run global pairwise alignment with a linear gap model.
"""
function _pairwise_align_global(left_bytes::AbstractVector{UInt8}, right_bytes::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap::Int)
    left_length = length(left_bytes)
    right_length = length(right_bytes)

    previous_scores = Vector{Int}(undef, left_length + 1)
    current_scores = similar(previous_scores)
    trace = Matrix{UInt8}(undef, left_length + 1, right_length + 1)

    @inbounds begin
        previous_scores[1] = 0
        trace[1, 1] = 0x00

        for i in 2:left_length + 1
            previous_scores[i] = previous_scores[i - 1] + gap
            trace[i, 1] = 0x02
        end

        for j in 2:right_length + 1
            current_scores[1] = previous_scores[1] + gap
            trace[1, j] = 0x03
            right_byte = right_bytes[j - 1]
            for i in 2:left_length + 1
                left_byte = left_bytes[i - 1]

                diag_score = previous_scores[i - 1] + _pairwise_score(scoring, left_byte, right_byte)
                up_score = previous_scores[i] + gap
                left_score = current_scores[i - 1] + gap

                cell_score = diag_score
                cell_trace = 0x01

                if up_score > cell_score
                    cell_score = up_score
                    cell_trace = 0x02
                end
                if left_score > cell_score
                    cell_score = left_score
                    cell_trace = 0x03
                end

                current_scores[i] = cell_score
                trace[i, j] = cell_trace
            end

            previous_scores, current_scores = current_scores, previous_scores
        end
    end

    return _pairwise_traceback(left_bytes, right_bytes, nothing, trace, left_length + 1, right_length + 1, previous_scores[left_length + 1])
end

"""
    _pairwise_align_local(left_bytes, right_bytes, scoring, gap)

Run local pairwise alignment with a linear gap model.
"""
function _pairwise_align_local(left_bytes::AbstractVector{UInt8}, right_bytes::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap::Int)
    left_length = length(left_bytes)
    right_length = length(right_bytes)

    scores = Matrix{Int}(undef, left_length + 1, right_length + 1)
    trace = Matrix{UInt8}(undef, left_length + 1, right_length + 1)

    best_score = 0
    best_i = 1
    best_j = 1

    @inbounds begin
        scores[1, 1] = 0
        trace[1, 1] = 0x00

        for i in 2:left_length + 1
            scores[i, 1] = 0
            trace[i, 1] = 0x00
        end

        for j in 2:right_length + 1
            scores[1, j] = 0
            trace[1, j] = 0x00
        end

        for j in 2:right_length + 1
            right_byte = right_bytes[j - 1]
            for i in 2:left_length + 1
                left_byte = left_bytes[i - 1]

                diag_score = scores[i - 1, j - 1] + _pairwise_score(scoring, left_byte, right_byte)
                up_score = scores[i - 1, j] + gap
                left_score = scores[i, j - 1] + gap

                cell_score = 0
                cell_trace = 0x00

                if diag_score > cell_score
                    cell_score = diag_score
                    cell_trace = 0x01
                end
                if up_score > cell_score
                    cell_score = up_score
                    cell_trace = 0x02
                end
                if left_score > cell_score
                    cell_score = left_score
                    cell_trace = 0x03
                end

                scores[i, j] = cell_score
                trace[i, j] = cell_trace

                if cell_score > best_score
                    best_score = cell_score
                    best_i = i
                    best_j = j
                end
            end
        end
    end

    return _pairwise_traceback(left_bytes, right_bytes, scores, trace, best_i, best_j, best_score)
end

"""
    _pairwise_align_affine_global(left_bytes, right_bytes, scoring, gap_open, gap_extend)

Run global pairwise alignment with an affine gap model.
"""
function _pairwise_align_affine_global(
    left_bytes::AbstractVector{UInt8},
    right_bytes::AbstractVector{UInt8},
    scoring::AbstractPairwiseScoring,
    gap_open::Int,
    gap_extend::Int,
)
    left_length = length(left_bytes)
    right_length = length(right_bytes)
    negative_infinity = typemin(Int) ÷ 4

    previous_match_scores = fill(negative_infinity, left_length + 1)
    previous_gap_left_scores = fill(negative_infinity, left_length + 1)
    previous_gap_right_scores = fill(negative_infinity, left_length + 1)
    current_match_scores = similar(previous_match_scores)
    current_gap_left_scores = similar(previous_gap_left_scores)
    current_gap_right_scores = similar(previous_gap_right_scores)

    match_trace = zeros(UInt8, left_length + 1, right_length + 1)
    gap_left_trace = zeros(UInt8, left_length + 1, right_length + 1)
    gap_right_trace = zeros(UInt8, left_length + 1, right_length + 1)

    previous_match_scores[1] = 0
    previous_gap_left_scores[1] = negative_infinity
    previous_gap_right_scores[1] = negative_infinity

    @inbounds begin
        for i in 2:left_length + 1
            from_match = previous_match_scores[i - 1] + gap_open
            from_gap = previous_gap_left_scores[i - 1] + gap_extend
            if from_match >= from_gap
                previous_gap_left_scores[i] = from_match
                gap_left_trace[i, 1] = _PAIRWISE_STATE_MATCH
            else
                previous_gap_left_scores[i] = from_gap
                gap_left_trace[i, 1] = _PAIRWISE_STATE_GAP_LEFT
            end
        end

        for j in 2:right_length + 1
            current_match_scores[1] = negative_infinity
            current_gap_left_scores[1] = negative_infinity

            from_match = previous_match_scores[1] + gap_open
            from_gap = previous_gap_right_scores[1] + gap_extend
            if from_match >= from_gap
                current_gap_right_scores[1] = from_match
                gap_right_trace[1, j] = _PAIRWISE_STATE_MATCH
            else
                current_gap_right_scores[1] = from_gap
                gap_right_trace[1, j] = _PAIRWISE_STATE_GAP_RIGHT
            end

            right_byte = right_bytes[j - 1]
            for i in 2:left_length + 1
                left_byte = left_bytes[i - 1]

                best_prev = previous_match_scores[i - 1]
                best_state = _PAIRWISE_STATE_MATCH
                if previous_gap_left_scores[i - 1] > best_prev
                    best_prev = previous_gap_left_scores[i - 1]
                    best_state = _PAIRWISE_STATE_GAP_LEFT
                end
                if previous_gap_right_scores[i - 1] > best_prev
                    best_prev = previous_gap_right_scores[i - 1]
                    best_state = _PAIRWISE_STATE_GAP_RIGHT
                end
                current_match_scores[i] = best_prev + _pairwise_score(scoring, left_byte, right_byte)
                match_trace[i, j] = best_state

                from_match = current_match_scores[i - 1] + gap_open
                from_gap = current_gap_left_scores[i - 1] + gap_extend
                if from_match >= from_gap
                    current_gap_left_scores[i] = from_match
                    gap_left_trace[i, j] = _PAIRWISE_STATE_MATCH
                else
                    current_gap_left_scores[i] = from_gap
                    gap_left_trace[i, j] = _PAIRWISE_STATE_GAP_LEFT
                end

                from_match = previous_match_scores[i] + gap_open
                from_gap = previous_gap_right_scores[i] + gap_extend
                if from_match >= from_gap
                    current_gap_right_scores[i] = from_match
                    gap_right_trace[i, j] = _PAIRWISE_STATE_MATCH
                else
                    current_gap_right_scores[i] = from_gap
                    gap_right_trace[i, j] = _PAIRWISE_STATE_GAP_RIGHT
                end
            end

            previous_match_scores, current_match_scores = current_match_scores, previous_match_scores
            previous_gap_left_scores, current_gap_left_scores = current_gap_left_scores, previous_gap_left_scores
            previous_gap_right_scores, current_gap_right_scores = current_gap_right_scores, previous_gap_right_scores
        end
    end

    best_score = previous_match_scores[left_length + 1]
    best_state = _PAIRWISE_STATE_MATCH
    if previous_gap_left_scores[left_length + 1] > best_score
        best_score = previous_gap_left_scores[left_length + 1]
        best_state = _PAIRWISE_STATE_GAP_LEFT
    end
    if previous_gap_right_scores[left_length + 1] > best_score
        best_score = previous_gap_right_scores[left_length + 1]
        best_state = _PAIRWISE_STATE_GAP_RIGHT
    end

    return _pairwise_traceback_affine(
        left_bytes,
        right_bytes,
        match_trace,
        gap_left_trace,
        gap_right_trace,
        left_length + 1,
        right_length + 1,
        best_state,
        best_score,
    )
end

"""
    _pairwise_align_affine_local(left_bytes, right_bytes, scoring, gap_open, gap_extend)

Run local pairwise alignment with an affine gap model.
"""
function _pairwise_align_affine_local(
    left_bytes::AbstractVector{UInt8},
    right_bytes::AbstractVector{UInt8},
    scoring::AbstractPairwiseScoring,
    gap_open::Int,
    gap_extend::Int,
)
    left_length = length(left_bytes)
    right_length = length(right_bytes)

    # O(N) memory buffers instead of O(N * M)
    previous_match_scores = zeros(Int, left_length + 1)
    previous_gap_left_scores = zeros(Int, left_length + 1)
    previous_gap_right_scores = zeros(Int, left_length + 1)
    
    current_match_scores = zeros(Int, left_length + 1)
    current_gap_left_scores = zeros(Int, left_length + 1)
    current_gap_right_scores = zeros(Int, left_length + 1)

    # Trace matrices still require O(N * M) but use only 1 byte per cell
    match_trace = zeros(UInt8, left_length + 1, right_length + 1)
    gap_left_trace = zeros(UInt8, left_length + 1, right_length + 1)
    gap_right_trace = zeros(UInt8, left_length + 1, right_length + 1)

    best_score = 0
    best_i = 1
    best_j = 1
    best_state = _PAIRWISE_TRACE_NONE

    @inbounds for j in 2:right_length + 1
        right_byte = right_bytes[j - 1]
        
        # Reset the first cell of the current column to 0
        current_match_scores[1] = 0
        current_gap_left_scores[1] = 0
        current_gap_right_scores[1] = 0

        for i in 2:left_length + 1
            left_byte = left_bytes[i - 1]

            best_prev = 0
            prev_state = _PAIRWISE_TRACE_NONE
            
            # Diagonal: previous_scores[i - 1] means scores[i-1, j-1]
            candidate = previous_match_scores[i - 1]
            if candidate > best_prev
                best_prev = candidate
                prev_state = _PAIRWISE_STATE_MATCH
            end
            candidate = previous_gap_left_scores[i - 1]
            if candidate > best_prev
                best_prev = candidate
                prev_state = _PAIRWISE_STATE_GAP_LEFT
            end
            candidate = previous_gap_right_scores[i - 1]
            if candidate > best_prev
                best_prev = candidate
                prev_state = _PAIRWISE_STATE_GAP_RIGHT
            end

            match_score = best_prev + _pairwise_score(scoring, left_byte, right_byte)
            if match_score > 0
                current_match_scores[i] = match_score
                match_trace[i, j] = prev_state
                if match_score > best_score
                    best_score = match_score
                    best_i = i
                    best_j = j
                    best_state = _PAIRWISE_STATE_MATCH
                end
            else
                current_match_scores[i] = 0
            end

            # Up: current_scores[i-1] means scores[i-1, j]
            from_match = current_match_scores[i - 1] + gap_open
            from_gap = current_gap_left_scores[i - 1] + gap_extend
            gap_score = 0
            gap_state = _PAIRWISE_TRACE_NONE
            if from_match > gap_score
                gap_score = from_match
                gap_state = _PAIRWISE_STATE_MATCH
            end
            if from_gap > gap_score
                gap_score = from_gap
                gap_state = _PAIRWISE_STATE_GAP_LEFT
            end
            if gap_score > 0
                current_gap_left_scores[i] = gap_score
                gap_left_trace[i, j] = gap_state
                if gap_score > best_score
                    best_score = gap_score
                    best_i = i
                    best_j = j
                    best_state = _PAIRWISE_STATE_GAP_LEFT
                end
            else
                current_gap_left_scores[i] = 0
            end

            # Left: previous_scores[i] means scores[i, j-1]
            from_match = previous_match_scores[i] + gap_open
            from_gap = previous_gap_right_scores[i] + gap_extend
            gap_score = 0
            gap_state = _PAIRWISE_TRACE_NONE
            if from_match > gap_score
                gap_score = from_match
                gap_state = _PAIRWISE_STATE_MATCH
            end
            if from_gap > gap_score
                gap_score = from_gap
                gap_state = _PAIRWISE_STATE_GAP_RIGHT
            end
            if gap_score > 0
                current_gap_right_scores[i] = gap_score
                gap_right_trace[i, j] = gap_state
                if gap_score > best_score
                    best_score = gap_score
                    best_i = i
                    best_j = j
                    best_state = _PAIRWISE_STATE_GAP_RIGHT
                end
            else
                current_gap_right_scores[i] = 0
            end
        end
        
        # Swap current memory to previous
        previous_match_scores, current_match_scores = current_match_scores, previous_match_scores
        previous_gap_left_scores, current_gap_left_scores = current_gap_left_scores, previous_gap_left_scores
        previous_gap_right_scores, current_gap_right_scores = current_gap_right_scores, previous_gap_right_scores
    end

    return _pairwise_traceback_affine(
        left_bytes,
        right_bytes,
        match_trace,
        gap_left_trace,
        gap_right_trace,
        best_i,
        best_j,
        best_state,
        best_score,
    )
end

"""
    _pairwise_traceback(...)

Reconstruct a nucleotide alignment from standard traceback data.
"""
function _pairwise_traceback(
    left_bytes::AbstractVector{UInt8},
    right_bytes::AbstractVector{UInt8},
    scores::Union{Nothing,Matrix{Int}},
    trace::Matrix{UInt8},
    i::Int,
    j::Int,
    best_score::Int,
)
    aligned_left = Vector{UInt8}(undef, length(left_bytes) + length(right_bytes))
    aligned_right = Vector{UInt8}(undef, length(left_bytes) + length(right_bytes))
    write_index = length(aligned_left)
    aligned_length = 0
    matches = 0

    current_i = i
    current_j = j

    while current_i > 1 || current_j > 1
        if scores !== nothing && scores[current_i, current_j] == 0 && trace[current_i, current_j] == 0x00
            break
        end

        direction = trace[current_i, current_j]
        direction == 0x00 && break

        if direction == 0x01
            left_byte = left_bytes[current_i - 1]
            right_byte = right_bytes[current_j - 1]
            aligned_left[write_index] = left_byte
            aligned_right[write_index] = right_byte
            matches += left_byte == right_byte ? 1 : 0
            current_i -= 1
            current_j -= 1
        elseif direction == 0x02
            aligned_left[write_index] = left_bytes[current_i - 1]
            aligned_right[write_index] = UInt8('-')
            current_i -= 1
        else
            aligned_left[write_index] = UInt8('-')
            aligned_right[write_index] = right_bytes[current_j - 1]
            current_j -= 1
        end

        aligned_length += 1
        write_index -= 1
    end

    start_index = write_index + 1
    left_string = String(aligned_left[start_index:end])
    right_string = String(aligned_right[start_index:end])
    identity = aligned_length == 0 ? 0.0 : matches / aligned_length

    return PairwiseAlignmentResult(left_string, right_string, best_score, matches, identity)
end

"""
    _pairwise_traceback_affine(...)

Reconstruct a nucleotide alignment from affine-gap traceback data.
"""
function _pairwise_traceback_affine(
    left_bytes::AbstractVector{UInt8},
    right_bytes::AbstractVector{UInt8},
    match_trace::Matrix{UInt8},
    gap_left_trace::Matrix{UInt8},
    gap_right_trace::Matrix{UInt8},
    i::Int,
    j::Int,
    start_state::UInt8,
    best_score::Int,
)
    aligned_left = Vector{UInt8}(undef, length(left_bytes) + length(right_bytes))
    aligned_right = Vector{UInt8}(undef, length(left_bytes) + length(right_bytes))
    write_index = length(aligned_left)
    aligned_length = 0
    matches = 0

    current_i = i
    current_j = j
    current_state = start_state

    while current_state != _PAIRWISE_TRACE_NONE && (current_i > 1 || current_j > 1)
        if current_state == _PAIRWISE_STATE_MATCH
            left_byte = left_bytes[current_i - 1]
            right_byte = right_bytes[current_j - 1]
            aligned_left[write_index] = left_byte
            aligned_right[write_index] = right_byte
            matches += left_byte == right_byte ? 1 : 0
            current_state = match_trace[current_i, current_j]
            current_i -= 1
            current_j -= 1
        elseif current_state == _PAIRWISE_STATE_GAP_LEFT
            aligned_left[write_index] = left_bytes[current_i - 1]
            aligned_right[write_index] = UInt8('-')
            current_state = gap_left_trace[current_i, current_j]
            current_i -= 1
        else
            aligned_left[write_index] = UInt8('-')
            aligned_right[write_index] = right_bytes[current_j - 1]
            current_state = gap_right_trace[current_i, current_j]
            current_j -= 1
        end

        aligned_length += 1
        write_index -= 1
    end

    start_index = write_index + 1
    left_string = String(aligned_left[start_index:end])
    right_string = String(aligned_right[start_index:end])
    identity = aligned_length == 0 ? 0.0 : matches / aligned_length

    return PairwiseAlignmentResult(left_string, right_string, best_score, matches, identity)
end
