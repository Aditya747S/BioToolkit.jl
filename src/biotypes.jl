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

# ---- Sequence type -----------------------------------------------------------

"""
    BioSequence{A <: BioAlphabet}

A typed biological sequence over alphabet `A`, stored as a contiguous byte
vector (one byte per symbol, uppercase-normalized).

Parametric on the alphabet so that functions can dispatch on the sequence
type without inspecting the data:

```julia
gc_content(seq::BioSequence{DNAAlphabet}) = ...   # only for DNA
translate(seq::BioSequence{DNAAlphabet})  = ...   # only for DNA
```

# Construction

```julia
DNASeq("ACGT")               # from a string
RNASeq("ACGU")
AASeq("MVLSPADKTNVK")
BioSequence{DNAAlphabet}(b"ACGT")  # from raw bytes
```
"""
struct BioSequence{A<:BioAlphabet}
    data::Vector{UInt8}

    function BioSequence{A}(data::Vector{UInt8}; validate::Bool=true) where {A<:BioAlphabet}
        if validate
            @inbounds for i in eachindex(data)
                isvalid_symbol(A, data[i]) || throw(ArgumentError(
                    "invalid symbol '$(Char(data[i]))' (byte $(data[i])) for $(A); " *
                    "expected one of: $(join(sort!([Char(b) for b in symbols(A)]), ", "))"))
            end
        end
        return new{A}(data)
    end
end

# Convenience aliases
const DNASeq = BioSequence{DNAAlphabet}
const RNASeq = BioSequence{RNAAlphabet}
const AASeq = BioSequence{AminoAcidAlphabet}

# ---- Constructors ------------------------------------------------------------

"""
    BioSequence{A}(s::String; validate=true)

Construct a `BioSequence{A}` from a string, normalising to uppercase.
"""
function BioSequence{A}(s::String; validate::Bool=true) where {A<:BioAlphabet}
    bytes = Vector{UInt8}(undef, ncodeunits(s))
    @inbounds for (i, byte) in enumerate(codeunits(s))
        bytes[i] = byte < 0x61 ? byte : (byte <= 0x7a ? byte - 0x20 : byte)
    end
    return BioSequence{A}(bytes; validate=validate)
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

function Base.show(io::IO, seq::BioSequence{A}) where {A}
    name = A === DNAAlphabet ? "DNASeq" :
           A === RNAAlphabet ? "RNASeq" :
           A === AminoAcidAlphabet ? "AASeq" : "BioSequence{$(A)}"
    n = length(seq.data)
    if n <= 60
        print(io, name, "(\"", String(copy(seq.data)), "\")")
    else
        print(io, name, "(\"", String(seq.data[1:30]), "…", String(seq.data[end-29:end]), "\") [", n, " nt]")
    end
end

"""
    String(seq::BioSequence) -> String

Convert a typed biological sequence back to a plain string.
"""
Base.String(seq::BioSequence) = String(copy(seq.data))

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
struct SummarizedExperiment
    assays::Dict{String,Matrix{Float64}}
    rowData::Dict{Symbol,Vector}
    colData::Dict{Symbol,Vector}
    metadata::Dict{Symbol,Any}

    function SummarizedExperiment(
        assays::Dict{String,Matrix{Float64}},
        rowData::Dict{Symbol,Vector},
        colData::Dict{Symbol,Vector},
        metadata::Dict{Symbol,Any}=Dict{Symbol,Any}()
    )
        isempty(assays) && throw(ArgumentError("at least one assay must be provided"))
        # Validate dimensions
        n_features, n_samples = size(first(values(assays)))
        for (name, matrix) in assays
            size(matrix) == (n_features, n_samples) || throw(DimensionMismatch("assay '$name' has dimensions $(size(matrix)), expected ($n_features, $n_samples)"))
        end
        for (key, vec) in rowData
            length(vec) == n_features || throw(DimensionMismatch("rowData column '$key' has length $(length(vec)), expected $n_features"))
        end
        for (key, vec) in colData
            length(vec) == n_samples || throw(DimensionMismatch("colData column '$key' has length $(length(vec)), expected $n_samples"))
        end
        return new(assays, rowData, colData, metadata)
    end
end

"""
    SummarizedExperiment(counts::Matrix; rowData=..., colData=..., metadata=...)

Main constructor for `SummarizedExperiment`.
"""
function SummarizedExperiment(
    counts::Matrix{<:Real};
    assay_name::String="counts",
    rowData::Dict{Symbol,Vector}=Dict{Symbol,Vector}(),
    colData::Dict{Symbol,Vector}=Dict{Symbol,Vector}(),
    metadata::Dict{Symbol,Any}=Dict{Symbol,Any}()
)
    return SummarizedExperiment(
        Dict(assay_name => Matrix{Float64}(counts)),
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
        if node.left !== nothing && left < node.left.left_endpoint
            return _itn_rotate_right(node)
        elseif node.left !== nothing
            node.left = _itn_rotate_left(node.left)
            return _itn_rotate_right(node)
        end
    end

    # Right-heavy
    if balance < -1
        if node.right !== nothing && (left > node.right.left_endpoint || (left == node.right.left_endpoint && right >= node.right.right_endpoint))
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
