# `biotypes.jl` - Native Biological Type System

## Overview

`biotypes.jl` is BioToolkit's zero-dependency biological type foundation. It provides **parametric type safety** for biological sequences and a small set of core genomics data structures implemented in pure Julia: typed biological sequences, run-length encoding, a `SummarizedExperiment` assay container, and an augmented interval tree.

### Purpose

When working with biological data, it is easy to accidentally apply a DNA-only operation such as translation, reverse-complement, or GC-content calculation to RNA or protein input. `biotypes.jl` prevents many of these mistakes through Julia's type system: every `BioSequence` carries its alphabet as a type parameter, and downstream functions dispatch on that parameter.

For example, a method declared for `BioSequence{DNAAlphabet}` will not accept `BioSequence{AminoAcidAlphabet}`. This gives BioToolkit a type-safe internal representation while preserving the ergonomic feel of string-like sequence manipulation.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Alphabets are singleton types** | `DNAAlphabet`, `RNAAlphabet`, and `AminoAcidAlphabet` carry no fields. They exist purely for dispatch, so the compiler can specialize methods without storing alphabet data in every sequence. |
| **One byte per symbol** | Sequences store raw `UInt8` values. A 2-bit DNA encoding would reduce memory, but it would complicate slicing, printing, validation, and interoperability with parsers and `String`-based code. BioToolkit targets interactive analysis more than whole-genome assembly storage. |
| **Uppercase normalization** | String constructors normalize ASCII lowercase input to uppercase. This removes case-sensitivity as a common source of bugs. |
| **Validation on construction** | Constructors validate by default and throw `ArgumentError` with the invalid symbol and expected alphabet. Trusted internal code can skip validation with `validate=false`. |
| **Gap is a valid symbol** | `'-'` is accepted in DNA, RNA, and amino-acid alphabets so aligned sequences can remain typed. |
| **String-like interface** | `BioSequence` supports length, indexing, iteration, conversion to `String`, equality with strings, and prefix checks so existing code can migrate gradually. |
| **Metadata containers are dependency-light** | `Rle`, `SummarizedExperiment`, and `IntervalTree` are implemented without bringing in heavy external packages. |

---

## Table of Contents

1. [Alphabet Types](#1-alphabet-types)
2. [Symbol Validation Functions](#2-symbol-validation-functions)
3. [BioSequence Type](#3-biosequence-type)
4. [BioSequence Constructors](#4-biosequence-constructors)
5. [String-Like Interface](#5-string-like-interface)
6. [Alphabet Query Functions](#6-alphabet-query-functions)
7. [Validation Helpers](#7-validation-helpers)
8. [Run-Length Encoding (`Rle`)](#8-run-length-encoding-rle)
9. [SummarizedExperiment](#9-summarizedexperiment)
10. [Interval Tree](#10-interval-tree)
11. [Implementation Notes and Caveats](#11-implementation-notes-and-caveats)
12. [Quick Reference](#12-quick-reference)
13. [Complete Usage Example](#13-complete-usage-example)

---

## 1. Alphabet Types

### `BioAlphabet`

```julia
abstract type BioAlphabet end
```

**Kind:** Abstract supertype

**Description:** Root of BioToolkit's biological alphabet hierarchy. Concrete alphabet types subtype `BioAlphabet` and are used as type parameters for `BioSequence{A}`.

`BioAlphabet` has no fields and is never instantiated directly. It exists to constrain type parameters and make dispatch signatures readable.

**Usage:**

```julia
BioSequence{DNAAlphabet}
BioSequence{RNAAlphabet}
BioSequence{AminoAcidAlphabet}
```

---

### `DNAAlphabet`

```julia
struct DNAAlphabet <: BioAlphabet end
```

**Kind:** Concrete singleton type

**Description:** Represents the standard IUPAC DNA alphabet, including ambiguity symbols and gaps. Lowercase input is accepted by the string constructor and stored as uppercase.

**Valid symbols:**

| Symbol | Meaning |
|---|---|
| `A` | Adenine |
| `C` | Cytosine |
| `G` | Guanine |
| `T` | Thymine |
| `N` | Any base |
| `R` | A or G, purine |
| `Y` | C or T, pyrimidine |
| `S` | G or C |
| `W` | A or T |
| `K` | G or T |
| `M` | A or C |
| `B` | C, G, or T; not A |
| `D` | A, G, or T; not C |
| `H` | A, C, or T; not G |
| `V` | A, C, or G; not T |
| `-` | Gap |

**Internal constant:**

```julia
const _DNA_VALID_BYTES = Set{UInt8}(UInt8.(collect("ACGTNRYSWKMBDHV-")))
```

---

### `RNAAlphabet`

```julia
struct RNAAlphabet <: BioAlphabet end
```

**Kind:** Concrete singleton type

**Description:** Represents the standard IUPAC RNA alphabet. It is identical to the DNA alphabet except that `U` replaces `T`.

**Valid symbols:**

| Symbol | Meaning |
|---|---|
| `A` | Adenine |
| `C` | Cytosine |
| `G` | Guanine |
| `U` | Uracil |
| `N` | Any base |
| `R` | A or G |
| `Y` | C or U |
| `S` | G or C |
| `W` | A or U |
| `K` | G or U |
| `M` | A or C |
| `B` | C, G, or U; not A |
| `D` | A, G, or U; not C |
| `H` | A, C, or U; not G |
| `V` | A, C, or G; not U |
| `-` | Gap |

**Internal constant:**

```julia
const _RNA_VALID_BYTES = Set{UInt8}(UInt8.(collect("ACGUNRYSWKMBDHV-")))
```

---

### `AminoAcidAlphabet`

```julia
struct AminoAcidAlphabet <: BioAlphabet end
```

**Kind:** Concrete singleton type

**Description:** Represents the IUPAC amino-acid alphabet plus common ambiguity, stop, rare amino-acid, unknown, and gap symbols.

**Valid symbols:**

| Symbol | Meaning |
|---|---|
| `A` | Alanine |
| `C` | Cysteine |
| `D` | Aspartic acid |
| `E` | Glutamic acid |
| `F` | Phenylalanine |
| `G` | Glycine |
| `H` | Histidine |
| `I` | Isoleucine |
| `K` | Lysine |
| `L` | Leucine |
| `M` | Methionine |
| `N` | Asparagine |
| `P` | Proline |
| `Q` | Glutamine |
| `R` | Arginine |
| `S` | Serine |
| `T` | Threonine |
| `V` | Valine |
| `W` | Tryptophan |
| `Y` | Tyrosine |
| `B` | Asparagine or aspartic acid |
| `Z` | Glutamine or glutamic acid |
| `X` | Unknown amino acid |
| `U` | Selenocysteine |
| `O` | Pyrrolysine |
| `*` | Translation stop |
| `-` | Gap |

**Internal constant:**

```julia
const _AA_VALID_BYTES = Set{UInt8}(UInt8.(collect("ACDEFGHIKLMNPQRSTVWYXBZUO*-")))
```

---

## 2. Symbol Validation Functions

### `symbols(::Type{A}) -> Set{UInt8}`

```julia
symbols(::Type{DNAAlphabet}) = _DNA_VALID_BYTES
symbols(::Type{RNAAlphabet}) = _RNA_VALID_BYTES
symbols(::Type{AminoAcidAlphabet}) = _AA_VALID_BYTES
```

**Description:** Returns the valid uppercase byte set for an alphabet type. This is the canonical validation source used by `isvalid_symbol`, `validate_sequence`, and the `BioSequence` constructor.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `A` | `Type{<:BioAlphabet}` | Alphabet type. Pass `DNAAlphabet`, not `DNAAlphabet()`. |

**Returns:** `Set{UInt8}`

**Example:**

```julia
julia> Char.(sort!(collect(symbols(RNAAlphabet))))
16-element Vector{Char}:
 '-'
 'A'
 'B'
 'C'
 'D'
 'G'
 'H'
 'K'
 'M'
 'N'
 'R'
 'S'
 'U'
 'V'
 'W'
 'Y'
```

---

### `gap(::Type{A}) -> UInt8`

```julia
gap(::Type{<:BioAlphabet}) = UInt8('-')
```

**Description:** Returns the gap character byte for any supported alphabet. All BioToolkit alphabets use `'-'`.

**Returns:** `UInt8`

**Example:**

```julia
julia> Char(gap(DNAAlphabet))
'-'
```

---

### `isvalid_symbol(::Type{A}, byte::UInt8) -> Bool`

```julia
@inline function isvalid_symbol(::Type{A}, byte::UInt8) where {A<:BioAlphabet}
    upper = byte < 0x61 ? byte : (byte <= 0x7a ? byte - 0x20 : byte)
    return upper in symbols(A)
end
```

**Description:** Checks whether a raw byte is valid for an alphabet. The byte is first normalized to uppercase if it is an ASCII lowercase letter.

**Uppercase logic:**

| Condition | Action |
|---|---|
| `byte < 0x61` | Use byte unchanged. |
| `0x61 <= byte <= 0x7a` | Subtract `0x20` to convert lowercase ASCII to uppercase. |
| `byte > 0x7a` | Use byte unchanged; normally fails validation. |

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `A` | `Type{<:BioAlphabet}` | Alphabet type. |
| `byte` | `UInt8` | Raw byte to validate. |

**Returns:** `Bool`

**Example:**

```julia
julia> isvalid_symbol(DNAAlphabet, UInt8('A'))
true

julia> isvalid_symbol(DNAAlphabet, UInt8('a'))
true

julia> isvalid_symbol(DNAAlphabet, UInt8('U'))
false

julia> isvalid_symbol(RNAAlphabet, UInt8('U'))
true
```

---

## 3. BioSequence Type

### `BioSequence{A <: BioAlphabet}`

```julia
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
```

**Kind:** Parametric concrete struct

**Description:** A typed biological sequence over alphabet `A`, stored as a contiguous byte vector. The alphabet is encoded in the type, not in a runtime field.

**Fields:**

| Field | Type | Description |
|---|---|---|
| `data` | `Vector{UInt8}` | Raw byte storage. For string construction, bytes are uppercase-normalized before validation. |

**Type parameter:**

| Parameter | Bound | Meaning |
|---|---|---|
| `A` | `<: BioAlphabet` | Biological alphabet for the sequence. |

**Inner constructor behavior:**

1. If `validate=true`, each byte is checked with `isvalid_symbol(A, data[i])`.
2. Invalid bytes throw `ArgumentError` with the invalid character, byte value, alphabet, and accepted symbols.
3. If `validate=false`, validation is skipped and the byte vector is stored directly.

**Type aliases:**

```julia
const DNASeq = BioSequence{DNAAlphabet}
const RNASeq = BioSequence{RNAAlphabet}
const AASeq  = BioSequence{AminoAcidAlphabet}
```

Because these are `const` aliases, not new wrapper types, constructors such as `DNASeq("ACGT")` resolve directly to `BioSequence{DNAAlphabet}(::String)`.

**Example:**

```julia
julia> DNASeq("ACGTNGATC")
DNASeq("ACGTNGATC")

julia> RNASeq("acgu")
RNASeq("ACGU")

julia> AASeq("MVLSPADKTNVK")
AASeq("MVLSPADKTNVK")

julia> DNASeq("ACGU")
ERROR: ArgumentError: invalid symbol 'U' (byte 85) for DNAAlphabet; expected one of: -, A, B, C, D, G, H, K, M, N, R, S, T, V, W, Y
```

---

## 4. BioSequence Constructors

### From raw bytes

```julia
BioSequence{A}(data::Vector{UInt8}; validate::Bool=true) where {A<:BioAlphabet}
```

**Description:** Constructs a typed sequence from raw byte storage. This is the inner constructor shown above.

**Parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data` | `Vector{UInt8}` | required | Byte vector containing sequence symbols. |
| `validate` | `Bool` | `true` | Whether to validate every byte against alphabet `A`. |

**Returns:** `BioSequence{A}`

**Example:**

```julia
julia> BioSequence{DNAAlphabet}(UInt8.(collect("ACGT")))
DNASeq("ACGT")

julia> DNASeq(UInt8.(collect("ACGT")); validate=false)
DNASeq("ACGT")
```

**Caution:** `validate=false` can create biologically invalid objects. It is intended for trusted internal paths.

---

### From a string

```julia
function BioSequence{A}(s::String; validate::Bool=true) where {A<:BioAlphabet}
    bytes = Vector{UInt8}(undef, ncodeunits(s))
    @inbounds for (i, byte) in enumerate(codeunits(s))
        bytes[i] = byte < 0x61 ? byte : (byte <= 0x7a ? byte - 0x20 : byte)
    end
    return BioSequence{A}(bytes; validate=validate)
end
```

**Description:** Constructs a typed sequence from a `String`. The constructor walks the string's code units, normalizes ASCII lowercase letters to uppercase, then delegates to the byte-vector constructor.

**Parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `s` | `String` | required | Input sequence text. |
| `validate` | `Bool` | `true` | Whether to validate normalized bytes. |

**Returns:** `BioSequence{A}`

**Example:**

```julia
julia> DNASeq("acgt")
DNASeq("ACGT")

julia> RNASeq("acgu")
RNASeq("ACGU")

julia> AASeq("mteyk*")
AASeq("MTEYK*")
```

**Important note:** Do not define methods such as:

```julia
DNASeq(s) = BioSequence{DNAAlphabet}(s)
```

Since `DNASeq` is already an alias for `BioSequence{DNAAlphabet}`, this kind of definition resolves recursively and is not needed.

---

## 5. String-Like Interface

`BioSequence` implements common `Base` methods so that typed biological sequences can be used in loops, comprehensions, indexing operations, comparisons, and display.

### Length and storage

| Method | Signature | Description | Complexity |
|---|---|---|---|
| `length` | `length(seq::BioSequence)` | Number of symbols. | O(1) |
| `sizeof` | `sizeof(seq::BioSequence)` | Number of bytes, equal to sequence length. | O(1) |
| `isempty` | `isempty(seq::BioSequence)` | Whether the sequence has zero symbols. | O(1) |
| `ncodeunits` | `ncodeunits(seq::BioSequence)` | Number of code units, equal to sequence length. | O(1) |
| `codeunits` | `codeunits(seq::BioSequence)` | Returns the underlying byte vector. | O(1) |

**Example:**

```julia
julia> seq = DNASeq("ACGT")
DNASeq("ACGT")

julia> length(seq)
4

julia> sizeof(seq)
4
```

**Caution:** `codeunits(seq)` returns `seq.data` directly. Treat it as read-only unless you are writing trusted internal code.

---

### Indexing

| Method | Signature | Description |
|---|---|---|
| `firstindex` | `firstindex(seq::BioSequence) -> Int` | Always returns `1`. |
| `lastindex` | `lastindex(seq::BioSequence) -> Int` | Returns `length(seq.data)`. |
| `getindex` scalar | `seq[i::Integer] -> Char` | Returns the symbol at one-based position `i`. |
| `getindex` range | `seq[r::UnitRange] -> BioSequence{A}` | Returns a same-alphabet subsequence. |

**Range slicing implementation:**

```julia
function Base.getindex(seq::BioSequence{A}, r::UnitRange{<:Integer}) where {A}
    return BioSequence{A}(seq.data[r]; validate=false)
end
```

Validation is skipped because the slice originates from an already validated sequence.

**Example:**

```julia
julia> seq = DNASeq("ACGTACGT")
DNASeq("ACGTACGT")

julia> seq[1]
'A'

julia> seq[3:6]
DNASeq("GTAC")
```

---

### Iteration

```julia
Base.iterate(seq::BioSequence)
Base.iterate(seq::BioSequence, i::Int)
Base.eltype(::Type{<:BioSequence}) = Char
```

**Description:** Iteration yields `Char` values.

**Example:**

```julia
julia> collect(DNASeq("ACG"))
3-element Vector{Char}:
 'A'
 'C'
 'G'
```

---

### Equality and prefix checks

```julia
Base.:(==)(a::BioSequence{A}, b::BioSequence{A}) where {A}
Base.:(==)(a::BioSequence, b::String)
Base.:(==)(a::String, b::BioSequence)
Base.startswith(sequence::BioSequence, prefix::AbstractString)
Base.startswith(sequence::BioSequence, prefix::BioSequence)
```

**Description:** Same-alphabet sequences compare by byte vector. Comparisons to strings convert the sequence to `String`. Prefix checks also use string conversion.

**Cross-alphabet note:** There is no specialized `==` method for `BioSequence{DNAAlphabet}` versus `BioSequence{RNAAlphabet}`. In ordinary Julia equality, different concrete sequence objects with the same bytes evaluate as not equal.

**Example:**

```julia
julia> DNASeq("ACGT") == DNASeq("ACGT")
true

julia> DNASeq("ACGT") == RNASeq("ACGT")
false

julia> DNASeq("ACGT") == "ACGT"
true

julia> startswith(DNASeq("ACGTACG"), "ACGT")
true
```

---

### Hashing

```julia
function Base.hash(seq::BioSequence, h::UInt)
    return hash(seq.data, hash(:BioSequence, h))
end
```

**Description:** Allows `BioSequence` values to be used in dictionaries and sets. The hash is based on a `:BioSequence` tag and the raw byte vector.

---

### Display

```julia
function Base.show(io::IO, seq::BioSequence{A}) where {A}
    name = A === DNAAlphabet ? "DNASeq" :
           A === RNAAlphabet ? "RNASeq" :
           A === AminoAcidAlphabet ? "AASeq" : "BioSequence{$(A)}"
    n = length(seq.data)
    if n <= 60
        print(io, name, "\"", String(copy(seq.data)), "\"")
    else
        print(io, name, "\"", String(seq.data[1:30]), "...", String(seq.data[end-29:end]), "\") [", n, " nt]")
    end
end
```

**Description:** Short sequences print in constructor-like form. Long sequences print the first 30 and last 30 symbols plus the total length.

---

### String conversion

```julia
Base.String(seq::BioSequence) = String(copy(seq.data))
Base.convert(::Type{String}, seq::BioSequence) = String(seq)
```

**Description:** Converts a typed sequence back to a plain Julia `String`. The byte vector is copied before conversion so the resulting immutable `String` does not alias mutable sequence storage.

---

## 6. Alphabet Query Functions

### `alphabet(seq::BioSequence{A}) -> Type{A}`

```julia
alphabet(::BioSequence{A}) where {A} = A
```

**Description:** Returns the alphabet type of a sequence.

**Returns:** `Type{A}`

**Example:**

```julia
julia> alphabet(DNASeq("ACGT"))
DNAAlphabet
```

---

### `isdna(seq) -> Bool`

```julia
isdna(::BioSequence{DNAAlphabet}) = true
isdna(::BioSequence) = false
isdna(::String) = false
```

**Description:** Returns `true` only for `BioSequence{DNAAlphabet}`.

---

### `isrna(seq) -> Bool`

```julia
isrna(::BioSequence{RNAAlphabet}) = true
isrna(::BioSequence) = false
isrna(::String) = false
```

**Description:** Returns `true` only for `BioSequence{RNAAlphabet}`.

---

### `isaminoacid(seq) -> Bool`

```julia
isaminoacid(::BioSequence{AminoAcidAlphabet}) = true
isaminoacid(::BioSequence) = false
isaminoacid(::String) = false
```

**Description:** Returns `true` only for `BioSequence{AminoAcidAlphabet}`.

---

## 7. Validation Helpers

### `validate_sequence(::Type{A}, s::String) -> Bool`

```julia
function validate_sequence(::Type{A}, s::String) where {A<:BioAlphabet}
    @inbounds for byte in codeunits(s)
        isvalid_symbol(A, byte) || return false
    end
    return true
end
```

**Description:** Pure predicate that checks whether every byte in a string is valid for alphabet `A`. Unlike constructors, this function does not throw for invalid symbols.

**Complexity:** O(n), where n is the number of code units in `s`.

**Example:**

```julia
julia> validate_sequence(DNAAlphabet, "ACGTN")
true

julia> validate_sequence(DNAAlphabet, "ACGU")
false
```

---

### `validate_aa_sequence(s::String) -> Bool`

```julia
validate_aa_sequence(s::String) = validate_sequence(AminoAcidAlphabet, s)
```

**Description:** Convenience wrapper for amino-acid validation.

---

## 8. Run-Length Encoding (`Rle`)

### `Rle{T}`

```julia
struct Rle{T}
    values::Vector{T}
    lengths::Vector{Int}

    function Rle{T}(values::Vector{T}, lengths::Vector{Int}) where {T}
        length(values) == length(lengths) || throw(ArgumentError("values and lengths must have same length"))
        all(lengths .> 0) || throw(ArgumentError("all lengths must be positive"))
        return new{T}(values, lengths)
    end
end
```

**Kind:** Parametric concrete struct

**Description:** Run-length encoded vector storing repeated values as `(value, length)` pairs. It is equivalent in spirit to Bioconductor's `IRanges::Rle`.

**Fields:**

| Field | Type | Description |
|---|---|---|
| `values` | `Vector{T}` | Run values in order. |
| `lengths` | `Vector{Int}` | Length of each run. Must be positive. |

**Validation:**

- `length(values) == length(lengths)`
- `all(lengths .> 0)`

---

### `Rle(data::AbstractVector{T})`

**Description:** Compresses adjacent identical values in one pass.

**Algorithm:**

1. Return an empty `Rle{T}` if input is empty.
2. Track the current value and run count.
3. Push a run whenever the value changes.
4. Push the final run after the loop.

**Complexity:** O(n)

**Example:**

```julia
julia> r = Rle([0, 0, 0, 1, 1, 1, 1, 1, 0, 0])
Rle{Int64}(3 runs, total length 10)

julia> r.values
3-element Vector{Int64}:
 0
 1
 0
```

---

### Rle interface methods

| Method | Description | Complexity |
|---|---|---|
| `length(rle::Rle)` | Returns decompressed length. | O(number of runs) because it sums lengths. |
| `isempty(rle::Rle)` | Returns whether there are no runs. | O(1) |
| `eltype(::Type{Rle{T}})` | Returns `T`. | O(1) |
| `rle[i]` | Returns the decoded value at one-based position `i`. | O(number of runs) |
| `decode(rle)` | Expands to a full vector. | O(decoded length) |
| `show(rle)` | Prints run count and total length. | O(number of runs) via `length(rle)` |

**Example:**

```julia
r = Rle([0,0,0,1,1,1,1,1,0,0])
r[1]       # 0
r[4]       # 1
r[10]      # 0
decode(r)  # original vector
```

---

## 9. SummarizedExperiment

### `SummarizedExperiment`

```julia
struct SummarizedExperiment
    assays::Dict{String,Matrix{Float64}}
    rowData::Dict{Symbol,Vector}
    colData::Dict{Symbol,Vector}
    metadata::Dict{Symbol,Any}
end
```

**Kind:** Concrete struct

**Description:** Matrix-like container that keeps assays aligned with row and column metadata. It is modeled after Bioconductor's `SummarizedExperiment`.

The central invariant is:

- every assay has the same `(n_features, n_samples)` dimensions;
- every `rowData` vector has length `n_features`;
- every `colData` vector has length `n_samples`.

**Fields:**

| Field | Type | Description |
|---|---|---|
| `assays` | `Dict{String,Matrix{Float64}}` | Named assay matrices. Rows are features; columns are samples. |
| `rowData` | `Dict{Symbol,Vector}` | Feature metadata. Each vector must match number of rows. |
| `colData` | `Dict{Symbol,Vector}` | Sample metadata. Each vector must match number of columns. |
| `metadata` | `Dict{Symbol,Any}` | Experiment-level metadata. |

---

### Full constructor

```julia
SummarizedExperiment(assays, rowData, colData, metadata=Dict{Symbol,Any}())
```

**Validation:**

1. At least one assay is required.
2. Every assay must have the same dimensions.
3. Every `rowData` vector must have length `n_features`.
4. Every `colData` vector must have length `n_samples`.

**Errors:**

| Error | Condition |
|---|---|
| `ArgumentError` | `assays` is empty. |
| `DimensionMismatch` | Assay, row metadata, or column metadata dimensions do not match. |

---

### Counts constructor

```julia
SummarizedExperiment(counts::Matrix{<:Real}; assay_name="counts", rowData=..., colData=..., metadata=...)
```

**Description:** Convenience constructor for a single assay. Converts `counts` to `Matrix{Float64}` and stores it under `assay_name`.

**Parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `counts` | `Matrix{<:Real}` | required | Primary assay matrix. |
| `assay_name` | `String` | `"counts"` | Name of the assay. |
| `rowData` | `Dict{Symbol,Vector}` | empty | Feature metadata. |
| `colData` | `Dict{Symbol,Vector}` | empty | Sample metadata. |
| `metadata` | `Dict{Symbol,Any}` | empty | Experiment metadata. |

**Example:**

```julia
counts = [10 20 30; 5 15 25; 8 12 18]

se = SummarizedExperiment(counts;
    rowData=Dict(:gene => ["BRCA1", "TP53", "MYC"]),
    colData=Dict(:condition => ["control", "treatment", "treatment"]))
```

---

### Accessor functions

| Function | Signature | Description |
|---|---|---|
| `assay` | `assay(se::SummarizedExperiment, name::String)` | Returns a named assay matrix. |
| `assay` | `assay(se::SummarizedExperiment)` | Returns the first assay in the dictionary. |
| `rowData` | `rowData(se::SummarizedExperiment)` | Returns feature metadata. |
| `colData` | `colData(se::SummarizedExperiment)` | Returns sample metadata. |
| `metadata` | `metadata(se::SummarizedExperiment)` | Returns experiment-level metadata. |

**Note:** `assay(se)` uses `first(values(se.assays))`. Dictionary order is not a semantic assay selection rule; prefer named access for reproducibility.

---

### `subset_features(se, feature_indices)`

**Description:** Returns a new experiment containing only selected rows/features. Assays and row metadata are subset together; column metadata and experiment metadata are copied.

**Parameters:**

| Parameter | Description |
|---|---|
| `feature_indices` | Any valid Julia row index: range, integer vector, boolean mask, etc. |

**Returns:** `SummarizedExperiment`

---

### `subset_samples(se, sample_indices)`

**Description:** Returns a new experiment containing only selected columns/samples. Assays and column metadata are subset together; row metadata and experiment metadata are copied.

---

### `Base.show(io::IO, se::SummarizedExperiment)`

Displays:

```text
SummarizedExperiment(n_features features, n_samples samples, n assay(s))
```

---

## 10. Interval Tree

The interval tree is a self-balancing AVL tree augmented with `max_end` for efficient overlap queries. It replaces a sorted-array approach that can degrade to O(n) for large genomic interval sets.

### `IntervalTreeNode{T}`

```julia
mutable struct IntervalTreeNode{T}
    left_endpoint::Int
    right_endpoint::Int
    max_end::Int
    payload::T
    left::Union{Nothing,IntervalTreeNode{T}}
    right::Union{Nothing,IntervalTreeNode{T}}
    height::Int
end
```

**Kind:** Mutable parametric struct

**Description:** Internal node storing one closed interval and its payload. `max_end` stores the largest right endpoint in the subtree, allowing overlap queries to skip subtrees that cannot overlap the query.

**Fields:**

| Field | Type | Description |
|---|---|---|
| `left_endpoint` | `Int` | Inclusive interval start. |
| `right_endpoint` | `Int` | Inclusive interval end. |
| `max_end` | `Int` | Maximum endpoint in this subtree. |
| `payload` | `T` | User data associated with the interval. |
| `left` | `Union{Nothing,IntervalTreeNode{T}}` | Left child. |
| `right` | `Union{Nothing,IntervalTreeNode{T}}` | Right child. |
| `height` | `Int` | AVL height. |

---

### `IntervalTree{T}`

```julia
struct IntervalTree{T}
    root::Base.RefValue{Union{Nothing,IntervalTreeNode{T}}}
end

IntervalTree{T}() where {T}
```

**Kind:** Parametric concrete struct

**Description:** Self-balancing augmented interval tree. The payload type `T` determines the type of values returned by `query_overlaps`.

**Fields:**

| Field | Type | Description |
|---|---|---|
| `root` | `Base.RefValue{Union{Nothing,IntervalTreeNode{T}}}` | Mutable root reference. Rotations can replace the root. |

---

### Internal AVL helper functions

| Function | Purpose |
|---|---|
| `_itn_height(node)` | Returns `node.height`, or `0` for `nothing`. |
| `_itn_max_end(node)` | Returns `node.max_end`, or `typemin(Int)` for `nothing`. |
| `_itn_update!(node)` | Recomputes `height` and `max_end` after insertion or rotation. |
| `_itn_balance(node)` | Returns AVL balance factor: left height minus right height. |
| `_itn_rotate_right(y)` | Performs a right AVL rotation and updates cached fields. |
| `_itn_rotate_left(x)` | Performs a left AVL rotation and updates cached fields. |
| `_itn_insert(node, left, right, payload)` | Recursive insert with AVL rebalancing. |
| `_itn_query(node, ql, qr, results)` | Recursive overlap query with pruning. |
| `_itn_count(node)` | Counts nodes in a subtree. |

---

### `Base.insert!(tree::IntervalTree{T}, left::Int, right::Int, payload::T)`

**Description:** Inserts a closed interval `[left, right]` and associated payload. The tree remains AVL-balanced and every affected node has its `max_end` recomputed.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `tree` | `IntervalTree{T}` | Tree to modify. |
| `left` | `Int` | Inclusive interval start. |
| `right` | `Int` | Inclusive interval end. |
| `payload` | `T` | User value associated with the interval. |

**Returns:** The same `tree`, allowing chaining.

**Ordering rule:** Intervals are ordered by `left_endpoint`, with ties broken by `right_endpoint`.

**Complexity:** O(log n)

---

### `query_overlaps(tree::IntervalTree{T}, query_left::Int, query_right::Int)`

**Description:** Returns all payloads whose intervals overlap the query interval `[query_left, query_right]`.

**Overlap rule:** Two intervals `[a, b]` and `[c, d]` overlap when `a <= d && c <= b`.

**Query algorithm:**

1. If the node is `nothing`, return.
2. If `node.max_end < query_left`, no interval in the subtree can overlap; prune.
3. Search the left subtree.
4. Test the current node.
5. If `node.left_endpoint > query_right`, all right-side intervals start too late; prune.
6. Search the right subtree.

**Complexity:** O(log n + k) for balanced/pruned queries, where k is the number of returned payloads.

**Ordering:** Results are returned in tree traversal order, not guaranteed genomic sort order.

---

### `Base.length(tree::IntervalTree) -> Int`

**Description:** Counts intervals by traversing the tree.

**Complexity:** O(n)

---

### `Base.isempty(tree::IntervalTree) -> Bool`

**Description:** Returns `true` if the tree contains no intervals.

**Complexity:** O(1)

---

## 11. Implementation Notes and Caveats

### Export behavior

`biotypes.jl` defines foundational types and functions used throughout BioToolkit. In the current source file, these definitions are available within the `BioToolkit` module because `biotypes.jl` is included by `BioToolkitBody.jl`. The file itself does not contain explicit `export` statements.

### `validate=false`

Use `validate=false` only when data are known to be valid. Examples include slicing an already validated `BioSequence`, reconstructing aligned sequences, parser internals that already performed validation, and performance-sensitive internal conversions.

### Mutable byte storage

`BioSequence` is an immutable struct, but its `data::Vector{UInt8}` field is mutable. Mutating `seq.data` or the object returned by `codeunits(seq)` can create invalid sequence objects. Treat sequence data as read-only unless writing trusted internal code.

### String conversion copies

`String(seq)` copies `seq.data` before constructing a `String`. This avoids ownership and aliasing problems with Julia's immutable string representation.

### `SummarizedExperiment` assay type

All assays are stored as `Matrix{Float64}`. The convenience constructor converts real-valued matrices to `Float64`, so integer count matrices do not preserve their integer element type inside the container.

### `Rle` random access

`rle[i]` scans the run list each time. This is efficient for few-run data but not ideal for many random queries over many runs. Decode once if repeated random access dominates.

### Interval coordinates

The interval tree stores plain integer intervals and does not enforce chromosome names, strand, or one-based indexing. BioToolkit convention is closed integer intervals; callers are responsible for using consistent coordinate systems.

---

## 12. Quick Reference

### Type hierarchy

```text
BioAlphabet (abstract)
├── DNAAlphabet
├── RNAAlphabet
└── AminoAcidAlphabet

BioSequence{A<:BioAlphabet}
├── DNASeq = BioSequence{DNAAlphabet}
├── RNASeq = BioSequence{RNAAlphabet}
└── AASeq  = BioSequence{AminoAcidAlphabet}

Rle{T}

SummarizedExperiment

IntervalTreeNode{T}
└── IntervalTree{T}
```

### Complete function reference

| Function | Signature | Returns | Complexity |
|---|---|---|---|
| `symbols` | `(::Type{A})` | `Set{UInt8}` | O(1) |
| `gap` | `(::Type{<:BioAlphabet})` | `UInt8` | O(1) |
| `isvalid_symbol` | `(::Type{A}, byte::UInt8)` | `Bool` | O(1) |
| `BioSequence{A}` | `(data::Vector{UInt8}; validate=true)` | `BioSequence{A}` | O(n) with validation |
| `BioSequence{A}` | `(s::String; validate=true)` | `BioSequence{A}` | O(n) |
| `length` | `(seq::BioSequence)` | `Int` | O(1) |
| `sizeof` | `(seq::BioSequence)` | `Int` | O(1) |
| `isempty` | `(seq::BioSequence)` | `Bool` | O(1) |
| `codeunits` | `(seq::BioSequence)` | `Vector{UInt8}` | O(1) |
| `getindex` | `(seq::BioSequence, i::Integer)` | `Char` | O(1) |
| `getindex` | `(seq::BioSequence{A}, r::UnitRange)` | `BioSequence{A}` | O(k) |
| `String` | `(seq::BioSequence)` | `String` | O(n) |
| `alphabet` | `(seq::BioSequence{A})` | `Type{A}` | O(1) |
| `isdna` | `(seq)` | `Bool` | O(1) |
| `isrna` | `(seq)` | `Bool` | O(1) |
| `isaminoacid` | `(seq)` | `Bool` | O(1) |
| `validate_sequence` | `(::Type{A}, s::String)` | `Bool` | O(n) |
| `validate_aa_sequence` | `(s::String)` | `Bool` | O(n) |
| `Rle` | `(data::AbstractVector)` | `Rle{T}` | O(n) |
| `decode` | `(rle::Rle{T})` | `Vector{T}` | O(n) |
| `assay` | `(se, name::String)` | `Matrix{Float64}` | O(1) |
| `assay` | `(se)` | `Matrix{Float64}` | O(1) |
| `rowData` | `(se)` | `Dict{Symbol,Vector}` | O(1) |
| `colData` | `(se)` | `Dict{Symbol,Vector}` | O(1) |
| `metadata` | `(se)` | `Dict{Symbol,Any}` | O(1) |
| `subset_features` | `(se, indices)` | `SummarizedExperiment` | O(a * selected rows + row metadata) |
| `subset_samples` | `(se, indices)` | `SummarizedExperiment` | O(a * selected columns + col metadata) |
| `insert!` | `(tree, left, right, payload)` | `IntervalTree` | O(log n) |
| `query_overlaps` | `(tree, ql, qr)` | `Vector{T}` | O(log n + k) |
| `length` | `(tree::IntervalTree)` | `Int` | O(n) |
| `isempty` | `(tree::IntervalTree)` | `Bool` | O(1) |

---

## 13. Complete Usage Example

```julia
using BioToolkit

# ---- Biological Sequences ----

dna = DNASeq("ATCGATCGATCG")
rna = RNASeq("AUCGAUCGAUCG")
protein = AASeq("MVLSPADKTNVKAAWGKVGA")

dna[1]                    # 'A'
dna[3:6]                  # DNASeq("CGAT")
String(protein)           # "MVLSPADKTNVKAAWGKVGA"

for sym in dna
    println(sym)
end

isdna(dna)                # true
isrna(dna)                # false
isaminoacid(protein)      # true
alphabet(rna)             # RNAAlphabet

validate_sequence(DNAAlphabet, "ATCGX")       # false
validate_sequence(AminoAcidAlphabet, "MKL*")  # true

# ---- Run-Length Encoding ----

coverage = [0,0,0,5,5,5,5,0,0,3,3,3,3,3]
rle = Rle(coverage)

rle.values                # [0, 5, 0, 3]
rle.lengths               # [3, 4, 2, 5]
rle[5]                    # 5
decode(rle)               # original vector

# ---- SummarizedExperiment ----

counts = rand(1:100, 5, 3)  # 5 genes x 3 samples
se = SummarizedExperiment(
    counts;
    rowData = Dict(:gene => ["BRCA1", "TP53", "MYC", "EGFR", "KRAS"]),
    colData = Dict(:condition => ["control", "treatment", "treatment"]),
    metadata = Dict(:organism => "human")
)

assay(se, "counts")
rowData(se)[:gene]
colData(se)[:condition]

se_small = subset_features(se, 1:3)
se_treated = subset_samples(se, [2, 3])

# ---- Interval Tree ----

tree = IntervalTree{String}()
insert!(tree, 100, 500, "gene1")
insert!(tree, 200, 800, "gene2")
insert!(tree, 50, 150, "gene3")

query_overlaps(tree, 120, 300)  # overlaps gene3, gene1, gene2
query_overlaps(tree, 600, 700)  # overlaps gene2
query_overlaps(tree, 10, 20)    # empty result
```
