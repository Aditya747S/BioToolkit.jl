# ==============================================================================
# sequence.jl — Nucleotide sequence operations
#
# Provides high-performance DNA/RNA sequence manipulation: translation,
# complementation, k-mer counting, GC analysis, and FASTA/FASTQ I/O.
#
# Functions use typed BioSequence inputs for biological sequence operations
# (DNASeq/RNASeq/AASeq) for compile-time alphabet safety.
# See biotypes.jl for the type hierarchy.
#
# Lookup tables use 256-element arrays indexed by raw byte value (+1 for
# Julia 1-based indexing). This gives O(1) per-base translation with
# excellent cache locality.
#
# References:
#   - IUPAC nucleotide codes: Cornish-Bowden (1985) NAR 13(9):3021-3030
#   - Codon adaptation index: Sharp & Li (1987) NAR 15(3):1281-1295
#   - Wallace rule for Tm: Wallace et al. (1979) NAR 6(11):3543-3557
#   - Molecular weights: Fasman (1975) Handbook of Biochemistry
# ==============================================================================

using Mmap

@inline function _register_sequence_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

const _DNA_CODE = fill(UInt8(255), 256)
const _DNA_COMPLEMENT = fill(UInt8('N'), 256)
const _CODON_TABLE = fill(UInt8('X'), 64)
const _CODON_AA = fill(UInt8('X'), 64)
const _CODON_TRIPLETS = Vector{String}(undef, 64)
const _DNA_MW_RESIDUE = fill(0.0, 256)
const _IUPAC_DNA_VALID = fill(false, 256)

for (byte, mass) in ((UInt8('A'), 331.2218), (UInt8('C'), 307.1971), (UInt8('G'), 347.2212), (UInt8('T'), 322.2085), (UInt8('U'), 306.17))
    _DNA_MW_RESIDUE[Int(byte) + 1] = mass
    _DNA_MW_RESIDUE[Int(byte | 0x20) + 1] = mass
end

for (byte, mass) in ((UInt8('R'), (331.2218 + 347.2212) / 2), (UInt8('Y'), (307.1971 + 322.2085) / 2), (UInt8('S'), (307.1971 + 347.2212) / 2), (UInt8('W'), (331.2218 + 322.2085) / 2), (UInt8('K'), (347.2212 + 322.2085) / 2), (UInt8('M'), (331.2218 + 307.1971) / 2), (UInt8('B'), (307.1971 + 347.2212 + 322.2085) / 3), (UInt8('D'), (331.2218 + 347.2212 + 322.2085) / 3), (UInt8('H'), (331.2218 + 307.1971 + 322.2085) / 3), (UInt8('V'), (331.2218 + 307.1971 + 347.2212) / 3), (UInt8('N'), (331.2218 + 307.1971 + 347.2212 + 322.2085) / 4))
    _DNA_MW_RESIDUE[Int(byte) + 1] = mass
    _DNA_MW_RESIDUE[Int(byte | 0x20) + 1] = mass
end

for byte in (UInt8('A'), UInt8('a'))
    _DNA_CODE[Int(byte) + 1] = UInt8(0)
    _DNA_COMPLEMENT[Int(byte) + 1] = UInt8('T')
    _IUPAC_DNA_VALID[Int(byte) + 1] = true
end

"""
    FastaIndexRecord

Compact metadata record for random access into an indexed FASTA file.
"""
struct FastaIndexRecord
    name::String
    sequence_length::Int
    byte_offset::Int64
    line_bases::Int
    line_bytes::Int
end

for byte in (UInt8('C'), UInt8('c'))
    _DNA_CODE[Int(byte) + 1] = UInt8(1)
    _DNA_COMPLEMENT[Int(byte) + 1] = UInt8('G')
    _IUPAC_DNA_VALID[Int(byte) + 1] = true
end
for byte in (UInt8('G'), UInt8('g'))
    _DNA_CODE[Int(byte) + 1] = UInt8(2)
    _DNA_COMPLEMENT[Int(byte) + 1] = UInt8('C')
    _IUPAC_DNA_VALID[Int(byte) + 1] = true
end
for byte in (UInt8('T'), UInt8('t'), UInt8('U'), UInt8('u'))
    _DNA_CODE[Int(byte) + 1] = UInt8(3)
    _DNA_COMPLEMENT[Int(byte) + 1] = UInt8('A')
    _IUPAC_DNA_VALID[Int(byte) + 1] = true
end
for byte in (UInt8('N'), UInt8('n'))
    _DNA_CODE[Int(byte) + 1] = UInt8(4)
    _DNA_COMPLEMENT[Int(byte) + 1] = UInt8('N')
    _IUPAC_DNA_VALID[Int(byte) + 1] = true
end

for (byte, complement) in ((UInt8('R'), UInt8('Y')), (UInt8('Y'), UInt8('R')), (UInt8('S'), UInt8('S')), (UInt8('W'), UInt8('W')), (UInt8('K'), UInt8('M')), (UInt8('M'), UInt8('K')), (UInt8('B'), UInt8('V')), (UInt8('V'), UInt8('B')), (UInt8('D'), UInt8('H')), (UInt8('H'), UInt8('D')))
    _DNA_CODE[Int(byte) + 1] = UInt8(4)
    _DNA_COMPLEMENT[Int(byte) + 1] = complement
    _IUPAC_DNA_VALID[Int(byte) + 1] = true
    lower = UInt8(byte | 0x20)
    _DNA_CODE[Int(lower) + 1] = UInt8(4)
    _DNA_COMPLEMENT[Int(lower) + 1] = UInt8(complement | 0x20)
    _IUPAC_DNA_VALID[Int(lower) + 1] = true
end

for (codon, amino_acid) in (
    ("ATA", 'I'), ("ATC", 'I'), ("ATT", 'I'), ("ATG", 'M'),
    ("ACA", 'T'), ("ACC", 'T'), ("ACG", 'T'), ("ACT", 'T'),
    ("AAC", 'N'), ("AAT", 'N'), ("AAA", 'K'), ("AAG", 'K'),
    ("AGC", 'S'), ("AGT", 'S'), ("AGA", 'R'), ("AGG", 'R'),
    ("CTA", 'L'), ("CTC", 'L'), ("CTG", 'L'), ("CTT", 'L'),
    ("CCA", 'P'), ("CCC", 'P'), ("CCG", 'P'), ("CCT", 'P'),
    ("CAC", 'H'), ("CAT", 'H'), ("CAA", 'Q'), ("CAG", 'Q'),
    ("CGA", 'R'), ("CGC", 'R'), ("CGG", 'R'), ("CGT", 'R'),
    ("GTA", 'V'), ("GTC", 'V'), ("GTG", 'V'), ("GTT", 'V'),
    ("GCA", 'A'), ("GCC", 'A'), ("GCG", 'A'), ("GCT", 'A'),
    ("GAC", 'D'), ("GAT", 'D'), ("GAA", 'E'), ("GAG", 'E'),
    ("GGA", 'G'), ("GGC", 'G'), ("GGG", 'G'), ("GGT", 'G'),
    ("TCA", 'S'), ("TCC", 'S'), ("TCG", 'S'), ("TCT", 'S'),
    ("TTC", 'F'), ("TTT", 'F'), ("TTA", 'L'), ("TTG", 'L'),
    ("TAC", 'Y'), ("TAT", 'Y'), ("TGC", 'C'), ("TGT", 'C'),
    ("TGG", 'W'),
    ("TAA", '*'), ("TAG", '*'), ("TGA", '*'))
    b1 = UInt8(codon[1])
    b2 = UInt8(codon[2])
    b3 = UInt8(codon[3])
    code1 = Int(_DNA_CODE[Int(b1) + 1])
    code2 = Int(_DNA_CODE[Int(b2) + 1])
    code3 = Int(_DNA_CODE[Int(b3) + 1])
    # Codon index: 6-bit value from three 2-bit base codes, plus 1 for Julia indexing.
    # Layout: [code1:2bits][code2:2bits][code3:2bits] → 0..63, stored at positions 1..64.
    codon_index = ((code1 << 4) | (code2 << 2) | code3) + 1
    _CODON_TABLE[codon_index] = UInt8(amino_acid)
    _CODON_AA[codon_index] = UInt8(amino_acid)
    _CODON_TRIPLETS[codon_index] = codon
end

# GPU acceleration is available when the CUDA extension is loaded;
# see lazy_gpu.jl for the dispatch mechanism.

"""
    _dna_code_byte(byte)

Map a nucleotide byte to its compact internal DNA code.
"""
@inline function _dna_code_byte(byte::UInt8)
    if byte == UInt8('A') || byte == UInt8('a')
        return UInt8(0)
    elseif byte == UInt8('C') || byte == UInt8('c')
        return UInt8(1)
    elseif byte == UInt8('G') || byte == UInt8('g')
        return UInt8(2)
    elseif byte == UInt8('T') || byte == UInt8('t') || byte == UInt8('U') || byte == UInt8('u')
        return UInt8(3)
    elseif byte == UInt8('N') || byte == UInt8('n')
        return UInt8(4)
    end

    return UInt8(255)
end

"""
    _kmer_base_from_code(code)

Map a compact k-mer code back to a nucleotide byte.
"""
@inline function _kmer_base_from_code(code::Int)
    code == 0 && return UInt8('A')
    code == 1 && return UInt8('C')
    code == 2 && return UInt8('G')
    code == 3 && return UInt8('T')
    return UInt8('N')
end

"""
    _decode_kmer_key(key, k)

Decode an integer k-mer key into its nucleotide string representation.
"""
function _decode_kmer_key(key::Int, k::Int)
    bytes = Vector{UInt8}(undef, k)
    value = key - 1

    @inbounds for index in k:-1:1
        value, code = divrem(value, 5)
        bytes[index] = _kmer_base_from_code(code)
    end

    return String(bytes)
end

"""
    validate_dna(sequence) -> Bool

Return `true` for typed DNA/RNA sequences and `false` for non-nucleotide
typed alphabets.
"""
validate_dna(::BioSequence{DNAAlphabet}) = true
validate_dna(::BioSequence{RNAAlphabet}) = true
validate_dna(::BioSequence{AminoAcidAlphabet}) = false
validate_dna(sequence::String) = isempty(sequence) || all(byte -> _IUPAC_DNA_VALID[Int(byte) + 1], codeunits(sequence))

@inline _sequence_to_dna(sequence::AbstractString) = DNASeq(String(sequence); validate=false)
@inline _sequence_to_aa(sequence::AbstractString) = AASeq(String(sequence); validate=false)

# Removed validate_dna(::AbstractString) - use BioSequence types instead

"""
    FastqRecordStream

Iterate over FASTQ records from an input stream.
"""
mutable struct FastqRecordStream{S}
    io::IO
    close_on_finish::Bool
    finished::Bool
    sequence_type::Type{S}
end

Base.IteratorSize(::Type{<:FastqRecordStream}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:FastqRecordStream}) = Base.HasEltype()

_fastq_record_alphabet(::Type{<:BioSequence{A}}) where {A <: BioAlphabet} = A
_fastq_record_alphabet(::Type) = DNAAlphabet

_fastq_sequence(::Type{S}, text::String) where {S<:BioSequence} = S(text)

Base.eltype(::Type{FastqRecordStream{S}}) where {S} = FastqRecord{_fastq_record_alphabet(S)}

"""
    Base.close(stream)

Close a FASTQ record stream and mark it as finished.
"""
function Base.close(stream::FastqRecordStream)
    stream.finished = true
    close(stream.io)
    return nothing
end

"""
    Base.iterate(stream)

Iterate over the next FASTQ record from a stream.
"""
function Base.iterate(stream::FastqRecordStream{S}) where {S}
    stream.finished && return nothing
    record = _read_fastq_record(stream.io, stream.sequence_type)
    if record === nothing
        stream.close_on_finish && !stream.finished && close(stream)
        return nothing
    end
    return record, nothing
end

"""
    Base.iterate(stream, state)

Continue FASTQ stream iteration.
"""
function Base.iterate(stream::FastqRecordStream{S}, ::Nothing) where {S}
    return iterate(stream)
end

"""
    each_fastq_record(path; sequence_type=DNASeq)

Create a FASTQ record iterator for a file path.
"""
function each_fastq_record(path::String; sequence_type::Type{S}=DNASeq) where {S <: BioSequence}
    return FastqRecordStream{S}(open(path, "r"), true, false, sequence_type)
end

"""
    each_fastq_record(io; sequence_type=DNASeq)

Create a FASTQ record iterator for an open IO stream.
"""
function each_fastq_record(io::IO; sequence_type::Type{S}=DNASeq) where {S <: BioSequence}
    return FastqRecordStream{S}(io, false, false, sequence_type)
end

"""
    _read_fastq_record(io)

Read a single FASTQ record from an IO stream.
"""
function _read_fastq_record(io::IO, sequence_type::Type{S}=DNASeq) where {S <: BioSequence}
    while !eof(io)
        header_line = strip(readline(io))
        isempty(header_line) && continue
        startswith(header_line, '@') || throw(ArgumentError("expected FASTQ header line starting with @"))

        description = strip(header_line[2:end])
        separator = findfirst(isspace, description)
        identifier = separator === nothing ? description : description[1:prevind(description, separator)]

        sequence_buffer = IOBuffer()
        while true
            eof(io) && throw(ArgumentError("truncated FASTQ record"))
            line = rstrip(readline(io))
            startswith(line, '+') && break
            write(sequence_buffer, line)
        end

        sequence = _fastq_sequence(sequence_type, uppercase(String(take!(sequence_buffer))))
        quality_buffer = IOBuffer()
        quality_length = 0

        while quality_length < length(sequence)
            eof(io) && throw(ArgumentError("truncated FASTQ quality block"))
            line = rstrip(readline(io))
            isempty(line) && continue
            write(quality_buffer, line)
            quality_length += ncodeunits(line)
        end

        quality = String(take!(quality_buffer))
        quality_length == length(sequence) || throw(ArgumentError("FASTQ sequence and quality lengths differ"))

        return FastqRecord(String(identifier), String(description), sequence, quality)
    end

    return nothing
end

"""
    count_nucleotides(sequence) -> NamedTuple{(:A,:C,:G,:T,:N), NTuple{5,Int}}

Count canonical nucleotides (A, C, G, T) and ambiguous/non-canonical
bases (N) in a nucleotide sequence.  Case-insensitive; IUPAC ambiguity
codes (R, Y, S, W, K, M, B, D, H, V) are counted as N.

Accepts typed `BioSequence` inputs.
"""
function count_nucleotides(sequence::BioSequence{A}) where {A <: BioAlphabet}
    a = 0; c = 0; g = 0; t = 0; n = 0

    @inbounds for byte in sequence.data
        code = _DNA_CODE[Int(byte) + 1]
        if code == 0
            a += 1
        elseif code == 1
            c += 1
        elseif code == 2
            g += 1
        elseif code == 3
            t += 1
        else
            n += 1
        end
    end

    return (A = a, C = c, G = g, T = t, N = n)
end

count_nucleotides(sequence::AbstractString) = count_nucleotides(_sequence_to_dna(sequence))

"""
    transcribe_dna(sequence)

Transcribe DNA to RNA by replacing thymine with uracil.
"""
function transcribe_dna(sequence::BioSequence{DNAAlphabet})
    length_sequence = length(sequence)
    length_sequence == 0 && return RNASeq(UInt8[]; validate=false)

    buffer = Vector{UInt8}(undef, length_sequence)

    @inbounds for index in eachindex(sequence.data)
        byte = sequence.data[index]
        if byte == UInt8('A')
            buffer[index] = UInt8('A')
        elseif byte == UInt8('C')
            buffer[index] = UInt8('C')
        elseif byte == UInt8('G')
            buffer[index] = UInt8('G')
        elseif byte == UInt8('T') || byte == UInt8('U')
            buffer[index] = UInt8('U')
        elseif byte == UInt8('N')
            buffer[index] = UInt8('N')
        else
            buffer[index] = byte
            if buffer[index] == UInt8('T')
                buffer[index] = UInt8('U')
            end
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, RNASeq(buffer; validate=false), "transcribe_dna")
end

transcribe_dna(sequence::AbstractString) = String(transcribe_dna(_sequence_to_dna(sequence)))

"""
    reverse_complement(sequence)

Return the reverse complement of a nucleotide string.  Uses the IUPAC
complement table: A↔T, C↔G, R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N.

Returns a typed `BioSequence`.
"""
function reverse_complement(sequence::BioSequence{DNAAlphabet})
    isempty(sequence) && return DNASeq(UInt8[]; validate=false)

    length_sequence = length(sequence)
    buffer = Vector{UInt8}(undef, length_sequence)

    complement = _DNA_COMPLEMENT
    source_index = length_sequence

    @inbounds for destination_index in 1:length_sequence
        buffer[destination_index] = complement[Int(sequence.data[source_index]) + 1]
        source_index -= 1
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, DNASeq(buffer; validate=false), "reverse_complement")
end

function reverse_complement(seq::BioSequence{RNAAlphabet})
    # Complement then swap T→U in result
    rc = reverse_complement(DNASeq(seq.data; validate=false))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, RNASeq(UInt8[b == UInt8('T') ? UInt8('U') : b for b in rc.data]; validate=false), "reverse_complement")
end

reverse_complement(sequence::AbstractString) = reverse_complement(codeunits(String(sequence)))

"""
    reverse_complement(bytes)

Return the reverse complement of a nucleotide byte vector.
"""
function reverse_complement(bytes::AbstractVector{UInt8})
    length_sequence = length(bytes)
    length_sequence == 0 && return ""

    buffer = Vector{UInt8}(undef, length_sequence)

    reverse_complement!(buffer, bytes)

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, unsafe_string(pointer(buffer), length_sequence), "reverse_complement")
end

"""
    reverse_complement!(buffer, bytes)

Write the reverse complement of `bytes` into `buffer`.
"""
function reverse_complement!(buffer::AbstractVector{UInt8}, bytes::AbstractVector{UInt8})
    length_sequence = length(bytes)
    length(buffer) >= length_sequence || throw(ArgumentError("buffer is too small"))

    complement = _DNA_COMPLEMENT
    source_index = length_sequence

    @inbounds for destination_index in 1:length_sequence
        buffer[destination_index] = complement[Int(bytes[source_index]) + 1]
        source_index -= 1
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, length_sequence, "reverse_complement!")
end

# Note: the Vector{UInt8} overload is intentionally removed — it was an exact duplicate
# of the AbstractVector{UInt8} overload above and caused method ambiguity warnings.
# Julia's dispatch already selects the concrete method when called with Vector{UInt8}.

"""
    gc_content(sequence) -> Float64

Compute the GC fraction of a nucleotide sequence.  Only unambiguous bases
(A, C, G, T/U) contribute to the denominator; the IUPAC ambiguity code
`S` (strong: C or G) is counted as GC.  Other ambiguity codes are skipped.

Returns 0.0 for empty sequences or sequences with no canonical bases.

    Accepts typed `BioSequence` inputs.
"""
function gc_content(sequence::BioSequence{A}) where {A <: BioAlphabet}
    isempty(sequence) && return 0.0

    gc = 0
    total = 0

    @inbounds for byte in sequence.data
        if byte == UInt8('C') || byte == UInt8('c') || byte == UInt8('G') || byte == UInt8('g') || byte == UInt8('S') || byte == UInt8('s')
            gc += 1
            total += 1
        elseif byte == UInt8('A') || byte == UInt8('a') || byte == UInt8('T') || byte == UInt8('t') || byte == UInt8('U') || byte == UInt8('u')
            total += 1
        elseif !_IUPAC_DNA_VALID[Int(byte) + 1]
            throw(ArgumentError("unsupported base '$(Char(byte))' in gc_content"))
        end
    end

    total == 0 && return 0.0
    gc_val = gc / total
    _ctx = active_provenance_context()
    return _register_sequence_result!(_ctx, gc_val, "gc_content"; parameters=(n_bases=total, gc_fraction=gc_val))
end

gc_content(sequence::AbstractString) = gc_content(_sequence_to_dna(sequence))

"""
    melting_temp(sequence; rna=false, method=:wallace)

Estimate oligonucleotide melting temperature in degrees Celsius.

The default `:wallace` method uses the classic 2/4 rule:
`2 * (A + T/U) + 4 * (G + C)`.
`rna=true` treats U as the thymine equivalent, but DNA sequences with U are
also accepted and counted as thymine-like bases.

`method=:basic` uses the longer-oligo approximation
`64.9 + 41 * (G + C - 16.4) / N` where `N` is the number of canonical bases.
"""
function melting_temp(sequence::BioSequence{A}; rna::Bool=false, method::Symbol=:wallace) where {A <: BioAlphabet}
    isempty(sequence) && return 0.0

    a = 0
    c = 0
    g = 0
    t = 0
    u = 0

    @inbounds for byte in sequence.data
        ch = byte | 0x20
        if ch == UInt8('a')
            a += 1
        elseif ch == UInt8('c')
            c += 1
        elseif ch == UInt8('g')
            g += 1
        elseif ch == UInt8('t')
            t += 1
        elseif ch == UInt8('u')
            u += 1
        elseif ch == UInt8('n')
            # Ambiguous bases are ignored in the wallace/basic estimate.
        else
            throw(ArgumentError("unsupported base '$(Char(byte))' in melting_temp"))
        end
    end

    method === :wallace || method === :basic || throw(ArgumentError("unsupported melting temperature method: $method"))

    if method === :wallace
        return 2.0 * (a + t + u) + 4.0 * (g + c)
    end

    canonical_bases = a + c + g + t + u
    canonical_bases == 0 && return 0.0
    return 64.9 + 41.0 * ((g + c) - 16.4) / canonical_bases
end

melting_temp(sequence::AbstractString; rna::Bool=false, method::Symbol=:wallace) = melting_temp(_sequence_to_dna(sequence); rna=rna, method=method)

"""
    dna_molecular_weight(sequence; stranded=:single)

Estimate the molecular weight of a DNA or RNA sequence in Daltons.
`stranded=:double` returns the weight of both strands combined.
Ambiguous IUPAC bases are averaged across their canonical possibilities.
"""
function dna_molecular_weight(sequence::BioSequence{A}; stranded::Symbol=:single) where {A <: BioAlphabet}
    stranded === :single || stranded === :double || throw(ArgumentError("stranded must be :single or :double"))

    weight = _dna_molecular_weight_single(sequence)
    stranded === :single && return weight
    return weight + _dna_molecular_weight_single(reverse_complement(sequence))
end

dna_molecular_weight(sequence::AbstractString; stranded::Symbol=:single) = dna_molecular_weight(_sequence_to_dna(sequence); stranded=stranded)

"""
    _dna_molecular_weight_single(sequence)

Compute the molecular weight of a single nucleic acid strand.
"""
function _dna_molecular_weight_single(sequence::BioSequence{A}) where {A <: BioAlphabet}
    length_sequence = length(sequence)
    length_sequence == 0 && return 0.0

    residue_mass = 0.0
    residue_count = 0

    @inbounds for byte in sequence.data
        mass = _DNA_MW_RESIDUE[Int(byte) + 1]
        mass == 0.0 && throw(ArgumentError("unsupported base '$(Char(byte))' in dna_molecular_weight"))
        residue_mass += mass
        residue_count += 1
    end

    residue_count == 0 && return 0.0
    return residue_mass - (residue_count - 1) * 18.01528
end

"""
    _codon_index(byte1, byte2, byte3)

Map three nucleotide bytes to a compact codon index.
"""
@inline function _codon_index(byte1::UInt8, byte2::UInt8, byte3::UInt8)
    code1 = _DNA_CODE[Int(byte1) + 1]
    code2 = _DNA_CODE[Int(byte2) + 1]
    code3 = _DNA_CODE[Int(byte3) + 1]
    (code1 > 3 || code2 > 3 || code3 > 3) && return 0
    return ((Int(code1) << 4) | (Int(code2) << 2) | Int(code3)) + 1
end

"""
    codon_usage(sequence; frame=1, include_stop=true)

Count codon usage in a single sequence.
"""
function codon_usage(sequence::BioSequence{DNAAlphabet}; frame::Integer=1, include_stop::Bool=true)
    1 <= frame <= 3 || throw(ArgumentError("frame must be 1, 2, or 3"))
    length_sequence = length(sequence)
    frame <= length_sequence || return Dict{String,Int}()

    counts = Dict{String,Int}()
    bytes = sequence.data

    @inbounds for index in frame:3:(length_sequence - 2)
        codon_index = _codon_index(bytes[index], bytes[index + 1], bytes[index + 2])
        codon_index == 0 && continue
        amino_acid = _CODON_AA[codon_index]
        !include_stop && amino_acid == UInt8('*') && continue
        codon = _CODON_TRIPLETS[codon_index]
        counts[codon] = get(counts, codon, 0) + 1
    end

    result = counts
    _ctx = active_provenance_context()
    return _register_sequence_result!(_ctx, result, "codon_usage"; parameters=(n_codons=sum(values(counts)), frame=Int(frame)))
end

"""
    codon_usage(sequences; frame=1, include_stop=true)

Count codon usage across a collection of sequences.
"""
function codon_usage(sequence::AbstractVector{<:BioSequence{DNAAlphabet}}; frame::Integer=1, include_stop::Bool=true)
    counts = Dict{String,Int}()
    for item in sequence
        for (codon, count) in codon_usage(item; frame=frame, include_stop=include_stop)
            counts[codon] = get(counts, codon, 0) + count
        end
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, counts, "codon_usage")
end

codon_usage(sequence::AbstractString; frame::Integer=1, include_stop::Bool=true) = codon_usage(_sequence_to_dna(sequence); frame=frame, include_stop=include_stop)

# Removed codon_usage(::AbstractVector{<:AbstractString}) - use Vector{DNASeq} instead

"""
    codon_usage_table(sequence; frame=1, include_stop=false)

Return codon usage as a table for a single sequence.
"""
function codon_usage_table(sequence::BioSequence{DNAAlphabet}; frame::Integer=1, include_stop::Bool=false)
    counts = codon_usage(sequence; frame=frame, include_stop=include_stop)
    total = sum(values(counts))
    total == 0 && return Dict{String,Float64}()

    table = Dict{String,Float64}()
    for (codon, count) in counts
        table[codon] = count / total
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, table, "codon_usage_table")
end

"""
    codon_usage_table(sequences; frame=1, include_stop=false)

Return codon usage as a table for multiple sequences.
"""
function codon_usage_table(sequences::AbstractVector{<:BioSequence{DNAAlphabet}}; frame::Integer=1, include_stop::Bool=false)
    counts = codon_usage(sequences; frame=frame, include_stop=include_stop)
    total = sum(values(counts))
    total == 0 && return Dict{String,Float64}()

    table = Dict{String,Float64}()
    for (codon, count) in counts
        table[codon] = count / total
    end
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, table, "codon_usage_table")
end

codon_usage_table(sequence::AbstractString; frame::Integer=1, include_stop::Bool=false) = codon_usage_table(_sequence_to_dna(sequence); frame=frame, include_stop=include_stop)

# Removed codon_usage_table(::AbstractVector{<:AbstractString}) - use Vector{DNASeq} instead

"""
    relative_codon_adaptiveness(reference; frame=1, pseudocount=0.5)

Compute relative codon adaptiveness values from a reference set.
"""
function relative_codon_adaptiveness(reference; frame::Integer=1, pseudocount::Real=0.5)
    pseudocount >= 0 || throw(ArgumentError("pseudocount must be nonnegative"))
    counts = codon_usage(reference; frame=frame, include_stop=false)
    weights = Dict{String,Float64}()
    amino_acid_max = Dict{Char,Float64}()

    for index in 1:64
        amino_acid = Char(_CODON_AA[index])
        amino_acid == '*' && continue
        codon = _CODON_TRIPLETS[index]
        count = get(counts, codon, 0) + pseudocount
        amino_acid_max[amino_acid] = max(get(amino_acid_max, amino_acid, 0.0), Float64(count))
    end

    for index in 1:64
        amino_acid = Char(_CODON_AA[index])
        amino_acid == '*' && continue
        codon = _CODON_TRIPLETS[index]
        denominator = get(amino_acid_max, amino_acid, 0.0)
        weights[codon] = denominator == 0.0 ? 0.0 : (get(counts, codon, 0) + pseudocount) / denominator
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, weights, "relative_codon_adaptiveness")
end

"""
    codon_adaptation_index(sequence; reference, frame=1, pseudocount=0.5)

Compute the codon adaptation index for a coding sequence.
"""
function codon_adaptation_index(sequence::BioSequence{DNAAlphabet}; reference, frame::Integer=1, pseudocount::Real=0.5)
    weights = relative_codon_adaptiveness(reference; frame=frame, pseudocount=pseudocount)
    length_sequence = length(sequence)
    length_sequence < 3 && return 0.0

    log_sum = 0.0
    used = 0
    bytes = sequence.data

    @inbounds for index in frame:3:(length_sequence - 2)
        codon_index = _codon_index(bytes[index], bytes[index + 1], bytes[index + 2])
        codon_index == 0 && continue
        amino_acid = _CODON_AA[codon_index]
        amino_acid == UInt8('*') && continue
        codon = _CODON_TRIPLETS[codon_index]
        weight = get(weights, codon, 0.0)
        weight <= 0 && return 0.0
        log_sum += log(weight)
        used += 1
    end

    used == 0 && return 0.0
    cai_val = exp(log_sum / used)
    _ctx = active_provenance_context()
    return _register_sequence_result!(_ctx, cai_val, "codon_adaptation_index"; parameters=(n_codons=used, frame=Int(frame), cai=cai_val))
end

codon_adaptation_index(sequence::AbstractString; reference, frame::Integer=1, pseudocount::Real=0.5) =
    codon_adaptation_index(_sequence_to_dna(sequence); reference=reference, frame=frame, pseudocount=pseudocount)

cai(sequence::BioSequence{DNAAlphabet}; reference, kwargs...) = codon_adaptation_index(sequence; reference=reference, kwargs...)
cai(sequence::AbstractString; reference, kwargs...) = codon_adaptation_index(sequence; reference=reference, kwargs...)

"""
    translate_dna(sequence; stop_at_stop=false)

Translate a DNA/RNA coding sequence into a single-letter amino acid string
using the standard genetic code (NCBI translation table 1).
using ..BioToolkit: ProvenanceParams, ThreadSafeProvenanceContext, new_provenance_id

The input length must be a multiple of 3.  Stop codons (TAA, TAG, TGA)
produce `'*'`.  When `stop_at_stop=true`, translation terminates at the
first stop codon; otherwise stop codons are omitted from the output.

Accepts typed `BioSequence{DNAAlphabet}` inputs and returns
`BioSequence{AminoAcidAlphabet}`.
"""
function translate_dna(sequence::BioSequence{DNAAlphabet}; stop_at_stop::Bool=false, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    length_sequence = length(sequence)
    length_sequence % 3 == 0 || throw(ArgumentError("DNA sequence length must be a multiple of 3 for translate_dna; got $(length_sequence)"))
    buffer = Vector{UInt8}(undef, length_sequence ÷ 3)
    written = translate_dna!(buffer, sequence.data; stop_at_stop=stop_at_stop)
    resize!(buffer, written)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, AASeq(buffer; validate=false), "translate_dna")
end

"""
    _translate_dna_bytes!(buffer, source, last; stop_at_stop=false)

Translate encoded DNA bytes into an amino-acid byte buffer.
"""
@inline function _translate_dna_bytes!(buffer::Vector{UInt8}, source, last::Int; stop_at_stop::Bool=false)
    code_table = _DNA_CODE
    codon_table = _CODON_TABLE

    index = 1
    out_index = 1

    # Codon index: 6-bit value from three 2-bit base codes + 1 for Julia indexing.
    # Ambiguous bases (code > 3) produce 'X' (unknown amino acid).
    @inbounds while index + 2 <= last
        code1 = code_table[Int(unsafe_load(source, index)) + 1]
        code2 = code_table[Int(unsafe_load(source, index + 1)) + 1]
        code3 = code_table[Int(unsafe_load(source, index + 2)) + 1]

        amino_acid = (code1 > 3 || code2 > 3 || code3 > 3) ? UInt8('X') : codon_table[((Int(code1) << 4) | (Int(code2) << 2) | Int(code3)) + 1]
        if stop_at_stop && amino_acid == UInt8('*')
            break
        elseif amino_acid != UInt8('*')
            buffer[out_index] = amino_acid
            out_index += 1
        end
        index += 3
    end

    return out_index - 1
end

"""
    translate_dna(bytes; stop_at_stop=false)

Translate a nucleotide byte vector into an amino-acid string.
"""
function translate_dna(bytes::AbstractVector{UInt8}; stop_at_stop::Bool=false)
    last = length(bytes)
    last % 3 == 0 || throw(ArgumentError("DNA sequence length must be a multiple of 3 for translate_dna; got $(last)"))
    buffer = Vector{UInt8}(undef, last ÷ 3)
    written = translate_dna!(buffer, bytes; stop_at_stop=stop_at_stop)
    resize!(buffer, written)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, unsafe_string(pointer(buffer), written), "translate_dna")
end

translate_dna(sequence::AbstractString; stop_at_stop::Bool=false) = translate_dna(codeunits(String(sequence)); stop_at_stop=stop_at_stop)

"""
    translate_dna!(buffer, bytes; stop_at_stop=false)

Translate nucleotide bytes into a preallocated amino-acid buffer.
"""
function translate_dna!(buffer::Vector{UInt8}, bytes::AbstractVector{UInt8}; stop_at_stop::Bool=false)
    last = length(bytes)
    last % 3 == 0 || throw(ArgumentError("DNA sequence length must be a multiple of 3 for translate_dna; got $(last)"))
    length(buffer) >= last ÷ 3 || throw(ArgumentError("buffer is too small"))
    code_table = _DNA_CODE
    codon_table = _CODON_TABLE

    index = 1
    out_index = 1

    @inbounds while index + 2 <= last
        code1 = code_table[Int(bytes[index]) + 1]
        code2 = code_table[Int(bytes[index + 1]) + 1]
        code3 = code_table[Int(bytes[index + 2]) + 1]

        amino_acid = (code1 > 3 || code2 > 3 || code3 > 3) ? UInt8('X') : codon_table[((Int(code1) << 4) | (Int(code2) << 2) | Int(code3)) + 1]
        if stop_at_stop && amino_acid == UInt8('*')
            break
        elseif amino_acid != UInt8('*')
            buffer[out_index] = amino_acid
            out_index += 1
        end
        index += 3
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, out_index - 1, "translate_dna!")
end

"""
    translate_dna!(buffer, bytes; stop_at_stop=false)

Translate nucleotide bytes into a preallocated amino-acid buffer.
"""
function translate_dna!(buffer::AbstractVector{UInt8}, bytes::AbstractVector{UInt8}; stop_at_stop::Bool=false, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    last = length(bytes)
    last % 3 == 0 || throw(ArgumentError("DNA sequence length must be a multiple of 3 for translate_dna; got $(last)"))
    length(buffer) >= last ÷ 3 || throw(ArgumentError("buffer is too small"))
    code_table = _DNA_CODE
    codon_table = _CODON_TABLE

    index = 1
    out_index = 1

    # When stop_at_stop=true, halt at the first stop codon without emitting it.
    # When stop_at_stop=false, skip stop codons and continue translating.
    @inbounds while index + 2 <= last
        code1 = code_table[Int(bytes[index]) + 1]
        code2 = code_table[Int(bytes[index + 1]) + 1]
        code3 = code_table[Int(bytes[index + 2]) + 1]

        amino_acid = (code1 > 3 || code2 > 3 || code3 > 3) ? UInt8('X') : codon_table[((Int(code1) << 4) | (Int(code2) << 2) | Int(code3)) + 1]
        if stop_at_stop && amino_acid == UInt8('*')
            break
        elseif amino_acid != UInt8('*')
            buffer[out_index] = amino_acid
            out_index += 1
        end

        index += 3
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, out_index - 1, "translate_dna!")
end

"""
    translate_dna!(buffer, bytes; stop_at_stop=false)

Translate nucleotide bytes into a preallocated amino-acid buffer.
"""
function translate_dna!(buffer::Vector{UInt8}, bytes::Vector{UInt8}; stop_at_stop::Bool=false, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    last = length(bytes)
    length(buffer) >= last ÷ 3 || throw(ArgumentError("buffer is too small"))
    GC.@preserve bytes begin

        return _translate_dna_bytes!(buffer, pointer(bytes), last; stop_at_stop=stop_at_stop)
    end
end

"""
    hamming_distance(left, right; use_cuda=false)

Compute the Hamming distance between two nucleotide byte vectors.
"""
function hamming_distance(left::AbstractVector{UInt8}, right::AbstractVector{UInt8}; use_cuda::Bool=false)
    if use_cuda
        _ensure_cuda_sequence!()
        left_cuda = _is_cuda_backed_array(left) ? left : CUDA.CuArray{UInt8}(Vector{UInt8}(left))
        right_cuda = _is_cuda_backed_array(right) ? right : CUDA.CuArray{UInt8}(Vector{UInt8}(right))
        return Base.invokelatest(_CUDA_SEQUENCE_HAMMING_IMPL[], left_cuda, right_cuda)
    end

    length(left) == length(right) || throw(ArgumentError("sequences must have the same length"))

    distance = 0
    left_length = length(left)
    word_count = left_length >>> 3          # how many full 8-byte words
    tail_start = (word_count << 3) + 1      # first byte after the last full word

    GC.@preserve left right begin
        left_words  = Ptr{UInt64}(pointer(left))
        right_words = Ptr{UInt64}(pointer(right))

        # Performance: xor two 64-bit words then count non-zero bytes via popcount.
        # A non-zero byte contributes at least one bit to xor; we need unequal-byte
        # count, so we detect each differing byte separately using the standard
        # "SWAR" (SIMD Within A Register) parallel trick:
        #   x = xor(L, R)
        #   a byte is non-zero iff:  (x - 0x0101...) & ~x & 0x8080... != 0
        # which Julia expresses as count_nonzero_bytes below.
        @inbounds for word_index in 1:word_count
            l = unsafe_load(left_words,  word_index)
            r = unsafe_load(right_words, word_index)
            xor_word = l ⊻ r
            # SWAR: detect ZERO bytes in xor_word (zero = matching bytes)
            # then subtract from 8 to get MISMATCH count
            zero_mask = (xor_word - 0x0101010101010101) & ~xor_word & 0x8080808080808080
            distance += 8 - count_ones(zero_mask)
        end

        @inbounds for index in tail_start:left_length
            distance += left[index] == right[index] ? 0 : 1
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, distance, "hamming_distance")
end

# Removed hamming_distance(::AbstractString, ::AbstractString) - use BioSequence types instead
function hamming_distance(left::AbstractString, right::AbstractString; use_cuda::Bool=false)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, hamming_distance(DNASeq(left; validate=false), DNASeq(right; validate=false); use_cuda=use_cuda), "hamming_distance")
end

"""
    find_orfs(sequence; min_aa=0)

Find open reading frames in a nucleotide sequence.
"""
function find_orfs(sequence::BioSequence{DNAAlphabet}; min_aa::Integer=0)
    bytes = sequence.data
    proteins = AASeq[]

    for frame in 1:3
        window_length = length(bytes) - frame + 1
        window_length <= 0 && continue
        complete_length = window_length - (window_length % 3)
        complete_length == 0 && continue
        frame_bytes = bytes[frame:frame + complete_length - 1]

        protein = AASeq(translate_dna(frame_bytes; stop_at_stop=false); validate=false)
        if length(protein) >= min_aa
            push!(proteins, protein)
        end

        reverse_protein = AASeq(translate_dna(reverse_complement(frame_bytes); stop_at_stop=false); validate=false)
        if length(reverse_protein) >= min_aa
            push!(proteins, reverse_protein)
        end
    end

    _ctx = active_provenance_context()
    return _register_sequence_result!(_ctx, proteins, "find_orfs"; parameters=(n_orfs=length(proteins), min_aa=Int(min_aa), sequence_length=length(sequence.data)))
end

find_orfs(sequence::AbstractString; min_aa::Integer=0) = find_orfs(_sequence_to_dna(sequence); min_aa=min_aa)

"""
    protein_search(sequence, motif)

Search for a protein motif in a translated nucleotide sequence.
"""
function protein_search(sequence::BioSequence{DNAAlphabet}, motif::BioSequence{AminoAcidAlphabet})
    proteins = find_orfs(sequence)
    hits = AASeq[]

    for protein in proteins
        occursin(String(motif), String(protein)) && push!(hits, protein)
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, hits, "protein_search")
end

protein_search(sequence::AbstractString, motif::BioSequence{AminoAcidAlphabet}) = protein_search(_sequence_to_dna(sequence), motif)
protein_search(sequence::BioSequence{DNAAlphabet}, motif::AbstractString) = protein_search(sequence, _sequence_to_aa(motif))
protein_search(sequence::AbstractString, motif::AbstractString) = protein_search(_sequence_to_dna(sequence), _sequence_to_aa(motif))

"""
    kmer_frequency(sequence, k; use_cuda=false)

Count k-mer frequencies in a nucleotide sequence.
"""
function kmer_frequency(sequence::BioSequence{A}, k::Integer; use_cuda::Bool=false) where {A <: BioAlphabet}
    if use_cuda
        bytes = Vector{UInt8}(sequence.data)
        return kmer_frequency(bytes, k; use_cuda=true)
    end

    k <= 0 && throw(ArgumentError("k must be positive"))
    lastindex(sequence) < k && return Dict{String,Int}()

    if k <= 31
        fast_dna = _kmer_frequency_unambiguous_dna(sequence, k)
        fast_dna !== nothing && return fast_dna
    end

    if k <= 55
        fast = _kmer_frequency_dna(sequence, k)
        fast !== nothing && return fast
    end

    counts = Dict{String,Int}()
    upper = String(sequence)

    for index in firstindex(upper):lastindex(upper) - k + 1
        kmer = upper[index:index+k-1]
        counts[kmer] = get(counts, kmer, 0) + 1
    end
    return counts
end

kmer_frequency(sequence::AbstractString, k::Integer; use_cuda::Bool=false) = kmer_frequency(_sequence_to_dna(sequence), k; use_cuda=use_cuda)

const _KMER_FAST_ALPHABET = (UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T'))

"""
    _decode_dna_kmer_key(key, k)

Decode a compact DNA k-mer key into its string form.
"""
@inline function _decode_dna_kmer_key(key::UInt64, k::Int)
    bytes = Vector{UInt8}(undef, k)
    value = key

    @inbounds for index in k:-1:1
        bytes[index] = _KMER_FAST_ALPHABET[Int(value & 0x03) + 1]
        value >>= 2
    end

    return String(bytes)
end

"""
    _kmer_frequency_unambiguous_dna(sequence, k)

Count k-mers for an unambiguous DNA sequence.
"""
function _kmer_frequency_unambiguous_dna(sequence::BioSequence{DNAAlphabet}, k::Integer)
    length_sequence = length(sequence)
    length_sequence < k && return Dict{String,Int}()

    bytes = sequence.data
    code_table = _DNA_CODE
    base_power = UInt64(1)
    for _ in 2:k
        base_power <<= 2
    end

    counts = Dict{UInt64,Int}()
    current_key = UInt64(0)

    @inbounds for index in 1:k
        byte = bytes[index]
        if byte == UInt8('U') || byte == UInt8('u')
            return nothing
        end
        code = code_table[Int(byte) + 1]
        code > 3 && return nothing
        current_key = (current_key << 2) | UInt64(code)
    end

    counts[current_key] = 1

    @inbounds for index in k + 1:length_sequence
        leading_byte = bytes[index - k]
        trailing_byte = bytes[index]
        if leading_byte == UInt8('U') || leading_byte == UInt8('u') || trailing_byte == UInt8('U') || trailing_byte == UInt8('u')
            return nothing
        end

        leading_code = code_table[Int(leading_byte) + 1]
        trailing_code = code_table[Int(trailing_byte) + 1]
        leading_code > 3 && return nothing
        trailing_code > 3 && return nothing

        current_key -= UInt64(leading_code) * base_power
        current_key = (current_key << 2) | UInt64(trailing_code)
        counts[current_key] = get(counts, current_key, 0) + 1
    end

    result = Dict{String,Int}()
    sizehint!(result, length(counts))
    for (key, count) in counts
        result[_decode_dna_kmer_key(key, k)] = count
    end

    return result
end

"""
    _kmer_frequency_dna(sequence, k)

Count k-mers for a DNA sequence that may contain ambiguous symbols.
"""
function _kmer_frequency_dna(sequence::BioSequence{DNAAlphabet}, k::Integer)
    upper = String(sequence)
    bytes = codeunits(upper)
    length(bytes) < k && return Dict{String,Int}()

    code_table = _DNA_CODE
    base_power = UInt128(1)
    for _ in 2:k
        base_power *= UInt128(5)
    end

    counts = Dict{UInt128,Int}()
    representatives = Dict{UInt128,String}()
    current_key = UInt128(0)

    @inbounds for index in 1:k
        code = code_table[Int(bytes[index]) + 1]
        code == 255 && return nothing
        current_key = current_key * 5 + UInt128(code)
    end

    counts[current_key] = 1
    representatives[current_key] = String(upper[firstindex(upper):firstindex(upper) + k - 1])

    @inbounds for index in k + 1:length(bytes)
        leading_code = code_table[Int(bytes[index - k]) + 1]
        trailing_code = code_table[Int(bytes[index]) + 1]
        leading_code == 255 && return nothing
        trailing_code == 255 && return nothing

        current_key -= UInt128(leading_code) * base_power
        current_key = current_key * 5 + UInt128(trailing_code)

        counts[current_key] = get(counts, current_key, 0) + 1
        if !haskey(representatives, current_key)
            start_index = index - k + 1
            representatives[current_key] = String(upper[start_index:index])
        end
    end

    result = Dict{String,Int}()
    sizehint!(result, length(counts))
    for (key, count) in counts
        result[representatives[key]] = count
    end

    return result
end

"""
    fasta_index(path)

Build a compact random-access index for a FASTA file.
"""
function fasta_index(path::String)
    records = Dict{String,FastaIndexRecord}()

    open(path, "r") do io
        current_name = ""
        sequence_length = 0
        byte_offset = 0
        line_bases = 0
        line_bytes = 0
        sequence_start = 0

        while !eof(io)
            line_start = position(io)
            raw_line = readline(io)
            stripped = strip(raw_line)

            if isempty(stripped)
                byte_offset = position(io)
                continue
            elseif startswith(stripped, '>')
                if !isempty(current_name)
                    records[current_name] = FastaIndexRecord(current_name, sequence_length, sequence_start, line_bases, line_bytes)
                end

                current_name = strip(stripped[2:end])
                sequence_length = 0
                sequence_start = position(io)
                line_bases = 0
                line_bytes = 0
            else
                if sequence_length == 0
                    sequence_start = line_start
                    line_bases = length(stripped)
                    line_bytes = Int(position(io) - line_start)
                end
                sequence_length += length(stripped)
            end

            byte_offset = position(io)
        end

        if !isempty(current_name)
            records[current_name] = FastaIndexRecord(current_name, sequence_length, sequence_start, line_bases, line_bytes)
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, records, "fasta_index")
end

"""
    fetch_fasta_sequence(filebytes, index, start_pos, stop_pos)

Fetch a sequence range from an indexed FASTA byte buffer.
"""
function fetch_fasta_sequence(filebytes::AbstractVector{UInt8}, index::FastaIndexRecord, start_pos::Integer, stop_pos::Integer)
    start_pos < 1 && throw(ArgumentError("start_pos must be at least 1"))
    stop_pos < start_pos && return ""
    stop_pos > index.sequence_length && throw(ArgumentError("stop_pos exceeds sequence length"))

    start_line = div(start_pos - 1, index.line_bases)
    stop_line = div(stop_pos - 1, index.line_bases)
    if start_line == stop_line
        src_index = index.byte_offset + start_line * index.line_bytes + mod(start_pos - 1, index.line_bases) + 1
        return unsafe_string(pointer(filebytes, src_index), stop_pos - start_pos + 1)
    end

    bytes = Vector{UInt8}(undef, stop_pos - start_pos + 1)
    out_index = 1
    current_pos = start_pos

    while current_pos <= stop_pos
        line_index = div(current_pos - 1, index.line_bases)
        line_offset = mod(current_pos - 1, index.line_bases)
        line_remaining = index.line_bases - line_offset
        read_length = min(line_remaining, stop_pos - current_pos + 1)

        src_index = index.byte_offset + line_index * index.line_bytes + line_offset + 1
        copyto!(bytes, out_index, filebytes, src_index, read_length)

        out_index += read_length
        current_pos += read_length
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, String(bytes), "fetch_fasta_sequence")
end

"""
    fetch_fasta_sequence(path, index, start_pos, stop_pos)

Fetch a sequence range from an indexed FASTA file path.
"""
function fetch_fasta_sequence(path::String, index::FastaIndexRecord, start_pos::Integer, stop_pos::Integer)

    open(path, "r") do io

        return fetch_fasta_sequence(Mmap.mmap(io), index, start_pos, stop_pos)
    end
end

"""
    hamming_distance(left, right; use_cuda=false)

Compute the Hamming distance between two nucleotide strings.
"""
function hamming_distance(left::BioSequence{A}, right::BioSequence{A}; use_cuda::Bool=false) where {A <: BioAlphabet}
    if use_cuda
        left_bytes = Vector{UInt8}(left.data)
        right_bytes = Vector{UInt8}(right.data)
        return hamming_distance(left_bytes, right_bytes; use_cuda=true)
    end
    return hamming_distance(left.data, right.data)
end

"""
    kmer_frequency(bytes, k; use_cuda=false)

Count k-mer frequencies in a nucleotide byte vector.
"""
function kmer_frequency(bytes::AbstractVector{UInt8}, k::Integer; use_cuda::Bool=false)
    if use_cuda
        _ensure_cuda_sequence!()
        bytes_cuda = _is_cuda_backed_array(bytes) ? bytes : CUDA.CuArray{UInt8}(Vector{UInt8}(bytes))
        return Base.invokelatest(_CUDA_SEQUENCE_KMER_IMPL[], bytes_cuda, k)
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, kmer_frequency(String(bytes), k), "kmer_frequency")
end

"""
    check_primer_dimers(fwd_seq, rev_seq; min_3prime_match=4)

Check primer pairs for potentially problematic 3' complementarity.
"""
function check_primer_dimers(fwd_seq::BioSequence{DNAAlphabet}, rev_seq::BioSequence{DNAAlphabet}; min_3prime_match::Int=4)
    fwd = String(fwd_seq)
    rev_rc = String(reverse_complement(rev_seq))
    best_match = 0
    best_window = ""

    max_window = min(length(fwd), length(rev_rc))
    for window in 3:max_window
        fwd_suffix = fwd[end-window+1:end]
        rev_suffix = rev_rc[end-window+1:end]
        if fwd_suffix == rev_suffix && window > best_match
            best_match = window
            best_window = fwd_suffix
        end
    end

    _ctx = active_provenance_context()
    result = (
        has_dimer = best_match >= min_3prime_match,
        max_3prime_match = best_match,
        score = max_window == 0 ? 0.0 : best_match / max_window,
        matched_sequence = best_window)
    return provenance_result!(_ctx, result, "check_primer_dimers")
end

check_primer_dimers(fwd_seq::AbstractString, rev_seq::AbstractString; min_3prime_match::Int=4) =
    check_primer_dimers(_sequence_to_dna(fwd_seq), _sequence_to_dna(rev_seq); min_3prime_match=min_3prime_match)

"""
    _find_all_occurrences(sequence, motif)

Find all exact occurrences of a motif in a sequence.
"""
function _find_all_occurrences(sequence::BioSequence{A}, motif::BioSequence{A}) where {A <: BioAlphabet}
    sequence_str = String(sequence)
    motif_str = String(motif)
    matches = Tuple{Int,Int}[]
    start_index = firstindex(sequence_str)
    motif_length = ncodeunits(motif_str)
    motif_length == 0 && return matches

    while start_index <= lastindex(sequence_str)
        found = findnext(motif_str, sequence_str, start_index)
        found === nothing && break
        push!(matches, (first(found), last(found)))
        start_index = nextind(sequence_str, last(found))
    end

    return matches
end

"""
    _resolve_genome_sequences(genome_index; genome_path=nothing)

Resolve genome sequence sources from either an index or a FASTA path.
"""
function _resolve_genome_sequences(genome_index; genome_path::Union{Nothing,String}=nothing)
    if genome_index isa AbstractDict
        sequence_values = collect(Base.values(genome_index))
        if isempty(sequence_values)
            return Dict{String,DNASeq}()
        elseif all(value -> value isa BioSequence{DNAAlphabet}, sequence_values)
            return Dict{String,DNASeq}(String(chrom) => DNASeq(sequence.data; validate=false) for (chrom, sequence) in genome_index)
        elseif all(value -> value isa AbstractString, sequence_values)
            return Dict{String,DNASeq}(String(chrom) => DNASeq(String(sequence); validate=false) for (chrom, sequence) in genome_index)
        elseif all(value -> value isa FastaIndexRecord, sequence_values)
            genome_path === nothing && throw(ArgumentError("genome_path is required when genome_index stores FastaIndexRecord values"))
            sequences = Dict{String,DNASeq}()
            for (chrom, index) in genome_index
                sequences[String(chrom)] = DNASeq(fetch_fasta_sequence(genome_path, index, 1, index.sequence_length))
            end
            return sequences
        end
    elseif genome_index isa String
        path = String(genome_index)
        return Dict{String,DNASeq}(name => DNASeq(fetch_fasta_sequence(path, index, 1, index.sequence_length)) for (name, index) in fasta_index(path))
    end

    throw(ArgumentError("genome_index must be a Dict of sequences or FastaIndexRecord objects, or a FASTA path"))
end

"""
    pcr_in_silico(primer_fwd, primer_rev, genome_index; kwargs...)

Simulate PCR in silico by locating primer hits and possible amplicons.
"""
function pcr_in_silico(primer_fwd::BioSequence{DNAAlphabet}, primer_rev::BioSequence{DNAAlphabet}, genome_index; genome_path::Union{Nothing,String}=nothing, max_product_length::Int=50_000)
    genome = _resolve_genome_sequences(genome_index; genome_path=genome_path)
    forward = DNASeq(String(primer_fwd); validate=false)
    reverse = reverse_complement(DNASeq(String(primer_rev); validate=false))
    amplicons = NamedTuple[]

    for (chrom, sequence) in genome
        forward_sites = _find_all_occurrences(sequence, forward)
        reverse_sites = _find_all_occurrences(sequence, reverse)
        isempty(forward_sites) && continue
        isempty(reverse_sites) && continue

        for (forward_start, forward_stop) in forward_sites
            for (reverse_start, reverse_stop) in reverse_sites
                reverse_stop <= forward_start && continue
                product_size = reverse_stop - forward_start + 1
                product_size <= 0 && continue
                product_size > max_product_length && continue
                push!(amplicons, (
                    chrom = chrom,
                    start = forward_start,
                    stop = reverse_stop,
                    amplicon_size = product_size,
                    forward_primer_start = forward_start,
                    reverse_primer_stop = reverse_stop))
            end
        end
    end

    sort!(amplicons; by = amplicon -> (amplicon.chrom, amplicon.start, amplicon.stop, amplicon.amplicon_size))
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, amplicons, "pcr_in_silico")
end

function pcr_in_silico(primer_fwd::String, primer_rev::String, genome_index; genome_path::Union{Nothing,String}=nothing, max_product_length::Int=50_000)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, pcr_in_silico(DNASeq(primer_fwd; validate=false), DNASeq(primer_rev; validate=false), genome_index; genome_path=genome_path, max_product_length=max_product_length), "pcr_in_silico")
end

# Removed pcr_in_silico(::AbstractString, ::AbstractString) - use DNASeq types instead

"""
    score_grna_efficiency(guide; pam="NGG")

Score a guide RNA sequence for simple CRISPR efficiency heuristics.
"""
function score_grna_efficiency(guide::BioSequence{DNAAlphabet}; pam::BioSequence{DNAAlphabet}=DNASeq("NGG"; validate=false))
    sequence = uppercase(strip(String(guide)))
    isempty(sequence) && return 0.0

    gc = gc_content(DNASeq(sequence; validate=false))
    gc_score = clamp(1.0 - abs(gc - 0.5) / 0.5, 0.0, 1.0)
    poly_penalty = occursin(r"A{4,}|C{4,}|G{4,}|T{4,}", sequence) ? 0.25 : 0.0
    seed_bonus = startswith(sequence, "G") ? 0.05 : 0.0
    pam_bonus = endswith(uppercase(String(pam)), "GG") ? (endswith(sequence, "G") ? 0.05 : 0.0) : 0.0

    score = 0.6 * gc_score + seed_bonus + pam_bonus - poly_penalty
    return clamp(score, 0.0, 1.0)
end

score_grna_efficiency(guide::AbstractString; pam::AbstractString="NGG") =
    score_grna_efficiency(_sequence_to_dna(guide); pam=_sequence_to_dna(pam))

# ─── GC Skew ──────────────────────────────────────────────────────────────────

"""
    gc_skew(sequence)

Compute cumulative GC skew for a sequence.
"""
function gc_skew(sequence::BioSequence{A}) where {A <: BioAlphabet}
    bytes = sequence.data
    len = length(bytes)
    skew_values = Vector{Int}(undef, len)
    current = 0

    @inbounds for i in 1:len
        ch = bytes[i] | 0x20  # lowercase
        current += (ch == UInt8('g')) - (ch == UInt8('c'))
        skew_values[i] = current
    end

    return skew_values
end

"""
    gc_skew(bytes)

Compute cumulative GC skew for a nucleotide byte vector.
"""
function gc_skew(bytes::AbstractVector{UInt8})
    len = length(bytes)
    skew_values = Vector{Int}(undef, len)
    current = 0

    @inbounds for i in 1:len
        ch = bytes[i] | 0x20
        current += (ch == UInt8('g')) - (ch == UInt8('c'))
        skew_values[i] = current
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, skew_values, "gc_skew")
end

gc_skew(sequence::AbstractString) = gc_skew(codeunits(String(sequence)))

"""
    minimum_skew(sequence)

Return positions where the cumulative GC skew reaches its minimum.
"""
function minimum_skew(sequence)
    skew_values = gc_skew(sequence)
    isempty(skew_values) && return Int[]
    min_val = minimum(skew_values)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, findall(==(min_val), skew_values), "minimum_skew")
end

# ─── Dotmatrix ────────────────────────────────────────────────────────────────

"""
    dotmatrix(seq1, seq2)

Compute a dot-plot matrix comparing two sequences.
Returns a Matrix{Int8} where 1 = match, 0 = mismatch.
Rows correspond to seq1, columns to seq2.
"""
function dotmatrix(seq1::BioSequence{A}, seq2::BioSequence{A}; use_cuda::Bool=false) where {A <: BioAlphabet}
    if use_cuda
        seq1_bytes = Vector{UInt8}(seq1.data)
        seq2_bytes = Vector{UInt8}(seq2.data)
        return dotmatrix(seq1_bytes, seq2_bytes; use_cuda=true)
    end
    return dotmatrix(seq1.data, seq2.data)
end

"""
    dotmatrix(seq1, seq2; use_cuda=false)

Compute a dot plot matrix for two nucleotide byte vectors.
"""
function dotmatrix(seq1::AbstractVector{UInt8}, seq2::AbstractVector{UInt8}; use_cuda::Bool=false)
    if use_cuda
        _ensure_cuda_sequence!()
        seq1_cuda = _is_cuda_backed_array(seq1) ? seq1 : CUDA.CuArray{UInt8}(Vector{UInt8}(seq1))
        seq2_cuda = _is_cuda_backed_array(seq2) ? seq2 : CUDA.CuArray{UInt8}(Vector{UInt8}(seq2))
        return Base.invokelatest(_CUDA_SEQUENCE_DOTMATRIX_IMPL[], seq1_cuda, seq2_cuda)
    end

    m = length(seq1)
    n = length(seq2)
    matrix = Matrix{Int8}(undef, m, n)

    @inbounds for j in 1:n
        b2 = seq2[j] | 0x20  # case-insensitive
        for i in 1:m
            matrix[i, j] = Int8((seq1[i] | 0x20) == b2)
        end
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, matrix, "dotmatrix")
end

# Removed dotmatrix(::AbstractString, ::AbstractString) - use BioSequence types instead
function dotmatrix(seq1::AbstractString, seq2::AbstractString; use_cuda::Bool=false)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, dotmatrix(DNASeq(seq1; validate=false), DNASeq(seq2; validate=false); use_cuda=use_cuda), "dotmatrix")
end

# ==============================================================================
# Missing Bioconductor-equivalent features
#
# The following functions mirror capabilities from R/Bioconductor's Biostrings,
# IRanges, and GenomicRanges packages that were absent from BioToolkit.
# ==============================================================================

# ---- dinucleotide_frequency (Biostrings::dinucleotideFrequency) ---------------

"""
    dinucleotide_frequency(sequence)

Count all 16 dinucleotide frequencies in a DNA/RNA sequence.
Returns a Dict{String,Int} with keys like "AA", "AC", ..., "TT".

Equivalent to R/Bioconductor's `Biostrings::dinucleotideFrequency`.
"""
function dinucleotide_frequency(sequence::BioSequence{A}) where {A <: BioAlphabet}
    counts = Dict{String,Int}()
    bytes = sequence.data
    len = length(bytes)
    len < 2 && return counts

    @inbounds for i in 1:(len - 1)
        b1 = bytes[i] < 0x61 ? bytes[i] : bytes[i] - 0x20
        b2 = bytes[i+1] < 0x61 ? bytes[i+1] : bytes[i+1] - 0x20
        key = String(UInt8[b1, b2])
        counts[key] = get(counts, key, 0) + 1
    end

    return counts
end

"""
    trinucleotide_frequency(sequence)

Count all 64 trinucleotide frequencies in a DNA/RNA sequence.
Returns a Dict{String,Int} with keys like "AAA", "AAC", etc.

Equivalent to R/Bioconductor's `Biostrings::trinucleotideFrequency`.
"""
function trinucleotide_frequency(sequence::BioSequence{A}) where {A <: BioAlphabet}
    counts = Dict{String,Int}()
    bytes = sequence.data
    len = length(bytes)
    len < 3 && return counts

    @inbounds for i in 1:(len - 2)
        b1 = bytes[i] < 0x61 ? bytes[i] : bytes[i] - 0x20
        b2 = bytes[i+1] < 0x61 ? bytes[i+1] : bytes[i+1] - 0x20
        b3 = bytes[i+2] < 0x61 ? bytes[i+2] : bytes[i+2] - 0x20
        key = String(UInt8[b1, b2, b3])
        counts[key] = get(counts, key, 0) + 1
    end

    return counts
end

# ---- letter_frequency (Biostrings::letterFrequency) --------------------------

"""
    letter_frequency(sequence)

Count every distinct character in a sequence. Returns a Dict{Char,Int}.
Case-insensitive: all characters are uppercased.

Equivalent to R/Bioconductor's `Biostrings::letterFrequency`.
"""
function letter_frequency(sequence::BioSequence{A}) where {A <: BioAlphabet}
    counts = Dict{Char,Int}()
    @inbounds for byte in sequence.data
        ch = Char(byte < 0x61 ? byte : (byte <= 0x7a ? byte - 0x20 : byte))
        counts[ch] = get(counts, ch, 0) + 1
    end
    return counts
end

# ---- sequence_complexity (linguistic complexity) -----------------------------

"""
    sequence_complexity(sequence; k::Int=3)

Compute the linguistic complexity of a sequence, defined as the ratio of
observed unique k-mers to the theoretical maximum for the given alphabet
and sequence length.  Values near 1.0 indicate high-complexity sequences;
low values indicate repetitive or low-complexity regions.

This is useful for masking low-complexity regions before BLAST searches,
similar to DUST/SEG algorithms.

Reference: Trifonov (1990) Bull Math Biol 52(1):35-40
"""
function sequence_complexity(sequence::BioSequence{A}; k::Int=3) where {A <: BioAlphabet}
    len = length(sequence)
    len < k && return 0.0

    observed = Set{UInt64}()
    bytes = sequence.data

    @inbounds for i in 1:(len - k + 1)
        h = UInt64(0)
        for j in 0:(k-1)
            b = bytes[i + j]
            b = b < 0x61 ? b : (b <= 0x7a ? b - 0x20 : b)
            h = (h << 8) | b
        end
        push!(observed, h)
    end

    # Theoretical maximum: min(4^k, L-k+1) for DNA
    max_possible = min(4^k, len - k + 1)
    return length(observed) / max_possible
end

# ---- generalized sequence utilities ----------------------------------------

"""
    oligonucleotide_frequency(sequence; k::Int=2, normalize::Bool=false)

Count all observed k-mer frequencies in a biological sequence.

When `normalize=true`, returns relative frequencies that sum to 1.0.
"""
function oligonucleotide_frequency(sequence::BioSequence{A}; k::Int=2, normalize::Bool=false) where {A <: BioAlphabet}
    counts = kmer_frequency(sequence, k)
    normalize || return counts

    total = sum(values(counts))
    total == 0 && return Dict{String,Float64}()
    return Dict(kmer => count / total for (kmer, count) in counts)
end

oligonucleotide_frequency(sequence::AbstractString; k::Int=2, normalize::Bool=false) = oligonucleotide_frequency(DNASeq(sequence; validate=false); k=k, normalize=normalize)

"""
    sequence_entropy(sequence; k::Int=1)

Compute the Shannon entropy of the observed k-mer distribution in a sequence.
For `k=1`, this is the single-letter entropy of the sequence.
"""
function sequence_entropy(sequence::BioSequence{A}; k::Int=1) where {A <: BioAlphabet}
    k <= 0 && throw(ArgumentError("k must be positive"))

    counts = k == 1 ? letter_frequency(sequence) : kmer_frequency(sequence, k)
    total = sum(values(counts))
    total == 0 && return 0.0

    entropy = 0.0
    for count in values(counts)
        probability = count / total
        entropy -= probability * log2(probability)
    end
    return entropy
end

sequence_entropy(sequence::AbstractString; k::Int=1) = sequence_entropy(DNASeq(sequence; validate=false); k=k)

"""
    sliding_window_gc_content(sequence; window::Int=100, step::Int=1)

Compute GC content in sliding windows across a DNA sequence.
Returns a named tuple with start positions, stop positions, and GC fractions.
"""
function sliding_window_gc_content(sequence::BioSequence{DNAAlphabet}; window::Int=100, step::Int=1)
    window <= 0 && throw(ArgumentError("window must be positive"))
    step <= 0 && throw(ArgumentError("step must be positive"))

    len = length(sequence)
    if len < window
        return (start_positions=Int[], stop_positions=Int[], gc_content=Float64[])
    end

    starts = Int[]
    stops = Int[]
    gc_content = Float64[]
    bytes = sequence.data

    @inbounds for start_pos in 1:step:(len - window + 1)
        stop_pos = start_pos + window - 1
        gc_count = 0
        for index in start_pos:stop_pos
            byte = bytes[index]
            byte = byte < 0x61 ? byte : (byte <= 0x7a ? byte - 0x20 : byte)
            gc_count += byte == UInt8('G') || byte == UInt8('C')
        end

        push!(starts, start_pos)
        push!(stops, stop_pos)
        push!(gc_content, gc_count / window)
    end

    _ctx = active_provenance_context()
    return provenance_result!(_ctx, (start_positions=starts, stop_positions=stops, gc_content=gc_content), "sliding_window_gc_content")
end

sliding_window_gc_content(sequence::AbstractString; window::Int=100, step::Int=1) = sliding_window_gc_content(DNASeq(sequence; validate=false); window=window, step=step)

