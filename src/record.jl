# ==============================================================================
# record.jl — Biological Sequence Record Types
#
# Provides parametric sequence records that bundle biological sequences with
# metadata, identifiers, and quality scores. These types are designed to
# provide parity with Bioconductor's DNAStringSet and SeqRecord structures.
#
# Records are parametric on the BioAlphabet {A} to ensure strict type safety:
# a FastqRecord{DNAAlphabet} cannot be accidentally treated as protein.
# ==============================================================================

using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_container_provenance!

"""
    SeqRecord{A <: BioAlphabet}

An annotated sequence record bundling a `BioSequence{A}` with metadata.
"""
struct SeqRecord{A <: BioAlphabet}
    sequence::BioSequence{A}
    identifier::String
    description::String
    metadata::Dict{Symbol, Any}
end

"""
    FastqRecord{A <: BioAlphabet}

A FASTQ record bundling a `BioSequence{A}` with quality scores.
"""
struct FastqRecord{A <: BioAlphabet}
    sequence::BioSequence{A}
    identifier::String
    description::String
    quality::String # Stored as raw Phred+33 ASCII string
    metadata::Dict{Symbol, Any}
end

"""
    SeqRecordLite{A <: BioAlphabet}

Lightweight record compatible with MSA/annotation/quality modules while using
typed `BioSequence{A}` storage.
"""
struct SeqRecordLite{A <: BioAlphabet}
    sequence::BioSequence{A}
    identifier::String
    name::String
    description::String
    annotations::Dict{Symbol, Any}
    letter_annotations::Dict{Symbol, Any}
end

@inline function _register_record_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

# ---- Constructors ------------------------------------------------------------

"""
    SeqRecord(sequence::BioSequence{A}; identifier="", description="", metadata=Dict())
"""
function SeqRecord(
    sequence::BioSequence{A};
    identifier::String="",
    description::String=identifier,
    metadata::AbstractDict=Dict{Symbol, Any}(),
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx)) where {A <: BioAlphabet}
    normalized_metadata = Dict{Symbol, Any}(metadata)
    ensure_provenance_id!(normalized_metadata)
    result = SeqRecord{A}(sequence, String(identifier), String(description), normalized_metadata)
    return _register_record_result!(_ctx, result, "SeqRecord"; parents=provenance_parent_ids(sequence), parameters=(identifier=String(identifier), description=String(description), length=length(sequence)))
end

"""
    FastqRecord(sequence::BioSequence{A}, quality; identifier="", description="")
"""
function FastqRecord(
    sequence::BioSequence{A},
    quality::String;
    identifier::String="",
    description::String=identifier,
    metadata::AbstractDict=Dict{Symbol, Any}(),
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx)) where {A <: BioAlphabet}
    length(sequence) == length(quality) || throw(ArgumentError("sequence and quality lengths must match"))
    normalized_metadata = Dict{Symbol, Any}(metadata)
    ensure_provenance_id!(normalized_metadata)
    result = FastqRecord{A}(sequence, String(identifier), String(description), String(quality), normalized_metadata)
    return _register_record_result!(_ctx, result, "FastqRecord"; parents=provenance_parent_ids(sequence), parameters=(identifier=String(identifier), description=String(description), length=length(sequence)))
end

# Backward-compatible positional constructor used in quality.jl.
function FastqRecord(
    identifier::String,
    description::String,
    sequence::BioSequence{A},
    quality::String) where {A <: BioAlphabet}
    return FastqRecord(sequence, quality; identifier=identifier, description=description)
end

function _infer_sequence_alphabet(sequence)
    sequence_text = String(sequence)
    if validate_sequence(DNAAlphabet, sequence_text)
        return DNAAlphabet
    elseif validate_sequence(RNAAlphabet, sequence_text)
        return RNAAlphabet
    elseif validate_sequence(AminoAcidAlphabet, sequence_text)
        return AminoAcidAlphabet
    end
    throw(ArgumentError("cannot infer alphabet for sequence"))
end

function _validate_letter_annotations!(sequence_length::Int, letter_annotations::Dict{Symbol, Any})
    for (key, value) in letter_annotations
        if value isa String
            length(value) == sequence_length || throw(ArgumentError("letter annotation '$key' length $(length(value)) must match sequence length $sequence_length"))
        elseif value isa AbstractVector
            length(value) == sequence_length || throw(ArgumentError("letter annotation '$key' length $(length(value)) must match sequence length $sequence_length"))
        end
    end
    return letter_annotations
end

function SeqRecordLite(
    sequence::BioSequence{A};
    identifier::String="",
    name::String=identifier,
    description::String=name,
    annotations::AbstractDict=Dict{Symbol, Any}(),
    letter_annotations::AbstractDict=Dict{Symbol, Any}(),
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx)) where {A <: BioAlphabet}
    normalized_letter_annotations = Dict{Symbol, Any}(letter_annotations)
    _validate_letter_annotations!(length(sequence), normalized_letter_annotations)
    annotations_copy = Dict{Symbol, Any}(annotations)
    ensure_provenance_id!(annotations_copy)
    result = SeqRecordLite{A}(
        sequence,
        String(identifier),
        String(name),
        String(description),
        annotations_copy,
        normalized_letter_annotations)
    return _register_record_result!(_ctx, result, "SeqRecordLite"; parents=provenance_parent_ids(sequence), parameters=(identifier=String(identifier), name=String(name), description=String(description), length=length(sequence)))
end

function SeqRecordLite(
    sequence::AbstractString;
    identifier::String="",
    name::String=identifier,
    description::String=name,
    annotations::AbstractDict=Dict{Symbol, Any}(),
    letter_annotations::AbstractDict=Dict{Symbol, Any}(),
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))
    sequence_text = String(sequence)
    alphabet = _infer_sequence_alphabet(sequence_text)
    typed_sequence = BioSequence{alphabet}(sequence_text)
    return SeqRecordLite(
        typed_sequence;
        identifier=identifier,
        name=name,
        description=description,
        annotations=annotations,
        letter_annotations=letter_annotations,
        _ctx=_ctx)
end

@inline _fastq_components(record::FastqRecord) = (record.identifier, record.description, String(record.sequence), record.quality)

@inline function _fastq_components(record::SeqRecordLite)
    quality = get(record.letter_annotations, :quality, nothing)
    quality === nothing && throw(ArgumentError("FASTQ records require a :quality letter annotation"))
    return (record.identifier, record.description, String(record.sequence), quality isa String ? String(quality) : String(UInt8[UInt8(character) for character in quality]))
end

# ---- Utility Functions -------------------------------------------------------

Base.length(record::Union{SeqRecord, FastqRecord, SeqRecordLite}) = length(record.sequence)
Base.String(record::SeqRecordLite) = String(record.sequence)
Base.convert(::Type{String}, record::SeqRecordLite) = String(record)

Base.:(==)(left::SeqRecord{A}, right::SeqRecord{A}) where {A} = left.sequence == right.sequence && left.identifier == right.identifier && left.description == right.description && _metadata_without_provenance(left.metadata) == _metadata_without_provenance(right.metadata)
Base.:(==)(left::FastqRecord{A}, right::FastqRecord{A}) where {A} = left.sequence == right.sequence && left.identifier == right.identifier && left.description == right.description && left.quality == right.quality && _metadata_without_provenance(left.metadata) == _metadata_without_provenance(right.metadata)
Base.:(==)(left::SeqRecordLite{A}, right::SeqRecordLite{A}) where {A} = left.sequence == right.sequence && left.identifier == right.identifier && left.name == right.name && left.description == right.description && _metadata_without_provenance(left.annotations) == _metadata_without_provenance(right.annotations) && left.letter_annotations == right.letter_annotations

function Base.show(io::IO, record::SeqRecord{A}) where {A}
    print(io, "SeqRecord{", A, "}(", record.identifier, ", ", length(record), " bp, ", container_provenance_summary(record), ")")
end

function Base.show(io::IO, record::FastqRecord{A}) where {A}
    print(io, "FastqRecord{", A, "}(", record.identifier, ", ", length(record), " bp, ", container_provenance_summary(record), ")")
end

function Base.show(io::IO, record::SeqRecordLite{A}) where {A}
    print(io, "SeqRecordLite{", A, "}(", record.identifier, ", ", length(record), " bp, ", container_provenance_summary(record), ")")
end

"""
    quality_scores(record::FastqRecord; offset=33) -> Vector{Int}

Convert the ASCII quality string into a vector of integer Phred scores.
"""
function quality_scores(record::FastqRecord; offset::Int=33)
    return [Int(c) - offset for c in record.quality]
end
