"""
    FastqRecord

Plain FASTQ record with identifier, description, sequence, and quality string.
"""
struct FastqRecord
    identifier::String
    description::String
    sequence::String
    quality::String
end

"""
    SeqRecordLite

Mutable lightweight sequence record with annotations and per-letter metadata.
"""
mutable struct SeqRecordLite
    sequence::String
    identifier::String
    name::String
    description::String
    annotations::Dict{Symbol,Any}
    letter_annotations::Dict{Symbol,Any}
end

"""
    SeqRecordLite(sequence; kwargs...)

Construct a lightweight sequence record from a raw sequence string.
"""
function SeqRecordLite(
    sequence::AbstractString;
    identifier::AbstractString="",
    name::AbstractString=identifier,
    description::AbstractString=name,
    annotations::AbstractDict=Dict{Symbol,Any}(),
    letter_annotations::AbstractDict=Dict{Symbol,Any}(),
)
    return SeqRecordLite(
        String(sequence),
        String(identifier),
        String(name),
        String(description),
        Dict{Symbol,Any}(annotations),
        Dict{Symbol,Any}(letter_annotations),
    )
end

"""
    SeqRecordLite(record::FastqRecord; annotations=...)

Convert a FASTQ record into a lightweight annotated record.
"""
SeqRecordLite(record::FastqRecord; annotations::AbstractDict=Dict{Symbol,Any}()) = SeqRecordLite(
    record.sequence;
    identifier=record.identifier,
    name=record.identifier,
    description=record.description,
    annotations=annotations,
    letter_annotations=Dict{Symbol,Any}(:quality => record.quality),
)

"""
    Base.length(record)

Return the sequence length of a lightweight record in code units.
"""
Base.length(record::SeqRecordLite) = ncodeunits(record.sequence)

"""
    Base.show(io, record)

Render a compact summary for a FASTQ record.
"""
function Base.show(io::IO, record::FastqRecord)
    print(io, "FastqRecord(", record.identifier, ", ", record.description, ", ", length(record.sequence), " bp)")
end

"""
    Base.show(io, record)

Render a compact summary for a lightweight sequence record.
"""
function Base.show(io::IO, record::SeqRecordLite)
    print(io, "SeqRecordLite(", record.identifier, ", ", length(record.sequence), " bp)")
end

"""
    _fastq_components(record)

Return the four FASTQ fields stored in a FASTQ record.
"""
function _fastq_components(record::FastqRecord)
    return record.identifier, record.description, record.sequence, record.quality
end

"""
    _fastq_components(record; quality_key=:quality)

Return the four FASTQ fields stored in a lightweight record.
"""
function _fastq_components(record::SeqRecordLite; quality_key::Symbol=:quality)
    quality = get(record.letter_annotations, quality_key, nothing)
    quality === nothing && throw(ArgumentError("missing $(quality_key) letter annotation"))
    quality isa AbstractString || throw(ArgumentError("FASTQ quality annotation must be a string"))
    return record.identifier, record.description, record.sequence, String(quality)
end

"""
    write_fastq(path, records; quality_key=:quality)

Write FASTQ records or lightweight records to disk in FASTQ format.
"""
function write_fastq(path::AbstractString, records; quality_key::Symbol=:quality)
    open(path, "w") do io
        for record in records
            identifier, description, sequence, quality = record isa FastqRecord ? _fastq_components(record) : _fastq_components(record; quality_key=quality_key)
            length(sequence) == length(quality) || throw(ArgumentError("sequence and quality lengths differ"))

            header = isempty(description) ? identifier : description
            write(io, "@")
            write(io, header)
            write(io, "\n")
            write(io, sequence)
            write(io, "\n+\n")
            write(io, quality)
            write(io, "\n")
        end
    end

    return path
end
