using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_container_provenance!, register_provenance!

@inline function _register_schema_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

@inline function _metadata_without_provenance(metadata::AbstractDict)
    return Dict(key => value for (key, value) in metadata if key != :provenance && key != "provenance" && key != PROVENANCE_ID_KEY && key != String(PROVENANCE_ID_KEY) && key != PROVENANCE_HASH_KEY && key != String(PROVENANCE_HASH_KEY))
end

"""
    VariantEvent

Compact Arrow-friendly representation of a variant record.
"""
struct VariantEvent
    chrom::UInt16
    pos::Int32
    id::UInt32
    ref::UInt8
    alt::UInt8
    qual::Float32
    hasqual::Bool
end

"""
    VariantTextRecord

Text-oriented variant record used as the parsing front end before Arrow encoding.
"""
struct VariantTextRecord
    chrom::String
    pos::Int32
    id::String
    ref::String
    alt::String
    qual::Union{Missing,Float32}
    filter::String
    info::String
    format::String
    samples::Vector{String}
    metadata::Dict{Symbol,Any}
end

function VariantTextRecord(
    chrom::AbstractString,
    pos::Integer,
    id::AbstractString,
    ref::AbstractString,
    alt::AbstractString,
    qual::Union{Missing,Real};
    filter::AbstractString="PASS",
    info::AbstractString=".",
    format::AbstractString="",
    samples::AbstractVector{<:AbstractString}=String[],
    metadata::AbstractDict=Dict{Symbol,Any}(),
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))
    qual_value = qual === missing ? missing : Float32(qual)
    metadata_copy = Dict{Symbol,Any}(metadata)
    ensure_provenance_id!(metadata_copy)
    result = VariantTextRecord(
        String(chrom),
        Int32(pos),
        String(id),
        String(ref),
        String(alt),
        qual_value,
        String(filter),
        String(info),
        String(format),
        String.(samples),
        metadata_copy)
    return _register_schema_result!(_ctx, result, "VariantTextRecord"; parents=String[], parameters=(chrom=String(chrom), pos=Int32(pos), id=String(id), ref=String(ref), alt=String(alt)))
end

"""
    VcfHeader

Structured VCF header metadata, including raw meta-lines and sample names.
"""
struct VcfHeader
    meta_lines::Vector{String}
    columns::Vector{String}
    sample_names::Vector{String}
end

VcfHeader() = VcfHeader(String[], String["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"], String[])

function VcfHeader(
    meta_lines::AbstractVector{<:AbstractString},
    columns::AbstractVector{<:AbstractString},
    sample_names::AbstractVector{<:AbstractString}=String[])
    return VcfHeader(String.(meta_lines), String.(columns), String.(sample_names))
end

"""
    VcfDocument

Full VCF document containing header metadata and parsed records.
"""
struct VcfDocument
    header::VcfHeader
    records::Vector{VariantTextRecord}
    metadata::Dict{Symbol,Any}
end

function VcfDocument(records::AbstractVector{<:VariantTextRecord}; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    metadata = Dict{Symbol,Any}()
    ensure_provenance_id!(metadata)
    result = VcfDocument(VcfHeader(), Vector{VariantTextRecord}(records), metadata)
    return _register_schema_result!(_ctx, result, "VcfDocument"; parents=provenance_parent_ids(records), parameters=(record_count=length(records), sample_count=length(result.header.sample_names)))
end

function VcfDocument(header::VcfHeader, records::AbstractVector{<:VariantTextRecord}; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    metadata = Dict{Symbol,Any}()
    ensure_provenance_id!(metadata)
    result = VcfDocument(header, Vector{VariantTextRecord}(records), metadata)
    return _register_schema_result!(_ctx, result, "VcfDocument"; parents=provenance_parent_ids(records), parameters=(record_count=length(records), sample_count=length(header.sample_names)))
end

function Base.:(==)(left::VariantTextRecord, right::VariantTextRecord)
    return isequal(left.chrom, right.chrom) &&
        isequal(left.pos, right.pos) &&
        isequal(left.id, right.id) &&
        isequal(left.ref, right.ref) &&
        isequal(left.alt, right.alt) &&
        isequal(left.qual, right.qual) &&
        isequal(left.filter, right.filter) &&
        isequal(left.info, right.info) &&
        isequal(left.format, right.format) &&
        isequal(left.samples, right.samples) &&
        isequal(_metadata_without_provenance(left.metadata), _metadata_without_provenance(right.metadata))
end

function Base.:(==)(left::VcfHeader, right::VcfHeader)
    return isequal(left.meta_lines, right.meta_lines) && isequal(left.columns, right.columns) && isequal(left.sample_names, right.sample_names)
end

Base.:(==)(left::VcfDocument, right::VcfDocument) = isequal(left.header, right.header) && isequal(left.records, right.records) && isequal(_metadata_without_provenance(left.metadata), _metadata_without_provenance(right.metadata))

Base.length(doc::VcfDocument) = length(doc.records)
Base.iterate(doc::VcfDocument, state...) = iterate(doc.records, state...)

function Base.show(io::IO, record::VariantTextRecord)
    print(io, "VariantTextRecord(", record.chrom, ":", record.pos, ", ", record.ref, "->", record.alt, ", ", container_provenance_summary(record), ")")
end

function Base.show(io::IO, document::VcfDocument)
    print(io, "VcfDocument(", length(document.records), " records, ", container_provenance_summary(document), ")")
end

"""
    arrow_schema(::Type{VariantEvent})

Return the Arrow schema for compact variant events.
"""
function arrow_schema(::Type{VariantEvent}; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    result = Tables.Schema(fieldnames(VariantEvent), fieldtypes(VariantEvent))
    _ctx !== nothing && register_provenance!(_ctx, "arrow_schema"; parameters=(field_count=length(fieldnames(VariantEvent)),))
    return result
end

"""
    encode_chromosome(chrom)

Encode a chromosome label into a compact numeric identifier.
"""
function encode_chromosome(chrom::String)
    normalized = startswith(chrom, "chr") ? chrom[4:end] : chrom

    if normalized == "X"
        return UInt16(23)
    elseif normalized == "Y"
        return UInt16(24)
    elseif normalized in ("M", "MT")
        return UInt16(25)
    end

    parsed = tryparse(Int, normalized)
    parsed === nothing && return UInt16(0)

    return UInt16(clamp(parsed, 0, typemax(UInt16)))
end

"""
    encode_base(base)

Encode a nucleotide base into a compact numeric value.
"""
function encode_base(base::String)
    upper = uppercase(base)
    upper == "A" && return UInt8(1)
    upper == "C" && return UInt8(2)
    upper == "G" && return UInt8(3)
    upper == "T" && return UInt8(4)
    upper == "N" && return UInt8(5)
    throw(ArgumentError("unsupported base sequence for compact VCF encoding: $(base)"))
end

"""
    encode_identifier(identifier)

Hash a variant identifier into a compact 32-bit value.
"""
function encode_identifier(identifier::String)
    hash_value = UInt32(2166136261)

    for byte in codeunits(identifier)
        hash_value ⊻= UInt32(byte)
        hash_value *= UInt32(16777619)
    end

    return hash_value
end

"""
    compact_variant_event(record)

Convert a text variant record into its compact Arrow-friendly form.
"""
function compact_variant_event(record::VariantTextRecord; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    qual_value = record.qual === missing ? Float32(NaN) : Float32(record.qual)
    result = VariantEvent(
        encode_chromosome(record.chrom),
        record.pos,
        encode_identifier(record.id),
        encode_base(record.ref),
        encode_base(record.alt),
        qual_value,
        record.qual !== missing)


    return _register_schema_result!(_ctx, result, "compact_variant_event"; parents=provenance_parent_ids(record), parameters=(chrom=record.chrom, pos=record.pos, id=record.id))
end
