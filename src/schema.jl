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
end

"""
    arrow_schema(::Type{VariantEvent})

Return the Arrow schema for compact variant events.
"""
function arrow_schema(::Type{VariantEvent})
    return Tables.Schema(fieldnames(VariantEvent), fieldtypes(VariantEvent))
end

"""
    encode_chromosome(chrom)

Encode a chromosome label into a compact numeric identifier.
"""
function encode_chromosome(chrom::AbstractString)
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
function encode_base(base::AbstractString)
    upper = uppercase(base)
    upper == "A" && return UInt8(1)
    upper == "C" && return UInt8(2)
    upper == "G" && return UInt8(3)
    upper == "T" && return UInt8(4)
    upper == "N" && return UInt8(5)
    return UInt8(0)
end

"""
    encode_identifier(identifier)

Hash a variant identifier into a compact 32-bit value.
"""
function encode_identifier(identifier::AbstractString)
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
function compact_variant_event(record::VariantTextRecord)
    qual_value = record.qual === missing ? Float32(NaN) : Float32(record.qual)
    return VariantEvent(
        encode_chromosome(record.chrom),
        record.pos,
        encode_identifier(record.id),
        encode_base(record.ref),
        encode_base(record.alt),
        qual_value,
        record.qual !== missing,
    )
end
