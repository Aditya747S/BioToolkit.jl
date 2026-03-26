using DataFrames

export AbstractFeatureLocation, FeatureLocationLite, CompoundFeatureLocation, SeqFeatureLite, AnnotatedSeqRecord, feature_spans, feature_bounds, feature_start, feature_stop, parse_feature_location, SangerTrace

"""
    AbstractFeatureLocation

Abstract parent type for all feature-location objects used by the annotation
model.
"""
abstract type AbstractFeatureLocation end

"""
    FeatureLocationLite

Simple contiguous feature interval with strand and partial-end flags.
"""
struct FeatureLocationLite <: AbstractFeatureLocation
    start::Int
    stop::Int
    strand::Int8
    partial_start::Bool
    partial_stop::Bool
end

"""
    CompoundFeatureLocation

Compound feature location composed of multiple parts such as joins or orders.
"""
struct CompoundFeatureLocation <: AbstractFeatureLocation
    operator::String
    parts::Vector{AbstractFeatureLocation}
    strand::Int8
end

"""
    SeqFeatureLite

Compact feature record storing the feature type, location, qualifiers, and id.
"""
struct SeqFeatureLite
    feature_type::String
    location::AbstractFeatureLocation
    qualifiers::Dict{String,Vector{String}}
    id::String
end

"""
    AnnotatedSeqRecord

Mutable annotated sequence record with metadata, per-letter annotations, and
feature storage.
"""
mutable struct AnnotatedSeqRecord
    sequence::String
    identifier::String
    name::String
    description::String
    annotations::Dict{Symbol,Any}
    letter_annotations::Dict{Symbol,Any}
    features::Vector{SeqFeatureLite}
end

"""
    SangerTrace

Chromatogram-style Sanger sequencing trace with signal arrays and qualities.
"""
struct SangerTrace
    sequence::String
    qualities::Vector{UInt8}
    trace_a::Vector{UInt16}
    trace_c::Vector{UInt16}
    trace_g::Vector{UInt16}
    trace_t::Vector{UInt16}
    annotations::Dict{String,Any}
end

"""
    FeatureLocationLite(start, stop; strand=1, partial_start=false, partial_stop=false)

Construct a contiguous feature location from integer coordinates.
"""
FeatureLocationLite(start::Integer, stop::Integer; strand::Integer=1, partial_start::Bool=false, partial_stop::Bool=false) = FeatureLocationLite(Int(start), Int(stop), Int8(strand), partial_start, partial_stop)

"""
    CompoundFeatureLocation(operator, parts; strand=1)

Construct a compound feature location from a list of child locations.
"""
CompoundFeatureLocation(operator::AbstractString, parts::AbstractVector{<:AbstractFeatureLocation}; strand::Integer=1) = CompoundFeatureLocation(String(operator), AbstractFeatureLocation[part for part in parts], Int8(strand))

"""
    SeqFeatureLite(feature_type, location; qualifiers=..., id="")

Construct a lightweight feature record from a type string and location.
"""
function SeqFeatureLite(feature_type::AbstractString, location::AbstractFeatureLocation; qualifiers::AbstractDict=Dict{String,Vector{String}}(), id::AbstractString="")
    return SeqFeatureLite(String(feature_type), location, Dict{String,Vector{String}}(qualifiers), String(id))
end

"""
    AnnotatedSeqRecord(sequence; kwargs...)

Construct a mutable annotated record from raw sequence and associated metadata.
"""
function AnnotatedSeqRecord(
    sequence::AbstractString;
    identifier::AbstractString="",
    name::AbstractString=identifier,
    description::AbstractString=name,
    annotations::AbstractDict=Dict{Symbol,Any}(),
    letter_annotations::AbstractDict=Dict{Symbol,Any}(),
    features::AbstractVector{SeqFeatureLite}=SeqFeatureLite[],
)
    feature_vector = features isa Vector{SeqFeatureLite} ? features : SeqFeatureLite[feature for feature in features]
    return AnnotatedSeqRecord(
        String(sequence),
        String(identifier),
        String(name),
        String(description),
        Dict{Symbol,Any}(annotations),
        Dict{Symbol,Any}(letter_annotations),
        feature_vector,
    )
end

"""
    Base.length(record)

Return the sequence length of an annotated record in code units.
"""
Base.length(record::AnnotatedSeqRecord) = ncodeunits(record.sequence)

"""
    Base.show(io, location)

Render a compact summary for a simple feature location.
"""
function Base.show(io::IO, location::FeatureLocationLite)
    print(io, "FeatureLocationLite(", location.start, ":", location.stop, ", strand=", location.strand, ")")
end

"""
    Base.show(io, location)

Render a compact summary for a compound feature location.
"""
function Base.show(io::IO, location::CompoundFeatureLocation)
    print(io, "CompoundFeatureLocation(", location.operator, ", parts=", length(location.parts), ", strand=", location.strand, ")")
end

"""
    Base.show(io, feature)

Render a compact summary for a lightweight feature record.
"""
function Base.show(io::IO, feature::SeqFeatureLite)
    print(io, "SeqFeatureLite(", feature.feature_type, ", id=", feature.id, ", location=")
    show(io, feature.location)
    print(io, ")")
end

"""
    Base.show(io, record)

Render a compact summary for an annotated sequence record.
"""
function Base.show(io::IO, record::AnnotatedSeqRecord)
    print(io, "AnnotatedSeqRecord(", record.identifier, ", ", length(record.sequence), " bp, features=", length(record.features), ")")
end

"""
    feature_spans(location)

Return the genomic span covered by a feature location.
"""
function feature_spans(location::FeatureLocationLite)
    return [(location.start, location.stop)]
end

"""
    feature_spans(location)

Return the spans covered by each part of a compound feature location.
"""
function feature_spans(location::CompoundFeatureLocation)
    return [feature_bounds(part) for part in location.parts]
end

"""
    feature_spans(feature)

Return the spans represented by a feature record.
"""
feature_spans(feature::SeqFeatureLite) = feature_spans(feature.location)

"""
    feature_start(location)

Return the first coordinate covered by a feature location.
"""
function feature_start(location::AbstractFeatureLocation)
    return feature_bounds(location)[1]
end

"""
    feature_start(feature)

Return the first coordinate covered by a feature record.
"""
feature_start(feature::SeqFeatureLite) = feature_start(feature.location)

"""
    feature_stop(location)

Return the last coordinate covered by a feature location.
"""
function feature_stop(location::AbstractFeatureLocation)
    return feature_bounds(location)[2]
end

"""
    feature_stop(feature)

Return the last coordinate covered by a feature record.
"""
feature_stop(feature::SeqFeatureLite) = feature_stop(feature.location)

feature_length(feature::SeqFeatureLite) = feature_length(feature.location)

"""
    feature_strand(location)

Return the strand encoded by a simple feature location.
"""
function feature_strand(location::FeatureLocationLite)
    return location.strand
end

"""
    feature_strand(location)

Return the strand encoded by a compound feature location.
"""
function feature_strand(location::CompoundFeatureLocation)
    return location.strand
end

"""
    feature_strand(feature)

Return the strand encoded by a feature record.
"""
feature_strand(feature::SeqFeatureLite) = feature_strand(feature.location)

"""
    feature_identifier(feature)

Return the best available identifier for a feature record.
"""
function feature_identifier(feature::SeqFeatureLite)
    return isempty(feature.id) ? _feature_identifier(feature.qualifiers) : feature.id
end

"""
    feature_annotations(feature)

Return the qualifier dictionary for a feature record.
"""
function feature_annotations(feature::SeqFeatureLite)
    return Dict{String,Vector{String}}(feature.qualifiers)
end

"""
    feature_annotation(feature, key; default=nothing)

Return the first qualifier value stored under `key`.
"""
function feature_annotation(feature::SeqFeatureLite, key::AbstractString; default=nothing)
    values = get(feature.qualifiers, String(key), nothing)
    values === nothing && return default
    isempty(values) && return default
    return values[1]
end

"""
    feature_annotation(feature, key; default=nothing)

Symbol-based convenience wrapper for `feature_annotation`.
"""
feature_annotation(feature::SeqFeatureLite, key::Symbol; default=nothing) = feature_annotation(feature, String(key); default=default)

"""
    DataFrames.DataFrame(features)

Convert a vector of features into a tabular DataFrame.
"""
function DataFrames.DataFrame(features::AbstractVector{<:SeqFeatureLite})
    rows = collect(features)
    return DataFrames.DataFrame(
        feature_type = [feature.feature_type for feature in rows],
        start = [feature_start(feature) for feature in rows],
        stop = [feature_stop(feature) for feature in rows],
        strand = [feature_strand(feature) for feature in rows],
        id = [feature_identifier(feature) for feature in rows],
    )
end

"""
    feature_contains(location, position)

Test whether a position falls inside a simple feature location.
"""
function feature_contains(location::FeatureLocationLite, position::Integer)
    start, stop = feature_bounds(location)
    lower = min(start, stop)
    upper = max(start, stop)
    return lower <= position <= upper
end

"""
    feature_contains(location, position)

Test whether a position falls inside any part of a compound feature location.
"""
function feature_contains(location::CompoundFeatureLocation, position::Integer)
    return any(part -> feature_contains(part, position), location.parts)
end

"""
    feature_contains(feature, position)

Test whether a position falls inside a feature record.
"""
feature_contains(feature::SeqFeatureLite, position::Integer) = feature_contains(feature.location, position)

"""
    feature_overlaps(left, right)

Test whether two feature locations overlap.
"""
function feature_overlaps(left::AbstractFeatureLocation, right::AbstractFeatureLocation)
    if left isa CompoundFeatureLocation
        return any(part -> feature_overlaps(part, right), left.parts)
    end

    if right isa CompoundFeatureLocation
        return any(part -> feature_overlaps(left, part), right.parts)
    end

    left_start, left_stop = feature_bounds(left::FeatureLocationLite)
    right_start, right_stop = feature_bounds(right::FeatureLocationLite)
    return max(left_start, right_start) <= min(left_stop, right_stop)
end

"""
    feature_overlaps(left, right)

Test whether two feature records or a record and location overlap.
"""
feature_overlaps(left::SeqFeatureLite, right::SeqFeatureLite) = feature_overlaps(left.location, right.location)
"""
    feature_overlaps(feature, location)

Test whether a feature record overlaps a feature location.
"""
feature_overlaps(feature::SeqFeatureLite, location::AbstractFeatureLocation) = feature_overlaps(feature.location, location)
"""
    feature_overlaps(location, feature)

Test whether a feature location overlaps a feature record.
"""
feature_overlaps(location::AbstractFeatureLocation, feature::SeqFeatureLite) = feature_overlaps(location, feature.location)

"""
    feature_extract(sequence, location)

Extract the subsequence covered by a feature location.
"""
function feature_extract(sequence::AbstractString, location::AbstractFeatureLocation)
    return feature_sequence(sequence, location)
end

"""
    feature_extract(record, feature)

Extract the subsequence covered by a feature record from an annotated record.
"""
feature_extract(record::AnnotatedSeqRecord, feature::SeqFeatureLite) = feature_sequence(record, feature)

"""
    feature_slice(record, slice_start, slice_stop; reverse_complemented=false)

Slice a sequence record while preserving annotation metadata where possible.
"""
function feature_slice(record::SeqRecordLite, slice_start::Integer, slice_stop::Integer; reverse_complemented::Bool=false)
    slice_start <= slice_stop || throw(ArgumentError("slice_start must be <= slice_stop"))
    start_index = max(1, Int(slice_start))
    stop_index = min(lastindex(record.sequence), Int(slice_stop))
    start_index > stop_index && return SeqRecordLite("")

    sliced_sequence = record.sequence[start_index:stop_index]
    if reverse_complemented
        sliced_sequence = reverse_complement(sliced_sequence)
    end

    sliced_letter_annotations = Dict{Symbol,Any}()
    for (key, value) in record.letter_annotations
        sliced_letter_annotations[key] = _slice_letter_annotation(value, start_index, stop_index; reverse_complemented=reverse_complemented)
    end

    return SeqRecordLite(
        sliced_sequence;
        identifier=record.identifier,
        name=record.name,
        description=record.description,
        annotations=copy(record.annotations),
        letter_annotations=sliced_letter_annotations,
    )
end

"""
    feature_slice(record, slice_start, slice_stop; reverse_complemented=false)

Slice an annotated record and remap its features into the sliced coordinate system.
"""
function feature_slice(record::AnnotatedSeqRecord, slice_start::Integer, slice_stop::Integer; reverse_complemented::Bool=false)
    return slice_annotated_record(record, slice_start, slice_stop; reverse_complemented=reverse_complemented)
end

"""
    feature_slice(record, slice_start, slice_stop; reverse_complemented=false)

Slice a GenBank record via its annotated-record representation.
"""
function feature_slice(record::GenBankRecord, slice_start::Integer, slice_stop::Integer; reverse_complemented::Bool=false)
    return feature_slice(annotate_genbank_record(record), slice_start, slice_stop; reverse_complemented=reverse_complemented)
end

feature_slice(record::SeqRecordLite, range::UnitRange{<:Integer}) = feature_slice(record, first(range), last(range))
feature_slice(record::AnnotatedSeqRecord, range::UnitRange{<:Integer}) = feature_slice(record, first(range), last(range))
feature_slice(record::GenBankRecord, range::UnitRange{<:Integer}) = feature_slice(record, first(range), last(range))

"""
    _span_bounds(span)

Extract start and stop coordinates from a named span.
"""
function _span_bounds(span::NamedTuple)
    hasproperty(span, :start) && hasproperty(span, :stop) || throw(ArgumentError("named spans must provide start and stop fields"))
    return Int(getproperty(span, :start)), Int(getproperty(span, :stop))
end

feature_slice(record::SeqRecordLite, span::NamedTuple; reverse_complemented::Bool=false) = feature_slice(record, _span_bounds(span)...; reverse_complemented=reverse_complemented)
feature_slice(record::AnnotatedSeqRecord, span::NamedTuple; reverse_complemented::Bool=false) = feature_slice(record, _span_bounds(span)...; reverse_complemented=reverse_complemented)
feature_slice(record::GenBankRecord, span::NamedTuple; reverse_complemented::Bool=false) = feature_slice(record, _span_bounds(span)...; reverse_complemented=reverse_complemented)

"""
    _slice_location_bounds(start, stop, slice_start, slice_stop)

Clip a feature location to a slice and return the overlapping bounds.
"""
function _slice_location_bounds(start::Int, stop::Int, slice_start::Int, slice_stop::Int)
    lower = max(min(start, stop), slice_start)
    upper = min(max(start, stop), slice_stop)
    lower > upper && return nothing
    return lower, upper
end

"""
    slice_feature_location(location, slice_start, slice_stop; reverse_complemented=false)

Remap a simple feature location into a sliced coordinate system.
"""
function slice_feature_location(location::FeatureLocationLite, slice_start::Integer, slice_stop::Integer; reverse_complemented::Bool=false)
    clipped = _slice_location_bounds(location.start, location.stop, Int(slice_start), Int(slice_stop))
    clipped === nothing && return nothing

    clipped_start, clipped_stop = clipped
    relative_start = clipped_start - Int(slice_start) + 1
    relative_stop = clipped_stop - Int(slice_start) + 1
    if reverse_complemented
        slice_length = Int(slice_stop) - Int(slice_start) + 1
        relative_start, relative_stop = slice_length - relative_stop + 1, slice_length - relative_start + 1
        return FeatureLocationLite(relative_start, relative_stop; strand=-location.strand, partial_start=location.partial_start || clipped_start > location.start, partial_stop=location.partial_stop || clipped_stop < location.stop)
    end

    return FeatureLocationLite(relative_start, relative_stop; strand=location.strand, partial_start=location.partial_start || clipped_start > location.start, partial_stop=location.partial_stop || clipped_stop < location.stop)
end

"""
    slice_feature_location(location, slice_start, slice_stop; reverse_complemented=false)

Remap a compound feature location into a sliced coordinate system.
"""
function slice_feature_location(location::CompoundFeatureLocation, slice_start::Integer, slice_stop::Integer; reverse_complemented::Bool=false)
    transformed = AbstractFeatureLocation[]
    parts = reverse_complemented ? Iterators.reverse(location.parts) : location.parts

    for part in parts
        sliced = slice_feature_location(part, slice_start, slice_stop; reverse_complemented=reverse_complemented)
        sliced === nothing && continue
        push!(transformed, sliced)
    end

    isempty(transformed) && return nothing
    return CompoundFeatureLocation(location.operator, transformed; strand=reverse_complemented ? -location.strand : location.strand)
end

"""
    slice_feature(feature, slice_start, slice_stop; reverse_complemented=false)

Remap a feature record into a sliced coordinate system.
"""
function slice_feature(feature::SeqFeatureLite, slice_start::Integer, slice_stop::Integer; reverse_complemented::Bool=false)
    sliced_location = slice_feature_location(feature.location, slice_start, slice_stop; reverse_complemented=reverse_complemented)
    sliced_location === nothing && return nothing
    return SeqFeatureLite(feature.feature_type, sliced_location, feature.qualifiers, feature.id)
end

"""
    _slice_letter_annotation(value, slice_start, slice_stop; reverse_complemented=false)

Slice per-letter annotations in the same span as the parent record.
"""
function _slice_letter_annotation(value, slice_start::Int, slice_stop::Int; reverse_complemented::Bool=false)
    value isa AbstractString || return value
    lastindex(value) < slice_stop && return value
    sliced = value[slice_start:slice_stop]
    return reverse_complemented ? reverse(sliced) : sliced
end

"""
    slice_annotated_record(record, slice_start, slice_stop; reverse_complemented=false)

Slice an annotated record while preserving metadata and remapping features.
"""
function slice_annotated_record(record::AnnotatedSeqRecord, slice_start::Integer, slice_stop::Integer; reverse_complemented::Bool=false)
    slice_start <= slice_stop || throw(ArgumentError("slice_start must be <= slice_stop"))
    start_index = max(1, Int(slice_start))
    stop_index = min(lastindex(record.sequence), Int(slice_stop))
    start_index > stop_index && return AnnotatedSeqRecord("")

    sliced_sequence = record.sequence[start_index:stop_index]
    if reverse_complemented
        sliced_sequence = reverse_complement(sliced_sequence)
    end

    sliced_features = SeqFeatureLite[]
    feature_iterator = reverse_complemented ? Iterators.reverse(record.features) : record.features
    for feature in feature_iterator
        sliced_feature = slice_feature(feature, start_index, stop_index; reverse_complemented=reverse_complemented)
        sliced_feature === nothing && continue
        push!(sliced_features, sliced_feature)
    end

    sliced_letter_annotations = Dict{Symbol,Any}()
    for (key, value) in record.letter_annotations
        sliced_letter_annotations[key] = _slice_letter_annotation(value, start_index, stop_index; reverse_complemented=reverse_complemented)
    end

    return AnnotatedSeqRecord(
        sliced_sequence;
        identifier=record.identifier,
        name=record.name,
        description=record.description,
        annotations=copy(record.annotations),
        letter_annotations=sliced_letter_annotations,
        features=sliced_features,
    )
end

Base.getindex(record::AnnotatedSeqRecord, range::UnitRange{<:Integer}) = slice_annotated_record(record, first(range), last(range))

"""
    reverse_complement(record)

Return a reverse-complemented annotated record with remapped features.
"""
function reverse_complement(record::AnnotatedSeqRecord)
    return slice_annotated_record(record, 1, length(record.sequence); reverse_complemented=true)
end

"""
    feature_summary(feature)

Return a short text summary of a feature record.
"""
function feature_summary(feature::SeqFeatureLite)
    return (
        feature_type = feature.feature_type,
        identifier = feature_identifier(feature),
        strand = feature_strand(feature),
        length = feature_length(feature),
        bounds = feature_bounds(feature.location),
        spans = feature_spans(feature),
        qualifiers = feature_annotations(feature),
    )
end

"""
    _feature_matches_region(feature, region)

Test whether a feature matches a query region.
"""
function _feature_matches_region(feature::SeqFeatureLite, region::AbstractFeatureLocation)
    return feature_overlaps(feature.location, region)
end

"""
    _feature_matches_region(feature, region)

Test whether a feature contains a single query position.
"""
function _feature_matches_region(feature::SeqFeatureLite, region::Integer)
    return feature_contains(feature.location, region)
end

"""
    _feature_matches_region(feature, region)

Test whether a feature overlaps a query interval.
"""
function _feature_matches_region(feature::SeqFeatureLite, region::Tuple{<:Integer,<:Integer})
    return feature_overlaps(feature.location, FeatureLocationLite(region[1], region[2]))
end

"""
    select_features(record; kwargs...)

Filter annotated features by type, region, strand, identifier, or overlap.
"""
function select_features(
    record::AnnotatedSeqRecord;
    feature_type=nothing,
    region=nothing,
    strand=nothing,
    qualifier_key=nothing,
    qualifier_value=nothing,
)
    selected = SeqFeatureLite[]
    for feature in record.features
        if feature_type !== nothing
            if feature_type isa AbstractVector
                any(candidate -> feature.feature_type == String(candidate), feature_type) || continue
            else
                feature.feature_type == String(feature_type) || continue
            end
        end

        if strand !== nothing && feature_strand(feature) != Int8(strand)
            continue
        end

        if region !== nothing && !_feature_matches_region(feature, region)
            continue
        end

        if qualifier_key !== nothing
            qualifier_values = get(feature.qualifiers, String(qualifier_key), nothing)
            qualifier_values === nothing && continue
            if qualifier_value !== nothing && all(value -> value != String(qualifier_value), qualifier_values)
                continue
            end
        elseif qualifier_value !== nothing
            any(qualifier_values -> any(value -> value == String(qualifier_value), qualifier_values), values(feature.qualifiers)) || continue
        end

        push!(selected, feature)
    end

    return selected
end

features_at(record::AnnotatedSeqRecord, position::Integer) = select_features(record; region=position)

"""
    features_overlapping(record, region)

Return the features in a record that overlap a query location.
"""
function features_overlapping(record::AnnotatedSeqRecord, region::AbstractFeatureLocation)
    return select_features(record; region=region)
end

"""
    features_overlapping(record, start, stop)

Return the features in a record that overlap a coordinate interval.
"""
function features_overlapping(record::AnnotatedSeqRecord, start::Integer, stop::Integer)
    return select_features(record; region=(start, stop))
end

"""
    feature_table(record)

Convert annotated features into a tabular DataFrame.
"""
function feature_table(record::AnnotatedSeqRecord)
    return [feature_summary(feature) for feature in record.features]
end

"""
    _split_top_level(text)

Split a nested feature-location expression at top-level separators.
"""
function _split_top_level(text::AbstractString)
    pieces = String[]
    depth = 0
    start_index = firstindex(text)

    for (index, char) in pairs(text)
        if char == '('
            depth += 1
        elseif char == ')'
            depth -= 1
        elseif char == ',' && depth == 0
            push!(pieces, strip(text[start_index:prevind(text, index)]))
            start_index = nextind(text, index)
        end
    end

    push!(pieces, strip(text[start_index:lastindex(text)]))
    filter!(piece -> !isempty(piece), pieces)
    return pieces
end

const _FEATURE_LOCATION_CACHE = Dict{String,AbstractFeatureLocation}()

"""
    _parse_feature_location_uncached(location)

Parse a feature location string without consulting the memoized cache.
"""
function _parse_feature_location_uncached(location::AbstractString)
    stripped = strip(location)
    isempty(stripped) && throw(ArgumentError("empty feature location"))

    if startswith(stripped, "complement(") && endswith(stripped, ")")
        inner = stripped[12:end-1]
        return _with_strand(_parse_feature_location_uncached(inner), Int8(-1))
    end

    if (startswith(stripped, "join(") || startswith(stripped, "order(")) && endswith(stripped, ")")
        operator = startswith(stripped, "join(") ? "join" : "order"
        inner = stripped[length(operator) + 2:end-1]
        parts = AbstractFeatureLocation[_parse_feature_location_uncached(part) for part in _split_top_level(inner)]
        return CompoundFeatureLocation(operator, parts; strand=1)
    end

    if occursin("^", stripped)
        position_text = replace(stripped, "^" => "")
        position = parse(Int, replace(position_text, r"[^0-9]" => ""))
        return FeatureLocationLite(position, position)
    end

    if occursin("..", stripped)
        start_part, stop_part = Base.split(stripped, "..", limit=2)
        start_text = strip(start_part)
        stop_text = strip(stop_part)
        partial_start = startswith(start_text, "<")
        partial_stop = startswith(stop_text, ">")
        start_value = parse(Int, replace(start_text, r"^[<>]" => ""))
        stop_value = parse(Int, replace(stop_text, r"^[<>]" => ""))
        return FeatureLocationLite(start_value, stop_value; partial_start=partial_start, partial_stop=partial_stop)
    end

    partial_start = startswith(stripped, "<")
    partial_stop = startswith(stripped, ">")
    position_text = strip(replace(stripped, r"^[<>]" => ""))
    position = parse(Int, position_text)
    return FeatureLocationLite(position, position; partial_start=partial_start, partial_stop=partial_stop)
end

"""
    _with_strand(location, strand)

Return a copy of a feature location with a new strand.
"""
function _with_strand(location::AbstractFeatureLocation, strand::Int8)
    if location isa FeatureLocationLite
        simple = location::FeatureLocationLite
        return FeatureLocationLite(simple.start, simple.stop; strand=strand, partial_start=simple.partial_start, partial_stop=simple.partial_stop)
    end

    compound = location::CompoundFeatureLocation
    return CompoundFeatureLocation(compound.operator, compound.parts; strand=strand)
end

"""
    parse_feature_location(location)

Parse a GenBank- or GFF-style feature location string into a typed location.
"""
function parse_feature_location(location::AbstractString)
    key = String(location)
    return get!(_FEATURE_LOCATION_CACHE, key) do
        _parse_feature_location_uncached(key)
    end
end

"""
    feature_bounds(location)

Return the outer coordinate bounds for a simple feature location.
"""
function feature_bounds(location::FeatureLocationLite)
    return location.start, location.stop
end

"""
    feature_bounds(location)

Return the outer coordinate bounds for a compound feature location.
"""
function feature_bounds(location::CompoundFeatureLocation)
    starts = Int[]
    stops = Int[]
    for part in location.parts
        start, stop = feature_bounds(part)
        push!(starts, start)
        push!(stops, stop)
    end
    return minimum(starts), maximum(stops)
end

"""
    feature_length(location)

Return the span length of a simple feature location.
"""
function feature_length(location::FeatureLocationLite)
    return abs(location.stop - location.start) + 1
end

"""
    feature_length(location)

Return the total span length of a compound feature location.
"""
function feature_length(location::CompoundFeatureLocation)
    return sum(feature_length(part) for part in location.parts)
end

"""
    feature_sequence(sequence, location)

Extract the sequence segment covered by a simple feature location.
"""
function feature_sequence(sequence::AbstractString, location::FeatureLocationLite)
    start = max(1, location.start)
    stop = min(lastindex(sequence), location.stop)
    start > stop && return ""
    subsequence = sequence[start:stop]
    return location.strand == -1 ? reverse_complement(subsequence) : subsequence
end

"""
    feature_sequence(sequence, location)

Extract the sequence segments covered by a compound feature location.
"""
function feature_sequence(sequence::AbstractString, location::CompoundFeatureLocation)
    extracted = join(feature_sequence(sequence, part) for part in location.parts)
    return location.strand == -1 ? reverse_complement(extracted) : extracted
end

"""
    _feature_identifier(qualifiers)

Derive a stable feature identifier from qualifier values.
"""
function _feature_identifier(qualifiers::Dict{String,Vector{String}})
    for key in ("gene", "locus_tag", "protein_id", "ID")
        haskey(qualifiers, key) && !isempty(qualifiers[key]) && return qualifiers[key][1]
    end
    return ""
end

"""
    SeqFeatureLite(feature::GenBankFeature)

Convert a parsed GenBank feature into a lightweight feature record.
"""
function SeqFeatureLite(feature::GenBankFeature)
    parsed_location = feature.parsed_location === nothing ? parse_feature_location(feature.location) : feature.parsed_location
    return SeqFeatureLite(
        feature.key,
        parsed_location,
        feature.qualifiers,
        _feature_identifier(feature.qualifiers),
    )
end

"""
    SeqFeatureLite(record::GffRecord)

Convert a parsed GFF record into a lightweight feature record.
"""
function SeqFeatureLite(record::GffRecord)
    qualifiers = Dict{String,Vector{String}}(record.attribute_map)
    strand = record.strand == "-" ? -1 : 1
    location = FeatureLocationLite(record.start, record.stop; strand=strand)
    return SeqFeatureLite(record.feature, location, qualifiers, _feature_identifier(qualifiers))
end

"""
    annotate_gff_records(records)

Convert parsed GFF records into lightweight feature records.
"""
function annotate_gff_records(records::AbstractVector{GffRecord})
    return AnnotatedSeqRecord[
        AnnotatedSeqRecord(
            "";
            identifier=record.chrom,
            name=record.source,
            description=record.feature,
            annotations=Dict{Symbol,Any}(
                :chrom => record.chrom,
                :source => record.source,
                :score => record.score,
                :strand => record.strand,
                :phase => record.phase,
                :attributes => record.attributes,
            ),
            features=[SeqFeatureLite(record)],
        ) for record in records
    ]
end

"""
    annotate_genbank_record(record)

Convert a parsed GenBank record into an annotated sequence record.
"""
function annotate_genbank_record(record::GenBankRecord)
    annotations = Dict{Symbol,Any}(
        :accessions => record.accession == "" ? String[] : [record.accession],
        :version => record.version,
        :source => record.source,
        :organism => record.organism,
        :definition => record.definition,
        :locus => record.locus,
    )
    features = Vector{SeqFeatureLite}(undef, length(record.features))
    @inbounds for index in eachindex(record.features)
        features[index] = SeqFeatureLite(record.features[index])
    end
    return AnnotatedSeqRecord(
        record.sequence;
        identifier=record.accession == "" ? record.locus : record.accession,
        name=record.locus,
        description=record.definition,
        annotations=annotations,
        features=features,
    )
end

"""
    annotate_genbank_records(records)

Convert parsed GenBank records into annotated sequence records.
"""
function annotate_genbank_records(records::AbstractVector{GenBankRecord})
    return AnnotatedSeqRecord[annotate_genbank_record(record) for record in records]
end

"""
    feature_sequence(record, feature)

Extract the sequence for a feature stored on an annotated record.
"""
function feature_sequence(record::AnnotatedSeqRecord, feature::SeqFeatureLite)
    return feature_sequence(record.sequence, feature.location)
end

"""
    feature_sequence(record, feature)

Extract the sequence for a feature on a parsed GenBank record.
"""
function feature_sequence(record::GenBankRecord, feature::GenBankFeature)
    parsed_location = feature.parsed_location === nothing ? parse_feature_location(feature.location) : feature.parsed_location
    return feature_sequence(record.sequence, parsed_location)
end

"""
    _feature_chrom(feature)

Read the chromosome or contig field from a GFF record.
"""
function _feature_chrom(feature::GffRecord)
    return feature.chrom
end

"""
    _feature_start(feature)

Read the start coordinate from a GFF record.
"""
function _feature_start(feature::GffRecord)
    return Int(feature.start)
end

"""
    _feature_stop(feature)

Read the stop coordinate from a GFF record.
"""
function _feature_stop(feature::GffRecord)
    return Int(feature.stop)
end

"""
    _feature_type(feature)

Read the feature type from a GFF record.
"""
function _feature_type(feature::GffRecord)
    return feature.feature
end

"""
    _feature_strand(feature)

Read the strand from a GFF record.
"""
function _feature_strand(feature::GffRecord)
    return feature.strand
end

"""
    _feature_identifier(feature)

Derive a stable identifier from a GFF record.
"""
function _feature_identifier(feature::GffRecord)
    for key in ("gene", "Name", "ID", "locus_tag")
        values = get(feature.attribute_map, key, nothing)
        values === nothing && continue
        isempty(values) && continue
        return values[1]
    end
    return ""
end

"""
    _variant_field(record, field)

Fetch a named field from a variant record with several fallback spellings.
"""
function _variant_field(record, field::Symbol)
    hasproperty(record, field) || throw(ArgumentError("variant record is missing $(field)"))
    return getproperty(record, field)
end

"""
    _variant_effect_rank(effect)

Rank a predicted variant effect by severity.
"""
function _variant_effect_rank(effect::AbstractString)
    lowered = lowercase(String(effect))
    return get(Dict(
        "stop-gain" => 6,
        "stop-loss" => 5,
        "missense" => 4,
        "synonymous" => 3,
        "coding" => 2,
        "utr" => 1,
        "intron" => 1,
        "intergenic" => 0,
    ), lowered, 2)
end

"""
    _complement_base(base)

Return the DNA complement of a nucleotide character.
"""
function _complement_base(base::Char)
    upper = uppercase(base)
    upper == 'A' && return 'T'
    upper == 'C' && return 'G'
    upper == 'G' && return 'C'
    upper == 'T' && return 'A'
    upper == 'U' && return 'A'
    return 'N'
end

"""
    _transcript_base(base, strand)

Return the strand-aware transcript base for a genomic nucleotide.
"""
function _transcript_base(base::Char, strand::AbstractString)
    return strand == "-" ? _complement_base(base) : uppercase(base)
end

"""
    _codon_effect(variant_ref, variant_alt, genomic_sequence, feature, pos)

Compute a codon-level consequence for a variant intersecting a coding feature.
"""
function _codon_effect(variant_ref::AbstractString, variant_alt::AbstractString, genomic_sequence::AbstractString, feature::GffRecord, pos::Integer)
    coding_sequence = feature_sequence(genomic_sequence, FeatureLocationLite(feature.start, feature.stop; strand=feature.strand == "-" ? -1 : 1))
    coding_sequence = uppercase(coding_sequence)
    coding_sequence = feature.strand == "-" ? reverse_complement(coding_sequence) : coding_sequence

    cdna_position = feature.strand == "-" ? feature.stop - Int(pos) + 1 : Int(pos) - feature.start + 1
    cdna_position < 1 && return (effect = "Coding", codon_ref = "", codon_alt = "")
    codon_start = 3 * div(cdna_position - 1, 3) + 1
    codon_start + 2 > lastindex(coding_sequence) && return (effect = "Coding", codon_ref = "", codon_alt = "")

    ref_codon = coding_sequence[codon_start:codon_start+2]
    codon_vector = collect(ref_codon)
    codon_index = cdna_position - codon_start + 1
    alt_base = _transcript_base(first(variant_alt), feature.strand)
    codon_vector[codon_index] = alt_base
    alt_codon = String(codon_vector)

    ref_aa = translate_dna(ref_codon; stop_at_stop=true)
    alt_aa = translate_dna(alt_codon; stop_at_stop=true)
    if ref_aa != alt_aa
        if alt_aa == "*"
            return (effect = "Stop-Gain", codon_ref = ref_codon, codon_alt = alt_codon)
        elseif ref_aa == "*"
            return (effect = "Stop-Loss", codon_ref = ref_codon, codon_alt = alt_codon)
        else
            return (effect = "Missense", codon_ref = ref_codon, codon_alt = alt_codon)
        end
    end

    return (effect = "Synonymous", codon_ref = ref_codon, codon_alt = alt_codon)
end

"""
    _feature_consequence(variant, feature; reference_sequences=nothing)

Compute the predicted consequence of a variant for a single annotated feature.
"""
function _feature_consequence(variant, feature::GffRecord; reference_sequences=nothing)
    chrom = String(_variant_field(variant, :chrom))
    position = Int(_variant_field(variant, :pos))
    variant_ref = uppercase(String(_variant_field(variant, :ref)))
    variant_alt = uppercase(String(_variant_field(variant, :alt)))

    (position < _feature_start(feature) || position > _feature_stop(feature)) && return nothing

    feature_type = lowercase(_feature_type(feature))
    gene = _feature_identifier(feature)

    if occursin("utr", feature_type)
        return (gene=gene, feature_type=_feature_type(feature), consequence="UTR", codon_ref="", codon_alt="")
    elseif occursin("intron", feature_type) || occursin("splice", feature_type)
        return (gene=gene, feature_type=_feature_type(feature), consequence="Intron", codon_ref="", codon_alt="")
    elseif occursin("cds", feature_type) || occursin("coding", feature_type)
        if reference_sequences !== nothing && haskey(reference_sequences, chrom)
            sequence = String(reference_sequences[chrom])
            codon_effect = _codon_effect(variant_ref, variant_alt, sequence, feature, position)
            return (gene=gene, feature_type=_feature_type(feature), consequence=codon_effect.effect, codon_ref=codon_effect.codon_ref, codon_alt=codon_effect.codon_alt)
        end
        return (gene=gene, feature_type=_feature_type(feature), consequence="Coding", codon_ref="", codon_alt="")
    else
        return (gene=gene, feature_type=_feature_type(feature), consequence="Gene", codon_ref="", codon_alt="")
    end
end

"""
    annotate_variants(variant_records, gene_features; reference_sequences=nothing)

Annotate variants against gene features and return their predicted consequences.
"""
function annotate_variants(variant_records::AbstractVector, gene_features::AbstractVector; reference_sequences=nothing)
    annotations = NamedTuple[]
    for variant in variant_records
        best = (gene="", feature_type="", consequence="Intergenic", codon_ref="", codon_alt="")
        best_rank = 0
        chrom = String(_variant_field(variant, :chrom))
        pos = Int(_variant_field(variant, :pos))

        for feature in gene_features
            feature isa GffRecord || continue
            _feature_chrom(feature) == chrom || continue
            pos < _feature_start(feature) || pos > _feature_stop(feature) && continue
            consequence = _feature_consequence(variant, feature; reference_sequences=reference_sequences)
            consequence === nothing && continue
            rank = _variant_effect_rank(consequence.consequence)
            if rank > best_rank
                best = consequence
                best_rank = rank
            end
        end

        push!(annotations, (
            chrom = chrom,
            pos = pos,
            ref = String(_variant_field(variant, :ref)),
            alt = String(_variant_field(variant, :alt)),
            gene = best.gene,
            feature_type = best.feature_type,
            consequence = best.consequence,
            codon_ref = best.codon_ref,
            codon_alt = best.codon_alt,
        ))
    end

    return annotations
end
