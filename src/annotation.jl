using DataFrames
using .BioToolkit: ProvenanceParams, ThreadSafeProvenanceContext, new_provenance_id

export AbstractFeatureLocation, FeatureLocationLite, CompoundFeatureLocation, SeqFeatureLite, AnnotatedSeqRecord, feature_spans, feature_bounds, feature_start, feature_stop, parse_feature_location, SangerTrace,
    AnnotatedSeqIndex, build_feature_index, indexed_features_overlapping, indexed_features_at,
    feature_overlaps_stranded, feature_distance, nearest_feature, preceding_features, following_features,
    feature_coverage, copy_feature_annotations,
    visualize_annotation_record_html, visualize_sanger_trace_html, visualize_variant_consequences_html

@inline function _register_annotation_result!(_explicit_ctx, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    _ctx = active_provenance_context(_explicit_ctx)
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

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
    # NOTE: parts is Vector{AbstractFeatureLocation} for generality.
    # In hot paths, use dispatch on FeatureLocationLite directly where possible.
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
mutable struct AnnotatedSeqRecord{A <: BioAlphabet}
    sequence::BioSequence{A}
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
    sequence::BioSequence{DNAAlphabet}
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
CompoundFeatureLocation(operator::String, parts::AbstractVector{<:AbstractFeatureLocation}; strand::Integer=1) = CompoundFeatureLocation(String(operator), AbstractFeatureLocation[part for part in parts], Int8(strand))

"""
    SeqFeatureLite(feature_type, location; qualifiers=..., id="")

Construct a lightweight feature record from a type string and location.
"""
function SeqFeatureLite(feature_type::String, location::AbstractFeatureLocation; qualifiers::AbstractDict=Dict{String,Vector{String}}(), id::String="")
    return SeqFeatureLite(String(feature_type), location, Dict{String,Vector{String}}(qualifiers), String(id))
end

"""
    AnnotatedSeqRecord(sequence; kwargs...)

Construct a mutable annotated record from raw sequence and associated metadata.
"""
function AnnotatedSeqRecord(
    sequence::BioSequence{A};
    identifier::String="",
    name::String=identifier,
    description::String=name,
    annotations::AbstractDict=Dict{Symbol,Any}(),
    letter_annotations::AbstractDict=Dict{Symbol,Any}(),
    features::AbstractVector{SeqFeatureLite}=SeqFeatureLite[],
) where {A <: BioAlphabet}
    feature_vector = features isa Vector{SeqFeatureLite} ? features : SeqFeatureLite[feature for feature in features]
    normalized_letter_annotations = Dict{Symbol,Any}(letter_annotations)
    for (key, value) in normalized_letter_annotations
        if value isa String
            length(value) == length(sequence) || throw(ArgumentError("letter annotation '$key' length $(length(value)) must match sequence length $(length(sequence))"))
        elseif value isa AbstractVector
            length(value) == length(sequence) || throw(ArgumentError("letter annotation '$key' length $(length(value)) must match sequence length $(length(sequence))"))
        end
    end
    annotations_copy = Dict{Symbol,Any}(annotations)
    ensure_provenance_id!(annotations_copy)
    return AnnotatedSeqRecord{A}(
        sequence,
        String(identifier),
        String(name),
        String(description),
        annotations_copy,
        normalized_letter_annotations,
        feature_vector,
    )
end

"""
    Base.length(record)

Return the sequence length of an annotated record in code units.
"""
Base.length(record::AnnotatedSeqRecord) = length(record.sequence)

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
    print(io, "AnnotatedSeqRecord(", record.identifier, ", ", length(record), " bp, features=", length(record.features), ", ", container_provenance_summary(record), ")")
end

"""
    feature_spans(location)

Return the genomic span covered by a feature location.
"""
function feature_spans(location::FeatureLocationLite; prov_ctx=nothing)
    result = [(location.start, location.stop)]
    _ctx = active_provenance_context(prov_ctx)


    return _register_annotation_result!(_ctx, result, "feature_spans"; parents=String[], parameters=(kind="FeatureLocationLite", span_count=length(result)))
end

"""
    feature_spans(location)

Return the spans covered by each part of a compound feature location.
"""
function feature_spans(location::CompoundFeatureLocation; prov_ctx=nothing)
    result = [feature_bounds(part) for part in location.parts]
    _ctx = active_provenance_context(prov_ctx)


    return _register_annotation_result!(_ctx, result, "feature_spans"; parents=String[], parameters=(kind="CompoundFeatureLocation", span_count=length(result)))
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
    result = feature_bounds(location)[1]
    _ctx = active_provenance_context()


    return _register_annotation_result!(_ctx, result, "feature_start"; parents=String[], parameters=(start=result))
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
    result = feature_bounds(location)[2]
    _ctx = active_provenance_context()


    return _register_annotation_result!(_ctx, result, "feature_stop"; parents=String[], parameters=(stop=result))
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

Return the qualifier dictionary for a feature record (read-only view).
Use `copy_feature_annotations(feature)` if you need a mutable copy.
"""
function feature_annotations(feature::SeqFeatureLite)
    return feature.qualifiers
end

"""
    copy_feature_annotations(feature)

Return a mutable deep copy of the qualifier dictionary for a feature record.
"""
function copy_feature_annotations(feature::SeqFeatureLite)
    return Dict{String,Vector{String}}(k => copy(v) for (k, v) in feature.qualifiers)
end

"""
    feature_annotation(feature, key; default=nothing)

Return the first qualifier value stored under `key`.
"""
function feature_annotation(feature::SeqFeatureLite, key::String; default=nothing)
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
feature_extract(sequence::BioSequence, location::AbstractFeatureLocation) = feature_sequence(sequence, location)
feature_extract(sequence, location::AbstractFeatureLocation) = feature_sequence(BioSequence{DNAAlphabet}(String(sequence)), location)

"""
    feature_extract(record, feature)

Extract the subsequence covered by a feature record from an annotated record.
"""
feature_extract(record::AnnotatedSeqRecord, feature::SeqFeatureLite) = feature_sequence(record, feature)

"""
    feature_slice(record, slice_start, slice_stop; reverse_complemented=false)

Slice a sequence record while preserving annotation metadata where possible.
"""
function feature_slice(record::SeqRecordLite{A}, slice_start::Integer, slice_stop::Integer; reverse_complemented::Bool=false) where {A <: BioAlphabet}
    slice_start <= slice_stop || throw(ArgumentError("slice_start must be <= slice_stop"))
    start_index = max(1, Int(slice_start))
    stop_index = min(lastindex(record.sequence), Int(slice_stop))
    _ctx = active_provenance_context()
    if start_index > stop_index
        result = SeqRecordLite(
            BioSequence{A}(UInt8[]; validate=false);
            identifier=record.identifier,
            name=record.name,
            description=record.description,
            annotations=copy(record.annotations),
            letter_annotations=Dict{Symbol,Any}(),
        )
        return _register_annotation_result!(_ctx, result, "feature_slice"; parents=provenance_parent_ids(record), parameters=(slice_start=slice_start, slice_stop=slice_stop, reverse_complemented=reverse_complemented, empty=true))
    end

    sliced_sequence = record.sequence[start_index:stop_index]
    if reverse_complemented
        sliced_sequence = reverse_complement(sliced_sequence)
    end

    sliced_letter_annotations = Dict{Symbol,Any}()
    for (key, value) in record.letter_annotations
        sliced_letter_annotations[key] = _slice_letter_annotation(value, start_index, stop_index; reverse_complemented=reverse_complemented)
    end

    result = SeqRecordLite(
        sliced_sequence;
        identifier=record.identifier,
        name=record.name,
        description=record.description,
        annotations=copy(record.annotations),
        letter_annotations=sliced_letter_annotations,
    )

    return _register_annotation_result!(_ctx, result, "feature_slice"; parents=provenance_parent_ids(record), parameters=(slice_start=slice_start, slice_stop=slice_stop, reverse_complemented=reverse_complemented, empty=false))
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
        return FeatureLocationLite(relative_start, relative_stop; strand=-location.strand, partial_start=location.partial_stop || clipped_stop < location.stop, partial_stop=location.partial_start || clipped_start > location.start)
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
    if value isa String
        lastindex(value) < slice_stop && throw(ArgumentError("letter annotation shorter than slice range"))
        sliced = value[slice_start:slice_stop]
        return reverse_complemented ? reverse(sliced) : sliced
    elseif value isa AbstractVector
        lastindex(value) < slice_stop && throw(ArgumentError("letter annotation shorter than slice range"))
        sliced = value[slice_start:slice_stop]
        return reverse_complemented ? reverse(sliced) : sliced
    end
    return value
end

"""
    slice_annotated_record(record, slice_start, slice_stop; reverse_complemented=false)

Slice an annotated record while preserving metadata and remapping features.
"""
function slice_annotated_record(record::AnnotatedSeqRecord{A}, slice_start::Integer, slice_stop::Integer; reverse_complemented::Bool=false) where {A <: BioAlphabet}
    slice_start <= slice_stop || throw(ArgumentError("slice_start must be <= slice_stop"))
    start_index = max(1, Int(slice_start))
    stop_index = min(length(record.sequence), Int(slice_stop))
    _ctx = active_provenance_context()
    if start_index > stop_index
        result = AnnotatedSeqRecord(
            BioSequence{A}(UInt8[]; validate=false);
            identifier=record.identifier,
            name=record.name,
            description=record.description,
            annotations=copy(record.annotations),
            letter_annotations=Dict{Symbol,Any}(),
            features=SeqFeatureLite[],
        )
        return _register_annotation_result!(_ctx, result, "slice_annotated_record"; parents=provenance_parent_ids(record), parameters=(slice_start=slice_start, slice_stop=slice_stop, reverse_complemented=reverse_complemented, empty=true))
    end

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

    result = AnnotatedSeqRecord(
        sliced_sequence;
        identifier=record.identifier,
        name=record.name,
        description=record.description,
        annotations=copy(record.annotations),
        letter_annotations=sliced_letter_annotations,
        features=sliced_features,
    )

    return _register_annotation_result!(_ctx, result, "slice_annotated_record"; parents=provenance_parent_ids(record), parameters=(slice_start=slice_start, slice_stop=slice_stop, reverse_complemented=reverse_complemented, empty=false))
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
    result = (
        feature_type = feature.feature_type,
        identifier = feature_identifier(feature),
        strand = feature_strand(feature),
        length = feature_length(feature),
        bounds = feature_bounds(feature.location),
        spans = feature_spans(feature),
        qualifiers = feature_annotations(feature),
    )
    _ctx = active_provenance_context()


    return _register_annotation_result!(_ctx, result, "feature_summary"; parents=provenance_parent_ids(feature), parameters=(feature_type=feature.feature_type, identifier=feature.id))
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
    prov_ctx=nothing,
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

    _ctx = active_provenance_context(prov_ctx)


    return _register_annotation_result!(_ctx, selected, "select_features"; parents=provenance_parent_ids(record), parameters=(feature_count=length(selected), feature_type=feature_type === nothing ? "all" : string(feature_type), has_region=region !== nothing))
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


# feature_table is defined below with extended as_dataframe support


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

const _FEATURE_LOCATION_CACHE_MAX_SIZE = 8192
const _FEATURE_LOCATION_CACHE = Dict{String,AbstractFeatureLocation}()
const _FEATURE_LOCATION_CACHE_LOCK = ReentrantLock()

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
        left_text, right_text = Base.split(stripped, "^"; limit=2)
        left_value = parse(Int, replace(strip(left_text), r"[^0-9]" => ""))
        right_value = parse(Int, replace(strip(right_text), r"[^0-9]" => ""))
        return FeatureLocationLite(left_value, right_value)
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
function parse_feature_location(location::AbstractString; prov_ctx=nothing)
    key = String(location)
    cached = lock(_FEATURE_LOCATION_CACHE_LOCK) do
        get(_FEATURE_LOCATION_CACHE, key, nothing)
    end
    if cached === nothing
        parsed = _parse_feature_location_uncached(key)
        lock(_FEATURE_LOCATION_CACHE_LOCK) do
            if length(_FEATURE_LOCATION_CACHE) >= _FEATURE_LOCATION_CACHE_MAX_SIZE
                empty!(_FEATURE_LOCATION_CACHE)
            end
            _FEATURE_LOCATION_CACHE[key] = parsed
        end
        cached = parsed
    end
    _ctx = active_provenance_context(prov_ctx)

    return _register_annotation_result!(_ctx, cached, "parse_feature_location"; parents=String[], parameters=(location=key, location_type=string(typeof(cached))))
end

"""
    feature_bounds(location)

Return the outer coordinate bounds for a simple feature location.
"""
@inline function feature_bounds(location::FeatureLocationLite)
    result = (location.start, location.stop)
    _ctx = active_provenance_context()


    return _register_annotation_result!(_ctx, result, "feature_bounds"; parents=String[], parameters=(kind="FeatureLocationLite", start=location.start, stop=location.stop))
end

"""
    feature_bounds(location)

Return the outer coordinate bounds for a compound feature location.
"""
function feature_bounds(location::CompoundFeatureLocation)
    isempty(location.parts) && throw(ArgumentError("CompoundFeatureLocation must have at least one part"))
    lo = typemax(Int)
    hi = typemin(Int)
    for part in location.parts
        s, e = feature_bounds(part)
        lo = min(lo, s)
        hi = max(hi, e)
    end
    result = (lo, hi)
    _ctx = active_provenance_context()


    return _register_annotation_result!(_ctx, result, "feature_bounds"; parents=String[], parameters=(kind="CompoundFeatureLocation", start=result[1], stop=result[2]))
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
function feature_sequence(sequence::BioSequence, location::FeatureLocationLite)
    lower = min(location.start, location.stop)
    upper = max(location.start, location.stop)
    start = max(1, lower)
    stop = min(length(sequence), upper)
    start > stop && return BioSequence{alphabet(sequence)}(UInt8[]; validate=false)
    subsequence = sequence[start:stop]
    result = location.strand == -1 ? reverse_complement(subsequence) : subsequence
    _ctx = active_provenance_context()


    return _register_annotation_result!(_ctx, result, "feature_sequence"; parents=provenance_parent_ids(sequence), parameters=(kind="FeatureLocationLite", length=length(result), strand=location.strand))
end

"""
    feature_sequence(sequence, location)

Extract the sequence segments covered by a compound feature location.
"""
function feature_sequence(sequence::BioSequence, location::CompoundFeatureLocation)
    alphabet_type = alphabet(sequence)
    isempty(location.parts) && return BioSequence{alphabet_type}(UInt8[]; validate=false)
    # Pre-allocate buffer with total length to avoid vcat + splat overhead
    total_len = 0
    for part in location.parts
        total_len += feature_length(part)
    end
    buffer = Vector{UInt8}(undef, total_len)
    offset = 0
    for part in location.parts
        part_seq = feature_sequence(sequence, part)
        n = length(part_seq)
        @inbounds copyto!(buffer, offset + 1, part_seq.data, 1, n)
        offset += n
    end
    # Trim to actual bytes written (may be less if sequence bounds were clamped)
    resize!(buffer, offset)
    concatenated = BioSequence{alphabet_type}(buffer; validate=false)
    result = location.strand == -1 ? reverse_complement(concatenated) : concatenated
    _ctx = active_provenance_context()


    return _register_annotation_result!(_ctx, result, "feature_sequence"; parents=provenance_parent_ids(sequence), parameters=(kind="CompoundFeatureLocation", length=length(result), strand=location.strand))
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
    result = AnnotatedSeqRecord[
        AnnotatedSeqRecord(
            BioSequence{DNAAlphabet}(UInt8[]; validate=false);
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
    _ctx = active_provenance_context()


    return _register_annotation_result!(_ctx, result, "annotate_gff_records"; parents=provenance_parent_ids(records), parameters=(record_count=length(result)))
end

"""
    annotate_genbank_record(record)

Convert a parsed GenBank record into an annotated sequence record.
"""
function annotate_genbank_record(record::GenBankRecord; prov_ctx=nothing)
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
    result = AnnotatedSeqRecord(
        record.sequence;
        identifier=record.accession == "" ? record.locus : record.accession,
        name=record.locus,
        description=record.definition,
        annotations=annotations,
        features=features,
    )
    _ctx = active_provenance_context(prov_ctx)


    return _register_annotation_result!(_ctx, result, "annotate_genbank_record"; parents=provenance_parent_ids(record), parameters=(feature_count=length(features), identifier=result.identifier))
end

"""
    annotate_genbank_records(records)

Convert parsed GenBank records into annotated sequence records.
"""
function annotate_genbank_records(records::AbstractVector{GenBankRecord}; prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    result = AnnotatedSeqRecord[annotate_genbank_record(record; prov_ctx=_ctx) for record in records]


    return _register_annotation_result!(_ctx, result, "annotate_genbank_records"; parents=provenance_parent_ids(records), parameters=(record_count=length(result)))
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
const _VARIANT_EFFECT_RANKS = Dict{String,Int}(
    "stop-gain" => 6,
    "stop-loss" => 5,
    "missense" => 4,
    "synonymous" => 3,
    "coding" => 2,
    "utr" => 1,
    "intron" => 1,
    "intergenic" => 0,
)

@inline function _variant_effect_rank(effect::String)
    return get(_VARIANT_EFFECT_RANKS, lowercase(effect), 2)
end

const _COMPLEMENT_TABLE = let
    table = fill('N', 128)
    table[Int('A') + 1] = 'T'; table[Int('a') + 1] = 'T'
    table[Int('C') + 1] = 'G'; table[Int('c') + 1] = 'G'
    table[Int('G') + 1] = 'C'; table[Int('g') + 1] = 'C'
    table[Int('T') + 1] = 'A'; table[Int('t') + 1] = 'A'
    table[Int('U') + 1] = 'A'; table[Int('u') + 1] = 'A'
    Tuple(table)
end

"""
    _complement_base(base)

Return the DNA complement of a nucleotide character using a branchless lookup table.
"""
@inline function _complement_base(base::Char)
    code = Int(base) + 1
    return (1 <= code <= 128) ? _COMPLEMENT_TABLE[code] : 'N'
end

"""
    _transcript_base(base, strand)

Return the strand-aware transcript base for a genomic nucleotide.
"""
function _transcript_base(base::Char, strand::String)
    return strand == "-" ? _complement_base(base) : uppercase(base)
end

"""
    _genomic_to_cds_position(genomic_pos, sorted_exons, strand)

Map a genomic coordinate into spliced CDS space, respecting exon boundaries
and strand. Returns -1 if the position falls outside all exons (e.g., intronic).

This is the core coordinate-mapping step that Bioconductor's
`VariantAnnotation::predictCoding` performs internally.
"""
function _genomic_to_cds_position(genomic_pos::Integer, sorted_exons::AbstractVector{GffRecord}, strand::String)
    pos = Int(genomic_pos)
    if strand == "-"
        offset = 0
        for exon in Iterators.reverse(sorted_exons)
            if pos >= exon.start && pos <= exon.stop
                return offset + (Int(exon.stop) - pos + 1)
            end
            offset += Int(exon.stop) - Int(exon.start) + 1
        end
    else
        offset = 0
        for exon in sorted_exons
            if pos >= exon.start && pos <= exon.stop
                return offset + (pos - Int(exon.start) + 1)
            end
            offset += Int(exon.stop) - Int(exon.start) + 1
        end
    end
    return -1
end

"""
    _group_cds_by_gene(features)

Group CDS/coding GFF features by gene identifier so multi-exon transcripts
can be properly spliced. Returns a Dict mapping gene ID → sorted exon list.
"""
function _group_cds_by_gene(features::AbstractVector{GffRecord})
    gene_cds = Dict{String,Vector{GffRecord}}()
    for feature in features
        ft = lowercase(_feature_type(feature))
        (occursin("cds", ft) || occursin("coding", ft)) || continue
        gene_id = _feature_identifier(feature)
        isempty(gene_id) && continue
        if !haskey(gene_cds, gene_id)
            gene_cds[gene_id] = GffRecord[]
        end
        push!(gene_cds[gene_id], feature)
    end
    for (_, exons) in gene_cds
        sort!(exons; by=f -> f.start)
    end
    return gene_cds
end

"""
    _codon_effect(variant_ref, variant_alt, genomic_sequence, feature, pos; cds_exons=nothing)

Compute a codon-level consequence for a variant intersecting a coding feature.

When `cds_exons` is provided (all CDS exons for the gene, sorted by start),
the coding sequence is properly spliced from exons only — introns are excluded.
The GFF3 `phase` field is respected to set the correct reading frame.

This mirrors Bioconductor's `VariantAnnotation::predictCoding` logic.
"""
function _codon_effect(variant_ref::BioSequence{DNAAlphabet}, variant_alt::BioSequence{DNAAlphabet}, genomic_sequence::BioSequence{DNAAlphabet}, feature::GffRecord, pos::Integer; cds_exons::Union{Nothing,AbstractVector{GffRecord}}=nothing)
    strand_int = feature.strand == "-" ? Int8(-1) : Int8(1)

    # Build the coding sequence: spliced from all exons if available,
    # otherwise fall back to single-exon extraction
    if cds_exons !== nothing && length(cds_exons) > 1
        sorted_exons = sort(cds_exons; by=e -> e.start)
        transcription_order_exons = feature.strand == "-" ? Iterators.reverse(sorted_exons) : sorted_exons
        parts = AbstractFeatureLocation[FeatureLocationLite(Int(e.start), Int(e.stop)) for e in transcription_order_exons]
        compound = CompoundFeatureLocation("join", parts; strand=strand_int)
        coding_sequence = feature_sequence(genomic_sequence, compound)
        cdna_position = _genomic_to_cds_position(pos, sorted_exons, feature.strand)
    else
        coding_sequence = feature_sequence(genomic_sequence, FeatureLocationLite(feature.start, feature.stop; strand=strand_int))
        cdna_position = feature.strand == "-" ? feature.stop - Int(pos) + 1 : Int(pos) - feature.start + 1
    end

    cdna_position < 1 && return (effect = "Coding", codon_ref = "", codon_alt = "")

    # Account for GFF3 phase: number of bases to skip to reach the first
    # complete codon in the first CDS exon of this gene
    phase_offset = 0
    if cds_exons !== nothing && !isempty(cds_exons)
        first_exon = feature.strand == "-" ? cds_exons[end] : cds_exons[1]
        phase_offset = first_exon.phase === missing ? 0 : Int(first_exon.phase)
    else
        phase_offset = feature.phase === missing ? 0 : Int(feature.phase)
    end
    cdna_position -= phase_offset
    cdna_position < 1 && return (effect = "Coding", codon_ref = "", codon_alt = "")

    codon_start = 3 * div(cdna_position - 1, 3) + 1
    codon_start + 2 > length(coding_sequence) && return (effect = "Coding", codon_ref = "", codon_alt = "")

    ref_codon = coding_sequence[codon_start:codon_start+2]
    codon_bytes = copy(ref_codon.data)
    codon_index = cdna_position - codon_start + 1
    (length(variant_ref) != length(variant_alt)) && return (effect = "Indel", codon_ref = "", codon_alt = "")
    length(variant_alt) == 0 && return (effect = "Coding", codon_ref = String(ref_codon), codon_alt = String(ref_codon))
    for k in 0:length(variant_alt)-1
        target_cdna = feature.strand == "-" ? cdna_position - k : cdna_position + k
        target_index = target_cdna - codon_start + 1
        if 1 <= target_index <= 3
            alt_base = _transcript_base(Char(variant_alt.data[k+1]), feature.strand)
            codon_bytes[target_index] = UInt8(alt_base)
        end
    end
    alt_codon = BioSequence{DNAAlphabet}(codon_bytes; validate=false)
    ref_codon_text = String(ref_codon)
    alt_codon_text = String(alt_codon)

    ref_aa = translate_dna(ref_codon; stop_at_stop=true)
    alt_aa = translate_dna(alt_codon; stop_at_stop=true)
    if ref_aa != alt_aa
        if String(alt_aa) == "*"
            return (effect = "Stop-Gain", codon_ref = ref_codon_text, codon_alt = alt_codon_text)
        elseif String(ref_aa) == "*"
            return (effect = "Stop-Loss", codon_ref = ref_codon_text, codon_alt = alt_codon_text)
        else
            return (effect = "Missense", codon_ref = ref_codon_text, codon_alt = alt_codon_text)
        end
    end

    return (effect = "Synonymous", codon_ref = ref_codon_text, codon_alt = alt_codon_text)
end

"""
    _feature_consequence(variant, feature; reference_sequences=nothing, cds_exons=nothing)

Compute the predicted consequence of a variant for a single annotated feature.
When `cds_exons` is provided, multi-exon CDS is properly spliced.
"""
function _feature_consequence(variant, feature::GffRecord; reference_sequences=nothing, cds_exons::Union{Nothing,AbstractVector{GffRecord}}=nothing)
    chrom = String(_variant_field(variant, :chrom))
    position = Int(_variant_field(variant, :pos))
    variant_ref = BioSequence{DNAAlphabet}(uppercase(String(_variant_field(variant, :ref))))
    variant_alt = BioSequence{DNAAlphabet}(uppercase(String(_variant_field(variant, :alt))))

    (position < _feature_start(feature) || position > _feature_stop(feature)) && return nothing

    feature_type = lowercase(_feature_type(feature))
    gene = _feature_identifier(feature)

    if occursin("utr", feature_type)
        return (gene=gene, feature_type=_feature_type(feature), consequence="UTR", codon_ref="", codon_alt="")
    elseif occursin("intron", feature_type) || occursin("splice", feature_type)
        return (gene=gene, feature_type=_feature_type(feature), consequence="Intron", codon_ref="", codon_alt="")
    elseif occursin("cds", feature_type) || occursin("coding", feature_type)
        if reference_sequences !== nothing && haskey(reference_sequences, chrom)
            sequence_entry = reference_sequences[chrom]
            sequence = sequence_entry isa BioSequence ? BioSequence{DNAAlphabet}(sequence_entry.data; validate=false) : BioSequence{DNAAlphabet}(String(sequence_entry))
            codon_effect = _codon_effect(variant_ref, variant_alt, sequence, feature, position; cds_exons=cds_exons)
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
function annotate_variants(variant_records::AbstractVector, gene_features::AbstractVector; reference_sequences=nothing, prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    annotations = NamedTuple[]

    # Build per-chromosome interval trees for O(log F + k) lookups
    chrom_trees = Dict{String,IntervalTree{Int}}()
    gff_features = GffRecord[]
    for (idx, feature) in enumerate(gene_features)
        feature isa GffRecord || continue
        push!(gff_features, feature)
        c = _feature_chrom(feature)
        if !haskey(chrom_trees, c)
            chrom_trees[c] = IntervalTree{Int}()
        end
        insert!(chrom_trees[c], Int(_feature_start(feature)), Int(_feature_stop(feature)), length(gff_features))
    end

    # Group CDS features by gene for proper multi-exon splicing
    gene_cds_map = _group_cds_by_gene(gff_features)

    for variant in variant_records
        best = (gene="", feature_type="", consequence="Intergenic", codon_ref="", codon_alt="")
        best_rank = 0
        chrom = String(_variant_field(variant, :chrom))
        pos = Int(_variant_field(variant, :pos))

        tree = get(chrom_trees, chrom, nothing)
        if tree !== nothing
            overlapping_indices = query_overlaps(tree, pos, pos)
            for feat_idx in overlapping_indices
                feature = gff_features[feat_idx]
                # Look up all CDS exons for this gene for proper splicing
                gene_id = _feature_identifier(feature)
                exons = get(gene_cds_map, gene_id, nothing)
                consequence = _feature_consequence(variant, feature; reference_sequences=reference_sequences, cds_exons=exons)
                consequence === nothing && continue
                rank = _variant_effect_rank(consequence.consequence)
                if rank > best_rank
                    best = consequence
                    best_rank = rank
                end
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

    return _register_annotation_result!(_ctx, annotations, "annotate_variants"; parents=provenance_parent_ids(variant_records, gene_features), parameters=(variant_count=length(variant_records), feature_count=length(gene_features), annotation_count=length(annotations)))
end

# ==============================================================================
# AnnotatedSeqIndex — Interval-tree-backed feature index
#
# Equivalent to Bioconductor's GRanges + IRanges findOverlaps / subsetByOverlaps,
# but with O(log n + k) queries using the IntervalTree from biotypes.jl.
# ==============================================================================

"""
    AnnotatedSeqIndex

Interval-tree-backed index over feature positions for O(log n + k) overlap,
containment, nearest, precede, and follow queries. Built lazily from an
`AnnotatedSeqRecord` and cached for reuse.

Equivalent to the indexed overlap operations in Bioconductor's
`GenomicRanges::findOverlaps`, `IRanges::nearest`, etc.
"""
struct AnnotatedSeqIndex
    tree::IntervalTree{Int}
    features::Vector{SeqFeatureLite}
    sorted_starts::Vector{Tuple{Int,Int}}  # (start, feature_index) sorted by start
    sorted_ends::Vector{Tuple{Int,Int}}    # (end, feature_index) sorted by end
end

"""
    build_feature_index(record)

Build an `AnnotatedSeqIndex` from an annotated record for fast spatial queries.
Maintains both start-sorted and end-sorted arrays for O(log N) nearest queries.
"""
function build_feature_index(record::AnnotatedSeqRecord)
    tree = IntervalTree{Int}()
    sorted_starts = Tuple{Int,Int}[]
    sorted_ends = Tuple{Int,Int}[]
    for (idx, feature) in enumerate(record.features)
        s, e = feature_bounds(feature.location)
        lo, hi = min(s, e), max(s, e)
        insert!(tree, lo, hi, idx)
        push!(sorted_starts, (lo, idx))
        push!(sorted_ends, (hi, idx))
    end
    sort!(sorted_starts; by=first)
    sort!(sorted_ends; by=first)
    return AnnotatedSeqIndex(tree, record.features, sorted_starts, sorted_ends)
end

"""
    indexed_features_overlapping(index, query_start, query_stop)

Return features overlapping the interval [query_start, query_stop]
in O(log n + k) time using the interval-tree index.
"""
function indexed_features_overlapping(index::AnnotatedSeqIndex, query_start::Integer, query_stop::Integer)
    indices = query_overlaps(index.tree, Int(query_start), Int(query_stop))
    return SeqFeatureLite[index.features[i] for i in indices]
end

"""
    indexed_features_overlapping(index, location)

Return features overlapping a feature location in O(log n + k) time.
"""
function indexed_features_overlapping(index::AnnotatedSeqIndex, location::AbstractFeatureLocation)
    s, e = feature_bounds(location)
    return indexed_features_overlapping(index, min(s, e), max(s, e))
end

"""
    indexed_features_at(index, position)

Return features containing a single position in O(log n + k) time.
"""
indexed_features_at(index::AnnotatedSeqIndex, position::Integer) = indexed_features_overlapping(index, Int(position), Int(position))

# ==============================================================================
# Strand-aware overlap queries (Bioconductor parity)
# ==============================================================================

"""
    feature_overlaps_stranded(left, right)

Test whether two feature locations overlap **and** share the same strand.
Equivalent to Bioconductor's `findOverlaps(..., ignore.strand=FALSE)`.
"""
function feature_overlaps_stranded(left::AbstractFeatureLocation, right::AbstractFeatureLocation)
    feature_strand(left) == feature_strand(right) || return false
    return feature_overlaps(left, right)
end

feature_overlaps_stranded(left::SeqFeatureLite, right::SeqFeatureLite) = feature_overlaps_stranded(left.location, right.location)
feature_overlaps_stranded(left::SeqFeatureLite, right::AbstractFeatureLocation) = feature_overlaps_stranded(left.location, right)
feature_overlaps_stranded(left::AbstractFeatureLocation, right::SeqFeatureLite) = feature_overlaps_stranded(left, right.location)

"""
    select_features_stranded(record; kwargs...)

Like `select_features` but only returns features on the specified strand
that overlap the query region. Strand-awareness is the default in Bioconductor.
"""
function select_features_stranded(
    record::AnnotatedSeqRecord;
    feature_type=nothing,
    region=nothing,
    strand=nothing,
    qualifier_key=nothing,
    qualifier_value=nothing,
    prov_ctx=nothing,
)
    return select_features(record; feature_type=feature_type, region=region,
        strand=strand, qualifier_key=qualifier_key, qualifier_value=qualifier_value, prov_ctx=prov_ctx)
end

# ==============================================================================
# Nearest / Precede / Follow operations (Bioconductor parity)
# ==============================================================================

"""
    feature_distance(location, position)

Return the minimum distance from a position to a feature location.
Returns 0 if the position is contained within the feature.
"""
function feature_distance(location::FeatureLocationLite, position::Integer)
    lo = min(location.start, location.stop)
    hi = max(location.start, location.stop)
    pos = Int(position)
    pos < lo && return lo - pos
    pos > hi && return pos - hi
    return 0
end

function feature_distance(location::CompoundFeatureLocation, position::Integer)
    d = typemax(Int)
    for part in location.parts
        d = min(d, feature_distance(part, Int(position)))
    end
    return d
end

feature_distance(feature::SeqFeatureLite, position::Integer) = feature_distance(feature.location, position)

"""
    nearest_feature(record, position; feature_type=nothing)

Return the feature nearest to a genomic position. Ties are broken by feature order.
Equivalent to Bioconductor's `IRanges::nearest`.
"""
function nearest_feature(record::AnnotatedSeqRecord, position::Integer; feature_type=nothing)
    best = nothing
    best_dist = typemax(Int)
    for feature in record.features
        if feature_type !== nothing && feature.feature_type != String(feature_type)
            continue
        end
        d = feature_distance(feature, position)
        if d < best_dist
            best = feature
            best_dist = d
        end
    end
    return best
end

"""
    nearest_feature(index, position; n=1)

Return the n nearest features to a position using O(log N) binary search
and interval-tree containment queries.

Uses `sorted_starts` for right-side candidates and `sorted_ends` for left-side
candidates, with early termination once remaining features cannot be closer
than the current n-th best distance.
"""
function nearest_feature(index::AnnotatedSeqIndex, position::Integer; n::Int=1)
    pos = Int(position)
    isempty(index.sorted_starts) && return SeqFeatureLite[]

    # O(log N + k): find all features containing pos (distance = 0)
    containing = query_overlaps(index.tree, pos, pos)
    if length(containing) >= n
        return SeqFeatureLite[index.features[i] for i in containing[1:n]]
    end

    seen = Set{Int}(containing)
    candidates = Tuple{Int,Int}[]  # (distance, feature_index)
    for idx in containing
        push!(candidates, (0, idx))
    end

    # Binary search in sorted_starts for right-side nearest: O(log N)
    # Features with start > pos → distance = start - pos
    right_ptr = searchsortedfirst(index.sorted_starts, (pos + 1, 0); by=first)

    # Binary search in sorted_ends for left-side nearest: O(log N)
    # Features with end < pos (and not containing pos) → distance = pos - end
    left_ptr = searchsortedlast(index.sorted_ends, (pos - 1, typemax(Int)); by=first)

    # Two-pointer expansion with early termination
    while left_ptr >= 1 || right_ptr <= length(index.sorted_starts)
        # Lower bounds on distance for next unchecked features
        right_lower = right_ptr <= length(index.sorted_starts) ? (index.sorted_starts[right_ptr][1] - pos) : typemax(Int)
        left_lower = left_ptr >= 1 ? (pos - index.sorted_ends[left_ptr][1]) : typemax(Int)

        # Early termination: if we have enough and remaining can't be closer
        if length(candidates) >= n
            nth_best = candidates[n][1]
            min(left_lower, right_lower) >= nth_best && break
        end

        # Expand toward the closer side
        if right_lower <= left_lower
            _, feat_idx = index.sorted_starts[right_ptr]
            if feat_idx ∉ seen
                dist = feature_distance(index.features[feat_idx], pos)
                ins_idx = searchsortedfirst(candidates, (dist, 0); by=first)
                insert!(candidates, ins_idx, (dist, feat_idx))
                push!(seen, feat_idx)
            end
            right_ptr += 1
        else
            _, feat_idx = index.sorted_ends[left_ptr]
            if feat_idx ∉ seen
                dist = feature_distance(index.features[feat_idx], pos)
                ins_idx = searchsortedfirst(candidates, (dist, 0); by=first)
                insert!(candidates, ins_idx, (dist, feat_idx))
                push!(seen, feat_idx)
            end
            left_ptr -= 1
        end
    end

    sort!(candidates; by=first)
    result_count = min(n, length(candidates))
    return SeqFeatureLite[index.features[candidates[i][2]] for i in 1:result_count]
end

"""
    preceding_features(record, position; feature_type=nothing)

Return features that end before the given position, sorted by proximity (nearest first).
Equivalent to Bioconductor's `IRanges::precede`.
"""
function preceding_features(record::AnnotatedSeqRecord, position::Integer; feature_type=nothing)
    pos = Int(position)
    candidates = Tuple{Int,SeqFeatureLite}[]
    for feature in record.features
        if feature_type !== nothing && feature.feature_type != String(feature_type)
            continue
        end
        s, e = feature_bounds(feature.location)
        hi = max(s, e)
        if hi < pos
            push!(candidates, (pos - hi, feature))
        end
    end
    sort!(candidates; by=first)
    return SeqFeatureLite[c[2] for c in candidates]
end

"""
    following_features(record, position; feature_type=nothing)

Return features that start after the given position, sorted by proximity (nearest first).
Equivalent to Bioconductor's `IRanges::follow`.
"""
function following_features(record::AnnotatedSeqRecord, position::Integer; feature_type=nothing)
    pos = Int(position)
    candidates = Tuple{Int,SeqFeatureLite}[]
    for feature in record.features
        if feature_type !== nothing && feature.feature_type != String(feature_type)
            continue
        end
        lo = min(feature_bounds(feature.location)...)
        if lo > pos
            push!(candidates, (lo - pos, feature))
        end
    end
    sort!(candidates; by=first)
    return SeqFeatureLite[c[2] for c in candidates]
end

"""
    feature_table(record; as_dataframe=false)

Convert annotated features into a summary table. When `as_dataframe=true`,
returns a `DataFrame` instead of a vector of named tuples.
"""
function feature_table(record::AnnotatedSeqRecord; prov_ctx=nothing, as_dataframe::Bool=false)
    result = [feature_summary(feature) for feature in record.features]
    if as_dataframe
        result = DataFrames.DataFrame(record.features)
    end
    _ctx = active_provenance_context(prov_ctx)

    return _register_annotation_result!(_ctx, result, "feature_table"; parents=provenance_parent_ids(record), parameters=(row_count=length(record.features), as_dataframe=as_dataframe))
end

"""
    feature_coverage(record; resolution=1)

Compute per-position feature coverage across the record sequence.
Returns a vector of Int counts. Useful for coverage plots and
finding regions with no annotation (gaps).
"""
function feature_coverage(record::AnnotatedSeqRecord; resolution::Int=1)
    seq_len = length(record.sequence)
    seq_len == 0 && return Int[]
    coverage_len = cld(seq_len, resolution)
    delta = zeros(Int, coverage_len + 1)
    for feature in record.features
        s, e = feature_bounds(feature.location)
        lo = cld(max(1, min(s, e)), resolution)
        hi = min(cld(min(seq_len, max(s, e)), resolution), coverage_len)
        lo > hi && continue
        delta[lo] += 1
        delta[hi + 1] -= 1
    end
    return cumsum(delta[1:coverage_len])
end

# ==============================================================================
# Interactive HTML5/Canvas Visualizers for Annotation Module
# ==============================================================================

function _features_to_json(features::Vector{SeqFeatureLite})
    json_items = String[]
    for (i, f) in enumerate(features)
        s, e = feature_bounds(f.location)
        st = feature_strand(f.location)
        id_str = _json_escape(_escape_html(f.id == "" ? "$(f.feature_type)_$i" : f.id))
        type_str = _json_escape(_escape_html(f.feature_type))
        qual_pairs = String[]
        for (k, v) in f.qualifiers
            push!(qual_pairs, _json_escape(_escape_html(k)) * "=" * _json_escape(_escape_html(join(v, ","))))
        end
        q_str = join(qual_pairs, "; ")
        push!(json_items, "{\"idx\":$i,\"type\":\"$type_str\",\"id\":\"$id_str\",\"start\":$s,\"stop\":$e,\"strand\":$st,\"qualifiers\":\"$q_str\"}")
    end
    return "[" * join(json_items, ",") * "]"
end

"""
    visualize_annotation_record_html(record::AnnotatedSeqRecord; title="Sequence Annotation Map") -> String

Generate an interactive HTML5/Canvas visualization for an annotated sequence record (GenBank/GFF features).
"""
function visualize_annotation_record_html(record::AnnotatedSeqRecord; title::String="Sequence Annotation Map")
    seq_len = length(record.sequence)
    feat_count = length(record.features)
    rec_id = _escape_html(record.identifier)
    rec_name = _escape_html(record.name)
    rec_desc = _escape_html(record.description)
    feats_json = _features_to_json(record.features)

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(rec_id) - Sequence Annotation Map</title>
    <style>
        :root {
            --bg: #0f172a; --panel: #1e293b; --border: #334155;
            --text: #f8fafc; --muted: #94a3b8; --accent: #38bdf8;
            --gene: #10b981; --cds: #059669; --mrna: #3b82f6;
            --promoter: #f59e0b; --rrna: #8b5cf6; --repeat: #f43f5e;
        }
        body { font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 24px; }
        .header { background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 20px 24px; margin-bottom: 20px; box-shadow: 0 4px 12px rgba(0,0,0,0.3); }
        .title { font-size: 1.5rem; font-weight: 700; color: var(--accent); margin: 0 0 8px 0; display: flex; align-items: center; gap: 10px; }
        .meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 16px; margin-top: 14px; }
        .meta-item { background: rgba(15,23,42,0.6); padding: 10px 14px; border-radius: 8px; border: 1px solid var(--border); }
        .meta-label { font-size: 0.75rem; text-transform: uppercase; color: var(--muted); letter-spacing: 0.5px; }
        .meta-val { font-size: 1.1rem; font-weight: 600; color: var(--text); margin-top: 2px; }
        .card { background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 20px; margin-bottom: 20px; }
        .controls { display: flex; flex-wrap: wrap; gap: 12px; align-items: center; margin-bottom: 16px; justify-content: space-between; }
        .search-box { background: #0f172a; border: 1px solid var(--border); color: #fff; padding: 8px 14px; border-radius: 6px; font-size: 0.9rem; min-width: 240px; }
        .btn { background: #3b82f6; color: #fff; border: none; padding: 8px 16px; border-radius: 6px; font-weight: 600; cursor: pointer; transition: 0.2s; }
        .btn:hover { background: #2563eb; }
        canvas { width: 100%; height: 380px; display: block; border-radius: 8px; background: #0b1329; }
        .table-container { max-height: 350px; overflow-y: auto; border: 1px solid var(--border); border-radius: 8px; }
        table { width: 100%; border-collapse: collapse; text-align: left; font-size: 0.88rem; }
        th { background: #0f172a; color: var(--accent); padding: 10px 14px; position: sticky; top: 0; }
        td { padding: 9px 14px; border-bottom: 1px solid var(--border); color: #cbd5e1; }
        tr:hover td { background: rgba(56, 189, 248, 0.08); }
        .badge { display: inline-block; padding: 3px 8px; border-radius: 4px; font-size: 0.75rem; font-weight: 600; text-transform: uppercase; }
    </style>
</head>
<body>
    <div class="header">
        <div class="title">🧬 $(_escape_html(title))</div>
        <div style="color: var(--muted); font-size: 0.95rem;">$(rec_id) — $(rec_desc)</div>
        <div class="meta-grid">
            <div class="meta-item"><div class="meta-label">Length</div><div class="meta-val">$(seq_len) bp</div></div>
            <div class="meta-item"><div class="meta-label">Features</div><div class="meta-val">$(feat_count)</div></div>
            <div class="meta-item"><div class="meta-label">Name</div><div class="meta-val">$(rec_name)</div></div>
        </div>
    </div>

    <div class="card">
        <div class="controls">
            <input type="text" id="searchInput" class="search-box" placeholder="🔍 Search features by ID or type..." oninput="filterFeatures()">
            <div style="display:flex; align-items:center; gap:10px;">
                <label style="font-size:0.85rem; color:var(--muted)">Zoom:</label>
                <input type="range" id="zoomRange" min="1" max="50" value="1" oninput="drawTrack()" style="width:140px">
                <button class="btn" onclick="resetZoom()">Reset View</button>
            </div>
        </div>
        <canvas id="mapCanvas"></canvas>
    </div>

    <div class="card">
        <h3 style="margin-top:0; color:var(--accent)">📋 Feature Table</h3>
        <div class="table-container">
            <table id="featureTable">
                <thead>
                    <tr><th>#</th><th>Type</th><th>ID</th><th>Start</th><th>Stop</th><th>Span</th><th>Strand</th><th>Qualifiers</th></tr>
                </thead>
                <tbody id="tableBody"></tbody>
            </table>
        </div>
    </div>

    <script>
        const features = $(feats_json);
        const seqLength = $(seq_len);
        const canvas = document.getElementById('mapCanvas');
        const ctx = canvas.getContext('2d');

        const colorMap = {
            'gene': '#10b981', 'cds': '#059669', 'mrna': '#3b82f6',
            'exon': '#60a5fa', 'promoter': '#f59e0b', 'rrna': '#8b5cf6',
            'trna': '#a855f7', 'repeat_region': '#f43f5e', 'misc_feature': '#06b6d4'
        };

        function getColor(type) {
            const t = type.toLowerCase();
            return colorMap[t] || '#0284c7';
        }

        function drawTrack() {
            canvas.width = canvas.parentElement.clientWidth * window.devicePixelRatio;
            canvas.height = 380 * window.devicePixelRatio;
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

            const w = canvas.parentElement.clientWidth;
            const h = 380;
            ctx.clearRect(0, 0, w, h);

            const zoom = parseFloat(document.getElementById('zoomRange').value);
            const viewWidth = seqLength / zoom;

            // Axis ruler
            ctx.strokeStyle = '#334155';
            ctx.lineWidth = 1;
            ctx.beginPath();
            ctx.moveTo(50, 40); ctx.lineTo(w - 50, 40);
            ctx.stroke();

            // Ticks
            ctx.fillStyle = '#94a3b8';
            ctx.font = '11px system-ui';
            const tickCount = 8;
            for (let i = 0; i <= tickCount; i++) {
                const pos = Math.round((i / tickCount) * viewWidth);
                const x = 50 + (i / tickCount) * (w - 100);
                ctx.beginPath(); ctx.moveTo(x, 35); ctx.lineTo(x, 45); ctx.stroke();
                ctx.fillText(pos + ' bp', x - 15, 28);
            }

            // Assign lanes
            const lanes = [];
            const filterTerm = document.getElementById('searchInput').value.toLowerCase();

            features.forEach(f => {
                if (filterTerm && !f.id.toLowerCase().includes(filterTerm) && !f.type.toLowerCase().includes(filterTerm)) {
                    return;
                }
                const x1 = 50 + (f.start / viewWidth) * (w - 100);
                const x2 = 50 + (f.stop / viewWidth) * (w - 100);
                const fw = Math.max(x2 - x1, 6);

                let laneIdx = 0;
                while (lanes[laneIdx] && lanes[laneIdx] > x1 - 10) {
                    laneIdx++;
                }
                lanes[laneIdx] = x1 + fw;

                const y = 80 + laneIdx * 34;
                if (y > h - 40) return;

                // Draw feature arrow/rect
                ctx.fillStyle = getColor(f.type);
                if (f.strand === 1) {
                    ctx.beginPath();
                    ctx.moveTo(x1, y);
                    ctx.lineTo(x1 + Math.max(fw - 8, 0), y);
                    ctx.lineTo(x1 + fw, y + 10);
                    ctx.lineTo(x1 + Math.max(fw - 8, 0), y + 20);
                    ctx.lineTo(x1, y + 20);
                    ctx.closePath();
                    ctx.fill();
                } else if (f.strand === -1) {
                    ctx.beginPath();
                    ctx.moveTo(x1 + fw, y);
                    ctx.lineTo(x1 + Math.min(8, fw), y);
                    ctx.lineTo(x1, y + 10);
                    ctx.lineTo(x1 + Math.min(8, fw), y + 20);
                    ctx.lineTo(x1 + fw, y + 20);
                    ctx.closePath();
                    ctx.fill();
                } else {
                    ctx.fillRect(x1, y, fw, 20);
                }

                // Label
                ctx.fillStyle = '#f8fafc';
                ctx.font = '11px system-ui';
                if (fw > 30) {
                    ctx.fillText(f.id, x1 + 6, y + 14);
                }
            });
        }

        function populateTable() {
            const tbody = document.getElementById('tableBody');
            tbody.innerHTML = '';
            features.forEach((f, idx) => {
                const tr = document.createElement('tr');
                const strandSymbol = f.strand === 1 ? '+' : f.strand === -1 ? '-' : '.';
                tr.innerHTML = `
                    <td>\${idx + 1}</td>
                    <td><span class="badge" style="background:\${getColor(f.type)}; color:#fff">\${f.type}</span></td>
                    <td><strong>\${f.id}</strong></td>
                    <td>\${f.start}</td>
                    <td>\${f.stop}</td>
                    <td>\${f.stop - f.start + 1} bp</td>
                    <td>\${strandSymbol}</td>
                    <td><small>\${f.qualifiers || '-'}</small></td>
                `;
                tbody.appendChild(tr);
            });
        }

        function filterFeatures() {
            drawTrack();
            const term = document.getElementById('searchInput').value.toLowerCase();
            const rows = document.querySelectorAll('#tableBody tr');
            rows.forEach(row => {
                const text = row.textContent.toLowerCase();
                row.style.display = text.includes(term) ? '' : 'none';
            });
        }

        function resetZoom() {
            document.getElementById('zoomRange').value = 1;
            document.getElementById('searchInput').value = '';
            filterFeatures();
        }

        window.addEventListener('resize', drawTrack);
        populateTable();
        setTimeout(drawTrack, 50);
    </script>
</body>
</html>
"""
end

function to_html(record::AnnotatedSeqRecord)
    return visualize_annotation_record_html(record)
end

function _sanger_trace_to_json(trace::SangerTrace)
    seq_str = _json_escape(String(trace.sequence))
    quals = Int.(trace.qualities)
    ta = Int.(trace.trace_a)
    tc = Int.(trace.trace_c)
    tg = Int.(trace.trace_g)
    tt = Int.(trace.trace_t)
    return "{\"seq\":\"$seq_str\",\"quals\":[$(join(quals, ","))],\"ta\":[$(join(ta, ","))],\"tc\":[$(join(tc, ","))],\"tg\":[$(join(tg, ","))],\"tt\":[$(join(tt, ","))]}"
end

"""
    visualize_sanger_trace_html(trace::SangerTrace; title="Sanger Chromatogram Trace") -> String

Generate an interactive HTML5/Canvas 4-channel chromatogram trace viewer with base calls and Phred quality scores.
"""
function visualize_sanger_trace_html(trace::SangerTrace; title::String="Sanger Chromatogram Trace")
    total_bases = length(trace.sequence)
    mean_q = isempty(trace.qualities) ? 0.0 : round(sum(trace.qualities) / length(trace.qualities), digits=1)
    q20_pct = isempty(trace.qualities) ? 0.0 : round(count(q -> q >= 20, trace.qualities) / length(trace.qualities) * 100, digits=1)
    q30_pct = isempty(trace.qualities) ? 0.0 : round(count(q -> q >= 30, trace.qualities) / length(trace.qualities) * 100, digits=1)
    trace_json = _sanger_trace_to_json(trace)

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(_escape_html(title))</title>
    <style>
        :root {
            --bg: #0f172a; --panel: #1e293b; --border: #334155;
            --text: #f8fafc; --muted: #94a3b8; --accent: #38bdf8;
            --a: #10b981; --c: #3b82f6; --g: #f59e0b; --t: #ef4444;
        }
        body { font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 24px; }
        .card { background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 20px; margin-bottom: 20px; }
        .title { font-size: 1.5rem; font-weight: 700; color: var(--accent); margin: 0 0 12px 0; }
        .meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 14px; margin-bottom: 20px; }
        .meta-item { background: rgba(15,23,42,0.6); padding: 10px 14px; border-radius: 8px; border: 1px solid var(--border); }
        .meta-label { font-size: 0.75rem; text-transform: uppercase; color: var(--muted); }
        .meta-val { font-size: 1.1rem; font-weight: 600; color: var(--text); margin-top: 2px; }
        canvas { width: 100%; height: 420px; display: block; border-radius: 8px; background: #0b1329; }
        .controls { display: flex; gap: 14px; align-items: center; margin-bottom: 12px; }
        .legend { display: flex; gap: 18px; font-weight: 600; font-size: 0.9rem; margin-top: 10px; }
        .leg-item { display: flex; align-items: center; gap: 6px; }
        .dot { width: 12px; height: 12px; border-radius: 50%; display: inline-block; }
    </style>
</head>
<body>
    <div class="card">
        <div class="title">📈 $(_escape_html(title))</div>
        <div class="meta-grid">
            <div class="meta-item"><div class="meta-label">Base Calls</div><div class="meta-val">$(total_bases) bp</div></div>
            <div class="meta-item"><div class="meta-label">Mean Quality</div><div class="meta-val">Q$(mean_q)</div></div>
            <div class="meta-item"><div class="meta-label">% Q ≥ 20</div><div class="meta-val">$(q20_pct)%</div></div>
            <div class="meta-item"><div class="meta-label">% Q ≥ 30</div><div class="meta-val">$(q30_pct)%</div></div>
        </div>

        <div class="controls">
            <label style="font-size:0.88rem; color:var(--muted)">Scroll Position:</label>
            <input type="range" id="posSlider" min="0" max="100" value="0" style="flex:1;" oninput="drawTrace()">
        </div>

        <canvas id="traceCanvas"></canvas>

        <div class="legend">
            <div class="leg-item"><span class="dot" style="background:var(--a)"></span> A</div>
            <div class="leg-item"><span class="dot" style="background:var(--c)"></span> C</div>
            <div class="leg-item"><span class="dot" style="background:var(--g)"></span> G</div>
            <div class="leg-item"><span class="dot" style="background:var(--t)"></span> T</div>
        </div>
    </div>

    <script>
        const trace = $(trace_json);
        const canvas = document.getElementById('traceCanvas');
        const ctx = canvas.getContext('2d');

        function drawTrace() {
            canvas.width = canvas.parentElement.clientWidth * window.devicePixelRatio;
            canvas.height = 420 * window.devicePixelRatio;
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

            const w = canvas.parentElement.clientWidth;
            const h = 420;
            ctx.clearRect(0, 0, w, h);

            const seqLen = trace.seq.length;
            if (seqLen === 0) return;

            const windowSize = Math.min(80, seqLen);
            const sliderVal = parseInt(document.getElementById('posSlider').value);
            const startIdx = Math.floor((sliderVal / 100) * Math.max(0, seqLen - windowSize));
            const endIdx = Math.min(startIdx + windowSize, seqLen);

            const step = (w - 60) / (endIdx - startIdx);

            // Draw Base Calls & Quality Bars
            for (let i = startIdx; i < endIdx; i++) {
                const x = 30 + (i - startIdx) * step + step / 2;
                const base = trace.seq[i];
                const q = trace.quals[i] || 0;

                // Base text
                ctx.fillStyle = base === 'A' ? '#10b981' : base === 'C' ? '#3b82f6' : base === 'G' ? '#f59e0b' : '#ef4444';
                ctx.font = 'bold 13px system-ui';
                ctx.textAlign = 'center';
                ctx.fillText(base, x, 30);

                // Quality bar (Phred 0-40)
                const barH = (q / 40) * 45;
                ctx.fillStyle = q < 20 ? '#ef4444' : q < 30 ? '#f59e0b' : '#10b981';
                ctx.fillRect(x - 3, 75 - barH, 6, barH);
            }

            // Draw Chromatogram Signals
            const channels = [
                { data: trace.ta, color: '#10b981' },
                { data: trace.tc, color: '#3b82f6' },
                { data: trace.tg, color: '#f59e0b' },
                { data: trace.tt, color: '#ef4444' }
            ];

            const traceLen = trace.ta.length;
            if (traceLen > 0) {
                const tracePointsPerBase = traceLen / seqLen;
                const traceStart = Math.floor(startIdx * tracePointsPerBase);
                const traceEnd = Math.min(Math.floor(endIdx * tracePointsPerBase), traceLen);

                let maxVal = 1;
                for (let i = traceStart; i < traceEnd; i++) {
                    maxVal = Math.max(maxVal, trace.ta[i]||0, trace.tc[i]||0, trace.tg[i]||0, trace.tt[i]||0);
                }

                channels.forEach(ch => {
                    ctx.strokeStyle = ch.color;
                    ctx.lineWidth = 1.8;
                    ctx.beginPath();
                    for (let i = traceStart; i < traceEnd; i++) {
                        const relPos = (i - traceStart) / (traceEnd - traceStart);
                        const x = 30 + relPos * (w - 60);
                        const val = ch.data[i] || 0;
                        const y = h - 20 - (val / maxVal) * (h - 120);
                        if (i === traceStart) ctx.moveTo(x, y);
                        else ctx.lineTo(x, y);
                    }
                    ctx.stroke();
                });
            }
        }

        window.addEventListener('resize', drawTrace);
        setTimeout(drawTrace, 50);
    </script>
</body>
</html>
"""
end

function to_html(trace::SangerTrace)
    return visualize_sanger_trace_html(trace)
end

function _variant_annotations_to_json(annotations)
    items = String[]
    for a in annotations
        c = _json_escape(_escape_html(String(get(a, :chrom, ""))))
        p = Int(get(a, :pos, 0))
        r = _json_escape(_escape_html(String(get(a, :ref, ""))))
        alt = _json_escape(_escape_html(String(get(a, :alt, ""))))
        g = _json_escape(_escape_html(String(get(a, :gene, ""))))
        ft = _json_escape(_escape_html(String(get(a, :feature_type, ""))))
        cq = _json_escape(_escape_html(String(get(a, :consequence, ""))))
        cref = _json_escape(_escape_html(String(get(a, :codon_ref, ""))))
        calt = _json_escape(_escape_html(String(get(a, :codon_alt, ""))))
        push!(items, "{\"chrom\":\"$c\",\"pos\":$p,\"ref\":\"$r\",\"alt\":\"$alt\",\"gene\":\"$g\",\"feature_type\":\"$ft\",\"consequence\":\"$cq\",\"codon_ref\":\"$cref\",\"codon_alt\":\"$calt\"}")
    end
    return "[" * join(items, ",") * "]"
end

"""
    visualize_variant_consequences_html(annotations; title="Variant Consequence Annotations") -> String

Generate an interactive HTML5 dashboard visualizing variant consequences and impact breakdown.
"""
function visualize_variant_consequences_html(annotations::AbstractVector; title::String="Variant Consequence Annotations")
    total_vars = length(annotations)
    missense_cnt = count(a -> get(a, :consequence, "") == "Missense", annotations)
    stop_cnt = count(a -> get(a, :consequence, "") in ("Stop-Gain", "Stop-Loss"), annotations)
    syn_cnt = count(a -> get(a, :consequence, "") == "Synonymous", annotations)
    noncoding_cnt = total_vars - missense_cnt - stop_cnt - syn_cnt
    vars_json = _variant_annotations_to_json(annotations)

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(_escape_html(title))</title>
    <style>
        :root {
            --bg: #0f172a; --panel: #1e293b; --border: #334155;
            --text: #f8fafc; --muted: #94a3b8; --accent: #38bdf8;
            --high: #ef4444; --mod: #f59e0b; --low: #10b981;
        }
        body { font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 24px; }
        .card { background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 20px; margin-bottom: 20px; }
        .title { font-size: 1.5rem; font-weight: 700; color: var(--accent); margin: 0 0 14px 0; }
        .meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 14px; margin-bottom: 20px; }
        .meta-item { background: rgba(15,23,42,0.6); padding: 10px 14px; border-radius: 8px; border: 1px solid var(--border); }
        .meta-label { font-size: 0.75rem; text-transform: uppercase; color: var(--muted); }
        .meta-val { font-size: 1.1rem; font-weight: 600; color: var(--text); margin-top: 2px; }
        .table-container { max-height: 380px; overflow-y: auto; border: 1px solid var(--border); border-radius: 8px; }
        table { width: 100%; border-collapse: collapse; text-align: left; font-size: 0.88rem; }
        th { background: #0f172a; color: var(--accent); padding: 10px 14px; position: sticky; top: 0; }
        td { padding: 9px 14px; border-bottom: 1px solid var(--border); color: #cbd5e1; }
        tr:hover td { background: rgba(56, 189, 248, 0.08); }
        .badge { display: inline-block; padding: 3px 8px; border-radius: 4px; font-size: 0.75rem; font-weight: 600; }
        .search-box { background: #0f172a; border: 1px solid var(--border); color: #fff; padding: 8px 14px; border-radius: 6px; font-size: 0.9rem; min-width: 240px; margin-bottom: 12px; }
    </style>
</head>
<body>
    <div class="card">
        <div class="title">🎯 $(_escape_html(title))</div>
        <div class="meta-grid">
            <div class="meta-item"><div class="meta-label">Total Variants</div><div class="meta-val">$(total_vars)</div></div>
            <div class="meta-item"><div class="meta-label">Stop-Gain/Loss</div><div class="meta-val" style="color:var(--high)">$(stop_cnt)</div></div>
            <div class="meta-item"><div class="meta-label">Missense</div><div class="meta-val" style="color:var(--mod)">$(missense_cnt)</div></div>
            <div class="meta-item"><div class="meta-label">Synonymous</div><div class="meta-val" style="color:var(--low)">$(syn_cnt)</div></div>
            <div class="meta-item"><div class="meta-label">Non-Coding</div><div class="meta-val">$(noncoding_cnt)</div></div>
        </div>

        <input type="text" id="searchInput" class="search-box" placeholder="🔍 Search variants by gene or consequence..." oninput="filterTable()">

        <div class="table-container">
            <table>
                <thead>
                    <tr><th>Chr</th><th>Position</th><th>Ref</th><th>Alt</th><th>Gene</th><th>Feature</th><th>Consequence</th><th>Codon</th></tr>
                </thead>
                <tbody id="tableBody"></tbody>
            </table>
        </div>
    </div>

    <script>
        const variants = $(vars_json);

        function getBadge(cq) {
            if (cq === 'Stop-Gain' || cq === 'Stop-Loss') return '<span class="badge" style="background:#ef4444; color:#fff">' + cq + '</span>';
            if (cq === 'Missense') return '<span class="badge" style="background:#f59e0b; color:#fff">' + cq + '</span>';
            if (cq === 'Synonymous') return '<span class="badge" style="background:#10b981; color:#fff">' + cq + '</span>';
            return '<span class="badge" style="background:#64748b; color:#fff">' + cq + '</span>';
        }

        function populateTable() {
            const tbody = document.getElementById('tableBody');
            tbody.innerHTML = '';
            variants.forEach(v => {
                const tr = document.createElement('tr');
                const codonStr = v.codon_ref ? (v.codon_ref + ' → ' + v.codon_alt) : '-';
                tr.innerHTML = `
                    <td>\${v.chrom}</td>
                    <td>\${v.pos}</td>
                    <td><code>\${v.ref}</code></td>
                    <td><code>\${v.alt}</code></td>
                    <td><strong>\${v.gene || '-'}</strong></td>
                    <td>\${v.feature_type || '-'}</td>
                    <td>\${getBadge(v.consequence)}</td>
                    <td><code>\${codonStr}</code></td>
                `;
                tbody.appendChild(tr);
            });
        }

        function filterTable() {
            const term = document.getElementById('searchInput').value.toLowerCase();
            const rows = document.querySelectorAll('#tableBody tr');
            rows.forEach(row => {
                row.style.display = row.textContent.toLowerCase().includes(term) ? '' : 'none';
            });
        }

        populateTable();
    </script>
</body>
</html>
"""
end
