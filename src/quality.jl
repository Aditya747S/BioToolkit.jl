using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!

@inline function _register_quality_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

"""
    phred_score(byte; offset=33)

Convert a single ASCII-encoded quality byte into a Phred score.
"""
@inline function phred_score(byte::UInt8; offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    score = Int(byte) - Int(offset)
    score < 0 && throw(ArgumentError("quality score below offset"))
    result = UInt8(score)
    return _register_quality_result!(_ctx, result, "phred_score"; parameters=(offset=Int(offset), score=result))
end

"""
    phred_scores(quality; offset=33)

Convert a FASTQ quality string into numeric Phred scores.
"""
function phred_scores(quality::String; offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))


    return phred_scores(codeunits(quality); offset=offset, _ctx=_ctx)
end

"""
    phred_scores(bytes; offset=33)

Convert a byte vector of encoded qualities into numeric Phred scores.
"""
function phred_scores(bytes::AbstractVector{UInt8}; offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    length_scores = length(bytes)
    scores = Vector{UInt8}(undef, length_scores)

    @inbounds for index in 1:length_scores
        scores[index] = phred_score(bytes[index]; offset=offset, _ctx=_ctx)
    end


    return _register_quality_result!(_ctx, scores, "phred_scores"; parents=provenance_parent_ids(bytes), parameters=(offset=Int(offset), score_count=length(scores)))
end

"""
    phred_scores(record; offset=33)

Convert the quality string stored in a FASTQ record into Phred scores.
"""
function phred_scores(record::FastqRecord; offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))


    return phred_scores(record.quality; offset=offset, _ctx=_ctx)
end

"""
    phred_scores(record; quality_key=:quality, offset=33)

Convert a quality annotation stored on a lightweight record into Phred scores.
"""
function phred_scores(record::SeqRecordLite; quality_key::Symbol=:quality, offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    quality = get(record.letter_annotations, quality_key, nothing)
    quality === nothing && throw(ArgumentError("missing $(quality_key) letter annotation"))
    quality isa String || throw(ArgumentError("FASTQ quality annotation must be a string"))
    length(quality) == length(record.sequence) || throw(ArgumentError("quality annotation length must match record sequence length"))


    return phred_scores(quality; offset=offset, _ctx=_ctx)
end

"""
    phred_string(scores; offset=33)

Convert numeric Phred scores back into an ASCII quality string.
"""
function phred_string(scores::AbstractVector{<:Integer}; offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    length_scores = length(scores)
    bytes = Vector{UInt8}(undef, length_scores)
    base_offset = UInt8(offset)

    @inbounds for index in 1:length_scores
        score = scores[index]
        score < 0 && throw(ArgumentError("quality scores must be nonnegative"))
        bytes[index] = base_offset + UInt8(score)
    end

    result = String(bytes)


    return _register_quality_result!(_ctx, result, "phred_string"; parents=provenance_parent_ids(scores), parameters=(offset=Int(offset), length=length_scores))
end

"""
    mean_quality(scores; offset=33)

Return the average Phred score for a numeric score vector.
"""
function mean_quality(scores::AbstractVector{<:Integer}; offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    isempty(scores) && return 0.0

    total = 0
    @inbounds for score in scores
        total += Int(score)
    end

    result = total / length(scores)


    return _register_quality_result!(_ctx, result, "mean_quality"; parents=provenance_parent_ids(scores), parameters=(offset=Int(offset), score=result))
end

"""
    mean_quality(quality; offset=33)

Return the average Phred score for a quality string.
"""
function mean_quality(quality::String; offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))


    return mean_quality(phred_scores(quality; offset=offset, _ctx=_ctx); offset=offset, _ctx=_ctx)
end

"""
    mean_quality(record; offset=33)

Return the average Phred score for a FASTQ record.
"""
function mean_quality(record::FastqRecord; offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))


    return mean_quality(record.quality; offset=offset, _ctx=_ctx)
end

"""
    mean_quality(record; quality_key=:quality, offset=33)

Return the average Phred score for a lightweight record's quality annotation.
"""
function mean_quality(record::SeqRecordLite; quality_key::Symbol=:quality, offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    quality = get(record.letter_annotations, quality_key, nothing)
    quality === nothing && throw(ArgumentError("missing $(quality_key) letter annotation"))
    quality isa String || throw(ArgumentError("FASTQ quality annotation must be a string"))


    return mean_quality(quality; offset=offset, _ctx=_ctx)
end

"""
    _quality_scores(record; offset=33)

Internal helper that extracts numeric quality scores from a FASTQ record.
"""
function _quality_scores(record::FastqRecord; offset::Integer=33, _ctx=active_provenance_context())
    return phred_scores(record.quality; offset=offset, _ctx=_ctx)
end

"""
    _quality_scores(record; quality_key=:quality, offset=33)

Internal helper that extracts numeric quality scores from a lightweight record.
"""
function _quality_scores(record::SeqRecordLite; quality_key::Symbol=:quality, offset::Integer=33, _ctx=active_provenance_context())
    return phred_scores(record; quality_key=quality_key, offset=offset, _ctx=_ctx)
end

"""
    quality_filter(record; kwargs...)

Apply read-level quality filters to a FASTQ or lightweight record.
"""
function quality_filter(
    record::Union{FastqRecord,SeqRecordLite};
    min_mean_quality::Real=20,
    min_base_quality::Real=0,
    max_low_quality_fraction::Real=1.0,
    quality_key::Symbol=:quality,
    offset::Integer=33,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))
    0 <= max_low_quality_fraction <= 1 || throw(ArgumentError("max_low_quality_fraction must be between 0 and 1"))
    scores = record isa FastqRecord ? _quality_scores(record; offset=offset, _ctx=_ctx) : _quality_scores(record; quality_key=quality_key, offset=offset, _ctx=_ctx)
    isempty(scores) && return _register_quality_result!(_ctx, false, "quality_filter"; parents=provenance_parent_ids(record), parameters=(min_mean_quality=Float64(min_mean_quality), min_base_quality=Float64(min_base_quality), max_low_quality_fraction=Float64(max_low_quality_fraction), result=false))

    mean_quality(scores; _ctx=_ctx) >= min_mean_quality || return false

    low_quality_count = 0
    @inbounds for score in scores
        Int(score) < min_base_quality && (low_quality_count += 1)
    end

    result = low_quality_count / length(scores) <= max_low_quality_fraction


    return _register_quality_result!(_ctx, result, "quality_filter"; parents=provenance_parent_ids(record), parameters=(min_mean_quality=Float64(min_mean_quality), min_base_quality=Float64(min_base_quality), max_low_quality_fraction=Float64(max_low_quality_fraction), result=result))
end

"""
    _adapter_trim_span(sequence, adapter; min_overlap=8, max_mismatches=1, from_end=:three_prime)

Find the span to keep after adapter trimming.
"""
@inline function _quality_upper_ascii(byte::UInt8)
    return byte >= 0x61 && byte <= 0x7a ? byte - 0x20 : byte
end

function _adapter_trim_span_bytes(sequence_bytes::AbstractVector{UInt8}, adapter_bytes::AbstractVector{UInt8}; min_overlap::Int=8, max_mismatches::Int=1, from_end::Symbol=:three_prime)
    sequence_length = length(sequence_bytes)
    adapter_length = length(adapter_bytes)
    min_overlap > 0 || throw(ArgumentError("min_overlap must be positive"))
    max_mismatches >= 0 || throw(ArgumentError("max_mismatches must be nonnegative"))
    sequence_length >= min_overlap || return nothing
    adapter_length >= min_overlap || return nothing

    max_overlap = min(sequence_length, adapter_length)
    if from_end === :three_prime
        for overlap in max_overlap:-1:min_overlap
            mismatches = 0
            seq_start = sequence_length - overlap
            @inbounds for j in 1:overlap
                left_byte = _quality_upper_ascii(sequence_bytes[seq_start + j])
                right_byte = _quality_upper_ascii(adapter_bytes[j])
                left_byte == right_byte && continue
                mismatches += 1
                mismatches > max_mismatches && break
            end
            mismatches <= max_mismatches && return 1, sequence_length - overlap
        end
    elseif from_end === :five_prime
        for overlap in max_overlap:-1:min_overlap
            mismatches = 0
            adapter_start = adapter_length - overlap
            @inbounds for j in 1:overlap
                left_byte = _quality_upper_ascii(sequence_bytes[j])
                right_byte = _quality_upper_ascii(adapter_bytes[adapter_start + j])
                left_byte == right_byte && continue
                mismatches += 1
                mismatches > max_mismatches && break
            end
            mismatches <= max_mismatches && return overlap + 1, sequence_length
        end
    else
        throw(ArgumentError("from_end must be :three_prime or :five_prime"))
    end

    return nothing
end

_adapter_trim_span(sequence::BioSequence, adapter::BioSequence; kwargs...) = _adapter_trim_span_bytes(sequence.data, adapter.data; kwargs...)

@inline _quality_adapter_sequence(adapter::AbstractString) = DNASeq(String(adapter))

"""
    adapter_trim(record; kwargs...)

Trim adapter contamination from a FASTQ record.
"""
function adapter_trim(
    record::FastqRecord;
    adapter::Union{BioSequence,AbstractString},
    min_overlap::Int=8,
    max_mismatches::Int=1,
    from_end::Symbol=:three_prime,
    quality_key::Symbol=:quality,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))
    adapter = adapter isa AbstractString ? _quality_adapter_sequence(adapter) : adapter
    trim_span = _adapter_trim_span(record.sequence, adapter; min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end)
    trim_span === nothing && return _register_quality_result!(_ctx, record, "adapter_trim"; parents=provenance_parent_ids(record), parameters=(min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end, trimmed=false))
    start_index, stop_index = trim_span
    start_index > stop_index && return _register_quality_result!(_ctx, nothing, "adapter_trim"; parents=provenance_parent_ids(record), parameters=(min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end, trimmed=true, kept=false))
    result = FastqRecord(record.identifier, record.description, record.sequence[start_index:stop_index], record.quality[start_index:stop_index])


    return _register_quality_result!(_ctx, result, "adapter_trim"; parents=provenance_parent_ids(record), parameters=(min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end, trimmed=true))
end

"""
    adapter_trim(record; kwargs...)

Trim adapter contamination from a lightweight record.
"""
function adapter_trim(
    record::SeqRecordLite;
    adapter::Union{BioSequence,AbstractString},
    min_overlap::Int=8,
    max_mismatches::Int=1,
    from_end::Symbol=:three_prime,
    quality_key::Symbol=:quality,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))
    adapter = adapter isa AbstractString ? _quality_adapter_sequence(adapter) : adapter
    trim_span = _adapter_trim_span(record.sequence, adapter; min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end)
    trim_span === nothing && return _register_quality_result!(_ctx, record, "adapter_trim"; parents=provenance_parent_ids(record), parameters=(min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end, trimmed=false))
    start_index, stop_index = trim_span
    start_index > stop_index && return _register_quality_result!(_ctx, nothing, "adapter_trim"; parents=provenance_parent_ids(record), parameters=(min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end, trimmed=true, kept=false))

    trimmed_annotations = copy(record.annotations)
    trimmed_letter_annotations = copy(record.letter_annotations)
    if haskey(trimmed_letter_annotations, quality_key)
        quality = trimmed_letter_annotations[quality_key]
        quality isa String && (trimmed_letter_annotations[quality_key] = quality[start_index:stop_index])
    end

    result = SeqRecordLite(
        record.sequence[start_index:stop_index];
        identifier=record.identifier,
        name=record.name,
        description=record.description,
        annotations=trimmed_annotations,
        letter_annotations=trimmed_letter_annotations)


    return _register_quality_result!(_ctx, result, "adapter_trim"; parents=provenance_parent_ids(record), parameters=(min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end, trimmed=true))
end

"""
    process_sequencing_record(record; kwargs...)

Run adapter trimming, low-quality trimming, length checks, and quality filters.
"""
function process_sequencing_record(
    record::Union{FastqRecord,SeqRecordLite};
    adapter::Union{Nothing,BioSequence,AbstractString}=nothing,
    min_mean_quality::Real=20,
    min_base_quality::Real=0,
    max_low_quality_fraction::Real=1.0,
    window::Int=4,
    threshold::Real=20,
    min_length::Int=0,
    min_overlap::Int=8,
    max_mismatches::Int=1,
    from_end::Symbol=:three_prime,
    quality_key::Symbol=:quality,
    offset::Integer=33,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    processed = record
    if adapter !== nothing
        processed = adapter_trim(processed; adapter=adapter, min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end, quality_key=quality_key, _ctx=_ctx)
        processed === nothing && return nothing
    end

    processed = trim_low_quality(processed; window=window, threshold=threshold, quality_key=quality_key, offset=offset, _ctx=_ctx)
    processed === nothing && return nothing
    ncodeunits(processed.sequence) >= min_length || return nothing
    quality_filter(processed; min_mean_quality=min_mean_quality, min_base_quality=min_base_quality, max_low_quality_fraction=max_low_quality_fraction, quality_key=quality_key, offset=offset, _ctx=_ctx) || return nothing
    return _register_quality_result!(_ctx, processed, "process_sequencing_record"; parents=provenance_parent_ids(record), parameters=(min_mean_quality=Float64(min_mean_quality), min_base_quality=Float64(min_base_quality), max_low_quality_fraction=Float64(max_low_quality_fraction), min_length=min_length))
end

"""
    sequencing_pipeline(records; kwargs...)

Apply sequencing cleanup to a batch of FASTQ or lightweight records.
"""
function sequencing_pipeline(
    records::AbstractVector{T};
    adapter::Union{Nothing,BioSequence,AbstractString}=nothing,
    min_mean_quality::Real=20,
    min_base_quality::Real=0,
    max_low_quality_fraction::Real=1.0,
    window::Int=4,
    threshold::Real=20,
    min_length::Int=0,
    min_overlap::Int=8,
    max_mismatches::Int=1,
    from_end::Symbol=:three_prime,
    quality_key::Symbol=:quality,
    offset::Integer=33,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx)) where {T<:Union{FastqRecord,SeqRecordLite}}
    processed = T[]
    for record in records
        filtered = process_sequencing_record(
            record;
            adapter=adapter,
            min_mean_quality=min_mean_quality,
            min_base_quality=min_base_quality,
            max_low_quality_fraction=max_low_quality_fraction,
            window=window,
            threshold=threshold,
            min_length=min_length,
            min_overlap=min_overlap,
            max_mismatches=max_mismatches,
            from_end=from_end,
            quality_key=quality_key,
            offset=offset,
            _ctx=_ctx)
        filtered === nothing && continue
        push!(processed, filtered)
    end


    return _register_quality_result!(_ctx, processed, "sequencing_pipeline"; parents=provenance_parent_ids(records), parameters=(input_count=length(records), output_count=length(processed), min_length=min_length))
end

"""
    _trim_window(scores, window, threshold)

Identify the retained span when trimming low-quality ends by sliding window.
"""
function _trim_window(scores::AbstractVector{UInt8}, window::Int, threshold::Real)
    window > 0 || throw(ArgumentError("window must be positive"))
    length(scores) >= window || return nothing

    prefix = Vector{Int}(undef, length(scores) + 1)
    prefix[1] = 0
    @inbounds for index in eachindex(scores)
        prefix[index + 1] = prefix[index] + Int(scores[index])
    end

    first_keep = nothing
    last_keep = nothing
    threshold_sum = threshold * window

    @inbounds for start_index in 1:(length(scores) - window + 1)
        window_sum = prefix[start_index + window] - prefix[start_index]
        if window_sum >= threshold_sum
            first_keep === nothing && (first_keep = start_index)
            last_keep = start_index + window - 1
        end
    end

    first_keep === nothing && return nothing
    return Int(first_keep)::Int, Int(last_keep)::Int
end

"""
    trim_low_quality(record; kwargs...)

Trim low-quality ends from a FASTQ record.
"""
function trim_low_quality(record::FastqRecord; window::Int=4, threshold::Real=20, quality_key::Symbol=:quality, offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    scores = phred_scores(record.quality; offset=offset, _ctx=_ctx)
    trimmed = _trim_window(scores, window, threshold)
    trimmed === nothing && return _register_quality_result!(_ctx, nothing, "trim_low_quality"; parents=provenance_parent_ids(record), parameters=(window=window, threshold=Float64(threshold), trimmed=false))

    start_index, stop_index = trimmed
    result = FastqRecord(
        record.identifier,
        record.description,
        record.sequence[start_index:stop_index],
        record.quality[start_index:stop_index])


    return _register_quality_result!(_ctx, result, "trim_low_quality"; parents=provenance_parent_ids(record), parameters=(window=window, threshold=Float64(threshold), trimmed=true))
end

"""
    trim_low_quality(record; kwargs...)

Trim low-quality ends from a lightweight record.
"""
function trim_low_quality(record::SeqRecordLite; quality_key::Symbol=:quality, window::Int=4, threshold::Real=20, offset::Integer=33, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    quality = get(record.letter_annotations, quality_key, nothing)
    quality === nothing && throw(ArgumentError("missing $(quality_key) letter annotation"))
    quality isa String || throw(ArgumentError("FASTQ quality annotation must be a string"))

    scores = phred_scores(quality; offset=offset, _ctx=_ctx)
    trimmed = _trim_window(scores, window, threshold)
    trimmed === nothing && return _register_quality_result!(_ctx, nothing, "trim_low_quality"; parents=provenance_parent_ids(record), parameters=(window=window, threshold=Float64(threshold), trimmed=false))

    start_index, stop_index = trimmed
    trimmed_annotations = copy(record.annotations)
    trimmed_letter_annotations = copy(record.letter_annotations)
    trimmed_letter_annotations[quality_key] = quality[start_index:stop_index]

    result = SeqRecordLite(
        record.sequence[start_index:stop_index];
        identifier=record.identifier,
        name=record.name,
        description=record.description,
        annotations=trimmed_annotations,
        letter_annotations=trimmed_letter_annotations)


    return _register_quality_result!(_ctx, result, "trim_low_quality"; parents=provenance_parent_ids(record), parameters=(window=window, threshold=Float64(threshold), trimmed=true))
end
