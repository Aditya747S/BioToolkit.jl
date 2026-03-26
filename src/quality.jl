"""
    phred_score(byte; offset=33)

Convert a single ASCII-encoded quality byte into a Phred score.
"""
@inline function phred_score(byte::UInt8; offset::Integer=33)
    score = Int(byte) - Int(offset)
    score < 0 && throw(ArgumentError("quality score below offset"))
    return UInt8(score)
end

"""
    phred_scores(quality; offset=33)

Convert a FASTQ quality string into numeric Phred scores.
"""
function phred_scores(quality::AbstractString; offset::Integer=33)
    return phred_scores(codeunits(quality); offset=offset)
end

"""
    phred_scores(bytes; offset=33)

Convert a byte vector of encoded qualities into numeric Phred scores.
"""
function phred_scores(bytes::AbstractVector{UInt8}; offset::Integer=33)
    length_scores = length(bytes)
    scores = Vector{UInt8}(undef, length_scores)

    @inbounds for index in 1:length_scores
        scores[index] = phred_score(bytes[index]; offset=offset)
    end

    return scores
end

"""
    phred_scores(record; offset=33)

Convert the quality string stored in a FASTQ record into Phred scores.
"""
function phred_scores(record::FastqRecord; offset::Integer=33)
    return phred_scores(record.quality; offset=offset)
end

"""
    phred_scores(record; quality_key=:quality, offset=33)

Convert a quality annotation stored on a lightweight record into Phred scores.
"""
function phred_scores(record::SeqRecordLite; quality_key::Symbol=:quality, offset::Integer=33)
    quality = get(record.letter_annotations, quality_key, nothing)
    quality === nothing && throw(ArgumentError("missing $(quality_key) letter annotation"))
    quality isa AbstractString || throw(ArgumentError("FASTQ quality annotation must be a string"))
    return phred_scores(quality; offset=offset)
end

"""
    phred_string(scores; offset=33)

Convert numeric Phred scores back into an ASCII quality string.
"""
function phred_string(scores::AbstractVector{<:Integer}; offset::Integer=33)
    length_scores = length(scores)
    bytes = Vector{UInt8}(undef, length_scores)
    base_offset = UInt8(offset)

    @inbounds for index in 1:length_scores
        score = scores[index]
        score < 0 && throw(ArgumentError("quality scores must be nonnegative"))
        bytes[index] = base_offset + UInt8(score)
    end

    return String(bytes)
end

"""
    mean_quality(scores; offset=33)

Return the average Phred score for a numeric score vector.
"""
function mean_quality(scores::AbstractVector{<:Integer}; offset::Integer=33)
    isempty(scores) && return 0.0

    total = 0
    @inbounds for score in scores
        total += Int(score)
    end

    return total / length(scores)
end

"""
    mean_quality(quality; offset=33)

Return the average Phred score for a quality string.
"""
function mean_quality(quality::AbstractString; offset::Integer=33)
    return mean_quality(phred_scores(quality; offset=offset); offset=offset)
end

"""
    mean_quality(record; offset=33)

Return the average Phred score for a FASTQ record.
"""
function mean_quality(record::FastqRecord; offset::Integer=33)
    return mean_quality(record.quality; offset=offset)
end

"""
    mean_quality(record; quality_key=:quality, offset=33)

Return the average Phred score for a lightweight record's quality annotation.
"""
function mean_quality(record::SeqRecordLite; quality_key::Symbol=:quality, offset::Integer=33)
    quality = get(record.letter_annotations, quality_key, nothing)
    quality === nothing && throw(ArgumentError("missing $(quality_key) letter annotation"))
    quality isa AbstractString || throw(ArgumentError("FASTQ quality annotation must be a string"))
    return mean_quality(quality; offset=offset)
end

"""
    _quality_scores(record; offset=33)

Internal helper that extracts numeric quality scores from a FASTQ record.
"""
function _quality_scores(record::FastqRecord; offset::Integer=33)
    return phred_scores(record.quality; offset=offset)
end

"""
    _quality_scores(record; quality_key=:quality, offset=33)

Internal helper that extracts numeric quality scores from a lightweight record.
"""
function _quality_scores(record::SeqRecordLite; quality_key::Symbol=:quality, offset::Integer=33)
    return phred_scores(record; quality_key=quality_key, offset=offset)
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
)
    0 <= max_low_quality_fraction <= 1 || throw(ArgumentError("max_low_quality_fraction must be between 0 and 1"))
    scores = record isa FastqRecord ? _quality_scores(record; offset=offset) : _quality_scores(record; quality_key=quality_key, offset=offset)
    isempty(scores) && return false

    mean_quality(scores) >= min_mean_quality || return false

    low_quality_count = 0
    @inbounds for score in scores
        Int(score) < min_base_quality && (low_quality_count += 1)
    end

    return low_quality_count / length(scores) <= max_low_quality_fraction
end

"""
    _adapter_trim_span(sequence, adapter; min_overlap=8, max_mismatches=1, from_end=:three_prime)

Find the span to keep after adapter trimming.
"""
function _adapter_trim_span(sequence::AbstractString, adapter::AbstractString; min_overlap::Int=8, max_mismatches::Int=1, from_end::Symbol=:three_prime)
    sequence_length = ncodeunits(sequence)
    adapter_length = ncodeunits(adapter)
    min_overlap > 0 || throw(ArgumentError("min_overlap must be positive"))
    max_mismatches >= 0 || throw(ArgumentError("max_mismatches must be nonnegative"))
    sequence_length >= min_overlap || return nothing
    adapter_length >= min_overlap || return nothing

    max_overlap = min(sequence_length, adapter_length)
    if from_end === :three_prime
        for overlap in max_overlap:-1:min_overlap
            sequence_window = sequence[sequence_length - overlap + 1:sequence_length]
            adapter_window = adapter[1:overlap]
            mismatches = 0
            @inbounds for (left_byte, right_byte) in zip(codeunits(sequence_window), codeunits(adapter_window))
                left_byte == right_byte && continue
                mismatches += 1
                mismatches > max_mismatches && break
            end
            mismatches <= max_mismatches && return 1, sequence_length - overlap
        end
    elseif from_end === :five_prime
        for overlap in max_overlap:-1:min_overlap
            sequence_window = sequence[1:overlap]
            adapter_window = adapter[adapter_length - overlap + 1:adapter_length]
            mismatches = 0
            @inbounds for (left_byte, right_byte) in zip(codeunits(sequence_window), codeunits(adapter_window))
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

"""
    adapter_trim(record; kwargs...)

Trim adapter contamination from a FASTQ record.
"""
function adapter_trim(
    record::FastqRecord;
    adapter::AbstractString,
    min_overlap::Int=8,
    max_mismatches::Int=1,
    from_end::Symbol=:three_prime,
    quality_key::Symbol=:quality,
)
    trim_span = _adapter_trim_span(record.sequence, adapter; min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end)
    trim_span === nothing && return record
    start_index, stop_index = trim_span
    start_index > stop_index && return nothing
    return FastqRecord(record.identifier, record.description, record.sequence[start_index:stop_index], record.quality[start_index:stop_index])
end

"""
    adapter_trim(record; kwargs...)

Trim adapter contamination from a lightweight record.
"""
function adapter_trim(
    record::SeqRecordLite;
    adapter::AbstractString,
    min_overlap::Int=8,
    max_mismatches::Int=1,
    from_end::Symbol=:three_prime,
    quality_key::Symbol=:quality,
)
    trim_span = _adapter_trim_span(record.sequence, adapter; min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end)
    trim_span === nothing && return record
    start_index, stop_index = trim_span
    start_index > stop_index && return nothing

    trimmed_annotations = copy(record.annotations)
    trimmed_letter_annotations = copy(record.letter_annotations)
    if haskey(trimmed_letter_annotations, quality_key)
        quality = trimmed_letter_annotations[quality_key]
        quality isa AbstractString && (trimmed_letter_annotations[quality_key] = quality[start_index:stop_index])
    end

    return SeqRecordLite(
        record.sequence[start_index:stop_index];
        identifier=record.identifier,
        name=record.name,
        description=record.description,
        annotations=trimmed_annotations,
        letter_annotations=trimmed_letter_annotations,
    )
end

"""
    process_sequencing_record(record; kwargs...)

Run adapter trimming, low-quality trimming, length checks, and quality filters.
"""
function process_sequencing_record(
    record::Union{FastqRecord,SeqRecordLite};
    adapter::Union{Nothing,AbstractString}=nothing,
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
)
    processed = record
    if adapter !== nothing
        processed = adapter_trim(processed; adapter=adapter, min_overlap=min_overlap, max_mismatches=max_mismatches, from_end=from_end, quality_key=quality_key)
        processed === nothing && return nothing
    end

    processed = trim_low_quality(processed; window=window, threshold=threshold, quality_key=quality_key, offset=offset)
    processed === nothing && return nothing
    ncodeunits(processed.sequence) >= min_length || return nothing
    quality_filter(processed; min_mean_quality=min_mean_quality, min_base_quality=min_base_quality, max_low_quality_fraction=max_low_quality_fraction, quality_key=quality_key, offset=offset) || return nothing
    return processed
end

"""
    sequencing_pipeline(records; kwargs...)

Apply sequencing cleanup to a batch of FASTQ or lightweight records.
"""
function sequencing_pipeline(
    records::AbstractVector{T};
    adapter::Union{Nothing,AbstractString}=nothing,
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
) where {T<:Union{FastqRecord,SeqRecordLite}}
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
        )
        filtered === nothing && continue
        push!(processed, filtered)
    end
    return processed
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
function trim_low_quality(record::FastqRecord; window::Int=4, threshold::Real=20, quality_key::Symbol=:quality, offset::Integer=33)
    scores = phred_scores(record.quality; offset=offset)
    trimmed = _trim_window(scores, window, threshold)
    trimmed === nothing && return nothing

    start_index, stop_index = trimmed
    return FastqRecord(
        record.identifier,
        record.description,
        record.sequence[start_index:stop_index],
        record.quality[start_index:stop_index],
    )
end

"""
    trim_low_quality(record; kwargs...)

Trim low-quality ends from a lightweight record.
"""
function trim_low_quality(record::SeqRecordLite; quality_key::Symbol=:quality, window::Int=4, threshold::Real=20, offset::Integer=33)
    quality = get(record.letter_annotations, quality_key, nothing)
    quality === nothing && throw(ArgumentError("missing $(quality_key) letter annotation"))
    quality isa AbstractString || throw(ArgumentError("FASTQ quality annotation must be a string"))

    scores = phred_scores(quality; offset=offset)
    trimmed = _trim_window(scores, window, threshold)
    trimmed === nothing && return nothing

    start_index, stop_index = trimmed
    trimmed_annotations = copy(record.annotations)
    trimmed_letter_annotations = copy(record.letter_annotations)
    trimmed_letter_annotations[quality_key] = quality[start_index:stop_index]

    return SeqRecordLite(
        record.sequence[start_index:stop_index];
        identifier=record.identifier,
        name=record.name,
        description=record.description,
        annotations=trimmed_annotations,
        letter_annotations=trimmed_letter_annotations,
    )
end