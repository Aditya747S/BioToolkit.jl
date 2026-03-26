module GenomicRanges

import ..BioToolkit: normalize_interval, BedRecord, GffRecord
using PooledArrays
using DataFrames

const PooledStringVector = PooledVector{String,UInt32,Vector{UInt32}}

export GenomicInterval, IntervalCollection, CoverageSegment, SeqInfo
export build_collection, read_intervals
export overlap, find_overlaps, nearest, find_nearest, follow, precede
export shift, flank, resize, promoters, narrow
export trim, gaps, complement, disjoin, pintersect, punion, psetdiff
export coverage

struct GenomicInterval
    chrom::String
    left::Int
    right::Int
    strand::Char
    metadata::Dict{String,Any}
end

struct IntervalCollection
    chroms::PooledStringVector
    starts::Vector{Int}
    ends::Vector{Int}
    strands::Vector{Char}
    metadata::Vector{Dict{String,Any}}
    intervals::Vector{GenomicInterval}
    chrom_indices::Dict{String,UnitRange{Int}}
    chrom_end_indices::Dict{String,Vector{Int}}
end

struct CoverageSegment
    chrom::String
    start::Int
    stop::Int
    depth::Int
end

struct SeqInfo
    chrom::String
    length::Int
    is_circular::Bool
end

GenomicInterval(chrom::AbstractString, left::Integer, right::Integer) = GenomicInterval(chrom, left, right, '.', Dict{String,Any}())
GenomicInterval(chrom::AbstractString, left::Integer, right::Integer, strand::Char) = GenomicInterval(chrom, left, right, strand, Dict{String,Any}())

function GenomicInterval(
    chrom::AbstractString,
    left::Integer,
    right::Integer,
    strand::Char,
    metadata::AbstractDict,
)
    strand in ('+', '-', '.') || throw(ArgumentError("strand must be '+', '-', or '.'"))
    left_int, right_int = normalize_interval(Int(left), Int(right))
    return GenomicInterval(String(chrom), left_int, right_int, strand, Dict{String,Any}(metadata))
end

GenomicInterval(record::BedRecord) = GenomicInterval(record.chrom, Int(record.start) + 1, Int(record.stop), '.')

function GenomicInterval(record::GffRecord)
    strand = isempty(record.strand) ? '.' : first(record.strand)
    strand in ('+', '-', '.') || (strand = '.')
    return GenomicInterval(record.chrom, Int(record.start), Int(record.stop), strand, record.attribute_map)
end

function Base.:(==)(left::GenomicInterval, right::GenomicInterval)
    return isequal(left.chrom, right.chrom) && isequal(left.left, right.left) && isequal(left.right, right.right) && isequal(left.strand, right.strand) && isequal(left.metadata, right.metadata)
end

function Base.hash(interval::GenomicInterval, state::UInt)
    state = hash(interval.chrom, state)
    state = hash(interval.left, state)
    state = hash(interval.right, state)
    state = hash(interval.strand, state)
    for (key, value) in sort(collect(pairs(interval.metadata)); by=first)
        state = hash(key, state)
        state = hash(value, state)
    end
    return state
end

function Base.show(io::IO, interval::GenomicInterval)
    print(io, "GenomicInterval(", interval.chrom, ":", interval.left, "-", interval.right, ", strand=", interval.strand, ")")
end

function Base.show(io::IO, collection::IntervalCollection)
    print(io, "IntervalCollection(", length(collection.intervals), " intervals, chromosomes=", length(collection.chrom_indices), ")")
end

function Base.show(io::IO, segment::CoverageSegment)
    print(io, "CoverageSegment(", segment.chrom, ":", segment.start, "-", segment.stop, ", depth=", segment.depth, ")")
end

function Base.show(io::IO, info::SeqInfo)
    print(io, "SeqInfo(", info.chrom, ", length=", info.length, ", circular=", info.is_circular, ")")
end

Base.length(collection::IntervalCollection) = length(collection.intervals)
Base.isempty(collection::IntervalCollection) = isempty(collection.intervals)
Base.iterate(collection::IntervalCollection, state...) = iterate(collection.intervals, state...)
Base.getindex(collection::IntervalCollection, index::Integer) = collection.intervals[index]

function _copy_metadata(interval::GenomicInterval)
    return Dict{String,Any}(interval.metadata)
end

function _new_interval(interval::GenomicInterval, left::Integer, right::Integer; strand::Char=interval.strand, metadata::AbstractDict=interval.metadata)
    return GenomicInterval(interval.chrom, left, right, strand, metadata)
end

function _map_intervals(transform, intervals::AbstractVector{<:GenomicInterval})
    return [transform(interval) for interval in intervals]
end

shift(intervals::AbstractVector{<:GenomicInterval}, delta::Integer) = _map_intervals(intervals) do interval
    _new_interval(interval, interval.left + Int(delta), interval.right + Int(delta))
end

shift(collection::IntervalCollection, delta::Integer) = build_collection(shift(collection.intervals, delta))

function flank(intervals::AbstractVector{<:GenomicInterval}, width::Integer; start::Bool=true)
    width > 0 || throw(ArgumentError("width must be positive"))
    result = GenomicInterval[]
    for interval in intervals
        if interval.strand == '-'
            if start
                push!(result, GenomicInterval(interval.chrom, interval.right + 1, interval.right + Int(width), interval.strand, _copy_metadata(interval)))
            else
                push!(result, GenomicInterval(interval.chrom, interval.left - Int(width), interval.left - 1, interval.strand, _copy_metadata(interval)))
            end
        else
            if start
                push!(result, GenomicInterval(interval.chrom, interval.left - Int(width), interval.left - 1, interval.strand, _copy_metadata(interval)))
            else
                push!(result, GenomicInterval(interval.chrom, interval.right + 1, interval.right + Int(width), interval.strand, _copy_metadata(interval)))
            end
        end
    end
    return result
end

flank(collection::IntervalCollection, width::Integer; start::Bool=true) = build_collection(flank(collection.intervals, width; start=start))

function resize(intervals::AbstractVector{<:GenomicInterval}, width::Integer; fix::Symbol=:start)
    width > 0 || throw(ArgumentError("width must be positive"))
    fix in (:start, :center, :end) || throw(ArgumentError("fix must be :start, :center, or :end"))
    result = GenomicInterval[]
    for interval in intervals
        new_left, new_right = if fix == :start
            interval.left, interval.left + Int(width) - 1
        elseif fix == :end
            interval.right - Int(width) + 1, interval.right
        else
            midpoint = (interval.left + interval.right) ÷ 2
            half_width = Int(width) ÷ 2
            if isodd(width)
                midpoint - half_width, midpoint + half_width
            else
                midpoint - half_width + 1, midpoint + half_width
            end
        end
        new_left <= new_right || throw(ArgumentError("width is too small for the requested fix"))
        push!(result, GenomicInterval(interval.chrom, new_left, new_right, interval.strand, _copy_metadata(interval)))
    end
    return result
end

resize(collection::IntervalCollection, width::Integer; fix::Symbol=:start) = build_collection(resize(collection.intervals, width; fix=fix))

function promoters(intervals::AbstractVector{<:GenomicInterval}, upstream::Integer, downstream::Integer)
    upstream >= 0 || throw(ArgumentError("upstream must be nonnegative"))
    downstream >= 0 || throw(ArgumentError("downstream must be nonnegative"))
    result = GenomicInterval[]
    for interval in intervals
        if interval.strand == '-'
            push!(result, GenomicInterval(interval.chrom, interval.right - Int(downstream), interval.right + Int(upstream) - 1, interval.strand, _copy_metadata(interval)))
        else
            push!(result, GenomicInterval(interval.chrom, interval.left - Int(upstream), interval.left + Int(downstream) - 1, interval.strand, _copy_metadata(interval)))
        end
    end
    return result
end

promoters(collection::IntervalCollection, upstream::Integer, downstream::Integer) = build_collection(promoters(collection.intervals, upstream, downstream))

function narrow(intervals::AbstractVector{<:GenomicInterval}, start::Integer=0, stop::Integer=0)
    start >= 0 || throw(ArgumentError("start must be nonnegative"))
    stop >= 0 || throw(ArgumentError("stop must be nonnegative"))
    result = GenomicInterval[]
    for interval in intervals
        new_left = interval.left + Int(start)
        new_right = interval.right - Int(stop)
        new_left <= new_right || continue
        push!(result, GenomicInterval(interval.chrom, new_left, new_right, interval.strand, _copy_metadata(interval)))
    end
    return result
end

narrow(collection::IntervalCollection, start::Integer=0, stop::Integer=0) = build_collection(narrow(collection.intervals, start, stop))

function _sort_key(interval::GenomicInterval)
    return (interval.chrom, interval.left, interval.right, interval.strand)
end

function _build_chrom_indices(intervals::Vector{GenomicInterval})
    chrom_indices = Dict{String,UnitRange{Int}}()
    chrom_end_indices = Dict{String,Vector{Int}}()
    isempty(intervals) && return chrom_indices, chrom_end_indices

    index = 1
    while index <= length(intervals)
        chrom = intervals[index].chrom
        stop = index
        while stop < length(intervals) && intervals[stop + 1].chrom == chrom
            stop += 1
        end
        chrom_indices[chrom] = index:stop

        end_order = collect(index:stop)
        sort!(end_order, by = value -> (intervals[value].right, intervals[value].left, intervals[value].strand))
        chrom_end_indices[chrom] = end_order

        index = stop + 1
    end

    return chrom_indices, chrom_end_indices
end

function build_collection(intervals::AbstractVector{<:GenomicInterval})
    sorted = GenomicInterval[interval for interval in intervals]
    sort!(sorted, by=_sort_key)
    chrom_indices, chrom_end_indices = _build_chrom_indices(sorted)
    chroms = PooledArray([interval.chrom for interval in sorted])
    starts = Int[interval.left for interval in sorted]
    ends = Int[interval.right for interval in sorted]
    strands = Char[interval.strand for interval in sorted]
    metadata = [Dict{String,Any}(interval.metadata) for interval in sorted]
    return IntervalCollection(chroms, starts, ends, strands, metadata, sorted, chrom_indices, chrom_end_indices)
end

build_collection(intervals) = build_collection(collect(intervals))

IntervalCollection(intervals::AbstractVector{<:GenomicInterval}) = build_collection(intervals)

function _find_last_start_leq(starts::Vector{Int}, range::UnitRange{Int}, bound::Int)
    lo = first(range)
    hi = last(range)
    lo > hi && return nothing
    answer = nothing

    while lo <= hi
        mid = (lo + hi) >>> 1
        if starts[mid] <= bound
            answer = mid
            lo = mid + 1
        else
            hi = mid - 1
        end
    end

    return answer
end

function _find_first_start_gt(starts::Vector{Int}, range::UnitRange{Int}, bound::Int)
    lo = first(range)
    hi = last(range)
    lo > hi && return nothing
    answer = nothing

    while lo <= hi
        mid = (lo + hi) >>> 1
        if starts[mid] > bound
            answer = mid
            hi = mid - 1
        else
            lo = mid + 1
        end
    end

    return answer
end

function _find_last_end_lt(ends::Vector{Int}, ordered_indices::Vector{Int}, bound::Int)
    lo = firstindex(ordered_indices)
    hi = lastindex(ordered_indices)
    lo > hi && return nothing
    answer = nothing

    while lo <= hi
        mid = (lo + hi) >>> 1
        index = ordered_indices[mid]
        if ends[index] < bound
            answer = index
            lo = mid + 1
        else
            hi = mid - 1
        end
    end

    return answer
end

function find_overlaps(query::GenomicInterval, subject::IntervalCollection)
    range = get(subject.chrom_indices, query.chrom, nothing)
    range === nothing && return GenomicInterval[]

    starts = subject.starts
    ends = subject.ends
    intervals = subject.intervals
    last_candidate = _find_last_start_leq(starts, range, query.right)
    last_candidate === nothing && return GenomicInterval[]

    hits = GenomicInterval[]
    @inbounds for index in first(range):last_candidate
        ends[index] < query.left && continue
        starts[index] > query.right && break
        push!(hits, intervals[index])
    end

    return hits
end

overlap(query::GenomicInterval, subject::IntervalCollection) = find_overlaps(query, subject)

function _distance(query::GenomicInterval, interval::GenomicInterval)
    interval.right < query.left && return query.left - interval.right - 1
    interval.left > query.right && return interval.left - query.right - 1
    return 0
end

function nearest(query::GenomicInterval, subject::IntervalCollection)
    range = get(subject.chrom_indices, query.chrom, nothing)
    range === nothing && return nothing
    intervals = subject.intervals
    starts = subject.starts
    ends = subject.ends
    best_interval = nothing
    best_key = nothing

    @inbounds for index in range
        interval = intervals[index]
        key = if starts[index] <= query.right && ends[index] >= query.left
            (0, starts[index], ends[index], Int(interval.strand))
        elseif ends[index] < query.left
            (1, _distance(query, interval), -ends[index], starts[index], ends[index])
        else
            (2, _distance(query, interval), starts[index], ends[index])
        end

        if best_key === nothing || key < best_key
            best_interval = interval
            best_key = key
        end
    end

    return best_interval
end

find_nearest(query::GenomicInterval, subject::IntervalCollection) = nearest(query, subject)

function precede(query::GenomicInterval, subject::IntervalCollection)
    range = get(subject.chrom_indices, query.chrom, nothing)
    range === nothing && return nothing
    intervals = subject.intervals
    upstream = _find_last_end_lt(subject.ends, subject.chrom_end_indices[query.chrom], query.left)
    upstream === nothing && return nothing
    return intervals[upstream]
end

function follow(query::GenomicInterval, subject::IntervalCollection)
    range = get(subject.chrom_indices, query.chrom, nothing)
    range === nothing && return nothing
    intervals = subject.intervals
    downstream = _find_first_start_gt(subject.starts, range, query.right)
    downstream === nothing && return nothing
    return intervals[downstream]
end

function _subtract_piece(piece::GenomicInterval, blocker::GenomicInterval)
    piece.chrom == blocker.chrom || return [piece]
    (piece.left <= blocker.right && piece.right >= blocker.left) || return [piece]

    pieces = GenomicInterval[]
    if blocker.left > piece.left
        push!(pieces, GenomicInterval(piece.chrom, piece.left, min(blocker.left - 1, piece.right), piece.strand, piece.metadata))
    end
    if blocker.right < piece.right
        push!(pieces, GenomicInterval(piece.chrom, max(blocker.right + 1, piece.left), piece.right, piece.strand, piece.metadata))
    end
    return pieces
end

function _subtract_many(piece::GenomicInterval, blockers::Vector{GenomicInterval})
    result = [piece]
    for blocker in blockers
        next_result = GenomicInterval[]
        for candidate in result
            append!(next_result, _subtract_piece(candidate, blocker))
        end
        result = next_result
        isempty(result) && break
    end
    return result
end

function _positions_by_chrom(intervals)
    positions = Dict{String,Set{Int}}()
    for interval in intervals
        chrom_positions = get!(positions, interval.chrom, Set{Int}())
        for position in interval.left:interval.right
            push!(chrom_positions, position)
        end
    end
    return positions
end

function _subtract_positions(left_positions::Dict{String,Set{Int}}, right_positions::Dict{String,Set{Int}})
    result = Dict{String,Set{Int}}()
    for (chrom, positions) in left_positions
        blocked = get(right_positions, chrom, Set{Int}())
        remaining = Set{Int}()
        for position in positions
            position in blocked && continue
            push!(remaining, position)
        end
        result[chrom] = remaining
    end
    return result
end

function _intervals_from_positions(positions::Dict{String,Set{Int}})
    intervals = GenomicInterval[]
    for chrom in sort(collect(keys(positions)))
        chrom_positions = sort(collect(positions[chrom]))
        isempty(chrom_positions) && continue

        start = chrom_positions[1]
        previous = start
        for position in chrom_positions[2:end]
            if position != previous + 1
                push!(intervals, GenomicInterval(chrom, start, previous, '.', Dict{String,Any}()))
                start = position
            end
            previous = position
        end
        push!(intervals, GenomicInterval(chrom, start, previous, '.', Dict{String,Any}()))
    end
    return intervals
end

function Base.setdiff(left::GenomicInterval, right::IntervalCollection)
    left_positions = _positions_by_chrom([left])
    right_positions = _positions_by_chrom(reduce(right).intervals)
    remaining_positions = _subtract_positions(left_positions, right_positions)
    return _intervals_from_positions(remaining_positions)
end

function Base.setdiff(left::IntervalCollection, right::IntervalCollection)
    left_positions = _positions_by_chrom(left.intervals)
    right_positions = _positions_by_chrom(reduce(right).intervals)
    remaining_positions = _subtract_positions(left_positions, right_positions)
    return build_collection(_intervals_from_positions(remaining_positions))
end

function _merge_metadata(metadata::Dict{String,Any})
    merged = Dict{String,Any}(metadata)
    merged["merged_from"] = get(merged, "merged_from", 1) + 1
    return merged
end

function Base.reduce(collection::IntervalCollection)
    isempty(collection) && return collection

    merged = GenomicInterval[]
    current = collection.intervals[1]

    @inbounds for index in 2:length(collection)
        candidate = collection.intervals[index]
        if candidate.chrom == current.chrom && candidate.left <= current.right + 1
            current = GenomicInterval(
                current.chrom,
                current.left,
                max(current.right, candidate.right),
                current.strand,
                _merge_metadata(current.metadata),
            )
        else
            push!(merged, current)
            current = candidate
        end
    end

    push!(merged, current)
    return build_collection(merged)
end

function _interval_collection_intersect(left::IntervalCollection, right::IntervalCollection)
    (isempty(left) || isempty(right)) && return build_collection(GenomicInterval[])

    pieces = GenomicInterval[]
    for left_interval in left.intervals
        for hit in find_overlaps(left_interval, right)
            push!(pieces, GenomicInterval(left_interval.chrom, max(left_interval.left, hit.left), min(left_interval.right, hit.right), left_interval.strand, _copy_metadata(left_interval)))
        end
    end
    return build_collection(pieces)
end

function Base.intersect(left::IntervalCollection, right::IntervalCollection)
    return _interval_collection_intersect(left, right)
end

function Base.union(left::IntervalCollection, right::IntervalCollection)
    return reduce(build_collection(vcat(left.intervals, right.intervals)))
end

function _coverage_segments_for_chrom(chrom::String, intervals::Vector{GenomicInterval}, range::UnitRange{Int})
    events = Dict{Int,Int}()

    @inbounds for index in range
        interval = intervals[index]
        events[interval.left] = get(events, interval.left, 0) + 1
        events[interval.right + 1] = get(events, interval.right + 1, 0) - 1
    end

    positions = sort!(collect(keys(events)))
    segments = CoverageSegment[]
    depth = 0
    previous_position = nothing

    for position in positions
        if previous_position !== nothing && depth > 0 && position > previous_position
            push!(segments, CoverageSegment(chrom, previous_position, position - 1, depth))
        end
        depth += events[position]
        previous_position = position
    end

    return segments
end

function find_overlaps_parallel(query::IntervalCollection, subject::IntervalCollection)
    results = Vector{Vector{GenomicInterval}}(undef, length(query))
    chrom_groups = Dict{String,Vector{Int}}()
    for (index, chrom) in enumerate(query.chroms)
        push!(get!(chrom_groups, chrom, Int[]), index)
    end

    chromosomes = collect(keys(chrom_groups))
    Threads.@threads for chrom_index in eachindex(chromosomes)
        chrom = chromosomes[chrom_index]
        for query_index in chrom_groups[chrom]
            results[query_index] = find_overlaps(query[query_index], subject)
        end
    end

    return results
end

function trim(intervals::AbstractVector{<:GenomicInterval}, seqlengths::Dict{String,Int})
    result = GenomicInterval[]
    for interval in intervals
        limit = get(seqlengths, interval.chrom, nothing)
        limit === nothing && continue
        new_left = max(1, interval.left)
        new_right = min(limit, interval.right)
        new_left <= new_right || continue
        push!(result, GenomicInterval(interval.chrom, new_left, new_right, interval.strand, _copy_metadata(interval)))
    end
    return result
end

trim(collection::IntervalCollection, seqlengths::Dict{String,Int}) = build_collection(trim(collection.intervals, seqlengths))
trim(intervals::AbstractVector{<:GenomicInterval}, seqinfo::Dict{String,SeqInfo}) = trim(intervals, Dict(chrom => info.length for (chrom, info) in seqinfo))
trim(collection::IntervalCollection, seqinfo::Dict{String,SeqInfo}) = trim(collection, Dict(chrom => info.length for (chrom, info) in seqinfo))

function gaps(intervals::AbstractVector{<:GenomicInterval}, seqlengths::Dict{String,Int})
    return gaps(build_collection(intervals), seqlengths)
end

function gaps(collection::IntervalCollection, seqlengths::Dict{String,Int})
    collection = reduce(collection)
    gap_intervals = GenomicInterval[]
    for (chrom, chrom_length) in sort(collect(pairs(seqlengths)); by=first)
        chrom_intervals = get(collection.chrom_indices, chrom, nothing)
        if chrom_intervals === nothing
            push!(gap_intervals, GenomicInterval(chrom, 1, chrom_length, '.', Dict{String,Any}()))
            continue
        end

        current = 1
        for interval in collection.intervals[chrom_intervals]
            if current < interval.left
                push!(gap_intervals, GenomicInterval(chrom, current, interval.left - 1, '.', Dict{String,Any}()))
            end
            current = max(current, interval.right + 1)
        end
        if current <= chrom_length
            push!(gap_intervals, GenomicInterval(chrom, current, chrom_length, '.', Dict{String,Any}()))
        end
    end
    return build_collection(gap_intervals)
end

gaps(intervals::AbstractVector{<:GenomicInterval}, seqinfo::Dict{String,SeqInfo}) = gaps(intervals, Dict(chrom => info.length for (chrom, info) in seqinfo))
gaps(collection::IntervalCollection, seqinfo::Dict{String,SeqInfo}) = gaps(collection, Dict(chrom => info.length for (chrom, info) in seqinfo))

complement(intervals::AbstractVector{<:GenomicInterval}, seqlengths::Dict{String,Int}) = gaps(intervals, seqlengths)
complement(collection::IntervalCollection, seqlengths::Dict{String,Int}) = gaps(collection, seqlengths)
complement(intervals::AbstractVector{<:GenomicInterval}, seqinfo::Dict{String,SeqInfo}) = gaps(intervals, seqinfo)
complement(collection::IntervalCollection, seqinfo::Dict{String,SeqInfo}) = gaps(collection, seqinfo)

function disjoin(collection::IntervalCollection)
    breakpoints = Dict{String,Set{Int}}()
    for interval in collection.intervals
        chrom_points = get!(breakpoints, interval.chrom, Set{Int}())
        push!(chrom_points, interval.left)
        push!(chrom_points, interval.right + 1)
    end

    pieces = GenomicInterval[]
    for chrom in sort(collect(keys(breakpoints)))
        points = sort(collect(breakpoints[chrom]))
        length(points) < 2 && continue
        for index in 1:length(points)-1
            left = points[index]
            right = points[index + 1] - 1
            left > right && continue
            candidate = GenomicInterval(chrom, left, right, '.', Dict{String,Any}())
            any(interval -> interval.left <= candidate.left && interval.right >= candidate.right, collection.intervals) && push!(pieces, candidate)
        end
    end

    return build_collection(pieces)
end

function _pairwise_zip(x::AbstractVector{<:GenomicInterval}, y::AbstractVector{<:GenomicInterval})
    length(x) == length(y) || throw(ArgumentError("interval vectors must have the same length"))
    return zip(x, y)
end

function pintersect(x::AbstractVector{<:GenomicInterval}, y::AbstractVector{<:GenomicInterval})
    results = Union{Nothing,GenomicInterval}[]
    for (left, right) in _pairwise_zip(x, y)
        if left.chrom == right.chrom && left.left <= right.right && left.right >= right.left
            push!(results, GenomicInterval(left.chrom, max(left.left, right.left), min(left.right, right.right), left.strand, _copy_metadata(left)))
        else
            push!(results, nothing)
        end
    end
    return results
end

function punion(x::AbstractVector{<:GenomicInterval}, y::AbstractVector{<:GenomicInterval})
    results = GenomicInterval[]
    for (left, right) in _pairwise_zip(x, y)
        left.chrom == right.chrom || throw(ArgumentError("intervals must be on the same chromosome"))
        push!(results, GenomicInterval(left.chrom, min(left.left, right.left), max(left.right, right.right), left.strand, _copy_metadata(left)))
    end
    return results
end

function psetdiff(x::AbstractVector{<:GenomicInterval}, y::AbstractVector{<:GenomicInterval})
    results = Vector{GenomicInterval}[]
    for (left, right) in _pairwise_zip(x, y)
        push!(results, setdiff(left, build_collection([right])))
    end
    return results
end

function DataFrames.DataFrame(intervals::AbstractVector{<:GenomicInterval})
    rows = collect(intervals)
    metadata_keys = Set{String}()
    for interval in rows
        union!(metadata_keys, keys(interval.metadata))
    end

    columns = Dict{Symbol,Vector{Any}}(
        :chrom => Any[interval.chrom for interval in rows],
        :start => Any[interval.left for interval in rows],
        :stop => Any[interval.right for interval in rows],
        :strand => Any[interval.strand for interval in rows],
    )

    for key in sort(collect(metadata_keys))
        columns[Symbol(key)] = Any[get(interval.metadata, key, missing) for interval in rows]
    end

    return DataFrames.DataFrame(columns)
end

function read_intervals(df::DataFrames.AbstractDataFrame)
    names_map = Dict(lowercase(String(name)) => name for name in DataFrames.names(df))
    chrom_name = get(names_map, "chrom", get(names_map, "chr", nothing))
    start_name = get(names_map, "start", get(names_map, "left", nothing))
    stop_name = get(names_map, "stop", get(names_map, "end", get(names_map, "right", nothing)))
    strand_name = get(names_map, "strand", nothing)

    chrom_name === nothing && throw(ArgumentError("DataFrame must include a Chrom or chr column"))
    start_name === nothing && throw(ArgumentError("DataFrame must include a Start or left column"))
    stop_name === nothing && throw(ArgumentError("DataFrame must include a Stop, End, or right column"))

    excluded = Set(filter(!isnothing, [chrom_name, start_name, stop_name, strand_name]))
    intervals = GenomicInterval[]

    for row in eachrow(df)
        chrom = String(row[chrom_name])
        left = Int(row[start_name])
        right = Int(row[stop_name])
        strand = strand_name === nothing || ismissing(row[strand_name]) ? '.' : Char(row[strand_name])
        metadata = Dict{String,Any}()

        for name in DataFrames.names(df)
            name in excluded && continue
            value = row[name]
            ismissing(value) && continue
            metadata[String(name)] = value
        end

        push!(intervals, GenomicInterval(chrom, left, right, strand, metadata))
    end

    return intervals
end

function coverage(collection::IntervalCollection)
    segments = CoverageSegment[]
    for chrom in sort(collect(keys(collection.chrom_indices)))
        append!(segments, _coverage_segments_for_chrom(chrom, collection.intervals, collection.chrom_indices[chrom]))
    end
    return segments
end

coverage(intervals::AbstractVector{<:GenomicInterval}) = coverage(build_collection(intervals))

function _parse_interval_row(headers::Vector{String}, fields::Vector{String})
    column_map = Dict(lowercase(strip(header)) => index for (index, header) in pairs(headers))

    chrom_index = get(column_map, "chr", get(column_map, "chrom", get(column_map, "chromosome", nothing)))
    start_index = get(column_map, "start", get(column_map, "left", get(column_map, "begin", nothing)))
    stop_index = get(column_map, "end", get(column_map, "stop", get(column_map, "right", nothing)))
    strand_index = get(column_map, "strand", nothing)

    chrom_index === nothing && throw(ArgumentError("interval header must include a chrom column"))
    start_index === nothing && throw(ArgumentError("interval header must include a start column"))
    stop_index === nothing && throw(ArgumentError("interval header must include an end column"))

    chrom = fields[chrom_index]
    left = parse(Int, fields[start_index])
    right = parse(Int, fields[stop_index])
    strand = strand_index === nothing || isempty(fields[strand_index]) ? '.' : first(fields[strand_index])

    metadata = Dict{String,Any}()
    for (index, header) in pairs(headers)
        index in (chrom_index, start_index, stop_index, strand_index) && continue
        index <= length(fields) || continue
        metadata[strip(header)] = fields[index]
    end

    return GenomicInterval(chrom, left, right, strand, metadata)
end

function read_intervals(io::IO)
    headers = nothing
    intervals = GenomicInterval[]

    for raw_line in eachline(io)
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, '#') && continue

        if headers === nothing
            headers = String[strip(field) for field in Base.split(line, ',')]
            continue
        end

        fields = String[strip(field) for field in Base.split(line, ',')]
        length(fields) >= length(headers) || throw(ArgumentError("interval row has fewer fields than header"))
        push!(intervals, _parse_interval_row(headers, fields))
    end

    headers === nothing && return GenomicInterval[]
    return intervals
end

function read_intervals(path::AbstractString)
    open(path, "r") do io
        return read_intervals(io)
    end
end

end