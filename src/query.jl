using BigWig

const _HAS_CUDA_QUERY = false

"""
    _window_coverage_serial(starts, stops, window_size)

Compute per-window coverage counts using the serial CPU path.
"""
function _window_coverage_serial(starts, stops, window_size::Integer)
    isempty(starts) && return Dict{Int,Int}()

    max_stop = Int(maximum(stops))
    max_window = fld(max_stop - 1, window_size)
    coverage = Vector{Int}(undef, max_window + 1)
    fill!(coverage, 0)

    @inbounds for index in eachindex(starts, stops)
        start = Int(starts[index])
        stop = Int(stops[index])
        stop <= start && continue

        first_window = fld(start, window_size)
        last_window = fld(stop - 1, window_size)

        for window_index in first_window:last_window
            window_start = window_index * window_size
            window_stop = window_start + window_size
            overlap_start = max(start, window_start)
            overlap_stop = min(stop, window_stop)
            overlap = overlap_stop - overlap_start
            overlap > 0 && (coverage[window_index + 1] += overlap)
        end
    end

    result = Dict{Int,Int}()
    @inbounds for (index, count) in pairs(coverage)
        count == 0 && continue
        result[index - 1] = count
    end

    return result
end

"""
    _filter_region_mask(columns, chrom, min_pos, max_pos)

Filter a tabular region using a boolean mask.
"""
function _filter_region_mask(columns, chrom::AbstractString, min_pos::Integer, max_pos::Integer)
    mask = (columns.chrom .== chrom) .& (columns.pos .>= min_pos) .& (columns.pos .<= max_pos)
    indices = findall(mask)

    return (
        chrom = columns.chrom[indices],
        pos = columns.pos[indices],
        id = columns.id[indices],
        ref = columns.ref[indices],
        alt = columns.alt[indices],
        qual = columns.qual[indices],
    )
end

"""
    _filter_region_sorted(columns, chrom, min_pos, max_pos)

Filter a sorted tabular region using binary search.
"""
function _filter_region_sorted(columns, chrom::AbstractString, min_pos::Integer, max_pos::Integer)
    chroms = columns.chrom
    start = searchsortedfirst(chroms, chrom)
    stop = searchsortedlast(chroms, chrom)
    start > stop && return _filter_region_mask(columns, chrom, min_pos, max_pos)

    positions = columns.pos[start:stop]
    pos_start = searchsortedfirst(positions, min_pos)
    pos_stop = searchsortedlast(positions, max_pos)
    pos_start > pos_stop && return (
        chrom = chroms[1:0],
        pos = columns.pos[1:0],
        id = columns.id[1:0],
        ref = columns.ref[1:0],
        alt = columns.alt[1:0],
        qual = columns.qual[1:0],
    )

    first_index = start + pos_start - 1
    last_index = start + pos_stop - 1
    return (
        chrom = chroms[first_index:last_index],
        pos = columns.pos[first_index:last_index],
        id = columns.id[first_index:last_index],
        ref = columns.ref[first_index:last_index],
        alt = columns.alt[first_index:last_index],
        qual = columns.qual[first_index:last_index],
    )
end

"""
    filter_region(table, chrom, min_pos, max_pos; sorted=false)

Filter rows in a genomic table by chromosome and coordinate range.
"""
function filter_region(
    table,
    chrom::AbstractString,
    min_pos::Integer,
    max_pos::Integer;
    sorted::Bool=false,
)
    columns = Tables.columntable(table)
    return sorted ? _filter_region_sorted(columns, chrom, min_pos, max_pos) : _filter_region_mask(columns, chrom, min_pos, max_pos)
end

"""
    _bin_positions_serial(positions, bin_size)

Count positions into integer bins on the CPU.
"""
function _bin_positions_serial(positions, bin_size::Integer)
    histogram = Dict{Int,Int}()

    @inbounds for position in positions
        bin_index = Int(fld(position, bin_size))
        histogram[bin_index] = get(histogram, bin_index, 0) + 1
    end

    return histogram
end

"""
    _bin_positions_threaded(positions, bin_size)

Count positions into integer bins using threaded CPU work.
"""
function _bin_positions_threaded(positions, bin_size::Integer)
    isempty(positions) && return Dict{Int,Int}()

    nthreads = Threads.nthreads()
    partials = [Dict{Int,Int}() for _ in 1:nthreads]
    len = length(positions)

    Threads.@threads :static for thread_index in 1:nthreads
        start_index = firstindex(positions) + div((thread_index - 1) * len, nthreads)
        stop_index  = firstindex(positions) + div(thread_index       * len, nthreads) - 1
        thread_index == nthreads && (stop_index = lastindex(positions))

        histogram = partials[thread_index]
        @inbounds for index in start_index:stop_index
            bin_index = Int(fld(Int(positions[index]), bin_size))
            histogram[bin_index] = get(histogram, bin_index, 0) + 1
        end
    end

    # Merge thread-local histograms
    histogram = Dict{Int,Int}()
    @inbounds for partial in partials
        for (bin_index, count) in partial
            histogram[bin_index] = get(histogram, bin_index, 0) + count
        end
    end

    return histogram
end

"""
    _window_coverage_threaded(starts, stops, window_size)

Compute per-window coverage counts using threaded CPU work.
"""
function _window_coverage_threaded(starts, stops, window_size::Integer)
    isempty(starts) && return Dict{Int,Int}()

    nthreads = Threads.nthreads()
    partials = [Dict{Int,Int}() for _ in 1:nthreads]
    len = length(starts)

    Threads.@threads :static for thread_index in 1:nthreads
        start_index = firstindex(starts) + div((thread_index - 1) * len, nthreads)
        stop_index = firstindex(starts) + div(thread_index * len, nthreads) - 1
        thread_index == nthreads && (stop_index = lastindex(starts))

        coverage = partials[thread_index]
        @inbounds for index in start_index:stop_index
            start = Int(starts[index])
            stop = Int(stops[index])
            stop <= start && continue

            first_window = fld(start, window_size)
            last_window = fld(stop - 1, window_size)

            for window_index in first_window:last_window
                window_start = window_index * window_size
                window_stop = window_start + window_size
                overlap_start = max(start, window_start)
                overlap_stop = min(stop, window_stop)
                overlap = overlap_stop - overlap_start
                window_key = window_index + 1
                overlap > 0 && (coverage[window_key] = get(coverage, window_key, 0) + overlap)
            end
        end
    end

    result = Dict{Int,Int}()
    @inbounds for partial in partials
        for (window_index, count) in partial
            result[window_index - 1] = get(result, window_index - 1, 0) + count
        end
    end

    return result
end

"""
    bin_positions(positions, bin_size; threaded=true, use_cuda=false)

Count positions into integer bins, optionally using CUDA.
"""
function bin_positions(positions, bin_size::Integer; threaded::Bool=true, use_cuda::Bool=false)
    if use_cuda
        _ensure_cuda_query!()
        positions_cuda = _is_cuda_backed_array(positions) ? positions : CUDA.CuArray{Int}(Int.(positions))
        return Base.invokelatest(_CUDA_QUERY_BIN_IMPL[], positions_cuda, bin_size)
    end

    return threaded ? _bin_positions_threaded(positions, bin_size) : _bin_positions_serial(positions, bin_size)
end

"""
    coverage_histogram(table, chrom, bin_size; threaded=true, use_cuda=false)

Build a binned coverage histogram from a table of genomic intervals.
"""
function coverage_histogram(table, chrom::AbstractString, bin_size::Integer; threaded::Bool=true, use_cuda::Bool=false)
    subset = filter_region(table, chrom, typemin(Int), typemax(Int))
    return bin_positions(subset.pos, bin_size; threaded=threaded, use_cuda=use_cuda)
end

"""
    window_coverage(starts, stops, window_size; threaded=true, use_cuda=false)

Compute per-window coverage across genomic intervals.
"""
function window_coverage(starts, stops, window_size::Integer; threaded::Bool=true, use_cuda::Bool=false)
    window_size <= 0 && throw(ArgumentError("window_size must be positive"))
    length(starts) == length(stops) || throw(ArgumentError("starts and stops must have the same length"))

    if use_cuda
        _ensure_cuda_query!()
        starts_cuda = _is_cuda_backed_array(starts) ? starts : CUDA.CuArray{Int}(Int.(starts))
        stops_cuda = _is_cuda_backed_array(stops) ? stops : CUDA.CuArray{Int}(Int.(stops))
        return Base.invokelatest(_CUDA_QUERY_WINDOW_IMPL[], starts_cuda, stops_cuda, window_size)
    end

    return threaded ? _window_coverage_threaded(starts, stops, window_size) : _window_coverage_serial(starts, stops, window_size)
end

"""
    window_coverage(table, chrom, window_size; sorted=false, threaded=true, use_cuda=false)

Compute per-window coverage for a filtered chromosome slice.
"""
function window_coverage(table, chrom::AbstractString, window_size::Integer; sorted::Bool=false, threaded::Bool=true, use_cuda::Bool=false)
    window_size <= 0 && throw(ArgumentError("window_size must be positive"))

    columns = Tables.columntable(table)

    if sorted
        chroms = columns.chrom
        start = searchsortedfirst(chroms, chrom)
        stop = searchsortedlast(chroms, chrom)
        start > stop && return Dict{Int,Int}()
        row_indices = start:stop
    else
        row_indices = findall(columns.chrom .== chrom)
        isempty(row_indices) && return Dict{Int,Int}()
    end

    return window_coverage(columns.start[row_indices], columns.stop[row_indices], window_size; threaded=threaded, use_cuda=use_cuda)
end

"""
    _write_bigwig_runs(writer, coverage_vector; chrom="chr1", start=1)

Write a coverage vector to BigWig using run-length encoded spans.
"""
function _write_bigwig_runs(writer, coverage_vector::AbstractVector{<:Real}; chrom::AbstractString="chr1", start::Integer=1)
    isempty(coverage_vector) && return nothing

    index = firstindex(coverage_vector)
    current_value = Float32(coverage_vector[index])
    run_start = index

    for cursor in (index + 1):lastindex(coverage_vector)
        value = Float32(coverage_vector[cursor])
        if value != current_value
            write(writer, (String(chrom), Int(start) + run_start - 1, Int(start) + cursor - 2, current_value))
            run_start = cursor
            current_value = value
        end
    end

    write(writer, (String(chrom), Int(start) + run_start - 1, Int(start) + lastindex(coverage_vector) - 1, current_value))

    return nothing
end

"""
    write_bigwig(coverage_vector, output_path; chrom="chr1", start=1)

Write a single coverage vector to a BigWig file.
"""
function write_bigwig(coverage_vector::AbstractVector{<:Real}, output_path::AbstractString; chrom::AbstractString="chr1", start::Integer=1)
    open(output_path, "w") do io
        writer = BigWig.Writer(io, [(String(chrom), Int(start) + length(coverage_vector) - 1)])
        try
            _write_bigwig_runs(writer, coverage_vector; chrom=chrom, start=start)
        finally
            close(writer)
        end
    end
    return output_path
end

"""
    _bigwig_start(starts, chrom)

Resolve a BigWig starting coordinate from either an integer or a mapping.
"""
function _bigwig_start(starts, chrom)
    starts isa Integer && return Int(starts)
    starts isa AbstractDict || throw(ArgumentError("starts must be an integer or a dictionary"))

    if haskey(starts, chrom)
        return Int(starts[chrom])
    end

    chrom_name = String(chrom)
    for key in keys(starts)
        String(key) == chrom_name && return Int(starts[key])
    end

    throw(KeyError(chrom))
end

"""
    _write_bigwig_coverage_map(io, coverage_vectors; starts=1)

Write a chromosome-to-coverage map into a BigWig writer.
"""
function _write_bigwig_coverage_map(io::IO, coverage_vectors::AbstractDict; starts=1)
    entries = collect(pairs(coverage_vectors))
    sort!(entries, by = entry -> String(first(entry)))

    chromlist = Tuple{String,Int}[]
    for (chrom, coverage_vector) in entries
        coverage_vector isa AbstractVector{<:Real} || throw(ArgumentError("coverage values must be vectors of real numbers"))
        chrom_start = _bigwig_start(starts, chrom)
        push!(chromlist, (String(chrom), chrom_start + length(coverage_vector) - 1))
    end

    writer = BigWig.Writer(io, chromlist)
    try
        for (chrom, coverage_vector) in entries
            chrom_start = _bigwig_start(starts, chrom)
            _write_bigwig_runs(writer, coverage_vector; chrom=String(chrom), start=chrom_start)
        end
    finally
        close(writer)
    end

    return nothing
end

"""
    write_bigwig(coverage_vectors, output_path; starts=1)

Write a chromosome-to-coverage map to a BigWig file.
"""
function write_bigwig(coverage_vectors::AbstractDict, output_path::AbstractString; starts=1)
    open(output_path, "w") do io
        _write_bigwig_coverage_map(io, coverage_vectors; starts=starts)
    end
    return output_path
end

"""
    normalize_interval(start, stop)

Return an interval as ordered integer bounds.
"""
function normalize_interval(start::Integer, stop::Integer)
    start <= stop ? (Int(start), Int(stop)) : (Int(stop), Int(start))
end

"""
    normalize_interval(interval)

Return a tuple interval with ordered integer bounds.
"""
normalize_interval(interval::Tuple{<:Integer,<:Integer}) = normalize_interval(interval[1], interval[2])

"""
    interval_length(start, stop)

Return the half-open length of an interval.
"""
function interval_length(start::Integer, stop::Integer)
    start, stop = normalize_interval(start, stop)
    return max(0, stop - start)
end

"""
    interval_length(interval)

Return the half-open length of a tuple interval.
"""
interval_length(interval::Tuple{<:Integer,<:Integer}) = interval_length(interval[1], interval[2])

"""
    interval_contains(start, stop, position)

Test whether a position falls inside a half-open interval.
"""
function interval_contains(start::Integer, stop::Integer, position::Integer)
    start, stop = normalize_interval(start, stop)
    return start <= position < stop
end

"""
    interval_contains(interval, position)

Test whether a position falls inside a tuple interval.
"""
interval_contains(interval::Tuple{<:Integer,<:Integer}, position::Integer) = interval_contains(interval[1], interval[2], position)

"""
    interval_overlaps(left_start, left_stop, right_start, right_stop)

Test whether two half-open intervals overlap.
"""
function interval_overlaps(left_start::Integer, left_stop::Integer, right_start::Integer, right_stop::Integer)
    left_start, left_stop = normalize_interval(left_start, left_stop)
    right_start, right_stop = normalize_interval(right_start, right_stop)
    return max(left_start, right_start) < min(left_stop, right_stop)
end

"""
    interval_overlaps(left, right)

Test whether two tuple intervals overlap.
"""
interval_overlaps(left::Tuple{<:Integer,<:Integer}, right::Tuple{<:Integer,<:Integer}) = interval_overlaps(left[1], left[2], right[1], right[2])

"""
    interval_intersection(left_start, left_stop, right_start, right_stop)

Return the overlap of two intervals, if any.
"""
function interval_intersection(left_start::Integer, left_stop::Integer, right_start::Integer, right_stop::Integer)
    interval_overlaps(left_start, left_stop, right_start, right_stop) || return nothing
    left_start, left_stop = normalize_interval(left_start, left_stop)
    right_start, right_stop = normalize_interval(right_start, right_stop)
    return max(left_start, right_start), min(left_stop, right_stop)
end

"""
    interval_intersection(left, right)

Return the overlap of two tuple intervals, if any.
"""
interval_intersection(left::Tuple{<:Integer,<:Integer}, right::Tuple{<:Integer,<:Integer}) = interval_intersection(left[1], left[2], right[1], right[2])

"""
    interval_union(left_start, left_stop, right_start, right_stop)

Return the union of two overlapping or touching intervals.
"""
function interval_union(left_start::Integer, left_stop::Integer, right_start::Integer, right_stop::Integer)
    left_start, left_stop = normalize_interval(left_start, left_stop)
    right_start, right_stop = normalize_interval(right_start, right_stop)
    left_stop < right_start && return nothing
    right_stop < left_start && return nothing
    return min(left_start, right_start), max(left_stop, right_stop)
end

"""
    interval_union(left, right)

Return the union of two tuple intervals when they touch or overlap.
"""
interval_union(left::Tuple{<:Integer,<:Integer}, right::Tuple{<:Integer,<:Integer}) = interval_union(left[1], left[2], right[1], right[2])

"""
    merge_intervals(intervals)

Merge a collection of overlapping intervals into disjoint spans.
"""
function merge_intervals(intervals::AbstractVector{<:Tuple{<:Integer,<:Integer}})
    isempty(intervals) && return Tuple{Int,Int}[]
    sorted_intervals = [normalize_interval(interval) for interval in intervals]
    sort!(sorted_intervals, by = interval -> (interval[1], interval[2]))

    merged = Tuple{Int,Int}[]
    current_start, current_stop = sorted_intervals[1]

    for interval in @view sorted_intervals[2:end]
        start, stop = interval
        if start <= current_stop
            current_stop = max(current_stop, stop)
        else
            push!(merged, (current_start, current_stop))
            current_start, current_stop = start, stop
        end
    end

    push!(merged, (current_start, current_stop))
    return merged
end

"""
    interval_difference(left_start, left_stop, right_start, right_stop)

Subtract one interval from another and return the remaining pieces.
"""
function interval_difference(left_start::Integer, left_stop::Integer, right_start::Integer, right_stop::Integer)
    left_start, left_stop = normalize_interval(left_start, left_stop)
    right_start, right_stop = normalize_interval(right_start, right_stop)
    pieces = Tuple{Int,Int}[]

    (right_stop <= left_start || right_start >= left_stop) && return [(left_start, left_stop)]

    if right_start > left_start
        push!(pieces, (left_start, min(right_start, left_stop)))
    end

    if right_stop < left_stop
        push!(pieces, (max(right_stop, left_start), left_stop))
    end

    return filter(piece -> piece[1] < piece[2], pieces)
end

"""
    interval_difference(left, right)

Subtract one tuple interval from another and return the remaining pieces.
"""
interval_difference(left::Tuple{<:Integer,<:Integer}, right::Tuple{<:Integer,<:Integer}) = interval_difference(left[1], left[2], right[1], right[2])
