using Test
using BioToolkit
using Random
using DataFrames

function _interval_signature(interval)
    return (interval.chrom, interval.left, interval.right, interval.strand, get(interval.metadata, "name", ""))
end

@testset "GenomicRanges Bioconductor-inspired accessors and interval utilities" begin
    intervals = [
        BioToolkit.GenomicInterval("chr1", 1, 10, '+', Dict("id" => 1)),
        BioToolkit.GenomicInterval("chr1", 8, 12, '+', Dict("id" => 2)),
        BioToolkit.GenomicInterval("chr2", 1, 3, '.', Dict("id" => 3)),
    ]
    collection = BioToolkit.GenomicRanges.build_collection(intervals)

    @test BioToolkit.GenomicRanges.start(intervals[1]) == 1
    @test BioToolkit.GenomicRanges.stop(intervals[1]) == 10
    @test BioToolkit.GenomicRanges.end_position(collection) == [10, 12, 3]
    @test BioToolkit.GenomicRanges.seqnames(intervals[1]) == "chr1"
    @test BioToolkit.GenomicRanges.strand(intervals[1]) == '+'
    @test BioToolkit.GenomicRanges.mcols(intervals[1])["id"] == 1
    @test length(BioToolkit.GenomicRanges.subset_by_overlaps(collection, BioToolkit.GenomicRanges.build_collection([BioToolkit.GenomicInterval("chr1", 9, 9)]))) == 2
    @test BioToolkit.GenomicRanges.overlaps_any(intervals, BioToolkit.GenomicRanges.build_collection([BioToolkit.GenomicInterval("chr1", 9, 9)])) == [true, true, false]
    @test BioToolkit.GenomicRanges.restrict(intervals, 3, 9)[1].left == 3
    terminator = BioToolkit.GenomicRanges.terminators(intervals[1:1], 2, 3)[1]
    @test (terminator.chrom, terminator.left, terminator.right, terminator.strand) == ("chr1", 9, 13, '+')
    @test length(BioToolkit.GenomicRanges.tile(intervals[1:1], 3)[1]) == 3
    @test length(BioToolkit.GenomicRanges.sliding_windows(intervals[1:1], 4)[1]) == 7
    @test !BioToolkit.GenomicRanges.is_disjoint(collection)
    @test BioToolkit.GenomicRanges.disjoint_bins(intervals) == [1, 2, 1]
    @test BioToolkit.GenomicRanges.duplicated_intervals(vcat(intervals, intervals[1:1])) == [false, false, false, true]
    @test BioToolkit.GenomicRanges.range_intervals(collection)[1].left == 1
    @test BioToolkit.GenomicRanges.poverlaps(intervals, intervals) == [true, true, true]
    @test length(BioToolkit.GenomicRanges.overlaps_ranges(intervals[1], collection)) == 2
    @test length(BioToolkit.GenomicRanges.find_overlap_pairs(intervals, collection)) == 5
    @test !BioToolkit.GenomicRanges.is_normal(intervals)
    @test length(BioToolkit.GenomicRanges.tile_genome(Dict("chr1" => 10, "chr2" => 5), 3)) == 3
end

_sorted_intervals(intervals) = sort(collect(intervals), by = _interval_signature)

function _brute_force_overlaps(query, intervals)
    result = BioToolkit.GenomicInterval[]
    for interval in intervals
        interval.chrom == query.chrom || continue
        interval.left <= query.right || continue
        interval.right >= query.left || continue
        push!(result, interval)
    end
    return result
end

function _distance(query, interval)
    interval.right < query.left && return query.left - interval.right - 1
    interval.left > query.right && return interval.left - query.right - 1
    return 0
end

function _brute_force_nearest(query, intervals)
    candidates = [interval for interval in intervals if interval.chrom == query.chrom]
    isempty(candidates) && return nothing

    best = nothing
    best_key = nothing

    for interval in candidates
        d = _distance(query, interval)
        key = if interval.left <= query.right && interval.right >= query.left
            (0, 1, interval.left, interval.right, Int(interval.strand))
        elseif interval.right < query.left
            (d, 1, -interval.right, interval.left, interval.right)
        else
            (d, 2, interval.left, interval.right)
        end

        if best_key === nothing || key < best_key
            best = interval
            best_key = key
        end
    end

    return best
end

function _brute_force_follow(query, intervals)
    candidates = [interval for interval in intervals if interval.chrom == query.chrom && interval.right < query.left]
    isempty(candidates) && return nothing
    best = candidates[1]
    best_key = (best.right, best.left, Int(best.strand))
    for interval in candidates[2:end]
        key = (interval.right, interval.left, Int(interval.strand))
        if key > best_key
            best = interval
            best_key = key
        end
    end
    return best
end

function _brute_force_precede(query, intervals)
    candidates = [interval for interval in intervals if interval.chrom == query.chrom && interval.left > query.right]
    isempty(candidates) && return nothing
    best = candidates[1]
    best_key = (best.left, best.right, Int(best.strand))
    for interval in candidates[2:end]
        key = (interval.left, interval.right, Int(interval.strand))
        if key < best_key
            best = interval
            best_key = key
        end
    end
    return best
end

function _brute_force_reduce(intervals)
    isempty(intervals) && return BioToolkit.IntervalCollection(BioToolkit.GenomicInterval[])
    sorted = sort(collect(intervals), by = interval -> (interval.chrom, interval.left, interval.right, interval.strand))
    merged = BioToolkit.GenomicInterval[]
    current = sorted[1]

    for interval in sorted[2:end]
        if interval.chrom == current.chrom && interval.left <= current.right + 1
            current = BioToolkit.GenomicInterval(current.chrom, current.left, max(current.right, interval.right), current.strand, current.metadata)
        else
            push!(merged, current)
            current = interval
        end
    end

    push!(merged, current)
    return BioToolkit.GenomicRanges.build_collection(merged)
end

function _positions_from_intervals(intervals)
    positions = Dict{String,Set{Int}}()
    for interval in intervals
        chrom_positions = get!(positions, interval.chrom, Set{Int}())
        for position in interval.left:interval.right
            push!(chrom_positions, position)
        end
    end
    return positions
end

function _positions_from_coverage(segments)
    positions = Dict{String,Set{Int}}()
    for segment in segments
        chrom_positions = get!(positions, segment.chrom, Set{Int}())
        for position in segment.start:segment.stop
            if segment.depth > 0
                push!(chrom_positions, position)
            end
        end
    end
    return positions
end

function _subtract_positions(a_intervals, b_intervals)
    a_positions = _positions_from_intervals(a_intervals)
    b_positions = _positions_from_intervals(b_intervals)
    result = Dict{String,Set{Int}}()

    for (chrom, positions) in a_positions
        remaining = Set{Int}()
        blocked = get(b_positions, chrom, Set{Int}())
        for position in positions
            position in blocked && continue
            push!(remaining, position)
        end
        result[chrom] = remaining
    end

    return result
end

function _intervals_from_positions(positions::Dict{String,Set{Int}})
    intervals = BioToolkit.GenomicInterval[]
    for chrom in sort(collect(keys(positions)))
        chrom_positions = sort(collect(positions[chrom]))
        isempty(chrom_positions) && continue

        start = chrom_positions[1]
        previous = start
        for position in chrom_positions[2:end]
            if position != previous + 1
                push!(intervals, BioToolkit.GenomicInterval(chrom, start, previous, '.', Dict{String,Any}()))
                start = position
            end
            previous = position
        end
        push!(intervals, BioToolkit.GenomicInterval(chrom, start, previous, '.', Dict{String,Any}()))
    end
    return intervals
end

function _coverage_positions(intervals)
    positions = Dict{String,Set{Int}}()
    for interval in intervals
        for position in interval.left:interval.right
            chrom_positions = get!(positions, interval.chrom, Set{Int}())
            push!(chrom_positions, position)
        end
    end
    return positions
end

@testset "GenomicRanges basics" begin
    intervals = [
        BioToolkit.GenomicInterval("chr2", 50, 60, '+', Dict("name" => "geneD")),
        BioToolkit.GenomicInterval("chr1", 100, 200, '+', Dict("name" => "geneA")),
        BioToolkit.GenomicInterval("chr1", 150, 250, '-', Dict("name" => "geneB")),
        BioToolkit.GenomicInterval("chr1", 300, 400, '.', Dict("name" => "geneC")),
    ]

    collection = BioToolkit.GenomicRanges.build_collection(intervals)
    @test length(collection) == 4
    @test collection[1].metadata["name"] == "geneA"
    @test collection[2].metadata["name"] == "geneB"
    @test collection.chrom_indices["chr1"] == 1:3
    @test collection.chrom_indices["chr2"] == 4:4

    tmp, io = mktemp()
    write(io, "chr,start,end,name,strand\n")
    write(io, "chr1,100,200,geneA,+\n")
    write(io, "chr1,150,250,geneB,-\n")
    write(io, "chr2,50,60,geneD,.\n")
    close(io)
    try
        parsed = BioToolkit.GenomicRanges.read_intervals(tmp)
        @test length(parsed) == 3
        @test parsed[1].metadata["name"] == "geneA"
        @test parsed[2].strand == '-'
        @test parsed[3].chrom == "chr2"
    finally
        rm(tmp, force=true)
    end
end

@testset "GenomicRanges queries" begin
    intervals = [
        BioToolkit.GenomicInterval("chr1", 100, 200, '+', Dict("name" => "geneA")),
        BioToolkit.GenomicInterval("chr1", 150, 250, '-', Dict("name" => "geneB")),
        BioToolkit.GenomicInterval("chr1", 300, 400, '.', Dict("name" => "geneC")),
        BioToolkit.GenomicInterval("chr2", 50, 60, '+', Dict("name" => "geneD")),
    ]
    collection = BioToolkit.GenomicRanges.build_collection(intervals)

    query = BioToolkit.GenomicInterval("chr1", 180, 180, '+', Dict{String,Any}())
    hits = BioToolkit.GenomicRanges.find_overlaps(query, collection)
    @test length(hits) == 2
    @test sort([hit.metadata["name"] for hit in hits]) == ["geneA", "geneB"]
    @test _sorted_intervals(_brute_force_overlaps(query, intervals)) == _sorted_intervals(hits)

    boundary_query = BioToolkit.GenomicInterval("chr1", 200, 200, '+', Dict{String,Any}())
    boundary_hits = BioToolkit.GenomicRanges.overlap(boundary_query, collection)
    @test sort([hit.metadata["name"] for hit in boundary_hits]) == ["geneA", "geneB"]

    nearest_query = BioToolkit.GenomicInterval("chr1", 1000, 1000, '+', Dict{String,Any}())
    nearest_hit = BioToolkit.GenomicRanges.find_nearest(nearest_query, collection)
    @test nearest_hit !== nothing
    @test nearest_hit.metadata["name"] == "geneC"
    @test _brute_force_nearest(nearest_query, intervals) == nearest_hit

    follow_query = BioToolkit.GenomicInterval("chr1", 260, 260, '+', Dict{String,Any}())
    follow_hit = BioToolkit.GenomicRanges.follow(follow_query, collection)
    @test follow_hit !== nothing
    @test follow_hit.metadata["name"] == "geneB"
    @test _brute_force_follow(follow_query, intervals) == follow_hit

    precede_hit = BioToolkit.GenomicRanges.precede(follow_query, collection)
    @test precede_hit !== nothing
    @test precede_hit.metadata["name"] == "geneC"
    @test _brute_force_precede(follow_query, intervals) == precede_hit

    parallel_hits = BioToolkit.GenomicRanges.find_overlaps_parallel(collection, collection)
    @test length(parallel_hits) == length(collection)
    @test all(index -> parallel_hits[index] == BioToolkit.GenomicRanges.find_overlaps(collection[index], collection), eachindex(parallel_hits))

    subtraction_query = BioToolkit.GenomicInterval("chr1", 100, 220, '+', Dict("name" => "peak"))
    subtraction_subject = BioToolkit.GenomicRanges.build_collection([
        BioToolkit.GenomicInterval("chr1", 120, 180, '+', Dict{String,Any}()),
        BioToolkit.GenomicInterval("chr1", 190, 210, '+', Dict{String,Any}()),
    ])
    subtraction = BioToolkit.GenomicRanges.setdiff(subtraction_query, subtraction_subject)
    @test [(piece.left, piece.right) for piece in subtraction] == [(100, 119), (181, 189), (211, 220)]

    left_collection = BioToolkit.GenomicRanges.build_collection([
        BioToolkit.GenomicInterval("chr1", 1, 5, '+', Dict{String,Any}()),
        BioToolkit.GenomicInterval("chr1", 8, 10, '+', Dict{String,Any}()),
    ])
    right_collection = BioToolkit.GenomicRanges.build_collection([
        BioToolkit.GenomicInterval("chr1", 4, 7, '+', Dict{String,Any}()),
        BioToolkit.GenomicInterval("chr1", 9, 12, '+', Dict{String,Any}()),
    ])
    collection_intersection = BioToolkit.GenomicRanges.intersect(left_collection, right_collection)
    @test [(interval.left, interval.right) for interval in collection_intersection] == [(4, 5), (9, 10)]
    collection_union = BioToolkit.GenomicRanges.union(left_collection, right_collection)
    @test [(interval.left, interval.right) for interval in collection_union] == [(1, 12)]

    reduced = BioToolkit.GenomicRanges.reduce(BioToolkit.GenomicRanges.build_collection([
        BioToolkit.GenomicInterval("chr1", 1, 5, '+', Dict("name" => "a")),
        BioToolkit.GenomicInterval("chr1", 6, 8, '+', Dict("name" => "b")),
        BioToolkit.GenomicInterval("chr1", 10, 12, '+', Dict("name" => "c")),
        BioToolkit.GenomicInterval("chr1", 11, 20, '+', Dict("name" => "d")),
        BioToolkit.GenomicInterval("chr2", 1, 2, '+', Dict("name" => "e")),
    ]))
    @test [(interval.chrom, interval.left, interval.right) for interval in reduced] == [("chr1", 1, 8), ("chr1", 10, 20), ("chr2", 1, 2)]

    coverage_segments = BioToolkit.GenomicRanges.coverage(BioToolkit.GenomicRanges.build_collection([
        BioToolkit.GenomicInterval("chr1", 1, 3, '+', Dict{String,Any}()),
        BioToolkit.GenomicInterval("chr1", 2, 4, '+', Dict{String,Any}()),
        BioToolkit.GenomicInterval("chr2", 5, 5, '+', Dict{String,Any}()),
    ]))
    @test [(segment.chrom, segment.start, segment.stop, segment.depth) for segment in coverage_segments] == [("chr1", 1, 1, 1), ("chr1", 2, 3, 2), ("chr1", 4, 4, 1), ("chr2", 5, 5, 1)]
end

@testset "GenomicRanges transforms" begin
    intervals = [
        BioToolkit.GenomicInterval("chr1", 10, 20, '+', Dict("name" => "geneA")),
        BioToolkit.GenomicInterval("chr1", 30, 40, '-', Dict("name" => "geneB")),
    ]

    shifted = BioToolkit.GenomicRanges.shift(intervals, 5)
    @test [(interval.left, interval.right) for interval in shifted] == [(15, 25), (35, 45)]

    flanked_start = BioToolkit.GenomicRanges.flank(intervals, 3)
    flanked_end = BioToolkit.GenomicRanges.flank(intervals, 3; start=false)
    @test [(interval.left, interval.right) for interval in flanked_start] == [(7, 9), (41, 43)]
    @test [(interval.left, interval.right) for interval in flanked_end] == [(21, 23), (27, 29)]

    resized_start = BioToolkit.GenomicRanges.resize(intervals, 4; fix=:start)
    resized_center = BioToolkit.GenomicRanges.resize(intervals, 5; fix=:center)
    resized_end = BioToolkit.GenomicRanges.resize(intervals, 3; fix=:end)
    @test [(interval.left, interval.right) for interval in resized_start] == [(10, 13), (30, 33)]
    @test [(interval.left, interval.right) for interval in resized_center] == [(13, 17), (33, 37)]
    @test [(interval.left, interval.right) for interval in resized_end] == [(18, 20), (38, 40)]

    promoter_regions = BioToolkit.GenomicRanges.promoters(intervals, 5, 2)
    @test [(interval.left, interval.right) for interval in promoter_regions] == [(5, 11), (38, 44)]

    narrowed = BioToolkit.GenomicRanges.narrow(intervals, 2, 3)
    @test [(interval.left, interval.right) for interval in narrowed] == [(12, 17), (32, 37)]

    seqinfo = Dict(
        "chr1" => BioToolkit.GenomicRanges.SeqInfo("chr1", 35, false),
        "chr2" => BioToolkit.GenomicRanges.SeqInfo("chr2", 10, false),
    )
    trimmed = BioToolkit.GenomicRanges.trim([
        BioToolkit.GenomicInterval("chr1", 1, 50, '+', Dict("name" => "a")),
        BioToolkit.GenomicInterval("chr2", 5, 15, '+', Dict("name" => "b")),
    ], seqinfo)
    @test [(interval.chrom, interval.left, interval.right) for interval in trimmed] == [("chr1", 1, 35), ("chr2", 5, 10)]

    gaps = BioToolkit.GenomicRanges.gaps(BioToolkit.GenomicRanges.build_collection([
        BioToolkit.GenomicInterval("chr1", 3, 4, '+', Dict{String,Any}()),
        BioToolkit.GenomicInterval("chr1", 6, 8, '+', Dict{String,Any}()),
    ]), Dict("chr1" => 10, "chr2" => 5))
    @test [(interval.chrom, interval.left, interval.right) for interval in gaps] == [("chr1", 1, 2), ("chr1", 5, 5), ("chr1", 9, 10), ("chr2", 1, 5)]

    disjoined = BioToolkit.GenomicRanges.disjoin(BioToolkit.GenomicRanges.build_collection([
        BioToolkit.GenomicInterval("chr1", 1, 4, '+', Dict{String,Any}()),
        BioToolkit.GenomicInterval("chr1", 3, 5, '+', Dict{String,Any}()),
        BioToolkit.GenomicInterval("chr1", 10, 12, '+', Dict{String,Any}()),
    ]))
    @test [(interval.left, interval.right) for interval in disjoined] == [(1, 2), (3, 4), (5, 5), (10, 12)]

    pairs_left = [
        BioToolkit.GenomicInterval("chr1", 1, 5, '+', Dict("name" => "a")),
        BioToolkit.GenomicInterval("chr1", 10, 15, '+', Dict("name" => "b")),
    ]
    pairs_right = [
        BioToolkit.GenomicInterval("chr1", 3, 7, '+', Dict("name" => "c")),
        BioToolkit.GenomicInterval("chr1", 20, 25, '+', Dict("name" => "d")),
    ]
    pairwise_intersection = BioToolkit.GenomicRanges.pintersect(pairs_left, pairs_right)
    pairwise_union = BioToolkit.GenomicRanges.punion(pairs_left, pairs_right)
    pairwise_difference = BioToolkit.GenomicRanges.psetdiff(pairs_left, pairs_right)
    @test pairwise_intersection[1] !== nothing
    @test pairwise_intersection[1].left == 3 && pairwise_intersection[1].right == 5
    @test pairwise_intersection[2] === nothing
    @test [(interval.left, interval.right) for interval in pairwise_union] == [(1, 7), (10, 25)]
    @test [[(piece.left, piece.right) for piece in pieces] for pieces in pairwise_difference] == [[(1, 2)], [(10, 15)]]

    df = DataFrames.DataFrame(intervals)
    roundtrip = BioToolkit.GenomicRanges.read_intervals(df)
    @test roundtrip == intervals
end

@testset "GenomicRanges transform fuzz" begin
    Random.seed!(20260324 + 1)
    for _ in 1:200
        count = rand(1:15)
        intervals = BioToolkit.GenomicInterval[]
        for _ in 1:count
            left = rand(1:200)
            right = left + rand(0:20)
            strand = rand(('+', '-', '.'))
            push!(intervals, BioToolkit.GenomicInterval("chr1", left, right, strand, Dict{String,Any}()))
        end

        delta = rand(-20:20)
        shifted = BioToolkit.GenomicRanges.shift(intervals, delta)
        @test [(s.left - o.left, s.right - o.right) for (o, s) in zip(intervals, shifted)] == [(delta, delta) for _ in 1:count]

        width = rand(1:20)
        start_flanks = BioToolkit.GenomicRanges.flank(intervals, width)
        end_flanks = BioToolkit.GenomicRanges.flank(intervals, width; start=false)
        @test all(interval -> (interval.right - interval.left + 1) == width, start_flanks)
        @test all(interval -> (interval.right - interval.left + 1) == width, end_flanks)

        if all(interval -> (interval.right - interval.left + 1) >= 3, intervals)
            narrowed = BioToolkit.GenomicRanges.narrow(intervals, 1, 1)
            @test all(pair -> begin
                n, o = pair
                (n.right - n.left + 1) == (o.right - o.left + 1) - 2
            end, zip(narrowed, intervals))
        end

        resized = BioToolkit.GenomicRanges.resize(intervals, 3; fix=:start)
        @test all(interval -> (interval.right - interval.left + 1) == 3, resized)
    end
end

@testset "GenomicRanges fuzz" begin
    Random.seed!(20260324)
    chroms = ["chr1", "chr2", "chrX"]

    for _ in 1:250
        intervals = BioToolkit.GenomicInterval[]
        for _ in 1:40
            chrom = chroms[rand(1:length(chroms))]
            left = rand(1:500)
            right = left + rand(0:40)
            strand = rand(('+', '-', '.'))
            push!(intervals, BioToolkit.GenomicInterval(chrom, left, right, strand, Dict("name" => string(chrom, ":", left, "-", right))))
        end

        collection = BioToolkit.GenomicRanges.build_collection(intervals)

        for _ in 1:20
            chrom = chroms[rand(1:length(chroms))]
            left = rand(1:500)
            right = left + rand(0:40)
            query = BioToolkit.GenomicInterval(chrom, left, right, '+', Dict{String,Any}())

            fast_hits = BioToolkit.GenomicRanges.find_overlaps(query, collection)
            slow_hits = _brute_force_overlaps(query, intervals)
            @test _sorted_intervals(fast_hits) == _sorted_intervals(slow_hits)

            fast_nearest = BioToolkit.GenomicRanges.find_nearest(query, collection)
            slow_nearest = _brute_force_nearest(query, intervals)
            fast_d = fast_nearest === nothing ? -1 : _distance(query, fast_nearest)
            slow_d = slow_nearest === nothing ? -1 : _distance(query, slow_nearest)
            @test fast_d == slow_d

            @test BioToolkit.GenomicRanges.follow(query, collection) == _brute_force_follow(query, intervals)
            @test BioToolkit.GenomicRanges.precede(query, collection) == _brute_force_precede(query, intervals)
        end

        reduced_fast = BioToolkit.GenomicRanges.reduce(collection)
        reduced_slow = _brute_force_reduce(intervals)
        @test [(interval.chrom, interval.left, interval.right) for interval in reduced_fast] == [(interval.chrom, interval.left, interval.right) for interval in reduced_slow]


        coverage_fast = BioToolkit.GenomicRanges.coverage(collection)
        fast_positions = _positions_from_coverage(coverage_fast)
        slow_positions = _coverage_positions(intervals)
        @test fast_positions == slow_positions
    end
end

@testset "GenomicRanges indexed bulk queries" begin
    empty_collection = BioToolkit.GenomicRanges.build_collection(BioToolkit.GenomicInterval[])
    @test isempty(empty_collection)
    @test BioToolkit.GenomicRanges.find_overlaps(
        BioToolkit.GenomicInterval("chr1", 1, 2), empty_collection) == BioToolkit.GenomicInterval[]
    @test BioToolkit.GenomicRanges.count_overlaps(
        BioToolkit.GenomicInterval("chr1", 1, 2), empty_collection) == 0

    subject = BioToolkit.GenomicRanges.build_collection([
        BioToolkit.GenomicInterval("chr1", 10, 20, '+', Dict("id" => 1)),
        BioToolkit.GenomicInterval("chr1", 15, 25, '-', Dict("id" => 2)),
        BioToolkit.GenomicInterval("chr1", 40, 50, '.', Dict("id" => 3)),
        BioToolkit.GenomicInterval("chr2", 1, 5, '.', Dict("id" => 4)),
    ])
    queries = [
        BioToolkit.GenomicInterval("chr1", 18, 18, '+'),
        BioToolkit.GenomicInterval("chr1", 30, 35, '+'),
        BioToolkit.GenomicInterval("chr2", 3, 3, '.'),
    ]

    bulk_hits = BioToolkit.GenomicRanges.find_overlaps(queries, subject)
    @test [sort(get.(getfield.(hits, :metadata), "id", 0)) for hits in bulk_hits] == [[1, 2], [], [4]]
    @test BioToolkit.GenomicRanges.count_overlaps(queries, subject) == [2, 0, 1]
    @test collect(BioToolkit.GenomicRanges.eachoverlap(queries[1], subject)) == bulk_hits[1]
    @test BioToolkit.GenomicRanges.hasintersection(queries[1], subject)
    @test !BioToolkit.GenomicRanges.hasintersection(queries[2], subject)
    @test length(collect(BioToolkit.GenomicRanges.eachoverlap(queries, subject))) == 3

    query_collection = BioToolkit.GenomicRanges.build_collection(queries)
    @test BioToolkit.GenomicRanges.find_overlaps(query_collection, subject) == bulk_hits
    @test BioToolkit.GenomicRanges.count_overlaps(query_collection, subject) == [2, 0, 1]

    tied_subject = BioToolkit.GenomicRanges.build_collection([
        BioToolkit.GenomicInterval("chr1", 10, 10),
        BioToolkit.GenomicInterval("chr1", 30, 30),
    ])
    tied = BioToolkit.GenomicRanges.nearest(BioToolkit.GenomicInterval("chr1", 20, 20), tied_subject; select=:all)
    @test length(tied) == 2
    @test BioToolkit.GenomicRanges.nearest(BioToolkit.GenomicInterval("chr1", 20, 20), tied_subject) isa BioToolkit.GenomicInterval
end

@testset "GenomicRanges storage and missing-row safety" begin
    empty_metadata_interval = BioToolkit.GenomicInterval("chr1", 1, 2)
    @test typeof(empty_metadata_interval) == BioToolkit.GenomicInterval{String,Int,Int,Char,NamedTuple{(),Tuple{}}}
    collection = BioToolkit.GenomicRanges.build_collection([empty_metadata_interval])
    @test eltype(collection.intervals) == BioToolkit.GenomicIntervalEmpty

    table = DataFrame(chrom=["chr1", missing, "chr2"], start=[1, 2, missing], stop=[3, 4, 8])
    parsed = BioToolkit.GenomicRanges.read_intervals(table)
    @test length(parsed) == 1
    @test parsed[1].chrom == "chr1"
end

@testset "Bioconductor API gap coverage" begin
    table = DataFrame(
        seq=["chr1", "chr1", missing],
        s=[0, 10, 20],
        e=[5, 15, 25],
        strand_col=['+', '-', '+'],
        feature=["a", "b", "ignored"])
    intervals = BioToolkit.make_granges_from_dataframe(
        table; chrom_col=:seq, start_col=:s, end_col=:e, strand_col=:strand_col,
        metadata_cols=[:feature], zero_based=true)
    @test length(intervals) == 2
    @test intervals[1].left == 1
    @test intervals[2].right == 15

    grouped = BioToolkit.make_granges_list_from_dataframe(
        DataFrame(group=["a", "a", "b"], chrom=["chr1", "chr1", "chr2"], start=[1, 5, 2], stop=[3, 8, 4]); group_col=:group)
    @test sort(collect(keys(grouped))) == ["a", "b"]
    @test length(grouped["a"]) == 2

    @test BioToolkit.which_as_iranges([false, true, true, false, true]) == [BioToolkit.IRange(2, 2), BioToolkit.IRange(5, 1)]
    @test BioToolkit.as_normal_iranges([BioToolkit.IRange(1, 3), BioToolkit.IRange(4, 2)]) == [BioToolkit.IRange(1, 5)]
    @test BioToolkit.break_in_chunks(10, 3) == [1:4, 5:8, 9:10]
    @test BioToolkit.is_small_genome(Dict("chr1" => 10, "chr2" => 20))
end

@testset "Expanded Bioconductor parity features" begin
    # GPos and IPos tests
    gp = BioToolkit.GPos("chr1", 100, '+')
    @test BioToolkit.is_gpos(gp)
    @test gp.left == 100 && gp.right == 100
    gps = BioToolkit.GPos("chr1", [10, 20, 30])
    @test length(gps) == 3 && all(BioToolkit.is_gpos, gps)

    ip = BioToolkit.IPos(50)
    @test BioToolkit.is_ipos(ip)
    @test ip.start == 50 && ip.width == 1
    ips = BioToolkit.IPos([1, 2, 3])
    @test length(ips) == 3 && all(BioToolkit.is_ipos, ips)

    # zoom tests
    iv = BioToolkit.GenomicInterval("chr1", 10, 20, '+', Dict{String,Any}())
    z_out = BioToolkit.zoom(iv, 2.0)
    @test BioToolkit.GenomicRanges.width(z_out) == 22
    z_in = BioToolkit.zoom_in(iv, 2.0)
    @test BioToolkit.GenomicRanges.width(z_in) == 6

    # seqlevels style mapping
    styles = ["chr1", "chrX", "chrMT"]
    mapped_ens = BioToolkit.map_seqlevels_style(styles, :Ensembl)
    @test mapped_ens == ["1", "X", "MT"]
    mapped_ucsc = BioToolkit.map_seqlevels_style(mapped_ens, :UCSC)
    @test mapped_ucsc == ["chr1", "chrX", "chrMT"]

    # coverage_by_transcript
    exons_tx1 = BioToolkit.build_collection([
        BioToolkit.GenomicInterval("chr1", 1, 5, '+'),
        BioToolkit.GenomicInterval("chr1", 10, 12, '+')
    ])
    transcripts = Dict("tx1" => exons_tx1)
    segs = [BioToolkit.CoverageSegment("chr1", 1, 15, 2)]
    tx_cov = BioToolkit.coverage_by_transcript(segs, transcripts)
    @test tx_cov["tx1"] == [2, 2, 2, 2, 2, 2, 2, 2] # 5 + 3 = 8 bases

    # extract_transcript_seqs and extract_upstream_seqs
    genome = Dict("chr1" => "ATCGATCGATCGATCGATCG")
    tx_seqs = BioToolkit.extract_transcript_seqs(genome, transcripts)
    @test tx_seqs["tx1"] == "ATCGATCG"

    up_seqs = BioToolkit.extract_upstream_seqs(genome, transcripts; width=3)
    @test up_seqs["tx1"] == ""
    
    exons_tx2 = BioToolkit.build_collection([
        BioToolkit.GenomicInterval("chr1", 10, 15, '-')
    ])
    transcripts_minus = Dict("tx2" => exons_tx2)
    tx_seqs_minus = BioToolkit.extract_transcript_seqs(genome, transcripts_minus)
    @test tx_seqs_minus["tx2"] == BioToolkit.GenomicRanges._reverse_complement_string("TCGATC")

    # extend_exons_into_introns & element_lengths
    ext = BioToolkit.extend_exons_into_introns(transcripts, 2)
    @test ext["tx1"][1].right == 7
    @test BioToolkit.element_lengths(transcripts) == Dict("tx1" => 2)

    # vcat for IntervalCollection
    c1 = BioToolkit.build_collection([BioToolkit.GenomicInterval("chr1", 1, 5)])
    c2 = BioToolkit.build_collection([BioToolkit.GenomicInterval("chr2", 10, 15)])
    c_cat = vcat(c1, c2)
    @test length(c_cat) == 2
end
