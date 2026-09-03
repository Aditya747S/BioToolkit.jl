# ==============================================================================
# genomicranges.jl — Genomic interval operations
#
# Provides an IntervalCollection backed by sorted arrays with chromosome
# partitioning for efficient overlap queries. The find_overlaps function
# uses binary search to find candidate intervals, giving O(log n + k)
# performance where k is the number of hits.
#
# The setdiff operation uses O(n log n) sweep-line interval subtraction
# rather than position-level enumeration, making it efficient even for
# megabase-scale genomic intervals.
#
# References:
#   - Interval overlap queries: based on UCSC binning scheme concepts
#   - Sweep-line interval arithmetic: Bentley & Ottmann (1979)
# ==============================================================================

module GenomicRanges

import ..BioToolkit: normalize_interval, BedRecord, GffRecord
using PooledArrays
using DataFrames
using Statistics
using ..BioToolkit: PROVENANCE_ID_KEY, PROVENANCE_METADATA_KEY, ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, metadata_provenance, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!, stamp_provenance!, with_provenance

@inline function _register_genomicranges_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

const PooledStringVector = PooledVector{String,UInt32,Vector{UInt32}}

export GenomicInterval, GenomicIntervalDict, GenomicIntervalEmpty, IntervalCollection, IRange, IRanges, CoverageSegment, SeqInfo
export build_collection, read_intervals
export overlap, find_overlaps, eachoverlap, hasintersection
export nearest, nearest_first, nearest_all, find_nearest, distance, follow, precede, find_overlaps_parallel
export shift, flank, resize, promoters, narrow
export trim, gaps, complement, disjoin, pintersect, punion, psetdiff
export coverage, count_overlaps
export OverlapProfileResult, profile_interval_overlaps
export width, mid, start, stop, end_position, seqnames, strand, ranges, mcols, metadata, seqinfo
export subset_by_overlaps, overlaps_any, restrict, terminators, tile, sliding_windows
export tile_genome, is_disjoint, disjoint_bins, order_intervals, duplicated_intervals
export range_intervals, poverlaps, overlaps_ranges, find_overlap_pairs, merge_by_overlaps, is_normal
export successive_iranges, sliding_iranges, centered_iranges, reflect, threebands
export seqlevels, seqlengths, isCircular, genome, keep_seqlevels, drop_seqlevels, rename_seqlevels, standard_chromosomes, seqlevels_in_use
export order_seqlevels, rank_seqlevels, sort_seqlevels, extract_seqlevels, extract_seqlevels_by_group, map_seqlevels, seqlevels_in_group, keep_standard_chromosomes, seqlevels_style
export orderSeqlevels, rankSeqlevels, sortSeqlevels, extractSeqlevels, extractSeqlevelsByGroup, mapSeqlevels, seqlevelsInGroup, keepStandardChromosomes, seqlevelsStyle
export make_granges_from_dataframe, make_granges_list_from_dataframe, make_granges_list_from_feature_fragments
export makeGRangesFromDataFrame, makeGRangesListFromDataFrame, makeGRangesListFromFeatureFragments
export make_iranges_from_dataframe, granges, grglist, subtract, is_small_genome, absolute_ranges, relative_ranges
export isSmallGenome, absoluteRanges, relativeRanges, tileGenome
export which_as_iranges, as_normal_iranges, break_in_chunks, heads, tails, distance_to_nearest
export features, transcripts, exons, cds, genes, transcripts_by, exons_by, cds_by
export transcripts_by_overlaps, exons_by_overlaps, cds_by_overlaps, introns_by_transcript, transcript_widths
export five_utrs_by_transcript, three_utrs_by_transcript, transcript_lengths, tidy_transcripts, tidy_exons, tidy_introns, exonic_parts, intronic_parts, id2name
export fiveUTRsByTranscript, threeUTRsByTranscript, transcriptLengths, tidyTranscripts, tidyExons, tidyIntrons, exonicParts, intronicParts
export GPos, IPos, is_gpos, is_ipos
export coverage_by_transcript, coverageByTranscript
export extract_transcript_seqs, extractTranscriptSeqs, extract_upstream_seqs, extractUpstreamSeqs
export extend_exons_into_introns, extendExonsIntoIntrons
export map_seqlevels_style, mapSeqlevelsStyle
export zoom, zoom_in, zoom_out, element_lengths, elementLengths

struct GenomicInterval{C<:AbstractString,L<:Integer,R<:Integer,S<:Char,M}
    chrom::C
    left::L
    right::R
    strand::S
    metadata::M

    function GenomicInterval(chrom::C, left::L, right::R, strand::S, metadata::M) where {C<:AbstractString,L<:Integer,R<:Integer,S<:Char,M}
        strand in ('+', '-', '.') || throw(ArgumentError("strand must be '+', '-', or '.'"))
        left_int, right_int = normalize_interval(left, right)
        metadata_copy = metadata isa AbstractDict ? Dict{String,Any}(string(key) => value for (key, value) in pairs(metadata)) : metadata
        if metadata_copy isa AbstractDict && metadata_provenance(metadata_copy) === nothing
            stamp_provenance!(metadata_copy; label="GenomicInterval", source="GenomicRanges/GenomicInterval", notes=["constructed genomic interval"], parameters=(chrom=String(chrom), left=left_int, right=right_int, strand=strand))
        end
        return new{C,Int,Int,S,typeof(metadata_copy)}(chrom, left_int, right_int, strand, metadata_copy)
    end
end

const GenomicIntervalDict = GenomicInterval{String,Int,Int,Char,Dict{String,Any}}
const GenomicIntervalEmpty = GenomicInterval{String,Int,Int,Char,NamedTuple{(),Tuple{}}}

"""A dependency-free IRanges equivalent using one-based closed coordinates."""
struct IRange
    start::Int
    width::Int
    function IRange(start::Integer, width::Integer)
        width >= 0 || throw(ArgumentError("width must be nonnegative"))
        new(Int(start), Int(width))
    end
end

IRange(; start::Integer=1, end_position::Union{Nothing,Integer}=nothing, width::Union{Nothing,Integer}=nothing) =
    end_position === nothing && width === nothing ? throw(ArgumentError("provide end_position or width")) :
    width === nothing ? IRange(start, Int(end_position) - Int(start) + 1) : IRange(start, width)

IRanges(starts::AbstractVector{<:Integer}, ends::AbstractVector{<:Integer}) = begin
    length(starts) == length(ends) || throw(ArgumentError("starts and ends must have equal length"))
    [IRange(first, last - first + 1) for (first, last) in zip(starts, ends)]
end
IRanges(; start=nothing, end_position=nothing, width=nothing) = begin
    start === nothing && end_position === nothing && width === nothing && return IRange[]
    start === nothing && throw(ArgumentError("start is required"))
    width === nothing && end_position === nothing && throw(ArgumentError("provide end_position or width"))
    starts = start isa AbstractVector ? collect(Int, start) : [Int(start)]
    values = width === nothing ? (end_position isa AbstractVector ? collect(Int, end_position) : [Int(end_position)]) :
             (width isa AbstractVector ? collect(Int, width) : [Int(width)])
    if width === nothing
        return IRanges(starts, values)
    end
    isempty(starts) && return IRange[]
    length(starts) == length(values) || (length(starts) == 1 || length(values) == 1) ||
        throw(ArgumentError("start and width must have compatible lengths"))
    [IRange(starts[mod1(i, length(starts))], values[mod1(i, length(values))]) for i in 1:max(length(starts), length(values))]
end

struct IntervalCollection{I<:GenomicInterval}
    chroms::PooledStringVector
    starts::Vector{Int}
    ends::Vector{Int}
    strands::Vector{Char}
    intervals::Vector{I}
    chrom_indices::Dict{String,UnitRange{Int}}
    chrom_end_indices::Dict{String,Vector{Int}}
    prefix_max_ends::Dict{String,Vector{Int}}
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

"""
    OverlapProfileResult

Reproducible timing and allocation measurements comparing the current indexed
overlap query with a direct coordinate scan. This is intentionally a profiling
baseline, not an alternative production index: it supplies the evidence needed
before replacing the interval representation or memory layout.
"""
struct OverlapProfileResult
    query_count::Int
    subject_count::Int
    repetitions::Int
    indexed_seconds::Float64
    naive_seconds::Float64
    indexed_allocated_bytes::Int
    naive_allocated_bytes::Int
    indexed_hit_count::Int
    naive_hit_count::Int
end

GenomicInterval(chrom::String, left::Integer, right::Integer) = GenomicInterval(chrom, left, right, '.', NamedTuple())
GenomicInterval(chrom::String, left::Integer, right::Integer, strand::Char) = GenomicInterval(chrom, left, right, strand, NamedTuple())

GenomicInterval(record::BedRecord) = GenomicInterval(record.chrom, Int(record.start) + 1, Int(record.stop), '.', NamedTuple())

function GenomicInterval(record::GffRecord)
    strand = isempty(record.strand) ? '.' : first(record.strand)
    strand in ('+', '-', '.') || (strand = '.')
    return GenomicInterval(record.chrom, Int(record.start), Int(record.stop), strand, record.attribute_map)
end

function _comparable_interval_metadata(metadata)
    ignored = Set((String(PROVENANCE_METADATA_KEY), String(PROVENANCE_ID_KEY)))
    return Dict{String,Any}(String(key) => value for (key, value) in pairs(metadata) if !(String(key) in ignored))
end

function Base.:(==)(left::GenomicInterval, right::GenomicInterval)
    return isequal(left.chrom, right.chrom) &&
           isequal(left.left, right.left) &&
           isequal(left.right, right.right) &&
           isequal(left.strand, right.strand) &&
           isequal(_comparable_interval_metadata(left.metadata), _comparable_interval_metadata(right.metadata))
end

function Base.hash(interval::GenomicInterval, state::UInt)
    state = hash(interval.chrom, state)
    state = hash(interval.left, state)
    state = hash(interval.right, state)
    state = hash(interval.strand, state)
    cm = _comparable_interval_metadata(interval.metadata)
    for key in sort!(collect(keys(cm)))
        state = hash(key, state)
        state = hash(cm[key], state)
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

"""
    width(interval) / mid(interval)

Return the width or midpoint of a genomic interval. Bioconductor parity.
"""
width(interval::GenomicInterval) = interval.right - interval.left + 1
mid(interval::GenomicInterval) = (interval.left + interval.right) ÷ 2
width(interval::IRange) = interval.width
mid(interval::IRange) = interval.start + (interval.width - 1) ÷ 2
width(intervals::AbstractVector{<:GenomicInterval}) = [width(i) for i in intervals]
mid(intervals::AbstractVector{<:GenomicInterval}) = [mid(i) for i in intervals]
start(interval::IRange) = interval.start
stop(interval::IRange) = interval.start + interval.width - 1
end_position(interval::IRange) = stop(interval)
start(intervals::AbstractVector{<:IRange}) = [x.start for x in intervals]
stop(intervals::AbstractVector{<:IRange}) = [stop(x) for x in intervals]
end_position(intervals::AbstractVector{<:IRange}) = stop(intervals)
width(intervals::AbstractVector{<:IRange}) = [width(x) for x in intervals]
mid(intervals::AbstractVector{<:IRange}) = [mid(x) for x in intervals]
width(collection::IntervalCollection) = width(collection.intervals)
mid(collection::IntervalCollection) = mid(collection.intervals)

# Bioconductor-compatible accessors. `end` is a Julia keyword, so the explicit
# `end_position` spelling is used for the right endpoint.
start(interval::GenomicInterval) = interval.left
stop(interval::GenomicInterval) = interval.right
end_position(interval::GenomicInterval) = interval.right
seqnames(interval::GenomicInterval) = interval.chrom
strand(interval::GenomicInterval) = interval.strand
ranges(interval::GenomicInterval) = (start=interval.left, stop=interval.right, width=width(interval))
mcols(interval::GenomicInterval) = interval.metadata
metadata(interval::GenomicInterval) = interval.metadata

start(intervals::AbstractVector{<:GenomicInterval}) = [interval.left for interval in intervals]
stop(intervals::AbstractVector{<:GenomicInterval}) = [interval.right for interval in intervals]
end_position(intervals::AbstractVector{<:GenomicInterval}) = stop(intervals)
seqnames(intervals::AbstractVector{<:GenomicInterval}) = [interval.chrom for interval in intervals]
strand(intervals::AbstractVector{<:GenomicInterval}) = [interval.strand for interval in intervals]
ranges(intervals::AbstractVector{<:GenomicInterval}) = [(start=interval.left, stop=interval.right, width=width(interval)) for interval in intervals]
mcols(intervals::AbstractVector{<:GenomicInterval}) = [interval.metadata for interval in intervals]
metadata(intervals::AbstractVector{<:GenomicInterval}) = mcols(intervals)

"""Construct consecutive ranges from widths and an optional gap."""
function successive_iranges(widths::AbstractVector{<:Integer}; gapwidth::Integer=0, from::Integer=1)
    gapwidth >= 0 || throw(ArgumentError("gapwidth must be nonnegative"))
    cursor = Int(from)
    result = IRange[]
    for value in widths
        value >= 0 || throw(ArgumentError("widths must be nonnegative"))
        push!(result, IRange(cursor, value))
        cursor += Int(value) + Int(gapwidth)
    end
    result
end

"""Construct all fixed-width windows in a sequence of length `len`."""
function sliding_iranges(len::Integer, window_width::Integer; shift::Integer=1)
    len >= 0 || throw(ArgumentError("len must be nonnegative"))
    window_width > 0 || throw(ArgumentError("window_width must be positive"))
    shift > 0 || throw(ArgumentError("shift must be positive"))
    len < window_width && return IRange[]
    [IRange(position, window_width) for position in 1:Int(shift):(Int(len) - Int(window_width) + 1)]
end

centered_iranges(center::Integer, flank::Integer) = begin
    flank >= 0 || throw(ArgumentError("flank must be nonnegative"))
    IRange(Int(center) - Int(flank), 2Int(flank) + 1)
end
centered_iranges(centers::AbstractVector{<:Integer}, flank::Integer) = [centered_iranges(center, flank) for center in centers]

"""Reflect ranges within an inclusive reference interval."""
function reflect(intervals::AbstractVector{<:IRange}, lower::Integer, upper::Integer)
    lower <= upper || throw(ArgumentError("lower must not exceed upper"))
    [IRange(Int(lower) + Int(upper) - stop(x), width(x)) for x in intervals]
end

"""Return left, interior, and right pieces around an interval."""
function threebands(x::IRange, lower::Integer, upper::Integer)
    lower <= upper || throw(ArgumentError("lower must not exceed upper"))
    left = IRange(start(x), max(0, Int(lower) - start(x)))
    middle_start = max(start(x), Int(lower))
    middle_end = min(stop(x), Int(upper))
    middle = middle_start <= middle_end ? IRange(middle_start, middle_end - middle_start + 1) : IRange(middle_start, 0)
    right_start = min(stop(x) + 1, Int(upper) + 1)
    right = IRange(right_start, max(0, stop(x) - right_start + 1))
    (left=left, middle=middle, right=right)
end

start(collection::IntervalCollection) = collection.starts
stop(collection::IntervalCollection) = collection.ends
end_position(collection::IntervalCollection) = collection.ends
seqnames(collection::IntervalCollection) = String.(collection.chroms)
strand(collection::IntervalCollection) = collection.strands
ranges(collection::IntervalCollection) = ranges(collection.intervals)
mcols(collection::IntervalCollection) = [interval.metadata for interval in collection.intervals]
metadata(collection::IntervalCollection) = mcols(collection)

"Return chromosome lengths and circularity metadata as a compact dictionary."
seqinfo(collection::IntervalCollection) = Dict(
    chrom => SeqInfo(chrom, maximum(collection.ends[collection.chrom_indices[chrom]]), false)
    for chrom in keys(collection.chrom_indices))
seqinfo(intervals::AbstractVector{<:GenomicInterval}) = seqinfo(build_collection(intervals))

seqlevels(collection::IntervalCollection) = sort!(collect(keys(collection.chrom_indices)))
seqlevels(intervals::AbstractVector{<:GenomicInterval}) = sort!(unique(String[interval.chrom for interval in intervals]))
seqlengths(collection::IntervalCollection) = Dict(chrom => info.length for (chrom, info) in seqinfo(collection))
seqlengths(intervals::AbstractVector{<:GenomicInterval}) = seqlengths(build_collection(intervals))
isCircular(collection::IntervalCollection) = Dict(chrom => info.is_circular for (chrom, info) in seqinfo(collection))
isCircular(intervals::AbstractVector{<:GenomicInterval}) = isCircular(build_collection(intervals))
genome(collection::IntervalCollection) = nothing
genome(intervals::AbstractVector{<:GenomicInterval}) = nothing

function keep_seqlevels(collection::IntervalCollection, levels)
    wanted = Set(String[string(level) for level in levels])
    return build_collection([interval for interval in collection.intervals if interval.chrom in wanted])
end
keep_seqlevels(intervals::AbstractVector{<:GenomicInterval}, levels) =
    [interval for interval in intervals if interval.chrom in Set(String[string(level) for level in levels])]
drop_seqlevels(collection::IntervalCollection, levels) = begin
    unwanted = Set(String[string(level) for level in levels])
    build_collection([interval for interval in collection.intervals if !(interval.chrom in unwanted)])
end
drop_seqlevels(intervals::AbstractVector{<:GenomicInterval}, levels) = begin
    unwanted = Set(String[string(level) for level in levels])
    [interval for interval in intervals if !(interval.chrom in unwanted)]
end
function rename_seqlevels(collection::IntervalCollection, mapping::AbstractDict)
    build_collection([GenomicInterval(String(get(mapping, interval.chrom, interval.chrom)), interval.left, interval.right, interval.strand, interval.metadata) for interval in collection.intervals])
end
rename_seqlevels(intervals::AbstractVector{<:GenomicInterval}, mapping::AbstractDict) =
    [GenomicInterval(String(get(mapping, interval.chrom, interval.chrom)), interval.left, interval.right, interval.strand, interval.metadata) for interval in intervals]
seqlevels_in_use(x) = seqlevels(x)
standard_chromosomes(x; prefixes=("chr",)) = [chrom for chrom in seqlevels(x) if any(startswith(chrom, prefix) for prefix in prefixes) && !occursin("_", chrom)]
keep_standard_chromosomes(x; prefixes=("chr",)) = keep_seqlevels(x, standard_chromosomes(x; prefixes=prefixes))

function _seqlevel_sort_key(level::AbstractString)
    value = replace(String(level), r"^chr"i => "")
    numeric = tryparse(Int, value)
    numeric === nothing ? (1, value == "X" ? 23 : value == "Y" ? 24 : value == "M" || value == "MT" ? 25 : 26, value) : (0, numeric, "")
end
order_seqlevels(levels) = sortperm(String[string(level) for level in levels], by=_seqlevel_sort_key)
rank_seqlevels(levels) = begin
    order = order_seqlevels(levels)
    ranks = zeros(Int, length(order))
    for (rank, index) in enumerate(order)
        ranks[index] = rank
    end
    ranks
end
sort_seqlevels(levels; rev::Bool=false) = sort(String[string(level) for level in levels], by=_seqlevel_sort_key, rev=rev)
extract_seqlevels(x) = seqlevels(x)
extract_seqlevels_by_group(x, group) = [level for level in seqlevels(x) if group === nothing || occursin(lowercase(String(group)), lowercase(level))]
seqlevels_in_group(levels, group) = extract_seqlevels_by_group(levels, group)
map_seqlevels(levels, mapping::AbstractDict) = [String(get(mapping, level, level)) for level in levels]
seqlevels_style(levels::Union{AbstractVector, AbstractSet, Tuple}) = all(startswith(String(level), "chr") for level in levels) ? :UCSC : :Ensembl
seqlevels_style(x) = seqlevels_style(seqlevels(x))

# Bioconductor spelling aliases. The implementation remains shared with the
# idiomatic Julia names above, so the aliases do not create a second code path.
const orderSeqlevels = order_seqlevels
const rankSeqlevels = rank_seqlevels
const sortSeqlevels = sort_seqlevels
const extractSeqlevels = extract_seqlevels
const extractSeqlevelsByGroup = extract_seqlevels_by_group
const mapSeqlevels = map_seqlevels
const seqlevelsInGroup = seqlevels_in_group
const keepStandardChromosomes = keep_standard_chromosomes
const seqlevelsStyle = seqlevels_style

"""Construct intervals from arbitrary DataFrame column names."""
function make_granges_from_dataframe(
    df::DataFrames.AbstractDataFrame;
    chrom_col=:chrom,
    start_col=:start,
    end_col=:stop,
    strand_col=nothing,
    metadata_cols=nothing,
    zero_based::Bool=false)
    chrom_col = Symbol(chrom_col)
    start_col = Symbol(start_col)
    end_col = Symbol(end_col)
    strand_col = strand_col === nothing ? nothing : Symbol(strand_col)
    names_set = Set(Symbol.(DataFrames.names(df)))
    for column in (chrom_col, start_col, end_col)
        column in names_set || throw(ArgumentError("missing required interval column: $(column)"))
    end
    selected = Dict{Symbol,Any}(
        :chrom => df[!, chrom_col],
        :start => zero_based ? [ismissing(value) ? missing : Int(value) + 1 for value in df[!, start_col]] : df[!, start_col],
        :stop => df[!, end_col])
    strand_col !== nothing && (selected[:strand] = df[!, strand_col])
    metadata_names = metadata_cols === nothing ?
        [column for column in DataFrames.names(df) if column ∉ (chrom_col, start_col, end_col, strand_col)] :
        (metadata_cols isa Symbol || metadata_cols isa AbstractString ? [Symbol(metadata_cols)] : Symbol[Symbol(column) for column in metadata_cols])
    for column in metadata_names
        selected[Symbol(column)] = df[!, column]
    end
    return read_intervals(DataFrames.DataFrame(selected))
end

"""Group a feature table into independent interval collections."""
function make_granges_list_from_dataframe(df::DataFrames.AbstractDataFrame; group_col=:group, kwargs...)
    group_col = Symbol(group_col)
    group_col in Set(Symbol.(DataFrames.names(df))) || throw(ArgumentError("missing group column: $(group_col)"))
    result = Dict{Any,IntervalCollection}()
    for group in unique(skipmissing(df[!, group_col]))
        row_indices = findall(value -> !ismissing(value) && value == group, df[!, group_col])
        rows = df[row_indices, :]
        result[group] = build_collection(make_granges_from_dataframe(rows; kwargs...))
    end
    return result
end

make_granges_list_from_feature_fragments(df::DataFrames.AbstractDataFrame; feature_col=:feature_id, kwargs...) =
    make_granges_list_from_dataframe(df; group_col=feature_col, kwargs...)

const makeGRangesFromDataFrame = make_granges_from_dataframe
const makeGRangesListFromDataFrame = make_granges_list_from_dataframe
const makeGRangesListFromFeatureFragments = make_granges_list_from_feature_fragments

function _feature_type_column(df::DataFrames.AbstractDataFrame, type_col)
    type_col in Set(DataFrames.names(df)) && return type_col
    for candidate in (:type, :feature, :feature_type)
        candidate in Set(DataFrames.names(df)) && return candidate
    end
    throw(ArgumentError("feature table requires a type/feature column"))
end

features(df::DataFrames.AbstractDataFrame; kwargs...) = make_granges_from_dataframe(df; kwargs...)

function _feature_subset(df::DataFrames.AbstractDataFrame, kind::AbstractString; type_col=:type, kwargs...)
    column = _feature_type_column(df, type_col)
    indices = findall(value -> !ismissing(value) && lowercase(String(value)) == lowercase(kind), df[!, column])
    make_granges_from_dataframe(df[indices, :]; kwargs...)
end

function _feature_dataframe(df::DataFrames.AbstractDataFrame, kind::AbstractString; type_col=:type)
    column = _feature_type_column(df, type_col)
    indices = findall(value -> !ismissing(value) && lowercase(String(value)) == lowercase(kind), df[!, column])
    df[indices, :]
end

transcripts(df::DataFrames.AbstractDataFrame; type_col=:type, kwargs...) = _feature_subset(df, "transcript"; type_col=type_col, kwargs...)
exons(df::DataFrames.AbstractDataFrame; type_col=:type, kwargs...) = _feature_subset(df, "exon"; type_col=type_col, kwargs...)
cds(df::DataFrames.AbstractDataFrame; type_col=:type, kwargs...) = _feature_subset(df, "cds"; type_col=type_col, kwargs...)
genes(df::DataFrames.AbstractDataFrame; type_col=:type, kwargs...) = _feature_subset(df, "gene"; type_col=type_col, kwargs...)

transcripts_by(df::DataFrames.AbstractDataFrame; group_col=:gene_id, type_col=:type, kwargs...) =
    make_granges_list_from_dataframe(_feature_dataframe(df, "transcript"; type_col=type_col); group_col=group_col, kwargs...)
exons_by(df::DataFrames.AbstractDataFrame; group_col=:transcript_id, type_col=:type, kwargs...) =
    make_granges_list_from_dataframe(_feature_dataframe(df, "exon"; type_col=type_col); group_col=group_col, kwargs...)
cds_by(df::DataFrames.AbstractDataFrame; group_col=:transcript_id, type_col=:type, kwargs...) =
    make_granges_list_from_dataframe(_feature_dataframe(df, "cds"; type_col=type_col); group_col=group_col, kwargs...)

function _features_by_overlaps(query, df::DataFrames.AbstractDataFrame, kind; type_col=:type, kwargs...)
    values = kind === nothing ? make_granges_from_dataframe(df; kwargs...) : _feature_subset(df, kind; type_col=type_col, kwargs...)
    query_intervals = query isa IntervalCollection ? query.intervals : query isa GenomicInterval ? [query] : collect(query)
    mask = [any(query.chrom == value.chrom && query.left <= value.right && query.right >= value.left
                for query in query_intervals) for value in values]
    return values[mask]
end

transcripts_by_overlaps(query, df::DataFrames.AbstractDataFrame; type_col=:type, kwargs...) = _features_by_overlaps(query, df, "transcript"; type_col=type_col, kwargs...)
exons_by_overlaps(query, df::DataFrames.AbstractDataFrame; type_col=:type, kwargs...) = _features_by_overlaps(query, df, "exon"; type_col=type_col, kwargs...)
cds_by_overlaps(query, df::DataFrames.AbstractDataFrame; type_col=:type, kwargs...) = _features_by_overlaps(query, df, "cds"; type_col=type_col, kwargs...)

function introns_by_transcript(df::DataFrames.AbstractDataFrame; transcript_col=:transcript_id, kwargs...)
    grouped = exons_by(df; group_col=transcript_col, kwargs...)
    result = Dict{Any,Vector{GenomicInterval}}()
    for (transcript, collection) in grouped
        introns = GenomicInterval[]
        for chrom in seqlevels(collection)
            values = collection.intervals[collection.chrom_indices[chrom]]
            for index in 2:length(values)
                values[index - 1].right + 1 <= values[index].left - 1 &&
                    push!(introns, GenomicInterval(chrom, values[index - 1].right + 1, values[index].left - 1, values[index].strand))
            end
        end
        result[transcript] = introns
    end
    result
end

function transcript_widths(df::DataFrames.AbstractDataFrame; transcript_col=:transcript_id, kwargs...)
    grouped = exons_by(df; group_col=transcript_col, kwargs...)
    Dict(transcript => sum(width(collection)) for (transcript, collection) in grouped)
end

function _feature_kind_subset(df::DataFrames.AbstractDataFrame, kinds; type_col=:type, kwargs...)
    column = _feature_type_column(df, type_col)
    wanted = Set(lowercase.(String[string(kind) for kind in kinds]))
    indices = findall(value -> !ismissing(value) && lowercase(String(value)) in wanted, df[!, column])
    make_granges_from_dataframe(df[indices, :]; kwargs...)
end

five_utrs_by_transcript(df::DataFrames.AbstractDataFrame; transcript_col=:transcript_id, type_col=:type, kwargs...) =
    make_granges_list_from_dataframe(_feature_kind_subset(df, ("five_prime_utr", "5utr", "5'utr"); type_col=type_col, kwargs...); group_col=transcript_col)
three_utrs_by_transcript(df::DataFrames.AbstractDataFrame; transcript_col=:transcript_id, type_col=:type, kwargs...) =
    make_granges_list_from_dataframe(_feature_kind_subset(df, ("three_prime_utr", "3utr", "3'utr"); type_col=type_col, kwargs...); group_col=transcript_col)

transcript_lengths(df::DataFrames.AbstractDataFrame; transcript_col=:transcript_id, kwargs...) =
    transcript_widths(df; transcript_col=transcript_col, kwargs...)

function _tidy_feature_table(df::DataFrames.AbstractDataFrame; kwargs...)
    values = make_granges_from_dataframe(df; kwargs...)
    order = sortperm(eachindex(values), by=index -> (values[index].chrom, values[index].left, values[index].right, values[index].strand))
    df[order, :]
end

tidy_transcripts(df::DataFrames.AbstractDataFrame; type_col=:type, kwargs...) = _tidy_feature_table(_feature_dataframe(df, "transcript"; type_col=type_col); kwargs...)
tidy_exons(df::DataFrames.AbstractDataFrame; type_col=:type, kwargs...) = _tidy_feature_table(_feature_dataframe(df, "exon"; type_col=type_col); kwargs...)
tidy_introns(df::DataFrames.AbstractDataFrame; type_col=:type, kwargs...) = _tidy_feature_table(_feature_dataframe(df, "intron"; type_col=type_col); kwargs...)

function exonic_parts(df::DataFrames.AbstractDataFrame; kwargs...)
    exons_collection = build_collection(exons(df; kwargs...))
    disjoin(exons_collection)
end

function intronic_parts(df::DataFrames.AbstractDataFrame; transcript_col=:transcript_id, kwargs...)
    introns = reduce(vcat, values(introns_by_transcript(df; transcript_col=transcript_col, kwargs...)); init=GenomicInterval[])
    disjoin(build_collection(introns))
end

function id2name(df::DataFrames.AbstractDataFrame; id_col=:id, name_col=:name)
    id_col = Symbol(id_col)
    name_col = Symbol(name_col)
    Dict(df[!, id_col][index] => df[!, name_col][index] for index in eachindex(df[!, id_col])
         if !ismissing(df[!, id_col][index]) && !ismissing(df[!, name_col][index]))
end

const fiveUTRsByTranscript = five_utrs_by_transcript
const threeUTRsByTranscript = three_utrs_by_transcript
const transcriptLengths = transcript_lengths
const tidyTranscripts = tidy_transcripts
const tidyExons = tidy_exons
const tidyIntrons = tidy_introns
const exonicParts = exonic_parts
const intronicParts = intronic_parts

function make_iranges_from_dataframe(df::DataFrames.AbstractDataFrame; start_col=:start, end_col=:stop, zero_based::Bool=false)
    start_col = Symbol(start_col)
    end_col = Symbol(end_col)
    starts = zero_based ? [ismissing(value) ? missing : Int(value) + 1 for value in df[!, start_col]] : df[!, start_col]
    ends = df[!, end_col]
    valid = [index for index in eachindex(starts) if !ismissing(starts[index]) && !ismissing(ends[index])]
    IRanges(start=starts[valid], end_position=ends[valid])
end

granges(interval::GenomicInterval) = [interval]
granges(collection::IntervalCollection) = collection.intervals
granges(intervals::AbstractVector{<:GenomicInterval}) = collect(intervals)

"""Group intervals by a metadata field, returning a dictionary of collections."""
function grglist(intervals::AbstractVector{<:GenomicInterval}, field::AbstractString)
    groups = Dict{Any,GenomicInterval[] }()
    for interval in intervals
        value = get(_comparable_interval_metadata(interval.metadata), field, missing)
        ismissing(value) && continue
        push!(get!(groups, value, GenomicInterval[]), interval)
    end
    Dict(value => build_collection(group) for (value, group) in groups)
end
grglist(collection::IntervalCollection, field::AbstractString) = grglist(collection.intervals, field)

subtract(x, y) = Base.setdiff(x, y)

is_small_genome(seqlengths::AbstractDict; threshold::Integer=100_000_000) =
    sum(Int(value) for value in values(seqlengths)) <= threshold

function _ordered_seqlengths(seqlengths::AbstractDict, chrom_order)
    chromosomes = chrom_order === nothing ? String[string(chrom) for chrom in keys(seqlengths)] : String[string(chrom) for chrom in chrom_order]
    all(haskey(seqlengths, chrom) for chrom in chromosomes) || throw(ArgumentError("chrom_order contains an unknown chromosome"))
    [(chrom, Int(seqlengths[chrom])) for chrom in chromosomes]
end

function absolute_ranges(intervals::AbstractVector{<:GenomicInterval}, seqlengths::AbstractDict; chrom_order=nothing)
    ordered = _ordered_seqlengths(seqlengths, chrom_order)
    offsets = Dict{String,Int}()
    cursor = 0
    for (chrom, length_chrom) in ordered
        offsets[chrom] = cursor
        cursor += length_chrom
    end
    [IRange(offsets[interval.chrom] + interval.left, width(interval)) for interval in intervals if haskey(offsets, interval.chrom)]
end

function relative_ranges(ranges::AbstractVector{<:IRange}, seqlengths::AbstractDict; chrom_order=nothing)
    ordered = _ordered_seqlengths(seqlengths, chrom_order)
    result = GenomicInterval[]
    offsets = 0
    for (chrom, length_chrom) in ordered
        chromosome_end = offsets + length_chrom
        for interval in ranges
            left = max(start(interval), offsets + 1)
            right = min(stop(interval), chromosome_end)
            left <= right || continue
            push!(result, GenomicInterval(chrom, left - offsets, right - offsets))
        end
        offsets = chromosome_end
    end
    result
end

function which_as_iranges(mask::AbstractVector{Bool})
    result = IRange[]
    index = 1
    while index <= length(mask)
        mask[index] || (index += 1; continue)
        first_index = index
        while index < length(mask) && mask[index + 1]
            index += 1
        end
        push!(result, IRange(first_index, index - first_index + 1))
        index += 1
    end
    result
end

function as_normal_iranges(ranges::AbstractVector{<:IRange})
    isempty(ranges) && return IRange[]
    ordered = sort(collect(ranges), by=x -> (start(x), stop(x)))
    result = IRange[ordered[1]]
    for current in Iterators.drop(ordered, 1)
        previous = result[end]
        if start(current) <= stop(previous) + 1
            result[end] = IRange(start(previous), max(stop(previous), stop(current)) - start(previous) + 1)
        else
            push!(result, current)
        end
    end
    result
end

function break_in_chunks(total_size::Integer, n_chunks::Integer, chunk_size::Integer=0)
    total_size >= 0 || throw(ArgumentError("total_size must be nonnegative"))
    n_chunks > 0 || throw(ArgumentError("n_chunks must be positive"))
    width = chunk_size > 0 ? Int(chunk_size) : max(1, cld(Int(total_size), Int(n_chunks)))
    result = UnitRange{Int}[]
    first_index = 1
    while first_index <= total_size
        last_index = min(Int(total_size), first_index + width - 1)
        push!(result, first_index:last_index)
        first_index = last_index + 1
    end
    result
end

heads(intervals::AbstractVector{<:GenomicInterval}, n::Integer=6) = intervals[1:min(Int(n), length(intervals))]
tails(intervals::AbstractVector{<:GenomicInterval}, n::Integer=6) = intervals[max(1, length(intervals) - Int(n) + 1):end]
heads(ranges::AbstractVector{<:IRange}, n::Integer=6) = ranges[1:min(Int(n), length(ranges))]
tails(ranges::AbstractVector{<:IRange}, n::Integer=6) = ranges[max(1, length(ranges) - Int(n) + 1):end]

shift(ranges::AbstractVector{<:IRange}, delta::Integer) = [IRange(start(interval) + Int(delta), width(interval)) for interval in ranges]

function resize(ranges::AbstractVector{<:IRange}, new_width::Integer; fix::Symbol=:start)
    new_width >= 0 || throw(ArgumentError("new_width must be nonnegative"))
    fix in (:start, :center, :end) || throw(ArgumentError("fix must be :start, :center, or :end"))
    [fix == :start ? IRange(start(interval), new_width) :
     fix == :end ? IRange(stop(interval) - Int(new_width) + 1, new_width) :
     IRange(mid(interval) - (Int(new_width) - 1) ÷ 2, new_width) for interval in ranges]
end

narrow(ranges::AbstractVector{<:IRange}, left::Integer=0, right::Integer=0) = [
    begin
        new_start = start(interval) + Int(left)
        new_width = width(interval) - Int(left) - Int(right)
        new_width >= 0 || throw(ArgumentError("narrowing removes the entire range"))
        IRange(new_start, new_width)
    end for interval in ranges]

function restrict(ranges::AbstractVector{<:IRange}, lower::Integer, upper::Integer)
    lower <= upper || throw(ArgumentError("lower must not exceed upper"))
    result = IRange[]
    for interval in ranges
        left = max(start(interval), Int(lower))
        right = min(stop(interval), Int(upper))
        left <= right && push!(result, IRange(left, right - left + 1))
    end
    result
end

function reduce_ranges(ranges::AbstractVector{<:IRange}; minoverlap::Integer=1)
    minoverlap >= 1 || throw(ArgumentError("minoverlap must be positive"))
    isempty(ranges) && return IRange[]
    ordered = sort(collect(ranges), by=x -> (start(x), stop(x)))
    result = IRange[ordered[1]]
    for current in Iterators.drop(ordered, 1)
        previous = result[end]
        if start(current) <= stop(previous) - Int(minoverlap) + 1
            result[end] = IRange(start(previous), max(stop(previous), stop(current)) - start(previous) + 1)
        else
            push!(result, current)
        end
    end
    result
end

function disjoin(ranges::AbstractVector{<:IRange})
    points = sort!(unique!(vcat([start(interval) for interval in ranges], [stop(interval) + 1 for interval in ranges])))
    result = IRange[]
    for index in 1:max(0, length(points) - 1)
        left, right = points[index], points[index + 1] - 1
        any(interval -> start(interval) <= left && stop(interval) >= right, ranges) && push!(result, IRange(left, right - left + 1))
    end
    result
end

function gaps(ranges::AbstractVector{<:IRange}, lower::Integer, upper::Integer)
    occupied = reduce_ranges(restrict(ranges, lower, upper))
    result = IRange[]
    cursor = Int(lower)
    for interval in occupied
        cursor < start(interval) && push!(result, IRange(cursor, start(interval) - cursor))
        cursor = max(cursor, stop(interval) + 1)
    end
    result
end

function pintersect(x::AbstractVector{<:IRange}, y::AbstractVector{<:IRange})
    length(x) == length(y) || throw(ArgumentError("range vectors must have equal length"))
    result = IRange[]
    for (a, b) in zip(x, y)
        if max(start(a), start(b)) <= min(stop(a), stop(b))
            push!(result, IRange(max(start(a), start(b)), min(stop(a), stop(b)) - max(start(a), start(b)) + 1))
        end
    end
    return result
end

function punion(x::AbstractVector{<:IRange}, y::AbstractVector{<:IRange})
    length(x) == length(y) || throw(ArgumentError("range vectors must have equal length"))
    [IRange(min(start(a), start(b)), max(stop(a), stop(b)) - min(start(a), start(b)) + 1) for (a, b) in zip(x, y)]
end

function psetdiff(x::AbstractVector{<:IRange}, y::AbstractVector{<:IRange})
    length(x) == length(y) || throw(ArgumentError("range vectors must have equal length"))
    [begin
        pieces = IRange[a]
        overlap = pintersect([a], [b])[1]
        overlap === nothing ? pieces : begin
            pieces = IRange[]
            start(a) < start(overlap) && push!(pieces, IRange(start(a), start(overlap) - start(a)))
            stop(overlap) < stop(a) && push!(pieces, IRange(stop(overlap) + 1, stop(a) - stop(overlap)))
            pieces
        end
    end for (a, b) in zip(x, y)]
end

distance(x::IRange, y::IRange) = max(max(start(x), start(y)) - min(stop(x), stop(y)) - 1, 0)

is_disjoint(ranges::AbstractVector{<:IRange}) = begin
    ordered = sort(collect(ranges), by=x -> (start(x), stop(x)))
    all(start(ordered[index]) > stop(ordered[index - 1]) for index in 2:length(ordered))
end

function disjoint_bins(ranges::AbstractVector{<:IRange})
    result = zeros(Int, length(ranges))
    bin_ends = Int[]
    for index in sortperm(ranges, by=x -> (start(x), stop(x)))
        interval = ranges[index]
        bin = findfirst(endpoint -> endpoint < start(interval), bin_ends)
        if bin === nothing
            push!(bin_ends, stop(interval))
            result[index] = length(bin_ends)
        else
            bin_ends[bin] = stop(interval)
            result[index] = bin
        end
    end
    result
end

function tile(ranges::AbstractVector{<:IRange}, n::Integer)
    n > 0 || throw(ArgumentError("n must be positive"))
    [begin
        q, remainder = divrem(width(interval), Int(n))
        cursor = start(interval)
        [begin
            tile_width = q + (tile_index <= remainder ? 1 : 0)
            current = IRange(cursor, tile_width)
            cursor += tile_width
            current
        end for tile_index in 1:Int(n) if q + (tile_index <= remainder ? 1 : 0) > 0]
    end for interval in ranges]
end

function sliding_windows(ranges::AbstractVector{<:IRange}, window_width::Integer; shift::Integer=1)
    window_width > 0 || throw(ArgumentError("window_width must be positive"))
    shift > 0 || throw(ArgumentError("shift must be positive"))
    [[IRange(position, window_width) for position in start(interval):Int(shift):(stop(interval) - Int(window_width) + 1)] for interval in ranges if width(interval) >= window_width]
end

function distance_to_nearest(queries::AbstractVector{<:GenomicInterval}, subject::IntervalCollection)
    result = NamedTuple{(:query_index, :subject, :distance),Tuple{Int,Union{Nothing,GenomicInterval},Int}}[]
    for (index, query) in Base.pairs(queries)
        hit = nearest(query, subject)
        push!(result, (query_index=index, subject=hit, distance=hit === nothing ? -1 : _distance(query, hit)))
    end
    result
end

function _copy_metadata(interval::GenomicInterval)
    return Dict{String,Any}(String(key) => value for (key, value) in pairs(interval.metadata))
end

function _new_interval(interval::GenomicInterval, left::Integer, right::Integer; strand::Char=interval.strand, metadata=interval.metadata)
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

"""Return the query elements having at least one indexed subject overlap."""
function subset_by_overlaps(
    queries::AbstractVector{<:GenomicInterval},
    subject::IntervalCollection;
    ignore_strand::Bool=true)
    return [query for query in queries if hasintersection(query, subject; ignore_strand=ignore_strand)]
end

subset_by_overlaps(query::IntervalCollection, subject::IntervalCollection; ignore_strand::Bool=true) =
    build_collection(subset_by_overlaps(query.intervals, subject; ignore_strand=ignore_strand))

overlaps_any(queries::AbstractVector{<:GenomicInterval}, subject::IntervalCollection; ignore_strand::Bool=true) =
    [hasintersection(query, subject; ignore_strand=ignore_strand) for query in queries]
overlaps_any(query::GenomicInterval, subject::IntervalCollection; ignore_strand::Bool=true) =
    hasintersection(query, subject; ignore_strand=ignore_strand)

"""Clip intervals to an inclusive coordinate range, dropping empty results."""
function restrict(
    intervals::AbstractVector{<:GenomicInterval},
    lower::Integer,
    upper::Integer)
    lower <= upper || throw(ArgumentError("lower must not exceed upper"))
    result = GenomicInterval[]
    for interval in intervals
        left = max(interval.left, Int(lower))
        right = min(interval.right, Int(upper))
        left <= right || continue
        push!(result, _new_interval(interval, left, right))
    end
    return result
end

restrict(collection::IntervalCollection, lower::Integer, upper::Integer) =
    build_collection(restrict(collection.intervals, lower, upper))

"""Return the terminal promoter-like window around each interval's end."""
function terminators(intervals::AbstractVector{<:GenomicInterval}, upstream::Integer, downstream::Integer)
    upstream >= 0 || throw(ArgumentError("upstream must be nonnegative"))
    downstream >= 0 || throw(ArgumentError("downstream must be nonnegative"))
    result = GenomicInterval[]
    for interval in intervals
        left, right = if interval.strand == '-'
            (interval.left - Int(downstream), interval.left + Int(upstream) - 1)
        else
            (interval.right - Int(upstream) + 1, interval.right + Int(downstream))
        end
        push!(result, _new_interval(interval, left, right))
    end
    return result
end

terminators(collection::IntervalCollection, upstream::Integer, downstream::Integer) =
    build_collection(terminators(collection.intervals, upstream, downstream))

"""Partition each interval into at most `n` nearly equal contiguous tiles."""
function tile(intervals::AbstractVector{<:GenomicInterval}, n::Integer)
    n > 0 || throw(ArgumentError("n must be positive"))
    result = Vector{Vector{GenomicInterval}}(undef, length(intervals))
    for (index, interval) in pairs(intervals)
        q, r = divrem(width(interval), Int(n))
        tiles = GenomicInterval[]
        cursor = interval.left
        for tile_index in 1:Int(n)
            tile_width = q + (tile_index <= r ? 1 : 0)
            tile_width == 0 && continue
            push!(tiles, _new_interval(interval, cursor, cursor + tile_width - 1))
            cursor += tile_width
        end
        result[index] = tiles
    end
    return result
end

tile(collection::IntervalCollection, n::Integer) = tile(collection.intervals, n)

"""Tile a genome described by chromosome lengths into contiguous bins."""
function tile_genome(seqlengths::AbstractDict, n::Integer)
    n > 0 || throw(ArgumentError("n must be positive"))
    chromosomes = sort!(String[string(chrom) for chrom in keys(seqlengths)])
    lengths = Dict(chrom => Int(seqlengths[chrom]) for chrom in chromosomes)
    total = sum(values(lengths))
    total > 0 || return GenomicInterval[]
    target = max(1, cld(total, Int(n)))
    result = GenomicInterval[]
    for chrom in chromosomes
        position = 1
        while position <= lengths[chrom]
            right = min(lengths[chrom], position + target - 1)
            push!(result, GenomicInterval(chrom, position, right))
            position = right + 1
        end
    end
    return result
end

const isSmallGenome = is_small_genome
const absoluteRanges = absolute_ranges
const relativeRanges = relative_ranges
const tileGenome = tile_genome

"""Generate fixed-width sliding windows over each interval."""
function sliding_windows(intervals::AbstractVector{<:GenomicInterval}, window_width::Integer; shift::Integer=1)
    window_width > 0 || throw(ArgumentError("window_width must be positive"))
    shift > 0 || throw(ArgumentError("shift must be positive"))
    result = Vector{Vector{GenomicInterval}}(undef, length(intervals))
    for (index, interval) in pairs(intervals)
        windows = GenomicInterval[]
        position = interval.left
        while position + Int(window_width) - 1 <= interval.right
            push!(windows, _new_interval(interval, position, position + Int(window_width) - 1))
            position += Int(shift)
        end
        result[index] = windows
    end
    return result
end

sliding_windows(collection::IntervalCollection, window_width::Integer; shift::Integer=1) =
    sliding_windows(collection.intervals, window_width; shift=shift)

"""Whether intervals are pairwise disjoint on each chromosome."""
function is_disjoint(intervals::AbstractVector{<:GenomicInterval})
    isempty(intervals) && return true
    ordered = sort(intervals, by=_sort_key)
    previous = ordered[1]
    for current in Iterators.drop(ordered, 1)
        current.chrom == previous.chrom && current.left <= previous.right && return false
        current.chrom == previous.chrom && (previous = current)
        current.chrom != previous.chrom && (previous = current)
    end
    return true
end

is_disjoint(collection::IntervalCollection) = is_disjoint(collection.intervals)

"""Assign the minimum greedy set of non-overlapping bins to intervals."""
function disjoint_bins(intervals::AbstractVector{<:GenomicInterval})
    result = zeros(Int, length(intervals))
    last_end = Dict{String,Vector{Int}}()
    for index in sortperm(intervals, by=_sort_key)
        interval = intervals[index]
        ends = get!(last_end, interval.chrom, Int[])
        chosen = findfirst(endpoint -> endpoint < interval.left, ends)
        if chosen === nothing
            push!(ends, interval.right)
            result[index] = length(ends)
        else
            ends[chosen] = interval.right
            result[index] = chosen
        end
    end
    return result
end

disjoint_bins(collection::IntervalCollection) = disjoint_bins(collection.intervals)

order_intervals(intervals::AbstractVector{<:GenomicInterval}; decreasing::Bool=false) =
    sortperm(intervals, by=_sort_key, rev=decreasing)
order_intervals(collection::IntervalCollection; decreasing::Bool=false) =
    order_intervals(collection.intervals; decreasing=decreasing)

duplicated_intervals(intervals::AbstractVector{<:GenomicInterval}) = begin
    seen = Set{GenomicInterval}()
    result = falses(length(intervals))
    for (index, interval) in pairs(intervals)
        result[index] = interval in seen
        push!(seen, interval)
    end
    result
end

duplicated_intervals(collection::IntervalCollection) = duplicated_intervals(collection.intervals)

"""Return one envelope interval per chromosome."""
function range_intervals(intervals::AbstractVector{<:GenomicInterval})
    isempty(intervals) && return GenomicInterval[]
    grouped = Dict{String,Vector{GenomicInterval}}()
    for interval in intervals
        push!(get!(grouped, interval.chrom, GenomicInterval[]), interval)
    end
    return [begin
        left = minimum(interval.left for interval in group)
        right = maximum(interval.right for interval in group)
        GenomicInterval(chrom, left, right, '.', Dict{String,Any}())
    end for (chrom, group) in sort!(collect(grouped); by=first)]
end

range_intervals(collection::IntervalCollection) = range_intervals(collection.intervals)

"""Pairwise overlap predicate, analogous to Bioconductor's `poverlaps`."""
function poverlaps(x::AbstractVector{<:GenomicInterval}, y::AbstractVector{<:GenomicInterval}; ignore_strand::Bool=true)
    length(x) == length(y) || throw(ArgumentError("interval vectors must have the same length"))
    return [a.chrom == b.chrom && a.left <= b.right && a.right >= b.left &&
            (ignore_strand || a.strand == '.' || b.strand == '.' || a.strand == b.strand)
            for (a, b) in zip(x, y)]
end

"""Return clipped overlap intervals for each query."""
function overlaps_ranges(query::GenomicInterval, subject::IntervalCollection; ignore_strand::Bool=true)
    return [GenomicInterval(query.chrom, max(query.left, hit.left), min(query.right, hit.right), query.strand, _copy_metadata(query))
            for hit in find_overlaps(query, subject; ignore_strand=ignore_strand)]
end

overlaps_ranges(queries::AbstractVector{<:GenomicInterval}, subject::IntervalCollection; ignore_strand::Bool=true) =
    [overlaps_ranges(query, subject; ignore_strand=ignore_strand) for query in queries]

"""Return zero-based-free `(query_index, subject_index)` overlap pairs."""
function find_overlap_pairs(queries::AbstractVector{<:GenomicInterval}, subject::IntervalCollection; ignore_strand::Bool=true)
    hits = Tuple{Int,Int}[]
    for (query_index, query) in Base.pairs(queries)
        subject_indices = _collection_overlap_indices(subject, query.chrom, query.left, query.right)
        if !ignore_strand && query.strand != '.'
            filter!(index -> begin
                hit_strand = subject.intervals[index].strand
                hit_strand == query.strand || hit_strand == '.'
            end, subject_indices)
        end
        for subject_index in subject_indices
            push!(hits, (query_index, subject_index))
        end
    end
    return hits
end

"""Merge overlapping query intervals and attach their source count."""
merge_by_overlaps(intervals::AbstractVector{<:GenomicInterval}) = reduce(build_collection(intervals))
merge_by_overlaps(collection::IntervalCollection) = reduce(collection)

is_normal(intervals::AbstractVector{<:GenomicInterval}) = is_disjoint(intervals) &&
    all(intervals[index].chrom != intervals[index - 1].chrom || intervals[index].left > intervals[index - 1].right + 1
        for index in 2:length(intervals))
is_normal(collection::IntervalCollection) = is_normal(collection.intervals)

function _sort_key(interval::GenomicInterval)
    return (interval.chrom, interval.left, interval.right, interval.strand)
end

function _build_chrom_indices(intervals::AbstractVector{<:GenomicInterval})
    chrom_indices = Dict{String,UnitRange{Int}}()
    chrom_end_indices = Dict{String,Vector{Int}}()
    prefix_max_ends = Dict{String,Vector{Int}}()
    isempty(intervals) && return chrom_indices, chrom_end_indices, prefix_max_ends

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

    for (chrom, range) in chrom_indices
        maxima = Vector{Int}(undef, length(range))
        running_max = typemin(Int)
        for (local_index, index) in enumerate(range)
            running_max = max(running_max, intervals[index].right)
            maxima[local_index] = running_max
        end
        prefix_max_ends[chrom] = maxima
    end

    return chrom_indices, chrom_end_indices, prefix_max_ends
end

function build_collection(intervals::AbstractVector{<:GenomicInterval})
    # Normalize annotation storage once while retaining a compact concrete
    # empty-metadata representation for large unannotated interval sets.
    annotated = any(!isempty(_comparable_interval_metadata(interval.metadata)) for interval in intervals)
    sorted = if annotated
        GenomicIntervalDict[
            GenomicInterval(String(interval.chrom), Int(interval.left), Int(interval.right), Char(interval.strand), _copy_metadata(interval))
            for interval in intervals]
    else
        GenomicIntervalEmpty[
            GenomicInterval(String(interval.chrom), Int(interval.left), Int(interval.right), Char(interval.strand), NamedTuple())
            for interval in intervals]
    end
    sort!(sorted, by=_sort_key)
    chrom_indices, chrom_end_indices, prefix_max_ends = _build_chrom_indices(sorted)
    chroms = PooledArray([interval.chrom for interval in sorted])
    starts = Int[interval.left for interval in sorted]
    ends = Int[interval.right for interval in sorted]
    strands = Char[interval.strand for interval in sorted]
    result = IntervalCollection(chroms, starts, ends, strands, sorted, chrom_indices, chrom_end_indices, prefix_max_ends)
    _ctx = active_provenance_context()

    return _register_genomicranges_result!(_ctx, result, "build_collection"; parameters=(n_intervals=length(sorted),))
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

"""Compact overlap index query over sorted starts and per-chromosome prefix maxima."""
function _collection_overlap_indices(subject::IntervalCollection, chrom::String, left::Int, right::Int)
    range = get(subject.chrom_indices, chrom, nothing)
    range === nothing && return Int[]
    last_index = _find_last_start_leq(subject.starts, range, right)
    last_index === nothing && return Int[]
    first_index = first(range)
    prefix = subject.prefix_max_ends[chrom]
    result = Int[]
    index = last_index
    while index >= first_index
        subject.ends[index] >= left && push!(result, index)
        index == first_index && break
        local_index = index - first_index + 1
        prefix[local_index - 1] < left && break
        index -= 1
    end
    sort!(result)
    return result
end

function _collection_overlap_count(subject::IntervalCollection, chrom::String, left::Int, right::Int)
    range = get(subject.chrom_indices, chrom, nothing)
    range === nothing && return 0
    last_index = _find_last_start_leq(subject.starts, range, right)
    last_index === nothing && return 0
    first_index = first(range)
    prefix = subject.prefix_max_ends[chrom]
    count = 0
    index = last_index
    while index >= first_index
        subject.ends[index] >= left && (count += 1)
        index == first_index && break
        local_index = index - first_index + 1
        prefix[local_index - 1] < left && break
        index -= 1
    end
    return count
end

function find_overlaps(query::GenomicInterval, subject::IntervalCollection; ignore_strand::Bool=true)
    indices = _collection_overlap_indices(subject, query.chrom, query.left, query.right)
    if !ignore_strand && query.strand != '.'
        filter!(i -> subject.intervals[i].strand == query.strand || subject.intervals[i].strand == '.', indices)
    end

    return with_provenance(subject.intervals[indices], "GenomicInterval", "GenomicRanges/find_overlaps"; notes=["overlap query against interval collection"], parameters=(chrom=query.chrom, left=query.left, right=query.right, hit_count=length(indices)))
end

"""
    count_overlaps(query, subject; ignore_strand=true)

Return the number of overlapping intervals without allocating the full result.
Equivalent to Bioconductor's `countOverlaps`.
"""
function count_overlaps(query::GenomicInterval, subject::IntervalCollection; ignore_strand::Bool=true)
    if !ignore_strand && query.strand != '.'
        indices = _collection_overlap_indices(subject, query.chrom, query.left, query.right)
        return count(i -> subject.intervals[i].strand == query.strand || subject.intervals[i].strand == '.', indices)
    end
    return _collection_overlap_count(subject, query.chrom, query.left, query.right)
end

"""Return all subject intervals overlapping every query interval."""
function find_overlaps(
    queries::AbstractVector{<:GenomicInterval},
    subject::IntervalCollection;
    ignore_strand::Bool=true)
    return [find_overlaps(query, subject; ignore_strand=ignore_strand) for query in queries]
end

function find_overlaps(query::IRange, subject::AbstractVector{<:IRange})
    [index for (index, interval) in Base.pairs(subject) if start(interval) <= stop(query) && stop(interval) >= start(query)]
end
find_overlaps(queries::AbstractVector{<:IRange}, subject::AbstractVector{<:IRange}) =
    [find_overlaps(query, subject) for query in queries]
count_overlaps(query::IRange, subject::AbstractVector{<:IRange}) = length(find_overlaps(query, subject))
count_overlaps(queries::AbstractVector{<:IRange}, subject::AbstractVector{<:IRange}) = [count_overlaps(query, subject) for query in queries]
overlaps_any(query::IRange, subject::AbstractVector{<:IRange}) = count_overlaps(query, subject) > 0
overlaps_any(queries::AbstractVector{<:IRange}, subject::AbstractVector{<:IRange}) = [overlaps_any(query, subject) for query in queries]
subset_by_overlaps(queries::AbstractVector{<:IRange}, subject::AbstractVector{<:IRange}) = [query for query in queries if overlaps_any(query, subject)]

find_overlaps(query::IntervalCollection, subject::IntervalCollection; ignore_strand::Bool=true) =
    find_overlaps(query.intervals, subject; ignore_strand=ignore_strand)
find_overlaps(query::GenomicInterval, subject::AbstractVector{<:GenomicInterval}; ignore_strand::Bool=true) =
    find_overlaps(query, build_collection(subject); ignore_strand=ignore_strand)

"""Count overlaps for a vector or collection without materializing hit intervals."""
function count_overlaps(
    queries::AbstractVector{<:GenomicInterval},
    subject::IntervalCollection;
    ignore_strand::Bool=true)
    return [count_overlaps(query, subject; ignore_strand=ignore_strand) for query in queries]
end

count_overlaps(query::IntervalCollection, subject::IntervalCollection; ignore_strand::Bool=true) =
    count_overlaps(query.intervals, subject; ignore_strand=ignore_strand)
count_overlaps(query::GenomicInterval, subject::AbstractVector{<:GenomicInterval}; ignore_strand::Bool=true) =
    count_overlaps(query, build_collection(subject); ignore_strand=ignore_strand)

"""Lazy overlap iterator for a single query, yielding subject intervals."""
function eachoverlap(query::GenomicInterval, subject::IntervalCollection; ignore_strand::Bool=true)
    return (hit for hit in find_overlaps(query, subject; ignore_strand=ignore_strand))
end

function eachoverlap(
    queries::AbstractVector{<:GenomicInterval},
    subject::IntervalCollection;
    ignore_strand::Bool=true)
    return ((query, hit) for query in queries for hit in eachoverlap(query, subject; ignore_strand=ignore_strand))
end

eachoverlap(query::IntervalCollection, subject::IntervalCollection; ignore_strand::Bool=true) =
    eachoverlap(query.intervals, subject; ignore_strand=ignore_strand)

"""Test for at least one overlap without allocating hit intervals."""
hasintersection(query::GenomicInterval, subject::IntervalCollection; ignore_strand::Bool=true) =
    count_overlaps(query, subject; ignore_strand=ignore_strand) > 0

@inline function _naive_overlap_count(query::GenomicInterval, subject::IntervalCollection)
    hits = 0
    for interval in subject.intervals
        interval.chrom == query.chrom || continue
        interval.left <= query.right && interval.right >= query.left && (hits += 1)
    end
    return hits
end

function _profile_overlap_kernel(kernel::Function, repetitions::Int)
    # Run once to compile specialized methods outside the measured samples.
    kernel()
    times = Vector{Float64}(undef, repetitions)
    allocations = Vector{Int}(undef, repetitions)
    hit_counts = Vector{Int}(undef, repetitions)
    for i in 1:repetitions
        measurement = @timed kernel()
        times[i] = measurement.time
        allocations[i] = measurement.bytes
        hit_counts[i] = measurement.value
    end
    return median(times), round(Int, median(allocations)), hit_counts[end]
end

"""
    profile_interval_overlaps(queries, subject; repetitions=5)

Measure the public indexed overlap path against an exact direct scan. Both
paths must produce the same total hit count; a mismatch raises an error rather
than reporting a misleading performance comparison. Timings exclude JIT
compilation and report median wall time and allocation volume across repeats.
"""
function profile_interval_overlaps(
    queries,
    subject::IntervalCollection;
    repetitions::Integer=5,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    repetitions > 0 || throw(ArgumentError("repetitions must be positive"))
    query_vector = GenomicInterval[query for query in queries]

    indexed_kernel = () -> sum(length(find_overlaps(query, subject)) for query in query_vector)
    naive_kernel = () -> sum(_naive_overlap_count(query, subject) for query in query_vector)
    indexed_seconds, indexed_bytes, indexed_hits = _profile_overlap_kernel(indexed_kernel, Int(repetitions))
    naive_seconds, naive_bytes, naive_hits = _profile_overlap_kernel(naive_kernel, Int(repetitions))
    indexed_hits == naive_hits || throw(ErrorException("indexed and naive overlap hit counts disagree: $(indexed_hits) != $(naive_hits)"))

    result = OverlapProfileResult(
        length(query_vector),
        length(subject),
        Int(repetitions),
        indexed_seconds,
        naive_seconds,
        indexed_bytes,
        naive_bytes,
        indexed_hits,
        naive_hits)
    return _register_genomicranges_result!(
        _ctx,
        result,
        "profile_interval_overlaps";
        parents=provenance_parent_ids(query_vector, subject),
        parameters=(
            query_count=result.query_count,
            subject_count=result.subject_count,
            repetitions=result.repetitions,
            indexed_hit_count=result.indexed_hit_count,
            naive_hit_count=result.naive_hit_count))
end

overlap(query::GenomicInterval, subject::IntervalCollection) = find_overlaps(query, subject)

function _distance(query::GenomicInterval, interval::GenomicInterval)
    interval.chrom == query.chrom || return -1
    interval.right < query.left && return query.left - interval.right - 1
    interval.left > query.right && return interval.left - query.right - 1
    return 0
end

"""
    distance(query::GenomicInterval, subject::IntervalCollection)

Calculate the distance to the nearest interval in the subject collection.
Returns the minimum distance as an integer, or -1 if no intervals exist on the same chromosome.
"""
function distance(query::GenomicInterval, subject::IntervalCollection)
    nearest_int = nearest(query, subject)
    nearest_int === nothing && return -1

    return _distance(query, nearest_int)
end

"""Find the nearest interval; `select=:all` retains every minimum-distance tie."""
function nearest(query::GenomicInterval, subject::IntervalCollection; select::Symbol=:first)
    select in (:first, :all) || throw(ArgumentError("select must be :first or :all"))
    range = get(subject.chrom_indices, query.chrom, nothing)
    range === nothing && return nothing

    # Fast path: check for overlapping intervals via interval tree (distance=0)
    containing = _collection_overlap_indices(subject, query.chrom, query.left, query.right)
    if !isempty(containing)
        selected = select == :all ? subject.intervals[containing] : subject.intervals[containing[1]]
        return with_provenance(selected, "GenomicInterval", "GenomicRanges/nearest"; parameters=(chrom=query.chrom, left=query.left, right=query.right, select=select))
    end

    best_dist = typemax(Int)
    candidates = Int[]

    # Every nearest right-side tie has the same smallest start coordinate.
    right_idx = searchsortedfirst(subject.starts, query.right + 1, first(range), last(range), Base.Order.Forward)
    if right_idx <= last(range)
        d = subject.starts[right_idx] - query.right - 1
        best_dist = d
        for index in right_idx:last(range)
            subject.starts[index] == subject.starts[right_idx] || break
            push!(candidates, index)
        end
    end

    # Every nearest left-side tie has the same largest end coordinate.
    end_indices = get(subject.chrom_end_indices, query.chrom, nothing)
    if end_indices !== nothing
        left_idx = _find_last_end_lt(subject.ends, end_indices, query.left)
        if left_idx !== nothing
            d = query.left - subject.ends[left_idx] - 1
            if d < best_dist
                best_dist = d
                empty!(candidates)
            end
            if d == best_dist
                end_position = findfirst(==(left_idx), end_indices)
                for position in end_position:-1:1
                    index = end_indices[position]
                    subject.ends[index] == subject.ends[left_idx] || break
                    push!(candidates, index)
                end
            end
        end
    end

    isempty(candidates) && return nothing
    sort!(unique!(candidates))
    selected = select == :all ? subject.intervals[candidates] : subject.intervals[candidates[1]]
    return with_provenance(selected, "GenomicInterval", "GenomicRanges/nearest"; parameters=(chrom=query.chrom, left=query.left, right=query.right, select=select))
end

find_nearest(query::GenomicInterval, subject::IntervalCollection; select::Symbol=:first) = nearest(query, subject; select=select)
nearest_all(query::GenomicInterval, subject::IntervalCollection) = begin
    result = nearest(query, subject; select=:all)
    result === nothing ? GenomicIntervalDict[] : GenomicIntervalDict[
        GenomicInterval(String(interval.chrom), Int(interval.left), Int(interval.right), Char(interval.strand), _copy_metadata(interval))
        for interval in result]
end
nearest_first(query::GenomicInterval, subject::IntervalCollection) = nearest(query, subject; select=:first)

"""
    precede(query::GenomicInterval, subject::IntervalCollection)

In Bioconductor, `precede(x, y)` is the index of the element in `y` that is preceded
by `x`. This is strand-aware:
- For `+` strand (or `.`): `y` is to the right of `x`.
- For `-` strand: `y` is to the left of `x`.
"""
function precede(query::GenomicInterval, subject::IntervalCollection)
    range = get(subject.chrom_indices, query.chrom, nothing)
    range === nothing && return nothing

    # Strand-aware: for '-' strand, "preceded by" means y is to the LEFT (ending before query.left)
    if query.strand == '-'
        end_indices = get(subject.chrom_end_indices, query.chrom, nothing)
        end_indices === nothing && return nothing
        idx = _find_last_end_lt(subject.ends, end_indices, query.left)
        return idx === nothing ? nothing : subject.intervals[idx]
    end

    # For '+' or '.': find the nearest interval starting after query.right
    idx = searchsortedfirst(subject.starts, query.right + 1, first(range), last(range), Base.Order.Forward)
    return idx <= last(range) ? subject.intervals[idx] : nothing
end

"""
    follow(query::GenomicInterval, subject::IntervalCollection)

Strand-aware: find the interval in `subject` that `query` follows.
- For `+` strand (or `.`): `y` is to the left of `x` (ending before query.left).
- For `-` strand: `y` is to the right of `x` (starting after query.right).
"""
function follow(query::GenomicInterval, subject::IntervalCollection)
    range = get(subject.chrom_indices, query.chrom, nothing)
    range === nothing && return nothing

    if query.strand == '-'
        # For '-' strand: find nearest interval starting after query.right
        idx = searchsortedfirst(subject.starts, query.right + 1, first(range), last(range), Base.Order.Forward)
        return idx <= last(range) ? subject.intervals[idx] : nothing
    end

    # For '+' or '.': find the nearest interval ending before query.left
    end_indices = get(subject.chrom_end_indices, query.chrom, nothing)
    end_indices === nothing && return nothing
    idx = _find_last_end_lt(subject.ends, end_indices, query.left)
    return idx === nothing ? nothing : subject.intervals[idx]
end

"""
    _subtract_piece(target, blocker)

Subtract a single blocker interval from a target, returning 0-2 remaining pieces.
"""
function _subtract_piece(target::GenomicInterval, blocker::GenomicInterval)
    target.chrom == blocker.chrom || return [target]
    blocker.right < target.left && return [target]
    blocker.left > target.right && return [target]
    result = GenomicInterval[]
    if blocker.left > target.left
        push!(result, GenomicInterval(target.chrom, target.left, blocker.left - 1, target.strand, target.metadata))
    end
    if blocker.right < target.right
        push!(result, GenomicInterval(target.chrom, blocker.right + 1, target.right, target.strand, target.metadata))
    end
    return result
end

function _subtract_many(piece::GenomicInterval, blockers::AbstractVector{<:GenomicInterval})
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

# ==============================================================================
# Interval subtraction via sweep-line arithmetic
#
# The previous implementation enumerated every individual position into
# Set{Int}, which is O(genome_length) in memory and time — catastrophic
# for megabase-scale intervals. This replacement operates on interval
# boundaries directly, achieving O(n log n) via sorting.
# ==============================================================================

"""
    _subtract_intervals_on_chrom(targets, blockers) -> Vector{GenomicInterval}

Subtract `blockers` from `targets` on a single chromosome using a sweep-line
algorithm. Both inputs must be sorted by left endpoint. Runs in O((m+n) log(m+n)).
"""
function _subtract_intervals_on_chrom(targets::AbstractVector{<:GenomicInterval}, blockers::AbstractVector{<:GenomicInterval})
    isempty(targets) && return GenomicInterval[]
    isempty(blockers) && return copy(targets)

    # Merge all blocker intervals into a non-overlapping sorted set
    merged_blockers = GenomicInterval[]
    current = blockers[1]
    for i in 2:length(blockers)
        b = blockers[i]
        if b.left <= current.right + 1
            current = GenomicInterval(current.chrom, current.left, max(current.right, b.right), '.', Dict{String,Any}())
        else
            push!(merged_blockers, current)
            current = b
        end
    end
    push!(merged_blockers, current)

    result = GenomicInterval[]
    bi = 1  # blocker index
    nb = length(merged_blockers)

    for target in targets
        # Advance past blockers that end before this target starts
        while bi <= nb && merged_blockers[bi].right < target.left
            bi += 1
        end

        cursor = target.left
        j = bi
        while j <= nb && merged_blockers[j].left <= target.right
            blocker = merged_blockers[j]
            if blocker.left > cursor
                # Gap before this blocker: emit [cursor, blocker.left - 1]
                push!(result, GenomicInterval(target.chrom, cursor, blocker.left - 1, target.strand, target.metadata))
            end
            cursor = max(cursor, blocker.right + 1)
            j += 1
        end
        # Remaining tail after all blockers
        if cursor <= target.right
            push!(result, GenomicInterval(target.chrom, cursor, target.right, target.strand, target.metadata))
        end
    end

    return result
end

function Base.setdiff(left::GenomicInterval, right::IntervalCollection)
    reduced = reduce(right)
    chrom = left.chrom
    range = get(reduced.chrom_indices, chrom, nothing)
    range === nothing && return [left]
    blockers = reduced.intervals[range]
    return _subtract_intervals_on_chrom([left], blockers)
end

function Base.setdiff(left::IntervalCollection, right::IntervalCollection)
    reduced_right = reduce(right)
    result = GenomicInterval[]
    for chrom in sort(collect(keys(left.chrom_indices)))
        range_left = left.chrom_indices[chrom]
        targets = left.intervals[range_left]
        range_right = get(reduced_right.chrom_indices, chrom, nothing)
        if range_right === nothing
            append!(result, targets)
        else
            blockers = reduced_right.intervals[range_right]
            append!(result, _subtract_intervals_on_chrom(targets, blockers))
        end
    end
    return build_collection(result)
end

function _merge_metadata(metadata)
    merged = Dict{String,Any}(String(key) => value for (key, value) in pairs(metadata))
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
                _merge_metadata(current.metadata))
        else
            push!(merged, current)
            current = candidate
        end
    end

    push!(merged, current)
    return build_collection(with_provenance(merged, "GenomicInterval", "GenomicRanges/reduce"; notes=["merged overlapping intervals"], parameters=(input_count=length(collection.intervals), output_count=length(merged))))
end

function _interval_collection_intersect(left::IntervalCollection, right::IntervalCollection)
    (isempty(left) || isempty(right)) && return build_collection(GenomicInterval[])

    pieces = GenomicInterval[]
    for left_interval in left.intervals
        for hit in find_overlaps(left_interval, right)
            push!(pieces, GenomicInterval(left_interval.chrom, max(left_interval.left, hit.left), min(left_interval.right, hit.right), left_interval.strand, _copy_metadata(left_interval)))
        end
    end
    return build_collection(with_provenance(pieces, "GenomicInterval", "GenomicRanges/intersect"; parameters=(left_count=length(left.intervals), right_count=length(right.intervals), output_count=length(pieces))))
end

function Base.intersect(left::IntervalCollection, right::IntervalCollection)
    return _interval_collection_intersect(left, right)
end

function Base.union(left::IntervalCollection, right::IntervalCollection)
    return reduce(build_collection(vcat(left.intervals, right.intervals)))
end

function _coverage_segments_for_chrom(chrom::String, intervals::AbstractVector{<:GenomicInterval}, range::UnitRange{Int})
    # Use sorted event vector instead of Dict for better cache performance
    events = Tuple{Int,Int}[]  # (position, delta)
    sizehint!(events, 2 * length(range))

    @inbounds for index in range
        interval = intervals[index]
        push!(events, (interval.left, 1))
        push!(events, (interval.right + 1, -1))
    end

    # Process closing events before opening events at the same coordinate.
    # The grouped-delta pass below preserves closed-interval semantics while
    # this ordering also makes the event stream deterministic.
    sort!(events; by=event -> (event[1], event[2]))

    segments = CoverageSegment[]
    depth = 0
    previous_position = nothing
    i = 1
    while i <= length(events)
        position = events[i][1]
        if previous_position !== nothing && depth > 0 && position > previous_position
            push!(segments, CoverageSegment(chrom, previous_position, position - 1, depth))
        end
        # Sum all deltas at the same position
        while i <= length(events) && events[i][1] == position
            depth += events[i][2]
            i += 1
        end
        previous_position = position
    end

    return segments
end

function find_overlaps_parallel(query::IntervalCollection, subject::IntervalCollection; multi_thread::Bool=true)
    results = Vector{Vector{GenomicInterval}}(undef, length(query))
    chrom_groups = Dict{String,Vector{Int}}()
    for (index, chrom) in enumerate(query.chroms)
        push!(get!(chrom_groups, chrom, Int[]), index)
    end

    chromosomes = collect(keys(chrom_groups))
    if multi_thread && length(chromosomes) > 1 && Threads.nthreads() > 1
        Threads.@threads for chrom_index in eachindex(chromosomes)
            chrom = chromosomes[chrom_index]
            for query_index in chrom_groups[chrom]
                results[query_index] = find_overlaps(query[query_index], subject)
            end
        end
    else
        for chrom in chromosomes
            for query_index in chrom_groups[chrom]
                results[query_index] = find_overlaps(query[query_index], subject)
            end
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

    return with_provenance(result, "GenomicInterval", "GenomicRanges/trim"; parameters=(input_count=length(intervals), output_count=length(result)))
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

    return build_collection(with_provenance(gap_intervals, "GenomicInterval", "GenomicRanges/gaps"; parameters=(input_count=length(collection.intervals), output_count=length(gap_intervals))))
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
        # Use interval tree for O(log n + k) containment checks instead of O(n)
        for index in 1:length(points)-1
            left = points[index]
            right = points[index + 1] - 1
            left > right && continue
            if _collection_overlap_count(collection, chrom, left, right) > 0
                push!(pieces, GenomicInterval(chrom, left, right, '.', Dict{String,Any}()))
            end
        end
    end

    return build_collection(with_provenance(pieces, "GenomicInterval", "GenomicRanges/disjoin"; parameters=(input_count=length(collection.intervals), output_count=length(pieces))))
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

    return with_provenance(results, "GenomicInterval", "GenomicRanges/pintersect"; parameters=(pair_count=length(results),))
end

function punion(x::AbstractVector{<:GenomicInterval}, y::AbstractVector{<:GenomicInterval})
    results = GenomicInterval[]
    for (left, right) in _pairwise_zip(x, y)
        left.chrom == right.chrom || throw(ArgumentError("intervals must be on the same chromosome"))
        push!(results, GenomicInterval(left.chrom, min(left.left, right.left), max(left.right, right.right), left.strand, _copy_metadata(left)))
    end

    return with_provenance(results, "GenomicInterval", "GenomicRanges/punion"; parameters=(pair_count=length(results),))
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
    n = length(rows)
    metadata_keys = Set{String}()
    for interval in rows
        union!(metadata_keys, keys(_comparable_interval_metadata(interval.metadata)))
    end

    # Use typed vectors instead of Any[] for performance
    columns = Dict{Symbol,AbstractVector}(
        :chrom => String[interval.chrom for interval in rows],
        :start => Int[interval.left for interval in rows],
        :stop => Int[interval.right for interval in rows],
        :strand => Char[interval.strand for interval in rows])

    for key in sort(collect(metadata_keys))
        columns[Symbol(key)] = Any[get(interval.metadata, key, missing) for interval in rows]
    end

    return with_provenance(DataFrames.DataFrame(columns), "GenomicIntervalTable", "GenomicRanges/DataFrame"; parameters=(row_count=n,))
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
        chrom_value = row[chrom_name]
        left_value = row[start_name]
        right_value = row[stop_name]
        (ismissing(chrom_value) || ismissing(left_value) || ismissing(right_value)) && continue
        chrom = String(chrom_value)
        left = Int(left_value)
        right = Int(right_value)
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

    return with_provenance(intervals, "GenomicInterval", "GenomicRanges/read_intervals"; notes=["parsed intervals from DataFrame"], parameters=(row_count=nrow(df), output_count=length(intervals)))
end

function coverage(collection::IntervalCollection)
    segments = CoverageSegment[]
    for chrom in sort(collect(keys(collection.chrom_indices)))
        append!(segments, _coverage_segments_for_chrom(chrom, collection.intervals, collection.chrom_indices[chrom]))
    end

    return with_provenance(segments, "CoverageSegment", "GenomicRanges/coverage"; parameters=(segment_count=length(segments), interval_count=length(collection.intervals)))
end

coverage(intervals::AbstractVector{<:GenomicInterval}) = coverage(build_collection(intervals))

function _parse_interval_row(headers::Vector{String}, fields::Vector{String}; column_map::Union{Nothing,Dict{String,Int}}=nothing)
    cmap = column_map === nothing ? Dict(lowercase(strip(header)) => index for (index, header) in pairs(headers)) : column_map

    chrom_index = get(cmap, "chr", get(cmap, "chrom", get(cmap, "chromosome", nothing)))
    start_index = get(cmap, "start", get(cmap, "left", get(cmap, "begin", nothing)))
    stop_index = get(cmap, "end", get(cmap, "stop", get(cmap, "right", nothing)))
    strand_index = get(cmap, "strand", nothing)

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
    column_map = nothing
    intervals = GenomicInterval[]

    for raw_line in eachline(io)
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, '#') && continue

        if headers === nothing
            headers = String[strip(field) for field in Base.split(line, ',')]
            # Pre-compute column map once instead of per-row
            column_map = Dict(lowercase(strip(header)) => index for (index, header) in pairs(headers))
            continue
        end

        fields = String[strip(field) for field in Base.split(line, ',')]
        length(fields) >= length(headers) || throw(ArgumentError("interval row has fewer fields than header"))
        push!(intervals, _parse_interval_row(headers, fields; column_map=column_map))
    end

    headers === nothing && return GenomicInterval[]

    return with_provenance(intervals, "GenomicInterval", "GenomicRanges/read_intervals"; notes=["parsed intervals from delimited text stream"], parameters=(output_count=length(intervals),))
end

function read_intervals(path::String)
    open(path, "r") do io
        intervals = read_intervals(io)
        return with_provenance(intervals, "GenomicInterval", "GenomicRanges/read_intervals"; notes=["parsed intervals from file"], parameters=(path=path, output_count=length(intervals)))
    end
end


# ==============================================================================
# Bioconductor Feature Expansions: GPos, IPos, Transcript & Seqlevels Utilities
# ==============================================================================

"""Construct a 1-bp GenomicInterval (or vector thereof) representing explicit genomic position(s)."""
function GPos(chrom::AbstractString, pos::Integer, strand::Char='.')
    return GenomicInterval(String(chrom), Int(pos), Int(pos), strand, NamedTuple())
end

function GPos(chrom::AbstractString, positions::AbstractVector{<:Integer}, strand::Char='.')
    return [GPos(chrom, p, strand) for p in positions]
end

"""Construct a 1-bp IRange (or vector thereof) representing integer position(s)."""
IPos(pos::Integer) = IRange(pos, 1)
IPos(positions::AbstractVector{<:Integer}) = [IPos(p) for p in positions]

"""Predicate checking if an interval is a 1-bp GPos position."""
is_gpos(interval::GenomicInterval) = width(interval) == 1
is_gpos(intervals::AbstractVector{<:GenomicInterval}) = [is_gpos(i) for i in intervals]

"""Predicate checking if a range is a 1-bp IPos position."""
is_ipos(interval::IRange) = width(interval) == 1
is_ipos(intervals::AbstractVector{<:IRange}) = [is_ipos(i) for i in intervals]

"""
    coverage_by_transcript(coverage_segments, transcripts)

Compute per-base positional coverage arrays across exon models for each transcript.
"""
function coverage_by_transcript(coverage_segments::AbstractVector{CoverageSegment}, transcripts::AbstractDict)
    chrom_cov = Dict{String, Vector{Int}}()
    for seg in coverage_segments
        cov_vec = get!(chrom_cov, seg.chrom, Int[])
        if length(cov_vec) < seg.stop
            append!(cov_vec, zeros(Int, seg.stop - length(cov_vec)))
        end
        for i in seg.start:seg.stop
            cov_vec[i] += seg.depth
        end
    end

    result = Dict{Any, Vector{Int}}()
    for (tx_id, collection) in transcripts
        intervals = collection isa IntervalCollection ? collection.intervals : collect(collection)
        tx_cov = Int[]
        for exon in sort(intervals, by=x -> x.left)
            cov_vec = get(chrom_cov, exon.chrom, Int[])
            for p in exon.left:exon.right
                depth = (p >= 1 && p <= length(cov_vec)) ? cov_vec[p] : 0
                push!(tx_cov, depth)
            end
        end
        result[tx_id] = tx_cov
    end
    return result
end

function coverage_by_transcript(chrom_cov::AbstractDict{<:AbstractString, <:AbstractVector{<:Integer}}, transcripts::AbstractDict)
    result = Dict{Any, Vector{Int}}()
    for (tx_id, collection) in transcripts
        intervals = collection isa IntervalCollection ? collection.intervals : collect(collection)
        tx_cov = Int[]
        for exon in sort(intervals, by=x -> x.left)
            cov_vec = get(chrom_cov, exon.chrom, Int[])
            for p in exon.left:exon.right
                depth = (p >= 1 && p <= length(cov_vec)) ? Int(cov_vec[p]) : 0
                push!(tx_cov, depth)
            end
        end
        result[tx_id] = tx_cov
    end
    return result
end

const coverageByTranscript = coverage_by_transcript

"""
    extract_transcript_seqs(genome_sequences, transcripts)

Extract and splice exonic sequences for each transcript according to its chromosomal coordinates and strand.
"""
function extract_transcript_seqs(genome_sequences::AbstractDict, transcripts::AbstractDict)
    result = Dict{Any, String}()
    for (tx_id, collection) in transcripts
        intervals = collection isa IntervalCollection ? collection.intervals : collect(collection)
        sorted_exons = sort(intervals, by=x -> x.left)
        isempty(sorted_exons) && continue
        
        tx_strand = sorted_exons[1].strand
        chrom = sorted_exons[1].chrom
        
        haskey(genome_sequences, chrom) || continue
        gen_seq_str = string(genome_sequences[chrom])
        
        parts = String[]
        for exon in sorted_exons
            start_pos = max(1, exon.left)
            stop_pos = min(length(gen_seq_str), exon.right)
            start_pos <= stop_pos && push!(parts, gen_seq_str[start_pos:stop_pos])
        end
        
        spliced = join(parts)
        if tx_strand == '-'
            spliced = _reverse_complement_string(spliced)
        end
        result[tx_id] = spliced
    end
    return result
end

function _reverse_complement_string(s::String)
    buf = Char[]
    sizehint!(buf, length(s))
    for c in reverse(s)
        comp = c == 'A' || c == 'a' ? 'T' :
               c == 'T' || c == 't' || c == 'U' || c == 'u' ? 'A' :
               c == 'C' || c == 'c' ? 'G' :
               c == 'G' || c == 'g' ? 'C' : c
        push!(buf, isuppercase(c) ? uppercase(comp) : lowercase(comp))
    end
    return String(buf)
end

const extractTranscriptSeqs = extract_transcript_seqs

"""
    extract_upstream_seqs(genome_sequences, transcripts; width=1000)

Extract upstream sequence of given width (promoter region) for each transcript.
"""
function extract_upstream_seqs(genome_sequences::AbstractDict, transcripts::AbstractDict; width::Integer=1000)
    width > 0 || throw(ArgumentError("width must be positive"))
    result = Dict{Any, String}()
    for (tx_id, collection) in transcripts
        intervals = collection isa IntervalCollection ? collection.intervals : collect(collection)
        sorted_exons = sort(intervals, by=x -> x.left)
        isempty(sorted_exons) && continue
        
        tx_strand = sorted_exons[1].strand
        chrom = sorted_exons[1].chrom
        
        haskey(genome_sequences, chrom) || continue
        gen_seq_str = string(genome_sequences[chrom])
        
        upstream_str = if tx_strand == '-'
            tss = sorted_exons[end].right
            start_pos = tss + 1
            stop_pos = min(length(gen_seq_str), tss + Int(width))
            start_pos <= stop_pos ? _reverse_complement_string(gen_seq_str[start_pos:stop_pos]) : ""
        else
            tss = sorted_exons[1].left
            start_pos = max(1, tss - Int(width))
            stop_pos = max(0, tss - 1)
            start_pos <= stop_pos ? gen_seq_str[start_pos:stop_pos] : ""
        end
        result[tx_id] = upstream_str
    end
    return result
end

const extractUpstreamSeqs = extract_upstream_seqs

"""
    extend_exons_into_introns(transcripts, extension)

Extend internal exon boundaries of each transcript by `extension` base pairs.
"""
function extend_exons_into_introns(transcripts::AbstractDict, extension::Integer)
    extension >= 0 || throw(ArgumentError("extension must be non-negative"))
    result = Dict{Any, Vector{GenomicInterval}}()
    for (tx_id, collection) in transcripts
        intervals = collection isa IntervalCollection ? collection.intervals : collect(collection)
        sorted = sort(intervals, by=x -> x.left)
        n = length(sorted)
        new_exons = GenomicInterval[]
        for i in 1:n
            exon = sorted[i]
            left = (i > 1) ? max(1, exon.left - Int(extension)) : exon.left
            right = (i < n) ? exon.right + Int(extension) : exon.right
            push!(new_exons, GenomicInterval(exon.chrom, left, right, exon.strand, _copy_metadata(exon)))
        end
        result[tx_id] = new_exons
    end
    return result
end

const extendExonsIntoIntrons = extend_exons_into_introns

"""
    map_seqlevels_style(levels, target_style::Symbol)

Map chromosome level names to `:UCSC` ("chr1"), `:Ensembl` ("1"), or `:NCBI` ("NC_000001.11").
"""
function map_seqlevels_style(levels::AbstractVector, target_style::Symbol)
    target_style in (:UCSC, :Ensembl, :NCBI) || throw(ArgumentError("target_style must be :UCSC, :Ensembl, or :NCBI"))
    result = String[]
    for level in levels
        str_level = String(level)
        clean = replace(str_level, r"^chr"i => "")
        mapped = if target_style == :UCSC
            startswith(str_level, "chr") ? str_level : "chr" * clean
        elseif target_style == :Ensembl
            clean
        else
            startswith(str_level, "NC_") ? str_level : "NC_" * clean
        end
        push!(result, mapped)
    end
    return result
end

const mapSeqlevelsStyle = map_seqlevels_style

"""
    zoom(interval, factor)

Expand (factor > 1.0) or shrink (factor < 1.0) interval width around its center.
"""
function zoom(interval::GenomicInterval, factor::Real)
    factor > 0 || throw(ArgumentError("zoom factor must be positive"))
    w = width(interval)
    new_w = max(1, round(Int, w * factor))
    m = mid(interval)
    half = (new_w - 1) ÷ 2
    new_left = max(1, m - half)
    new_right = new_left + new_w - 1
    return GenomicInterval(interval.chrom, new_left, new_right, interval.strand, _copy_metadata(interval))
end

zoom(intervals::AbstractVector{<:GenomicInterval}, factor::Real) = [zoom(i, factor) for i in intervals]
zoom(collection::IntervalCollection, factor::Real) = build_collection(zoom(collection.intervals, factor))

zoom_in(x, factor::Real=2.0) = zoom(x, 1.0 / factor)
zoom_out(x, factor::Real=2.0) = zoom(x, factor)

"""Return element lengths of grouped interval collections in a dictionary."""
element_lengths(grouped::AbstractDict) = Dict(k => length(v) for (k, v) in grouped)
const elementLengths = element_lengths

"""Concatenate multiple IntervalCollection instances into a single IntervalCollection."""
function Base.vcat(collections::IntervalCollection...)
    all_intervals = reduce(vcat, [c.intervals for c in collections]; init=GenomicInterval[])
    return build_collection(all_intervals)
end

end
