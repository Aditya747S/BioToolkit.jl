using Statistics
using Random

export GenomeViewport, AbstractTrack, GeneTrack, CoverageTrack, AlignmentTrack, GenomeBrowser
export GeneSegment, GenePlacement, GeneRenderPlan, CoverageBin, CoverageRenderPlan
export ReadMismatch, ReadPlacement, ReadConnector, AlignmentRenderPlan
export genome_lod, render_track, render_browser, export_figure

const _GENOME_BROWSER_CHROMOSOME_LOD = 10_000.0
const _GENOME_BROWSER_GENE_LOD = 20.0
const _GENOME_BROWSER_DEFAULT_PIXEL_WIDTH = 1200

abstract type AbstractTrack end

struct GenomeViewport
    chrom::String
    range::UnitRange{Int}
    pixel_width::Int
    resolution::Float64

    function GenomeViewport(chrom::AbstractString, range::UnitRange{<:Integer}, pixel_width::Integer)
        pixel_width > 0 || throw(ArgumentError("pixel_width must be positive"))
        first_position = Int(first(range))
        last_position = Int(last(range))
        last_position >= first_position || throw(ArgumentError("range must be non-empty"))
        resolution = (last_position - first_position) / Float64(pixel_width)
        return new(String(chrom), first_position:last_position, Int(pixel_width), Float64(resolution))
    end
end

GenomeViewport(chrom::AbstractString, start::Integer, stop::Integer; width::Integer=_GENOME_BROWSER_DEFAULT_PIXEL_WIDTH, pixel_width::Union{Nothing,Integer}=nothing) = GenomeViewport(chrom, Int(start):Int(stop), Int(pixel_width === nothing ? width : pixel_width))
GenomeViewport(chrom::AbstractString, start::Integer, stop::Integer, pixel_width::Integer) = GenomeViewport(chrom, Int(start):Int(stop), Int(pixel_width))

struct GeneTrack <: AbstractTrack
    intervals::Vector{GenomicInterval}
    style::Symbol
    height::Int
    color::Any
    label_field::Union{Nothing,Symbol}
    max_labels::Int
end

struct CoverageTrack <: AbstractTrack
    source::Any
    chrom::Union{Nothing,String}
    height::Int
    color::Any
    style::Symbol
    window_size::Union{Nothing,Int}
    sorted::Bool
    start::Int
end

struct AlignmentTrack <: AbstractTrack
    source::Any
    show_mismatches::Bool
    height::Int
    color::Any
    max_reads::Int
    rng_seed::Int
    show_pairs::Bool
end

struct GenomeBrowser
    tracks::Vector{AbstractTrack}
    viewport::GenomeViewport
    gap::Int
    title::String
end

struct GeneSegment
    kind::Symbol
    left::Int
    right::Int
end

struct GenePlacement
    interval::GenomicInterval
    row::Int
    label::Union{Nothing,String}
    segments::Vector{GeneSegment}
end

struct GeneRenderPlan
    lod::Symbol
    placements::Vector{GenePlacement}
    density::Vector{Float64}
    bin_edges::Vector{Int}
    row_count::Int
    labels::Vector{String}
end

struct CoverageBin
    left::Int
    right::Int
    mean_value::Float64
    min_value::Float64
    max_value::Float64
end

struct CoverageRenderPlan
    lod::Symbol
    bins::Vector{CoverageBin}
    values::Vector{Float64}
    bin_size::Int
end

struct ReadMismatch
    position::Int
    reference_base::Union{Missing,Char}
    read_base::Union{Missing,Char}
end

struct ReadPlacement
    qname::String
    row::Int
    left::Int
    right::Int
    blocks::Vector{Tuple{Int,Int}}
    mismatches::Vector{ReadMismatch}
    paired::Bool
    strand::Char
end

struct ReadConnector
    qname::String
    left::Int
    right::Int
    mate_left::Int
    mate_right::Int
end

struct AlignmentRenderPlan
    lod::Symbol
    reads::Vector{ReadPlacement}
    connectors::Vector{ReadConnector}
    downsampled::Bool
end

Base.length(browser::GenomeBrowser) = length(browser.tracks)
Base.iterate(browser::GenomeBrowser, state...) = iterate(browser.tracks, state...)

function Base.show(io::IO, viewport::GenomeViewport)
    print(io, "GenomeViewport(", viewport.chrom, ":", first(viewport.range), "-", last(viewport.range), ", pixel_width=", viewport.pixel_width, ", resolution=", round(viewport.resolution; digits=3), ")")
end

function Base.show(io::IO, browser::GenomeBrowser)
    print(io, "GenomeBrowser(", length(browser.tracks), " tracks, ", browser.viewport.chrom, ":", first(browser.viewport.range), "-", last(browser.viewport.range), ", title=", repr(browser.title), ")")
end

function Base.show(io::IO, ::MIME"text/plain", viewport::GenomeViewport)
    show(io, viewport)
end

function Base.show(io::IO, ::MIME"text/plain", browser::GenomeBrowser)
    show(io, browser)
end

GenomeBrowser(tracks::AbstractVector{<:AbstractTrack}, viewport::GenomeViewport; gap::Integer=12, title::AbstractString="Genome Browser") = GenomeBrowser(AbstractTrack[track for track in tracks], viewport, Int(gap), String(title))
GenomeBrowser(tracks::AbstractVector{<:AbstractTrack}, chrom::AbstractString, start::Integer, stop::Integer; width::Integer=_GENOME_BROWSER_DEFAULT_PIXEL_WIDTH, gap::Integer=12, title::AbstractString="Genome Browser") = GenomeBrowser(tracks, GenomeViewport(chrom, start, stop; width=width); gap=gap, title=title)

GeneTrack(intervals::AbstractVector{<:GenomicInterval}; style::Symbol=:squish, height::Integer=100, color=:black, label_field::Union{Nothing,Symbol}=nothing, max_labels::Integer=12) = GeneTrack(GenomicInterval[interval for interval in intervals], style, Int(height), color, label_field, Int(max_labels))
GeneTrack(collection::IntervalCollection; kwargs...) = GeneTrack(collection.intervals; kwargs...)

CoverageTrack(source; chrom::Union{Nothing,AbstractString}=nothing, height::Integer=150, color=:blue, style::Symbol=:line, window_size::Union{Nothing,Integer}=nothing, sorted::Bool=false, start::Integer=1) = CoverageTrack(source, chrom === nothing ? nothing : String(chrom), Int(height), color, style, window_size === nothing ? nothing : Int(window_size), sorted, Int(start))

AlignmentTrack(source; show_mismatches::Bool=true, height::Integer=300, color=:steelblue, max_reads::Integer=500, rng_seed::Integer=0, show_pairs::Bool=true) = AlignmentTrack(source, show_mismatches, Int(height), color, Int(max_reads), Int(rng_seed), show_pairs)

genome_lod(viewport::GenomeViewport) = viewport.resolution > _GENOME_BROWSER_CHROMOSOME_LOD ? :chromosome : viewport.resolution > _GENOME_BROWSER_GENE_LOD ? :gene : :basepair

function _viewport_start(viewport::GenomeViewport)
    return first(viewport.range)
end

function _viewport_stop(viewport::GenomeViewport)
    return last(viewport.range)
end

function _viewport_span(viewport::GenomeViewport)
    return _viewport_stop(viewport) - _viewport_start(viewport) + 1
end

function _visible_interval(interval::GenomicInterval, viewport::GenomeViewport)
    interval.chrom == viewport.chrom || return false
    interval.right < _viewport_start(viewport) && return false
    interval.left > _viewport_stop(viewport) && return false
    return true
end

function _visible_intervals(intervals::AbstractVector{<:GenomicInterval}, viewport::GenomeViewport)
    visible = [interval for interval in intervals if _visible_interval(interval, viewport)]
    sort!(visible, by = interval -> (interval.left, interval.right, interval.strand))
    return visible
end

function _metadata_value(metadata::AbstractDict, key::AbstractString, default=nothing)
    haskey(metadata, key) && return metadata[key]
    symbol_key = Symbol(key)
    haskey(metadata, symbol_key) && return metadata[symbol_key]
    return default
end

function _gene_label(interval::GenomicInterval, label_field::Union{Nothing,Symbol})
    if label_field !== nothing
        value = _metadata_value(interval.metadata, String(label_field), nothing)
        value !== nothing && return String(value)
    end
    for key in ("name", "gene_id", "gene", "ID")
        value = _metadata_value(interval.metadata, key, nothing)
        value !== nothing && return String(value)
    end
    return nothing
end

function _coerce_segment(value, kind::Symbol)
    if value isa GenomicInterval
        return GeneSegment(kind, value.left, value.right)
    elseif value isa AbstractRange
        return GeneSegment(kind, Int(first(value)), Int(last(value)))
    elseif value isa Tuple && length(value) >= 2
        return GeneSegment(kind, Int(value[1]), Int(value[2]))
    elseif value isa Pair
        return GeneSegment(kind, Int(first(value)), Int(last(value)))
    else
        throw(ArgumentError("unsupported segment type $(typeof(value))"))
    end
end

function _gene_segments(interval::GenomicInterval, lod::Symbol)
    if lod != :basepair
        return GeneSegment[GeneSegment(:gene, interval.left, interval.right)]
    end

    segments = GeneSegment[]
    for (key, kind) in (("exons", :exon), ("cds", :cds), ("utr5", :utr5), ("utr3", :utr3))
        values = _metadata_value(interval.metadata, key, nothing)
        values === nothing && continue
        if values isa AbstractVector
            append!(segments, (_coerce_segment(value, kind) for value in values))
        else
            push!(segments, _coerce_segment(values, kind))
        end
    end

    isempty(segments) && push!(segments, GeneSegment(:gene, interval.left, interval.right))
    sort!(segments, by = segment -> (segment.left, segment.right))
    return segments
end

function _pack_genes(intervals::AbstractVector{<:GenomicInterval}, viewport::GenomeViewport, lod::Symbol, label_field::Union{Nothing,Symbol}, max_labels::Integer)
    placements = GenePlacement[]
    row_ends = Int[]
    labels = String[]

    for interval in _visible_intervals(intervals, viewport)
        row = 1
        while row <= length(row_ends) && interval.left <= row_ends[row]
            row += 1
        end
        if row > length(row_ends)
            push!(row_ends, interval.right)
        else
            row_ends[row] = interval.right
        end

        label = _gene_label(interval, label_field)
        label !== nothing && length(labels) < Int(max_labels) && push!(labels, label)
        push!(placements, GenePlacement(interval, row, label, _gene_segments(interval, lod)))
    end

    return placements, labels, length(row_ends)
end

function _gene_density(intervals::AbstractVector{<:GenomicInterval}, viewport::GenomeViewport, bin_size::Int)
    span = _viewport_span(viewport)
    bins = max(1, cld(span, bin_size))
    density = zeros(Float64, bins)
    start_position = _viewport_start(viewport)
    stop_position = _viewport_stop(viewport)

    for interval in _visible_intervals(intervals, viewport)
        left = max(interval.left, start_position) - start_position
        right = min(interval.right, stop_position) - start_position
        first_bin = fld(left, bin_size) + 1
        last_bin = fld(right, bin_size) + 1
        for bin in first_bin:last_bin
            density[bin] += 1.0
        end
    end

    bin_edges = [start_position + (bin_index - 1) * bin_size for bin_index in 1:bins]
    return density, bin_edges
end

function _dense_coverage(source::AbstractVector{<:Real}, start::Int)
    dense = Float64.(source)
    return dense, start
end

function _dense_coverage(source::SparseCoverageVector, start::Int)
    dense = zeros(Float64, max(source.span, maximum(source.positions; init=0)))
    for (position, depth) in zip(source.positions, source.depths)
        1 <= position <= length(dense) || continue
        dense[position] = Float64(depth)
    end
    return dense, start
end

function _coverage_source(track::CoverageTrack, viewport::GenomeViewport)
    if track.source isa AbstractDict
        chrom = track.chrom === nothing ? viewport.chrom : track.chrom
        haskey(track.source, chrom) && return track.source[chrom]
        for key in keys(track.source)
            String(key) == chrom && return track.source[key]
        end
        return nothing
    end
    return track.source
end

function _coverage_bin_size(track::CoverageTrack, viewport::GenomeViewport, lod::Symbol)
    if track.window_size !== nothing
        return max(1, track.window_size)
    end
    if lod == :chromosome
        return max(1, ceil(Int, viewport.resolution) * 10)
    elseif lod == :gene
        return max(1, ceil(Int, viewport.resolution))
    else
        return 1
    end
end

function _coverage_bins_from_vector(values::AbstractVector{<:Real}, start::Int, viewport::GenomeViewport, bin_size::Int)
    bins = CoverageBin[]
    viewport_start = _viewport_start(viewport)
    viewport_stop = _viewport_stop(viewport)
    source_start = start
    source_stop = start + length(values) - 1
    overlap_start = max(viewport_start, source_start)
    overlap_stop = min(viewport_stop, source_stop)
    overlap_start > overlap_stop && return bins, Float64[]

    first_bin_start = overlap_start
    while first_bin_start <= overlap_stop
        bin_stop = min(first_bin_start + bin_size - 1, overlap_stop)
        local_start = first_bin_start - source_start + 1
        local_stop = bin_stop - source_start + 1
        slice = Float64.(values[local_start:local_stop])
        push!(bins, CoverageBin(first_bin_start, bin_stop, mean(slice), minimum(slice), maximum(slice)))
        first_bin_start = bin_stop + 1
    end

    return bins, [bin.mean_value for bin in bins]
end

function _coverage_bins_from_table(table, viewport::GenomeViewport, bin_size::Int; sorted::Bool=false)
    coverage = BioToolkit.window_coverage(table, viewport.chrom, bin_size; sorted=sorted)
    bins = CoverageBin[]
    for bin_index in sort(collect(keys(coverage)))
        left = _viewport_start(viewport) + bin_index * bin_size
        right = min(left + bin_size - 1, _viewport_stop(viewport))
        value = Float64(coverage[bin_index]) / Float64(bin_size)
        push!(bins, CoverageBin(left, right, value, value, value))
    end
    return bins, [bin.mean_value for bin in bins]
end

function _coverage_plan(track::CoverageTrack, viewport::GenomeViewport, lod::Symbol)
    source = _coverage_source(track, viewport)
    source === nothing && return CoverageRenderPlan(lod, CoverageBin[], Float64[], _coverage_bin_size(track, viewport, lod))

    bin_size = _coverage_bin_size(track, viewport, lod)
    if Tables.istable(source)
        bins, values = _coverage_bins_from_table(source, viewport, bin_size; sorted=track.sorted)
        return CoverageRenderPlan(lod, bins, values, bin_size)
    elseif source isa SparseCoverageVector
        values, start = _dense_coverage(source, track.start)
        bins, summarized = _coverage_bins_from_vector(values, start, viewport, bin_size)
        return CoverageRenderPlan(lod, bins, summarized, bin_size)
    elseif source isa AbstractVector{<:Real}
        values, start = _dense_coverage(source, track.start)
        bins, summarized = _coverage_bins_from_vector(values, start, viewport, bin_size)
        return CoverageRenderPlan(lod, bins, summarized, bin_size)
    elseif source isa AbstractDict
        chrom = track.chrom === nothing ? viewport.chrom : track.chrom
        haskey(source, chrom) || return CoverageRenderPlan(lod, CoverageBin[], Float64[], bin_size)
        nested = source[chrom]
        nested isa SparseCoverageVector && return _coverage_plan(CoverageTrack(nested; chrom=chrom, height=track.height, color=track.color, style=track.style, window_size=track.window_size, sorted=track.sorted, start=track.start), viewport, lod)
        nested isa AbstractVector{<:Real} && return _coverage_plan(CoverageTrack(nested; chrom=chrom, height=track.height, color=track.color, style=track.style, window_size=track.window_size, sorted=track.sorted, start=track.start), viewport, lod)
        Tables.istable(nested) && return _coverage_plan(CoverageTrack(nested; chrom=chrom, height=track.height, color=track.color, style=track.style, window_size=track.window_size, sorted=track.sorted, start=track.start), viewport, lod)
        throw(ArgumentError("unsupported coverage source in dictionary for chromosome $(chrom)"))
    else
        throw(ArgumentError("unsupported coverage source $(typeof(source))"))
    end
end

function _read_reference_span(record::BamRecord)
    return Int(record.pos) + 1, Int(record.pos) + _bam_reference_span(record.cigar)
end

function _read_blocks(record::BamRecord)
    blocks = Tuple{Int,Int}[]
    cursor = Int(record.pos) + 1
    for op in record.cigar
        if op.op in ('M', '=', 'X')
            push!(blocks, (cursor, cursor + op.length - 1))
            cursor += op.length
        elseif op.op in ('D', 'N')
            cursor += op.length
        end
    end
    return blocks
end

function _parse_md_number(md::AbstractString, index::Int)
    cursor = index
    while cursor <= lastindex(md) && isdigit(md[cursor])
        cursor += 1
    end
    return parse(Int, md[index:cursor - 1]), cursor
end

function _read_mismatches(record::BamRecord)
    md = get(record.tags, "MD", nothing)
    md isa AbstractString || return ReadMismatch[]

    mismatches = ReadMismatch[]
    reference_cursor = Int(record.pos) + 1
    index = firstindex(md)
    while index <= lastindex(md)
        character = md[index]
        if isdigit(character)
            consumed, index = _parse_md_number(md, index)
            reference_cursor += consumed
        elseif character == '^'
            index += 1
            while index <= lastindex(md) && !isdigit(md[index])
                reference_cursor += 1
                index += 1
            end
        else
            push!(mismatches, ReadMismatch(reference_cursor, character, missing))
            reference_cursor += 1
            index += 1
        end
    end

    return mismatches
end

function _read_connector(record::BamRecord)
    record.mate_refname === nothing && return nothing
    record.mate_refname != record.refname && return nothing
    record.mate_pos < 0 && return nothing
    mate_left = Int(record.mate_pos) + 1
    mate_right = mate_left + max(0, abs(Int(record.template_length)))
    left, right = _read_reference_span(record)
    return ReadConnector(record.qname, left, right, mate_left, mate_right)
end

function _read_records(source, viewport::GenomeViewport)
    region = GenomicInterval(viewport.chrom, _viewport_start(viewport), _viewport_stop(viewport))
    if source isa AbstractString
        return read_bam(source, region).records
    elseif source isa BamFile
        return [record for record in source.records if _bam_overlaps(record, region)]
    elseif source isa AbstractVector{<:BamRecord}
        return [record for record in source if _bam_overlaps(record, region)]
    else
        throw(ArgumentError("unsupported alignment source $(typeof(source))"))
    end
end

function _reservoir_sample(records::Vector{BamRecord}, max_reads::Int, rng::AbstractRNG)
    length(records) <= max_reads && return records, false
    sampled = records[1:max_reads]
    for index in max_reads+1:length(records)
        chosen = rand(rng, 1:index)
        if chosen <= max_reads
            sampled[chosen] = records[index]
        end
    end
    return sampled, true
end

function _pack_reads(records::Vector{BamRecord})
    placements = ReadPlacement[]
    connectors = ReadConnector[]
    row_ends = Int[]

    for record in records
        left, right = _read_reference_span(record)
        row = 1
        while row <= length(row_ends) && left <= row_ends[row]
            row += 1
        end
        if row > length(row_ends)
            push!(row_ends, right)
        else
            row_ends[row] = right
        end

        blocks = _read_blocks(record)
        mismatches = _read_mismatches(record)
        paired = record.mate_refname !== nothing && record.mate_refname == record.refname && record.mate_pos >= 0
        paired && push!(connectors, _read_connector(record))
        strand = (record.flag & 0x10) == 0x10 ? '-' : '+'
        push!(placements, ReadPlacement(record.qname, row, left, right, blocks, mismatches, paired, strand))
    end

    return placements, connectors
end

function _alignment_plan(track::AlignmentTrack, viewport::GenomeViewport, lod::Symbol)
    records = _read_records(track.source, viewport)
    sampled, downsampled = _reservoir_sample(records, track.max_reads, MersenneTwister(track.rng_seed))
    placements, connectors = _pack_reads(sampled)
    if !track.show_mismatches
        placements = [ReadPlacement(read.qname, read.row, read.left, read.right, read.blocks, ReadMismatch[], read.paired, read.strand) for read in placements]
    end
    if !track.show_pairs
        connectors = ReadConnector[]
    end
    return AlignmentRenderPlan(lod, placements, connectors, downsampled)
end

function render_track(track::GeneTrack, viewport::GenomeViewport)
    lod = genome_lod(viewport)
    if lod == :chromosome
        bin_size = max(1, ceil(Int, viewport.resolution) * 10)
        density, bin_edges = _gene_density(track.intervals, viewport, bin_size)
        return GeneRenderPlan(lod, GenePlacement[], density, bin_edges, 0, String[])
    end

    placements, labels, row_count = _pack_genes(track.intervals, viewport, lod, track.label_field, track.max_labels)
    return GeneRenderPlan(lod, placements, Float64[], Int[], row_count, labels)
end

function render_track(track::CoverageTrack, viewport::GenomeViewport)
    lod = genome_lod(viewport)
    return _coverage_plan(track, viewport, lod)
end

function render_track(track::AlignmentTrack, viewport::GenomeViewport)
    lod = genome_lod(viewport)
    return _alignment_plan(track, viewport, lod)
end

render_track(track::AbstractTrack, viewport::GenomeViewport) = throw(ArgumentError("unsupported track type $(typeof(track))"))

function render_browser(browser::GenomeBrowser)
    return [(track = track, plan = render_track(track, browser.viewport)) for track in browser.tracks]
end

function _render_axis!(backend::Module, axis, plan::GeneRenderPlan, track::GeneTrack)
    if plan.lod == :chromosome
        scatter = getfield(backend, :scatter!)
        line = getfield(backend, :lines!)
        x = collect(plan.bin_edges)
        y = plan.density
        isempty(x) || scatter(axis, x, y; color=track.color, markersize=8)
        isempty(x) || line(axis, x, y; color=track.color)
    else
        line = getfield(backend, :lines!)
        text = getfield(backend, :text!)
        for placement in plan.placements
            for segment in placement.segments
                line(axis, [segment.left, segment.right], fill(placement.row, 2); color=track.color, linewidth=2)
            end
            placement.label === nothing || text(axis, [placement.interval.left], [placement.row + 0.2]; text=[placement.label], color=track.color, fontsize=10)
        end
    end
    return axis
end

function _render_axis!(backend::Module, axis, plan::CoverageRenderPlan, track::CoverageTrack)
    line = getfield(backend, :lines!)
    scatter = getfield(backend, :scatter!)
    if isempty(plan.bins)
        return axis
    end
    x = Float64[]
    y = Float64[]
    for bin in plan.bins
        push!(x, bin.left)
        push!(x, bin.right)
        push!(y, bin.mean_value)
        push!(y, bin.mean_value)
    end
    line(axis, x, y; color=track.color)
    scatter(axis, [bin.left for bin in plan.bins], [bin.mean_value for bin in plan.bins]; color=track.color, markersize=4)
    return axis
end

function _render_axis!(backend::Module, axis, plan::AlignmentRenderPlan, track::AlignmentTrack)
    line = getfield(backend, :lines!)
    scatter = getfield(backend, :scatter!)
    text = getfield(backend, :text!)
    for read in plan.reads
        for (left, right) in read.blocks
            line(axis, [left, right], fill(read.row, 2); color=track.color, linewidth=2)
        end
        for mismatch in read.mismatches
            scatter(axis, [mismatch.position], [read.row]; color=:red, markersize=10)
        end
        text(axis, [read.left], [read.row + 0.2]; text=[read.qname], color=track.color, fontsize=9)
    end
    for connector in plan.connectors
        line(axis, [connector.left, connector.mate_left], [1, 1]; color=:gray, linewidth=1)
    end
    return axis
end

function render!(scene, track::AbstractTrack, viewport::GenomeViewport)
    throw(ArgumentError("Makie extension required to render $(typeof(track))"))
end

function export_figure(browser::GenomeBrowser, path::AbstractString; dpi::Integer=300)
    throw(ArgumentError("CairoMakie extension required to export browser figures"))
end