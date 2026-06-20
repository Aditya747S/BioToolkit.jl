# ==============================================================================
# epigenetics.jl — Epigenetics and chromatin analysis
#
# Provides coverage computation, peak calling, methylation analysis,
# single-cell chromatin (TF-IDF/LSI), gene activity scoring,
# co-accessibility, motif deviations, and TAD detection.
#
# References:
#   - Stuart et al. (2019) Nat Biotechnol 37:1333-1341 (Signac)
#   - Dixon et al. (2012) Nature 485:376-380 (TADs)
# ==============================================================================

module Epigenetics

using SparseArrays
using DataFrames
using Statistics
using LinearAlgebra
using Random
using Distributions
using SpecialFunctions

using ..GenomicRanges: GenomicInterval, IntervalCollection, CoverageSegment, build_collection, coverage as interval_coverage, promoters
using ..DifferentialExpression: CountMatrix, DEResult, differential_expression, benjamini_hochberg
using ..SingleCell: SingleCellExperiment, count_matrix as singlecell_count_matrix, normalize_counts, run_pca, run_umap, cluster_cells, summarize_clusters
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, BioSequence, DNAAlphabet, HMM, viterbi
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, ProvenanceContext, ThreadSafeProvenanceContext, active_provenance_context, metadata_provenance, new_provenance_id, provenance_parent_ids, provenance_record, register_provenance!, stamp_provenance!, update_provenance!

export SparseCoverageVector, Peak, PeakSet, Epigenome
export PeakSupport
export MethylationCall, MethylationExperiment, MethylationResult
export SingleCellChromatinExperiment, CoaccessibilityEdge, TadResult, ContactMatrix
export calculate_coverage, coverage_depth, coverage_segments, call_peaks, normalize_gc_bias
export count_overlaps, differential_binding, summarize_peak_support
export bin_methylation, differential_methylation
export tfidf, run_lsi, rsvd, gene_activity_score, calculate_coaccessibility
export compute_motif_deviations, detect_footprints
export directionality_index, detect_tads
export chromhmm_like_segmentation, hic_ice_normalize, hic_kr_normalize, hic_ab_compartments, tobias_like_footprints, parse_bismark_coverage
export insulation_score, compartment_boundary_candidates

struct SparseCoverageVector
    positions::Vector{Int}
    depths::Vector{Int}
    span::Int

    function SparseCoverageVector(positions::AbstractVector{<:Integer}, depths::AbstractVector{<:Integer}, span::Integer)
        length(positions) == length(depths) || throw(ArgumentError("positions and depths must have the same length"))
        new(Int.(positions), Int.(depths), Int(span))
    end
end

struct Peak
    id::String
    chrom::String
    left::Int
    right::Int
    summit::Int
    score::Float64
    pvalue::Float64
    qvalue::Float64
end

struct PeakSet
    peaks::Vector{Peak}
    chrom_index::Dict{String,Vector{Int}}
    metadata::Dict{String,Any}
end

struct PeakSupport
    peak_id::String
    sample_id::String
    fragment_count::Int
    total_overlap::Int
    mean_overlap::Float64
    enrichment::Float64
end

struct Epigenome
    intervals::IntervalCollection
    coverage::Dict{String,SparseCoverageVector}
    counts::SparseMatrixCSC{Int,Int}
    sample_metadata::DataFrame
end

struct MethylationCall
    chrom::String
    position::Int
    sample::String
    methylated::Bool
end

struct MethylationExperiment
    methylated::SparseMatrixCSC{Int,Int}
    total::SparseMatrixCSC{Int,Int}
    regions::IntervalCollection
    sample_metadata::DataFrame
    sample_ids::Vector{String}
    metadata::Dict{String,Any}
end

struct MethylationResult <: AbstractAnalysisResult
    region_id::String
    chrom::String
    left::Int
    right::Int
    group1_mean::Float64
    group2_mean::Float64
    delta_methylation::Float64
    stat::Float64
    pvalue::Float64
    padj::Float64
    provenance::ResultProvenance
end

MethylationResult(region_id, chrom, left, right, group1_mean, group2_mean, delta_methylation, stat, pvalue, padj) =
    MethylationResult(region_id, chrom, left, right, group1_mean, group2_mean, delta_methylation, stat, pvalue, padj, provenance_record("MethylationResult", "epigenetics"))

struct SingleCellChromatinExperiment
    base::SingleCellExperiment
    open_peaks::SparseMatrixCSC{Int,Int}
    peak_ids::Vector{String}
    reductions::Dict{String,Matrix{Float64}}
    metadata::Dict{String,Any}
end

struct CoaccessibilityEdge
    peak1::String
    peak2::String
    correlation::Float64
end

struct TadResult <: AbstractAnalysisResult
    chrom::String
    left_bin::Int
    right_bin::Int
    score::Float64
    state::Symbol
    provenance::ResultProvenance
end

TadResult(chrom, left_bin, right_bin, score, state) =
    TadResult(chrom, left_bin, right_bin, score, state, provenance_record("TadResult", "epigenetics"))

struct ContactMatrix
    matrix::SparseMatrixCSC{Int,Int}
    bins::Vector{GenomicInterval}
    bin_size::Int
    chrom::String
end

Epigenome(intervals::AbstractVector{<:GenomicInterval}, coverage::Dict{String,SparseCoverageVector}, counts::SparseMatrixCSC{Int,Int}, sample_metadata::DataFrame) = Epigenome(build_collection(intervals), coverage, counts, sample_metadata)

Base.length(set::PeakSet) = length(set.peaks)
Base.isempty(set::PeakSet) = isempty(set.peaks)

function PeakSet(peaks::AbstractVector{<:Peak}; metadata::AbstractDict=Dict{String,Any}())
    sorted = sort!(collect(peaks); by = peak -> (peak.chrom, peak.left, peak.right, peak.id))
    chrom_index = Dict{String,Vector{Int}}()
    for (index, peak) in enumerate(sorted)
        push!(get!(chrom_index, peak.chrom, Int[]), index)
    end
    metadata_copy = Dict{String,Any}(string(key) => value for (key, value) in metadata)
    metadata_provenance(metadata_copy) === nothing && stamp_provenance!(
        metadata_copy;
        label="PeakSet",
        source="PeakSet",
        notes=["constructed from peak vector"],
        parameters=(peak_count=length(sorted), chrom_count=length(chrom_index)))
    return PeakSet(sorted, chrom_index, metadata_copy)
end

function MethylationExperiment(methylated::SparseMatrixCSC{Int,Int}, total::SparseMatrixCSC{Int,Int}, regions::IntervalCollection, sample_metadata::DataFrame, sample_ids::Vector{String}; metadata::AbstractDict=Dict{String,Any}())
    metadata_copy = Dict{String,Any}(string(key) => value for (key, value) in metadata)
    metadata_provenance(metadata_copy) === nothing && stamp_provenance!(
        metadata_copy;
        label="MethylationExperiment",
        source="MethylationExperiment",
        notes=["constructed from binned methylation counts"],
        parameters=(region_count=length(regions.intervals), sample_count=length(sample_ids)))
    return MethylationExperiment(methylated, total, regions, sample_metadata, sample_ids, metadata_copy)
end

function DataFrames.DataFrame(peaks::AbstractVector{<:Peak})
    rows = collect(peaks)
    return DataFrames.DataFrame(
        id = [peak.id for peak in rows],
        chrom = [peak.chrom for peak in rows],
        start = [peak.left for peak in rows],
        stop = [peak.right for peak in rows],
        summit = [peak.summit for peak in rows],
        score = [peak.score for peak in rows],
        pvalue = [peak.pvalue for peak in rows],
        qvalue = [peak.qvalue for peak in rows])
end

function DataFrames.DataFrame(set::PeakSet)
    return DataFrames.DataFrame(set.peaks)
end

function DataFrames.DataFrame(support::AbstractVector{<:PeakSupport})
    rows = collect(support)
    return DataFrames.DataFrame(
        peak_id = [row.peak_id for row in rows],
        sample_id = [row.sample_id for row in rows],
        fragment_count = [row.fragment_count for row in rows],
        total_overlap = [row.total_overlap for row in rows],
        mean_overlap = [row.mean_overlap for row in rows],
        enrichment = [row.enrichment for row in rows])
end

function DataFrames.DataFrame(results::AbstractVector{<:MethylationResult})
    rows = collect(results)
    return DataFrames.DataFrame(
        region_id = [row.region_id for row in rows],
        chrom = [row.chrom for row in rows],
        start = [row.left for row in rows],
        stop = [row.right for row in rows],
        group1_mean = [row.group1_mean for row in rows],
        group2_mean = [row.group2_mean for row in rows],
        delta_methylation = [row.delta_methylation for row in rows],
        stat = [row.stat for row in rows],
        pvalue = [row.pvalue for row in rows],
        padj = [row.padj for row in rows])
end

function DataFrames.DataFrame(results::AbstractVector{<:TadResult})
    rows = collect(results)
    return DataFrames.DataFrame(
        chrom = [row.chrom for row in rows],
        left_bin = [row.left_bin for row in rows],
        right_bin = [row.right_bin for row in rows],
        score = [row.score for row in rows],
        state = [row.state for row in rows])
end

function coverage_depth(coverage::SparseCoverageVector, position::Integer)
    position < 1 && return 0
    idx = searchsortedlast(coverage.positions, Int(position))
    idx == 0 && return 0

    return coverage.depths[idx]
end

function coverage_segments(coverage::SparseCoverageVector; chrom::String="unknown")
    segments = CoverageSegment[]
    isempty(coverage.positions) && return segments
    for index in 1:length(coverage.positions)-1
        depth = coverage.depths[index]
        depth <= 0 && continue
        start_pos = coverage.positions[index]
        stop_pos = coverage.positions[index + 1] - 1
        start_pos <= stop_pos || continue
        push!(segments, CoverageSegment(String(chrom), start_pos, stop_pos, depth))
    end

    return segments
end

function _trimmed_mean(values::AbstractVector{<:Real}; trim::Real=0.2)
    isempty(values) && return 0.0
    trim < 0 && throw(ArgumentError("trim must be nonnegative"))
    trim <= 0 && return mean(Float64.(values))
    sorted = sort!(Float64.(values))
    count_values = length(sorted)
    drop = min(floor(Int, count_values * trim), count_values ÷ 2)
    keep = sorted[drop + 1:count_values - drop]
    isempty(keep) && return mean(sorted)
    return mean(keep)
end

function _median_absolute_deviation(values::AbstractVector{<:Real})
    isempty(values) && return 0.0
    center = median(Float64.(values))
    return median(abs.(Float64.(values) .- center))
end

function _tricube(weight::Real)
    abs(weight) >= 1 && return 0.0
    value = 1 - abs(weight)^3
    return value^3
end

function _weighted_local_linear(x_values::Vector{Float64}, y_values::Vector{Float64}, target::Real, weights::Vector{Float64})
    design = hcat(ones(length(x_values)), x_values .- Float64(target))
    weighted_design = design .* sqrt.(weights)
    weighted_response = y_values .* sqrt.(weights)
    coefficients = weighted_design 
        weighted_response
    return coefficients[1]
end

function _robust_lowess(x_values::AbstractVector{<:Real}, y_values::AbstractVector{<:Real}; span::Real=0.5, iterations::Int=2)
    n = length(x_values)
    n == length(y_values) || throw(ArgumentError("x_values and y_values must have the same length"))
    n == 0 && return Float64[]
    n == 1 && return Float64[y_values[1]]

    x = Float64.(x_values)
    y = Float64.(y_values)
    span < 0 && throw(ArgumentError("span must be nonnegative"))
    k = clamp(ceil(Int, span * n), 2, n)
    fitted = similar(y)
    robust = ones(Float64, n)

    for _ in 1:max(iterations, 1)
        for index in 1:n
            distances = abs.(x .- x[index])
            bandwidth = sort(distances)[k]
            if bandwidth <= eps(Float64)
                fitted[index] = y[index]
                continue
            end
            local_weights = [robust[pos] * _tricube(distances[pos] / bandwidth) for pos in 1:n]
            if sum(local_weights) <= 0
                fitted[index] = y[index]
                continue
            end
            fitted[index] = _weighted_local_linear(x, y, x[index], local_weights)
        end

        residuals = y .- fitted
        scale = 6.0 * max(_median_absolute_deviation(residuals), eps(Float64))
        if scale <= eps(Float64)
            break
        end
        robust = Float64[]
        sizehint!(robust, n)
        for residual in residuals
            u = abs(residual) / scale
            push!(robust, u < 1 ? (1 - u^2)^2 : 0.0)
        end
    end

    return fitted
end

function _gc_binned_fit(gc_values::AbstractVector{<:Real}, observed::AbstractVector{<:Real})
    n = length(gc_values)
    n == length(observed) || throw(ArgumentError("gc_values and observed must have the same length"))
    n == 0 && return Float64[]
    n == 1 && return Float64[observed[1]]

    order = sortperm(Float64.(gc_values))
    sorted_gc = Float64.(gc_values)[order]
    sorted_observed = Float64.(observed)[order]
    bin_count = clamp(round(Int, sqrt(n)), 2, min(10, n))
    bin_edges = round.(Int, range(1, n + 1, length=bin_count + 1))
    anchor_gc = Float64[]
    anchor_observed = Float64[]

    for bin_index in 1:bin_count
        left = bin_edges[bin_index]
        right = max(left, bin_edges[bin_index + 1] - 1)
        indices = left:right
        isempty(indices) && continue
        push!(anchor_gc, median(sorted_gc[indices]))
        push!(anchor_observed, median(sorted_observed[indices]))
    end

    isempty(anchor_gc) && return fill(median(Float64.(observed)), n)
    if length(anchor_gc) == 1
        return fill(anchor_observed[1], n)
    end

    ordering = sortperm(anchor_gc)
    anchor_gc = anchor_gc[ordering]
    anchor_observed = anchor_observed[ordering]

    fitted = similar(Float64.(gc_values))
    for (index, gc) in enumerate(Float64.(gc_values))
        if gc <= first(anchor_gc)
            fitted[index] = first(anchor_observed)
            continue
        elseif gc >= last(anchor_gc)
            fitted[index] = last(anchor_observed)
            continue
        end
        upper = searchsortedfirst(anchor_gc, gc)
        lower = upper - 1
        span = anchor_gc[upper] - anchor_gc[lower]
        weight = span <= eps(Float64) ? 0.0 : (gc - anchor_gc[lower]) / span
        fitted[index] = (1 - weight) * anchor_observed[lower] + weight * anchor_observed[upper]
    end
    return fitted
end

function _segment_length(segment::CoverageSegment)
    return segment.stop - segment.start + 1
end

function _finalize_peak(chrom::String, segments::Vector{CoverageSegment}, lambda::Real, peaks::Vector{Peak})
    isempty(segments) && return
    start_pos = first(segments).start
    stop_pos = last(segments).stop
    summit_segment = argmax(getfield.(segments, :depth))
    summit = (segments[summit_segment].start + segments[summit_segment].stop) ÷ 2
    depth_sum = sum(segment.depth * _segment_length(segment) for segment in segments)
    total_bases = sum(_segment_length(segment) for segment in segments)
    average_depth = total_bases > 0 ? depth_sum / total_bases : 0.0
    best_depth = maximum(segment.depth for segment in segments)
    score = max(average_depth, Float64(best_depth))
    pvalue = _poisson_pvalue(best_depth, lambda)
    push!(peaks, Peak("$(chrom):$(start_pos)-$(stop_pos)", chrom, start_pos, stop_pos, summit, score, pvalue, pvalue))
end

function _coverage_span(fragments::AbstractVector{<:GenomicInterval}, chrom::String)
    maxima = [fragment.right for fragment in fragments if fragment.chrom == chrom]
    isempty(maxima) && return 0
    return maximum(maxima)
end

function calculate_coverage(fragments::AbstractVector{<:GenomicInterval}; chrom_lengths::AbstractDict=Dict{String,Int}())
    events = Dict{String,Vector{Tuple{Int,Int}}}()
    for fragment in fragments
        push!(get!(events, fragment.chrom, Tuple{Int,Int}[]), (fragment.left, 1))
        push!(events[fragment.chrom], (fragment.right + 1, -1))
    end

    coverage = Dict{String,SparseCoverageVector}()
    for (chrom, chrom_events) in events
        sort!(chrom_events; by = event -> (event[1], -event[2]))
        positions = Int[]
        depths = Int[]
        current_depth = 0
        current_position = first(chrom_events)[1]
        if current_position > 1
            push!(positions, 1)
            push!(depths, 0)
        end
        index = 1
        while index <= length(chrom_events)
            position = chrom_events[index][1]
            delta = 0
            while index <= length(chrom_events) && chrom_events[index][1] == position
                delta += chrom_events[index][2]
                index += 1
            end
            current_depth += delta
            push!(positions, position)
            push!(depths, current_depth)
            current_position = position
        end
        span = get(chrom_lengths, chrom, max(current_position, positions[end]))
        if isempty(depths) || depths[end] != 0
            push!(positions, span + 1)
            push!(depths, 0)
        elseif positions[end] <= span
            push!(positions, span + 1)
            push!(depths, 0)
        end
        coverage[chrom] = SparseCoverageVector(positions, depths, span)
    end

    for (chrom, span) in chrom_lengths
        haskey(coverage, chrom) && continue
        coverage[String(chrom)] = SparseCoverageVector([1, Int(span) + 1], [0, 0], Int(span))
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "calculate_coverage";
        parameters=(fragment_count=length(fragments), chrom_count=length(coverage)))
    end
    return coverage
end

function _coverage_mean(coverage::SparseCoverageVector, start_pos::Int, stop_pos::Int)
    start_pos > stop_pos && return 0.0
    total = 0.0
    for index in 1:length(coverage.positions)-1
        depth = coverage.depths[index]
        depth <= 0 && continue
        segment_start = coverage.positions[index]
        segment_stop = coverage.positions[index + 1] - 1
        overlap_start = max(start_pos, segment_start)
        overlap_stop = min(stop_pos, segment_stop)
        overlap_start <= overlap_stop || continue
        total += depth * (overlap_stop - overlap_start + 1)
    end
    return total / (stop_pos - start_pos + 1)
end

function _poisson_pvalue(depth::Real, background::Real)
    background <= 0 && return 0.0
    return ccdf(Poisson(max(background, eps(Float64))), max(Int(round(depth)) - 1, 0))
end

function call_peaks(coverage::Dict{String,SparseCoverageVector}; pvalue_threshold::Real=0.01, min_depth::Int=1, merge_gap::Int=0)
    peaks = Peak[]
    for (chrom, chrom_coverage) in sort!(collect(coverage); by = first)
        segments = coverage_segments(chrom_coverage; chrom=chrom)
        isempty(segments) && continue
        positive_segments = [segment for segment in segments if segment.depth > 0]
        isempty(positive_segments) && continue
        background = max(_trimmed_mean([segment.depth for segment in positive_segments]; trim=0.25), 1.0)
        lambda_floor = max(background, Float64(min_depth), eps(Float64))

        current_segments = CoverageSegment[]
        current_background = lambda_floor
        best_segment_depth = 0
        best_segment_pvalue = 1.0

        for (index, segment) in enumerate(positive_segments)
            local_left = max(1, index - 3)
            local_right = min(length(positive_segments), index + 3)
            local_neighbours = Float64[positive_segments[pos].depth for pos in local_left:local_right if pos != index]
            local_background = isempty(local_neighbours) ? lambda_floor : max(_trimmed_mean(local_neighbours; trim=0.2), lambda_floor)
            lambda = max(lambda_floor, local_background)
            pvalue = _poisson_pvalue(segment.depth, lambda)
            keep = segment.depth >= min_depth && pvalue <= pvalue_threshold

            if keep
                if isempty(current_segments)
                    current_segments = [segment]
                    current_background = lambda
                    best_segment_depth = segment.depth
                    best_segment_pvalue = pvalue
                else
                    gap = segment.start - last(current_segments).stop - 1
                    if gap <= merge_gap
                        push!(current_segments, segment)
                        current_background = max(current_background, lambda)
                        if segment.depth > best_segment_depth
                            best_segment_depth = segment.depth
                            best_segment_pvalue = pvalue
                        end
                    else
                        _finalize_peak(chrom, current_segments, current_background, peaks)
                        current_segments = [segment]
                        current_background = lambda
                        best_segment_depth = segment.depth
                        best_segment_pvalue = pvalue
                    end
                end
            elseif !isempty(current_segments)
                gap = segment.start - last(current_segments).stop - 1
                gap <= merge_gap || begin
                    _finalize_peak(chrom, current_segments, current_background, peaks)
                    current_segments = CoverageSegment[]
                    current_background = lambda_floor
                    best_segment_depth = 0
                    best_segment_pvalue = 1.0
                end
            end
        end

        isempty(current_segments) || _finalize_peak(chrom, current_segments, current_background, peaks)
    end

    qvalues = isempty(peaks) ? Float64[] : benjamini_hochberg([peak.pvalue for peak in peaks])
    adjusted = Peak[Peak(peak.id, peak.chrom, peak.left, peak.right, peak.summit, peak.score, peak.pvalue, qvalues[index]) for (index, peak) in enumerate(peaks)]
    result = PeakSet(
        adjusted;
        metadata=Dict{String,Any}(
            "peak_count" => length(adjusted),
            "chrom_count" => length(coverage),
            "provenance" => provenance_record(
                "PeakSet",
                "Epigenetics/call_peaks";
                notes=["Poisson enrichment over positive coverage segments"],
                parameters=(pvalue_threshold=Float64(pvalue_threshold), min_depth=Int(min_depth), merge_gap=Int(merge_gap), chrom_count=length(coverage), peak_count=length(adjusted)))))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "call_peaks";
        parameters=(pvalue_threshold=Float64(pvalue_threshold), min_depth=Int(min_depth),
        merge_gap=Int(merge_gap), chrom_count=length(coverage), peak_count=length(adjusted)))
    end
    return result
end

function _window_gc(sequence::BioSequence{DNAAlphabet}, start_pos::Int, stop_pos::Int)
    start_pos > stop_pos && return 0.0
    left = max(start_pos, 1)
    right = min(stop_pos, length(sequence))
    left > right && return 0.0
    gc = 0
    @inbounds for idx in left:right
        byte = sequence.data[idx]
        if byte == UInt8('G') || byte == UInt8('g') || byte == UInt8('C') || byte == UInt8('c')
            gc += 1
        end
    end
    return gc / (right - left + 1)
end

function normalize_gc_bias(sequence::BioSequence{DNAAlphabet}, coverage::SparseCoverageVector; window_size::Int=50)
    sequence_length = min(length(sequence), coverage.span)
    windows = GenomicInterval[]
    gc_values = Float64[]
    observed = Float64[]
    index = 1
    while index <= sequence_length
        stop_pos = min(index + window_size - 1, sequence_length)
        push!(windows, GenomicInterval("sequence", index, stop_pos, '.'))
        push!(gc_values, _window_gc(sequence, index, stop_pos))
        push!(observed, _coverage_mean(coverage, index, stop_pos))
        index = stop_pos + 1
    end

    if isempty(gc_values)
        return (windows=windows, gc=gc_values, observed=observed, fitted=Float64[], corrected_coverage=coverage)
    end

    fitted = length(gc_values) <= 8 || length(unique(gc_values)) <= 3 ? _gc_binned_fit(gc_values, observed) : _robust_lowess(gc_values, observed; span=min(0.6, max(0.35, 3 / length(gc_values))), iterations=2)
    if any(isnan, fitted) || any(fitted .<= eps(Float64))
        fitted = _gc_binned_fit(gc_values, observed)
    end
    baseline = max(median(observed), eps(Float64))
    corrected = baseline .* (observed ./ max.(fitted, eps(Float64)))
    corrected = max.(corrected, 0.0)
    corrected_positions = Int[windows[1].left]
    corrected_depths = Int[round(Int, corrected[1])]
    for (window, depth) in zip(windows[2:end], corrected[2:end])
        push!(corrected_positions, window.left)
        push!(corrected_depths, round(Int, depth))
    end
    push!(corrected_positions, sequence_length + 1)
    push!(corrected_depths, 0)

    return (
        windows=windows,
        gc=gc_values,
        observed=observed,
        fitted=fitted,
        corrected_coverage=SparseCoverageVector(corrected_positions, corrected_depths, coverage.span))
end

normalize_gc_bias(sequence::AbstractString, coverage::SparseCoverageVector; window_size::Int=50) =
    normalize_gc_bias(BioSequence{DNAAlphabet}(String(sequence); validate=false), coverage; window_size=window_size)

function _peak_lookup(peaks::PeakSet)
    lookup = Dict{String,Vector{Int}}()
    for (index, peak) in enumerate(peaks.peaks)
        push!(get!(lookup, peak.chrom, Int[]), index)
    end
    return lookup
end

function _midpoint(fragment::GenomicInterval)
    return (fragment.left + fragment.right) ÷ 2
end

function _overlap_length(left1::Int, right1::Int, left2::Int, right2::Int)
    overlap_left = max(left1, left2)
    overlap_right = min(right1, right2)
    overlap_left <= overlap_right || return 0
    return overlap_right - overlap_left + 1
end

function _best_peak_index(fragment::GenomicInterval, peaks::PeakSet, chrom_peaks::Vector{Int})
    best_index = 0
    best_overlap = 0
    best_distance = typemax(Int)
    fragment_midpoint = _midpoint(fragment)
    for peak_index in chrom_peaks
        peak = peaks.peaks[peak_index]
        overlap = _overlap_length(fragment.left, fragment.right, peak.left, peak.right)
        overlap <= 0 && continue
        peak_midpoint = (peak.left + peak.right) ÷ 2
        distance = abs(fragment_midpoint - peak_midpoint)
        if overlap > best_overlap || (overlap == best_overlap && distance < best_distance)
            best_index = peak_index
            best_overlap = overlap
            best_distance = distance
        end
    end
    return best_index, best_overlap
end

function summarize_peak_support(fragments_by_sample::AbstractDict, peaks::PeakSet)
    sample_ids = sort!(collect(String.(keys(fragments_by_sample))))
    lookup = _peak_lookup(peaks)
    support = PeakSupport[]

    for peak in peaks.peaks
        peak_length = max(peak.right - peak.left + 1, 1)
        background_density = peak.score / peak_length
        for sample_id in sample_ids
            fragments = fragments_by_sample[sample_id]
            fragment_count = 0
            total_overlap = 0
            for fragment in fragments
                fragment.chrom != peak.chrom && continue
                overlap = _overlap_length(fragment.left, fragment.right, peak.left, peak.right)
                overlap <= 0 && continue
                fragment_count += 1
                total_overlap += overlap
            end
            mean_overlap = fragment_count > 0 ? total_overlap / fragment_count : 0.0
            enrichment = background_density > 0 ? (total_overlap / peak_length) / background_density : float(total_overlap)
            push!(support, PeakSupport(peak.id, sample_id, fragment_count, total_overlap, mean_overlap, enrichment))
        end
    end

    return support
end

function count_overlaps(fragments_by_sample::AbstractDict, peaks::PeakSet; multi_thread::Bool=true)
    sample_ids = sort!(collect(String.(keys(fragments_by_sample))))
    peak_ids = [peak.id for peak in peaks.peaks]
    counts = zeros(Int, length(peak_ids), length(sample_ids))
    lookup = _peak_lookup(peaks)

    if multi_thread && length(sample_ids) > 1 && Base.Threads.nthreads() > 1
        Base.Threads.@threads for sample_index in eachindex(sample_ids)
            sample_id = sample_ids[sample_index]
            fragments = fragments_by_sample[sample_id]
            column = zeros(Int, length(peak_ids))
            for fragment in fragments
                chrom_peaks = get(lookup, fragment.chrom, Int[])
                isempty(chrom_peaks) && continue
                best_index, best_overlap = _best_peak_index(fragment, peaks, chrom_peaks)
                best_index == 0 && continue
                if best_overlap > 0
                    column[best_index] += 1
                end
            end
            counts[:, sample_index] = column
        end
    else
        for sample_index in eachindex(sample_ids)
            sample_id = sample_ids[sample_index]
            fragments = fragments_by_sample[sample_id]
            column = zeros(Int, length(peak_ids))
            for fragment in fragments
                chrom_peaks = get(lookup, fragment.chrom, Int[])
                isempty(chrom_peaks) && continue
                best_index, best_overlap = _best_peak_index(fragment, peaks, chrom_peaks)
                best_index == 0 && continue
                if best_overlap > 0
                    column[best_index] += 1
                end
            end
            counts[:, sample_index] = column
        end
    end

    return CountMatrix(sparse(counts), peak_ids, sample_ids)
end

function differential_binding(fragments_by_sample::AbstractDict, peaks::PeakSet, design::AbstractVector{<:Symbol}; min_total::Int=10, shrink::Bool=false, multi_thread::Bool=true, sfType::Symbol=:poscounts)
    count_matrix = count_overlaps(fragments_by_sample, peaks; multi_thread=multi_thread)
    length(design) == length(count_matrix.sample_ids) || throw(ArgumentError("design must match samples"))
    result = differential_expression(count_matrix, design; min_total=min_total, shrink=shrink, sfType=sfType)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "differential_binding";
        parameters=(n_samples=length(count_matrix.sample_ids), n_peaks=length(peaks),
        min_total=min_total, shrink=shrink))
    end
    return result
end

function _bin_key(chrom::String, position::Int, bin_width::Int)
    bin_start = ((position - 1) ÷ bin_width) * bin_width + 1
    return (chrom, bin_start)
end

function _estimate_precision(methylated::AbstractVector{<:Integer}, total::AbstractVector{<:Integer}, mu::Real)
    informative = total .> 0
    count(informative) <= 1 && return 10.0
    proportions = Float64[methylated[index] / total[index] for index in eachindex(total) if informative[index]]
    weights = Float64[total[index] for index in eachindex(total) if informative[index]]
    weighted_mean = sum(weights .* proportions) / sum(weights)
    centered = proportions .- weighted_mean
    weighted_var = sum(weights .* centered .^ 2) / max(sum(weights) - 1, 1)
    binomial_var = max(weighted_mean * (1 - weighted_mean) / max(mean(weights), 1.0), eps(Float64))
    if weighted_var <= binomial_var + eps(Float64)
        return 1000.0
    end
    precision = (weighted_mean * (1 - weighted_mean) - binomial_var) / max(weighted_var - binomial_var, eps(Float64))
    return clamp(precision, 1.0, 500.0)
end

function bin_methylation(calls::AbstractVector{<:MethylationCall}; bin_width::Int=500, sample_metadata::DataFrame=DataFrame())
    sample_ids = sort!(unique(String.(call.sample for call in calls)))
    sample_index = Dict(sample_id => index for (index, sample_id) in enumerate(sample_ids))
    bins = Dict{Tuple{String,Int},Int}()
    bin_intervals = GenomicInterval[]
    methylated = Int[]
    total = Int[]

    function ensure_bin(chrom::String, bin_start::Int)
        key = (chrom, bin_start)
        existing = get(bins, key, 0)
        existing != 0 && return existing
        push!(bin_intervals, GenomicInterval(chrom, bin_start, bin_start + bin_width - 1, '.'))
        push!(methylated, 0)
        push!(total, 0)
        bins[key] = length(bin_intervals)
        return bins[key]
    end

    methylated_matrix = Dict{Tuple{Int,Int},Int}()
    total_matrix = Dict{Tuple{Int,Int},Int}()
    chrom_spans = Dict{String,Int}()

    for call in calls
        bin_index = ensure_bin(call.chrom, _bin_key(call.chrom, call.position, bin_width)[2])
        sample_idx = sample_index[String(call.sample)]
        key = (bin_index, sample_idx)
        total_matrix[key] = get(total_matrix, key, 0) + 1
        call.methylated && (methylated_matrix[key] = get(methylated_matrix, key, 0) + 1)
        chrom_spans[call.chrom] = max(get(chrom_spans, call.chrom, 0), call.position)
    end

    m = zeros(Int, length(bin_intervals), length(sample_ids))
    t = zeros(Int, length(bin_intervals), length(sample_ids))
    for ((row, col), value) in methylated_matrix
        m[row, col] = value
    end
    for ((row, col), value) in total_matrix
        t[row, col] = value
    end

    result = MethylationExperiment(
        sparse(m),
        sparse(t),
        build_collection(bin_intervals),
        sample_metadata,
        sample_ids;
        metadata=Dict{String,Any}(
            "bin_width" => Int(bin_width),
            "provenance" => provenance_record(
                "MethylationExperiment",
                "Epigenetics/bin_methylation";
                notes=["per-call aggregation into fixed genomic bins"],
                parameters=(bin_width=Int(bin_width), region_count=length(bin_intervals), sample_count=length(sample_ids), call_count=length(calls)))))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "bin_methylation";
        parameters=(bin_width=Int(bin_width), call_count=length(calls),
        sample_count=length(sample_ids), region_count=length(bin_intervals)))
    end
    return result
end

function _beta_binomial_loglik(methylated::AbstractVector{<:Integer}, total::AbstractVector{<:Integer}, mu::Real, precision::Real)
    μ = clamp(Float64(mu), eps(Float64), 1 - eps(Float64))
    ϕ = max(Float64(precision), eps(Float64))
    alpha = μ * ϕ
    beta = (1 - μ) * ϕ
    loglik = 0.0
    for (m, n) in zip(methylated, total)
        n == 0 && continue
        loglik += loggamma(n + 1) - loggamma(m + 1) - loggamma(n - m + 1)
        loglik += loggamma(m + alpha) + loggamma(n - m + beta) - loggamma(n + alpha + beta)
        loglik += loggamma(alpha + beta) - loggamma(alpha) - loggamma(beta)
    end
    return loglik
end

function _fit_beta_binomial(methylated::AbstractVector{<:Integer}, total::AbstractVector{<:Integer}; precision_grid=[1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0])
    informative = total .> 0
    any(informative) || return (mu=0.5, precision=1.0, loglik=0.0)
    μ = clamp(sum(methylated[informative]) / sum(total[informative]), eps(Float64), 1 - eps(Float64))
    estimated_precision = _estimate_precision(methylated[informative], total[informative], μ)
    candidate_grid = sort(unique(vcat(precision_grid, [estimated_precision / 4, estimated_precision / 2, estimated_precision, estimated_precision * 2, estimated_precision * 4])))
    best_precision = first(candidate_grid)
    best_loglik = -Inf
    for precision in candidate_grid
        loglik = _beta_binomial_loglik(methylated[informative], total[informative], μ, precision)
        if loglik > best_loglik
            best_loglik = loglik
            best_precision = precision
        end
    end
    return (mu=μ, precision=best_precision, loglik=best_loglik)
end

function differential_methylation(experiment::MethylationExperiment, design::AbstractVector{<:Symbol}; min_total::Int=10)
    groups = unique(design)
    length(groups) == 2 || throw(ArgumentError("design must contain exactly two conditions"))
    length(design) == length(experiment.sample_ids) || throw(ArgumentError("design must match sample count"))

    methylated = Matrix{Int}(experiment.methylated)
    total = Matrix{Int}(experiment.total)
    results = Vector{MethylationResult}(undef, size(total, 1))

    for row in 1:size(total, 1)
        row_total = total[row, :]
        row_methylated = methylated[row, :]
        if sum(row_total) < min_total
            region = experiment.regions.intervals[row]
            results[row] = MethylationResult(string(row), region.chrom, region.left, region.right, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0)
            continue
        end
        group1 = design .== groups[1]
        group2 = design .== groups[2]
        fit_null = _fit_beta_binomial(row_methylated, row_total)
        fit_group1 = _fit_beta_binomial(row_methylated[group1], row_total[group1])
        fit_group2 = _fit_beta_binomial(row_methylated[group2], row_total[group2])
        stat = max(2 * ((fit_group1.loglik + fit_group2.loglik) - fit_null.loglik), 0.0)
        pvalue = ccdf(Chisq(1), stat)
        group1_mean = (sum(row_methylated[group1]) + fit_null.mu * fit_null.precision) / (sum(row_total[group1]) + fit_null.precision)
        group2_mean = (sum(row_methylated[group2]) + fit_null.mu * fit_null.precision) / (sum(row_total[group2]) + fit_null.precision)
        region = experiment.regions.intervals[row]
        results[row] = MethylationResult(string(row), region.chrom, region.left, region.right, group1_mean, group2_mean, group2_mean - group1_mean, stat, pvalue, 1.0)
    end

    qvalues = benjamini_hochberg([result.pvalue for result in results])
    final = [MethylationResult(result.region_id, result.chrom, result.left, result.right, result.group1_mean, result.group2_mean, result.delta_methylation, result.stat, result.pvalue, qvalues[index]) for (index, result) in enumerate(results)]
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "differential_methylation";
        parameters=(n_regions=size(total,1), min_total=min_total, n_sig=count(r -> r.padj < 0.05, final)))
    end
    return final
end

function _binary_accessibility(matrix::SparseMatrixCSC{Int,Int})
    rows, columns, values = findnz(matrix)
    return sparse(rows, columns, ones(Float64, length(values)), size(matrix, 1), size(matrix, 2))
end

function tfidf(chromatin::SingleCellChromatinExperiment)
    binary = _binary_accessibility(chromatin.open_peaks)
    rows, columns, values = findnz(binary)
    nrows, ncols = size(binary)
    row_counts = zeros(Float64, nrows)
    column_counts = zeros(Float64, ncols)

    for position in eachindex(values)
        row_counts[rows[position]] += 1
        column_counts[columns[position]] += 1
    end

    idf = log.(1 .+ ncols ./ (1 .+ row_counts))
    transformed = Vector{Float64}(undef, length(values))
    for position in eachindex(values)
        column = columns[position]
        row = rows[position]
        tf = 1 / max(column_counts[column], eps(Float64))
        transformed[position] = log1p(tf * idf[row])
    end
    result = sparse(rows, columns, transformed, nrows, ncols)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "tfidf";
        parameters=(n_peaks=nrows, n_cells=ncols, nnz=length(values)))
    end
    return result
end

function rsvd(matrix::AbstractMatrix{<:Real}; rank::Int=50, oversample::Int=10, n_iter::Int=2)
    a = Matrix{Float64}(matrix)
    rows, cols = size(a)
    target = min(rank + oversample, min(rows, cols))
    random_matrix = randn(cols, target)
    y = a * random_matrix
    for _ in 1:n_iter
        y = a * (a' * y)
    end
    q = qr(y).Q[:, 1:target]
    b = q' * a
    svd_small = svd(b; full=false)
    components = min(rank, size(svd_small.V, 2))

    return q * svd_small.U[:, 1:components], svd_small.S[1:components], svd_small.V[:, 1:components]
end

function run_lsi(chromatin::SingleCellChromatinExperiment; n_components::Int=30, randomized::Bool=false)
    matrix = Matrix{Float64}(tfidf(chromatin))
    if randomized
        _, singular_values, right_vectors = rsvd(matrix; rank=n_components)
        embedding = right_vectors .* reshape(singular_values, 1, :) 
    else
        svd_result = svd(matrix; full=false)
        components = min(n_components, size(svd_result.V, 2))
        embedding = svd_result.V[:, 1:components] * Diagonal(svd_result.S[1:components])
    end
    chromatin.reductions["lsi"] = Matrix{Float64}(embedding)

    return chromatin.reductions["lsi"]
end

function _overlap_interval(left1::Int, right1::Int, left2::Int, right2::Int)
    return max(left1, left2) <= min(right1, right2)
end

function gene_activity_score(chromatin::SingleCellChromatinExperiment, genes::AbstractVector{<:GenomicInterval}, peak_intervals::AbstractVector{<:GenomicInterval}; promoter_upstream::Int=2000, promoter_downstream::Int=500)
    peak_matrix = _binary_accessibility(chromatin.open_peaks)
    gene_matrix = zeros(Float64, length(genes), size(peak_matrix, 2))
    promoters_only = promoters(genes, promoter_upstream, promoter_downstream)
    for (gene_index, gene) in enumerate(genes)
        for (peak_index, peak) in enumerate(peak_intervals)
            peak.chrom != gene.chrom && continue
            (_overlap_interval(gene.left, gene.right, peak.left, peak.right) || _overlap_interval(promoters_only[gene_index].left, promoters_only[gene_index].right, peak.left, peak.right)) || continue
            gene_matrix[gene_index, :] .+= peak_matrix[peak_index, :]
        end
    end

    return gene_matrix
end

function calculate_coaccessibility(chromatin::SingleCellChromatinExperiment; min_correlation::Real=0.3, max_pairs::Int=5000)
    matrix = Matrix{Float64}(chromatin.open_peaks .> 0)
    edges = CoaccessibilityEdge[]
    for left in 1:size(matrix, 1)-1
        for right in left+1:size(matrix, 1)
            correlation = cor(matrix[left, :], matrix[right, :])
            isnan(correlation) && continue
            correlation < min_correlation && continue
            push!(edges, CoaccessibilityEdge(chromatin.peak_ids[left], chromatin.peak_ids[right], correlation))
            length(edges) >= max_pairs && return edges
        end
    end

    return edges
end

function compute_motif_deviations(chromatin::SingleCellChromatinExperiment, motif_peak_membership::AbstractDict)
    matrix = Matrix{Float64}(chromatin.open_peaks .> 0)
    background = vec(mean(matrix, dims=1))
    deviations = Dict{String,Vector{Float64}}()
    for (motif_name, peak_indices) in motif_peak_membership
        isempty(peak_indices) && continue
        observed = vec(mean(matrix[collect(peak_indices), :], dims=1))
        baseline = mean(background)
        spread = std(background)
        deviations[String(motif_name)] = spread > 0 ? (observed .- baseline) ./ spread : observed .- baseline
    end

    return deviations
end

struct FootprintResult <: AbstractAnalysisResult
    motif_id::String
    position::Int
    centrality::Float64
    flank_mean::Float64
    center_mean::Float64
    provenance::ResultProvenance
end

FootprintResult(motif_id, position, centrality, flank_mean, center_mean) =
    FootprintResult(motif_id, position, centrality, flank_mean, center_mean, provenance_record("FootprintResult", "epigenetics"))

function detect_footprints(signal::AbstractVector{<:Real}, motif_positions::AbstractVector{<:Integer}; flank::Int=10, center::Int=2)
    results = FootprintResult[]
    n = length(signal)
    for position in motif_positions
        left = max(1, Int(position) - flank)
        right = min(n, Int(position) + flank)
        center_left = max(left, Int(position) - center)
        center_right = min(right, Int(position) + center)
        flank_values = Float64[signal[index] for index in left:right if index < center_left || index > center_right]
        center_values = Float64[signal[index] for index in center_left:center_right]
        isempty(center_values) && continue
        flank_mean = isempty(flank_values) ? 0.0 : mean(flank_values)
        center_mean = mean(center_values)
        push!(results, FootprintResult("motif_$(position)", Int(position), flank_mean - center_mean, flank_mean, center_mean))
    end

    return results
end

function _weighted_contact_sum(matrix::AbstractMatrix{<:Real}, index::Int, start_index::Int, stop_index::Int; decay::Real=1.0)
    start_index > stop_index && return 0.0
    total = 0.0
    for partner in start_index:stop_index
        distance = abs(partner - index)
        weight = 1.0 / (1 + distance)^decay
        total += Float64(matrix[partner, index]) * weight
    end
    return total
end

function _moving_average(values::AbstractVector{<:Real}, window::Int)
    window <= 1 && return Float64.(values)
    radius = max(window ÷ 2, 1)
    smoothed = zeros(Float64, length(values))
    for index in eachindex(values)
        left = max(firstindex(values), index - radius)
        right = min(lastindex(values), index + radius)
        smoothed[index] = mean(Float64[values[pos] for pos in left:right])
    end
    return smoothed
end

function directionality_index(contact::ContactMatrix; window::Int=5)
    matrix = Matrix{Float64}(contact.matrix)
    n = size(matrix, 1)
    di = zeros(Float64, n)
    for index in 1:n
        upstream_start = max(1, index - window)
        upstream_stop = index - 1
        downstream_start = index + 1
        downstream_stop = min(n, index + window)
        upstream = upstream_stop >= upstream_start ? _weighted_contact_sum(matrix, index, upstream_start, upstream_stop; decay=1.0) : 0.0
        downstream = downstream_stop >= downstream_start ? _weighted_contact_sum(matrix, index, downstream_start, downstream_stop; decay=1.0) : 0.0
        total = upstream + downstream
        if total > 0
            di[index] = sign(downstream - upstream) * ((downstream - upstream)^2) / (total + 1.0)
        end
    end

    return _moving_average(di, max(3, window))
end

function _merge_tad_candidates(candidates::Vector{TadResult})
    isempty(candidates) && return TadResult[]
    sorted = sort(candidates; by = tad -> (tad.chrom, tad.left_bin, tad.right_bin, tad.score))
    merged = TadResult[]
    current = sorted[1]
    for candidate in sorted[2:end]
        if candidate.chrom == current.chrom && candidate.left_bin <= current.right_bin + 1
            current = TadResult(current.chrom, current.left_bin, max(current.right_bin, candidate.right_bin), max(current.score, candidate.score), current.state)
        else
            push!(merged, current)
            current = candidate
        end
    end
    push!(merged, current)
    return merged
end

function detect_tads(contact::ContactMatrix; window::Int=5, threshold::Real=1.0, min_bins::Int=2)
    di = directionality_index(contact; window=window)
    center = median(di)
    spread = max(_median_absolute_deviation(di), std(di), eps(Float64))
    standardized = (di .- center) ./ spread
    observations = UInt8[]
    for value in standardized
        if value > threshold
            push!(observations, UInt8('R'))
        elseif value < -threshold
            push!(observations, UInt8('L'))
        else
            push!(observations, UInt8('N'))
        end
    end

    hmm = HMM(
        ["left", "domain", "right"],
        [UInt8('L'), UInt8('N'), UInt8('R')],
        [0.05, 0.9, 0.05],
        [0.92 0.06 0.02; 0.03 0.94 0.03; 0.02 0.06 0.92],
        [0.85 0.1 0.05; 0.15 0.7 0.15; 0.05 0.1 0.85])
    path, _ = viterbi(hmm, observations)

    tads = TadResult[]
    index = 1
    while index <= length(path)
        path[index] == 2 || (index += 1; continue)
        start_index = index
        while index <= length(path) && path[index] == 2
            index += 1
        end
        stop_index = index - 1
        stop_index - start_index + 1 < min_bins && continue
        local_score = mean(abs.(di[start_index:stop_index]))
        local_score >= threshold / 4 || continue
        push!(tads, TadResult(contact.chrom, start_index, stop_index, local_score, :domain))
    end

    isempty(tads) || return tads

    fallback = TadResult[]
    if length(standardized) >= min_bins
        for index in 2:length(standardized)-1
            left_sign = sign(standardized[index - 1])
            right_sign = sign(standardized[index + 1])
            if left_sign != 0 && right_sign != 0 && left_sign != right_sign && abs(standardized[index]) <= max(threshold * 2, 1.0)
                left_bin = max(1, index - min_bins + 1)
                right_bin = min(length(standardized), index + min_bins - 1)
                push!(fallback, TadResult(contact.chrom, left_bin, right_bin, mean(abs.(di[left_bin:right_bin])), :domain))
            end
        end
    end
    if isempty(fallback) && !isempty(standardized)
        best = argmin(abs.(standardized))
        left_bin = max(1, best - min_bins + 1)
        right_bin = min(length(standardized), best + min_bins - 1)
        push!(fallback, TadResult(contact.chrom, left_bin, right_bin, mean(abs.(di[left_bin:right_bin])), :domain))
    end

    return _merge_tad_candidates(fallback)
end

function _epi_zscore_columns(x::AbstractMatrix{<:Real})
    data = Matrix{Float64}(x)
    for j in axes(data, 2)
        col = @view data[:, j]
        μ = mean(col)
        σ = std(col)
        if isfinite(σ) && σ > 0
            col .= (col .- μ) ./ σ
        else
            col .= 0.0
        end
    end
    return data
end

function _epi_kmeans_lloyd(x::AbstractMatrix{<:Real}, k::Int; max_iter::Int=50, seed::Int=1)
    n = size(x, 1)
    n >= k || throw(ArgumentError("k must be <= number of rows"))
    rng = MersenneTwister(seed)
    centers = Matrix{Float64}(x[rand(rng, 1:n, k), :])
    labels = ones(Int, n)

    for _ in 1:max_iter
        changed = false
        for i in 1:n
            xi = @view x[i, :]
            best = 1
            bestd = Inf
            for c in 1:k
                d = sum(abs2, xi .- @view(centers[c, :]))
                if d < bestd
                    bestd = d
                    best = c
                end
            end
            if labels[i] != best
                labels[i] = best
                changed = true
            end
        end

        for c in 1:k
            idx = findall(==(c), labels)
            if isempty(idx)
                centers[c, :] .= x[rand(rng, 1:n), :]
            else
                centers[c, :] .= vec(mean(x[idx, :], dims=1))
            end
        end
        !changed && break
    end

    return labels, centers
end

"""
    chromhmm_like_segmentation(mark_matrix; n_states=15)

ChromHMM-style state segmentation over bins using k-means state emission
estimates and empirical transition probabilities.
"""
function chromhmm_like_segmentation(mark_matrix::AbstractMatrix{<:Real}; n_states::Int=15, seed::Int=1)
    X = _epi_zscore_columns(mark_matrix)
    labels, centers = _epi_kmeans_lloyd(X, n_states; seed=seed)

    T = zeros(Float64, n_states, n_states)
    for i in 1:(length(labels) - 1)
        T[labels[i], labels[i + 1]] += 1
    end
    for i in 1:n_states
        row_sum = sum(T[i, :])
        if row_sum > 0
            T[i, :] ./= row_sum
        else
            T[i, :] .= 1 / n_states
        end
    end

    result_tuple = (state=labels, emissions=centers, transition=T)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "chromhmm_like_segmentation";
        parameters=(n_states=n_states, n_bins=size(mark_matrix, 1), n_marks=size(mark_matrix, 2), seed=seed))
    end
    return result_tuple
end

"""
    hic_ice_normalize(contact; max_iter=200, tol=1e-6)

ICE-like iterative matrix balancing for Hi-C contact maps.
"""
function hic_ice_normalize(contact::AbstractMatrix{<:Real}; max_iter::Int=200, tol::Real=1e-6)
    M = max.(Float64.(contact), 0.0)
    n = size(M, 1)
    size(M, 2) == n || throw(DimensionMismatch("contact matrix must be square"))
    bias = ones(Float64, n)

    for _ in 1:max_iter
        row_sums = vec(sum(M, dims=2))
        valid = row_sums .> 0
        scale = ones(Float64, n)
        scale[valid] .= 1.0 ./ row_sums[valid]
        M = Diagonal(scale) * M * Diagonal(scale)
        bias .*= scale
        if any(valid) && maximum(abs.(row_sums[valid] .- mean(row_sums[valid]))) < Float64(tol)
            break
        end
    end

    return (normalized=M, bias=bias)
end

"""
    hic_kr_normalize(contact; max_iter=200, tol=1e-6)

Knight-Ruiz style matrix balancing for Hi-C contact maps.
"""
function hic_kr_normalize(contact::AbstractMatrix{<:Real}; max_iter::Int=200, tol::Real=1e-6)
    M = max.(Float64.(contact), 0.0)
    n = size(M, 1)
    size(M, 2) == n || throw(DimensionMismatch("contact matrix must be square"))

    bias = ones(Float64, n)
    valid = vec(sum(M, dims=2)) .> 0

    for _ in 1:max_iter
        previous = copy(bias)
        Mb = M * bias
        for i in 1:n
            if valid[i] && Mb[i] > 0
                bias[i] /= sqrt(Mb[i])
            end
        end

        s = bias[valid]
        if !isempty(s)
            gmean = exp(mean(log.(max.(s, eps(Float64)))))
            gmean > 0 && (bias ./= gmean)
        end

        if maximum(abs.(bias .- previous)) < Float64(tol)
            break
        end
    end

    normalized = Diagonal(bias) * M * Diagonal(bias)

    return (normalized=normalized, bias=bias)
end

"""
    hic_ab_compartments(contact)

A/B compartment calling from observed/expected transformed Hi-C matrix.
"""
function hic_ab_compartments(contact::AbstractMatrix{<:Real})
    M = max.(Float64.(contact), 0.0)
    n = size(M, 1)
    size(M, 2) == n || throw(DimensionMismatch("contact matrix must be square"))

    expected = zeros(Float64, n, n)
    for d in 0:(n - 1)
        vals = [M[i, i + d] for i in 1:(n - d)]
        μ = isempty(vals) ? 1.0 : max(mean(vals), eps(Float64))
        for i in 1:(n - d)
            expected[i, i + d] = μ
            expected[i + d, i] = μ
        end
    end

    oe = log2.(M .+ 1.0) .- log2.(expected .+ 1.0)
    c = cor(oe)
    c[.!isfinite.(c)] .= 0.0
    ev = eigen(Symmetric(c))
    pc1 = ev.vectors[:, argmax(ev.values)]

    return (pc1=pc1, compartment=ifelse.(pc1 .>= 0, "A", "B"))
end

"""
    tobias_like_footprints(cut_profile, motif_centers)

TOBIAS-style footprint depletion summary around motif centers.
"""
function tobias_like_footprints(cut_profile::AbstractVector{<:Real}, motif_centers::AbstractVector{<:Integer}; flank::Int=20, center_width::Int=10)
    x = Float64.(cut_profile)
    n = length(x)
    flank >= 1 || throw(ArgumentError("flank must be >= 1"))
    center_width >= 1 || throw(ArgumentError("center_width must be >= 1"))

    out = DataFrame(center=Int[], flank_mean=Float64[], center_mean=Float64[], depletion_score=Float64[])
    for c in motif_centers
        left = max(1, c - flank)
        right = min(n, c + flank)
        center_left = max(1, c - div(center_width, 2))
        center_right = min(n, c + div(center_width, 2))

        flank_vals = vcat(x[left:(center_left - 1)], x[(center_right + 1):right])
        center_vals = x[center_left:center_right]
        fμ = isempty(flank_vals) ? 0.0 : mean(flank_vals)
        cμ = isempty(center_vals) ? 0.0 : mean(center_vals)
        push!(out, (Int(c), fμ, cμ, fμ - cμ))
    end
    sort!(out, :depletion_score, rev=true)

    return out
end

"""
    parse_bismark_coverage(lines)

Parse Bismark-style coverage lines and summarize CpG/non-CpG methylation.
"""
function parse_bismark_coverage(lines::AbstractVector{<:AbstractString})
    chrom = String[]
    start = Int[]
    stop = Int[]
    methylation_pct = Float64[]
    meth_count = Int[]
    unmeth_count = Int[]
    context = String[]

    for line in lines
        s = split(strip(String(line)), '\t')
        length(s) < 6 && continue
        push!(chrom, s[1])
        push!(start, parse(Int, s[2]))
        push!(stop, parse(Int, s[3]))
        push!(methylation_pct, parse(Float64, s[4]))
        push!(meth_count, parse(Int, s[5]))
        push!(unmeth_count, parse(Int, s[6]))
        push!(context, iseven(parse(Int, s[2])) ? "CpG" : "CHH")
    end

    df = DataFrame(
        chrom=chrom,
        start=start,
        stop=stop,
        methylation_pct=methylation_pct,
        meth_count=meth_count,
        unmeth_count=unmeth_count,
        context=context)

    summary = combine(groupby(df, :context),
        :meth_count => sum => :meth,
        :unmeth_count => sum => :unmeth)
    summary[!, :beta] = summary.meth ./ max.(summary.meth .+ summary.unmeth, 1)

    return (calls=df, summary=summary)
end

"""
    insulation_score(contact; window=5)

Compute insulation profile from local cross-boundary contact intensity.
"""
function insulation_score(contact::AbstractMatrix{<:Real}; window::Int=5)
    M = max.(Float64.(contact), 0.0)
    n = size(M, 1)
    size(M, 2) == n || throw(DimensionMismatch("contact matrix must be square"))
    window >= 1 || throw(ArgumentError("window must be >= 1"))

    ins = fill(NaN, n)
    for i in 1:n
        left = max(1, i - window + 1):i
        right_start = i + 1
        right_stop = min(n, i + window)
        if right_start > right_stop
            continue
        end
        block = @view M[left, right_start:right_stop]
        μ = mean(block)
        ins[i] = log2(μ + 1.0)
    end

    finite = isfinite.(ins)
    z = zeros(Float64, n)
    if any(finite)
        μ = mean(ins[finite])
        σ = std(ins[finite])
        if isfinite(σ) && σ > 0
            z[finite] .= (ins[finite] .- μ) ./ σ
        end
    end

    return DataFrame(bin=1:n, insulation=ins, zscore=z)
end

"""
    compartment_boundary_candidates(insulation; q=0.2)

Identify candidate boundaries as local minima in insulation tracks.
"""
function compartment_boundary_candidates(insulation; q::Real=0.2)
    0.0 < q < 1.0 || throw(ArgumentError("q must be in (0, 1)"))
    vals = insulation isa DataFrame ? Float64.(insulation.insulation) : Float64.(insulation)
    n = length(vals)
    n >= 3 || return DataFrame(bin=Int[], insulation=Float64[], is_boundary=Bool[])

    finite = vals[isfinite.(vals)]
    isempty(finite) && return DataFrame(bin=Int[], insulation=Float64[], is_boundary=Bool[])
    cutoff = quantile(finite, Float64(q))

    out = DataFrame(bin=Int[], insulation=Float64[], is_boundary=Bool[])
    for i in 2:(n - 1)
        v = vals[i]
        isfinite(v) || continue
        is_min = v <= vals[i - 1] && v <= vals[i + 1] && v <= cutoff
        push!(out, (i, v, is_min))
    end
    sort!(out, :insulation)

    return out
end

end
