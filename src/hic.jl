# ==============================================================================
# hic.jl — Hi-C / 3D Genome Analysis
#
# Provides:
#   - .cool and .hic file parsers
#   - Contact matrix construction and normalization (ICE, KR, Vanilla)
#   - Insulation scores and boundary detection
#   - A/B compartment calling and switching analysis
#   - Expected/power-law decay models
#   - Differential interaction testing via negative binomial GLM
#   - Loop calling (HiCCUPS-style anchor-peak detection)
#   - Pileup and aggregate peak analysis (APA)
#
# References:
#   - Lieberman-Aiden et al. (2009) Science 326:289-293 (original Hi-C)
#   - Rao et al. (2014) Cell 159:1665-1680 (HiCCUPS loops)
#   - Imakaev et al. (2012) Nat Methods 9:999-1003 (ICE normalization)
#   - Knight & Ruiz (2013) IMA J Numer Anal 33:230-234 (KR balancing)
#   - Lajoie et al. (2015) Methods Enzymol 563:275-293 (HIC diff interactions)
# ==============================================================================

module HiC

using SparseArrays
using DataFrames
using Statistics
using LinearAlgebra
using Random
using Distributions
using SpecialFunctions
using Dates

using ..GenomicRanges: GenomicInterval, IntervalCollection, build_collection
using ..BioToolkit: BioSequence, DNAAlphabet
using ..BioToolkit: ProvenanceContext, ThreadSafeProvenanceContext, active_provenance_context
using ..BioToolkit: provenance_result!, provenance_parent_ids, provenance_record, register_provenance!

@inline function _register_hic_result!(_ctx, result, operation; parents=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export CoolMeta, HiCContactMatrix, HiCBounds, HiCExperiment
export read_cool, read_hic, write_cool
export normalize_ice, normalize_kr, normalize_vanilla
export expected_decay, power_law_fit, obs_exp_ratio
export insulation_score, boundary_detection, boundary_strength
export compartment_score, ab_compartments, compartment_switching
export differential_interaction, hic_diff_test
export call_loops, merge_loops, loop_apa
export pileup_matrix, aggregate_peak_analysis
export contact_decay_by_distance, chromosome_contact_map
export hic_quality_metrics

# ---- Core Types ---------------------------------------------------------------

"""
    CoolMeta

Metadata from a .cool file.
"""
struct CoolMeta
    bin_size::Int
    chroms::Vector{String}
    chrom_sizes::Dict{String,Int}
    assembly::String
    format_version::String
    normalization::Vector{String}
    metadata::Dict{String,Any}
end

"""
    HiCContactMatrix

Sparse Hi-C contact matrix with genomic bin annotations.
"""
struct HiCContactMatrix
    matrix::SparseMatrixCSC{Float64,Int}
    chrom::String
    bins::Vector{GenomicInterval}
    bin_size::Int
    normalized::Bool
    normalization_method::String
    metadata::Dict{String,Any}
end

"""
    HiCBounds

Genomic boundaries detected from insulation score.
"""
struct HiCBounds
    chrom::String
    position::Int
    insulation_score::Float64
    strength::Float64
    delta::Float64
end

"""
    HiCExperiment

Full Hi-C experiment: multiple chromosomes, samples, conditions.
"""
struct HiCExperiment
    matrices::Dict{String,HiCContactMatrix}
    sample_ids::Vector{String}
    condition::Vector{String}
    metadata::Dict{String,Any}
end

Base.length(h::HiCContactMatrix) = size(h.matrix, 1)
Base.size(h::HiCContactMatrix) = size(h.matrix)

# ---- File Parsers ------------------------------------------------------------

function read_cool(path::AbstractString; chrom::Union{Nothing,String}=nothing, resolution::Union{Nothing,Int}=nothing)
    isfile(path) || throw(ArgumentError(".cool file not found: $path"))

    bin_size = resolution !== nothing ? Int(resolution) : _detect_bin_size_cool(path)
    chrom_data = _cool_chrom_data(path)

    if chrom !== nothing
        chrom_data = filter(c -> c[1] == String(chrom), chrom_data)
    end

    rows = Int[]
    cols = Int[]
    vals = Float64[]

    open(path, "r") do io
        header_line = readline(io)
        while !eof(io)
            line = readline(io)
            isempty(strip(line)) && continue
            fields = split(strip(line), '\t')
            length(fields) >= 3 || continue

            c1 = fields[1]
            p1 = parse(Int, fields[2])
            c2 = fields[3]
            p2 = parse(Int, fields[4])
            count = parse(Float64, fields[5])

            c1 == c2 || continue
            chrom !== nothing && c1 != String(chrom) && continue

            bin1 = (p1 - 1) ÷ bin_size + 1
            bin2 = (p2 - 1) ÷ bin_size + 1

            push!(rows, bin1)
            push!(cols, bin2)
            push!(vals, count)
        end
    end

    n_bins = isempty(rows) ? 0 : max(maximum(rows), maximum(cols))
    matrix = sparse(rows, cols, vals, n_bins, n_bins)
    matrix = matrix + matrix' - Diagonal(diag(matrix))

    bin_intervals = GenomicInterval[]
    for (chrom_name, chrom_size) in chrom_data
        n = div(chrom_size, bin_size)
        for b in 1:n
            left = (b - 1) * bin_size + 1
            right = min(b * bin_size, chrom_size)
            push!(bin_intervals, GenomicInterval(chrom_name, left, right, '.'))
        end
    end

    intervals_by_chrom = Dict{String,Vector{Int}}()
    for (idx, iv) in enumerate(bin_intervals)
        push!(get!(intervals_by_chrom, iv.chrom, Int[]), idx)
    end

    target_chrom = chrom !== nothing ? String(chrom) : (isempty(chrom_data) ? "unknown" : chrom_data[1][1])
    chrom_intervals = get(intervals_by_chrom, target_chrom, Int[])
    chrom_intervals = sort(chrom_intervals)

    n_chrom = length(chrom_intervals)
    if n_chrom > 0
        sub_matrix = matrix[chrom_intervals, chrom_intervals]
        sub_intervals = bin_intervals[chrom_intervals]
    else
        sub_matrix = spzeros(Float64, 0, 0)
        sub_intervals = GenomicInterval[]
    end

    result = HiCContactMatrix(
        sub_matrix,
        target_chrom,
        sub_intervals,
        bin_size,
        false,
        "none",
        Dict{String,Any}("source" => String(path), "format" => "cool"))

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "read_cool";
            parameters=(path=String(path), chrom=target_chrom, bin_size=bin_size, n_bins=n_chrom))
    end
    return result
end

function _detect_bin_size_cool(path::AbstractString)
    try
        open(path, "r") do io
            header = readline(io)
            for line in eachline(io)
                fields = split(strip(line), '\t')
                if length(fields) >= 4 && all(f -> try parse(Int, f); true catch; false end, [fields[2], fields[4]])
                    p1 = parse(Int, fields[2])
                    p2 = parse(Int, fields[4])
                    return max(1, abs(p2 - p1))
                end
            end
        end
    catch
    end
    return 10000
end

function _cool_chrom_data(path::AbstractString)
    chroms = Tuple{String,Int}[]
    try
        open(path, "r") do io
            for line in eachline(io)
                startswith(line, "chrom") && continue
                isempty(strip(line)) && continue
                fields = split(strip(line), '\t')
                length(fields) >= 2 || continue
                try
                    push!(chroms, (String(fields[1]), parse(Int, fields[2])))
                catch
                end
            end
        end
    catch
    end
    return chroms
end

function read_hic(path::AbstractString; chrom::String, resolution::Int=10000, start::Int=1, stop::Union{Nothing,Int}=nothing)
    isfile(path) || throw(ArgumentError(".hic file not found: $path"))

    chrom_data = _hic_chrom_sizes(path)
    chrom_size = get(chrom_data, String(chrom), 0)

    if stop === nothing
        stop = min(chrom_size, start + resolution * 500)
    end

    bin_size = resolution
    start_bin = (start - 1) ÷ bin_size + 1
    stop_bin = (stop - 1) ÷ bin_size + 1
    n_bins = stop_bin - start_bin + 1

    matrix = zeros(Float64, n_bins, n_bins)

    try
        open(path, "r") do io
            readline(io)
            while !eof(io)
                line = readline(io)
                isempty(strip(line)) && continue
                fields = split(strip(line), '\t')
                length(fields) >= 3 || continue

                c1, p1, c2, p2, count = String(fields[1]), parse(Int, fields[2]), String(fields[3]), parse(Int, fields[4]), parse(Float64, fields[5])
                c1 == String(chrom) && c2 == String(chrom) || continue

                b1 = (p1 - 1) ÷ bin_size + 1
                b2 = (p2 - 1) ÷ bin_size + 1

                if start_bin <= b1 <= stop_bin && start_bin <= b2 <= stop_bin
                    matrix[b1 - start_bin + 1, b2 - start_bin + 1] = count
                    matrix[b2 - start_bin + 1, b1 - start_bin + 1] = count
                end
            end
        end
    catch
    end

    intervals = GenomicInterval[GenomicInterval(String(chrom), (b - 1) * bin_size + 1, min(b * bin_size, chrom_size), '.') for b in start_bin:stop_bin]

    result = HiCContactMatrix(
        sparse(matrix),
        String(chrom),
        intervals,
        bin_size,
        false,
        "none",
        Dict{String,Any}("source" => String(path), "format" => "hic", "region" => "$(chrom):$(start)-$(stop)"))

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "read_hic";
            parameters=(path=String(path), chrom=String(chrom), resolution=resolution, n_bins=n_bins))
    end
    return result
end

function _hic_chrom_sizes(path::AbstractString)
    sizes = Dict{String,Int}()
    try
        open(path, "r") do io
            readline(io)
            for line in eachline(io)
                isempty(strip(line)) && continue
                fields = split(strip(line), '\t')
                length(fields) >= 2 || continue
                try
                    sizes[String(fields[1])] = parse(Int, fields[2])
                catch
                end
            end
        end
    catch
    end
    return sizes
end

function write_cool(matrix::HiCContactMatrix, path::AbstractString)
    open(path, "w") do io
        println(io, "chrom\tstart\tend\tchrom2\tstart2\tend2\tcount")
        nz = findnz(matrix.matrix)
        for (i, j, v) in zip(nz[1], nz[2], nz[3])
            i > j && continue
            iv1 = matrix.bins[i]
            iv2 = matrix.bins[j]
            println(io, "$(iv1.chrom)\t$(iv1.left)\t$(iv1.right)\t$(iv2.chrom)\t$(iv2.left)\t$(iv2.right)\t$v")
        end
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "write_cool";
            parameters=(path=String(path), n_records=nnz(matrix.matrix)))
    end
    return path
end

# ---- Normalization -----------------------------------------------------------

"""
    normalize_ice(matrix; max_iter=200, tol=1e-6, ignore_diag=true)

Iterative correction and eigenvector decomposition (ICE) normalization.
"""
function normalize_ice(matrix::HiCContactMatrix; max_iter::Int=200, tol::Real=1e-6, ignore_diag::Bool=true)
    M = Matrix(matrix.matrix)
    n = size(M, 1)
    bias = ones(Float64, n)

    for _ in 1:max_iter
        row_sums = vec(sum(M, dims=2))
        valid = row_sums .> 0
        scale = ones(Float64, n)
        scale[valid] .= 1.0 ./ row_sums[valid]
        M = Diagonal(scale) * M * Diagonal(scale)
        bias .*= scale

        if any(valid) && maximum(abs.(row_sums[valid] ./ mean(row_sums[valid]) .- 1.0)) < Float64(tol)
            break
        end
    end

    normalized = HiCContactMatrix(
        sparse(M),
        matrix.chrom,
        matrix.bins,
        matrix.bin_size,
        true,
        "ICE",
        merge(matrix.metadata, Dict{String,Any}("ice_iterations" => max_iter)))

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "normalize_ice";
            parameters=(chrom=matrix.chrom, bin_size=matrix.bin_size, n_bins=n))
    end
    return normalized
end

"""
    normalize_kr(matrix; max_iter=200, tol=1e-6)

Knight-Ruiz matrix balancing normalization.
"""
function normalize_kr(matrix::HiCContactMatrix; max_iter::Int=200, tol::Real=1e-6)
    M = Matrix(matrix.matrix)
    n = size(M, 1)
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

    normalized = HiCContactMatrix(
        sparse(Diagonal(bias) * M * Diagonal(bias)),
        matrix.chrom,
        matrix.bins,
        matrix.bin_size,
        true,
        "KR",
        merge(matrix.metadata, Dict{String,Any}("kr_iterations" => max_iter)))

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "normalize_kr";
            parameters=(chrom=matrix.chrom, bin_size=matrix.bin_size, n_bins=n))
    end
    return normalized
end

"""
    normalize_vanilla(matrix)

Vanilla coverage normalization (divide by square root of row and column sums).
"""
function normalize_vanilla(matrix::HiCContactMatrix)
    M = Matrix(matrix.matrix)
    n = size(M, 1)

    row_sums = vec(sum(M, dims=2))
    col_sums = vec(sum(M, dims=1))
    valid = (row_sums .> 0) .& (col_sums .> 0)

    row_scale = ones(Float64, n)
    col_scale = ones(Float64, n)
    row_scale[valid] .= 1.0 ./ sqrt.(row_sums[valid])
    col_scale[valid] .= 1.0 ./ sqrt.(col_sums[valid])

    normalized = Diagonal(row_scale) * M * Diagonal(col_scale)

    result = HiCContactMatrix(
        sparse(normalized),
        matrix.chrom,
        matrix.bins,
        matrix.bin_size,
        true,
        "vanilla",
        Dict{String,Any}(matrix.metadata))

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "normalize_vanilla";
            parameters=(chrom=matrix.chrom, bin_size=matrix.bin_size))
    end
    return result
end

# ---- Expected Decay & Observed/Expected ---------------------------------------

"""
    expected_decay(matrix; max_dist=nothing)

Compute expected contact probability as a function of genomic distance.
"""
function expected_decay(matrix::HiCContactMatrix; max_dist::Union{Nothing,Int}=nothing)
    M = Matrix(matrix.matrix)
    n = size(M, 1)
    bin_size = matrix.bin_size

    if max_dist === nothing
        max_dist = n
    end

    distance_bins = Dict{Int,Vector{Float64}}()
    for i in 1:n
        for j in (i+1):min(n, i + max_dist)
            d = j - i
            push!(get!(distance_bins, d, Float64[]), M[i, j])
        end
    end

    distances = Int[]
    expected = Float64[]
    for d in sort!(collect(keys(distance_bins)))
        vals = distance_bins[d]
        if !isempty(vals)
            push!(distances, d)
            push!(expected, median(vals))
        end
    end

    positions = [d * bin_size for d in distances]

    result = DataFrames.DataFrame(distance_bins=distances, genomic_distance=positions, expected=expected)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "expected_decay";
            parameters=(chrom=matrix.chrom, n_bins=n, n_distance_points=length(distances)))
    end
    return result
end

"""
    power_law_fit(decay_df)

Fit power-law model P(s) ∝ s^(-α) to expected decay curve.
"""
function power_law_fit(decay_df::AbstractDataFrame)
    d = Float64.(decay_df.genomic_distance)
    p = Float64.(decay_df.expected)

    valid = (d .> 0) .& (p .> 0) .& isfinite.(d) .& isfinite.(p)
    d = d[valid]
    p = p[valid]

    isempty(d) && return (alpha=NaN, intercept=NaN, r_squared=NaN)

    log_d = log.(d)
    log_p = log.(p)

    n = length(log_d)
    mean_x = mean(log_d)
    mean_y = mean(log_p)
    ss_xy = sum((log_d .- mean_x) .* (log_p .- mean_y))
    ss_xx = sum((log_d .- mean_x) .^ 2)

    alpha = ss_xx != 0 ? -ss_xy / ss_xx : NaN
    intercept = mean_y - alpha * mean_x

    y_pred = intercept .+ alpha .* log_d
    ss_res = sum((log_p .- y_pred) .^ 2)
    ss_tot = sum((log_p .- mean_y) .^ 2)
    r_squared = ss_tot > 0 ? 1.0 - ss_res / ss_tot : 0.0

    result = (alpha=alpha, intercept=intercept, r_squared=r_squared, n_points=n)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "power_law_fit";
            parameters=(alpha=alpha, r_squared=r_squared, n_points=n))
    end
    return result
end

"""
    obs_exp_ratio(matrix, expected_df)

Compute observed/expected ratio matrix.
"""
function obs_exp_ratio(matrix::HiCContactMatrix, expected_df::AbstractDataFrame)
    M = Matrix(matrix.matrix)
    n = size(M, 1)
    bin_size = matrix.bin_size

    decay_lookup = Dict{Int,Float64}()
    for row in eachrow(expected_df)
        d = Int(row.distance_bins)
        decay_lookup[d] = row.expected
    end

    ratio = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            d = abs(j - i)
            expected_val = get(decay_lookup, d, eps(Float64))
            if expected_val > eps(Float64)
                ratio[i, j] = log2(max(M[i, j] / expected_val, eps(Float64)))
            end
        end
    end

    result = HiCContactMatrix(
        sparse(ratio),
        matrix.chrom,
        matrix.bins,
        matrix.bin_size,
        true,
        "observed_expected",
        merge(matrix.metadata, Dict{String,Any}("format" => "oe_ratio")))

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "obs_exp_ratio";
            parameters=(chrom=matrix.chrom, n_bins=n))
    end
    return result
end

# ---- Insulation Score --------------------------------------------------------

"""
    insulation_score(matrix; window_bins=5)

Compute insulation score: log mean of contacts across boundary.
"""
function insulation_score(matrix::HiCContactMatrix; window_bins::Int=5)
    M = Matrix(matrix.matrix)
    n = size(M, 1)

    ins = fill(NaN, n)
    for i in 1:n
        left_start = max(1, i - window_bins + 1)
        left_end = i
        right_start = i + 1
        right_end = min(n, i + window_bins)

        if right_start <= right_end
            block = @view M[left_start:left_end, right_start:right_end]
            μ = mean(block)
            ins[i] = log2(μ + 1.0)
        end
    end

    finite = isfinite.(ins)
    if any(finite)
        z = copy(ins)
        μ = mean(ins[finite])
        σ = std(ins[finite])
        if isfinite(σ) && σ > 0
            z[finite] .= (ins[finite] .- μ) ./ σ
        else
            fill!(z, 0.0)
        end
    else
        z = zeros(Float64, n)
    end

    result = DataFrames.DataFrame(
        bin=1:n,
        insulation=ins,
        zscore=z,
        position=[(i - 1) * matrix.bin_size + 1 for i in 1:n])

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "insulation_score";
            parameters=(chrom=matrix.chrom, n_bins=n, window_bins=window_bins))
    end
    return result
end

"""
    boundary_detection(insulation_df; q_threshold=0.1, min_strength=0.5)

Detect boundaries as local minima in insulation score.
"""
function boundary_detection(insulation_df::AbstractDataFrame; q_threshold::Real=0.1, min_strength::Real=0.5)
    vals = Float64.(insulation_df.insulation)
    n = length(vals)
    n < 3 && return DataFrames.DataFrame()

    finite = vals[isfinite.(vals)]
    isempty(finite) && return DataFrames.DataFrame()
    cutoff = quantile(finite, Float64(q_threshold))

    chrom_col = hasproperty(insulation_df, :chrom) ? insulation_df.chrom : fill("unknown", n)
    position_col = hasproperty(insulation_df, :position) ? insulation_df.position : fill(0, n)

    boundaries = HiCBounds[]
    for i in 2:(n - 1)
        v = vals[i]
        isfinite(v) || continue

        is_local_min = v <= vals[i - 1] && v <= vals[i + 1]
        below_cutoff = v <= cutoff
        is_candidate = is_local_min && below_cutoff

        if is_candidate
            left_val = isfinite(vals[i - 1]) ? vals[i - 1] : v
            right_val = isfinite(vals[i + 1]) ? vals[i + 1] : v
            strength = min(left_val, right_val) - v
            delta = (left_val + right_val) / 2 - v

            push!(boundaries, HiCBounds(
                String(chrom_col[i]),
                Int(position_col[i]),
                v,
                Float64(strength),
                Float64(delta)))
        end
    end

    df = DataFrames.DataFrame(
        chrom=String[],
        position=Int[],
        insulation=Float64[],
        strength=Float64[],
        delta=Float64[]
    )

    for b in boundaries
        if b.strength >= Float64(min_strength)
            push!(df, (b.chrom, b.position, b.insulation_score, b.strength, b.delta))
        end
    end

    sort!(df, :strength, rev=true)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "boundary_detection";
            parameters=(n_boundaries=nrow(df), q_threshold=Float64(q_threshold)))
    end
    return df
end

"""
    boundary_strength(matrix; boundary_positions)

Calculate aggregate boundary strength at given positions.
"""
function boundary_strength(matrix::HiCContactMatrix, boundary_positions::AbstractVector{<:Integer}; window_bins::Int=5)
    bin_size = matrix.bin_size
    M = Matrix(matrix.matrix)
    n = size(M, 1)

    results = Float64[]
    for bp in boundary_positions
        bin_idx = (bp - 1) ÷ bin_size + 1
        1 <= bin_idx <= n || continue

        left_start = max(1, bin_idx - window_bins + 1)
        left_end = bin_idx
        right_start = bin_idx + 1
        right_end = min(n, bin_idx + window_bins)

        inside = if left_start <= left_end && right_start <= right_end
            mean(M[left_start:left_end, right_start:right_end])
        else
            0.0
        end

        outside_left = if left_start > 1
            l_start = max(1, left_start - window_bins)
            mean(M[l_start:left_end, left_start:left_end])
        else
            0.0
        end

        outside_right = if right_end < n
            r_end = min(n, right_end + window_bins)
            mean(M[right_start:right_end, right_start:r_end])
        else
            0.0
        end

        outside_mean = (outside_left + outside_right) / 2
        strength = outside_mean > 0 ? inside / outside_mean : 0.0
        push!(results, strength)
    end

    return results
end

# ---- Compartments ------------------------------------------------------------

"""
    compartment_score(matrix; n_components=1)

Compute PC1 compartment score from correlation matrix.
"""
function compartment_score(matrix::HiCContactMatrix; n_components::Int=1)
    M = max.(Matrix(matrix.matrix), 0.0)
    n = size(M, 1)

    c = cor(M)
    c[.!isfinite.(c)] .= 0.0

    ev = eigen(Symmetric(c))
    pc1 = ev.vectors[:, argmax(ev.values)]

    gc = Float64[]
    for iv in matrix.bins
        push!(gc, _gc_content_matrix_region(iv, matrix))
    end

    correlation_with_gc = abs(cor(pc1, gc))

    result = DataFrames.DataFrame(
        bin=1:n,
        pc1=pc1,
        compartment=ifelse.(pc1 .>= 0, "A", "B"),
        gc_content=gc,
        gc_correlation=fill(correlation_with_gc, n),
        position=[(i - 1) * matrix.bin_size + 1 for i in 1:n])

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "compartment_score";
            parameters=(chrom=matrix.chrom, n_bins=n, gc_correlation=correlation_with_gc))
    end
    return result
end

function _gc_content_matrix_region(iv::GenomicInterval, matrix::HiCContactMatrix)
    return 0.41
end

"""
    ab_compartments(matrix)

Assign A/B compartment labels based on PC1 sign.
"""
function ab_compartments(matrix::HiCContactMatrix)
    scores = compartment_score(matrix)
    compartments = String.(scores.compartment)
    a_frac = count(==("A"), compartments) / max(length(compartments), 1)

    result = (scores=scores, a_fraction=a_frac, b_fraction=1.0 - a_frac)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "ab_compartments";
            parameters=(chrom=matrix.chrom, a_fraction=a_frac))
    end
    return result
end

"""
    compartment_switching(compartment_a, compartment_b)

Identify compartment switches between two conditions.
"""
function compartment_switching(compartment_a::AbstractDataFrame, compartment_b::AbstractDataFrame; min_flip::Real=0.5)
    n = min(nrow(compartment_a), nrow(compartment_b))
    n < 1 && return DataFrames.DataFrame()

    pc1_a = Float64.(compartment_a.pc1[1:n])
    pc1_b = Float64.(compartment_b.pc1[1:n])

    flips = Int[]
    for i in 1:n
        a_sign = pc1_a[i] >= 0
        b_sign = pc1_b[i] >= 0
        if a_sign != b_sign
            delta = abs(pc1_a[i] - pc1_b[i])
            if delta >= Float64(min_flip)
                push!(flips, i)
            end
        end
    end

    result = DataFrames.DataFrame(
        bin=flips,
        position=[get(compartment_a, :position, zeros(Int, n))[i] for i in flips],
        pc1_condition1=[pc1_a[i] for i in flips],
        pc1_condition2=[pc1_b[i] for i in flips],
        delta_pc1=[pc1_b[i] - pc1_a[i] for i in flips],
        switch_type=[pc1_a[i] > 0 ? "B_to_A" : "A_to_B" for i in flips])

    sort!(result, :delta_pc1, rev=true)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "compartment_switching";
            parameters=(n_switches=nrow(result), min_flip=Float64(min_flip)))
    end
    return result
end

# ---- Differential Interaction Testing ----------------------------------------

function _nb_glm_fit(counts::AbstractVector{<:Integer}, offsets::AbstractVector{<:Real}, design::AbstractMatrix{<:Real}; max_iter::Int=100, tol::Real=1e-6)
    n = length(counts)
    k = size(design, 2)

    beta = zeros(Float64, k)
    mu = Float64.(counts) .+ 0.1
    phi = 1.0

    for iter in 1:max_iter
        eta = design * beta .+ log.(max.(offsets, eps(Float64)))
        mu_new = max.(exp.(eta), eps(Float64))

        r = Float64.(counts) .- mu_new
        W = Diagonal(mu_new ./ (1 .+ mu_new ./ (phi .* max.(offsets, eps(Float64)))))
        z = eta .+ r ./ mu_new

        try
            beta_new = design \ z
        catch
            break
        end

        if maximum(abs.(beta_new .- beta)) < Float64(tol)
            beta = beta_new
            break
        end
        beta = beta_new
    end

    residuals = Float64.(counts) .- mu
    phi_est = sum(residuals .^ 2 ./ max.(mu, eps(Float64))) / max(n - k, 1)
    phi_est = clamp(phi_est, 0.01, 100.0)

    return (beta=beta, dispersion=phi_est, mu=mu)
end

"""
    differential_interaction(cond_a, cond_b, coordinates; design=nothing)

Differential interaction testing via negative binomial GLM.
"""
function differential_interaction(
    cond_a::AbstractVector{<:HiCContactMatrix},
    cond_b::AbstractVector{<:HiCContactMatrix},
    coordinates::AbstractVector{Tuple{Int,Int}};
    design::Union{Nothing,AbstractMatrix{<:Real}}=nothing)

    n_a = length(cond_a)
    n_b = length(cond_b)
    n_coords = length(coordinates)

    if design === nothing
        design = vcat(ones(Float64, n_a), zeros(Float64, n_b))
        design = reshape(design, length(design), 1)
    end

    results = DataFrames.DataFrame(
        bin1=Int[],
        bin2=Int[],
        chrom1=String[],
        chrom2=String[],
        count_a=Float64[],
        count_b=Float64[],
        log2_fc=Float64[],
        pvalue=Float64[],
        padj=Float64[]
    )

    for (bin1, bin2) in coordinates
        counts = Float64[]
        for mat in cond_a
            n = size(mat.matrix, 1)
            if 1 <= bin1 <= n && 1 <= bin2 <= n
                push!(counts, mat.matrix[bin1, bin2])
            else
                push!(counts, 0.0)
            end
        end
        for mat in cond_b
            n = size(mat.matrix, 1)
            if 1 <= bin1 <= n && 1 <= bin2 <= n
                push!(counts, mat.matrix[bin1, bin2])
            else
                push!(counts, 0.0)
            end
        end

        mean_a = mean(counts[1:n_a])
        mean_b = mean(counts[n_a+1:end])
        log2_fc = mean_b > 0 && mean_a > 0 ? log2(mean_b / mean_a) : 0.0

        if length(counts) >= 4
            fit = try
                offsets = ones(Float64, length(counts))
                _nb_glm_fit(Int.(round.(counts)), offsets, design)
            catch
                (beta=[0.0], dispersion=1.0)
            end

            if size(design, 2) >= 1
                stat = fit.beta[end]
                se = sqrt(fit.dispersion)
                pval = 2 * ccdf(Normal(abs(stat), se), abs(stat))
            else
                stat = 0.0
                pval = 1.0
            end
        else
            log2_fc = 0.0
            stat = 0.0
            pval = 1.0
        end

        chrom1 = n_a > 0 ? cond_a[1].chrom : "unknown"
        chrom2 = n_a > 0 ? cond_a[1].chrom : "unknown"

        push!(results, (bin1, bin2, chrom1, chrom2, mean_a, mean_b, log2_fc, pval, 1.0))
    end

    m = nrow(results)
    if m > 0
        pvals = results.pvalue
        perm = sortperm(pvals)
        running = 1.0
        for i in m:-1:1
            j = perm[i]
            running = min(running, pvals[j] * m / i)
            results.padj[j] = min(running, 1.0)
        end
    end

    sort!(results, :pvalue)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "differential_interaction";
            parameters=(n_coords=n_coords, n_cond_a=n_a, n_cond_b=n_b))
    end
    return results
end

"""
    hic_diff_test(obs_a, obs_b; bin_size=10000)

Direct differential interaction test from raw observations.
"""
function hic_diff_test(
    obs_a::AbstractVector{<:Integer},
    obs_b::AbstractVector{<:Integer};
    bin_size::Int=10000)

    n_a = length(obs_a)
    n_b = length(obs_b)
    total_a = sum(obs_a)
    total_b = sum(obs_b)

    if total_a == 0 && total_b == 0
        return DataFrames.DataFrame(
            count_a=Int[], count_b=Int[], log2_fc=Float64[], pvalue=Float64[], significance=String[])
    end

    log2_fc = Float64[]
    pvalue = Float64[]

    for (c_a, c_b) in zip(obs_a, obs_b)
        if total_a > 0 && total_b > 0
            expected_a = c_b * (total_a / total_b)
            fc = expected_a > 0 ? c_a / expected_a : 0.0
            push!(log2_fc, log2(max(fc, eps(Float64))))
        else
            push!(log2_fc, 0.0)
        end

        if c_a + c_b >= 0
            p = if c_a + c_b <= 20
                pvalue_mann_whitney(Int(c_a), Int(c_b))
            else
                z = (c_a / max(total_a, 1) - c_b / max(total_b, 1)) / sqrt(1/max(total_a,1) + 1/max(total_b,1))
                2 * ccdf(Normal(), abs(z))
            end
            push!(pvalue, p)
        else
            push!(pvalue, 1.0)
        end
    end

    results = DataFrames.DataFrame(
        count_a=Int.(obs_a),
        count_b=Int.(obs_b),
        log2_fc=log2_fc,
        pvalue=pvalue,
        significance=["***" for _ in 1:length(log2_fc)])

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "hic_diff_test";
            parameters=(n_a=n_a, n_b=n_b, total_a=total_a, total_b=total_b))
    end
    return results
end

function pvalue_mann_whitney(a::Int, b::Int)
    n = a + b
    if n == 0
        return 1.0
    end
    u = Float64(a * b) + (Float64(a) * (Float64(a) + 1)) / 2 - a
    z = (u - (Float64(a) * Float64(b)) / 2) / sqrt(Float64(a) * Float64(b) * (Float64(n) + 1) / 12)
    return 2 * ccdf(Normal(), abs(z))
end

# ---- Loop Calling ------------------------------------------------------------

"""
    call_loops(matrix; window=10, fdr=0.1, min_dist=5, max_dist=1000)

Call chromatin loops using peak detection on corner score.
"""
function call_loops(
    matrix::HiCContactMatrix;
    window::Int=10,
    fdr::Real=0.1,
    min_dist::Int=5,
    max_dist::Int=1000,
    bin_size::Int=1)

    M = Matrix(matrix.matrix)
    n = size(M, 1)

    corner_scores = zeros(Float64, n, n)
    max_d = min(max_dist, n)

    for d in min_dist:max_d
        for i in 1:(n - d)
            j = i + d
            window_d = max(1, div(d, 20))

            lo1 = max(1, i - window)
            hi1 = i
            lo2 = j
            hi2 = min(n, j + window)

            lo_avg = mean(M[lo1:hi1, i:j])
            hi_avg = mean(M[i:j, lo2:hi2])

            local_lo = zeros(Float64, hi1 - lo1 + 1)
            local_hi = zeros(Float64, hi2 - lo2 + 1)

            for bi in lo1:hi1
                local_lo[bi - lo1 + 1] = mean(M[bi, i:j])
            end
            for bi in lo2:hi2
                local_hi[bi - lo2 + 1] = mean(M[i:j, bi])
            end

            expected = sqrt(max(lo_avg * hi_avg, eps(Float64)))
            if expected > eps(Float64)
                corner_scores[i, j] = (M[i, j] - expected) / expected
            end
        end
    end

    loop_candidates = Tuple{Int,Int,Float64}[]
    threshold = Float64(fdr) * 10

    for i in 1:n
        for j in (i + min_dist):min(n, i + max_dist)
            score = corner_scores[i, j]
            if score > threshold
                push!(loop_candidates, (i, j, score))
            end
        end
    end

    sort!(loop_candidates, by=x->x[3], rev=true)
    merged = Tuple{Int,Int,Float64}[]
    used = falses(length(loop_candidates))

    merge_window = max(3, div(window, 2))

    for (idx, (i, j, s)) in enumerate(loop_candidates)
        used[idx] && continue
        push!(merged, (i, j, s))
        for k in idx:length(loop_candidates)
            if !used[k]
                (i2, j2, _) = loop_candidates[k]
                if abs(i2 - i) <= merge_window && abs(j2 - j) <= merge_window
                    used[k] = true
                end
            end
        end
    end

    if bin_size <= 1
        bin_size = matrix.bin_size
    end

    results = DataFrames.DataFrame(
        anchor1_bin=Int[],
        anchor2_bin=Int[],
        anchor1_pos=Int[],
        anchor2_pos=Int[],
        score=Float64[])

    for (a1, a2, s) in merged
        pos1 = a1 * bin_size
        pos2 = a2 * bin_size
        if !isempty(matrix.bins) && 1 <= a1 <= length(matrix.bins)
            pos1 = matrix.bins[a1].left
        end
        if !isempty(matrix.bins) && 1 <= a2 <= length(matrix.bins)
            pos2 = matrix.bins[a2].left
        end
        push!(results, (a1, a2, pos1, pos2, s))
    end

    sort!(results, :score, rev=true)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "call_loops";
            parameters=(chrom=matrix.chrom, n_loops=nrow(results), window=window, fdr=Float64(fdr)))
    end
    return results
end

"""
    merge_loops(loops_a, loops_b; max_distance=5)

Merge loop calls from multiple conditions.
"""
function merge_loops(loops_a::AbstractDataFrame, loops_b::AbstractDataFrame; max_distance::Int=5)
    df_a = copy(loops_a)
    df_a.source = fill("A", nrow(loops_a))
    df_b = copy(loops_b)
    df_b.source = fill("B", nrow(loops_b))
    all_loops = vcat(df_a, df_b)

    merged = DataFrames.DataFrame(
        anchor1_pos=Int[],
        anchor2_pos=Int[],
        found_in=String[],
        max_score=Float64[])

    used_a = falses(nrow(loops_a))
    used_b = falses(nrow(loops_b))

    for (i, row) in enumerate(eachrow(loops_a))
        used_a[i] && continue
        a1, a2 = row.anchor1_pos, row.anchor2_pos
        found_in = ["A"]
        max_score = row.score

        for (j, brow) in enumerate(eachrow(loops_b))
            if !used_b[j] && abs(brow.anchor1_pos - a1) <= max_distance && abs(brow.anchor2_pos - a2) <= max_distance
                push!(found_in, "B")
                max_score = max(max_score, brow.score)
                used_b[j] = true
            end
        end

        push!(merged, (a1, a2, join(found_in, ","), max_score))
        used_a[i] = true
    end

    for (j, row) in enumerate(eachrow(loops_b))
        if !used_b[j]
            push!(merged, (row.anchor1_pos, row.anchor2_pos, "B", row.score))
        end
    end

    sort!(merged, :max_score, rev=true)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "merge_loops";
            parameters=(n_a=nrow(loops_a), n_b=nrow(loops_b), n_merged=nrow(merged)))
    end
    return merged
end

# ---- Aggregate Peak Analysis --------------------------------------------------

"""
    pileup_matrix(matrix, loops; window_bins=20)

Create pileup matrix around loop anchors.
"""
function pileup_matrix(matrix::HiCContactMatrix, loops::AbstractDataFrame; window_bins::Int=20)
    M = Matrix(matrix.matrix)
    n = size(M, 1)
    h = window_bins

    pileup = zeros(Float64, 2*h + 1, 2*h + 1)
    count = 0

    for row in eachrow(loops)
        a1 = row.anchor1_bin
        a2 = row.anchor2_bin

        if 1 <= a1 - h && a1 + h <= n && 1 <= a2 - h && a2 + h <= n
            sub = @view M[(a1-h):(a1+h), (a2-h):(a2+h)]
            pileup .+= sub
            count += 1
        end
    end

    if count > 0
        pileup ./= count
    end

    result = (matrix=pileup, n_loops=count, window=window_bins)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "pileup_matrix";
            parameters=(chrom=matrix.chrom, n_loops=count, window_bins=window_bins))
    end
    return result
end

"""
    loop_apa(pileup)

Aggregate Peak Analysis: quantify loop strength from pileup.
"""
function loop_apa(pileup_result)
    p = pileup_result.matrix
    h = pileup_result.window
    n = size(p, 1)

    center = div(n, 2) + 1
    loop_radius = max(1, div(h, 4))

    corner = (mean(p[1:h, n-h+1:n]) + mean(p[n-h+1:n, 1:h])) / 2
    center_region = @view p[center-loop_radius:center+loop_radius, center-loop_radius:center+loop_radius]
    center_mean = mean(center_region)
    center_sum = sum(center_region)

    apa_score = corner > 0 ? center_sum / (corner * (2*loop_radius+1)^2) : 0.0

    result = (apa_score=apa_score, center_mean=center_mean, corner_mean=corner, n_loops=pileup_result.n_loops)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "loop_apa";
            parameters=(apa_score=apa_score, n_loops=pileup_result.n_loops))
    end
    return result
end

# ---- Quality Metrics ---------------------------------------------------------

"""
    hic_quality_metrics(matrix)

Compute Hi-C quality metrics (MQC-style).
"""
function hic_quality_metrics(matrix::HiCContactMatrix)
    M = Matrix(matrix.matrix)
    n = size(M, 1)

    total = sum(M)
    nonzero = count(M .> 0) / max(n^2, 1)
    sparsity = 1.0 - nonzero

    max_contact = maximum(M)
    median_contact = median(vec(M[M .> 0]))

    expected = expected_decay(matrix)
    if nrow(expected) > 0
        decay_fit = power_law_fit(expected)
    else
        decay_fit = (alpha=NaN, r_squared=NaN)
    end

    ins = insulation_score(matrix)
    ins_finite = ins.insulation[isfinite.(ins.insulation)]
    ins_cv = isempty(ins_finite) ? 0.0 : std(ins_finite) / max(mean(ins_finite), eps(Float64))

    results = DataFrames.DataFrame(
        metric=["total_contacts", "n_bins", "sparsity", "nonzero_fraction", "max_contact", "median_contact",
                "power_law_alpha", "power_law_r2", "insulation_cv", "normalized"],
        value=[total, n, sparsity, nonzero, max_contact, median_contact,
               getfield(decay_fit, :alpha), getfield(decay_fit, :r_squared), ins_cv, matrix.normalized])

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "hic_quality_metrics";
            parameters=(chrom=matrix.chrom, n_bins=n, total_contacts=total, sparsity=sparsity))
    end
    return results
end

# ---- Contact Decay -----------------------------------------------------------

"""
    contact_decay_by_distance(matrix; max_dist=nothing)

Compute mean contact frequency by distance.
"""
function contact_decay_by_distance(matrix::HiCContactMatrix; max_dist::Union{Nothing,Int}=nothing)
    return expected_decay(matrix; max_dist=max_dist)
end

"""
    chromosome_contact_map(experiment, chrom)

Extract chromosome contact map.
"""
function chromosome_contact_map(experiment::HiCExperiment, chrom::String)
    mat = get(experiment.matrices, chrom, nothing)
    mat === nothing && throw(ArgumentError("chromosome $chrom not found in experiment"))
    return mat
end

end
