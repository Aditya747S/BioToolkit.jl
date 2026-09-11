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
using HDF5
using CodecZlib
using Dates

using ..GenomicRanges: GenomicInterval, IntervalCollection, build_collection
using ..BioToolkit: BioSequence, DNAAlphabet, BlenderIntegrator
using ..BlenderIntegrator: BlenderHiCPayload, BlenderMaterial, to_blender_payload

using ..BioToolkit: ProvenanceContext, ThreadSafeProvenanceContext, active_provenance_context
using ..BioToolkit: provenance_result!, provenance_parent_ids, provenance_record, register_provenance!

@inline function _register_hic_result!(_ctx, result, operation; parents=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export read_cool, read_mcool, read_hic, write_hic, write_cool, write_contacts_tsv, read_contacts_tsv
export read_cool, read_mcool, read_hic, write_cool, write_contacts_tsv, read_contacts_tsv
export normalize_ice, normalize_kr, normalize_vanilla
export expected_decay, power_law_fit, obs_exp_ratio
export insulation_score, boundary_detection, boundary_strength
export compartment_score, ab_compartments, compartment_switching
export differential_interaction, hic_diff_test
export call_loops, merge_loops, loop_apa
export GInteraction, InteractionSet, interactions_from_matrix, interactions_from_pairs
export interaction_counts, anchors1, anchors2, swap_anchors, cis_interactions, trans_interactions
export find_interactions, coarsen, write_mcool
export read_pairs, write_pairs, read_hicpro, write_hicpro
export saddle_plot, compartment_strength, call_tads, virtual_4c, aggregate_at_features, regions
export hic_loess_normalize
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

function BlenderIntegrator.to_blender_payload(hic::HiCContactMatrix; tube_radius::Float64=0.2, name::String="HiC_3D_Spline")
    n_bins = size(hic.matrix, 1)
    # 3D layout from the data: classical MDS on a 1 - O/E contact similarity
    # (previously the payload placed bins on a parametric helix unrelated to
    # the contact matrix). Large maps are subsampled to <= 400 bins for the
    # eigendecomposition and interpolation is omitted (sparse representative
    # bins are used directly).
    pts = zeros(Float64, n_bins, 3)
    if n_bins >= 3
        idx = n_bins <= 400 ? collect(1:n_bins) : round.(Int, range(1, n_bins; length=400))
        k = length(idx)
        sub = Matrix(hic.matrix[idx, idx])
        dmax = maximum(sub)
        dmax > 0 || (dmax = 1.0)
        D = 1.0 .- sub ./ dmax
        D[diagind(D)] .= 0.0
        D .= (D .+ D') ./ 2.0
        J = Matrix{Float64}(I, k, k) .- 1.0 / k
        B = -0.5 * J * (D .^ 2) * J
        ev = eigen(Symmetric(B))
        order = sortperm(ev.values; rev=true)
        for (dim, oi) in enumerate(order[1:min(3, k)])
            vals = sqrt.(max.(ev.values[oi], 0.0))
            for (row, bi) in enumerate(idx)
                pts[bi, dim] = vals[row]
            end
        end
    end
    scores = vec(sum(hic.matrix, dims=2))
    anchors = [(1, min(10, n_bins))]
    mat = BlenderMaterial(name=name * "_mat", color=(0.7, 0.2, 0.9, 1.0), roughness=0.3)
    return BlenderHiCPayload(name, pts, scores, anchors, tube_radius, mat)
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

# Cooler single-resolution (.cool) and multi-resolution (.mcool) HDF5 readers
# per the cooler specification (/chroms, /bins, /pixels, /indexes), plus a
# documented TSV fallback for generic 5-column contact tables (the previous
# behaviour, which mis-presented itself as a .cool reader).
_is_hdf5_file(path::AbstractString) = let bytes = read(path, 8); length(bytes) >= 8 && bytes == UInt8[0x89, UInt8('H'), UInt8('D'), UInt8('F'), 0x0d, 0x0a, 0x1a, 0x0a]; end

function read_cool(path::AbstractString; chrom::Union{Nothing,String}=nothing, resolution::Union{Nothing,Int}=nothing)
    isfile(path) || throw(ArgumentError(".cool file not found: $path"))
    if _is_hdf5_file(path)
        return _read_cool_hdf5(path; chrom=chrom, resolution=resolution, group="")
    end
    # TSV fallback (documented): generic chrom/pos/chrom2/pos/count table.
    @warn "read_cool: $path is not HDF5; falling back to the generic 5-column TSV contact reader (real .cool files are HDF5)"
    return _read_contacts_tsv(path; chrom=chrom, resolution=resolution)
end

"""
    read_mcool(path; resolution=10000, chrom=nothing)

Read one resolution of a multi-resolution `.mcool` (HDF5) file.
"""
function read_mcool(path::AbstractString; resolution::Int=10000, chrom::Union{Nothing,String}=nothing)
    isfile(path) || throw(ArgumentError(".mcool file not found: $path"))
    _is_hdf5_file(path) || throw(ArgumentError(".mcool must be HDF5: $path"))
    return _read_cool_hdf5(path; chrom=chrom, resolution=resolution, group="resolutions/$(resolution)/")
end

function _read_cool_hdf5(path::AbstractString; chrom::Union{Nothing,String}=nothing, resolution::Union{Nothing,Int}=nothing, group::String="")
    HDF5.h5open(path, "r") do h5
        g = isempty(group) ? h5 : h5[group]
        chrom_names = Vector{String}(read(g["chroms/names"]))
        chrom_lengths = Vector{Int}(read(g["chroms/lengths"]))
        bin_chrom = Vector{String}(read(g["bins/chrom"]))
        bin_start = Vector{Int}(read(g["bins/start"]))
        bin_end = Vector{Int}(read(g["bins/end"]))
        bin1_id = Vector{Int}(read(g["pixels/bin1_id"]))
        bin2_id = Vector{Int}(read(g["pixels/bin2_id"]))
        counts = Vector{Float64}(read(g["pixels/count"]))

        bin_size = 0
        if haskey(HDF5.attrs(g), "bin-size")
            bin_size = Int(HDF5.attrs(g)["bin-size"])
        end

        # Restrict to one chromosome and build a within-chromosome matrix.
        target = chrom === nothing ? String(bin_chrom[1]) : String(chrom)
        keep = findall(==(target), bin_chrom)
        isempty(keep) && throw(ArgumentError("chromosome $target not present in $path"))
        local_index = Dict{Int,Int}(gid => i for (i, gid) in enumerate(keep))
        nb = length(keep)

        rows = Int[]; cols = Int[]; vals = Float64[]
        for (b1, b2, v) in zip(bin1_id, bin2_id, counts)
            i1 = get(local_index, b1 + 1, 0)   # cooler bin ids are 0-based
            i2 = get(local_index, b2 + 1, 0)
            (i1 > 0 && i2 > 0) || continue
            push!(rows, i1); push!(cols, i2); push!(vals, Float64(v))
        end
        matrix = sparse(rows, cols, vals, nb, nb)
        matrix = matrix + matrix' - Diagonal(diag(matrix))

        intervals = GenomicInterval[GenomicInterval(target, bin_start[i] + 1, bin_end[i], '.') for i in keep]
        if bin_size == 0 && nb > 1
            bin_size = Int(bin_start[keep[2]] - bin_start[keep[1]])
        end
        bin_size = resolution !== nothing ? Int(resolution) : max(bin_size, 1)

        chrom_size = get(Dict(zip(chrom_names, chrom_lengths)), target, maximum(bin_end[keep]))

        result = HiCContactMatrix(
            matrix,
            target,
            intervals,
            bin_size,
            false,
            "none",
            Dict{String,Any}("source" => String(path), "format" => "cool",
                             "chrom_sizes" => Dict(zip(chrom_names, chrom_lengths)),
                             "chromosome_size" => chrom_size))

        _ctx = active_provenance_context()
        if _ctx !== nothing
            register_provenance!(_ctx, "read_cool";
                parameters=(path=String(path), chrom=target, bin_size=bin_size, n_bins=nb, format="hdf5"))
        end
        return result
    end
end

"""
    read_contacts_tsv(path; chrom=nothing, resolution=nothing)

Read a generic 5-column contact TSV (chrom pos chrom2 pos count) — the format
the old `read_cool` actually parsed.
"""
function read_contacts_tsv(path::AbstractString; chrom::Union{Nothing,String}=nothing, resolution::Union{Nothing,Int}=nothing)
    isfile(path) || throw(ArgumentError("contact TSV not found: $path"))
    return _read_contacts_tsv(path; chrom=chrom, resolution=resolution)
end

function _read_contacts_tsv(path::AbstractString; chrom::Union{Nothing,String}=nothing, resolution::Union{Nothing,Int}=nothing)
    bin_size = resolution !== nothing ? Int(resolution) : _detect_bin_size_cool(path)
    chrom_data = _cool_chrom_data(path)

    if chrom !== nothing
        chrom_data = filter(c -> c[1] == String(chrom), chrom_data)
    end

    rows = Int[]
    cols = Int[]
    vals = Float64[]
    first_chrom = nothing

    open(path, "r") do io
        header_line = readline(io)
        while !eof(io)
            line = readline(io)
            isempty(strip(line)) && continue
            fields = split(strip(line), '\t')
            length(fields) >= 3 || continue

            # Support both the 5-column generic layout (chrom pos chrom2 pos
            # count) and the 7-column layout the old writer produced (chrom
            # start end chrom2 start2 end2 count) — previously these were
            # incompatible, so round trips crashed. Column 4 is numeric in
            # the 5-column layout (pos2) and a chromosome name in the
            # 7-column layout.
            local c1, p1, c2, p2, count
            if length(fields) >= 7 && !occursin(r"^\d+$", fields[4])
                c1 = fields[1]; p1 = parse(Int, fields[2]); c2 = fields[4]; p2 = parse(Int, fields[5]); count = parse(Float64, fields[7])
            else
                c1 = fields[1]; p1 = parse(Int, fields[2]); c2 = fields[3]; p2 = parse(Int, fields[4]); count = parse(Float64, fields[5])
            end

            c1 == c2 || continue
            chrom !== nothing && c1 != String(chrom) && continue

            bin1 = (p1 - 1) ÷ bin_size + 1
            bin2 = (p2 - 1) ÷ bin_size + 1

            first_chrom === nothing && (first_chrom = String(c1))
            push!(rows, bin1)
            push!(cols, bin2)
            push!(vals, count)
        end
    end

    n_bins = isempty(rows) ? 0 : max(maximum(rows), maximum(cols))
    matrix = sparse(rows, cols, vals, n_bins, n_bins)
    matrix = matrix + matrix' - Diagonal(diag(matrix))

    # Build bins from the observed bin extent (one chromosome per TSV). The
    # previous chrom-size-file machinery zeroed out the matrix whenever the
    # TSV lacked a usable chrom-size table, silently dropping every contact.
    target_chrom = first_chrom !== nothing ? first_chrom :
                   (chrom !== nothing ? String(chrom) : "unknown")
    bin_intervals = GenomicInterval[GenomicInterval(
        target_chrom, (b - 1) * bin_size + 1, b * bin_size, '.') for b in 1:n_bins]
    sub_matrix = matrix
    sub_intervals = bin_intervals

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
    # NOTE: `return` inside an `open(...) do` block only exits the do-closure
    # and is discarded by the caller — the original version therefore always
    # returned the default 10000. The result is captured explicitly instead.
    detected = 0
    try
        open(path, "r") do io
            header = readline(io)
            for line in eachline(io)
                fields = split(strip(line), '\t')
                isnum(f) = try parse(Int, f); true catch; false end
                # 7-column layout carries the exact bin width (end - start + 1).
                # 5-column layouts only expose positions, so the smallest
                # nonzero observed |pos2 - pos1| is used (approximate).
                if length(fields) >= 7 && isnum(fields[2]) && !isnum(fields[4]) && isnum(fields[3])
                    detected = max(1, parse(Int, fields[3]) - parse(Int, fields[2]) + 1)
                    break
                elseif length(fields) >= 5 && isnum(fields[2]) && isnum(fields[4]) && fields[4] != fields[2]
                    detected = max(1, abs(parse(Int, fields[4]) - parse(Int, fields[2])))
                    break
                end
            end
        end
    catch
    end
    return detected > 0 ? detected : 10000
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

# Generic 5-column TSV fallback (the original read_hic behaviour), kept under
# its own name and used by the binary reader's dispatcher for non-HIC files.
function _read_hic_tsv(path::AbstractString; chrom::String, resolution::Int=10000, start::Int=1, stop::Union{Nothing,Int}=nothing)
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
    # Accepts chrom-size lines (exactly 2 fields: name + integer length).
    # Contact lines (>= 3 fields) are ignored, so a mixed file works without
    # assuming a header — the unconditional first-line skip previously ate
    # the chrom-size line when no header was present.
    sizes = Dict{String,Int}()
    try
        open(path, "r") do io
            for line in eachline(io)
                isempty(strip(line)) && continue
                fields = split(strip(line), '\t')
                length(fields) == 2 || continue
                ischrom = !occursin(r"^\d+$", fields[1])
                isnum = occursin(r"^\d+$", fields[2])
                (ischrom && isnum) || continue
                sizes[String(fields[1])] = parse(Int, fields[2])
            end
        end
    catch
    end
    return sizes
end

"""
    write_cool(matrix, path)

Write a cooler-specification HDF5 `.cool` file (/chroms, /bins, /pixels,
/indexes with chrom_offset and bin1_offset). The previous TSV output is still
available as `write_contacts_tsv`.
"""
function write_cool(matrix::HiCContactMatrix, path::AbstractString)
    nb = size(matrix.matrix, 1)
    chrom_size = isempty(matrix.bins) ? nb * matrix.bin_size : maximum(iv.right for iv in matrix.bins)
    chrom_name = matrix.chrom

    nz = findnz(matrix.matrix)
    # Cooler stores each contact once: pixels with bin1_id <= bin2_id only.
    keep_px = nz[1] .<= nz[2]
    # Global 0-based bin ids (single-chromosome matrix: ids == within-chrom ids).
    bin1_ids = Int.(nz[1][keep_px]) .- 1
    bin2_ids = Int.(nz[2][keep_px]) .- 1
    counts = Float64.(nz[3][keep_px])

    # bin1_offset: cumulative non-zero pixels per bin (cooler index).
    bin1_offset = zeros(Int, nb + 1)
    for b in bin1_ids
        bin1_offset[min(b + 2, nb + 1)] += 1
    end
    cumsum!(bin1_offset, bin1_offset)

    HDF5.h5open(path, "w") do h5
        HDF5.attrs(h5)["format"] = "HDF5::Cooler"
        HDF5.attrs(h5)["format-version"] = "0.1.1"
        HDF5.attrs(h5)["bin-size"] = matrix.bin_size
        HDF5.attrs(h5)["nbins"] = nb
        HDF5.attrs(h5)["nchroms"] = 1
        HDF5.attrs(h5)["npixels"] = length(counts)
        HDF5.attrs(h5)["genome-assembly"] = "unknown"

        g = HDF5.create_group(h5, "chroms")
        write(g, "names", String[chrom_name])
        write(g, "lengths", Int[chrom_size])

        g = HDF5.create_group(h5, "bins")
        write(g, "chrom", String[chrom_name for _ in 1:nb])
        write(g, "start", Int[(i - 1) * matrix.bin_size for i in 1:nb])
        write(g, "end", Int[min(i * matrix.bin_size, chrom_size) for i in 1:nb])

        g = HDF5.create_group(h5, "pixels")
        write(g, "bin1_id", bin1_ids)
        write(g, "bin2_id", bin2_ids)
        write(g, "count", counts)

        g = HDF5.create_group(h5, "indexes")
        write(g, "chrom_offset", Int[0, nb])
        write(g, "bin1_offset", bin1_offset)
    end

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "write_cool";
            parameters=(path=String(path), n_records=length(counts), format="hdf5"))
    end
    return path
end

"""
    write_contacts_tsv(matrix, path)

Write contacts as a generic 5-column TSV (chrom pos chrom2 pos count) — the
format the old `write_cool` actually produced.
"""
function write_contacts_tsv(matrix::HiCContactMatrix, path::AbstractString)
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
        register_provenance!(_ctx, "write_contacts_tsv";
            parameters=(path=String(path), n_records=nnz(matrix.matrix)))
    end
    return path
end

# ---- Normalization -----------------------------------------------------------

"""
    normalize_ice(matrix; max_iter=200, tol=1e-6, ignore_diag=true)

Iterative correction and eigenvector decomposition (ICE) normalization.
"""
function normalize_ice(matrix::HiCContactMatrix; max_iter::Int=200, tol::Real=1e-6, ignore_diag::Bool=true,
                       blacklist_bins::AbstractVector{<:Integer}=Int[],
                       weights::Union{Nothing,AbstractVector{<:Real}}=nothing)
    M = Matrix(matrix.matrix)
    n = size(M, 1)
    bias = ones(Float64, n)
    original_diag = diag(M)
    # ignore_diag=true (ICE standard): the diagonal is masked out of the
    # balancing iterations and zeroed in the balanced matrix (previously the
    # kwarg was accepted and silently ignored).
    ignore_diag && (M[diagind(M)] .= 0.0)

    # Blacklisted bins and per-bin weights (weighted ICE): excluded bins are
    # masked out of the balancing system and receive a zero bias, so their
    # contacts contribute to neither their own nor their partners' scaling.
    excluded = falses(n)
    for b in blacklist_bins
        1 <= b <= n || throw(ArgumentError("blacklist bin $b out of range 1:$n"))
        excluded[b] = true
    end
    w = ones(Float64, n)
    if weights !== nothing
        length(weights) == n || throw(DimensionMismatch("weights length must equal the number of bins"))
        all(x -> isfinite(x) && x >= 0, weights) || throw(ArgumentError("weights must be finite and non-negative"))
        w .= Float64.(weights)
    end
    excluded .|= (w .== 0.0)
    M[excluded, :] .= 0.0
    M[:, excluded] .= 0.0

    for _ in 1:max_iter
        row_sums = vec(sum(M, dims=2))
        valid = (row_sums .> 0) .& .!excluded
        scale = zeros(Float64, n)
        scale[valid] .= 1.0 ./ sqrt.(row_sums[valid] ./ max.(w[valid], eps(Float64)))
        M = Diagonal(scale) * M * Diagonal(scale)
        bias .*= scale

        if any(valid) && maximum(abs.(row_sums[valid] ./ mean(row_sums[valid]) .- 1.0)) < Float64(tol)
            break
        end
    end
    ignore_diag && (M[diagind(M)] .= 0.0)
    bias[excluded] .= 0.0

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
# Knight-Ruiz-family matrix balancing: find the diagonal scaling D such that
# every row of D A D sums to 1, i.e. solve the nonlinear system A x = 1 ./ x
# (x = bias). This is solved with the Osborne (1960) fixed-point iteration
# x_i <- x_i / sqrt(x_i (A x)_i), the log-domain counterpart of the KR
# balancing problem, which converges monotonically for irreducible
# non-negative matrices (the published KR accelerates the same system with
# CG; this variant trades speed for guaranteed convergence).
function _kr_balance!(A::AbstractMatrix{Float64}; max_iter::Int=200, tol::Real=1e-6)
    nv = size(A, 1)
    x = ones(Float64, nv)
    for _ in 1:max_iter
        Ax = A * x
        r = x .* Ax                       # current row sums of D A D
        all(isfinite, r) || break
        maximum(abs.(r .- 1.0)) < Float64(tol) && break
        x ./= sqrt.(max.(r, eps(Float64)))
    end
    # The equation A x = 1 ./ x fixes the overall scale, so no post-hoc
    # renormalisation is applied (renormalising inside the loop shifts the
    # target row sum away from 1).
    return x
end

function normalize_kr(matrix::HiCContactMatrix; max_iter::Int=200, tol::Real=1e-6,
                      blacklist_bins::AbstractVector{<:Integer}=Int[],
                      weights::Union{Nothing,AbstractVector{<:Real}}=nothing)
    M = Matrix(matrix.matrix)
    n = size(M, 1)
    excluded = falses(n)
    for b in blacklist_bins
        1 <= b <= n || throw(ArgumentError("blacklist bin $b out of range 1:$n"))
        excluded[b] = true
    end
    if weights !== nothing
        length(weights) == n || throw(DimensionMismatch("weights length must equal the number of bins"))
        all(x -> isfinite(x) && x >= 0, weights) || throw(ArgumentError("weights must be finite and non-negative"))
        excluded .|= (Float64.(weights) .== 0.0)
        M = M .* reshape(Float64.(weights), :, 1) .* reshape(Float64.(weights), 1, :)
    end
    M[excluded, :] .= 0.0
    M[:, excluded] .= 0.0
    valid = (vec(sum(M, dims=2)) .> 0) .& .!excluded
    bias = zeros(Float64, n)

    if count(valid) > 1
        A = Matrix(M[valid, valid])
        x = _kr_balance!(A; max_iter=max_iter, tol=tol)
        bias[valid] = x
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
function compartment_score(matrix::HiCContactMatrix; n_components::Int=1, genome::Union{Nothing,AbstractDict{<:AbstractString}}=nothing)
    M0 = max.(Matrix(matrix.matrix), 0.0)
    n = size(M0, 1)

    # Observed/expected by genomic distance (cooltools-style) before the
    # correlation step: raw contact maps are dominated by the distance-decay,
    # which swamps the compartment signal (the previous version correlated the
    # raw matrix).
    oe = Matrix{Float64}(M0)
    for d in 1:(n - 1)
        vals = [M0[i, i + d] for i in 1:(n - d)]
        vals = vals[isfinite.(vals)]
        expected = isempty(vals) ? 0.0 : mean(vals)
        expected > 0 || continue
        for i in 1:(n - d)
            oe[i, i + d] = M0[i, i + d] / expected
            oe[i + d, i] = oe[i, i + d]
        end
    end

    c = cor(oe)
    c[.!isfinite.(c)] .= 0.0

    ev = eigen(Symmetric(c))
    pc1 = ev.vectors[:, argmax(ev.values)]

    gc = Float64[]
    for iv in matrix.bins
        push!(gc, _gc_content_matrix_region(iv, matrix; genome=genome))
    end

    # Sign phasing: the eigenvector sign is arbitrary, so when a genome is
    # provided the axis is oriented such that the high-GC (gene-rich) bins are
    # the A compartment — the standard convention.
    if all(isfinite, gc) && any(gc .> 0)
        cor_gc = cor(pc1, gc)
        if isfinite(cor_gc) && cor_gc < 0
            pc1 = -pc1
        end
    end
    correlation_with_gc = all(isfinite, gc) ? abs(cor(pc1, gc)) : NaN

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

# GC content of a bin computed from a real sequence. The genome argument is a
# Dict of chromosome name => sequence (BioSequence or AbstractString); without
# it the function refuses to run instead of returning the previous fabricated
# constant 0.41.
function _gc_content_matrix_region(iv::GenomicInterval, matrix::HiCContactMatrix; genome::Union{Nothing,AbstractDict{<:AbstractString}}=nothing)
    genome === nothing && return NaN
    seq = get(genome, iv.chrom, nothing)
    seq === nothing && return NaN
    lo = clamp(iv.left, 1, length(seq))
    hi = clamp(iv.right, 1, length(seq))
    hi < lo && return NaN
    sub = seq[lo:hi]
    total = length(sub)
    total == 0 && return NaN
    gc = 0
    for b in Vector{UInt8}(codeunits(sub))
        u = b < 0x61 ? b : (b <= 0x7a ? b - 0x20 : b)
        gc += (u == UInt8('G') || u == UInt8('C')) ? 1 : 0
    end
    return gc / total
end

"""
    ab_compartments(matrix)

Assign A/B compartment labels based on PC1 sign.
"""
function ab_compartments(matrix::HiCContactMatrix; genome::Union{Nothing,AbstractDict{<:AbstractString}}=nothing)
    scores = compartment_score(matrix; genome=genome)
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

# Negative-binomial (NB2) GLM with log link fitted by IRLS. The IRLS weights
# enter the weighted normal equations (X'WX) beta = X'Wz, the NB2 dispersion is
# re-estimated from Pearson residuals, and the returned mu is the converged
# mean (the previous version computed weights it never used and returned the
# pre-iteration initialisation).
function _nb_glm_fit(counts::AbstractVector{<:Integer}, offsets::AbstractVector{<:Real}, design::AbstractMatrix{<:Real}; max_iter::Int=100, tol::Real=1e-6)
    n = length(counts)
    k = size(design, 2)
    y = Float64.(counts)
    lo = max.(Float64.(offsets), eps(Float64))

    beta = zeros(Float64, k)
    mu = max.(y, 0.1)
    phi = 1.0
    converged = false

    for iter in 1:max_iter
        eta = design * beta .+ log.(lo)
        mu = max.(exp.(eta), eps(Float64))

        # NB2 working weights and working response.
        W = Diagonal(mu ./ (1.0 .+ mu ./ (phi .* lo)))
        z = eta .+ (y .- mu) ./ mu

        XTW = design' * W * design
        beta_new = try
            (XTW + Diagonal(1e-10 .* ones(k))) \ (design' * (W * z))
        catch
            break
        end

        step = maximum(abs.(beta_new .- beta))
        beta = beta_new

        # Pearson-based NB2 dispersion update (damped for stability).
        pearson = (y .- mu) .^ 2 ./ (mu .+ mu .^ 2 ./ (phi .* lo))
        phi_new = clamp(sum(pearson) / max(n - k, 1), 0.01, 100.0)
        phi = 0.5 * phi + 0.5 * phi_new

        if step < Float64(tol)
            converged = true
            break
        end
    end

    eta = design * beta .+ log.(lo)
    mu = max.(exp.(eta), eps(Float64))
    W = Diagonal(mu ./ (1.0 .+ mu ./ (phi .* lo)))
    cov_unscaled_diag = try
        diag(inv(Symmetric(Matrix(design' * W * design))))
    catch
        diag(pinv(Matrix(design' * W * design)))
    end
    se = sqrt.(max.(phi .* cov_unscaled_diag, 0.0))

    return (beta=beta, se=se, dispersion=phi, cov_unscaled_diag=cov_unscaled_diag, mu=mu, converged=converged)
end

"""
    differential_interaction(cond_a, cond_b, coordinates; design=nothing, dispersion=:eb)

Differential interaction testing via negative binomial GLM with replicate
counts in the design. `dispersion=:eb` shrinks per-bin NB2 dispersions toward
the common dispersion with 10 prior degrees of freedom (edgeR-style
moderation); `dispersion=:per_bin` tests with unmoderated per-bin estimates.
"""
# Prior degrees of freedom for the empirical-Bayes dispersion shrinkage
# (edgeR's default prior.df is 10).
const _EB_PRIOR_DF = 10.0

function differential_interaction(
    cond_a::AbstractVector{<:HiCContactMatrix},
    cond_b::AbstractVector{<:HiCContactMatrix},
    coordinates::AbstractVector{Tuple{Int,Int}};
    design::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
    dispersion::Symbol=:eb)

    n_a = length(cond_a)
    n_b = length(cond_b)
    n_coords = length(coordinates)

    if design === nothing
        # Default design: intercept + condition indicator (A = 1, B = 0). The
        # previous single-column design had no intercept, which forced the B
        # mean to exp(0) = 1 and made every bin look "significant".
        design = hcat(ones(Float64, n_a + n_b),
                      vcat(ones(Float64, n_a), zeros(Float64, n_b)))
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

    fits = Vector{NamedTuple}(undef, n_coords)
    counts_all = Vector{Vector{Float64}}(undef, n_coords)
    mean_a = zeros(Float64, n_coords)
    mean_b = zeros(Float64, n_coords)
    log2_fc = zeros(Float64, n_coords)

    for (ci, (bin1, bin2)) in enumerate(coordinates)
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
        counts_all[ci] = counts

        mean_a[ci] = mean(counts[1:n_a])
        mean_b[ci] = mean(counts[n_a+1:end])
        log2_fc[ci] = mean_b[ci] > 0 && mean_a[ci] > 0 ? log2(mean_b[ci] / mean_a[ci]) : 0.0

        if length(counts) >= 4
            fits[ci] = try
                offsets = ones(Float64, length(counts))
                _nb_glm_fit(Int.(round.(counts)), offsets, design)
            catch
                (beta=[0.0], dispersion=1.0, cov_unscaled_diag=[0.0], converged=false)
            end
        else
            fits[ci] = (beta=[0.0], dispersion=NaN, cov_unscaled_diag=[NaN], converged=false)
        end
    end

    # Empirical-Bayes dispersion sharing (edgeR-style): per-bin NB2
    # dispersions are shrunk toward the common (df-weighted mean) dispersion
    # so that low-count bins borrow strength from the genome-wide trend
    # instead of producing unstable Wald tests.
    # Common dispersion: df-weighted geometric mean of the per-bin estimates
    # (edgeR averages on the log scale, where dispersion estimates are
    # approximately symmetric).
    log_phi_sum = 0.0
    df_total = 0.0
    for fit in fits
        isfinite(get(fit, :dispersion, NaN)) || continue
        df_here = max(length(counts_all[1]) - size(design, 2), 1)
        log_phi_sum += df_here * log(clamp(fit.dispersion, 0.01, 100.0))
        df_total += df_here
    end
    phi_common = df_total > 0 ? clamp(exp(log_phi_sum / df_total), 0.01, 100.0) : 1.0

    for (ci, (bin1, bin2)) in enumerate(coordinates)
        counts = counts_all[ci]
        fit = fits[ci]
        if length(counts) >= 4 && size(design, 2) >= 1 && isfinite(get(fit, :dispersion, NaN))
            stat = fit.beta[end]
            df_here = max(length(counts) - size(design, 2), 1)
            phi_bin = fit.dispersion
            phi_used = dispersion === :eb ?
                       exp((df_here * log(phi_bin) + _EB_PRIOR_DF * log(phi_common)) / (df_here + _EB_PRIOR_DF)) :
                       phi_bin
            cov_diag = get(fit, :cov_unscaled_diag, [NaN])[end]
            se = sqrt(max(phi_used * cov_diag, 0.0))
            pval = se > 0 ? 2 * ccdf(Normal(), abs(stat / se)) : 1.0
        else
            pval = 1.0
        end

        chrom1 = n_a > 0 ? cond_a[1].chrom : "unknown"
        chrom2 = n_a > 0 ? cond_a[1].chrom : "unknown"

        push!(results, (bin1, bin2, chrom1, chrom2, mean_a[ci], mean_b[ci], log2_fc[ci], pval, 1.0))
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

        p = _twocount_pvalue(Int(c_a), Int(c_b), total_a, total_b)
        push!(pvalue, p)
    end

    results = DataFrames.DataFrame(
        count_a=Int.(obs_a),
        count_b=Int.(obs_b),
        log2_fc=log2_fc,
        pvalue=pvalue,
        significance=[_significance_stars(p) for p in pvalue])

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "hic_diff_test";
            parameters=(n_a=n_a, n_b=n_b, total_a=total_a, total_b=total_b))
    end
    return results
end

# Exact two-sided conditional test for a pair of counts: given the total,
# the split is Binomial(n, 0.5) under equal enrichment (depth-normalised
# totals are compared by the caller via the z path). The previous
# "pvalue_mann_whitney" applied a U-statistic to two scalars, which is not a
# valid test; the same name now performs the correct conditional test.
function pvalue_mann_whitney(a::Int, b::Int)
    return _twocount_pvalue(a, b, a, b)
end

function _twocount_pvalue(c_a::Int, c_b::Int, total_a::Int, total_b::Int)
    n = c_a + c_b
    n == 0 && return 1.0
    if total_a > 0 && total_b > 0 && (total_a != total_b)
        # Depth-adjusted conditional binomial: expected split ratio r = tb/ta.
        r = total_b / total_a
        p0 = r / (1.0 + r)
        k = max(c_a, c_b)
        tail = cdf(Binomial(n, p0), min(c_a, c_b)) + ccdf(Binomial(n, p0), k - 1)
        return clamp(min(tail, 1.0), 0.0, 1.0)
    end
    # Equal-depth exact two-sided binomial mid-p style test.
    k = max(c_a, c_b)
    p = min(1.0, 2.0 * ccdf(Binomial(n, 0.5), k - 1))
    return clamp(p, 0.0, 1.0)
end

_significance_stars(p::Real) = p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : p < 0.1 ? "." : ""

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
    bin_size::Int=1,
    min_expected::Real=1.0)

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
            # A floor on the local expected contact level keeps the ratio
            # score finite and comparable across positions (both for
            # candidates and for the Monte-Carlo null below).
            if expected >= Float64(min_expected)
                corner_scores[i, j] = (M[i, j] - expected) / expected
            end
        end
    end

    loop_candidates = Tuple{Int,Int,Float64}[]

    # Empirical null for the corner score: sample matched-distance off-diagonal
    # positions and score them with the same corner statistic (HiCCUPS-style
    # Monte-Carlo FDR instead of the previous `threshold = fdr * 10`).
    rng = MersenneTwister(2026)
    null_scores = Float64[]
    n_null = min(20000, n * n)
    for _ in 1:n_null
        i0 = rand(rng, 1:n)
        d0 = rand(rng, min_dist:min(max_dist, n - 1))
        j0 = i0 + d0
        j0 <= n || continue
        s0 = corner_scores[i0, j0]
        isfinite(s0) && push!(null_scores, s0)
    end

    # BH q-values from the empirical null.
    function _mc_qvalue(score::Float64)
        isempty(null_scores) && return 1.0
        tail = count(>=(score), null_scores) / length(null_scores)
        return clamp(tail, 0.0, 1.0)
    end

    # The empirical tail IS the estimated FDR at that score threshold, so
    # candidates are selected directly on it (applying BH on top would double
    # the correction).
    for i in 1:n
        for j in (i + min_dist):min(n, i + max_dist)
            score = corner_scores[i, j]
            isfinite(score) || continue
            score <= 0 && continue
            q = _mc_qvalue(score)
            q <= Float64(fdr) && push!(loop_candidates, (i, j, score))
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

    has_bins = hasproperty(loops_a, :anchor1_bin) && hasproperty(loops_b, :anchor1_bin)
    merged = DataFrames.DataFrame(
        anchor1_pos=Int[],
        anchor2_pos=Int[],
        found_in=String[],
        max_score=Float64[])
    has_bins && (merged.anchor1_bin = Int[]; merged.anchor2_bin = Int[])

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

        if has_bins
            push!(merged, (Int(row.anchor1_pos), Int(row.anchor2_pos), join(found_in, ","), max_score, Int(row.anchor1_bin), Int(row.anchor2_bin)))
        else
            push!(merged, (Int(row.anchor1_pos), Int(row.anchor2_pos), join(found_in, ","), max_score))
        end
        used_a[i] = true
    end

    for (j, row) in enumerate(eachrow(loops_b))
        if !used_b[j]
            if has_bins
                push!(merged, (Int(row.anchor1_pos), Int(row.anchor2_pos), "B", row.score, Int(row.anchor1_bin), Int(row.anchor2_bin)))
            else
                push!(merged, (Int(row.anchor1_pos), Int(row.anchor2_pos), "B", row.score))
            end
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

    bin_size = max(matrix.bin_size, 1)
    for row in eachrow(loops)
        # Accept loop tables with bin columns (call_loops output) or position
        # columns (merge_loops output); positions are converted with the
        # matrix bin size.
        a1 = hasproperty(loops, :anchor1_bin) ? Int(row.anchor1_bin) :
             (hasproperty(loops, :anchor1_pos) ? Int((row.anchor1_pos - 1) ÷ bin_size + 1) : 0)
        a2 = hasproperty(loops, :anchor2_bin) ? Int(row.anchor2_bin) :
             (hasproperty(loops, :anchor2_pos) ? Int((row.anchor2_pos - 1) ÷ bin_size + 1) : 0)
        (a1 > 0 && a2 > 0) || continue

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

    max_contact = n == 0 || total == 0 ? 0.0 : maximum(M)
    median_contact = count(M .> 0) == 0 ? 0.0 : median(vec(M[M .> 0]))

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

# =============================================================================
# InteractionSet / GInteractions data model (InteractionSet-parity and beyond)
# =============================================================================
# Region-pair representation with assay matrices: the Bioconductor
# InteractionSet/GInteractions data model, plus operations InteractionSet does
# not ship (cooler/mcool round-trip built in, saddle plots, MC-FDR loops,
# provenance on every constructor).

"""
    GInteraction

A single interaction between two genomic anchors (region-pair representation,
as in Bioconductor GInteractions). `anchor1` is always the earlier anchor on
the genome; use `swap_anchors` to flip orientation.
"""
struct GInteraction
    anchor1::GenomicInterval
    anchor2::GenomicInterval
    metadata::Dict{Symbol,Any}
end

function GInteraction(a1::GenomicInterval, a2::GenomicInterval)
    if (a1.chrom, a1.left, a1.right) <= (a2.chrom, a2.left, a2.right)
        return GInteraction(a1, a2, Dict{Symbol,Any}())
    else
        return GInteraction(a2, a1, Dict{Symbol,Any}())
    end
end

"""
    InteractionSet

A vector of `GInteraction`s together with named assay matrices
(`n_interactions × n_samples`) — the Bioconductor InteractionSet model with
cooler-native construction.
"""
struct InteractionSet
    interactions::Vector{GInteraction}
    assays::Dict{String,AbstractMatrix{Float64}}
    sample_names::Vector{String}
    metadata::Dict{String,Any}
end

"""
    interactions_from_matrix(matrix; assay="counts", min_value=0) → InteractionSet

Convert a per-chromosome `HiCContactMatrix` into an InteractionSet with one
GInteraction per non-zero contact pair and the counts as the first assay.
"""
function interactions_from_matrix(matrix::HiCContactMatrix; assay::String="counts", min_value::Real=0.0,
                                  prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    interactions = GInteraction[]
    values = Float64[]
    nz = findnz(matrix.matrix)
    for (i, j, v) in zip(nz[1], nz[2], nz[3])
        v <= min_value && continue
        i <= j || continue
        iv1 = matrix.bins[i]
        iv2 = matrix.bins[j]
        push!(interactions, GInteraction(iv1, iv2))
        push!(values, Float64(v))
    end
    assays = Dict{String,AbstractMatrix{Float64}}(assay => reshape(values, :, 1))
    result = InteractionSet(interactions, assays, ["default"],
                            Dict{String,Any}("source" => "HiCContactMatrix", "chrom" => matrix.chrom,
                                             "bin_size" => matrix.bin_size))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "interactions_from_matrix";
            parameters=(chrom=matrix.chrom, n_interactions=length(interactions), min_value=Float64(min_value)))
    end
    return result
end

"""
    interactions_from_pairs(df; bin_size, chrom_sizes=nothing) → InteractionSet

Build an InteractionSet from a contact/pairs DataFrame with columns
`chrom1, pos1, chrom2, pos2` (positions are 1-based). Anchors are single-base
intervals at each end of every contact.
"""
function interactions_from_pairs(df::DataFrame; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    hasproperty(df, :chrom1) || throw(ArgumentError("pairs table must contain chrom1/pos1/chrom2/pos2"))
    interactions = GInteraction[]
    for row in eachrow(df)
        a1 = GenomicInterval(String(row.chrom1), Int(row.pos1), Int(row.pos1), '.')
        a2 = GenomicInterval(String(row.chrom2), Int(row.pos2), Int(row.pos2), '.')
        push!(interactions, GInteraction(a1, a2))
    end
    result = InteractionSet(interactions, Dict{String,AbstractMatrix{Float64}}(), String[],
                            Dict{String,Any}("source" => "pairs"))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "interactions_from_pairs"; parameters=(n_interactions=length(interactions),))
    end
    return result
end

Base.length(iset::InteractionSet) = length(iset.interactions)
Base.getindex(iset::InteractionSet, i::Integer) = iset.interactions[i]

"""
    counts(iset; assay="counts") → Matrix

Assay matrix (`n_interactions × n_samples`) for the given assay name.
"""
interaction_counts(iset::InteractionSet; assay::String="counts") = iset.assays[assay]

anchors1(iset::InteractionSet) = [g.anchor1 for g in iset.interactions]
anchors2(iset::InteractionSet) = [g.anchor2 for g in iset.interactions]

"""
    regions(iset) → Vector{GenomicInterval}

Deduplicated, sorted registry of all anchor intervals in an InteractionSet
(InteractionSet `regions()` parity).
"""
function regions(iset::InteractionSet)
    seen = Set{Tuple{String,Int,Int}}()
    out = GenomicInterval[]
    for g in iset.interactions
        for iv in (g.anchor1, g.anchor2)
            key = (iv.chrom, iv.left, iv.right)
            if !in(key, seen)
                push!(seen, key)
                push!(out, iv)
            end
        end
    end
    sort!(out; by=iv -> (iv.chrom, iv.left, iv.right))
    return out
end

"""
    swap_anchors(iset) → InteractionSet

Swap the two anchors of every interaction (InteractionSet `swapAnchors`
parity).
"""
function swap_anchors(iset::InteractionSet; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    flipped = [GInteraction(g.anchor2, g.anchor1, g.metadata) for g in iset.interactions]
    result = InteractionSet(flipped, iset.assays, iset.sample_names, iset.metadata)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "swap_anchors"; parameters=(n_interactions=length(flipped),))
    end
    return result
end

"""
    cis_interactions(iset) / trans_interactions(iset) → InteractionSet

Split an InteractionSet into intra-chromosomal (cis) and inter-chromosomal
(trans) parts.
"""
function cis_interactions(iset::InteractionSet)
    keep = [g.anchor1.chrom == g.anchor2.chrom for g in iset.interactions]
    return InteractionSet(iset.interactions[keep], Dict(k => v[keep, :] for (k, v) in iset.assays),
                          iset.sample_names, iset.metadata)
end

function trans_interactions(iset::InteractionSet)
    keep = [g.anchor1.chrom != g.anchor2.chrom for g in iset.interactions]
    return InteractionSet(iset.interactions[keep], Dict(k => v[keep, :] for (k, v) in iset.assays),
                          iset.sample_names, iset.metadata)
end

"""
    find_interactions(iset, interval; anchor=:any) → Vector{Int}

Indices of interactions whose anchors overlap `interval`. `anchor` selects
`:any`, `:both`, `:first`, or `:second`.
"""
# Overlap predicate with maxgap/minoverlap semantics (GenomicRanges-style):
# two 1-based inclusive intervals overlap when the overlap span is at least
# `minoverlap`; a gap of up to `maxgap` bases between them still counts.
function _interval_overlaps(a::GenomicInterval, b::GenomicInterval; maxgap::Int=0, minoverlap::Int=1)
    a.chrom == b.chrom || return false
    overlap = min(a.right, b.right) - max(a.left, b.left) + 1
    overlap >= minoverlap && return true
    gap = max(a.left, b.left) - min(a.right, b.right) - 1
    return gap >= 0 && gap <= maxgap
end

function find_interactions(iset::InteractionSet, interval::GenomicInterval; anchor::Symbol=:any,
                           maxgap::Int=0, minoverlap::Int=1)
    anchor in (:any, :both, :first, :second) || throw(ArgumentError("anchor must be :any, :both, :first or :second"))
    (maxgap >= 0 && minoverlap >= 1) || throw(ArgumentError("maxgap must be >= 0 and minoverlap >= 1"))
    hits = Int[]
    for (i, g) in enumerate(iset.interactions)
        o1 = _interval_overlaps(g.anchor1, interval; maxgap=maxgap, minoverlap=minoverlap)
        o2 = _interval_overlaps(g.anchor2, interval; maxgap=maxgap, minoverlap=minoverlap)
        hit = anchor === :any ? (o1 || o2) :
              anchor === :both ? (o1 && o2) :
              anchor === :first ? o1 : o2
        hit && push!(hits, i)
    end
    return hits
end

# =============================================================================
# Matrix operations: coarsening (zooming)
# =============================================================================

"""
    coarsen(matrix, factor) → HiCContactMatrix

Zoom a contact map out by an integer `factor`: blocks of `factor × factor`
bins are summed (the cooler re-binning semantics). The result has
`bin_size × factor` bins.
"""
function coarsen(matrix::HiCContactMatrix, factor::Integer;
                 prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    factor >= 1 || throw(ArgumentError("factor must be >= 1"))
    factor == 1 && return matrix
    n = size(matrix.matrix, 1)
    nb = cld(n, factor)
    M = Matrix(matrix.matrix)
    out = zeros(Float64, nb, nb)
    for bi in 1:nb
        ilo = (bi - 1) * factor + 1
        ihi = min(bi * factor, n)
        for bj in bi:nb
            jlo = (bj - 1) * factor + 1
            jhi = min(bj * factor, n)
            s = sum(M[ilo:ihi, jlo:jhi])
            out[bi, bj] = s
            out[bj, bi] = s
        end
    end
    new_bins = GenomicInterval[]
    for b in 1:nb
        lo = (b - 1) * factor
        hi = min(b * factor, n)
        left = isempty(matrix.bins) ? lo * matrix.bin_size + 1 : matrix.bins[lo + 1].left
        right = isempty(matrix.bins) ? hi * matrix.bin_size : matrix.bins[hi].right
        push!(new_bins, GenomicInterval(matrix.chrom, left, right, '.'))
    end
    result = HiCContactMatrix(sparse(out), matrix.chrom, new_bins, matrix.bin_size * factor,
                              matrix.normalized, matrix.normalization_method,
                              merge(matrix.metadata, Dict{String,Any}("coarsened_from" => n, "factor" => factor)))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "coarsen";
            parameters=(chrom=matrix.chrom, factor=factor, n_bins_from=n, n_bins_to=nb))
    end
    return result
end

# =============================================================================
# Multi-resolution writing (.mcool)
# =============================================================================

function _write_cool_group(parent, groupname::String, matrix::HiCContactMatrix)
    nb = size(matrix.matrix, 1)
    chrom_size = isempty(matrix.bins) ? nb * matrix.bin_size : maximum(iv.right for iv in matrix.bins)
    chrom_name = matrix.chrom

    nz = findnz(matrix.matrix)
    keep_px = nz[1] .<= nz[2]
    bin1_ids = Int.(nz[1][keep_px]) .- 1
    bin2_ids = Int.(nz[2][keep_px]) .- 1
    counts = Float64.(nz[3][keep_px])

    bin1_offset = zeros(Int, nb + 1)
    for b in bin1_ids
        bin1_offset[min(b + 2, nb + 1)] += 1
    end
    cumsum!(bin1_offset, bin1_offset)

    g = isempty(groupname) ? parent : HDF5.create_group(parent, groupname)
    HDF5.attrs(g)["format"] = "HDF5::Cooler"
    HDF5.attrs(g)["format-version"] = "0.1.1"
    HDF5.attrs(g)["bin-size"] = matrix.bin_size
    HDF5.attrs(g)["nbins"] = nb
    HDF5.attrs(g)["nchroms"] = 1
    HDF5.attrs(g)["npixels"] = length(counts)
    HDF5.attrs(g)["genome-assembly"] = "unknown"

    gc = HDF5.create_group(g, "chroms")
    write(gc, "names", String[chrom_name])
    write(gc, "lengths", Int[chrom_size])

    gb = HDF5.create_group(g, "bins")
    write(gb, "chrom", String[chrom_name for _ in 1:nb])
    write(gb, "start", Int[(i - 1) * matrix.bin_size for i in 1:nb])
    write(gb, "end", Int[min(i * matrix.bin_size, chrom_size) for i in 1:nb])

    gp = HDF5.create_group(g, "pixels")
    write(gp, "bin1_id", bin1_ids)
    write(gp, "bin2_id", bin2_ids)
    write(gp, "count", counts)

    gi = HDF5.create_group(g, "indexes")
    write(gi, "chrom_offset", Int[0, nb])
    write(gi, "bin1_offset", bin1_offset)
    return g
end

"""
    multi_resolution_map(matrix; factors=[1, 2, 4, 5, 10]) → Dict{Int, HiCContactMatrix}

In-memory multi-resolution contact maps: one coarsened matrix per factor,
keyed by resolution (bin_size × factor) — the in-memory counterpart of
`write_mcool`.
"""
function multi_resolution_map(matrix::HiCContactMatrix; factors::AbstractVector{<:Integer}=[1, 2, 4, 5, 10],
                              prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    out = Dict{Int,HiCContactMatrix}()
    for f in factors
        f >= 1 || throw(ArgumentError("factors must be >= 1"))
        m = coarsen(matrix, Int(f); _ctx=nothing)
        out[matrix.bin_size * Int(f)] = m
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "multi_resolution_map";
            parameters=(chrom=matrix.chrom, resolutions=Tuple(keys(out))))
    end
    return out
end

"""
    write_mcool(matrix, path; factors=[1, 2, 4, 5, 10])

Write a multi-resolution `.mcool` file: one cooler group per requested
coarsening factor of `matrix` (resolution = bin_size × factor), read back
with `read_mcool(path; resolution=...)`.
"""
function write_mcool(matrix::HiCContactMatrix, path::AbstractString; factors::AbstractVector{<:Integer}=[1, 2, 4, 5, 10],
                     prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    all(f -> f >= 1, factors) || throw(ArgumentError("factors must be >= 1"))
    HDF5.h5open(path, "w") do h5
        HDF5.attrs(h5)["format"] = "HDF5::MCOOL"
        for f in factors
            res_matrix = coarsen(matrix, Int(f); _ctx=nothing)
            groupname = "resolutions/$(matrix.bin_size * Int(f))"   # cooler convention: groups named by resolution
            _write_cool_group(h5, groupname, res_matrix)
        end
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "write_mcool";
            parameters=(path=String(path), factors=Tuple(Int.(factors)), bin_size=matrix.bin_size))
    end
    return path
end

# =============================================================================
# 4DN .pairs format
# =============================================================================

"""
    read_pairs(path) → DataFrame

Read a 4DN `.pairs` file (optionally gzip-compressed): header lines starting
with `#` are parsed for chrom sizes and the column layout, data lines are
returned as a DataFrame with columns `read_id, chrom1, pos1, chrom2, pos2,
strand1, strand2` (missing columns become `missing` values).
"""
function read_pairs(path::AbstractString;
                    prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    isfile(path) || throw(ArgumentError(".pairs file not found: $path"))
    raw = _is_gzip_file(path) ? readlines(CodecZlib.GzipDecompressorStream(open(path))) : readlines(path)

    chrom_sizes = Dict{String,Int}()
    columns = String[]
    data_rows = Vector{Vector{String}}()   # NB: Vector{Vector{String}}[] would be a triple-nested literal
    for line in raw
        isempty(strip(line)) && continue
        if startswith(line, '#')
            if (m = match(r"^#chromsize:\s+(\S+)\s+(\d+)", line)) !== nothing
                chrom_sizes[String(m[1])] = parse(Int, m[2])
            elseif (m = match(r"^#columns:\s*(.+)$", line)) !== nothing
                columns = String.(split(strip(m[1])))
            end
            continue
        end
        push!(data_rows, String.(split(strip(line))))
    end
    isempty(columns) && (columns = ["read_id", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"])

    result = DataFrame()
    for (ci, col) in enumerate(columns)
        vals = Vector{Union{String,Missing}}(missing, length(data_rows))
        for (ri, fields) in enumerate(data_rows)
            ci <= length(fields) && (vals[ri] = fields[ci])
        end
        result[!, Symbol(col)] = vals
    end

    # Coerce the standard positional columns to proper types.
    for c in (:pos1, :pos2)
        hasproperty(result, c) && (result[!, c] = [v === missing ? missing : parse(Int, v) for v in result[!, c]])
    end

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "read_pairs";
            parameters=(path=String(path), n_rows=nrow(result), n_chroms=length(chrom_sizes)))
    end
    return result
end

_is_gzip_file(path::AbstractString) = let bytes = read(path, 2); length(bytes) >= 2 && bytes[1] == 0x1f && bytes[2] == 0x8b; end

"""
    write_pairs(path, df; chrom_sizes=Dict{String,Int}(), format_version="1.0.0")

Write a DataFrame with columns `chrom1, pos1, chrom2, pos2` (plus optional
`read_id`, `strand1`, `strand2`) as a 4DN `.pairs` file.
"""
function write_pairs(path::AbstractString, df::DataFrame; chrom_sizes::Dict{String,Int}=Dict{String,Int}(),
                     format_version::String="1.0.0",
                     prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    hasproperty(df, :chrom1) || throw(ArgumentError("pairs table must contain chrom1/pos1/chrom2/pos2"))
    open(path, "w") do io
        println(io, "## pairs format $(format_version)")
        for (c, len) in sort(collect(chrom_sizes); by=x -> x[1])
            println(io, "#chromsize: $(c) $(len)")
        end
        cols = ["read_id", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]
        println(io, "#columns: " * join(cols, ' '))
        n = nrow(df)
        for i in 1:n
            read_id = hasproperty(df, :read_id) ? String(df.read_id[i]) : "."
            s1 = hasproperty(df, :strand1) ? String(df.strand1[i]) : "+"
            s2 = hasproperty(df, :strand2) ? String(df.strand2[i]) : "+"
            println(io, join([read_id, String(df.chrom1[i]), Int(df.pos1[i]), String(df.chrom2[i]), Int(df.pos2[i]), s1, s2], '\t'))
        end
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "write_pairs"; parameters=(path=String(path), n_rows=nrow(df)))
    end
    return path
end

# =============================================================================
# HiC-Pro format (genome-wide matrix + abs bed)
# =============================================================================

"""
    read_hicpro(matrix_path, abs_path) → Dict{String, HiCContactMatrix}

Read a HiC-Pro contact map: `<prefix>_matrix.matrix` (bin1 bin2 count, 1-based
global bin ids) plus `<prefix>_abs.bed` (bin_id chrom start end; 0-based
start). Returns one `HiCContactMatrix` per chromosome.
"""
function read_hicpro(matrix_path::AbstractString, abs_path::AbstractString;
                     prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    isfile(matrix_path) || throw(ArgumentError("HiC-Pro matrix not found: $matrix_path"))
    isfile(abs_path) || throw(ArgumentError("HiC-Pro abs bed not found: $abs_path"))

    bin_ids = Int[]
    bin_chroms = String[]
    bin_lefts = Int[]
    bin_rights = Int[]
    for line in eachline(abs_path)
        isempty(strip(line)) && continue
        f = split(strip(line))
        length(f) >= 4 || continue
        push!(bin_ids, parse(Int, f[1]))
        push!(bin_chroms, String(f[2]))
        push!(bin_lefts, parse(Int, f[3]) + 1)   # abs.bed start is 0-based
        push!(bin_rights, parse(Int, f[4]))
    end
    id_index = Dict(b => i for (i, b) in enumerate(bin_ids))

    per_chrom = Dict{String,Dict{Tuple{Int,Int},Float64}}()
    for line in eachline(matrix_path)
        isempty(strip(line)) && continue
        f = split(strip(line))
        length(f) >= 3 || continue
        b1 = get(id_index, parse(Int, f[1]), 0)
        b2 = get(id_index, parse(Int, f[2]), 0)
        (b1 > 0 && b2 > 0) || continue
        c = bin_chroms[b1]
        c == bin_chroms[b2] || continue          # cis maps only
        d = get!(per_chrom, c, Dict{Tuple{Int,Int},Float64}())
        d[(b1, b2)] = parse(Float64, f[3])
    end

    result = Dict{String,HiCContactMatrix}()
    for (chrom, entries) in per_chrom
        keep = [i for (i, c) in enumerate(bin_chroms) if c == chrom]
        remap = Dict(b => i for (i, b) in enumerate(keep))
        nb = length(keep)
        M = spzeros(Float64, nb, nb)
        for ((b1, b2), v) in entries
            M[remap[b1], remap[b2]] = v
        end
        intervals = GenomicInterval[GenomicInterval(chrom, bin_lefts[i], bin_rights[i], '.') for i in keep]
        result[chrom] = HiCContactMatrix(M, chrom, intervals, 0, false, "none",
                                         Dict{String,Any}("source" => String(matrix_path), "format" => "hicpro"))
    end

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "read_hicpro";
            parameters=(matrix=String(matrix_path), n_chroms=length(result)))
    end
    return result
end

"""
    write_hicpro(matrix, prefix)

Write a HiC-Pro pair of files: `prefix_abs.bed` and `prefix_matrix.matrix`
(1-based global bin ids, upper-triangle pixels).
"""
function write_hicpro(matrix::HiCContactMatrix, prefix::AbstractString;
                      prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    nb = size(matrix.matrix, 1)
    open("$(prefix)_abs.bed", "w") do io
        for (i, iv) in enumerate(matrix.bins)
            println(io, "$(i)\t$(iv.chrom)\t$(iv.left - 1)\t$(iv.right)")
        end
    end
    open("$(prefix)_matrix.matrix", "w") do io
        # HiC-Pro matrices list both triangle entries for every contact.
        nz = findnz(matrix.matrix)
        for (i, j, v) in zip(nz[1], nz[2], nz[3])
            println(io, "$(i)\t$(j)\t$(v)")
            i != j && println(io, "$(j)\t$(i)\t$(v)")
        end
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "write_hicpro"; parameters=(prefix=String(prefix), n_bins=nb))
    end
    return ("$(prefix)_abs.bed", "$(prefix)_matrix.matrix")
end

# =============================================================================
# Saddle plots and compartment strength (cooltools parity)
# =============================================================================

"""
    saddle_plot(matrix, pc1; n_bins=25) → NamedTuple

Compute a compartment saddle: the average distance-normalised contact
frequency (O/E) between bins stratified by `pc1` quantile. Returns
`(saddle=K×K matrix, strength, quantile_edges, bin_assignments)`.
`strength = ⟨AA + BB⟩ − ⟨AB + BA⟩`, where ⟨AA⟩ is the mean O/E in the
high-high corner blocks — the cooltools compartment-strength statistic.
"""
function saddle_plot(matrix::HiCContactMatrix, pc1::AbstractVector{<:Real}; n_bins::Int=25, min_dist_bins::Int=5,
                     prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    n = size(matrix.matrix, 1)
    length(pc1) == n || throw(DimensionMismatch("pc1 must have one value per bin"))
    M0 = max.(Matrix(matrix.matrix), 0.0)

    # O/E by genomic distance.
    oe = zeros(Float64, n, n)
    for d in 1:(n - 1)
        vals = Float64[]
        for i in 1:(n - d)
            v = M0[i, i + d]
            isfinite(v) && push!(vals, v)
        end
        expected = isempty(vals) ? 0.0 : mean(vals)
        expected > 0 || continue
        for i in 1:(n - d)
            oe[i, i + d] = M0[i, i + d] / expected
            oe[i + d, i] = oe[i, i + d]
        end
    end
    # Long-range only: near-diagonal pairs have O/E ~ 1 regardless of
    # compartment and would dilute the saddle contrast (cooltools default).
    for d in 1:min(min_dist_bins, n - 1)
        for i in 1:(n - d)
            oe[i, i + d] = 0.0
            oe[i + d, i] = 0.0
        end
    end

    order = sortperm(Float64.(pc1))
    qedges = round.(Int, range(0, n; length=n_bins + 1))
    qbin = zeros(Int, n)
    for (qi, range) in enumerate(1:length(qedges) - 1)
        for r in qedges[range]:(qedges[range + 1])
            r >= 1 && r <= n && (qbin[order[r]] = qi)
        end
    end

    saddle = fill(NaN, n_bins, n_bins)
    for k1 in 1:n_bins, k2 in 1:n_bins
        s = 0.0
        c = 0
        for i in 1:n
            qbin[i] == k1 || continue
            for j in 1:n
                i == j && continue
                qbin[j] == k2 || continue
                v = oe[i, j]
                (isfinite(v) && v > 0) || continue
                s += v
                c += 1
            end
        end
        c > 0 && (saddle[k1, k2] = s / c)
    end

    corner = max(1, n_bins ÷ 4)   # standard: extreme quantile block per side
    aa = _safe_mean(saddle[1:corner, 1:corner])          # low-low (B compartment)
    bb = _safe_mean(saddle[(end - corner + 1):end, (end - corner + 1):end])  # high-high (A)
    ab = (_safe_mean(saddle[1:corner, (end - corner + 1):end]) + _safe_mean(saddle[(end - corner + 1):end, 1:corner])) / 2
    strength = (isnan(aa) || isnan(bb) || isnan(ab)) ? NaN : ((aa + bb) / 2 - ab)

    result = (saddle=saddle, strength=strength, quantile_edges=qedges, bin_assignments=qbin)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "saddle_plot";
            parameters=(chrom=matrix.chrom, n_bins=n_bins, strength=strength))
    end
    return result
end

# Mean over finite values only (NaN-safe; avoids extending Statistics.mean,
# which would shadow the imported function inside this module).
function _safe_mean(values)
    finite = Float64[]
    for v in values
        v isa Real && isfinite(v) && push!(finite, Float64(v))
    end
    isempty(finite) && return NaN
    return sum(finite) / length(finite)
end

"""
    compartment_strength(matrix, pc1; n_bins=25)

Convenience wrapper returning just the saddle compartment strength.
"""
function compartment_strength(matrix::HiCContactMatrix, pc1::AbstractVector{<:Real}; n_bins::Int=25)
    return saddle_plot(matrix, pc1; n_bins=n_bins).strength
end

# =============================================================================
# TAD calling (TopDom-style)
# =============================================================================

"""
    call_tads(matrix; window_bins=10, min_size_bins=3) → DataFrame

TopDom-style TAD caller: compute the mean contact level in a sliding
`(2w+1) × (2w+1)` window, take local minima of the resulting profile as TAD
boundaries, and emit TADs as the intervals between consecutive boundaries
with a strength score (mean within-TAD contact over mean flanking contact).
"""
function call_tads(matrix::HiCContactMatrix; window_bins::Int=10, min_size_bins::Int=3,
                   prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    M = max.(Matrix(matrix.matrix), 0.0)
    n = size(M, 1)
    w = max(window_bins, 1)

    # TopDom gap statistic: at a TAD boundary, contacts between the two
    # flanking windows are depleted relative to the windows' internal
    # contacts. gap(i) = <up,up> + <down,down> - 2·<up,down>; boundaries are
    # local maxima of the gap.
    profile = zeros(Float64, n)
    for i in 2:(n - 1)
        ulo = max(1, i - w); uhi = i - 1
        dlo = i + 1; dhi = min(n, i + w)
        (uhi >= ulo && dhi >= dlo) || continue
        up = M[ulo:uhi, ulo:uhi]
        down = M[dlo:dhi, dlo:dhi]
        cross = M[ulo:uhi, dlo:dhi]
        mean_up = up[up .> 0] |> (x -> isempty(x) ? 0.0 : mean(x))
        mean_down = down[down .> 0] |> (x -> isempty(x) ? 0.0 : mean(x))
        mean_cross = cross[cross .> 0] |> (x -> isempty(x) ? 0.0 : mean(x))
        profile[i] = mean_up + mean_down - 2.0 * mean_cross
    end

    # Local maxima of the gap = TAD boundaries.
    boundaries = Int[1]
    for i in 2:(n - 1)
        profile[i] > 0 && profile[i] >= profile[i - 1] && profile[i] >= profile[i + 1] && push!(boundaries, i)
    end
    push!(boundaries, n)
    unique!(sort!(boundaries))

    rows = NamedTuple[]
    for k in 1:(length(boundaries) - 1)
        lo = boundaries[k]
        hi = boundaries[k + 1]
        (hi - lo + 1) >= min_size_bins || continue
        within_vals = M[lo:hi, lo:hi][M[lo:hi, lo:hi] .> 0]
        within = isempty(within_vals) ? 0.0 : mean(within_vals)
        flo = max(1, lo - w)
        fhi = min(n, hi + w)
        flank_mask = trues(n, n)
        flank_mask[lo:hi, lo:hi] .= false
        flank_vals = M[flo:fhi, flo:fhi][flank_mask[flo:fhi, flo:fhi]]
        flanking = isempty(flank_vals) ? 0.0 : mean(flank_vals)
        strength = flanking > 0 ? within / flanking : NaN
        left = isempty(matrix.bins) ? (lo - 1) * matrix.bin_size + 1 : matrix.bins[lo].left
        right = isempty(matrix.bins) ? hi * matrix.bin_size : matrix.bins[hi].right
        push!(rows, (tad_index=length(rows) + 1, start_bin=lo, end_bin=hi,
                     chrom=matrix.chrom, left=left, right=right,
                     strength=strength, mean_contact=within))
    end

    result = DataFrame(rows)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "call_tads";
            parameters=(chrom=matrix.chrom, window_bins=window_bins, n_tads=nrow(result)))
    end
    return result
end

# =============================================================================
# Virtual 4C and feature aggregation
# =============================================================================

"""
    virtual_4c(matrix, viewpoint_bin; min_dist_bins=2) → DataFrame

One-dimensional contact profile between `viewpoint_bin` and every other bin
(cis restriction window `min_dist_bins` around the viewpoint is excluded),
with genomic positions — virtual 4C.
"""
function virtual_4c(matrix::HiCContactMatrix, viewpoint_bin::Integer; min_dist_bins::Int=2,
                    prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    n = size(matrix.matrix, 1)
    1 <= viewpoint_bin <= n || throw(ArgumentError("viewpoint_bin out of range"))
    rows = NamedTuple[]
    for b in 1:n
        abs(b - viewpoint_bin) < min_dist_bins && continue
        v = matrix.matrix[viewpoint_bin, b]
        pos = isempty(matrix.bins) ? (b - 1) * matrix.bin_size + 1 : matrix.bins[b].left
        push!(rows, (bin=b, position=Int(pos), contact=Float64(v), distance_bins=abs(b - viewpoint_bin)))
    end
    result = DataFrame(rows)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "virtual_4c";
            parameters=(chrom=matrix.chrom, viewpoint=viewpoint_bin, min_dist_bins=min_dist_bins))
    end
    return result
end

"""
    aggregate_at_features(matrix, features; window_bins=10) → NamedTuple

Average the contact sub-matrix around an arbitrary set of genomic features
(peaks, promoters, …): the generalisation of loop APA to any feature set.
Features are converted to central bins by midpoint.
"""
# Features accept any interval-like objects exposing chrom/left/right
# (GenomicInterval and GenomicIntervalEmpty both qualify).
function aggregate_at_features(matrix::HiCContactMatrix, features::AbstractVector; window_bins::Int=10,
                               prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    n = size(matrix.matrix, 1)
    h = max(window_bins, 1)
    M = Matrix(matrix.matrix)
    pileup = zeros(Float64, 2 * h + 1, 2 * h + 1)
    count = 0

    for f in features
        f.chrom == matrix.chrom || continue
        center = clamp(div(f.left + f.right, 2) ÷ max(matrix.bin_size, 1) + 1, 1, n)
        (center - h >= 1 && center + h <= n) || continue
        pileup .+= M[(center - h):(center + h), (center - h):(center + h)]
        count += 1
    end
    count > 0 && (pileup ./= count)

    result = (matrix=pileup, n_features=count, window_bins=h)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "aggregate_at_features";
            parameters=(chrom=matrix.chrom, n_features=count, window_bins=h))
    end
    return result
end

# =============================================================================
# HiCcompare-style loess normalisation between conditions
# =============================================================================

"""
    hic_loess_normalize(cond_a, cond_b, coordinates) → DataFrame

HiCcompare (Stansfield et al. 2018) style normalisation between two
conditions at matched bin pairs: for each pair, `M = log2(mean_b / mean_a)`
and `A = log2(sqrt(mean_a * mean_b))`; a robust binned local-linear trend is
fitted on the (A, M) plane and subtracted. Returns
`(bin1, bin2, mean_a, mean_b, A, M, M_adj, pvalue, padj)` with a Wald z-test
on the adjusted M per pair — the HiCcompare differential test.
"""
function hic_loess_normalize(
    cond_a::AbstractVector{<:HiCContactMatrix},
    cond_b::AbstractVector{<:HiCContactMatrix},
    coordinates::AbstractVector{Tuple{Int,Int}};
    n_bins::Int=100,
    prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    n_a = length(cond_a)
    table = NamedTuple[]
    for (bin1, bin2) in coordinates
        va = Float64[]
        vb = Float64[]
        for mat in cond_a
            (1 <= bin1 <= size(mat.matrix, 1) && 1 <= bin2 <= size(mat.matrix, 1)) && push!(va, mat.matrix[bin1, bin2])
        end
        for mat in cond_b
            (1 <= bin1 <= size(mat.matrix, 1) && 1 <= bin2 <= size(mat.matrix, 1)) && push!(vb, mat.matrix[bin1, bin2])
        end
        mean_a = isempty(va) ? 0.0 : mean(va)
        mean_b = isempty(vb) ? 0.0 : mean(vb)
        (mean_a > 0 && mean_b > 0) || continue
        A = 0.5 * log2(mean_a * mean_b)
        M = log2(mean_b / mean_a)
        push!(table, (bin1=bin1, bin2=bin2, mean_a=mean_a, mean_b=mean_b, A=A, M=M))
    end
    isempty(table) && return DataFrame(bin1=Int[], bin2=Int[], mean_a=Float64[], mean_b=Float64[],
                                       A=Float64[], M=Float64[], M_adj=Float64[], pvalue=Float64[], padj=Float64[])

    df = DataFrame(table)
    fit_A, fit_M = _robust_binned_trend(df.A, df.M; n_bins=n_bins)
    M_adj = df.M .- [InterpLinear(fit_A, fit_M, a) for a in df.A]
    resid_sd = std(M_adj .- mean(M_adj); corrected=true)
    resid_sd = resid_sd > 0 ? resid_sd : 1.0

    pvals = [2 * ccdf(Normal(), abs(m) / resid_sd) for m in M_adj]
    padj = _bh_adjust(pvals)

    result = DataFrames.DataFrame(
        bin1=df.bin1, bin2=df.bin2, mean_a=df.mean_a, mean_b=df.mean_b,
        A=df.A, M=df.M, M_adj=M_adj, pvalue=pvals, padj=padj)
    sort!(result, :pvalue)

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "hic_loess_normalize";
            parameters=(n_pairs=nrow(result), n_cond_a=n_a, n_cond_b=length(cond_b), resid_sd=resid_sd))
    end
    return result
end

# Robust binned local-linear trend: quantile-bin x, take medians, iterate with
# outlier down-weighting (6·MAD), then linear interpolation between bin centers.
function _robust_binned_trend(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; n_bins::Int=100, n_robust::Int=3)
    ord = sortperm(collect(x))
    xs = x[ord]
    ys = y[ord]
    m = length(xs)
    nb = clamp(n_bins, 1, max(m ÷ 20, 1))
    centers = Float64[]
    medians = Float64[]
    edges = round.(Int, range(1, m + 1; length=nb + 1))
    for b in 1:nb
        lo = edges[b]
        hi = edges[b + 1] - 1
        hi < lo && continue
        push!(centers, median(xs[lo:hi]))
        push!(medians, median(ys[lo:hi]))
    end
    weights = ones(Float64, length(centers))
    for _ in 1:n_robust
        fitted = [InterpLinear(centers, medians, a) for a in centers]
        resid = medians .- fitted
        mad = median(abs.(resid .- median(resid))) * 1.4826
        mad > 0 || break
        weights = [abs(r) > 6 * mad ? 0.0 : 1.0 for r in resid]
        for b in eachindex(centers)
            weights[b] == 1.0 || (medians[b] = fitted[b])   # keep outliers on the curve
        end
    end
    return centers, medians
end

function InterpLinear(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real}, q::Real)
    isempty(xs) && return 0.0
    q <= first(xs) && return first(ys)
    q >= last(xs) && return last(ys)
    idx = searchsortedfirst(xs, q)
    idx <= 1 && return first(ys)
    x0, x1 = xs[idx - 1], xs[idx]
    y0, y1 = ys[idx - 1], ys[idx]
    x1 == x0 && return y0
    t = (q - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)
end

function _bh_adjust(pvals::AbstractVector{<:Real})
    n = length(pvals)
    n == 0 && return Float64[]
    order = sortperm(collect(pvals))
    out = Vector{Float64}(undef, n)
    running = 1.0
    for i in n:-1:1
        running = min(running, pvals[order[i]] * n / i)
        out[order[i]] = min(running, 1.0)
    end
    return out
end


# =============================================================================#
# Binary Juicer .hic format — writer and reader                               #
#                                                                              #
# Layout verified against the reference implementation (aidenlab/straw,       #
# C++/straw.cpp: readHeader, readResolutionsFromHeader, readFooter,           #
# readMatrix, readMatrixZoomData, populateBlockMap, readBlock).               #
# =============================================================================

const _HIC_VERSION = 9
const _HIC_BLOCK_DIM = 1000     # bins per block dimension (Juicer default)

_is_hic_magic(path::AbstractString) = let m = read(path, 3); length(m) == 3 && m == UInt8[UInt8('H'), UInt8('I'), UInt8('C')]; end

function _write_i32_le(io::IO, v::Integer)
    raw = UInt32(v & 0xffffffff)
    write(io, UInt8(raw & 0xff))
    write(io, UInt8((raw >> 8) & 0xff))
    write(io, UInt8((raw >> 16) & 0xff))
    write(io, UInt8((raw >> 24) & 0xff))
end

function _read_i32_le(io::IO)
    b1 = read(io, UInt8); b2 = read(io, UInt8); b3 = read(io, UInt8); b4 = read(io, UInt8)
    return reinterpret(Int32, UInt32(b1) | UInt32(b2) << 8 | UInt32(b3) << 16 | UInt32(b4) << 24)
end

function _write_i16_le(io::IO, v::Integer)
    raw = UInt16(v & 0xffff)
    write(io, UInt8(raw & 0xff))
    write(io, UInt8((raw >> 8) & 0xff))
end

function _read_i16_le(io::IO)
    b1 = read(io, UInt8); b2 = read(io, UInt8)
    return reinterpret(Int16, UInt16(b1) | UInt16(b2) << 8)
end

function _read_i64_le(io::IO)
    raw = UInt64(0)
    for k in 0:7
        raw |= UInt64(read(io, UInt8)) << (8 * k)
    end
    return reinterpret(Int64, raw)
end

function _write_f32_le(io::IO, v::Real)
    write(io, reinterpret(UInt32, Float32(v)))
end

function _read_f32_le(io::IO)
    return reinterpret(Float32, read(io, UInt32))
end

# NUL-terminated strings, as used throughout the .hic format
function _hic_write_cstring(io::IO, str::AbstractString)
    write(io, str)
    write(io, UInt8(0))
end

function _hic_read_cstring(io::IO)
    bytes = UInt8[]
    while !eof(io)
        b = read(io, UInt8)
        b == 0x00 && return String(bytes)
        push!(bytes, b)
    end
    throw(ArgumentError("corrupt .hic: unterminated string at end of file"))
end

"""
    write_hic(matrix_or_experiment, path; genome_id="unknown", block_dim=1000)

Write a binary Juicer `.hic` file (version 9, zlib-compressed blocks) in the
exact layout of the aidenlab straw reference implementation: header with
attributes dictionary and chromosome table, `numBpResolutions` list, observed/
NONE matrix sections per chromosome pair (`c1Idx`, `c2Idx`, one MatrixZoomData
with unit "BP", sum/occupancy/stddev/percent95 statistics, block bin count and
column count, block map), zlib blocks holding
`nRecords, binXOffset, binYOffset, useShort, useShortBinX, useShortBinY, type=1`
records with float counts, and a footer master index with `"c1_c2"` keys.
Accepts a single-chromosome `HiCContactMatrix` or an `HiCExperiment` (cis pairs
are written; trans pairs are omitted — straw treats missing pairs as empty).
Reads back with `read_hic(path; chrom=..., resolution=...)`.
"""
function write_hic(matrix::HiCContactMatrix, path::AbstractString; genome_id::String="unknown",
                   block_dim::Int=_HIC_BLOCK_DIM,
                   prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    chroms = String[matrix.chrom]
    lengths = Int[isempty(matrix.bins) ? size(matrix.matrix, 1) * matrix.bin_size :
                  maximum(iv.right for iv in matrix.bins)]
    return _write_hic_file(Dict(matrix.chrom => matrix), chroms, lengths, path, genome_id, block_dim)
end

function write_hic(experiment::HiCExperiment, path::AbstractString; genome_id::String="unknown",
                   block_dim::Int=_HIC_BLOCK_DIM,
                   prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    chroms = sort(collect(keys(experiment.matrices)))
    lengths = Int[maximum(iv.right for iv in experiment.matrices[c].bins; init=0) for c in chroms]
    return _write_hic_file(experiment.matrices, chroms, lengths, path, genome_id, block_dim)
end

function _write_hic_file(matrices::Dict{String,HiCContactMatrix}, chroms::Vector{String},
                         lengths::Vector{Int}, path::AbstractString, genome_id::String, block_dim::Int)
    isempty(chroms) && throw(ArgumentError("no chromosomes to write"))
    idx = Dict(c => i - 1 for (i, c) in enumerate(chroms))   # 0-based chr indices
    resolutions = sort(unique(Int(m.bin_size) for m in values(matrices)))

    function header_bytes(master_position::Int64)
        hb = IOBuffer()
        _hic_write_cstring(hb, "HIC")                # magic "HIC\0", as Juicer writes
        _write_i32_le(hb, _HIC_VERSION)
        write(hb, master_position)
        _hic_write_cstring(hb, genome_id)
        write(hb, Int64(0))                      # NVI position (v9, i64 per reference)
        write(hb, Int64(0))                      # NVI length (v9, i64)
        _write_i32_le(hb, 0)                     # nattributes (empty dictionary)
        _write_i32_le(hb, length(chroms))
        for (c, l) in zip(chroms, lengths)
            _hic_write_cstring(hb, c)
            write(hb, Int64(l))                  # v9 stores i64 chromosome lengths
        end
        _write_i32_le(hb, length(resolutions))   # numBpResolutions
        for r in resolutions
            _write_i32_le(hb, r)
        end
        _write_i32_le(hb, 0)                     # numStringResolutions
        return take!(hb)
    end

    header = header_bytes(0)
    header_size = Int64(length(header))

    # Body: for every cis pair, the matrix section (c1, c2, nRes, MZD with the
    # block map) followed by the zlib blocks. Block entries hold ABSOLUTE file
    # positions; the block index trails the MatrixZoomData inside the same
    # section, per readMatrix/readMatrixZoomData in the reference.
    body = IOBuffer()
    footer_entries = Tuple{String,Int64,Int}[]   # (key, fpos, sizeinbytes)

    for c1 in chroms
        haskey(matrices, c1) || continue
        key = "$(idx[c1])_$(idx[c1])"
        M = Matrix(matrices[c1].matrix)
        nb = size(M, 1)
        nblocks = cld(nb, block_dim)

        # Build zlib-compressed block payloads: type-1 row-major records with
        # float counts, upper-triangle only, absolute 0-based bin ids.
        block_entries = Tuple{Int,Int64,Int}[]               # (blockNumber, 0, compressedSize)
        payloads = Vector{Vector{UInt8}}()
        for br in 0:(nblocks - 1)
            for bc in br:(nblocks - 1)
                xlo = bc * block_dim + 1
                xhi = min((bc + 1) * block_dim, nb)
                ylo = br * block_dim + 1
                yhi = min((br + 1) * block_dim, nb)
                rows_data = Tuple{Int,Vector{Tuple{Int,Float64}}}[]
                nrec = 0
                for i in ylo:yhi
                    row_recs = Tuple{Int,Float64}[]
                    for j in xlo:xhi
                        v = M[i, j]
                        (j - 1 >= i - 1) && v > 0 || continue
                        push!(row_recs, (j - 1, v))
                    end
                    isempty(row_recs) && continue
                    push!(rows_data, (i - 1, row_recs))
                    nrec += length(row_recs)
                end
                nrec == 0 && continue
                raw_block = IOBuffer()
                _write_i32_le(raw_block, nrec)
                _write_i32_le(raw_block, bc * block_dim)      # binXOffset
                _write_i32_le(raw_block, br * block_dim)      # binYOffset
                write(raw_block, UInt8(1))                    # useShort = false (float counts)
                write(raw_block, UInt8(1))                    # useShortBinX = false (i32)
                write(raw_block, UInt8(1))                    # useShortBinY = false (i32)
                write(raw_block, UInt8(1))                    # type 1: row-major lists
                _write_i32_le(raw_block, length(rows_data))
                for (bin_y, row_recs) in rows_data
                    _write_i32_le(raw_block, bin_y)              # binY (absolute, 0-based)
                    _write_i32_le(raw_block, length(row_recs))   # colCount
                    for (bx_idx, v) in row_recs
                        _write_i32_le(raw_block, bx_idx)
                        _write_f32_le(raw_block, v)
                    end
                end
                compressed = transcode(CodecZlib.ZlibCompressor, take!(raw_block))
                push!(block_entries, (br * nblocks + bc, 0, length(compressed)))
                push!(payloads, compressed)
            end
        end

        # Build the section in memory first, so the block positions written
        # into the block map are the final absolute ones (single pass, no
        # draft copy — the earlier two-pass draft left the body duplicated
        # and the footer pointing at the wrong offset).
        section_start = header_size + position(body)
        total = sum(M)
        section = IOBuffer()
        _write_i32_le(section, idx[c1])
        _write_i32_le(section, idx[c1])
        _write_i32_le(section, 1)                            # one resolution
        _hic_write_cstring(section, "BP")
        _write_i32_le(section, 0)                            # old zoom index
        _write_f32_le(section, Float32(total))
        _write_f32_le(section, Float32(count(M .> 0)))
        _write_f32_le(section, Float32(std(M[M .> 0]; corrected=false)))
        _write_f32_le(section, Float32(total > 0 ? total * 0.05 : 0.0))
        _write_i32_le(section, matrices[c1].bin_size)
        _write_i32_le(section, block_dim)                    # blockBinCount
        _write_i32_le(section, nblocks)                      # blockColumnCount
        _write_i32_le(section, length(block_entries))        # nBlocks
        # block entries: positions computed against the final section layout.
        # Layout: header(12) + cstrlen("BP") + 4 + 16 + 12 + 4
        #         + nBlocks * 16 + sum_k (4 + size_k)
        mzd_head = 12 + 3 + 4 + 16 + 12 + 4
        blocks_region_start = mzd_head + length(block_entries) * 16
        block_positions = Int64[]
        offset = blocks_region_start
        for (_, _, sz) in block_entries
            push!(block_positions, section_start + offset)
            offset += sz   # straw: raw zlib bytes at position; size from the index
        end
        for (k, (bn, _, sz)) in enumerate(block_entries)
            _write_i32_le(section, bn)
            write(section, Int64(block_positions[k]))
            _write_i32_le(section, sz)
        end
        for (k, (_, _, _)) in enumerate(block_entries)
            write(section, payloads[k])
        end
        section_bytes = take!(section)
        size_in_bytes = length(section_bytes)
        write(body, section_bytes)
        push!(footer_entries, (key, section_start, size_in_bytes))
    end

    # Footer (master index): nBytes (i64 for v9), nEntries, entries of
    # (key NUL-terminated, fpos i64, size i32) — per readFooter in the
    # reference. nBytes covers the index entries that follow it.
    index_bytes = IOBuffer()
    _write_i32_le(index_bytes, length(footer_entries))
    for (key, fpos, sz) in footer_entries
        _hic_write_cstring(index_bytes, key)
        write(index_bytes, Int64(fpos))
        _write_i32_le(index_bytes, sz)
    end
    index_data = take!(index_bytes)

    master_pos = Int64(header_size + position(body))
    open(path, "w") do io
        write(io, header_bytes(master_pos))
        write(io, take!(body))
        write(io, Int64(length(index_data)))
        write(io, index_data)
    end

    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "write_hic";
            parameters=(path=String(path), n_chroms=length(chroms), version=_HIC_VERSION))
    end
    return path
end


"""
    read_hic(path; chrom, resolution, start=1, stop=nothing, unit=:BP)

Read a binary Juicer `.hic` file — observed contacts, NONE normalization —
following the aidenlab straw reference implementation (versions 6-9 for the
header; v7+ block encoding with per-block offsets, compression flags and
row-major/dense-lower-triangle record types; zlib blocks; v9 i64 chromosome
lengths). Returns the `HiCContactMatrix` for `chrom` at `resolution` over
`[start, stop]`. Non-BP units, unsupported versions, and unknown block
record layouts throw explicit errors. Generic 5-column TSV contact files
fall back to the TSV reader.
"""
function read_hic(path::AbstractString; chrom::String, resolution::Int=10000, start::Int=1, stop::Union{Nothing,Int}=nothing, unit::Symbol=:BP)
    isfile(path) || throw(ArgumentError(".hic file not found: $path"))
    if !_is_hic_magic(path)
        @warn "read_hic: $path is not a binary .hic file; falling back to the generic 5-column TSV contact reader"
        return _read_hic_tsv(path; chrom=chrom, resolution=resolution, start=start, stop=stop)
    end
    unit in (:BP, :FRAG) || throw(ArgumentError("unit must be :BP or :FRAG"))
    io = open(path, "r")
    try
        magic = _hic_read_cstring(io)
        magic == "HIC" || throw(ArgumentError("not a binary .hic file (bad magic): $path"))
        version = _read_i32_le(io)
        version >= 6 || throw(ArgumentError("unsupported .hic version $version (reader supports 6-9)"))
        master_pos = _read_i64_le(io)
        genome_id = _hic_read_cstring(io)
        nvi_pos, nvi_len = Int64(0), Int64(0)
        if version > 8
            nvi_pos = _read_i64_le(io)   # NVI position
            nvi_len = _read_i64_le(io)   # NVI length
        end
        n_attributes = _read_i32_le(io)
        for _ in 1:n_attributes
            _hic_read_cstring(io)   # attribute key
            _hic_read_cstring(io)   # attribute value
        end
        n_chroms = _read_i32_le(io)
        chrom_index = Dict{String,Int}()
        chrom_lengths = Dict{String,Int}()
        for i in 0:(n_chroms - 1)
            name = _hic_read_cstring(io)
            len = version > 8 ? _read_i64_le(io) : Int64(_read_i32_le(io))
            haskey(chrom_index, name) || (chrom_index[name] = i)
            chrom_lengths[name] = len
        end
        haskey(chrom_index, chrom) || throw(ArgumentError("chromosome $chrom not present in $path"))
        c1_idx = chrom_index[chrom]

        # Footer (master index): nBytes (i64 v9 / i32 v6-8), nEntries, entries
        # of (key "c1_c2" NUL-terminated, fpos i64, size i32). Indices in the
        # keys are 0-based chromosome indices.
        seek(io, master_pos)
        version > 8 ? _read_i64_le(io) : _read_i32_le(io)   # nBytes of index
        n_entries = _read_i32_le(io)
        matrix_pos = nothing
        target_key = "$(c1_idx)_$(c1_idx)"
        for _ in 1:n_entries
            key = _hic_read_cstring(io)
            fpos = _read_i64_le(io)
            _read_i32_le(io)   # size in bytes of that matrix section
            if key == target_key
                matrix_pos = fpos
                break
            end
        end

        chrom_size = chrom_lengths[chrom]
        resolution >= 1 || throw(ArgumentError("resolution must be >= 1"))
        nb_total = cld(chrom_size, resolution)
        start_bin = (start - 1) ÷ resolution + 1
        stop_bin = stop === nothing ? min(nb_total, start_bin + 499) : (stop - 1) ÷ resolution + 1
        stop_bin = clamp(stop_bin, 1, nb_total)
        n_bins = max(stop_bin - start_bin + 1, 0)
        matrix = spzeros(Float64, n_bins, n_bins)

        if matrix_pos !== nothing
            # Matrix section: c1, c2, nRes, then one MatrixZoomData per
            # resolution until the requested one is found.
            seek(io, matrix_pos)
            _read_i32_le(io)   # c1
            _read_i32_le(io)   # c2
            n_res = _read_i32_le(io)
            found_mzd = false
            block_map = Dict{Int,Tuple{Int64,Int}}()
            block_bin_count = _HIC_BLOCK_DIM
            block_column_count = cld(nb_total, _HIC_BLOCK_DIM)
            for _ in 1:n_res
                mz_unit = _hic_read_cstring(io)
                _read_i32_le(io)          # old zoom index
                _read_f32_le(io)          # sum counts
                _read_f32_le(io)          # occupied cell count
                _read_f32_le(io)          # std dev
                _read_f32_le(io)          # percent 95
                mz_binsize = _read_i32_le(io)
                mz_block_bin_count = _read_i32_le(io)
                mz_block_column_count = _read_i32_le(io)
                n_blocks = _read_i32_le(io)
                if String(mz_unit) == String(unit) && mz_binsize == resolution
                    block_bin_count = mz_block_bin_count
                    block_column_count = mz_block_column_count
                    for _ in 1:n_blocks
                        bn = _read_i32_le(io)
                        bpos = _read_i64_le(io)
                        bsize = _read_i32_le(io)
                        block_map[bn] = (bpos, bsize)
                    end
                    found_mzd = true
                    break
                else
                    skip(io, n_blocks * (4 + 8 + 4))
                end
            end
            found_mzd || @warn "read_hic: resolution $resolution not present for $chrom in $path"

            # Needed blocks (upper-triangle grid): block (row, col) with
            # row = bin ÷ blockBinCount. Requested records may live in upper
            # blocks even when the query is lower-triangle (data symmetric).
            needed = Set{Int}()
            for r in (start_bin - 1) ÷ block_bin_count : (stop_bin - 1) ÷ block_bin_count
                for c in (start_bin - 1) ÷ block_bin_count : (stop_bin - 1) ÷ block_bin_count
                    push!(needed, r * block_column_count + c)
                    push!(needed, c * block_column_count + r)
                end
            end
            for bn in needed
                haskey(block_map, bn) || continue
                (bpos, bsize) = block_map[bn]
                seek(io, bpos)
                compressed = read(io, bsize)
                raw = IOBuffer(transcode(CodecZlib.ZlibDecompressor, compressed))
                n_records = _read_i32_le(raw)
                bin_x_offset = _read_i32_le(raw)
                bin_y_offset = _read_i32_le(raw)
                use_short = read(raw, UInt8) == 0x00        # 0 means useShort = true
                use_short_bin_x = version > 8 ? read(raw, UInt8) == 0x00 : true
                use_short_bin_y = version > 8 ? read(raw, UInt8) == 0x00 : true
                rec_type = read(raw, UInt8)
                if rec_type == 0x01
                    # row-major: rows of (binY, colCount, [(binX, counts)])
                    n_rows = use_short_bin_y ? _read_i16_le(raw) : _read_i32_le(raw)
                    for _ in 1:n_rows
                        bin_y = bin_y_offset + (use_short_bin_y ? _read_i16_le(raw) : _read_i32_le(raw))
                        col_count = use_short_bin_x ? _read_i16_le(raw) : _read_i32_le(raw)
                        for _ in 1:col_count
                            bin_x = bin_x_offset + (use_short_bin_x ? _read_i16_le(raw) : _read_i32_le(raw))
                            counts = use_short ? Float64(_read_i16_le(raw)) : Float64(_read_f32_le(raw))
                            _hic_store_record!(matrix, bin_x, bin_y, counts, start_bin, stop_bin)
                        end
                    end
                elseif rec_type == 0x02
                    # dense lower-triangle: nPts, w; point i at (row, col) =
                    # (i ÷ w, i % w); bin1 = xoff + col, bin2 = yoff + row
                    n_pts = _read_i32_le(raw)
                    w = _read_i16_le(raw)
                    w == 0 && continue
                    for i in 0:(n_pts - 1)
                        row = i ÷ w
                        col = i - row * w
                        bin_x = bin_x_offset + col
                        bin_y = bin_y_offset + row
                        if use_short
                            c = _read_i16_le(raw)
                            c != -32768 && _hic_store_record!(matrix, bin_x, bin_y, Float64(c), start_bin, stop_bin)
                        else
                            v = _read_f32_le(raw)
                            isfinite(v) && _hic_store_record!(matrix, bin_x, bin_y, Float64(v), start_bin, stop_bin)
                        end
                    end
                else
                    throw(ArgumentError("corrupt .hic: unknown block record type $rec_type"))
                end
            end
        end

        intervals = GenomicInterval[GenomicInterval(chrom, (b - 1) * resolution + 1,
                                                    min(b * resolution, chrom_size), '.') for b in start_bin:stop_bin]
        result = HiCContactMatrix(matrix, chrom, intervals, resolution, false, "NONE",
                                  Dict{String,Any}("source" => String(path), "format" => "hic",
                                                   "genome_id" => genome_id, "nvi" => (nvi_pos, nvi_len),
                                                   "region" => "$(chrom):$(start)-$(stop === nothing ? "end" : stop)"))
        _ctx = active_provenance_context()
        if _ctx !== nothing
            register_provenance!(_ctx, "read_hic";
                parameters=(path=String(path), chrom=chrom, resolution=resolution, n_bins=n_bins, version=version))
        end
        return result
    finally
        close(io)
    end
end

function _hic_store_record!(matrix::AbstractMatrix{Float64}, bin_x::Integer, bin_y::Integer, counts::Real, start_bin::Integer, stop_bin::Integer)
    for (bx, by) in ((bin_x, bin_y), (bin_y, bin_x))   # symmetric fill
        (start_bin - 1) <= bx < stop_bin || continue
        (start_bin - 1) <= by < stop_bin || continue
        matrix[bx - (start_bin - 1) + 1, by - (start_bin - 1) + 1] = counts
    end
    return matrix
end

end
