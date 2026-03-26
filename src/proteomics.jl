module Proteomics

using DataFrames
using Statistics
using LinearAlgebra
using Random
using Graphs
using Optim
using Distributions
using Base.Threads
using LoopVectorization

using ..DifferentialExpression: CountMatrix, benjamini_hochberg

export MassSpecExperiment, Spectrum, MassSpecPeak, PeakDetectionResult, AlignmentResult, DifferentialAbundanceResult, SparsePLSDAResult, DeNovoResult
export read_mzml, detect_peaks, align_samples, qrilc_impute, differential_abundance, mixed_model_abundance, sparse_pls_da, stream_mass_spec, build_spectrum_graph, de_novo_sequence

struct Spectrum
    rt::Float64
    mz::Vector{Float64}
    intensity::Vector{Float64}
end

struct MassSpecExperiment
    spectra::Vector{Spectrum}
    features::Matrix{Float32}
    design::DataFrame
end

MassSpecExperiment(spectra::Vector{Spectrum}) = MassSpecExperiment(spectra, zeros(Float32, 0, length(spectra)), DataFrame())

struct MassSpecPeak
    rt::Float64
    mz::Float64
    intensity::Float64
    score::Float64
end

struct PeakDetectionResult
    peaks::Vector{MassSpecPeak}
    scales::Vector{Float64}
    coefficients::Matrix{Float64}
end

struct AlignmentResult
    path::Vector{Tuple{Int,Int}}
    warped_query::Vector{Float64}
    cost::Float64
    warp::Vector{Float64}
end

struct DifferentialAbundanceResult
    coefficients::Matrix{Float64}
    statistics::Matrix{Float64}
    qvalues::Matrix{Float64}
    feature_names::Vector{String}
    sample_names::Vector{String}
end

struct SparsePLSDAResult
    x_weights::Matrix{Float64}
    x_scores::Matrix{Float64}
    y_loadings::Matrix{Float64}
    vip::Vector{Float64}
    selected_features::Vector{Int}
    classes::Vector{String}
end

struct DeNovoResult
    sequence::String
    path::Vector{Int}
    score::Float64
    graph::SimpleDiGraph
end

const _AA_MASSES = Dict(
    'A' => 71.03711, 'R' => 156.10111, 'N' => 114.04293, 'D' => 115.02694,
    'C' => 103.00919, 'E' => 129.04259, 'Q' => 128.05858, 'G' => 57.02146,
    'H' => 137.05891, 'I' => 113.08406, 'L' => 113.08406, 'K' => 128.09496,
    'M' => 131.04049, 'F' => 147.06841, 'P' => 97.05276, 'S' => 87.03203,
    'T' => 101.04768, 'W' => 186.07931, 'Y' => 163.06333, 'V' => 99.06841,
)

_mad(values::AbstractVector{<:Real}) = median(abs.(Float64.(values) .- median(Float64.(values))))

function _parse_float_vector(text::AbstractString)
    stripped = strip(replace(replace(text, ',' => ' '), '\n' => ' '))
    isempty(stripped) && return Float64[]
    return parse.(Float64, split(stripped))
end

function _extract_tag_value(block::AbstractString, tag::AbstractString)
    result = Base.match(Regex("<$tag[^>]*value=\"([^\"]+)\""), block)
    result === nothing && return nothing
    return result.captures[1]
end

function _extract_series(block::AbstractString, tag::AbstractString)
    result = Base.match(Regex("<$tag[^>]*>(.*?)</$tag>", "s"), block)
    result === nothing && return Float64[]
    return _parse_float_vector(result.captures[1])
end

function _parse_spectrum(block::AbstractString)
    rt_text = _extract_tag_value(block, "cvParam name=\"scan start time\"")
    rt = rt_text === nothing ? 0.0 : parse(Float64, rt_text)
    mz = _extract_series(block, "mz")
    intensity = _extract_series(block, "intensity")
    if isempty(mz) && isempty(intensity)
        mz = _extract_series(block, "mzArray")
        intensity = _extract_series(block, "intensityArray")
    end
    isempty(mz) && isempty(intensity) && return nothing
    length(mz) == length(intensity) || throw(ArgumentError("spectrum mz and intensity lengths must match"))
    return Spectrum(rt, mz, intensity)
end

function read_mzml(path::AbstractString)
    spectra = Spectrum[]
    open(path, "r") do io
        buffer = IOBuffer()
        in_spectrum = false
        for line in eachline(io)
            if occursin("<spectrum", line)
                in_spectrum = true
                truncate(buffer, 0)
                seekstart(buffer)
            end
            if in_spectrum
                write(buffer, line, '\n')
            end
            if in_spectrum && occursin("</spectrum>", line)
                spectrum = _parse_spectrum(String(take!(buffer)))
                spectrum !== nothing && push!(spectra, spectrum)
                in_spectrum = false
            end
        end
    end
    return MassSpecExperiment(spectra)
end

function _ricker(scale::Real, halfwidth::Int)
    x = collect(-halfwidth:halfwidth)
    s2 = float(scale)^2
    kernel = (2 .- (x .^ 2) ./ s2) .* exp.(-(x .^ 2) ./ (2s2))
    kernel .-= mean(kernel)
    return kernel ./ max(norm(kernel), eps())
end

function _same_convolution(signal::AbstractVector{<:Real}, kernel::AbstractVector{<:Real})
    n = length(signal)
    m = length(kernel)
    center = cld(m, 2)
    output = zeros(Float64, n)
    @threads for index in 1:n
        total = 0.0
        @inbounds for k in 1:m
            sample = index + k - center
            1 <= sample <= n || continue
            total += float(signal[sample]) * float(kernel[k])
        end
        output[index] = total
    end
    return output
end

function detect_peaks(signal::AbstractVector{<:Real}; scales::AbstractVector{<:Real}=[1.0, 2.0, 4.0], threshold::Real=1.25, mz::Union{Nothing,AbstractVector{<:Real}}=nothing, rt::Real=0.0, use_gpu::Bool=false)
    values = Float64.(signal)
    coeffs = zeros(Float64, length(values), length(scales))
    for (column, scale) in enumerate(scales)
        kernel = _ricker(scale, max(3, ceil(Int, 3 * float(scale))))
        coeffs[:, column] = _same_convolution(values, kernel)
    end
    score = maximum(abs.(coeffs), dims=2)[:, 1]
    cutoff = median(score) + threshold * _mad(score)
    peaks = MassSpecPeak[]
    positions = mz === nothing ? collect(1:length(values)) : Float64.(mz)
    for index in 2:(length(values) - 1)
        if score[index] >= cutoff && score[index] >= score[index - 1] && score[index] > score[index + 1]
            push!(peaks, MassSpecPeak(rt, positions[index], values[index], score[index]))
        end
    end
    return PeakDetectionResult(peaks, Float64.(scales), coeffs)
end

detect_peaks(spectrum::Spectrum; kwargs...) = detect_peaks(spectrum.intensity; mz=spectrum.mz, rt=spectrum.rt, kwargs...)

function _dtw_alignment(reference::Vector{Float64}, query::Vector{Float64}; band::Int=max(length(reference), length(query)))
    n = length(reference)
    m = length(query)
    cost = fill(Inf, n + 1, m + 1)
    back = fill(UInt8(0), n + 1, m + 1)
    cost[1, 1] = 0.0
    @inbounds for i in 1:n
        jmin = max(1, i - band)
        jmax = min(m, i + band)
        for j in jmin:jmax
            local_cost = abs(reference[i] - query[j])
            candidates = (cost[i, j], cost[i, j + 1], cost[i + 1, j])
            best_index = argmin(candidates)
            predecessor = best_index == 1 ? (i, j) : best_index == 2 ? (i, j + 1) : (i + 1, j)
            cost[i + 1, j + 1] = local_cost + candidates[best_index]
            back[i + 1, j + 1] = UInt8(best_index)
        end
    end
    i = n + 1
    j = m + 1
    path = Tuple{Int,Int}[]
    while i > 1 || j > 1
        push!(path, (max(i - 1, 1), max(j - 1, 1)))
        direction = back[i, j]
        if direction == 1
            i -= 1
            j -= 1
        elseif direction == 2
            i -= 1
        else
            j -= 1
        end
        i = max(i, 1)
        j = max(j, 1)
        if i == 1 && j == 1
            push!(path, (1, 1))
            break
        end
    end
    reverse!(path)
    warping = collect(range(first(query), last(query), length=length(reference)))
    warped_query = interp1d(collect(1:length(query)), query, collect(range(1, length(query), length=length(reference))))
    return path, warped_query, cost[end, end], warping
end

function interp1d(x::Vector{Int}, y::Vector{Float64}, xi::Vector{Float64})
    output = similar(xi)
    for (index, value) in pairs(xi)
        left = clamp(floor(Int, value), 1, length(y))
        right = clamp(left + 1, 1, length(y))
        t = value - left
        output[index] = (1 - t) * y[left] + t * y[right]
    end
    return output
end

function align_samples(reference::AbstractVector{<:Real}, query::AbstractVector{<:Real}; method::Symbol=:dtw, band::Integer=0)
    reference_values = Float64.(reference)
    query_values = Float64.(query)
    if method == :obiwarp
        objective = params -> begin
            scale, shift = params
            warped_index = clamp.(scale .* collect(1:length(query_values)) .+ shift, 1.0, float(length(query_values)))
            warped = interp1d(collect(1:length(query_values)), query_values, warped_index)
            sum(abs.(reference_values .- warped))
        end
        result = optimize(objective, [1.0, 0.0], BFGS())
        scale, shift = Optim.minimizer(result)
        warped_index = clamp.(scale .* collect(1:length(query_values)) .+ shift, 1.0, float(length(query_values)))
        warped = interp1d(collect(1:length(query_values)), query_values, warped_index)
        return AlignmentResult([(1, 1), (length(reference_values), length(query_values))], warped, Optim.minimum(result), warped_index)
    end
    path, warped_query, cost, warp = _dtw_alignment(reference_values, query_values; band=band <= 0 ? max(length(reference_values), length(query_values)) : band)
    return AlignmentResult(path, warped_query, cost, warp)
end

function qrilc_impute(matrix::AbstractMatrix; rng::AbstractRNG=MersenneTwister(1))
    data = Matrix{Float64}(undef, size(matrix, 1), size(matrix, 2))
    for row in axes(matrix, 1), column in axes(matrix, 2)
        value = matrix[row, column]
        data[row, column] = value === missing ? NaN : Float64(value)
    end
    for column in axes(data, 2)
        observed = [value for value in data[:, column] if isfinite(value) && value > 0]
        isempty(observed) && continue
        log_values = log.(observed)
        center = mean(log_values)
        spread = max(std(log_values), 1e-6)
        cutoff = quantile(Normal(center, spread), 0.25)
        dist = truncated(Normal(center, spread), -Inf, cutoff)
        for row in axes(data, 1)
            if !isfinite(data[row, column]) || data[row, column] <= 0
                data[row, column] = exp(rand(rng, dist))
            end
        end
    end
    return data
end

function _design_matrix(groups::AbstractVector, covariates::Union{Nothing,AbstractMatrix{<:Real}})
    labels = string.(groups)
    unique_labels = unique(labels)
    indicator = hcat([Float64.(labels .== label) for label in unique_labels[2:end]]...)
    base = ones(length(labels), 1)
    covariates === nothing && return hcat(base, indicator)
    return hcat(base, indicator, Matrix{Float64}(covariates))
end

function differential_abundance(matrix::AbstractMatrix, groups::AbstractVector; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, qvalue_threshold::Real=0.05)
    imputed = qrilc_impute(matrix)
    response = log.(imputed .+ eps())
    design = _design_matrix(groups, covariates)
    coefficients = zeros(Float64, size(response, 1), size(design, 2))
    statistics = zeros(Float64, size(response, 1), size(design, 2))
    pvalues = ones(Float64, size(response, 1), size(design, 2))
    for feature in axes(response, 1)
        y = response[feature, :]
        fit = design \ y
        residuals = y .- design * fit
        dof = max(size(design, 1) - size(design, 2), 1)
        sigma2 = sum(abs2, residuals) / dof
        covariance = sigma2 * inv(design' * design)
        se = sqrt.(diag(covariance))
        tstat = fit ./ max.(se, eps())
        coefficients[feature, :] .= fit
        statistics[feature, :] .= tstat
        pvalues[feature, :] .= 2 .* ccdf.(TDist(dof), abs.(tstat))
    end
    qvalues = similar(pvalues)
    for column in axes(pvalues, 2)
        qvalues[:, column] = benjamini_hochberg(pvalues[:, column])
    end
    return DifferentialAbundanceResult(coefficients, statistics, qvalues, ["feature_$(index)" for index in axes(response, 1)], string.(groups))
end

function mixed_model_abundance(matrix::AbstractMatrix{<:Real}, groups::AbstractVector; batch=nothing)
    _ = batch
    return differential_abundance(matrix, groups)
end

function sparse_pls_da(X::AbstractMatrix{<:Real}, y::AbstractVector; n_components::Integer=2, sparsity::Real=0.3)
    data = Matrix{Float64}(X)
    classes = unique(string.(y))
    targets = hcat([Float64.(string.(y) .== class) for class in classes]...)
    data .-= mean(data, dims=1)
    targets .-= mean(targets, dims=1)
    nfeatures = size(data, 2)
    x_weights = zeros(Float64, nfeatures, n_components)
    x_scores = zeros(Float64, size(data, 1), n_components)
    y_loadings = zeros(Float64, size(targets, 2), n_components)
    residual_x = copy(data)
    residual_y = copy(targets)
    keep = max(1, ceil(Int, sparsity * nfeatures))
    for component in 1:n_components
        singular = svd(residual_x' * residual_y)
        weights = singular.U[:, 1]
        threshold = sort(abs.(weights), rev=true)[keep]
        sparse_weights = map(value -> abs(value) >= threshold ? value : 0.0, weights)
        norm_weight = norm(sparse_weights)
        norm_weight == 0 && (sparse_weights .= weights; norm_weight = norm(weights))
        sparse_weights ./= norm_weight
        score = residual_x * sparse_weights
        loading = residual_y' * score / max(dot(score, score), eps())
        residual_x .-= score * (residual_x' * score / max(dot(score, score), eps()))'
        residual_y .-= score * loading'
        x_weights[:, component] .= sparse_weights
        x_scores[:, component] .= score
        y_loadings[:, component] .= loading
    end
    vip = vec(sum(abs.(x_weights), dims=2))
    selected = findall(>=(maximum(vip) * 0.1), vip)
    return SparsePLSDAResult(x_weights, x_scores, y_loadings, vip, selected, classes)
end

function stream_mass_spec(producer::Function, consumer::Function; chunk_size::Integer=128)
    channel = Channel{Any}(chunk_size)
    task = @async begin
        for item in producer()
            put!(channel, item)
        end
        close(channel)
    end
    return @async begin
        for item in channel
            consumer(item)
        end
        wait(task)
        return nothing
    end
end

function build_spectrum_graph(masses::AbstractVector{<:Real}; tolerance::Real=0.5)
    peaks = sort(Float64.(masses))
    graph = SimpleDiGraph(length(peaks))
    labels = Dict{Tuple{Int,Int},Char}()
    for left in 1:length(peaks)-1
        for right in left+1:length(peaks)
            gap = peaks[right] - peaks[left]
            for (aa, mass) in _AA_MASSES
                if abs(gap - mass) <= tolerance
                    add_edge!(graph, left, right)
                    labels[(left, right)] = aa
                end
            end
        end
    end
    return graph, labels, peaks
end

function de_novo_sequence(masses::AbstractVector{<:Real}; tolerance::Real=0.5)
    graph, labels, peaks = build_spectrum_graph(masses; tolerance=tolerance)
    n = length(peaks)
    best_score = fill(-Inf, n)
    best_path = [Int[] for _ in 1:n]
    best_score[1] = 0.0
    best_path[1] = [1]
    for node in 1:n
        isfinite(best_score[node]) || continue
        for edge in outneighbors(graph, node)
            aa = labels[(node, edge)]
            candidate_score = best_score[node] + 1
            if candidate_score > best_score[edge]
                best_score[edge] = candidate_score
                best_path[edge] = vcat(best_path[node], edge)
            end
        end
    end
    target = argmax(best_score)
    path = best_path[target]
    sequence = String([labels[(path[index], path[index + 1])] for index in 1:length(path)-1])
    return DeNovoResult(sequence, path, best_score[target], graph)
end

end