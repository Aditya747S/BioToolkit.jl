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
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, AASeq, AminoAcidAlphabet, BioSequence
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_container_provenance!, register_provenance!

export MassSpecExperiment, Spectrum, MassSpecPeak, PeakDetectionResult, AlignmentResult, DifferentialAbundanceResult, SparsePLSDAResult, DeNovoResult
export read_mzml, detect_peaks, align_samples, qrilc_impute, differential_abundance, mixed_model_abundance, sparse_pls_da, stream_mass_spec, build_spectrum_graph, de_novo_sequence
export dia_like_quantification, phosphosite_localization, glycoproteomics_motif_table, project_peptides_to_structure
export protein_inference_top3, ptm_site_enrichment

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

struct PeakDetectionResult <: AbstractAnalysisResult
    peaks::Vector{MassSpecPeak}
    scales::Vector{Float64}
    coefficients::Matrix{Float64}
    provenance::ResultProvenance
end

PeakDetectionResult(peaks, scales, coefficients) =
    PeakDetectionResult(peaks, scales, coefficients, provenance_record("PeakDetectionResult", "proteomics"))

struct AlignmentResult <: AbstractAnalysisResult
    path::Vector{Tuple{Int,Int}}
    warped_query::Vector{Float64}
    cost::Float64
    warp::Vector{Float64}
    provenance::ResultProvenance
end

AlignmentResult(path, warped_query, cost, warp) =
    AlignmentResult(path, warped_query, cost, warp, provenance_record("AlignmentResult", "proteomics"))

struct DifferentialAbundanceResult <: AbstractAnalysisResult
    coefficients::Matrix{Float64}
    statistics::Matrix{Float64}
    qvalues::Matrix{Float64}
    feature_names::Vector{String}
    sample_names::Vector{String}
    provenance::ResultProvenance
end

DifferentialAbundanceResult(coefficients, statistics, qvalues, feature_names, sample_names) =
    DifferentialAbundanceResult(coefficients, statistics, qvalues, feature_names, sample_names, provenance_record("DifferentialAbundanceResult", "proteomics"))

struct SparsePLSDAResult <: AbstractAnalysisResult
    x_weights::Matrix{Float64}
    x_scores::Matrix{Float64}
    y_loadings::Matrix{Float64}
    vip::Vector{Float64}
    selected_features::Vector{Int}
    classes::Vector{String}
    provenance::ResultProvenance
end

SparsePLSDAResult(x_weights, x_scores, y_loadings, vip, selected_features, classes) =
    SparsePLSDAResult(x_weights, x_scores, y_loadings, vip, selected_features, classes, provenance_record("SparsePLSDAResult", "proteomics"))

struct DeNovoResult <: AbstractAnalysisResult
    sequence::BioSequence{AminoAcidAlphabet}
    path::Vector{Int}
    score::Float64
    graph::SimpleDiGraph
    provenance::ResultProvenance
end

DeNovoResult(sequence, path, score, graph) =
    DeNovoResult(sequence, path, score, graph, provenance_record("DeNovoResult", "proteomics"))

const _AA_MASSES = Dict(
    'A' => 71.03711, 'R' => 156.10111, 'N' => 114.04293, 'D' => 115.02694,
    'C' => 103.00919, 'E' => 129.04259, 'Q' => 128.05858, 'G' => 57.02146,
    'H' => 137.05891, 'I' => 113.08406, 'L' => 113.08406, 'K' => 128.09496,
    'M' => 131.04049, 'F' => 147.06841, 'P' => 97.05276, 'S' => 87.03203,
    'T' => 101.04768, 'W' => 186.07931, 'Y' => 163.06333, 'V' => 99.06841)

@inline function _register_proteomics_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

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

function read_mzml(path::String)
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
    result = MassSpecExperiment(spectra)
    _ctx = active_provenance_context()


    return _register_proteomics_result!(_ctx, result, "read_mzml"; parents=String[], parameters=(path=path, spectrum_count=length(spectra)))
end

function _ricker(scale::Real, halfwidth::Int)
    x = collect(-halfwidth:halfwidth)
    s2 = float(scale)^2
    kernel = (2 .- (x .^ 2) ./ s2) .* exp.(-(x .^ 2) ./ (2s2))
    kernel .-= mean(kernel)
    return kernel ./ max(norm(kernel), eps())
end

function _same_convolution(signal::AbstractVector{<:Real}, kernel::AbstractVector{<:Real}; multi_thread::Bool=true)
    n = length(signal)
    m = length(kernel)
    center = cld(m, 2)
    output = zeros(Float64, n)
    if multi_thread && n > 1 && Threads.nthreads() > 1
        @threads for index in 1:n
            total = 0.0
            @inbounds for k in 1:m
                sample = index + k - center
                1 <= sample <= n || continue
                total += float(signal[sample]) * float(kernel[k])
            end
            output[index] = total
        end
    else
        for index in 1:n
            total = 0.0
            @inbounds for k in 1:m
                sample = index + k - center
                1 <= sample <= n || continue
                total += float(signal[sample]) * float(kernel[k])
            end
            output[index] = total
        end
    end
    return output
end

function detect_peaks(signal::AbstractVector{<:Real}; scales::AbstractVector{<:Real}=[1.0, 2.0, 4.0], threshold::Real=1.25, mz::Union{Nothing,AbstractVector{<:Real}}=nothing, rt::Real=0.0, use_gpu::Bool=false, use_cuda::Bool=use_gpu, multi_thread::Bool=true)
    if use_cuda
        @warn "CUDA peak detection is not enabled in this build; falling back to CPU implementation."
    end
    values = Float64.(signal)
    coeffs = zeros(Float64, length(values), length(scales))
    for (column, scale) in enumerate(scales)
        kernel = _ricker(scale, max(3, ceil(Int, 3 * float(scale))))
        coeffs[:, column] = _same_convolution(values, kernel; multi_thread=multi_thread)
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
    result = PeakDetectionResult(peaks, Float64.(scales), coeffs)
    _ctx = active_provenance_context()


    return _register_proteomics_result!(_ctx, result, "detect_peaks"; parents=String[], parameters=(peak_count=length(peaks), scale_count=length(scales), threshold=float(threshold), use_cuda=use_cuda, multi_thread=multi_thread))
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
        alignment = AlignmentResult([(1, 1), (length(reference_values), length(query_values))], warped, Optim.minimum(result), warped_index)
        _ctx = active_provenance_context()
        return _register_proteomics_result!(_ctx, alignment, "align_samples"; parents=String[], parameters=(method=method, band=band, cost=alignment.cost))
    end
    path, warped_query, cost, warp = _dtw_alignment(reference_values, query_values; band=band <= 0 ? max(length(reference_values), length(query_values)) : band)
    alignment = AlignmentResult(path, warped_query, cost, warp)
    _ctx = active_provenance_context()


    return _register_proteomics_result!(_ctx, alignment, "align_samples"; parents=String[], parameters=(method=method, band=band, cost=alignment.cost))
end

function qrilc_impute(matrix::AbstractMatrix; rng::AbstractRNG=MersenneTwister(1), prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
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
    return _register_proteomics_result!(_ctx, data, "qrilc_impute"; parents=provenance_parent_ids(matrix), parameters=(row_count=size(matrix, 1), column_count=size(matrix, 2)))
end

function _design_matrix(groups::AbstractVector, covariates::Union{Nothing,AbstractMatrix{<:Real}})
    labels = string.(groups)
    unique_labels = unique(labels)
    indicator = hcat([Float64.(labels .== label) for label in unique_labels[2:end]]...)
    base = ones(length(labels), 1)
    covariates === nothing && return hcat(base, indicator)
    return hcat(base, indicator, Matrix{Float64}(covariates))
end

function differential_abundance(matrix::AbstractMatrix, groups::AbstractVector; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing, qvalue_threshold::Real=0.05, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    imputed = qrilc_impute(matrix; _ctx=_ctx)
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
    result = DifferentialAbundanceResult(coefficients, statistics, qvalues, ["feature_$(index)" for index in axes(response, 1)], string.(groups))


    return _register_proteomics_result!(_ctx, result, "differential_abundance"; parents=provenance_parent_ids(matrix), parameters=(feature_count=size(response, 1), sample_count=size(response, 2), qvalue_threshold=float(qvalue_threshold)))
end

function mixed_model_abundance(matrix::AbstractMatrix{<:Real}, groups::AbstractVector; batch=nothing)
    _ = batch
    _ctx = active_provenance_context()


    return differential_abundance(matrix, groups; _ctx=_ctx)
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
    result = SparsePLSDAResult(x_weights, x_scores, y_loadings, vip, selected, classes)
    _ctx = active_provenance_context()


    return _register_proteomics_result!(_ctx, result, "sparse_pls_da"; parents=provenance_parent_ids(X), parameters=(component_count=n_components, sparsity=float(sparsity), feature_count=size(X, 2)))
end

function stream_mass_spec(producer::Function, consumer::Function; chunk_size::Integer=128)
    channel = Channel{Any}(chunk_size)
    task = @async begin
        for item in producer()
            put!(channel, item)
        end
        close(channel)
    end
    task_handle = @async begin
        for item in channel
            consumer(item)
        end
        wait(task)
        return nothing
    end
    _ctx === nothing || register_provenance!(_ctx, "stream_mass_spec"; parents=String[], parameters=(chunk_size=chunk_size))
    return task_handle
end

function build_spectrum_graph(masses::AbstractVector{<:Real}; tolerance::Real=0.5, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
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
    result = (graph, labels, peaks)


    return _register_proteomics_result!(_ctx, result, "build_spectrum_graph"; parents=provenance_parent_ids(masses), parameters=(mass_count=length(masses), tolerance=float(tolerance)))
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
    sequence = AASeq(String([labels[(path[index], path[index + 1])] for index in 1:length(path)-1]))
    result = DeNovoResult(sequence, path, best_score[target], graph)
    _ctx = active_provenance_context()


    return _register_proteomics_result!(_ctx, result, "de_novo_sequence"; parents=provenance_parent_ids(masses), parameters=(tolerance=float(tolerance), path_length=length(path), score=best_score[target]))
end

@inline _prot_bh(p::Vector{Float64}) = benjamini_hochberg(clamp.(p, 0.0, 1.0))

"""
    dia_like_quantification(intensity, groups)

DIA-like differential quantification across two groups.
"""
function dia_like_quantification(intensity::AbstractMatrix{<:Real}, groups::AbstractVector; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    X = log1p.(Float64.(intensity))
    g = Symbol.(groups)
    labels = unique(g)
    length(labels) == 2 || throw(ArgumentError("exactly two groups are currently supported"))
    a = findall(==(labels[1]), g)
    b = findall(==(labels[2]), g)

    lfc = vec(mean(X[:, b], dims=2) .- mean(X[:, a], dims=2))
    p = Float64[]
    for i in 1:size(X, 1)
        xa = @view X[i, a]
        xb = @view X[i, b]
        s2 = var(xa) / max(length(a), 1) + var(xb) / max(length(b), 1)
        if s2 <= eps(Float64)
            push!(p, 1.0)
        else
            z = (mean(xb) - mean(xa)) / sqrt(s2)
            push!(p, 2 * ccdf(Normal(), abs(z)))
        end
    end

    result = DataFrame(feature_id=["feature_$(i)" for i in 1:size(X, 1)], log2_fc=lfc, pvalue=p, padj=_prot_bh(p))


    return _register_proteomics_result!(_ctx, result, "dia_like_quantification"; parents=provenance_parent_ids(intensity), parameters=(feature_count=size(X, 1), sample_count=size(X, 2), group_count=length(labels)))
end

"""
    phosphosite_localization(peptide, site_probabilities)

Localize phosphosite with confidence classes from posterior probabilities.
"""
function phosphosite_localization(peptide::AbstractString, site_probabilities::AbstractDict{<:Integer,<:Real})
    isempty(site_probabilities) && throw(ArgumentError("site_probabilities must not be empty"))
    prob, site = findmax(Dict(Int(k) => Float64(v) for (k, v) in site_probabilities))
    residue = site <= ncodeunits(peptide) ? peptide[site] : '?'
    confidence = prob >= 0.75 ? "class-I" : prob >= 0.5 ? "class-II" : "class-III"
    result = (peptide=String(peptide), site=site, residue=string(residue), probability=prob, confidence=confidence)
    _ctx = active_provenance_context()


    return _register_proteomics_result!(_ctx, result, "phosphosite_localization"; parents=String[], parameters=(site=site, probability=prob, confidence=confidence))
end

phosphosite_localization(peptide::BioSequence{AminoAcidAlphabet}, site_probabilities::AbstractDict{<:Integer,<:Real}) = phosphosite_localization(String(peptide), site_probabilities)

"""
    glycoproteomics_motif_table(peptides)

Detect N-X-S/T glycosylation motifs in peptide sequences.
"""
function glycoproteomics_motif_table(peptides::AbstractVector{<:AbstractString}; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    out = DataFrame(peptide=String[], n_glyco_motifs=Int[], motif_positions=String[])
    for pep in peptides
        p = uppercase(String(pep))
        pos = Int[]
        for i in 1:(ncodeunits(p) - 2)
            if p[i] == 'N' && p[i + 1] != 'P' && (p[i + 2] == 'S' || p[i + 2] == 'T')
                push!(pos, i)
            end
        end
        push!(out, (String(pep), length(pos), isempty(pos) ? "" : join(pos, ",")))
    end


    return _register_proteomics_result!(_ctx, out, "glycoproteomics_motif_table"; parents=provenance_parent_ids(peptides), parameters=(peptide_count=length(peptides), motif_count=sum(out.n_glyco_motifs)))
end

glycoproteomics_motif_table(peptides::AbstractVector{<:BioSequence{AminoAcidAlphabet}}) = glycoproteomics_motif_table(String.(peptides))

"""
    project_peptides_to_structure(peptides, residue_embeddings)

Project peptide sequences to structure-aware embedding space by averaging
residue-level vectors.
"""
function project_peptides_to_structure(peptides::AbstractVector{<:AbstractString}, residue_embeddings::AbstractDict{<:AbstractString,<:AbstractVector{<:Real}})
    dim = isempty(residue_embeddings) ? 0 : length(first(values(residue_embeddings)))
    dim > 0 || throw(ArgumentError("residue_embeddings must contain at least one entry"))

    out = zeros(Float64, length(peptides), dim)
    for (i, pep) in enumerate(peptides)
        tokens = [string(c) for c in collect(uppercase(String(pep)))]
        vecs = [Float64.(residue_embeddings[t]) for t in tokens if haskey(residue_embeddings, t)]
        if !isempty(vecs)
            out[i, :] .= vec(mean(reduce(hcat, vecs), dims=2))
        end
    end
    _ctx = active_provenance_context()


    return _register_proteomics_result!(_ctx, out, "project_peptides_to_structure"; parents=provenance_parent_ids(peptides), parameters=(peptide_count=length(peptides), dimension=dim))
end

project_peptides_to_structure(peptides::AbstractVector{<:BioSequence{AminoAcidAlphabet}}, residue_embeddings::AbstractDict{<:AbstractString,<:AbstractVector{<:Real}}) = project_peptides_to_structure(String.(peptides), residue_embeddings)

"""
    protein_inference_top3(peptide_table; protein_col=:protein, intensity_col=:intensity)

Top3 protein inference from peptide intensities.
"""
function protein_inference_top3(peptide_table::DataFrame; protein_col::Symbol=:protein, intensity_col::Symbol=:intensity, sample_col::Union{Nothing,Symbol}=nothing, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    hasproperty(peptide_table, protein_col) || throw(ArgumentError("missing protein column"))
    hasproperty(peptide_table, intensity_col) || throw(ArgumentError("missing intensity column"))

    tbl = copy(peptide_table)
    tbl[!, intensity_col] = Float64.(tbl[!, intensity_col])

    if sample_col === nothing || !hasproperty(tbl, sample_col)
        out = DataFrame(protein=String[], n_peptides=Int[], top3_intensity=Float64[])
        for g in groupby(tbl, protein_col)
            vals = sort(Float64.(g[!, intensity_col]); rev=true)
            push!(out, (String(g[1, protein_col]), nrow(g), sum(vals[1:min(3, length(vals))])))
        end
        sort!(out, :top3_intensity, rev=true)
        return _register_proteomics_result!(_ctx, out, "protein_inference_top3"; parents=provenance_parent_ids(peptide_table), parameters=(protein_col=protein_col, intensity_col=intensity_col, sample_col="none", protein_count=nrow(out)))
    end

    out = DataFrame(sample=String[], protein=String[], n_peptides=Int[], top3_intensity=Float64[])
    for g in groupby(tbl, [sample_col, protein_col])
        vals = sort(Float64.(g[!, intensity_col]); rev=true)
        push!(out, (String(g[1, sample_col]), String(g[1, protein_col]), nrow(g), sum(vals[1:min(3, length(vals))])))
    end
    sort!(out, [:sample, :top3_intensity], rev=[false, true])
    _ctx = active_provenance_context()


    return _register_proteomics_result!(_ctx, out, "protein_inference_top3"; parents=provenance_parent_ids(peptide_table), parameters=(protein_col=protein_col, intensity_col=intensity_col, sample_col=String(sample_col), protein_count=nrow(out)))
end

"""
    ptm_site_enrichment(modified_sites, background_sites)

Site-wise hypergeometric enrichment for PTM-bearing residues.
"""
function ptm_site_enrichment(modified_sites::AbstractVector{<:AbstractString}, background_sites::AbstractVector{<:AbstractString}; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    mod_counts = Dict{String,Int}()
    bg_counts = Dict{String,Int}()
    for site in modified_sites
        s = String(site)
        mod_counts[s] = get(mod_counts, s, 0) + 1
    end
    for site in background_sites
        s = String(site)
        bg_counts[s] = get(bg_counts, s, 0) + 1
    end

    n = length(modified_sites)
    N = max(length(background_sites), n)
    n > 0 || begin
        result = DataFrame(site=String[], modified_count=Int[], background_count=Int[], enrichment=Float64[], pvalue=Float64[], padj=Float64[])
        return _register_proteomics_result!(_ctx, result, "ptm_site_enrichment"; parents=provenance_parent_ids(modified_sites, background_sites), parameters=(modified_count=0, background_count=length(background_sites)))
    end

    sites = sort!(collect(keys(mod_counts)))
    k_mod = Int[]
    k_bg = Int[]
    enrich = Float64[]
    pvals = Float64[]

    for site in sites
        k = get(mod_counts, site, 0)
        K = get(bg_counts, site, 0)
        p = K > 0 ? ccdf(Hypergeometric(N, K, n), k - 1) : 1.0
        e = (k / max(n, 1)) / max(K / max(N, 1), eps(Float64))
        push!(k_mod, k)
        push!(k_bg, K)
        push!(enrich, e)
        push!(pvals, p)
    end

    out = DataFrame(site=sites, modified_count=k_mod, background_count=k_bg, enrichment=enrich, pvalue=pvals)
    out[!, :padj] = _prot_bh(pvals)
    sort!(out, [:padj, :enrichment], rev=[false, true])


    return _register_proteomics_result!(_ctx, out, "ptm_site_enrichment"; parents=provenance_parent_ids(modified_sites, background_sites), parameters=(modified_count=length(modified_sites), background_count=length(background_sites)))
end

end
