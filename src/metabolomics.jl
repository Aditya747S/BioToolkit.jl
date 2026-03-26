module Metabolomics

using DataFrames
using Statistics
using LinearAlgebra
using Random

using ..Proteomics: MassSpecExperiment, Spectrum, qrilc_impute, differential_abundance, detect_peaks, align_samples
using ..Proteomics: _mad

export MassSpecExperiment, Spectrum, MetabolomicsSourceTrackingResult, metabolomics_source_tracking_model, metabolomics_source_tracking, metabolomics_source_tracking_posterior_summary
export metabolomics_differential_abundance, metabolomics_streaming_analysis, annotate_metabolite_features, quantify_metabolite_variation

struct MetabolomicsSourceTrackingResult
    chain
    mean_proportions::Vector{Float64}
    median_proportions::Vector{Float64}
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
end

function metabolomics_source_tracking_posterior_summary(chain)
    samples = Matrix{Float64}(Array(chain))
    mean_proportions = vec(mean(samples, dims=1))
    median_proportions = vec(median(samples, dims=1))
    lower_bounds = [quantile(view(samples, :, column), 0.025) for column in axes(samples, 2)]
    upper_bounds = [quantile(view(samples, :, column), 0.975) for column in axes(samples, 2)]
    return MetabolomicsSourceTrackingResult(chain, mean_proportions, median_proportions, Float64.(lower_bounds), Float64.(upper_bounds))
end

const _METABOLOMICS_TURING_LOADED = Ref(false)
const _METABOLOMICS_TURING_MODULE = Ref{Any}(nothing)
const _METABOLOMICS_MODEL_IMPL = Ref{Any}(nothing)

function _metabolomics_turing_module()
    turing = _METABOLOMICS_TURING_MODULE[]
    turing === nothing && throw(ArgumentError("Turing is required for metabolomics_source_tracking"))
    return turing
end

function metabolomics_source_tracking_model(observed::AbstractVector{<:Integer}, source_profiles::AbstractMatrix{<:Real})
    model_impl = _METABOLOMICS_MODEL_IMPL[]
    model_impl === nothing && throw(ArgumentError("Turing is required for metabolomics_source_tracking"))
    return Base.invokelatest(model_impl, observed, source_profiles)
end

function metabolomics_source_tracking(observed::AbstractVector{<:Integer}, source_profiles::AbstractMatrix{<:Real}; draws::Integer=250, rng::AbstractRNG=MersenneTwister(1))
    draws > 0 || throw(ArgumentError("draws must be positive"))
    turing = _metabolomics_turing_module()
    model_impl = _METABOLOMICS_MODEL_IMPL[]
    model_impl === nothing && throw(ArgumentError("Turing is required for metabolomics_source_tracking"))
    model = Base.invokelatest(model_impl, observed, source_profiles)
    chain = Base.invokelatest(turing.sample, rng, model, turing.NUTS(), draws)
    return metabolomics_source_tracking_posterior_summary(chain)
end

function metabolomics_differential_abundance(matrix::AbstractMatrix{<:Real}, groups::AbstractVector; covariates::Union{Nothing,AbstractMatrix{<:Real}}=nothing)
    return differential_abundance(matrix, groups; covariates=covariates)
end

function metabolomics_streaming_analysis(source::Function, sink::Function; threshold::Real=1.0)
    task = @async begin
        for experiment in source()
            peaks = map(detect_peaks, experiment.spectra)
            sink(experiment, peaks)
            if any(result -> length(result.peaks) >= threshold, peaks)
                break
            end
        end
    end
    return task
end

function annotate_metabolite_features(features::AbstractMatrix{<:Real}; labels::AbstractVector{<:AbstractString}=String[])
    names = isempty(labels) ? ["feature_$(index)" for index in axes(features, 1)] : String.(labels)
    return DataFrame(feature=names, mean_intensity=vec(mean(features; dims=2)), variance=vec(var(features; dims=2)))
end

function quantify_metabolite_variation(features::AbstractMatrix{<:Real}; robust::Bool=true)
    data = Matrix{Float64}(features)
    summary = robust ? median(data, dims=2) : mean(data, dims=2)
    dispersion = robust ? map(row -> _mad(collect(row)), eachrow(data)) : vec(std(data; dims=2))
    return DataFrame(center=vec(summary), dispersion=Float64.(dispersion))
end

end