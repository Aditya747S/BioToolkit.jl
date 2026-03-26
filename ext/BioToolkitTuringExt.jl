module BioToolkitTuringExt

using BioToolkit
using Turing

const METABOLOMICS_SOURCE_TRACKING_MODEL_IMPL = Turing.@model function _metabolomics_source_tracking_model_impl(observed::AbstractVector{<:Integer}, source_profiles::AbstractMatrix{<:Real})
	nsources = size(source_profiles, 2)
	proportions ~ Dirichlet(fill(1.0, nsources))
	mixture = Matrix{Float64}(source_profiles) * proportions
	mixture = mixture ./ sum(mixture)
	observed ~ Multinomial(sum(observed), mixture)
end

const SOURCE_TRACKING_MODEL_IMPL = Turing.@model function _source_tracking_model_impl(observed::AbstractVector{<:Integer}, source_profiles::AbstractMatrix{<:Real})
	nsources = size(source_profiles, 2)
	proportions ~ Dirichlet(fill(1.0, nsources))
	mixture = Matrix{Float64}(source_profiles) * proportions
	mixture = mixture ./ sum(mixture)
	observed ~ Multinomial(sum(observed), mixture)
end

function __init__()
	BioToolkit.Metabolomics._METABOLOMICS_MODEL_IMPL[] = METABOLOMICS_SOURCE_TRACKING_MODEL_IMPL
	BioToolkit.Metabolomics._METABOLOMICS_TURING_MODULE[] = Turing
	BioToolkit.Metabolomics._METABOLOMICS_TURING_LOADED[] = true

	BioToolkit.Microbiome._SOURCE_TRACKING_MODEL_IMPL[] = SOURCE_TRACKING_MODEL_IMPL
	BioToolkit.Microbiome._SOURCE_TRACKING_TURING_MODULE[] = Turing
	BioToolkit.Microbiome._SOURCE_TRACKING_TURING_LOADED[] = true
end

end
