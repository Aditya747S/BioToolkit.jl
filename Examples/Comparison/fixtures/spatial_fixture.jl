module SpatialBenchmarkFixtures

using Random
using Statistics
using BioToolkit

export spatial_deconvolution_fixture, align_cell_type_order, evaluate_fraction_recovery

function _normalize_rows(matrix::AbstractMatrix{<:Real})
    output = Matrix{Float64}(matrix)
    for row_index in 1:size(output, 1)
        row_sum = sum(output[row_index, :])
        if row_sum > 0
            output[row_index, :] ./= row_sum
        else
            output[row_index, :] .= 1.0 / size(output, 2)
        end
    end
    return output
end

function _correlation_or_zero(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    if std(a) <= eps(Float64) || std(b) <= eps(Float64)
        return 0.0
    end
    return cor(a, b)
end

function _jsd(p::AbstractVector{<:Real}, q::AbstractVector{<:Real})
    p_safe = max.(Float64.(p), eps(Float64))
    q_safe = max.(Float64.(q), eps(Float64))
    p_safe ./= sum(p_safe)
    q_safe ./= sum(q_safe)
    midpoint = 0.5 .* (p_safe .+ q_safe)
    kl_pm = sum(p_safe .* log2.(p_safe ./ midpoint))
    kl_qm = sum(q_safe .* log2.(q_safe ./ midpoint))
    return sqrt(max(0.0, 0.5 * (kl_pm + kl_qm)))
end

function evaluate_fraction_recovery(predicted_fractions::AbstractMatrix{<:Real}, truth_fractions::AbstractMatrix{<:Real})
    size(predicted_fractions) == size(truth_fractions) || throw(DimensionMismatch("predicted and truth fraction matrices must have identical dimensions"))

    predicted = _normalize_rows(predicted_fractions)
    truth = _normalize_rows(truth_fractions)

    rmse = sqrt(mean((predicted .- truth) .^ 2))
    mae = mean(abs.(predicted .- truth))

    correlations = [_correlation_or_zero(predicted[:, col], truth[:, col]) for col in axes(truth, 2)]
    dominant_accuracy = mean(argmax(predicted[row, :]) == argmax(truth[row, :]) for row in axes(truth, 1))
    mean_jsd = mean(_jsd(predicted[row, :], truth[row, :]) for row in axes(truth, 1))

    return (
        rmse=rmse,
        mae=mae,
        mean_correlation=mean(correlations),
        dominant_accuracy=dominant_accuracy,
        mean_jsd=mean_jsd,
    )
end

function align_cell_type_order(result::DeconvolutionResult, truth_ids::Vector{String})
    mapping = Dict(id => idx for (idx, id) in enumerate(result.cell_type_ids))
    reordered = zeros(Float64, size(result.cell_type_fractions, 1), length(truth_ids))

    for (truth_col, cell_type_id) in enumerate(truth_ids)
        haskey(mapping, cell_type_id) || throw(ArgumentError("missing cell type '$cell_type_id' in deconvolution result"))
        reordered[:, truth_col] = result.cell_type_fractions[:, mapping[cell_type_id]]
    end

    return reordered
end

function spatial_deconvolution_fixture(; seed::Integer=2026)
    Random.seed!(seed)

    n_genes = 90
    n_types = 3
    cells_per_type = 60

    gene_ids = ["gene_$(idx)" for idx in 1:n_genes]
    cell_type_ids = ["TypeA", "TypeB", "TypeC"]

    type_profiles = fill(2.0, n_genes, n_types)
    marker_block = 20

    for type_idx in 1:n_types
        marker_start = (type_idx - 1) * marker_block + 1
        marker_stop = type_idx * marker_block
        type_profiles[marker_start:marker_stop, type_idx] .= 22.0
    end

    type_profiles[61:75, :] .= 6.0

    n_reference_cells = n_types * cells_per_type
    reference_counts = zeros(Int, n_genes, n_reference_cells)
    reference_labels = String[]
    reference_cell_ids = String[]

    cell_cursor = 1
    for type_idx in 1:n_types
        base = type_profiles[:, type_idx]
        base ./= sum(base)
        for local_idx in 1:cells_per_type
            depth = 380 + rand(-40:40)
            expected = depth .* base
            noisy = expected .+ sqrt.(expected) .* randn(n_genes)
            reference_counts[:, cell_cursor] = round.(Int, max.(noisy, 0.0))
            push!(reference_labels, cell_type_ids[type_idx])
            push!(reference_cell_ids, "ref_$(type_idx)_$(local_idx)")
            cell_cursor += 1
        end
    end

    reference = SingleCellExperiment(reference_counts, gene_ids, reference_cell_ids)

    truth_fractions = [
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0;
        0.75 0.25 0.0;
        0.25 0.75 0.0;
        0.0 0.65 0.35;
        0.0 0.35 0.65;
        0.55 0.0 0.45;
        0.45 0.0 0.55;
        0.35 0.35 0.30;
        0.25 0.25 0.50;
        0.50 0.30 0.20;
    ]

    n_spots = size(truth_fractions, 1)
    spot_counts = zeros(Int, n_genes, n_spots)

    for spot_idx in 1:n_spots
        depth = 520 + rand(-50:50)
        mixed_profile = type_profiles * truth_fractions[spot_idx, :]
        mixed_profile ./= sum(mixed_profile)
        expected = depth .* mixed_profile
        noisy = expected .+ sqrt.(expected) .* randn(n_genes)
        spot_counts[:, spot_idx] = round.(Int, max.(noisy, 0.0))
    end

    spot_ids = ["spot_$(idx)" for idx in 1:n_spots]
    theta = range(0.0, 2.0 * pi; length=n_spots + 1)[1:end-1]
    coords = hcat(cos.(theta), sin.(theta))

    spatial_sce = SingleCellExperiment(spot_counts, gene_ids, spot_ids; spatial_coords=coords)
    spatial = SpatialExperiment(spatial_sce)

    return (
        reference=reference,
        reference_labels=reference_labels,
        spatial=spatial,
        truth_fractions=truth_fractions,
        cell_type_ids=cell_type_ids,
    )
end

end
