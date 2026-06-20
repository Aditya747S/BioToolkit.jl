using BioToolkit
using Random
using Statistics
using LinearAlgebra
using SparseArrays
using Plots

Random.seed!(42)

function synthetic_counts(n_genes::Int, n_cells::Int; latent_dim::Int=4)
    gene_loadings = zeros(Float64, n_genes, latent_dim)
    gene_groups = rand(1:latent_dim, n_genes)
    for gene_index in 1:n_genes
        dominant_factor = gene_groups[gene_index]
        gene_loadings[gene_index, dominant_factor] = 2.5 + rand()
        for factor_index in 1:latent_dim
            factor_index == dominant_factor && continue
            gene_loadings[gene_index, factor_index] = 0.1 * randn()
        end
    end

    cell_scores = randn(latent_dim, n_cells)
    signal = gene_loadings * cell_scores / sqrt(latent_dim)
    baseline = 35.0 .+ 4.0 .* randn(n_genes, 1)
    noisy_counts = max.(round.(Int, baseline .+ 14.0 .* signal .+ 0.35 .* randn(n_genes, n_cells)), 0)
    noisy_counts[noisy_counts .< 2] .= 0
    return sparse(noisy_counts)
end

function procrustes_align(source::AbstractMatrix{<:Real}, target::AbstractMatrix{<:Real})
    source_matrix = Matrix{Float64}(source)
    target_matrix = Matrix{Float64}(target)

    source_mean = vec(mean(source_matrix, dims=1))
    target_mean = vec(mean(target_matrix, dims=1))
    source_centered = source_matrix .- reshape(source_mean, 1, :)
    target_centered = target_matrix .- reshape(target_mean, 1, :)

    u, singular_values, vt = svd(transpose(source_centered) * target_centered)
    rotation = u * vt
    scale = sum(singular_values) / max(sum(abs2, source_centered), eps(Float64))
    aligned = scale .* (source_centered * rotation) .+ reshape(target_mean, 1, :)
    return aligned, rotation, scale, source_mean, target_mean
end

function projection_consistency_report(; n_genes::Int=1000, n_cells::Int=500, n_components::Int=30, reference_fraction::Real=0.8, save_plot::Bool=true, plot_path::AbstractString="projection_consistency.png")
    reference_cells = round(Int, reference_fraction * n_cells)
    query_cells = n_cells - reference_cells

    counts = synthetic_counts(n_genes, n_cells)
    gene_ids = ["gene_$(i)" for i in 1:n_genes]
    cell_ids = ["cell_$(i)" for i in 1:n_cells]

    println("Calculating ground-truth PCA on the full synthetic dataset...")
    sce_full = SingleCellExperiment(counts, gene_ids, cell_ids)
    full_normalized = normalize_counts(sce_full)
    full_truth = run_pca(sce_full; normalized=full_normalized, n_components=n_components, use_variable_features=false)

    reference_idx = 1:reference_cells
    query_idx = (reference_cells + 1):n_cells

    sce_reference = SingleCellExperiment(counts[:, reference_idx], gene_ids, cell_ids[reference_idx])
    sce_query = SingleCellExperiment(counts[:, query_idx], gene_ids, cell_ids[query_idx])

    println("Fitting the projection model on the 80 percent reference split...")
    reference_projection = fit_singlecell_projection_model(
        sce_reference;
        n_components=n_components,
        use_variable_features=false,
    )

    projected_reference = project_singlecell(sce_reference, reference_projection; reduction_name="pca_projection_reference")
    projected_query = project_singlecell(sce_query, reference_projection; reduction_name="pca_projection_query")

    true_reference = full_truth[reference_idx, :]
    true_query = full_truth[query_idx, :]

    println("Aligning the reference projection back onto the full-data PCA frame...")
    aligned_reference, rotation, scale, projected_mean, truth_mean = procrustes_align(projected_reference, true_reference)
    aligned_query = scale .* ((Matrix{Float64}(projected_query) .- reshape(projected_mean, 1, :)) * rotation) .+ reshape(truth_mean, 1, :)

    pc1_correlation = cor(aligned_query[:, 1], true_query[:, 1])
    pc2_correlation = cor(aligned_query[:, 2], true_query[:, 2])
    rmse = sqrt(mean((aligned_query .- true_query) .^ 2))

    println("------------------------------------------------")
    println("PROJECTION CONSISTENCY REPORT")
    println("------------------------------------------------")
    println("Reference cells: $(reference_cells)")
    println("Query cells: $(query_cells)")
    println("PC1 Correlation (Projected vs True): $(round(pc1_correlation, digits=4))")
    println("PC2 Correlation (Projected vs True): $(round(pc2_correlation, digits=4))")
    println("Aligned RMSE: $(round(rmse, digits=4))")

    if save_plot
        p = scatter(
            true_query[:, 1],
            true_query[:, 2];
            label="True Query",
            marker=:circle,
            alpha=0.75,
            title="Projection Consistency",
            xlabel="PC1",
            ylabel="PC2",
        )
        scatter!(
            p,
            aligned_query[:, 1],
            aligned_query[:, 2];
            label="Projected Query",
            marker=:cross,
            alpha=0.75,
        )
        savefig(p, plot_path)
        println("Plot saved to $(plot_path)")
        println("------------------------------------------------")
    end

    return (
        pc1_correlation=pc1_correlation,
        pc2_correlation=pc2_correlation,
        rmse=rmse,
        projected_query=aligned_query,
        true_query=true_query,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    report = projection_consistency_report()
    @assert report.pc1_correlation > 0.99 "Projection consistency failed: PC1 correlation too low"
    @assert report.pc2_correlation > 0.99 "Projection consistency failed: PC2 correlation too low"
    println("SUCCESS: projection model is internally consistent on the held-out query set.")
end