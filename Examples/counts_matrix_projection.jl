using BioToolkit
using CSV
using DataFrames
using Statistics
using LinearAlgebra
using SparseArrays
using Plots

function load_counts_matrix(csv_path::AbstractString)
    counts_df = CSV.read(csv_path, DataFrame)
    counts_matrix = Matrix{Int}(counts_df[:, 2:end])
    return sparse(counts_matrix)
end

function procrustes_julia(source::AbstractMatrix{<:Real}, target::AbstractMatrix{<:Real})
    source_matrix = Matrix{Float64}(source)
    target_matrix = Matrix{Float64}(target)

    source_mean = vec(mean(source_matrix, dims=1))
    target_mean = vec(mean(target_matrix, dims=1))
    source_centered = source_matrix .- reshape(source_mean, 1, :)
    target_centered = target_matrix .- reshape(target_mean, 1, :)

    factorization = svd(transpose(source_centered) * target_centered)
    rotation = factorization.U * factorization.Vt
    scale = sum(factorization.S) / max(sum(abs2, source_centered), eps(Float64))
    return scale, rotation, source_mean, target_mean
end

function counts_matrix_projection_report(; csv_path::AbstractString=joinpath(@__DIR__, "counts_matrix.csv"), reference_cells::Int=400, n_components::Int=30, save_plot::Bool=true, plot_path::AbstractString=joinpath(@__DIR__, "counts_matrix_projection.png"))
    counts_sparse = load_counts_matrix(csv_path)
    n_genes, n_cells = size(counts_sparse)
    reference_cells < n_cells || throw(ArgumentError("reference_cells must be smaller than the number of cells"))

    gene_ids = ["gene_$(i)" for i in 1:n_genes]
    cell_ids = ["cell_$(i)" for i in 1:n_cells]
    reference_idx = 1:reference_cells
    query_idx = (reference_cells + 1):n_cells

    sce_full = SingleCellExperiment(counts_sparse, gene_ids, cell_ids)
    full_normalized = normalize_counts(sce_full)
    truth_coords = run_pca(sce_full; normalized=full_normalized, n_components=n_components, use_variable_features=false)

    sce_reference = SingleCellExperiment(counts_sparse[:, reference_idx], gene_ids, cell_ids[reference_idx])
    sce_query = SingleCellExperiment(counts_sparse[:, query_idx], gene_ids, cell_ids[query_idx])

    projection_model = fit_singlecell_projection_model(
        sce_reference;
        n_components=n_components,
        use_variable_features=false,
        scale_features=false,
    )
    projected_reference = project_singlecell(sce_reference, projection_model; reduction_name="counts_matrix_projection_reference")
    projected_query = project_singlecell(sce_query, projection_model; reduction_name="counts_matrix_projection_query")

    true_reference = truth_coords[reference_idx, :]
    true_query = truth_coords[query_idx, :]

    scale, rotation, projected_mean, truth_mean = procrustes_julia(projected_reference, true_reference)
    aligned_query = scale .* ((Matrix{Float64}(projected_query) .- reshape(projected_mean, 1, :)) * rotation) .+ reshape(truth_mean, 1, :)

    pc1_corr = cor(aligned_query[:, 1], true_query[:, 1])
    rmse = sqrt(mean((aligned_query .- true_query) .^ 2))

    println("------------------------------------------------")
    println("COUNTS MATRIX PROJECTION REPORT")
    println("------------------------------------------------")
    println("Genes: $(n_genes)")
    println("Cells: $(n_cells)")
    println("Reference cells: $(reference_cells)")
    println("Query cells: $(length(query_idx))")
    println("Julia RMSE (on R data): $(round(rmse, digits=4))")
    println("Julia Cor PC1: $(round(pc1_corr, digits=4))")

    if save_plot
        p = scatter(
            true_query[:, 1],
            true_query[:, 2];
            label="True Query",
            marker=:circle,
            alpha=0.75,
            title="Counts Matrix Projection",
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
        rmse=rmse,
        pc1_corr=pc1_corr,
        projected_query=aligned_query,
        true_query=true_query,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    report = counts_matrix_projection_report()
    @assert report.pc1_corr > 0.99 "Projection consistency failed: PC1 correlation too low"
    @assert report.rmse < 1.0 "Projection consistency failed: RMSE too high"
    println("SUCCESS: counts_matrix.csv projection benchmark passed.")
end