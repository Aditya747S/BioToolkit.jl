using DataFrames

@testset "BioPlotting" begin
    results = [
        BioToolkit.DEResult("gene1", 10.0, 2.5, 0.5, 12.0, 1e-6, 1e-4),
        BioToolkit.DEResult("gene2", 12.0, -2.1, 0.6, 11.0, 2e-4, 0.01),
        BioToolkit.DEResult("gene3", 8.0, 0.2, 0.4, 1.0, 0.8, 0.9),
    ]

    points, summary = BioToolkit.volcano_data(results; lfc_cutoff=1.0, fdr_cutoff=0.05, label_top=1)
    @test length(points) == 3
    @test summary.significant == 2
    @test any(point -> point.labeled, points)

    ma_points, ma_summary = BioToolkit.ma_data(results)
    @test length(ma_points) == 3
    @test ma_summary.significant == 2

    heatmap = BioToolkit.clustered_heatmap([1.0 2.0 3.0; 4.0 5.0 6.0]; row_labels=["r1", "r2"], column_labels=["c1", "c2", "c3"])
    @test size(heatmap.matrix) == (2, 3)
    @test length(heatmap.row_order) == 2
    @test length(heatmap.column_order) == 3
end

@testset "Enrichment" begin
    database = BioToolkit.build_annotation_database([
        BioToolkit.EnrichmentTerm("GO:0001", "Signal transduction", "GO", ["gene1", "gene2"], String[]),
        BioToolkit.EnrichmentTerm("GO:0002", "Cell cycle", "GO", ["gene2", "gene3"], ["GO:0001"]),
        BioToolkit.EnrichmentTerm("KEGG:0001", "Metabolism", "KEGG", ["gene3", "gene4"], String[]),
    ])

    builtin = BioToolkit.builtin_annotation_database()
    @test !isempty(builtin.terms)
    @test any(term -> term.namespace == "GO", values(builtin.terms))
    @test any(term -> term.namespace == "KEGG", values(builtin.terms))
    @test builtin.terms["GO:0006281"].parents == ["GO:0006974"]
    @test builtin.terms["GO:0006302"].parents == ["GO:0006281"]
    @test builtin.terms["GO:0000082"].parents == ["GO:0007049"]
    @test builtin.terms["KEGG:hsa05200"].parents == ["KEGG:hsa00001"]
    @test builtin.terms["KEGG:hsa04110"].parents == ["KEGG:hsa05200"]
    @test builtin.terms["GO:0006302"].genes ⊆ builtin.terms["GO:0006281"].genes
    @test builtin.terms["GO:0006281"].genes ⊆ builtin.terms["GO:0006974"].genes
    @test builtin.terms["GO:0007049"].genes ⊆ builtin.terms["GO:0008150"].genes
    @test builtin.terms["KEGG:hsa04110"].genes ⊆ builtin.terms["KEGG:hsa05200"].genes
    @test builtin.terms["KEGG:hsa05200"].genes ⊆ builtin.terms["KEGG:hsa00001"].genes

    builtin_go = BioToolkit.go_enrichment(["BRCA1", "BRCA2", "RAD51"], builtin)
    @test !isempty(builtin_go)
    @test any(result -> result.term_id == "GO:0006302", builtin_go)
    @test any(result -> result.term_id == "GO:0006281", builtin_go)

    mapped = BioToolkit.map_ids(database.mapper, ["gene1", "gene5"])
    @test mapped == ["gene1", "gene5"]

    results = BioToolkit.enrichment_test(["gene1", "gene2"], database)
    @test !isempty(results)
    @test results[1].term_id in ("GO:0001", "GO:0002")
    @test results[1].padj <= 1.0

    go_results = BioToolkit.go_enrichment(["gene1", "gene2"], database)
    @test all(result -> result.namespace == "GO", go_results)

    kegg_results = BioToolkit.kegg_enrichment(["gene3"], database)
    @test all(result -> result.namespace == "KEGG", kegg_results)

    dotplot_result = BioToolkit.dotplot(results; top_n=2)
    @test length(dotplot_result.results) <= 2
end

@testset "SingleCell" begin
    counts = [10 0 1 2; 12 1 0 3; 0 8 2 0; 1 1 7 9]
    spatial_coords = Float64[0.0 0.0; 1.0 0.0; 2.0 0.0; 3.0 0.0]
    experiment = BioToolkit.SingleCellExperiment(counts, ["gene1", "gene2", "gene3", "gene4"], ["cell1", "cell2", "cell3", "cell4"]; spatial_coords=spatial_coords)
    @test length(experiment.counts.nzval) < length(experiment.gene_ids) * length(experiment.cell_ids)
    @test experiment.spatial_coords == spatial_coords

    spatial_counts = [10 10 1 1; 10 0 10 0; 1 1 1 1; 2 2 2 2]
    spatial_experiment = BioToolkit.SingleCellExperiment(spatial_counts, ["gene1", "gene2", "gene3", "gene4"], ["cell1", "cell2", "cell3", "cell4"]; spatial_coords=spatial_coords)
    spatial_neighbors = BioToolkit.find_spatial_neighbors(spatial_experiment; k=1)
    @test length(spatial_neighbors) == 4
    @test haskey(spatial_experiment.neighbors, "spatial_neighbors")
    spatial_test = BioToolkit.moran_i_test(spatial_experiment, "gene1"; neighbors=spatial_neighbors, normalize=false)
    @test spatial_test isa BioToolkit.SpatialMoranResult
    @test spatial_test.moran_i > 0
    @test BioToolkit.moran_i_test(spatial_experiment, "gene2"; neighbors=spatial_neighbors, normalize=false).moran_i < 0
    spatial_genes = BioToolkit.find_spatially_variable_genes(spatial_experiment; neighbors=spatial_neighbors, normalize=false, top_n=2)
    @test spatial_genes[1].gene == "gene1"

    matrix = BioToolkit.normalize_counts(experiment)
    @test size(matrix) == (4, 4)

    saved = mktempdir() do tempdir
        prefix = joinpath(tempdir, "singlecell")
        BioToolkit.save_singlecell_experiment(prefix, experiment)
        loaded = BioToolkit.load_singlecell_experiment(prefix)
        @test Array(loaded.counts) == Array(experiment.counts)
        @test loaded.gene_ids == experiment.gene_ids
        @test loaded.cell_ids == experiment.cell_ids
        @test loaded.spatial_coords == experiment.spatial_coords
        loaded
    end
    @test size(saved.counts) == size(experiment.counts)

    transformed, genes, cells = BioToolkit.sctransform(experiment; min_total=0)
    @test size(transformed, 2) == 4
    @test length(genes) == size(transformed, 1)
    @test length(cells) == 4

    experiment.metadata["obs"] = DataFrames.DataFrame(_index=copy(experiment.cell_ids), batch=["a", "a", "b", "b"])
    experiment.metadata["var"] = DataFrames.DataFrame(_index=copy(experiment.gene_ids), family=["housekeeping", "marker", "marker", "marker"])
    experiment.metadata["layers"] = Dict("spliced" => counts .+ 1)
    experiment.reductions["pca"] = Float64[1.0 0.0; 0.8 0.2; 0.2 0.8; 0.0 1.0]

    h5ad_path = mktempdir() do tempdir
        path = joinpath(tempdir, "singlecell.h5ad")
        BioToolkit.write_h5ad(experiment, path)
        loaded = BioToolkit.read_h5ad(path)
        @test Array(loaded.counts) == Array(experiment.counts)
        @test loaded.cell_ids == experiment.cell_ids
        @test loaded.gene_ids == experiment.gene_ids
        @test loaded.spatial_coords == experiment.spatial_coords
        @test haskey(loaded.metadata, "layers")
        @test haskey(loaded.reductions, "pca")
        loaded
    end
    @test size(h5ad_path.counts) == size(experiment.counts)

    annotation_experiment = BioToolkit.SingleCellExperiment(
        [20 1 1 0; 18 1 0 0; 0 19 1 1; 1 17 0 0],
        ["MS4A1", "CD3D", "LYZ", "EPCAM"],
        ["cell1", "cell2", "cell3", "cell4"]
    )
    annotation_labels = [1, 1, 2, 2]
    annotation = BioToolkit.annotate_cell_types(annotation_experiment; labels=annotation_labels, method=:marker, top_n=2)
    @test annotation isa BioToolkit.CellTypeAnnotationResult
    @test length(annotation.assigned_labels) == 2
    @test annotation.cell_labels[1] == annotation.cell_labels[2]
    @test all(label -> label in annotation.reference_labels, annotation.cell_labels)

    reference = BioToolkit.SingleCellExperiment(
        [30 28 1 0; 1 0 25 27; 0 1 1 0; 0 0 0 0],
        ["MS4A1", "CD3D", "LYZ", "EPCAM"],
        ["ref1", "ref2", "ref3", "ref4"]
    )
    reference.metadata["cell_type_labels"] = ["B cell", "B cell", "T cell", "T cell"]
    reference_annotation = BioToolkit.annotate_cell_types(annotation_experiment; labels=annotation_labels, reference=reference, method=:reference)
    @test reference_annotation isa BioToolkit.CellTypeAnnotationResult
    @test length(reference_annotation.assigned_labels) == 2

    ambient = BioToolkit.remove_ambient_rna(experiment; empty_droplets=[3, 4])
    @test ambient isa BioToolkit.AmbientRNARemovalResult
    @test length(ambient.contamination_fraction) == length(experiment.cell_ids)
    @test length(ambient.background_cells) == 2
    @test size(ambient.corrected_counts) == size(experiment.counts)
    @test all(ambient.corrected_counts .<= experiment.counts)

    viewer = BioToolkit.interactive_singlecell_viewer(experiment; embedding=Float64[0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0], labels=[1, 1, 2, 2], k=2)
    @test viewer isa BioToolkit.SingleCellViewer
    @test all(pair -> pair.first in experiment.gene_ids, BioToolkit.top_expressed_genes(experiment, 1; top_n=2))
    @test occursin("cell1", BioToolkit.cell_hover_text(viewer, 1))
    @test BioToolkit.lasso_select_cells(viewer, [( -0.5, -0.5 ), ( 1.5, -0.5 ), ( 1.5, 0.5 ), ( -0.5, 0.5 )]) == [1, 2]
    @test length(BioToolkit.recluster_singlecell_viewer!(viewer; k=2)) == 4

    pca = BioToolkit.run_pca(experiment; normalized=matrix, n_components=2)
    @test size(pca, 2) == 2
    @test haskey(experiment.reductions, "pca")

    experiment.variable_features["vst"] = ["gene1", "gene2", "gene3"]
    projection_model = BioToolkit.fit_singlecell_projection_model(experiment; n_components=2, use_variable_features=true)
    @test projection_model isa BioToolkit.SingleCellProjectionModel
    projected_reference = BioToolkit.project_singlecell(experiment, projection_model)
    @test size(projected_reference) == (4, 2)

    reordered_experiment = BioToolkit.SingleCellExperiment(counts[[3, 1, 4, 2], :], ["gene3", "gene1", "gene4", "gene2"], copy(experiment.cell_ids))
    projected_reordered = BioToolkit.project_singlecell(reordered_experiment, projection_model)
    @test projected_reordered ≈ projected_reference

    projected_via_pca = BioToolkit.run_pca(reordered_experiment; projection_model=projection_model)
    @test projected_via_pca ≈ projected_reference

    pseudotime = BioToolkit.calculate_pseudotime(
        experiment;
        root_cell="cell1",
        embedding=Float64[0.0 0.0; 1.0 0.0; 2.0 0.0; 3.0 0.0],
        k=2,
    )
    @test length(pseudotime.pseudotime) == 4
    @test pseudotime.root_cell == "cell1"
    @test pseudotime.pseudotime[1] == 0.0
    @test pseudotime.pseudotime[end] ≈ 1.0
    @test all(diff(pseudotime.pseudotime) .>= 0)
    @test haskey(experiment.reductions, "pseudotime")

    trajectory_plot = BioToolkit.plot_trajectory(experiment; root_cell="cell1", embedding=Float64[0.0 0.0; 1.0 0.0; 2.0 0.0; 3.0 0.0], k=2)
    @test trajectory_plot isa Plots.Plot

    spliced = counts .+ 1
    unspliced = 2 .* spliced
    BioToolkit.attach_velocity_layers!(experiment; spliced=spliced, unspliced=unspliced)
    velocity = BioToolkit.calculate_rna_velocity(experiment; normalize=false)
    @test size(velocity.velocity) == size(counts)
    @test length(velocity.cell_scores) == size(counts, 2)
    @test length(velocity.latent_time) == size(counts, 2)

    velocity_plot = BioToolkit.plot_rna_velocity(experiment; velocity_result=velocity, embedding=Float64[0.0 0.0; 1.0 0.0; 2.0 0.0; 3.0 0.0], root_cell="cell1")
    @test velocity_plot isa Plots.Plot
    quiver_plot = BioToolkit.plot_rna_velocity_quiver(experiment; velocity_result=velocity, embedding=Float64[0.0 0.0; 1.0 0.0; 2.0 0.0; 3.0 0.0], root_cell="cell1")
    @test quiver_plot isa Plots.Plot

    dynamical_velocity = BioToolkit.calculate_dynamical_rna_velocity(experiment; spliced=spliced, unspliced=unspliced, normalize=false)
    @test dynamical_velocity isa BioToolkit.DynamicalRNAVelocityResult
    @test size(dynamical_velocity.velocity) == size(counts)
    @test length(dynamical_velocity.alpha) == size(counts, 1)
    @test length(dynamical_velocity.latent_time) == size(counts, 2)
    dynamical_velocity_plot = BioToolkit.plot_rna_velocity(experiment; velocity_result=dynamical_velocity, embedding=Float64[0.0 0.0; 1.0 0.0; 2.0 0.0; 3.0 0.0], root_cell="cell1")
    @test dynamical_velocity_plot isa Plots.Plot

    protein_counts = [8 1 2 0; 7 2 1 1; 1 7 2 8; 2 2 6 7]
    protein_experiment = BioToolkit.SingleCellExperiment(protein_counts, ["prot1", "prot2", "prot3", "prot4"], ["cell1", "cell2", "cell3", "cell4"])
    wnn = BioToolkit.weighted_nearest_neighbors(Dict("rna" => experiment, "protein" => protein_experiment); k=2, n_components=2)
    @test wnn isa BioToolkit.WNNResult
    @test length(wnn.combined_neighbors) == 4
    @test all(length(neighborhood) <= 2 for neighborhood in wnn.combined_neighbors)
    @test haskey(wnn.modality_embeddings, "rna")

    perturbation = BioToolkit.predict_perturbation(experiment, "gene1"; normalize=false, top_n=2, diffusion_steps=1)
    @test perturbation isa BioToolkit.PerturbationPredictionResult
    @test perturbation.target_gene == "gene1"
    @test length(perturbation.ranked_genes) == 2
    @test all(gene -> gene in experiment.gene_ids, perturbation.ranked_genes)

    communication_labels = [1, 1, 2, 2]
    communication_pairs = [BioToolkit.LigandReceptorPair("gene1", "gene2", "toy", "manual")]
    communication = BioToolkit.find_cell_communication(experiment, communication_labels; pairs=communication_pairs, normalize=false, min_expression=0.0)
    @test !isempty(communication)
    @test all(result -> result isa BioToolkit.CellCommunicationResult, communication)
    @test all(result -> result.pathway == "toy", communication)
    network = BioToolkit.communication_network(experiment, communication_labels; pairs=communication_pairs, normalize=false, min_expression=0.0)
    @test size(network.score_matrix) == (2, 2)
    @test all(result -> result isa BioToolkit.CellCommunicationResult, network.interactions)
    @test BioToolkit.communication_pathway_summary(network)[1].pathway == "toy"
    communication_plot = BioToolkit.plot_communication_network(network)
    @test communication_plot isa Plots.Plot
    pathway_plot = BioToolkit.plot_communication_pathways(network)
    @test pathway_plot isa Plots.Plot
    report = BioToolkit.rank_ligand_receptor_report(network; pathway="toy", top_n=5)
    @test report isa BioToolkit.LigandReceptorReport
    @test !isempty(report.interactions)
    report_plot = BioToolkit.plot_ligand_receptor_report(report)
    @test report_plot isa Plots.Plot

    archive = mktempdir() do tempdir
        archive_path = joinpath(tempdir, "singlecell.sca")
        saved_archive = BioToolkit.save_singlecell_archive(archive_path, experiment; layers=Dict("spliced" => spliced, "unspliced" => unspliced), embeddings=Dict("pca" => pca))
        loaded_archive = BioToolkit.load_singlecell_archive(archive_path)
        @test haskey(loaded_archive.layers, "spliced")
        @test haskey(loaded_archive.embeddings, "pca")
        @test saved_archive.gene_ids == loaded_archive.gene_ids
        @test saved_archive.spatial_coords == loaded_archive.spatial_coords
        restored = BioToolkit.archive_to_singlecell_experiment(loaded_archive)
        @test Array(restored.counts) == Array(experiment.counts)
        @test haskey(restored.metadata, "layers")
        @test haskey(restored.reductions, "pca")
        @test restored.spatial_coords == experiment.spatial_coords
        browser = BioToolkit.browse_singlecell_archive(loaded_archive; layer="spliced", embedding="pca")
        @test browser isa BioToolkit.SingleCellArchiveBrowser
        @test BioToolkit.archive_layer(browser, "spliced") == spliced
        @test BioToolkit.archive_embedding(browser, "pca") == pca
        loaded_archive
    end
    @test archive isa BioToolkit.SingleCellArchive

    if CUDA.functional()
        gpu_experiment = BioToolkit.gpu_singlecell_experiment(experiment)
        @test gpu_experiment.counts isa CUDA.CuArray
        gpu_normalized = BioToolkit.normalize_counts(gpu_experiment; use_cuda=true)
        @test gpu_normalized isa CUDA.CuArray
        gpu_velocity = BioToolkit.calculate_rna_velocity(gpu_experiment; spliced=spliced, unspliced=unspliced, normalize=false, use_cuda=true)
        @test length(gpu_velocity.cell_scores) == size(counts, 2)
    end

    @test length(BioToolkit.default_ligand_receptor_pairs()) >= 20

    variable_features = BioToolkit.find_variable_features(experiment; n_features=2, method=:variance)
    @test length(variable_features) == 2
    @test haskey(experiment.variable_features, "variance")

    neighbors = BioToolkit.find_neighbors(experiment; embedding=pca, k=2)
    @test length(neighbors) == 4
    @test haskey(experiment.neighbors, "neighbors")

    seurat_labels = BioToolkit.find_clusters(experiment; neighbors=neighbors, method=:leiden)
    @test length(seurat_labels) == 4
    @test sum(values(BioToolkit.summarize_clusters(seurat_labels))) == 4

    louvain_labels = BioToolkit.find_clusters(experiment; neighbors=neighbors, method=:louvain)
    @test length(louvain_labels) == 4

    umap = BioToolkit.run_umap(experiment; embedding=pca)
    @test size(umap, 2) <= 2

    labels = BioToolkit.cluster_cells(experiment; embedding=pca, method=:graph, k=2, random_seed=1)
    @test length(labels) == 4
    @test sum(values(BioToolkit.summarize_clusters(labels))) == 4

    marker_labels = BioToolkit.cluster_cells(experiment; embedding=pca, method=:kmeans, n_clusters=2, random_seed=1)
    @test length(unique(marker_labels)) == 2

    markers = BioToolkit.find_cluster_markers(experiment, marker_labels; cluster_id=marker_labels[1], min_total=0, shrink=false)
    @test !isempty(markers)
    @test all(result -> result isa BioToolkit.DEResult, markers)

    wilcox_markers = BioToolkit.find_markers(experiment, marker_labels; ident_1=marker_labels[1], ident_2=marker_labels[2], min_total=0)
    @test !isempty(wilcox_markers)
    @test all(result -> result isa BioToolkit.DEResult, wilcox_markers)

    pb_counts = [10 11 12 13 30 29 31 32; 5 4 6 5 6 7 5 6; 2 3 2 3 8 7 8 7; 1 1 1 1 1 1 1 1]
    pb_experiment = BioToolkit.SingleCellExperiment(pb_counts, ["gene1", "gene2", "gene3", "gene4"], ["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"])
    pb_labels = [1, 1, 1, 1, 2, 2, 2, 2]
    pb_samples = ["rep1", "rep1", "rep2", "rep2", "rep3", "rep3", "rep4", "rep4"]
    pseudobulk_markers = BioToolkit.find_markers(pb_experiment, pb_labels; ident_1=1, ident_2=2, test=:deseq2, pseudobulk=true, sample_labels=pb_samples, min_total=0, shrink=false)
    @test !isempty(pseudobulk_markers)
    @test all(result -> result isa BioToolkit.DEResult, pseudobulk_markers)

    marker_summary = BioToolkit.cluster_marker_summary(experiment, marker_labels; min_total=0, shrink=false, top_n=2)
    @test !isempty(marker_summary)
    @test all(value -> isa(value, Vector{BioToolkit.DEResult}), values(marker_summary))
    @test all(length(value) <= 2 for value in values(marker_summary))

    integrated = BioToolkit.integrate_data([experiment, experiment]; n_components=2)
    @test size(integrated.corrected_matrix) == (4, 8)
    @test length(integrated.gene_ids) == 4
    @test length(integrated.cell_ids) == 8
end