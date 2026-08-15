using Test
using Random
using Statistics
using DataFrames
using BioToolkit

@testset "Deep Learning Module Comprehensive Tests" begin
    Random.seed!(42)

    # 1. Single-Cell Matrix Fixtures
    n_genes = 20
    n_cells = 30
    counts = abs.(randn(n_genes, n_cells) .* 10.0)
    gene_ids = ["Gene_$i" for i in 1:n_genes]
    cell_ids = ["Cell_$i" for i in 1:n_cells]

    @testset "scVI Latent Model & Embeddings" begin
        emb = scvi_like_embedding(counts; n_latent=5, backend=:cpu)
        @test size(emb.latent) == (n_cells, 5)
        @test size(emb.loadings) == (n_genes, 5)

        model = fit_scvi_model(counts; n_latent=5, backend=:cpu, max_epochs=2)
        @test size(model.latent) == (n_cells, 5)
        @test model.training_history.epochs == 2

        trans = transform_scvi(model, counts)
        @test size(trans) == (n_cells, 5)

        pred_lat = predict_scvi_latent(model, counts)
        @test size(pred_lat) == (n_cells, 5)
    end

    @testset "CellAssign Model & Mapping" begin
        marker_sets = Dict(
            "T_cell" => ["Gene_1", "Gene_2", "Gene_3"],
            "B_cell" => ["Gene_4", "Gene_5", "Gene_6"]
        )

        mapping = cellassign_like_mapping(counts, gene_ids, marker_sets)
        @test length(mapping.predicted_label) == n_cells
        @test size(mapping.score_matrix) == (n_cells, 2)
        @test sort(mapping.label_order) == ["B_cell", "T_cell"]

        fit_ca = fit_cellassign_model(counts, gene_ids, marker_sets)
        pred_ca = predict_cell_types(fit_ca, counts)
        @test length(pred_ca.predicted_label) == n_cells
        @test size(pred_ca.probability) == (n_cells, 2)
    end

    @testset "Sequence Transformer Embeddings" begin
        dna_strings = ["ACGTACGTACGT", "TGCATGCATGCA", "GGCCGGCCGGCC"]
        dna_seqs = BioToolkit.DNASeq.(dna_strings)

        # Geneformer
        gf_str = geneformer_like_embedding(dna_strings; k=3, dim=16, threaded=true)
        gf_seq = geneformer_like_embedding(dna_seqs; k=3, dim=16, threaded=true)
        @test size(gf_str) == (3, 16)
        @test size(gf_seq) == (3, 16)

        # scGPT
        gpt_str = scgpt_like_embedding(dna_strings; token_dim=24, k=3, threaded=true)
        gpt_seq = scgpt_like_embedding(dna_seqs; token_dim=24, k=3, threaded=true)
        @test size(gpt_str) == (3, 24)
        @test size(gpt_seq) == (3, 24)

        # scBERT
        bert_str = scbert_like_embedding(dna_strings; dim=20, max_len=32, threaded=true)
        bert_seq = scbert_like_embedding(dna_seqs; dim=20, max_len=32, threaded=true)
        @test size(bert_str) == (3, 20)
        @test size(bert_seq) == (3, 20)

        # High-level sequence transformer fit
        st_fit = fit_sequence_transformer_embedding(dna_strings; token_dim=16, k=3)
        @test size(st_fit.latent) == (3, 16)
    end

    @testset "Attention GRN & GNN GRN" begin
        grn_def = attention_grn(counts; top_k=5, temperature=1.0)
        @test "source" in names(grn_def)
        @test "target" in names(grn_def)
        @test nrow(grn_def) > 0

        grn_genes = attention_grn(counts; gene_ids=gene_ids, top_k=5, temperature=1.0)
        @test grn_genes.source[1] in gene_ids

        gnn_grn = gene_regulatory_network_gnn(counts; gene_ids=gene_ids, k_neighbors=5, n_propagation=2, top_k_edges=20)
        @test "source_gene" in names(gnn_grn.edges)
        @test gnn_grn.edges.source_gene[1] in gene_ids
        @test nrow(gnn_grn.edges) <= 20
    end

    @testset "Batch Correction & Contrastive Embeddings" begin
        batches = vcat(fill("Batch1", 15), fill("Batch2", 15))
        bc = batch_corrected_latent(counts, batches; n_latent=4, backend=:cpu)
        @test size(bc.latent) == (n_cells, 4)
        @test length(bc.batch_levels) == 2

        con = contrastive_cell_embedding(counts; n_latent=4, dropout=0.1, seed=123)
        @test size(con.latent) == (n_cells, 4)
        @test isfinite(con.alignment_loss)
    end

    @testset "scGen Perturbation Response" begin
        control_expr = abs.(randn(n_genes, 15))
        treated_expr = abs.(randn(n_genes, 15))
        query_ctrl   = abs.(randn(n_genes, 10))

        scg = scgen_like_perturbation(control_expr, treated_expr, query_ctrl; n_latent=4)
        @test size(scg.predicted_treated) == (n_genes, 10)

        fit_scg = fit_scgen_model(control_expr, treated_expr)
        pred_scg = predict_perturbation_response(fit_scg, query_ctrl)
        @test size(pred_scg.predicted_response) == (n_genes, 10)
    end

    @testset "GraphSCA Label Transfer" begin
        adjacency = rand(n_cells, n_cells)
        adjacency = 0.5 .* (adjacency .+ adjacency')
        known_labels = vcat(fill("T_cell", 10), fill("B_cell", 10), fill("unknown", 10))

        transfer = graphsca_label_transfer(counts, adjacency, known_labels; alpha=0.7, n_iter=10)
        @test length(transfer.predicted_label) == n_cells
        @test all(l -> l in ["B_cell", "T_cell"], transfer.predicted_label)
    end

    @testset "Cell Type Denoising (MAGIC)" begin
        denoised = cell_type_denoising(counts; k=5, t=2, n_pcs=4)
        @test size(denoised.denoised) == (n_genes, n_cells)
        @test size(denoised.diffusion_operator) == (n_cells, n_cells)
    end

    @testset "Sparse Autoencoder Features" begin
        sae = sparse_autoencoder_features(counts; n_features=8, sparsity_k=2, n_iter=20, lr=1e-2)
        @test size(sae.features) == (n_genes, 8)
        @test size(sae.feature_activations) == (n_cells, 8)
        @test isfinite(sae.reconstruction_loss)
    end

    @testset "Trajectory Neural ODE" begin
        spliced = abs.(randn(n_genes, n_cells) .* 5.0)
        unspliced = abs.(randn(n_genes, n_cells) .* 2.0)

        traj = trajectory_neural_ode(spliced, unspliced; n_latent=4, n_steps=5, dt=0.1)
        @test length(traj.pseudotime) == n_cells
        @test size(traj.latent_velocity) == (n_cells, 4)
        @test size(traj.trajectory) == (n_cells, 4, 6)
    end

    @testset "Weighted Nearest Neighbour Multi-Modal" begin
        mod1 = abs.(randn(15, n_cells))
        mod2 = abs.(randn(10, n_cells))
        wnn = multimodal_wnn_embedding([mod1, mod2]; n_latent=5, n_neighbors=5)
        @test size(wnn.latent, 1) == n_cells
        @test size(wnn.modality_weights) == (n_cells, 2)
        @test size(wnn.nn_graph) == (n_cells, n_cells)
    end

    @testset "Protein Sequence Embeddings" begin
        peptides = ["MKTIIALSYIFCLVFA", "ACDEFGHIKLMNPQRSTVWY"]
        aa_seqs = BioToolkit.AASeq.(peptides)

        p_emb1 = protein_sequence_embedding(peptides; dim=32, k=3, include_properties=true)
        p_emb2 = protein_sequence_embedding(aa_seqs; dim=32, k=3, include_properties=true)
        @test size(p_emb1) == (2, 32)
        @test size(p_emb2) == (2, 32)
    end

    @testset "Zero-Shot Cell Annotation" begin
        archetypes = Dict(
            "T_cell" => ["Gene_1", "Gene_2"],
            "B_cell" => ["Gene_3", "Gene_4"]
        )
        zs = zero_shot_cell_annotation(counts, gene_ids, archetypes)
        @test nrow(zs.cells) == n_cells
        @test size(zs.score_matrix) == (n_cells, 2)
    end

    @testset "Self-Supervised Pretraining & Cell Cycle & Deep Factorization" begin
        ss = self_supervised_pretraining(counts; mask_fraction=0.1, n_latent=4, n_iter=10)
        @test size(ss.latent) == (n_cells, 4)
        @test length(ss.reconstruction_loss_history) == 10

        cc = cell_cycle_regression(counts, gene_ids; s_genes=["Gene_1", "Gene_2"], g2m_genes=["Gene_3", "Gene_4"])
        @test size(cc.corrected_counts) == (n_genes, n_cells)
        @test length(cc.phase) == n_cells

        mod1 = abs.(randn(10, n_cells))
        mod2 = abs.(randn(12, n_cells))
        dfe = deep_factorization_embedding([mod1, mod2]; n_factors=4, n_iter=10)
        @test size(dfe.factors) == (n_cells, 4)
        @test length(dfe.loadings) == 2
    end

    @testset "Interactive HTML Visualizers" begin
        emb = scvi_like_embedding(counts; n_latent=3)
        html_emb = embedding_to_html(emb; labels=vcat(fill("TypeA", 15), fill("TypeB", 15)), title="Test Latent Plot")
        @test contains(html_emb, "plotly")
        @test contains(html_emb, "Test Latent Plot")

        tmp_emb = tempname() * ".html"
        export_embedding_html(emb, tmp_emb)
        @test isfile(tmp_emb)
        rm(tmp_emb, force=true)

        grn = attention_grn(counts; top_k=5)
        html_grn = grn_to_html(grn; title="Test GRN Plot")
        @test contains(html_grn, "vis-network")
        @test contains(html_grn, "Test GRN Plot")

        tmp_grn = tempname() * ".html"
        export_grn_html(grn, tmp_grn)
        @test isfile(tmp_grn)
        rm(tmp_grn, force=true)

        spliced = abs.(randn(n_genes, n_cells))
        unspliced = abs.(randn(n_genes, n_cells))
        traj = trajectory_neural_ode(spliced, unspliced; n_latent=3, n_steps=3)
        html_traj = trajectory_to_html(traj; title="Test Trajectory Plot")
        @test contains(html_traj, "plotly")
        @test contains(html_traj, "Test Trajectory Plot")

        archetypes = Dict("T" => ["Gene_1"], "B" => ["Gene_2"])
        zs = zero_shot_cell_annotation(counts, gene_ids, archetypes)
        html_annot = cell_type_annotation_html(zs; title="Test Cell Type Report")
        @test contains(html_annot, "plotly")
        @test contains(html_annot, "Test Cell Type Report")
    end
end
