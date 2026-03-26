using DataFrames
using SparseArrays

@testset "SystemsBio" begin
    counts = BioToolkit.CountMatrix(
        sparse(Int[
            10 12 11  9;
             4  5  4  6;
            20 19 21 18;
             8  7  9  8;
        ]),
        ["g1", "g2", "g3", "g4"],
        ["s1", "s2", "s3", "s4"],
    )

    soft = BioToolkit.pick_soft_threshold(counts; powers=2:4)
    @test soft isa BioToolkit.SoftThresholdResult
    @test soft.best_power in soft.powers

    network = BioToolkit.find_modules(counts; power=2, tom_threshold=0.05, edge_threshold=0.05)
    @test network isa BioToolkit.GeneNetwork
    @test length(network.node_to_gene) == 4
    @test length(network.modules) == 4
    @test length(network.connectivity) == 4

    eigengenes = BioToolkit.module_eigengenes(counts, network.modules)
    @test size(eigengenes, 1) == 4
    @test size(eigengenes, 2) == length(unique(network.modules))

    tree = BioToolkit.module_dendrogram(network)
    @test tree isa BioToolkit.PhyloTree
    @test !isempty(BioToolkit.get_terminals(tree))

    design = [1.0 0.0;
              1.0 0.0;
              1.0 1.0;
              1.0 1.0]
    limma = BioToolkit.limma_fit(counts, design; coefficient_names=["intercept", "group"])
    @test limma isa BioToolkit.LimmaResult
    @test size(limma.coefficients, 1) == 4
    @test size(limma.coefficients, 2) == 2

    deresults = BioToolkit.limma_deresults(limma; coefficient_index=2)
    @test length(deresults) == 4
    @test all(result -> result isa BioToolkit.DEResult, deresults)

    db = BioToolkit.build_annotation_database([
        BioToolkit.EnrichmentTerm("TERM:1", "toy pathway", "GO", ["g1", "g2"], String[]),
        BioToolkit.EnrichmentTerm("TERM:2", "background", "GO", ["g3", "g4"], String[]),
    ])
    gsea_results = BioToolkit.gsea(deresults, db; permutations=25, min_size=1)
    @test all(result -> result isa BioToolkit.EnrichmentResult, gsea_results)
    @test !isempty(gsea_results)

    dag = BioToolkit.infer_network(permutedims(Matrix{Float64}(counts.counts)); gene_ids=counts.gene_ids)
    @test dag isa BioToolkit.NetworkInferenceResult
    @test length(dag.node_to_gene) == 4

    assays = [
        [1.0 2.0 3.0; 2.0 3.0 4.0; 3.0 4.0 5.0; 4.0 5.0 6.0],
        [2.0 1.0 0.0; 3.0 2.0 1.0; 4.0 3.0 2.0; 5.0 4.0 3.0],
    ]
    mofa = BioToolkit.multi_omics_factor_analysis(assays; n_factors=2)
    @test mofa isa BioToolkit.MultiOmicsFactorAnalysisResult
    @test size(mofa.factors, 1) == 4
    @test size(mofa.factors, 2) == 2
    @test length(mofa.assay_loadings) == 2
end
