using DataFrames
using Turing

@testset "Microbiome" begin
    counts = [10 0 3;
              0 5 1;
              2 1 0]
    taxonomy = DataFrame(
        taxon_id=["tax1", "tax2", "tax3"],
        kingdom=["Bacteria", "Bacteria", "Bacteria"],
        phylum=["P1", "P1", "P2"],
        class=["C1", "C1", "C2"],
        order=["O1", "O1", "O2"],
        family=["F1", "F1", "F2"],
        genus=["G1", "G2", "G3"],
        species=["S1", "S2", "S3"],
    )
    metadata = DataFrame(sample_id=["s1", "s2", "s3"], treatment=["case", "case", "ctrl"])
    tree = BioToolkit.PhyloTree([
        BioToolkit.PhyloTree("tax1"; branch_length=0.1),
        BioToolkit.PhyloTree("tax2"; branch_length=0.2),
        BioToolkit.PhyloTree("tax3"; branch_length=0.3),
    ]; name="root")

    profile = BioToolkit.CommunityProfile(counts, ["tax1", "tax2", "tax3"], ["s1", "s2", "s3"], taxonomy, tree, metadata)
    subset = profile[1:2, 1:2]
    @test subset.counts.gene_ids == ["tax1", "tax2"]
    @test subset.counts.sample_ids == ["s1", "s2"]
    @test subset.taxonomy.taxon_id == ["tax1", "tax2"]
    @test subset.metadata.sample_id == ["s1", "s2"]

    clr = BioToolkit.clr_transform(profile)
    @test size(clr) == (3, 3)
    @test all(abs.(sum(clr, dims=1)) .< 1e-8)

    ilr = BioToolkit.ilr_transform(profile)
    @test size(ilr) == (2, 3)

    bray_curtis = BioToolkit.bray_curtis(profile)
    @test size(bray_curtis) == (3, 3)
    @test bray_curtis ≈ bray_curtis'
    @test all(iszero, diag(bray_curtis))
    @test BioToolkit.pairwise_bray_curtis(profile) ≈ bray_curtis

    unifrac = BioToolkit.unifrac(profile)
    @test size(unifrac) == (3, 3)
    @test unifrac ≈ unifrac'
    @test all(iszero, diag(unifrac))
    @test BioToolkit.pairwise_unifrac(profile) ≈ unifrac

    weighted_unifrac = BioToolkit.weighted_unifrac(profile)
    @test size(weighted_unifrac) == (3, 3)
    @test weighted_unifrac ≈ weighted_unifrac'
    @test BioToolkit.pairwise_unifrac(profile; weighted=true) ≈ weighted_unifrac

    @test length(BioToolkit.shannon_entropy(profile)) == 3
    @test length(BioToolkit.simpson_index(profile)) == 3
    @test length(BioToolkit.faith_pd(profile)) == 3

    pcoa = BioToolkit.pcoa(bray_curtis)
    @test size(pcoa.coordinates, 1) == 3
    @test size(pcoa.coordinates, 2) == 2

    nmds = BioToolkit.nmds(bray_curtis; dimensions=2, n_starts=3, maxiters=50, random_seed=1)
    @test size(nmds.coordinates) == (3, 2)
    @test nmds.stress >= 0.0

    ancom = BioToolkit.ancom(profile, metadata.treatment)
    @test length(ancom.taxon_ids) == 3
    @test length(ancom.w_stat) == 3
    @test all((0 .<= ancom.qvalue) .& (ancom.qvalue .<= 1))

    songbird = BioToolkit.songbird(profile, hcat(ones(3), [1.0, 1.0, 0.0]))
    @test size(songbird.coefficients, 1) == 2
    @test size(songbird.coefficients, 2) == 2
    @test songbird.iterations >= 0

    network = BioToolkit.cooccurrence_network(profile; threshold=0.1)
    @test length(network.taxa) == 3
    @test length(network.taxa) == 3
    @test !isempty(network.weights)
    network_plot = BioToolkit.network_plot(network; layout=:spring, edge_scale=3.0, show_labels=false)
    @test network_plot !== nothing

    source_profiles = [0.7 0.1; 0.2 0.7; 0.1 0.2]
    source_model = BioToolkit.source_tracking_model([10, 6, 4], source_profiles)
    @test source_model !== nothing
    source_result = BioToolkit.source_tracking([10, 6, 4], source_profiles; draws=2)
    @test source_result isa BioToolkit.SourceTrackingResult
    @test length(source_result.mean_proportions) == 2
    @test length(source_result.median_proportions) == 2
    @test all((0 .<= source_result.mean_proportions) .& (source_result.mean_proportions .<= 1))
end
