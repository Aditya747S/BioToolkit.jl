using Test
using BioToolkit

@testset "Spatial transcriptomics deconvolution" begin
    genes = ["G1", "G2", "G3", "G4", "G5", "G6"]

    reference_counts = [
        80 75 90  5  4  6;
        60 55 65  6  5  4;
        10  8  9 12 10 11;
        12 10 11 10 11  9;
         3  2  4 70 80 75;
         2  3  2 65 70 60
    ]
    reference_cells = ["r1", "r2", "r3", "r4", "r5", "r6"]
    reference_labels = [1, 1, 1, 2, 2, 2]

    spatial_counts = [
        85  5 40;
        62  5 32;
         9 11 10;
        11 10 10;
         3 78 35;
         2 67 30
    ]
    spot_ids = ["s1", "s2", "s3"]
    coords = [
        0.0 0.0;
        1.0 0.0;
        0.5 0.8
    ]

    reference = SingleCellExperiment(reference_counts, genes, reference_cells)
    spatial_sce = SingleCellExperiment(spatial_counts, genes, spot_ids; spatial_coords=coords)

    @testset "Reference profile construction" begin
        reference_profiles = build_reference_matrix(reference, reference_labels; marker_top_n=4, min_total=1)

        @test length(reference_profiles.cell_type_ids) == 2
        @test size(reference_profiles.profiles, 2) == 2
        @test !isempty(reference_profiles.genes)
        @test all(sum(reference_profiles.profiles[:, idx]) ≈ 1.0 for idx in 1:size(reference_profiles.profiles, 2))

        spatial = SpatialExperiment(spatial_sce)
        @test size(spatial.spatial_coords) == (3, 2)
    end

    @testset "RCTD-like deconvolution" begin
        spatial = SpatialExperiment(spatial_sce)
        result = rctd_deconvolution(spatial, reference, reference_labels; marker_top_n=4, min_total=1, maxiter=120)

        @test result.method == :rctd
        @test size(result.cell_type_fractions) == (3, 2)
        @test all(abs(sum(result.cell_type_fractions[i, :]) - 1.0) < 1e-6 for i in 1:3)
        @test all(isfinite, result.residuals)

        type1 = findfirst(==("1"), result.cell_type_ids)
        type2 = findfirst(==("2"), result.cell_type_ids)
        @test result.cell_type_fractions[1, type1] > 0.65
        @test result.cell_type_fractions[2, type2] > 0.65
        @test result.cell_type_fractions[3, type1] > 0.25
        @test result.cell_type_fractions[3, type2] > 0.25
    end

    @testset "Cell2location-like MAP deconvolution" begin
        spatial = SpatialExperiment(spatial_sce, coords)
        result = cell2location_deconvolution(spatial, reference, reference_labels; marker_top_n=4, min_total=1, maxiter=120, dirichlet_alpha=2.5)

        @test result.method == :cell2location
        @test size(result.cell_type_fractions) == (3, 2)
        @test all(abs(sum(result.cell_type_fractions[i, :]) - 1.0) < 1e-6 for i in 1:3)
        @test all(isfinite, result.residuals)

        type1 = findfirst(==("1"), result.cell_type_ids)
        type2 = findfirst(==("2"), result.cell_type_ids)
        @test result.cell_type_fractions[1, type1] > result.cell_type_fractions[1, type2]
        @test result.cell_type_fractions[2, type2] > result.cell_type_fractions[2, type1]
    end
end
