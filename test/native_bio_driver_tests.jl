using Test
using BigWig

@testset "Native bio driver helpers" begin
    @testset "Sequence utilities" begin
        dimers = BioToolkit.check_primer_dimers("AAAAAA", "TTTTTT")
        @test dimers.has_dimer
        @test dimers.max_3prime_match >= 4

        genome = Dict("chr1" => "TTTACGATGCCTT")
        amplicons = BioToolkit.pcr_in_silico("ACGA", "AGGC", genome)
        @test !isempty(amplicons)
        @test amplicons[1].chrom == "chr1"

        guide_score = BioToolkit.score_grna_efficiency("GAGTCCGAGCAGAAGAAGAATGG")
        @test 0.0 <= guide_score <= 1.0
    end

    @testset "Coverage export" begin
        mktempdir() do dir
            output_path = joinpath(dir, "coverage.bw")
            written = BioToolkit.write_bigwig([0, 0, 5, 5, 1, 1], output_path; chrom="chr1")
            @test isfile(written)

            open(written, "r") do io
                reader = BigWig.Reader(io)
                @test BigWig.chromlist(reader)[1][1] == "chr1"
                @test BigWig.chromlist(reader)[1][2] == 6
                @test BigWig.mean(reader, "chr1", 1, 6; usezoom=false) == 2.0
                @test BigWig.maximum(reader, "chr1", 1, 6; usezoom=false) == 5.0
                @test BigWig.minimum(reader, "chr1", 1, 6; usezoom=false) == 0.0
            end
        end

        mktempdir() do dir
            output_path = joinpath(dir, "multi.bw")
            written = BioToolkit.write_bigwig(Dict("chr2" => [2, 2, 2], "chr1" => [0, 4, 4, 0]), output_path)
            @test isfile(written)

            open(written, "r") do io
                reader = BigWig.Reader(io)
                @test BigWig.chromlist(reader) == [("chr1", 4), ("chr2", 3)]
                @test BigWig.mean(reader, "chr1", 1, 4; usezoom=false) == 2.0
                @test BigWig.maximum(reader, "chr2", 1, 3; usezoom=false) == 2.0
            end
        end
    end

    @testset "Single-cell QC" begin
        counts = [10 0 2 0; 0 10 1 2; 4 0 8 1]
        genes = ["MCM5", "HMGB2", "PCNA"]
        cells = ["cell1", "cell2", "cell3", "cell4"]
        experiment = BioToolkit.SingleCellExperiment(counts, genes, cells)

        doublets = BioToolkit.detect_doublets(experiment; n_simulated=10, k=3)
        @test length(doublets.scores) == length(cells)
        @test length(doublets.doublets) == length(cells)

        integrated = BioToolkit.integrate_batches(counts, ["batch1", "batch1", "batch2", "batch2"])
        @test size(integrated.corrected_matrix) == size(counts)

        cycle = BioToolkit.score_cell_cycle(experiment)
        @test length(cycle.phase) == length(cells)
    end

    @testset "Annotation bridge" begin
        variant = BioToolkit.VariantTextRecord("chr1", 4, "var1", "G", "A", missing)
        feature = BioToolkit.GffRecord("chr1", "src", "CDS", 1, 6, missing, "+", missing, "ID=gene1", Dict("ID" => ["gene1"]))
        annotations = BioToolkit.annotate_variants([variant], [feature]; reference_sequences=Dict("chr1" => "ATGGAA"))
        @test annotations[1].consequence in ("Missense", "Synonymous", "Stop-Gain", "Stop-Loss")
        @test annotations[1].gene == "gene1"
    end

    @testset "Interval and structure helpers" begin
        intervals = [BioToolkit.GenomicInterval("chr1", 5, 10)]
        gaps = BioToolkit.complement(intervals, Dict("chr1" => 20))
        @test length(gaps.intervals) == 2

        residue1 = BioToolkit.Residue("ALA", 1, ' ', [BioToolkit.Atom(1, "CA", 0.0, 0.0, 0.0), BioToolkit.Atom(2, "CB", 1.0, 0.0, 0.0)])
        residue2 = BioToolkit.Residue("GLY", 2, ' ', [BioToolkit.Atom(3, "CA", 0.0, 0.0, 4.0)])
        ligand = BioToolkit.Residue("LIG", 3, ' ', [BioToolkit.Atom(4, "C1", 0.0, 0.0, 4.5; hetatm=true)])
        chain_a = BioToolkit.Chain("A", [residue1, residue2])
        chain_b = BioToolkit.Chain("B", [ligand])

        interface = BioToolkit.calculate_interface_residues(chain_a, chain_b)
        @test !isempty(interface)

        nearby = BioToolkit.residues_within_radius(chain_a, ligand; radius=5.0)
        @test !isempty(nearby)
    end
end