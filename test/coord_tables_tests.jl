using Test
using DataFrames

@testset "Coordinate and table integration" begin
    @testset "BED, GFF, and FASTA coordinates" begin
        mktempdir() do dir
            bed_path = joinpath(dir, "region.bed")
            gff_path = joinpath(dir, "region.gff")

            open(bed_path, "w") do io
                write(io, "chr1\t9\t20\n")
            end

            open(gff_path, "w") do io
                write(io, "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tID=gene1\n")
            end

            bed = first(BioToolkit.read_bed(bed_path))
            gff = first(BioToolkit.read_gff(gff_path))

            bed_interval = BioToolkit.GenomicRanges.GenomicInterval(bed)
            gff_interval = BioToolkit.GenomicRanges.GenomicInterval(gff)

            @test bed_interval.chrom == gff_interval.chrom
            @test bed_interval.left == 10
            @test bed_interval.right == 20

            interval_table = DataFrame([bed_interval])
            @test Symbol.(names(interval_table)) == [:chrom, :start, :stop, :strand]
            @test interval_table.chrom[1] == "chr1"
            @test interval_table.start[1] == 10
            @test interval_table.stop[1] == 20

            feature = BioToolkit.SeqFeatureLite("gene", BioToolkit.FeatureLocationLite(10, 20; strand=1); id="gene1")
            feature_table = DataFrame([feature])
            @test Symbol.(names(feature_table)) == [:feature_type, :start, :stop, :strand, :id]
            @test feature_table.start[1] == 10
            @test feature_table.stop[1] == 20

            fasta_path = joinpath(dir, "slice.fasta")
            open(fasta_path, "w") do io
                write(io, ">chr1\n")
                write(io, "ACGTACGTACGTACGTACGTACGTACGTAC\n")
            end

            index = BioToolkit.fasta_index(fasta_path)
            slice = BioToolkit.fetch_fasta_sequence(fasta_path, index["chr1"], 10, 20)
            @test length(slice) == 11
            @test slice == "CGTACGTACGT"
        end
    end

    @testset "Result table adapters" begin
        de = BioToolkit.DEResult("gene1", 12.0, 1.5, 0.2, 3.0, 0.01, 0.02)
        de_table = DataFrame([de])
        @test de_table.gene_id[1] == "gene1"
        @test de_table.log2_fold_change[1] == 1.5

        gwas = BioToolkit.GWASResult(
            ["rs1", "rs2"],
            ["chr1", "chr2"],
            [10, 20],
            [("A", "G"), ("C", "T")],
            ["gene1", "gene2"],
            [0.1, 0.2],
            [0.01, 0.02],
            [10.0, 20.0],
            [0.001, 0.002],
            50,
            ["cov1"],
            "phenotype",
            "linear",
        )
        gwas_table = DataFrame(gwas)
        @test Symbol.(names(gwas_table)) == [:snp_id, :chromosome, :position, :allele1, :allele2, :gene_id, :beta, :standard_error, :zscore, :pvalue, :sample_size, :phenotype_name, :method]
        @test gwas_table.snp_id == ["rs1", "rs2"]

        peak_table = DataFrame([BioToolkit.Peak("peak1", "chr1", 10, 20, 15, 9.5, 0.001, 0.01)])
        @test peak_table.chrom[1] == "chr1"

        motif_table = DataFrame([BioToolkit.MotifHit(4, Int8(1), 7.5, "ACGT")])
        @test motif_table.start[1] == 4
        @test motif_table.score[1] == 7.5
    end

    @testset "Parser hardening" begin
        @test BioToolkit.parse_bed_record("chr1\t1") === nothing
        @test BioToolkit.parse_gff_record("chr1\tsrc\tgene") === nothing
        @test_throws ArgumentError BioToolkit.parse_bed_record("chr1\tnot_an_int\t20")
        @test_throws ArgumentError BioToolkit.parse_gff_record("chr1\tsrc\tgene\tfoo\t20\t.\t+\t.\tID=x")

        mktempdir() do dir
            bed_path = joinpath(dir, "mixed.bed")
            open(bed_path, "w") do io
                write(io, "chr1\t9\t20\n")
                write(io, "chr1\tbad\t20\n")
                write(io, "chr1\t30\t40\n")
            end
            @test_throws ArgumentError BioToolkit.read_bed(bed_path)

            gff_path = joinpath(dir, "mixed.gff")
            open(gff_path, "w") do io
                write(io, "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tID=gene1\n")
                write(io, "chr1\tsrc\tgene\tfoo\t20\t.\t+\t.\tID=bad\n")
            end
            @test_throws ArgumentError BioToolkit.read_gff(gff_path)
        end
    end
end