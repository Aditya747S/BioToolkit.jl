using Test

@testset "BAM support" begin
    mktempdir() do dir
        bam_path = joinpath(dir, "sample.bam")

        header = BioToolkit.BamHeader([
            BioToolkit.BamReference("chr1", 1_000),
            BioToolkit.BamReference("chr2", 500),
        ])

        records = [
            BioToolkit.BamRecord("read1", "chr1", 9, [BioToolkit.BamCigarOp(10, 'M')], "ACGTACGTAA"; quality="IIIIIIIIII"),
            BioToolkit.BamRecord("read2", "chr2", 19, [BioToolkit.BamCigarOp(6, 'M')], "GATTAC"; quality="JJJJJJ"),
        ]

        bam = BioToolkit.BamFile(records; header=header)
        BioToolkit.write_bam(bam_path, bam)

        @test isfile(bam_path)
        @test isfile(string(bam_path, ".bai"))

        roundtrip = BioToolkit.read_bam(bam_path)
        @test roundtrip.header == header
        @test roundtrip.records == records

        stream_reader = BioToolkit.read_bam(bam_path; materialize=false)
        @test stream_reader isa BioToolkit.AbstractBamRecordReader
        @test stream_reader.header == header
        @test collect(stream_reader) == records

        region = BioToolkit.GenomicInterval("chr1", 10, 15)
        region_hits = BioToolkit.read_bam(bam_path, region)
        @test length(region_hits) == 1
        @test region_hits.records[1].qname == "read1"
        @test region_hits.records[1].refname == "chr1"

        region_stream = BioToolkit.read_bam(bam_path, region; materialize=false)
        @test region_stream isa BioToolkit.AbstractBamRecordReader
        @test collect(region_stream) == [records[1]]

        missing_region_stream = BioToolkit.read_bam(bam_path, BioToolkit.GenomicInterval("chrX", 1, 50); materialize=false)
        @test missing_region_stream isa BioToolkit.AbstractBamRecordReader
        @test isempty(collect(missing_region_stream))

        @test_throws ArgumentError BioToolkit.read_bam(joinpath(dir, "unsupported.cram"))
        @test_throws ArgumentError BioToolkit.write_bam(joinpath(dir, "unsupported.cram"), bam)
    end
end