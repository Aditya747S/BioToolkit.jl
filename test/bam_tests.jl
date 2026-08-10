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

        empty_region = BioToolkit.GenomicInterval("chr1", 800, 900)
        empty_reader = BioToolkit.stream_bam(bam_path, empty_region)
        @test empty_reader isa BioToolkit.AbstractBamRecordReader
        @test isempty(collect(empty_reader))

        header_3 = BioToolkit.BamHeader([
            BioToolkit.BamReference("chr1", 1_000),
            BioToolkit.BamReference("chr2", 500),
            BioToolkit.BamReference("chr3", 500),
        ])
        bam_3 = BioToolkit.BamFile([records[1]]; header=header_3)
        bam3_path = joinpath(dir, "sample3.bam")
        BioToolkit.write_bam(bam3_path, bam_3)
        no_chunks_reader = BioToolkit.stream_bam(bam3_path, BioToolkit.GenomicInterval("chr3", 1, 100))
        @test no_chunks_reader isa BioToolkit.BamEmptyReader
        @test isempty(collect(no_chunks_reader))

        unmapped_rec = BioToolkit.BamRecord("unmapped_read", "chr1", 100, BioToolkit.BamCigarOp[], "ACGT"; flag=BioToolkit.SAM_FLAG_UNMAPPED)
        @test BioToolkit._bam_query_span(unmapped_rec) == (0, -1)
        ic = BioToolkit.granges([unmapped_rec, records[1]])
        @test length(ic) == 1

        rec_copy = BioToolkit.BamRecord("read1", "chr1", 9, [BioToolkit.BamCigarOp(10, 'M')], "ACGTACGTAA"; quality="IIIIIIIIII")
        @test records[1] == rec_copy
        @test hash(records[1]) == hash(rec_copy)
        s = Set([records[1]])
        @test rec_copy in s

        sam_path = joinpath(dir, "sample.sam")
        rec_with_tags = BioToolkit.BamRecord("tagged_read", "chr1", 50, [BioToolkit.BamCigarOp(4, 'M')], "ACGT"; quality="IIII", tags=Dict("AS"=>Int32(42), "FZ"=>Int32[1, 2, 3], "MD"=>"4M"))
        sam_file = BioToolkit.BamFile([rec_with_tags]; header=header)
        BioToolkit.write_sam(sam_path, sam_file)

        sam_content = read(sam_path, String)
        @test contains(sam_content, "AS:i:42")
        @test contains(sam_content, "FZ:B:i,1,2,3")

        sam_roundtrip = BioToolkit.read_sam(sam_path)
        rt_rec = sam_roundtrip.records[1]
        @test rt_rec.tags["AS"] == Int32(42)
        @test rt_rec.tags["FZ"] == Int32[1, 2, 3]
        @test rt_rec.tags["MD"] == "4M"

        invalid_rec = BioToolkit.BamRecord("bad_read", "chrUnknown", 10, [BioToolkit.BamCigarOp(4, 'M')], "ACGT")
        bad_bam = BioToolkit.BamFile([invalid_rec]; header=header)
        @test_throws ArgumentError BioToolkit.write_bam(joinpath(dir, "bad.bam"), bad_bam)

        dup_rec = BioToolkit.BamRecord("dup_read", "chr1", 9, [BioToolkit.BamCigarOp(10, 'M')], "ACGTACGTAA"; flag=BioToolkit.SAM_FLAG_DUPLICATE)
        cov_all = BioToolkit.bam_coverage([records[1], dup_rec, unmapped_rec], header.references; exclude_flags=0)
        @test cov_all["chr1"][10] == 2
        cov_filtered = BioToolkit.bam_coverage([records[1], dup_rec, unmapped_rec], header.references)
        @test cov_filtered["chr1"][10] == 1

        @test_throws ArgumentError BioToolkit.read_bam(joinpath(dir, "unsupported.cram"))
        @test_throws ArgumentError BioToolkit.write_bam(joinpath(dir, "unsupported.cram"), bam)
    end
end