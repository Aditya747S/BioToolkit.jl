using Test
using BioToolkit
import BioToolkit: BamRecord, BamHeader, BamFile, BamCigarOp, BamReference
import BioToolkit: read_bam, write_bam, read_sam, write_sam
import BioToolkit: leftposition, rightposition, alignlength, readname, cigar_rle, ismapped
import BioToolkit: is_paired, is_proper_pair, is_unmapped, is_mate_unmapped
import BioToolkit: is_reverse_strand, is_mate_reverse_strand, is_read1, is_read2
import BioToolkit: is_secondary, is_qc_failed, is_duplicate, is_supplementary
import BioToolkit: SAM_FLAG_PAIRED, SAM_FLAG_PROPER_PAIR, SAM_FLAG_UNMAPPED
import BioToolkit: SAM_FLAG_MATE_UNMAPPED, SAM_FLAG_REVERSE, SAM_FLAG_MATE_REVERSE
import BioToolkit: SAM_FLAG_READ1, SAM_FLAG_READ2, SAM_FLAG_SECONDARY
import BioToolkit: SAM_FLAG_QC_FAIL, SAM_FLAG_DUPLICATE, SAM_FLAG_SUPPLEMENTARY

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

        long_cigar_ops = [BioToolkit.BamCigarOp(1, 'M') for _ in 1:65537]
        long_cigar_rec = BioToolkit.BamRecord("long_read", "chr1", 1, long_cigar_ops, repeat("A", 65537))
        long_bam = BioToolkit.BamFile([long_cigar_rec]; header=header)
        long_bam_path = joinpath(dir, "long_cigar.bam")
        BioToolkit.write_bam(long_bam_path, long_bam)
        long_rt = BioToolkit.read_bam(long_bam_path)
        @test length(long_rt.records[1].cigar) == 65537
        @test BioToolkit.leftposition(records[1]) == 10
        @test BioToolkit.rightposition(records[1]) == 19
        @test BioToolkit.alignlength(records[1]) == 10
        @test BioToolkit.readname(records[1]) == "read1"
        @test BioToolkit.ismapped(records[1]) == true
        @test BioToolkit.ismapped(unmapped_rec) == false
        @test BioToolkit.cigar_rle(records[1]) == (['M'], [10])

        @test_throws ArgumentError BioToolkit.read_bam(joinpath(dir, "unsupported.cram"))
        @test_throws ArgumentError BioToolkit.write_bam(joinpath(dir, "unsupported.cram"), bam)
    end

    @testset "XAM.jl - Auxiliary Data Tests" begin
        rec_empty = BamRecord("read_empty", "chr1", 100, [BamCigarOp(10, 'M')], BioToolkit.BioSequence{BioToolkit.DNAAlphabet}("ACGTACGTAA"); tags=Dict{String,Any}())
        @test isempty(rec_empty.tags)

        rich_tags = Dict{String,Any}(
            "AS" => Int8(-18),
            "NM" => Int16(1),
            "XA" => Float32(3.14),
            "XB" => "some text",
            "XC" => Int32[10, -5, 8],
            "FZ" => UInt8[0x01, 0x02, 0x03],
            "MD" => "10M"
        )
        rec_rich = BamRecord("read_rich", "chr1", 100, [BamCigarOp(10, 'M')], BioToolkit.BioSequence{BioToolkit.DNAAlphabet}("ACGTACGTAA"); tags=rich_tags)
        @test length(rec_rich.tags) == 7
        @test rec_rich.tags["AS"] === Int8(-18)
        @test rec_rich.tags["NM"] === Int16(1)
        @test rec_rich.tags["XA"] === Float32(3.14)
        @test rec_rich.tags["XB"] == "some text"
        @test rec_rich.tags["XC"] == Int32[10, -5, 8]
        @test rec_rich.tags["FZ"] == UInt8[0x01, 0x02, 0x03]

        mktempdir() do dir
            hdr = BamHeader([BamReference("chr1", 1000)])
            bam_in = BamFile([rec_rich]; header=hdr)
            bam_path = joinpath(dir, "rich_tags.bam")
            write_bam(bam_path, bam_in)
            bam_out = read_bam(bam_path)
            rt_rec = bam_out.records[1]
            @test rt_rec.tags["AS"] == Int8(-18)
            @test rt_rec.tags["NM"] == Int16(1)
            @test rt_rec.tags["XA"] == Float32(3.14)
            @test rt_rec.tags["XB"] == "some text"
            @test rt_rec.tags["XC"] == Int32[10, -5, 8]
            @test rt_rec.tags["FZ"] == UInt8[0x01, 0x02, 0x03]
        end
    end

    @testset "XAM.jl - SAM Flags Suite" begin
        flag_tests = [
            (is_paired, SAM_FLAG_PAIRED, true),
            (is_proper_pair, SAM_FLAG_PROPER_PAIR, true),
            (is_unmapped, SAM_FLAG_UNMAPPED, true),
            (ismapped, SAM_FLAG_UNMAPPED, false),
            (is_mate_unmapped, SAM_FLAG_MATE_UNMAPPED, true),
            (is_reverse_strand, SAM_FLAG_REVERSE, true),
            (is_mate_reverse_strand, SAM_FLAG_MATE_REVERSE, true),
            (is_read1, SAM_FLAG_READ1, true),
            (is_read2, SAM_FLAG_READ2, true),
            (is_secondary, SAM_FLAG_SECONDARY, true),
            (is_qc_failed, SAM_FLAG_QC_FAIL, true),
            (is_duplicate, SAM_FLAG_DUPLICATE, true),
            (is_supplementary, SAM_FLAG_SUPPLEMENTARY, true),
        ]

        for (func, flag_bit, expected_when_set) in flag_tests
            rec_set = BamRecord("test_flag", "chr1", 1, [BamCigarOp(4, 'M')], BioToolkit.BioSequence{BioToolkit.DNAAlphabet}("ACGT"); flag=flag_bit)
            @test func(rec_set) == expected_when_set

            rec_unset = BamRecord("test_flag", "chr1", 1, [BamCigarOp(4, 'M')], BioToolkit.BioSequence{BioToolkit.DNAAlphabet}("ACGT"); flag=UInt16(0))
            @test func(rec_unset) == !expected_when_set
        end
    end

    @testset "XAM.jl - Record Accessors & Alignment Interval Tests" begin
        ops = [BamCigarOp(8, 'M'), BamCigarOp(2, 'I'), BamCigarOp(4, 'M'), BamCigarOp(1, 'D'), BamCigarOp(3, 'M')]
        seq = BioToolkit.BioSequence{BioToolkit.DNAAlphabet}("TTAGATAAAGGATACTG")
        rec = BamRecord("r001", "chr1", 6, ops, seq; flag=UInt16(99), mapq=30)

        @test readname(rec) == "r001"
        @test ismapped(rec) == true
        @test rec.pos == 6
        @test leftposition(rec) == 7
        @test alignlength(rec) == 16  # 8 + 4 + 1 + 3 = 16
        @test rightposition(rec) == 22 # 7 + 16 - 1 = 22
        @test cigar_rle(rec) == (['M', 'I', 'M', 'D', 'M'], [8, 2, 4, 1, 3])
    end

    @testset "XAM.jl - Crosscheck SAM vs BAM Records" begin
        mktempdir() do dir
            hdr = BamHeader([BamReference("chr1", 10000), BamReference("chr2", 5000)])
            cigar1 = [BamCigarOp(27, 'M'), BamCigarOp(1, 'D'), BamCigarOp(73, 'M')]
            seq1 = BioToolkit.BioSequence{BioToolkit.DNAAlphabet}(repeat("A", 100))
            rec1 = BamRecord("read_cross_1", "chr1", 1, cigar1, seq1; mapq=60, quality=repeat("#", 100), tags=Dict{String,Any}("AS" => Int32(42), "MD" => "27M1D73M"))

            cigar2 = [BamCigarOp(50, 'M')]
            seq2 = BioToolkit.BioSequence{BioToolkit.DNAAlphabet}(repeat("C", 50))
            rec2 = BamRecord("read_cross_2", "chr2", 500, cigar2, seq2; mapq=40, quality=repeat("I", 50), tags=Dict{String,Any}("NM" => Int32(0)))

            bam_file = BamFile([rec1, rec2]; header=hdr)

            bam_path = joinpath(dir, "cross.bam")
            sam_path = joinpath(dir, "cross.sam")
            write_bam(bam_path, bam_file)
            write_sam(sam_path, bam_file)

            bam_read = read_bam(bam_path)
            sam_read = read_sam(sam_path)

            @test length(bam_read.records) == length(sam_read.records) == 2

            for (b_rec, s_rec) in zip(bam_read.records, sam_read.records)
                @test readname(b_rec) == readname(s_rec)
                @test b_rec.refname == s_rec.refname
                @test leftposition(b_rec) == leftposition(s_rec)
                @test rightposition(b_rec) == rightposition(s_rec)
                @test alignlength(b_rec) == alignlength(s_rec)
                @test b_rec.mapq == s_rec.mapq
                @test b_rec.flag == s_rec.flag
                @test String(b_rec.sequence) == String(s_rec.sequence)
                @test b_rec.quality == s_rec.quality
                @test b_rec.cigar == s_rec.cigar
            end
        end
    end

    @testset "XAM.jl - Convert SAM to BAM Roundtrip" begin
        mktempdir() do dir
            sam_str = "seq1\t81\tchr1\t1051\t60\t70M\t=\t1821\t702\tTCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGAC\t*\tNM:i:0\tMD:Z:70\tMC:Z:70M\tAS:i:70\tXS:i:0\n"
            sam_io = IOBuffer(sam_str)
            sam_file = read_sam(sam_io)
            @test length(sam_file.records) == 1
            rec = sam_file.records[1]
            @test rec.qname == "seq1"
            @test rec.flag == UInt16(81)
            @test rec.refname == "chr1"
            @test leftposition(rec) == 1051

            bam_path = joinpath(dir, "converted.bam")
            write_bam(bam_path, sam_file)

            bam_read = read_bam(bam_path)
            @test length(bam_read.records) == 1
            b_rec = bam_read.records[1]
            @test b_rec.qname == rec.qname
            @test b_rec.flag == rec.flag
            @test b_rec.refname == rec.refname
            @test leftposition(b_rec) == leftposition(rec)
            @test alignlength(b_rec) == alignlength(rec)
            @test String(b_rec.sequence) == String(rec.sequence)
        end
    end

    @testset "XAM.jl - Read Ultra Long CIGAR (>65535 ops)" begin
        long_ops = [BamCigarOp(1, 'M') for _ in 1:72000]
        long_seq = BioToolkit.BioSequence{BioToolkit.DNAAlphabet}(repeat("A", 72000))
        rec_long = BamRecord("long_read_xam", "chr1", 0, long_ops, long_seq)

        hdr = BamHeader([BamReference("chr1", 100000)])
        bam_in = BamFile([rec_long]; header=hdr)

        mktempdir() do dir
            bam_path = joinpath(dir, "long_xam.bam")
            write_bam(bam_path, bam_in)

            bam_out = read_bam(bam_path)
            rt_long = bam_out.records[1]

            @test length(rt_long.cigar) == 72000
            @test rt_long.cigar == long_ops
            @test alignlength(rt_long) == 72000
            @test rightposition(rt_long) == 72000
        end
    end
end