using Test
using BioToolkit

@testset "Long-read genomics and SV calling" begin
    @testset "Read wrappers and minimizers" begin
        fastq = FastqRecord(DNASeq("ACGTTGCAACGT"), "IIIIIIIIIIII"; identifier="read1")
        nanopore = NanoporeRead(fastq; signal=fill(0.25f0, 16), channel_info=Dict("channel" => 42))
        pacbio = PacBioRead(fastq; signal=fill(0.5f0, 16), channel_info=Dict("movie" => "m540"))

        @test String(nanopore.sequence) == "ACGTTGCAACGT"
        @test nanopore.channel_info["channel"] == 42
        @test length(pacbio.signal) == 16

        hash_a = canonical_kmer_hash("ACGTAC")
        hash_b = canonical_kmer_hash(String(reverse_complement("ACGTAC")))
        @test hash_a == hash_b

        sketch = minimizer_sketch("TTGACCTGATCGTAGCTAGGATCCTAGATCGATG"; k=5, w=4)
        @test !isempty(sketch)
        @test all(seed.position >= 1 for seed in sketch)
        @test all(seed.hash isa UInt64 for seed in sketch)

        reference = Dict("chr1" => "TTGACCTGATCGTAGCTAGGATCCTAGATCGATG")
        index = build_minimizer_index(reference; k=5, w=4)
        candidates = find_candidate_regions("GATCGTAGCTAGGATC" , index; k=5, w=4)
        @test !isempty(candidates)
        @test first(candidates).chrom == "chr1"
        @test 6 <= first(candidates).start <= 16
    end

    @testset "Banded semiglobal alignment" begin
        query = "GATCGTAGCTA"
        reference = "TTTTTGATCGTAGCTAGGG"
        aln = banded_semiglobal_alignment(query, reference; bandwidth=6)

        @test aln.score > 0
        @test aln.query_start == 1
        @test aln.query_end == length(query)
        @test 1 <= aln.ref_start <= aln.ref_end <= length(reference)
        @test length(aln.aligned_query) == length(aln.aligned_reference)
        @test !isempty(aln.cigar)

        query_with_insertion = "GATCGTTAGCTA"
        aln2 = banded_semiglobal_alignment(query_with_insertion, reference; bandwidth=8)
        @test length(aln2.cigar) >= 1
        @test any(op.op in ('I', 'D') for op in aln2.cigar)

        long_query = "ACGT"^25
        long_reference = "T"^120 * long_query * "G"^40
        aln3 = banded_semiglobal_alignment(long_query, long_reference)
        @test replace(aln3.aligned_query, "-" => "") == long_query
        @test replace(aln3.aligned_reference, "-" => "") == long_query
    end

    @testset "SV evidence extraction and clustering" begin
        split_record = BamRecord(
            "split1",
            "chr1",
            100,
            [BamCigarOp(20, 'M'), BamCigarOp(60, 'D'), BamCigarOp(20, 'M')],
            "A"^40;
            mapq=60,
            mate_refname="chr1",
            mate_pos=260,
            template_length=400,
        )
        split_events = split_read_evidence(split_record; min_sv_size=50)
        @test any(event.sv_type == :DEL for event in split_events)

        trans_record = BamRecord(
            "pair1",
            "chr1",
            200,
            [BamCigarOp(50, 'M')],
            "C"^50;
            flag=0,
            mate_refname="chr2",
            mate_pos=1200,
            template_length=0,
        )
        trans_events = discordant_pair_evidence(trans_record)
        @test any(event.sv_type == :TRA for event in trans_events)

        inv_record = BamRecord(
            "pair2",
            "chr1",
            300,
            [BamCigarOp(50, 'M')],
            "G"^50;
            flag=0x30,
            mate_refname="chr1",
            mate_pos=450,
            template_length=900,
        )
        inv_events = discordant_pair_evidence(inv_record; insert_size_mean=500.0, insert_size_std=100.0, z_threshold=3.0)
        @test any(event.sv_type == :INV for event in inv_events)
        @test any(event.sv_type == :DEL for event in inv_events)

        record_a = BamRecord(
            "cluster_a",
            "chr1",
            100,
            [BamCigarOp(25, 'M'), BamCigarOp(70, 'D'), BamCigarOp(25, 'M')],
            "A"^50;
            mate_refname="chr1",
            mate_pos=260,
            template_length=420,
        )
        record_b = BamRecord(
            "cluster_b",
            "chr1",
            120,
            [BamCigarOp(25, 'M'), BamCigarOp(68, 'D'), BamCigarOp(25, 'M')],
            "A"^50;
            mate_refname="chr1",
            mate_pos=280,
            template_length=430,
        )

        calls = call_structural_variants([record_a, record_b]; min_sv_size=50, cluster_window=80, min_support=2)
        @test !isempty(calls)
        @test any(call.sv_type == :DEL for call in calls)
        del_call = first(filter(call -> call.sv_type == :DEL, calls))
        @test length(del_call.supporting_reads) >= 2
        @test del_call.genotype in ("0/1", "1/1")
    end
end
