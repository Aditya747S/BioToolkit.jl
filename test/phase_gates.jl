# ==========================================================================
# Phase gates 
#
# G1: the package loads at all (implicit — this file only runs after include).
# G2: a fast smoke pass over core entry points that must keep working across
#     correction phases. This file is intentionally lightweight: no CUDA, no
#     plotting, no network. Detailed correctness tests live in the per-domain
#     test files; this catches phase-to-phase regressions.
# ==========================================================================
using Test
using BioToolkit

@testset "Phase gate: core sequence surface" begin
    seq = BioToolkit.PackedDNASeq("ACGTAGCGTA")
    @test BioToolkit.gc_content(seq) ≈ 0.5
    @test String(BioToolkit.reverse_complement(seq)) == "TACGCTACGT"
    @test String(BioToolkit.translate_dna("ATGAAATAA")) == "MK*"
    @test BioToolkit.melting_temp("ACGT") == 12.0
    kmers = BioToolkit.kmer_frequency(BioToolkit.DNASeq("ACGTACGT"), 2)
    @test kmers["AC"] == 2
end

@testset "Phase gate: parsing round trips" begin
    mktempdir() do dir
        fasta_path = joinpath(dir, "t.fasta")
        rec1 = BioToolkit.SeqRecord(BioToolkit.DNASeq("ATGC"); identifier="s1")
        rec2 = BioToolkit.SeqRecord(BioToolkit.DNASeq("GGTT"); identifier="s2")
        BioToolkit.write_fasta(fasta_path, [rec1, rec2])
        records = BioToolkit.read_fasta(fasta_path)
        @test [r.identifier for r in records] == ["s1", "s2"]
        @test String(records[1].sequence) == "ATGC"

        fastq_path = joinpath(dir, "t.fastq")
        fq1 = BioToolkit.FastqRecord(BioToolkit.DNASeq("ACGT"), "IIII"; identifier="r1")
        BioToolkit.write_fastq(fastq_path, [fq1])
        fq = BioToolkit.read_fastq(fastq_path)
        @test String(fq[1].sequence) == "ACGT"
        @test BioToolkit.phred_scores(fq[1].quality) == [40, 40, 40, 40]

        vcf_path = joinpath(dir, "t.vcf")
        variant = BioToolkit.VariantTextRecord("chr1", 100, "rs1", "A", "G", 50.0f0)
        BioToolkit.write_vcf(vcf_path, [variant])
        vcf = BioToolkit.read_vcf(vcf_path)
        @test length(vcf) == 1
        @test vcf[1].chrom == "chr1"
    end
end

@testset "Phase gate: intervals and alignment" begin
    intervals = [BioToolkit.GenomicInterval("chr1", 10, 20), BioToolkit.GenomicInterval("chr1", 15, 25)]
    collection = BioToolkit.build_collection(intervals)
    @test length(BioToolkit.find_overlaps(BioToolkit.GenomicInterval("chr1", 18, 22), collection)) == 2

    aln = BioToolkit.pairwise_align("ACGTACGT", "ACGTTCGT")
    @test aln.score > 0
    @test BioToolkit.melting_temp("ACGCGTACGCGT"; method=:nearest_neighbor) > 0
end

@testset "Phase gate: enrichment and DE entry points" begin
    terms = [
        BioToolkit.EnrichmentTerm("GO:1", "pathway one", "BP", ["g1", "g2", "g3", "g4"], String[]),
        BioToolkit.EnrichmentTerm("GO:2", "pathway two", "BP", ["g5", "g6", "g7", "g8"], String[]),
    ]
    db = BioToolkit.build_annotation_database(terms)
    result = BioToolkit.enrichment_test(["g1", "g2"], db)
    @test result[1].pvalue <= 1.0

    ranked = [("g1", 3.2), ("g2", 2.1), ("g5", 1.4), ("g6", 1.1), ("g7", 0.6), ("g3", 0.2), ("g4", -1.5), ("g8", -2.0)]
    gsea = BioToolkit.fgsea_like(ranked, Dict("GO:1" => ["g1", "g2", "g3", "g4"]); n_permutations=200, min_size=2, max_size=500)
    @test all(r -> r.pvalue <= 1.0, gsea)
end
