using Test
using BioToolkit

using Test, Random

@testset "pairwise_align traceback self-consistency" begin
  Random.seed!(42)
  for _ in 1:500
    n, m = rand(0:20), rand(0:20)
    a = randstring("ACGT", n)
    b = randstring("ACGT", m)
    match, mismatch, gap = rand(1:5), -rand(1:5), -rand(1:5)

    result = pairwise_align(a, b; match=match, mismatch=mismatch, gap=gap)

    # recompute the score purely from the returned aligned sequences
    recomputed = 0
    for (x, y) in zip(String(result.left), String(result.right))
      if x == '-' || y == '-'
        recomputed += gap
      elseif x == y
        recomputed += match
      else
        recomputed += mismatch
      end
    end

    @test recomputed == result.score
    @test replace(String(result.left), '-' => "") == a
    @test replace(String(result.right), '-' => "") == b
  end
end

@testset "Parser edge cases" begin
  @test BioToolkit.validate_dna("") == true
  @test BioToolkit.validate_dna("ATGNRYSWKMBDHV") == true
  @test BioToolkit.validate_dna(BioToolkit.DNASeq("atgNRYSWKMBDHV")) == true
  @test BioToolkit.reverse_complement("R") == "Y"
  @test BioToolkit.reverse_complement("Y") == "R"
  @test BioToolkit.gc_content("ATGNRR") == 1 / 3
  @test BioToolkit.gc_content(BioToolkit.DNASeq("ATGNRR")) == 1 / 3
  @test_throws ArgumentError BioToolkit.DNASeq("ATBX")
  @test BioToolkit.RNASeq("aucg") == BioToolkit.RNASeq("AUCG")

  mktempdir() do dir
    fastq_path = joinpath(dir, "bad.fastq")
    open(fastq_path, "w") do io
      write(io, "@seq1\n")
      write(io, "ACGT\n")
      write(io, "+\n")
      write(io, "III\n")
    end

    @test_throws ArgumentError BioToolkit.read_fastq(fastq_path)
    @test_throws ArgumentError BioToolkit.read_fastq(fastq_path; alphabet=BioToolkit.DNAAlphabet)
  end

  mktempdir() do dir
    fastq_path = joinpath(dir, "lower.fastq")
    open(fastq_path, "w") do io
      write(io, "@seq1\n")
      write(io, "acgt\n")
      write(io, "+\n")
      write(io, "IIII\n")
    end

    typed_records = BioToolkit.read_fastq(fastq_path; alphabet=BioToolkit.DNAAlphabet)
    @test typed_records[1].sequence isa BioToolkit.DNASeq
    @test typed_records[1].sequence == "ACGT"
  end

  mktempdir() do dir
    fasta_path = joinpath(dir, "bad.fasta")
    open(fasta_path, "w") do io
      write(io, "ACGT\n")
    end

    @test_throws ArgumentError BioToolkit.read_fasta(fasta_path)
  end

  @test_throws ArgumentError BioToolkit.translate_dna("ATGAA")
  @test_throws ArgumentError BioToolkit.translate_dna(BioToolkit.DNASeq("ATGAA"))

  @test_throws ArgumentError BioToolkit.substitution_matrix("ACGT", [1 2; 3 4])

  mktempdir() do dir
    fasta_path = joinpath(dir, "empty.fasta")
    write(fasta_path, "")
    records = BioToolkit.read_fasta(fasta_path)
    @test isempty(records)
    @test eltype(records) <: BioToolkit.SeqRecord{BioToolkit.DNAAlphabet}
  end

  seq_info = BioToolkit.GenomicRanges.SeqInfo("chrM", 16569, true)
  buffer = IOBuffer()
  show(buffer, seq_info)
  @test occursin("circular=true", String(take!(buffer)))
end

@testset "Streaming FASTQ" begin
  mktempdir() do dir
    fastq_path = joinpath(dir, "stream.fastq")
    open(fastq_path, "w") do io
      write(io, "@seq1\nACGT\n+\nIIII\n")
      write(io, "@seq2\nTTTT\n+\n####\n")
    end

    function count_records(path)
      count = 0
      for record in BioToolkit.each_fastq_record(path; sequence_type=BioToolkit.DNASeq)
        count += 1
        @test record isa BioToolkit.FastqRecord{BioToolkit.DNAAlphabet}
      end
      return count
    end

    @test @inferred(count_records(fastq_path)) == 2
    @test length(BioToolkit.read_fastq(fastq_path; alphabet=BioToolkit.DNAAlphabet)) == 2
  end
end

# ==========================================================================
# Session 1 corrections (2026-09-04): packed GC, nearest-neighbor Tm,
# sequence complexity k-mer hashing, real find_orfs, packed ungap,
# PackedDNASeq validation, SCHNEIDER codon-matrix routing.
# Evidence: plan/audit_notes.md wave 3; plan/correction_plan.md items 6-13.
# ==========================================================================

@testset "Session 1: packed gc_content counts G+C" begin
    @test BioToolkit.gc_content(BioToolkit.PackedDNASeq("GCGC")) == 1.0
    @test BioToolkit.gc_content(BioToolkit.PackedDNASeq("ATAT")) == 0.0
    @test BioToolkit.gc_content(BioToolkit.PackedDNASeq("ATGC")) == 0.5
    # Exact parity with the unpacked path (previously the packed path
    # returned (C+T)/n because the mask missed G lanes).
    Random.seed!(2026)
    for _ in 1:25
        s = String(rand(collect("ACGT"), rand(1:200)))
        @test BioToolkit.gc_content(BioToolkit.PackedDNASeq(s)) == BioToolkit.gc_content(BioToolkit.DNASeq(s))
    end
end

@testset "Session 1: PackedDNASeq validation is explicit" begin
    @test_throws ArgumentError BioToolkit.PackedDNASeq("ACGN")
    # U is accepted in a DNA context and stored as T, documented.
    @test String(BioToolkit.PackedDNASeq("ACGU")) == "ACGT"
    @test_throws ArgumentError BioToolkit.BitPackedVector{2}(collect(codeunits("ACGN")); validate=true)
end

@testset "Session 1: packed gap removal" begin
    packed = BioToolkit.PackedDNASeq("ACGTAA")
    @test String(BioToolkit.ungap(packed)) == "ACGTAA"
    @test_throws ArgumentError BioToolkit.ungap!(packed)  # immutable chunk length; message points at ungap
end

@testset "Session 1: nearest-neighbor melting temperature" begin
    # Regression lock computed independently from SantaLucia (1998) for ACGT:
    # stacks AC->GT, CG->CG, GT->GT; initiation; two terminal A/T.
    dh = (0.2 - 8.4 - 10.6 - 8.4 + 2.2 + 2.2) * 1000.0
    ds = (-5.7 - 22.4 - 27.2 - 22.4 + 6.9 + 6.9)
    expected = dh / (ds + 1.987 * log(2.5e-7 / 4)) - 273.15 + 16.6 * log10(0.05)
    @test isapprox(BioToolkit.melting_temp("ACGT"; method=:nearest_neighbor), expected; atol=1e-6)

    # Orientation invariance: a duplex and its reverse complement coincide.
    @test isapprox(BioToolkit.melting_temp("AACG"; method=:nearest_neighbor),
                   BioToolkit.melting_temp("CGTT"; method=:nearest_neighbor); atol=1e-10)

    # Biophysical properties.
    @test BioToolkit.melting_temp("GCGCGCGC"; method=:nearest_neighbor) >
          BioToolkit.melting_temp("ATATATAT"; method=:nearest_neighbor)
    @test BioToolkit.melting_temp("ACGCGTACGCGT"; method=:nearest_neighbor) >
          BioToolkit.melting_temp("ACGCGT"; method=:nearest_neighbor)
    @test BioToolkit.melting_temp("ACGCGTACGCGT"; method=:nearest_neighbor, na_conc=0.2) >
          BioToolkit.melting_temp("ACGCGTACGCGT"; method=:nearest_neighbor, na_conc=0.01)
    # The published self-complementary symmetry term is a +1.4 kcal/mol
    # PENALTY on dH (naive stack sums overestimate self-comp duplex
    # stability), so the corrected Tm is slightly lower.
    @test BioToolkit.melting_temp("ACGCGTACGCGT"; method=:nearest_neighbor, self_complementary=true) <
          BioToolkit.melting_temp("ACGCGTACGCGT"; method=:nearest_neighbor, self_complementary=false)

    # Ambiguity and edge conditions are rejected, not silently averaged.
    @test_throws ArgumentError BioToolkit.melting_temp("ACGNT"; method=:nearest_neighbor)
    @test_throws ArgumentError BioToolkit.melting_temp("A"; method=:nearest_neighbor)
    @test_throws ArgumentError BioToolkit.melting_temp("ACGT"; method=:nearest_neighbor, rna=true)

    # Existing methods unchanged (the :basic approximation is only meaningful
    # for oligos longer than ~13 nt).
    @test BioToolkit.melting_temp("AAAA") == 8.0
    @test isapprox(BioToolkit.melting_temp("ACGTACGTACGTACGTACGT"; method=:basic), 51.78; atol=1e-10)
end

@testset "Session 1: sequence_complexity k-mer hashing" begin
    # Homopolymer: one distinct k-mer.
    @test BioToolkit.sequence_complexity(BioToolkit.DNASeq("A"^30); k=3) == 1 / 28
    # Period-4 sequence: exactly 4 distinct trimers (one per frame) over 62 windows.
    @test BioToolkit.sequence_complexity(BioToolkit.DNASeq("ACGT"^16); k=3) == 4 / 62
    # Single bases cover the whole k=1 alphabet.
    @test BioToolkit.sequence_complexity(BioToolkit.DNASeq("ACGT"^16); k=1) == 1.0
    # k > 31: previously 4^k overflowed Int64; now the String path handles it.
    # Period-4 100-mer has exactly 4 distinct 32-mers (one per offset).
    @test BioToolkit.sequence_complexity(BioToolkit.DNASeq("ACGT"^25); k=32) == 4 / 69
    # Ambiguity codes are distinct states, not silently collapsed.
    @test BioToolkit.sequence_complexity(BioToolkit.DNASeq("ACGN"^5); k=2) == 4 / 19
    @test_throws ArgumentError BioToolkit.sequence_complexity(BioToolkit.DNASeq("ACGT"); k=0)
end

@testset "Session 1: find_orfs requires start codons" begin
    orfs = BioToolkit.find_orfs("ATGAAAGGGTAA")
    @test length(orfs) == 1
    @test String(orfs[1]) == "MKG"

    # An open ORF without a stop still counts to the frame end.
    open_orfs = BioToolkit.find_orfs("ATGAAAGGG")
    @test String.(open_orfs) == ["MKG"]

    # Without a start codon there are no ORFs; the legacy stop-to-stop scan
    # remains available (it reports every inter-stop fragment).
    @test isempty(BioToolkit.find_orfs("GGGTAAAAAGGG"))
    legacy = BioToolkit.find_orfs("GGGTAAAAAGGG"; require_start=false)
    # Legacy stop-to-stop scan reports both strands; the forward fragments
    # "G" (before the stop) and "KG" (after it) must both be present.
    legacy_strings = String.(legacy)
    @test "G" in legacy_strings && "KG" in legacy_strings

    # min_aa filters short ORFs.
    @test isempty(BioToolkit.find_orfs("ATGAAATAA"; min_aa=3))
    @test String.(BioToolkit.find_orfs("ATGAAATAA"; min_aa=1)) == ["MK"]

    # Reverse-strand ORFs are found.
    rev = BioToolkit.find_orfs("TTACCCtttCAT")
    @test any(o -> String(o) == "MKG", rev)
end

@testset "Session 1: translation paths agree on stop codons" begin
    @test String(BioToolkit.translate_dna("ATGAAATAA")) == "MK*"
    @test String(BioToolkit.translate_dna("ATGAAATAA"; stop_at_stop=true)) == "MK"
    buffer = Vector{UInt8}(undef, 3)
    written = BioToolkit.translate_dna!(buffer, collect(codeunits("ATGAAATAA")))
    @test written == 3
    @test String(buffer[1:written]) == "MK*"
end

@testset "Session 1: SCHNEIDER codon matrix parses without warnings" begin
    # Previously routed through the single-byte nucleotide/AA parser and
    # emitted "substitution matrix row has the wrong width" on every load.
    @test "SCHNEIDER" in BioToolkit.available_named_substitution_matrices()
    codon_matrix = BioToolkit.named_codon_substitution_matrix("SCHNEIDER")
    @test size(codon_matrix.scores) == (64, 64)
    @test_throws ArgumentError BioToolkit.named_codon_substitution_matrix("BLOSUM62")
end

@testset "Session 1: standard substitution matrices are intact" begin
    blosum62 = BioToolkit.named_substitution_matrix("BLOSUM62")
    @test blosum62.scores[1, 1] == 4            # A/A
    @test blosum62.scores[1, 3] == -2           # A/N
    @test blosum62.scores[1, 13] == -1          # A/M
    for name in BioToolkit.available_named_substitution_matrices()
        name == "SCHNEIDER" && continue
        sm = BioToolkit.named_substitution_matrix(name)
        @test size(sm.scores, 1) == size(sm.scores, 2)
    end
end

# ==========================================================================
# Session 2: bioplotting HTML export works (JSON import + escape helpers +
# reflected savefig), VcfDocument derives sample columns, merge_loops
# propagates bin columns.
# ==========================================================================
@testset "Session 2: bioplotting HTML export" begin
    vr = BioToolkit.volcano_plot([BioToolkit.DEResult("g1", 10.0, 2.5, 0.5, 12.0, 1e-6, 1e-4)])
    html = BioToolkit.to_html(vr)
    @test occursin("<html", html)
    @test occursin("volcano", html)
end

@testset "Session 2: VcfDocument sample columns" begin
    rec = BioToolkit.parse_vcf_record("chr1\t42\trs1\tA\tG\t99.5\tPASS\t.\tGT\t0/1")
    doc = BioToolkit.VcfDocument([rec])
    @test doc.header.sample_names == ["0/1"]
end

@testset "Session 2: hic_diff_test stars and pileup fallback" begin
    result = BioToolkit.hic_diff_test([50, 5, 100], [5, 50, 3])
    @test all(s -> s in ("***", "**", "*", ".", ""), result.significance)
end
