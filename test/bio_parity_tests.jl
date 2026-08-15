using Test
using BioToolkit

@testset "Sequence Analysis & Alignment Utilities" begin

    @testset "CIGAR Engine & Count Helpers" begin
        seq1 = DNASeq("ACGTACGTAC")
        seq2 = DNASeq("AC--ACGTAC")
        aln = needleman_wunsch(seq1, seq2)
        
        c = cigar(aln)
        @test !isempty(c)
        @test occursin("M", c) || occursin("I", c) || occursin("D", c)
        
        parsed = parse_cigar("10M2I5M3D")
        @test length(parsed) == 4
        @test parsed[1] == ('M', 10)
        @test parsed[2] == ('I', 2)
        @test parsed[3] == ('M', 5)
        @test parsed[4] == ('D', 3)
        
        @test count_matches(aln) >= 0
        @test count_mismatches(aln) >= 0
        @test count_insertions(aln) >= 0
        @test count_deletions(aln) >= 0
        @test count_aligned(aln) == length(aln.left)
    end

    @testset "Coordinate Mapping Infrastructure" begin
        # Perfect alignment
        seq1 = DNASeq("ACGTACGT")
        seq2 = DNASeq("ACGTACGT")
        aln = needleman_wunsch(seq1, seq2)
        
        @test seq2ref(aln, 1) == 1
        @test seq2ref(aln, 4) == 4
        @test ref2seq(aln, 3) == 3
        @test seq2aln(aln, 5) == 5
        @test ref2aln(aln, 5) == 5
        @test aln2seq(aln, 5) == 5
        @test aln2ref(aln, 5) == 5
        
        # Alignment with gap
        seq_gap1 = DNASeq("ACGTACGT")
        seq_gap2 = DNASeq("AC--ACGT")
        aln_gap = needleman_wunsch(seq_gap1, seq_gap2)
        @test seq2aln(aln_gap, 1) == 1
        @test aln2seq(aln_gap, 1) == 1
    end

    @testset "Banded, Semi-Global, Overlap & Unified pairalign" begin
        s1 = DNASeq("ACGTACGTACGTACGT")
        s2 = DNASeq("ACGTACCTACGTACGT")
        
        # Banded NW
        res_banded = banded_needleman_wunsch(s1, s2; k=5)
        @test res_banded.score > 0
        
        # Semi-Global & Overlap
        res_sg = semi_global_align(s1, s2)
        @test res_sg.score > 0
        
        res_ov = overlap_align(s1, s2)
        @test res_ov.score > 0
        
        # Unified pairalign dispatcher
        @test pairalign(:global, s1, s2).score == needleman_wunsch(s1, s2).score
        @test pairalign(:local, s1, s2).score == smith_waterman(s1, s2).score
        @test pairalign(:banded, s1, s2; k=5).score > 0
        @test pairalign(:semi_global, s1, s2).score > 0
        @test pairalign(:overlap, s1, s2).score > 0
    end

    @testset "FOLDSEEK3DI & GRANTHAM1974 Matrices" begin
        mat_foldseek = substitution_matrix("FOLDSEEK3DI")
        @test size(mat_foldseek.scores) == (20, 20)
        
        mat_grantham = substitution_matrix("GRANTHAM1974")
        @test size(mat_grantham.scores) == (20, 20)
        
        mat_blosum = substitution_matrix(:BLOSUM62)
        @test size(mat_blosum.scores) == (24, 24)
    end

    @testset "NCBI Genetic Translation Tables (1-33)" begin
        # Table 1: Standard
        tbl1 = ncbi_trans_table(1)
        @test tbl1["ATG"] == 'M'
        @test tbl1["TAA"] == '*'
        
        # Table 2: Vertebrate Mitochondrial (AGA, AGG -> *)
        tbl2 = ncbi_trans_table(2)
        @test tbl2["AGA"] == '*'
        @test tbl2["ATA"] == 'M'
        
        # Table 5: Invertebrate Mitochondrial (AGA, AGG -> S)
        tbl5 = ncbi_trans_table(5)
        @test tbl5["AGA"] == 'S'
        
        # Translation with table parameter
        dna = DNASeq("ATGAGA")
        @test String(translate_dna(dna; table=1)) == "MR"
        @test String(translate_dna(dna; table=2, stop_at_stop=true)) == "M"
    end

    @testset "Myers' Bit-Vector Approximate Search & BioRegex" begin
        # Approximate search
        matches = approximate_search("ACGT", "ACGCACGTACGA", 1; overlap=true)
        @test !isempty(matches)
        @test matches[1] isa UnitRange{Int}
        
        query = ApproximateSearchQuery("ACGT", 1)
        @test !isempty(approximate_search(query, "ACGCACGTACGA"))
        
        # BioRegex
        biore = biore"RGY" # R=[AG], G, Y=[CTU]
        @test occursin(biore, "AGCT")
        @test occursin(biore, "GGCT")
        @test !occursin(biore, "TTTT")
        
        found = findall(biore, "ACAGCTACGGCT")
        @test !isempty(found)
    end

    @testset "String Macros, In-place Transformations & Predicates" begin
        # String macros
        d = dna"ACGT"
        r = rna"ACGU"
        a = aa"MKST"
        @test d isa DNASeq
        @test r isa RNASeq
        @test a isa AASeq
        
        # In-place complement & reverse_complement!
        buf = Vector{UInt8}(undef, 4)
        complement!(buf, codeunits("ACGT"))
        @test String(buf) == "TGCA"
        
        # In-place mutating BioSequence functions
        seq = DNASeq("ACGT")
        complement!(seq)
        @test String(seq) == "TGCA"
        reverse_complement!(seq)
        @test String(seq) == "TGCA" # reverse of TGCA is ACGT, complement of ACGT is TGCA, reverse_complement of TGCA is TGCA
        
        seq2 = DNASeq("A-C-G-T")
        ungap!(seq2)
        @test String(seq2) == "ACGT"
        
        seq3 = DNASeq("TTTT")
        canonical!(seq3)
        @test String(seq3) == "AAAA"
    end

    @testset "Semi-Global, Overlap & Banded Alignment Modes" begin
        # Semi-Global Alignment (free end gaps)
        res_sg = semi_global_align("ACGTACGT", "ACGT")
        @test res_sg.score > 0
        @test res_sg.matches == 4

        # Overlap Alignment
        res_ov = overlap_align("ACGTACGT", "CGTACGTT")
        @test res_ov.score > 0
        @test res_ov.matches >= 6

        # Banded Needleman-Wunsch Alignment
        res_bnw = banded_needleman_wunsch("ACGTACGT", "ACGTACGT"; k=5)
        @test res_bnw.score == 16
        @test res_bnw.matches == 8
    end

    @testset "2-Bit Packed Sequence Infrastructure Across Operations" begin
        pseq1 = BioToolkit.PackedDNASeq("ACGTACGTACGTACGTACGTACGTACGTACGT")
        pseq2 = BioToolkit.PackedDNASeq("ACGTACGTACGTACGTACGTACGTACGTACGA")
        @test pseq1 isa BioToolkit.PackedDNASeq
        @test length(pseq1) == 32
        @test String(pseq1) == "ACGTACGTACGTACGTACGTACGTACGTACGT"
        @test pseq1[1] == 'A'
        @test pseq1[2] == 'C'
        @test pseq1[3] == 'G'
        @test pseq1[4] == 'T'
        @test pseq1[32] == 'T'
        @test pseq1 == DNASeq("ACGTACGTACGTACGTACGTACGTACGTACGT")
        
        # 1. Fast Bitwise GC Content
        @test gc_content(pseq1) == 0.5
        
        # 2. Alignment Algorithms
        res_nw = needleman_wunsch(pseq1, pseq2)
        @test res_nw.score > 0
        @test res_nw.matches == 31

        res_sw = smith_waterman(pseq1, pseq2)
        @test res_sw.score > 0

        res_sg = semi_global_align(pseq1, pseq2)
        @test res_sg.score > 0

        res_ov = overlap_align(pseq1, pseq2)
        @test res_ov.score > 0

        res_bnw = banded_needleman_wunsch(pseq1, pseq2; k=5)
        @test res_bnw.score > 0

        # 3. Distance & K-mer Analysis
        @test hamming_distance(pseq1, pseq2) == 1
        kmers = kmer_frequency(pseq1, 4)
        @test kmers["ACGT"] == 8

        # 4. In-place Mutating Functions
        pseq_mut = BioToolkit.PackedDNASeq("ACGT")
        complement!(pseq_mut)
        @test String(pseq_mut) == "TGCA"
        reverse_complement!(pseq_mut)
        @test String(pseq_mut) == "TGCA"

        # 5. Motif Analysis
        mcounts = motif_counts([pseq1, pseq2])
        @test size(mcounts.counts, 2) == 32
    end

    @testset "Biological Symbol System & 3-Letter Amino Acid Parser" begin
        # 1. Nucleotide symbol constants and 1-hot compatibility
        @test DNA_A isa BioToolkit.DNA
        @test Char(DNA_A) == 'A'
        @test Char(DNA_T) == 'T'
        @test Char(RNA_U) == 'U'
        @test iscompatible(DNA_A, DNA_R) # R is A or G
        @test !iscompatible(DNA_C, DNA_R)
        @test iscompatible(DNA_C, DNA_N) # N is any base
        
        # 2. Symbol Predicates & Complement
        @test isGC(DNA_C)
        @test isGC(DNA_G)
        @test !isGC(DNA_A)
        @test ispurine(DNA_A)
        @test ispurine(DNA_G)
        @test ispyrimidine(DNA_C)
        @test ispyrimidine(DNA_T)
        @test isambiguous(DNA_N)
        @test iscertain(DNA_A)
        @test isgap(DNA_Gap)
        @test BioToolkit.complement(DNA_A) == DNA_T
        @test BioToolkit.complement(RNA_U) == RNA_A

        # 3. 3-Letter Amino Acid Parser
        @test parse_amino_acid("ALA") == 'A'
        @test parse_amino_acid("CYS") == 'C'
        @test parse_amino_acid("MET") == 'M'
        @test parse_amino_acid("SEC") == 'U'
        @test parse_amino_acid("PYL") == 'O'
        @test parse_amino_acid("ASX") == 'B'

        # 4. 3-Letter Sequence Construction
        aa_seq = AASeq("ALA-CYS-MET-VAL")
        @test String(aa_seq) == "ACMV"

        # 5. Pairwise Alignment Dispatcher Modes
        s1 = DNASeq("ATGCGATCG")
        s2 = DNASeq("ATGCGATCG")
        @test pairalign(:global, s1, s2).matches == 9
        @test pairalign(:local, s1, s2).matches == 9
        @test pairalign(:semi_global, s1, s2).matches == 9
        @test pairalign(:overlap, s1, s2).matches == 9
        @test pairalign(:banded, s1, s2).matches == 9
        @test pairalign(:pair_hmm, s1, s2).consensus_alignment.matches == 9
        @test pairalign(:differentiable, s1, s2) isa BioToolkit.PairwiseAlignmentResult
        @test pairalign(:codon, s1, s2) isa BioToolkit.PairwiseAlignmentResult
    end

end





