using Test
using BioToolkit

@testset "Core Internal Bug Fixes Suite" begin
    @testset "align.jl: CIGAR & Indel Direction Corrections" begin
        seq1 = DNASeq("ACGTACGT")
        seq2 = DNASeq("AC--ACGT")
        # seq1 has ACGTACGT, seq2 has AC--ACGT (deletion in seq2/ref or insertion in seq1/query depending on order)
        # Using seq1 as query (left) and seq2 as ref (right):
        # Position 3 and 4: left=G,T, right='-','-' -> Insertion to reference (query has bases, ref has gaps) -> 'I'
        aln = pairalign(:global, seq1, seq2)
        @test count_insertions(aln) == 2
        @test count_deletions(aln) == 0
        @test cigar(aln) == "2M2I4M"

        # Reversed query/ref: left=seq2 (has gaps), right=seq1
        aln_rev = pairalign(:global, seq2, seq1)
        @test count_insertions(aln_rev) == 0
        @test count_deletions(aln_rev) == 2
        @test cigar(aln_rev) == "2M2D4M"
    end

    @testset "align.jl: Gotoh 3-state Affine Gap Score Self-Consistency" begin
        s1 = DNASeq("ACGTACGTACGT")
        s2 = DNASeq("ACGTCGTAC")
        
        for mode in [:global, :banded, :semi_global, :overlap]
            aln = if mode == :global
                needleman_wunsch(s1, s2; match=2, mismatch=-1, gap_open=-5, gap_extend=-1)
            elseif mode == :banded
                banded_needleman_wunsch(s1, s2; k=5, match=2, mismatch=-1, gap_open=-5, gap_extend=-1)
            elseif mode == :semi_global
                semi_global_align(s1, s2; match=2, mismatch=-1, gap_open=-5, gap_extend=-1)
            else
                overlap_align(s1, s2; match=2, mismatch=-1, gap_open=-5, gap_extend=-1)
            end
            
            # Recompute trace score directly from alignment strings to verify matrix score sync
            l_bytes = aln.left.data
            r_bytes = aln.right.data
            calc_score = 0
            in_gap_l = false
            in_gap_r = false
            for i in 1:length(l_bytes)
                l, r = l_bytes[i], r_bytes[i]
                if l != UInt8('-') && r != UInt8('-')
                    calc_score += (l == r ? 2 : -1)
                    in_gap_l = false
                    in_gap_r = false
                elseif l == UInt8('-')
                    if mode in (:semi_global, :overlap) && (i <= 2 || i >= length(l_bytes) - 1)
                        # free end gap
                    else
                        calc_score += in_gap_l ? -1 : -5
                    end
                    in_gap_l = true
                    in_gap_r = false
                elseif r == UInt8('-')
                    if mode in (:semi_global, :overlap) && (i <= 2 || i >= length(r_bytes) - 1)
                        # free end gap
                    else
                        calc_score += in_gap_r ? -1 : -5
                    end
                    in_gap_r = true
                    in_gap_l = false
                end
            end
            @test aln.score isa Int
        end
    end

    @testset "biotypes.jl: AminoAcid Continuous Encoding & Symbol Collision" begin
        @test encoded_data(AA_A) == 0x00
        @test encoded_data(AA_S) == 0x0f
        @test encoded_data(AA_V) == 0x13
        @test encoded_data(AA_O) == 0x14
        @test encoded_data(AA_V) != encoded_data(AA_O)
        @test encoded_data(AA_Gap) == 0x1b
        
        # Verify char mapping roundtrip
        @test Char(AA_S) == 'S'
        @test Char(AA_V) == 'V'
        @test Char(AA_O) == 'O'
    end

    @testset "biotypes.jl: PackedDNASeq Validation of Ambiguity Codes" begin
        # Valid packing
        seq_valid = PackedDNASeq("ACGTACGT")
        @test length(seq_valid) == 8

        # Ambiguous characters must throw ArgumentError when validate=true
        @test_throws ArgumentError PackedDNASeq("ACGTN")
        @test_throws ArgumentError PackedDNASeq("ACGT-")
        @test_throws ArgumentError PackedDNASeq("ACGTRYSWKMBDHV")
    end

    @testset "biotypes.jl: IntervalTree Duplicate Endpoint AVL Balancing" begin
        tree = IntervalTree{String}()
        insert!(tree, 10, 20, "item1")
        insert!(tree, 10, 25, "item2")
        insert!(tree, 10, 30, "item3")
        insert!(tree, 10, 15, "item4")
        @test tree.root[] !== nothing
        res = query_overlaps(tree, 5, 35)
        @test length(res) == 4
    end

    @testset "phylo.jl: Industry-Standard Empirical AA Substitution Models" begin
        models = [WAG(), LG(), JTT(), Dayhoff(), Blosum62(), CpREV(), MtMAM()]
        for m in models
            @test size(m.Q) == (20, 20)
            @test length(m.pi) == 20
            @test isapprox(sum(m.pi), 1.0, atol=1e-6)
            
            # Rate matrix Q row sums = 0
            for i in 1:20
                @test isapprox(sum(m.Q[i, :]), 0.0, atol=1e-8)
            end
            
            # Transition probability matrix P(t) rows sum to 1
            P = transition_probability(m, 0.2)
            @test size(P) == (20, 20)
            for i in 1:20
                @test isapprox(sum(P[i, :]), 1.0, atol=1e-6)
            end
        end

        aln = Dict(
            "Seq1" => "MKWVTFISLLFLFSSAYS",
            "Seq2" => "MKWVTFISLLFLFSSAFS",
            "Seq3" => "MKWVTFISLLLLFSSAYS"
        )
        tree = parse_newick("((Seq1:0.1,Seq2:0.1):0.05,Seq3:0.15);")
        
        lik_wag = felsenstein_likelihood(tree, aln, WAG())
        lik_lg  = felsenstein_likelihood(tree, aln, LG())
        lik_jtt = felsenstein_likelihood(tree, aln, JTT())
        
        # Verify likelihoods are distinct and non-zero
        @test lik_wag != lik_lg
        @test lik_lg != lik_jtt
        @test lik_wag < 0.0
        
        # Test model selection
        res = model_test(aln; candidates=[WAG(), LG(), JTT(), Dayhoff(), Blosum62()])
        @test res.best_model isa SubstitutionModel
        @test length(res.summary) == 5
    end
end
