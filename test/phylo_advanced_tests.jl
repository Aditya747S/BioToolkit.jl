using Test
using BioToolkit

@testset "Phylogenetics & MSA Advanced Capabilities" begin
    @testset "Custom CTMC & Spectral Decomposition" begin
        # Test 4x4 DNA custom rate matrix
        Q = [-1.0 0.33 0.33 0.34;
              0.33 -1.0 0.34 0.33;
              0.33 0.34 -1.0 0.33;
              0.34 0.33 0.33 -1.0]
        pi = [0.25, 0.25, 0.25, 0.25]
        model = CustomCTMC(Q, pi)
        
        P = transition_probability(model, 0.5)
        @test size(P) == (4, 4)
        @test all(P .>= 0.0)
        @test all(P .<= 1.0)
        @test isapprox(sum(P, dims=2), ones(4, 1), atol=1e-4)
    end

    @testset "Codon Substitution Models (GY94 & MG94)" begin
        gy = GY94(kappa=2.0, omega=1.5)
        @test gy isa CustomCTMC
        @test size(gy.Q) == (61, 61)

        mg = MG94(kappa=2.0, omega=1.0)
        @test mg isa CustomCTMC
        @test size(mg.Q) == (61, 61)

        P_gy = transition_probability(gy, 0.1)
        @test size(P_gy) == (61, 61)
        @test isapprox(sum(P_gy, dims=2), ones(61, 1), atol=1e-3)
    end

    @testset "Fixed Effects Likelihood Site Selection Test (FEL)" begin
        seq1 = "ATGGCCATTGTAATGGCCATTGTA"
        seq2 = "ATGGCCATTGTAATGGCCATTGTA"
        seq3 = "ATGGCCATTGTAATGGCCATTGTA"
        aln = Dict("s1" => seq1, "s2" => seq2, "s3" => seq3)

        t1 = PhyloTree("s1"; branch_length=0.1)
        t2 = PhyloTree("s2"; branch_length=0.1)
        t3 = PhyloTree("s3"; branch_length=0.1)
        tree = PhyloTree([t1, t2, t3])

        fel_res = fel_selection_test(aln, tree)
        @test haskey(fel_res, :site)
        @test haskey(fel_res, :omega)
        @test haskey(fel_res, :p_value)
        @test length(fel_res.site) == 8 # 24 / 3 = 8 codons
    end

    @testset "Subtree Pruning & Regrafting ML Search (SPR)" begin
        seqs = Dict(
            "TaxonA" => "ATGCATGCATGC",
            "TaxonB" => "ATGCATGCATGT",
            "TaxonC" => "ATGCATGGATGC",
            "TaxonD" => "ATGCATGGATGT"
        )
        best_tree = spr_ml_search(seqs; model=JC69(), max_iters=5)
        @test best_tree isa PhyloTree
        @test count_terminals(best_tree) == 4
    end

    @testset "Native Progressive MSA & Gap Trimming" begin
        seqs = ["ATGCATGC", "ATGCACGC", "ATGCAAGC"]
        ids = ["s1", "s2", "s3"]
        msa = progressive_msa(seqs; identifiers=ids)
        @test msa isa BioToolkit.MultipleSequenceAlignment
        @test length(msa) == 3
        @test BioToolkit.get_alignment_length(msa) >= 8

        # Gap trimming test
        rec1 = BioToolkit.SeqRecordLite(BioToolkit._msa_sequence_from_string("ATGC--AT"); identifier="s1", name="s1")
        rec2 = BioToolkit.SeqRecordLite(BioToolkit._msa_sequence_from_string("ATGC--AC"); identifier="s2", name="s2")
        rec3 = BioToolkit.SeqRecordLite(BioToolkit._msa_sequence_from_string("ATGCGGAT"); identifier="s3", name="s3")
        gapped_msa = BioToolkit.MultipleSequenceAlignment([rec1, rec2, rec3])
        
        @test BioToolkit.get_alignment_length(gapped_msa) == 8
        trim_gaps!(gapped_msa; max_gap_fraction=0.5)
        @test BioToolkit.get_alignment_length(gapped_msa) == 6
    end

    @testset "Free-Rate Models (+R) & SIMMAP" begin
        frm = DiscreteFreeRateModel([0.2, 1.0, 2.5], [0.2, 0.5, 0.3])
        @test frm isa DiscreteFreeRateModel
        @test length(frm.rates) == 3
        @test isapprox(sum(frm.weights .* frm.rates), 1.0, atol=1e-4)

        t1 = PhyloTree("A"; branch_length=0.2)
        t2 = PhyloTree("B"; branch_length=0.2)
        tree = PhyloTree([t1, t2])
        Q = [-1.0 1.0; 1.0 -1.0]
        tips = Dict("A" => 1, "B" => 2)
        
        sims = simmap(tree, tips, Q; n_sims=3)
        @test length(sims) == 3
    end

    @testset "Neighbor-Net Split Networks & Joint Ancestral Reconstruction" begin
        dm = [0.0 0.1 0.4;
              0.1 0.0 0.5;
              0.4 0.5 0.0]
        taxa = ["A", "B", "C"]
        net = neighbor_net(dm, taxa)
        @test length(net) > 0

        t1 = PhyloTree("A"; branch_length=0.1)
        t2 = PhyloTree("B"; branch_length=0.1)
        tree = PhyloTree([t1, t2])
        tip_states = Dict("A" => "A", "B" => "G")
        joint_states = joint_ancestral_reconstruction(tree, tip_states, JC69())
        @test joint_states isa Dict
    end

    @testset "Advanced MSA: Profile Alignment, End Trimming, Clustering & PSSM" begin
        # 1. Profile-Profile Alignment
        msa1 = progressive_msa(["ATGCATGC", "ATGCACGC"]; identifiers=["s1", "s2"])
        msa2 = progressive_msa(["ATGCAAGC", "ATGCAGGC"]; identifiers=["s3", "s4"])
        prof_aln = profile_profile_align(msa1, msa2)
        @test prof_aln isa BioToolkit.MultipleSequenceAlignment
        @test length(prof_aln) == 4

        # 2. End Trimming
        rec1 = BioToolkit.SeqRecordLite(BioToolkit._msa_sequence_from_string("---ATGCATGC---"); identifier="s1", name="s1")
        rec2 = BioToolkit.SeqRecordLite(BioToolkit._msa_sequence_from_string("---ATGCACGC---"); identifier="s2", name="s2")
        gapped_ends = BioToolkit.MultipleSequenceAlignment([rec1, rec2])
        trim_ends!(gapped_ends; min_coverage=0.5)
        @test BioToolkit.get_alignment_length(gapped_ends) == 8

        # 3. Clustering
        clustered = cluster_alignment(prof_aln; min_identity=0.8)
        @test clustered isa BioToolkit.MultipleSequenceAlignment
        @test length(clustered) >= 1

        # 4. PSSM & Information Content
        pssm_data = msa_pssm(msa1)
        @test haskey(pssm_data, :pssm)
        @test haskey(pssm_data, :information_content)
        @test length(pssm_data.information_content) == BioToolkit.get_alignment_length(msa1)
    end

    @testset "R Package Parity (phangorn / ape / DECIPHER)" begin
        # 1. Consistency & Retention Indices
        t1 = PhyloTree("s1"; branch_length=0.1)
        t2 = PhyloTree("s2"; branch_length=0.1)
        tree = PhyloTree([t1, t2])
        aln = Dict("s1" => "ATGC", "s2" => "ATGC")
        ci = consistency_index(tree, aln)
        ri = retention_index(tree, aln)
        @test ci == 1.0
        @test ri == 1.0

        # 2. Weighted Robinson-Foulds Distance
        wrf = weighted_robinson_foulds(tree, tree)
        @test wrf == 0.0

        # 3. Maximum Clade Credibility
        mcc = max_clade_credibility([tree, tree])
        @test mcc isa PhyloTree

        # 4. Alignment Distance Matrix
        msa1 = progressive_msa(["ATGCATGC", "ATGCACGC"]; identifiers=["s1", "s2"])
        dist_res = alignment_distance_matrix(msa1; model=:kimura)
        @test haskey(dist_res, :matrix)
        @test size(dist_res.matrix) == (2, 2)

        # 5. Mask Alignment
        masked = mask_alignment(msa1; max_entropy=1.5)
        @test masked isa BioToolkit.MultipleSequenceAlignment
    end

    @testset "Likelihood +G/+I Heterogeneity & Unrooted Tree Architecture" begin
        # 1. UnrootedPhyloTree Architecture
        u1 = UnrootedPhyloTree("A")
        u2 = UnrootedPhyloTree("B")
        BioToolkit.add_edge!(u1, u2, 0.2)
        @test u1 isa AbstractPhyloTree
        @test u2 isa AbstractPhyloTree
        @test length(u1.neighbors) == 1
        @test u1.edge_lengths[1] == 0.2

        # 2. felsenstein_likelihood with +G (Gamma Rates) and +I (Invariable Sites)
        t1 = PhyloTree("s1"; branch_length=0.1)
        t2 = PhyloTree("s2"; branch_length=0.1)
        tree = PhyloTree([t1, t2])
        aln = Dict("s1" => "ATGCATGC", "s2" => "ATGCACGC")
        
        gamma_model = DiscreteGammaRate(0.5; n_categories=4)
        inv_model = InvariableSites(0.2)
        
        lik_base = felsenstein_likelihood(tree, aln, JC69())
        lik_gamma = felsenstein_likelihood(tree, aln, JC69(); rates=gamma_model)
        lik_inv = felsenstein_likelihood(tree, aln, JC69(); p_inv=inv_model)
        lik_combo = felsenstein_likelihood(tree, aln, JC69(); rates=gamma_model, p_inv=inv_model)
        
        @test lik_base isa Float64
        @test lik_gamma isa Float64
        @test lik_inv isa Float64
        @test lik_combo isa Float64
        @test lik_base != lik_gamma
        @test lik_base != lik_inv
    end
end
