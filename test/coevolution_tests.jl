using Test
using Random
using BioToolkit

@testset "Co-evolutionary contact inference and folding" begin
    @testset "Alignment filtering and sequence reweighting" begin
        msa_with_gaps = MultipleSequenceAlignment([
            "ACD-EFGH",
            "ACDGEFGH",
            "ACD-EYGH",
            "ACDGEYGH",
        ])

        filtered = filter_alignment_for_dca(msa_with_gaps; max_gap_fraction=0.6, min_sequence_coverage=0.7)
        @test length(filtered) == 4
        @test size(filtered)[2] >= 7

        encoded = [1 2 3 4; 1 2 3 4; 2 2 3 4; 2 3 3 4]
        weights = sequence_reweighting(encoded; identity_threshold=1.0)
        @test length(weights) == size(encoded, 1)
        @test all(weight > 0 for weight in weights)
        @test weights[1] ≈ weights[2]
    end

    @testset "Pseudo-likelihood model fitting and contact scoring" begin
        rng = MersenneTwister(2026)
        amino_acids = collect("ACDEFGHIKLMNPQRSTVWY")
        sequences = String[]

        for _ in 1:80
            state = rand(rng, Bool)
            chars = [rand(rng, amino_acids) for _ in 1:10]
            chars[1] = state ? 'A' : 'G'
            chars[6] = state ? 'V' : 'L'
            push!(sequences, String(chars))
        end

        msa = MultipleSequenceAlignment(sequences)
        model = fit_pseudolikelihood_model(msa; max_gap_fraction=0.8, identity_threshold=0.9, regularization=0.05)

        @test size(model.fields, 1) == size(model.raw_scores, 1)
        @test size(model.couplings, 1) == size(model.couplings, 2)
        @test size(model.raw_scores) == size(model.apc_scores)
        @test model.effective_sequences > 1.0

        scores = compute_contact_scores(model; apc=true, min_separation=1)
        @test scores ≈ scores'
        @test all(scores[i, i] == 0.0 for i in axes(scores, 1))

        contact_map, _ = predict_contact_map(msa; return_model=true, min_separation=1, top_l=12, regularization=0.05)
        @test size(contact_map.scores, 1) == 10
        @test maximum(contact_map.scores) <= 1.0 + 1e-8
        @test minimum(contact_map.scores) >= -1e-8

        top_pairs = top_contact_pairs(contact_map; top_n=5, min_separation=1)
        @test !isempty(top_pairs)
        @test any((pair[1] == 1 && pair[2] == 6) || (pair[1] == 6 && pair[2] == 1) for pair in top_pairs)
    end

    @testset "Contact-guided folding" begin
        scores = zeros(Float64, 8, 8)
        scores[1, 6] = scores[6, 1] = 1.0
        scores[2, 7] = scores[7, 2] = 0.8
        scores[3, 8] = scores[8, 3] = 0.7
        cmap = ContactMap(scores, collect(1:8))

        structure = fold_from_contacts(cmap; top_n=6, iterations=500, learning_rate=0.02, seed=11)
        @test length(structure.models) == 1
        @test length(structure.models[1].chains) == 1
        @test length(structure.models[1].chains[1].residues) == 8

        atoms = [atom for residue in structure.models[1].chains[1].residues for atom in residue.atoms]
        @test !isempty(atoms)
        @test all(isfinite(atom.x) && isfinite(atom.y) && isfinite(atom.z) for atom in atoms)

        first_ca = structure.models[1].chains[1].residues[1].atoms[1]
        sixth_ca = structure.models[1].chains[1].residues[6].atoms[1]
        distance = sqrt((first_ca.x - sixth_ca.x)^2 + (first_ca.y - sixth_ca.y)^2 + (first_ca.z - sixth_ca.z)^2)
        @test distance < 12.0
    end
end
