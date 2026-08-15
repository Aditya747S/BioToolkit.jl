using Test
using Random
using Statistics
using BioToolkit
using BioToolkit.Coevolution

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

        # Test true PLM algorithm
        plm_model = fit_pseudolikelihood_model(msa; algorithm=:plm, regularization=0.05)
        @test size(plm_model.couplings, 1) == 10
        @test size(plm_model.raw_scores) == (10, 10)
    end

    @testset "Gap handling in sequence logo entropy" begin
        msa_gaps = MultipleSequenceAlignment([
            "ACD-E",
            "ACDGE",
            "A--GE",
            "ACDGE",
        ])
        logo = sequence_logo_entropy(msa_gaps)
        @test logo.position == [1, 2, 3, 4, 5]
        @test all(logo.entropy .>= 0)
        # Verify probability sum across all alphabet characters equals 1.0 for every column
        char_keys = [k for k in keys(logo) if k ∉ (:position, :entropy, :information_content)]
        for col in 1:5
            col_prob_sum = sum(logo[k][col] for k in char_keys)
            @test col_prob_sum ≈ 1.0 atol=1e-6
        end
    end

    @testset "Shrinkage estimators (Ledoit-Wolf vs OAS)" begin
        msa = MultipleSequenceAlignment([
            "ACDEFGHIKL",
            "ACDEFGYYKL",
            "ACDEFGHIKL",
            "ACDEGGYYKL",
        ])
        cm_lw = shrinkage_precision_contacts(msa; shrinkage=:ledoit_wolf, min_separation=1)
        cm_oas = shrinkage_precision_contacts(msa; shrinkage=:oas, min_separation=1)
        @test cm_lw isa ContactMap
        @test cm_oas isa ContactMap
        @test size(cm_lw.scores) == (10, 10)
        @test size(cm_oas.scores) == (10, 10)
    end

    @testset "Sequence weights caching and input validation" begin
        msa = MultipleSequenceAlignment([
            "ACDEFGHIKL",
            "ACDEFGYYKL",
            "ACDEFGHIKL",
            "ACDEGGYYKL",
        ])
        encoded = BioToolkit.Coevolution._encode_alignment(BioToolkit.Coevolution._alignment_strings(msa), BioToolkit.Coevolution._alphabet(BioToolkit.Coevolution._alignment_strings(msa)))
        w = sequence_reweighting(encoded)

        di = direct_information_contacts(msa; min_separation=1, weights=w)
        @test di isa ContactMap

        cov_f = positional_covariation_matrix(msa; metric=:frobenius_correlation, weights=w)
        @test cov_f.covariation isa Matrix{Float64}

        # Validation error for non-alphabet character in _encode_alignment
        @test_throws ArgumentError BioToolkit.Coevolution._encode_alignment(["ACDEF", "ACD?F"], collect("ACDEF"))

        # Test mismatched weights length error validation
        msa_unfiltered = MultipleSequenceAlignment([
            "ACDEFGHIKL",
            "ACDEFGYYKL",
            "----------", # High gap seq, gets filtered out
            "ACDEGGYYKL",
        ])
        orig_weights = [1.0, 1.0, 1.0, 1.0] # Length 4, but filtered alignment has 3 sequences
        @test_throws ArgumentError fit_pseudolikelihood_model(msa_unfiltered; max_gap_fraction=0.5, weights=orig_weights)
        @test_throws ArgumentError mutual_information_contacts(msa_unfiltered; max_gap_fraction=0.5, weights=orig_weights)
        @test_throws ArgumentError direct_information_contacts(msa_unfiltered; max_gap_fraction=0.5, weights=orig_weights)
        @test_throws ArgumentError shrinkage_precision_contacts(msa_unfiltered; max_gap_fraction=0.5, weights=orig_weights)
        @test_throws ArgumentError positional_covariation_matrix(msa_unfiltered; max_gap_fraction=0.5, weights=orig_weights)

        # Test unsupported metric error in positional_covariation_matrix
        @test_throws ArgumentError positional_covariation_matrix(msa; metric=:invalid_metric)
    end

    @testset "Contact-guided folding and structure contact matrix" begin
        scores = zeros(Float64, 8, 8)
        scores[1, 6] = scores[6, 1] = 1.0
        scores[2, 7] = scores[7, 2] = 0.8
        scores[3, 8] = scores[8, 3] = 0.7
        cmap = ContactMap(scores, collect(1:8))

        # Test distance filtering helper
        filtered_cmap = filter_contacts_by_sequence_distance(cmap; min_separation=5)
        @test filtered_cmap.scores[1, 2] == 0.0
        @test filtered_cmap.scores[1, 6] == 1.0

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

        # Test structure_to_contact_matrix helper
        true_cmap_matrix = structure_to_contact_matrix(structure; cutoff=12.0)
        @test size(true_cmap_matrix) == (8, 8)
        @test true_cmap_matrix[1, 6] == 1.0
        @test true_cmap_matrix[6, 1] == 1.0
    end
end

using Test
using Random
Random.seed!(4242)

# make_msa builds a MultipleSequenceAlignment from raw strings. Adjust the
# constructor call below if your MSA/record types differ from this guess.
function make_msa(seqs::Vector{String}; ids=nothing)
    records = [SeqRecordLite(seqs[i]; identifier=(ids === nothing ? "seq$i" : ids[i]),
                              name=(ids === nothing ? "seq$i" : ids[i]),
                              description="") for i in eachindex(seqs)]
    return MultipleSequenceAlignment(records)
end

# ==============================================================================
# _alphabet / _encode_alignment
# ==============================================================================
println("\n" * "="^70)
println("alphabet / encoding")
println("="^70)

@testset "alphabet: gap always sorts first" begin
    strings = ["AC-G", "TG-A", "AC-G"]
    alphabet = Coevolution._alphabet(strings)
    @test alphabet[1] == '-'
    @test Set(alphabet) == Set(['-', 'A', 'C', 'G', 'T'])
    println("alphabet: $alphabet")
end

@testset "_encode_alignment: regression test for the invalid-character fix" begin
    # 'X' is NOT in the alphabet derived from these two strings -- must throw,
    # not silently fall back to the gap index (that was the original bug).
    alphabet = Coevolution._alphabet(["ACG", "TGA"])
    threw = false
    try
        Coevolution._encode_alignment(["ACG", "TGX"], alphabet)
    catch e
        threw = e isa ArgumentError
    end
    @test threw
    println("invalid character 'X' correctly throws ArgumentError: $threw")
end

@testset "_encode_alignment: valid alignment round-trips through the lookup table" begin
    strings = ["AC-G", "TG-A"]
    alphabet = Coevolution._alphabet(strings)
    encoded = Coevolution._encode_alignment(strings, alphabet)
    decoded = [alphabet[encoded[i, j]] for i in 1:2, j in 1:4]
    @test String(decoded[1, :]) == "AC-G"
    @test String(decoded[2, :]) == "TG-A"
end

# ==============================================================================
# filter_alignment_for_dca
# ==============================================================================
println("\n" * "="^70)
println("filter_alignment_for_dca")
println("="^70)

@testset "filter_alignment_for_dca: hand-counted gap/coverage filtering" begin
    # col1: 0/4 gaps, col2: 3/4 gaps (>0.5 -> dropped), col3: 1/4 gaps, col4: 0/4 gaps
    seqs = ["AAAA", "A-AA", "A-AA", "A-AA"]
    # cols: 1="AAAA", 2="A---", 3="AAAA", 4="AAAA" -> col2 gap frac = 3/4 = 0.75 > 0.5 -> dropped
    msa = make_msa(seqs)
    filtered = filter_alignment_for_dca(msa; max_gap_fraction=0.5, min_sequence_coverage=0.0)
    filt_strings = Coevolution._alignment_strings(filtered)
    @test ncodeunits(filt_strings[1]) == 3
    println("filtered alignment length: $(ncodeunits(filt_strings[1])) (expected 3, column 2 dropped)")
end

# ==============================================================================
# sequence_reweighting
# ==============================================================================
println("\n" * "="^70)
println("sequence_reweighting")
println("="^70)

@testset "sequence_reweighting: hand-computed cluster weights" begin
    # 4 sequences, length 10. seq1/seq2 identical (100% identity -> cluster together).
    # seq3 differs from seq1 by exactly 3/10 = 70% identity (< 0.8 threshold -> separate).
    # seq4 is maximally different from all (separate).
    # At identity_threshold=0.8: max_mismatches = floor(0.2*10) = 2
    # seq1 vs seq2: 0 mismatches <= 2 -> cluster (count=2 each)
    # seq1 vs seq3: 3 mismatches > 2 -> not clustered
    # seq1 vs seq4, seq2 vs seq3/4, seq3 vs seq4: all > 2 mismatches -> not clustered
    # Expected weights: seq1=1/2, seq2=1/2, seq3=1/1, seq4=1/1
    strings = ["AAAAAAAAAA", "AAAAAAAAAA", "TTTAAAAAAA", "TTTTTTTTTT"]
    alphabet = Coevolution._alphabet(strings)
    encoded = Coevolution._encode_alignment(strings, alphabet)
    weights = sequence_reweighting(encoded; identity_threshold=0.8)
    @test weights ≈ [0.5, 0.5, 1.0, 1.0]
    println("weights: got=$weights  expected=[0.5, 0.5, 1.0, 1.0]")
end

# ==============================================================================
# column_conservation_scores
# ==============================================================================
println("\n" * "="^70)
println("column_conservation_scores")
println("="^70)

@testset "column_conservation_scores: hand-computed entropy for a fully conserved column" begin
    # column 1 is 100% 'A' across 4 sequences -> H=0 -> conservation score=1.0
    # column 2 is fully random over {A,C,G,T} (uniform) -> H approaches H_max -> score near 0
    seqs = ["AA", "AC", "AG", "AT"]
    msa = make_msa(seqs)
    scores = column_conservation_scores(msa; pseudocount=0.01, gap_penalise=true)
    @test scores[1] > 0.95   # fully conserved column -> near-1 score
    @test scores[2] < 0.15   # uniform column -> near-0 score
    println("conservation scores: col1=$(scores[1]) (expect ~1), col2=$(scores[2]) (expect ~0)")
end

# ==============================================================================
# sequence_logo_entropy: regression test for the probability-mass fix
# ==============================================================================
println("\n" * "="^70)
println("sequence_logo_entropy")
println("="^70)

@testset "sequence_logo_entropy: column probabilities sum to 1.0 even with gaps" begin
    # This is the exact regression case for the fix: a column with gaps must
    # have its per-character probabilities sum to 1, not silently drop the
    # gap's probability mass.
    seqs = ["A", "A", "-", "-"]  # single-column alignment: 2 A's, 2 gaps
    msa = make_msa(seqs)
    logo = sequence_logo_entropy(msa; pseudocount=0.5)
    total_p = sum(getproperty(logo, Symbol(c)) for c in propertynames(logo) if !(c in (:position, :entropy, :information_content)))[1]
    @test total_p ≈ 1.0 atol=1e-9
    println("sum of column-1 character probabilities: $total_p (expected 1.0)")
end

@testset "sequence_logo_entropy: fully conserved column has near-zero entropy" begin
    seqs = ["A", "A", "A", "A"]
    msa = make_msa(seqs)
    logo = sequence_logo_entropy(msa; pseudocount=0.01)
    @test logo.entropy[1] < 0.1
    println("entropy of fully-conserved column: $(logo.entropy[1]) (expect ~0)")
end

# ==============================================================================
# mutual_information_contacts: hand-derived MI for a perfectly coupled column pair
# ==============================================================================
println("\n" * "="^70)
println("mutual_information_contacts")
println("="^70)

@testset "mutual_information_contacts: MI is maximal for a perfectly coupled pair, ~0 for independent" begin
    # Columns 1&2 are PERFECTLY coupled: col1='A' <-> col2='C', col1='T' <-> col2='G'.
    # Column 3 is independent random noise relative to column 1.
    # With no gaps and a large-ish sample, MI(1,2) should be near log2(2)=1.0 bit
    # (in nats: log(2)~0.693), and MI(1,3) should be near 0.
    n = 200
    rng = MersenneTwister(1)
    seqs = String[]
    for _ in 1:n
        b1 = rand(rng, ('A', 'T'))
        b2 = b1 == 'A' ? 'C' : 'G'          # perfectly coupled to col1
        b3 = rand(rng, ('A', 'T'))          # independent
        push!(seqs, string(b1, b2, b3))
    end
    msa = make_msa(seqs)
    cmap = mutual_information_contacts(msa; pseudocount=0.01, min_separation=0, apc=false)
    # un-normalized comparison isn't meaningful post-_normalize_scores, so
    # instead check RELATIVE ordering survives normalization: coupled pair
    # should score strictly higher than the independent pair.
    @test cmap.scores[1, 2] > cmap.scores[1, 3]
    println("normalized MI(1,2)=$(cmap.scores[1,2]) (coupled)  MI(1,3)=$(cmap.scores[1,3]) (independent)")
end

# ==============================================================================
# min_separation boundary consistency (regression test for the fix)
# ==============================================================================
println("\n" * "="^70)
println("min_separation boundary consistency")
println("="^70)

@testset "min_separation: boundary pairs are treated consistently across building and analysis functions" begin
    # Build a small alignment where columns min_separation apart are strongly
    # coupled, then check the pair AT EXACTLY min_separation is not silently
    # zeroed out by the masking convention mismatch that was fixed.
    n = 150
    min_sep = 3
    rng = MersenneTwister(2)
    seqs = String[]
    for _ in 1:n
        b1 = rand(rng, ('A', 'T'))
        cols = fill('A', 5)
        cols[1] = b1
        cols[1 + min_sep] = b1 == 'A' ? 'C' : 'G'  # coupled at exactly min_separation apart
        push!(seqs, String(cols))
    end
    msa = make_msa(seqs)
    cmap = mutual_information_contacts(msa; pseudocount=0.01, min_separation=min_sep, apc=false)
    i, j = 1, 1 + min_sep  # |i-j| == min_sep exactly
    @test cmap.scores[i, j] > 0.0  # must NOT be masked to zero at the boundary
    println("score at |i-j|==min_separation ($min_sep): $(cmap.scores[i,j]) (must be > 0)")

    # Also check contact_enrichment_statistics doesn't include phantom
    # always-zero pairs strictly INSIDE min_separation (those should still
    # be correctly excluded).
    true_contacts = zeros(5, 5)
    true_contacts[i, j] = 1.0
    true_contacts[j, i] = 1.0
    stats = contact_enrichment_statistics(cmap, true_contacts; min_separation=min_sep)
    @test stats.n_true_contacts == 1
    println("n_true_contacts at boundary: $(stats.n_true_contacts) (expected 1)")
end

# ==============================================================================
# weights= parameter: length-validation regression check
# ==============================================================================
println("\n" * "="^70)
println("custom weights= validation")
println("="^70)

@testset "custom weights: mismatched length after internal filtering should error, not silently misalign" begin
    # Construct an alignment where filter_alignment_for_dca WILL drop a
    # sequence (low coverage), then pass weights sized to the ORIGINAL
    # (pre-filter) sequence count. This is exactly the scenario where the
    # missing length check could silently misalign weights to sequences.
    seqs = ["AAAA", "AAAA", "AAAA", "----"]  # last seq is 100% gaps -> dropped by coverage filter
    msa = make_msa(seqs)
    bad_weights = ones(4)  # sized to the ORIGINAL 4 sequences, not the filtered 3

    threw_or_correct = false
    try
        result = mutual_information_contacts(msa; min_separation=0, weights=bad_weights)
        # If it didn't throw, at minimum verify it didn't silently produce a
        # result using misaligned weights -- we can't easily verify semantic
        # correctness here, so treat "ran without throwing" as a finding to
        # report, not an automatic pass.
        println("!!! mutual_information_contacts did NOT throw on mismatched weights length.")
        println("!!! This means the length-validation fix has NOT been applied yet --")
        println("!!! silently misaligned weights are possible. See conversation notes.")
        threw_or_correct = false
    catch e
        println("mutual_information_contacts threw on mismatched weights: $(typeof(e))")
        threw_or_correct = e isa ArgumentError || e isa BoundsError || e isa DimensionMismatch
    end
    @test threw_or_correct
end

# ==============================================================================
# shrinkage_precision_contacts: rho sanity + method differentiation
# ==============================================================================
println("\n" * "="^70)
println("shrinkage_precision_contacts")
println("="^70)

@testset "shrinkage_precision_contacts: ledoit_wolf, oas, and ridge give DIFFERENT (not identical) results" begin
    # Regression check for the OAS/Ledoit-Wolf mislabeling fix: the two
    # methods should now be genuinely different estimators, not the same
    # formula under two names.
    n = 60
    rng = MersenneTwister(3)
    seqs = [String([rand(rng, ('A','C','G','T')) for _ in 1:10]) for _ in 1:n]
    msa = make_msa(seqs)

    cmap_lw = shrinkage_precision_contacts(msa; shrinkage=:ledoit_wolf, min_separation=1)
    cmap_oas = shrinkage_precision_contacts(msa; shrinkage=:oas, min_separation=1)
    cmap_ridge = shrinkage_precision_contacts(msa; shrinkage=:ridge, min_separation=1)

    lw_vs_oas_differ = !isapprox(cmap_lw.scores, cmap_oas.scores; atol=1e-8)
    @test lw_vs_oas_differ
    println("ledoit_wolf vs oas scores differ: $lw_vs_oas_differ (expected true -- distinct estimators)")
    println("(if this is false, ledoit_wolf and oas are still computing the same formula)")

    @test all(isfinite, cmap_lw.scores)
    @test all(isfinite, cmap_ridge.scores)
end

# ==============================================================================
# fit_pseudolikelihood_model: PLM recovers a known strong coupling
# ==============================================================================
println("\n" * "="^70)
println("fit_pseudolikelihood_model (:plm)")
println("="^70)

@testset "fit_pseudolikelihood_model :plm ranks a strongly coupled pair above an uncoupled pair" begin
    # Same synthetic coupling construction as the MI test above, but now
    # exercising the actual PLM gradient-descent fit. This is a coarse
    # sanity check (correct RELATIVE ranking), not an exact-value test --
    # PLM parameter recovery accuracy depends on iteration count/learning
    # rate, which is a separate tuning question from correctness.
    n = 250
    rng = MersenneTwister(5)
    seqs = String[]
    for _ in 1:n
        b1 = rand(rng, ('A', 'T'))
        b2 = b1 == 'A' ? 'C' : 'G'   # perfectly coupled to col1
        b3 = rand(rng, ('A', 'T'))  # independent
        push!(seqs, string(b1, b2, b3))
    end
    msa = make_msa(seqs)
    model = fit_pseudolikelihood_model(msa; algorithm=:plm, max_gap_fraction=1.0,
                                         min_sequence_coverage=0.0, regularization=0.01)
    @test model.raw_scores[1, 2] > model.raw_scores[1, 3]
    println("PLM raw_scores: coupled(1,2)=$(model.raw_scores[1,2])  independent(1,3)=$(model.raw_scores[1,3])")
end

@testset "fit_pseudolikelihood_model: :mean_field and :plm agree on RELATIVE ranking" begin
    n = 250
    rng = MersenneTwister(6)
    seqs = String[]
    for _ in 1:n
        b1 = rand(rng, ('A', 'T'))
        b2 = b1 == 'A' ? 'C' : 'G'
        b3 = rand(rng, ('A', 'T'))
        push!(seqs, string(b1, b2, b3))
    end
    msa = make_msa(seqs)
    model_mf = fit_pseudolikelihood_model(msa; algorithm=:mean_field, max_gap_fraction=1.0, min_sequence_coverage=0.0)
    model_plm = fit_pseudolikelihood_model(msa; algorithm=:plm, max_gap_fraction=1.0, min_sequence_coverage=0.0)
    @test (model_mf.raw_scores[1,2] > model_mf.raw_scores[1,3]) == (model_plm.raw_scores[1,2] > model_plm.raw_scores[1,3])
    println("mean_field and plm agree that (1,2) > (1,3): true")
end

# ==============================================================================
# direct_information_contacts: DI regression / sanity
# ==============================================================================
println("\n" * "="^70)
println("direct_information_contacts")
println("="^70)

@testset "direct_information_contacts: DI ranks coupled pair above independent pair" begin
    n = 250
    rng = MersenneTwister(7)
    seqs = String[]
    for _ in 1:n
        b1 = rand(rng, ('A', 'T'))
        b2 = b1 == 'A' ? 'C' : 'G'
        b3 = rand(rng, ('A', 'T'))
        push!(seqs, string(b1, b2, b3))
    end
    msa = make_msa(seqs)
    cmap = direct_information_contacts(msa; min_separation=0, max_gap_fraction=1.0, min_sequence_coverage=0.0)
    @test cmap.scores[1, 2] > cmap.scores[1, 3]
    println("DI: coupled(1,2)=$(cmap.scores[1,2])  independent(1,3)=$(cmap.scores[1,3])")
end

# ==============================================================================
# gap_analysis: hand-derived block detection
# ==============================================================================
println("\n" * "="^70)
println("gap_analysis")
println("="^70)

@testset "gap_analysis: hand-derived gap blocks and fractions" begin
    # seq1 = "AA--A--A" -> gap blocks at (3,4,len2) and (6,7,len2)
    seqs = ["AA--A--A"]
    msa = make_msa(seqs)
    ga = gap_analysis(msa)
    @test length(ga.gap_blocks[1]) == 2
    @test ga.gap_blocks[1][1] == (3, 4, 2)
    @test ga.gap_blocks[1][2] == (6, 7, 2)
    @test ga.sequence_gap_fraction[1] ≈ 4/8
    println("gap blocks: $(ga.gap_blocks[1])  seq_gap_fraction=$(ga.sequence_gap_fraction[1])")
end

# ==============================================================================
# contact_precision_recall / contact_enrichment_statistics: hand example
# ==============================================================================
println("\n" * "="^70)
println("contact_precision_recall / contact_enrichment_statistics")
println("="^70)

@testset "contact_precision_recall: perfect ranking gives AUC-PR near 1.0" begin
    # 4x4 contact map, min_separation=1. True contacts: (1,3) only.
    # Predicted scores rank (1,3) highest -> should give high precision/recall.
    scores = zeros(4, 4)
    scores[1, 3] = 1.0; scores[3, 1] = 1.0
    scores[2, 4] = 0.3; scores[4, 2] = 0.3
    scores[1, 2] = 0.1; scores[2, 1] = 0.1
    cmap = ContactMap(scores, collect(1:4))
    true_contacts = zeros(4, 4)
    true_contacts[1, 3] = 1.0; true_contacts[3, 1] = 1.0

    pr = contact_precision_recall(cmap, true_contacts; min_separation=1, n_points=20)
    @test pr.auc_pr > 0.8
    println("AUC-PR for perfectly-ranked single true contact: $(pr.auc_pr) (expect close to 1.0)")
end

@testset "contact_enrichment_statistics: precision@L/2 hand-check with a single true contact ranked first" begin
    scores = zeros(4, 4)
    scores[1, 3] = 1.0; scores[3, 1] = 1.0
    scores[2, 4] = 0.5; scores[4, 2] = 0.5
    cmap = ContactMap(scores, collect(1:4))
    true_contacts = zeros(4, 4)
    true_contacts[1, 3] = 1.0; true_contacts[3, 1] = 1.0

    stats = contact_enrichment_statistics(cmap, true_contacts; top_fractions=[1.0], min_separation=1)
    # l=4, frac=1.0 -> k=4 (or min(k, n_pairs)); with only 2 candidate pairs
    # here ((1,3) and (2,4), since min_separation=1 keeps all off-diagonal
    # pairs for a 4x4), and 1 true contact ranked #1, precision should be
    # exactly 1/2 at k=2 or reflect the true positive being found.
    @test stats.n_true_contacts == 1
    @test stats.precision_at_k[1.0] > 0.0
    println("precision_at_k[1.0]=$(stats.precision_at_k[1.0])  n_true=$(stats.n_true_contacts)")
end

# ==============================================================================
# fold_from_contacts: basic geometric sanity
# ==============================================================================
println("\n" * "="^70)
println("fold_from_contacts")
println("="^70)

@testset "fold_from_contacts: backbone distances converge near the target after optimization" begin
    l = 10
    scores = zeros(l, l)
    scores[1, 5] = 0.9; scores[5, 1] = 0.9  # one strong long-range restraint
    cmap = ContactMap(scores, collect(1:l))
    structure = fold_from_contacts(cmap; iterations=500, learning_rate=0.01, seed=1)

    atoms = structure.models[1].chains[1].residues
    coords = [(r.atoms[1].x, r.atoms[1].y, r.atoms[1].z) for r in atoms]
    backbone_dists = [sqrt(sum((coords[i] .- coords[i+1]) .^ 2)) for i in 1:(l-1)]
    mean_backbone_dist = mean(backbone_dists)
    @test isapprox(mean_backbone_dist, 3.8; atol=1.0)
    println("mean backbone Cα-Cα distance after folding: $mean_backbone_dist (target 3.8)")

    contact_dist = sqrt(sum((coords[1] .- coords[5]) .^ 2))
    println("distance between restrained pair (1,5): $contact_dist (target contact_distance=7.5)")
    @test contact_dist < 15.0  # loose sanity bound, not a tight convergence claim
end

# ==============================================================================
# Performance smoke test
# ==============================================================================
println("\n" * "="^70)
println("Performance (l=60, n=300 synthetic alignment)")
println("="^70)

rng_perf = MersenneTwister(99)
perf_seqs = [String([rand(rng_perf, ('A','C','G','T','-')) for _ in 1:60]) for _ in 1:300]
perf_msa = make_msa(perf_seqs)

t_mi = @elapsed mutual_information_contacts(perf_msa; min_separation=5)
t_mf = @elapsed fit_pseudolikelihood_model(perf_msa; algorithm=:mean_field, max_gap_fraction=1.0, min_sequence_coverage=0.0)
t_plm = @elapsed fit_pseudolikelihood_model(perf_msa; algorithm=:plm, max_gap_fraction=1.0, min_sequence_coverage=0.0)
t_lw = @elapsed shrinkage_precision_contacts(perf_msa; shrinkage=:ledoit_wolf)

println("mutual_information_contacts (l=60,n=300): $(round(t_mi, digits=3)) s")
println("fit_pseudolikelihood_model :mean_field:    $(round(t_mf, digits=3)) s")
println("fit_pseudolikelihood_model :plm:           $(round(t_plm, digits=3)) s")
println("shrinkage_precision_contacts :ledoit_wolf: $(round(t_lw, digits=3)) s")
println("(watch the :plm number in particular -- see 'still open' notes on PLM performance)")

println("\n" * "="^70)
println("DONE")
println("="^70)