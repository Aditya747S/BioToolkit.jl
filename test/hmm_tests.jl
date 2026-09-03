using Test
using BioToolkit

@testset "Hidden Markov Models" begin
    @testset "Fair / Loaded Casino Die" begin
        states = ["Fair", "Loaded"]
        alphabet = UInt8[UInt8('1'), UInt8('2'), UInt8('3'), UInt8('4'), UInt8('5'), UInt8('6')]
        
        initial = [0.5, 0.5]
        transitions = [
            0.95  0.05;
            0.10  0.90
        ]
        
        emissions = [
            1/6  1/6  1/6  1/6  1/6  1/6;
            0.1  0.1  0.1  0.1  0.1  0.5
        ]
        
        hmm = HMM(states, alphabet, initial, transitions, emissions; log_space=false)
        
        # Test 1: Mostly normal rolls (Fair)
        seq_fair = "123145213612"
        path_f, p_f = viterbi(hmm, seq_fair)
        @test length(path_f) == length(seq_fair)
        @test all(==(1), path_f) # Should all be Fair
        
        # Test 2: Sequence of 6s (Loaded)
        seq_loaded = "666666666666"
        path_l, p_l = viterbi(hmm, seq_loaded)
        @test all(==(2), path_l) # Should all be Loaded
        
        # Test 3: Transitions
        seq_mixed = "12341236666666662312"
        path_m, p_m = viterbi(hmm, seq_mixed)
        # Start Fair
        @test path_m[1:4] == [1, 1, 1, 1]
        # Middle Loaded
        @test all(==(2), path_m[8:16])
        # End Fair
        @test path_m[19:20] == [1, 1]
        
        # Forward Probability check
        f_mixed = forward(hmm, seq_mixed)
        # Total forward probability should sum all paths; should be strictly >= Viterbi path
        @test f_mixed >= p_m
        if !isinf(f_mixed)
            @test isfinite(f_mixed)
        end
        
        # Backward Probability array check
        b_array = backward(hmm, seq_mixed)
        @test size(b_array) == (2, length(seq_mixed))
        # End state should be 0.0 (log(1.0))
        @test b_array[1, end] == 0.0
        @test b_array[2, end] == 0.0
    end
end

@testset "HMM numerical and edge-case contracts" begin
    alphabet = UInt8['A', 'B']
    hmm = HMM(["s1", "s2"], alphabet,
        [1.0, 0.0],
        [1.0 0.0; 0.25 0.75],
        [1.0 0.0; 0.5 0.5])

    # Structural zeros must remain impossible; they must not be replaced by
    # eps-sized probability mass.
    @test hmm.initial[2] == -Inf
    @test hmm.transitions[1, 2] == -Inf
    @test hmm.emissions[1, 2] == -Inf
    @test forward(hmm, "AAAA") == 0.0
    @test_throws ArgumentError viterbi(hmm, "B")
    @test_throws ArgumentError posterior_state_probabilities(hmm, "B")

    # Log-space inputs are normalized row-wise without losing -Inf entries.
    log_hmm = HMM(["s1", "s2"], alphabet,
        log.([2.0, 1.0]),
        log.([2.0 0.0; 1.0 3.0]),
        log.([3.0 1.0; 1.0 1.0]); log_space=true)
    @test all(isapprox.(exp.(log_hmm.initial), [2 / 3, 1 / 3]))
    @test all(isapprox.(vec(sum(exp.(log_hmm.transitions), dims=2)), [1.0, 1.0]))

    # Long sequences remain finite in log-space and posterior columns sum to 1.
    long_hmm = HMM(["s1", "s2"], alphabet, [0.5, 0.5],
        [0.99 0.01; 0.02 0.98], [0.8 0.2; 0.3 0.7])
    ll = forward(long_hmm, repeat("AB", 5000))
    @test isfinite(ll)
    gamma = posterior_state_probabilities(long_hmm, "ABAB")
    @test all(isapprox.(vec(sum(gamma, dims=1)), 1.0; atol=1e-10))
    gamma_fb, ll_fb = forward_backward(long_hmm, "ABAB")
    @test gamma_fb ≈ gamma
    @test ll_fb == forward(long_hmm, "ABAB")
    @test isempty(viterbi(long_hmm, "")[1])

    # HMM utilities & interface methods: size, copy, nparams, statdists, rand, joint_loglikelihood
    @test size(long_hmm) == (2, 2)
    @test size(long_hmm, 1) == 2
    @test nparams(long_hmm) == 1 + 2*(1) + 2*(1) # (2-1) + 2*1 + 2*1 = 5
    
    hmm_copied = copy(long_hmm)
    @test hmm_copied.states == long_hmm.states
    @test hmm_copied.initial == long_hmm.initial
    
    sd = statdists(long_hmm)
    @test !isempty(sd)
    @test sum(sd[1]) ≈ 1.0

    # Simulation and Joint Likelihood
    z_sim, obs_sim = rand(long_hmm, 50; seq=true)
    @test length(z_sim) == 50
    @test length(obs_sim) == 50
    @test all(b -> b in long_hmm.alphabet, obs_sim)

    jll = joint_loglikelihood(long_hmm, obs_sim, z_sim)
    @test isfinite(jll)
    @test jll <= forward(long_hmm, obs_sim) + 1e-10

    # Model Selection & Entropy
    seqs = ["ABAB", "BABA", "AAAA"]
    @test isfinite(aic(long_hmm, seqs))
    @test isfinite(bic(long_hmm, seqs))
    ent = state_entropy(gamma)
    @test length(ent) == 4
    @test all(e -> e >= 0, ent)

    # StatsAPI Aliases
    @test logdensityof(long_hmm, "ABAB") == forward(long_hmm, "ABAB")
    fit_lls = fit!(copy(long_hmm), seqs; max_iter=5)
    @test length(fit_lls) >= 1
end

@testset "GaussianHMM continuous observation tests" begin
    ghmm = GaussianHMM(["Background", "Peak"], [0.8, 0.2],
                       [0.9 0.1; 0.2 0.8], [0.0, 10.0], [1.0, 2.0])
    obs_signal = [0.1, -0.2, 0.3, 9.8, 10.2, 9.5, 0.0]
    path, logp = viterbi(ghmm, obs_signal)
    @test length(path) == 7
    @test path[4:6] == [2, 2, 2] # Peak state
    @test isfinite(forward(ghmm, obs_signal))
    
    post = posterior_state_probabilities(ghmm, obs_signal)
    @test size(post) == (2, 7)
    
    bw_lls = baum_welch!(ghmm, [obs_signal]; max_iter=5)
    @test length(bw_lls) >= 1
end

@testset "ChromHMM multi-track epigenomic state tests" begin
    chmm = ChromHMM(["Promoter", "Enhancer", "Heterochromatin"],
                    ["H3K4me3", "H3K27ac", "H3K27me3"],
                    [0.33, 0.33, 0.34],
                    [0.8 0.1 0.1; 0.1 0.8 0.1; 0.1 0.1 0.8],
                    [0.9 0.8 0.01;   # Promoter: high H3K4me3, H3K27ac
                     0.1 0.9 0.01;   # Enhancer: high H3K27ac
                     0.01 0.01 0.9]) # Heterochromatin: high H3K27me3

    # Mark matrix: 3 marks x 5 genomic bins
    marks = [1.0 1.0 0.0 0.0 0.0;
             1.0 1.0 1.0 0.0 0.0;
             0.0 0.0 0.0 1.0 1.0]

    path, logp = viterbi(chmm, marks)
    @test length(path) == 5
    @test path[1] == 1 # Promoter
    @test path[4] == 3 # Heterochromatin

    segs = segment_chromatin(chmm, marks)
    @test !isempty(segs)
    @test segs[1].state_name == "Promoter"
end

@testset "GenericEmissionHMM & Blender payload tests" begin
    # Custom emission functions (e.g. Poisson log-pdf for sequencing counts)
    log_factorial(k) = k <= 1 ? 0.0 : sum(log(i) for i in 1:k)
    poisson_logpdf(λ, k) = k * log(λ) - λ - log_factorial(k)
    
    gem = GenericEmissionHMM(["LowReadCount", "HighReadCount"],
                             [0.5, 0.5],
                             [0.9 0.1; 0.1 0.9],
                             [k -> poisson_logpdf(2.0, k),
                              k -> poisson_logpdf(50.0, k)])
    
    @test n_states(gem) == 2
    @test states(gem) == ["LowReadCount", "HighReadCount"]
    
    read_counts = [1, 2, 0, 52, 48, 55, 1]
    post = posterior_state_probabilities(gem, read_counts)
    @test size(post) == (2, 7)
    @test post[2, 4] > 0.9 # Position 4 (count=52) is HighReadCount state
    
    # 3D Blender Payload creation test
    long_hmm = HMM(["s1", "s2"], UInt8['A', 'B'], [0.5, 0.5],
                   [0.99 0.01; 0.02 0.98], [0.8 0.2; 0.3 0.7])
    payload = to_blender_payload(long_hmm, "ABAB")
    @test payload isa BioToolkit.BlenderIntegrator.BlenderHMMPayload
    @test payload.name == "HMM_Landscape"
    @test size(payload.posterior_matrix) == (2, 4)
end
