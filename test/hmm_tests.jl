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
