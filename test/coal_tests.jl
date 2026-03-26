using BioToolkit
using Test

@testset "Coalescent Simulation" begin
    n = 5
    ne = 1000
    mu = 0.001
    seq_len = 50
    
    tree, haps = BioToolkit.simulate_coalescent(n; ne=ne, mu=mu, seq_len=seq_len)
    
    @test count_terminals(tree) == 5
    @test length(haps) == 5
    @test all(length(h) == seq_len for h in values(haps))
    
    # Check if we have some variation (it's random but with mu=0.001 is likely)
    all_identical = true
    first_hap = first(values(haps))
    for h in values(haps)
        if h != first_hap
            all_identical = false
            break
        end
    end
    println("Variation produced: ", !all_identical)
    
    # Test large simulation
    println("Running larger simulation (n=50)...")
    tree2, haps2 = BioToolkit.simulate_coalescent(50; ne=10000, mu=1e-8, seq_len=1000)
    @test count_terminals(tree2) == 50
    println("Done.")
end
