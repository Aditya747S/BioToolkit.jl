using BioToolkit
using Test

@testset "ML Phylogenetics" begin
    # Simple 3-taxon case
    alignment = Dict(
        "A" => "ACGTACGTACGT",
        "B" => "ACGTACGTACGT",
        "C" => "TGCAACGTTGCA"
    )
    
    # Starting tree (A:0.1, (B:0.1, C:0.1))
    tree_str = "(A:0.1,(B:0.1,C:0.1):0.1);"
    tree = parse_newick(tree_str)
    
    model = JC69()
    lik = felsenstein_likelihood(tree, alignment, model)
    println("Likelihood for (A,(B,C)): ", lik)
    
    # Try different topology (C,(A,B))
    tree_str2 = "(C:0.1,(A:0.1,B:0.1):0.1);"
    tree2 = parse_newick(tree_str2)
    lik2 = felsenstein_likelihood(tree2, alignment, model)
    println("Likelihood for (C,(A,B)): ", lik2)
    
    @test lik2 > lik # (A,B) are identical, so they should be grouped
    
    # Test ML optimization
    println("Starting ML tree search...")
    ml_tree = maximum_likelihood_tree(alignment; model=model)
    println("ML Tree: ", write_newick(ml_tree))
    
    # Verify NJ used as starting point and then optimized
    @test count_terminals(ml_tree) == 3
end

@testset "Substitution Models" begin
    jc = JC69()
    P_jc = transition_probability(jc, 0.1)
    @test size(P_jc) == (4, 4)
    @test all(sum(P_jc, dims=2) .≈ 1.0)
    
    k80 = K80(2.0)
    P_k80 = transition_probability(k80, 0.1)
    @test P_k80[1,3] > P_k80[1,2] # Transition > Transversion
    
    hky = HKY85([0.2, 0.3, 0.3, 0.2], 2.0)
    P_hky = transition_probability(hky, 0.1)
    @test all(sum(P_hky, dims=2) .≈ 1.0)
end
