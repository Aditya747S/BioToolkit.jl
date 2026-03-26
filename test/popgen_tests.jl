using BioToolkit
using Test
using Statistics
using Random
using LinearAlgebra

@testset "Basic PopGen Features" begin
    # Setup
    l1 = Locus{Int}((1, 2))
    l2 = Locus{Int}((1, 1))
    ind1 = PopGenIndividual{Int}("i1", [l1, l2])
    ind2 = PopGenIndividual{Int}("i2", [l1, l1])
    pop = Population{Int}("TestPop", [ind1, ind2])

    @testset "Frequencies" begin
        f1 = allele_frequencies(pop, 1)
        @test f1[1] == 0.5
        @test f1[2] == 0.5
        
        g1 = genotype_frequencies(pop, 1)
        @test g1[Set([1, 2])] == 1.0
    end

    @testset "Heterozygosity" begin
        @test heterozygosity_observed(pop, 1) == 1.0
        @test heterozygosity_expected(pop, 1) == 0.5
    end

    @testset "HWE" begin
        # Perfect HWE: 25 AA, 50 AB, 25 BB
        inds = []
        for _ in 1:25; push!(inds, PopGenIndividual{Int}("A", [Locus{Int}((1, 1))])); end
        for _ in 1:50; push!(inds, PopGenIndividual{Int}("B", [Locus{Int}((1, 2))])); end
        for _ in 1:25; push!(inds, PopGenIndividual{Int}("C", [Locus{Int}((2, 2))])); end
        pop_hwe = Population{Int}("HWE", inds)
        
        p_val = hardy_weinberg_test(pop_hwe, 1)
        @test p_val > 0.9
    end
end

@testset "AMOVA & Differentiation" begin
    # 3 individuals in 2 populations
    # pop1: (1,1), (1,1)
    # pop2: (2,2)
    p1 = Population{Int}("p1", [PopGenIndividual{Int}("i1", [Locus{Int}((1, 1))]), 
                               PopGenIndividual{Int}("i2", [Locus{Int}((1, 1))])])
    p2 = Population{Int}("p2", [PopGenIndividual{Int}("i3", [Locus{Int}((2, 2))])])
    
    # Distance matrix (arbitrary for test)
    # i1-i2: 0, i1-i3: 1, i2-i3: 1
    dist = [0.0 0.0 1.0; 0.0 0.0 1.0; 1.0 1.0 0.0]
    phi, vw, va = amova(dist, [2, 1])
    @test phi > 0.5
    
    # FST
    f_stats = f_statistics([p1, p2], 1)
    @test f_stats[2] > 0.5 # High FST
end

@testset "Advanced PopGen Features" begin
    # 1. HWE Exact
    l_hwe = Locus{Int}((1, 2))
    ind_hwe = PopGenIndividual{Int}("h1", [l_hwe])
    pop_hwe = Population{Int}("PHWE", [ind_hwe, ind_hwe])
    @test hardy_weinberg_exact(pop_hwe, 1) > 0.0
    
    # Ne LD
    l1 = Locus{Int}((1, 1))
    l2 = Locus{Int}((1, 2))
    ind1 = PopGenIndividual{Int}("i1", [l1, l2])
    pop_ne = Population{Int}("PNE", [ind1, ind1, ind1, ind1])
    ne_ld = estimate_ne_ld(pop_ne, [(1, 2)])
    @test ne_ld > 0.0
    
    # LD Mapping
    ld_map = ld_mapping(pop_ne, [1, 2], 1)
    @test length(ld_map) == 1
    @test haskey(ld_map, (1, 2))
    
    # 2. Migration & G-Stats
    ind1 = PopGenIndividual{Int}("p1", [Locus{Int}((1, 1))])
    ind2 = PopGenIndividual{Int}("p2", [Locus{Int}((2, 2))])
    pop1 = Population{Int}("P1", [ind1, ind1])
    pop2 = Population{Int}("P2", [ind2, ind2])
    @test migration_rate([pop1, pop2], 1) == 0.0
    
    g_stats = g_statistics([pop1, pop2], 1)
    @test g_stats[1] == 1.0 # GST = 1.0 for completely fixed differences
    
    # 3. Genetic Distance
    d_rogers = genetic_distance(pop1, pop2, 1, method=:rogers)
    @test d_rogers == 1.0 # sqrt( (1-0)^2 + (0-1)^2 ) / sqrt(2) = sqrt(2)/sqrt(2) = 1.0
    
    # 4. Spatial & PCA/PCoA
    d1 = [0.0 1.0 2.0; 1.0 0.0 1.5; 2.0 1.5 0.0]
    m1 = [0.0 1.0 2.0; 1.0 0.0 1.5; 2.0 1.5 0.0]
    m2 = [0.0 2.0 4.0; 2.0 0.0 3.0; 4.0 3.0 0.0]
    corr, p = mantel_test(m1, m2, permutations=19)
    @test corr > 0.9
    @test p <= 1.0 # Just test it returns a valid p-value
    
    # PCoA
    coords, evals = population_pcoa(d1)
    @test size(coords, 1) == 3
    @test length(evals) >= 1
    
    # 5. Mismatch Distribution
    s1 = BioToolkit.SeqRecordLite("AAAAA", identifier="s1")
    s2 = BioToolkit.SeqRecordLite("AAATT", identifier="s2")
    msa = BioToolkit.MultipleSequenceAlignment([s1, s2])
    mm = mismatch_distribution(msa)
    @test mm == [2]
end

@testset "GenePop Record Helpers" begin
    tmp_path, tmp_io = mktemp()
    close(tmp_io)
    try
        open(tmp_path, "w") do io
            println(io, "BioToolkit test")
            println(io, "Locus_1")
            println(io, "Locus_2")
            println(io, "Pop")
            println(io, "Ind1, 001001 001002")
            println(io, "Ind2, 002002 002001")
            println(io, "Pop")
            println(io, "Ind3, 001002 001001")
        end

        record = read_genepop_record(tmp_path)
        @test record.comment_line == "BioToolkit test"
        @test record.loci_list == ["Locus_1", "Locus_2"]
        @test length(record.populations) == 2

        pop_splits = split_in_pops(record, ["A", "B"])
        @test sort(collect(keys(pop_splits))) == ["A", "B"]
        @test length(pop_splits["A"].populations) == 1

        locus_splits = split_in_loci(record)
        @test haskey(locus_splits, "Locus_1")
        @test length(locus_splits["Locus_1"].loci_list) == 1

        working = deepcopy(record)
        remove_population!(working, 0)
        @test length(working.populations) == 1

        working = deepcopy(record)
        remove_locus_by_position!(working, 0)
        @test working.loci_list == ["Locus_2"]

        working = deepcopy(record)
        remove_locus_by_name!(working, "Locus_2")
        @test working.loci_list == ["Locus_1"]

        @test occursin("Pop", sprint(show, record))
    finally
        isfile(tmp_path) && rm(tmp_path)
    end
end

@testset "Advanced Neutrality & Sweeps" begin
    s1 = BioToolkit.SeqRecordLite("AAAAA", identifier="s1")
    s2 = BioToolkit.SeqRecordLite("AAAAT", identifier="s2")
    s3 = BioToolkit.SeqRecordLite("AAATT", identifier="s3")
    msa = BioToolkit.MultipleSequenceAlignment([s1, s2, s3])
    
    # Ewens-Watterson
    l1 = Locus{Int}((1, 1))
    ind1 = PopGenIndividual{Int}("i1", [l1])
    pop = Population{Int}("P", [ind1, ind1, ind1])
    obs, exp_f = ewens_watterson_test(pop, 1)
    @test obs == 1.0
    
    # SweepFinder CLR
    clrs = sweepfinder_clr(msa, 10)
    @test length(clrs) == 10
    @test all(clrs .>= 0.0)
end

@testset "Kinship & GWAS" begin
    # Mock data for GRM
    # ind1: 1/1 (AA -> code 2), ind2: 1/2 (AB -> code 1)
    l1 = Locus{Int}((1, 1))
    l2 = Locus{Int}((1, 2))
    ind1 = PopGenIndividual{Int}("i1", [l1])
    ind2 = PopGenIndividual{Int}("i2", [l2])
    pop = Population{Int}("P", [ind1, ind2])
    
    G = genetic_relationship_matrix([pop])
    @test size(G) == (2, 2)
    @test G[1, 1] > 0.0
    
    # inbreeding
    f = inbreeding_coefficient(G, 1)
    @test f == G[1, 1] - 1.0
    
    # GWAS
    # 2 inds, 1 locus
    genos = [2.0; 1.0;;] # 2x1 matrix
    phenos = [10.0, 5.0]
    p_vals = linear_mixed_model_scan(genos, phenos, G)
    @test length(p_vals) == 1
    @test 0.0 <= p_vals[1] <= 1.0
end

println("All PopGen tests passed successfully!")
