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
        @test length(g1) == 1
        @test only(values(g1)) == 1.0
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

@testset "PopGen correctness regressions" begin
    missing_pop = Population{Int}("missing", [
        PopGenIndividual{Int}("a", [Locus{Int}((1, 2))]),
        PopGenIndividual{Int}("b", [Locus{Int}((0, 0))]),
        PopGenIndividual{Int}("c", [Locus{Int}((2, 1))])
    ])
    frequencies = allele_frequencies(missing_pop, 1)
    @test frequencies == Dict(1 => 0.5, 2 => 0.5)
    @test heterozygosity_observed(missing_pop, 1) == 1.0
    @test_throws ArgumentError allele_frequencies(missing_pop, 0)

    # Exact HWE remains finite for samples where direct probability products underflow.
    large_hwe = Population{Int}("large", vcat(
        [PopGenIndividual{Int}("aa", [Locus{Int}((1, 1))]) for _ in 1:600],
        [PopGenIndividual{Int}("bb", [Locus{Int}((2, 2))]) for _ in 1:600]))
    @test 0.0 <= hardy_weinberg_exact(large_hwe, 1) <= 1.0

    unphased = Population{Int}("unphased", [
        PopGenIndividual{Int}("u$(i)", [Locus{Int}((1, 2)), Locus{Int}((1, 2))]) for i in 1:20])
    @test linkage_disequilibrium(unphased, 1, 2)[3] < 1e-10

    p_a = Population{Int}("a", [PopGenIndividual{Int}("a", [Locus{Int}((1, 1))])])
    p_b = Population{Int}("b", [PopGenIndividual{Int}("b", [Locus{Int}((2, 2))])])
    @test genetic_distance(p_a, p_b, 1, method=:rogers) == 1.0
    @test genetic_distance(p_a, p_b, 1, method=:nei) == Inf
    @test_throws ArgumentError genetic_distance(p_a, p_b, 1, method=:unknown)
    @test_throws ArgumentError population_pcoa([0.0 1.0; 2.0 0.0])
end

@testset "Population workflow parity" begin
    pa = Population{Int}("A", [
        PopGenIndividual{Int}("a1", [Locus{Int}((1, 1)), Locus{Int}((1, 2))]),
        PopGenIndividual{Int}("a2", [Locus{Int}((1, 1)), Locus{Int}((0, 0))])
    ])
    pb = Population{Int}("B", [
        PopGenIndividual{Int}("b1", [Locus{Int}((2, 2)), Locus{Int}((2, 2))]),
        PopGenIndividual{Int}("b2", [Locus{Int}((2, 2)), Locus{Int}((2, 2))])
    ])
    data = [pa, pb]
    @test populations(data) == ["A", "B"]
    @test populations(data; counts=true) == Dict("A" => 2, "B" => 2)
    @test loci(data) == [1, 2]
    @test samplenames(data) == ["a1", "a2", "b1", "b2"]
    @test missingdata(data; by=:sample)[2].missing == 1
    @test missingdata(data; by=:population)[1].missing == 1
    @test richness(data; by=:locus)[1].richness == 2
    @test richness(data; by=:population)[1].richness == 1
    @test alleleaverage(data).mean == 2.0
    fst = pairwise_fst(data)
    @test fst.populations == ["A", "B"]
    @test fst.estimates[1, 2] > 0.5
    @test fst.loci_used[1, 2] == 2
    stats = summary_statistics(data)
    @test stats.n_loci == 2
    @test 0.0 <= stats.fst <= 1.0
    @test summary(data) == stats
end

@testset "PopGen.jl Parity Features" begin
    pa = Population{Int}("A", [
        PopGenIndividual{Int}("a1", [Locus{Int}((1, 1)), Locus{Int}((1, 2))]),
        PopGenIndividual{Int}("a2", [Locus{Int}((1, 1)), Locus{Int}((1, 1))])
    ])
    pb = Population{Int}("B", [
        PopGenIndividual{Int}("b1", [Locus{Int}((2, 2)), Locus{Int}((2, 2))]),
        PopGenIndividual{Int}("b2", [Locus{Int}((2, 2)), Locus{Int}((2, 2))])
    ])
    data = [pa, pb]

    # 1. Hudson & Weir-Cockerham single locus
    h_fst = hudson_fst(pa, pb, 1)
    @test !isnan(h_fst) && h_fst >= 0.0

    wc_fst = weir_cockerham_fst(pa, pb, 1)
    @test !isnan(wc_fst) && wc_fst >= 0.0

    # 2. Pairwise FST with all methods
    fst_hudson = pairwise_fst(data; method=:hudson)
    @test fst_hudson.method === :hudson
    @test fst_hudson.estimates[1, 2] > 0.5

    fst_wc = pairwise_fst(data; method=:weir_cockerham)
    @test fst_wc.method === :weir_cockerham
    @test fst_wc.estimates[1, 2] > 0.5

    fst_amova = pairwise_fst(data; method=:amova)
    @test fst_amova.method === :amova

    # 3. Permutation test
    perm_res = fst_permutation_test(data; method=:nei, permutations=19)
    @test size(perm_res.pvalues) == (2, 2)
    @test perm_res.pvalues[1, 1] == 0.0

    # 4. Sample heterozygosity
    sample_ho = sample_heterozygosity(data)
    @test length(sample_ho) == 4
    @test sample_ho[1].sample == "a1"
    @test sample_ho[1].ho == 0.5

    # 5. Kinship estimators
    for m in [:queller_goodnight, :blouin, :li_horvitz, :ritland, :lynch_li, :lynch_ritland, :moran, :loiselle]
        kin = pairwise_kinship(data; method=m)
        @test size(kin.matrix) == (4, 4)
        @test kin.matrix[1, 1] == 1.0
    end

    # 6. Clustering suite
    km = population_kmeans(data; k=2)
    @test length(km.assignments) == 2

    hc = population_hclust(data)
    @test length(hc.linkage) == 1

    kmed = population_kmedoids(data; k=2)
    @test length(kmed.assignments) == 2

    fcm = population_fuzzycmeans(data; k=2)
    @test size(fcm.membership) == (2, 2)

    db = population_dbscan(data; eps=1.0)
    @test length(db.assignments) == 2

    for clus_m in [:kmeans, :kmedoids, :hclust, :fuzzycmeans, :dbscan]
        res = population_cluster(data; method=clus_m)
        @test res !== nothing
    end

    # 7. Pairwise identical genotypes
    pi_res = pairwise_identical(data)
    @test size(pi_res.matrix) == (4, 4)
    @test pi_res.matrix[1, 1] == 1.0

    # 8. Multi-pop HWE
    hwe_dict = hardy_weinberg_test(data, 1)
    @test haskey(hwe_dict, "A") && haskey(hwe_dict, "B")

    hwe_exact_dict = hardy_weinberg_exact(data, 1)
    @test haskey(hwe_exact_dict, "A") && haskey(hwe_exact_dict, "B")

    # 9. STRUCTURE format read/write
    tmp_str, str_io = mktemp()
    close(str_io)
    try
        write_structure(tmp_str, data)
        read_pops = read_structure_record(tmp_str)
        @test length(read_pops) >= 1
    finally
        isfile(tmp_str) && rm(tmp_str)
    end
end

println("All PopGen tests passed successfully!")
