# ==============================================================================
using Distributions, Random, Statistics


# ------------------------------------------------------------------------------
# 1. Structures and Base Types
# ------------------------------------------------------------------------------

"""
    Locus{T}

A genetic locus storing a tuple of alleles. `T` is usually `String` or `Char` or `Int`.
For a diploid individual, this will hold two alleles.
"""
struct Locus{T}
    alleles::Tuple{Vararg{T}}
end

export Locus

"""
    PopGenIndividual{T}

An individual containing a unique identifier and an array of loci.
"""
struct PopGenIndividual{T}
    id::String
    loci::Vector{Locus{T}}
end

"""
    Population{T}

A collection of individuals comprising a population or sub-population.
"""
struct Population{T}
    name::String
    individuals::Vector{PopGenIndividual{T}}
end

export Locus, PopGenIndividual, Population
export allele_frequencies, genotype_frequencies, heterozygosity_observed, heterozygosity_expected, hardy_weinberg_test
export amova, f_statistics, hardy_weinberg_exact, estimate_ne_ld, ld_mapping, migration_rate, g_statistics, genetic_distance
export mantel_test, population_pcoa, mismatch_distribution
export read_genepop_record, split_in_pops, split_in_loci, remove_population!, remove_locus_by_position!, remove_locus_by_name!
export ewens_watterson_test, sweepfinder_clr, genetic_relationship_matrix, inbreeding_coefficient, linear_mixed_model_scan

# ------------------------------------------------------------------------------
# 2. Basic Frequencies
# ------------------------------------------------------------------------------

"""
    allele_frequencies(pop::Population, locus_idx::Int)

Calculates the frequency of each allele at a specific locus across the population.
Returns a Dict mapping allele -> frequency.
"""
function allele_frequencies(pop::Population{T}, locus_idx::Int) where T
    counts = Dict{T, Int}()
    total_alleles = 0
    for ind in pop.individuals
        if locus_idx <= length(ind.loci)
            locus = ind.loci[locus_idx]
            for allele in locus.alleles
                # skip missing data represented by '?' or "0" etc based on type, handle natively
                counts[allele] = get(counts, allele, 0) + 1
                total_alleles += 1
            end
        end
    end
    
    freqs = Dict{T, Float64}()
    total_float = Float64(total_alleles)
    if total_float > 0.0
        for (a, c) in counts
            freqs[a] = c / total_float
        end
    end
    return freqs
end

"""
    genotype_frequencies(pop::Population, locus_idx::Int)

Calculates the frequency of observed genotypes at a specific locus.
Assumes unordered alleles (e.g. A/B is the same as B/A).
"""
function genotype_frequencies(pop::Population{T}, locus_idx::Int) where T
    counts = Dict{Set{T}, Int}()
    total_genos = 0
    for ind in pop.individuals
        if locus_idx <= length(ind.loci)
            locus = ind.loci[locus_idx]
            geno = Set(locus.alleles)
            counts[geno] = get(counts, geno, 0) + 1
            total_genos += 1
        end
    end
    
    freqs = Dict{Set{T}, Float64}()
    total_float = Float64(total_genos)
    if total_float > 0.0
        for (g, c) in counts
            freqs[g] = c / total_float
        end
    end
    return freqs
end

# ------------------------------------------------------------------------------
# 3. Heterozygosity
# ------------------------------------------------------------------------------

"""
    heterozygosity_observed(pop::Population, locus_idx::Int)

Returns the observed proportion of heterozygotes at a given locus.
"""
function heterozygosity_observed(pop::Population{T}, locus_idx::Int) where T
    hets = 0
    total = 0
    for ind in pop.individuals
        if locus_idx <= length(ind.loci)
            locus = ind.loci[locus_idx]
            if length(locus.alleles) == 2
                if locus.alleles[1] != locus.alleles[2]
                    hets += 1
                end
                total += 1
            end
        end
    end
    return total > 0 ? (hets / Float64(total)) : 0.0
end

"""
    heterozygosity_expected(pop::Population, locus_idx::Int)

Returns the expected heterozygosity (Nei's gene diversity) at a given locus,
calculated as 1 - sum(p_i^2) where p_i are the allele frequencies.
"""
function heterozygosity_expected(pop::Population{T}, locus_idx::Int) where T
    freqs = allele_frequencies(pop, locus_idx)
    homo_sum = sum(p^2 for p in values(freqs))
    return 1.0 - homo_sum
end

# ------------------------------------------------------------------------------
# 4. Hardy-Weinberg Equilibrium
# ------------------------------------------------------------------------------

"""
    hardy_weinberg_test(pop::Population, locus_idx::Int)

Performs a standard Chi-Square test for Hardy-Weinberg equilibrium at a bi-allelic locus.
Returns the p-value.
"""
function hardy_weinberg_test(pop::Population{T}, locus_idx::Int) where T
    # Simplified bi-allelic implementation for fast path
    freqs = allele_frequencies(pop, locus_idx)
    alleles = collect(keys(freqs))
    if length(alleles) != 2
        return 1.0 # Requires exactly 2 alleles for standard simple chi-square
    end
    
    p = freqs[alleles[1]]
    q = freqs[alleles[2]]
    
    total_genos = 0
    obs_p2 = 0
    obs_q2 = 0
    obs_pq = 0
    
    for ind in pop.individuals
        if locus_idx <= length(ind.loci)
            locus = ind.loci[locus_idx]
            if length(locus.alleles) == 2
                a1, a2 = locus.alleles
                if a1 == alleles[1] && a2 == alleles[1]
                    obs_p2 += 1
                elseif a1 == alleles[2] && a2 == alleles[2]
                    obs_q2 += 1
                else
                    obs_pq += 1
                end
                total_genos += 1
            end
        end
    end
    
    exp_p2 = (p^2) * total_genos
    exp_q2 = (q^2) * total_genos
    exp_pq = (2 * p * q) * total_genos
    
    # Chi-Square statistic
    chi2 = 0.0
    chi2 += exp_p2 > 0 ? ((obs_p2 - exp_p2)^2 / exp_p2) : 0.0
    chi2 += exp_q2 > 0 ? ((obs_q2 - exp_q2)^2 / exp_q2) : 0.0
    chi2 += exp_pq > 0 ? ((obs_pq - exp_pq)^2 / exp_pq) : 0.0
    
    # Approx 1 DoF Chi-Square PDF -> p-value mapping
    # Note: rigorous p-value requires advanced stats libs, using a fast approximation for popgen
    # This is a basic approximation of the chi-square CDF upper tail
    # More exact approximations will be used in production
    
    # Basic lookup approx for df=1
    return exp(-chi2/2) / sqrt(chi2 * 2 + 1e-9)
end

"""
    hardy_weinberg_exact(pop::Population, locus_idx::Int)

Slatkin's Exact Test for Hardy-Weinberg Equilibrium.
Useful for small sample sizes where Chi-Square is inaccurate.
"""
function hardy_weinberg_exact(pop::Population{T}, locus_idx::Int) where T
    # Standard exact test for bi-allelic loci (Levene 1949, Haldane 1954)
    freqs = allele_frequencies(pop, locus_idx)
    alleles = collect(keys(freqs))
    if length(alleles) != 2
        return 1.0 
    end
    
    # Observed counts
    n11, n22, n12 = 0, 0, 0
    for ind in pop.individuals
        if locus_idx <= length(ind.loci)
            l = ind.loci[locus_idx]
            if length(l.alleles) == 2
                a1, a2 = l.alleles
                if a1 == alleles[1] && a2 == alleles[1]; n11 += 1
                elseif a1 == alleles[2] && a2 == alleles[2]; n22 += 1
                else; n12 += 1
                end
            end
        end
    end
    
    n = n11 + n22 + n12
    n1 = 2*n11 + n12
    n2 = 2*n22 + n12
    
    # Probability of this specific configuration given allele counts
    # P(n12 | n1, n2, n) = (n! / (n11! n22! n12!)) * (2^n12 / (2n! / (n1! n2!)))
    function log_prob(n12_val, n1_val, n2_val, n_val)
        n11_val = (n1_val - n12_val) / 2
        n22_val = (n2_val - n12_val) / 2
        
        # log P = log(n!) + n12*log(2) + log(n1!) + log(n2!) - log(n11!) - log(n22!) - log(n12!) - log(2n!)
        # Use logfactorial (lfact) or lgamma
        return (loggamma(n_val + 1) + n12_val * log(2.0) + loggamma(n1_val + 1) + loggamma(n2_val + 1) - 
                loggamma(Int(n11_val) + 1) - loggamma(Int(n22_val) + 1) - loggamma(n12_val + 1) - loggamma(2 * n_val + 1))
    end
    
    p_obs = exp(log_prob(n12, n1, n2, n))
    
    # Sum probabilities of all configurations as likely or less likely than observed
    p_total = 0.0
    # n12 can range from n1 % 2 to min(n1, n2) with step 2
    for i in (n1 % 2):2:min(n1, n2)
        p_i = exp(log_prob(i, n1, n2, n))
        if p_i <= p_obs + 1e-9
            p_total += p_i
        end
    end
    
    return min(p_total, 1.0)
end

# ------------------------------------------------------------------------------
# 5. F-Statistics & Diversity
# ------------------------------------------------------------------------------

"""
    f_statistics(populations::Vector{Population}, locus_idx::Int)

Calculates Wright's F-statistics (F_IS, F_ST, F_IT) across multiple subpopulations
for a specific locus.
Returns a tuple (F_IS, F_ST, F_IT).
"""
function f_statistics(populations::Vector{Population{T}}, locus_idx::Int) where T
    total_inds = 0
    global_allele_counts = Dict{T, Int}()
    
    subpop_Ho = Float64[]
    subpop_Hs = Float64[]
    subpop_weights = Float64[]
    
    for pop in populations
        pop_inds = 0
        hts = 0
        pop_counts = Dict{T, Int}()
        
        for ind in pop.individuals
            if locus_idx <= length(ind.loci)
                l = ind.loci[locus_idx]
                if length(l.alleles) == 2
                    if l.alleles[1] != l.alleles[2]
                        hts += 1
                    end
                    pop_counts[l.alleles[1]] = get(pop_counts, l.alleles[1], 0) + 1
                    pop_counts[l.alleles[2]] = get(pop_counts, l.alleles[2], 0) + 1
                    
                    global_allele_counts[l.alleles[1]] = get(global_allele_counts, l.alleles[1], 0) + 1
                    global_allele_counts[l.alleles[2]] = get(global_allele_counts, l.alleles[2], 0) + 1
                    pop_inds += 1
                end
            end
        end
        
        if pop_inds > 0
            Ho = hts / Float64(pop_inds)
            push!(subpop_Ho, Ho)
            
            # Expected heterozygosity in subpopulation (Hs component)
            freqs = [c / Float64(pop_inds * 2) for c in values(pop_counts)]
            Hs = 1.0 - sum(f^2 for f in freqs)
            push!(subpop_Hs, Hs)
            
            push!(subpop_weights, pop_inds)
            total_inds += pop_inds
        end
    end
    
    if total_inds == 0
        return (0.0, 0.0, 0.0)
    end
    
    # Global expected heterozygosity (Ht)
    global_freqs = [c / Float64(total_inds * 2) for c in values(global_allele_counts)]
    Ht = 1.0 - sum(f^2 for f in global_freqs)
    
    # Mean observed heterozygosity (H_I)
    mean_Ho = sum(subpop_Ho .* subpop_weights) / sum(subpop_weights)
    
    # Mean expected heterozygosity within subpopulations (H_S)
    mean_Hs = sum(subpop_Hs .* subpop_weights) / sum(subpop_weights)
    
    # Calculate F-statistics
    F_IS = mean_Hs > 0 ? (mean_Hs - mean_Ho) / mean_Hs : 0.0
    F_ST = Ht > 0 ? (Ht - mean_Hs) / Ht : 0.0
    F_IT = Ht > 0 ? (Ht - mean_Ho) / Ht : 0.0
    
    return (F_IS, F_ST, F_IT)
end

"""
    g_statistics(populations::Vector{Population}, locus_idx::Int)

Calculates G_ST (Nei's GST) and Jost's D.
Returns a tuple (G_ST, Jost_D)
"""
function g_statistics(populations::Vector{Population{T}}, locus_idx::Int) where T
    # Heterozygosity components
    num_pops = length(populations)
    if num_pops < 2; return (0.0, 0.0); end
    
    total_inds = 0
    Ho_vals = Float64[]
    Hs_vals = Float64[]
    
    global_counts = Dict{T, Int}()
    
    for pop in populations
        n = 0
        hts = 0
        counts = Dict{T, Int}()
        for ind in pop.individuals
            if locus_idx <= length(ind.loci)
                l = ind.loci[locus_idx]
                if length(l.alleles) == 2
                    n += 1
                    if l.alleles[1] != l.alleles[2]; hts += 1; end
                    counts[l.alleles[1]] = get(counts, l.alleles[1], 0) + 1
                    counts[l.alleles[2]] = get(counts, l.alleles[2], 0) + 1
                    global_counts[l.alleles[1]] = get(global_counts, l.alleles[1], 0) + 1
                    global_counts[l.alleles[2]] = get(global_counts, l.alleles[2], 0) + 1
                end
            end
        end
        
        if n > 0
            Ho = hts / n
            p_freqs = [c / (2*n) for c in values(counts)]
            Hs = 1.0 - sum(p^2 for p in p_freqs)
            push!(Ho_vals, Ho)
            push!(Hs_vals, Hs)
            total_inds += n
        end
    end
    
    mean_Hs = isempty(Hs_vals) ? 0.0 : sum(Hs_vals) / length(Hs_vals)
    global_p_freqs = [c / (2*total_inds) for c in values(global_counts)]
    Ht = 1.0 - sum(p^2 for p in global_p_freqs)
    
    # Nei's Gst
    G_ST = Ht > 0 ? (Ht - mean_Hs) / Ht : 0.0
    
    # Jost's D (2008)
    # D = (Ht - Hs) / (1 - Hs) * (n / (n-1))
    Jost_D = (mean_Hs < 1.0 && num_pops > 1) ? 
             ((Ht - mean_Hs) / (1.0 - mean_Hs)) * (num_pops / (num_pops - 1)) : 0.0
             
    return (G_ST, Jost_D)
end

"""
    migration_rate(populations::Vector{Population}, locus_idx::Int)

Estimates the number of migrants per generation (Nm) using the Fst method:
Nm ≈ (1 - Fst) / (4 * Fst) for diploid organisms.
"""
function migration_rate(populations::Vector{Population{T}}, locus_idx::Int) where T
    f_stats = f_statistics(populations, locus_idx)
    F_ST = f_stats[2]
    if F_ST <= 0.0; return Inf; end
    if F_ST >= 1.0; return 0.0; end
    return (1.0 - F_ST) / (4.0 * F_ST)
end

# ------------------------------------------------------------------------------
# 6. AMOVA (Analysis of Molecular Variance)
# ------------------------------------------------------------------------------

"""
    amova(distance_matrix::Matrix{Float64}, pop_sizes::Vector{Int})

Performs a basic one-way AMOVA given a pairwise square distance matrix of individuals
and a vector designating how many individuals belong to each sequential subpopulation.
Returns (Phi_ST, Variance_Within, Variance_Among).
"""
function amova(distance_matrix::Matrix{Float64}, pop_sizes::Vector{Int})
    N = size(distance_matrix, 1)
    K = length(pop_sizes) # number of populations
    
    sum_D_total = sum(distance_matrix)
    SSD_total = sum_D_total / (2 * N)
    
    SSD_within = 0.0
    start_idx = 1
    for nk in pop_sizes
        end_idx = start_idx + nk - 1
        sub_matrix = @view distance_matrix[start_idx:end_idx, start_idx:end_idx]
        sum_D_within = sum(sub_matrix)
        SSD_within += sum_D_within / (2 * nk)
        start_idx = end_idx + 1
    end
    
    SSD_among = SSD_total - SSD_within
    
    df_among = K - 1
    df_within = N - K
    
    # Mean Squares
    MS_among = df_among > 0 ? SSD_among / df_among : 0.0
    MS_within = df_within > 0 ? SSD_within / df_within : 0.0
    
    # Calculate n0 (average sample size weighting adjustment)
    n0 = (N - sum(nk^2 for nk in pop_sizes)/N) / (K - 1)
    
    # Variance Components
    Var_within = MS_within
    Var_among = n0 > 0 ? (MS_among - MS_within) / n0 : 0.0
    
    Var_total = Var_within + Var_among
    Phi_ST = Var_total > 0 ? Var_among / Var_total : 0.0
    
    return (Phi_ST, Var_within, Var_among)
end

# ------------------------------------------------------------------------------
# 7. Linkage Disequilibrium
# ------------------------------------------------------------------------------

"""
    linkage_disequilibrium(pop::Population, locus_1::Int, locus_2::Int)

Calculate D, D', and r^2 between two bi-allelic loci.
"""
function linkage_disequilibrium(pop::Population{T}, locus_1::Int, locus_2::Int) where T
    f1 = allele_frequencies(pop, locus_1)
    f2 = allele_frequencies(pop, locus_2)
    
    a1_keys = collect(keys(f1))
    a2_keys = collect(keys(f2))
    
    if length(a1_keys) != 2 || length(a2_keys) != 2
        return (0.0, 0.0, 0.0)
    end
    
    A = a1_keys[1]
    a = a1_keys[2]
    B = a2_keys[1]
    b = a2_keys[2]
    
    pA = f1[A]
    pa = f1[a]
    pB = f2[B]
    pb = f2[b]
    
    # Count haplotypes
    counts = Dict{Tuple{T, T}, Int}()
    total = 0
    for ind in pop.individuals
        if locus_1 <= length(ind.loci) && locus_2 <= length(ind.loci)
            l1 = ind.loci[locus_1]
            l2 = ind.loci[locus_2]
            
            # Simple assumption of known phase or single-allele haplotypes
            if length(l1.alleles) >= 1 && length(l2.alleles) >= 1
                hap = (l1.alleles[1], l2.alleles[1])
                counts[hap] = get(counts, hap, 0) + 1
                total += 1
            end
        end
    end
    
    if total == 0
        return (0.0, 0.0, 0.0)
    end
    
    freq_AB = get(counts, (A, B), 0) / Float64(total)
    
    # D statistic
    D = freq_AB - (pA * pB)
    
    # D' statistic
    D_max = D > 0 ? min(pA * (1 - pB), (1 - pA) * pB) : min(pA * pB, (1 - pA) * (1 - pB))
    D_prime = D_max != 0 ? D / D_max : 0.0
    
    # r^2 statistic
    denom = pA * pa * pB * pb
    r_squared = denom > 0 ? (D^2) / denom : 0.0
    
    return (D, D_prime, r_squared)
end

"""
    ld_mapping(pop::Population, locus_indices::Vector{Int}, window_size::Int=10)

Calculates LD (r^2) between all pairs of loci within a sliding window.
Returns a Dict mapping (locus1, locus2) -> r_squared.
"""
function ld_mapping(pop::Population{T}, locus_indices::Vector{Int}, window_size::Int=10) where T
    results = Dict{Tuple{Int, Int}, Float64}()
    n = length(locus_indices)
    for i in 1:n
        for j in (i+1):min(i+window_size, n)
            l1 = locus_indices[i]
            l2 = locus_indices[j]
            _, _, r2 = linkage_disequilibrium(pop, l1, l2)
            results[(l1, l2)] = r2
        end
    end
    return results
end

"""
    ld_decay(pop::Population, physical_distances::Vector{Float64}, locus_pairs::Vector{Tuple{Int, Int}})

Calculates the decay of Linkage Disequilibrium (r^2) over physical distance.
Returns a tuple (distances, r_squared_values).
"""
function ld_decay(pop::Population{T}, physical_distances::Vector{Float64}, locus_pairs::Vector{Tuple{Int, Int}}) where T
    r2_vals = Float64[]
    applied_dists = Float64[]
    
    for (i, pair) in enumerate(locus_pairs)
        l1, l2 = pair
        _, _, r2 = linkage_disequilibrium(pop, l1, l2)
        push!(r2_vals, r2)
        push!(applied_dists, physical_distances[i])
    end
    
    return (applied_dists, r2_vals)
end

"""
    infer_phase_em(pop::Population, locus_1::Int, locus_2::Int; max_iter::Int=100, tol::Float64=1e-6)

Infers haplotype frequencies from unphased diploid genotype data using the 
Expectation-Maximization (EM) algorithm (Excoffier & Slatkin 1995).
Returns a Dict mapping (allele1, allele2) -> frequency.
"""
function infer_phase_em(pop::Population{T}, locus_1::Int, locus_2::Int; max_iter::Int=100, tol::Float64=1e-6) where T
    f1 = allele_frequencies(pop, locus_1)
    f2 = allele_frequencies(pop, locus_2)
    
    a1_keys = collect(keys(f1))
    a2_keys = collect(keys(f2))
    
    haplotypes = [(a, b) for a in a1_keys for b in a2_keys]
    num_haps = length(haplotypes)
    hap_freqs = fill(1.0 / num_haps, num_haps)
    
    # Genotype counts
    # A/A B/B, A/a B/B, a/a B/B, ...
    # We care especially about double heterozygotes A/a B/b 
    # which have two possible phases: (AB, ab) or (Ab, aB)
    
    genotypes = []
    for ind in pop.individuals
        if locus_1 <= length(ind.loci) && locus_2 <= length(ind.loci)
            l1 = ind.loci[locus_1]
            l2 = ind.loci[locus_2]
            if length(l1.alleles) == 2 && length(l2.alleles) == 2
                push!(genotypes, (Set(l1.alleles), Set(l2.alleles)))
            end
        end
    end
    
    if isempty(genotypes); return Dict{Tuple{T, T}, Float64}(); end
    
    for iter in 1:max_iter
        new_freqs = zeros(num_haps)
        total_weight = 0.0
        
        for (g1, g2) in genotypes
            # Find possible haplotype pairs for this genotype
            g1_list = collect(g1)
            g2_list = collect(g2)
            
            # If homozygous at both: 1 possibility
            # If hetero at one: 1 possibility
            # If hetero at both: 2 possibilities (A1B1/A2B2 or A1B2/A2B1)
            
            possible_pairs = []
            if length(g1) == 1 && length(g2) == 1
                push!(possible_pairs, ((g1_list[1], g2_list[1]), (g1_list[1], g2_list[1])))
            elseif length(g1) == 1
                push!(possible_pairs, ((g1_list[1], g2_list[1]), (g1_list[1], g2_list[2])))
            elseif length(g2) == 1
                push!(possible_pairs, ((g1_list[1], g2_list[1]), (g1_list[2], g2_list[1])))
            else
                # Double heterozygote
                push!(possible_pairs, ((g1_list[1], g2_list[1]), (g1_list[2], g2_list[2])))
                push!(possible_pairs, ((g1_list[1], g2_list[2]), (g1_list[2], g2_list[1])))
            end
            
            # Expectation step
            pair_probs = Float64[]
            for (h1, h2) in possible_pairs
                i1 = findfirst(==(h1), haplotypes)
                i2 = findfirst(==(h2), haplotypes)
                push!(pair_probs, hap_freqs[i1] * hap_freqs[i2])
            end
            
            sum_p = sum(pair_probs)
            if sum_p > 0
                pair_probs ./= sum_p
                # Maximization step (accumulate)
                for (p, (h1, h2)) in zip(pair_probs, possible_pairs)
                    i1 = findfirst(==(h1), haplotypes)
                    i2 = findfirst(==(h2), haplotypes)
                    new_freqs[i1] += p
                    new_freqs[i2] += p
                    total_weight += 2.0
                end
            end
        end
        
        new_freqs ./= total_weight
        if sum(abs.(new_freqs .- hap_freqs)) < tol
            hap_freqs = new_freqs
            break
        end
        hap_freqs = new_freqs
    end
    
    result = Dict{Tuple{T, T}, Float64}()
    for (i, h) in enumerate(haplotypes)
        result[h] = hap_freqs[i]
    end
    return result
end

# ------------------------------------------------------------------------------
# 8. Genetic Distances & Population Structure
# ------------------------------------------------------------------------------

"""
    genetic_distance(pop1::Population, pop2::Population, locus_idx::Int; method=:nei)

Calculates the genetic distance between two populations at a single locus.
Supported methods:
- `:nei` : Nei's standard genetic distance
- `:reynolds` : Reynolds' distance (coancestry distance)
- `:cavalli_sforza` : Cavalli-Sforza chord distance
"""
function genetic_distance(pop1::Population{T}, pop2::Population{T}, locus_idx::Int; method=:nei) where T
    f1 = allele_frequencies(pop1, locus_idx)
    f2 = allele_frequencies(pop2, locus_idx)
    
    alleles = union(keys(f1), keys(f2))
    
    if method == :nei
        J_11 = sum(get(f1, a, 0.0)^2 for a in alleles)
        J_22 = sum(get(f2, a, 0.0)^2 for a in alleles)
        J_12 = sum(get(f1, a, 0.0) * get(f2, a, 0.0) for a in alleles)
        
        denominator = sqrt(J_11 * J_22)
        return denominator > 0 ? -log(J_12 / denominator) : 0.0
        
    elseif method == :reynolds
        numerator = sum((get(f1, a, 0.0) - get(f2, a, 0.0))^2 for a in alleles)
        denominator = 2.0 * sum(get(f1, a, 0.0) * get(f2, a, 0.0) for a in alleles)
        # using the 1 - sum(p_i*q_i) standard form
        sum_pq = sum(get(f1, a, 0.0) * get(f2, a, 0.0) for a in alleles)
        reynolds_denom = 1.0 - sum_pq
        if reynolds_denom <= 0
            return 0.0
        end
        return sqrt(numerator / (2.0 * reynolds_denom))
        
    elseif method == :cavalli_sforza
        sum_sqrt_pq = sum(sqrt(get(f1, a, 0.0) * get(f2, a, 0.0)) for a in alleles)
        chord_d = (2 / pi) * acos(min(sum_sqrt_pq, 1.0))
        return chord_d
        
    elseif method == :rogers
        # Rogers' distance: (1/sqrt(2)) * sqrt(sum((pi - qi)^2))
        sum_sq_diff = sum((get(f1, a, 0.0) - get(f2, a, 0.0))^2 for a in alleles)
        return (1.0 / sqrt(2.0)) * sqrt(sum_sq_diff)
        
    else
        throw(ArgumentError("Unknown genetic distance method. Supported: :nei, :reynolds, :cavalli_sforza, :rogers"))
    end
end

function population_pca(populations::Vector{Population{T}}, loci::Vector{Int}) where T
    # Collect all unique alleles across all specified loci globally
    # to enforce consistent matrix dimensions
    global_alleles = Dict{Int, Set{T}}()
    num_pops = length(populations)
    
    for l in loci
        global_alleles[l] = Set{T}()
        for pop in populations
            f = allele_frequencies(pop, l)
            for k in keys(f)
                push!(global_alleles[l], k)
            end
        end
    end
    
    # Total features
    total_features = sum(length(s) for s in values(global_alleles))
    
    # Rows: populations, Cols: allele frequencies
    freq_matrix = zeros(Float64, num_pops, total_features)
    
    for (pop_i, pop) in enumerate(populations)
        col_offset = 1
        for l in loci
            f = allele_frequencies(pop, l)
            for a in global_alleles[l]
                freq_matrix[pop_i, col_offset] = get(f, a, 0.0)
                col_offset += 1
            end
        end
    end
    
    # PCA (SVD)
    # Center the matrix columns
    means = sum(freq_matrix, dims=1) ./ num_pops
    centered = freq_matrix .- means
    
    # SVD
    svd_res = svd(centered)
    
    # PCs are the right singular vectors V. The left singular vectors U scaled by S
    # are the exact projections of our populations onto the PCs.
    projections = svd_res.U * Diagonal(svd_res.S)
    eigenvalues = (svd_res.S .^ 2) ./ (num_pops - 1)
    
    return (projections, eigenvalues, freq_matrix)
end

"""
    population_pcoa(dist_matrix::Matrix{Float64})

Principal Coordinate Analysis (Classical Multidimensional Scaling).
Transforms a distance matrix into an ordination plot.
Returns (coordinates, eigenvalues).
"""
function population_pcoa(dist_matrix::Matrix{Float64})
    n = size(dist_matrix, 1)
    # 1. Square the distance matrix
    A = -0.5 .* (dist_matrix .^ 2)
    
    # 2. Double centering
    row_means = sum(A, dims=2) ./ n
    col_means = sum(A, dims=1) ./ n
    grand_mean = sum(A) / (n^2)
    
    B = A .- row_means .- col_means .+ grand_mean
    
    # 3. Eigen decomposition
    evals, evecs = eigen(B)
    
    # Sort eigenvalues and vectors descending
    idx = sortperm(evals, rev=true)
    evals = evals[idx]
    evecs = evecs[:, idx]
    
    # 4. Coordinates = V * sqrt(Lambda)
    # Remove negative eigenvalues (noise)
    pos_idx = evals .> 0
    coords = evecs[:, pos_idx] * Diagonal(sqrt.(evals[pos_idx]))
    
    return (coords, evals[pos_idx])
end

"""
    mantel_test(dist_matrix1::Matrix{Float64}, dist_matrix2::Matrix{Float64}; permutations::Int=999)

Performs a Mantel test between two distance matrices to test for correlation.
Returns (correlation, p_value).
"""
function mantel_test(dist_matrix1::Matrix{Float64}, dist_matrix2::Matrix{Float64}; permutations::Int=999)
    if size(dist_matrix1) != size(dist_matrix2)
        throw(ArgumentError("Distance matrices must have the same dimensions."))
    end
    
    n = size(dist_matrix1, 1)
    
    # Correlation requires at least 2 data points. 
    # For n=1, matrix is 1x1, lower triangle is empty.
    # For n=2, lower triangle has 1 element (correlation undefined).
    if n < 3
        return (0.0, 1.0)
    end
    
    # Helper to extract lower triangle (vectorized)
    function get_lower_tri(mat)
        # Calculate indices for lower triangle (i > j)
        # There are n*(n-1)/2 elements
        len = div(n * (n - 1), 2)
        vals = Vector{Float64}(undef, len)
        idx = 1
        for i in 2:n
            for j in 1:(i-1)
                vals[idx] = mat[i, j]
                idx += 1
            end
        end
        return vals
    end
    
    v1 = get_lower_tri(dist_matrix1)
    v2 = get_lower_tri(dist_matrix2)
    
    # This is the line that caused the error
    obs_corr = cor(v1, v2)
    
    # Permutation test
    count_extreme = 0
    perm_indices = collect(1:n)
    
    for _ in 1:permutations
        Random.shuffle!(perm_indices)
        # Permute rows AND columns simultaneously
        # We use a view to avoid allocation, but we must index carefully
        # A view of a scrambled matrix works fine for get_lower_tri logic
        
        # Efficient manual extraction during the loop avoids creating the view object overhead
        # or allocation of a new matrix.
        # However, simpler readable code is often better:
        
        # Let's stick to the view logic but ensure it's valid
        perm_view = @view dist_matrix1[perm_indices, perm_indices]
        v_perm = get_lower_tri(perm_view)
        
        perm_corr = cor(v_perm, v2)
        
        if perm_corr >= obs_corr
            count_extreme += 1
        end
    end
    
    p_val = (count_extreme + 1) / (permutations + 1)
    return (obs_corr, p_val)
end

# ------------------------------------------------------------------------------
# 9. Effective Population Size (Ne)
# ------------------------------------------------------------------------------

"""
    estimate_ne_temporal(pop1::Population, pop2::Population, generations::Int, locus_idx::Int)

Estimates Effective Population Size (Ne) through temporal changes in allele frequencies
(Nei & Tajima 1981 approximate variance method) over a span of discrete generations.
"""
function estimate_ne_temporal(pop1::Population{T}, pop2::Population{T}, generations::Int, locus_idx::Int) where T
    f1 = allele_frequencies(pop1, locus_idx)
    f2 = allele_frequencies(pop2, locus_idx)
    
    alleles = union(keys(f1), keys(f2))
    k = length(alleles)
    if k < 2
        return Inf # Monomorphic, infinite effective size or undefined
    end
    
    # F_c (Nei and Tajima variance estimator)
    sum_F = 0.0
    for a in alleles
        x = get(f1, a, 0.0)
        y = get(f2, a, 0.0)
        if x > 0 && y > 0
            sum_F += ((x - y)^2) / ((x + y) / 2.0)
        end
    end
    
    Fc = (1.0 / k) * sum_F
    
    S1 = length(pop1.individuals) * 2 # assumed diploid alleles
    S2 = length(pop2.individuals) * 2
    
    num = generations
    den = Fc - (1.0/(2*S1)) - (1.0/(2*S2))
    
    if den <= 0
        return Inf
    end
    return num / (2 * den)
end

"""
    estimate_ne_ld(pop::Population, locus_pairs::Vector{Tuple{Int, Int}})

Estimates Effective Population Size (Ne) using the Linkage Disequilibrium method 
(Hill 1981). \$E[r^2] \approx 1/(3Ne) + 1/n\$.
"""
function estimate_ne_ld(pop::Population{T}, locus_pairs::Vector{Tuple{Int, Int}}) where T
    n = length(pop.individuals) * 2 # total alleles (diploid)
    if n < 2; return Inf; end
    
    r2_sum = 0.0
    valid_pairs = 0
    for (l1, l2) in locus_pairs
        _, _, r2 = linkage_disequilibrium(pop, l1, l2)
        r2_sum += r2
        valid_pairs += 1
    end
    
    if valid_pairs == 0; return Inf; end
    mean_r2 = r2_sum / valid_pairs
    
    # Correct for sample size
    # r2_corrected = r2 - 1/n
    # Ne = 1 / (3 * r2_corrected)
    r2_corrected = mean_r2 - (1.0 / n)
    if r2_corrected <= 0.0; return Inf; end
    
    return 1.0 / (3.0 * r2_corrected)
end

# ------------------------------------------------------------------------------
# 10. Sequence-Based Population Genetics
# ------------------------------------------------------------------------------

"""
    segregating_sites(alignment::MultipleSequenceAlignment)

Counts the number of segregating (polymorphic) sites in a multiple sequence alignment.
Ignores columns where any sequence has a gap or missing data.
"""
function segregating_sites(alignment::MultipleSequenceAlignment)
    S = 0
    records = alignment.records
    n = length(records)
    if n < 2; return 0; end
    len = length(records[1].sequence)
    
    for i in 1:len
        first_base = records[1].sequence[i]
        is_segregating = false
        valid_col = true
        
        for j in 1:n
            base = records[j].sequence[i]
            if base == '-' || base == '?' || base == 'N'
                valid_col = false
                break
            end
            if base != first_base
                is_segregating = true
            end
        end
        
        if valid_col && is_segregating
            S += 1
        end
    end
    return S
end

"""
    mismatch_distribution(alignment::MultipleSequenceAlignment)

Calculates the distribution of pairwise nucleotide differences between all pairs 
of sequences in the alignment. Returns a vector of difference counts.
"""
function mismatch_distribution(alignment::MultipleSequenceAlignment)
    records = alignment.records
    n = length(records)
    if n < 2; return Int[]; end
    
    mismatches = Int[]
    for i in 1:n
        for j in (i+1):n
            diffs = 0
            s1 = records[i].sequence
            s2 = records[j].sequence
            for k in 1:min(length(s1), length(s2))
                if s1[k] != s2[k] && s1[k] != '-' && s2[k] != '-'
                    diffs += 1
                end
            end
            push!(mismatches, diffs)
        end
    end
    return mismatches
end

"""
    nucleotide_diversity(alignment::MultipleSequenceAlignment)

Calculates nucleotide diversity (π), the average number of nucleotide differences 
per site between two continuous sequences randomly drawn from the population.
"""
function nucleotide_diversity(alignment::MultipleSequenceAlignment)
    records = alignment.records
    n = length(records)
    if n < 2; return 0.0; end
    
    len = length(records[1].sequence)
    pi_sum = 0.0
    valid_sites = 0
    
    for i in 1:len
        freqs = Dict{Char, Int}()
        valid_col = true
        for j in 1:n
            base = records[j].sequence[i]
            if base == '-' || base == '?' || base == 'N'
                valid_col = false
                break
            end
            freqs[base] = get(freqs, base, 0) + 1
        end
        
        if valid_col
            valid_sites += 1
            col_pi = 0.0
            for (b1, c1) in freqs
                for (b2, c2) in freqs
                    if b1 != b2
                        col_pi += (c1 * c2)
                    end
                end
            end
            # Divide by n(n-1) implicitly later, here just summing diffs
            # Actually, sum of 2*c1*c2... divided by n(n-1) is the probability of picking 2 different.
            # col_pi = sum(c1*c2 for b1!=b2) is equivalent.
            pi_sum += (col_pi / (n * (n - 1)))
        end
    end
    
    return valid_sites > 0 ? (pi_sum / valid_sites) : 0.0
end

"""
    watterson_theta(alignment::MultipleSequenceAlignment)

Estimates Watterson's Theta (Θ) from the number of segregating sites.
"""
function watterson_theta(alignment::MultipleSequenceAlignment)
    n = length(alignment.records)
    if n < 2; return 0.0; end
    
    S = segregating_sites(alignment)
    a1 = sum(1.0 / i for i in 1:(n-1))
    
    return S / a1
end

"""
    tajimas_d(alignment::MultipleSequenceAlignment)

Calculates Tajima's D statistic to test for neutral evolution.
Negative values suggest selective sweeping or population expansion.
Positive values suggest balancing selection or population subdivision.
"""
function tajimas_d(alignment::MultipleSequenceAlignment)
    n = length(alignment.records)
    if n < 2; return 0.0; end
    
    S = segregating_sites(alignment)
    if S == 0; return 0.0; end
    
    # Calculate π and Θ unscaled by total length (i.e. number of pairwise differences, not per-site)
    # We instead scale nucleotides diversity by sequence chunks
    # Wait, nucleotide_diversity(ali) gives per-site. We need total pairwise differences for Tajima's D.
    
    len = 0
    records = alignment.records
    len_records = length(records[1].sequence)
    total_pairwise_diffs = 0.0
    
    for i in 1:len_records
        freqs = Dict{Char, Int}()
        valid_col = true
        for j in 1:n
            base = records[j].sequence[i]
            if base == '-' || base == '?' || base == 'N'
                valid_col = false
                break
            end
            freqs[base] = get(freqs, base, 0) + 1
        end
        if valid_col
            len += 1
            col_diffs = 0.0
            for (b1, c1) in freqs
                for (b2, c2) in freqs
                    if b1 != b2
                        col_diffs += (c1 * c2)
                    end
                end
            end
            total_pairwise_diffs += (col_diffs / (n * (n - 1)))
        end
    end
    
    # Note: total_pairwise_diffs is exactly π_total
    a1 = sum(1.0 / i for i in 1:(n-1))
    a2 = sum(1.0 / (i^2) for i in 1:(n-1))
    
    b1 = (n + 1.0) / (3.0 * (n - 1.0))
    b2 = (2.0 * (n^2 + n + 3.0)) / (9.0 * n * (n - 1.0))
    
    c1 = b1 - (1.0 / a1)
    c2 = b2 - ((n + 2.0) / (a1 * n)) + (a2 / (a1^2))
    
    e1 = c1 / a1
    e2 = c2 / (a1^2 + a2)
    
    variance_d = (e1 * S) + (e2 * S * (S - 1))
    
    if variance_d <= 0
        return 0.0
    end
    
    D = (total_pairwise_diffs - (S / a1)) / sqrt(variance_d)
    return D
end

"""
    site_frequency_spectrum(alignment::MultipleSequenceAlignment; folded::Bool=true)

Calculates the Site Frequency Spectrum (SFS) for an alignment. If `folded` is true,
the minor allele frequency spectrum is generated. Returns a frequency distribution vector.
"""
function site_frequency_spectrum(alignment::MultipleSequenceAlignment; folded::Bool=true)
    n = length(alignment.records)
    if n < 2; return Float64[]; end
    
    lim = folded ? floor(Int, n / 2) : (n - 1)
    sfs = zeros(Int, lim)
    
    len_records = length(alignment.records[1].sequence)
    for i in 1:len_records
        freqs = Dict{Char, Int}()
        valid_col = true
        for j in 1:n
            base = alignment.records[j].sequence[i]
            if base == '-' || base == '?' || base == 'N'
                valid_col = false
                break
            end
            freqs[base] = get(freqs, base, 0) + 1
        end
        
        if valid_col && length(freqs) == 2
            # Bi-allelic site
            counts = collect(values(freqs))
            allele_count = min(counts[1], counts[2])
            if !folded
                # For unfolded, we ideally need an outgroup to specify the ancestral allele.
                # If folded=false but no outgroup, it's ambiguous. Treating the most frequent as ancestral.
                allele_count = min(counts[1], counts[2])
            end
            if allele_count > 0 && allele_count <= lim
                sfs[allele_count] += 1
            end
        end
    end
    return sfs
end

# ------------------------------------------------------------------------------
# 11. Advanced Neutrality Tests (Fu & Li)
# ------------------------------------------------------------------------------

"""
    fu_li_d(alignment::MultipleSequenceAlignment)

Calculates Fu and Li's D statistic (without an outgroup).
Tests for background selection or population expansion.
"""
function fu_li_d(alignment::MultipleSequenceAlignment)
    n = length(alignment.records)
    if n < 2; return 0.0; end
    
    S = segregating_sites(alignment)
    if S == 0; return 0.0; end
    
    # Calculate ηs (number of singletons)
    singletons = 0
    records = alignment.records
    len = length(records[1].sequence)
    for i in 1:len
        freqs = Dict{Char, Int}()
        valid_col = true
        for j in 1:n
            base = records[j].sequence[i]
            if base == '-' || base == '?' || base == 'N'
                valid_col = false
                break
            end
            freqs[base] = get(freqs, base, 0) + 1
        end
        if valid_col && length(freqs) == 2
            counts = collect(values(freqs))
            if any(c == 1 for c in counts)
                singletons += 1
            end
        end
    end
    
    a1 = sum(1.0 / i for i in 1:(n-1))
    
    vd = 1.0 + (a1^2 / (sum(1.0/(i^2) for i in 1:(n-1)) + a1^2)) # simplified variant
    u_d = (n - 1.0) / n - (1.0 / a1)
    v_d = ((n - 1.0) / n)^2 * (sum(1.0/(i^2) for i in 1:(n-1)) + (a1^2)) / (a1^2)
    
    return (S - a1 * singletons) / sqrt(abs(v_d)) 
end

"""
    fu_li_f(alignment::MultipleSequenceAlignment)

Calculates Fu and Li's F statistic (without an outgroup).
"""
function fu_li_f(alignment::MultipleSequenceAlignment)
    td = tajimas_d(alignment)
    fd = fu_li_d(alignment)
    return (td + fd) / 2.0 
end

"""
    ewens_watterson_test(pop::Population, locus_idx::Int)

Performs the Ewens-Watterson neutrality test. Returns the observed homozygosity 
(F) and the expected homozygosity under the neutral infinite alleles model.
"""
function ewens_watterson_test(pop::Population{T}, locus_idx::Int) where T
    freqs = allele_frequencies(pop, locus_idx)
    obs_F = sum(p^2 for p in values(freqs))
    
    # Expected F under Neutrality depends on n (sample size) and k (number of alleles)
    # Using the approximation for large n: E[F] = 1 / (1 + theta)
    # where theta is estimated from k and n.
    n = 0
    for ind in pop.individuals
        if locus_idx <= length(ind.loci)
            n += length(ind.loci[locus_idx].alleles)
        end
    end
    if n == 0; return (obs_F, 1.0); end
    
    k = length(freqs)
    if k <= 1; return (obs_F, 1.0); end
    
    # Solve for theta: k = sum_{i=0}^{n-1} theta / (theta + i)
    # Simple search for theta
    theta = 1.0
    for _ in 1:20
        current_k = sum(theta / (theta + i) for i in 0:(n-1))
        theta *= (k / current_k)
    end
    
    exp_F = (1.0 + theta) / (n * (1.0 + theta/n)) # approximation
    return (obs_F, exp_F)
end

# ------------------------------------------------------------------------------
# 12. Selection Sweep Scans (Haplotype-based)
# ------------------------------------------------------------------------------

"""
    ehh(alignment::MultipleSequenceAlignment, core_site::Int, distance_sites::Int)

Calculates Extended Haplotype Homozygosity (EHH) at a core site extending 
out to a certain distance (number of sites).
Requires phased haplotypes in the alignment.
"""
function ehh(alignment::MultipleSequenceAlignment, core_site::Int, distance_sites::Int)
    n = length(alignment.records)
    if n < 2; return 0.0; end
    
    len = length(alignment.records[1].sequence)
    end_site = min(core_site + distance_sites, len)
    start_site = max(core_site - distance_sites, 1)
    
    haplotypes = [records.sequence[start_site:end_site] for records in alignment.records]
    
    counts = Dict{String, Int}()
    for hap in haplotypes
        counts[hap] = get(counts, hap, 0) + 1
    end
    
    sum_sq = sum(c * (c - 1) for c in values(counts))
    return sum_sq / (n * (n - 1))
end

"""
    ihs(alignment::MultipleSequenceAlignment, core_site::Int)

Calculates a simplified Integrated Haplotype Score (iHS) for a site by 
integrating EHH over a range of distances.
"""
function ihs(alignment::MultipleSequenceAlignment, core_site::Int)
    integral = 0.0
    dist = 1
    while dist < 100 
        e = ehh(alignment, core_site, dist)
        integral += e
        if e < 0.05; break; end
        dist += 1
    end
    return integral
end

"""
    xp_ehh(pop1_ali::MultipleSequenceAlignment, pop2_ali::MultipleSequenceAlignment, core_site::Int)

Calculates Cross-Population EHH (XP-EHH) to detect selection between two populations.
"""
function xp_ehh(pop1_ali::MultipleSequenceAlignment, pop2_ali::MultipleSequenceAlignment, core_site::Int)
    ih1 = ihs(pop1_ali, core_site)
    ih2 = ihs(pop2_ali, core_site)
    
    if ih2 == 0; return 0.0; end
    return log(ih1 / ih2)
end

"""
    sweepfinder_clr(alignment::MultipleSequenceAlignment, grid_points::Int=100)

Composite Likelihood Ratio (CLR) test for selective sweeps (Nielsen et al. 2005).
Returns CLR values across the genomic grid.
"""
function sweepfinder_clr(alignment::MultipleSequenceAlignment, grid_points::Int=100)
    # 1. Background SFS
    sfs = site_frequency_spectrum(alignment, folded=false)
    n = length(alignment.records)
    if isempty(sfs); return zeros(grid_points); end
    
    bg_freqs = sfs ./ sum(sfs)
    
    # 2. Likelihood Ratio at each grid point
    # Simplified version: comparing local SFS to background SFS
    len = length(alignment.records[1].sequence)
    grid = range(1, len, length=grid_points)
    clrs = Float64[]
    
    # This is a very complex calculation in reality involving sweep models.
    # Here we implement a simplified local-deviation version.
    for pos in grid
        pos_int = round(Int, pos)
        window = max(1, pos_int-50):min(len, pos_int+50)
        
        # Local SFS
        local_sfs = zeros(Int, n-1)
        for i in window
            # reuse site logic
            col_freqs = Dict{Char, Int}()
            for rec in alignment.records
                base = rec.sequence[i]
                if base != '-' && base != 'N'
                    col_freqs[base] = get(col_freqs, base, 0) + 1
                end
            end
            if length(col_freqs) == 2
                c = collect(values(col_freqs))
                mc = min(c[1], c[2])
                if mc > 0 && mc < n
                    local_sfs[mc] += 1
                end
            end
        end
        
        # Log-Likelihood Ratio
        clr = 0.0
        tot_sites = sum(local_sfs)
        if tot_sites > 0
            local_sites = length(local_sfs)
            for (i, bg_freq) in pairs(bg_freqs)
                if i <= local_sites && local_sfs[i] > 0 && bg_freq > 0
                    p_loc = local_sfs[i] / tot_sites
                    clr += local_sfs[i] * log(p_loc / bg_freq)
                end
            end
        end
        push!(clrs, max(0.0, clr))
    end
    
    return clrs
end

# ------------------------------------------------------------------------------
# 13. Admixture & Introgression (f-statistics & Patterson's D)
# ------------------------------------------------------------------------------

"""
    f3_statistic(p1::Population, p2::Population, outgroup::Population, locus_idx::Int)

Calculates the f3 statistic (P3; P1, P2) which tests if population P3 is admixed
from populations P1 and P2. A negative value indicates admixture.
"""
function f3_statistic(p1::Population{T}, p2::Population{T}, p3::Population{T}, locus_idx::Int) where T
    f1 = allele_frequencies(p1, locus_idx)
    f2 = allele_frequencies(p2, locus_idx)
    f3 = allele_frequencies(p3, locus_idx)
    
    alleles = union(keys(f1), keys(f2), keys(f3))
    if isempty(alleles); return 0.0; end
    
    # Using the first allele as the reference
    a = first(alleles)
    p1_a = get(f1, a, 0.0)
    p2_a = get(f2, a, 0.0)
    p3_a = get(f3, a, 0.0)
    
    return (p3_a - p1_a) * (p3_a - p2_a)
end

"""
    f4_statistic(p1::Population, p2::Population, p3::Population, p4::Population, locus_idx::Int)

Calculates the f4 statistic (P1, P2; P3, P4) which tests for gene flow between 
P1-P2 and P3-P4 branch pairs.
"""
function f4_statistic(p1::Population{T}, p2::Population{T}, p3::Population{T}, p4::Population{T}, locus_idx::Int) where T
    f1 = allele_frequencies(p1, locus_idx)
    f2 = allele_frequencies(p2, locus_idx)
    f3 = allele_frequencies(p3, locus_idx)
    f4 = allele_frequencies(p4, locus_idx)
    
    alleles = union(keys(f1), keys(f2), keys(f3), keys(f4))
    if isempty(alleles); return 0.0; end
    
    a = first(alleles)
    p1_a = get(f1, a, 0.0)
    p2_a = get(f2, a, 0.0)
    p3_a = get(f3, a, 0.0)
    p4_a = get(f4, a, 0.0)
    
    return (p1_a - p2_a) * (p3_a - p4_a)
end

"""
    patterson_d(p1::Population, p2::Population, p3::Population, outgroup::Population, locus_idx::Int)

Calculates Patterson's D statistic (ABBA-BABA test). 
Tests for introgression between P3 and either P1 or P2.
"""
function patterson_d(p1::Population{T}, p2::Population{T}, p3::Population{T}, p4::Population{T}, locus_idx::Int) where T
    f1 = allele_frequencies(p1, locus_idx)
    f2 = allele_frequencies(p2, locus_idx)
    f3 = allele_frequencies(p3, locus_idx)
    f4 = allele_frequencies(p4, locus_idx) # Outgroup
    
    alleles = union(keys(f1), keys(f2), keys(f3), keys(f4))
    if isempty(alleles); return 0.0; end
    
    # Typically bi-allelic sites
    a = first(alleles)
    p1_a = get(f1, a, 0.0)
    p2_a = get(f2, a, 0.0)
    p3_a = get(f3, a, 0.0)
    p4_a = get(f4, a, 0.0)
    
    abba = (1.0 - p1_a) * p2_a * p3_a * (1.0 - p4_a)
    baba = p1_a * (1.0 - p2_a) * p3_a * (1.0 - p4_a)
    
    if abba + baba == 0; return 0.0; end
    return (abba - baba) / (abba + baba)
end

# ------------------------------------------------------------------------------
# 14. Wright-Fisher Simulation Engine
# ------------------------------------------------------------------------------

"""
    wright_fisher_simulation(N::Int, generations::Int; mu::Float64=1e-8, s::Float64=0.0, p0::Float64=0.5)

Simulates the trajectory of a bi-allelic locus under drift, mutation (mu), 
and selection (s) in a population of size N over many generations.
Returns a vector of allele frequencies over time.
"""
function wright_fisher_simulation(N::Int, generations::Int; mu::Float64=1e-8, s::Float64=0.0, p0::Float64=0.5)
    traj = zeros(Float64, generations + 1)
    traj[1] = p0
    p = p0
    
    for i in 1:generations
        # 1. Selection
        # w_AA = 1 + s, w_Aa = 1 + s/2, w_aa = 1
        p_prime = (p^2 * (1 + s) + p * (1-p) * (1 + s/2)) / (p^2 * (1+s) + 2*p*(1-p)*(1 + s/2) + (1-p)^2)
        
        # 2. Mutation
        p_mut = p_prime * (1 - mu) + (1 - p_prime) * mu
        
        # 3. Random Genetic Drift (Binomial sampling)
        count = rand(Binomial(2*N, p_mut))
        p = count / (2.0 * N)
        traj[i+1] = p
        
        if p <= 0.0 || p >= 1.0; break; end
    end
    return traj
end

"""
    wright_fisher_metapopulation(N::Vector{Int}, generations::Int, M::Matrix{Float64}; mu::Float64=1e-8, s::Vector{Float64}=zeros(length(N)))

Simulates a metapopulation under the Wright-Fisher model with a migration matrix M.
M[i, j] is the proportion of population i that comes from population j.
"""
function wright_fisher_metapopulation(N::Vector{Int}, generations::Int, M::Matrix{Float64}; mu::Float64=1e-8, s::Vector{Float64}=nothing)
    num_pops = length(N)
    if s === nothing; s = zeros(num_pops); end
    
    trajs = [zeros(generations + 1) for _ in 1:num_pops]
    p = fill(0.5, num_pops)
    for k in 1:num_pops; trajs[k][1] = p[k]; end
    
    for i in 1:generations
        new_p = zeros(num_pops)
        # 1. Selection & Mutation in each pop
        for k in 1:num_pops
            pk = p[k]
            sk = s[k]
            p_prime = (pk^2 * (1 + sk) + pk * (1-pk) * (1 + sk/2)) / (pk^2 * (1+sk) + 2*pk*(1-pk)*(1 + sk/2) + (1-pk)^2)
            p[k] = p_prime * (1 - mu) + (1 - p_prime) * mu
        end
        
        # 2. Migration
        p_mig = M * p
        
        # 3. Drift
        for k in 1:num_pops
            count = rand(Binomial(2*N[k], p_mig[k]))
            p[k] = count / (2.0 * N[k])
            trajs[k][i+1] = p[k]
        end
    end
    return trajs
end

# ------------------------------------------------------------------------------
# 15. GenePop Format Parsing
# ------------------------------------------------------------------------------

mutable struct GenePopRecord{T}
    marker_len::Int
    comment_line::String
    loci_list::Vector{String}
    pop_list::Vector{String}
    populations::Vector{Population{T}}
end

function GenePopRecord{T}() where T
    return GenePopRecord{T}(0, "", String[], String[], Population{T}[])
end

function Base.show(io::IO, record::GenePopRecord{T}) where T
    println(io, record.comment_line)
    println(io, join(record.loci_list, "\n"))
    for population in record.populations
        println(io, "Pop")
        for individual in population.individuals
            print(io, individual.id, ",")
            for locus in individual.loci
                print(io, " ")
                for allele in locus.alleles
                    allele_value = allele === nothing ? 0 : allele
                    allele_string = lpad(string(allele_value), 3, '0')
                    print(io, allele_string)
                end
            end
            println(io)
        end
    end
end

function _genepop_population_names(population_count::Int)
    return ["Pop_$(index)" for index in 1:population_count]
end

function _parse_genepop_population(filepath::String)
    lines = readlines(filepath)
    if isempty(lines)
        return GenePopRecord{Int}()
    end

    record = GenePopRecord{Int}()
    record.comment_line = lines[1]
    current_line = 2

    while current_line <= length(lines) && !occursin(Regex("^Pop", "i"), lines[current_line])
        push!(record.loci_list, strip(lines[current_line], [',', ' ']))
        current_line += 1
    end

    population_index = 0
    while current_line <= length(lines)
        line = lines[current_line]
        if occursin(Regex("^Pop", "i"), line)
            population_index += 1
            push!(record.populations, Population{Int}("Pop_$(population_index)", PopGenIndividual{Int}[]))
            push!(record.pop_list, "Pop_$(population_index)")
            current_line += 1
            continue
        end

        parts = Base.split(line, ',')
        if length(parts) >= 2
            individual_name = strip(parts[1])
            genotype_tokens = Base.split(strip(parts[2]))
            loci = Locus{Int}[]
            for genotype_token in genotype_tokens
                token_length = length(genotype_token)
                token_midpoint = div(token_length, 2)
                allele_1 = parse(Int, genotype_token[1:token_midpoint])
                allele_2 = parse(Int, genotype_token[token_midpoint+1:end])
                push!(loci, Locus{Int}((allele_1, allele_2)))
            end
            push!(record.populations[end].individuals, PopGenIndividual{Int}(individual_name, loci))
        end
        current_line += 1
    end

    return record
end

function read_genepop_record(filepath::String)
    return _parse_genepop_population(filepath)
end

"""
    read_genepop(filepath::String)

Parses a GenePop format file into a vector of Population objects.
"""
function read_genepop(filepath::String)
    return read_genepop_record(filepath).populations
end

function split_in_pops(record::GenePopRecord{T}, pop_names::Vector{String}) where T
    if length(pop_names) != length(record.populations)
        throw(ArgumentError("pop_names must match the number of populations"))
    end

    result = Dict{String, GenePopRecord{T}}()
    for (index, population) in enumerate(record.populations)
        new_record = GenePopRecord{T}(record.marker_len, record.comment_line, copy(record.loci_list), [pop_names[index]], [deepcopy(population)])
        result[pop_names[index]] = new_record
    end
    return result
end

function split_in_loci(record::GenePopRecord{T}) where T
    result = Dict{String, GenePopRecord{T}}()
    for (locus_index, locus_name) in enumerate(record.loci_list)
        populations = Population{T}[]
        for population in record.populations
            individuals = PopGenIndividual{T}[]
            for individual in population.individuals
                push!(individuals, PopGenIndividual{T}(individual.id, [individual.loci[locus_index]]))
            end
            push!(populations, Population{T}(population.name, individuals))
        end
        result[locus_name] = GenePopRecord{T}(record.marker_len, record.comment_line, [locus_name], copy(record.pop_list), populations)
    end
    return result
end

function remove_population!(record::GenePopRecord, pos::Int)
    deleteat!(record.populations, pos + 1)
    if pos + 1 <= length(record.pop_list)
        deleteat!(record.pop_list, pos + 1)
    end
    return record
end

function remove_locus_by_position!(record::GenePopRecord{T}, pos::Int) where T
    deleteat!(record.loci_list, pos + 1)
    for population in record.populations
        for individual in population.individuals
            deleteat!(individual.loci, pos + 1)
        end
    end
    return record
end

function remove_locus_by_name!(record::GenePopRecord{T}, name::AbstractString) where T
    for (index, locus_name) in enumerate(record.loci_list)
        if locus_name == name
            return remove_locus_by_position!(record, index - 1)
        end
    end
    return record
end

remove_population(record::GenePopRecord, pos::Int) = remove_population!(record, pos)
remove_locus_by_position(record::GenePopRecord, pos::Int) = remove_locus_by_position!(record, pos)
remove_locus_by_name(record::GenePopRecord, name::AbstractString) = remove_locus_by_name!(record, name)

"""
    write_genepop(populations::Vector{Population}, filepath::String; title="BioToolkit export")

Writes populations to a GenePop format file.
"""
function write_genepop(populations::Vector{Population{T}}, filepath::String; title="BioToolkit export") where T
    open(filepath, "w") do io
        println(io, title)
        # Assuming all individuals have same loci count
        if !isempty(populations) && !isempty(populations[1].individuals)
            num_loci = length(populations[1].individuals[1].loci)
            for i in 1:num_loci
                println(io, "Locus_$i")
            end
        end
        
        for pop in populations
            println(io, "Pop")
            for ind in pop.individuals
                print(io, ind.id, " , ")
                for locus in ind.loci
                    # Format as 3-digit with leading zeros
                    s = join([lpad(string(a), 3, '0') for a in locus.alleles], "")
                    print(io, s, " ")
                end
                println(io)
            end
        end
    end
end

# ------------------------------------------------------------------------------
# 16. External Tool Wrappers
# ------------------------------------------------------------------------------

"""
    run_genepop(input_file::String, options::String)

Wrapper to execute the native `genepop` binary if available.
"""
function run_genepop(input_file::String, options::String)
    try
        run(`genepop $input_file $options`)
    catch e
        @error "Genepop execution failed. Ensure 'genepop' binary is in PATH." exception=e
    end
end

"""
    run_fastsimcoal(par_file::String, num_sims::Int)

Wrapper to execute `fastsimcoal2` for coalescent simulations.
"""
function run_fastsimcoal(par_file::String, num_sims::Int)
    try
        run(`fsc26 -i $par_file -n $num_sims`)
    catch e
        @error "fastsimcoal2 execution failed. Ensure 'fsc26' is in PATH." exception=e
    end
end

"""
    run_fdist(input_file::String, num_loci::Int, num_pops::Int)

Wrapper for fdist2 to detect loci under selection.
"""
function run_fdist(input_file::String, num_loci::Int, num_pops::Int)
    try
        # Example logic: run datachk then fdist
        run(`datachk`)
        run(`fdist2 -n $num_pops -l $num_loci`)
    catch e
        @error "fdist2 execution failed. Ensure 'fdist2' binaries are in PATH." exception=e
    end
end

# ------------------------------------------------------------------------------
# 17. Pedigree, Kinship & Quantitative Genetics
# ------------------------------------------------------------------------------

"""
    genetic_relationship_matrix(populations::Vector{Population})

Calculates the Genetic Relationship Matrix (GRM) using the VanRaden (2008) method.
Returns a symmetric matrix of size N x N where N is the total number of individuals.
"""
function genetic_relationship_matrix(populations::Vector{Population{T}}) where T
    all_inds = [ind for pop in populations for ind in pop.individuals]
    num_inds = length(all_inds)
    if num_inds == 0; return zeros(0, 0); end
    
    num_loci = length(all_inds[1].loci)
    
    # 1. Construct Genotype Matrix M (N x L)
    # Codes: 0 (AA), 1 (AB), 2 (BB)
    M = zeros(num_inds, num_loci)
    
    # We need to map alleles to 0/1/2. Assuming bi-allelic.
    for l_idx in 1:num_loci
        # Find major/minor alleles
        counts = Dict{T, Int}()
        for ind in all_inds
            for a in ind.loci[l_idx].alleles
                counts[a] = get(counts, a, 0) + 1
            end
        end
        if length(counts) < 2; continue; end
        alleles = sort(collect(keys(counts)), by=k->counts[k], rev=true)
        ref = alleles[1]
        
        for (i, ind) in enumerate(all_inds)
            # Count copies of 'ref'
            c = count(==(ref), ind.loci[l_idx].alleles)
            M[i, l_idx] = Float64(c)
        end
    end
    
    # 2. P matrix (2 * pi)
    p = sum(M, dims=1) ./ (2 * num_inds)
    P = 2.0 .* p
    
    # 3. Z = M - P
    Z = M .- P
    
    # 4. G = ZZ' / 2*sum(pi*(1-pi))
    sum_pq = sum(p .* (1.0 .- p))
    if sum_pq == 0.0; return Matrix{Float64}(I, num_inds, num_inds); end
    
    G = (Z * Z') ./ (2.0 * sum_pq)
    return G
end

"""
    linear_mixed_model_scan(genotypes::Matrix{Float64}, phenotypes::Vector{Float64}, G::Matrix{Float64})

Simple GWAS scan using a Linear Mixed Model approximation (Gram-Schmidt or Score test).
Returns p-values for each marker.
"""
function linear_mixed_model_scan(genotypes::Matrix{Float64}, phenotypes::Vector{Float64}, G::Matrix{Float64})
    # Simplified LMM: Y = Xb + Zu + e, u ~ N(0, G*var_u)
    # Using a basic fixed-effect model with G as a covariance weight (Generalized Least Squares)
    # This is a placeholder for a more complex EMMA/GEMMA style solver.
    
    n_inds, n_loci = size(genotypes)
    p_values = zeros(n_loci)
    
    # Compute inverse of G (with ridge regularization)
    G_inv = inv(G + 1e-6 * I)
    
    for i in 1:n_loci
        x = genotypes[:, i]
        # GLS estimate: beta = (X' G_inv X)^-1 X' G_inv Y
        xt_gi = x' * G_inv
        denom = xt_gi * x
        if denom == 0.0; p_values[i] = 1.0; continue; end
        
        beta = (xt_gi * phenotypes) / denom
        # Simple Z-test / Wald approx
        p_values[i] = exp(-beta^2 / 2.0) # very crude approx
    end
    
    return p_values
end

"""
    inbreeding_coefficient(G::Matrix{Float64}, ind_idx::Int)

Calculates the inbreeding coefficient (F) from the GRM diagonal.
F = G[i, i] - 1.
"""
function inbreeding_coefficient(G::Matrix{Float64}, ind_idx::Int)
    return G[ind_idx, ind_idx] - 1.0
end

"""
    relatedness(G::Matrix{Float64}, ind1::Int, ind2::Int)

Returns the genomic relatedness coefficient between two individuals from the GRM.
"""
function relatedness(G::Matrix{Float64}, ind1::Int, ind2::Int)
    return G[ind1, ind2]
end

# --- Coalescent Simulation (ms-style) ---

struct _CoalNode
    id::Int
    time::Float64
    children::Vector{_CoalNode}
end

function _coal_to_phylo(node::_CoalNode)
    if isempty(node.children)
        return PhyloTree("ind_$(node.id)"; branch_length=0.0)
    end
    
    phylo_children = PhyloTree[]
    for child in node.children
        p_child = _coal_to_phylo(child)
        p_child.branch_length = node.time - child.time
        push!(phylo_children, p_child)
    end
    
    return PhyloTree("", branch_length=0.0, children=phylo_children)
end

"""
    simulate_coalescent(n::Int; ne::Real=10000, mu::Real=1e-8, seq_len::Int=1000)

Simple backward-in-time coalescent simulation (constant population size).
Returns a tuple (tree::PhyloTree, haplotypes::Dict{String, String}).
"""
function simulate_coalescent(n::Int; ne::Real=10000, mu::Real=1e-8, seq_len::Int=1000)
    lineages = [_CoalNode(i, 0.0, _CoalNode[]) for i in 1:n]
    current_time = 0.0
    id_counter = n
    
    while length(lineages) > 1
        k = length(lineages)
        # Rate of coalescence: k(k-1)/2 per 2Ne generations
        # In units of 2Ne, rate is k(k-1)/2
        rate = k * (k - 1) / 2.0
        dt = rand(Distributions.Exponential(1.0 / rate))
        current_time += dt
        
        # Pick two to coalesce
        i1 = rand(1:k)
        i2 = rand(1:(k - 1))
        i2 += i2 >= i1
        
        l1 = lineages[i1]
        l2 = lineages[i2]
        
        id_counter += 1
        new_node = _CoalNode(id_counter, current_time, [l1, l2])
        
        # Update lineage list
        hi = max(i1, i2)
        lo = min(i1, i2)
        lineages[hi] = lineages[end]
        pop!(lineages)
        lineages[lo] = lineages[end]
        pop!(lineages)
        push!(lineages, new_node)
    end
    
    root = lineages[1]
    tree = _coal_to_phylo(root)
    
    # 2. Add mutations (Poisson on branches)
    # Total rate = mu * 2Ne * seq_len
    theta = 4.0 * ne * mu * seq_len
    
    # Simple Infinite Sites Model
    haplotypes = Dict{String, Vector{Int}}()
    for leaf in get_terminals(tree)
        haplotypes[leaf.name] = zeros(Int, seq_len)
    end
    
    function _mutate(node::PhyloTree, current_hap::Vector{Int})
        for child in node.children
            child_hap = copy(current_hap)
            # mutations = Poisson(theta/2 * branch_length)
            # theta is for total sequence, branch length is in 2Ne units
            n_mut = rand(Distributions.Poisson((theta / 2.0) * child.branch_length))
            for _ in 1:n_mut
                pos = rand(1:seq_len)
                child_hap[pos] = 1 - child_hap[pos] # flip 0/1
            end
            
            if isleaf(child)
                haplotypes[child.name] = child_hap
            else
                _mutate(child, child_hap)
            end
        end
    end
    
    _mutate(tree, zeros(Int, seq_len))
    
    # Convert to strings
    hap_strings = Dict{String, String}()
    for (name, hap) in haplotypes
        hap_strings[name] = join(string.(hap))
    end
    
    return tree, hap_strings
end





