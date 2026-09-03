# ==============================================================================
# popgen.jl — Population genetics analysis
#
# Provides allele/genotype frequency estimation, Hardy-Weinberg equilibrium
# testing (chi-squared and exact), F-statistics, AMOVA, linkage
# disequilibrium, genetic distance matrices, effective population size
# estimation, and migration rate inference.
#
# References:
#   - Weir & Cockerham (1984) Evolution 38(6):1358-1370 (F-statistics)
#   - Excoffier et al. (1992) Genetics 131(2):479-491 (AMOVA)
#   - Hill (1981) Heredity 47(2):229-239 (LD-based Ne estimation)
#   - Wigginton et al. (2005) AJHG 76(5):887-893 (exact HW test)
# ==============================================================================
using Distributions, Random, Statistics
using ..BioToolkit: ProvenanceParams, ThreadSafeProvenanceContext, new_provenance_id

@inline function _register_popgen_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
  return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end


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
export amova, f_statistics, hardy_weinberg_exact, linkage_disequilibrium, infer_phase_em, estimate_ne_ld, ld_mapping, ld_decay, migration_rate, g_statistics, genetic_distance
export mantel_test, population_pcoa, mismatch_distribution
export read_genepop_record, split_in_pops, split_in_loci, remove_population!, remove_locus_by_position!, remove_locus_by_name!
export ewens_watterson_test, sweepfinder_clr, genetic_relationship_matrix, inbreeding_coefficient, relatedness, linear_mixed_model_scan
export populations, loci, samplenames, missingdata, richness, alleleaverage, pairwise_fst, PairwiseFSTResult, summary_statistics
# --- PopGen.jl-parity exports ---
export weir_cockerham_fst, hudson_fst, fst_permutation_test
export kinship_queller_goodnight, kinship_ritland, kinship_lynch_li, kinship_lynch_ritland
export kinship_li_horvitz, kinship_moran, kinship_blouin, kinship_loiselle
export pairwise_kinship, sample_heterozygosity
export population_kmeans, population_hclust, population_kmedoids, population_fuzzycmeans, population_dbscan, population_cluster
export pairwise_identical, read_structure_record, write_structure

const DEFAULT_MISSING_ALLELES = (missing, 0, "0", '?', "?", '.', "NA", "N", '-')

@inline _is_missing_allele(allele, missing_alleles) =
  ismissing(allele) || any(marker -> isequal(allele, marker), missing_alleles)

@inline function _checked_locus(locus_idx::Integer)
  locus_idx >= 1 || throw(ArgumentError("locus_idx must be positive; got $locus_idx"))
  return Int(locus_idx)
end

@inline function _valid_diploid(locus::Locus, missing_alleles)
  return length(locus.alleles) == 2 &&
         !_is_missing_allele(locus.alleles[1], missing_alleles) &&
         !_is_missing_allele(locus.alleles[2], missing_alleles)
end

# ------------------------------------------------------------------------------
# 2. Basic Frequencies
# ------------------------------------------------------------------------------

"""
    allele_frequencies(pop::Population, locus_idx::Int)

Calculates the frequency of each allele at a specific locus across the population.
Returns a Dict mapping allele -> frequency.
"""
function allele_frequencies(pop::Population{T}, locus_idx::Integer;
                            missing_alleles=DEFAULT_MISSING_ALLELES) where T
  locus_idx = _checked_locus(locus_idx)
  counts = Dict{T,Int}()
  total_alleles = 0
  for ind in pop.individuals
    locus_idx <= length(ind.loci) || continue
    for allele in ind.loci[locus_idx].alleles
      _is_missing_allele(allele, missing_alleles) && continue
      counts[allele] = get(counts, allele, 0) + 1
      total_alleles += 1
    end
  end
  total_alleles == 0 && return Dict{T,Float64}()
  return Dict(a => c / Float64(total_alleles) for (a, c) in counts)
end

"""
    genotype_frequencies(pop::Population, locus_idx::Int)

Calculates the frequency of observed genotypes at a specific locus.
Assumes unordered alleles (e.g. A/B is the same as B/A).
"""
function genotype_frequencies(pop::Population{T}, locus_idx::Integer;
                              missing_alleles=DEFAULT_MISSING_ALLELES) where T
  locus_idx = _checked_locus(locus_idx)
  counts = Dict{Tuple{Vararg{T}},Int}()
  total_genos = 0
  for ind in pop.individuals
    locus_idx <= length(ind.loci) || continue
    alleles = ind.loci[locus_idx].alleles
    any(a -> _is_missing_allele(a, missing_alleles), alleles) && continue
    # Preserve allele dosage; sort by string representation to avoid hash collisions
    genotype = Tuple(sort(collect(alleles), by=string))
    counts[genotype] = get(counts, genotype, 0) + 1
    total_genos += 1
  end
  total_genos == 0 && return Dict{Tuple{Vararg{T}},Float64}()
  return Dict(g => c / Float64(total_genos) for (g, c) in counts)
end

# ------------------------------------------------------------------------------
# 3. Heterozygosity
# ------------------------------------------------------------------------------

"""
    heterozygosity_observed(pop::Population, locus_idx::Int)

Returns the observed proportion of heterozygotes at a given locus.
"""
function heterozygosity_observed(pop::Population{T}, locus_idx::Integer;
                                 missing_alleles=DEFAULT_MISSING_ALLELES) where T
  locus_idx = _checked_locus(locus_idx)
  heterozygotes = total = 0
  for ind in pop.individuals
    locus_idx <= length(ind.loci) || continue
    locus = ind.loci[locus_idx]
    _valid_diploid(locus, missing_alleles) || continue
    heterozygotes += locus.alleles[1] != locus.alleles[2]
    total += 1
  end
  return total == 0 ? 0.0 : heterozygotes / Float64(total)
end

"""
    heterozygosity_expected(pop::Population, locus_idx::Int)

Returns the expected heterozygosity (Nei's gene diversity) at a given locus,
calculated as 1 - sum(p_i^2) where p_i are the allele frequencies.
"""
function heterozygosity_expected(pop::Population{T}, locus_idx::Integer;
                                 missing_alleles=DEFAULT_MISSING_ALLELES,
                                 unbiased::Bool=false) where T
  freqs = allele_frequencies(pop, locus_idx; missing_alleles=missing_alleles)
  he = 1.0 - sum(p^2 for p in values(freqs))
  if !unbiased
    return he
  end
  n = sum(1 for ind in pop.individuals if locus_idx <= length(ind.loci) &&
          _valid_diploid(ind.loci[locus_idx], missing_alleles)) * 2
  return n > 1 ? he * n / (n - 1) : 0.0
end

# ------------------------------------------------------------------------------
# 4. Hardy-Weinberg Equilibrium
# ------------------------------------------------------------------------------

"""
    hardy_weinberg_test(pop::Population, locus_idx::Int)

Performs a standard Chi-Square test for Hardy-Weinberg equilibrium at a bi-allelic locus.
Returns the p-value.
"""
function _biallelic_genotype_counts(pop::Population{T}, locus_idx::Integer, missing_alleles) where T
  locus_idx = _checked_locus(locus_idx)
  counts = Dict{T,Int}()
  genotypes = Vector{Tuple{T,T}}()
  for ind in pop.individuals
    locus_idx <= length(ind.loci) || continue
    locus = ind.loci[locus_idx]
    _valid_diploid(locus, missing_alleles) || continue
    a, b = locus.alleles
    counts[a] = get(counts, a, 0) + 1
    counts[b] = get(counts, b, 0) + 1
    push!(genotypes, (a, b))
  end
  length(counts) == 2 || return nothing
  alleles = collect(keys(counts))
  n11 = n22 = n12 = 0
  for (a, b) in genotypes
    if a == alleles[1] && b == alleles[1]
      n11 += 1
    elseif a == alleles[2] && b == alleles[2]
      n22 += 1
    else
      n12 += 1
    end
  end
  return n11, n22, n12
end

function hardy_weinberg_test(pop::Population{T}, locus_idx::Integer;
                             missing_alleles=DEFAULT_MISSING_ALLELES) where T
  observed = _biallelic_genotype_counts(pop, locus_idx, missing_alleles)
  observed === nothing && return 1.0
  obs_p2, obs_q2, obs_pq = observed
  total_genos = obs_p2 + obs_q2 + obs_pq
  total_genos == 0 && return 1.0
  p = (2obs_p2 + obs_pq) / (2.0 * total_genos)
  q = 1.0 - p
  exp_p2, exp_q2, exp_pq = p^2 * total_genos, q^2 * total_genos, 2p * q * total_genos
  if exp_p2 == 0 || exp_q2 == 0 || exp_pq == 0
    return 1.0  # chi-square undefined for zero expected counts; fail to reject H0
  end
  chi2 = (obs_p2 - exp_p2)^2 / exp_p2 +
         (obs_q2 - exp_q2)^2 / exp_q2 +
         (obs_pq - exp_pq)^2 / exp_pq
  return ccdf(Chisq(1), chi2)
end

"""
    hardy_weinberg_exact(pop::Population, locus_idx::Int)

Slatkin's Exact Test for Hardy-Weinberg Equilibrium.
Useful for small sample sizes where Chi-Square is inaccurate.
"""
function hardy_weinberg_exact(pop::Population{T}, locus_idx::Integer;
                              missing_alleles=DEFAULT_MISSING_ALLELES) where T
  observed = _biallelic_genotype_counts(pop, locus_idx, missing_alleles)
  observed === nothing && return 1.0
  n11, n22, n12 = observed
  n = n11 + n22 + n12
  n == 0 && return 1.0
  n1, n2 = 2n11 + n12, 2n22 + n12
  log_prob(h) = begin
    aa, bb = (n1 - h) ÷ 2, (n2 - h) ÷ 2
    loggamma(n + 1) + h * log(2.0) + loggamma(n1 + 1) + loggamma(n2 + 1) -
      loggamma(aa + 1) - loggamma(bb + 1) - loggamma(h + 1) - loggamma(2n + 1)
  end
  configurations = collect((n1 % 2):2:min(n1, n2))
  log_probs = log_prob.(configurations)
  max_logp = maximum(log_probs)
  weights = exp.(log_probs .- max_logp)
  observed_logp = log_prob(n12)
  numerator = sum(w for (lp, w) in zip(log_probs, weights) if lp <= observed_logp + 1e-12)
  return min(numerator / sum(weights), 1.0)
end

"""
    hardy_weinberg_test(data::Vector{Population}, locus_idx::Int)

Run Hardy-Weinberg Chi-square test across all populations for a given locus.
Returns a Dict mapping population name -> p-value.
"""
function hardy_weinberg_test(data::AbstractVector{<:Population}, locus_idx::Integer; missing_alleles=DEFAULT_MISSING_ALLELES)
  Dict(pop.name => hardy_weinberg_test(pop, locus_idx; missing_alleles=missing_alleles) for pop in data)
end

"""
    hardy_weinberg_exact(data::Vector{Population}, locus_idx::Int)

Run Hardy-Weinberg Exact test across all populations for a given locus.
Returns a Dict mapping population name -> p-value.
"""
function hardy_weinberg_exact(data::AbstractVector{<:Population}, locus_idx::Integer; missing_alleles=DEFAULT_MISSING_ALLELES)
  Dict(pop.name => hardy_weinberg_exact(pop, locus_idx; missing_alleles=missing_alleles) for pop in data)
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
  global_allele_counts = Dict{T,Int}()

  subpop_Ho = Float64[]
  subpop_Hs = Float64[]
  subpop_weights = Float64[]

  for pop in populations
    pop_inds = 0
    hts = 0
    pop_counts = Dict{T,Int}()

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

  result = (F_IS, F_ST, F_IT)
  _ctx = active_provenance_context()
  return _register_popgen_result!(_ctx, result, "f_statistics"; parameters=(n_populations=length(populations), locus_idx=locus_idx, F_ST=F_ST))
end

"""
    g_statistics(populations::Vector{Population}, locus_idx::Int)

Calculates G_ST (Nei's GST) and Jost's D.
Returns a tuple (G_ST, Jost_D)
"""
function g_statistics(populations::Vector{Population{T}}, locus_idx::Int) where T
  # Heterozygosity components
  num_pops = length(populations)
  if num_pops < 2
    ;
    return (0.0, 0.0);
  end

  total_inds = 0
  Ho_vals = Float64[]
  Hs_vals = Float64[]

  global_counts = Dict{T,Int}()

  for pop in populations
    n = 0
    hts = 0
    counts = Dict{T,Int}()
    for ind in pop.individuals
      if locus_idx <= length(ind.loci)
        l = ind.loci[locus_idx]
        if length(l.alleles) == 2
          n += 1
          if l.alleles[1] != l.alleles[2]
            ;
            hts += 1;
          end
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
  if F_ST <= 0.0
    ;
    return Inf;
  end
  if F_ST >= 1.0
    ;
    return 0.0;
  end
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

  result = (Phi_ST, Var_within, Var_among)
  _ctx = active_provenance_context()
  return _register_popgen_result!(_ctx, result, "amova"; parameters=(n_populations=length(pop_sizes), n_individuals=N, Phi_ST=Phi_ST))
end

# ------------------------------------------------------------------------------
# 7. Linkage Disequilibrium
# ------------------------------------------------------------------------------

"""
    linkage_disequilibrium(pop::Population, locus_1::Int, locus_2::Int)

Calculate D, D', and r^2 between two bi-allelic loci.
"""
function linkage_disequilibrium(pop::Population{T}, locus_1::Integer, locus_2::Integer;
                               missing_alleles=DEFAULT_MISSING_ALLELES,
                               max_iter::Int=200, tol::Real=1e-8) where T
  locus_1 = _checked_locus(locus_1)
  locus_2 = _checked_locus(locus_2)
  locus_1 == locus_2 && throw(ArgumentError("LD requires two distinct loci"))
  haplotypes = infer_phase_em(pop, locus_1, locus_2; max_iter=max_iter, tol=Float64(tol), missing_alleles=missing_alleles)
  isempty(haplotypes) && return (0.0, 0.0, 0.0)
  alleles_1 = unique(first(h) for h in keys(haplotypes))
  alleles_2 = unique(last(h) for h in keys(haplotypes))
  (length(alleles_1) == 2 && length(alleles_2) == 2) || return (0.0, 0.0, 0.0)
  A, B = alleles_1[1], alleles_2[1]
  pA = sum(freq for ((a, _), freq) in haplotypes if a == A)
  pB = sum(freq for ((_, b), freq) in haplotypes if b == B)
  pAB = get(haplotypes, (A, B), 0.0)
  D = pAB - pA * pB
  dmax = D >= 0.0 ? min(pA * (1.0 - pB), (1.0 - pA) * pB) :
                        min(pA * pB, (1.0 - pA) * (1.0 - pB))
  denominator = pA * (1.0 - pA) * pB * (1.0 - pB)
  return (D, dmax == 0.0 ? 0.0 : D / dmax, denominator == 0.0 ? 0.0 : D^2 / denominator)
end

"""
    ld_mapping(pop::Population, locus_indices::Vector{Int}, window_size::Int=10)

Calculates LD (r^2) between all pairs of loci within a sliding window.
Returns a Dict mapping (locus1, locus2) -> r_squared.
"""
function ld_mapping(pop::Population{T}, locus_indices::Vector{Int}, window_size::Int=10) where T
  results = Dict{Tuple{Int,Int},Float64}()
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
function ld_decay(pop::Population{T}, physical_distances::Vector{Float64}, locus_pairs::Vector{Tuple{Int,Int}}) where T
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
function infer_phase_em(pop::Population{T}, locus_1::Integer, locus_2::Integer;
                        max_iter::Int=100, tol::Float64=1e-6,
                        missing_alleles=DEFAULT_MISSING_ALLELES) where T
  max_iter > 0 || throw(ArgumentError("max_iter must be positive"))
  tol > 0 || throw(ArgumentError("tol must be positive"))
  locus_1, locus_2 = _checked_locus(locus_1), _checked_locus(locus_2)
  alleles_1, alleles_2 = Set{T}(), Set{T}()
  raw_genotypes = Vector{NTuple{4,T}}()
  for ind in pop.individuals
    (locus_1 <= length(ind.loci) && locus_2 <= length(ind.loci)) || continue
    g1, g2 = ind.loci[locus_1], ind.loci[locus_2]
    (_valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)) || continue
    a, b = g1.alleles
    c, d = g2.alleles
    push!(raw_genotypes, (a, b, c, d))
    union!(alleles_1, (a, b)); union!(alleles_2, (c, d))
  end
  isempty(raw_genotypes) && return Dict{Tuple{T,T},Float64}()
  haplotypes = [(a, b) for a in alleles_1 for b in alleles_2]
  hap_index = Dict(h => i for (i, h) in enumerate(haplotypes))
  possible = Vector{Vector{Tuple{Int,Int}}}(undef, length(raw_genotypes))
  for (k, (a, b, c, d)) in enumerate(raw_genotypes)
    pairs = Tuple{Int,Int}[]
    for (h1, h2) in (((a, c), (b, d)), ((a, d), (b, c)))
      i, j = hap_index[h1], hap_index[h2]
      pair = i <= j ? (i, j) : (j, i)
      pair in pairs || push!(pairs, pair)
    end
    possible[k] = pairs
  end
  frequencies = fill(inv(Float64(length(haplotypes))), length(haplotypes))
  next_frequencies = similar(frequencies)
  for _ in 1:max_iter
    fill!(next_frequencies, 0.0)
    for pairs in possible
      weights = [frequencies[i] * frequencies[j] * (i == j ? 1.0 : 2.0) for (i, j) in pairs]
      normalizer = sum(weights)
      normalizer > 0.0 || continue
      for ((i, j), weight) in zip(pairs, weights)
        posterior = weight / normalizer
        next_frequencies[i] += posterior
        next_frequencies[j] += posterior
      end
    end
    next_frequencies ./= 2.0 * length(possible)
    maximum(abs.(next_frequencies .- frequencies)) < tol && (frequencies .= next_frequencies; break)
    frequencies .= next_frequencies
  end
  return Dict(h => frequencies[i] for (i, h) in enumerate(haplotypes))
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
function genetic_distance(pop1::Population{T}, pop2::Population{T}, locus_idx::Integer;
                          method::Symbol=:nei, missing_alleles=DEFAULT_MISSING_ALLELES) where T
  f1 = allele_frequencies(pop1, locus_idx; missing_alleles=missing_alleles)
  f2 = allele_frequencies(pop2, locus_idx; missing_alleles=missing_alleles)
  alleles = union(keys(f1), keys(f2))
  isempty(alleles) && return NaN
  p = (get(f1, a, 0.0) for a in alleles)
  q = (get(f2, a, 0.0) for a in alleles)
  if method === :rogers
    return sqrt(sum((x - y)^2 for (x, y) in zip(p, q)) / 2.0)
  elseif method === :nei
    j11 = sum(x^2 for x in p)
    j22 = sum(y^2 for y in q)
    j12 = sum(x * y for (x, y) in zip(p, q))
    denom = sqrt(j11 * j22)
    denom == 0.0 && return Inf  # monomorphic populations
    identity = j12 / denom
    return identity <= 0.0 ? Inf : -log(identity)
  elseif method === :cavalli_sforza
    return sqrt(max(0.0, 2.0 * (1.0 - sum(sqrt(x * y) for (x, y) in zip(p, q)))))
  elseif method === :reynolds
    numerator = sum((x - y)^2 for (x, y) in zip(p, q))
    denominator = 2.0 * sum(1.0 - x * y for (x, y) in zip(p, q))
    return denominator == 0.0 ? 0.0 : -log1p(-min(numerator / denominator, 1.0))
  end
  throw(ArgumentError("Unknown genetic distance method: $method. Supported methods are :nei, :reynolds, :cavalli_sforza, and :rogers."))
end

function population_pca(populations::AbstractVector{<:Population{T}}, locus_indices::AbstractVector{<:Integer};
                        missing_alleles=DEFAULT_MISSING_ALLELES) where T
  npopulations = length(populations)
  npopulations > 0 || return (zeros(Float64, 0, 0), Float64[], zeros(Float64, 0, 0))
  selected_loci = unique(_checked_locus.(locus_indices))
  features = Tuple{Int,T}[]
  for locus_idx in selected_loci
    alleles = Set{T}()
    for pop in populations
      union!(alleles, keys(allele_frequencies(pop, locus_idx; missing_alleles=missing_alleles)))
    end
    append!(features, ((locus_idx, allele) for allele in sort!(collect(alleles), by=hash)))
  end
  nfeatures = length(features)
  frequency_matrix = zeros(Float64, npopulations, nfeatures)
  for (row, pop) in enumerate(populations)
    by_locus = Dict{Int,Dict{T,Float64}}()
    for locus_idx in selected_loci
      by_locus[locus_idx] = allele_frequencies(pop, locus_idx; missing_alleles=missing_alleles)
    end
    for (column, (locus_idx, allele)) in enumerate(features)
      frequency_matrix[row, column] = get(by_locus[locus_idx], allele, 0.0)
    end
  end
  npopulations == 1 && return (zeros(Float64, 1, 0), Float64[], frequency_matrix)
  centered = frequency_matrix .- mean(frequency_matrix, dims=1)
  decomposition = svd(centered; full=false)
  projections = decomposition.U * Diagonal(decomposition.S)
  eigenvalues = decomposition.S .^ 2 ./ (npopulations - 1)
  result = (projections, eigenvalues, frequency_matrix)
  _ctx = active_provenance_context()
  return _register_popgen_result!(_ctx, result, "population_pca"; parameters=(n_populations=npopulations, n_loci=length(selected_loci)))
end

"""
    population_pcoa(dist_matrix::Matrix{Float64})

Principal Coordinate Analysis (Classical Multidimensional Scaling).
Transforms a distance matrix into an ordination plot.
Returns (coordinates, eigenvalues).
"""
function population_pcoa(dist_matrix::AbstractMatrix{<:Real})
  n, m = size(dist_matrix)
  n == m || throw(ArgumentError("dist_matrix must be square"))
  n > 0 || return (zeros(Float64, 0, 0), Float64[])
  all(isfinite, dist_matrix) || throw(ArgumentError("dist_matrix must contain only finite values"))
  isapprox(dist_matrix, transpose(dist_matrix); atol=sqrt(eps(Float64)), rtol=sqrt(eps(Float64))) ||
    throw(ArgumentError("dist_matrix must be symmetric"))
  all(iszero, diag(dist_matrix)) || throw(ArgumentError("dist_matrix diagonal must be zero"))
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

  _ctx = active_provenance_context()
  return provenance_result!(_ctx, (coords, evals[pos_idx]), "population_pcoa")
end

"""
    mantel_test(dist_matrix1::Matrix{Float64}, dist_matrix2::Matrix{Float64}; permutations::Int=999)

Performs a Mantel test between two distance matrices to test for correlation.
Returns (correlation, p_value).
"""
function mantel_test(dist_matrix1::AbstractMatrix{<:Real}, dist_matrix2::AbstractMatrix{<:Real};
                     permutations::Integer=999, rng::AbstractRNG=Random.default_rng(), two_sided::Bool=false)
  size(dist_matrix1) == size(dist_matrix2) || throw(ArgumentError("distance matrices must have the same dimensions"))
  n, m = size(dist_matrix1)
  n == m || throw(ArgumentError("distance matrices must be square"))
  permutations >= 0 || throw(ArgumentError("permutations must be non-negative"))
  n < 3 && return (0.0, 1.0)
  nvalues = n * (n - 1) ÷ 2
  observed_1, observed_2, permuted_1 = Vector{Float64}(undef, nvalues), Vector{Float64}(undef, nvalues), Vector{Float64}(undef, nvalues)
  index = 1
  for row in 2:n, column in 1:(row - 1)
    observed_1[index] = dist_matrix1[row, column]
    observed_2[index] = dist_matrix2[row, column]
    index += 1
  end
  observed = cor(observed_1, observed_2)
  isfinite(observed) || return (0.0, 1.0)
  permutation = collect(1:n)
  extreme = 0
  for _ in 1:permutations
    Random.shuffle!(rng, permutation)
    index = 1
    for row in 2:n, column in 1:(row - 1)
      permuted_1[index] = dist_matrix1[permutation[row], permutation[column]]
      index += 1
    end
    statistic = cor(permuted_1, observed_2)
    if two_sided ? abs(statistic) >= abs(observed) : statistic >= observed
      extreme += 1
    end
  end
  result = (observed, (extreme + 1) / (permutations + 1))
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, result, "mantel_test")
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
function estimate_ne_ld(pop::Population{T}, locus_pairs::Vector{Tuple{Int,Int}}) where T
  n = length(pop.individuals) * 2 # total alleles (diploid)
  if n < 2
    ;
    return Inf;
  end

  r2_sum = 0.0
  valid_pairs = 0
  for (l1, l2) in locus_pairs
    _, _, r2 = linkage_disequilibrium(pop, l1, l2)
    r2_sum += r2
    valid_pairs += 1
  end

  if valid_pairs == 0
    ;
    return Inf;
  end
  mean_r2 = r2_sum / valid_pairs

  # Correct for sample size
  # r2_corrected = r2 - 1/n
  # Ne = 1 / (3 * r2_corrected)
  r2_corrected = mean_r2 - (1.0 / n)
  if r2_corrected <= 0.0
    ;
    return Inf;
  end

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
  if n < 2
    ;
    return 0;
  end
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
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, S, "segregating_sites")
end

"""
    mismatch_distribution(alignment::MultipleSequenceAlignment)

Calculates the distribution of pairwise nucleotide differences between all pairs
of sequences in the alignment. Returns a vector of difference counts.
"""
function mismatch_distribution(alignment::MultipleSequenceAlignment)
  records = alignment.records
  n = length(records)
  if n < 2
    ;
    return Int[];
  end

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
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, mismatches, "mismatch_distribution")
end

"""
    nucleotide_diversity(alignment::MultipleSequenceAlignment)

Calculates nucleotide diversity (π), the average number of nucleotide differences
per site between two continuous sequences randomly drawn from the population.
"""
function nucleotide_diversity(alignment::MultipleSequenceAlignment)
  records = alignment.records
  n = length(records)
  if n < 2
    ;
    return 0.0;
  end

  len = length(records[1].sequence)
  pi_sum = 0.0
  valid_sites = 0

  for i in 1:len
    freqs = Dict{Char,Int}()
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

  _ctx = active_provenance_context()
  return provenance_result!(_ctx, valid_sites > 0 ? (pi_sum / valid_sites) : 0.0, "nucleotide_diversity")
end

"""
    watterson_theta(alignment::MultipleSequenceAlignment)

Estimates Watterson's Theta (Θ) from the number of segregating sites.
"""
function watterson_theta(alignment::MultipleSequenceAlignment)
  n = length(alignment.records)
  if n < 2
    ;
    return 0.0;
  end

  S = segregating_sites(alignment)
  a1 = sum(1.0 / i for i in 1:(n-1))

  _ctx = active_provenance_context()
  return provenance_result!(_ctx, S / a1, "watterson_theta")
end

"""
    tajimas_d(alignment::MultipleSequenceAlignment)

Calculates Tajima's D statistic to test for neutral evolution.
Negative values suggest selective sweeping or population expansion.
Positive values suggest balancing selection or population subdivision.
"""
function tajimas_d(alignment::MultipleSequenceAlignment)
  n = length(alignment.records)
  if n < 2
    ;
    return 0.0;
  end

  S = segregating_sites(alignment)
  if S == 0
    ;
    return 0.0;
  end

  # Calculate π and Θ unscaled by total length (i.e. number of pairwise differences, not per-site)
  # We instead scale nucleotides diversity by sequence chunks
  # Wait, nucleotide_diversity(ali) gives per-site. We need total pairwise differences for Tajima's D.

  len = 0
  records = alignment.records
  len_records = length(records[1].sequence)
  total_pairwise_diffs = 0.0

  for i in 1:len_records
    freqs = Dict{Char,Int}()
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
  _ctx = active_provenance_context()
  return _register_popgen_result!(_ctx, D, "tajimas_d"; parameters=(n_sequences=n, n_segregating_sites=S, tajimas_d=D))
end

"""
    site_frequency_spectrum(alignment::MultipleSequenceAlignment; folded::Bool=true)

Calculates the Site Frequency Spectrum (SFS) for an alignment. If `folded` is true,
the minor allele frequency spectrum is generated. Returns a frequency distribution vector.
"""
function site_frequency_spectrum(alignment::MultipleSequenceAlignment; folded::Bool=true)
  n = length(alignment.records)
  if n < 2
    ;
    return Float64[];
  end

  lim = folded ? floor(Int, n / 2) : (n - 1)
  sfs = zeros(Int, lim)

  len_records = length(alignment.records[1].sequence)
  for i in 1:len_records
    freqs = Dict{Char,Int}()
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
      if folded
        allele_count = min(counts[1], counts[2])
      else
        # unfolded mode requires an outgroup to determine ancestral allele.
        # Without an outgroup, folded=false behaves identically to folded=true.
        # To use true unfolded SFS, provide an outgroup sequence.
        allele_count = min(counts[1], counts[2])
      end
      if allele_count > 0 && allele_count <= lim
        sfs[allele_count] += 1
      end
    end
  end
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, sfs, "site_frequency_spectrum")
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
  n < 2 && return 0.0

  S = segregating_sites(alignment)
  S == 0 && return 0.0

  # Calculate ηs (number of singletons)
  singletons = 0
  records = alignment.records
  len = length(records[1].sequence)
  for i in 1:len
    freqs = Dict{Char,Int}()
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
  a2 = sum(1.0 / (i^2) for i in 1:(n-1))

  # Exact Fu & Li D variance components (without outgroup)
  u_D = n / (n - 1.0) - 1.0 / a1
  v_D = 1.0 + (n / (n - 1.0))^2 * (a2 - 1.0) / (a1^2) - 1.0 / a1

  var_D = u_D * S + v_D * S^2
  if var_D <= 0
    return 0.0
  end

  D = (S - a1 * singletons) / sqrt(var_D)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, D, "fu_li_d")
end

"""
    fu_li_f(alignment::MultipleSequenceAlignment)

Calculates Fu and Li's F statistic (without an outgroup).
"""
function fu_li_f(alignment::MultipleSequenceAlignment)
  n = length(alignment.records)
  n < 2 && return 0.0

  S = segregating_sites(alignment)
  S == 0 && return 0.0

  # Count singletons and total pairwise differences (pi)
  singletons = 0
  records = alignment.records
  len = length(records[1].sequence)
  total_pairwise_diffs = 0.0

  for i in 1:len
    freqs = Dict{Char,Int}()
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
      if length(freqs) == 2
        counts = collect(values(freqs))
        if any(c == 1 for c in counts)
          singletons += 1
        end
      end
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

  a1 = sum(1.0 / i for i in 1:(n-1))
  a2 = sum(1.0 / (i^2) for i in 1:(n-1))

  # Exact Fu & Li F variance components (without outgroup)
  u_F = (4.0 * n - 6.0) / (n - 1.0) - 1.0 / a1
  v_F = 1.0 + (n / (n - 1.0))^2 * a2 - 2.0 / (n - 1.0) - 1.0 / a1

  var_F = u_F * S + v_F * S^2
  if var_F <= 0
    return 0.0
  end

  F = (total_pairwise_diffs - singletons) / sqrt(var_F)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, F, "fu_li_f")
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
  if n == 0
    ;
    return (obs_F, 1.0);
  end

  k = length(freqs)
  if k <= 1
    ;
    return (obs_F, 1.0);
  end

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
  if n < 2
    ;
    return 0.0;
  end

  len = length(alignment.records[1].sequence)
  end_site = min(core_site + distance_sites, len)
  start_site = max(core_site - distance_sites, 1)

  haplotypes = [records.sequence[start_site:end_site] for records in alignment.records]

  counts = Dict{String,Int}()
  for hap in haplotypes
    counts[hap] = get(counts, hap, 0) + 1
  end

  sum_sq = sum(c * (c - 1) for c in values(counts))
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, sum_sq / (n * (n - 1)), "ehh")
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
    if e < 0.05
      ;
      break;
    end
    dist += 1
  end
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, integral, "ihs")
end

"""
    xp_ehh(pop1_ali::MultipleSequenceAlignment, pop2_ali::MultipleSequenceAlignment, core_site::Int)

Calculates Cross-Population EHH (XP-EHH) to detect selection between two populations.
"""
function xp_ehh(pop1_ali::MultipleSequenceAlignment, pop2_ali::MultipleSequenceAlignment, core_site::Int)
  ih1 = ihs(pop1_ali, core_site)
  ih2 = ihs(pop2_ali, core_site)

  if ih2 == 0
    ;
    return 0.0;
  end
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, log(ih1 / ih2), "xp_ehh")
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
  if isempty(sfs)
    ;
    return zeros(grid_points);
  end

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
      col_freqs = Dict{Char,Int}()
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

  _ctx = active_provenance_context()
  return provenance_result!(_ctx, clrs, "sweepfinder_clr")
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
  isempty(alleles) && return 0.0

  # Sum over all alleles: f3 = sum((p3 - p1) * (p3 - p2))
  val = sum(a -> (get(f3, a, 0.0) - get(f1, a, 0.0)) * (get(f3, a, 0.0) - get(f2, a, 0.0)), alleles)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, val, "f3_statistic")
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
  isempty(alleles) && return 0.0

  # Sum over all alleles: f4 = sum((p1 - p2) * (p3 - p4))
  val = sum(a -> (get(f1, a, 0.0) - get(f2, a, 0.0)) * (get(f3, a, 0.0) - get(f4, a, 0.0)), alleles)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, val, "f4_statistic")
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
  length(alleles) != 2 && return 0.0 # ABBA-BABA requires exactly 2 alleles

  a = first(alleles)
  p1_a = get(f1, a, 0.0)
  p2_a = get(f2, a, 0.0)
  p3_a = get(f3, a, 0.0)
  p4_a = get(f4, a, 0.0)

  abba = (1.0 - p1_a) * p2_a * p3_a * (1.0 - p4_a)
  baba = p1_a * (1.0 - p2_a) * p3_a * (1.0 - p4_a)

  if abba + baba == 0
    return 0.0
  end
  val = (abba - baba) / (abba + baba)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, val, "patterson_d")
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

    if p <= 0.0 || p >= 1.0
      ;
      break;
    end
  end
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, traj, "wright_fisher_simulation")
end

"""
    wright_fisher_metapopulation(N::Vector{Int}, generations::Int, M::Matrix{Float64}; mu::Float64=1e-8, s::Vector{Float64}=zeros(length(N)))

Simulates a metapopulation under the Wright-Fisher model with a migration matrix M.
M[i, j] is the proportion of population i that comes from population j.
"""
function wright_fisher_metapopulation(N::Vector{Int}, generations::Int, M::Matrix{Float64}; mu::Float64=1e-8, s::Vector{Float64}=nothing)
  num_pops = length(N)
  if s === nothing
    ;
    s = zeros(num_pops);
  end

  trajs = [zeros(generations + 1) for _ in 1:num_pops]
  p = fill(0.5, num_pops)
  for k in 1:num_pops
    ;
    trajs[k][1] = p[k];
  end

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
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, trajs, "wright_fisher_metapopulation")
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
        allele_2 = parse(Int, genotype_token[(token_midpoint+1):end])
        push!(loci, Locus{Int}((allele_1, allele_2)))
      end
      push!(record.populations[end].individuals, PopGenIndividual{Int}(individual_name, loci))
    end
    current_line += 1
  end

  return record
end

function read_genepop_record(filepath::String)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, _parse_genepop_population(filepath), "read_genepop_record")
end

"""
    read_genepop(filepath::String)

Parses a GenePop format file into a vector of Population objects.
"""
function read_genepop(filepath::String)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, read_genepop_record(filepath).populations, "read_genepop")
end

function split_in_pops(record::GenePopRecord{T}, pop_names::Vector{String}) where T
  if length(pop_names) != length(record.populations)
    throw(ArgumentError("pop_names must match the number of populations"))
  end

  result = Dict{String,GenePopRecord{T}}()
  for (index, population) in enumerate(record.populations)
    new_record = GenePopRecord{T}(record.marker_len, record.comment_line, copy(record.loci_list), [pop_names[index]], [deepcopy(population)])
    result[pop_names[index]] = new_record
  end
  return result
end

function split_in_loci(record::GenePopRecord{T}) where T
  result = Dict{String,GenePopRecord{T}}()
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
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, record, "remove_population!")
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

function remove_locus_by_name!(record::GenePopRecord{T}, name::String) where T
  for (index, locus_name) in enumerate(record.loci_list)
    if locus_name == name
      return remove_locus_by_position!(record, index - 1)
    end
  end
  return record
end

remove_population(record::GenePopRecord, pos::Int) = remove_population!(record, pos)
remove_locus_by_position(record::GenePopRecord, pos::Int) = remove_locus_by_position!(record, pos)
remove_locus_by_name(record::GenePopRecord, name::String) = remove_locus_by_name!(record, name)

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
function genetic_relationship_matrix(populations::AbstractVector{<:Population{T}};
                                     missing_alleles=DEFAULT_MISSING_ALLELES) where T
  individuals = [individual for population in populations for individual in population.individuals]
  nindividuals = length(individuals)
  nindividuals == 0 && return zeros(Float64, 0, 0)
  nloci = maximum((length(individual.loci) for individual in individuals); init=0)
  centered = zeros(Float64, nindividuals, nloci)
  nused = 0
  denominator = 0.0
  for locus_idx in 1:nloci
    counts = Dict{T,Int}()
    observed = falses(nindividuals)
    for (row, individual) in enumerate(individuals)
      locus_idx <= length(individual.loci) || continue
      genotype = individual.loci[locus_idx]
      _valid_diploid(genotype, missing_alleles) || continue
      a, b = genotype.alleles
      counts[a] = get(counts, a, 0) + 1
      counts[b] = get(counts, b, 0) + 1
      observed[row] = true
    end
    length(counts) == 2 || continue
    reference = first(sort!(collect(keys(counts)), by=allele -> (-counts[allele], hash(allele))))
    ncalled = count(observed)
    p = counts[reference] / (2.0 * ncalled)
    (p == 0.0 || p == 1.0) && continue
    nused += 1
    denominator += 2.0 * p * (1.0 - p)
    for (row, individual) in enumerate(individuals)
      if observed[row]
        centered[row, nused] = count(==(reference), individual.loci[locus_idx].alleles) - 2.0 * p
      end
    end
  end
  (nused == 0 || denominator == 0.0) && return Matrix{Float64}(I, nindividuals, nindividuals)
  return view(centered, :, 1:nused) * transpose(view(centered, :, 1:nused)) / denominator
end

"""
    linear_mixed_model_scan(genotypes::Matrix{Float64}, phenotypes::Vector{Float64}, G::Matrix{Float64})

Simple GWAS scan using a Linear Mixed Model approximation (Gram-Schmidt or Score test).
Returns p-values for each marker.
"""
function linear_mixed_model_scan(genotypes::AbstractMatrix{<:Real}, phenotypes::AbstractVector{<:Real},
                                 G::AbstractMatrix{<:Real}; ridge::Real=1e-6)
  nindividuals, nloci = size(genotypes)
  nindividuals == length(phenotypes) || throw(DimensionMismatch("genotypes and phenotypes have different sample counts"))
  size(G) == (nindividuals, nindividuals) || throw(DimensionMismatch("G must be n_samples by n_samples"))
  nindividuals > 0 || throw(ArgumentError("at least one sample is required"))
  if nindividuals <= 2
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, ones(Float64, nloci), "linear_mixed_model_scan")
  end
  ridge > 0 || throw(ArgumentError("ridge must be positive"))
  covariance = Symmetric(Matrix{Float64}(G) + ridge * I)
  factor = cholesky(covariance)
  intercept = factor.L \ ones(Float64, nindividuals)
  response = factor.L \ Float64.(phenotypes)
  intercept_norm = dot(intercept, intercept)
  response .-= intercept .* (dot(intercept, response) / intercept_norm)
  pvalues = ones(Float64, nloci)
  degrees_of_freedom = nindividuals - 2
  for locus_idx in 1:nloci
    predictor = factor.L \ Float64.(view(genotypes, :, locus_idx))
    predictor .-= intercept .* (dot(intercept, predictor) / intercept_norm)
    ssx = dot(predictor, predictor)
    ssx > eps(Float64) || continue
    beta = dot(predictor, response) / ssx
    residual = response .- beta .* predictor
    variance = dot(residual, residual) / degrees_of_freedom
    variance > 0.0 || continue
    z = beta / sqrt(variance / ssx)
    pvalues[locus_idx] = 2.0 * ccdf(Normal(), abs(z))
  end
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, pvalues, "linear_mixed_model_scan")
end

"""
    inbreeding_coefficient(G::Matrix{Float64}, ind_idx::Int)

Calculates the inbreeding coefficient (F) from the GRM diagonal.
F = G[i, i] - 1.
"""
function inbreeding_coefficient(G::Matrix{Float64}, ind_idx::Int)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, G[ind_idx, ind_idx] - 1.0, "inbreeding_coefficient")
end

"""
    relatedness(G::Matrix{Float64}, ind1::Int, ind2::Int)

Returns the genomic relatedness coefficient between two individuals from the GRM.
"""
function relatedness(G::Matrix{Float64}, ind1::Int, ind2::Int)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, G[ind1, ind2], "relatedness")
end

# ------------------------------------------------------------------------------
# 18. Population data exploration and multi-locus summaries
# ------------------------------------------------------------------------------

"""A labeled symmetric matrix of pairwise population FST estimates."""
struct PairwiseFSTResult
  populations::Vector{String}
  estimates::Matrix{Float64}
  method::Symbol
  loci_used::Matrix{Int}
end

function Base.show(io::IO, result::PairwiseFSTResult)
  print(io, "PairwiseFSTResult(method=:$(result.method), populations=",
        join(result.populations, ", "), ")")
end

"""Return population names, or a name-to-sample-count map when `counts=true`."""
function populations(data::AbstractVector{<:Population}; counts::Bool=false)
  counts || return [pop.name for pop in data]
  return Dict(pop.name => length(pop.individuals) for pop in data)
end
populations(pop::Population; counts::Bool=false) = counts ? Dict(pop.name => length(pop.individuals)) : [pop.name]

"""Return one-based locus indices available in a population collection."""
function loci(data::AbstractVector{<:Population})
  maximum((length(ind.loci) for pop in data for ind in pop.individuals); init=0) |> x -> collect(1:x)
end
loci(record::GenePopRecord) = copy(record.loci_list)

"""Return individual identifiers in population order."""
samplenames(data::AbstractVector{<:Population}) = [ind.id for pop in data for ind in pop.individuals]
samplenames(pop::Population) = [ind.id for ind in pop.individuals]

 function _population_locus_stats(pop::Population{T}, locus_idx::Int, missing_alleles) where T
  counts = Dict{T,Int}()
  n = heterozygotes = 0
  for ind in pop.individuals
    locus_idx <= length(ind.loci) || continue
    genotype = ind.loci[locus_idx]
    _valid_diploid(genotype, missing_alleles) || continue
    a, b = genotype.alleles
    counts[a] = get(counts, a, 0) + 1
    counts[b] = get(counts, b, 0) + 1
    heterozygotes += a != b
    n += 1
  end
  n == 0 && return (n=0, ho=NaN, he=NaN, counts=counts)
  he = 1.0 - sum((count / (2.0 * n))^2 for count in values(counts))
  return (n=n, ho=heterozygotes / n, he=he, counts=counts)
end

"""
    missingdata(data; by=:sample, missing_alleles=DEFAULT_MISSING_ALLELES)

Return non-destructive missing-call summaries by `:sample`, `:population`,
`:locus`, or `:locusxpopulation`. An absent locus is counted as missing.
"""
function missingdata(data::AbstractVector{<:Population}; by::Union{Symbol,AbstractString}=:sample,
                     missing_alleles=DEFAULT_MISSING_ALLELES)
  mode = Symbol(lowercase(String(by)))
  mode in (:sample, :population, :locus, :locusxpopulation) ||
    throw(ArgumentError("by must be :sample, :population, :locus, or :locusxpopulation"))
  nloci = length(loci(data))
  incomplete(ind, locus_idx) = locus_idx > length(ind.loci) || !_valid_diploid(ind.loci[locus_idx], missing_alleles)
  rows = NamedTuple[]
  if mode === :sample
    for pop in data, ind in pop.individuals
      nmissing = count(locus_idx -> incomplete(ind, locus_idx), 1:nloci)
      push!(rows, (population=pop.name, sample=ind.id, missing=nmissing, total=nloci, percent=nloci == 0 ? 0.0 : nmissing / nloci))
    end
  elseif mode === :population
    for pop in data
      total = length(pop.individuals) * nloci
      nmissing = sum(count(locus_idx -> incomplete(ind, locus_idx), 1:nloci) for ind in pop.individuals)
      push!(rows, (population=pop.name, missing=nmissing, total=total, percent=total == 0 ? 0.0 : nmissing / total))
    end
  elseif mode === :locus
    for locus_idx in 1:nloci
      total = sum(length(pop.individuals) for pop in data)
      nmissing = sum(incomplete(ind, locus_idx) for pop in data for ind in pop.individuals)
      push!(rows, (locus=locus_idx, missing=nmissing, total=total, percent=total == 0 ? 0.0 : nmissing / total))
    end
  else
    for pop in data, locus_idx in 1:nloci
      total = length(pop.individuals)
      nmissing = count(ind -> incomplete(ind, locus_idx), pop.individuals)
      push!(rows, (population=pop.name, locus=locus_idx, missing=nmissing, total=total, percent=total == 0 ? 0.0 : nmissing / total))
    end
  end
  return rows
end

"""
    richness(data; by=:locus)

Count distinct observed alleles per locus, optionally separately per population.
"""
function richness(data::AbstractVector{<:Population}; by::Union{Symbol,AbstractString}=:locus,
                  missing_alleles=DEFAULT_MISSING_ALLELES)
  mode = Symbol(lowercase(String(by)))
  mode in (:locus, :population) || throw(ArgumentError("by must be :locus or :population"))
  rows = NamedTuple[]
  if mode === :locus
    for locus_idx in loci(data)
      alleles = Set{Any}()
      for pop in data, ind in pop.individuals
        locus_idx <= length(ind.loci) || continue
        for allele in ind.loci[locus_idx].alleles
          _is_missing_allele(allele, missing_alleles) || push!(alleles, allele)
        end
      end
      push!(rows, (locus=locus_idx, richness=length(alleles)))
    end
  else
    for pop in data, locus_idx in loci(data)
      alleles = Set{Any}()
      for ind in pop.individuals
        locus_idx <= length(ind.loci) || continue
        for allele in ind.loci[locus_idx].alleles
          _is_missing_allele(allele, missing_alleles) || push!(alleles, allele)
        end
      end
      push!(rows, (population=pop.name, locus=locus_idx, richness=length(alleles)))
    end
  end
  return rows
end

"""Return the mean and sample standard deviation of per-locus allelic richness."""
function alleleaverage(data::AbstractVector{<:Population}; missing_alleles=DEFAULT_MISSING_ALLELES)
  values = [row.richness for row in richness(data; missing_alleles=missing_alleles)]
  isempty(values) && return (mean=NaN, stdev=NaN)
  return (mean=mean(values), stdev=length(values) > 1 ? std(values) : 0.0)
end

# --- FST Estimator Helpers ---

function _nei_pairwise_fst(pop1::Population{T}, pop2::Population{T}, locus_idx::Int, missing_alleles) where T
  s1 = _population_locus_stats(pop1, locus_idx, missing_alleles)
  s2 = _population_locus_stats(pop2, locus_idx, missing_alleles)
  (s1.n == 0 || s2.n == 0) && return (NaN, 0)
  total = 2 * (s1.n + s2.n)
  pooled = Dict{T,Int}(s1.counts)
  for (allele, count) in s2.counts
    pooled[allele] = get(pooled, allele, 0) + count
  end
  ht = 1.0 - sum((count / total)^2 for count in values(pooled))
  hs = (s1.he + s2.he) / 2.0
  return (ht == 0.0 ? NaN : (ht - hs) / ht, 1)
end

function _hudson_pairwise_fst(pop1::Population{T}, pop2::Population{T}, locus_idx::Int, missing_alleles) where T
  s1 = _population_locus_stats(pop1, locus_idx, missing_alleles)
  s2 = _population_locus_stats(pop2, locus_idx, missing_alleles)
  (s1.n == 0 || s2.n == 0) && return (NaN, NaN, 0)
  alleles = union(keys(s1.counts), keys(s2.counts))
  isempty(alleles) && return (NaN, NaN, 0)
  a = first(alleles)
  p1 = get(s1.counts, a, 0) / (2.0 * s1.n)
  p2 = get(s2.counts, a, 0) / (2.0 * s2.n)
  n1, n2 = 2.0 * s1.n, 2.0 * s2.n
  num = (p1 - p2)^2 - (p1 * (1.0 - p1) / (n1 - 1.0)) - (p2 * (1.0 - p2) / (n2 - 1.0))
  den = p1 * (1.0 - p2) + p2 * (1.0 - p1)
  return (num, den, 1)
end

function _weir_cockerham_components(pop1::Population{T}, pop2::Population{T}, locus_idx::Int, missing_alleles) where T
  s1 = _population_locus_stats(pop1, locus_idx, missing_alleles)
  s2 = _population_locus_stats(pop2, locus_idx, missing_alleles)
  (s1.n == 0 || s2.n == 0) && return (0.0, 0.0, 0.0, 0)
  n1, n2 = Float64(s1.n), Float64(s2.n)
  n_bar = (n1 + n2) / 2.0
  n_c = 2.0 * n_bar - (n1^2 + n2^2) / (2.0 * n_bar)
  alleles = union(keys(s1.counts), keys(s2.counts))
  a_tot = b_tot = c_tot = 0.0
  r = 2.0 # 2 populations
  for a in alleles
    p1 = get(s1.counts, a, 0) / (2.0 * n1)
    p2 = get(s2.counts, a, 0) / (2.0 * n2)
    p_bar = (n1 * p1 + n2 * p2) / (n1 + n2)
    s2_val = (n1 * (p1 - p_bar)^2 + n2 * (p2 - p_bar)^2) / ((r - 1.0) * n_bar)
    h1 = 0.0
    for ind in pop1.individuals
      locus_idx <= length(ind.loci) || continue
      _valid_diploid(ind.loci[locus_idx], missing_alleles) || continue
      al = ind.loci[locus_idx].alleles
      h1 += (al[1] == a && al[2] != a) || (al[1] != a && al[2] == a)
    end
    h2 = 0.0
    for ind in pop2.individuals
      locus_idx <= length(ind.loci) || continue
      _valid_diploid(ind.loci[locus_idx], missing_alleles) || continue
      al = ind.loci[locus_idx].alleles
      h2 += (al[1] == a && al[2] != a) || (al[1] != a && al[2] == a)
    end
    h_bar = (h1 + h2) / (n1 + n2)
    a_comp = (n_bar / n_c) * (s2_val - (1.0 / (n_bar - 1.0)) * (p_bar * (1.0 - p_bar) - ((r - 1.0) / r) * s2_val - 0.25 * h_bar))
    b_comp = (n_bar / (n_bar - 1.0)) * (p_bar * (1.0 - p_bar) - ((r - 1.0) / r) * s2_val - ((2.0 * n_bar - 1.0) / (4.0 * n_bar)) * h_bar)
    c_comp = 0.5 * h_bar
    a_tot += a_comp; b_tot += b_comp; c_tot += c_comp
  end
  return (a_tot, b_tot, c_tot, 1)
end

"""
    hudson_fst(pop1::Population, pop2::Population, locus_idx::Int)

Calculate Hudson et al. (1992) FST estimator for a single locus.
"""
function hudson_fst(pop1::Population{T}, pop2::Population{T}, locus_idx::Int; missing_alleles=DEFAULT_MISSING_ALLELES) where T
  num, den, valid = _hudson_pairwise_fst(pop1, pop2, locus_idx, missing_alleles)
  valid == 0 && return NaN
  return den == 0.0 ? 0.0 : num / den
end

"""
    weir_cockerham_fst(pop1::Population, pop2::Population, locus_idx::Int)

Calculate Weir & Cockerham (1984) theta (FST) for a single locus.
"""
function weir_cockerham_fst(pop1::Population{T}, pop2::Population{T}, locus_idx::Int; missing_alleles=DEFAULT_MISSING_ALLELES) where T
  a, b, c, valid = _weir_cockerham_components(pop1, pop2, locus_idx, missing_alleles)
  valid == 0 && return NaN
  denom = a + b + c
  return denom == 0.0 ? 0.0 : a / denom
end

"""
    pairwise_fst(data; method=:nei, parallel=true)

Compute pairwise FST across all populations using the specified method (`:nei`, `:hudson`, `:weir_cockerham`, or `:amova`).
Parallelized using `Threads.@threads` for multi-core performance.
"""
function pairwise_fst(data::AbstractVector{<:Population}; method::Symbol=:nei,
                      missing_alleles=DEFAULT_MISSING_ALLELES, parallel::Bool=true)
  method in (:nei, :hudson, :weir_cockerham, :amova) ||
    throw(ArgumentError("Method $method not recognized. Must be :nei, :hudson, :weir_cockerham, or :amova"))
  names = populations(data)
  n = length(data)
  estimates, loci_used = fill(NaN, n, n), zeros(Int, n, n)
  for i in 1:n
    estimates[i, i] = 0.0
  end

  pairs = [(i, j) for i in 1:n for j in (i+1):n]
  all_loci = loci(data)

  if method === :nei
    if parallel
      Threads.@threads for k in 1:length(pairs)
        i, j = pairs[k]
        numerator = denominator = 0.0
        used = 0
        for locus_idx in all_loci
          fst, valid = _nei_pairwise_fst(data[i], data[j], locus_idx, missing_alleles)
          valid == 0 || (numerator += fst; denominator += 1.0; used += 1)
        end
        val = denominator == 0.0 ? NaN : numerator / denominator
        estimates[i, j] = estimates[j, i] = val
        loci_used[i, j] = loci_used[j, i] = used
      end
    else
      for (i, j) in pairs
        numerator = denominator = 0.0
        used = 0
        for locus_idx in all_loci
          fst, valid = _nei_pairwise_fst(data[i], data[j], locus_idx, missing_alleles)
          valid == 0 || (numerator += fst; denominator += 1.0; used += 1)
        end
        val = denominator == 0.0 ? NaN : numerator / denominator
        estimates[i, j] = estimates[j, i] = val
        loci_used[i, j] = loci_used[j, i] = used
      end
    end
  elseif method === :hudson
    for (i, j) in pairs
      num_sum = den_sum = 0.0
      used = 0
      for locus_idx in all_loci
        num, den, valid = _hudson_pairwise_fst(data[i], data[j], locus_idx, missing_alleles)
        if valid != 0 && !isnan(num) && !isnan(den)
          num_sum += num; den_sum += den; used += 1
        end
      end
      val = den_sum == 0.0 ? NaN : num_sum / den_sum
      estimates[i, j] = estimates[j, i] = val
      loci_used[i, j] = loci_used[j, i] = used
    end
  elseif method === :weir_cockerham
    for (i, j) in pairs
      a_sum = b_sum = c_sum = 0.0
      used = 0
      for locus_idx in all_loci
        a, b, c, valid = _weir_cockerham_components(data[i], data[j], locus_idx, missing_alleles)
        if valid != 0
          a_sum += a; b_sum += b; c_sum += c; used += 1
        end
      end
      tot = a_sum + b_sum + c_sum
      val = tot == 0.0 ? NaN : a_sum / tot
      estimates[i, j] = estimates[j, i] = val
      loci_used[i, j] = loci_used[j, i] = used
    end
  elseif method === :amova
    # Individual pairwise genetic distance matrix AMOVA
    for (i, j) in pairs
      sub_pops = [data[i], data[j]]
      pop_sizes = [length(data[i].individuals), length(data[j].individuals)]
      N = sum(pop_sizes)
      inds = [data[i].individuals; data[j].individuals]
      dist_mat = zeros(Float64, N, N)
      for x in 1:N, y in (x+1):N
        d = 0.0
        valid_l = 0
        for l_idx in all_loci
          g1 = l_idx <= length(inds[x].loci) ? inds[x].loci[l_idx] : nothing
          g2 = l_idx <= length(inds[y].loci) ? inds[y].loci[l_idx] : nothing
          if g1 !== nothing && g2 !== nothing && _valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)
            d += (g1.alleles[1] != g2.alleles[1]) + (g1.alleles[2] != g2.alleles[2])
            valid_l += 1
          end
        end
        dist_mat[x, y] = dist_mat[y, x] = valid_l > 0 ? d / valid_l : 0.0
      end
      phi_st, _, _ = amova(dist_mat, pop_sizes)
      estimates[i, j] = estimates[j, i] = phi_st
      loci_used[i, j] = loci_used[j, i] = length(all_loci)
    end
  end

  _ctx = active_provenance_context()
  res = PairwiseFSTResult(names, estimates, method, loci_used)
  return _register_popgen_result!(_ctx, res, "pairwise_fst"; parameters=(method=method, n_populations=n))
end

"""
    fst_permutation_test(data::Vector{Population}; method=:nei, permutations=999)

Perform permutation testing on pairwise FST by randomly shuffling individual population labels.
Returns a NamedTuple `(observed=PairwiseFSTResult, pvalues=Matrix{Float64})`.
"""
function fst_permutation_test(data::AbstractVector{<:Population}; method::Symbol=:nei,
                              permutations::Int=999, rng::AbstractRNG=Random.default_rng())
  obs = pairwise_fst(data; method=method, parallel=false)
  n_pops = length(data)
  pop_sizes = [length(p.individuals) for p in data]
  all_inds = [ind for p in data for ind in p.individuals]

  greater_counts = zeros(Int, n_pops, n_pops)

  for _ in 1:permutations
    shuffled = Random.shuffle(rng, all_inds)
    perm_pops = Population[]
    start_i = 1
    for (idx, psize) in enumerate(pop_sizes)
      pinds = shuffled[start_i:(start_i + psize - 1)]
      push!(perm_pops, Population(data[idx].name, pinds))
      start_i += psize
    end
    perm_fst = pairwise_fst(perm_pops; method=method, parallel=false)
    for i in 1:n_pops, j in (i+1):n_pops
      if perm_fst.estimates[i, j] >= obs.estimates[i, j]
        greater_counts[i, j] += 1
        greater_counts[j, i] += 1
      end
    end
  end

  pvals = (greater_counts .+ 1.0) ./ (permutations + 1.0)
  for i in 1:n_pops
    pvals[i, i] = 0.0
  end

  _ctx = active_provenance_context()
  res = (observed=obs, pvalues=pvals)
  return _register_popgen_result!(_ctx, res, "fst_permutation_test"; parameters=(method=method, permutations=permutations))
end

"""
    summary_statistics(data)

Return pooled multi-locus observed/expected heterozygosity and Nei FST.
The FST is the ratio of summed among-population diversity to summed total diversity.
"""
function summary_statistics(data::AbstractVector{<:Population}; missing_alleles=DEFAULT_MISSING_ALLELES)
  isempty(data) && return (n_populations=0, n_loci=0, observed_heterozygosity=NaN, expected_heterozygosity=NaN, fst=NaN)
  ho_sum = he_sum = weight = fst_num = fst_den = 0.0
  used_loci = 0
  for locus_idx in loci(data)
    stats = [_population_locus_stats(pop, locus_idx, missing_alleles) for pop in data]
    valid = filter(s -> s.n > 0, stats)
    isempty(valid) && continue
    local_weight = sum(s.n for s in valid)
    ho_sum += sum(s.ho * s.n for s in valid)
    he_sum += sum(s.he * s.n for s in valid)
    weight += local_weight
    pooled = Dict{Any,Int}()
    for s in valid, (allele, count) in s.counts
      pooled[allele] = get(pooled, allele, 0) + count
    end
    ht = 1.0 - sum((count / (2.0 * local_weight))^2 for count in values(pooled))
    hs = sum(s.he * s.n for s in valid) / local_weight
    ht > 0.0 && (fst_num += ht - hs; fst_den += ht)
    used_loci += 1
  end
  return (n_populations=length(data), n_loci=used_loci,
          observed_heterozygosity=weight == 0.0 ? NaN : ho_sum / weight,
          expected_heterozygosity=weight == 0.0 ? NaN : he_sum / weight,
          fst=fst_den == 0.0 ? NaN : fst_num / fst_den)
end

Base.summary(data::AbstractVector{<:Population}) = summary_statistics(data)

# --- Per-sample Heterozygosity ---

"""
    sample_heterozygosity(data::Vector{Population})

Calculate per-individual observed heterozygosity across all loci.
Returns a vector of NamedTuples: `(population=String, sample=String, ho=Float64, total_loci=Int)`.
"""
function sample_heterozygosity(data::AbstractVector{<:Population}; missing_alleles=DEFAULT_MISSING_ALLELES)
  results = NamedTuple{(:population, :sample, :ho, :total_loci), Tuple{String, String, Float64, Int}}[]
  for pop in data
    for ind in pop.individuals
      n_het = 0
      valid = 0
      for g in ind.loci
        if _valid_diploid(g, missing_alleles)
          valid += 1
          n_het += (g.alleles[1] != g.alleles[2])
        end
      end
      ho = valid > 0 ? n_het / valid : NaN
      push!(results, (population=pop.name, sample=ind.id, ho=ho, total_loci=valid))
    end
  end
  return results
end

# --- Kinship Moment Estimators ---

function _get_allele_freq_dict(data::AbstractVector{<:Population}, missing_alleles)
  freq_map = Dict{Int, Dict{Any, Float64}}()
  all_l = loci(data)
  for l in all_l
    counts = Dict{Any, Int}()
    total = 0
    for pop in data
      for ind in pop.individuals
        l <= length(ind.loci) || continue
        for a in ind.loci[l].alleles
          _is_missing_allele(a, missing_alleles) && continue
          counts[a] = get(counts, a, 0) + 1
          total += 1
        end
      end
    end
    if total > 0
      freq_map[l] = Dict{Any, Float64}(a => c / Float64(total) for (a, c) in counts)
    end
  end
  return freq_map
end

"""
    kinship_queller_goodnight(ind1, ind2, freq_map)
Queller & Goodnight (1989) pairwise relatedness estimator.
"""
function kinship_queller_goodnight(ind1::PopGenIndividual, ind2::PopGenIndividual, freq_map; missing_alleles=DEFAULT_MISSING_ALLELES)
  num1 = num2 = den1 = den2 = 0.0
  n_loci = min(length(ind1.loci), length(ind2.loci))
  for l in 1:n_loci
    g1, g2 = ind1.loci[l], ind2.loci[l]
    (_valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)) || continue
    f_dict = get(freq_map, l, nothing)
    f_dict === nothing && continue
    a, b = g1.alleles[1], g1.alleles[2]
    c, d = g2.alleles[1], g2.alleles[2]
    fa, fb = get(f_dict, a, 0.0), get(f_dict, b, 0.0)
    fc, fd = get(f_dict, c, 0.0), get(f_dict, d, 0.0)
    ident = Float64((a == c) + (a == d) + (b == c) + (b == d))
    num1 += ident - 2.0 * (fa + fb)
    num2 += ident - 2.0 * (fc + fd)
    den1 += 2.0 * (1.0 + (a == b) - fa - fb)
    den2 += 2.0 * (1.0 + (c == d) - fc - fd)
  end
  r1 = den1 == 0.0 ? NaN : num1 / den1
  r2 = den2 == 0.0 ? NaN : num2 / den2
  return (r1 + r2) / 2.0
end

"""
    kinship_blouin(ind1, ind2)
Blouin et al. (1996) identity-by-state similarity score.
"""
function kinship_blouin(ind1::PopGenIndividual, ind2::PopGenIndividual; missing_alleles=DEFAULT_MISSING_ALLELES)
  res = 0.0
  valid = 0
  n_loci = min(length(ind1.loci), length(ind2.loci))
  for l in 1:n_loci
    g1, g2 = ind1.loci[l], ind2.loci[l]
    (_valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)) || continue
    a, b = g1.alleles[1], g1.alleles[2]
    c, d = g2.alleles[1], g2.alleles[2]
    res += Float64((a == c || a == d) + (b == c || b == d))
    valid += 1
  end
  return valid == 0 ? NaN : res / (2.0 * valid)
end

"""
    kinship_li_horvitz(ind1, ind2)
Li & Horvitz (1953) allele-sharing score.
"""
function kinship_li_horvitz(ind1::PopGenIndividual, ind2::PopGenIndividual; missing_alleles=DEFAULT_MISSING_ALLELES)
  res = 0.0
  valid = 0
  n_loci = min(length(ind1.loci), length(ind2.loci))
  for l in 1:n_loci
    g1, g2 = ind1.loci[l], ind2.loci[l]
    (_valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)) || continue
    a, b = g1.alleles[1], g1.alleles[2]
    c, d = g2.alleles[1], g2.alleles[2]
    res += Float64((a == c) + (a == d) + (b == c) + (b == d))
    valid += 1
  end
  return valid == 0 ? NaN : res / (4.0 * valid)
end

"""
    kinship_ritland(ind1, ind2, freq_map)
Ritland (1996) pairwise relatedness estimator.
"""
function kinship_ritland(ind1::PopGenIndividual, ind2::PopGenIndividual, freq_map; missing_alleles=DEFAULT_MISSING_ALLELES)
  numer = denom = 0.0
  n_loci = min(length(ind1.loci), length(ind2.loci))
  for l in 1:n_loci
    g1, g2 = ind1.loci[l], ind2.loci[l]
    (_valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)) || continue
    f_dict = get(freq_map, l, nothing)
    f_dict === nothing && continue
    K = length(f_dict)
    K <= 1 && continue
    A = K - 1.0
    a, b = g1.alleles[1], g1.alleles[2]
    c, d = g2.alleles[1], g2.alleles[2]
    R = 0.0
    for (allele, frq) in f_dict
      frq > 0.0 || continue
      s1 = Float64((a == allele) + (b == allele))
      s2 = Float64((c == allele) + (d == allele))
      R += (s1 * s2) / (4.0 * frq)
    end
    r_loc = (2.0 / A) * (R - 1.0)
    numer += r_loc * A
    denom += A
  end
  return denom == 0.0 ? NaN : numer / denom
end

"""
    kinship_lynch_li(ind1, ind2, freq_map)
Lynch & Li (1993) relatedness estimator.
"""
function kinship_lynch_li(ind1::PopGenIndividual, ind2::PopGenIndividual, freq_map; missing_alleles=DEFAULT_MISSING_ALLELES)
  num = den = 0.0
  n_loci = min(length(ind1.loci), length(ind2.loci))
  for l in 1:n_loci
    g1, g2 = ind1.loci[l], ind2.loci[l]
    (_valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)) || continue
    f_dict = get(freq_map, l, nothing)
    f_dict === nothing && continue
    a, b = g1.alleles[1], g1.alleles[2]
    c, d = g2.alleles[1], g2.alleles[2]
    s0 = sum(p^2 for p in values(f_dict)) - 0.5 * sum(p^3 for p in values(f_dict))
    s_xy = 0.5 * (((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (a == b))) +
                  ((a == c) + (a == d) + (b == c) + (b == d)) / (2.0 * (1.0 + (c == d))))
    num += s_xy - s0
    den += 1.0 - s0
  end
  return den == 0.0 ? NaN : num / den
end

"""
    kinship_lynch_ritland(ind1, ind2, freq_map)
Lynch & Ritland (1999) regression estimator of relatedness.
"""
function kinship_lynch_ritland(ind1::PopGenIndividual, ind2::PopGenIndividual, freq_map; missing_alleles=DEFAULT_MISSING_ALLELES)
  num_sum = den_sum = 0.0
  n_loci = min(length(ind1.loci), length(ind2.loci))
  for l in 1:n_loci
    g1, g2 = ind1.loci[l], ind2.loci[l]
    (_valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)) || continue
    f_dict = get(freq_map, l, nothing)
    f_dict === nothing && continue
    a, b = g1.alleles[1], g1.alleles[2]
    c, d = g2.alleles[1], g2.alleles[2]
    fa, fb = get(f_dict, a, 0.0), get(f_dict, b, 0.0)
    fc, fd = get(f_dict, c, 0.0), get(f_dict, d, 0.0)
    (fa > 0 && fb > 0 && fc > 0 && fd > 0) || continue
    n1 = fa * Float64((b == c) + (b == d)) + fb * Float64((a == c) + (a == d)) - 4.0 * fa * fb
    d1 = 2.0 * (1.0 + (a == b)) * (fa + fb) - 8.0 * fa * fb
    w1 = ((1.0 + (a == b)) * (fa + fb) - 4.0 * fa * fb) / (2.0 * fa * fb)
    d1 > 0.0 && w1 > 0.0 || continue
    num_sum += (n1 / d1) * w1
    den_sum += w1
  end
  return den_sum == 0.0 ? NaN : num_sum / den_sum
end

"""
    kinship_moran(ind1, ind2, freq_map)
Moran-style genetic covariance kinship coefficient.
"""
function kinship_moran(ind1::PopGenIndividual, ind2::PopGenIndividual, freq_map; missing_alleles=DEFAULT_MISSING_ALLELES)
  num = den = 0.0
  n_loci = min(length(ind1.loci), length(ind2.loci))
  for l in 1:n_loci
    g1, g2 = ind1.loci[l], ind2.loci[l]
    (_valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)) || continue
    f_dict = get(freq_map, l, nothing)
    f_dict === nothing && continue
    a, b = g1.alleles[1], g1.alleles[2]
    c, d = g2.alleles[1], g2.alleles[2]
    for (allele, fq) in f_dict
      g1_val = (Float64(a == allele) + Float64(b == allele)) / 2.0 - fq
      g2_val = (Float64(c == allele) + Float64(d == allele)) / 2.0 - fq
      num += g1_val * g2_val
      den += (g1_val^2 + g2_val^2)
    end
  end
  return den == 0.0 ? NaN : num / (den / 2.0)
end

"""
    kinship_loiselle(ind1, ind2, freq_map)
Loiselle et al. (1987) kinship estimator.
"""
function kinship_loiselle(ind1::PopGenIndividual, ind2::PopGenIndividual, freq_map; missing_alleles=DEFAULT_MISSING_ALLELES)
  num = den = 0.0
  n_loci = min(length(ind1.loci), length(ind2.loci))
  for l in 1:n_loci
    g1, g2 = ind1.loci[l], ind2.loci[l]
    (_valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)) || continue
    f_dict = get(freq_map, l, nothing)
    f_dict === nothing && continue
    a, b = g1.alleles[1], g1.alleles[2]
    c, d = g2.alleles[1], g2.alleles[2]
    for (allele, fq) in f_dict
      p1 = (Float64(a == allele) + Float64(b == allele)) / 2.0 - fq
      p2 = (Float64(c == allele) + Float64(d == allele)) / 2.0 - fq
      num += p1 * p2
      den += fq * (1.0 - fq)
    end
  end
  return den == 0.0 ? NaN : num / den
end

"""
    pairwise_kinship(data::Vector{Population}; method=:queller_goodnight)

Calculate pairwise kinship/relatedness between all individuals across populations.
Returns `(sample_names=Vector{String}, matrix=Matrix{Float64}, method=Symbol)`.
"""
function pairwise_kinship(data::AbstractVector{<:Population}; method::Symbol=:queller_goodnight,
                          missing_alleles=DEFAULT_MISSING_ALLELES)
  all_inds = [ind for pop in data for ind in pop.individuals]
  names = [ind.id for ind in all_inds]
  N = length(all_inds)
  matrix = fill(NaN, N, N)
  freq_map = _get_allele_freq_dict(data, missing_alleles)

  for i in 1:N
    matrix[i, i] = 1.0
    for j in (i+1):N
      ind1, ind2 = all_inds[i], all_inds[j]
      val = if method === :queller_goodnight
        kinship_queller_goodnight(ind1, ind2, freq_map; missing_alleles=missing_alleles)
      elseif method === :blouin
        kinship_blouin(ind1, ind2; missing_alleles=missing_alleles)
      elseif method === :li_horvitz
        kinship_li_horvitz(ind1, ind2; missing_alleles=missing_alleles)
      elseif method === :ritland
        kinship_ritland(ind1, ind2, freq_map; missing_alleles=missing_alleles)
      elseif method === :lynch_li
        kinship_lynch_li(ind1, ind2, freq_map; missing_alleles=missing_alleles)
      elseif method === :lynch_ritland
        kinship_lynch_ritland(ind1, ind2, freq_map; missing_alleles=missing_alleles)
      elseif method === :moran
        kinship_moran(ind1, ind2, freq_map; missing_alleles=missing_alleles)
      elseif method === :loiselle
        kinship_loiselle(ind1, ind2, freq_map; missing_alleles=missing_alleles)
      else
        throw(ArgumentError("Unknown kinship method: $method"))
      end
      matrix[i, j] = matrix[j, i] = val
    end
  end

  _ctx = active_provenance_context()
  res = (sample_names=names, matrix=matrix, method=method)
  return _register_popgen_result!(_ctx, res, "pairwise_kinship"; parameters=(method=method, n_samples=N))
end

# --- Population Clustering Wrappers ---

"""
    population_kmeans(data::Vector{Population}; k::Int)

Perform K-means clustering on population allele frequency profiles.
"""
function population_kmeans(data::AbstractVector{<:Population}; k::Int=2, missing_alleles=DEFAULT_MISSING_ALLELES)
  pops = populations(data)
  n_pops = length(data)
  all_l = loci(data)
  
  locus_alleles = Dict{Int, Vector{Any}}()
  for l in all_l
    al_set = Set{Any}()
    for pop in data
      af = allele_frequencies(pop, l; missing_alleles=missing_alleles)
      for a in keys(af)
        push!(al_set, a)
      end
    end
    locus_alleles[l] = collect(al_set)
  end

  features = Tuple{Int, Any}[]
  for l in all_l
    for a in locus_alleles[l]
      push!(features, (l, a))
    end
  end
  
  dim = length(features)
  dim == 0 && return (populations=pops, assignments=ones(Int, n_pops), centers=zeros(Float64, k, 0))

  mat = zeros(Float64, n_pops, dim)
  for (p_i, pop) in enumerate(data)
    for (f_i, (l, a)) in enumerate(features)
      af = allele_frequencies(pop, l; missing_alleles=missing_alleles)
      mat[p_i, f_i] = get(af, a, 0.0)
    end
  end

  k_actual = min(k, n_pops)
  centers = copy(mat[1:k_actual, :])
  assignments = zeros(Int, n_pops)

  for iter in 1:50
    for i in 1:n_pops
      dists = [sum((mat[i, :] .- centers[c, :]).^2) for c in 1:k_actual]
      assignments[i] = argmin(dists)
    end
    for c in 1:k_actual
      members = findall(==(c), assignments)
      if !isempty(members)
        centers[c, :] .= mean(mat[members, :], dims=1)[:]
      end
    end
  end

  return (populations=pops, assignments=assignments, centers=centers)
end

"""
    population_hclust(data::Vector{Population})

Perform hierarchical clustering (UPGMA) on population genetic distances.
Returns linkage matrix and ordered population names.
"""
function population_hclust(data::AbstractVector{<:Population}; missing_alleles=DEFAULT_MISSING_ALLELES)
  fst_res = pairwise_fst(data; method=:nei, missing_alleles=missing_alleles, parallel=false)
  dist_mat = copy(fst_res.estimates)
  n = length(fst_res.populations)

  clusters = [[i] for i in 1:n]
  heights = Float64[]
  linkage = Tuple{Int, Int, Float64}[]

  active = fill(true, n)
  curr_dist = copy(dist_mat)

  for step in 1:(n-1)
    min_d = Inf
    best_pair = (0, 0)
    for i in 1:n
      active[i] || continue
      for j in (i+1):n
        active[j] || continue
        if curr_dist[i, j] < min_d
          min_d = curr_dist[i, j]
          best_pair = (i, j)
        end
      end
    end
    u, v = best_pair
    push!(linkage, (u, v, min_d))
    active[v] = false
    for k in 1:n
      if active[k] && k != u
        curr_dist[u, k] = curr_dist[k, u] = (curr_dist[u, k] + curr_dist[v, k]) / 2.0
      end
    end
  end

  return (populations=fst_res.populations, linkage=linkage)
end

"""
    population_kmedoids(data::Vector{Population}; k::Int=2)

Perform K-medoids (PAM) clustering on population pairwise genetic distances.
"""
function population_kmedoids(data::AbstractVector{<:Population}; k::Int=2, missing_alleles=DEFAULT_MISSING_ALLELES)
  fst_res = pairwise_fst(data; method=:nei, missing_alleles=missing_alleles, parallel=false)
  dist_mat = fst_res.estimates
  pops = fst_res.populations
  n = length(pops)
  k_actual = min(k, n)

  # PAM initialization: pick initial medoids evenly
  medoids = collect(1:k_actual)
  assignments = zeros(Int, n)

  for iter in 1:50
    # Assign points to nearest medoid
    for i in 1:n
      dists = [dist_mat[i, m] for m in medoids]
      assignments[i] = medoids[argmin(dists)]
    end
    # Update medoids
    for c_idx in 1:k_actual
      m_curr = medoids[c_idx]
      members = findall(==(m_curr), assignments)
      isempty(members) && continue
      best_m = m_curr
      min_cost = sum(dist_mat[i, m_curr] for i in members)
      for candidate in members
        cost = sum(dist_mat[i, candidate] for i in members)
        if cost < min_cost
          min_cost = cost
          best_m = candidate
        end
      end
      medoids[c_idx] = best_m
    end
  end

  medoid_names = [pops[m] for m in medoids]
  return (populations=pops, medoids=medoid_names, assignments=assignments)
end

"""
    population_fuzzycmeans(data::Vector{Population}; k::Int=2, m::Float64=2.0)

Perform Fuzzy C-means clustering on population allele frequency profiles.
Returns population names and membership matrix `U` of size `(k, N)`.
"""
function population_fuzzycmeans(data::AbstractVector{<:Population}; k::Int=2, m::Float64=2.0, missing_alleles=DEFAULT_MISSING_ALLELES)
  pops = populations(data)
  n = length(data)
  all_l = loci(data)

  locus_alleles = Dict{Int, Vector{Any}}()
  for l in all_l
    al_set = Set{Any}()
    for pop in data
      af = allele_frequencies(pop, l; missing_alleles=missing_alleles)
      for a in keys(af)
        push!(al_set, a)
      end
    end
    locus_alleles[l] = collect(al_set)
  end

  features = Tuple{Int, Any}[]
  for l in all_l
    for a in locus_alleles[l]
      push!(features, (l, a))
    end
  end
  dim = length(features)
  dim == 0 && return (populations=pops, membership=ones(Float64, k, n) ./ k)

  X = zeros(Float64, n, dim)
  for (p_i, pop) in enumerate(data)
    for (f_i, (l, a)) in enumerate(features)
      af = allele_frequencies(pop, l; missing_alleles=missing_alleles)
      X[p_i, f_i] = get(af, a, 0.0)
    end
  end

  k_actual = min(k, n)
  # Initialize membership randomly & normalize columns
  rng = Random.default_rng()
  U = Random.rand(rng, k_actual, n)
  U ./= sum(U, dims=1)

  centers = zeros(Float64, k_actual, dim)
  for iter in 1:50
    # Update cluster centers
    for c in 1:k_actual
      um = U[c, :].^m
      centers[c, :] .= sum(X .* um, dims=1)[:] ./ (sum(um) + 1e-12)
    end
    # Update membership U
    for i in 1:n
      dists = [sqrt(sum((X[i, :] .- centers[c, :]).^2)) + 1e-12 for c in 1:k_actual]
      for c in 1:k_actual
        U[c, i] = 1.0 / sum((dists[c] ./ dists).^(2.0 / (m - 1.0)))
      end
    end
  end

  return (populations=pops, membership=U)
end

"""
    population_dbscan(data::Vector{Population}; eps::Float64=0.5, min_neighbors::Int=1)

Perform DBSCAN density-based clustering on population pairwise genetic distances.
"""
function population_dbscan(data::AbstractVector{<:Population}; eps::Float64=0.5, min_neighbors::Int=1, missing_alleles=DEFAULT_MISSING_ALLELES)
  fst_res = pairwise_fst(data; method=:nei, missing_alleles=missing_alleles, parallel=false)
  dist_mat = fst_res.estimates
  pops = fst_res.populations
  n = length(pops)

  labels = zeros(Int, n) # 0 = unvisited/noise
  cluster_id = 0

  for i in 1:n
    labels[i] == 0 || continue
    neighbors = findall(k -> dist_mat[i, k] <= eps, 1:n)
    if length(neighbors) < min_neighbors
      labels[i] = -1 # Noise
    else
      cluster_id += 1
      labels[i] = cluster_id
      queue = copy(neighbors)
      idx = 1
      while idx <= length(queue)
        q = queue[idx]
        if labels[q] == -1
          labels[q] = cluster_id
        elseif labels[q] == 0
          labels[q] = cluster_id
          q_neighbors = findall(k -> dist_mat[q, k] <= eps, 1:n)
          if length(q_neighbors) >= min_neighbors
            append!(queue, filter(x -> !(x in queue), q_neighbors))
          end
        end
        idx += 1
      end
    end
  end

  return (populations=pops, assignments=labels)
end

"""
    population_cluster(data::Vector{Population}; method::Symbol=:kmeans, kwargs...)

Unified clustering function interfacing `:kmeans`, `:kmedoids`, `:hclust`, `:fuzzycmeans`, and `:dbscan`.
"""
function population_cluster(data::AbstractVector{<:Population}; method::Symbol=:kmeans, kwargs...)
  if method === :kmeans
    return population_kmeans(data; kwargs...)
  elseif method === :kmedoids
    return population_kmedoids(data; kwargs...)
  elseif method === :hclust
    return population_hclust(data; kwargs...)
  elseif method === :fuzzycmeans
    return population_fuzzycmeans(data; kwargs...)
  elseif method === :dbscan
    return population_dbscan(data; kwargs...)
  else
    throw(ArgumentError("Unknown clustering method: $method. Options are :kmeans, :kmedoids, :hclust, :fuzzycmeans, :dbscan"))
  end
end

"""
    pairwise_identical(data::Vector{Population})

Calculate pairwise identical genotype matching fraction between all individuals.
Returns `(sample_names=Vector{String}, matrix=Matrix{Float64})`.
"""
function pairwise_identical(data::AbstractVector{<:Population}; missing_alleles=DEFAULT_MISSING_ALLELES)
  all_inds = [ind for pop in data for ind in pop.individuals]
  names = [ind.id for ind in all_inds]
  N = length(all_inds)
  matrix = zeros(Float64, N, N)

  all_l = loci(data)

  for i in 1:N
    matrix[i, i] = 1.0
    for j in (i+1):N
      ind1, ind2 = all_inds[i], all_inds[j]
      same_genotypes = 0
      valid_loci = 0
      n_loci = min(length(ind1.loci), length(ind2.loci))
      for l in 1:n_loci
        g1, g2 = ind1.loci[l], ind2.loci[l]
        (_valid_diploid(g1, missing_alleles) && _valid_diploid(g2, missing_alleles)) || continue
        valid_loci += 1
        a1, b1 = g1.alleles[1], g1.alleles[2]
        a2, b2 = g2.alleles[1], g2.alleles[2]
        # Match unordered genotypes
        if (a1 == a2 && b1 == b2) || (a1 == b2 && b1 == a2)
          same_genotypes += 1
        end
      end
      val = valid_loci > 0 ? same_genotypes / Float64(valid_loci) : NaN
      matrix[i, j] = matrix[j, i] = val
    end
  end

  return (sample_names=names, matrix=matrix)
end

# --- STRUCTURE Format I/O ---

"""
    read_structure_record(filepath::String; ploidy::Int=2, missing_code="-9")

Read a STRUCTURE format genotype file (.str or .structure) into a Vector of `Population`.
"""
function read_structure_record(filepath::String; ploidy::Int=2, missing_code::String="-9")
  lines = readlines(filepath)
  filter!(l -> !isempty(strip(l)) && !startswith(strip(l), "#"), lines)
  isempty(lines) && return Population[]

  # Check if first line is locus header
  tokens1 = Base.split(strip(lines[1]))
  has_header = false
  loc_names = String[]
  start_line = 1

  if length(tokens1) > 0 && tryparse(Int, tokens1[1]) === nothing && tryparse(Int, tokens1[end]) === nothing
    loc_names = [String(t) for t in tokens1]
    start_line = 2
  end

  pop_dict = Dict{String, Vector{PopGenIndividual{Int}}}()

  for i in start_line:length(lines)
    line = strip(lines[i])
    isempty(line) && continue
    toks = Base.split(line)
    length(toks) >= 2 || continue
    sample_id = String(toks[1])
    pop_id = String(toks[2])
    raw_alleles = toks[3:end]
    loci_vec = Locus{Int}[]

    for chunk in Iterators.partition(raw_alleles, ploidy)
      length(chunk) == ploidy || continue
      al1 = tryparse(Int, chunk[1])
      al2 = tryparse(Int, chunk[2])
      a1 = (al1 === nothing || String(chunk[1]) == missing_code) ? 0 : al1
      a2 = (al2 === nothing || String(chunk[2]) == missing_code) ? 0 : al2
      push!(loci_vec, Locus{Int}((a1, a2)))
    end

    ind = PopGenIndividual{Int}(sample_id, loci_vec)
    pop_list = get!(pop_dict, pop_id, PopGenIndividual{Int}[])
    push!(pop_list, ind)
  end

  return [Population{Int}(pop_name, inds) for (pop_name, inds) in pop_dict]
end

"""
    write_structure(filepath::String, data::Vector{Population}; missing_code="-9")

Write population dataset to a STRUCTURE format genotype file.
"""
function write_structure(filepath::String, data::AbstractVector{<:Population}; missing_code::String="-9")
  open(filepath, "w") do io
    all_l = loci(data)
    # Header line
    println(io, join(["Locus_$l" for l in all_l], "\t"))
    for (pop_idx, pop) in enumerate(data)
      for ind in pop.individuals
        toks = String[ind.id, string(pop.name)]
        for l in all_l
          if l <= length(ind.loci)
            al = ind.loci[l].alleles
            push!(toks, string(al[1] == 0 ? missing_code : al[1]))
            push!(toks, string(al[2] == 0 ? missing_code : al[2]))
          else
            push!(toks, missing_code)
            push!(toks, missing_code)
          end
        end
        println(io, join(toks, "\t"))
      end
    end
  end
  return filepath
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
    i2 = rand(1:(k-1))
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
  haplotypes = Dict{String,Vector{Int}}()
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
  hap_strings = Dict{String,String}()
  for (name, hap) in haplotypes
    hap_strings[name] = join(string.(hap))
  end

  _ctx = active_provenance_context()
  return provenance_result!(_ctx, (tree, hap_strings), "simulate_coalescent")
end


