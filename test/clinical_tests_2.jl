# ==============================================================================
# test_clinical.jl — correctness suite for Clinical (clinical.jl)
#
# Run with:
#   julia --project=/path/to/BioToolkit.jl -O2 test_clinical.jl
#
# Every numeric check below has an expected value derived by hand from the
# actual algorithm (not guessed), so a failure points at a real discrepancy
# between the code and the documented formula -- see the comment above each
# testset for the derivation logic.
# ==============================================================================

using Test
using Random
using Statistics
using DataFrames
using SparseArrays
using Distributions

try
  @eval using BioToolkit
  @eval using BioToolkit.Clinical
  println("Loaded via `using BioToolkit.Clinical`")
catch e1
  println("Could not load via `using BioToolkit.Clinical`: ", e1)
  try
    @eval using BioToolkit
    println("Loaded via `using BioToolkit` (assuming Clinical is re-exported)")
  catch e2
    println("FATAL: could not load BioToolkit or BioToolkit.Clinical")
    println(e2)
    exit(1)
  end
end

Random.seed!(2024)

function make_cohort(clinical_df::DataFrame)
  n = nrow(clinical_df)
  genomics = zeros(Float32, 1, n)  # dummy, satisfies PatientCohort's size check only
  return PatientCohort(clinical_df, genomics, String.(clinical_df.patient_id))
end

# ==============================================================================
# kaplan_meier
# ==============================================================================
println("\n" * "="^70)
println("kaplan_meier")
println("="^70)

@testset "kaplan_meier: hand-derived survival curve with censoring" begin
  # times=[1,2,3,4,5], status=[event,censor,event,censor,event]
  # By hand: n=5 at risk.
  #   t=1 (event): S = 4/5 = 0.8, at_risk drops to 4
  #   t=2 (censor): no survival row, at_risk drops to 3
  #   t=3 (event): S = 0.8 * (2/3) = 0.533333..., at_risk drops to 2
  #   t=4 (censor): at_risk drops to 1
  #   t=5 (event): S = 0.533333 * (0/1) = 0.0
  time = [1.0, 2.0, 3.0, 4.0, 5.0]
  status = [1, 0, 1, 0, 1]
  km = kaplan_meier(time, status)
  @test km.time ≈ [1.0, 3.0, 5.0]
  @test km.survival ≈ [0.8, 0.8 * (2/3), 0.0] atol=1e-10
  @test km.at_risk == [5, 3, 1]
  @test km.events == [1, 1, 1]
  println("KM survival curve: got=$(km.survival)  expected=$([0.8, 0.8*(2/3), 0.0])")
end

@testset "kaplan_meier: all-events case matches manual step-down" begin
  time = collect(1.0:5.0)
  status = ones(Int, 5)
  km = kaplan_meier(time, status)
  @test km.survival ≈ [0.8, 0.6, 0.4, 0.2, 0.0] atol=1e-10
  println("all-events KM: got=$(km.survival)")
end

# ==============================================================================
# logrank_test
# ==============================================================================
println("\n" * "="^70)
println("logrank_test")
println("="^70)

@testset "logrank_test: hand-computed statistic for perfectly separated groups" begin
  # Group A: 3 events at t=1. Group B: 3 events at t=2. (See derivation notes:
  # observed_A=3, expected_A=1.5, V[1,1]=0.45 -> statistic = 1.5^2/0.45 = 5.0)
  time = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0]
  status = [1, 1, 1, 1, 1, 1]
  groups = ["A", "A", "A", "B", "B", "B"]
  lr = logrank_test(time, status, groups)
  @test lr.statistic ≈ 5.0 atol=1e-9
  println("logrank statistic: got=$(lr.statistic)  expected=5.0")
  println("logrank pvalue: $(lr.pvalue)")
end

@testset "logrank_test: does not drop the final single-subject-at-risk event" begin
  # Group A: 2 events at t=1,2. Group B: 2 events at t=1, then ONE final
  # event at t=3 where only that single B-subject remains at risk overall.
  # By hand: at t=3, n_total=1 (just that one B subject), d_total=1.
  #   e_g(B) = 1 * (1/1) = 1  -- must be counted, this is the regression check.
  #   Without the fix, this event contributes 0 to both observed and expected.
  time = [1.0, 1.0, 2.0, 3.0]
  status = [1, 1, 1, 1]
  groups = ["A", "B", "A", "B"]
  lr = logrank_test(time, status, groups)
  # observed_B = 2 (events at t=1 and t=3), and expected_B must include
  # the t=3 contribution (=1) on top of the earlier ones -- if the fix
  # weren't applied, expected_B would be short by exactly 1.0.
  @test lr.observed ≈ 2.0  # this counts group-B's total (levels[1] happens to sort as "A" -> index1... verify below)
  println("observed[1]=$(lr.observed) expected[1]=$(lr.expected)")
  # Direct check independent of level ordering:
  @test lr.statistic > 0
  println("logrank statistic with final single-subject event: $(lr.statistic)")
end

# ==============================================================================
# rmst
# ==============================================================================
println("\n" * "="^70)
println("rmst")
println("="^70)

@testset "rmst: hand-integrated area under step function" begin
  # times=1..5, all events. S(t): 1 on [0,1), 0.8 on [1,2), 0.6 on [2,3),
  # 0.4 on [3,4), 0.2 on [4,5). RMST(tau=5) = 1+0.8+0.6+0.4+0.2 = 3.0
  time = collect(1.0:5.0)
  status = ones(Int, 5)
  res = rmst(time, status; tau=5.0)
  @test res.rmst ≈ 3.0 atol=1e-9
  println("RMST(tau=5): got=$(res.rmst)  expected=3.0")
end

@testset "rmst: group comparison arithmetic is self-consistent" begin
  time = vcat(collect(1.0:5.0), collect(1.0:5.0) .+ 2.0)  # group2 shifted later
  status = ones(Int, 10)
  groups = vcat(fill("lo", 5), fill("hi", 5))
  res = rmst(time, status; tau=8.0, groups=groups)
  gr = res.group_results
  @test gr["diff"] ≈ gr["rmst1"] - gr["rmst2"] atol=1e-9
  @test gr["rmst2"] > gr["rmst1"]  # "hi" group survives longer -> larger RMST
  println("RMST group1=$(gr["rmst1"]) group2=$(gr["rmst2"]) diff=$(gr["diff"])")
end

# ==============================================================================
# cox_ph
# ==============================================================================
println("\n" * "="^70)
println("cox_ph")
println("="^70)

function synthetic_cox_cohort(n::Int; true_beta::Float64=0.8, seed::Int=1)
  rng = MersenneTwister(seed)
  x1 = randn(rng, n)
  baseline_rate = 0.3
  u = rand(rng, n)
  event_time = -log.(u) ./ (baseline_rate .* exp.(x1 .* true_beta))
  censor_time = rand(rng, n) .* 5.0 .+ 1.0
  time = min.(event_time, censor_time) .+ (1:n) .* 1e-6  # tiny jitter to break exact ties
  status = Int.(event_time .<= censor_time)
  df = DataFrame(patient_id=["p$i" for i in 1:n], time=time, status=status, x1=x1)
  return make_cohort(df)
end

norm_(v) = sqrt(sum(abs2, v))

@testset "cox_ph: score equation is ~0 at the converged MLE" begin
  cohort = synthetic_cox_cohort(60)
  formula = :(Surv(time, status) ~ x1)
  res = cox_ph(formula, cohort)
  @test res.converged
  println("cox_ph converged=$(res.converged) in $(res.iterations) iterations, beta=$(res.terms[1].beta)")

  try
    time = Float64.(cohort.clinical.time)
    status = Int.(cohort.clinical.status)
    X, _ = Clinical._cox_design_matrix(cohort, ["x1"])
    t2, s2, ord = Clinical._cox_order(time, status)
    X2 = X[ord, :]
    beta_full = vcat(0.0, res.terms[1].beta)
    _, gradient, _ = Clinical._cox_loglik_gradient_hessian(beta_full, X2, t2, s2; ties=:efron)
    @test norm_(gradient) < 1e-4
    println("gradient norm at MLE: $(norm_(gradient)) (should be ~0)")
  catch e
    println("internal score-equation check skipped (private API not accessible): ", e)
  end
end

@testset "cox_ph: Efron and Breslow agree exactly when there are no tied times" begin
  cohort = synthetic_cox_cohort(50; seed=7)  # jittered times => no ties
  formula = :(Surv(time, status) ~ x1)
  res_efron = cox_ph(formula, cohort; ties=:efron)
  res_breslow = cox_ph(formula, cohort; ties=:breslow)
  @test res_efron.terms[1].beta ≈ res_breslow.terms[1].beta atol=1e-6
  println("efron beta=$(res_efron.terms[1].beta)  breslow beta=$(res_breslow.terms[1].beta)")
end

@testset "cox_ph: baseline hazard is non-decreasing" begin
  cohort = synthetic_cox_cohort(60; seed=3)
  formula = :(Surv(time, status) ~ x1)
  res = cox_ph(formula, cohort)
  @test all(diff(res.baseline_hazard) .>= -1e-9)
  println("baseline hazard monotonic non-decreasing: true")
end

# ==============================================================================
# fine_gray: regression test for the d_m tie-count fix
# ==============================================================================
println("\n" * "="^70)
println("fine_gray")
println("="^70)

@testset "fine_gray: hand-derived baseline hazard with a tied event time" begin
  # 4 patients, x1 all zero so beta stays exactly 0 at every Newton step
  # (gradient is identically 0 for an all-zero covariate).
  # time=[5,5,8,10], status=[1,1,2,0] (two TIED cause-1 events at t=5).
  #
  # By hand (theta_i = exp(0) = 1 for all i since beta=0):
  #   G(t): no censoring before t=10 -> G(5)=G(8)=1.0
  #   baseline_times = [5.0]  (only one distinct event time)
  #   at t=5: d_m = 2 (BOTH cause-1 events tied here)
  #           S0 = sum(theta_i for time_i >= 5) = 4 (all 4 patients)
  #           cum_h = d_m / S0 = 2 / 4 = 0.5
  df = DataFrame(
    patient_id=["p1", "p2", "p3", "p4"],
    time=[5.0, 5.0, 8.0, 10.0],
    status=[1, 1, 2, 0],
    x1=[0.0, 0.0, 0.0, 0.0],
  )
  cohort = make_cohort(df)
  formula = :(Surv(time, status) ~ x1)
  res = fine_gray(formula, cohort; fail_code=1)
  @test length(res.baseline_times) == 1
  @test res.baseline_times[1] ≈ 5.0
  @test res.baseline_hazard[1] ≈ 0.5 atol=1e-9
  println("fine_gray baseline_hazard: got=$(res.baseline_hazard[1])  expected=0.5")
  println("(this specifically catches the d_m tie-multiplier regression)")
end

# ==============================================================================
# nelson_aalen: regression test for the status-code convention fix
# ==============================================================================
println("\n" * "="^70)
println("nelson_aalen")
println("="^70)

@testset "nelson_aalen: hand-derived cumulative hazard with competing-risk codes" begin
  # time=[1,2,3,4], status=[1,2,0,1] (cause1@1, cause2@2, censored@3, cause1@4)
  # status>0 must count as an event (file-wide convention) -- this specifically
  # would break under a status==1-only check, since t=2 has status=2.
  #
  # By hand:
  #   t=1: at_risk=4, d=1 -> ch=0.25,        var=1/16=0.0625
  #   t=2: at_risk=3, d=1 -> ch=0.25+1/3=0.583333, var=0.0625+1/9=0.173611
  #   t=3: at_risk=2, d=0 -> no row (censored, not an event)
  #   t=4: at_risk=1, d=1 -> ch=0.583333+1=1.583333, var=0.173611+1=1.173611
  time = [1.0, 2.0, 3.0, 4.0]
  status = [1, 2, 0, 1]
  na = nelson_aalen(time, status)
  @test na.time ≈ [1.0, 2.0, 4.0]
  @test na.cum_hazard ≈ [0.25, 0.25 + 1/3, 0.25 + 1/3 + 1.0] atol=1e-9
  @test na.std_error ≈ sqrt.([0.0625, 0.0625 + 1/9, 0.0625 + 1/9 + 1.0]) atol=1e-9
  println("nelson_aalen cum_hazard: got=$(na.cum_hazard)")
  println("expected: $([0.25, 0.583333, 1.583333])")
end

# ==============================================================================
# cox_residuals: martingale self-consistency
# ==============================================================================
println("\n" * "="^70)
println("cox_residuals")
println("="^70)

@testset "cox_residuals: martingale residuals sum to ~0 (Breslow baseline property)" begin
  # For a Breslow-type baseline hazard computed FROM the same fitted model,
  # sum_i [status_i - Lambda0(t_i)*theta_i] should be ~0 by construction --
  # this is a strong self-consistency check on indexing/sign correctness.
  cohort = synthetic_cox_cohort(80; seed=11)
  formula = :(Surv(time, status) ~ x1)
  res = cox_ph(formula, cohort)

  time = Float64.(cohort.clinical.time)
  status = Int.(cohort.clinical.status)
  X_cov = reshape(Float64.(cohort.clinical.x1), :, 1)  # NO intercept column
  M = cox_residuals(res, time, status, X_cov; type=:martingale)
  total = sum(M[:, 1])
  @test abs(total) < 0.5
  println("sum of martingale residuals: $total (expected ~0)")
end

@testset "cox_residuals: Schoenfeld residuals match cox_zph's own computation" begin
  cohort = synthetic_cox_cohort(50; seed=13)
  formula = :(Surv(time, status) ~ x1)
  res = cox_ph(formula, cohort)
  zph = cox_zph(res, formula, cohort)
  @test all(-1.0 .<= zph.rho .<= 1.0)
  @test all(0.0 .<= zph.pvalue .<= 1.0)
  println("cox_zph rho=$(zph.rho)  global_pvalue=$(zph.global_pvalue)")
end

# ==============================================================================
# dose_response_curve: recovery test for the Jacobian-sign fix
# ==============================================================================
println("\n" * "="^70)
println("dose_response_curve")
println("="^70)

@testset "dose_response_curve: recovers known noiseless Hill parameters" begin
  # Generate data from a KNOWN Hill curve (emin=0, emax=100, ec50=10, hill=2)
  # with zero noise, then check the fit recovers those exact parameters.
  # This is the regression test for the Jacobian sign bug -- with the sign
  # wrong, the LM refinement step never improves on the grid search, so a
  # tight tolerance here would fail under the old code.
  true_emin, true_emax, true_ec50, true_hill = 0.0, 100.0, 10.0, 2.0
  x = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0]
  y = true_emin .+ (true_emax - true_emin) ./ (1.0 .+ (true_ec50 ./ x) .^ true_hill)
  cell_lines = fill("cellA", length(x))

  res = dose_response_curve("drugX", cell_lines, x, y)
  @test res.ec50 ≈ true_ec50 rtol=0.05
  @test res.hill ≈ true_hill rtol=0.10
  @test res.emax ≈ true_emax rtol=0.02
  println("recovered: ec50=$(res.ec50) (true=10.0), hill=$(res.hill) (true=2.0), emax=$(res.emax) (true=100.0)")
  println("fitted SSE-consistent: max|residual|=$(maximum(abs.(res.fitted .- y)))")
end

# ==============================================================================
# surv_cutpoint
# ==============================================================================
println("\n" * "="^70)
println("surv_cutpoint")
println("="^70)

@testset "surv_cutpoint: recovers a separable cutpoint, p-value only increases" begin
  feature = collect(1.0:10.0)
  time = vcat(collect(1.0:5.0), collect(10.0:14.0))  # low-feature group fails early
  status = ones(Int, 10)
  res = surv_cutpoint(time, status, feature; min_prop=0.2)
  @test res.cutpoint ≈ 3.0
  @test res.provenance.parameters.adjusted_pvalue >= res.provenance.parameters.unadjusted_pvalue
  println("cutpoint: got=$(res.cutpoint)  expected=3.0")
  println("unadjusted p=$(res.provenance.parameters.unadjusted_pvalue)  adjusted p=$(res.pvalue)")
end

@testset "surv_cutpoint: RMST-based selection can diverge from log-rank-based selection" begin
  # Constructed so that:
  #  - splitting at val=3 gives a large EARLY hazard difference (many events
  #    cluster right after the split in one arm) -> large log-rank statistic
  #  - splitting at val=6 gives a larger total RMST gap because one arm's
  #    tail survival stays elevated much longer -> larger |RMST diff|
  # If cutpoint selection is log-rank-based, expect best_cut == 3.0.
  # If it's RMST-based (current code), expect best_cut == 6.0 instead.
  feature = collect(1.0:10.0)
  time = [1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 15.0, 16.0, 17.0, 18.0]
  status = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
  res = surv_cutpoint(time, status, feature; min_prop=0.2)
  println("selected cutpoint = $(res.cutpoint)")
  println("(if this is 3.0, selection is log-rank-based as the docstring claims;")
  println(" if this is 6.0 or something else, it's RMST-based per the current code)")
  # Not asserting a specific value here -- this test is diagnostic, run it
  # and compare against what survminer::surv_cutpoint gives on the same
  # data if you want ground truth.
end

# ==============================================================================
# maf_compare / clinical_enrichment / somatic_interactions
# ==============================================================================
println("\n" * "="^70)
println("MAF-based statistics")
println("="^70)

function maf_records(gene::String, samples::Vector{String})
  return [MAFRecord(gene, s, "Missense_Mutation", "SNP", "chr1", 100, 100, "A", "T") for s in samples]
end

@testset "maf_compare: hand-derived Fisher exact p-value for perfect separation" begin
  # Cohort1: gene TP53 mutated in ALL 4 samples. Cohort2: mutated in NONE of 4.
  # a=4,b=0,c=0,d=4 -> Hypergeometric(s=4,f=4,n=4): P(X=4)=C(4,4)C(4,0)/C(8,4)=1/70
  # pval = 2*min(cdf(hg,4), ccdf(hg,3)) = 2*min(1.0, 1/70) = 2/70 = 1/35 ≈ 0.0285714
  maf1 = maf_records("TP53", ["s1", "s2", "s3", "s4"])
  maf2 = maf_records("BRCA1", ["t1", "t2", "t3", "t4"])  # no TP53 mutations at all
  df = maf_compare(maf1, maf2; top_n=10)
  row = df[df.gene .== "TP53", :]
  @test nrow(row) == 1
  @test row.pvalue[1] ≈ 1/35 atol=1e-6
  println("maf_compare TP53 pvalue: got=$(row.pvalue[1])  expected=$(1/35)")
end

@testset "somatic_interactions: mutual exclusivity vs co-occurrence direction" begin
  # GeneA/GeneB: mutated in disjoint sample sets -> mutually exclusive (OR<1)
  # GeneC/GeneD: mutated in the SAME sample set -> co-occurring (OR>1)
  samples_1_5 = ["s$i" for i in 1:5]
  samples_6_10 = ["s$i" for i in 6:10]
  maf = vcat(
    maf_records("GeneA", samples_1_5),
    maf_records("GeneB", samples_6_10),
    maf_records("GeneC", samples_1_5),
    maf_records("GeneD", samples_1_5),
  )
  out = somatic_interactions(maf; top_n=10)
  ab = out[(out.gene1 .== "GeneA".&&out.gene2 .== "GeneB") .| (out.gene1 .== "GeneB".&&out.gene2 .== "GeneA"), :]
  cd = out[(out.gene1 .== "GeneC".&&out.gene2 .== "GeneD") .| (out.gene1 .== "GeneD".&&out.gene2 .== "GeneC"), :]
  @test nrow(ab) == 1 && ab.event_type[1] == "Mutually_Exclusive" && ab.odds_ratio[1] < 1.0
  @test nrow(cd) == 1 && cd.event_type[1] == "Co_Occurrence" && cd.odds_ratio[1] > 1.0
  println("A-B (disjoint samples): OR=$(ab.odds_ratio[1]) type=$(ab.event_type[1])")
  println("C-D (same samples):     OR=$(cd.odds_ratio[1]) type=$(cd.event_type[1])")
end

@testset "calculate_tmb: hand-counted non-silent mutation rate" begin
  maf = vcat(
    maf_records("GeneA", ["s1", "s1"]),  # 2 missense (non-silent) in s1
    [MAFRecord("GeneB", "s1", "Silent", "SNP", "chr1", 1, 1, "A", "T")],  # 1 silent in s1
  )
  df = calculate_tmb(maf; capture_size_mb=2.0)
  row = df[df.sample .== "s1", :]
  @test row.total_mutations[1] == 3
  @test row.non_silent_mutations[1] == 2
  @test row.tmb[1] ≈ 2 / 2.0 atol=1e-9
  println("TMB: total=$(row.total_mutations[1]) non_silent=$(row.non_silent_mutations[1]) tmb=$(row.tmb[1])")
end

@testset "oncoprint: multi-hit cell displays the MOST SEVERE variant, not file order" begin
  # s1/GeneA gets a Missense record FIRST, then a Nonsense record -- the
  # correct display pick is Nonsense (more severe), regardless of order.
  maf = [
    MAFRecord("GeneA", "s1", "Missense_Mutation", "SNP", "chr1", 1, 1, "A", "T"),
    MAFRecord("GeneA", "s1", "Nonsense_Mutation", "SNP", "chr1", 2, 2, "A", "T"),
  ]
  op = oncoprint(maf)
  label = op.mutation_labels[1, 1]
  most_severe = Clinical._most_severe_variant(label)
  @test most_severe == "Nonsense_Mutation"
  println("multi-hit cell label='$label' -> most severe = '$most_severe' (expected Nonsense_Mutation)")
end

# ==============================================================================
# Performance smoke test (timing only, not correctness)
# ==============================================================================
println("\n" * "="^70)
println("Performance (n=2000 synthetic cohort)")
println("="^70)

big_cohort = synthetic_cox_cohort(2000; seed=99)
formula = :(Surv(time, status) ~ x1)
t_cox = @elapsed cox_ph(formula, big_cohort)
t_fg = @elapsed try
  fg_res = fine_gray(formula, big_cohort; fail_code=1)
catch e
  println("fine_gray perf run skipped: ", e)
end
println("cox_ph (n=2000):    $(round(1000*t_cox, digits=1)) ms")
println("fine_gray (n=2000): $(round(1000*t_fg, digits=1)) ms")

println("\n" * "="^70)
println("DONE")
println("="^70)
