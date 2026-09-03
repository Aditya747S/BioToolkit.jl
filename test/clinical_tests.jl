using SparseArrays
using DataFrames
using Plots
using Test
using BioToolkit

@testset "Clinical" begin
    clinical = DataFrame(
        patient_id=["P1", "P2", "P3", "P4"],
        age=[60.0, 55.0, 70.0, 65.0],
        stage=["IV", "II", "IV", "I"],
        time=[10.0, 8.0, 12.0, 15.0],
        status=[1, 0, 1, 1],
        site=["A", "A", "B", "B"]
    )
    genomics = BioToolkit.CountMatrix(sparse([10 20 30 40; 1 0 2 3]), ["G1", "G2"], ["P1", "P2", "P3", "P4"])
    cohort = BioToolkit.PatientCohort(clinical, genomics, ["P1", "P2", "P3", "P4"])

    clinical_row, genomics_column = cohort["P1"]
    @test clinical_row.patient_id == "P1"
    @test length(genomics_column) == 2

    subset = cohort[cohort.clinical.stage .== "IV"]
    @test length(subset.patient_ids) == 2
    @test subset.genomics.sample_ids == ["P1", "P3"]

    km = BioToolkit.kaplan_meier(clinical.time, clinical.status)
    @test !isempty(km.time)
    @test all(0 .<= km.survival .<= 1)
    km_plot = BioToolkit.kaplan_meier_plot(km)
    @test km_plot isa Plots.Plot

    logrank = BioToolkit.logrank_test(clinical.time, clinical.status, clinical.stage)
    @test isfinite(logrank.statistic)
    @test 0.0 <= logrank.pvalue <= 1.0

    # Cox PH with ties & cluster test
    cox = BioToolkit.cox_ph(:(Surv(time, status) ~ age), cohort; ties=:efron, cluster=:site)
    @test !isempty(cox.terms)
    @test all(term -> isfinite(term.beta) && isfinite(term.hazard_ratio), cox.terms)
    @test isfinite(cox.loglik)
    forest = BioToolkit.forest_plot(cox)
    @test forest isa Plots.Plot

    # Cox ZPH diagnostics
    zph = BioToolkit.cox_zph(cox, :(Surv(time, status) ~ age), cohort)
    @test length(zph.rho) == 1
    @test isfinite(zph.global_chisq)

    # MAF reading, summarizing, TMB, somatic interactions
    maf_path = Base.tempname()
    open(maf_path, "w") do io
        println(io, join(["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "Variant_Type", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"], '\t'))
        println(io, join(["TP53", "P1", "Missense_Mutation", "SNP", "17", "7579472", "7579472", "C", "T"], '\t'))
        println(io, join(["TP53", "P3", "Nonsense_Mutation", "SNP", "17", "7579472", "7579472", "C", "A"], '\t'))
        println(io, join(["EGFR", "P1", "Silent", "SNP", "7", "5500000", "5500000", "A", "G"], '\t'))
        println(io, join(["KRAS", "P2", "Missense_Mutation", "SNP", "12", "25398284", "25398284", "C", "T"], '\t'))
    end
    maf = BioToolkit.read_maf(maf_path)
    @test length(maf) == 4
    summary = BioToolkit.summarize_maf(maf)
    @test summary.total_mutations == 4
    @test summary.per_sample["P1"] == 2
    @test summary.per_gene["TP53"] == 2

    # TMB & Somatic Interactions
    tmb = BioToolkit.calculate_tmb(maf; capture_size_mb=38.0)
    @test nrow(tmb) == 3
    @test "tmb" in names(tmb)

    interactions = BioToolkit.somatic_interactions(maf; top_n=5)
    @test nrow(interactions) >= 1
    @test "event_type" in names(interactions)

    # Oncoprint & plot with kwargs syntax fix
    oncoprint = BioToolkit.oncoprint(maf)
    @test size(oncoprint.matrix) == (3, 3)
    @test oncoprint.mutation_labels[findfirst(==("TP53"), oncoprint.genes), findfirst(==("P1"), oncoprint.samples)] != ""
    op_plot = BioToolkit.oncoprint_plot(oncoprint; title="Custom Oncoprint Title")
    @test op_plot isa Plots.Plot

    # Stored XSS Protection Test in HTML Output
    malicious_cox = BioToolkit.CoxResult(
        [BioToolkit.CoxTermResult("<script>alert('xss')</script>", 0.5, 1.65, 0.2, 2.5, 0.01, 1.1, 2.4)],
        [10.0], [0.1], -12.5, 5, true
    )
    cox_html = BioToolkit.to_html(malicious_cox)
    @test contains(cox_html, "escapeHtml")
    @test contains(cox_html, "escapeHtml(t.term)")

    # TCGA count file ingestion test
    tcga_a = Base.tempname()
    tcga_b = Base.tempname()
    open(tcga_a, "w") do io
        println(io, "gene\tcount")
        println(io, "TP53\t10")
        println(io, "EGFR\t5")
    end
    open(tcga_b, "w") do io
        println(io, "gene\tcount")
        println(io, "TP53\t7")
        println(io, "KRAS\t2")
    end
    tcga_counts = BioToolkit.tcga_ingest([tcga_a, tcga_b], ["S1", "S2"])
    @test tcga_counts isa BioToolkit.CountMatrix
    @test tcga_counts.sample_ids == ["S1", "S2"]
    @test sort(tcga_counts.gene_ids) == ["EGFR", "KRAS", "TP53"]
    @test sum(tcga_counts.counts) == 24

    # Survival ROC (IPCW-weighted)
    roc = BioToolkit.survival_roc(clinical.time, clinical.status, clinical.age, 9.0)
    @test 0.0 <= roc.auc <= 1.0

    # Competing Risks CIF & Fine-Gray
    cif = BioToolkit.cif_curve(clinical.time, clinical.status, [1, 0, 2, 1])
    @test !isempty(cif.time)
    @test haskey(cif.cumulative_incidence, 1)

    fg = BioToolkit.fine_gray(:(Surv(time, status) ~ age), cohort; fail_code=1)
    @test fg.converged

    # Dose-Response LM Optimization
    response = BioToolkit.dose_response_curve("drugA", ["CL1", "CL2", "CL3", "CL4"], [0.1, 0.5, 1.0, 10.0], [0.9, 0.7, 0.5, 0.2])
    @test length(response.fitted) == 4
    @test response.ic50 > 0

    # Neural Cox fast O(N log N)
    neural = BioToolkit.neural_cox([1.0 0.2; 0.4 1.2; 0.7 0.9; 0.5 0.5], [4.0, 3.0, 2.0, 5.0], [1, 0, 1, 1]; epochs=25)
    @test length(neural.risk_scores) == 4

    # Propensity Score Matching with missing data handling
    ps_data = DataFrame(
        treatment=[1, 1, 0, 0, 1, 0],
        age=[50.0, 60.0, 52.0, 58.0, missing, 55.0],
        bmi=[22.0, 28.0, 23.0, 27.0, 25.0, 24.0]
    )
    ps_drop = BioToolkit.propensity_score_match(ps_data; missing_handling=:drop)
    @test nrow(ps_drop.matches) >= 1

    ps_imp = BioToolkit.propensity_score_match(ps_data; missing_handling=:impute)
    @test nrow(ps_imp.matches) >= 1

    # RMST Test
    rmst_res = BioToolkit.rmst(clinical.time, clinical.status; tau=12.0, groups=clinical.stage .== "IV")
    @test rmst_res.rmst > 0
    @test rmst_res.group_results !== nothing
    @test haskey(rmst_res.group_results, "diff")

    # Bioconductor Parity Extensions Tests
    # 1. Surv Cutpoint
    cut_res = BioToolkit.surv_cutpoint([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [1, 1, 1, 0, 1, 0], [10.5, 12.0, 25.0, 28.0, 30.0, 35.0])
    @test cut_res.cutpoint > 0
    @test cut_res.high_group_count + cut_res.low_group_count == 6

    # 2. MAF Compare
    m1 = [BioToolkit.MAFRecord("TP53", "S1", "Missense_Mutation", "SNP", "chr17", 7577120, 7577120, "A", "T")]
    m2 = [BioToolkit.MAFRecord("KRAS", "S2", "Missense_Mutation", "SNP", "chr12", 25398284, 25398284, "G", "A")]
    comp_df = BioToolkit.maf_compare(m1, m2; cohort1_name="CohortA", cohort2_name="CohortB")
    @test nrow(comp_df) >= 1
    @test hasproperty(comp_df, :odds_ratio)

    # 3. Clinical Enrichment
    clin_df = DataFrame(patient_id=["S1", "S2"], stage=["I", "IV"])
    enrich_df = BioToolkit.clinical_enrichment(m1, clin_df, :stage)
    @test nrow(enrich_df) >= 1
    @test hasproperty(enrich_df, :odds_ratio)
    @test hasproperty(enrich_df, :pvalue)
    @test hasproperty(enrich_df, :fdr)
    @test all(0.0 .<= enrich_df.pvalue .<= 1.0)

    # 4. Nelson-Aalen Cumulative Hazard (including multi-cause event codes)
    na_res = BioToolkit.nelson_aalen([1.0, 2.0, 3.0, 4.0], [1, 0, 2, 1])
    @test length(na_res.cum_hazard) >= 1
    @test all(na_res.cum_hazard .>= 0.0)

    # 5. Cox Residuals (Schoenfeld, Martingale, Deviance)
    sch_resids = BioToolkit.cox_residuals(cox, clinical.time, clinical.status, [clinical.age clinical.time]; type=:schoenfeld)
    @test size(sch_resids, 1) == length(clinical.time)
    
    mart_resids = BioToolkit.cox_residuals(cox, clinical.time, clinical.status, [clinical.age clinical.time]; type=:martingale)
    @test size(mart_resids, 1) == length(clinical.time)
    
    dev_resids = BioToolkit.cox_residuals(cox, clinical.time, clinical.status, [clinical.age clinical.time]; type=:deviance)
    @test size(dev_resids, 1) == length(clinical.time)
end