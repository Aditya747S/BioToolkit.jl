using SparseArrays
using DataFrames
using Plots

@testset "Clinical" begin
    clinical = DataFrame(
        patient_id=["P1", "P2", "P3"],
        age=[60.0, 55.0, 70.0],
        stage=["IV", "II", "IV"],
        time=[10.0, 8.0, 12.0],
        status=[1, 0, 1],
    )
    genomics = BioToolkit.CountMatrix(sparse([10 20 30; 1 0 2]), ["G1", "G2"], ["P1", "P2", "P3"])
    cohort = BioToolkit.PatientCohort(clinical, genomics, ["P1", "P2", "P3"])

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

    cox = BioToolkit.cox_ph(:(Surv(time, status) ~ age), cohort)
    @test !isempty(cox.terms)
    @test all(term -> isfinite(term.beta) && isfinite(term.hazard_ratio), cox.terms)
    forest = BioToolkit.forest_plot(cox)
    @test forest isa Plots.Plot

    maf_path = tempname()
    open(maf_path, "w") do io
        println(io, join(["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "Variant_Type", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"], '\t'))
        println(io, join(["TP53", "P1", "Missense_Mutation", "SNP", "17", "7579472", "7579472", "C", "T"], '\t'))
        println(io, join(["TP53", "P3", "Nonsense_Mutation", "SNP", "17", "7579472", "7579472", "C", "A"], '\t'))
    end
    maf = BioToolkit.read_maf(maf_path)
    @test length(maf) == 2
    summary = BioToolkit.summarize_maf(maf)
    @test summary.total_mutations == 2
    @test summary.per_sample["P1"] == 1
    @test summary.per_gene["TP53"] == 2
    oncoprint = BioToolkit.oncoprint(maf)
    @test size(oncoprint.matrix) == (1, 2)
    @test oncoprint.mutation_labels[1, 1] != ""
    oncoprint_plot = plot(oncoprint)
    @test oncoprint_plot isa Plots.Plot

    tcga_a = tempname()
    tcga_b = tempname()
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

    roc = BioToolkit.survival_roc(clinical.time, clinical.status, clinical.age, 9.0)
    @test 0.0 <= roc.auc <= 1.0

    cif = BioToolkit.cif_curve(clinical.time, clinical.status, [1, 0, 2])
    @test !isempty(cif.time)
    @test haskey(cif.cumulative_incidence, 1)

    response = BioToolkit.dose_response_curve("drugA", ["CL1", "CL2", "CL3", "CL4"], [0.1, 0.5, 1.0, 10.0], [0.9, 0.7, 0.5, 0.2])
    @test length(response.fitted) == 4
    @test response.ic50 > 0

    neural = BioToolkit.neural_cox([1.0 0.2; 0.4 1.2; 0.7 0.9], [4.0, 3.0, 2.0], [1, 0, 1]; epochs=25)
    @test length(neural.risk_scores) == 3
end