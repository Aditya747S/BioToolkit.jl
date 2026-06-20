using SparseArrays
using LinearAlgebra
using Random
using Statistics
using Logging: with_logger, NullLogger
using DataFrames

function _build_synthetic_count_matrix()
    gene_ids = ["gene$(index)" for index in 1:51]
    sample_ids = ["sample$(index)" for index in 1:4]
    counts = zeros(Int, length(gene_ids), length(sample_ids))

    for gene in 2:length(gene_ids)
        counts[gene, 1] = 120 + mod(gene, 7)
        counts[gene, 2] = 118 + mod(gene, 5)
        counts[gene, 3] = 122 + mod(gene, 3)
        counts[gene, 4] = 121 + mod(gene, 11)
    end

    counts[1, 1] = 2
    counts[1, 2] = 3
    counts[1, 3] = 300
    counts[1, 4] = 320

    return BioToolkit.CountMatrix(sparse(counts), gene_ids, sample_ids)
end

@testset "DifferentialExpression basics" begin
    count_matrix = _build_synthetic_count_matrix()
    design = [:Control, :Control, :Treat, :Treat]

    # DESeq2-faithful default normalization uses median-ratio factors.
    dds_default = BioToolkit.DESeqDataSet(count_matrix, DataFrame(condition=design), design)
    @test dds_default.sf_type == :ratio

    norm_factors = BioToolkit.calc_norm_factors(count_matrix)
    @test length(norm_factors) == 4
    @test all(factor -> factor > 0, norm_factors)

    dispersions = BioToolkit.estimate_dispersions(count_matrix, norm_factors)
    @test length(dispersions) == size(count_matrix.counts, 1)
    @test all(dispersion -> dispersion > 0, dispersions)

    bh = BioToolkit.benjamini_hochberg([0.01, 0.04, 0.2, 0.5])
    @test bh[1] <= bh[2] <= bh[3] <= bh[4]
    @test all(value -> 0.0 <= value <= 1.0, bh)

    results = BioToolkit.differential_expression(count_matrix, design)
    @test length(results) == length(count_matrix.gene_ids)

    gene1 = results[1]
    gene2 = results[2]
    @test gene1.gene_id == "gene1"
    @test gene1.log2_fold_change > 0.5
    @test gene1.padj < 0.05
    @test gene2.padj > 0.05
    @test gene1.base_mean > 0
    @test gene1.lfc_se > 0
    @test gene1.stat >= 0

    filtered_counts = BioToolkit.DifferentialExpression.filter_low_counts(count_matrix; min_total = 500)
    @test size(filtered_counts.counts, 1) == 1
    @test filtered_counts.gene_ids == ["gene1"]

    design_matrix = hcat(
        ones(Float64, length(design)),
        Float64[condition == :Treat ? 1.0 : 0.0 for condition in design],
    )
    solver = BioToolkit.DifferentialExpression.GLMSolver(design_matrix)
    beta = zeros(size(design_matrix, 2))
    fitted_beta, fitted_se, fitted_stat = BioToolkit.DifferentialExpression.fit_gene_fast!(
        solver,
        collect(count_matrix.counts[1, :]),
        log.(norm_factors),
        0.1,
    )
    @test all(isfinite, fitted_beta)
    @test isfinite(fitted_se)
    @test isfinite(fitted_stat)

    shrunk = BioToolkit.DifferentialExpression.shrink_lfc(results)
    @test length(shrunk) == length(results)
    @test abs(shrunk[1].log2_fold_change) <= abs(results[1].log2_fold_change)

    transformed = BioToolkit.DifferentialExpression.vst(count_matrix, dispersions; norm_factors = norm_factors)
    @test size(transformed) == size(Matrix(count_matrix.counts))
    @test all(isfinite, transformed)
end

@testset "DifferentialExpression batch correction" begin
    rng = MersenneTwister(42)
    genes = 24
    samples = 8
    truth = randn(rng, genes, 1) * ones(1, samples)
    biology = zeros(Float64, genes, samples)
    biology[:, 1:4] .+= 3.0
    biology[:, 5:8] .-= 3.0
    batch_effect = zeros(Float64, genes, samples)
    batch_effect[:, 1:2:8] .+= 2.5
    batch_effect[:, 2:2:8] .-= 2.5
    data = truth .+ biology .+ batch_effect

    bio_design = reshape([fill(0.0, 4); fill(1.0, 4)], :, 1)
    batch_labels = [isodd(index) ? "batch1" : "batch2" for index in 1:samples]

    corrected_linear = BioToolkit.remove_batch_effect(data, batch_labels; bio_design=bio_design)
    @test corrected_linear ≈ truth .+ biology atol=1e-8

    corrected_combat = BioToolkit.combat_correction(data, batch_labels; bio_design=bio_design)
    batch_gap_before = abs(mean(data[:, isodd.(1:samples)]) - mean(data[:, iseven.(1:samples)]))
    batch_gap_after = abs(mean(corrected_combat[:, isodd.(1:samples)]) - mean(corrected_combat[:, iseven.(1:samples)]))
    @test size(corrected_combat) == size(data)
    @test batch_gap_after < batch_gap_before

    hidden_factor = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
    hidden_loading = randn(rng, genes)
    surrogate_data = data .+ hidden_loading * hidden_factor'
    surrogates = BioToolkit.estimate_surrogates(surrogate_data; bio_design=bio_design, nperm=0, max_components=1, rng=rng)
    @test size(surrogates) == (samples, 1)

    corrected_surrogate = BioToolkit.remove_batch_effect(surrogate_data, surrogates; bio_design=bio_design)
    @test norm(corrected_surrogate * hidden_factor) < 0.1 * norm(surrogate_data * hidden_factor)
end

@testset "DifferentialExpression validation" begin
    count_matrix = _build_synthetic_count_matrix()
    @test_throws ArgumentError BioToolkit.CountMatrix(sparse(zeros(Int, 2, 3)), ["g1"], ["s1", "s2", "s3"])
    @test_throws ArgumentError BioToolkit.differential_expression(count_matrix, [:Control, :Treat])

    single_gene = BioToolkit.CountMatrix(sparse([0 0 0 0]), ["geneX"], ["s1", "s2", "s3", "s4"])
    @test_throws ArgumentError BioToolkit.differential_expression(single_gene, [:Control, :Control, :Treat, :Treat])
    single_result = BioToolkit.differential_expression(single_gene, [:Control, :Control, :Treat, :Treat]; sfType=:poscounts)
    @test single_result[1].pvalue == 1.0
    @test single_result[1].padj == 1.0
end

@testset "DifferentialExpression edgeR and IHW" begin
    count_matrix = _build_synthetic_count_matrix()
    design = [:Control, :Control, :Treat, :Treat]

    disp_fit = BioToolkit.estimate_dispersion_edgeR(count_matrix, design)
    @test haskey(Dict(pairs(disp_fit)), :tagwise)
    @test length(disp_fit.tagwise) == length(count_matrix.gene_ids)

    et = BioToolkit.exact_test_edgeR(
        count_matrix,
        design;
        norm_factors=disp_fit.norm_factors,
        dispersion=disp_fit.tagwise,
    )

    @test size(et, 1) == length(count_matrix.gene_ids)
    @test all(v -> 0.0 <= v <= 1.0, et.pvalue)
    @test all(v -> 0.0 <= v <= 1.0, et.padj)
    @test et[1, :logFC] > 1.0
    @test abs(et[1, :logFC]) == maximum(abs.(et.logFC))

    pvals = vcat(fill(0.8, 150), fill(1e-6, 30), fill(0.2, 20))
    fstat = collect(1:length(pvals))

    q_bh = BioToolkit.benjamini_hochberg(pvals)
    q_ihw = BioToolkit.DifferentialExpression.ihw_qvalue(pvals, fstat; n_bins=10, n_folds=5)

    @test length(q_ihw) == length(pvals)
    @test all(v -> 0.0 <= v <= 1.0, q_ihw)
    @test count(<(0.1), q_ihw) >= count(<(0.1), q_bh)
end

@testset "DifferentialExpression math regressions" begin
    exact_left = BioToolkit.DifferentialExpression._edgeR_exact_twosided_pvalue(10, 2, 2.0, 8.0)
    exact_right = BioToolkit.DifferentialExpression._edgeR_exact_twosided_pvalue(10, 2, 8.0, 2.0)
    @test exact_left ≈ exact_right atol=1e-12

    spline_edge = BioToolkit.DifferentialExpression._natural_cubic_spline_eval([0.5, 0.95], [0.8, 1.2], 1.0)
    @test spline_edge ≈ 1.2444444444444445 atol=1e-12

    beta_matrix = [0.0 0.0; 10.0 1.0; -10.0 -1.0]
    covariance_matrices = [
        [0.0 0.0; 0.0 0.0],
        [1.0 0.0; 0.0 0.25],
        [1.0 0.0; 0.0 0.25],
    ]
    prior_var = BioToolkit.DifferentialExpression.estimate_beta_prior_var(beta_matrix, covariance_matrices; coef_idx=2)
    @test prior_var ≈ 1.75 atol=1e-10

    design_matrix = [1.0 0.0; 1.0 0.0; 1.0 1.0; 1.0 1.0]
    solver = BioToolkit.DifferentialExpression.GLMSolver(design_matrix)
    fitted_beta, fitted_cov, _, _ = BioToolkit.DifferentialExpression.fit_gene_fast!(
        solver,
        [0.0, 0.0, 10.0, 10.0],
        zeros(4),
        0.1;
        minmu=5.0,
        return_covariance=true,
    )
    @test minimum(exp.(design_matrix * fitted_beta)) < 5.0
    @test solver.eta ≈ design_matrix * fitted_beta atol=1e-8
    @test size(fitted_cov) == (2, 2)

    squeezed_ql, df_prior, ql_prior = BioToolkit.DifferentialExpression._squeeze_quasi_dispersions([0.1, 0.2, 5.0, 0.15, 0.18], 3)
    @test df_prior >= 0.0
    @test length(ql_prior) == 5
    @test squeezed_ql[3] < 5.0

    trend_values, trend_type = BioToolkit.DifferentialExpression._fit_dispersion_trend(
        [1.0, 2.0, 4.0, 8.0, 16.0, 32.0],
        [10.0, 0.2, 9.0, 0.25, 8.0, 0.3];
        fit_type=:parametric,
    )
    @test trend_type == :local
    @test length(trend_values) == 6

    count_matrix = _build_synthetic_count_matrix()
    design = [:Control, :Control, :Treat, :Treat]
    dds_prior = BioToolkit.DESeqDataSet(count_matrix, DataFrame(condition=design), design)
    with_logger(NullLogger()) do
        BioToolkit.DESeq(dds_prior; betaPrior=true, betaPriorVar=0.25, maxit=50, reference_level=:Control, target_level=:Treat)
    end
    @test dds_prior.metadata["betaPriorVarScale"] == "log2"
    @test dds_prior.metadata["betaPriorVarNatural"] ≈ 0.25 * log(2)^2 atol=1e-12

    beta_input = reshape([0.0, 3.0], 1, 2)
    cov_input = [[0.05 0.0; 0.0 0.2]]
    res_input = [BioToolkit.DEResult("g1", 10.0, 3.0 / log(2), sqrt(0.2) / log(2), 0.0, 0.01, NaN, false, true)]
    shrunk_beta, shrunk_cov, shrunk_res = BioToolkit.DifferentialExpression._apply_beta_prior_results(
        copy(beta_input),
        copy(cov_input),
        res_input,
        ["Intercept", "condition_Treat_vs_Control"],
        BioToolkit.DifferentialExpression._beta_prior_var_log2_to_natural(0.25),
    )
    @test abs(shrunk_beta[1, 2]) < abs(beta_input[1, 2])
    @test shrunk_cov[1][2, 2] < cov_input[1][2, 2]
    @test abs(shrunk_res[1].log2_fold_change) < abs(res_input[1].log2_fold_change)

    cook_design = [:Control, :Control, :Control, :Treat, :Treat, :Treat]
    cook_design_matrix = [1.0 0.0; 1.0 0.0; 1.0 0.0; 1.0 1.0; 1.0 1.0; 1.0 1.0]
    dummy_counts = BioToolkit.CountMatrix(sparse(zeros(Int, 3, 6)), ["g1", "g2", "g3"], ["s1", "s2", "s3", "s4", "s5", "s6"])
    dds = BioToolkit.DESeqDataSet(dummy_counts, DataFrame(condition=cook_design), cook_design)
    dds.test = :Wald
    dds.model_matrix = cook_design_matrix
    dds.wald_results = [
        BioToolkit.DEResult("g1", 10.0, 0.0, 1.0, 0.0, 0.9, NaN, false, true),
        BioToolkit.DEResult("g2", 10.0, 0.0, 1.0, 0.0, 0.04, NaN, false, true),
        BioToolkit.DEResult("g3", 10.0, 0.0, 1.0, 0.0, 0.06, NaN, false, true),
    ]
    dds.cooks = [10.0 10.0 10.0 10.0 10.0 10.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
    dds.metadata["cooksCutoff"] = 1.0

    tab = BioToolkit.results(dds; independent_filter=false)
    @test ismissing(tab[1, :pvalue])
    @test ismissing(tab[1, :padj])
    @test tab[2, :padj] ≈ 0.06 atol=1e-12
    @test tab[3, :padj] ≈ 0.06 atol=1e-12

    # Regression: contrast results must honor independent filtering and keep
    # filtered genes at missing padj (DESeq2-compatible semantics).
    contrast_counts = vcat(fill(120, 5, 6), fill(1, 5, 6))
    contrast_cm = BioToolkit.CountMatrix(
        sparse(contrast_counts),
        ["g$(i)" for i in 1:10],
        ["s$(i)" for i in 1:6],
    )
    contrast_design = [:Control, :Control, :Control, :Treat, :Treat, :Treat]
    contrast_dds = BioToolkit.DESeqDataSet(
        contrast_cm,
        DataFrame(condition=contrast_design),
        contrast_design,
    )
    contrast_dds.test = :Wald
    contrast_dds.size_factors = ones(6)
    contrast_dds.coefficient_names = ["Intercept", "condition_Treat_vs_Control"]

    beta = zeros(10, 2)
    # Deterministic p-value pattern chosen to produce a non-trivial
    # independent-filter cutoff (> low-count baseMean).
    beta[:, 2] .= [1.47, 2.48, 1.66, 2.10, 2.56, 0.45, 0.35, 0.56, 0.52, 0.34]
    contrast_dds.wald_beta_matrix = beta
    contrast_dds.wald_covariances = [[1.0 0.0; 0.0 1.0] for _ in 1:10]

    contrast_tab = BioToolkit.results(
        contrast_dds,
        [:condition, :Treat, :Control];
        independent_filter=true,
    )

    cutoff = get(contrast_dds.metadata, "indepFilterCutoff", NaN)
    @test isfinite(cutoff)
    filtered_rows = findall(<=(cutoff), contrast_tab[!, :base_mean])
    @test !isempty(filtered_rows)
    @test all(ismissing(contrast_tab[row, :padj]) for row in filtered_rows)
end
