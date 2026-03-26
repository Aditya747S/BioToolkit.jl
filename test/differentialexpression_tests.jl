using SparseArrays
using LinearAlgebra
using Random

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
    single_result = BioToolkit.differential_expression(single_gene, [:Control, :Control, :Treat, :Treat])
    @test single_result[1].pvalue == 1.0
    @test single_result[1].padj == 1.0
end