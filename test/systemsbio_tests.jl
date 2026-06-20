using DataFrames
using Graphs
using SparseArrays
using Statistics

@testset "SystemsBio" begin
    counts = BioToolkit.CountMatrix(
        sparse(Int[
            10 12 11  9;
             4  5  4  6;
            20 19 21 18;
             8  7  9  8;
        ]),
        ["g1", "g2", "g3", "g4"],
        ["s1", "s2", "s3", "s4"],
    )

    soft = BioToolkit.pick_soft_threshold(counts; powers=2:4)
    @test soft isa BioToolkit.SoftThresholdResult
    @test soft.best_power in soft.powers

    network = BioToolkit.find_modules(counts; power=2, tom_threshold=0.05, edge_threshold=0.05)
    @test network isa BioToolkit.GeneNetwork
    @test length(network.node_to_gene) == 4
    @test length(network.modules) == 4
    @test length(network.connectivity) == 4

    @testset "weighted TOM formula" begin
        adjacency = Float64[
            0.0 0.3 0.5 0.1;
            0.3 0.0 0.2 0.4;
            0.5 0.2 0.0 0.6;
            0.1 0.4 0.6 0.0;
        ]
        observed = BioToolkit.SystemsBio._tom_matrix(adjacency)
        shared = adjacency * adjacency
        degrees = vec(sum(adjacency, dims=2))
        expected = (shared .+ adjacency) ./ max.(min.(degrees, transpose(degrees)) .+ 1 .- adjacency, eps(Float64))
        for i in axes(expected, 1)
            expected[i, i] = 1.0
        end
        @test isapprox(observed, expected; atol=1e-12)
        @test all(observed[i, i] == 1.0 for i in axes(observed, 1))
    end

    eigengenes = BioToolkit.module_eigengenes(counts, network.modules)
    @test size(eigengenes, 1) == 4
    @test size(eigengenes, 2) == length(unique(network.modules))

    tree = BioToolkit.module_dendrogram(network)
    @test tree isa BioToolkit.PhyloTree
    @test !isempty(BioToolkit.get_terminals(tree))

    design = [1.0 0.0;
              1.0 0.0;
              1.0 1.0;
              1.0 1.0]
    limma = BioToolkit.limma_fit(counts, design; coefficient_names=["intercept", "group"], use_voom=false)
    @test limma isa BioToolkit.LimmaResult
    @test size(limma.coefficients, 1) == 4
    @test size(limma.coefficients, 2) == 2

    @testset "voom transform" begin
        voom = BioToolkit.voom_transform(counts, design)
        @test voom isa BioToolkit.VoomResult
        @test size(voom.log_cpm) == size(Matrix{Float64}(counts.counts))
        @test size(voom.weights) == size(Matrix{Float64}(counts.counts))
        @test all(isfinite, voom.weights)
        @test all(>(0), voom.weights)

        limma_voom = BioToolkit.limma_fit(counts, design; coefficient_names=["intercept", "group"], use_voom=true)
        @test all(isfinite, limma_voom.pvalues)
        @test all(p -> 0.0 <= p <= 1.0, limma_voom.pvalues)

        voom_qw = BioToolkit.voom_with_quality_weights(counts, design)
        @test voom_qw.voom isa BioToolkit.VoomResult
        @test size(voom_qw.voom.weights) == size(Matrix{Float64}(counts.counts))
        @test length(voom_qw.sample_weights) == size(counts.counts, 2)
        @test all(isfinite, voom_qw.sample_weights)
        @test all(>(0), voom_qw.sample_weights)
        @test isapprox(exp(mean(log.(voom_qw.sample_weights))), 1.0; atol=1e-10)
    end

    @testset "limma empirical Bayes" begin
        s2 = exp.(range(log(0.05), log(2.0); length=40))
        covariate = collect(range(1.0, 5.0; length=40))
        s0_sq, d0, trend = BioToolkit.SystemsBio.estimate_limma_hyperparameters(s2, covariate; return_trend=true, df_residual=4.0)
        @test s0_sq > 0
        @test d0 >= 0
        @test length(trend) == length(s2)
        @test all(isfinite, trend)
        @test all(>(0), trend)

        s2_const = fill(0.2, 20)
        cov_const = collect(1.0:20.0)
        _, d0_const, trend_const = BioToolkit.SystemsBio.estimate_limma_hyperparameters(s2_const, cov_const; return_trend=true, df_residual=4.0)
        @test isfinite(d0_const) || isinf(d0_const)
        @test all(isfinite, trend_const)
        @test all(>(0), trend_const)

        cov_negative = collect(range(-4.0, 2.0; length=40))
        s0_neg, d0_neg, trend_neg = BioToolkit.SystemsBio.estimate_limma_hyperparameters(s2, cov_negative; return_trend=true, df_residual=4.0)
        @test s0_neg > 0
        @test d0_neg >= 0
        @test length(trend_neg) == length(s2)
        @test all(isfinite, trend_neg)
        @test all(>(0), trend_neg)

        raw_counts = Matrix{Float64}(counts.counts)
        lib_sizes = vec(sum(raw_counts, dims=1))
        norm_factors = BioToolkit.calc_norm_factors(counts; method=:tmm)
        effective_lib_sizes = lib_sizes .* norm_factors
        pseudocount = 0.5
        log_cpm = log2.(raw_counts .+ pseudocount) .- log2.(reshape(effective_lib_sizes ./ 1e6 .+ pseudocount, 1, :))

        X = Matrix{Float64}(design)
        xtx_inv = BioToolkit.pinv(X' * X)
        response = permutedims(log_cpm)
        beta = xtx_inv * (X' * response)
        residuals = response .- X * beta
        df_residual = max(size(X, 1) - BioToolkit.rank(X), 1)
        s2_all = vec(sum(abs2, residuals; dims=1)) ./ df_residual
        gene_means = vec(BioToolkit.mean(log_cpm, dims=2))
        _, d0_fit, prior_trend = BioToolkit.SystemsBio.estimate_limma_hyperparameters(s2_all, gene_means; return_trend=true, df_residual=df_residual)

        expected_post = if isfinite(d0_fit)
            (df_residual .* s2_all .+ d0_fit .* prior_trend) ./ (df_residual + d0_fit)
        else
            prior_trend
        end
        @test isapprox(limma.moderated_variance, expected_post; atol=1e-10)

        df_total = min(df_residual + d0_fit, length(s2_all) * df_residual)
        scale = max(BioToolkit.diag(xtx_inv)[2], eps(Float64))
        se_group = sqrt.(max.(expected_post .* scale, eps(Float64)))
        t_group = limma.coefficients[:, 2] ./ se_group
        p_expected = 2 .* BioToolkit.ccdf.(BioToolkit.TDist(max(df_total, 1.0)), abs.(t_group))
        @test isapprox(limma.pvalues[:, 2], p_expected; atol=1e-10)

        raw_fit = (
            coefficients=copy(limma.coefficients),
            stdev_unscaled=limma.standard_errors ./ reshape(sqrt.(max.(limma.moderated_variance, eps(Float64))), :, 1),
            sigma2=copy(limma.moderated_variance),
            df_residual=2.0,
            gene_ids=copy(limma.gene_ids),
            coefficient_names=copy(limma.coefficient_names),
            base_mean=copy(limma.base_mean),
        )
        eb_small = BioToolkit.eBayes(raw_fit; robust=true, return_hyperparameters=true)
        @test eb_small.fit isa BioToolkit.LimmaResult
        @test !eb_small.robust_used
        @test all(isfinite, eb_small.fit.pvalues)

        expanded_raw = (
            coefficients=vcat(raw_fit.coefficients, raw_fit.coefficients, raw_fit.coefficients),
            stdev_unscaled=vcat(raw_fit.stdev_unscaled, raw_fit.stdev_unscaled, raw_fit.stdev_unscaled),
            sigma2=vcat(raw_fit.sigma2, raw_fit.sigma2, raw_fit.sigma2),
            df_residual=4.0,
            gene_ids=["g$(index)" for index in 1:12],
            coefficient_names=copy(raw_fit.coefficient_names),
            base_mean=vcat(raw_fit.base_mean, raw_fit.base_mean, raw_fit.base_mean),
        )
        eb_robust = BioToolkit.eBayes(expanded_raw; robust=true, return_hyperparameters=true)
        @test eb_robust.robust_used
        @test eb_robust.d0 >= 0
        @test all(isfinite, eb_robust.fit.t_statistics)
    end

    @testset "limma contrast, correlation, and batch APIs" begin
        contrast_matrix = [1.0 2.0; 0.0 0.0]
        contrast_fit = BioToolkit.contrasts_fit(limma, contrast_matrix; coefficient_names=["c1", "c2"])
        @test contrast_fit isa BioToolkit.LimmaResult
        @test size(contrast_fit.coefficients) == (size(limma.coefficients, 1), 2)
        finite_columns = [all(isfinite, contrast_fit.pvalues[:, index]) for index in 1:2]
        nan_columns = [all(isnan, contrast_fit.pvalues[:, index]) for index in 1:2]
        @test count(identity, finite_columns) == 1
        @test count(identity, nan_columns) == 1

        expression = Matrix{Float64}(counts.counts)
        block = [:b1, :b1, :b2, :b2]
        correlation = BioToolkit.duplicateCorrelation(expression, design, block)
        @test isfinite(correlation.correlation)
        @test abs(correlation.correlation) < 1.0

        block_singleton = [:b1, :b1, :b2, :s]
        correlation_singleton = BioToolkit.duplicateCorrelation(expression, design, block_singleton)
        @test isfinite(correlation_singleton.correlation)

        batch = [:A, :A, :B, :B]
        corrected = BioToolkit.remove_batch_effect_limma(expression, batch)
        @test size(corrected) == size(expression)
        mean_a = vec(mean(corrected[:, 1:2], dims=2))
        mean_b = vec(mean(corrected[:, 3:4], dims=2))
        @test maximum(abs.(mean_a .- mean_b)) < 1e-8

        covariates = [0.0, 0.0, 1.0, 1.0]
        corrected_cov = BioToolkit.remove_batch_effect_limma(expression, batch; covariates=covariates)
        @test all(isfinite, corrected_cov)
    end

    deresults = BioToolkit.limma_deresults(limma; coefficient_index=2)
    @test length(deresults) == 4
    @test all(result -> result isa BioToolkit.DEResult, deresults)

    db = BioToolkit.build_annotation_database([
        BioToolkit.EnrichmentTerm("TERM:1", "toy pathway", "GO", ["g1", "g2"], String[]),
        BioToolkit.EnrichmentTerm("TERM:2", "background", "GO", ["g3", "g4"], String[]),
    ])
    gsea_results = BioToolkit.gsea(deresults, db; permutations=25, min_size=1)
    @test all(result -> result isa BioToolkit.EnrichmentResult, gsea_results)
    @test !isempty(gsea_results)

    dag = BioToolkit.infer_network(permutedims(Matrix{Float64}(counts.counts)); gene_ids=counts.gene_ids)
    @test dag isa BioToolkit.NetworkInferenceResult
    @test length(dag.node_to_gene) == 4

    constrained = BioToolkit.infer_network(permutedims(Matrix{Float64}(counts.counts)); gene_ids=counts.gene_ids, max_parents=1, max_iterations=20)
    @test !is_cyclic(constrained.graph)
    @test all(indegree(constrained.graph, node) <= 1 for node in 1:nv(constrained.graph))

    assays = [
        [1.0 2.0 3.0; 2.0 3.0 4.0; 3.0 4.0 5.0; 4.0 5.0 6.0],
        [2.0 1.0 0.0; 3.0 2.0 1.0; 4.0 3.0 2.0; 5.0 4.0 3.0],
    ]
    mofa = BioToolkit.multi_omics_factor_analysis(assays; n_factors=2)
    @test mofa isa BioToolkit.MultiOmicsFactorAnalysisResult
    @test size(mofa.factors, 1) == 4
    @test size(mofa.factors, 2) == 2
    @test length(mofa.assay_loadings) == 2
end
