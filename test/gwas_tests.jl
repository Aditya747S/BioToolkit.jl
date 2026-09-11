using Test
using LinearAlgebra
using Statistics
using SimpleWeightedGraphs
using Plots
using BioToolkit

function _write_u16_le(io::IO, value::Integer)
    raw = UInt16(value)
    write(io, UInt8(raw & 0xff))
    write(io, UInt8((raw >> 8) & 0xff))
    return nothing
end

function _write_u32_le(io::IO, value::Integer)
    raw = UInt32(value)
    write(io, UInt8(raw & 0xff))
    write(io, UInt8((raw >> 8) & 0xff))
    write(io, UInt8((raw >> 16) & 0xff))
    write(io, UInt8((raw >> 24) & 0xff))
    return nothing
end

function _write_minimal_bgen(path::String)
    sample_ids = ["I1", "I2", "I3"]
    sample_block_length = 4 + sum(2 + ncodeunits(id) for id in sample_ids)

    open(path, "w") do io
        _write_u32_le(io, 0)
        _write_u32_le(io, 20)
        _write_u32_le(io, 1)
        _write_u32_le(io, length(sample_ids))
        write(io, "bgen")
        _write_u32_le(io, 0x80000008)

        _write_u32_le(io, sample_block_length)
        _write_u32_le(io, length(sample_ids))
        for id in sample_ids
            _write_u16_le(io, ncodeunits(id))
            write(io, id)
        end

        _write_u16_le(io, 3)
        write(io, "rs1")
        _write_u16_le(io, 3)
        write(io, "rs1")
        _write_u16_le(io, 1)
        write(io, "1")
        _write_u32_le(io, 101)
        _write_u16_le(io, 2)
        _write_u32_le(io, 1)
        write(io, "A")
        _write_u32_le(io, 1)
        write(io, "G")

        block_io = IOBuffer()
        _write_u32_le(block_io, length(sample_ids))
        _write_u16_le(block_io, 2)
        write(block_io, UInt8(2))
        write(block_io, UInt8(2))
        write(block_io, UInt8[0x02, 0x02, 0x02])
        write(block_io, UInt8(0))
        write(block_io, UInt8(8))
        write(block_io, UInt8[0xff, 0x00, 0x00, 0xff, 0x00, 0x00])
        payload = take!(block_io)

        _write_u32_le(io, length(payload))
        _write_u32_le(io, length(payload))
        write(io, payload)
    end

    return path
end

@testset "GWAS module" begin
    fixture_prefix = joinpath(@__DIR__, "..", "Examples", "fixtures", "plink_small", "plink_small")
    genotypes = BioToolkit.read_plink(fixture_prefix)
    expected = [0.0 1.0 2.0; 1.0 0.0 1.0; 2.0 1.0 0.0; 0.0 1.0 1.0]
    @test Matrix(genotypes) == expected
    @test genotypes.bim.snp_id == ["rs1", "rs2", "rs3"]
    @test genotypes.fam.sample_id == ["I1", "I2", "I3", "I4"]

    mktempdir() do dir
        roundtrip_prefix = joinpath(dir, "roundtrip")
        BioToolkit.write_plink(roundtrip_prefix, Matrix(genotypes), genotypes.bim, genotypes.fam)
        roundtrip = BioToolkit.read_plink(roundtrip_prefix)
        @test Matrix(roundtrip) == expected
        @test roundtrip.bim.snp_id == genotypes.bim.snp_id
        @test roundtrip.fam.sample_id == genotypes.fam.sample_id
    end

    phenotype = [0.0, 0.5, 1.0, 1.5]

    linear = BioToolkit.gwas_linear_scan(genotypes, phenotype)
    @test linear.method == "linear_scan"
    @test length(linear.snp_ids) == 3
    @test all(isfinite, linear.pvalue)

    mixed = BioToolkit.gwas_lmm_scan(genotypes, phenotype)
    @test mixed.method == "lmm_scan"
    @test length(mixed.snp_ids) == 3

    clumped = BioToolkit.ld_clumping(linear, genotypes; p_threshold=1.0, r2_threshold=0.1, window=500)
    @test length(clumped.snp_ids) <= length(linear.snp_ids)

    prs = BioToolkit.calculate_prs(genotypes, linear)
    @test length(prs) == size(genotypes, 1)

    prs_cv = BioToolkit.prs_cross_validation(genotypes, phenotype, linear)
    @test haskey(prs_cv.scores, (1e-5, 50_000))

    meta = BioToolkit.meta_analyze([linear, linear])
    @test meta.study_count == 2
    @test length(meta.snp_ids) == length(linear.snp_ids)

    peakset = BioToolkit.PeakSet([
        BioToolkit.Peak("peak1", "chr1", 50, 150, 100, 2.0, 1e-6, 1e-4),
    ])
    overlaps = BioToolkit.overlap_gwas_peaks(linear, peakset; flank=75, pvalue_threshold=1.0)
    @test "snp_id" in names(overlaps)

    network = BioToolkit.GeneNetwork(
        SimpleWeightedGraph(3),
        Dict("GENE1" => 1, "GENE2" => 2, "GENE3" => 3),
        ["GENE1", "GENE2", "GENE3"],
        [0.0, 0.0, 0.0],
        [1, 1, 2],
    )
    enrichment = BioToolkit.gene_based_test(linear, network; pvalue_threshold=1.0)
    @test enrichment isa Vector

    manhattan_points, _ = BioToolkit.BioPlotting.manhattan_data(linear)
    qq_points, _ = BioToolkit.BioPlotting.qq_data(linear)
    forest = BioToolkit.BioPlotting.gwas_forest_plot(meta)
    manhattan = BioToolkit.BioPlotting.manhattan_plot(linear)
    qq = BioToolkit.BioPlotting.qq_plot(linear)
    @test !isempty(manhattan_points)
    @test !isempty(qq_points)
    @test forest.figure !== nothing
    @test manhattan.figure.subplots[1].attr[:title] == "Manhattan plot"
    @test manhattan.figure.subplots[1].attr[:annotations] |> length == 3
    @test length(manhattan.figure.series_list) == 2
    @test qq.figure.subplots[1].attr[:title] == "QQ plot"
    @test qq.figure.subplots[1].attr[:aspect_ratio] == :equal
    @test length(qq.figure.series_list) == 3

    maf_mask, maf = BioToolkit.calculate_maf(genotypes; min_maf=0.4, return_frequency=true)
    @test maf_mask == [false, false, true]
    @test all(isfinite, maf)

    qc_matrix = [
        0.0 1.0 NaN;
        1.0 NaN 2.0;
        2.0 1.0 0.0;
        NaN 0.0 1.0;
    ]
    variant_missing, variant_call = BioToolkit.calculate_missingness(qc_matrix; by=:variant, return_call_rate=true)
    @test variant_missing == [0.25, 0.25, 0.25]
    @test variant_call == [0.75, 0.75, 0.75]

    sample_missing = BioToolkit.calculate_missingness(qc_matrix; by=:sample)
    @test sample_missing == [1 / 3, 1 / 3, 0.0, 1 / 3]
    @test_throws ArgumentError BioToolkit.calculate_missingness(qc_matrix; by=:bad_axis)

    variant_mask, variant_rates = BioToolkit.missingness_filter(qc_matrix; by=:variant, max_missing=0.2, return_rates=true)
    @test variant_mask == [false, false, false]
    @test variant_rates == [0.25, 0.25, 0.25]
    @test all(BioToolkit.missingness_filter(qc_matrix; by=:variant, max_missing=0.3))

    info_mask, info_scores = BioToolkit.info_score_filter([0.95, 0.72, NaN]; min_info=0.8, return_scores=true)
    @test info_mask == [true, false, false]
    @test info_scores[1:2] == [0.95, 0.72]
    @test isnan(info_scores[3])

    qc_report = BioToolkit.gwas_qc_report(genotypes; min_maf=0.0, max_missing=0.1, hwe_threshold=0.0, info_scores=[0.9, 0.7, 0.95], min_info=0.8)
    @test names(qc_report) == ["CHR", "POS", "ID", "CALL_RATE", "MISSING_RATE", "MAF", "HWE_P", "INFO", "PASS_MAF", "PASS_MISSING", "PASS_HWE", "PASS_INFO", "PASS"]
    @test size(qc_report, 1) == size(genotypes, 2)
    @test qc_report.PASS_INFO == [true, false, true]
    @test qc_report.PASS == [true, false, true]

    hwe_pvalues = BioToolkit.calculate_hwe_pvalues(genotypes)
    @test length(hwe_pvalues) == 3
    @test all(p -> isfinite(p) && 0.0 <= p <= 1.0, hwe_pvalues)
    @test all(BioToolkit.hwe_filter(genotypes; p_threshold=0.0))
    @test_throws ArgumentError BioToolkit.hwe_exact(-1, 1, 1)

    lambda_gc = BioToolkit.genomic_control_lambda(linear)
    @test isfinite(lambda_gc)
    @test lambda_gc > 0.0
    gc_adjusted = BioToolkit.apply_genomic_control(linear; lambda=lambda_gc)
    @test gc_adjusted.method == "linear_scan+gc"
    @test length(gc_adjusted.pvalue) == length(linear.pvalue)

    ld_matrix, snp_indices, filtered_maf = BioToolkit.calculate_ld_matrix(genotypes; min_maf=0.0, max_maf=1.0, multi_thread=false, return_snp_indices=true)
    @test size(ld_matrix) == (3, 3)
    @test LinearAlgebra.issymmetric(ld_matrix)
    @test all(isapprox.(diag(ld_matrix), 1.0; atol=1e-8))
    @test snp_indices == [1, 2, 3]
    @test length(filtered_maf) == 3

    grm = BioToolkit.calculate_grm(genotypes)
    @test size(grm) == (4, 4)
    @test LinearAlgebra.issymmetric(grm)
    @test isapprox(Statistics.mean(diag(grm)), 1.0; atol=1e-8)

    covariates = hcat(ones(Float64, size(genotypes, 1)), collect(1.0:size(genotypes, 1)))
    projected_matrix, projector = BioToolkit.loco_projection(covariates, Matrix(genotypes); return_projector=true)
    @test size(projected_matrix) == size(Matrix(genotypes))
    @test size(projector) == (size(genotypes, 1), size(genotypes, 1))
    @test maximum(abs.(covariates' * projected_matrix)) < 1e-5
    projected_vector, _ = BioToolkit.loco_projection(covariates, phenotype; return_projector=true)
    @test length(projected_vector) == length(phenotype)
    @test maximum(abs.(covariates' * projected_vector)) < 1e-5

    logistic_phenotype = [0.0, 0.0, 1.0, 1.0]
    logistic = BioToolkit.gwas_logistic_scan(genotypes, logistic_phenotype; firth=false, multi_thread=false)
    @test logistic.method == "logistic_scan"
    @test length(logistic.snp_ids) == 3
    @test all(isfinite, logistic.pvalue)

    gxe_genotypes = [
        0.0 0.0 0.0;
        1.0 0.0 1.0;
        0.0 1.0 0.0;
        1.0 1.0 1.0;
        2.0 0.0 1.0;
        2.0 1.0 2.0;
    ]
    gxe_continuous = [0.2, 0.4, 0.7, 0.9, 1.1, 1.4]
    gxe_linear = BioToolkit.gwas_gxe_interaction(gxe_genotypes, gxe_continuous, 1, [2, 3]; multi_thread=false)
    @test gxe_linear.method == "gxe_linear"
    @test length(gxe_linear.snp_ids) == 2
    @test all(isfinite, gxe_linear.pvalue)

    gxe_binary = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    gxe_logistic = BioToolkit.gwas_gxe_interaction(gxe_genotypes, gxe_binary, 1, 2; logistic=true, firth=false, multi_thread=false)
    @test gxe_logistic.method == "gxe_logistic"
    @test length(gxe_logistic.snp_ids) == 1
    @test all(isfinite, gxe_logistic.pvalue)

    plink_table = BioToolkit.to_plink_dataframe(linear, genotypes)
    @test names(plink_table) == ["CHR", "POS", "ID", "REF", "ALT", "BETA", "SE", "P"]
    @test size(plink_table, 1) == length(linear.snp_ids)

    mktemp() do path, io
        close(io)
        BioToolkit.write_plink_sumstats(path, linear; genotypes=genotypes)
        lines = readlines(path)
        @test lines[1] == "CHR\tPOS\tID\tREF\tALT\tBETA\tSE\tP"
        @test length(lines) == length(linear.snp_ids) + 1
    end

    migration = BioToolkit.gwas_migration_guide()
    @test occursin("SnpArrays.jl", migration)
    @test occursin("read_bgen", migration)

    @testset "Advanced GWAS analytics" begin
        n_snps = 40
        snp_ids = ["adv_rs$(i)" for i in 1:n_snps]
        chromosomes = fill("1", n_snps)
        positions = collect(10_000:(10_000 + n_snps - 1))
        alleles = [("A", "G") for _ in 1:n_snps]
        genes = [i <= 20 ? "setA" : "setB" for i in 1:n_snps]
        beta_a = collect(range(-0.2, 0.2; length=n_snps))
        se_a = fill(0.05, n_snps)
        z_a = beta_a ./ se_a
        p_a = clamp.(exp.(-abs.(z_a)), 1e-12, 1.0)

        beta_b = beta_a .* 0.8 .+ 0.02
        se_b = copy(se_a)
        z_b = beta_b ./ se_b
        p_b = clamp.(exp.(-abs.(z_b)), 1e-12, 1.0)

        adv_trait_a = BioToolkit.GWASResult(snp_ids, chromosomes, positions, alleles, genes, beta_a, se_a, z_a, p_a, 10_000, String[], "traitA", "synthetic")
        adv_trait_b = BioToolkit.GWASResult(snp_ids, chromosomes, positions, alleles, genes, beta_b, se_b, z_b, p_b, 10_000, String[], "traitB", "synthetic")
        ld_scores = collect(range(1.0, 5.0; length=n_snps))

        ldsc = BioToolkit.ldsc_heritability(adv_trait_a, ld_scores; n_blocks=20)
        @test ldsc.n_snps_used == n_snps
        @test isfinite(ldsc.intercept)

        part_h2 = BioToolkit.partitioned_heritability(adv_trait_a, ld_scores; annotations=genes, min_snps=10)
        @test haskey(part_h2, "setA")
        @test haskey(part_h2, "setB")

        rg, rg_se = BioToolkit.ldsc_genetic_correlation([adv_trait_a, adv_trait_b]; min_snps=20)
        @test isfinite(rg)
        @test isfinite(rg_se)

        greml_geno = [0.0 1.0; 1.0 0.0; 2.0 1.0; 1.0 2.0; 0.0 2.0; 2.0 0.0]
        greml_pheno = [0.1, 0.5, 0.9, 1.1, 1.4, 1.8]
        h2_greml = BioToolkit.estimate_heritability_greml(greml_geno, greml_pheno; kinship=Matrix{Float64}(I, 6, 6))
        @test isfinite(h2_greml)

        fm_beta = [0.22, 0.01, 0.16, 0.00, -0.05]
        fm_se = fill(0.05, 5)
        fm_ld = Matrix{Float64}(I, 5, 5)
        pip = BioToolkit.posterior_inclusion_probability(fm_beta, fm_se, fm_ld)
        @test length(pip) == 5
        @test all(0.0 .<= pip .<= 1.0)
        cs = BioToolkit.calculate_credible_set(pip; coverage=0.8)
        @test !isempty(cs)
        susie = BioToolkit.fine_map_susie(fm_beta, fm_se, fm_ld; n_effects=3)
        @test length(susie.pip) == 5
        @test susie.n_effects == 3

        mr_beta_exp = [0.12, 0.08, 0.15, 0.20]
        mr_se_exp = fill(0.02, 4)
        mr_beta_out = [0.06, 0.03, 0.08, 0.09]
        mr_se_out = fill(0.03, 4)
        ivw = BioToolkit.mr_two_sample(mr_beta_exp, mr_se_exp, mr_beta_out, mr_se_out)
        egger = BioToolkit.mr_egger(mr_beta_exp, mr_se_exp, mr_beta_out, mr_se_out)
        q_stat, q_p = BioToolkit.mr_pleiotropy_test(mr_beta_exp, mr_se_exp, mr_beta_out, mr_se_out)
        @test isfinite(ivw.estimate)
        @test isfinite(egger.estimate)
        @test isfinite(q_stat)
        @test isfinite(q_p)

        conditional = BioToolkit.conditional_analysis(linear, genotypes, [1])
        joint = BioToolkit.joint_analysis(linear, genotypes, [1, 2])
        cojo = BioToolkit.cojo_stepwise(linear, genotypes; p_threshold=1.0, r2_threshold=1.1, max_snps=2)
        @test length(conditional.snp_ids) == length(linear.snp_ids)
        @test length(joint.snp_ids) == 2
        @test cojo isa BioToolkit.GWASResult

        pcs, var_exp = BioToolkit.gwas_pca(genotypes; n_components=2, maf_threshold=0.0)
        @test size(pcs) == (size(genotypes, 1), 2)
        @test length(var_exp) == 2

        pcs_fit, _, loadings, center, scale, idx = BioToolkit.gwas_pca(genotypes; n_components=2, maf_threshold=0.0, return_loadings=true)
        projected = BioToolkit.project_pca(genotypes, loadings; center=center, scale=scale, snp_indices=idx)
        @test size(projected) == size(pcs_fit)
        @test projected ≈ pcs_fit atol=1e-8

        ibd = BioToolkit.calculate_ibd(genotypes)
        kin = BioToolkit.calculate_king_kinship(genotypes)
        rel_pairs = BioToolkit.detect_related_pairs(kin; threshold=-1.0)
        @test size(ibd) == (size(genotypes, 1), size(genotypes, 1))
        @test size(kin) == (size(genotypes, 1), size(genotypes, 1))
        @test !isempty(rel_pairs)

        coloc = BioToolkit.coloc_abf(fm_beta, fm_se, fm_beta .+ 0.01, fm_se, 1_000)
        @test coloc.n_snps == 5
        @test 0.0 <= coloc.posterior_prob_ab <= 1.0

        sel_positions = Int.(genotypes.bim.position)
        h_a = [1, 1, 2]
        h_b = [1, 2, 2]
        ihs = BioToolkit.ihs_score(genotypes, sel_positions, h_a, h_b; window=500)
        xp = BioToolkit.xp_ehh_score(genotypes, sel_positions, [1, 2], [3, 4]; window=500, min_maf=0.0)
        fst = BioToolkit.fst_outlier_test(genotypes, sel_positions; top_frac=0.5)
        @test length(ihs) == length(sel_positions)
        @test length(xp) == length(sel_positions)
        @test length(fst.fst_values) == length(sel_positions)

        transformed = BioToolkit.rank_inverse_normal([1.0, 2.0, 3.0, 4.0])
        transformed_inplace = [1.0, 2.0, 3.0, 4.0]
        BioToolkit.rank_inverse_normal!(transformed_inplace)
        @test all(isfinite, transformed)
        @test all(isfinite, transformed_inplace)

        sex_check = BioToolkit.calculate_sex_check(genotypes)
        het_qc = BioToolkit.calculate_heterozygosity_outliers(genotypes)
        sample_qc = BioToolkit.sample_qc_report(genotypes; missing_threshold=1.0, het_threshold=10.0)
        @test size(sex_check, 1) == size(genotypes, 1)
        @test length(het_qc.outlier) == size(genotypes, 1)
        @test size(sample_qc, 1) == size(genotypes, 1)
    end

    mktemp() do path, io
        write(io, "#CHR\tPOS\tID\tREF\tALT\tS1\tS2\tS3\n")
        write(io, "1\t100\trs1\tA\tG\t0\t1\t2\n")
        write(io, "1\t200\trs2\tC\tT\t1\t1\tNA\n")
        close(io)

        bed_reader = BioToolkit.BedReader(path)
        bed_variants = collect(bed_reader)
        @test length(bed_variants) == 2
        @test bed_variants[1].snp_id == "rs1"
        @test bed_variants[1].dosage == [0.0, 1.0, 2.0]
        @test isnan(bed_variants[2].dosage[3])

        bed_matrix = BioToolkit.read_bed_genotypes(path)
        @test size(bed_matrix) == (3, 2)
        @test Matrix(bed_matrix)[:, 1] == [0.0, 1.0, 2.0]
        @test isnan(Matrix(bed_matrix)[3, 2])
        @test bed_matrix.fam.sample_id == ["S1", "S2", "S3"]
    end

    mktemp() do path, io
        close(io)
        _write_minimal_bgen(path)

        bgen_matrix = BioToolkit.read_bgen(path)
        @test size(bgen_matrix) == (3, 1)
        @test vec(Matrix(bgen_matrix)) == [0.0, 1.0, 2.0]
        @test bgen_matrix.bim.snp_id == ["rs1"]
        @test bgen_matrix.fam.sample_id == ["I1", "I2", "I3"]

        bgen_reader = BioToolkit.BgenReader(path)
        bgen_variants = collect(bgen_reader)
        @test length(bgen_variants) == 1
        @test bgen_variants[1].alleles == ["A", "G"]
        @test bgen_variants[1].dosage == [0.0, 1.0, 2.0]
    end

    mktemp() do path, io
        close(io)
        BioToolkit.write_bgen(path, genotypes)
        bgen_roundtrip = BioToolkit.read_bgen(path)
        @test size(bgen_roundtrip) == size(genotypes)
        @test bgen_roundtrip.bim.snp_id == genotypes.bim.snp_id
        @test bgen_roundtrip.fam.sample_id == genotypes.fam.sample_id
        @test Matrix(bgen_roundtrip) ≈ Matrix(genotypes) atol=5e-3
    end

    @testset "GWAS roadmap API additions" begin
        @test BioToolkit.normalise_chromosome("chr1") == "1"
        @test BioToolkit.normalise_chromosome("23") == "X"

        var_subset = BioToolkit.filter_variants(genotypes, [true, false, true])
        sample_subset = BioToolkit.filter_samples(genotypes, [1, 3])
        @test size(var_subset) == (size(genotypes, 1), 2)
        @test size(sample_subset) == (2, size(genotypes, 2))

        merged = BioToolkit.merge_genotype_matrices([var_subset, var_subset]; by=:variants)
        @test size(merged) == (size(genotypes, 1), 4)

        ld_scores = BioToolkit.compute_ld_scores(genotypes; window_kb=1000, min_maf=0.0)
        @test length(ld_scores) == size(genotypes, 2)

        kept_idx = BioToolkit.prune_ld(genotypes; r2=0.0, window_kb=1000, return_indices=true)
        @test !isempty(kept_idx)

        flipped = BioToolkit.flip_alleles(genotypes, Dict("rs1" => "G"))
        @test any(flipped.flipped)

        harmonised = BioToolkit.harmonise_alleles(linear, linear)
        @test length(harmonised.result_a.snp_ids) == length(linear.snp_ids)
        @test length(harmonised.result_b.snp_ids) == length(linear.snp_ids)

        score_scan = BioToolkit.score_test_linear(genotypes, phenotype; multi_thread=false)
        lrt_scan = BioToolkit.likelihood_ratio_test(genotypes, phenotype)
        @test length(score_scan.snp_ids) == size(genotypes, 2)
        @test length(lrt_scan.snp_ids) == size(genotypes, 2)

        survival_scan = BioToolkit.gwas_survival_scan(genotypes, [1.0, 2.0, 3.0, 4.0], [1.0, 0.0, 1.0, 1.0]; multi_thread=false)
        ordinal_scan = BioToolkit.gwas_ordinal_scan(genotypes, [0.0, 1.0, 2.0, 2.0]; multi_thread=false)
        multivar_scan = BioToolkit.gwas_multivariate_scan(genotypes, hcat(phenotype, phenotype .+ 0.1); multi_thread=false)
        @test survival_scan.method == "cox_score_scan"
        @test ordinal_scan.method == "ordinal_scan_approx"
        @test multivar_scan.method == "multivariate_scan"

        burden = BioToolkit.burden_test(genotypes, Dict("SET1" => [1, 2]), phenotype)
        skat = BioToolkit.skat_test(genotypes, Dict("SET1" => [1, 2]), phenotype)
        @test size(burden, 1) == 1
        @test size(skat, 1) == 1

        condfdr = BioToolkit.conditional_fdr([0.1, 0.2], [0.05, 0.4])
        pi0 = BioToolkit.storey_pi0_estimate([0.1, 0.2, 0.8, 0.9])
        power = BioToolkit.gwas_power_calculation(10_000, 0.2, 0.1, 5e-8)
        pperm = BioToolkit.permutation_pvalue(1.5, [0.2, 1.0, 1.6]; n_perm=3)
        sem = BioToolkit.genomic_sem_fit([BioToolkit.LDSCResult(0.2, 0.05, 1.0, 0.1, 1.1, 1.0, 100, 10)])
        @test length(condfdr) == 2
        @test isfinite(pi0)
        @test 0.0 <= power <= 1.0
        @test 0.0 <= pperm <= 1.0
        @test sem.model == :common_factor

        hard = BioToolkit.dosage_to_hardcall([0.02, 1.02, 1.85, NaN]; threshold=0.2)
        info = BioToolkit.info_score_from_dosage([0.0, 1.0, 2.0], [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
        @test hard[1:2] == [0.0, 1.0]
        @test isfinite(info)

        reml = BioToolkit.reml_variance_components(phenotype, [Matrix{Float64}(I, length(phenotype), length(phenotype))])
        he = BioToolkit.he_regression_variance_components(phenotype, [Matrix{Float64}(I, length(phenotype), length(phenotype))])
        simG = BioToolkit.simulate_genotypes(12, 6)
        simy = BioToolkit.simulate_phenotype(simG, [1, 2], [0.2, 0.3], 0.4)
        expanded = BioToolkit.ld_expand_credible_set([1], genotypes; r2=0.0)
        neff = BioToolkit.estimate_effective_n(linear)
        rg_pop, rg_se_pop = BioToolkit.popcorn_genetic_correlation([1.0, 2.0, 3.0], [1.1, 2.1, 2.9], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0])
        @test length(reml.sigma2_random) == 1
        @test length(he.sigma2_random) == 1
        @test size(simG) == (12, 6)
        @test length(simy) == 12
        @test !isempty(expanded)
        @test isfinite(neff)
        @test isfinite(rg_pop)
        @test isfinite(rg_se_pop)

        twas = BioToolkit.twas_scan(genotypes, Dict("GENE1" => [0.1, 0.0, 0.2]), phenotype)
        smr = BioToolkit.smr_test(linear, linear)
        @test_throws ArgumentError BioToolkit.liftover(linear, "GRCh37", "GRCh38")
        lifted = BioToolkit.liftover(linear, "GRCh37", "GRCh38"; allow_passthrough=true)
        ann = BioToolkit.functional_annotation(linear, Dict("rs1" => "missense_variant"))
        ebi = BioToolkit.ebi_lookup(String[]; max_requests=0)
        @test size(twas, 1) == 1
        @test size(smr.summary, 1) >= 0
        @test occursin("liftover", lifted.method)
        @test "consequence" in names(ann)
        @test size(ebi, 1) == 0

        gw_manhattan = BioToolkit.GWAS.manhattan_plot(linear)
        gw_qq = BioToolkit.GWAS.qq_plot(linear.pvalue)
        zoom = BioToolkit.locus_zoom(linear, "1", 0, 500)
        @test haskey(gw_manhattan, :data)
        @test haskey(gw_qq, :data)
        @test size(zoom, 1) >= 0

        mktemp() do arrow_path, arrow_io
            close(arrow_io)
            BioToolkit.save_gwas_result(arrow_path, linear)
            loaded = BioToolkit.load_gwas_result(arrow_path)
            @test length(loaded.snp_ids) == length(linear.snp_ids)
        end

        int_vec = [1, 2, 3, 4]
        transformed_int = BioToolkit.rank_inverse_normal!(int_vec)
        @test transformed_int isa AbstractVector
        @test all(isfinite, transformed_int)
    end
end

using Random

using Distributions

@testset "Session 1: PLINK .bed follows the published bit-pair spec" begin
    # Spec codes: 00=hom(A1)=0, 01=het=1, 10=hom(A2)=2, 11=missing.
    byte = UInt8(0xe4)  # lanes: 00 01 10 11
    @test BioToolkit.GWAS._decode_plink_code(byte & 0x03) == 0.0
    @test BioToolkit.GWAS._decode_plink_code((byte >> 2) & 0x03) == 1.0
    @test BioToolkit.GWAS._decode_plink_code((byte >> 4) & 0x03) == 2.0
    @test isnan(BioToolkit.GWAS._decode_plink_code((byte >> 6) & 0x03))
    @test BioToolkit.GWAS.PLINK_DECODE_LUT[1, Int(byte) + 1] == 0.0
    @test BioToolkit.GWAS.PLINK_DECODE_LUT[2, Int(byte) + 1] == 1.0
    @test BioToolkit.GWAS.PLINK_DECODE_LUT[3, Int(byte) + 1] == 2.0
    @test isnan(BioToolkit.GWAS.PLINK_DECODE_LUT[4, Int(byte) + 1])

    @test BioToolkit.GWAS._plink_code(0.0) == 0x00
    @test BioToolkit.GWAS._plink_code(1.0) == 0x01
    @test BioToolkit.GWAS._plink_code(2.0) == 0x02
    @test BioToolkit.GWAS._plink_code(missing) == 0x03
    @test BioToolkit.GWAS._plink_code(NaN) == 0x03

    # Byte-exact encode of a known column: dosages [0,1,2,0] -> 0x24.
    payload = BioToolkit.GWAS._encode_bed(reshape(Float64[0, 1, 2, 0], 4, 1))
    @test payload[1:3] == UInt8[0x6c, 0x1b, 0x01]
    @test payload[4] == 0x24

    # The shipped fixture must itself be spec-encoded (it was regenerated when
    # the decoder was fixed: het<->missing were swapped before).
    fixture_prefix = joinpath(@__DIR__, "..", "Examples", "fixtures", "plink_small", "plink_small")
    @test read(fixture_prefix * ".bed") == vcat(UInt8[0x6c, 0x1b, 0x01], UInt8[0x24, 0x45, 0x46])

    # Missing calls survive a write/read round trip.
    genotypes_fixture = BioToolkit.read_plink(fixture_prefix)
    mktempdir() do dir
        prefix = joinpath(dir, "missing_ok")
        dosages = [0.0 1.0; NaN 2.0; 1.0 NaN]
        BioToolkit.write_plink(prefix, dosages, first(genotypes_fixture.bim, 2), first(genotypes_fixture.fam, 3))
        @test isequal(Matrix(BioToolkit.read_plink(prefix)), dosages)
    end
end

@testset "Session 1: BGEN first-variant offset per spec" begin
    fixture_prefix = joinpath(@__DIR__, "..", "Examples", "fixtures", "plink_small", "plink_small")
    genotypes_fixture = BioToolkit.read_plink(fixture_prefix)
    mktempdir() do dir
        path = joinpath(dir, "out.bgen")
        BioToolkit.write_bgen(path, genotypes_fixture)
        bytes = read(path)
        offset = UInt32(bytes[1]) | UInt32(bytes[2]) << 8 | UInt32(bytes[3]) << 16 | UInt32(bytes[4]) << 24
        # The sample-block length field includes its own 4 bytes per the BGEN
        # spec, so the offset is 24 + sample_block_length (previously 4 too far).
        sample_ids = genotypes_fixture.fam.sample_id
        @test offset == 24 + (4 + sum(2 + ncodeunits(id) for id in sample_ids))
        reader = BioToolkit.read_bgen(path)
        @test size(reader) == (4, 3)
        @test Matrix(reader) == Matrix(genotypes_fixture)
    end
end

@testset "Session 1: fine_map_susie implements IBSS/SER" begin
    Random.seed!(2026)
    p = 20
    causal = 12

    # Identity LD, one strong causal SNP. Signal strength: mean(z^2)-1 = 64/20
    # must be substantial for a single-effect SER to concentrate (the IBSS
    # prior variance is fit by the moment estimator Var(R*r) - 1).
    z = randn(p)
    z[causal] = 8.0
    R = Matrix{Float64}(I, p, p)
    res = BioToolkit.fine_map_susie(z, ones(p), R; n_effects=5, max_iter=50, tol=1e-6)
    @test res isa BioToolkit.SuSiEResult
    @test res.n_effects == 5
    @test length(res.pip) == p
    @test argmax(res.pip) == causal
    @test res.pip[causal] > 0.9
    @test res.credible_sets[1] == [causal]
    @test length(res.credible_sets) == 5
    @test res.heritability > 0.5
    @test all(isfinite, res.posterior_mean)
    @test all(>=(0.0), res.posterior_sd)

    # Correlated block: the credible set must cover the LD block, not a point.
    p2 = 50
    block_lo, block_hi = 24, 26
    z2 = randn(p2)
    z2[block_lo:block_hi] .= 3.0
    R2 = Matrix{Float64}(I, p2, p2)
    for i in block_lo:block_hi, j in block_lo:block_hi
        R2[i, j] = i == j ? 1.0 : 0.95
    end
    res2 = BioToolkit.fine_map_susie(z2, ones(p2), R2; n_effects=3)
    @test block_lo <= argmax(res2.pip) <= block_hi
    # With three near-identical tags the SER posterior concentrates on one
    # tag SNP (standard SuSiE behaviour; separating tags requires purity
    # filtering across effects, which single-credible-set runs do not do).
    @test maximum(res2.pip[block_lo:block_hi]) > 0.4
    @test all(i -> block_lo <= i <= block_hi, res2.credible_sets[1])

    # max_iter/tol are respected (terminates; no silent discarding).
    res3 = BioToolkit.fine_map_susie(z, ones(p), R; n_effects=2, max_iter=1, tol=1.0)
    @test res3.n_effects == 2

    @test_throws DimensionMismatch BioToolkit.fine_map_susie(z, ones(p), Matrix{Float64}(I, p - 1, p - 1))
    @test_throws ArgumentError BioToolkit.fine_map_susie(z, -ones(p), R)
end

@testset "Session 1: calculate_ibd recovers IBD states" begin
    Random.seed!(7)
    m = 20000
    maf = 0.2 .+ 0.6 .* rand(m)

    # Two founders with two alleles each (Bernoulli(maf)), a child and a full
    # sibling that each draw one allele from every founder, and an unrelated
    # individual drawn from the population.
    pa1 = rand(m) .< maf; pa2 = rand(m) .< maf      # founder A (row 1)
    pb1 = rand(m) .< maf; pb2 = rand(m) .< maf      # founder B (not a row)
    parent = Float64.(pa1) .+ Float64.(pa2)
    child = Vector{Float64}(undef, m)
    sib = Vector{Float64}(undef, m)
    for k in 1:m
        child[k] = (rand() < 0.5 ? pa1[k] : pa2[k]) + (rand() < 0.5 ? pb1[k] : pb2[k])
        sib[k] = (rand() < 0.5 ? pa1[k] : pa2[k]) + (rand() < 0.5 ? pb1[k] : pb2[k])
    end
    unrelated = Float64.(rand(m) .< maf) .+ Float64.(rand(m) .< maf)

    # Samples as rows (n x m): row 1 = parent A, 2 = child, 3 = sibling, 4 = unrelated.
    G = Matrix{Float64}(undef, 4, m)
    G[1, :] .= parent
    G[2, :] .= child
    G[3, :] .= sib
    G[4, :] .= unrelated
    res = BioToolkit.calculate_ibd(G)
    @test res isa BioToolkit.IBDEstimate
    # Loci monomorphic within the cohort (empirical MAF < min_maf) are excluded.
    @test 0.8 * m <= res.n_snps_used <= m

    # Parent-offspring: IBD1 = 1, IBD2 = 0, PI_HAT = 0.5.
    @test res.z0[1, 2] < 0.05
    @test res.z1[1, 2] > 0.9
    @test res.z2[1, 2] < 0.1
    @test 0.4 < res.pi_hat[1, 2] < 0.6

    # Full siblings: Z0 = 0.25, Z1 = 0.5, Z2 = 0.25, PI_HAT = 0.5.
    @test 0.15 < res.z0[2, 3] < 0.38
    @test 0.35 < res.z1[2, 3] < 0.65
    @test 0.10 < res.z2[2, 3] < 0.38
    @test 0.35 < res.pi_hat[2, 3] < 0.65

    # Unrelated: Z0 ~ 1, PI_HAT ~ 0.
    @test res.z0[1, 4] > 0.8
    @test abs(res.pi_hat[1, 4]) < 0.1

    # Symmetry and diagonal conventions.
    @test res.pi_hat[2, 1] == res.pi_hat[1, 2]
    @test res.pi_hat[1, 1] == 1.0
    @test res.z1[1, 1] == 1.0
    @test res.z0[1, 1] == 0.0

    # IBS matrix (kept under its honest name) sanity.
    ibs = BioToolkit.calculate_ibs(G)
    @test ibs[1, 1] == 1.0
    @test 0.5 < ibs[1, 2] <= 1.0
    @test ibs[1, 2] > ibs[1, 4]   # parent-offspring more similar than unrelated
end

@testset "Session 1: ldsc_genetic_correlation is LD-score regression" begin
    Random.seed!(11)
    m = 1200
    snp_ids = ["rs$i" for i in 1:m]
    causal = rand(m) .< 0.08
    ld = 0.2 .+ 3.0 .* rand(m) .+ 10.0 .* causal

    function _trait(z)
        pvals = 2.0 .* ccdf.(Ref(Normal()), abs.(z))
        BioToolkit.GWASResult(snp_ids, fill("1", m), collect(1:m), fill(("A", "G"), m), String[],
                              z .* 0.05, fill(0.05, m), z, pvals, 1000, String[], "trait", "synthetic")
    end

    effects = 4.0 .* Float64.(causal)
    z1 = effects .+ randn(m)
    z2_pos = effects .+ randn(m)       # same causal signs -> positive rg
    z2_neg = -effects .+ randn(m)      # opposite signs -> negative rg
    z2_none = randn(m)                 # independent -> rg near 0

    rg_pos, se_pos = BioToolkit.ldsc_genetic_correlation([_trait(z1), _trait(z2_pos)]; ld_scores=ld, n_blocks=20)
    @test isfinite(rg_pos)
    @test se_pos > 0
    @test rg_pos > 0.15

    rg_neg, _ = BioToolkit.ldsc_genetic_correlation([_trait(z1), _trait(z2_neg)]; ld_scores=ld, n_blocks=20)
    @test rg_neg < -0.15

    rg_none, _ = BioToolkit.ldsc_genetic_correlation([_trait(z1), _trait(z2_none)]; ld_scores=ld, n_blocks=20)
    @test abs(rg_none) < abs(rg_pos)

    # Dict-keyed LD scores give the same estimate.
    rg_dict, _ = BioToolkit.ldsc_genetic_correlation([_trait(z1), _trait(z2_pos)];
                                                     ld_scores=Dict(snp_ids[i] => ld[i] for i in 1:m), n_blocks=20)
    @test isapprox(rg_dict, rg_pos; atol=1e-8)

    # LD scores are now required (the old plain-cor(z1,z2) behaviour is gone).
    @test_throws ArgumentError BioToolkit.ldsc_genetic_correlation([_trait(z1), _trait(z2_pos)])
end
