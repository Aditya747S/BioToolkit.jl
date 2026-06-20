using Test
using LinearAlgebra
using Statistics
using SimpleWeightedGraphs
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