using DataFrames
if Base.find_package("Turing") !== nothing
    using Turing
end

@testset "Proteomics and Metabolomics" begin
    mzml_text = """<mzML>
<spectrum id=\"scan1\">
<cvParam name=\"scan start time\" value=\"1.5\"/>
<mz>100.0 150.0 200.0 260.0</mz>
<intensity>10.0 50.0 5.0 60.0</intensity>
</spectrum>
<spectrum id=\"scan2\">
<cvParam name=\"scan start time\" value=\"2.0\"/>
<mz>100.0 155.0 210.0 260.0</mz>
<intensity>8.0 42.0 6.0 58.0</intensity>
</spectrum>
</mzML>"""

    mktempdir() do dir
        path = joinpath(dir, "synthetic.mzml")
        open(path, "w") do io
            write(io, mzml_text)
        end

        experiment = BioToolkit.read_mzml(path)
        @test experiment isa BioToolkit.MassSpecExperiment
        @test length(experiment.spectra) == 2
        @test experiment.spectra[1].rt == 1.5

        peak_result = BioToolkit.detect_peaks(experiment.spectra[1]; threshold=0.1)
        @test peak_result isa BioToolkit.PeakDetectionResult
        @test !isempty(peak_result.peaks)

        alignment = BioToolkit.align_samples(experiment.spectra[1].intensity, experiment.spectra[2].intensity)
        @test alignment isa BioToolkit.AlignmentResult
        @test !isempty(alignment.path)

        masses = [0.0, 57.02146, 128.05857]
        denovo = BioToolkit.de_novo_sequence(masses; tolerance=0.2)
        @test denovo isa BioToolkit.DeNovoResult
        @test denovo.sequence == "GA"
    end

    @testset "Typed protein biotypes" begin
        protein = BioToolkit.AASeq("ACDEFGHIKLMNPQRSTVWY")
        short = BioToolkit.AASeq("YYY")

        @test BioToolkit.protein_mass(protein; type="monoisotopic") > 2000.0
        @test BioToolkit.protein_mass(protein; type="average") > BioToolkit.protein_mass(protein; type="monoisotopic")
        @test BioToolkit.extinction_coefficient(short) == 4470
        @test BioToolkit.instability_index(protein) > 0.0
        @test BioToolkit.gravy(protein) isa Float64
        @test BioToolkit.aliphatic_index(protein) > 0.0

        params = BioToolkit.protparam(protein)
        @test params.length == 20
        @test params.extinction_coefficient == 7115
        @test params.isoelectric_point > 0.0
    end

    matrix = Matrix{Union{Missing,Float64}}([
        1.0 2.0 0.0 4.0;
        2.0 missing 1.0 3.0;
        4.0 1.0 2.0 missing
    ])
    imputed = BioToolkit.qrilc_impute(matrix)
    @test size(imputed) == size(matrix)
    @test all(isfinite, imputed)

    groups = ["case", "case", "ctrl", "ctrl"]
    abund = BioToolkit.differential_abundance(matrix, groups)
    @test abund isa BioToolkit.DifferentialAbundanceResult
    @test size(abund.coefficients, 1) == 3

    pls = BioToolkit.sparse_pls_da([1.0 2.0 3.0; 2.0 2.5 3.5; 3.0 3.5 4.0; 4.0 4.5 5.0], ["A", "A", "B", "B"])
    @test pls isa BioToolkit.SparsePLSDAResult
    @test !isempty(pls.selected_features)

    source_profiles = [0.7 0.3;
                       0.2 0.8;
                       0.1 0.1]
    observed = [7, 2, 1]
    if Base.get_extension(BioToolkit, :BioToolkitTuringExt) !== nothing || isdefined(Main, :Turing)
        tracking = BioToolkit.metabolomics_source_tracking(observed, source_profiles; draws=10)
        @test tracking isa BioToolkit.MetabolomicsSourceTrackingResult
        @test length(tracking.mean_proportions) == 2
    end

    annotated = BioToolkit.annotate_metabolite_features([1.0 2.0; 2.0 4.0])
    @test nrow(annotated) == 2
    @test "mean_intensity" in names(annotated)
end
# ==========================================================================
# Session 3 corrections (2026-09-06): hypergeometric tail, real incomplete
# beta, charge-correct adduct masses, isotope enrichment denominator,
# batch-correction provenance binding.
# Evidence: plan/audit_notes.md wave 6; plan/correction_plan.md items 28-32.
# ==========================================================================
using Test
using BioToolkit
using Random
using Statistics
using SpecialFunctions

@testset "Session 3: hypergeometric tail (not point probability)" begin
    # Deep-enrichment case: a = overlap 8, pathway 10, query 10, bg 100.
    # Point probability P(X=8) < tail P(X>=8); the old code returned the point.
    a, b, c, d = 8, 2, 2, 88
    p = BioToolkit.Metabolomics._hypergeometric_pvalue(a, b, c, d)
    @test 0.0 < p <= 1.0
    # tail must be >= point probability
    point = exp((loggamma(10) - loggamma(8) - loggamma(2)) +
                (loggamma(90) - loggamma(2) - loggamma(88)) -
                (loggamma(100) - loggamma(10) - loggamma(90)))
    @test p >= point
    @test p < 0.05                       # 8 of 10 hits in a 10/100 pathway: clearly enriched
    # Non-enriched: overlap equals random expectation
    p_null = BioToolkit.Metabolomics._hypergeometric_pvalue(1, 9, 9, 81)
    @test p_null > 0.3
end

@testset "Session 3: incomplete beta equals known values" begin
    # I_0.5(0.5, 0.5) = 0.5 by symmetry
    @test isapprox(BioToolkit.Metabolomics._ibeta(0.5, 0.5, 0.5), 0.5; atol=1e-10)
    # I_x(1, 1) = x
    @test isapprox(BioToolkit.Metabolomics._ibeta(1.0, 1.0, 0.3), 0.3; atol=1e-12)
    @test isapprox(BioToolkit.Metabolomics._ibeta(1.0, 1.0, 0.9), 0.9; atol=1e-12)
    # I_x(a,b) + I_{1-x}(b,a) = 1
    @test isapprox(BioToolkit.Metabolomics._ibeta(2.0, 3.0, 0.4) +
                   BioToolkit.Metabolomics._ibeta(3.0, 2.0, 0.6), 1.0; atol=1e-12)
    # t-test p-value: t = 2, df = 5 → p ≈ 0.1019 (textbook value)
    df = 5.0; t = 2.0
    x = df / (df + t^2)
    p = BioToolkit.Metabolomics._ibeta(df/2, 0.5, x)
    @test isapprox(p, 0.10193947882985832; atol=1e-10)   # matches Distributions.ccdf(TDist(5), 2)*2
end

@testset "Session 3: charge-correct adduct m/z" begin
    db = DataFrame(name=["glucose"], formula=["C6H12O6"], monoisotopic_mass=[180.0634], hmdb_id=[""], kegg_id=[""], class=[""])
    mzs = [180.0634 + 1.007276,            # [M+H]+
           (180.0634 + 2*1.007276)/2,      # [M+2H]2+
           2*180.0634 + 1.007276]          # [2M+H]+
    result = BioToolkit.annotate_adducts(mzs, db; ppm_tol=5.0, mode=:positive)
    adducts_found = Set(result.adduct)
    @test "[M+H]+" in adducts_found
    @test "[M+2H]2+" in adducts_found      # previously impossible (delta wrong by M/2)
    @test "[2M+H]+" in adducts_found       # previously off by M
    # each observed mz annotates to the right adduct at ~0 ppm
    for row in eachrow(result)
        @test abs(row.delta_ppm) < 2.0
    end
end

@testset "Session 3: isotope enrichment denominator" begin
    Random.seed!(4)
    # 6-carbon metabolite, 4 channels (M+0..M+3), one sample
    labeled = reshape(Float64[40, 30, 20, 10], 1, 4)     # fractions .4/.3/.2/.1
    unlabeled = reshape(Float64[98, 1, 0.7, 0.3], 1, 4)
    traces = BioToolkit.isotope_tracer_analysis(labeled, unlabeled, ["met6"]; n_carbons=[6])
    t = traces[1]
    @test t.n_carbons == 6
    # mean enrichment = (0.3*1 + 0.2*2 + 0.1*3)/6 = 10/6 ≈ 0.1667
    @test isapprox(t.mean_enrichment, (0.3 + 0.4 + 0.3) / 6; atol=1e-12)
    @test t.excess_enrichment ≈ t.mean_enrichment - 0.011 atol=1e-12
end

@testset "Session 3: batch correction runs with provenance" begin
    Random.seed!(5)
    mat = rand(20, 12) .* 100
    batches = vcat(fill("A", 6), fill("B", 6))
    corrected = BioToolkit.batch_effect_correction_metabolomics(mat, batches; method=:combat_mean)
    @test size(corrected) == size(mat)
    # batch means now aligned to the grand mean
    gm = vec(mean(mat, dims=2))
    a_mean = vec(mean(corrected[:, 1:6], dims=2))
    b_mean = vec(mean(corrected[:, 7:12], dims=2))
    @test maximum(abs.(a_mean .- gm)) < 1e-9
    @test maximum(abs.(b_mean .- gm)) < 1e-9
end
