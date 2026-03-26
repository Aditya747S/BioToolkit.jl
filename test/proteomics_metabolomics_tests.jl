using DataFrames
using Turing

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
    tracking = BioToolkit.metabolomics_source_tracking(observed, source_profiles; draws=10)
    @test tracking isa BioToolkit.MetabolomicsSourceTrackingResult
    @test length(tracking.mean_proportions) == 2

    annotated = BioToolkit.annotate_metabolite_features([1.0 2.0; 2.0 4.0])
    @test nrow(annotated) == 2
    @test "mean_intensity" in names(annotated)
end