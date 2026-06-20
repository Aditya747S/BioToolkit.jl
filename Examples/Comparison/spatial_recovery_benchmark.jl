using BioToolkit
using Printf

include(joinpath(@__DIR__, "fixtures", "spatial_fixture.jl"))
using .SpatialBenchmarkFixtures

function _fraction_metric_table(results)
    io = IOBuffer()
    println(io, "| Method | RMSE | MAE | Mean Correlation | Dominant Accuracy | Mean JSD |")
    println(io, "|---|---:|---:|---:|---:|---:|")

    for row in results
        println(io,
            "| ", row.method,
            " | ", @sprintf("%.4f", row.rmse),
            " | ", @sprintf("%.4f", row.mae),
            " | ", @sprintf("%.4f", row.mean_correlation),
            " | ", @sprintf("%.4f", row.dominant_accuracy),
            " | ", @sprintf("%.4f", row.mean_jsd),
            " |"
        )
    end

    return String(take!(io))
end

function run_spatial_recovery_benchmark(; output_path::AbstractString=joinpath(@__DIR__, "results", "spatial_recovery_metrics.md"))
    fixture = spatial_deconvolution_fixture()

    rctd_result = rctd_deconvolution(
        fixture.spatial,
        fixture.reference,
        fixture.reference_labels,
        marker_top_n=25,
        maxiter=350,
    )

    c2l_result = cell2location_deconvolution(
        fixture.spatial,
        fixture.reference,
        fixture.reference_labels,
        marker_top_n=25,
        maxiter=350,
    )

    methods = NamedTuple[
        (method="RCTD-like", result=rctd_result),
        (method="Cell2Location-like", result=c2l_result),
    ]

    rows = NamedTuple[]
    for method_result in methods
        aligned_predicted = align_cell_type_order(method_result.result, fixture.cell_type_ids)
        metrics = evaluate_fraction_recovery(aligned_predicted, fixture.truth_fractions)
        push!(rows, merge((method=method_result.method,), metrics))
    end

    table = _fraction_metric_table(rows)
    mkpath(dirname(output_path))

    open(output_path, "w") do io
        println(io, "# Spatial Deconvolution Recovery")
        println(io)
        println(io, "Synthetic benchmark with known cell-type fractions per spot. Lower RMSE/MAE/JSD and higher correlation/dominant accuracy indicate better recovery.")
        println(io)
        println(io, "- Genes: $(size(fixture.reference.counts, 1))")
        println(io, "- Reference cells: $(size(fixture.reference.counts, 2))")
        println(io, "- Spatial spots: $(size(fixture.spatial.experiment.counts, 2))")
        println(io)
        print(io, table)
        println(io)
    end

    println(table)
    return (output_path=String(output_path), rows=rows)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    report = run_spatial_recovery_benchmark()

    for row in report.rows
        @assert row.rmse <= 0.30 "Spatial benchmark failed RMSE threshold for $(row.method)"
        @assert row.mean_correlation >= 0.80 "Spatial benchmark failed correlation threshold for $(row.method)"
    end

    println("Spatial recovery benchmark complete: ", report.output_path)
end
