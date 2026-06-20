using BioToolkit
using Printf

include(joinpath(@__DIR__, "fixtures", "longread_fixture.jl"))
using .LongReadBenchmarkFixtures

function _sv_recovery_table(rows)
    io = IOBuffer()
    println(io, "| SV Type | TP | FP | FN | Precision | Recall | F1 |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|")

    for row in rows
        println(io,
            "| ", row.sv_type,
            " | ", row.tp,
            " | ", row.fp,
            " | ", row.fn,
            " | ", @sprintf("%.3f", row.precision),
            " | ", @sprintf("%.3f", row.recall),
            " | ", @sprintf("%.3f", row.f1),
            " |"
        )
    end

    return String(take!(io))
end

function run_longread_recovery_benchmark(; output_path::AbstractString=joinpath(@__DIR__, "results", "longread_recovery_metrics.md"))
    fixture = longread_fixture()
    predicted_calls = call_structural_variants(fixture.records; fixture.call_kwargs...)
    rows = evaluate_sv_recovery(predicted_calls, fixture.truth; position_tolerance=70)

    table = _sv_recovery_table(rows)
    mkpath(dirname(output_path))

    open(output_path, "w") do io
        println(io, "# Long-read Structural Variant Recovery")
        println(io)
        println(io, "Synthetic benchmark with known DEL/INV/TRA truth events. Metrics are computed by event-type aware matching with +/-70 bp positional tolerance.")
        println(io)
        println(io, "- Truth events: $(length(fixture.truth))")
        println(io, "- Predicted calls: $(length(predicted_calls))")
        println(io)
        print(io, table)
        println(io)
    end

    println(table)
    return (
        output_path=String(output_path),
        rows=rows,
        predicted_calls=predicted_calls,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    report = run_longread_recovery_benchmark()
    overall_index = findfirst(row -> row.sv_type == "ALL", report.rows)
    @assert !isnothing(overall_index) "Long-read benchmark missing ALL summary row"
    overall = report.rows[overall_index]

    @assert overall.precision >= 0.70 "Long-read benchmark failed precision recovery threshold"
    @assert overall.recall >= 0.95 "Long-read benchmark failed recall recovery threshold"
    println("Long-read recovery benchmark complete: ", report.output_path)
end
