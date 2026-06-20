using BioToolkit
using Printf

include(joinpath(@__DIR__, "fixtures", "pipeline_fixture.jl"))
using .PipelineBenchmarkFixtures

function _pipeline_metric_table(metrics)
    io = IOBuffer()
    println(io, "| Metric | Observed | Target | Pass |")
    println(io, "|---|---:|---:|:---:|")

    for row in metrics
        observed_str = isa(row.observed, Bool) ? string(row.observed) : @sprintf("%.4f", Float64(row.observed))
        target_str = isa(row.target, Bool) ? string(row.target) : @sprintf("%.4f", Float64(row.target))
        pass_str = row.pass ? "yes" : "no"
        println(io, "| $(row.metric) | $(observed_str) | $(target_str) | $(pass_str) |")
    end

    return String(take!(io))
end

function run_pipeline_recovery_benchmark(; output_path::AbstractString=joinpath(@__DIR__, "results", "pipeline_recovery_metrics.md"))
    fixture = pipeline_fixture()
    report = evaluate_pipeline_recovery(fixture)
    table = _pipeline_metric_table(report.rows)

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, "# Pipeline DAG Recovery")
        println(io)
        println(io, "Synthetic deterministic DAG benchmark that evaluates structural validity, execution equivalence, parallel consistency, and cache reuse.")
        println(io)
        println(io, "- Nodes: $(length(fixture.graph.nodes))")
        println(io, "- Initial cache primed: yes")
        println(io)
        print(io, table)
        println(io)
    end

    println(table)
    return (
        output_path=String(output_path),
        rows=report.rows,
        outputs_serial=report.outputs_serial,
        outputs_parallel=report.outputs_parallel,
        cache_hit_rate=report.cache_hit_rate,
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    report = run_pipeline_recovery_benchmark()
    for row in report.rows
        @assert row.pass "Pipeline benchmark failed metric: $(row.metric)"
    end

    println("Pipeline recovery benchmark complete: ", report.output_path)
end
