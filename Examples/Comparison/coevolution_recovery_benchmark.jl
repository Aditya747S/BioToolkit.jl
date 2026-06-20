using BioToolkit
using Printf

include(joinpath(@__DIR__, "fixtures", "coevolution_fixture.jl"))
using .CoevolutionBenchmarkFixtures

function _contact_metric_table(rows)
    io = IOBuffer()
    println(io, "| Evaluation | Top-N | Precision | Recall | F1 | Mean Signal Gap | Mean Rank Recovery | Mean True-Pair Distance |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")

    for row in rows
        println(io,
            "| ", row.evaluation,
            " | ", row.top_n,
            " | ", @sprintf("%.3f", row.precision),
            " | ", @sprintf("%.3f", row.recall),
            " | ", @sprintf("%.3f", row.f1),
            " | ", @sprintf("%.3f", row.signal_gap),
            " | ", @sprintf("%.3f", row.rank_recovery),
            " | ", @sprintf("%.3f", row.mean_true_pair_distance),
            " |"
        )
    end

    return String(take!(io))
end

function run_coevolution_recovery_benchmark(; output_path::AbstractString=joinpath(@__DIR__, "results", "coevolution_recovery_metrics.md"))
    fixture = coevolution_fixture()
    alignment_size = size(fixture.alignment)
    predicted_scores = predict_contact_map(fixture.alignment; pseudocount=1e-3)

    strict_top_n = length(fixture.truth_pairs)
    relaxed_top_n = 2 * strict_top_n

    strict_metrics = evaluate_contact_recovery(predicted_scores, fixture.truth_pairs; top_n=strict_top_n)
    relaxed_metrics = evaluate_contact_recovery(predicted_scores, fixture.truth_pairs; top_n=relaxed_top_n)

    rows = [
        merge((evaluation="Top-Truth", top_n=strict_top_n), strict_metrics),
        merge((evaluation="Top-2xTruth", top_n=relaxed_top_n), relaxed_metrics),
    ]

    table = _contact_metric_table(rows)
    mkpath(dirname(output_path))

    open(output_path, "w") do io
        println(io, "# Coevolution Contact Recovery")
        println(io)
        println(io, "Synthetic MSA benchmark with implanted coevolving residue pairs and known 3D-contact distances. Metrics capture precision-recall and score separation quality.")
        println(io)
        println(io, "- Alignment depth: $(alignment_size[1])")
        println(io, "- Alignment length: $(alignment_size[2])")
        println(io, "- Implanted contact pairs: $(length(fixture.truth_pairs))")
        println(io)
        print(io, table)
        println(io)
    end

    println(table)
    return (output_path=String(output_path), rows=rows)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    report = run_coevolution_recovery_benchmark()

    strict_row = report.rows[1]
    @assert strict_row.precision >= 0.45 "Coevolution benchmark failed strict precision threshold"
    @assert strict_row.recall >= 0.45 "Coevolution benchmark failed strict recall threshold"

    println("Coevolution recovery benchmark complete: ", report.output_path)
end
