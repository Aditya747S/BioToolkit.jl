pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using Random
using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_records(record_count::Integer=120, width::Integer=80)
    Random.seed!(2024)
    template = repeat("ACGT", cld(width, 4))[1:width]
    bases = ('A', 'C', 'G', 'T')
    records = SeqRecordLite[]

    for record_index in 1:record_count
        chars = collect(template)
        for position in 1:width
            if mod(record_index + position, 17) == 0
                chars[position] = '-'
            elseif mod(record_index * position, 19) == 0
                chars[position] = bases[mod(record_index + position - 2, length(bases)) + 1]
            end
        end
        push!(records, SeqRecordLite(String(chars); identifier="seq$(record_index)", description="bootstrap_input"))
    end

    return records
end

function main()
    alignment = MultipleSequenceAlignment(build_records())
    bootstrap_ms, consensus = repeat_elapsed_ms(() -> BioToolkit.bootstrap_consensus_tree(alignment; replicates=25, threshold=0.5, method=:hamming, constructor=:nj), 5)

    newick = BioToolkit.write_newick(consensus)
    @assert occursin("0.", newick) || occursin("1.0", newick)
    @assert length(BioToolkit.get_terminals(consensus)) == length(alignment.records)

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "consensus_support_newick_compare.py"))`, String)
    python_consensus_ms = parse(Float64, match(r"python_consensus_support_newick_ms=([0-9.]+)", python_output).captures[1])
    python_consensus_ok = match(r"python_consensus_support_newick_ok=(true|false)", python_output).captures[1] == "true"

    @assert python_consensus_ok

    println("Consensus support Newick benchmark")
    println("  record_count=120")
    println("  julia_consensus_support_newick_ms=$(round(bootstrap_ms / 5; digits=4))")
    println("  python_consensus_support_newick_ms=$(round(python_consensus_ms; digits=4))")
    println("  julia_newick_length=$(length(newick))")
    print(python_output)
end

main()
