pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using Random
using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_sequences(sequence_count::Integer=60, width::Integer=160)
    Random.seed!(1234)
    template = repeat("ACGT", cld(width, 4))[1:width]
    bases = ('A', 'C', 'G', 'T')
    sequences = String[]

    for sequence_index in 1:sequence_count
        chars = collect(template)
        for position in 1:width
            if mod(sequence_index + position, 13) == 0
                chars[position] = bases[mod(sequence_index + position - 2, length(bases)) + 1]
            elseif mod(sequence_index * position, 29) == 0
                chars[position] = bases[mod(sequence_index + 2 * position - 3, length(bases)) + 1]
            end
        end
        push!(sequences, String(chars))
    end

    return sequences
end

function build_tree()
    sequences = build_sequences()
    labels = ["seq$(index)" for index in eachindex(sequences)]
    distances = BioToolkit.distance_matrix(sequences; method=:hamming)
    return BioToolkit.neighbor_joining_tree(distances, labels)
end

function main()
    tree = build_tree()
    midpoint_ms, midpoint_tree = repeat_elapsed_ms(() -> BioToolkit.midpoint_root(tree), 100)

    @assert length(BioToolkit.get_terminals(midpoint_tree)) == length(BioToolkit.get_terminals(tree))
    @assert occursin("(", BioToolkit.write_newick(midpoint_tree))

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "midpoint_root_compare.py"))`, String)
    python_midpoint_ms = parse(Float64, match(r"python_midpoint_root_ms=([0-9.]+)", python_output).captures[1])
    python_midpoint_ok = match(r"python_midpoint_root_ok=(true|false)", python_output).captures[1] == "true"

    @assert python_midpoint_ok

    println("Midpoint rooting benchmark")
    println("  sequence_count=60")
    println("  julia_midpoint_root_ms=$(round(midpoint_ms / 100; digits=4))")
    println("  python_midpoint_root_ms=$(round(python_midpoint_ms; digits=4))")
    print(python_output)
end

main()
