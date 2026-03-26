pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using Random
using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_tree(sequence_count::Integer=40, width::Integer=100)
    Random.seed!(4242)
    template = repeat("ACGT", cld(width, 4))[1:width]
    bases = ('A', 'C', 'G', 'T')
    sequences = String[]

    for sequence_index in 1:sequence_count
        chars = collect(template)
        for position in 1:width
            if mod(sequence_index + position, 13) == 0
                chars[position] = bases[mod(sequence_index + position - 2, length(bases)) + 1]
            elseif mod(sequence_index * position, 23) == 0
                chars[position] = bases[mod(sequence_index + 2 * position - 3, length(bases)) + 1]
            end
        end
        push!(sequences, String(chars))
    end

    labels = ["seq$(index)" for index in eachindex(sequences)]
    distances = BioToolkit.distance_matrix(sequences; method=:hamming)
    return BioToolkit.neighbor_joining_tree(distances, labels)
end

function main()
    tree = build_tree()
    dot_ms, dot_output = repeat_elapsed_ms(() -> BioToolkit.tree_to_dot(tree), 200)
    mermaid_ms, mermaid_output = repeat_elapsed_ms(() -> BioToolkit.tree_to_mermaid(tree), 200)

    @assert occursin("digraph", dot_output)
    @assert occursin("graph TD", mermaid_output)

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "graph_export_compare.py"))`, String)
    python_ascii_ms = parse(Float64, match(r"python_tree_ascii_ms=([0-9.]+)", python_output).captures[1])
    python_newick_ms = parse(Float64, match(r"python_tree_newick_ms=([0-9.]+)", python_output).captures[1])
    python_export_ok = match(r"python_tree_export_ok=(true|false)", python_output).captures[1] == "true"

    @assert python_export_ok

    println("Phylogenetics export benchmark")
    println("  julia_tree_dot_ms=$(round(dot_ms / 200; digits=4))")
    println("  julia_tree_mermaid_ms=$(round(mermaid_ms / 200; digits=4))")
    println("  python_tree_ascii_ms=$(round(python_ascii_ms; digits=4))")
    println("  python_tree_newick_ms=$(round(python_newick_ms; digits=4))")
    print(python_output)
end

main()
