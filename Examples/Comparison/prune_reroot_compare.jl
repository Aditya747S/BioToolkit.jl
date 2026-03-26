pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using Random
using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_tree(sequence_count::Integer=80, width::Integer=140)
    Random.seed!(31415)
    template = repeat("ACGT", cld(width, 4))[1:width]
    bases = ('A', 'C', 'G', 'T')
    sequences = String[]

    for sequence_index in 1:sequence_count
        chars = collect(template)
        for position in 1:width
            if mod(sequence_index + position, 11) == 0
                chars[position] = bases[mod(sequence_index + position - 2, length(bases)) + 1]
            elseif mod(sequence_index * position, 17) == 0
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
    keep_names = ["seq$(index)" for index in 1:2:40]

    prune_ms, pruned = repeat_elapsed_ms(() -> BioToolkit.prune(tree, keep_names), 100)
    reroot_ms, rerooted = repeat_elapsed_ms(() -> BioToolkit.root_with_outgroup(tree, "seq1"), 100)

    @assert length(BioToolkit.get_terminals(pruned)) == length(keep_names)
    @assert occursin("seq1", BioToolkit.write_newick(rerooted))

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "prune_reroot_compare.py"))`, String)
    python_prune_ms = parse(Float64, match(r"python_prune_ms=([0-9.]+)", python_output).captures[1])
    python_reroot_ms = parse(Float64, match(r"python_reroot_ms=([0-9.]+)", python_output).captures[1])
    python_prune_ok = match(r"python_prune_ok=(true|false)", python_output).captures[1] == "true"
    python_reroot_ok = match(r"python_reroot_ok=(true|false)", python_output).captures[1] == "true"

    @assert python_prune_ok
    @assert python_reroot_ok

    println("Pruning and rerooting benchmark")
    println("  sequence_count=80")
    println("  julia_prune_ms=$(round(prune_ms / 100; digits=4))")
    println("  julia_reroot_ms=$(round(reroot_ms / 100; digits=4))")
    println("  python_prune_ms=$(round(python_prune_ms; digits=4))")
    println("  python_reroot_ms=$(round(python_reroot_ms; digits=4))")
    print(python_output)
end

main()