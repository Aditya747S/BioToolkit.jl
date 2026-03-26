pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using Random
using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_distance_matrix(n::Integer=200)
    Random.seed!(2024)
    D = zeros(Float64, n, n)
    for i in 1:n
        for j in (i+1):n
            val = rand()
            D[i, j] = val
            D[j, i] = val
        end
    end
    names = ["Taxon_$i" for i in 1:n]
    return D, names
end

function main()
    D, names = build_distance_matrix(200)
    
    # Warmup
    BioToolkit.neighbor_joining_tree(D, names)
    
    # Benchmark
    nj_ms, _ = repeat_elapsed_ms(() -> BioToolkit.neighbor_joining_tree(D, names), 10)
    
    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "nj_compare.py"))`, String)
    python_nj_ms = parse(Float64, match(r"python_nj_ms=([0-9.]+)", python_output).captures[1])
    
    println("Neighbor Joining Tree Construction benchmark")
    println("  taxa=200")
    println("  julia_nj_ms=$(round(nj_ms / 10; digits=4))")
    println("  python_nj_ms=$(round(python_nj_ms; digits=4))")
    print(python_output)
end

main()
