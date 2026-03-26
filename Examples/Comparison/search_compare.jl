pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))
using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

println("Generating sequences...")
# Target: 100,000 bp
target_base = repeat("ACGT", 25000)
target = target_base[1:50000] * "GGCCGATCGATCGATCGATCGA" * target_base[50023:end]
targets = [target]

# Query: 1,000 bp, containing the homologous sequence
query = repeat("T", 500) * "GGCCGATCGATCGATCGATCGA" * repeat("A", 478)

# Setup scoring
matrix = BioToolkit.substitution_matrix("ACGT"; match=2, mismatch=-1)
scoring = BioToolkit.MatrixPairwiseScoring(matrix)

println("Building K-mer Index...")
index_time = @elapsed begin
    index = BioToolkit.build_index(targets, ["chr1"]; k=6)
end
println("Index built in $(round(index_time * 1000, digits=2)) ms")

function do_search()
    BioToolkit.local_search(query, index; scoring=scoring, x_drop=8, min_score=20)
end

function do_smith_waterman()
    BioToolkit.local_align(query, target; match=2, mismatch=-1, gap=-2)
end

function main()
    println("Warming up JIT...")
    # Warmup
    _ = do_search()
    _ = do_smith_waterman()
    
    println("Running Heuristic Search (100 iterations)...")
    search_ms, hsps = repeat_elapsed_ms(do_search, 100)
    
    println("Running Smith-Waterman (1 iteration)...")
    sw_ms, sw_res = repeat_elapsed_ms(do_smith_waterman, 1) # Only 1 rep because it's O(N*M)
    
    println("\n--- BLAST-like Heuristic vs Smith-Waterman ---")
    println("Target size: 100,000 bp")
    println("Query size:  1,000 bp")
    println("Result: Heuristic Search: $(round(search_ms / 100, digits=3)) ms per iter")
    println("Result: Smith-Waterman:   $(round(sw_ms / 1, digits=3)) ms per iter")
    
    speedup = (sw_ms / 1) / (search_ms / 100)
    println("Heuristic Speedup:        $(round(speedup, digits=1))x")
    
    println("\nValidating outputs:")
    println("  Smith-Waterman score: ", sw_res.score)
    println("  Heuristic best score: ", isempty(hsps) ? 0 : hsps[1].score)
end

main()
