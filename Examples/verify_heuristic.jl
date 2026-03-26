pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..")))
using BioToolkit
include(joinpath(@__DIR__, "..", "scripts", "benchmark_helpers.jl"))

println("Generating sequences with GAPS for robust verification...")
# Target: 5,000 bp
target_base = repeat("ACGT", 1250)

# Plant a homologous sequence with an insertion and deletion
# True sequence: "GGAATTCC" * "TT" (insertion) * "GG" * "C" (deletion) * "AATT"
# Query will just look for "GGAATTCCGGAATT"
# Let's just use something simple
homology_target = "GGCCGATCGATC" * "A" * "GATCGATCGA" # inserted 'A'
target = target_base[1:2500] * homology_target * target_base[2524:end]
targets = [target]

# Query: 500 bp, containing the homologous sequence WITHOUT the insertion
homology_query  = "GGCCGATCGATC" * "GATCGATCGA"
query = repeat("T", 200) * homology_query * repeat("A", 278)

matrix = BioToolkit.substitution_matrix("ACGT"; match=2, mismatch=-1)
scoring = BioToolkit.MatrixPairwiseScoring(matrix)

index = BioToolkit.build_index(targets; k=6)

println("\n--- Heuristic Search (Ungapped only) ---")
hsps = BioToolkit.local_search(query, index; scoring=scoring, x_drop=5, min_score=10)
if isempty(hsps)
    println("No HSP found!")
else
    best_hsp = hsps[1]
    println("Found HSP natively! Score: ", best_hsp.score)
    println("Query bounds: $(best_hsp.query_start) to $(best_hsp.query_end)")
    if best_hsp.score < 44 # Full match is 22 * 2 = 44. Ungapped will break at the 'A' insertion.
        println("As expected, ungapped extension broke at the insertion!")
    end
end

println("\n--- Exact Smith-Waterman ---")
sw_res = BioToolkit.local_align(query, target; match=2, mismatch=-1, gap=-2)
println("SW Score: $(sw_res.score)")
