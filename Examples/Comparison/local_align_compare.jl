pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))
using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

seq1 = repeat("ACGTGCCG", 125) # 1000 bp
seq2 = repeat("ACGTGCG", 100)  # 700 bp

function do_local()
    BioToolkit.local_align(seq1, seq2; match=2, mismatch=-1, gap_open=-5, gap_extend=-1)
end

function main()
    do_local() # heat
    reps = 10
    total_ms, _ = repeat_elapsed_ms(do_local, reps)
    ms_per_rep = total_ms / reps

    println("Julia Local Align benchmark")
    println("  reps=$reps")
    println("  julia_local_align_ms=$(round(ms_per_rep; digits=4))")

    run(`conda run -n general python $(joinpath(@__DIR__, "local_align_compare.py"))`)
end

main()
