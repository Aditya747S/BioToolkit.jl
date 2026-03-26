pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

protein_seq = repeat("ACDEFGHIKLMNPQRSTVWY", 1000) # 20k residues
dna_seq = repeat("ACGTGCGCAGCT", 5000) # 60k bp

function calc_protein_basic()
    BioToolkit.protein_mass(protein_seq)
    BioToolkit.extinction_coefficient(protein_seq)
    BioToolkit.instability_index(protein_seq)
    BioToolkit.isoelectric_point(protein_seq)
    BioToolkit.gravy(protein_seq)
end

function calc_protein_bundle()
    BioToolkit.aliphatic_index(protein_seq)
    BioToolkit.protparam(protein_seq)
end

function calc_skew()
    BioToolkit.gc_skew(dna_seq)
end

function main()
    # verify logic works
    _ = calc_protein_basic()
    _ = calc_protein_bundle()
    _ = calc_skew()

    reps = 100

    protein_ms, _ = repeat_elapsed_ms(calc_protein_basic, reps)
    bundle_ms, _ = repeat_elapsed_ms(calc_protein_bundle, reps)
    skew_ms, _ = repeat_elapsed_ms(calc_skew, reps)

    protein_ms /= reps
    bundle_ms /= reps
    skew_ms /= reps


    println("Julia Protein Stats benchmark")
    println("  reps=$reps")
    println("  julia_protein_ms=$(round(protein_ms; digits=4))")
    println("  julia_protein_bundle_ms=$(round(bundle_ms; digits=4))")
    println("  julia_skew_ms=$(round(skew_ms; digits=4))")
    
    # Run Python comparisons
    run(`conda run -n general python $(joinpath(@__DIR__, "protein_stats_compare.py"))`)
end

main()
