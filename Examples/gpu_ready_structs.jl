using BioToolkit

records = [
    BioToolkit.parse_vcf_record("chr1\t10\trs1\tA\tG\t99.5"),
    BioToolkit.parse_vcf_record("chr1\t20\trs2\tC\tT\t80.0"),
    BioToolkit.parse_vcf_record("chr2\t30\trs3\tG\tA\t70.0"),
]

events = BioToolkit.compact_variant_event.(records)

println("GPU-ready example")
println("  element type isbits: ", isbitstype(eltype(events)))
println("  number of events: ", length(events))
println("  first event: ", events[1])
println("  text records stay on CPU: ", !isbitstype(BioToolkit.VariantTextRecord))
println("  BED records stay on CPU: ", !isbitstype(BioToolkit.BedRecord))

println()
println("Optional CUDA workflow if you have CUDA.jl installed:")
println("  1. `using CUDA`")
println("  2. `gpu_events = CuArray(events)`")
println("  3. Run your own kernel or broadcast over `gpu_events`")
println("  The package itself does not depend on CUDA.jl, so this step is optional.")
println("  CPU threaded paths remain available in the package APIs.")
