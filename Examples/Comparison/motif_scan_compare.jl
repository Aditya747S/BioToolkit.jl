pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_motif_inputs()
    instances = [
        SeqRecordLite("TTTACGTTTT"; identifier="s1"),
        SeqRecordLite("GGGACGTGGG"; identifier="s2"),
        SeqRecordLite("CCCACGTCCC"; identifier="s3"),
        SeqRecordLite("TTTACGTTTT"; identifier="s4"),
    ]
    sequence = repeat("TTTACGTTTTGGGACGTGGGCCCACGTCCC", 1000)
    return instances, sequence
end

function main()
    instances, sequence = build_motif_inputs()
    counts = BioToolkit.motif_counts(instances)
    pwm = BioToolkit.motif_pwm(counts; pseudocount=0.5)

    JuliaScanner = () -> BioToolkit.motif_scan_both_strands(sequence, pwm; threshold=0.0)
    JuliaScanner()
    julia_ms, hits = elapsed_ms(JuliaScanner)

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "motif_scan_compare.py"))`, String)

    python_ms = parse(Float64, match(r"python_motif_scan_ms=([0-9.]+)", python_output).captures[1])
    python_hit_count = parse(Int, match(r"python_motif_hit_count=(\d+)", python_output).captures[1])
    python_first_score = parse(Float64, match(r"python_motif_first_score=([-0-9.]+)", python_output).captures[1])
    python_ok = match(r"python_motif_scan_ok=(true|false)", python_output).captures[1] == "true"

    @assert python_ok
    @assert length(hits) == python_hit_count
    @assert !isempty(hits)
    @assert isapprox(hits[1].score, python_first_score; atol=1e-6)

    println("Motif scan benchmark")
    println("  sequence_length=", length(sequence))
    println("  julia_scan_ms=", round(julia_ms, digits=4))
    println("  python_scan_ms=", round(python_ms, digits=4))
    println("  julia_hit_count=", length(hits))
    println("  julia_information_content=", round(BioToolkit.motif_information_content(pwm), digits=6))
    print(python_output)
end

main()
