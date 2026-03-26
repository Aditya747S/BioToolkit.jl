pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_sequences()
    left = repeat("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 30)
    right = replace(left, 'A' => 'T'; count=180)
    return left, right
end

function main()
    left, right = build_sequences()
    left_record = SeqRecordLite(left; identifier="left", description="left sequence")
    right_record = SeqRecordLite(right; identifier="right", description="right sequence")

    BioToolkit.pairwise_align(left_record, right_record; match=2, mismatch=-1, gap=-2)

    alignment_ms, alignment = elapsed_ms() do
        BioToolkit.pairwise_align(left_record, right_record; match=2, mismatch=-1, gap=-2)
    end

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "pairwise_align_compare.py"))`, String)
    python_score = parse(Float64, match(r"python_pairwise_score=(-?[0-9.]+)", python_output).captures[1])
    python_matches = parse(Int, match(r"python_pairwise_matches=(\d+)", python_output).captures[1])
    python_identity = parse(Float64, match(r"python_pairwise_identity=([0-9.]+)", python_output).captures[1])
    python_length = parse(Int, match(r"python_pairwise_length=(\d+)", python_output).captures[1])

    @assert isapprox(Float64(alignment.score), python_score; atol=1e-8)
    @assert alignment.matches == python_matches
    @assert isapprox(round(alignment.identity, digits=6), python_identity; atol=1e-12)
    @assert length(alignment.left) == python_length

    println("Pairwise alignment benchmark")
    println("  left_length=", length(left))
    println("  right_length=", length(right))
    println("  julia_pairwise_ms=", round(alignment_ms, digits=4))
    println("  julia_pairwise_score=", alignment.score)
    println("  julia_pairwise_identity=", round(alignment.identity, digits=6))
    print(python_output)
end

main()
