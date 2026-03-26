pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_sequences()
    left = join(("ACGT"[(mod(index * 7 + 3, 4) + 1)] for index in 1:1200))
    right = left[1:500] * left[571:end]
    return left, right
end

function main()
    left, right = build_sequences()
    left_record = SeqRecordLite(left; identifier="left", description="affine benchmark left sequence")
    right_record = SeqRecordLite(right; identifier="right", description="affine benchmark right sequence")

    BioToolkit.pairwise_align(left_record, right_record; match=2, mismatch=-1, gap_open=-5, gap_extend=-1)

    alignment_ms, alignment = elapsed_ms() do
        BioToolkit.pairwise_align(left_record, right_record; match=2, mismatch=-1, gap_open=-5, gap_extend=-1)
    end

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "pairwise_affine_compare.py"))`, String)
    python_score = parse(Float64, match(r"python_pairwise_affine_score=(-?[0-9.]+)", python_output).captures[1])
    python_matches = parse(Int, match(r"python_pairwise_affine_matches=(\d+)", python_output).captures[1])
    python_identity = parse(Float64, match(r"python_pairwise_affine_identity=([0-9.]+)", python_output).captures[1])
    python_length = parse(Int, match(r"python_pairwise_affine_length=(\d+)", python_output).captures[1])

    @assert isapprox(Float64(alignment.score), python_score; atol=1e-8)
    @assert alignment.matches == python_matches
    @assert isapprox(round(alignment.identity, digits=6), python_identity; atol=1e-12)
    @assert length(alignment.left) == python_length

    println("Affine-gap alignment benchmark")
    println("  left_length=", length(left))
    println("  right_length=", length(right))
    println("  julia_affine_ms=", round(alignment_ms, digits=4))
    println("  julia_affine_score=", alignment.score)
    println("  julia_affine_identity=", round(alignment.identity, digits=6))
    print(python_output)
end

main()
