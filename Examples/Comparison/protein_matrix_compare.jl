pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_sequences()
    left = repeat("MTEITAAMVKELRESTGAGMMDCKNALSETQHEWAY", 12)
    right_chars = collect(left)
    for index in 1:3:length(right_chars)
        if right_chars[index] == 'A'
            right_chars[index] = 'G'
        elseif right_chars[index] == 'G'
            right_chars[index] = 'A'
        elseif right_chars[index] == 'M'
            right_chars[index] = 'L'
        end
    end
    right = String(right_chars)
    return left, right
end

function parse_inventory(output::AbstractString)
    match_result = match(r"python_named_matrix_inventory=([A-Z0-9.,_]+)", output)
    match_result === nothing && throw(ArgumentError("missing Python named matrix inventory"))
    return split(match_result.captures[1], ',')
end

function main()
    left, right = build_sequences()
    left_record = SeqRecordLite(left; identifier="left", description="protein benchmark left sequence")
    right_record = SeqRecordLite(right; identifier="right", description="protein benchmark right sequence")

    names = BioToolkit.available_named_substitution_matrices()
    alignable_names = filter(name -> name != "GENETIC" && name != "SCHNEIDER" && name != "TRANS", names)
    load_ms, _ = elapsed_ms() do
        for name in alignable_names
            BioToolkit.named_substitution_matrix(name)
        end
    end

    blosum62_global_ms, blosum62_global = elapsed_ms() do
        BioToolkit.needleman_wunsch(left_record, right_record; substitution_matrix=BioToolkit.named_substitution_matrix(:BLOSUM62), gap_open=-10, gap_extend=-1)
    end
    blosum62_local_ms, blosum62_local = elapsed_ms() do
        BioToolkit.smith_waterman(left_record, right_record; substitution_matrix=BioToolkit.named_substitution_matrix(:BLOSUM62), gap_open=-10, gap_extend=-1)
    end
    pam250_global_ms, pam250_global = elapsed_ms() do
        BioToolkit.needleman_wunsch(left_record, right_record; substitution_matrix=BioToolkit.named_substitution_matrix(:PAM250), gap_open=-10, gap_extend=-1)
    end
    pam250_local_ms, pam250_local = elapsed_ms() do
        BioToolkit.smith_waterman(left_record, right_record; substitution_matrix=BioToolkit.named_substitution_matrix(:PAM250), gap_open=-10, gap_extend=-1)
    end

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "protein_matrix_compare.py"))`, String)
    python_inventory = parse_inventory(python_output)
    @assert sort(python_inventory) == names

    python_load_ms = parse(Float64, match(r"python_named_matrix_load_ms=([0-9.]+)", python_output).captures[1])
    python_blosum62_global_ms = parse(Float64, match(r"python_blosum62_global_ms=([0-9.]+)", python_output).captures[1])
    python_blosum62_global_score = parse(Float64, match(r"python_blosum62_global_score=(-?[0-9.]+)", python_output).captures[1])
    python_blosum62_global_identity = parse(Float64, match(r"python_blosum62_global_identity=([0-9.]+)", python_output).captures[1])
    python_blosum62_global_length = parse(Int, match(r"python_blosum62_global_length=(\d+)", python_output).captures[1])

    python_blosum62_local_ms = parse(Float64, match(r"python_blosum62_local_ms=([0-9.]+)", python_output).captures[1])
    python_blosum62_local_score = parse(Float64, match(r"python_blosum62_local_score=(-?[0-9.]+)", python_output).captures[1])
    python_blosum62_local_identity = parse(Float64, match(r"python_blosum62_local_identity=([0-9.]+)", python_output).captures[1])
    python_blosum62_local_length = parse(Int, match(r"python_blosum62_local_length=(\d+)", python_output).captures[1])

    python_pam250_global_ms = parse(Float64, match(r"python_pam250_global_ms=([0-9.]+)", python_output).captures[1])
    python_pam250_global_score = parse(Float64, match(r"python_pam250_global_score=(-?[0-9.]+)", python_output).captures[1])
    python_pam250_global_identity = parse(Float64, match(r"python_pam250_global_identity=([0-9.]+)", python_output).captures[1])
    python_pam250_global_length = parse(Int, match(r"python_pam250_global_length=(\d+)", python_output).captures[1])

    python_pam250_local_ms = parse(Float64, match(r"python_pam250_local_ms=([0-9.]+)", python_output).captures[1])
    python_pam250_local_score = parse(Float64, match(r"python_pam250_local_score=(-?[0-9.]+)", python_output).captures[1])
    python_pam250_local_identity = parse(Float64, match(r"python_pam250_local_identity=([0-9.]+)", python_output).captures[1])
    python_pam250_local_length = parse(Int, match(r"python_pam250_local_length=(\d+)", python_output).captures[1])

    @assert isapprox(Float64(blosum62_global.score), python_blosum62_global_score; atol=1e-8)
    @assert isapprox(round(blosum62_global.identity, digits=6), python_blosum62_global_identity; atol=1e-12)
    @assert length(blosum62_global.left) == python_blosum62_global_length
    @assert isapprox(Float64(blosum62_local.score), python_blosum62_local_score; atol=1e-8)
    @assert isapprox(round(blosum62_local.identity, digits=6), python_blosum62_local_identity; atol=1e-12)
    @assert length(blosum62_local.left) == python_blosum62_local_length
    @assert isapprox(Float64(pam250_global.score), python_pam250_global_score; atol=1e-8)
    @assert isapprox(round(pam250_global.identity, digits=6), python_pam250_global_identity; atol=1e-12)
    @assert length(pam250_global.left) == python_pam250_global_length
    @assert isapprox(Float64(pam250_local.score), python_pam250_local_score; atol=1e-8)
    @assert isapprox(round(pam250_local.identity, digits=6), python_pam250_local_identity; atol=1e-12)
    @assert length(pam250_local.left) == python_pam250_local_length

    println("Named substitution matrix protein benchmark")
    println("  matrix_count=", length(names))
    println("  alignable_matrix_count=", length(alignable_names))
    println("  julia_named_matrix_load_ms=", round(load_ms; digits=4))
    println("  julia_blosum62_global_ms=", round(blosum62_global_ms; digits=4))
    println("  julia_blosum62_local_ms=", round(blosum62_local_ms; digits=4))
    println("  julia_pam250_global_ms=", round(pam250_global_ms; digits=4))
    println("  julia_pam250_local_ms=", round(pam250_local_ms; digits=4))
    println("  julia_blosum62_global_score=", blosum62_global.score)
    println("  julia_blosum62_global_identity=", round(blosum62_global.identity; digits=6))
    println("  julia_blosum62_local_score=", blosum62_local.score)
    println("  julia_blosum62_local_identity=", round(blosum62_local.identity; digits=6))
    println("  julia_pam250_global_score=", pam250_global.score)
    println("  julia_pam250_global_identity=", round(pam250_global.identity; digits=6))
    println("  julia_pam250_local_score=", pam250_local.score)
    println("  julia_pam250_local_identity=", round(pam250_local.identity; digits=6))
    print(python_output)
end

main()
