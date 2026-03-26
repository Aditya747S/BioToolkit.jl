pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_sequences(record_count::Integer=8, width::Integer=180)
    template = repeat("ACGT", cld(width, 4))[1:width]
    bases = ('A', 'C', 'G', 'T')
    records = SeqRecordLite[]

    for record_index in 1:record_count
        chars = collect(template)
        for position in 1:width
            if mod(record_index + position, 19) == 0
                chars[position] = bases[mod(record_index + position - 2, length(bases)) + 1]
            elseif mod(record_index + position, 13) == 0
                chars[position] = bases[mod(record_index + position, length(bases)) + 1]
            end
        end
        sequence = String(chars)
        push!(records, SeqRecordLite(sequence; identifier="seq$(record_index)", description="synthetic_sequence"))
    end

    return records
end

function resolve_binary(candidate::AbstractString)
    resolved = Sys.which(candidate)
    resolved !== nothing && return resolved

    home = homedir()
    for root in (
        get(ENV, "CONDA_PREFIX", ""),
        joinpath(home, "miniconda3", "envs", "general"),
        joinpath(home, ".conda", "envs", "general"),
        joinpath(home, "anaconda3", "envs", "general"),
    )
        isempty(root) && continue
        path = joinpath(root, "bin", candidate)
        isfile(path) && return path
    end

    throw(ArgumentError("unable to locate binary: $(candidate)"))
end

function assert_alignment(alignment, expected_count::Integer)
    @assert length(alignment) == expected_count
    @assert get_alignment_length(alignment) > 0
    @assert all(!isempty(record.sequence) for record in alignment.records)
    return nothing
end

function main()
    records = build_sequences()
    clustalo = resolve_binary("clustalo")
    muscle = resolve_binary("muscle")
    temp_dir = mktempdir()
    try
        julia_clustal_ms, julia_clustal = elapsed_ms() do
            BioToolkit.clustal_msa(records; executable=clustalo, output_format="fasta")
        end

        julia_muscle_ms, julia_muscle = elapsed_ms() do
            BioToolkit.muscle_msa(records; executable=muscle, output_format="fasta")
        end

        assert_alignment(julia_clustal, length(records))
        assert_alignment(julia_muscle, length(records))

        python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "clustal_muscle_compare.py")) $(clustalo) $(muscle)`, String)

        python_clustal_ms = parse(Float64, match(r"python_clustal_msa_ms=([0-9.]+)", python_output).captures[1])
        python_muscle_ms = parse(Float64, match(r"python_muscle_msa_ms=([0-9.]+)", python_output).captures[1])
        python_records = parse(Int, match(r"python_external_msa_records=(\d+)", python_output).captures[1])
        python_columns = parse(Int, match(r"python_external_msa_columns=(\d+)", python_output).captures[1])

        @assert python_records == length(records)
        @assert python_columns > 0
        @assert occursin("python_clustal_msa_ok=1", python_output)
        @assert occursin("python_muscle_msa_ok=1", python_output)

        println("External MSA benchmark")
        println("  record_count=", length(records))
        println("  sequence_length=180")
        println("  clustal_binary=", clustalo)
        println("  muscle_binary=", muscle)
        println("  julia_clustal_msa_ms=", round(julia_clustal_ms, digits=4))
        println("  julia_muscle_msa_ms=", round(julia_muscle_ms, digits=4))
        println("  julia_external_msa_records=", length(julia_clustal))
        println("  julia_external_msa_columns=", get_alignment_length(julia_clustal))
        print(python_output)
    finally
        rm(temp_dir; recursive=true, force=true)
    end
end

main()
