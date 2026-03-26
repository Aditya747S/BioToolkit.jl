pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_aligned_records(record_count::Integer=500, width::Integer=120)
    template = repeat("ACGT", cld(width, 4))[1:width]
    bases = ('A', 'C', 'G', 'T')
    records = SeqRecordLite[]

    for record_index in 1:record_count
        chars = collect(template)
        for position in 1:width
            if mod(record_index + position, 17) == 0
                chars[position] = '-'
            elseif mod(record_index + position, 11) == 0
                chars[position] = bases[mod(record_index + position - 2, length(bases)) + 1]
            end
        end
        push!(records, SeqRecordLite(String(chars); identifier="seq$(record_index)", description="aligned_sequence"))
    end

    return records
end

function assert_record_sequences(expected, actual)
    @assert length(expected) == length(actual)
    for (expected_record, actual_record) in zip(expected, actual)
        @assert expected_record.sequence == actual_record.sequence
        @assert expected_record.identifier == actual_record.identifier
    end
end

function main()
    records = build_aligned_records()
    temp_dir = mktempdir()
    try
        alignment = MultipleSequenceAlignment(records)

        roundtrip_ms, roundtrip = elapsed_ms() do
            path = joinpath(temp_dir, "alignment.phy")
            BioToolkit.write_alignment(path, alignment, "phylip")
            BioToolkit.read_alignment(path, "phylip")
        end

        @assert length(alignment) == 500
        @assert get_alignment_length(alignment) == 120
        assert_record_sequences(alignment.records, roundtrip.records)

        python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "phylip_compare.py"))`, String)

        python_roundtrip_ms = parse(Float64, match(r"python_phylip_roundtrip_ms=([0-9.]+)", python_output).captures[1])
        python_roundtrip_ok = match(r"python_phylip_roundtrip_ok=(true|false)", python_output).captures[1] == "true"

        @assert python_roundtrip_ok

        println("PHYLIP round-trip benchmark")
        println("  record_count=", length(records))
        println("  julia_roundtrip_ms=", round(roundtrip_ms, digits=4))
        println("  python_roundtrip_ms=", round(python_roundtrip_ms, digits=4))
        print(python_output)
    finally
        rm(temp_dir; recursive=true, force=true)
    end
end

main()
