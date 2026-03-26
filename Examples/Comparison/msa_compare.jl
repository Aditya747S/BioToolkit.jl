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
    end
end

function main()
    records = build_aligned_records()
    temp_dir = mktempdir()
    try

    construct_ms, alignment = elapsed_ms() do
        MultipleSequenceAlignment(records; annotations=Dict(:source => "synthetic"), column_annotations=Dict(:clustal_consensus => repeat("*", 120), :emboss_consensus => repeat("|", 120), :secondary_structure => repeat(".", 120)))
    end

    column_ms, column = elapsed_ms() do
        alignment[:, 25]
    end

    row_ms, row = elapsed_ms() do
        alignment[250]
    end

    appended = copy(alignment)
    append_ms, _ = elapsed_ms() do
        push!(appended, SeqRecordLite(repeat("ACGT", 30); identifier="seq501", description="aligned sequence"))
    end

    frequency_ms, frequency = elapsed_ms() do
        BioToolkit.alignment_column_frequencies(alignment, 25)
    end

    fasta_roundtrip_ms, fasta_roundtrip = elapsed_ms() do
        path = joinpath(temp_dir, "alignment.fasta")
        BioToolkit.write_alignment(path, alignment, "fasta")
        BioToolkit.read_alignment(path, "fasta")
    end

    clustal_roundtrip_ms, clustal_roundtrip = elapsed_ms() do
        path = joinpath(temp_dir, "alignment.aln")
        BioToolkit.write_alignment(path, alignment, "clustal")
        BioToolkit.read_alignment(path, "clustal")
    end

    stockholm_roundtrip_ms, stockholm_roundtrip = elapsed_ms() do
        path = joinpath(temp_dir, "alignment.sto")
        BioToolkit.write_alignment(path, alignment, "stockholm")
        BioToolkit.read_alignment(path, "stockholm")
    end

    @assert length(alignment) == 500
    @assert get_alignment_length(alignment) == 120
    @assert length(column) == 500
    @assert row.identifier == "seq250"
    @assert length(appended) == 501
    @assert frequency['A'] > 0.0
    assert_record_sequences(alignment.records, fasta_roundtrip.records)
    assert_record_sequences(alignment.records, clustal_roundtrip.records)
    assert_record_sequences(alignment.records, stockholm_roundtrip.records)
    @assert all(expected.identifier == actual.identifier for (expected, actual) in zip(alignment.records, fasta_roundtrip.records))
    @assert all(expected.identifier == actual.identifier for (expected, actual) in zip(alignment.records, clustal_roundtrip.records))
    @assert all(expected.identifier == actual.identifier for (expected, actual) in zip(alignment.records, stockholm_roundtrip.records))
    @assert isempty(fasta_roundtrip.annotations)
    @assert isempty(clustal_roundtrip.annotations)
    @assert isempty(stockholm_roundtrip.annotations)
    @assert !haskey(fasta_roundtrip.column_annotations, :clustal_consensus)
    @assert clustal_roundtrip.column_annotations[:clustal_consensus] == repeat("*", 120)
    @assert length(clustal_roundtrip.column_annotations) == 1
    @assert stockholm_roundtrip.column_annotations[:clustal_consensus] == repeat("*", 120)
    @assert stockholm_roundtrip.column_annotations[:emboss_consensus] == repeat("|", 120)
    @assert stockholm_roundtrip.column_annotations[:secondary_structure] == repeat(".", 120)

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "msa_compare.py"))`, String)

    python_construct = parse(Float64, match(r"python_msa_construct_ms=([0-9.]+)", python_output).captures[1])
    python_column = parse(Float64, match(r"python_msa_column_ms=([0-9.]+)", python_output).captures[1])
    python_row = parse(Float64, match(r"python_msa_row_ms=([0-9.]+)", python_output).captures[1])
    python_append = parse(Float64, match(r"python_msa_append_ms=([0-9.]+)", python_output).captures[1])
    python_frequency = parse(Float64, match(r"python_msa_frequency_ms=([0-9.]+)", python_output).captures[1])
    python_fasta_roundtrip = parse(Float64, match(r"python_msa_fasta_roundtrip_ms=([0-9.]+)", python_output).captures[1])
    python_clustal_roundtrip = parse(Float64, match(r"python_msa_clustal_roundtrip_ms=([0-9.]+)", python_output).captures[1])
    python_stockholm_roundtrip = parse(Float64, match(r"python_msa_stockholm_roundtrip_ms=([0-9.]+)", python_output).captures[1])
    python_rows = parse(Int, match(r"python_msa_rows=(\d+)", python_output).captures[1])
    python_columns = parse(Int, match(r"python_msa_columns=(\d+)", python_output).captures[1])

    @assert python_rows == 500
    @assert python_columns == 120
    @assert occursin("python_msa_fasta_roundtrip_ok=1", python_output)
    @assert occursin("python_msa_clustal_roundtrip_ok=1", python_output)
    @assert occursin("python_msa_stockholm_roundtrip_ok=1", python_output)

    println("Multiple sequence alignment benchmark")
    println("  record_count=", length(records))
    println("  column_index=25")
    println("  julia_msa_construct_ms=", round(construct_ms, digits=4))
    println("  julia_msa_column_ms=", round(column_ms, digits=4))
    println("  julia_msa_row_ms=", round(row_ms, digits=4))
    println("  julia_msa_append_ms=", round(append_ms, digits=4))
    println("  julia_msa_frequency_ms=", round(frequency_ms, digits=4))
    println("  julia_msa_fasta_roundtrip_ms=", round(fasta_roundtrip_ms, digits=4))
    println("  julia_msa_clustal_roundtrip_ms=", round(clustal_roundtrip_ms, digits=4))
    println("  julia_msa_stockholm_roundtrip_ms=", round(stockholm_roundtrip_ms, digits=4))
    println("  julia_msa_rows=", length(alignment))
    println("  julia_msa_columns=", get_alignment_length(alignment))
    print(python_output)
    finally
        rm(temp_dir; recursive=true, force=true)
    end
end

main()
