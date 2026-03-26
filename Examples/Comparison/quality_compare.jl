pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_records(record_count::Integer, read_length::Integer)
    bases = ('A', 'C', 'G', 'T')
    high_quality = repeat("I", read_length - 12)
    low_quality = repeat("!", 12)
    records = FastqRecord[]

    for index in 1:record_count
        sequence = Vector{Char}(undef, read_length)
        for position in 1:read_length
            sequence[position] = bases[mod(index + position - 2, length(bases)) + 1]
        end
        quality = index % 5 == 0 ? high_quality * low_quality : repeat("I", read_length)
        push!(records, FastqRecord("read$(index)", "read$(index)", String(sequence), quality))
    end

    return records
end

function main()
    record_count = 20_000
    read_length = 100
    records = build_records(record_count, read_length)
    score_vector = collect(0:39)

    BioToolkit.phred_scores(first(records))
    BioToolkit.mean_quality(first(records))
    BioToolkit.phred_string(score_vector)
    BioToolkit.trim_low_quality(first(records); window=8, threshold=25)

    decode_ms, decode_total = elapsed_ms() do
        total = 0
        for record_index in eachindex(records)
            total += length(BioToolkit.phred_scores(records[record_index]))
        end
        total
    end

    mean_ms, mean_total = elapsed_ms() do
        total = 0.0
        for record_index in eachindex(records)
            total += BioToolkit.mean_quality(records[record_index])
        end
        total
    end

    string_ms, string_value = elapsed_ms() do
        BioToolkit.phred_string(score_vector)
    end

    trim_ms, trim_total = elapsed_ms() do
        total = 0
        for record_index in eachindex(records)
            trimmed = BioToolkit.trim_low_quality(records[record_index]; window=8, threshold=25)
            total += trimmed === nothing ? 0 : length(trimmed.sequence)
        end
        total
    end

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "quality_compare.py"))`, String)

    python_decode_total = parse(Int, match(r"python_decode_total=(\d+)", python_output).captures[1])
    python_mean_total = parse(Float64, match(r"python_mean_total=([0-9.]+)", python_output).captures[1])
    python_phred_string_len = parse(Int, match(r"python_phred_string_len=(\d+)", python_output).captures[1])
    python_trim_total = parse(Int, match(r"python_trim_total=(\d+)", python_output).captures[1])
    biopython_decode_total = parse(Int, match(r"biopython_decode_total=(\d+)", python_output).captures[1])
    biopython_mean_total = parse(Float64, match(r"biopython_mean_total=([0-9.]+)", python_output).captures[1])
    biopython_trim_total = parse(Int, match(r"biopython_trim_total=(\d+)", python_output).captures[1])

    @assert decode_total == python_decode_total == biopython_decode_total
    @assert isapprox(round(mean_total, digits=4), python_mean_total; atol=1e-12)
    @assert isapprox(round(mean_total, digits=4), biopython_mean_total; atol=1e-12)
    @assert length(string_value) == python_phred_string_len
    @assert trim_total == python_trim_total == biopython_trim_total
    @assert length(string_value) == length(score_vector)

    println("Quality toolkit benchmark")
    println("  records=", record_count)
    println("  read_length=", read_length)
    println("  julia_decode_ms=", round(decode_ms, digits=4))
    println("  julia_decode_total=", decode_total)
    println("  julia_mean_ms=", round(mean_ms, digits=4))
    println("  julia_mean_total=", round(mean_total, digits=4))
    println("  julia_phred_string_ms=", round(string_ms, digits=4))
    println("  julia_phred_string_len=", length(string_value))
    println("  julia_trim_ms=", round(trim_ms, digits=4))
    println("  julia_trim_total=", trim_total)
    print(python_output)
end

main()
