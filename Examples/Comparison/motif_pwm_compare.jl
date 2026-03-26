pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function build_sequences(count::Integer, motif_length::Integer)
    bases = ('A', 'C', 'G', 'T')
    template = collect("ACGTACGTACGT")
    motif_length == length(template) || throw(ArgumentError("motif_length must equal $(length(template))"))
    records = SeqRecordLite[]

    for index in 1:count
        buffer = copy(template)
        if mod(index, 7) == 0
            buffer[4] = bases[mod(index, length(bases)) + 1]
        end
        if mod(index, 11) == 0
            buffer[9] = bases[mod(index + 1, length(bases)) + 1]
        end
        sequence = String(buffer)
        push!(records, SeqRecordLite(sequence; identifier="motif$(index)"))
    end

    return records
end

function main()
    record_count = 20_000
    motif_length = 12
    records = build_sequences(record_count, motif_length)
    jaspar_data = ">MA0001.1 Example motif\nA [ 3 1 0 0 ]\nC [ 0 2 4 0 ]\nG [ 1 1 0 5 ]\nT [ 0 0 0 0 ]\n"
    jaspar_path = tempname()
    open(jaspar_path, "w") do io
        write(io, jaspar_data)
    end

    BioToolkit.motif_counts(records)
    BioToolkit.motif_pwm(records; pseudocount=0.0)
    BioToolkit.read_jaspar(jaspar_path)

    counts_ms, counts = elapsed_ms() do
        BioToolkit.motif_counts(records)
    end

    pwm_ms, pwm = elapsed_ms() do
        BioToolkit.motif_pwm(counts; pseudocount=0.0)
    end

    jaspar_ms, jaspar_profiles = elapsed_ms() do
        BioToolkit.read_jaspar(jaspar_path)
    end

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "motif_pwm_compare.py"))`, String)

    python_consensus = match(r"python_consensus=(.+)", python_output).captures[1]
    python_pwm_first = parse(Float64, match(r"python_pwm_first=([0-9.]+)", python_output).captures[1])
    python_jaspar_name = match(r"python_jaspar_name=(.+)", python_output).captures[1]
    python_jaspar_consensus = match(r"python_jaspar_consensus=(.+)", python_output).captures[1]
    python_jaspar_pwm_first = parse(Float64, match(r"python_jaspar_pwm_first=([0-9.]+)", python_output).captures[1])
    @assert motif_consensus(counts) == python_consensus
    @assert isapprox(pwm.values[1, 1], python_pwm_first; atol=1e-8)
    @assert startswith(jaspar_profiles[1].name, python_jaspar_name)
    @assert motif_consensus(jaspar_profiles[1]) == python_jaspar_consensus
    @assert isapprox(jaspar_profiles[1].pwm.values[1, 1], python_jaspar_pwm_first; atol=1e-8)

    println("Motif PWM benchmark")
    println("  records=", record_count)
    println("  motif_length=", motif_length)
    println("  julia_counts_ms=", round(counts_ms, digits=4))
    println("  julia_pwm_ms=", round(pwm_ms, digits=4))
    println("  julia_jaspar_ms=", round(jaspar_ms, digits=4))
    println("  julia_consensus=", motif_consensus(counts))
    println("  julia_pwm_first=", round(pwm.values[1, 1], digits=6))
    println("  julia_jaspar_name=", jaspar_profiles[1].name)
    println("  julia_jaspar_consensus=", motif_consensus(jaspar_profiles[1]))
    println("  julia_jaspar_pwm_first=", round(jaspar_profiles[1].pwm.values[1, 1], digits=6))
    print(python_output)
end

main()
