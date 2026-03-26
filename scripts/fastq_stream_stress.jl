#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using BioToolkit

function write_dummy_fastq(path::AbstractString, record_count::Integer; sequence::AbstractString="ACGTACGTACGTACGT", quality::AbstractString="IIIIIIIIIIIIIIII")
    open(path, "w") do io
        for index in 1:record_count
            write(io, "@seq$(index)\n")
            write(io, sequence, "\n")
            write(io, "+\n")
            write(io, quality, "\n")
        end
    end
    return path
end

function count_fastq_records(path::AbstractString)
    total = 0
    stream = BioToolkit.each_fastq_record(path)
    try
        for record in stream
            total += 1
            _ = record.sequence
        end
    finally
        close(stream)
    end
    return total
end

function main()
    output_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(pwd(), "dummy_10m.fastq")
    record_count = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10_000_000

    println("writing=", output_path)
    println("records=", record_count)
    write_dummy_fastq(output_path, record_count)
    GC.gc()
    counted = count_fastq_records(output_path)
    println("counted=", counted)
    return nothing
end

main()