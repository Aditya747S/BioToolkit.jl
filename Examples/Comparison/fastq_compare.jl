pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function write_sample_fastq(path::AbstractString, record_count::Integer, read_length::Integer)
    bases = ('A', 'C', 'G', 'T')

    open(path, "w") do io
        for record_index in 1:record_count
            write(io, "@read$(record_index)\n")

            for base_index in 1:read_length
                base = bases[mod(record_index + base_index - 2, length(bases)) + 1]
                write(io, base)
            end
            write(io, "\n")

            write(io, "+\n")
            for _ in 1:read_length
                write(io, 'I')
            end
            write(io, "\n")
        end
    end
end

function main()
    record_count = 30_000
    read_length = 100

    mktempdir() do dir
        fastq_path = joinpath(dir, "sample.fastq")
        python_script = joinpath(@__DIR__, "fastq_compare.py")
        written_fastq_path = joinpath(dir, "written.fastq")

        write_sample_fastq(fastq_path, record_count, read_length)

        BioToolkit.read_fastq(fastq_path)

        julia_ms, julia_records = elapsed_ms() do
            BioToolkit.read_fastq(fastq_path)
        end

        write_ms, _ = elapsed_ms() do
            BioToolkit.write_fastq(written_fastq_path, julia_records)
        end

        python_output = read(`conda run -n general --no-capture-output python $python_script $fastq_path`, String)

        total_bases = sum(length(record.sequence) for record in julia_records)

        println("FASTQ parse/write benchmark")
        println("  records: ", record_count)
        println("  read_length: ", read_length)
        println("  julia_fastq_ms=", round(julia_ms, digits=4))
        println("  julia_fastq_records=", length(julia_records))
        println("  julia_fastq_bases=", total_bases)
        println("  julia_fastq_write_ms=", round(write_ms, digits=4))
        print(python_output)
    end
end

main()
