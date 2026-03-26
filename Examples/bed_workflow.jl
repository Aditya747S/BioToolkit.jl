using BioToolkit
using Arrow
using Tables

function write_sample_bed(path::AbstractString)
    open(path, "w") do io
        write(io, "chr1\t0\t10\n")
        write(io, "chr1\t10\t20\n")
        write(io, "chr1\t20\t40\n")
        write(io, "chr2\t5\t15\n")
    end
end

mktempdir() do dir
    input_path = joinpath(dir, "sample.bed")
    output_path = joinpath(dir, "sample.arrow")

    write_sample_bed(input_path)
    BioToolkit.ingest_bed(input_path, output_path; chunk_size=2)

    table = BioToolkit.load_arrow_table(output_path)
    total_interval_length = sum(table.stop .- table.start)
    window_coverage = BioToolkit.window_coverage(table, "chr1", 10)
    partitions = length(collect(Tables.partitions(Arrow.Stream(output_path))))

    println("BED example")
    println("  loaded rows: ", Tables.rowcount(table))
    println("  total interval length: ", total_interval_length)
    println("  chr1 window coverage: ", sort(collect(window_coverage)))
    println("  Arrow partitions written: ", partitions)
    println("  first interval: ", (table.chrom[1], table.start[1], table.stop[1]))
end
