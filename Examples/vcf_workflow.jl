using BioToolkit
using Arrow
using Tables

function write_sample_vcf(path::AbstractString)
    open(path, "w") do io
        write(io, "##fileformat=VCFv4.2\n")
        write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")
        write(io, "chr1\t10\trs1\tA\tG\t99.5\n")
        write(io, "chr1\t20\trs2\tC\tT\t80.0\n")
        write(io, "chr2\t30\trs3\tG\tA\t70.0\n")
        write(io, "chr1\t40\trs4\tT\tC\t60.0\n")
    end
end

mktempdir() do dir
    input_path = joinpath(dir, "sample.vcf")
    output_path = joinpath(dir, "sample.arrow")

    write_sample_vcf(input_path)
    BioToolkit.ingest_vcf(input_path, output_path; chunk_size=2)

    table = BioToolkit.load_arrow_table(output_path)
    subset = BioToolkit.filter_region(table, "chr1", 0, 100)
    hist = BioToolkit.coverage_histogram(table, "chr1", 10)

    sorted_table = (
        chrom = ["chr1", "chr1", "chr1", "chr2"],
        pos = Int32[10, 20, 40, 30],
        id = ["rs1", "rs2", "rs4", "rs3"],
        ref = ["A", "C", "T", "G"],
        alt = ["G", "T", "C", "A"],
        qual = Union{Missing,Float32}[Float32(99.5), Float32(80.0), Float32(60.0), Float32(70.0)],
    )
    sorted_subset = BioToolkit.filter_region(sorted_table, "chr1", 15, 45; sorted=true)
    partitions = length(collect(Tables.partitions(Arrow.Stream(output_path))))

    println("VCF example")
    println("  loaded rows: ", Tables.rowcount(table))
    println("  chr1 subset rows: ", length(subset.pos))
    println("  sorted chr1 subset rows: ", length(sorted_subset.pos))
    println("  Arrow partitions written: ", partitions)
    println("  chr1 histogram bins: ", sort(collect(keys(hist))))
end
