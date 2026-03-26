pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
using Arrow
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function write_sample_genbank(path::AbstractString)
    open(path, "w") do io
        for record_index in 1:500
            locus = record_index == 1 ? "SCU49845" : "SCU$(lpad(49845 + record_index - 1, 5, '0'))"
            accession = record_index == 1 ? "U49845" : "U$(lpad(49845 + record_index - 1, 5, '0'))"
            write(io, "LOCUS       $(locus)       50 bp    DNA             PLN       21-JUN-1999\n")
            write(io, "DEFINITION  Saccharomyces cerevisiae TCP1-beta gene.\n")
            write(io, "ACCESSION   $(accession)\n")
            write(io, "VERSION     $(accession).1\n")
            write(io, "SOURCE      Saccharomyces cerevisiae (baker's yeast)\n")
            write(io, "  ORGANISM  Saccharomyces cerevisiae\n")
            write(io, "            Eukaryota; Fungi; Ascomycota; Saccharomycetes.\n")
            write(io, "FEATURES             Location/Qualifiers\n")
            write(io, "     source          1..50\n")
            write(io, "                     /organism=\"Saccharomyces cerevisiae\"\n")
            write(io, "     gene            1..206\n")
            write(io, "                     /gene=\"TCP1\"\n")
            write(io, "ORIGIN\n")
            write(io, "        1 atgaccaatg cctgctgctg ctgctgctgc tgctgctgct gctgctgctg\n")
            write(io, "//\n")
        end
    end
end

function main()
    mktempdir() do dir
        input_path = joinpath(dir, "sample.gb")
        arrow_path = joinpath(dir, "sample.arrow")
        warm_arrow_path = joinpath(dir, "warm.arrow")
        write_sample_genbank(input_path)

        BioToolkit.read_genbank(input_path)
        BioToolkit.ingest_genbank(input_path, warm_arrow_path; chunk_size=500)

        parse_ms, records = elapsed_ms() do
            BioToolkit.read_genbank(input_path)
        end

        ingest_ms, _ = elapsed_ms() do
            BioToolkit.ingest_genbank(input_path, arrow_path; chunk_size=500)
        end

        python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "genbank_compare.py")) $input_path`, String)

        python_locus = match(r"python_genbank_locus=(.+)", python_output).captures[1]
        python_accession = match(r"python_genbank_accession=(.+)", python_output).captures[1]
        python_version = match(r"python_genbank_version=(.+)", python_output).captures[1]
        python_feature_count = parse(Int, match(r"python_genbank_feature_count=(\d+)", python_output).captures[1])
        python_sequence_length = parse(Int, match(r"python_genbank_sequence_length=(\d+)", python_output).captures[1])
        python_gene = match(r"python_genbank_gene=(.+)", python_output).captures[1]

        @assert length(records) == 500
        @assert records[1].locus == python_locus
        @assert records[1].accession == python_accession
        @assert records[1].version == python_version
        @assert length(records[1].features) == python_feature_count
        @assert length(records[1].sequence) == python_sequence_length
        @assert records[1].features[2].qualifiers["gene"][1] == python_gene

        table = Arrow.Table(arrow_path)
        @assert length(table.locus) == 500
        @assert table.locus[1] == records[1].locus
        @assert table.feature_count[1] == Int32(length(records[1].features))

        println("GenBank benchmark")
        println("  records=", length(records))
        println("  julia_parse_ms=", round(parse_ms, digits=4))
        println("  julia_ingest_ms=", round(ingest_ms, digits=4))
        println("  julia_locus=", records[1].locus)
        println("  julia_feature_count=", length(records[1].features))
        println("  julia_sequence_length=", length(records[1].sequence))
        print(python_output)
    end
end

main()
