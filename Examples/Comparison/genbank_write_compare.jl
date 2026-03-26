pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function write_sample_genbank(path::AbstractString, n::Integer=500)
    sequence = repeat("ACGT", 20)
    open(path, "w") do io
        for record_index in 1:n
            locus = record_index == 1 ? "WRIT001" : "WRIT$(lpad(1000 + record_index - 1, 6, '0'))"
            accession = record_index == 1 ? "WRIT001" : "WRIT$(lpad(1000 + record_index - 1, 6, '0'))"
            write(io, "LOCUS       $(locus)       80 bp    DNA             PLN       21-JUN-1999\n")
            write(io, "DEFINITION  Synthetic GenBank writer benchmark record with wrapped metadata fields.\n")
            write(io, "ACCESSION   $(accession)\n")
            write(io, "VERSION     $(accession).1\n")
            write(io, "KEYWORDS    synthetic benchmark record with canonical writer formatting\n")
            write(io, "SOURCE      Synthetic construct\n")
            write(io, "  ORGANISM  Synthetic construct\n")
            write(io, "            Artificial sequences.\n")
            write(io, "COMMENT     This benchmark exercises stricter GenBank serialization and wrapping.\n")
            write(io, "FEATURES             Location/Qualifiers\n")
            write(io, "     source          1..80\n")
            write(io, "                     /organism=\"Synthetic construct\"\n")
            write(io, "     gene            complement(5..20)\n")
            write(io, "                     /gene=\"geneA\"\n")
            write(io, "     CDS             join(1..8,41..48)\n")
            write(io, "                     /gene=\"geneB\"\n")
            write(io, "ORIGIN\n")
            write(io, "        1 $(sequence[1:10]) $(sequence[11:20]) $(sequence[21:30]) $(sequence[31:40]) $(sequence[41:50])\n")
            write(io, "       51 $(sequence[51:60]) $(sequence[61:70]) $(sequence[71:80])\n")
            write(io, "//\n")
        end
    end
end

function main()
    mktempdir() do dir
        input_path = joinpath(dir, "sample.gb")
        output_path = joinpath(dir, "written.gb")
        write_sample_genbank(input_path)

        records = BioToolkit.read_genbank(input_path)
        @assert length(records) == 500

        write_ms, _ = elapsed_ms() do
            BioToolkit.write_genbank(output_path, records)
        end

        roundtrip = BioToolkit.read_genbank(output_path)
        output_text = read(output_path, String)
        @assert roundtrip == records
        @assert occursin("KEYWORDS    synthetic benchmark record", output_text)
        @assert occursin("COMMENT     This benchmark exercises", output_text)

        python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "genbank_write_compare.py")) $input_path`, String)

        python_write_ms = parse(Float64, match(r"python_genbank_write_ms=(\d+(?:\.\d+)?)", python_output).captures[1])
        python_roundtrip_ok = match(r"python_genbank_roundtrip_ok=(true|false)", python_output).captures[1] == "true"
        python_output_contains_keywords = match(r"python_genbank_output_keywords=(true|false)", python_output).captures[1] == "true"

        @assert python_roundtrip_ok
        @assert python_output_contains_keywords

        println("GenBank writer benchmark")
        println("  records=", length(records))
        println("  julia_write_ms=", round(write_ms, digits=4))
        println("  python_write_ms=", round(python_write_ms, digits=4))
        print(python_output)
    end
end

main()
