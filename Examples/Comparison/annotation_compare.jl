pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function write_sample_genbank(path::AbstractString)
    sequence = repeat("ACGT", 20)
    open(path, "w") do io
        write(io, "LOCUS       ANN00001       80 bp    DNA             PLN       21-JUN-1999\n")
        write(io, "DEFINITION  Synthetic annotated sequence.\n")
        write(io, "ACCESSION   ANN00001\n")
        write(io, "VERSION     ANN00001.1\n")
        write(io, "SOURCE      Synthetic construct\n")
        write(io, "  ORGANISM  Synthetic construct\n")
        write(io, "            Artificial sequences.\n")
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

function main()
    mktempdir() do dir
        input_path = joinpath(dir, "sample.gb")
        write_sample_genbank(input_path)

        BioToolkit.read_genbank(input_path)
        BioToolkit.annotate_genbank_record(first(BioToolkit.read_genbank(input_path)))

        parse_ms, annotated = elapsed_ms() do
            BioToolkit.annotate_genbank_record(first(BioToolkit.read_genbank(input_path)))
        end

        python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "annotation_compare.py")) $input_path`, String)

        python_feature_count = parse(Int, match(r"python_annotation_feature_count=(\d+)", python_output).captures[1])
        python_source_type = match(r"python_annotation_source_type=(.+)", python_output).captures[1]
        python_source_extract = match(r"python_annotation_source_extract=(.+)", python_output).captures[1]
        python_feature1_type = match(r"python_annotation_feature1_type=(.+)", python_output).captures[1]
        python_feature1_strand = parse(Int, match(r"python_annotation_feature1_strand=(-?\d+)", python_output).captures[1])
        python_feature1_extract = match(r"python_annotation_feature1_extract=(.+)", python_output).captures[1]
        python_feature2_type = match(r"python_annotation_feature2_type=(.+)", python_output).captures[1]
        python_feature2_strand = parse(Int, match(r"python_annotation_feature2_strand=(-?\d+)", python_output).captures[1])
        python_feature2_extract = match(r"python_annotation_feature2_extract=(.+)", python_output).captures[1]
        python_gene = match(r"python_annotation_gene=(.+)", python_output).captures[1]

        @assert annotated.annotations[:locus] == "ANN00001"
        @assert length(annotated.features) == python_feature_count
        @assert annotated.features[1].feature_type == python_source_type
        @assert feature_sequence(annotated, annotated.features[1]) == python_source_extract
        @assert annotated.features[2].feature_type == python_feature1_type
        @assert Int(annotated.features[2].location.strand) == python_feature1_strand
        @assert feature_sequence(annotated, annotated.features[2]) == python_feature1_extract
        @assert annotated.features[3].feature_type == python_feature2_type
        @assert Int(annotated.features[3].location.strand) == python_feature2_strand
        @assert feature_sequence(annotated, annotated.features[3]) == python_feature2_extract
        @assert annotated.features[2].id == python_gene

        println("Sequence annotation benchmark")
        println("  parse_ms=", round(parse_ms, digits=4))
        println("  feature_count=", length(annotated.features))
        println("  source_type=", annotated.features[1].feature_type)
        println("  source_extract=", feature_sequence(annotated, annotated.features[1]))
        println("  feature1_type=", annotated.features[2].feature_type)
        println("  feature1_extract=", feature_sequence(annotated, annotated.features[2]))
        println("  feature2_type=", annotated.features[3].feature_type)
        println("  feature2_extract=", feature_sequence(annotated, annotated.features[3]))
        print(python_output)
    end
end

main()
