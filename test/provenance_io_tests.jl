using Test
using BioToolkit

function _write_text_file(path::AbstractString, text::AbstractString)
    open(path, "w") do io
        write(io, text)
    end
    return path
end

function _node_by_operation(ctx::ProvenanceContext, operation::AbstractString)
    matches = [node for node in values(ctx.nodes) if node.operation == operation]
    @test !isempty(matches)
    return first(matches)
end

@testset "Provenance I/O" begin
    @testset "FASTA and FASTQ" begin
        ctx = BioToolkit.ProvenanceContext()
        fasta_dir = mktempdir()
        fasta_input = _write_text_file(joinpath(fasta_dir, "input.fasta"), ">seq1\nACGTACGT\n")
        fasta_records = BioToolkit.read_fasta(fasta_input; prov_ctx=ctx)
        @test length(fasta_records) == 1
        @test BioToolkit.container_provenance_id(fasta_records[1]) !== nothing
        @test occursin("provenance=", BioToolkit.container_provenance_summary(fasta_records[1]))

        fasta_output = joinpath(mktempdir(), "roundtrip.fasta")
        BioToolkit.write_fasta(fasta_output, fasta_records; prov_ctx=ctx)

        fasta_read_node = _node_by_operation(ctx, "read_fasta")
        fasta_write_node = _node_by_operation(ctx, "write_fasta")
        @test fasta_read_node.parameters["record_count"] == 1
        @test fasta_write_node.parameters["record_count"] == 1
        @test haskey(fasta_write_node.parameters, "hash")
        @test fasta_write_node.parent_ids == [BioToolkit.container_provenance_id(fasta_records[1])]
        @test isfile(fasta_output)

        fastq_dir = mktempdir()
        fastq_input = _write_text_file(joinpath(fastq_dir, "input.fastq"), "@fq1\nACGT\n+\nIIII\n")
        fastq_records = BioToolkit.read_fastq(fastq_input; prov_ctx=ctx)
        @test length(fastq_records) == 1
        @test BioToolkit.container_provenance_id(fastq_records[1]) !== nothing

        fastq_output = joinpath(mktempdir(), "roundtrip.fastq")
        BioToolkit.write_fastq(fastq_output, fastq_records; prov_ctx=ctx)

        fastq_write_node = _node_by_operation(ctx, "write_fastq")
        @test fastq_write_node.parent_ids == [BioToolkit.container_provenance_id(fastq_records[1])]
        @test haskey(fastq_write_node.parameters, "hash")
        @test isfile(fastq_output)
    end

    @testset "VCF and GenBank" begin
        ctx = BioToolkit.ProvenanceContext()
        vcf_dir = mktempdir()
        vcf_input = _write_text_file(joinpath(vcf_dir, "input.vcf"), "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\nchr1\t42\trs1\tA\tG\t99.5\n")
        vcf_document = BioToolkit.read_vcf_document(vcf_input; prov_ctx=ctx)
        @test BioToolkit.container_provenance_id(vcf_document) !== nothing
        @test BioToolkit.container_provenance_id(vcf_document.records[1]) !== nothing

        vcf_output = joinpath(mktempdir(), "roundtrip.vcf")
        BioToolkit.write_vcf_document(vcf_output, vcf_document; prov_ctx=ctx)
        vcf_write_node = _node_by_operation(ctx, "write_vcf_document")
        @test length(vcf_write_node.parent_ids) == 2
        @test haskey(vcf_write_node.parameters, "hash")
        @test isfile(vcf_output)

        vcf_alias_output = joinpath(mktempdir(), "roundtrip_alias.vcf")
        BioToolkit.write_vcf(vcf_alias_output, vcf_document; prov_ctx=ctx)
        @test isfile(vcf_alias_output)

        gb_dir = mktempdir()
        gb_input = _write_text_file(joinpath(gb_dir, "input.gb"), "LOCUS       SCU49845       50 bp    DNA             PLN       21-JUN-1999\nDEFINITION  Synthetic GenBank record for provenance testing.\nACCESSION   U49845\nVERSION     U49845.1\nKEYWORDS    synthetic; test record\nSOURCE      Saccharomyces cerevisiae\n            Eukaryota; Fungi; Ascomycota.\nCOMMENT     This record exercises provenance round-tripping.\nFEATURES             Location/Qualifiers\n     source          1..50\n                     /organism=\"Saccharomyces cerevisiae\"\nORIGIN\n        1 atgaccaatg cctgctgctg ctgctgctgc tgctgctgct gctgctgctg\n//\n")
        gb_records = BioToolkit.read_genbank(gb_input; prov_ctx=ctx)
        @test length(gb_records) == 1
        @test BioToolkit.container_provenance_id(gb_records[1]) !== nothing

        gb_output = joinpath(mktempdir(), "roundtrip.gb")
        BioToolkit.write_genbank(gb_output, gb_records; prov_ctx=ctx)
        gb_write_node = _node_by_operation(ctx, "write_genbank")
        @test gb_write_node.parent_ids == [BioToolkit.container_provenance_id(gb_records[1])]
        @test haskey(gb_write_node.parameters, "hash")
        @test isfile(gb_output)
    end

    @testset "BED and GFF" begin
        ctx = BioToolkit.ProvenanceContext()
        bed_dir = mktempdir()
        bed_input = _write_text_file(joinpath(bed_dir, "input.bed"), "chr1\t10\t20\nchr1\t15\t25\n")
        bed_records = BioToolkit.read_bed(bed_input; prov_ctx=ctx)
        @test !isempty(bed_records)
        bed_output = joinpath(mktempdir(), "roundtrip.bed")
        BioToolkit.write_bed(bed_output, bed_records; prov_ctx=ctx)
        bed_read_node = _node_by_operation(ctx, "read_bed")
        bed_write_node = _node_by_operation(ctx, "write_bed")
        @test bed_read_node.parameters["record_count"] == length(bed_records)
        @test isempty(bed_write_node.parent_ids)
        @test haskey(bed_write_node.parameters, "hash")
        @test isfile(bed_output)

        gff_dir = mktempdir()
        gff_input = _write_text_file(joinpath(gff_dir, "input.gff"), "chr1\tBioToolkit\tgene\t1\t50\t.\t+\t.\tID=gene1;Name=gene1\n")
        gff_records = BioToolkit.read_gff(gff_input; prov_ctx=ctx)
        @test !isempty(gff_records)
        gff_output = joinpath(mktempdir(), "roundtrip.gff")
        BioToolkit.write_gff(gff_output, gff_records; prov_ctx=ctx)
        gff_read_node = _node_by_operation(ctx, "read_gff")
        gff_write_node = _node_by_operation(ctx, "write_gff")
        @test gff_read_node.parameters["record_count"] == length(gff_records)
        @test isempty(gff_write_node.parent_ids)
        @test haskey(gff_write_node.parameters, "hash")
        @test isfile(gff_output)
    end

    @testset "BAM and gene prediction" begin
        ctx = BioToolkit.ProvenanceContext()

        bam_dir = mktempdir()
        bam_path = joinpath(bam_dir, "sample.bam")
        header = BioToolkit.BamHeader([
            BioToolkit.BamReference("chr1", 1_000),
        ])
        records = [BioToolkit.BamRecord("read1", "chr1", 9, [BioToolkit.BamCigarOp(10, 'M')], "ACGTACGTAA"; quality="IIIIIIIIII")]
        bam = BioToolkit.BamFile(records; header=header)
        BioToolkit.write_bam(bam_path, bam; prov_ctx=ctx)

        bam_write_node = _node_by_operation(ctx, "write_bam")
        @test bam_write_node.parameters["record_count"] == 1
        @test bam_write_node.parameters["write_index"]
        @test isfile(bam_path)

        roundtrip = BioToolkit.read_bam(bam_path; prov_ctx=ctx)
        @test roundtrip.records == records
        bam_read_node = _node_by_operation(ctx, "read_bam")
        @test bam_read_node.parameters["record_count"] == 1

        seq_ctx = BioToolkit.ProvenanceContext()
        genes = BioToolkit.predict_genes_hmm(BioToolkit.DNASeq("ATATAAATTTTAAATATATATAAATTTTAAATATATATAAT"); p_coding=0.05, prov_ctx=seq_ctx)
        @test all(gene -> gene isa BioToolkit.GeneInterval, genes)
        @test _node_by_operation(seq_ctx, "predict_genes_hmm").parameters["gene_count"] == length(genes)

        starts = BioToolkit.find_start_codons(BioToolkit.DNASeq("AAATGAAATG"); prov_ctx=seq_ctx)
        @test !isempty(starts)
        @test _node_by_operation(seq_ctx, "find_start_codons").parameters["site_count"] == length(starts)
    end

    @testset "Annotation provenance" begin
        ctx = BioToolkit.ProvenanceContext()

        location = BioToolkit.parse_feature_location("join(1..4,9..12)"; prov_ctx=ctx)
        @test BioToolkit.feature_spans(location; prov_ctx=ctx) == [(1, 4), (9, 12)]

        gb_path = joinpath(mktempdir(), "input.gb")
        open(gb_path, "w") do io
            write(io, "LOCUS       LOCUS1       12 bp    DNA     PLN       01-JAN-2000\nDEFINITION  Synthetic record.\nACCESSION   ACC1\nVERSION     ACC1.1\nSOURCE      Synthetic\nFEATURES             Location/Qualifiers\n     CDS             1..4\n                     /gene=\"gene1\"\n     CDS             join(1..4,9..12)\n                     /gene=\"gene2\"\nORIGIN\n        1 acgtacgtacgt\n//\n")
        end
        record = first(BioToolkit.read_genbank(gb_path; prov_ctx=ctx))
        annotated = BioToolkit.annotate_genbank_record(record; prov_ctx=ctx)
        @test length(annotated.features) == 2
        @test BioToolkit.feature_table(annotated; prov_ctx=ctx) isa AbstractVector
        @test BioToolkit.select_features(annotated; feature_type="CDS", prov_ctx=ctx) |> length == 2

        variants = [BioToolkit.VariantTextRecord("chr1", 4, "var1", "G", "A", missing)]
        feature = BioToolkit.GffRecord("chr1", "src", "CDS", 1, 6, missing, "+", missing, "ID=gene1", Dict("ID" => ["gene1"]))
        annotations = BioToolkit.annotate_variants(variants, [feature]; reference_sequences=Dict("chr1" => "ATGGAA"), prov_ctx=ctx)
        @test annotations[1].gene == "gene1"
        @test _node_by_operation(ctx, "annotate_variants").parameters["annotation_count"] == length(annotations)
    end

    @testset "Provenance export" begin
        ctx = BioToolkit.ProvenanceContext()
        export_dir = mktempdir()
        input_path = _write_text_file(joinpath(export_dir, "input.fasta"), ">seq1\nACGT\n")
        records = BioToolkit.read_fasta(input_path; prov_ctx=ctx)
        output_path = joinpath(mktempdir(), "out.fasta")
        BioToolkit.write_fasta(output_path, records; prov_ctx=ctx)
        json_text = BioToolkit.export_provenance_json(ctx)
        @test occursin("write_fasta", json_text)
        @test occursin("read_fasta", json_text)
        @test occursin("wasGeneratedBy", json_text)
    end
end
