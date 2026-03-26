using Test
using BioToolkit
using Arrow
using Plots
using Turing
using Tables

include(joinpath(@__DIR__, "..", "ext", "BioToolkitTuringExt.jl"))
BioToolkitTuringExt.__init__()

const test_files = sort(filter(file -> endswith(file, ".jl") && file != "runtests.jl", readdir(@__DIR__)))

for file in test_files
    include(joinpath(@__DIR__, file))
end

@testset "BioToolkit" begin
    @testset "VCF parsing" begin
        record = BioToolkit.parse_vcf_record("chr1\t42\trs1\tA\tG\t99.5")
        @test record.chrom == "chr1"
        @test record.pos == 42
        @test record.id == "rs1"
        @test record.alt == "G"
        @test record.qual == Float32(99.5)

        compact = BioToolkit.compact_variant_event(record)
        @test compact.chrom == UInt16(1)
        @test compact.pos == 42
        @test compact.hasqual == true
        @test isbitstype(typeof(compact))

        mktempdir() do dir
            input_path = joinpath(dir, "sample.vcf")
            output_path = joinpath(dir, "written.vcf")

            open(input_path, "w") do io
                write(io, "##fileformat=VCFv4.2\n")
                write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")
                write(io, "chr1\t42\trs1\tA\tG\t99.5\n")
            end

            parsed = BioToolkit.read_vcf(input_path)
            @test parsed == [record]

            BioToolkit.write_vcf(output_path, parsed)
            @test BioToolkit.read_vcf(output_path) == parsed
        end
    end

    @testset "Sequence toolkit" begin
        sequence = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
        sequence_bytes = Vector{UInt8}(codeunits(sequence))

        @test BioToolkit.validate_dna(sequence) == true
        @test BioToolkit.count_nucleotides(sequence) == (A = 9, C = 8, G = 14, T = 8, N = 0)
        @test BioToolkit.transcribe_dna("ATGTTT") == "AUGUUU"
        @test BioToolkit.melting_temp("ATGC") == 12.0
        @test BioToolkit.melting_temp("AUGC"; rna=true) == 12.0
        @test BioToolkit.dna_molecular_weight("ATGC") ≈ 1253.80 atol=0.1
        @test BioToolkit.dna_molecular_weight("ATGC"; stranded=:double) > BioToolkit.dna_molecular_weight("ATGC")
        @test BioToolkit.reverse_complement("ATGC") == "GCAT"
        @test BioToolkit.reverse_complement(sequence_bytes) == BioToolkit.reverse_complement(sequence)
        @test round(BioToolkit.gc_content(sequence), digits=3) == 0.564
        @test BioToolkit.translate_dna(sequence; stop_at_stop=true) == "MAIVMGR"
        @test BioToolkit.translate_dna(sequence_bytes; stop_at_stop=true) == BioToolkit.translate_dna(sequence; stop_at_stop=true)
        @test !isempty(BioToolkit.find_orfs(sequence; min_aa=5))
        @test BioToolkit.kmer_frequency("ATAT", 2) == Dict("AT" => 2, "TA" => 1)
        @test BioToolkit.codon_usage("ATGGCCATG") == Dict("ATG" => 2, "GCC" => 1)
        @test BioToolkit.codon_usage_table("ATGGCCATG") == Dict("ATG" => 2 / 3, "GCC" => 1 / 3)
        @test BioToolkit.relative_codon_adaptiveness("ATGGCCATG")["ATG"] == 1.0
        @test BioToolkit.codon_adaptation_index("ATGGCCATG"; reference="ATGGCCATG") == 1.0
        @test BioToolkit.cai("ATGGCCATG"; reference="ATGGCCATG") == 1.0
    end

        @testset "BLAST parsing" begin
            blast_xml = """<?xml version="1.0"?>
<BlastOutput>
    <BlastOutput_program>blastn</BlastOutput_program>
    <BlastOutput_version>BLASTN 2.15.0+</BlastOutput_version>
    <BlastOutput_db>synthetic_db</BlastOutput_db>
    <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
    <BlastOutput_query-def>synthetic query</BlastOutput_query-def>
    <BlastOutput_query-len>12</BlastOutput_query-len>
    <BlastOutput_iterations>
        <Iteration>
            <Iteration_query-ID>Query_1</Iteration_query-ID>
            <Iteration_query-def>synthetic query</Iteration_query-def>
            <Iteration_query-len>12</Iteration_query-len>
            <Iteration_hits>
                <Hit>
                    <Hit_id>subject1</Hit_id>
                    <Hit_def>subject one</Hit_def>
                    <Hit_len>12</Hit_len>
                    <Hit_hsps>
                        <Hsp>
                            <Hsp_bit-score>42.0</Hsp_bit-score>
                            <Hsp_evalue>1e-10</Hsp_evalue>
                            <Hsp_query-from>2</Hsp_query-from>
                            <Hsp_query-to>9</Hsp_query-to>
                            <Hsp_hit-from>3</Hsp_hit-from>
                            <Hsp_hit-to>10</Hsp_hit-to>
                            <Hsp_identity>8</Hsp_identity>
                            <Hsp_positive>8</Hsp_positive>
                            <Hsp_align-len>8</Hsp_align-len>
                            <Hsp_qseq>ACGTACGT</Hsp_qseq>
                            <Hsp_hseq>ACGTACGT</Hsp_hseq>
                            <Hsp_midline>||||||||</Hsp_midline>
                        </Hsp>
                    </Hit_hsps>
                </Hit>
            </Iteration_hits>
        </Iteration>
    </BlastOutput_iterations>
</BlastOutput>
"""

                records = BioToolkit.parse_blast_xml(IOBuffer(blast_xml))
                @test length(records) == 1
                @test records[1].program == "blastn"
                @test records[1].database == "synthetic_db"
                @test records[1].query_id == "Query_1"
                @test records[1].query_description == "synthetic query"
                @test records[1].query_length == 12
                @test length(records[1].hits) == 1
                @test records[1].hits[1].id == "subject1"
                @test records[1].hits[1].description == "subject one"
                @test records[1].hits[1].length == 12
                @test length(records[1].hits[1].hsps) == 1
                @test records[1].hits[1].hsps[1].bit_score == 42.0
                @test records[1].hits[1].hsps[1].evalue == 1e-10
                @test records[1].hits[1].hsps[1].query_start == 2
                @test records[1].hits[1].hsps[1].hit_end == 10

                tabular = """query1\tsubject1\t100.0\t8\t0\t0\t2\t9\t3\t10\t1e-10\t42.0
query1\tsubject2\t87.5\t8\t1\t0\t5\t12\t1\t8\t2e-5\t31.0
query2\tsubject3\t75.0\t4\t1\t0\t1\t4\t4\t7\t0.01\t20.0
"""

                tabular_records = BioToolkit.parse_blast_tabular(IOBuffer(tabular))
                @test length(tabular_records) == 2
                @test tabular_records[1].query_id == "query1"
                @test length(tabular_records[1].hits) == 2
                @test tabular_records[1].hits[1].target_id in ("subject1", "subject2")
                @test tabular_records[2].query_id == "query2"
                @test tabular_records[2].hits[1].hsps[1].alignment_length == 4
        end

        @testset "MEDLINE parsing" begin
                medline_text = """PMID- 12345678
TI  - Synthetic PubMed title
AB  - First line of abstract.
            Second line of abstract.
FAU - Ada Lovelace
FAU - Alan Turing
JT  - Journal of Synthetic Biology
DP  - 2024 Jan
MH  - Bioinformatics
AID - 10.1234/example [doi]

PMID- 87654321
TI  - Another record
AB  - Short abstract.
JT  - Another Journal
DP  - 2023
"""

                medline_records = BioToolkit.parse_medline_text(medline_text)
                @test length(medline_records) == 2
                @test medline_records[1].pmid == "12345678"
                @test occursin("Second line", medline_records[1].abstract)
                @test medline_records[1].year == 2024
                @test medline_records[1].doi == "10.1234/example"
                @test medline_records[1].authors[1] == "Ada Lovelace"

                medline_xml = """<PubmedArticleSet>
    <PubmedArticle>
        <MedlineCitation>
            <PMID Version=\"1\">12345678</PMID>
            <Article>
                <ArticleTitle>Synthetic PubMed title</ArticleTitle>
                <Abstract><AbstractText>XML abstract text.</AbstractText></Abstract>
                <Journal><Title>Journal of Synthetic Biology</Title><JournalIssue><PubDate><Year>2024</Year></PubDate></JournalIssue></Journal>
                <AuthorList>
                    <Author><LastName>Lovelace</LastName><ForeName>Ada</ForeName></Author>
                </AuthorList>
            </Article>
            <MeshHeadingList><MeshHeading><DescriptorName>Bioinformatics</DescriptorName></MeshHeading></MeshHeadingList>
        </MedlineCitation>
        <PubmedData>
            <ArticleIdList><ArticleId IdType=\"doi\">10.1234/example</ArticleId></ArticleIdList>
        </PubmedData>
    </PubmedArticle>
</PubmedArticleSet>"""

                medline_xml_records = BioToolkit.parse_medline_xml(medline_xml)
                @test length(medline_xml_records) == 1
                @test medline_xml_records[1].pmid == "12345678"
                @test medline_xml_records[1].journal == "Journal of Synthetic Biology"
                @test medline_xml_records[1].doi == "10.1234/example"
                @test medline_xml_records[1].authors[1] == "Ada Lovelace"
        end

    @testset "FASTA indexing and distance" begin
        mktempdir() do dir
            fasta_path = joinpath(dir, "sample.fasta")

            open(fasta_path, "w") do io
                write(io, ">seq1\n")
                write(io, "ATGCATGC\n")
                write(io, "ATGCATGC\n")
            end

            index = BioToolkit.fasta_index(fasta_path)
            @test haskey(index, "seq1")
            @test index["seq1"].sequence_length == 16
            @test BioToolkit.fetch_fasta_sequence(fasta_path, index["seq1"], 5, 12) == "ATGCATGC"
            @test BioToolkit.hamming_distance("ATGC", "ATGT") == 1
        end
    end

    @testset "Filtering and histogramming" begin
        table = (
            chrom = ["chr1", "chr1", "chr2"],
            pos = Int32[10, 120, 30],
            id = ["a", "b", "c"],
            ref = ["A", "C", "G"],
            alt = ["T", "G", "A"],
            qual = Union{Missing,Float32}[missing, Float32(10), Float32(20)],
        )

        subset = BioToolkit.filter_region(table, "chr1", 0, 100)
        @test subset.pos == Int32[10]
        @test subset.chrom == ["chr1"]

        sorted_table = (
            chrom = ["chr1", "chr1", "chr1", "chr2"],
            pos = Int32[5, 10, 15, 30],
            id = ["a", "b", "c", "d"],
            ref = ["A", "C", "G", "T"],
            alt = ["T", "G", "A", "C"],
            qual = Union{Missing,Float32}[Float32(1), Float32(2), Float32(3), Float32(4)],
        )

        sorted_subset = BioToolkit.filter_region(sorted_table, "chr1", 6, 14; sorted=true)
        @test sorted_subset.pos == Int32[10]
        @test sorted_subset.id == ["b"]

        hist = BioToolkit.bin_positions(Int32[1, 2, 11, 12], 10)
        @test hist[0] == 2
        @test hist[1] == 2

        threaded_hist = BioToolkit.bin_positions(Int32[1, 2, 11, 12], 10; threaded=true)
        @test threaded_hist == hist
    end

    @testset "FASTA parsing" begin
        mktempdir() do dir
            fasta_path = joinpath(dir, "sample.fasta")

            open(fasta_path, "w") do io
                write(io, ">seq1\n")
                write(io, "ATGC\n")
                write(io, ">seq2\n")
                write(io, "GATTACA\n")
            end

            records = BioToolkit.read_fasta(fasta_path)
            @test records == [("seq1", "ATGC"), ("seq2", "GATTACA")]
        end
    end

    @testset "FASTQ parsing" begin
        mktempdir() do dir
            fastq_path = joinpath(dir, "sample.fastq")

            open(fastq_path, "w") do io
                write(io, "@seq1 description\n")
                write(io, "ACGT\n")
                write(io, "+\n")
                write(io, "!!!!\n")
                write(io, "@seq2\n")
                write(io, "GATTACA\n")
                write(io, "+\n")
                write(io, "HHHHHHH\n")
            end

            records = BioToolkit.read_fastq(fastq_path)
            @test length(records) == 2
            @test records[1].identifier == "seq1"
            @test records[1].description == "seq1 description"
            @test records[1].sequence == "ACGT"
            @test records[1].quality == "!!!!"
            @test records[2].identifier == "seq2"
            @test records[2].sequence == "GATTACA"
        end
    end

    @testset "FASTQ writing and records" begin
        mktempdir() do dir
            output_path = joinpath(dir, "written.fastq")
            fastq_record = BioToolkit.FastqRecord("seq1", "seq1 description", "ACGT", "IIII")
            seq_record = BioToolkit.SeqRecordLite(
                "GATTACA";
                identifier="seq2",
                description="seq2 description",
                annotations=Dict(:organism => "synthetic"),
                letter_annotations=Dict(:quality => "HHHHHHH"),
            )

            @test length(seq_record) == 7
            @test seq_record.annotations[:organism] == "synthetic"

            BioToolkit.write_fastq(output_path, [fastq_record, seq_record])
            read_back = BioToolkit.read_fastq(output_path)

            @test read_back[1] == fastq_record
            @test read_back[2].sequence == "GATTACA"
            @test read_back[2].quality == "HHHHHHH"

            @test BioToolkit._fastq_components(fastq_record) == ("seq1", "seq1 description", "ACGT", "IIII")
            @test BioToolkit._fastq_components(seq_record) == ("seq2", "seq2 description", "GATTACA", "HHHHHHH")
            @test_throws ArgumentError BioToolkit._fastq_components(BioToolkit.SeqRecordLite("ACGT"; identifier="bad"))
        end
    end

    @testset "Quality pipeline" begin
        high_quality = BioToolkit.FastqRecord("high", "high quality", "ACGTAGATCGGAAGAGC", "IIIIIIIIIIIIIIIII")
        low_quality = BioToolkit.FastqRecord("low", "low quality", "ACGT", "!!!!")

        @test BioToolkit.quality_filter(high_quality; min_mean_quality=30, min_base_quality=30)
        @test !BioToolkit.quality_filter(low_quality; min_mean_quality=30, min_base_quality=30)

        trimmed_fastq = BioToolkit.adapter_trim(high_quality; adapter="AGATCGGAAGAGC", min_overlap=8)
        @test trimmed_fastq.sequence == "ACGT"
        @test trimmed_fastq.quality == "IIII"

        high_seq = BioToolkit.SeqRecordLite(
            "ACGTAGATCGGAAGAGC";
            identifier="seq1",
            letter_annotations=Dict(:quality => "IIIIIIIIIIIIIIIII"),
        )
        trimmed_seq = BioToolkit.adapter_trim(high_seq; adapter="AGATCGGAAGAGC", min_overlap=8)
        @test trimmed_seq.sequence == "ACGT"
        @test trimmed_seq.letter_annotations[:quality] == "IIII"

        pipeline_input = [high_quality, low_quality]
        pipeline_output = BioToolkit.sequencing_pipeline(pipeline_input; adapter="AGATCGGAAGAGC", min_mean_quality=30, min_base_quality=30, min_length=4)
        @test length(pipeline_output) == 1
        @test pipeline_output[1].sequence == "ACGT"

        seq_pipeline_input = [high_seq]
        seq_pipeline_output = BioToolkit.sequencing_pipeline(seq_pipeline_input; adapter="AGATCGGAAGAGC", min_mean_quality=30, min_base_quality=30, min_length=4)
        @test length(seq_pipeline_output) == 1
        @test seq_pipeline_output[1].sequence == "ACGT"
    end

    @testset "Intermediate helpers" begin
        @test BioToolkit._normalize_external_alignment_format("fa") == "fasta"
        @test BioToolkit._normalize_external_alignment_format("CLUSTAL") == "clustal"
        @test BioToolkit._normalize_external_alignment_format("sto") == "stockholm"
        @test BioToolkit._alignment_output_extension("fasta") == "fasta"
        @test BioToolkit._alignment_output_extension("clustal") == "aln"
        @test BioToolkit._alignment_output_extension("stockholm") == "sto"
        @test_throws ArgumentError BioToolkit._normalize_external_alignment_format("invalid")

        empty_msa = BioToolkit.MultipleSequenceAlignment(BioToolkit.SeqRecordLite[])
        @test BioToolkit.get_alignment_length(empty_msa) == 0
        @test isempty(BioToolkit.alignment_column_counts(empty_msa, 1))
        @test isempty(BioToolkit.alignment_column_frequencies(empty_msa, 1))
        @test BioToolkit.consensus_sequence(empty_msa) == ""

        @test_throws ArgumentError BioToolkit.MultipleSequenceAlignment([
            BioToolkit.SeqRecordLite("ACG"; identifier="short"),
            BioToolkit.SeqRecordLite("ACGT"; identifier="long"),
        ])
        @test_throws ArgumentError BioToolkit.MultipleSequenceAlignment([
            BioToolkit.SeqRecordLite("ACGT"; identifier="seq1"),
            BioToolkit.SeqRecordLite("ACGT"; identifier="seq2"),
        ]; column_annotations=Dict(:bad => "***"))
        @test_throws ArgumentError BioToolkit.format_alignment(BioToolkit.MultipleSequenceAlignment(BioToolkit.SeqRecordLite[]), "unsupported")
    end

    @testset "Pairwise alignment" begin
        left = BioToolkit.SeqRecordLite("ACGT"; identifier="left")
        right = BioToolkit.SeqRecordLite("ACGT"; identifier="right")
        alignment = BioToolkit.pairwise_align(left, right; match=2, mismatch=-1, gap=-2)
        affine_alignment = BioToolkit.pairwise_align(left, right; match=2, mismatch=-1, gap_open=-2, gap_extend=-2)
        matrix = BioToolkit.substitution_matrix("ACGT"; match=3, mismatch=-2)
        named_matrix = BioToolkit.named_substitution_matrix(:BLOSUM62)
        matrix_alignment = BioToolkit.pairwise_align(left, right; substitution_matrix=matrix, gap=-2)

        @test alignment.left == "ACGT"
        @test alignment.right == "ACGT"
        @test alignment.score == 8
        @test alignment.identity == 1.0
        @test affine_alignment.score == alignment.score
        @test affine_alignment.identity == alignment.identity
        @test matrix_alignment.score == 12
        @test matrix_alignment.identity == 1.0
        @test BioToolkit.needleman_wunsch(left, right; match=2, mismatch=-1, gap=-2).score == alignment.score
        @test BioToolkit.smith_waterman("TTTAGCCTT", "AGC"; match=2, mismatch=-1, gap=-2).score == 6
        @test BioToolkit.pairwise_align("AA", "AA"; substitution_matrix=named_matrix, gap=-4).score == 8

        codon_matrix = BioToolkit.named_codon_substitution_matrix(:SCHNEIDER)
        codon_global = BioToolkit.pairwise_align_codons("ATGATG", "ATGATG"; substitution_matrix=codon_matrix, gap=-3)
        codon_local = BioToolkit.smith_waterman_codons("ATGATG", "ATGATG"; substitution_matrix=codon_matrix, gap=-3)
        @test codon_global.left == "ATGATG"
        @test codon_global.right == "ATGATG"
        @test codon_global.identity == 1.0
        @test codon_local.score == codon_global.score

        gapped_left = BioToolkit.SeqRecordLite("AAAAAA"; identifier="gapped_left")
        gapped_right = BioToolkit.SeqRecordLite("AA"; identifier="gapped_right")
        linear_gap = BioToolkit.pairwise_align(gapped_left, gapped_right; match=2, mismatch=-1, gap=-2)
        affine_gap = BioToolkit.pairwise_align(gapped_left, gapped_right; match=2, mismatch=-1, gap_open=-2, gap_extend=-2)

        @test affine_gap.score == linear_gap.score
        @test affine_gap.matches == linear_gap.matches
        @test affine_gap.identity == linear_gap.identity
    end

    @testset "Local alignment (Smith-Waterman)" begin
        left = "TTTAGCCTT"
        right = "AGC"
        local_res = BioToolkit.local_align(left, right; match=2, mismatch=-1, gap=-2)
        @test local_res.score == 6
        @test local_res.left == "AGC"
        @test local_res.right == "AGC"
        
        left3 = "AAAAAAGGGGGAAAAAA"
        right3 = "AAAAAAAAAAAA"
        # Match=2, mismatch=-3, gap_open=-5, gap_extend=-1
        # GGGGG is 5 chars gap. Affine gap = -5 + (-1*4) = -9
        # Mismatch = -3 * 5 = -15
        local_affine = BioToolkit.local_align(left3, right3; match=2, mismatch=-3, gap_open=-5, gap_extend=-1)
        @test local_affine.score == 15  # 12 A's match = 24. Minus 9 for gap = 15.
        
        # Linear gap = -3. 5 gaps = -15.
        # Local alignment drops below 0 during the 5 gaps (12 - 15 = -3 -> 0)
        # So it simply matches one of the 6-A blocks for a score of 12.
        local_linear = BioToolkit.local_align(left3, right3; match=2, mismatch=-3, gap=-3)
        @test local_linear.score == 12   
        
        @test local_affine.score > local_linear.score
    end

    @testset "Multiple sequence alignment" begin
        alignment = BioToolkit.MultipleSequenceAlignment([
            BioToolkit.SeqRecordLite("ACGT"; identifier="seq1"),
            BioToolkit.SeqRecordLite("A-GT"; identifier="seq2"),
            BioToolkit.SeqRecordLite("AC-T"; identifier="seq3"),
        ]; annotations=Dict(:source => "synthetic"), column_annotations=Dict(:clustal_consensus => "****", :emboss_consensus => "||||", :column_scores => [1, 2, 3, 4]))

        @test length(alignment) == 3
        @test BioToolkit.get_alignment_length(alignment) == 4
        @test alignment[1].identifier == "seq1"
        @test alignment[1, 1] == 'A'
        @test alignment[1, :].identifier == "seq1"
        @test alignment[:, 2] == "C-C"
        @test alignment[1:2].records[2].identifier == "seq2"
        @test alignment[1:2, 2] == "C-"
        @test BioToolkit.consensus_sequence(alignment) == "ACGT"
        @test alignment.column_annotations[:clustal_consensus] == "****"
        @test alignment.column_annotations[:emboss_consensus] == "||||"
        @test alignment.column_annotations[:column_scores] == [1, 2, 3, 4]
        @test BioToolkit.alignment_column_counts(alignment, 2) == Dict('C' => 2)
        @test BioToolkit.alignment_column_frequencies(alignment, 2) == Dict('C' => 1.0)
        @test BioToolkit.alignment_column_symbol(alignment, 1) == '*'
        @test BioToolkit.alignment_symbol_line(alignment; style=:emboss) == "||||"
        @test BioToolkit.clustal_consensus(alignment) == "****"
        @test BioToolkit.alignment_symbol_line(alignment; style=:clustal) == "****"
        @test startswith(BioToolkit.format_alignment(alignment, "fasta"), ">seq1\nACGT")
        @test occursin("CLUSTAL multiple sequence alignment", BioToolkit.format_alignment(alignment, "clustal"))
        @test occursin("|", BioToolkit.format_alignment(alignment, "clustal"; symbol_style=:emboss))
        @test occursin("emboss_consensus ||||", BioToolkit.format_alignment(alignment, "stockholm"))

        protein_alignment = BioToolkit.MultipleSequenceAlignment([
            BioToolkit.SeqRecordLite("AB"; identifier="p1"),
            BioToolkit.SeqRecordLite("BB"; identifier="p2"),
        ])
        scoring = BioToolkit.SubstitutionMatrix("AB"; match=1, mismatch=1)
        @test BioToolkit.alignment_column_symbol(protein_alignment, 1; scoring=scoring) == ':'
        @test BioToolkit.clustal_consensus(protein_alignment; scoring=scoring) == ":*"

        sliced = alignment[:, 2:3]
        @test sliced.column_annotations[:clustal_consensus] == "**"
        @test sliced.column_annotations[:column_scores] == [2, 3]

        extended = copy(alignment)
        push!(extended, BioToolkit.SeqRecordLite("ACGT"; identifier="seq4"))
        @test length(extended) == 4
        @test extended[4].identifier == "seq4"

        combined = alignment + BioToolkit.MultipleSequenceAlignment([
            BioToolkit.SeqRecordLite("GG"; identifier="seq1"),
            BioToolkit.SeqRecordLite("TT"; identifier="seq2"),
            BioToolkit.SeqRecordLite("AA"; identifier="seq3"),
        ]; annotations=Dict(:source => "synthetic"), column_annotations=Dict(:clustal_consensus => "**", :emboss_consensus => "||", :column_scores => [5, 6]))
        @test combined[1].sequence == "ACGTGG"
        @test combined.annotations[:source] == "synthetic"
        @test combined.column_annotations[:clustal_consensus] == "******"
        @test combined.column_annotations[:emboss_consensus] == "||||||"
        @test combined.column_annotations[:column_scores] == [1, 2, 3, 4, 5, 6]

        stockholm = BioToolkit.format_alignment(alignment, "stockholm")
        @test startswith(stockholm, "# STOCKHOLM 1.0")
        @test occursin("#=GC clustal_consensus ****", stockholm)
        @test occursin("#=GC emboss_consensus ||||", stockholm)
        @test occursin("#=GC column_scores 1234", stockholm)
        @test endswith(stockholm, "//")

        fasta_roundtrip = BioToolkit.read_alignment(IOBuffer(BioToolkit.format_alignment(alignment, "fasta")), "fasta")
        @test fasta_roundtrip.records[1].sequence == "ACGT"
        @test fasta_roundtrip.records[1].identifier == "seq1"

        clustal_roundtrip = BioToolkit.read_alignment(IOBuffer(BioToolkit.format_alignment(alignment, "clustal")), "clustal")
        @test clustal_roundtrip.records[2].sequence == "A-GT"
        @test clustal_roundtrip.column_annotations[:clustal_consensus] == "****"

        stockholm_roundtrip = BioToolkit.read_alignment(IOBuffer(stockholm), "stockholm")
        @test stockholm_roundtrip.records[3].sequence == "AC-T"
        @test stockholm_roundtrip.column_annotations[:clustal_consensus] == "****"
        @test stockholm_roundtrip.column_annotations[:emboss_consensus] == "||||"
        @test stockholm_roundtrip.column_annotations[:column_scores] == "1234"

        msf_text = BioToolkit.format_alignment(alignment, "msf")
        @test startswith(msf_text, "!!NA_MULTIPLE_ALIGNMENT")
        msf_roundtrip = BioToolkit.read_alignment(IOBuffer(msf_text), "msf")
        @test msf_roundtrip.records[1].sequence == "ACGT"
        @test msf_roundtrip.records[2].sequence == "A-GT"

        pir_text = BioToolkit.format_alignment(alignment, "pir")
        @test startswith(pir_text, ">P1;")
        pir_roundtrip = BioToolkit.read_alignment(IOBuffer(pir_text), "pir")
        @test pir_roundtrip.records[1].sequence == "ACGT"
        @test pir_roundtrip.records[2].sequence == "A-GT"

        nexus_text = BioToolkit.format_alignment(alignment, "nexus")
        @test startswith(nexus_text, "#NEXUS")
        @test occursin("Begin data;", nexus_text)
        nexus_roundtrip = BioToolkit.read_alignment(IOBuffer(nexus_text), "nexus")
        @test nexus_roundtrip.records[1].sequence == "ACGT"
        @test nexus_roundtrip.records[2].sequence == "A-GT"

        phylip_alignment = BioToolkit.MultipleSequenceAlignment([
            BioToolkit.SeqRecordLite("ACGT"; identifier="seq1"),
            BioToolkit.SeqRecordLite("A-GT"; identifier="seq2"),
        ])
        phylip_text = BioToolkit.format_alignment(phylip_alignment, "phylip")
        @test startswith(phylip_text, "2 4")
        phylip_roundtrip = BioToolkit.read_alignment(IOBuffer(phylip_text), "phylip")
        @test phylip_roundtrip.records[1].sequence == "ACGT"
        @test phylip_roundtrip.records[2].sequence == "A-GT"

        relaxed_phylip = BioToolkit.MultipleSequenceAlignment([
            BioToolkit.SeqRecordLite("ACGT"; identifier="long_identifier_1"),
            BioToolkit.SeqRecordLite("A-GT"; identifier="long_identifier_2"),
        ])
        relaxed_text = BioToolkit.format_alignment(relaxed_phylip, "phylip-relaxed")
        @test occursin("long_identifier_1", relaxed_text)
        relaxed_roundtrip = BioToolkit.read_alignment(IOBuffer(relaxed_text), "auto")
        @test relaxed_roundtrip.records[1].identifier == "long_identifier_1"
        @test relaxed_roundtrip.records[2].sequence == "A-GT"

        fasta_buffer = IOBuffer()
        BioToolkit.write_alignment(fasta_buffer, alignment, "fasta")
        @test startswith(String(take!(fasta_buffer)), ">seq1\nACGT")

        fasta_auto = BioToolkit.read_alignment(IOBuffer(BioToolkit.format_alignment(alignment, "fasta")), "auto")
        clustal_auto = BioToolkit.read_alignment(IOBuffer(BioToolkit.format_alignment(alignment, "clustal")), "auto")
        stockholm_auto = BioToolkit.read_alignment(IOBuffer(stockholm), "auto")
        @test fasta_auto.records[1].sequence == "ACGT"
        @test clustal_auto.column_annotations[:clustal_consensus] == "****"
        @test stockholm_auto.column_annotations[:column_scores] == "1234"

        empty_alignment = BioToolkit.MultipleSequenceAlignment(BioToolkit.SeqRecordLite[])
        @test_throws ArgumentError BioToolkit.format_alignment(empty_alignment, "fasta")
        @test_throws ArgumentError BioToolkit.format_alignment(empty_alignment, "clustal")
        @test_throws ArgumentError BioToolkit.format_alignment(empty_alignment, "stockholm")

        mktempdir() do dir
            stockholm_path = joinpath(dir, "alignment.sto")
            BioToolkit.write_alignment(stockholm_path, alignment, "stockholm")
            written = read(stockholm_path, String)
            @test startswith(written, "# STOCKHOLM 1.0")
            @test endswith(written, "//")
        end

        rendered = sprint(io -> show(io, MIME("text/plain"), alignment))
        @test occursin("Alignment with 3 rows and 4 columns", rendered)
        @test occursin("clustal_consensus ****", rendered)
        @test occursin("seq1 ACGT", rendered)

        built = BioToolkit.progressive_multiple_sequence_alignment([
            BioToolkit.SeqRecordLite("ACGT"; identifier="seq1"),
            BioToolkit.SeqRecordLite("ACGT"; identifier="seq2"),
        ])
        @test length(built) == 2
        @test built[2].sequence == "ACGT"
        @test BioToolkit.sort!(copy(built); key = record -> record.identifier)[1].identifier == "seq1"

        mutated = copy(alignment)
        deleteat!(mutated, 2)
        @test length(mutated) == 2
        @test mutated[2].identifier == "seq3"

        @test BioToolkit.copy(alignment).column_annotations[:column_scores] == [1, 2, 3, 4]

        if Sys.which("clustalo") !== nothing || Sys.which("clustalw") !== nothing || Sys.which("clustalw2") !== nothing
            clustal_executable = Sys.which("clustalo")
            clustal_executable === nothing && (clustal_executable = Sys.which("clustalw2"))
            clustal_executable === nothing && (clustal_executable = Sys.which("clustalw"))
            clustal_external = BioToolkit.clustal_msa([
                BioToolkit.SeqRecordLite("ACGT"; identifier="seq1"),
                BioToolkit.SeqRecordLite("ACGT"; identifier="seq2"),
            ]; executable=clustal_executable, output_format="fasta")
            @test length(clustal_external) == 2
            @test clustal_external[1].sequence == "ACGT"
        end

        if Sys.which("muscle") !== nothing
            muscle_external = BioToolkit.muscle_msa([
                BioToolkit.SeqRecordLite("ACGT"; identifier="seq1"),
                BioToolkit.SeqRecordLite("ACGT"; identifier="seq2"),
            ]; executable=Sys.which("muscle"), output_format="fasta")
            @test length(muscle_external) == 2
            @test muscle_external[1].sequence == "ACGT"
        end
    end

    @testset "GenBank parsing" begin
        mktempdir() do dir
            input_path = joinpath(dir, "sample.gb")

            open(input_path, "w") do io
                write(io, "LOCUS       SCU49845       50 bp    DNA             PLN       21-JUN-1999\n")
                write(io, "DEFINITION  Saccharomyces cerevisiae TCP1-beta gene.\n")
                write(io, "ACCESSION   U49845\n")
                write(io, "VERSION     U49845.1\n")
                write(io, "KEYWORDS    synthetic; test record\n")
                write(io, "SOURCE      Saccharomyces cerevisiae (baker's yeast)\n")
                write(io, "            Eukaryota; Fungi; Ascomycota; Saccharomycetes.\n")
                write(io, "COMMENT     This record exercises metadata round-tripping.\n")
                write(io, "FEATURES             Location/Qualifiers\n")
                write(io, "     source          1..50\n")
                write(io, "                     /organism=\"Saccharomyces cerevisiae\"\n")
                write(io, "     gene            1..206\n")
                write(io, "                     /gene=\"TCP1\"\n")
                write(io, "ORIGIN\n")
                write(io, "        1 atgaccaatg cctgctgctg ctgctgctgc tgctgctgct gctgctgctg\n")
                write(io, "//\n")
            end

            records = BioToolkit.read_genbank(input_path)
            @test length(records) == 1
            @test records[1].locus == "SCU49845"
            @test records[1].accession == "U49845"
            @test records[1].version == "U49845.1"
            @test occursin("synthetic; test record", records[1].keywords)
            @test occursin("metadata round-tripping", records[1].comment)
            @test startswith(records[1].sequence, "ATGACCAATGCCTGCTGCTG")
            @test length(records[1].features) == 2
            @test records[1].features[2].key == "gene"
            @test records[1].features[2].qualifiers["gene"][1] == "TCP1"

            roundtrip_path = joinpath(dir, "roundtrip.gb")
            BioToolkit.write_genbank(roundtrip_path, records)
            written_lines = readlines(roundtrip_path)
            @test any(line -> startswith(line, "KEYWORDS    synthetic; test record."), written_lines)
            @test any(line -> startswith(line, "COMMENT     This record exercises metadata round-tripping."), written_lines)
            @test any(line -> startswith(line, "FEATURES             Location/Qualifiers"), written_lines)
            roundtrip = BioToolkit.read_genbank(roundtrip_path)
            @test roundtrip[1].keywords == records[1].keywords
            @test roundtrip[1].comment == records[1].comment
            @test roundtrip[1].locus_line == records[1].locus_line

            arrow_path = joinpath(dir, "sample.arrow")
            BioToolkit.ingest_genbank(input_path, arrow_path; chunk_size=1)

            table = Arrow.Table(arrow_path)
            @test table.locus == ["SCU49845"]
            @test table.locus_line[1] == strip("LOCUS       SCU49845       50 bp    DNA             PLN       21-JUN-1999")
            @test table.accession == ["U49845"]
            @test table.version == ["U49845.1"]
            @test occursin("synthetic; test record", table.keywords[1])
            @test occursin("metadata round-tripping", table.comment[1])
            @test table.feature_count == Int32[2]
            @test occursin("gene", table.feature_keys[1])
            @test occursin("1..206", table.feature_locations[1])
        end
    end

    @testset "Sequence annotation objects" begin
        simple = BioToolkit.parse_feature_location("12..20")
        @test simple isa BioToolkit.FeatureLocationLite
        @test BioToolkit.feature_bounds(simple) == (12, 20)
        @test BioToolkit.feature_length(simple) == 9
        @test BioToolkit.feature_spans(simple) == [(12, 20)]
        @test BioToolkit.feature_start(simple) == 12
        @test BioToolkit.feature_stop(simple) == 20
        @test BioToolkit.feature_strand(simple) == Int8(1)
        @test BioToolkit.feature_contains(simple, 15)
        @test !BioToolkit.feature_contains(simple, 21)

        compound = BioToolkit.parse_feature_location("join(1..4,10..12)")
        @test compound isa BioToolkit.CompoundFeatureLocation
        @test BioToolkit.feature_bounds(compound) == (1, 12)
        @test BioToolkit.feature_length(compound) == 7
        @test BioToolkit.feature_spans(compound) == [(1, 4), (10, 12)]

        complement = BioToolkit.parse_feature_location("complement(5..8)")
        @test complement isa BioToolkit.FeatureLocationLite
        @test complement.strand == -1

        record = BioToolkit.GenBankRecord(
            "TEST",
            "Synthetic annotated record",
            "ACC1",
            "ACC1.1",
            "Synthetic source",
            "Synthetic organism",
            "ACGTACGTACGTACGT",
            [
                BioToolkit.GenBankFeature("gene", "complement(5..8)", Dict("gene" => ["gene1"])),
                BioToolkit.GenBankFeature("CDS", "join(1..4,9..12)", Dict("gene" => ["gene2"])),
            ],
        )

        annotated = BioToolkit.annotate_genbank_record(record)
        @test annotated.identifier == "ACC1"
        @test annotated.annotations[:locus] == "TEST"
        @test length(annotated.features) == 2
        @test BioToolkit.feature_sequence(annotated, annotated.features[1]) == "ACGT"
        @test BioToolkit.feature_sequence(annotated, annotated.features[2]) == "ACGTACGT"
        @test BioToolkit.feature_extract(annotated, annotated.features[1]) == "ACGT"
        @test annotated.features[1].id == "gene1"
        @test annotated.features[2].id == "gene2"
        @test BioToolkit.feature_identifier(annotated.features[1]) == "gene1"
        @test BioToolkit.feature_annotation(annotated.features[1], :gene) == "gene1"
        @test BioToolkit.feature_summary(annotated.features[1]).identifier == "gene1"
        @test BioToolkit.select_features(annotated; feature_type="gene")[1].id == "gene1"
        @test BioToolkit.features_at(annotated, 6)[1].id == "gene1"
        @test BioToolkit.features_overlapping(annotated, 1, 4)[1].id == "gene2"
        @test length(BioToolkit.feature_table(annotated)) == 2

        fasta_like = BioToolkit.SeqRecordLite(
            "ACGTACGTACGT";
            identifier="plain1",
            name="plain1",
            description="plain record",
            annotations=Dict(:source => "synthetic"),
            letter_annotations=Dict(:quality => "IIIIIIIIIIII"),
        )
        sliced_fasta = BioToolkit.feature_slice(fasta_like, 3:8)
        @test sliced_fasta.sequence == "GTACGT"
        @test sliced_fasta.identifier == "plain1"
        @test sliced_fasta.letter_annotations[:quality] == "IIIIII"

        named_slice = BioToolkit.feature_slice(fasta_like, (start=4, stop=9))
        @test named_slice.sequence == "TACGTA"

        sliced_gb = BioToolkit.feature_slice(record, 2:10)
        @test sliced_gb isa BioToolkit.AnnotatedSeqRecord
        @test sliced_gb.sequence == "CGTACGTAC"
        @test length(sliced_gb.features) == 2

        sliced = annotated[2:10]
        @test sliced.sequence == "CGTACGTAC"
        @test BioToolkit.feature_spans(sliced.features[1]) == [(4, 7)]
        @test sliced.features[1].location.strand == -1
        @test BioToolkit.feature_spans(sliced.features[2]) == [(1, 3), (8, 9)]

        rc_record = BioToolkit.reverse_complement(annotated)
        @test rc_record.sequence == BioToolkit.reverse_complement(annotated.sequence)
        @test rc_record.features[1].location.strand == -1
        @test rc_record.features[2].location.strand == 1

        gff_feature = BioToolkit.SeqFeatureLite(BioToolkit.parse_gff_record("chr1\tsource\tgene\t10\t20\t.\t-\t.\tID=gene3;Name=Gene3"))
        @test gff_feature.feature_type == "gene"
        @test BioToolkit.feature_identifier(gff_feature) == "gene3"
        @test BioToolkit.feature_strand(gff_feature) == Int8(-1)
        @test BioToolkit.feature_annotation(gff_feature, :Name) == "Gene3"
        @test gff_feature.qualifiers["Name"] == ["Gene3"]
        @test BioToolkit.feature_contains(gff_feature, 15)
        @test BioToolkit.feature_overlaps(gff_feature, annotated.features[2])
    end

    @testset "Genomic interval operations" begin
        @test BioToolkit.normalize_interval((10, 5)) == (5, 10)
        @test BioToolkit.interval_length((5, 10)) == 5
        @test BioToolkit.interval_contains((5, 10), 5)
        @test !BioToolkit.interval_contains((5, 10), 10)
        @test BioToolkit.interval_overlaps((5, 10), (8, 12))
        @test !BioToolkit.interval_overlaps((5, 10), (10, 12))
        @test BioToolkit.interval_intersection((5, 10), (8, 12)) == (8, 10)
        @test BioToolkit.interval_union((5, 10), (8, 12)) == (5, 12)
        @test BioToolkit.merge_intervals([(10, 15), (1, 5), (5, 8)]) == [(1, 8), (10, 15)]
        @test BioToolkit.interval_difference((1, 10), (3, 6)) == [(1, 3), (6, 10)]
    end

    @testset "GFF parsing" begin
        record = BioToolkit.parse_gff_record("chr1\tsource\tgene\t10\t20\t.\t+\t.\tID=gene1;Name=Gene1")
        @test record.chrom == "chr1"
        @test record.source == "source"
        @test record.feature == "gene"
        @test record.start == 10
        @test record.stop == 20
        @test record.score === missing
        @test record.strand == "+"
        @test record.phase === missing
        @test record.attributes == "ID=gene1;Name=Gene1"
        @test record.attribute_map["ID"] == ["gene1"]

        mktempdir() do dir
            input_path = joinpath(dir, "sample.gff")
            output_path = joinpath(dir, "written.gff")

            open(input_path, "w") do io
                write(io, "##gff-version 3\n")
                write(io, "chr1\tsource\tgene\t10\t20\t.\t+\t.\tID=gene1;Name=Gene1\n")
            end

            parsed = BioToolkit.read_gff(input_path)
            @test parsed == [record]
            @test parsed[1].attribute_map["Name"] == ["Gene1"]

            BioToolkit.write_gff(output_path, parsed)
            @test BioToolkit.read_gff(output_path) == parsed
        end

        mktempdir() do dir
            input_path = joinpath(dir, "sample.gff")
            output_path = joinpath(dir, "sample.arrow")

            open(input_path, "w") do io
                write(io, "##gff-version 3\n")
                write(io, "chr1\tsource\tgene\t10\t20\t.\t+\t.\tID=gene1\n")
                write(io, "chr1\tsource\tmRNA\t10\t20\t5.5\t+\t0\tID=tx1;Parent=gene1\n")
            end

            BioToolkit.ingest_gff(input_path, output_path; chunk_size=1)

            table = Arrow.Table(output_path)
            @test table.chrom == ["chr1", "chr1"]
            @test table.feature == ["gene", "mRNA"]
            @test table.start == Int32[10, 10]
            @test table.score[1] === missing
            @test table.score[2] == Float32(5.5)
            @test table.phase[1] === missing
            @test table.phase[2] == Int8(0)
        end
    end

    @testset "Motif utilities" begin
        records = [
            BioToolkit.SeqRecordLite("ACGT"; identifier="s1"),
            BioToolkit.SeqRecordLite("ACGA"; identifier="s2"),
            BioToolkit.SeqRecordLite("ACGT"; identifier="s3"),
        ]

        counts = BioToolkit.motif_counts(records)
        @test counts.counts[:, 1] == [3, 0, 0, 0]
        @test BioToolkit.motif_consensus(counts) == "ACGT"
        @test BioToolkit.motif_consensus(counts; threshold=0.6) == "ACGT"
        @test BioToolkit.motif_consensus(counts; threshold=0.9) == "ACGN"

        @test_throws ArgumentError BioToolkit.motif_counts(["ACGT", "ACG"])
        @test_throws ArgumentError BioToolkit.motif_counts(["AXGT"]) 

        pwm = BioToolkit.motif_pwm(records; pseudocount=0.5)
        @test size(pwm.values) == (4, 4)
        @test pwm.values[1, 1] > 0

        pwm_background = BioToolkit.motif_pwm(records; pseudocount=0.5, background=Dict('A' => 0.4, 'C' => 0.2, 'G' => 0.2, 'T' => 0.2))
        @test size(pwm_background.values) == (4, 4)
        @test pwm_background.values[1, 1] > pwm.values[1, 1] - 0.1
        @test_throws ArgumentError BioToolkit.motif_pwm(counts; pseudocount=-1)
        @test BioToolkit.motif_information_content(pwm) > 0.0

        frequencies = BioToolkit.motif_frequency_matrix(counts; pseudocount=0.5)
        @test size(frequencies.values) == (4, 4)
        @test isapprox(sum(frequencies.values[:, 1]), 1.0; atol=1e-12)
        @test BioToolkit.motif_entropy(frequencies) >= 0.0
        @test BioToolkit.motif_relative_entropy(frequencies) > 0.0

        discovery_sequences = ["TTTACGTTTT", "GGGACGTGGG", "CCCACGTCCC"]
        discoveries = BioToolkit.discover_motifs(discovery_sequences; k=4, top_n=1, min_support=2, max_mismatches=0)
        @test length(discoveries) == 1
        @test discoveries[1].seed == "ACGT"
        @test discoveries[1].support == 3
        @test discoveries[1].information_content > 0.0
        @test discoveries[1].counts.counts[:, 1] == [3, 0, 0, 0]

        scan_hits = BioToolkit.motif_scan(discovery_sequences[1], discoveries[1].pwm; threshold=0.0)
        @test any(hit -> hit.start == 4, scan_hits)
        @test BioToolkit.motif_scan(discovery_sequences[1], discoveries[1].pwm; threshold=100.0) == BioToolkit.MotifHit[]

        both_strand_hits = BioToolkit.motif_scan_both_strands("TTTACGTTTTGGGACGTGGG", discoveries[1].pwm; threshold=0.0)
        @test any(hit -> hit.strand == Int8(1), both_strand_hits)
    end

    @testset "Motif parsers and logos" begin
        meme_data = """MEME version 5.5.0
ALPHABET= ACGT

MOTIF MotifOne
letter-probability matrix: alength= 4 w= 4 nsites= 3 E= 1e-5
0.80 0.10 0.05 0.05
0.10 0.70 0.10 0.10
0.05 0.10 0.80 0.05
0.25 0.25 0.25 0.25
Motif 1 sites sorted by position p-value
seq1  2  1e-4  ACGT
seq2  8  2e-4  ACGT
"""
        meme_path, meme_io = mktemp()
        write(meme_io, meme_data)
        close(meme_io)

        alignace_data = """ALIGNACE 3.0
MOTIF AlignOne
seq1  1  ACGT
seq2  4  ACGT
seq3  7  ACGT
"""
        alignace_path, alignace_io = mktemp()
        write(alignace_io, alignace_data)
        close(alignace_io)

        try
            meme_profiles = BioToolkit.read_meme(meme_path)
            @test length(meme_profiles) == 1
            @test meme_profiles[1].name == "MotifOne"
            @test length(meme_profiles[1].occurrences) == 2
            @test BioToolkit.motif_consensus(meme_profiles[1]) == "ACGN"
            @test occursin("<svg", BioToolkit.sequence_logo_svg(meme_profiles[1]; width=360, height=160, title="MotifOne"))

            alignace_profiles = BioToolkit.read_alignace(alignace_path)
            @test length(alignace_profiles) == 1
            @test alignace_profiles[1].name == "AlignOne"
            @test length(alignace_profiles[1].occurrences) == 3
            @test alignace_profiles[1].counts.counts[:, 1] == [3, 0, 0, 0]
            @test occursin("<svg", BioToolkit.sequence_logo_svg(alignace_profiles[1].pwm; width=360, height=160))

            jaspar_data = ">MA0001.1 Example motif\nA [ 3 1 0 0 ]\nC [ 0 2 4 0 ]\nG [ 1 1 0 5 ]\nT [ 0 0 0 0 ]\n"
            jaspar_path, jaspar_io = mktemp()
            write(jaspar_io, jaspar_data)
            close(jaspar_io)
            try
                jaspar_profiles = BioToolkit.read_jaspar(jaspar_path)
                @test length(jaspar_profiles) == 1
                @test jaspar_profiles[1].metadata["source"] == "JASPAR"
                @test jaspar_profiles[1].metadata["jaspar_id"] == "MA0001.1"
                @test jaspar_profiles[1].name == "Example motif"
                @test size(jaspar_profiles[1].counts.counts) == (4, 4)
            finally
                rm(jaspar_path, force=true)
            end
        finally
            rm(meme_path, force=true)
            rm(alignace_path, force=true)
        end
    end

    @testset "Chunked ingestion" begin
        mktempdir() do dir
            input_path = joinpath(dir, "sample.vcf")
            output_path = joinpath(dir, "sample.arrow")

            open(input_path, "w") do io
                write(io, "##fileformat=VCFv4.2\n")
                write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")
                write(io, "chr1\t10\trs1\tA\tG\t99.5\n")
                write(io, "chr1\t20\trs2\tC\tT\t80.0\n")
                write(io, "chr2\t30\trs3\tG\tA\t70.0\n")
            end

            BioToolkit.ingest_vcf(input_path, output_path; chunk_size=2)

            stream = Arrow.Stream(output_path)
            @test length(collect(Tables.partitions(stream))) == 2

            table = Arrow.Table(output_path)
            @test table.pos == Int32[10, 20, 30]
        end
    end

    @testset "BED ingestion" begin
        mktempdir() do dir
            input_path = joinpath(dir, "sample.bed")
            output_path = joinpath(dir, "sample-bed.arrow")
            written_path = joinpath(dir, "written.bed")

            open(input_path, "w") do io
                write(io, "chr1\t0\t10\n")
                write(io, "chr1\t10\t20\n")
                write(io, "chr2\t20\t30\n")
            end

            parsed = BioToolkit.read_bed(input_path)
            @test parsed[1].chrom == "chr1"
            @test parsed[1].start == Int32(0)
            @test parsed[1].stop == Int32(10)

            BioToolkit.write_bed(written_path, parsed)
            @test BioToolkit.read_bed(written_path) == parsed

            BioToolkit.ingest_bed(input_path, output_path; chunk_size=2)

            table = Arrow.Table(output_path)
            @test table.start == Int32[0, 10, 20]
            @test table.stop == Int32[10, 20, 30]
            @test table.chrom == ["chr1", "chr1", "chr2"]

            coverage = BioToolkit.window_coverage(table, "chr1", 10)
            @test coverage == Dict(0 => 10, 1 => 10)
        end
    end

    @testset "GC skew and dotmatrix" begin
        @test BioToolkit.gc_skew("GCGCG") == [1, 0, 1, 0, 1]
        @test BioToolkit.gc_skew("CCCC") == [-1, -2, -3, -4]
        @test BioToolkit.gc_skew("GGGG") == [1, 2, 3, 4]
        @test BioToolkit.gc_skew("AAAA") == [0, 0, 0, 0]

        min_pos = BioToolkit.minimum_skew("CCCC")
        @test min_pos == [4]

        min_pos2 = BioToolkit.minimum_skew("CCGGCC")
        @test min_pos2 == [2, 6]  # cumulative: -1,-2,-1,0,-1,-2 → min=-2 at 2,6

        dm = BioToolkit.dotmatrix("ACGT", "ACGT")
        @test size(dm) == (4, 4)
        @test dm[1, 1] == Int8(1)  # A-A match
        @test dm[1, 2] == Int8(0)  # A-C mismatch
        @test dm[2, 2] == Int8(1)  # C-C match
        @test sum(dm) == 4         # identity → 4 matches on diagonal

        dm2 = BioToolkit.dotmatrix("AA", "AAA")
        @test size(dm2) == (2, 3)
        @test all(dm2 .== Int8(1))  # all A's match
    end

    @testset "Protein statistics" begin
        # Small known protein: "ACDEFGHIKLMNPQRSTVWY" (all 20 standard AAs)
        all_aa = "ACDEFGHIKLMNPQRSTVWY"

        # Protein mass (monoisotopic) — sum of all 20 residue masses + water
        mass_mono = BioToolkit.protein_mass(all_aa; type="monoisotopic")
        @test mass_mono > 2000.0  # all 20 AAs should be > 2 kDa
        @test mass_mono < 3000.0

        mass_avg = BioToolkit.protein_mass(all_aa; type="average")
        @test mass_avg > mass_mono  # average masses are slightly larger

        # Extinction coefficient: Y=1 → 1490, W=1 → 5500, C=1 → 125
        ec = BioToolkit.extinction_coefficient(all_aa)
        @test ec == 1490 + 5500 + 125  # 7115

        # Single residue tests
        @test BioToolkit.extinction_coefficient("YYY") == 3 * 1490
        @test BioToolkit.extinction_coefficient("WWW") == 3 * 5500
        @test BioToolkit.extinction_coefficient("CCC") == 3 * 125

        # Instability index — needs at least 2 residues
        ii = BioToolkit.instability_index(all_aa)
        @test ii isa Float64
        @test_throws ArgumentError BioToolkit.instability_index("A")

        # GRAVY
        gravy_val = BioToolkit.gravy(all_aa)
        @test gravy_val isa Float64
        @test BioToolkit.gravy("III") ≈ 4.5  # I has hydropathicity 4.5

        # Aliphatic index: for "AVILAAAA" → high due to many A's
        ai = BioToolkit.aliphatic_index("AVIL")
        @test ai > 0.0
        @test ai isa Float64

        # Isoelectric point — should be between 0 and 14
        pi_val = BioToolkit.isoelectric_point(all_aa)
        @test 0.0 < pi_val < 14.0

        # Known pI: pure Lysine (K) → pI ≈ 10.5
        pi_k = BioToolkit.isoelectric_point("KKKKK")
        @test 9.0 < pi_k < 12.0

        # Known pI: pure Glutamic acid (E) → pI ≈ 3.2
        pi_e = BioToolkit.isoelectric_point("EEEEE")
        @test 2.0 < pi_e < 5.0

        # protparam returns a NamedTuple with all fields
        pp = BioToolkit.protparam(all_aa)
        @test pp.length == 20
        @test pp.molecular_weight_mono ≈ mass_mono
        @test pp.molecular_weight_avg ≈ mass_avg
        @test pp.extinction_coefficient == ec
        @test pp.gravy ≈ gravy_val
        @test 0.0 < pp.isoelectric_point < 14.0
    end

    @testset "Search heuristics (BLAST-like)" begin
        targets = ["ACGTGCCGATCGATCGATCGA", "ACGTGCGTAGCTAGCTG", "ATGGCCATTGTAATGGGCCGCTGAA"]
        index = BioToolkit.build_index(targets, ["seq1", "seq2", "seq3"]; k=4)
        @test index.k == 4
        
        query = "GCCGATCG"
        matrix = BioToolkit.substitution_matrix("ACGT"; match=2, mismatch=-1)
        scoring = BioToolkit.MatrixPairwiseScoring(matrix)
        
        hsps = BioToolkit.local_search(query, index; scoring=scoring, x_drop=5, min_score=10)
        
        @test !isempty(hsps)
        best_hsp = hsps[1]
        @test best_hsp.target_idx == 1
        @test best_hsp.score == 16 
        @test best_hsp.query_start == 1
        @test best_hsp.query_end == 8
    end

end
