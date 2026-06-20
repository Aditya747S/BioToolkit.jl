using Test
using BioToolkit

@testset "Parser edge cases" begin
    @test BioToolkit.validate_dna("") == true
    @test BioToolkit.validate_dna("ATGNRYSWKMBDHV") == true
    @test BioToolkit.validate_dna(BioToolkit.DNASeq("atgNRYSWKMBDHV")) == true
    @test BioToolkit.reverse_complement("R") == "Y"
    @test BioToolkit.reverse_complement("Y") == "R"
    @test BioToolkit.gc_content("ATGNRR") == 1 / 3
    @test BioToolkit.gc_content(BioToolkit.DNASeq("ATGNRR")) == 1 / 3
    @test_throws ArgumentError BioToolkit.DNASeq("ATBX")
    @test BioToolkit.RNASeq("aucg") == BioToolkit.RNASeq("AUCG")

    mktempdir() do dir
        fastq_path = joinpath(dir, "bad.fastq")
        open(fastq_path, "w") do io
            write(io, "@seq1\n")
            write(io, "ACGT\n")
            write(io, "+\n")
            write(io, "III\n")
        end

        @test_throws ArgumentError BioToolkit.read_fastq(fastq_path)
        @test_throws ArgumentError BioToolkit.read_fastq(fastq_path; alphabet=BioToolkit.DNAAlphabet)
    end

    mktempdir() do dir
        fastq_path = joinpath(dir, "lower.fastq")
        open(fastq_path, "w") do io
            write(io, "@seq1\n")
            write(io, "acgt\n")
            write(io, "+\n")
            write(io, "IIII\n")
        end

        typed_records = BioToolkit.read_fastq(fastq_path; alphabet=BioToolkit.DNAAlphabet)
        @test typed_records[1].sequence isa BioToolkit.DNASeq
        @test typed_records[1].sequence == "ACGT"
    end

    mktempdir() do dir
        fasta_path = joinpath(dir, "bad.fasta")
        open(fasta_path, "w") do io
            write(io, "ACGT\n")
        end

        @test_throws ArgumentError BioToolkit.read_fasta(fasta_path)
    end

    @test_throws ArgumentError BioToolkit.translate_dna("ATGAA")
    @test_throws ArgumentError BioToolkit.translate_dna(BioToolkit.DNASeq("ATGAA"))

    @test_throws ArgumentError BioToolkit.substitution_matrix("ACGT", [1 2; 3 4])

    mktempdir() do dir
        fasta_path = joinpath(dir, "empty.fasta")
        write(fasta_path, "")
        records = BioToolkit.read_fasta(fasta_path)
        @test isempty(records)
        @test eltype(records) <: BioToolkit.SeqRecord{BioToolkit.DNAAlphabet}
    end

    seq_info = BioToolkit.GenomicRanges.SeqInfo("chrM", 16569, true)
    buffer = IOBuffer()
    show(buffer, seq_info)
    @test occursin("circular=true", String(take!(buffer)))
end

@testset "Streaming FASTQ" begin
    mktempdir() do dir
        fastq_path = joinpath(dir, "stream.fastq")
        open(fastq_path, "w") do io
            write(io, "@seq1\nACGT\n+\nIIII\n")
            write(io, "@seq2\nTTTT\n+\n####\n")
        end

        function count_records(path)
            count = 0
            for record in BioToolkit.each_fastq_record(path; sequence_type=BioToolkit.DNASeq)
                count += 1
                @test record isa BioToolkit.FastqRecord{BioToolkit.DNAAlphabet}
            end
            return count
        end

        @test @inferred(count_records(fastq_path)) == 2
        @test length(BioToolkit.read_fastq(fastq_path; alphabet=BioToolkit.DNAAlphabet)) == 2
    end
end