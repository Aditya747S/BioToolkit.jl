using Test
using Printf
using LinearAlgebra
using SparseArrays
using DataFrames
using BioToolkit

@testset "Dirty Real-World Stress Suite" begin
    @testset "I/O & Parsing" begin
        # Dirty Data Pattern: CRLF FASTA, wrapped sequence lines, whitespace in headers,
        # embedded spaces/dashes in sequence text, and explicit empty records.
        # Expected Behavior: FASTA should parse permissively, preserve header text,
        # uppercase the sequence payload, and keep empty records as zero-length sequences.
        @testset "FASTA" begin
            fasta_text = ">read one with spaces\r\nac gt-\r\nnn\r\n>empty record\r\n>single line header\r\nACgtry\r\n"
            records = BioToolkit.read_fasta(IOBuffer(fasta_text))

            @test length(records) == 3
            @test records[1].identifier == "read one with spaces"
            @test String(records[1].sequence) == "AC GT-NN"
            @test records[2].identifier == "empty record"
            @test length(records[2].sequence) == 0
            @test records[3].identifier == "single line header"
            @test String(records[3].sequence) == "ACGTRY"
        end

        # Dirty Data Pattern: FASTQ identifiers with spaces, lower-case sequence text,
        # and truncated quality strings that are shorter than the sequence.
        # Expected Behavior: headers should be split into identifier/description, valid
        # records should uppercase the sequence, and length mismatches must throw cleanly.
        @testset "FASTQ" begin
            fastq_text = "@read with spaces\r\nacgtn\r\n+\r\n!!!!!\r\n"
            records = BioToolkit.read_fastq(IOBuffer(fastq_text))

            @test length(records) == 1
            @test records[1].identifier == "read"
            @test records[1].description == "read with spaces"
            @test String(records[1].sequence) == "ACGTN"
            @test records[1].quality == "!!!!!"

            empty_fastq = "@empty\n\n+\n\n"
            empty_records = BioToolkit.read_fastq(IOBuffer(empty_fastq))
            @test length(empty_records) == 1
            @test length(empty_records[1].sequence) == 0
            @test empty_records[1].quality == ""

            short_quality = "@short\nACGT\n+\n!!!\n"
            @test_throws ArgumentError BioToolkit.read_fastq(IOBuffer(short_quality))
        end

        # Dirty Data Pattern: a truncated BAM file boundary plus synthetic in-memory
        # records with unmapped positions and CIGAR strings that over-consume reference.
        # Expected Behavior: corrupt BAM input should fail loudly; direct record handling
        # should not segfault and overlap logic should stay defined on pathological records.
        @testset "BAM" begin
            mktempdir() do dir
                bam_path = joinpath(dir, "truncated.bam")
                open(bam_path, "w") do io
                    write(io, "BAM")
                end
                @test_throws Exception BioToolkit.read_bam(bam_path)
            end

            unmapped = BioToolkit.BamRecord(
                "unmapped_read",
                "chr1",
                -1,
                [BioToolkit.BamCigarOp(10, 'M')],
                BioToolkit.DNASeq("ACGTACGTAA");
                flag=4,
                mapq=0)
            overlong = BioToolkit.BamRecord(
                "overlong_cigar",
                "chr1",
                100,
                [BioToolkit.BamCigarOp(100, 'M')],
                BioToolkit.DNASeq("ACGT");
                flag=0,
                mapq=60)

            region = BioToolkit.GenomicInterval("chr1", 150, 160)
            @test BioToolkit._bam_overlaps(unmapped, region) == false
            @test BioToolkit._bam_overlaps(overlong, region) == true
        end

        # Dirty Data Pattern: missing INFO/FORMAT fields, missing sample columns,
        # genotype placeholders like ./., and path-based VCF decoding.
        # Expected Behavior: text parsers should default missing fields cleanly and
        # the GWAS reader should convert missing genotypes to NaN rather than crashing.
        @testset "VCF" begin
            minimal_line = "chr1\t42\trs1\tA\tG\t99.5"
            parsed = BioToolkit.parse_vcf_record(minimal_line)
            @test parsed.chrom == "chr1"
            @test parsed.pos == 42
            @test parsed.filter == "PASS"
            @test parsed.info == "."
            @test parsed.format == ""
            @test isempty(parsed.samples)

            document_text = "##fileformat=VCFv4.2\n" * minimal_line * "\n"
            document = BioToolkit.read_vcf_document(IOBuffer(document_text))
            @test length(document.records) == 1
            @test document.header.sample_names == String[]
            @test document.records[1].format == ""
            @test document.records[1].info == "."

            mktempdir() do dir
                vcf_path = joinpath(dir, "dirty.vcf")
                open(vcf_path, "w") do io
                    write(io, "##fileformat=VCFv4.2\n")
                    write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n")
                    write(io, "chr1\t10\trs_missing\tA\tG\t.\tPASS\t.\tGT\t0/1\t./.\t1/1\n")
                    write(io, "chr1\t20\trs_mono\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\n")
                    write(io, "chr1\t30\trs_signal\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n")
                end

                genotypes = BioToolkit.GWAS.read_vcf(vcf_path)
                @test size(genotypes) == (3, 3)
                @test isfinite(genotypes[1, 1])
                @test isnan(genotypes[2, 1])

                variant_missingness = BioToolkit.GWAS.calculate_missingness(genotypes; by=:variant)
                sample_missingness = BioToolkit.GWAS.calculate_missingness(genotypes; by=:sample)
                maf_mask, maf = BioToolkit.GWAS.calculate_maf(genotypes; return_frequency=true)

                @test variant_missingness ≈ [1 / 3, 0.0, 0.0]
                @test sample_missingness ≈ [0.0, 1 / 3, 0.0]
                @test maf_mask == Bool[true, false, true]
                @test all(isfinite, maf[maf_mask])
            end
        end

        # Dirty Data Pattern: BED 0-based intervals, GFF 1-based intervals, malformed
        # attributes, large chromosome names, and inverted coordinates.
        # Expected Behavior: BED/GFF records should parse valid records, reject inverted
        # coordinates, and normalize into closed 1-based GenomicInterval coordinates.
        @testset "BED/GFF" begin
            bed_records = BioToolkit.read_bed(IOBuffer("chr1_gl000220_random\t0\t10\n"))
            @test length(bed_records) == 1
            @test bed_records[1].chrom == "chr1_gl000220_random"
            bed_interval = BioToolkit.GenomicInterval(bed_records[1])
            @test bed_interval.left == 1
            @test bed_interval.right == 10

            gff_text = "chr1_gl000220_random\tsource\tgene\t1\t10\t.\t+\t.\tID=gene 1;Name=bad quote;Note=unescaped spaces\n"
            gff_records = BioToolkit.read_gff(IOBuffer(gff_text))
            @test length(gff_records) == 1
            @test gff_records[1].chrom == "chr1_gl000220_random"
            @test gff_records[1].attribute_map["ID"] == ["gene 1"]
            @test gff_records[1].attribute_map["Note"] == ["unescaped spaces"]

            gff_interval = BioToolkit.GenomicInterval(gff_records[1])
            @test gff_interval.left == 1
            @test gff_interval.right == 10

            @test_throws ArgumentError BioToolkit.read_bed(IOBuffer("chr1\t10\t1\n"))
            @test_throws ArgumentError BioToolkit.read_bed(IOBuffer("chr1\t-1\t10\n"))
            @test_throws ArgumentError BioToolkit.read_gff(IOBuffer("chr1\tsource\tgene\t10\t1\t.\t+\t.\tID=bad\n"))
        end
    end

    @testset "Sequence Typing & Ambiguity" begin
        # Dirty Data Pattern: mixed case symbols, IUPAC ambiguity codes, all-gap strings,
        # zero-length sequences, and invalid whitespace in typed sequence constructors.
        # Expected Behavior: valid ambiguity codes should uppercase and survive typing;
        # invalid whitespace should fail at construction; empty and all-gap sequences
        # should be representable without special casing.
        @testset "DNA / AA typing" begin
            dna = BioToolkit.DNASeq("aCgTnryw--")
            aa = BioToolkit.AASeq("xbz*uo")

            @test String(dna) == "ACGTNRYW--"
            @test BioToolkit.validate_dna(dna) == true
            @test length(BioToolkit.DNASeq("")) == 0
            @test String(BioToolkit.DNASeq("----")) == "----"
            @test String(aa) == "XBZ*UO"
            @test_throws ArgumentError BioToolkit.DNASeq("AC GT")
        end

        # Dirty Data Pattern: cross-alphabet coercion between typed DNA and protein
        # sequences in pairwise alignment.
        # Expected Behavior: method dispatch should reject incompatible alphabets cleanly
        # instead of trying to coerce them into a shared runtime representation.
        @testset "Type coercion" begin
            @test_throws MethodError BioToolkit.pairwise_align(BioToolkit.DNASeq("ACGT"), BioToolkit.AASeq("ACGT"))
        end
    end

    @testset "Genomic Intervals & Arithmetic" begin
        # Dirty Data Pattern: intervals that touch at exactly one base, zero-length spans,
        # flipped start/stop coordinates, and high-cardinality interval collections.
        # Expected Behavior: closed intervals should count a shared boundary as overlap,
        # zero-length intervals should remain valid, and large collections should still
        # resolve overlaps deterministically.
        @testset "Boundary arithmetic" begin
            left = BioToolkit.GenomicInterval("chr1", 100, 200)
            right = BioToolkit.GenomicInterval("chr1", 200, 300)
            overlap_query = BioToolkit.GenomicInterval("chr1", 200, 200)

            left_collection = BioToolkit.build_collection([left])
            right_collection = BioToolkit.build_collection([right])
            merged = BioToolkit.intersect(left_collection, right_collection)
            reduced = BioToolkit.setdiff(left, right_collection)

            @test length(BioToolkit.find_overlaps(overlap_query, BioToolkit.build_collection([left, right]))) == 2
            @test length(merged) == 1
            @test merged[1].left == 200
            @test merged[1].right == 200
            @test length(reduced) == 1
            @test reduced[1].left == 100
            @test reduced[1].right == 199

            zero_length = BioToolkit.GenomicInterval("chr1", 50, 50)
            @test length(BioToolkit.find_overlaps(zero_length, BioToolkit.build_collection([zero_length]))) == 1

            flipped = BioToolkit.GenomicInterval("chr1", 300, 100)
            @test flipped.left == 100
            @test flipped.right == 300
        end

        # Dirty Data Pattern: a very large set of tiny intervals against a single wide
        # query interval, scaled down from a million-interval production load.
        # Expected Behavior: the overlap index should stay stable and return every hit
        # without blowing up on a quadratic scan.
        @testset "Scale stress" begin
            interval_count = 25_000
            tiny_intervals = [BioToolkit.GenomicInterval("chrStress", index, index) for index in 1:interval_count]
            collection = BioToolkit.build_collection(tiny_intervals)
            hits = BioToolkit.find_overlaps(BioToolkit.GenomicInterval("chrStress", 1, interval_count), collection)

            @test length(hits) == interval_count
        end
    end

    @testset "Single-Cell & Omics Data" begin
        # Dirty Data Pattern: extreme sparsity, dead genes/cells, and mismatched row
        # metadata lengths.
        # Expected Behavior: sparse matrices should normalize without NaNs/Infs; shape
        # mismatches should fail loudly at construction time.
        @testset "Sparse matrices" begin
            counts = spzeros(Int, 12, 8)
            counts[1, 1] = 5
            counts[4, 4] = 7
            counts[10, 8] = 1
            counts[12, 2] = 4

            gene_ids = ["gene_$(index)" for index in 1:12]
            cell_ids = ["cell_$(index)" for index in 1:8]
            sce = BioToolkit.SingleCellExperiment(counts, gene_ids, cell_ids)
            normalized = BioToolkit.normalize_counts(sce)

            @test size(normalized) == size(counts)
            @test all(isfinite, Array(normalized))
            @test BioToolkit.count_matrix(sce).sample_ids == cell_ids

            @test_throws ArgumentError BioToolkit.SingleCellExperiment(counts, gene_ids[1:end-1], cell_ids)

            se = BioToolkit.SummarizedExperiment(
                Matrix{Float64}(counts);
                rowData=Dict{Symbol,Vector}(:gene_id => Vector{String}(gene_ids)),
                colData=Dict{Symbol,Vector}(:cell_id => Vector{String}(cell_ids)))
            @test size(BioToolkit.assay(se)) == size(counts)
            @test_throws DimensionMismatch BioToolkit.SummarizedExperiment(
                Matrix{Float64}(counts);
                rowData=Dict{Symbol,Vector}(:gene_id => Any[gene_ids[1:end-1]...]),
                colData=Dict{Symbol,Vector}(:cell_id => Any[cell_ids...]))
        end

        # Dirty Data Pattern: omics workflows with singular statistical models are
        # handled in the GWAS linear-model core; here we keep the single-cell layer on
        # the sparse/count/metadata edge cases that are most likely to appear in practice.
        # Expected Behavior: the single-cell constructor and normalization path should
        # stay robust even when the matrix is nearly empty.
        @testset "Near-empty normalization" begin
            counts = spzeros(Int, 5, 4)
            counts[2, 3] = 11
            sce = BioToolkit.SingleCellExperiment(counts, ["g$(i)" for i in 1:5], ["c$(i)" for i in 1:4])
            normalized = BioToolkit.normalize_counts(sce; scale_factor=1e4, log_transform=true)

            @test all(isfinite, Array(normalized))
            @test any(x -> x > 0, Array(normalized))
        end
    end

    @testset "Alignment & Search" begin
        # Dirty Data Pattern: identical sequences, completely unrelated sequences,
        # and a 10 bp query against a 50,000 bp target to stress numerical stability.
        # Expected Behavior: exact matches should score exactly, anti-matches should
        # stay negative, and long/short asymmetry should not destabilize the scorer.
        @testset "Pairwise alignment" begin
            ident = BioToolkit.DNASeq("ACGTAC")
            ident_res = BioToolkit.pairwise_align(ident, ident; match=2, mismatch=-1, gap=-2)
            @test ident_res.score == 12
            @test ident_res.matches == length(ident)
            @test ident_res.identity == 1.0

            mismatch_res = BioToolkit.pairwise_align(BioToolkit.DNASeq("AAAA"), BioToolkit.DNASeq("TTTT"); match=1, mismatch=-1, gap=-2)
            @test mismatch_res.score < 0
            @test mismatch_res.matches == 0

            long_target = BioToolkit.DNASeq(repeat("T", 50_000))
            short_query = BioToolkit.DNASeq(repeat("T", 10))
            long_res = BioToolkit.pairwise_align(short_query, long_target; is_local=true, match=2, mismatch=-1, gap=-2)
            @test long_res.score >= 20
        end

        # Dirty Data Pattern: k-mer search queries containing ambiguity codes that are
        # absent from the target index.
        # Expected Behavior: the index should simply return no hits instead of trying to
        # reinterpret ambiguous bytes as valid seeds.
        @testset "K-mer search" begin
            index = BioToolkit.build_index([BioToolkit.DNASeq("ACGTACGTAC")]; k=4)
            hits = BioToolkit.local_search(
                BioToolkit.DNASeq("ACNN"),
                index;
                scoring=BioToolkit.LinearPairwiseScoring(2, -1),
                x_drop=5,
                min_score=4)

            @test isempty(hits)
        end
    end

    @testset "GWAS & Population Genetics" begin
        # Dirty Data Pattern: NaN missing genotypes, monomorphic variants, tiny sample
        # sizes, and a singular/collinear covariate matrix.
        # Expected Behavior: QC helpers should ignore NaN safely, monomorphic SNPs should
        # be filtered or returned as null signal, and the linear model should stay finite.
        @testset "GWAS QC and scan" begin
            mktempdir() do dir
                vcf_path = joinpath(dir, "dirty_gwas.vcf")
                open(vcf_path, "w") do io
                    write(io, "##fileformat=VCFv4.2\n")
                    write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n")
                    write(io, "chr1\t10\trs_missing\tA\tG\t.\tPASS\t.\tGT\t0/1\t./.\t1/1\n")
                    write(io, "chr1\t20\trs_mono\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\n")
                    write(io, "chr1\t30\trs_signal\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n")
                end

                genotypes = BioToolkit.GWAS.read_vcf(vcf_path)
                @test size(genotypes) == (3, 3)

                variant_missingness = BioToolkit.GWAS.calculate_missingness(genotypes; by=:variant)
                sample_missingness = BioToolkit.GWAS.calculate_missingness(genotypes; by=:sample)
                maf_mask, maf = BioToolkit.GWAS.calculate_maf(genotypes; return_frequency=true)

                @test variant_missingness ≈ [1 / 3, 0.0, 0.0]
                @test sample_missingness ≈ [0.0, 1 / 3, 0.0]
                @test maf_mask == Bool[true, false, true]
                @test all(isfinite, maf[maf_mask])

                phenotype = [0.0, 1.0, 0.0]
                underdetermined_covariates = hcat(
                    ones(3),
                    [0.0, 1.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [1.0, 1.0, 1.0],
                )
                underdetermined_failed = false
                try
                    BioToolkit.GWAS.gwas_linear_scan(genotypes, phenotype; covariates=underdetermined_covariates)
                catch err
                    underdetermined_failed = err isa DimensionMismatch
                end
                @test underdetermined_failed

                square_singular_covariates = hcat(
                    ones(3),
                    [0.0, 1.0, 0.0],
                    [0.0, 1.0, 0.0],
                )
                singular_failed = false
                try
                    BioToolkit.GWAS.gwas_linear_scan(genotypes, phenotype; covariates=square_singular_covariates)
                catch err
                    singular_failed = err isa LinearAlgebra.SingularException
                end
                @test singular_failed

                valid_covariates = hcat(ones(3), [0.0, 1.0, 0.0])
                scan = BioToolkit.GWAS.gwas_linear_scan(genotypes, phenotype; covariates=valid_covariates)

                @test scan.sample_size == 3
                @test length(scan.pvalue) == 3
                @test all(isfinite, scan.pvalue)
            end

            @test BioToolkit.GWAS.hwe_exact(4, 0, 0) == 1.0
            @test BioToolkit.GWAS.calculate_maf([0.0 0.0; 0.0 0.0]; return_frequency=true)[1] == Bool[false, false]
        end

        # Dirty Data Pattern: monomorphic loci, tiny populations, and degenerate allele
        # counts.
        # Expected Behavior: population-genetic summaries should return finite values and
        # treat the locus as non-informative rather than blowing up on divide-by-zero.
        @testset "Population genetics" begin
            population = BioToolkit.Population("tiny", [
                BioToolkit.PopGenIndividual("ind1", [BioToolkit.Locus(("A", "A"))]),
            ])

            @test BioToolkit.allele_frequencies(population, 1) == Dict("A" => 1.0)
            @test BioToolkit.genotype_frequencies(population, 1) == Dict(Set(["A"]) => 1.0)
            @test BioToolkit.heterozygosity_observed(population, 1) == 0.0
            @test BioToolkit.heterozygosity_expected(population, 1) == 0.0
            @test BioToolkit.hardy_weinberg_test(population, 1) == 1.0
        end
    end

    @testset "Structure & 3D" begin
        # Dirty Data Pattern: alternate locations, missing backbone atoms, insertion codes,
        # and heterogens/waters mixed into the same coordinate file.
        # Expected Behavior: parsing should keep the structure usable, collapse altlocs
        # deterministically, and preserve HETATM vs ATOM annotations.
        @testset "PDB parsing and selection" begin
            pdb_text = join([
                "HEADER    DIRTY TEST STRUCTURE",
                @sprintf("%-6s%5d %-4s%c%3s %1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s", "ATOM", 1, "N", ' ', "ALA", "A", 100, 'A', 11.0, 12.0, 13.0, 1.00, 10.00, "N", ""),
                @sprintf("%-6s%5d %-4s%c%3s %1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s", "ATOM", 2, "CA", 'A', "ALA", "A", 100, 'A', 12.0, 12.5, 13.5, 0.60, 10.00, "C", ""),
                @sprintf("%-6s%5d %-4s%c%3s %1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s", "ATOM", 3, "CA", 'B', "ALA", "A", 100, 'A', 12.1, 12.6, 13.6, 0.40, 10.00, "C", ""),
                @sprintf("%-6s%5d %-4s%c%3s %1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s", "ATOM", 4, "C", ' ', "ALA", "A", 100, 'A', 13.0, 13.0, 14.0, 1.00, 10.00, "C", ""),
                @sprintf("%-6s%5d %-4s%c%3s %1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s", "ATOM", 5, "CA", ' ', "GLY", "A", 100, 'B', 14.0, 14.0, 15.0, 1.00, 10.00, "C", ""),
                @sprintf("%-6s%5d %-4s%c%3s %1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s", "ATOM", 6, "C", ' ', "GLY", "A", 100, 'B', 15.0, 15.0, 16.0, 1.00, 10.00, "C", ""),
                @sprintf("%-6s%5d %-4s%c%3s %1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s", "HETATM", 7, "O", ' ', "HOH", "A", 200, ' ', 16.0, 16.0, 17.0, 1.00, 10.00, "O", ""),
                "END",
            ], "\n")

            structure = BioToolkit.read_pdb(IOBuffer(pdb_text))
            chain = structure.models[1].chains[1]

            @test length(structure.models) == 1
            @test length(structure.models[1].chains) == 1
            @test length(chain.residues) == 3
            @test chain.residues[1].insertion_code == 'A'
            @test chain.residues[2].insertion_code == 'B'
            @test count(atom -> atom.hetatm, BioToolkit.structure_atoms(structure)) == 1
            @test length(BioToolkit.select_residues(structure; property=:protein)) == 2
            @test length(BioToolkit.select_residues(structure; property=:water)) == 1

            collapsed = BioToolkit.collapse_altlocs(structure)
            collapsed_residue = collapsed.models[1].chains[1].residues[1]
            @test count(atom -> strip(atom.name) == "CA", collapsed_residue.atoms) == 1
            @test any(atom -> strip(atom.name) == "CA" && atom.altloc == 'A', collapsed_residue.atoms)
            @test String(BioToolkit.sequence_from_structure(chain)) == "AGX"
        end
    end
end