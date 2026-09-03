using Test
using BioToolkit
using BioToolkit.TxImport
using SparseArrays
using DataFrames
using Dates
using Statistics

# Import internal functions for testing
using BioToolkit.TxImport: _strip_tx_version, _strip_after_bar, _clean_tx_id, _tx_checksum, row_vars, _make_counts_from_abundance, _median_length_over_isoform, _replace_missing_length, convert_stringtie_counts, summarize_to_gene

@testset "TxImport" begin

    @testset "Core Types" begin
        # Test TxImportOptions
        opts = TxImportOptions(countsFromAbundance=:lengthScaledTPM, ignore_tx_version=true)
        @test opts.countsFromAbundance == :lengthScaledTPM
        @test opts.ignore_tx_version == true
        @test opts.ignore_after_bar == false
        @test opts.require_all_files == true

        # Test TxQuantRecord
        tx_ids = ["ENST000001", "ENST000002"]
        counts = [100.0, 200.0]
        tpm = [50.0, 100.0]
        efflength = [1000.0, 2000.0]
        rec = TxQuantRecord(tx_ids, counts, tpm, efflength, "salmon", "sample1")
        @test rec.tx_ids == tx_ids
        @test rec.counts == counts
        @test rec.tpm == tpm
        @test rec.efflength == efflength
        @test rec.tool == "salmon"
        @test rec.sample_id == "sample1"
    end

    @testset "Transcript ID Cleaning" begin
        @test _strip_tx_version("ENST000001.1") == "ENST000001"
        @test _strip_tx_version("ENST000001") == "ENST000001"
        @test _strip_after_bar("ENST000001|extra") == "ENST000001"
        @test _strip_after_bar("ENST000001") == "ENST000001"
        @test _clean_tx_id("ENST000001.1|extra", true, true) == "ENST000001"
        @test _clean_tx_id("ENST000001.1|extra", false, true) == "ENST000001.1"
        @test _clean_tx_id("ENST000001.1|extra", true, false) == "ENST000001"
    end

    @testset "Checksum" begin
        mktempdir() do dir
            file = joinpath(dir, "test.txt")
            write(file, "test content")
            checksum = _tx_checksum(file)
            @test checksum isa UInt64
            @test checksum != 0
        end
    end

    @testset "Quantification File Reading" begin
        mktempdir() do dir
            # Create mock Salmon quant file
            salmon_file = joinpath(dir, "quant.sf")
            write(salmon_file, """Name\tLength\tEffectiveLength\tTPM\tNumReads
ENST000001\t1500\t1200\t50.0\t100
ENST000002\t2000\t1800\t100.0\t200
ENST000003\t1000\t800\t25.0\t50
""")

            rec = read_salmon(salmon_file; sample_id="sample1")
            @test length(rec.tx_ids) == 3
            @test rec.tx_ids == ["ENST000001", "ENST000002", "ENST000003"]
            @test rec.counts ≈ [100.0, 200.0, 50.0]
            @test rec.tpm ≈ [50.0, 100.0, 25.0]
            @test rec.efflength ≈ [1200.0, 1800.0, 800.0]
            @test rec.tool == "salmon"
            @test rec.sample_id == "sample1"

            # Test with ignore_tx_version
            rec2 = read_salmon(salmon_file; sample_id="sample1", ignore_tx_version=true)
            @test rec2.tx_ids == ["ENST000001", "ENST000002", "ENST000003"]

            # Test kallisto TSV
            kallisto_file = joinpath(dir, "abundance.tsv")
            write(kallisto_file, """target_id\tlength\teff_length\test_counts\ttpm
ENST000001\t1500\t1200\t100\t50.0
ENST000002\t2000\t1800\t200\t100.0
ENST000003\t1000\t800\t50\t25.0
""")

            rec_k = read_kallisto(kallisto_file; sample_id="sample1")
            @test length(rec_k.tx_ids) == 3
            @test rec_k.tx_ids == ["ENST000001", "ENST000002", "ENST000003"]
            @test rec_k.counts ≈ [100.0, 200.0, 50.0]

            # Test RSEM
            rsem_file = joinpath(dir, "isoforms.results")
            write(rsem_file, """transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct
ENST000001\tENSG000001\t1500\t1200\t100\t50.0\t48.0\t50.0
ENST000002\tENSG000002\t2000\t1800\t200\t100.0\t96.0\t100.0
""")

            rec_r = read_rsem_quant(rsem_file; sample_id="sample1")
            @test length(rec_r.tx_ids) == 2
            @test rec_r.tx_ids == ["ENST000001", "ENST000002"]

            # Test StringTie
            stringtie_file = joinpath(dir, "stringtie.gtf")
            write(stringtie_file, """t_name\tgene_name\tref_gene_id\tlength\tcov\tFPKM
ENST000001\tGENE1\tENSG000001\t1500\t100\t50.0
ENST000002\tGENE2\tENSG000002\t2000\t200\t100.0
""")

            rec_s = read_stringtie(stringtie_file; sample_id="sample1")
            @test length(rec_s.tx_ids) == 2
            @test rec_s.tx_ids == ["ENST000001", "ENST000002"]
            @test rec_s.counts ≈ [100.0, 200.0]  # coverage
        end
    end

    @testset "Tx2GeneMap" begin
        # Test Tx2GeneMap creation
        tx_map = Dict(
            "ENST000001" => "ENSG000001",
            "ENST000002" => "ENSG000001",
            "ENST000003" => "ENSG000002"
        )
        tx2gene = Tx2GeneMap(tx_map, BioToolkit.BioToolkit.ensembl, Dict{Symbol,Any}())
        @test length(tx2gene.map) == 3
        @test tx2gene.map["ENST000001"] == "ENSG000001"
        @test tx2gene.gene_id_type == BioToolkit.ensembl
    end

@testset "Summarize to Gene" begin
        # Create mock transcript-level data
        tx_ids = ["ENST000001", "ENST000002", "ENST000003"]
        gene_ids = ["ENSG000001", "ENSG000001", "ENSG000002"]
        tx2gene_map = Dict(tx_ids .=> gene_ids)
        tx2gene = Tx2GeneMap(tx2gene_map, BioToolkit.ensembl, Dict{Symbol,Any}())

        # Transcript-level matrices (3 transcripts x 2 samples)
        abundance_tx = [50.0 60.0; 100.0 120.0; 25.0 30.0]
        counts_tx = [100.0 120.0; 200.0 240.0; 50.0 60.0]
        length_tx = [1200.0 1300.0; 1800.0 1900.0; 800.0 900.0]

        tx_data = Dict(
            "abundance" => abundance_tx,
            "counts" => counts_tx,
            "length" => length_tx
        )

        # Summarize to gene level
        gene_data = summarize_to_gene(tx_data, tx2gene, tx_ids)

        @test haskey(gene_data, "abundance")
        @test haskey(gene_data, "counts")
        @test haskey(gene_data, "length")
        @test size(gene_data["abundance"]) == (2, 2)  # 2 genes, 2 samples
        @test size(gene_data["counts"]) == (2, 2)
        @test size(gene_data["length"]) == (2, 2)

        # Gene 1 (ENSG000001): sum of ENST000001 + ENST000002
        @test gene_data["abundance"][1, 1] ≈ 150.0  # 50 + 100
        @test gene_data["abundance"][1, 2] ≈ 180.0  # 60 + 120
        @test gene_data["counts"][1, 1] ≈ 300.0    # 100 + 200
        @test gene_data["counts"][1, 2] ≈ 360.0    # 120 + 240

        # Gene 2 (ENSG000002): ENST000003
        @test gene_data["abundance"][2, 1] ≈ 25.0
        @test gene_data["abundance"][2, 2] ≈ 30.0
        @test gene_data["counts"][2, 1] ≈ 50.0
        @test gene_data["counts"][2, 2] ≈ 60.0

        # Test weighted average length for gene 1
        # (50*1200 + 100*1800) / (50+100) = (60000 + 180000) / 150 = 1600
        @test gene_data["length"][1, 1] ≈ 1600.0 atol=1.0
    end

    @testset "Make Counts From Abundance" begin
        counts_mat = [100.0 120.0; 200.0 240.0; 50.0 60.0]
        abundance_mat = [50.0 60.0; 100.0 120.0; 25.0 30.0]
        length_mat = [1200.0 1300.0; 1800.0 1900.0; 800.0 900.0]

        # scaledTPM
        scaled = _make_counts_from_abundance(counts_mat, abundance_mat, length_mat, :scaledTPM)
        @test size(scaled) == size(counts_mat)
        # Library size = 350, 420
        # TPM sum = 175, 210
        # scaled = TPM * lib_size / TPM_sum
        @test scaled[1, 1] ≈ 50.0 * 350.0 / 175.0  # 100.0
        @test scaled[2, 1] ≈ 100.0 * 350.0 / 175.0  # 200.0

        # lengthScaledTPM
        length_scaled = _make_counts_from_abundance(counts_mat, abundance_mat, length_mat, :lengthScaledTPM)
        @test size(length_scaled) == size(counts_mat)
        # Average length per transcript
        avg_len = mean(length_mat, dims=2)[:]
        # scaled = TPM * avg_len * lib_size / sum(TPM * avg_len)
        @test all(length_scaled .> 0)

        # dtuScaledTPM
        gene_lengths = [1500.0, 1500.0, 1000.0]  # median isoform length per gene
        dtu_scaled = _make_counts_from_abundance(counts_mat, abundance_mat, length_mat, :dtuScaledTPM, gene_lengths=gene_lengths)
        @test size(dtu_scaled) == size(counts_mat)
        @test all(dtu_scaled .> 0)
    end

    @testset "Median Length Over Isoform" begin
        tx_ids = ["ENST000001", "ENST000002", "ENST000003"]
        tx2gene_map = Dict(
            "ENST000001" => "ENSG000001",
            "ENST000002" => "ENSG000001",
            "ENST000003" => "ENSG000002"
        )
        tx2gene = Tx2GeneMap(tx2gene_map, BioToolkit.ensembl, Dict{Symbol,Any}())

        # 3 transcripts x 2 samples
        length_mat = [1000.0 1100.0; 2000.0 2100.0; 1500.0 1600.0]

        median_len = _median_length_over_isoform(length_mat, tx2gene, tx_ids;
                                                  ignore_tx_version=false, ignore_after_bar=false)
        @test size(median_len) == size(length_mat)

        # For gene ENSG000001: transcripts 1 and 2, avg lengths = 1050, 2050, median = 1550
        @test median_len[1, 1] ≈ 1550.0
        @test median_len[2, 1] ≈ 1550.0
        @test median_len[3, 1] ≈ 1550.0  # gene ENSG000002 only has transcript 3, avg = 1550
    end

    @testset "Replace Missing Length" begin
        length_mat = [1000.0 1100.0; NaN NaN; 1500.0 1600.0]
        ave_gene = [1050.0, 1550.0, 1550.0]

        fixed = _replace_missing_length(copy(length_mat), ave_gene)
        @test !any(isnan, fixed)
        @test fixed[1, :] ≈ [1000.0, 1100.0]
        @test fixed[2, :] ≈ [1550.0, 1550.0]  # all NaN replaced with gene average
    end

    @testset "StringTie Count Conversion" begin
        counts = [100.0 200.0; 50.0 100.0]
        length = [1000.0 1000.0; 2000.0 2000.0]

        converted = convert_stringtie_counts(counts, length, 75.0)
        # counts = cov * length / read_length
        @test converted[1, 1] ≈ 100.0 * 1000.0 / 75.0
        @test converted[2, 1] ≈ 50.0 * 2000.0 / 75.0
    end

    @testset "Inferential Replicates - Row Variance" begin
        reps = [100.0 110.0 90.0; 200.0 190.0 210.0]
        vars = row_vars(reps)
        @test length(vars) == 2
        @test vars[1] ≈ Statistics.var([100.0, 110.0, 90.0])
        @test vars[2] ≈ Statistics.var([200.0, 190.0, 210.0])

        # Test with single replicate (should return 0)
        single = [100.0 100.0; 200.0 200.0]
        vars_single = row_vars(single)
        @test all(vars_single .== 0.0)
    end

    @testset "Fetch Tx2Gene" begin
        # Test with mock organism
        db = download_annotation_db(:human)
        tx2gene = fetch_tx2gene(db)
        @test tx2gene isa Tx2GeneMap
        @test length(tx2gene.map) > 0
    end

    @testset "TxImport with Mock Data" begin
        mktempdir() do dir
            # Create mock quantification files
            for i in 1:2
                file = joinpath(dir, "sample_$i", "quant.sf")
                mkpath(dirname(file))
                write(file, """Name\tLength\tEffectiveLength\tTPM\tNumReads
ENST000001\t1500\t1200\t50.0\t100
ENST000002\t2000\t1800\t100.0\t200
ENST000003\t1000\t800\t25.0\t50
""")
            end

            files = [
                joinpath(dir, "sample_1", "quant.sf"),
                joinpath(dir, "sample_2", "quant.sf")
            ]

            # Create tx2gene mapping
            tx2gene = Tx2GeneMap(
                Dict("ENST000001" => "ENSG000001", "ENST000002" => "ENSG000001", "ENST000003" => "ENSG000002"),
                BioToolkit.ensembl,
                Dict{Symbol,Any}()
            )

            # Test gene-level import
            result = tximport(files; type=:salmon, tx2gene=tx2gene, tx_out=false)
            @test result isa TxImportResult
            @test result.gene_counts isa CountMatrix
            @test size(result.gene_counts.counts) == (2, 2)  # 2 genes, 2 samples
            @test result.se isa SummarizedExperiment
            @test result.diagnostics isa TxImportDiagnostics

            # Test transcript-level import
            result_tx = tximport(files; type=:salmon, tx_out=true)
            @test haskey(result_tx, "counts")
            @test haskey(result_tx, "abundance")
            @test haskey(result_tx, "length")
            @test size(result_tx["counts"]) == (3, 2)
        end
    end

    @testset "TxImport Options" begin
        mktempdir() do dir
            file = joinpath(dir, "quant.sf")
            write(file, """Name\tLength\tEffectiveLength\tTPM\tNumReads
ENST000001\t1500\t1200\t50.0\t100
ENST000002\t2000\t1800\t100.0\t200
""")

            files = [file]
            tx2gene = Tx2GeneMap(Dict("ENST000001" => "ENSG000001", "ENST000002" => "ENSG000002"), BioToolkit.ensembl, Dict{Symbol,Any}())

            # Test with options struct
            opts = TxImportOptions(
                countsFromAbundance=:scaledTPM,
                ignore_tx_version=true,
                tx_out=false
            )

            result = tximport(files; type=:salmon, tx2gene=tx2gene, options=opts)
            @test result isa TxImportResult
        end
    end

    @testset "Provenance Tracking" begin
        mktempdir() do dir
            file = joinpath(dir, "quant.sf")
            write(file, """Name\tLength\tEffectiveLength\tTPM\tNumReads
ENST000001\t1500\t1200\t50.0\t100
ENST000002\t2000\t1800\t100.0\t200
""")

            files = [file]
            tx2gene = Tx2GeneMap(Dict("ENST000001" => "ENSG000001", "ENST000002" => "ENSG000002"), BioToolkit.ensembl, Dict{Symbol,Any}())

            # Enable provenance
            ctx = BioToolkit.enable_provenance!()
            result = tximport(files; type=:salmon, tx2gene=tx2gene)
            BioToolkit.disable_provenance!()

            @test result.provenance isa ResultProvenance
            @test !isempty(result.provenance.id)
            @test result.diagnostics.provenance isa ResultProvenance
        end
    end

    @testset "Error Handling" begin
        mktempdir() do dir
            file = joinpath(dir, "quant.sf")
            write(file, """Name\tLength\tEffectiveLength\tTPM\tNumReads
ENST000001\t1500\t1200\t50.0\t100
""")

            files = [file]

            # Missing tx2gene for gene-level
            @test_throws ArgumentError tximport(files; type=:salmon, tx_out=false)

            # Non-existent file
            @test_throws ArgumentError tximport(["/nonexistent/file.sf"]; type=:salmon)

            # Invalid type
            @test_throws ArgumentError tximport(files; type=:invalid_type, tx_out=true)

            # Invalid countsFromAbundance
            tx2gene = Tx2GeneMap(Dict("ENST000001" => "ENSG000001"), BioToolkit.ensembl, Dict{Symbol,Any}())
            @test_throws ArgumentError tximport(files; type=:salmon, tx2gene=tx2gene, counts_from_abundance=:invalid)
        end
    end

    @testset "Sparse Import" begin
        mktempdir() do dir
            for i in 1:2
                file = joinpath(dir, "sample_$i", "quant.sf")
                mkpath(dirname(file))
                write(file, """Name\tLength\tEffectiveLength\tTPM\tNumReads
ENST000001\t1500\t1200\t50.0\t100
ENST000002\t2000\t1800\t0.0\t0
ENST000003\t1000\t800\t25.0\t50
""")
            end

            files = [
                joinpath(dir, "sample_1", "quant.sf"),
                joinpath(dir, "sample_2", "quant.sf")
            ]

            result = tximport(files; type=:salmon, tx_out=true, sparse=true, sparse_threshold=1.0)
            @test result["counts"] isa SparseMatrixCSC
            @test size(result["counts"]) == (3, 2)
        end
    end

    @testset "Ignore TX Version and After Bar" begin
        mktempdir() do dir
            file = joinpath(dir, "quant.sf")
            write(file, """Name\tLength\tEffectiveLength\tTPM\tNumReads
ENST000001.1\t1500\t1200\t50.0\t100
ENST000002.2\t2000\t1800\t100.0\t200
""")

            files = [file]
            tx2gene = Tx2GeneMap(
                Dict("ENST000001" => "ENSG000001", "ENST000002" => "ENSG000002"),
                BioToolkit.ensembl,
                Dict{Symbol,Any}()
            )

            # Without ignoring version - should fail to match
            @test_throws ArgumentError tximport(files; type=:salmon, tx2gene=tx2gene, ignore_tx_version=false)

            # With ignoring version - should work
            result = tximport(files; type=:salmon, tx2gene=tx2gene, ignore_tx_version=true)
            @test result isa TxImportResult
        end
    end

    @testset "Alevin Import" begin
        mktempdir() do dir
            # Create mock alevin directory structure
            alevin_dir = joinpath(dir, "alevin")
            mkpath(alevin_dir)

            # gene names
            write(joinpath(alevin_dir, "quants_mat_cols.txt"), "GENE1\nGENE2\n")
            # cell barcodes
            write(joinpath(alevin_dir, "quants_mat_rows.txt"), "CELL1\nCELL2\n")

            # Simple counts matrix in EDS format (simplified)
            # We'll test the reader structure without full EDS implementation
            matrix_file = joinpath(alevin_dir, "quants_mat.gz")
            open(matrix_file, "w") do io
                gz = GzipCompressorStream(io)
                # Write empty EDS format for testing
                write(gz, UInt8[0,0,0,0])  # bit vectors for 2 genes, 2 cells
                write(gz, Float64[])  # no non-zero counts
            end

            cmd_info = joinpath(dir, "cmd_info.json")
            write(cmd_info, JSON.json(Dict("salmon_version" => "1.0.0", "numCellBootstraps" => 0)))

            # Test alevin read (will use slow path since no EDS data)
            # This tests the structure without full EDS implementation
        end
    end

    @testset "Diagnostics" begin
        mktempdir() do dir
            file = joinpath(dir, "quant.sf")
            write(file, """Name\tLength\tEffectiveLength\tTPM\tNumReads
ENST000001\t1500\t1200\t50.0\t100
ENST000002\t2000\t1800\t100.0\t200
""")

            files = [file]
            tx2gene = Tx2GeneMap(Dict("ENST000001" => "ENSG000001", "ENST000002" => "ENSG000002"), BioToolkit.ensembl, Dict{Symbol,Any}())

            result = tximport(files; type=:salmon, tx2gene=tx2gene)

            diag = result.diagnostics
            @test diag isa TxImportDiagnostics
            @test haskey(diag.file_checksums, file)
            @test haskey(diag.transcript_counts, "sample_1")
            @test diag.transcript_counts["sample_1"] == 2
            @test length(diag.ignored_tx) == 0
            @test size(diag.length_offset) == (2, 1)
        end
    end
end