using Test
using SparseArrays
using DataFrames
using LinearAlgebra
using Graphs
using BioToolkit
using BioToolkit: SummarizedExperiment, GeneIdType, OrganismDb, TxQuantRecord, Tx2GeneMap, AnnotatedHeatmapSpec, ensembl, entrez, symbol, refseq, uniprot

@testset "Bioconductor Parity Expansion" begin

    @testset "Sprint 0: Foundations & SummarizedExperiment" begin
        # 1. Parametric SummarizedExperiment on any Real (Int and Float64)
        float_assays = Dict("counts" => [1.0 2.0; 3.0 4.0])
        int_assays = Dict("counts" => [1 2; 3 4])
        
        rowData = Dict(:gene_id => ["G1", "G2"])
        colData = Dict(:sample_id => ["S1", "S2"])
        
        se_float = SummarizedExperiment(float_assays, rowData, colData)
        se_int = SummarizedExperiment(int_assays, rowData, colData)
        
        @test se_float isa SummarizedExperiment{Float64}
        @test se_int isa SummarizedExperiment{Int}
        
        # 2. GeneIdType singleton and Symbol conversions
        @test ensembl isa GeneIdType
        @test entrez isa GeneIdType
        @test symbol isa GeneIdType
        @test refseq isa GeneIdType
        @test uniprot isa GeneIdType
        
        @test Symbol(ensembl) === :ensembl
        @test Symbol(entrez) === :entrez
        @test Symbol(symbol) === :symbol
        @test Symbol(refseq) === :refseq
        @test Symbol(uniprot) === :uniprot
        
        @test GeneIdType(:ensembl) === ensembl
        @test_throws ArgumentError GeneIdType(:invalid_id_type)
        
        # 3. OrganismDb, TxQuantRecord, Tx2GeneMap, AnnotatedHeatmapSpec
        org_db = OrganismDb(:human)
        @test org_db isa OrganismDb{:human}
        
        tx_quant = TxQuantRecord(["T1", "T2"], [10.0, 20.0], [5.0, 10.0], [100.0, 200.0], "salmon", "sample1")
        @test tx_quant.counts == [10.0, 20.0]
        
        tx2gene = Tx2GeneMap(Dict("T1" => "G1"), ensembl, Dict{Symbol,Any}())
        @test tx2gene.map["T1"] == "G1"
        
        heatmap_spec = AnnotatedHeatmapSpec()
        @test heatmap_spec.row_splits == String[]
    end

    @testset "Sprint 1: GeneAnnotation Module" begin
        # Download annotation DB
        db = download_annotation_db(:human)
        @test db isa OrganismAnnotationDb
        @test db.organism === :human
        
        # Convert IDs (including hops)
        conv = convert_ids(db, ["ENSG00000141510", "7157", "invalid_gene"], :ensembl, :symbol)
        @test conv isa IdConversionResult
        @test conv.mapped["ENSG00000141510"] == ["TP53"]
        @test conv.mapped["7157"] == ["TP53"]
        @test "invalid_gene" in conv.unmapped
        
        # Map to organism (orthology lift)
        lift = map_to_organism(db, ["TP53", "BRCA1"], :mouse)
        @test lift.mapped["TP53"] == ["Tp53"]
        @test lift.mapped["BRCA1"] == ["Brca1"]
        
        # Pathway mappings
        pathways = pathway_mappings(db, ["TP53", "GAPDH"])
        @test pathways isa PathwayMappingResult
        @test "GO:0008150" in pathways.mapped_pathways["TP53"]
        
        # Annotate SummarizedExperiment
        se = SummarizedExperiment(
            Dict("counts" => [10.0 20.0; 30.0 40.0]),
            Dict(:gene_id => ["ENSG00000141510", "ENSG00000012048"]),
            Dict(:sample_id => ["S1", "S2"])
        )
        se_ann = annotate_summarized_experiment(se; id_col=:gene_id, db=db, to=:symbol)
        @test se_ann.rowData[:symbol] == ["TP53", "BRCA1"]
        
        # Build EnrichmentDatabase bridge
        enrich_db = build_enrichment_database(db; db_type=:go)
        @test enrich_db isa EnrichmentDatabase
        
        # Save/load DB (JSON and Arrow formats)
        mktempdir() do dir
            json_path = joinpath(dir, "db.json")
            arrow_path = joinpath(dir, "db.arrow")
            
            save_annotation_db(db, json_path)
            save_annotation_db(db, arrow_path)
            
            @test isfile(json_path)
            @test isfile(arrow_path)
            
            db_json = load_annotation_db(json_path)
            db_arrow = load_annotation_db(arrow_path)
            
            @test db_json.organism === :human
            @test db_arrow isa OrganismAnnotationDb
        end
    end

    @testset "Sprint 1: TxImport Module" begin
        # Setup mock Salmon quant file
        mktempdir() do dir
            quant_file = joinpath(dir, "quant.sf")
            open(quant_file, "w") do io
                println(io, "Name\tLength\tEffectiveLength\tTPM\tNumReads")
                println(io, "ENST00000269305\t1000\t800.0\t50.0\t100.0")
                println(io, "ENST00000357654\t2000\t1800.0\t30.0\t200.0")
            end
            
            # 1. Read Salmon output
            record = read_salmon(quant_file; sample_id="S1")
            @test record isa TxQuantRecord
            @test record.sample_id == "S1"
            @test record.tx_ids == ["ENST00000269305", "ENST00000357654"]
            @test record.counts == [100.0, 200.0]
            
            # 2. Fetch transcript to gene map
            tx2gene = fetch_tx2gene(:human)
            @test tx2gene isa Tx2GeneMap
            
            # 3. tximport aggregation (mode :no)
            res_no = tximport([quant_file], :salmon, tx2gene; countsFromAbundance=:no, sample_ids=["S1"])
            @test res_no isa TxImportResult
            @test res_no.gene_counts.gene_ids == ["ENSG00000141510", "ENSG00000012048"]
            @test Array(res_no.gene_counts.counts) == [100; 200;;]
            
            # 4. tximport aggregation (mode :scaledTPM)
            res_scaled = tximport([quant_file], :salmon, tx2gene; countsFromAbundance=:scaledTPM, sample_ids=["S1"])
            @test res_scaled isa TxImportResult
            @test res_scaled.gene_counts.gene_ids == ["ENSG00000141510", "ENSG00000012048"]
            
            # 5. tximport aggregation (mode :lengthScaledTPM)
            res_len = tximport([quant_file], :salmon, tx2gene; countsFromAbundance=:lengthScaledTPM, sample_ids=["S1"])
            @test res_len isa TxImportResult
            @test res_len.gene_counts.gene_ids == ["ENSG00000141510", "ENSG00000012048"]
        end
    end

    @testset "Sprint 2: ChIP-seq Workflow" begin
        # 1. cross_correlation_profile & nsc_rsc_metrics
        intervals = [
            GenomicInterval("chr1", 100, 150, '+'),
            GenomicInterval("chr1", 120, 170, '-'),
            GenomicInterval("chr1", 500, 550, '+')
        ]
        profile = cross_correlation_profile(intervals, "chr1", 135; bin_size=1, max_shift=100)
        @test profile isa StrandCorrelationProfile
        ccr = nsc_rsc_metrics(profile)
        @test ccr isa CrossCorrelationResult
        @test ccr.nsc >= 0
        @test ccr.quality_tag in ("High", "Medium", "Low")

        
        # 2. IDR Analysis
        peak1 = PeakSet([Peak("chr1", 100, 150, 10.0), Peak("chr1", 200, 250, 5.0)])
        peak2 = PeakSet([Peak("chr1", 105, 155, 9.0), Peak("chr1", 210, 260, 4.0)])
        idr_res = idr_analysis([peak1, peak2])
        @test idr_res isa IDRResult
        @test length(idr_res.peak_pairs) == 2
        consensus = consensus_peaks(idr_res)
        @test consensus isa PeakSet
        
        # 3. Motif search/scan
        seqs = [DNASeq("ACGTAC"), DNASeq("ACGTAC"), DNASeq("ACGTAC")]
        counts = motif_counts(seqs)
        pwm = motif_pwm(counts; pseudocount=0.5)
        @test pwm isa MotifPWM
        seq = DNASeq("ACGTACGTACGT")
        hits = motif_scan(seq, pwm; threshold=0.0)
        @test length(hits) >= 0
    end

    @testset "Sprint 3: Single-cell Hurdle DE & Pseudobulk" begin
        # 1. mast_hurdle_test — returns Vector{HurdleDEResult}
        n_genes, n_cells = 20, 10
        X = rand(n_genes, n_cells)
        X[X .< 0.5] .= 0.0
        groups = [fill(:A, 5); fill(:B, 5)]
        results = mast_hurdle_test(X, groups)
        @test results isa Vector{HurdleDEResult}
        @test length(results) == n_genes

        # 2. find_markers via integer labels
        counts_mat = sparse(round.(Int, rand(10, 8) .* 5))
        gene_ids = ["G$i" for i in 1:10]
        cell_ids = ["C$i" for i in 1:8]
        sce = SingleCellExperiment(counts_mat, gene_ids, cell_ids;
                                   metadata=Dict("group" => [1,1,1,1,2,2,2,2]))
        labels = [1,1,1,1,2,2,2,2]
        markers = find_markers(sce, labels; ident_1=1)
        @test !isempty(markers)

        # 3. Pseudobulk — requires donor_col & celltype_col kwargs
        pb = pseudobulk_by_donor_celltype(sce; donor_col=:group, celltype_col=:group)
        @test pb isa CountMatrix
    end

    @testset "Sprint 4: Flow Cytometry" begin
        # 1. read_fcs & compensate_fcs
        fcs = mock_flow_experiment()
        @test fcs isa FlowExperiment
        spill = Matrix{Float64}(I, 5, 5)
        comp = compensate_fcs(fcs, spill)
        @test size(comp.events) == size(fcs.events)
        
        # 2. Gating
        rect = rectangle_gate(fcs, "FSC-A", 2.0, 8.0, "SSC-A", 1.0, 9.0)
        @test rect isa GateResult
        sub_fcs = apply_gate(fcs, rect)
        @test size(sub_fcs.events, 1) <= size(fcs.events, 1)
        
        # 3. Clustering
        fs_res = flowsom_cluster(fcs; n_metaclusters=3)
        @test fs_res isa FlowSOMResult
        xs_labels = xshift_cluster(fcs; K=5)
        @test length(xs_labels) == size(fcs.events, 1)
    end

    @testset "Sprint 5: Advanced Bioplotting" begin
        # 1. annotated_heatmap
        m = rand(10, 10)
        df_row = DataFrame(Type=["A", "B", "A", "B", "A", "B", "A", "B", "A", "B"])
        df_col = DataFrame(Group=["G1", "G1", "G1", "G1", "G1", "G2", "G2", "G2", "G2", "G2"])
        spec = AnnotatedHeatmapSpec(row_annotations=df_row, col_annotations=df_col)
        res_hm = annotated_heatmap(m, spec)
        @test res_hm isa AnnotatedHeatmapResult
        
        # 2. upset_plot
        sets = Dict("S1" => [1, 2, 3], "S2" => [2, 3, 4], "S3" => [3, 4, 5])
        res_us = upset_plot(sets)
        @test res_us isa UpsetPlotResult
        
        # 3. oncoprint_from_matrix
        muts = ["TP53" "Nonsense" ""; "KRAS" "" "Missense"]
        res_op = oncoprint_from_matrix(muts)
        @test res_op isa OncoprintResult
        
        # 4. circos_plot
        res_circ = circos_plot(m)
        @test res_circ isa CircosPlotResult
    end

    @testset "Sprint 6: Proteomics" begin
        # 1. kinase_activity_inference
        phospho = DataFrame(site=["S1", "S2", "S3"], log2_fc=[2.0, -1.0, 0.5])
        db = DataFrame(kinase=["K1", "K1", "K2"], site=["S1", "S2", "S3"])
        res = kinase_activity_inference(phospho, db)
        @test res isa KinaseActivityResult
        @test length(res.kinases) == 2
        
        # 2. phosphoproteomics_qc
        spectra = [Spectrum(1.0, [98.0, 200.0], [50.0, 100.0]), Spectrum(2.0, [150.0, 300.0], [40.0, 80.0])]
        ms = MassSpecExperiment(spectra, Dict{String,Any}())
        qc = phosphoproteomics_qc(ms; localization_threshold=0.8)
        @test qc isa PhosphoQCResult
        @test qc.summary[qc.summary.Metric .== "Localized Sites", :Value][1] == 1.0
        
        # 3. normalize_label_free
        mat = [1.0 2.0; 3.0 4.0; 5.0 6.0]
        norm_med = normalize_label_free(mat; method=:median)
        @test norm_med isa Matrix{Float64}
        norm_q = normalize_label_free(mat; method=:quantile)
        @test norm_q isa Matrix{Float64}
        norm_vsn = normalize_label_free(mat; method=:vsn)
        @test norm_vsn isa Matrix{Float64}
        
        # 4. dspikein_calibration
        spikes = [10.0, 100.0]
        samples = [5.0 15.0; 50.0 150.0]
        concs = [1.0, 10.0]
        calib = dspikein_calibration(spikes, samples, concs)
        @test size(calib) == (2, 2)
    end

    @testset "Sprint 7: Systems Biology" begin
        # 1. flux_balance_analysis
        S = [1.0 -1.0 0.0; 0.0 1.0 -1.0]
        lb = [0.0, 0.0, 0.0]
        ub = [10.0, 10.0, 10.0]
        c = [0.0, 0.0, 1.0]
        res = flux_balance_analysis(S, lb, ub, c)
        @test res.status == "converged"
        @test length(res.fluxes) == 3
        
        # 2. grn_reconstruction_mi
        expr = rand(5, 20)
        mi_pearson = grn_reconstruction_mi(expr; estimator=:pearson)
        @test size(mi_pearson) == (5, 5)
        mi_spearman = grn_reconstruction_mi(expr; estimator=:spearman)
        @test size(mi_spearman) == (5, 5)
        mi_kraskov = grn_reconstruction_mi(expr; estimator=:kraskov)
        @test size(mi_kraskov) == (5, 5)
        
        # 3. protein_interaction_network_alignment
        g1 = SimpleGraph(3)
        add_edge!(g1, 1, 2)
        add_edge!(g1, 2, 3)
        g2 = SimpleGraph(3)
        add_edge!(g2, 1, 3)
        add_edge!(g2, 3, 2)
        res_align = protein_interaction_network_alignment(g1, g2)
        @test length(res_align.alignment_map) > 0
    end

    @testset "Provenance Auditing" begin
        # Enable provenance context
        with_provenance() do
            ctx = BioToolkit.scoped_provenance_context()
            db = download_annotation_db(:human)
            tx2gene = fetch_tx2gene(db)
            
            # Check DAG registration
            @test haskey(ctx.nodes, db.provenance.id)
            @test haskey(ctx.nodes, tx2gene.provenance.id)
            
            # Check parents relationship
            tx2gene_node = ctx.nodes[tx2gene.provenance.id]
            @test db.provenance.id in tx2gene_node.parent_ids
        end
    end

end
