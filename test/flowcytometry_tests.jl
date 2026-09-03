using Test
using BioToolkit
using LinearAlgebra
using Statistics

@testset "FlowCytometry Comprehensive Test Suite" begin
    @testset "FCS Data Reshaping (List Mode Parity)" begin
        n_events = 4
        n_channels = 3
        # FCS List Mode Binary Data:
        # Event 1: [1.0, 2.0, 3.0]
        # Event 2: [4.0, 5.0, 6.0]
        # Event 3: [7.0, 8.0, 9.0]
        # Event 4: [10.0, 11.0, 12.0]
        data_float32 = Float32[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
        raw_bytes = reinterpret(UInt8, data_float32)

        meta = Dict{String, Any}(
            "\$DATATYPE" => "F",
            "\$BYTEORD" => "1,2,3,4",
            "\$TOT" => "4",
            "\$PAR" => "3"
        )

        events = BioToolkit.FlowCytometry._read_fcs_data(Vector{UInt8}(raw_bytes), meta, n_events, n_channels)
        
        @test size(events) == (4, 3)
        @test events[1, :] == [1.0, 2.0, 3.0]
        @test events[2, :] == [4.0, 5.0, 6.0]
        @test events[3, :] == [7.0, 8.0, 9.0]
        @test events[4, :] == [10.0, 11.0, 12.0]

        # Test Float64 mode
        data_float64 = Float64[10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
        raw64 = reinterpret(UInt8, data_float64)
        meta64 = Dict{String, Any}("\$DATATYPE" => "D", "\$BYTEORD" => "1,2,3,4", "\$TOT" => "2", "\$PAR" => "3")
        events64 = BioToolkit.FlowCytometry._read_fcs_data(Vector{UInt8}(raw64), meta64, 2, 3)
        @test events64[1, :] == [10.0, 20.0, 30.0]
        @test events64[2, :] == [40.0, 50.0, 60.0]
    end

    @testset "Spillover Matrix Parsing & Auto-Compensation" begin
        chans = ["FITC-A", "PE-A", "APC-A"]
        meta = Dict{String, Any}(
            "\$SPILLOVER" => "3,FITC-A,PE-A,APC-A,1.0,0.02,0.001,0.05,1.0,0.01,0.01,0.03,1.0"
        )
        spill, parsed_chans = parse_spillover(meta)
        @test parsed_chans == chans
        @test size(spill) == (3, 3)
        @test spill[1, 1] == 1.0
        @test spill[1, 2] == 0.02
        @test spill[2, 1] == 0.05

        fcs = mock_flow_experiment(n_events=50, channels=chans)
        fcs.metadata["\$SPILLOVER"] = meta["\$SPILLOVER"]

        # Explicit compensation
        comp_fcs1 = compensate_fcs(fcs, spill)
        @test size(comp_fcs1.events) == (50, 3)
        @test all(isfinite, comp_fcs1.events)

        # Automatic zero-arg compensation from metadata
        comp_fcs2 = compensate_fcs(fcs)
        @test comp_fcs2.events ≈ comp_fcs1.events
    end

    @testset "Ellipsoid Gate (Triangular Solver & Regularization)" begin
        fcs = mock_flow_experiment(n_events=200, seed=42)
        cov_mat = [2.0 0.5; 0.5 1.5]
        center_vec = [15.0, 13.0]
        
        gate = ellipsoid_gate(fcs, ["FSC-A", "SSC-A"], center_vec, cov_mat; confidence=0.95)
        @test gate isa GateResult
        @test gate.gate_type === :ellipsoid
        @test length(gate.indices) > 0

        sub_fcs = apply_gate(fcs, gate)
        @test size(sub_fcs.events, 1) == length(gate.indices)
    end

    @testset "FlowSOM Metaclustering (High Grid Size Scalability)" begin
        fcs = mock_flow_experiment(n_events=500, seed=123)
        
        # Test standard grid size (10x10 = 100 codes)
        res10 = flowsom_cluster(fcs; n_metaclusters=5, grid_size=10, n_epochs=3)
        @test res10 isa FlowSOMResult
        @test length(res10.metaclusters) == 500
        @test length(unique(res10.metaclusters)) <= 5
        @test size(res10.cluster_centers, 1) <= 5

        # Test large grid size (30x30 = 900 codes) for O(N^2) optimized metaclustering
        t_start = time()
        res30 = flowsom_cluster(fcs; n_metaclusters=8, grid_size=30, n_epochs=1)
        t_elapsed = time() - t_start
        @test length(res30.metaclusters) == 500
        @test length(unique(res30.metaclusters)) <= 8
        @test t_elapsed < 5.0 # Fast performance guaranteed
    end

    @testset "XShift Clustering (RP-Tree & Refined Graph kNN)" begin
        fcs = mock_flow_experiment(n_events=1200, seed=777)
        
        # Exact kNN
        exact_clusters = xshift_cluster(fcs; K=15, approximate=false)
        @test length(exact_clusters) == 1200
        @test maximum(exact_clusters) >= 1

        # Approximate multi-tree RP-Tree + local join kNN
        approx_clusters = xshift_cluster(fcs; K=15, approximate=true)
        @test length(approx_clusters) == 1200
        @test maximum(approx_clusters) >= 1
    end

    @testset "Logicle Auto-Estimation & Transformation Utilities" begin
        fcs = mock_flow_experiment(n_events=100)
        
        # Single channel estimation
        params = estimate_logicle_params(fcs, "FITC-A")
        @test params.T >= 1000.0
        @test params.W > 0

        # All channels estimation dictionary
        all_params = estimate_logicle_params(fcs)
        @test length(all_params) == length(fcs.channels)
        @test haskey(all_params, "FITC-A")

        # Auto-estimate Logicle transform
        transformed = logicle_transform(fcs; auto_estimate=true)
        @test all(isfinite, transformed.events)

        # Standard Logicle transform & roundtrip
        trans_manual = logicle_transform(fcs; channels=["FITC-A"], T=params.T, W=params.W)
        inv_transformed = inverse_logicle_transform(trans_manual; channels=["FITC-A"], T=params.T, W=params.W)
        @test inv_transformed.events[:, 1] ≈ fcs.events[:, 1] atol=1e-3

        # Arcsinh transform
        asinh_fcs = arcsinh_transform(fcs; cofactor=5.0)
        @test all(isfinite, asinh_fcs.events)

        # Hlog transform
        hlog_fcs = hlog_transform(fcs; b=1.0)
        @test all(isfinite, hlog_fcs.events)
    end
end
