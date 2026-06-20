using Test
using DataFrames
using BioToolkit

function _node_by_operation(ctx::BioToolkit.ProvenanceContext, operation::AbstractString)
    matches = [node for node in values(ctx.nodes) if node.operation == operation]
    @test !isempty(matches)
    return first(matches)
end

@testset "Pipeline orchestration DAG" begin
    @testset "Graph validation and scheduling" begin
        graph = PipelineGraph(cache_dir=mktempdir())
        add_node!(graph, pipeline_node(:a, () -> 2; cache=false))
        add_node!(graph, pipeline_node(:b, x -> x + 1; dependencies=[:a], cache=false))
        add_node!(graph, pipeline_node(:c, x -> 3x; dependencies=[:a], cache=false))
        add_node!(graph, pipeline_node(:d, (b, c) -> b + c; dependencies=[:b, :c], cache=false))

        @test validate_pipeline(graph)
        levels = execution_levels(graph)
        @test length(levels) == 3
        @test levels[1] == [:a]
        @test Set(levels[2]) == Set([:b, :c])
        @test levels[3] == [:d]

        outputs = execute_pipeline(graph; parallel=true)
        @test outputs[:d] == 9
    end

    @testset "Cycle detection" begin
        cyc = PipelineGraph(cache_dir=mktempdir())
        add_node!(cyc, pipeline_node(:x, () -> 1; cache=false))
        add_node!(cyc, pipeline_node(:y, x -> x + 1; dependencies=[:x], cache=false))
        add_dependency!(cyc, :x, :y)
        @test_throws ArgumentError validate_pipeline(cyc)
    end

    @testset "Caching behavior and file artifacts" begin
        cache_dir = mktempdir()
        cache_graph = PipelineGraph(cache_dir=cache_dir)
        eval_count = Ref(0)

        add_node!(cache_graph, pipeline_node(:expensive, () -> begin
            eval_count[] += 1
            return 42
        end; cache=true))

        first_run = execute_pipeline(cache_graph; parallel=false)
        second_run = execute_pipeline(cache_graph; parallel=false)

        @test first_run[:expensive] == 42
        @test second_run[:expensive] == 42
        @test eval_count[] == 1

        temp_file = joinpath(cache_dir, "artifact.txt")
        write(temp_file, "alpha")
        artifact_a = file_artifact(temp_file)
        write(temp_file, "alpha-beta")
        artifact_b = file_artifact(temp_file)

        @test artifact_a.path == temp_file
        @test artifact_a.md5 != ""
        @test artifact_a.md5 != artifact_b.md5
        @test artifact_b.size > artifact_a.size
    end

    @testset "Provenance tracking" begin
        ctx = BioToolkit.ProvenanceContext()
        graph = PipelineGraph(cache_dir=mktempdir(); prov_ctx=ctx)
        node = pipeline_node(:root, () -> 7; cache=false, prov_ctx=ctx)
        add_node!(graph, node; prov_ctx=ctx)

        outputs = execute_pipeline(graph; parallel=false, prov_ctx=ctx)
        @test outputs[:root] == 7

        temp_file = joinpath(mktempdir(), "artifact.txt")
        write(temp_file, "alpha")
        artifact = file_artifact(temp_file; prov_ctx=ctx)
        @test artifact.md5 != ""

        plan = slurm_array_plan(DataFrame(sample_id=["s1"], value=[1]); prov_ctx=ctx)
        @test nrow(plan) == 1

        exec_node = _node_by_operation(ctx, "execute_pipeline")
        artifact_node = _node_by_operation(ctx, "file_artifact")
        plan_node = _node_by_operation(ctx, "slurm_array_plan")

        @test exec_node.parameters["node_id"] == "root"
        @test artifact_node.parameters["status"] == "ok"
        @test plan_node.parameters["row_count"] == 1
    end

    @testset "Template wrappers" begin
        template_graph = PipelineGraph(cache_dir=mktempdir())

        add_node!(template_graph, pipeline_node(:reads, () -> ["GATCGTAGCTA", "ATCGTAGCTAA"]; cache=false))
        add_node!(template_graph, pipeline_node(:reference, () -> "TTTTTGATCGTAGCTAGGG"; cache=false))
        add_node!(template_graph, template_align_reads_node(:alignments; dependencies=[:reads, :reference], cache=false, bandwidth=8))

        add_node!(template_graph, pipeline_node(:features, () -> [("geneA", 6, 14), ("geneB", 15, 20)]; cache=false))
        add_node!(template_graph, template_count_features_node(:feature_counts; dependencies=[:alignments, :features], cache=false))

        outputs = execute_pipeline(template_graph; parallel=false)
        @test length(outputs[:alignments]) == 2
        @test outputs[:feature_counts]["geneA"] >= 1

        de_graph = PipelineGraph(cache_dir=mktempdir())
        counts = [30 28 10 12; 8 7 25 23; 18 17 19 20; 40 38 12 11]
        design = ["A", "A", "B", "B"]
        add_node!(de_graph, pipeline_node(:count_matrix, () -> counts; cache=false))
        add_node!(de_graph, pipeline_node(:design, () -> design; cache=false))
        add_node!(de_graph, template_differential_expression_node(:de; dependencies=[:count_matrix, :design], cache=false, min_total=0, fitType=:local))

        de_outputs = execute_pipeline(de_graph; parallel=false)
        @test de_outputs[:de] isa Vector
    end
end
