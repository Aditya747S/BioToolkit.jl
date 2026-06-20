module PipelineBenchmarkFixtures

using BioToolkit

export pipeline_fixture, evaluate_pipeline_recovery

function _counted(counters::Dict{Symbol,Int}, node_id::Symbol, f::Function)
    return function(args...; kwargs...)
        counters[node_id] = get(counters, node_id, 0) + 1
        return f(args...; kwargs...)
    end
end

function pipeline_fixture()
    counters = Dict{Symbol,Int}()
    cache_dir = mktempdir()

    reads = [
        "GATCGTAGCTA",
        "ATCGTAGCTAA",
        "TTGATCGTAGC",
        "GATCGTAGCTG",
    ]
    reference = "TTTTTGATCGTAGCTAGGGTTTGATCGTAGCTA"
    features = [
        (id="geneA", start=6, stop=14),
        (id="geneB", start=15, stop=22),
        (id="geneC", start=23, stop=30),
    ]

    expected_alignments = align_reads(reads, reference; bandwidth=8)
    expected_counts = count_features(expected_alignments, features)
    expected_summary = sort(collect(expected_counts); by=pair -> pair.first)

    graph = PipelineGraph(cache_dir=cache_dir)

    add_node!(graph, pipeline_node(:reads, _counted(counters, :reads, () -> copy(reads)); cache=true))
    add_node!(graph, pipeline_node(:reference, _counted(counters, :reference, () -> reference); cache=true))
    add_node!(graph, pipeline_node(:features, _counted(counters, :features, () -> copy(features)); cache=true))

    add_node!(graph, pipeline_node(
        :alignments,
        _counted(counters, :alignments, (read_batch, ref_seq; bandwidth=8) -> align_reads(read_batch, ref_seq; bandwidth=bandwidth));
        dependencies=[:reads, :reference],
        inputs=Dict{Symbol,Any}(:bandwidth => 8),
        cache=true,
    ))

    add_node!(graph, pipeline_node(
        :feature_counts,
        _counted(counters, :feature_counts, (aln, feature_set) -> count_features(aln, feature_set));
        dependencies=[:alignments, :features],
        cache=true,
    ))

    add_node!(graph, pipeline_node(
        :qc,
        _counted(counters, :qc, read_batch -> length(read_batch));
        dependencies=[:reads],
        cache=true,
    ))

    add_node!(graph, pipeline_node(
        :final,
        _counted(counters, :final, (feature_counts, n_reads) -> (
            summary=sort(collect(feature_counts); by=pair -> pair.first),
            n_reads=n_reads,
            total_feature_counts=sum(values(feature_counts)),
        ));
        dependencies=[:feature_counts, :qc],
        cache=true,
    ))

    return (
        graph=graph,
        counters=counters,
        cached_nodes=[:reads, :reference, :features, :alignments, :feature_counts, :qc, :final],
        expected_levels=[
            [:reads, :reference, :features],
            [:alignments, :qc],
            [:feature_counts],
            [:final],
        ],
        expected_final=(
            summary=expected_summary,
            n_reads=length(reads),
            total_feature_counts=sum(values(expected_counts)),
        ),
    )
end

@inline _bool_score(flag::Bool) = flag ? 1.0 : 0.0

function evaluate_pipeline_recovery(fixture)
    validate_ok = validate_pipeline(fixture.graph)
    levels = execution_levels(fixture.graph)
    level_ok = levels == fixture.expected_levels

    outputs_serial_first = execute_pipeline(fixture.graph; parallel=false)
    counters_after_first = Dict{Symbol,Int}(fixture.counters)

    outputs_serial_second = execute_pipeline(fixture.graph; parallel=false)
    counters_after_second = Dict{Symbol,Int}(fixture.counters)

    outputs_parallel = execute_pipeline(fixture.graph; parallel=true)

    final_expected = fixture.expected_final
    final_serial_ok = outputs_serial_first[:final] == final_expected
    final_repeat_ok = outputs_serial_second[:final] == final_expected
    final_parallel_ok = outputs_parallel[:final] == final_expected

    cache_hits = 0
    for node_id in fixture.cached_nodes
        first_count = get(counters_after_first, node_id, 0)
        second_count = get(counters_after_second, node_id, 0)
        if second_count == first_count
            cache_hits += 1
        end
    end
    cache_hit_rate = cache_hits / length(fixture.cached_nodes)

    rows = [
        (metric="DAG validity", observed=_bool_score(validate_ok), target=1.0, pass=validate_ok),
        (metric="Execution level recovery", observed=_bool_score(level_ok), target=1.0, pass=level_ok),
        (metric="Serial output recovery", observed=_bool_score(final_serial_ok), target=1.0, pass=final_serial_ok),
        (metric="Repeat-run output recovery", observed=_bool_score(final_repeat_ok), target=1.0, pass=final_repeat_ok),
        (metric="Parallel consistency recovery", observed=_bool_score(final_parallel_ok), target=1.0, pass=final_parallel_ok),
        (metric="Cache hit rate", observed=cache_hit_rate, target=1.0, pass=cache_hit_rate >= 0.999),
    ]

    return (
        rows=rows,
        outputs_serial=outputs_serial_first,
        outputs_parallel=outputs_parallel,
        cache_hit_rate=cache_hit_rate,
        levels=levels,
    )
end

end
