using Test
using DataFrames
using DataAPI
using BioToolkit

@testset "Alignment depth APIs" begin
    posterior = BioToolkit.pairhmm_align("ACGT", "ACGT")
    @test posterior isa BioToolkit.PosteriorAlignmentResult
    @test size(posterior.posterior_matrix) == (5, 5, 3)
    @test posterior.consensus_alignment.identity == 1.0
    @test isfinite(posterior.log_likelihood)

    soft_score = BioToolkit.soft_alignment_score("ACGT", "ACGT", BioToolkit.DifferentiableScoring([2.0, -1.0, -2.0]; temperature=0.5))
    @test isfinite(soft_score)
    @test soft_score > 0

    graph = BioToolkit.SequenceGraph([BioToolkit.DNASeq("AC"), BioToolkit.DNASeq("GT")], [(1, 2)])
    graph_alignment = BioToolkit.align_to_graph(BioToolkit.DNASeq("ACGT"), graph)
    @test graph_alignment isa BioToolkit.GraphAlignmentResult
    @test graph_alignment.graph_path == [1, 2]
    @test graph_alignment.query_alignment.identity == 1.0

    profile_a = BioToolkit.AlignmentProfileHMM(UInt8.(collect("ACGT")), Float32[1 0 0 0; 0 1 0 0])
    profile_b = BioToolkit.AlignmentProfileHMM(UInt8.(collect("ACGT")), Float32[1 0 0 0; 0 0 1 0])
    profile_alignment = BioToolkit.align_profiles(profile_a, profile_b)
    @test profile_alignment isa BioToolkit.ProfileAlignmentResult
    @test length(profile_alignment.path) >= 2
end

@testset "Packed genotype storage" begin
    pg = BioToolkit.PackedGenotypes(UInt8[0, 1, 2, 3, 1])
    @test length(pg) == 5
    @test collect(pg) == UInt8[0, 1, 2, 3, 1]
    pg[1] = 2
    @test pg[1] == 2
    @test isnan(BioToolkit.genotype_dosage(pg)[4])
    @test BioToolkit.genotype_missingness(pg) == 0.2

    from_gt = BioToolkit.packed_genotypes_from_gt(["0/0", "0/1", "1/1", "./."])
    @test collect(from_gt) == UInt8[0, 1, 2, 3]
    both = BioToolkit.genotype_and(from_gt, BioToolkit.PackedGenotypes(UInt8[1, 1, 0, 3]))
    @test collect(both) == UInt8[0, 1, 0, 3]
end

@testset "Why provenance primitives" begin
    df = DataFrame(gene=["A", "B", "C"], score=[0.1, 0.9, 0.8])
    filtered = BioToolkit.trace_filter(df, row -> row.score > 0.5; expression="score > 0.5")
    @test nrow(filtered) == 2
    why = BioToolkit.why_in_result(filtered, 1)
    @test why isa BioToolkit.WhyProvenanceNode
    @test why.kept
    @test why.values["gene"] == "B"

    backend = BioToolkit.InMemoryProvenanceBackend()
    node = BioToolkit.ProvenanceNode("n1", "operation", Dict{String,Any}("x" => 1), String[], "now")
    BioToolkit.register_provenance_backend!(backend, node)
    @test haskey(BioToolkit.provenance_backend_nodes(backend), "n1")
end


@testset "Profiled genomic overlap baseline" begin
    subject = BioToolkit.build_collection([
        BioToolkit.GenomicInterval("chr1", 10, 20),
        BioToolkit.GenomicInterval("chr1", 30, 50),
        BioToolkit.GenomicInterval("chr2", 5, 15),
    ])
    profile = BioToolkit.profile_interval_overlaps([
        BioToolkit.GenomicInterval("chr1", 15, 35),
        BioToolkit.GenomicInterval("chr2", 8, 9),
    ], subject; repetitions=2)
    @test profile isa BioToolkit.OverlapProfileResult
    @test profile.indexed_hit_count == 3
    @test profile.naive_hit_count == profile.indexed_hit_count
    @test profile.indexed_seconds >= 0
    @test profile.naive_seconds >= 0
end

@testset "Typed provenance-aware workflow" begin
    cache_dir = mktempdir()
    workflow = BioToolkit.BioWorkflow(cache_dir=cache_dir)
    add_one = BioToolkit.bio_workflow_node(:add_one, x -> x + 1, Tuple{Int}, Int; dependencies=[:input], version="1")
    double = BioToolkit.bio_workflow_node(:double, x -> 2 * x, Tuple{Int}, Int; dependencies=[:add_one], version="1")
    BioToolkit.add_workflow_node!(workflow, add_one)
    BioToolkit.add_workflow_node!(workflow, double)
    @test BioToolkit.validate_workflow(workflow)
    @test BioToolkit.workflow_execution_levels(workflow) == [[:add_one], [:double]]

    first_run = run(workflow, Dict{Symbol,Any}(:input => 4))
    @test first_run.outputs[:double] == 10
    @test first_run.node_status[:add_one] == :computed
    second_run = run(workflow, Dict{Symbol,Any}(:input => 4))
    @test second_run.node_status[:add_one] == :replayed_from_cache
    @test :double in second_run.replayed_nodes

    bad = BioToolkit.bio_workflow_node(:bad, x -> string(x), Tuple{Int}, Int; dependencies=[:input], cache=false)
    invalid_workflow = BioToolkit.BioWorkflow(cache_dir=mktempdir())
    BioToolkit.add_workflow_node!(invalid_workflow, bad)
    @test_throws ArgumentError run(invalid_workflow, Dict{Symbol,Any}(:input => 1))
end


@testset "Provenance migrations and uncertainty propagation" begin
    migrator = BioToolkit.ProvenanceMigrator(current_version=v"1.0")
    BioToolkit.register_migration!(migrator, v"0.9" => v"1.0") do payload
        payload["entity"] = get(payload, "entities", Dict{String,Any}())
        delete!(payload, "entities")
        return payload
    end
    legacy = Dict{String,Any}(
        "schema_version" => "0.9",
        "entities" => Dict{String,Any}("input" => Dict{String,Any}("prov:label" => "input", "prov:value" => Dict{String,Any}())),
        "activity" => Dict{String,Any}(),
        "wasDerivedFrom" => Dict{String,Any}(),
    )
    migrated = BioToolkit.migrate_provenance_payload(migrator, legacy)
    @test migrated["schema_version"] == "1.0"
    @test haskey(migrated, "entity")

    ctx = BioToolkit.ProvenanceContext()
    BioToolkit.register_provenance!(ctx, "input", "measure")
    BioToolkit.register_provenance!(ctx, "model", "fit"; parents=["input"])
    propagated = BioToolkit.propagate_uncertainty!(
        ctx;
        local_variances=Dict("input" => 4.0, "model" => 1.0),
        sensitivities=Dict("model" => Dict("input" => 2.0)))
    @test propagated isa BioToolkit.UncertaintyPropagationResult
    @test propagated.node_variances["input"] == 4.0
    @test propagated.node_variances["model"] == 17.0
    @test ctx.nodes["model"].parameters["uncertainty_variance"] == 17.0
end


@testset "NB-Wald CRISPR screen with RRA ranking" begin
    control = [500.0 520.0 480.0; 510.0 490.0 505.0; 450.0 460.0 455.0; 440.0 435.0 445.0]
    treatment = [45.0 55.0 50.0; 470.0 495.0 480.0; 420.0 430.0 425.0; 435.0 420.0 430.0]
    guide_map = DataFrame(guide=["A1", "B1", "A2", "B2"], gene=["GENE_A", "GENE_B", "GENE_A", "GENE_B"])
    result = BioToolkit.crispr_screen_nb(treatment, control, guide_map; min_mean_count=1.0)
    @test result isa BioToolkit.CRISPRScreenResult
    @test nrow(result.guide_results) == 4
    @test nrow(result.gene_results) == 2
    @test result.global_dispersion >= 1e-8
    @test result.gene_results.gene[1] == "GENE_A"
    @test result.gene_results.direction[1] == "depleted"
    @test BioToolkit.crispr_screen_analysis(treatment, control, guide_map; min_reads=1).gene[1] == "GENE_A"
end


@testset "Deterministic distributed provenance merge" begin
    distributed = BioToolkit.DistributedProvenanceContext([2, 1])
    worker_one = BioToolkit.worker_provenance_context(distributed, 1)
    worker_two = BioToolkit.worker_provenance_context(distributed, 2)
    BioToolkit.register_provenance!(worker_one, "input", "load")
    BioToolkit.register_provenance!(worker_two, "model", "fit"; parents=["input"])
    master = BioToolkit.ProvenanceContext()
    BioToolkit.merge_distributed_provenance!(master, distributed)
    @test length(master.nodes) == 2
    @test master.nodes["model"].parent_ids == ["input"]
    @test [node.id for node in BioToolkit.provenance_ancestors(master, "model")] == ["input"]

    conflicting = BioToolkit.DistributedProvenanceContext([1, 2])
    BioToolkit.register_provenance!(BioToolkit.worker_provenance_context(conflicting, 1), "same", "first")
    BioToolkit.register_provenance!(BioToolkit.worker_provenance_context(conflicting, 2), "same", "second")
    @test_throws ArgumentError BioToolkit.merge_distributed_provenance!(BioToolkit.ProvenanceContext(), conflicting)
end


@testset "FM-index CRISPR off-target enumeration" begin
    guide = BioToolkit.DNASeq("ACGT")
    reference = BioToolkit.DNASeq("ACGTAGGTTTACGAAGG")
    index = BioToolkit.build_offtarget_index(reference; chromosome="chrTest", checkpoint_interval=2, suffix_array_sample_rate=2)
    @test index isa BioToolkit.OffTargetIndex
    indexed = BioToolkit.enumerate_off_targets(guide, index; max_mismatches=1)
    @test nrow(indexed) == 2
    @test indexed.chromosome == ["chrTest", "chrTest"]
    @test indexed.position == [1, 11]
    @test indexed.mismatches == [0, 1]
    @test all(indexed.bulges .== 0)
    diagnostics = DataAPI.metadata(indexed, "off_target_search_diagnostics")
    @test diagnostics isa BioToolkit.OffTargetSearchDiagnostics
    @test diagnostics.verified_sites == 2

    repeated = BioToolkit.DNASeq(repeat("AAAA", 16))
    repeated_index = BioToolkit.build_offtarget_index(repeated; checkpoint_interval=4, suffix_array_sample_rate=2)
    @test_throws ArgumentError BioToolkit.enumerate_off_targets(BioToolkit.DNASeq("AAAA"), repeated_index; max_mismatches=0, max_candidates=1)

    function brute_force_sites(guide_string, contigs, max_mismatches)
        complement = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')
        reverse_complement(sequence) = String(reverse([complement[base] for base in sequence]))
        hits = Tuple{String,Int,Int8,Int}[]
        for (chromosome, sequence) in sort!(collect(contigs); by=first)
            for strand in (Int8(1), Int8(-1))
                working = strand == 1 ? sequence : reverse_complement(sequence)
                for start in 1:(length(working) - length(guide_string) - 3 + 1)
                    site = working[start:start + length(guide_string) - 1]
                    mismatches = count(index -> guide_string[index] != site[index], eachindex(guide_string))
                    pam = working[start + length(guide_string):start + length(guide_string) + 2]
                    mismatches <= max_mismatches && BioToolkit.pam_matches("NGG", pam) || continue
                    position = strand == 1 ? start : length(sequence) - (start + length(guide_string) - 1) + 1
                    push!(hits, (chromosome, position, strand, mismatches))
                end
            end
        end
        return sort!(hits)
    end

    contigs = Dict("chr2" => "TTTACGAAGG", "chr1" => "ACGTAGGACGAAGG")
    multi_index = BioToolkit.build_offtarget_index(contigs; checkpoint_interval=2, suffix_array_sample_rate=2)
    multi_hits = BioToolkit.enumerate_off_targets(guide, multi_index; max_mismatches=1)
    observed = sort!([(row.chromosome, row.position, row.strand, row.mismatches) for row in eachrow(multi_hits)])
    @test observed == brute_force_sites("ACGT", contigs, 1)
end
