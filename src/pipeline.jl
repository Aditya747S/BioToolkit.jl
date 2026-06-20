module Pipeline

using Base.Threads
using DataFrames
using Graphs
using Serialization
using SHA

using ..LongRead: BandedAlignmentResult, NanoporeRead, PacBioRead, banded_semiglobal_alignment
using ..DifferentialExpression: differential_expression
using ..BioToolkit: BioSequence, FastqRecord, read_fasta
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_record, provenance_result!, provenance_summary, register_provenance!

export FileArtifact, PipelineNode, PipelineGraph
export file_artifact, pipeline_node, add_node!, add_dependency!, validate_pipeline, execution_levels, execute_pipeline
export align_reads, count_features
export template_read_fasta_node, template_align_reads_node, template_count_features_node, template_differential_expression_node
export terra_workspace_plan, anvil_workspace_manifest, distributed_pipeline_map
export retry_execute_pipeline, containerized_node
export slurm_array_plan

struct FileArtifact
    path::String
    md5::String
    size::Int
    mtime::Float64
end

struct PipelineNode
    id::Symbol
    func::Function
    dependencies::Vector{Symbol}
    inputs::Dict{Symbol,Any}
    cache::Bool
end

mutable struct PipelineGraph
    graph::SimpleDiGraph
    node_index::Dict{Symbol,Int}
    index_node::Vector{Symbol}
    nodes::Dict{Symbol,PipelineNode}
    cache_dir::String
end

function PipelineGraph(; cache_dir::AbstractString=".biotoolkit_cache", prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    mkpath(cache_dir)
    graph = PipelineGraph(SimpleDiGraph(0), Dict{Symbol,Int}(), Symbol[], Dict{Symbol,PipelineNode}(), String(cache_dir))
    return _register_pipeline_result!(_ctx, graph, "PipelineGraph"; parameters=(cache_dir=String(cache_dir), node_count=0))
end

function _register_pipeline_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, value, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, value, operation; parents=parents, parameters=parameters)
end

function file_artifact(path::AbstractString; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    file_path = String(path)
    if !isfile(file_path)
        artifact = FileArtifact(file_path, "", 0, 0.0)
        return _register_pipeline_result!(_ctx, artifact, "file_artifact"; parameters=(path=file_path, size=artifact.size, mtime=artifact.mtime, md5=artifact.md5, status="missing"))
    end

    data = read(file_path)
    digest = bytes2hex(SHA.sha1(data))
    status = stat(file_path)
    artifact = FileArtifact(file_path, digest, Int(status.size), Float64(status.mtime))
    return _register_pipeline_result!(_ctx, artifact, "file_artifact"; parameters=(path=file_path, size=artifact.size, mtime=artifact.mtime, md5=artifact.md5, status="ok"))
end

function file_artifact_provenance(artifact::FileArtifact)
    status = isempty(artifact.md5) ? :warn : :ok
    notes = isempty(artifact.md5) ? ["file is missing on disk"] : String[]
    return provenance_record(
        "FileArtifact",
        artifact.path;
        status=status,
        notes=notes,
        parameters=(size=artifact.size, mtime=artifact.mtime, md5=artifact.md5))
end

function pipeline_node_provenance(node::PipelineNode)
    return provenance_record(
        string(node.id),
        "pipeline_node";
        status=:ok,
        notes=node.cache ? ["node output may be cached"] : ["node output is uncached"],
        parameters=(dependencies=copy(node.dependencies), input_keys=sort!(collect(string.(keys(node.inputs)))), cache=node.cache))
end

function _register_pipeline_node!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, node::PipelineNode, operation::AbstractString)
    return _register_pipeline_result!(_ctx, node, operation; parameters=(node_id=String(node.id), dependencies=copy(node.dependencies), input_keys=sort!(collect(string.(keys(node.inputs)))), cache=node.cache))
end

function _show_file_artifact(io::IO, artifact::FileArtifact)
    provenance = file_artifact_provenance(artifact)
    print(io, "FileArtifact(", artifact.path, ", size=", artifact.size, ", md5=")
    print(io, isempty(artifact.md5) ? "missing" : artifact.md5)
    print(io, ", mtime=", round(artifact.mtime; digits=3), ") | provenance=", provenance_summary(provenance))
end

function Base.show(io::IO, artifact::FileArtifact)
    _show_file_artifact(io, artifact)
end

function Base.show(io::IO, ::MIME"text/plain", artifact::FileArtifact)
    _show_file_artifact(io, artifact)
end

function Base.show(io::IO, node::PipelineNode)
    print(io, "PipelineNode(", node.id, ", deps=", join(string.(node.dependencies), ", "), ", cache=", node.cache, ")")
end

function Base.show(io::IO, ::MIME"text/plain", node::PipelineNode)
    print(io, "PipelineNode(", node.id, ", deps=", join(string.(node.dependencies), ", "), ", cache=", node.cache, ", provenance=", provenance_summary(pipeline_node_provenance(node)), ")")
end

function Base.show(io::IO, pipeline::PipelineGraph)
    print(io, "PipelineGraph(", nv(pipeline.graph), " nodes, cache_dir=", pipeline.cache_dir, ")")
end

function Base.show(io::IO, ::MIME"text/plain", pipeline::PipelineGraph)
    print(io, "PipelineGraph(", nv(pipeline.graph), " nodes, cache_dir=", pipeline.cache_dir, ")")
end

function pipeline_node(id::Symbol, func::Function; dependencies::AbstractVector{Symbol}=Symbol[], inputs::AbstractDict{Symbol,Any}=Dict{Symbol,Any}(), cache::Bool=true, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    node = PipelineNode(id, func, collect(dependencies), Dict{Symbol,Any}(inputs), cache)
    return _register_pipeline_node!(_ctx, node, "pipeline_node")
end

function add_node!(pipeline::PipelineGraph, node::PipelineNode; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    haskey(pipeline.nodes, node.id) && throw(ArgumentError("node $(node.id) already exists"))

    add_vertex!(pipeline.graph)
    index = nv(pipeline.graph)
    pipeline.node_index[node.id] = index
    push!(pipeline.index_node, node.id)
    pipeline.nodes[node.id] = node

    for dep in node.dependencies
        haskey(pipeline.node_index, dep) || throw(ArgumentError("dependency $(dep) must be added before node $(node.id)"))
        add_edge!(pipeline.graph, pipeline.node_index[dep], index)
    end

    _ctx === nothing || register_provenance!(_ctx, "add_node!"; parents=String[], parameters=(node_id=String(node.id), dependency_count=length(node.dependencies), cache=node.cache))
    return pipeline
end

function add_dependency!(pipeline::PipelineGraph, node_id::Symbol, dependency_id::Symbol; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    haskey(pipeline.nodes, node_id) || throw(ArgumentError("unknown node: $(node_id)"))
    haskey(pipeline.nodes, dependency_id) || throw(ArgumentError("unknown dependency: $(dependency_id)"))

    node = pipeline.nodes[node_id]
    dependency_id in node.dependencies || push!(node.dependencies, dependency_id)
    add_edge!(pipeline.graph, pipeline.node_index[dependency_id], pipeline.node_index[node_id])
    _ctx === nothing || register_provenance!(_ctx, "add_dependency!"; parents=String[], parameters=(node_id=String(node_id), dependency_id=String(dependency_id)))
    return pipeline
end

function _cycle_dfs(graph::SimpleDiGraph)
    state = fill(UInt8(0), nv(graph))

    function visit(vertex::Int)
        if state[vertex] == UInt8(1)
            return true
        elseif state[vertex] == UInt8(2)
            return false
        end

        state[vertex] = UInt8(1)
        for neighbor in outneighbors(graph, vertex)
            visit(neighbor) && return true
        end
        state[vertex] = UInt8(2)
        return false
    end

    for vertex in 1:nv(graph)
        if state[vertex] == UInt8(0) && visit(vertex)
            return true
        end
    end
    return false
end

function validate_pipeline(pipeline::PipelineGraph; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    for node in values(pipeline.nodes)
        for dep in node.dependencies
            haskey(pipeline.nodes, dep) || throw(ArgumentError("node $(node.id) has missing dependency $(dep)"))
        end
    end

    _cycle_dfs(pipeline.graph) && throw(ArgumentError("pipeline graph contains a cycle"))
    _ctx === nothing || register_provenance!(_ctx, "validate_pipeline"; parents=String[], parameters=(node_count=length(pipeline.nodes), edge_count=ne(pipeline.graph)))
    return true
end

function execution_levels(pipeline::PipelineGraph; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    validate_pipeline(pipeline; _ctx=_ctx)

    indegrees = indegree(pipeline.graph)
    ready = sort!(findall(==(0), indegrees))
    levels = Vector{Vector{Symbol}}()

    while !isempty(ready)
        push!(levels, [pipeline.index_node[index] for index in ready])
        next_ready = Int[]
        for node_index in ready
            for neighbor in outneighbors(pipeline.graph, node_index)
                indegrees[neighbor] -= 1
                if indegrees[neighbor] == 0
                    push!(next_ready, neighbor)
                end
            end
        end
        ready = sort!(next_ready)
    end

    _ctx === nothing || register_provenance!(_ctx, "execution_levels"; parents=String[], parameters=(node_count=length(pipeline.nodes), level_count=length(levels)))
    return levels
end

function _state_hash(node::PipelineNode, dependency_values::Vector{Any})
    io = IOBuffer()
    serialize(io, node.id)
    serialize(io, dependency_values)
    sorted_inputs = sort!(collect(node.inputs); by=pair -> string(pair.first))
    serialize(io, sorted_inputs)
    return bytes2hex(SHA.sha1(take!(io)))
end

function _cache_path(pipeline::PipelineGraph, node::PipelineNode, digest::String)
    return joinpath(pipeline.cache_dir, "$(node.id)_$(digest).jls")
end

function _execute_node(pipeline::PipelineGraph, node::PipelineNode, outputs::Dict{Symbol,Any})
    dependency_values = Any[]
    for dep in node.dependencies
        haskey(outputs, dep) || throw(ArgumentError("missing dependency output $(dep) for node $(node.id)"))
        push!(dependency_values, outputs[dep])
    end

    digest = _state_hash(node, dependency_values)
    cache_path = _cache_path(pipeline, node, digest)

    if node.cache && isfile(cache_path)
        cached = open(cache_path, "r") do io
            deserialize(io)
        end
        return cached
    end

    result = node.func(dependency_values...; node.inputs...)
    if node.cache
        open(cache_path, "w") do io
            serialize(io, result)
        end
    end
    return result
end

function execute_pipeline(pipeline::PipelineGraph; context::AbstractDict{Symbol,Any}=Dict{Symbol,Any}(), parallel::Bool=true, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    levels = execution_levels(pipeline; _ctx=_ctx)
    outputs = Dict{Symbol,Any}(context)

    for level in levels
        if parallel && length(level) > 1
            level_results = Vector{Pair{Symbol,Any}}(undef, length(level))
            Threads.@threads for idx in eachindex(level)
                node_id = level[idx]
                node = pipeline.nodes[node_id]
                value = _execute_node(pipeline, node, outputs)
                level_results[idx] = node_id => value
            end
            for pair in level_results
                outputs[pair.first] = pair.second
                _ctx === nothing || register_provenance!(_ctx, "execute_pipeline"; parents=provenance_parent_ids([outputs[dep] for dep in pipeline.nodes[pair.first].dependencies]...), parameters=(node_id=String(pair.first), dependencies=copy(pipeline.nodes[pair.first].dependencies), input_keys=sort!(collect(string.(keys(pipeline.nodes[pair.first].inputs)))), cache=pipeline.nodes[pair.first].cache))
            end
        else
            for node_id in level
                node = pipeline.nodes[node_id]
                outputs[node_id] = _execute_node(pipeline, node, outputs)
                _ctx === nothing || register_provenance!(_ctx, "execute_pipeline"; parents=provenance_parent_ids([outputs[dep] for dep in node.dependencies]...), parameters=(node_id=String(node_id), dependencies=copy(node.dependencies), input_keys=sort!(collect(string.(keys(node.inputs)))), cache=node.cache))
            end
        end
    end

    return outputs
end

@inline function _sequence_string(read)
    if read isa NanoporeRead || read isa PacBioRead
        return String(read.sequence)
    elseif read isa FastqRecord
        return String(read.sequence)
    elseif read isa BioSequence
        return String(read)
    elseif read isa AbstractString
        return String(read)
    end
    throw(ArgumentError("unsupported read type for alignment: $(typeof(read))"))
end

function align_reads(reads::AbstractVector, reference; bandwidth::Integer=64)
    reference_sequence = _sequence_string(reference)
    result = [banded_semiglobal_alignment(_sequence_string(read), reference_sequence; bandwidth=bandwidth) for read in reads]
    _ctx = active_provenance_context()
    return _register_pipeline_result!(_ctx, result, "align_reads"; parents=provenance_parent_ids(reads, reference), parameters=(bandwidth=Int(bandwidth), read_count=length(reads), alignment_count=length(result)))
end

function count_features(alignments::AbstractVector{<:BandedAlignmentResult}, features::AbstractVector)
    counts = Dict{String,Int}()
    for feature in features
        if feature isa NamedTuple
            feature_id = String(feature.id)
            start = Int(feature.start)
            stop = Int(feature.stop)
        elseif feature isa Tuple && length(feature) == 3
            feature_id = String(feature[1])
            start = Int(feature[2])
            stop = Int(feature[3])
        else
            throw(ArgumentError("feature entries must be NamedTuple(id, start, stop) or (id, start, stop) tuples"))
        end

        overlap_count = 0
        for aln in alignments
            if max(start, aln.ref_start) <= min(stop, aln.ref_end)
                overlap_count += 1
            end
        end
        counts[feature_id] = overlap_count
    end
    _ctx = active_provenance_context()
    return _register_pipeline_result!(_ctx, counts, "count_features"; parents=provenance_parent_ids(alignments, features), parameters=(alignment_count=length(alignments), feature_count=length(features), result_count=length(counts)))
end

function template_read_fasta_node(id::Symbol, fasta_path::AbstractString; cache::Bool=true)
    return pipeline_node(id, () -> read_fasta(fasta_path); cache=cache)
end

function template_align_reads_node(id::Symbol; dependencies::AbstractVector{Symbol}=[:reads, :reference], cache::Bool=true, bandwidth::Integer=64)
    inputs = Dict{Symbol,Any}(:bandwidth => bandwidth)
    return pipeline_node(id, (reads, reference; bandwidth=64) -> align_reads(reads, reference; bandwidth=bandwidth); dependencies=dependencies, inputs=inputs, cache=cache)
end

function template_count_features_node(id::Symbol; dependencies::AbstractVector{Symbol}=[:alignments, :features], cache::Bool=true)
    return pipeline_node(id, (alignments, features) -> count_features(alignments, features); dependencies=dependencies, cache=cache)
end

function template_differential_expression_node(id::Symbol; dependencies::AbstractVector{Symbol}=[:count_matrix, :design], cache::Bool=true, kwargs...)
    return pipeline_node(id, (count_matrix, design; kwargs...) -> differential_expression(count_matrix, design; kwargs...); dependencies=dependencies, inputs=Dict{Symbol,Any}(kwargs), cache=cache)
end

"""
    terra_workspace_plan(samples; workspace_name, bucket="")

Build a Terra-style sample entity table for workspace import.
"""
function terra_workspace_plan(samples::DataFrame; workspace_name::AbstractString, bucket::AbstractString="")
    sample_ids = "sample_id" in names(samples) ? String.(samples.sample_id) : ["sample_$(i)" for i in 1:nrow(samples)]
    result = DataFrame(
        workspace=fill(String(workspace_name), nrow(samples)),
        bucket=fill(String(bucket), nrow(samples)),
        sample_id=sample_ids,
        entity_type=fill("sample", nrow(samples)))
    _ctx = active_provenance_context()
    return _register_pipeline_result!(_ctx, result, "terra_workspace_plan"; parents=provenance_parent_ids(samples), parameters=(workspace_name=String(workspace_name), bucket=String(bucket), row_count=nrow(result)))
end

"""
    anvil_workspace_manifest(file_paths; entity_type="sample")

Create an AnVIL-ready URI manifest for entity ingestion.
"""
function anvil_workspace_manifest(file_paths::AbstractVector{<:AbstractString}; entity_type::AbstractString="sample")
    result = DataFrame(entity_type=fill(String(entity_type), length(file_paths)), uri=String.(file_paths), basename=basename.(String.(file_paths)))
    _ctx = active_provenance_context()
    return _register_pipeline_result!(_ctx, result, "anvil_workspace_manifest"; parents=provenance_parent_ids(file_paths), parameters=(entity_type=String(entity_type), row_count=nrow(result)))
end

"""
    distributed_pipeline_map(f, xs; nworkers=Threads.nthreads())

Threaded async map helper for embarrassingly parallel pipeline stages.
"""
function distributed_pipeline_map(f::Function, xs; nworkers::Int=Threads.nthreads())
    items = collect(xs)
    result = collect(asyncmap(f, items; ntasks=max(1, Int(nworkers))))
    _ctx = active_provenance_context()
    return _register_pipeline_result!(_ctx, result, "distributed_pipeline_map"; parents=provenance_parent_ids(items), parameters=(worker_count=Int(nworkers), input_count=length(items), output_count=length(result)))
end

"""
    containerized_node(id, func; image, command="")

Attach container execution metadata to a pipeline node specification.
"""
function containerized_node(id::Symbol, func::Function; image::AbstractString, command::AbstractString="", dependencies::AbstractVector{Symbol}=Symbol[], inputs::AbstractDict{Symbol,Any}=Dict{Symbol,Any}(), cache::Bool=true)
    merged_inputs = Dict{Symbol,Any}(inputs)
    merged_inputs[:container_image] = String(image)
    isempty(command) || (merged_inputs[:container_command] = String(command))
    return pipeline_node(id, func; dependencies=dependencies, inputs=merged_inputs, cache=cache)
end

"""
    retry_execute_pipeline(pipeline; max_retries=2)

Execute pipeline with bounded retries, returning attempts and status.
"""
function retry_execute_pipeline(pipeline::PipelineGraph; max_retries::Int=2, kwargs...)
    _ctx = active_provenance_context()
    max_retries >= 0 || throw(ArgumentError("max_retries must be >= 0"))
    attempts = 0
    last_error = nothing
    while attempts <= max_retries
        attempts += 1
        try
            result = execute_pipeline(pipeline; kwargs...)
            outcome = (status=:ok, attempts=attempts, result=result)
            return _register_pipeline_result!(_ctx, outcome, "retry_execute_pipeline"; parameters=(max_retries=max_retries, attempts=attempts, status="ok"))
        catch err
            last_error = err
            attempts > max_retries && break
        end
    end
    outcome = (status=:failed, attempts=attempts, error=last_error)
    error_type = last_error === nothing ? "" : String(typeof(last_error))
    error_text = last_error === nothing ? "" : sprint(showerror, last_error)
    return _register_pipeline_result!(_ctx, outcome, "retry_execute_pipeline"; parameters=(max_retries=max_retries, attempts=attempts, status="failed", error_type=error_type, error_message=error_text))
end

"""
    slurm_array_plan(samples; job_name="biotoolkit", cpus_per_task=8)

Create a SLURM job-array plan table for scalable HPC execution.
"""
function slurm_array_plan(samples::DataFrame; job_name::AbstractString="biotoolkit", cpus_per_task::Int=8, gpus_per_task::Int=0, mem_gb::Int=32, time_limit::AbstractString="12:00:00", prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    nrow(samples) > 0 || throw(ArgumentError("samples must contain at least one row"))
    cpus_per_task >= 1 || throw(ArgumentError("cpus_per_task must be >= 1"))
    gpus_per_task >= 0 || throw(ArgumentError("gpus_per_task must be >= 0"))
    mem_gb >= 1 || throw(ArgumentError("mem_gb must be >= 1"))

    sid = hasproperty(samples, :sample_id) ? String.(samples.sample_id) : ["sample_$(i)" for i in 1:nrow(samples)]
    arr = 1:nrow(samples)
    result = DataFrame(
        array_task_id=collect(arr),
        sample_id=sid,
        job_name=fill(String(job_name), nrow(samples)),
        cpus_per_task=fill(Int(cpus_per_task), nrow(samples)),
        gpus_per_task=fill(Int(gpus_per_task), nrow(samples)),
        mem_gb=fill(Int(mem_gb), nrow(samples)),
        time_limit=fill(String(time_limit), nrow(samples)),
        slurm_directive=["sbatch --array=$(first(arr))-$(last(arr)) --cpus-per-task=$(cpus_per_task) --mem=$(mem_gb)G --time=$(time_limit)" for _ in arr])
    return _register_pipeline_result!(_ctx, result, "slurm_array_plan"; parents=provenance_parent_ids(samples), parameters=(job_name=String(job_name), cpus_per_task=Int(cpus_per_task), gpus_per_task=Int(gpus_per_task), mem_gb=Int(mem_gb), time_limit=String(time_limit), row_count=nrow(result)))
end

end
