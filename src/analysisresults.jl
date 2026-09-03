using DataFrames
using DataAPI
using Dates
using JSON
using DataStructures: OrderedDict
using Pkg
using Random
using SHA
using UUIDs

# ==============================================================================
# analysisresults.jl usage guide
# ==============================================================================
# This file provides two related facilities:
#
# 1. A common analysis-result protocol
#    - Subtype `AbstractAnalysisResult` for domain result structs.
#    - Store a `provenance::ResultProvenance` field when the result should carry
#      compact embedded provenance.
#    - Override `analysis_result_type(::Type{YourResult})` or
#      `analysis_result_fields(::Type{YourResult})` only when the default summary
#      and one-row `DataFrame(result)` view are not sufficient.
#    - Use `analysis_result_summary(result)`,
#      `analysis_result_provenance(result)`, `summary(result)`, and
#      `DataFrame(result)` for uniform display/export behavior.
#
# 2. Data provenance tracking
#    There are three supported user-facing ways to record provenance:
#
#    A. Global interactive tracking:
#       ```julia
#       ctx = enable_provenance!()
#       result = some_biotoolkit_function(data)
#       export_provenance_json(ctx)
#       disable_provenance!()
#       ```
#       This is convenient in notebooks and scripts. Avoid sharing a plain
#       `ProvenanceContext` across threaded workloads; use
#       `ThreadSafeProvenanceContext()` for concurrent writers.
#
#    B. Scoped advanced-user tracking:
#       ```julia
#       ctx = ProvenanceContext()
#       result = with_provenance(ctx) do
#           step1 = some_biotoolkit_function(data)
#           some_other_function(step1)
#       end
#       ```
#       This is the preferred per-analysis/per-experiment API. It routes all
#       nested BioToolkit calls that use `active_provenance_context()` into `ctx`
#       without forcing every public function to expose a `prov_ctx=` keyword.
#       The scope is task-local and restores the previous state automatically.
#
#    C. Explicit call-site tracking:
#       ```julia
#       ctx = ProvenanceContext()
#       result = function_that_accepts_it(data; prov_ctx=ctx)
#       ```
#       Use this where a public function already exposes `prov_ctx=`, or inside
#       library code that needs to override any ambient scoped/global context.
#
#    For BioToolkit implementers:
#    - Resolve context once near the public boundary with
#      `_ctx = active_provenance_context(prov_ctx)` when a function has
#      `prov_ctx=`, or `_ctx = active_provenance_context()` otherwise.
#    - Pass `_ctx` through internal helper calls when they record provenance.
#    - Finish by calling `provenance_result!(_ctx, result, "operation";
#      parents=provenance_parent_ids(inputs...), parameters=(...))`.
#    - Do not manually inspect global refs or task-local storage from domain
#      modules; `active_provenance_context` is the only resolver they should use.
#
# Public provenance utilities:
# - `ProvenanceContext(; max_nodes, max_depth)` stores a DAG of operations.
# - `ThreadSafeProvenanceContext()` wraps a context with a lock for threaded use.
# - `register_provenance!` manually adds DAG nodes.
# - `provenance_result!` stamps supported result containers and registers a DAG
#   node; `provenance_result!(nothing, ...)` is a no-op.
# - `with_provenance(result, label, source; ...)` stamps metadata on an existing
#   result/table and, when a scoped/global context is active, also registers it.
# - `stamp_provenance!`, `update_provenance!`, `metadata_provenance`,
#   `has_provenance`, and `requires_provenance` manage/query embedded metadata.
# - `container_provenance_id`, `container_provenance_hash`,
#   `container_provenance_summary`, and `provenance_parent_ids` connect result
#   containers into parent/child lineages.
# - `export_provenance_json` and `import_provenance_json` serialize/restore the
#   PROV-JSON representation; `provenance_lineage_table` gives a tabular view.
# - `find_nodes`, `provenance_chain`, `provenance_ancestors`,
#   `provenance_descendants`, `provenance_diff`, `merge_provenance_contexts`,
#   `provenance_to_dot`, and `provenance_to_mermaid` inspect or visualize DAGs.
# - `new_provenance_id()` creates a random SHA-256-style identifier;
#   `new_provenance_id(content)` creates a deterministic content-addressed ID.
# - `provenance_structural_hash`, `provenance_content_hash`, and
#   `file_provenance_hash` provide cheap shape hashes, full in-memory content
#   hashes, and streaming file hashes respectively.
# - `stamp_provenance_with_hash!` records compact hash references without copying
#   large matrices/sequences into provenance payloads.
# - `@provenance` is for explicit low-level scoped registration around an
#   expression; `@pure_provenance` is for deterministic pure-function nodes.
# - `generate_methods_section(ctx)` converts the DAG into methods-style prose for
#   manuscripts/reports. Review the prose before publication.
#
# Scientific-software cautions:
# - Provenance must never change numerical results. Keep instrumentation outside
#   the algorithmic core and pass `_ctx` explicitly through helper boundaries.
# - Prefer small, serializable parameter summaries over storing full datasets.
# - Use content/file hashes for integrity claims; structural hashes prove only
#   shape/type stability.
# - Treat `max_nodes`/`max_depth` as audit safeguards. Limit violations throw
#   rather than silently dropping lineage.
# ==============================================================================


# ==============================================================================
# Global (opt-in) provenance context
# ==============================================================================
# Beginners: just call `enable_provenance!()` once before running their pipeline.
# Every BioToolkit function will then record itself automatically — no prov_ctx
# argument ever needed.
#
# Power users / library authors: pass prov_ctx= explicitly to route into a
# specific context (e.g. per-thread or per-experiment), overriding the global.
# ==============================================================================

# NOTE: ProvenanceContext is defined later in this file; these containers are
# untyped to avoid a forward-reference error. Type safety is enforced at the
# public API boundary.
const _GLOBAL_PROV_CTX = Ref{Any}(nothing)
const _GLOBAL_PROV_LOCK = ReentrantLock()
const _PROV_ENABLED = Ref{Bool}(false)  # guarded by _GLOBAL_PROV_LOCK for all writes
const _SCOPED_PROV_STACK_KEY = :_biotoolkit_scoped_prov_stack

@inline function _is_valid_provenance_context(ctx)
  return ctx === nothing || ctx isa ProvenanceContext || ctx isa ThreadSafeProvenanceContext
end

function _check_provenance_context(ctx, name::AbstractString)
  _is_valid_provenance_context(ctx) ||
    throw(ArgumentError("$name must be a ProvenanceContext, ThreadSafeProvenanceContext, or nothing"))
  return ctx
end

function _scoped_provenance_stack()
  stack = get(task_local_storage(), _SCOPED_PROV_STACK_KEY, nothing)
  if stack === nothing
    stack = Any[]
    task_local_storage(_SCOPED_PROV_STACK_KEY, stack)
  end
  return stack::Vector{Any}
end

@inline function _active_scoped_provenance_context()
  stack = get(task_local_storage(), _SCOPED_PROV_STACK_KEY, nothing)
  stack === nothing && return nothing
  isempty(stack) && return nothing
  return stack[end]
end

"""
    enable_provenance!([ctx]) → ProvenanceContext

Activate **global automatic provenance tracking**. After this call every
BioToolkit function records itself in `ctx` without the caller writing a single
`prov_ctx=` argument.

```julia
enable_provenance!()          # create a fresh context automatically
result = pairwise_align(a, b) # recorded silently
println(generate_methods_section(get_provenance_context()))
```

Pass an existing `ProvenanceContext` to resume a session or merge with prior work:

```julia
ctx = ProvenanceContext()
enable_provenance!(ctx)
```

!!! note "Thread safety"
    The global `Ref` storing the context is read atomically, but the
    `ProvenanceContext` it points to is **not** thread-safe for concurrent
    reads and writes. Concurrent calls from multiple threads that all write
    to the same shared context will race on its internal `OrderedDict`.

    For any multi-threaded workload use `ThreadSafeProvenanceContext`, which
    wraps every mutation in a `ReentrantLock`. The shared global is intentionally
    optimised for the common single-threaded interactive case.
"""
function enable_provenance!()
  # ProvenanceContext is defined later in this file; we invoke it by name
  # once the module is fully loaded — this function is only called at runtime.
  Base.lock(_GLOBAL_PROV_LOCK) do
    _GLOBAL_PROV_CTX[] = ProvenanceContext()
    _PROV_ENABLED[] = true  # safe: called inside _GLOBAL_PROV_LOCK
    return _GLOBAL_PROV_CTX[]
  end
end
function enable_provenance!(ctx)
  _check_provenance_context(ctx, "enable_provenance! ctx")
  Base.lock(_GLOBAL_PROV_LOCK) do
    _GLOBAL_PROV_CTX[] = ctx
    _PROV_ENABLED[] = ctx !== nothing
    return ctx
  end
end
export enable_provenance!

"""
    disable_provenance!()

Deactivate global provenance tracking. Functions revert to zero-cost no-op mode.
"""
function disable_provenance!()
  Base.lock(_GLOBAL_PROV_LOCK) do
    _GLOBAL_PROV_CTX[] = nothing
    _PROV_ENABLED[] = false
  end
  return nothing
end
export disable_provenance!

"""
    get_provenance_context() → Union{Nothing, ProvenanceContext}

Return the active global `ProvenanceContext`, or `nothing` if provenance is
disabled. Use this to inspect or export the recorded lineage after a pipeline:

```julia
enable_provenance!()
# ... run analysis ...
ctx = get_provenance_context()
println(generate_methods_section(ctx))
json = export_provenance_json(ctx)
```
"""
function get_provenance_context()
  Base.lock(_GLOBAL_PROV_LOCK) do
    return _GLOBAL_PROV_CTX[]
  end
end
export get_provenance_context

"""
    with_global_provenance(f, [ctx]) → result

Run `f()` with `ctx` (or a fresh context) as the global provenance context,
then restore the previous context. Thread-safe for scoped per-experiment use.

```julia
result = with_global_provenance() do
    normalize(data)   # recorded in a temporary context
end
```
"""
function with_global_provenance(f)
  previous = Base.lock(_GLOBAL_PROV_LOCK) do
    prev = _GLOBAL_PROV_CTX[]
    _GLOBAL_PROV_CTX[] = ProvenanceContext()
    _PROV_ENABLED[] = true
    prev
  end
  try
    return f()
  finally
    Base.lock(_GLOBAL_PROV_LOCK) do
      _GLOBAL_PROV_CTX[] = previous
      _PROV_ENABLED[] = previous !== nothing
    end
  end
end
function with_global_provenance(f, ctx)
  previous = Base.lock(_GLOBAL_PROV_LOCK) do
    prev = _GLOBAL_PROV_CTX[]
    _GLOBAL_PROV_CTX[] = ctx
    _PROV_ENABLED[] = true
    prev
  end
  try
    return f()
  finally
    Base.lock(_GLOBAL_PROV_LOCK) do
      _GLOBAL_PROV_CTX[] = previous
      _PROV_ENABLED[] = previous !== nothing
    end
  end
end
export with_global_provenance

"""
    active_provenance_context([explicit])

**Internal helper used by BioToolkit functions.**

Returns `explicit` if it is not `nothing`, otherwise falls back to a task-local
`with_provenance(ctx) do ... end` scope and then to the global context set by
`enable_provenance!()`. This is the single call every public function makes to
resolve which context (if any) to record into.
"""
@inline function active_provenance_context(explicit=nothing)
  explicit !== nothing && return explicit
  scoped = _active_scoped_provenance_context()
  scoped !== nothing && return scoped
  _PROV_ENABLED[] || return nothing
  Base.lock(_GLOBAL_PROV_LOCK) do
    return _PROV_ENABLED[] ? _GLOBAL_PROV_CTX[] : nothing
  end
end

"""
    AbstractAnalysisResult

Common supertype for BioToolkit analysis outputs.

Concrete result types can use this to participate in a shared display and
tabular-conversion protocol. The default `summary(result)` method prints a
compact field preview, `show(result)` uses that summary, and
`DataFrame(result)` produces a one-row table when a domain-specific adapter is
not more specific.
"""
abstract type AbstractAnalysisResult end

"""
    ResultProvenance

Compact provenance record for a BioToolkit analysis result.

The fields are intentionally lightweight so they can be surfaced in compact
`show` output and table exports without requiring every domain object to store
its own provenance payload.
"""
struct ResultProvenance
  id::String
  label::String
  source::String
  status::Symbol
  warnings::Vector{String}
  errors::Vector{String}
  fallbacks::Vector{String}
  notes::Vector{String}
  parameters::NamedTuple
  parent_ids::Vector{String}
  timestamp::String
end

const PROVENANCE_ID_KEY = :provenance_id
const PROVENANCE_HASH_KEY = :provenance_hash
const PROVENANCE_SCHEMA_VERSION = "1.0"

@inline _metadata_key(metadata::AbstractDict, key::Symbol) = keytype(metadata) <: Symbol ? key : String(key)
@inline _metadata_key(metadata::AbstractDict, key::AbstractString) = keytype(metadata) <: Symbol ? Symbol(key) : String(key)

struct ProvenanceNode
  id::String
  operation::String
  parameters::Dict{String,Any}
  parent_ids::Vector{String}
  timestamp::String
end

struct EnvironmentSnapshot
  julia_version::String
  manifest_hash::Union{String,Nothing}
  package_versions::Dict{String,String}
  cpu_threads::Int
  git_commit::Union{String,Nothing}
  random_state_hash::String
  timestamp::String
end

mutable struct ProvenanceContext
  nodes::OrderedDict{String,ProvenanceNode}
  environment::EnvironmentSnapshot
  max_nodes::Union{Int,Nothing}
  max_depth::Union{Int,Nothing}
  schema_version::String
  placeholder_ids::Set{String}
  node_depths::Dict{String,Int}
end

const _GIT_COMMIT_CACHE = Ref{Tuple{Union{Nothing,String},Union{Nothing,String}}}((nothing, nothing))

function _git_commit()
  proj = Base.active_project()
  # When there's no active project, fall back to pwd() as cache key so that
  # different directories don't share a cached git commit hash.
  cache_key = proj === nothing ? pwd() : proj
  cached_key, cached_commit = _GIT_COMMIT_CACHE[]
  cached_key == cache_key && return cached_commit
  # Walk up from active project to find .git
  dir = proj === nothing ? pwd() : dirname(proj)
  while dir != dirname(dir)  # stop at root
    git_path = joinpath(dir, ".git")
    if isdir(git_path) || isfile(git_path)
      try
        commit = readchomp(`git -C $dir rev-parse HEAD`)
        result = isempty(commit) ? nothing : commit
        _GIT_COMMIT_CACHE[] = (cache_key, result)
        return result
      catch
        _GIT_COMMIT_CACHE[] = (cache_key, nothing)
        return nothing
      end
    end
    dir = dirname(dir)
  end
  _GIT_COMMIT_CACHE[] = (cache_key, nothing)
  return nothing
end

function _manifest_hash()
  path = Base.active_project()
  path === nothing && return nothing
  manifest_path = joinpath(dirname(path), "Manifest.toml")
  isfile(manifest_path) || return nothing
  return bytes2hex(sha256(read(manifest_path)))
end

function _package_versions()
  versions = Dict{String,String}()
  try
    for (_, dep) in Pkg.dependencies()
      # dep.name is always a String; dep.version can be nothing (stdlib, path deps)
      dep.version === nothing && continue
      versions[dep.name] = string(dep.version)
    end
  catch
    return Dict{String,String}()
  end
  return versions
end

function capture_environment_snapshot()
  rng_hash = try
    bytes2hex(sha256(codeunits(repr(copy(Random.default_rng())))))
  catch
    bytes2hex(sha256(codeunits("unavailable")))
  end
  return EnvironmentSnapshot(
    string(VERSION),
    _manifest_hash(),
    _package_versions(),
    Threads.nthreads(),
    _git_commit(),
    rng_hash,
    _provenance_timestamp(),
  )
end

ProvenanceContext(; max_nodes::Union{Int,Nothing}=nothing, max_depth::Union{Int,Nothing}=nothing, schema_version::AbstractString=PROVENANCE_SCHEMA_VERSION) =
  ProvenanceContext(OrderedDict{String,ProvenanceNode}(), capture_environment_snapshot(), max_nodes, max_depth, String(schema_version), Set{String}(), Dict{String,Int}())
ProvenanceContext(nodes::OrderedDict{String,ProvenanceNode}) =
  ProvenanceContext(nodes, capture_environment_snapshot(), nothing, nothing, PROVENANCE_SCHEMA_VERSION, Set{String}(), Dict{String,Int}())
ProvenanceContext(nodes::OrderedDict{String,ProvenanceNode}, environment::EnvironmentSnapshot, max_nodes::Union{Int,Nothing}, max_depth::Union{Int,Nothing}) =
  ProvenanceContext(nodes, environment, max_nodes, max_depth, PROVENANCE_SCHEMA_VERSION, Set{String}(), Dict{String,Int}())
ProvenanceContext(nodes::OrderedDict{String,ProvenanceNode}, environment::EnvironmentSnapshot, max_nodes::Union{Int,Nothing}, max_depth::Union{Int,Nothing}, schema_version::String) =
  ProvenanceContext(nodes, environment, max_nodes, max_depth, schema_version, Set{String}(), Dict{String,Int}())

"""
    ThreadSafeProvenanceContext

A thread-safe wrapper around `ProvenanceContext` that serialises all
node registrations through a `ReentrantLock`.  Use this whenever multiple
Julia tasks or `@threads` loops share a single provenance context.

!!! note
    `ProvenanceContext` is **not** thread-safe by itself; concurrent writes to
    its `OrderedDict` will race. Use `ThreadSafeProvenanceContext` for any
    multi-threaded workload.
"""
mutable struct ThreadSafeProvenanceContext
  ctx::ProvenanceContext
  lock::ReentrantLock
end

ThreadSafeProvenanceContext() = ThreadSafeProvenanceContext(ProvenanceContext(), ReentrantLock())

# Delegate read-only access to the inner context
function Base.getproperty(ts::ThreadSafeProvenanceContext, s::Symbol)
  s === :lock && return getfield(ts, :lock)
  s === :ctx && return getfield(ts, :ctx)
  return Base.lock(getfield(ts, :lock)) do
    val = getproperty(getfield(ts, :ctx), s)
    s === :nodes ? copy(val) : val
  end
end
Base.propertynames(ts::ThreadSafeProvenanceContext, private::Bool=false) =
  (:ctx, :lock, propertynames(getfield(ts, :ctx), private)...)

Base.summary(ts::ThreadSafeProvenanceContext) = "ThreadSafe" * summary(ts.ctx)
Base.show(io::IO, ts::ThreadSafeProvenanceContext) = print(io, summary(ts))

# Thread-safe registration — acquires the lock before writing
function register_provenance!(ts::ThreadSafeProvenanceContext, node::ProvenanceNode)
  Base.lock(ts.lock) do
    register_provenance!(ts.ctx, node)
  end
  return node
end
function register_provenance!(ts::ThreadSafeProvenanceContext, operation::AbstractString; kwargs...)
  node = ProvenanceNode(
    new_provenance_id(),
    String(operation),
    _provenance_parameters_dict(get(kwargs, :parameters, NamedTuple())),
    String.(get(kwargs, :parents, String[])),
    _provenance_timestamp(),
  )
  Base.lock(ts.lock) do
    register_provenance!(ts.ctx, node)
  end
  return node
end
function register_provenance!(ts::ThreadSafeProvenanceContext, node_id::AbstractString, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
  node = ProvenanceNode(String(node_id), String(operation), _provenance_parameters_dict(parameters), String.(parents), _provenance_timestamp())
  Base.lock(ts.lock) do
    register_provenance!(ts.ctx, node)
  end
  return node
end

# Keep the whole stamp-and-register operation under one lock. Calling the inner
# context directly here is safe because the lock is reentrant for nested helpers.
function provenance_result!(ts::ThreadSafeProvenanceContext, result, operation::AbstractString; kwargs...)
  Base.lock(ts.lock) do
    return provenance_result!(ts.ctx, result, operation; kwargs...)
  end
end

export ThreadSafeProvenanceContext

@inline function _provenance_timestamp()
  return Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sss\Z")
end

"""    new_provenance_id() → String

Generate a unique, non-deterministic provenance ID using time, randomness,
and thread ID. Suitable for tracking operations whose output depends on
non-deterministic state (I/O, random initialization, etc.).
"""
function new_provenance_id()
  return bytes2hex(sha256(codeunits(string(UUIDs.uuid4(), ':', time_ns(), ':', rand(UInt)))))
end

"""
    new_provenance_id(content) → String

Generate a **deterministic** provenance ID by SHA-256 hashing `content`.
Two calls with identical content always return the same ID — enabling
content-addressable, reproducible lineage graphs for pure functions.

Use this whenever the output is a deterministic function of the inputs:

```julia
id = new_provenance_id(string(method, ':', size(matrix)))
```
"""
function new_provenance_id(content::AbstractString)
  return bytes2hex(sha256(codeunits(String(content))))
end
function new_provenance_id(content::AbstractVector{UInt8})
  return bytes2hex(sha256(content))
end

const _PROV_MAX_ARRAY_LEN = Ref{Int}(64)
const _PROV_MAX_STRING_LEN = Ref{Int}(256)

function set_provenance_limits!(; max_array_len::Int, max_string_len::Int)
  _PROV_MAX_ARRAY_LEN[] = max_array_len
  _PROV_MAX_STRING_LEN[] = max_string_len
end
export set_provenance_limits!

function _provenance_array_summary_hash(value::AbstractArray)
  ctx_sha = SHA.SHA2_256_CTX()
  SHA.update!(ctx_sha, codeunits(string(size(value), ':', eltype(value), ':', length(value))))
  # Restrict pointer fast path to concrete Array to avoid MethodError on SubArray etc.
  if isbitstype(eltype(value)) && value isa Array && IndexStyle(value) === IndexLinear()
    GC.@preserve value begin
      p = convert(Ptr{UInt8}, pointer(value))
      sz = sizeof(value)
      if p != C_NULL && sz > 0
        buf = unsafe_wrap(Array, p, sz)
        SHA.update!(ctx_sha, buf)
      end
    end
  else
    for item in value
      SHA.update!(ctx_sha, codeunits(repr(item)))
      SHA.update!(ctx_sha, UInt8[0])
    end
  end
  return bytes2hex(SHA.digest!(ctx_sha))
end

function _provenance_json_value(value)
  if value isa AbstractDict
    return OrderedDict(string(key) => _provenance_json_value(item) for (key, item) in value)
  elseif value isa AbstractArray && !(value isa AbstractString)
    n = length(value)
    if n > _PROV_MAX_ARRAY_LEN[]
      # Store a compact summary instead of the full array to avoid RAM doubling
      return OrderedDict(
        "__truncated__" => true,
        "length" => n,
        "eltype" => string(eltype(value)),
        "hash" => _provenance_array_summary_hash(value)
      )
    end
    return [_provenance_json_value(item) for item in value]
  elseif value isa Tuple
    return [_provenance_json_value(item) for item in value]
  elseif value isa NamedTuple
    return OrderedDict(string(key) => _provenance_json_value(getproperty(value, key)) for key in keys(value))
  elseif value isa Symbol
    return string(value)
  elseif value isa AbstractString
    s = String(value)
    nb = ncodeunits(s)
    if nb > _PROV_MAX_STRING_LEN[]
      # Truncate by codeunit count (bytes), cutting on a valid UTF-8 boundary
      limit = _PROV_MAX_STRING_LEN[]
      cut = min(limit, nb)
      # Step back from the cut point to a valid character boundary
      while cut > 0 && cut < nb && !isvalid(s, cut + 1)
        cut -= 1
      end
      trunc = s[1:min(cut, nb)]
      return string(trunc, "…[hash:", bytes2hex(sha256(codeunits(s)))[1:16], "]")
    end
    return s
  elseif value === nothing || value isa Number || value isa Bool
    return value
  else
    return string(value)
  end
end

@inline function _provenance_parameters_dict(parameters)
  if parameters isa NamedTuple
    return Dict{String,Any}(string(key) => _provenance_json_value(getproperty(parameters, key)) for key in keys(parameters))
  elseif parameters isa AbstractDict
    return Dict{String,Any}(string(key) => _provenance_json_value(value) for (key, value) in parameters)
  elseif parameters isa Base.Pairs
    return Dict{String,Any}(string(key) => _provenance_json_value(value) for (key, value) in parameters)
  elseif parameters isa AbstractVector{<:Pair}
    return Dict{String,Any}(string(key) => _provenance_json_value(value) for (key, value) in parameters)
  else
    return Dict{String,Any}()
  end
end

function _can_store_provenance_metadata(metadata::AbstractDict)
  key_type = keytype(metadata)
  value_type = valtype(metadata)
  can_key = key_type === Any || String <: key_type || Symbol <: key_type
  can_value = value_type === Any || ResultProvenance <: value_type || AbstractString <: value_type
  return can_key && can_value
end

function container_provenance_dict(container)
  container isa AbstractDict && return _can_store_provenance_metadata(container) ? container : nothing
  if hasproperty(container, :metadata)
    metadata = getproperty(container, :metadata)
    metadata isa AbstractDict && return metadata
  end
  if hasproperty(container, :annotations)
    annotations = getproperty(container, :annotations)
    annotations isa AbstractDict && return annotations
  end
  return nothing
end

function container_provenance_id(result::AbstractAnalysisResult)
  if hasproperty(result, :provenance)
    return result.provenance.id
  end
  return nothing
end

function container_provenance_id(container)
  metadata = container_provenance_dict(container)
  metadata === nothing && return nothing
  if haskey(metadata, PROVENANCE_ID_KEY)
    value = metadata[PROVENANCE_ID_KEY]
    return value === nothing ? nothing : String(value)
  elseif haskey(metadata, String(PROVENANCE_ID_KEY))
    value = metadata[String(PROVENANCE_ID_KEY)]
    return value === nothing ? nothing : String(value)
  end
  return nothing
end

function container_provenance_id(table::DataFrames.AbstractDataFrame)
  for key in (String(PROVENANCE_ID_KEY), string(PROVENANCE_ID_KEY))
    value = try
      DataAPI.metadata(table, key)
    catch
      nothing
    end
    value === nothing || return String(value)
  end
  return nothing
end

function ensure_provenance_id!(metadata::AbstractDict)
  if haskey(metadata, PROVENANCE_ID_KEY)
    metadata[PROVENANCE_ID_KEY] = String(metadata[PROVENANCE_ID_KEY])
    return String(metadata[PROVENANCE_ID_KEY])
  elseif haskey(metadata, String(PROVENANCE_ID_KEY))
    identifier = String(metadata[String(PROVENANCE_ID_KEY)])
    if keytype(metadata) <: Symbol
      metadata[PROVENANCE_ID_KEY] = identifier
      delete!(metadata, String(PROVENANCE_ID_KEY))
    else
      metadata[String(PROVENANCE_ID_KEY)] = identifier
    end
    return identifier
  else
    identifier = new_provenance_id()
    metadata[_metadata_key(metadata, PROVENANCE_ID_KEY)] = identifier
    return identifier
  end
end

function ensure_provenance_id!(table::DataFrames.AbstractDataFrame)
  identifier = container_provenance_id(table)
  if identifier === nothing
    identifier = new_provenance_id()
    DataAPI.metadata!(table, String(PROVENANCE_ID_KEY), identifier; style=:note)
  end
  return identifier
end

function _append_provenance_ids!(ids::Vector{String}, value)
  if value isa AbstractVector
    for item in value
      _append_provenance_ids!(ids, item)
    end
    return ids
  elseif value isa Tuple
    for item in value
      _append_provenance_ids!(ids, item)
    end
    return ids
  elseif value isa AbstractSet
    for item in value
      _append_provenance_ids!(ids, item)
    end
    return ids
  else
    identifier = container_provenance_id(value)
    identifier === nothing || push!(ids, identifier)
    return ids
  end
end

function provenance_parent_ids(values...)
  ids = String[]
  for value in values
    _append_provenance_ids!(ids, value)
  end
  return unique(ids)
end

function _provenance_node_depth(ctx::ProvenanceContext, node::ProvenanceNode)
  isempty(node.parent_ids) && return 1
  d = 1
  for pid in node.parent_ids
    pd = get(ctx.node_depths, pid, nothing)
    if pd === nothing
      parent = get(ctx.nodes, pid, nothing)
      pd = parent === nothing ? 1 : _provenance_node_depth(ctx, parent)
    end
    d = max(d, 1 + pd)
  end
  return d
end

# Use a dedicated sentinel key to mark auto-registered placeholder nodes.
# This avoids the fragile heuristic of checking parameter names.
const _PROV_PLACEHOLDER_KEY = "__placeholder__"

function _is_generic_self_registration(node::ProvenanceNode)
  get(node.parameters, _PROV_PLACEHOLDER_KEY, false) === true && return true
  isempty(node.parameters) && return true
  return length(node.parameters) == 1 &&
         get(node.parameters, "operation", nothing) == node.operation
end

function _prune_generic_self_registrations!(ctx::ProvenanceContext, node::ProvenanceNode)
  _is_generic_self_registration(node) && return nothing
  isempty(ctx.placeholder_ids) && return nothing

  stale_ids = String[]
  replacement = Dict{String,String}()
  for id in ctx.placeholder_ids
    existing = get(ctx.nodes, id, nothing)
    if existing !== nothing && existing.operation == node.operation
      push!(stale_ids, id)
      replacement[id] = node.id
    end
  end
  isempty(stale_ids) && return nothing

  for id in stale_ids
    delete!(ctx.nodes, id)
    delete!(ctx.placeholder_ids, id)
    delete!(ctx.node_depths, id)
  end
  for (nid, n) in collect(ctx.nodes)
    new_parent_ids = copy(n.parent_ids)
    changed = false
    for (i, pid) in enumerate(new_parent_ids)
      if haskey(replacement, pid)
        new_parent_ids[i] = replacement[pid]
        changed = true
      end
    end
    changed && (ctx.nodes[nid] = ProvenanceNode(n.id, n.operation, n.parameters, new_parent_ids, n.timestamp))
  end
  return nothing
end

function register_provenance!(ctx::ProvenanceContext, node::ProvenanceNode)
  # Compute read-only state BEFORE any mutation so that a limit-check throw
  # leaves ctx completely untouched (no dangling placeholder_ids, no
  # irreversible prune/rewire of existing nodes).
  is_generic = _is_generic_self_registration(node)
  depth = _provenance_node_depth(ctx, node)

  if ctx.max_nodes !== nothing && !haskey(ctx.nodes, node.id) && length(ctx.nodes) >= ctx.max_nodes
    throw(ArgumentError("ProvenanceContext node limit reached (max_nodes=$(ctx.max_nodes)); refusing to drop node '$(node.operation)'. Increase max_nodes= to record all steps."))
  end
  if ctx.max_depth !== nothing && depth > ctx.max_depth
    throw(ArgumentError("ProvenanceContext depth limit reached (max_depth=$(ctx.max_depth)); refusing to drop node '$(node.operation)'. Increase max_depth= to record full lineage."))
  end

  # All checks passed — now safe to mutate ctx.
  if is_generic
    push!(ctx.placeholder_ids, node.id)
  else
    _prune_generic_self_registrations!(ctx, node)
  end
  ctx.nodes[node.id] = node
  ctx.node_depths[node.id] = depth
  return node
end

function register_provenance!(ctx::ProvenanceContext, node_id::AbstractString, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
  node = ProvenanceNode(String(node_id), String(operation), _provenance_parameters_dict(parameters), String.(parents), _provenance_timestamp())
  return register_provenance!(ctx, node)
end

function register_provenance!(ctx::ProvenanceContext, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
  return register_provenance!(ctx, new_provenance_id(), operation; parents=parents, parameters=parameters)
end

function register_container_provenance!(ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, container, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple(), provenance_hash::Union{Nothing,AbstractString}=nothing)
  metadata = container_provenance_dict(container)
  metadata === nothing && return nothing
  identifier = ensure_provenance_id!(metadata)
  provenance_hash === nothing || (metadata[_metadata_key(metadata, PROVENANCE_HASH_KEY)] = String(provenance_hash))
  ctx === nothing || register_provenance!(ctx, identifier, operation; parents=parents, parameters=parameters)
  return identifier
end

function register_container_provenance!(ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, table::DataFrames.AbstractDataFrame, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple(), provenance_hash::Union{Nothing,AbstractString}=nothing)
  identifier = ensure_provenance_id!(table)
  provenance_hash === nothing || _set_container_provenance_hash!(table, provenance_hash)
  ctx === nothing || register_provenance!(ctx, identifier, operation; parents=parents, parameters=parameters)
  return identifier
end

function _container_provenance_hash(container)
  metadata = container_provenance_dict(container)
  metadata === nothing && return nothing
  if haskey(metadata, PROVENANCE_HASH_KEY)
    value = metadata[PROVENANCE_HASH_KEY]
    return value === nothing ? nothing : String(value)
  elseif haskey(metadata, String(PROVENANCE_HASH_KEY))
    value = metadata[String(PROVENANCE_HASH_KEY)]
    return value === nothing ? nothing : String(value)
  end
  return nothing
end

function _container_provenance_hash(table::DataFrames.AbstractDataFrame)
  value = try
    DataAPI.metadata(table, String(PROVENANCE_HASH_KEY))
  catch
    nothing
  end
  return value === nothing ? nothing : String(value)
end

function _set_container_provenance_id!(container, identifier::AbstractString)
  metadata = container_provenance_dict(container)
  metadata === nothing && return nothing
  metadata[_metadata_key(metadata, PROVENANCE_ID_KEY)] = String(identifier)
  return String(identifier)
end

function _set_container_provenance_id!(table::DataFrames.AbstractDataFrame, identifier::AbstractString)
  DataAPI.metadata!(table, String(PROVENANCE_ID_KEY), String(identifier); style=:note)
  return String(identifier)
end

function _set_container_provenance_hash!(container, provenance_hash::AbstractString)
  metadata = container_provenance_dict(container)
  metadata === nothing && return nothing
  metadata[_metadata_key(metadata, PROVENANCE_HASH_KEY)] = String(provenance_hash)
  return String(provenance_hash)
end

function _set_container_provenance_hash!(table::DataFrames.AbstractDataFrame, provenance_hash::AbstractString)
  DataAPI.metadata!(table, String(PROVENANCE_HASH_KEY), String(provenance_hash); style=:note)
  return String(provenance_hash)
end

function ensure_container_provenance_id!(result::AbstractAnalysisResult)
  if hasproperty(result, :provenance)
    return result.provenance.id
  end
  return nothing
end

function ensure_container_provenance_id!(container)
  metadata = container_provenance_dict(container)
  metadata === nothing && return nothing
  return ensure_provenance_id!(metadata)
end

ensure_container_provenance_id!(table::DataFrames.AbstractDataFrame) = ensure_provenance_id!(table)

"""
    container_provenance_hash(container) → Union{String,Nothing}

Return the SHA-256 content hash attached to `container` by a previous
`stamp_provenance_with_hash!` call, or `nothing` if absent.
"""
container_provenance_hash(container) = _container_provenance_hash(container)
export container_provenance_hash

function container_provenance_summary(container)
  metadata = container_provenance_dict(container)
  metadata === nothing && return "provenance=untracked"
  parts = String[]
  identifier = container_provenance_id(container)
  identifier === nothing || push!(parts, string("id=", identifier))
  if haskey(metadata, PROVENANCE_HASH_KEY)
    push!(parts, string("hash=", metadata[PROVENANCE_HASH_KEY]))
  end
  isempty(parts) && return "provenance=untracked"
  return string("provenance=", join(parts, ", "))
end

function container_provenance_summary(table::DataFrames.AbstractDataFrame)
  parts = String[]
  identifier = container_provenance_id(table)
  provenance_hash = _container_provenance_hash(table)
  identifier === nothing || push!(parts, string("id=", identifier))
  provenance_hash === nothing || push!(parts, string("hash=", provenance_hash))
  isempty(parts) && return "provenance=untracked"
  return string("provenance=", join(parts, ", "))
end

function Base.summary(node::ProvenanceNode)
  return string("ProvenanceNode(", node.operation, ", parents=", length(node.parent_ids), ")")
end

Base.show(io::IO, node::ProvenanceNode) = print(io, summary(node))
Base.show(io::IO, ::MIME"text/plain", node::ProvenanceNode) = print(io, summary(node))

Base.hash(node::ProvenanceNode, h::UInt) = hash(node.id, h)
Base.:(==)(a::ProvenanceNode, b::ProvenanceNode) = a.id == b.id

function Base.summary(ctx::ProvenanceContext)
  return string("ProvenanceContext(", length(ctx.nodes), " nodes)")
end

Base.show(io::IO, ctx::ProvenanceContext) = print(io, summary(ctx))
Base.show(io::IO, ::MIME"text/plain", ctx::ProvenanceContext) = print(io, summary(ctx))

"""
    find_nodes(ctx, operation_pattern) → Vector{ProvenanceNode}
    
Find all nodes whose operation matches `operation_pattern` (regex or exact string).
"""
function find_nodes(ctx::ProvenanceContext, pattern::AbstractString)
  regex = try
    Regex(pattern)
  catch
    ;
    nothing
  end
  return ProvenanceNode[
    node for node in values(ctx.nodes)
             if regex !== nothing ? occursin(regex, node.operation) : node.operation == pattern
  ]
end

function find_nodes(ts::ThreadSafeProvenanceContext, pattern::AbstractString)
  Base.lock(ts.lock) do
    return find_nodes(ts.ctx, pattern)
  end
end
export find_nodes

function export_provenance_json(ctx::ProvenanceContext)
  entity = OrderedDict{String,Any}()
  activity = OrderedDict{String,Any}()
  was_generated_by = OrderedDict{String,Any}()
  was_derived_from = OrderedDict{String,Any}()

  for node in values(ctx.nodes)
    entity[node.id] = OrderedDict(
      "prov:type" => "prov:Entity",
      "prov:label" => node.operation,
      "prov:value" => _provenance_json_value(node.parameters),
    )
    activity[node.id] = OrderedDict(
      "prov:type" => "prov:Activity",
      "prov:label" => node.operation,
      "prov:startedAtTime" => node.timestamp,
      "prov:endedAtTime" => node.timestamp,
    )
    if isempty(node.parent_ids)
      was_generated_by[node.id] = OrderedDict(
        "prov:entity" => node.id,
        "prov:activity" => node.id,
        "prov:time" => node.timestamp,
      )
    else
      for parent_id in node.parent_ids
        was_generated_by[string(parent_id, "->", node.id)] = OrderedDict(
          "prov:entity" => node.id,
          "prov:activity" => node.id,
          "prov:time" => node.timestamp,
        )
        was_derived_from[string(parent_id, "=>", node.id)] = OrderedDict(
          "prov:generatedEntity" => node.id,
          "prov:usedEntity" => parent_id,
        )
      end
    end
  end

  payload = OrderedDict(
    "schema_version" => ctx.schema_version,
    "prefix" => OrderedDict("prov" => "http://www.w3.org/ns/prov#"),
    "environment" => _environment_json(ctx.environment),
    "entity" => entity,
    "activity" => activity,
    "wasGeneratedBy" => was_generated_by,
    "wasDerivedFrom" => was_derived_from,
  )
  return JSON.json(payload)
end

function export_provenance_json(ts::ThreadSafeProvenanceContext)
  snapshot = Base.lock(ts.lock) do
    copy(ts.ctx)
  end
  return export_provenance_json(snapshot)
end

function _environment_json(snapshot::EnvironmentSnapshot)
  return OrderedDict(
    "julia_version" => snapshot.julia_version,
    "manifest_hash" => snapshot.manifest_hash,
    "package_versions" => snapshot.package_versions,
    "cpu_threads" => snapshot.cpu_threads,
    "git_commit" => snapshot.git_commit,
    "random_state_hash" => snapshot.random_state_hash,
    "timestamp" => snapshot.timestamp,
  )
end

function _environment_snapshot_from_payload(payload)
  payload isa AbstractDict || return capture_environment_snapshot()
  return EnvironmentSnapshot(
    string(get(payload, "julia_version", string(VERSION))),
    get(payload, "manifest_hash", nothing) === nothing ? nothing : string(get(payload, "manifest_hash", nothing)),
    Dict{String,String}(string(k) => string(v) for (k, v) in get(payload, "package_versions", Dict{String,String}())),
    Int(get(payload, "cpu_threads", Threads.nthreads())),
    get(payload, "git_commit", nothing) === nothing ? nothing : string(get(payload, "git_commit", nothing)),
    string(get(payload, "random_state_hash", "")),
    string(get(payload, "timestamp", _provenance_timestamp())),
  )
end

mutable struct ProvenanceMigrator
  current_version::VersionNumber
  migrations::Dict{Tuple{VersionNumber,VersionNumber},Function}
end

ProvenanceMigrator(; current_version::VersionNumber=VersionNumber(PROVENANCE_SCHEMA_VERSION)) =
  ProvenanceMigrator(current_version, Dict{Tuple{VersionNumber,VersionNumber},Function}())

const DEFAULT_PROVENANCE_MIGRATOR = ProvenanceMigrator()

@inline function _schema_version_string(version::VersionNumber)
  return version.patch == 0 && isempty(version.prerelease) && isempty(version.build) ?
         "$(version.major).$(version.minor)" : string(version)
end

function register_migration!(migrator::ProvenanceMigrator, transition::Pair{VersionNumber,VersionNumber}, migration::Function)
  source, destination = transition.first, transition.second
  source < destination || throw(ArgumentError("provenance migrations must advance schema versions"))
  haskey(migrator.migrations, (source, destination)) && throw(ArgumentError("migration $(source) => $(destination) is already registered"))
  migrator.migrations[(source, destination)] = migration
  return migrator
end

register_migration!(migration::Function, migrator::ProvenanceMigrator, transition::Pair{VersionNumber,VersionNumber}) =
  register_migration!(migrator, transition, migration)

function _provenance_migration_path(migrator::ProvenanceMigrator, source::VersionNumber, target::VersionNumber)
  source == target && return Tuple{VersionNumber,VersionNumber}[]
  frontier = [(source, Tuple{VersionNumber,VersionNumber}[])]
  visited = Set([source])
  while !isempty(frontier)
    version, path = popfirst!(frontier)
    edges = sort!([edge for edge in keys(migrator.migrations) if edge[1] == version]; by=edge -> edge[2])
    for edge in edges
      next_version = edge[2]
      next_version in visited && continue
      next_path = [path; edge]
      next_version == target && return next_path
      push!(visited, next_version)
      push!(frontier, (next_version, next_path))
    end
  end
  throw(ArgumentError("no provenance migration path from schema $(source) to $(target)"))
end

"""
    migrate_provenance_payload(migrator, payload; target_version=migrator.current_version)

Upgrade a parsed PROV-JSON object through registered, forward-only migrations.
Each migration receives and returns a dictionary-like payload. The framework
sets `schema_version` after every step, preventing a migration from silently
claiming an incorrect schema revision.
"""
function migrate_provenance_payload(
  migrator::ProvenanceMigrator,
  payload::AbstractDict;
  target_version::VersionNumber=migrator.current_version)

  schema = get(payload, "schema_version", nothing)
  schema === nothing && throw(ArgumentError("provenance payload is missing schema_version"))
  source_version = try
    VersionNumber(string(schema))
  catch
    throw(ArgumentError("invalid provenance schema_version=$(schema)"))
  end
  migrated = Dict{String,Any}(string(key) => value for (key, value) in deepcopy(payload))
  for edge in _provenance_migration_path(migrator, source_version, target_version)
    updated = migrator.migrations[edge](migrated)
    updated isa AbstractDict || throw(ArgumentError("provenance migration $(edge[1]) => $(edge[2]) must return an AbstractDict"))
    migrated = Dict{String,Any}(string(key) => value for (key, value) in updated)
    migrated["schema_version"] = _schema_version_string(edge[2])
  end
  return migrated
end

# Topologically sort node IDs so that parents are always registered before
# children — otherwise _provenance_node_depth gets permanently wrong depths
# cached in ctx.node_depths for children that are visited before their parents
# (which is likely when iterating Dict key order from JSON parsing).
function _import_topological_order(all_node_ids::Vector{String}, parents_by_child::Dict{String,Vector{String}})
  id_set = Set(all_node_ids)
  children_of = Dict{String,Vector{String}}()
  indegree = Dict{String,Int}(id => 0 for id in all_node_ids)
  for (child, parents) in parents_by_child
    haskey(indegree, child) || continue
    for parent in parents
      parent in id_set || continue
      indegree[child] += 1
      push!(get!(children_of, parent, String[]), child)
    end
  end
  ready = sort!([id for (id, d) in indegree if d == 0])
  ordered = String[]
  while !isempty(ready)
    id = popfirst!(ready)
    push!(ordered, id)
    for child in sort!(get(children_of, id, String[]))
      indegree[child] -= 1
      indegree[child] == 0 && push!(ready, child)
    end
  end
  # leftover ids (cycles / dangling refs) appended at the end, deterministic order
  append!(ordered, sort!(setdiff(all_node_ids, ordered)))
  return ordered
end

function import_provenance_json(json_text::AbstractString; migrator::ProvenanceMigrator=DEFAULT_PROVENANCE_MIGRATOR)
  payload = try
    JSON.parse(String(json_text); dicttype=Dict)
  catch e
    throw(ArgumentError("import_provenance_json: invalid JSON — $(sprint(showerror, e))"))
  end
  payload isa AbstractDict || throw(ArgumentError("import_provenance_json: expected a JSON object"))
  payload = migrate_provenance_payload(migrator, payload)
  ctx = ProvenanceContext(schema_version=_schema_version_string(migrator.current_version))
  if haskey(payload, "environment")
    ctx.environment = _environment_snapshot_from_payload(payload["environment"])
  end
  entities = get(payload, "entity", Dict())
  activities = get(payload, "activity", Dict())
  parents_by_child = Dict{String,Vector{String}}()
  for rel in values(get(payload, "wasDerivedFrom", Dict()))
    child = string(get(rel, "prov:generatedEntity", ""))
    parent = string(get(rel, "prov:usedEntity", ""))
    isempty(child) || isempty(parent) || push!(get!(parents_by_child, child, String[]), parent)
  end
  all_node_ids = String[string(id) for id in unique(vcat(collect(keys(entities)), collect(keys(activities))))]
  # Topological order ensures parents are registered before children,
  # so _provenance_node_depth correctly looks up already-cached parent depths.
  sorted_ids = _import_topological_order(all_node_ids, parents_by_child)
  for id in sorted_ids
    sid = string(id)
    entity = get(entities, sid, Dict())
    activity = get(activities, sid, Dict())
    operation = string(get(entity, "prov:label", get(activity, "prov:label", "unknown")))
    params = _provenance_parameters_dict(get(entity, "prov:value", Dict()))
    timestamp = string(get(activity, "prov:startedAtTime", _provenance_timestamp()))
    parents = get(parents_by_child, sid, String[])
    register_provenance!(ctx, ProvenanceNode(sid, operation, params, parents, timestamp))
  end
  return ctx
end

function import_provenance_json(::Type{ThreadSafeProvenanceContext}, json_text::AbstractString; kwargs...)
  return ThreadSafeProvenanceContext(import_provenance_json(json_text; kwargs...), ReentrantLock())
end

"""
    Base.empty!(ctx::ProvenanceContext) → ctx

Reset a `ProvenanceContext` in-place: clears all recorded nodes and refreshes
the environment snapshot. The same `ctx` object is returned, so existing
references remain valid.

Useful for long-running services or REPL sessions that reuse a single context:

```julia
enable_provenance!()
# ... run pipeline A ...
println(generate_methods_section(get_provenance_context()))

Base.empty!(get_provenance_context())   # reset for next pipeline
# ... run pipeline B ...
```
"""
function Base.empty!(ctx::ProvenanceContext)
  empty!(ctx.nodes)
  empty!(ctx.placeholder_ids)
  empty!(ctx.node_depths)
  ctx.environment = capture_environment_snapshot()
  return ctx
end

function Base.empty!(ts::ThreadSafeProvenanceContext)
  Base.lock(ts.lock) do
    empty!(ts.ctx)
  end
  return ts
end

const PROVENANCE_METADATA_KEY = "provenance"

@inline function _provenance_field(payload, field::AbstractString, default)
  if payload isa AbstractDict
    haskey(payload, field) && return payload[field]
    symbol_field = Symbol(field)
    haskey(payload, symbol_field) && return payload[symbol_field]
  elseif payload isa NamedTuple
    symbol_field = Symbol(field)
    hasproperty(payload, symbol_field) && return getproperty(payload, symbol_field)
  end
  return default
end

@inline function _provenance_string_vector(value)
  value === nothing && return String[]
  value isa AbstractString && return [String(value)]
  value isa AbstractVector && return String.(value)
  value isa Tuple && return [string(item) for item in value]
  value isa AbstractSet && return [string(item) for item in value]
  return [string(value)]
end

@inline function _provenance_parameters(value)
  value isa NamedTuple && return value
  value isa AbstractDict || return NamedTuple()
  isempty(value) && return NamedTuple()
  names = Tuple(Symbol.(collect(keys(value))))
  vals = Tuple(collect(Base.values(value)))
  return NamedTuple{names}(vals)
end

"""
    provenance_record(label, source; kwargs...)

Construct a provenance record that can be stored in metadata or returned by
data-loading helpers.
"""
function provenance_record(
  label::AbstractString,
  source::AbstractString;
  id::Union{Nothing,AbstractString}=nothing,
  status::Symbol=:ok,
  warnings::AbstractVector{<:AbstractString}=String[],
  errors::AbstractVector{<:AbstractString}=String[],
  fallbacks::AbstractVector{<:AbstractString}=String[],
  notes::AbstractVector{<:AbstractString}=String[],
  parameters::NamedTuple=NamedTuple(),
  parent_ids::AbstractVector{<:AbstractString}=String[],
  timestamp::Union{Nothing,AbstractString}=nothing)
  return ResultProvenance(
    id === nothing ? new_provenance_id() : String(id),
    String(label),
    String(source),
    Symbol(status),
    String.(warnings),
    String.(errors),
    String.(fallbacks),
    String.(notes),
    parameters,
    String.(parent_ids),
    timestamp === nothing ? _provenance_timestamp() : String(timestamp),
  )
end

# Adapter for callers that accept raw metadata payloads or an already
# normalised ResultProvenance without branching at each call site.
function provenance_record(provenance::ResultProvenance;
  id::AbstractString=provenance.id,
  label::AbstractString=provenance.label,
  source::AbstractString=provenance.source,
  status::Symbol=provenance.status,
  warnings::AbstractVector{<:AbstractString}=provenance.warnings,
  errors::AbstractVector{<:AbstractString}=provenance.errors,
  fallbacks::AbstractVector{<:AbstractString}=provenance.fallbacks,
  notes::AbstractVector{<:AbstractString}=provenance.notes,
  parameters::NamedTuple=provenance.parameters,
  parent_ids::AbstractVector{<:AbstractString}=provenance.parent_ids,
  timestamp::AbstractString=provenance.timestamp)
  return provenance_record(label, source; id=id, status=status, warnings=warnings,
    errors=errors, fallbacks=fallbacks, notes=notes, parameters=parameters,
    parent_ids=parent_ids, timestamp=timestamp)
end

function provenance_record(payload::Union{AbstractDict,NamedTuple})
  label = _provenance_field(payload, "label", "untracked")
  source = _provenance_field(payload, "source", "unknown")
  status = _provenance_field(payload, "status", :ok)
  warnings = _provenance_string_vector(_provenance_field(payload, "warnings", String[]))
  errors = _provenance_string_vector(_provenance_field(payload, "errors", String[]))
  fallbacks = _provenance_string_vector(_provenance_field(payload, "fallbacks", String[]))
  notes = _provenance_string_vector(_provenance_field(payload, "notes", String[]))
  parameters = _provenance_parameters(_provenance_field(payload, "parameters", NamedTuple()))
  id = _provenance_field(payload, "id", "")
  parent_ids = _provenance_string_vector(_provenance_field(payload, "parent_ids", String[]))
  timestamp = _provenance_field(payload, "timestamp", _provenance_timestamp())
  return provenance_record(string(label), string(source); id=string(id), status=Symbol(status), warnings=warnings, errors=errors, fallbacks=fallbacks, notes=notes, parameters=parameters, parent_ids=parent_ids, timestamp=string(timestamp))
end

function provenance_summary(provenance::ResultProvenance)
  parts = String[
    provenance.label,
    isempty(provenance.id) ? "id=untracked" : string("id=", provenance.id),
    string("source=", provenance.source),
    string("status=", provenance.status),
    string("timestamp=", provenance.timestamp),
  ]
  !isempty(provenance.parent_ids) && push!(parts, string("parents=", join(provenance.parent_ids, ",")))
  !isempty(provenance.warnings) && push!(parts, string("warnings=", join(provenance.warnings, "; ")))
  !isempty(provenance.errors) && push!(parts, string("errors=", join(provenance.errors, "; ")))
  !isempty(provenance.fallbacks) && push!(parts, string("fallbacks=", join(provenance.fallbacks, "; ")))
  !isempty(provenance.notes) && push!(parts, string("notes=", join(provenance.notes, "; ")))
  length(provenance.parameters) > 0 && push!(parts, string("parameters=", repr(provenance.parameters)))
  return join(parts, " | ")
end

provenance_summary(provenance::Union{AbstractDict,NamedTuple}) = provenance_summary(provenance_record(provenance))
provenance_summary(::Nothing) = "provenance=untracked"

struct WhyProvenanceNode <: AbstractAnalysisResult
  row_index::Int
  expression::String
  values::Dict{String,Any}
  kept::Bool
  provenance::ResultProvenance
end

struct LineageTracer <: AbstractAnalysisResult
  expression::String
  nodes::Vector{WhyProvenanceNode}
  provenance::ResultProvenance
end

function trace_filter(df::DataFrame, predicate::Function; expression::AbstractString="anonymous predicate")
  kept = Bool[]
  nodes = WhyProvenanceNode[]
  for (i, row) in enumerate(eachrow(df))
    row_values = Dict{String,Any}(String(name) => row[name] for name in propertynames(df))
    decision = Bool(predicate(row))
    push!(kept, decision)
    push!(nodes, WhyProvenanceNode(i, String(expression), row_values, decision, provenance_record("WhyProvenanceNode", "analysisresults/trace_filter"; parameters=(row_index=i, kept=decision, expression=String(expression)))))
  end
  result = df[kept, :]
  tracer = LineageTracer(String(expression), nodes, provenance_record("LineageTracer", "analysisresults/trace_filter"; parameters=(row_count=nrow(df), kept_count=nrow(result), expression=String(expression))))
  DataAPI.metadata!(result, "lineage_tracer", tracer; style=:note)
  return result
end

function why_in_result(df::DataFrame, row_index::Integer)
  meta = DataAPI.metadata(df)
  haskey(meta, "lineage_tracer") || throw(ArgumentError("DataFrame does not contain lineage_tracer metadata; create it with trace_filter"))
  tracer = meta["lineage_tracer"]
  tracer isa LineageTracer || throw(ArgumentError("lineage_tracer metadata is not a LineageTracer"))
  kept_nodes = [node for node in tracer.nodes if node.kept]
  1 <= row_index <= length(kept_nodes) || throw(BoundsError(df, row_index))
  return kept_nodes[Int(row_index)]
end




# ==============================================================================
# Provenance storage backends and row-level why-provenance
# ==============================================================================

abstract type AbstractProvenanceBackend end

struct InMemoryProvenanceBackend <: AbstractProvenanceBackend
  nodes::OrderedDict{String,ProvenanceNode}
end

InMemoryProvenanceBackend() = InMemoryProvenanceBackend(OrderedDict{String,ProvenanceNode}())

struct SQLiteProvenanceBackend <: AbstractProvenanceBackend
  path::String
  wal::Bool
  initialized::Bool
end

SQLiteProvenanceBackend(path::AbstractString; wal::Bool=true) = SQLiteProvenanceBackend(String(path), Bool(wal), false)

function register_provenance_backend!(backend::InMemoryProvenanceBackend, node::ProvenanceNode)
  backend.nodes[node.id] = node
  return node
end

function register_provenance_backend!(backend::SQLiteProvenanceBackend, node::ProvenanceNode)
  throw(ArgumentError("SQLiteProvenanceBackend is declared but not active because SQLite.jl is not a hard dependency; use InMemoryProvenanceBackend or enable a future SQLite extension"))
end

function provenance_backend_nodes(backend::InMemoryProvenanceBackend)
  return backend.nodes
end

function provenance_backend_nodes(backend::SQLiteProvenanceBackend)
  throw(ArgumentError("SQLiteProvenanceBackend query support requires the future SQLite extension"))
end

"""
    UncertaintyPropagationResult

First-order propagated variances over a provenance DAG. Each child-parent
sensitivity is a local derivative, and independent input uncertainty is assumed.
This explicit contract avoids presenting unavailable Hessian calculations as a
certainty estimate.
"""
struct UncertaintyPropagationResult <: AbstractAnalysisResult
  node_variances::Dict{String,Float64}
  local_variances::Dict{String,Float64}
  sensitivities::Dict{String,Dict{String,Float64}}
  provenance::ResultProvenance
end

function _normalise_uncertainty_variances(local_variances::AbstractDict)
  result = Dict{String,Float64}()
  for (node_id, variance) in local_variances
    value = Float64(variance)
    isfinite(value) && value >= 0 || throw(ArgumentError("local variance for $(node_id) must be finite and nonnegative"))
    result[String(node_id)] = value
  end
  return result
end

function _normalise_uncertainty_sensitivities(sensitivities::AbstractDict)
  result = Dict{String,Dict{String,Float64}}()
  for (child_id, weights) in sensitivities
    weights isa AbstractDict || throw(ArgumentError("sensitivities for $(child_id) must be an AbstractDict of parent derivatives"))
    child_weights = Dict{String,Float64}()
    for (parent_id, derivative) in weights
      value = Float64(derivative)
      isfinite(value) || throw(ArgumentError("sensitivity $(child_id) <- $(parent_id) must be finite"))
      child_weights[String(parent_id)] = value
    end
    result[String(child_id)] = child_weights
  end
  return result
end

"""
    propagate_uncertainty!(ctx; local_variances, sensitivities=Dict())

Apply the first-order delta method to every node in `ctx`: `Var(y) = Var_local
+ sum((dy/dx)^2 Var(x))`. Parent contributions default to unit sensitivity when
no derivative is supplied. The inferred variance and assumptions are stamped
on each provenance node, while the full result retains the derivative map.
"""
function propagate_uncertainty!(
  ctx::ProvenanceContext;
  local_variances::AbstractDict,
  sensitivities::AbstractDict=Dict{String,Dict{String,Float64}}())

  locals = _normalise_uncertainty_variances(local_variances)
  derivatives = _normalise_uncertainty_sensitivities(sensitivities)
  variances = Dict{String,Float64}()
  visiting = Set{String}()

  function variance_for(node_id::String)
    haskey(variances, node_id) && return variances[node_id]
    node = get(ctx.nodes, node_id, nothing)
    node === nothing && return 0.0
    node_id in visiting && throw(ArgumentError("cannot propagate uncertainty through a cyclic provenance graph"))
    push!(visiting, node_id)
    variance = get(locals, node_id, 0.0)
    weights = get(derivatives, node_id, Dict{String,Float64}())
    for parent_id in node.parent_ids
      derivative = get(weights, parent_id, 1.0)
      variance += derivative^2 * variance_for(parent_id)
    end
    delete!(visiting, node_id)
    variances[node_id] = variance
    parameters = copy(node.parameters)
    parameters["uncertainty_variance"] = variance
    parameters["uncertainty_method"] = "first_order_delta_independent_inputs"
    ctx.nodes[node_id] = ProvenanceNode(node.id, node.operation, parameters, node.parent_ids, node.timestamp)
    return variance
  end

  for node_id in keys(ctx.nodes)
    variance_for(node_id)
  end
  provenance = provenance_record(
    "UncertaintyPropagationResult",
    "analysisresults/propagate_uncertainty!";
    parameters=(node_count=length(variances), local_variance_count=length(locals), sensitivity_count=length(derivatives)),
    notes=["first-order delta method; independent local uncertainty assumption"])
  return UncertaintyPropagationResult(variances, locals, derivatives, provenance)
end

propagate_uncertainty!(ts::ThreadSafeProvenanceContext; kwargs...) = Base.lock(ts.lock) do
  propagate_uncertainty!(ts.ctx; kwargs...)
end


# ==============================================================================
# Shared model-depth result types
# ==============================================================================

struct BackendConfig <: AbstractAnalysisResult
  backend::Symbol
  device::Symbol
  deterministic::Bool
  seed::Union{Nothing,Int}
  parameters::Dict{Symbol,Any}
  provenance::ResultProvenance
end

BackendConfig(backend::Symbol; device::Symbol=:cpu, deterministic::Bool=true, seed::Union{Nothing,Int}=nothing, parameters::Dict{Symbol,Any}=Dict{Symbol,Any}()) =
  BackendConfig(backend, device, deterministic, seed, parameters, provenance_record("BackendConfig", "analysisresults/BackendConfig"; parameters=(backend=backend, device=device, deterministic=deterministic, seed=seed)))

struct ConvergenceReport <: AbstractAnalysisResult
  converged::Bool
  iterations::Int
  objective::Float64
  delta::Float64
  message::String
  provenance::ResultProvenance
end

ConvergenceReport(converged::Bool, iterations::Integer, objective::Real, delta::Real, message::AbstractString="") =
  ConvergenceReport(converged, Int(iterations), Float64(objective), Float64(delta), String(message), provenance_record("ConvergenceReport", "analysisresults/ConvergenceReport"; parameters=(converged=converged, iterations=Int(iterations))))

struct ModelFitDiagnostics <: AbstractAnalysisResult
  method::Symbol
  backend::BackendConfig
  convergence::ConvergenceReport
  hyperparameters::Dict{Symbol,Any}
  training_loss::Vector{Float64}
  warnings::Vector{String}
  provenance::ResultProvenance
end

ModelFitDiagnostics(method::Symbol, backend::BackendConfig, convergence::ConvergenceReport; hyperparameters::Dict{Symbol,Any}=Dict{Symbol,Any}(), training_loss::Vector{Float64}=Float64[], warnings::Vector{String}=String[]) =
  ModelFitDiagnostics(method, backend, convergence, hyperparameters, training_loss, warnings, provenance_record("ModelFitDiagnostics", "analysisresults/ModelFitDiagnostics"; notes=warnings, parameters=(method=method, loss_count=length(training_loss))))

struct ValidationBenchmarkResult <: AbstractAnalysisResult
  benchmark_name::String
  method::Symbol
  metrics::Dict{Symbol,Float64}
  passed::Bool
  provenance::ResultProvenance
end

ValidationBenchmarkResult(benchmark_name::AbstractString, method::Symbol, metrics::Dict{Symbol,<:Real}, passed::Bool) =
  ValidationBenchmarkResult(String(benchmark_name), method, Dict{Symbol,Float64}(k => Float64(v) for (k, v) in metrics), passed, provenance_record("ValidationBenchmarkResult", "analysisresults/ValidationBenchmarkResult"; parameters=(benchmark_name=String(benchmark_name), method=method, passed=passed)))

struct CalibrationCurveResult <: AbstractAnalysisResult
  predicted::Vector{Float64}
  observed::Vector{Float64}
  bin_edges::Vector{Float64}
  method::Symbol
  provenance::ResultProvenance
end

CalibrationCurveResult(predicted::AbstractVector{<:Real}, observed::AbstractVector{<:Real}, bin_edges::AbstractVector{<:Real}, method::Symbol) =
  CalibrationCurveResult(Float64.(predicted), Float64.(observed), Float64.(bin_edges), method, provenance_record("CalibrationCurveResult", "analysisresults/CalibrationCurveResult"; parameters=(method=method, bin_count=max(length(bin_edges)-1, 0))))


@inline _provenance_metadata_key(metadata::AbstractDict) = _metadata_key(metadata, PROVENANCE_METADATA_KEY)

@inline function _metadata_value(metadata::AbstractDict, value)
  value isa valtype(metadata) && return value
  valtype(metadata) <: AbstractString && return provenance_summary(value)
  return value
end

function metadata_provenance(metadata::AbstractDict)
  for key in (PROVENANCE_METADATA_KEY, Symbol(PROVENANCE_METADATA_KEY))
    haskey(metadata, key) || continue
    value = metadata[key]
    value isa ResultProvenance && return value
    value isa Union{AbstractDict,NamedTuple} && return provenance_record(value)
  end
  return nothing
end

function metadata_provenance(table::DataFrames.AbstractDataFrame)
  keys_iter = try
    DataAPI.metadatakeys(table)
  catch
    ()
  end
  for key in keys_iter
    key == PROVENANCE_METADATA_KEY || key == String(PROVENANCE_METADATA_KEY) || continue
    value = try
      DataAPI.metadata(table, String(key))
    catch
      nothing
    end
    value === nothing && continue
    value isa ResultProvenance && return value
    value isa Union{AbstractDict,NamedTuple} && return provenance_record(value)
  end
  return nothing
end

function stamp_provenance!(metadata::AbstractDict, provenance::ResultProvenance)
  identifier = if isempty(provenance.id)
    ensure_provenance_id!(metadata)
  else
    metadata[_metadata_key(metadata, PROVENANCE_ID_KEY)] = provenance.id
    provenance.id
  end
  metadata[_provenance_metadata_key(metadata)] =
    _metadata_value(metadata, isempty(provenance.id) ? provenance_record(provenance; id=identifier) : provenance)
  return metadata
end

function stamp_provenance!(table::DataFrames.AbstractDataFrame, provenance::ResultProvenance)
  identifier = if isempty(provenance.id)
    ensure_provenance_id!(table)
  else
    _set_container_provenance_id!(table, provenance.id)
    provenance.id
  end
  DataAPI.metadata!(table, PROVENANCE_METADATA_KEY,
    isempty(provenance.id) ? provenance_record(provenance; id=identifier) : provenance;
    style=:note)
  return table
end

function stamp_provenance!(metadata::AbstractDict;
  label::AbstractString="data",
  source::AbstractString="unknown",
  status::Symbol=:ok,
  warnings::AbstractVector{<:AbstractString}=String[],
  errors::AbstractVector{<:AbstractString}=String[],
  fallbacks::AbstractVector{<:AbstractString}=String[],
  notes::AbstractVector{<:AbstractString}=String[],
  parameters::NamedTuple=NamedTuple(),
  id::Union{Nothing,AbstractString}=nothing,
  parent_ids::AbstractVector{<:AbstractString}=String[],
  timestamp::Union{Nothing,AbstractString}=nothing)
  return stamp_provenance!(metadata, provenance_record(label, source; id=id, status=status, warnings=warnings, errors=errors, fallbacks=fallbacks, notes=notes, parameters=parameters, parent_ids=parent_ids, timestamp=timestamp))
end

function stamp_provenance!(table::DataFrames.AbstractDataFrame;
  label::AbstractString="data",
  source::AbstractString="unknown",
  status::Symbol=:ok,
  warnings::AbstractVector{<:AbstractString}=String[],
  errors::AbstractVector{<:AbstractString}=String[],
  fallbacks::AbstractVector{<:AbstractString}=String[],
  notes::AbstractVector{<:AbstractString}=String[],
  parameters::NamedTuple=NamedTuple(),
  id::Union{Nothing,AbstractString}=nothing,
  parent_ids::AbstractVector{<:AbstractString}=String[],
  timestamp::Union{Nothing,AbstractString}=nothing)
  return stamp_provenance!(table, provenance_record(label, source; id=id, status=status, warnings=warnings, errors=errors, fallbacks=fallbacks, notes=notes, parameters=parameters, parent_ids=parent_ids, timestamp=timestamp))
end

function update_provenance!(metadata::AbstractDict;
  label::Union{Nothing,AbstractString}=nothing,
  source::Union{Nothing,AbstractString}=nothing,
  status::Union{Nothing,Symbol}=nothing,
  warnings::AbstractVector{<:AbstractString}=String[],
  errors::AbstractVector{<:AbstractString}=String[],
  fallbacks::AbstractVector{<:AbstractString}=String[],
  notes::AbstractVector{<:AbstractString}=String[],
  parameters::NamedTuple=NamedTuple(),
  id::Union{Nothing,AbstractString}=nothing,
  parent_ids::AbstractVector{<:AbstractString}=String[],
  timestamp::Union{Nothing,AbstractString}=nothing)
  current = metadata_provenance(metadata)
  if current === nothing
    return stamp_provenance!(metadata;
      label=label === nothing ? "data" : label,
      source=source === nothing ? "unknown" : source,
      status=status === nothing ? :ok : status,
      warnings=warnings,
      errors=errors,
      fallbacks=fallbacks,
      notes=notes,
      parameters=parameters,
      id=id,
      parent_ids=parent_ids,
      timestamp=timestamp,
    )
  end

  merged = provenance_record(
    label === nothing ? current.label : label,
    source === nothing ? current.source : source;
    id=id === nothing ? current.id : id,
    status=status === nothing ? current.status : status,
    warnings=vcat(current.warnings, String.(warnings)),
    errors=vcat(current.errors, String.(errors)),
    fallbacks=vcat(current.fallbacks, String.(fallbacks)),
    notes=vcat(current.notes, String.(notes)),
    parameters=merge(current.parameters, parameters),
    parent_ids=isempty(parent_ids) ? current.parent_ids : unique(vcat(current.parent_ids, String.(parent_ids))),
    timestamp=timestamp === nothing ? current.timestamp : timestamp,
  )
  return stamp_provenance!(metadata, merged)
end

function update_provenance!(table::DataFrames.AbstractDataFrame;
  label::Union{Nothing,AbstractString}=nothing,
  source::Union{Nothing,AbstractString}=nothing,
  status::Union{Nothing,Symbol}=nothing,
  warnings::AbstractVector{<:AbstractString}=String[],
  errors::AbstractVector{<:AbstractString}=String[],
  fallbacks::AbstractVector{<:AbstractString}=String[],
  notes::AbstractVector{<:AbstractString}=String[],
  parameters::NamedTuple=NamedTuple(),
  id::Union{Nothing,AbstractString}=nothing,
  parent_ids::AbstractVector{<:AbstractString}=String[],
  timestamp::Union{Nothing,AbstractString}=nothing)
  current = metadata_provenance(table)
  if current === nothing
    return stamp_provenance!(table;
      label=label === nothing ? "data" : label,
      source=source === nothing ? "unknown" : source,
      status=status === nothing ? :ok : status,
      warnings=warnings,
      errors=errors,
      fallbacks=fallbacks,
      notes=notes,
      parameters=parameters,
      id=id,
      parent_ids=parent_ids,
      timestamp=timestamp,
    )
  end

  merged = provenance_record(
    label === nothing ? current.label : label,
    source === nothing ? current.source : source;
    id=id === nothing ? current.id : id,
    status=status === nothing ? current.status : status,
    warnings=vcat(current.warnings, String.(warnings)),
    errors=vcat(current.errors, String.(errors)),
    fallbacks=vcat(current.fallbacks, String.(fallbacks)),
    notes=vcat(current.notes, String.(notes)),
    parameters=merge(current.parameters, parameters),
    parent_ids=isempty(parent_ids) ? current.parent_ids : unique(vcat(current.parent_ids, String.(parent_ids))),
    timestamp=timestamp === nothing ? current.timestamp : timestamp,
  )
  return stamp_provenance!(table, merged)
end

function with_provenance(f::Function, ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext})
  _check_provenance_context(ctx, "with_provenance ctx")
  stack = _scoped_provenance_stack()
  push!(stack, ctx)
  try
    return f()
  finally
    pop!(stack)
  end
end

function with_provenance(f::Function)
  ctx = ProvenanceContext()
  return with_provenance(f, ctx)
end

function scoped_provenance_context()
  return _active_scoped_provenance_context()
end

function with_provenance(result::AbstractAnalysisResult, label::AbstractString, source::AbstractString;
  prov_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}=nothing,
  operation::AbstractString=source,
  parents::AbstractVector{<:AbstractString}=String[],
  id::Union{Nothing,AbstractString}=nothing,
  status::Symbol=:ok,
  warnings::AbstractVector{<:AbstractString}=String[],
  errors::AbstractVector{<:AbstractString}=String[],
  fallbacks::AbstractVector{<:AbstractString}=String[],
  notes::AbstractVector{<:AbstractString}=String[],
  parameters::NamedTuple=NamedTuple(),
  parent_ids::AbstractVector{<:AbstractString}=String[],
  timestamp::Union{Nothing,AbstractString}=nothing,
  provenance_hash::Union{Nothing,AbstractString}=nothing,
  _register_graph::Bool=true,
)
  resolved_ctx = active_provenance_context(prov_ctx)
  if _register_graph && resolved_ctx !== nothing
    graph_parents = isempty(parent_ids) ? parents : parent_ids
    return provenance_result!(resolved_ctx, result, operation; parents=graph_parents, parameters=parameters, label=label, source=source, status=status, warnings=warnings, errors=errors, fallbacks=fallbacks, notes=notes, provenance_hash=provenance_hash)
  end

  # Build the new ResultProvenance
  provenance_id = id === nothing ? (hasproperty(result, :provenance) ? result.provenance.id : new_provenance_id()) : id
  ts = timestamp === nothing ? _provenance_timestamp() : timestamp
  new_prov = ResultProvenance(
    String(provenance_id),
    String(label),
    String(source),
    status,
    String.(warnings),
    String.(errors),
    String.(fallbacks),
    String.(notes),
    parameters,
    String.(parent_ids),
    String(ts)
  )

  # Reconstruct the struct with the new provenance field
  T = typeof(result)
  new_fields = ntuple(Val(fieldcount(T))) do i
    f = fieldname(T, i)
    f === :provenance ? new_prov : getfield(result, f)
  end
  return T(new_fields...)
end

function with_provenance(result, label::AbstractString, source::AbstractString;
  prov_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}=nothing,
  operation::AbstractString=source,
  parents::AbstractVector{<:AbstractString}=String[],
  id::Union{Nothing,AbstractString}=nothing,
  status::Symbol=:ok,
  warnings::AbstractVector{<:AbstractString}=String[],
  errors::AbstractVector{<:AbstractString}=String[],
  fallbacks::AbstractVector{<:AbstractString}=String[],
  notes::AbstractVector{<:AbstractString}=String[],
  parameters::NamedTuple=NamedTuple(),
  parent_ids::AbstractVector{<:AbstractString}=String[],
  timestamp::Union{Nothing,AbstractString}=nothing,
  provenance_hash::Union{Nothing,AbstractString}=nothing,
  _register_graph::Bool=true,
)
  resolved_ctx = active_provenance_context(prov_ctx)
  if _register_graph && resolved_ctx !== nothing
    graph_parents = isempty(parent_ids) ? parents : parent_ids
    return provenance_result!(resolved_ctx, result, operation; parents=graph_parents, parameters=parameters, label=label, source=source, status=status, warnings=warnings, errors=errors, fallbacks=fallbacks, notes=notes, provenance_hash=provenance_hash)
  end

  if result isa DataFrames.AbstractDataFrame
    provenance_hash === nothing || _set_container_provenance_hash!(result, provenance_hash)
    update_provenance!(result; label=label, source=source, id=id, status=status, warnings=warnings, errors=errors, fallbacks=fallbacks, notes=notes, parameters=parameters, parent_ids=parent_ids, timestamp=timestamp)
    return result
  elseif result isa AbstractVector && !isempty(result)
    # Intentional design: a Vector is the *output* of a batch operation.
    # We stamp provenance on each element only if it already carries its own
    # provenance dict (opt-in containers). Elements that do NOT implement
    # container_provenance_dict are left untouched — the caller is responsible
    # for attaching batch-level provenance at a higher layer.
    # The batch id is NOT shared across all elements to avoid false equivalence;
    # each element gets its own fresh id so it can be individually traced later.
    for item in result
      metadata = container_provenance_dict(item)
      metadata === nothing && continue
      provenance_hash === nothing || _set_container_provenance_hash!(item, provenance_hash)
      # id=nothing so each element gets its own distinct provenance id
      update_provenance!(metadata; label=label, source=source, id=nothing,
        status=status, warnings=warnings, errors=errors, fallbacks=fallbacks,
        notes=notes, parameters=parameters, parent_ids=parent_ids, timestamp=timestamp)
    end
    return result
  else
    metadata = container_provenance_dict(result)
    if metadata !== nothing
      provenance_hash === nothing || _set_container_provenance_hash!(result, provenance_hash)
      update_provenance!(metadata; label=label, source=source, id=id, status=status, warnings=warnings, errors=errors, fallbacks=fallbacks, notes=notes, parameters=parameters, parent_ids=parent_ids, timestamp=timestamp)
    end
    return result
  end
end

@inline function provenance_result!(::Nothing, result, operation::AbstractString; kwargs...)
  return result
end

function provenance_result!(prov_ctx::ProvenanceContext, result, operation::AbstractString;
  parents::AbstractVector{<:AbstractString}=String[],
  parameters=NamedTuple(),
  id::Union{Nothing,AbstractString}=nothing,
  label::Union{Nothing,AbstractString}=nothing,
  source::Union{Nothing,AbstractString}=nothing,
  status::Symbol=:ok,
  warnings::AbstractVector{<:AbstractString}=String[],
  errors::AbstractVector{<:AbstractString}=String[],
  fallbacks::AbstractVector{<:AbstractString}=String[],
  notes::AbstractVector{<:AbstractString}=String[],
  provenance_hash::Union{Nothing,AbstractString}=nothing,
)
  parent_ids = String.(parents)
  identifier = if id === nothing
    ensure_container_provenance_id!(result)
  else
    requested_id = String(id)
    _set_container_provenance_id!(result, requested_id)
    requested_id
  end
  if provenance_hash !== nothing
    _set_container_provenance_hash!(result, provenance_hash)
  end

  timestamp = _provenance_timestamp()
  if prov_ctx !== nothing
    node_id = identifier === nothing ? new_provenance_id() : identifier
    node = ProvenanceNode(String(node_id), String(operation), _provenance_parameters_dict(parameters), parent_ids, timestamp)
    register_provenance!(prov_ctx, node)
    identifier = node.id
  end

  if identifier !== nothing || result isa DataFrames.AbstractDataFrame || container_provenance_dict(result) !== nothing || (result isa AbstractVector && any(item -> container_provenance_dict(item) !== nothing, result))
    return with_provenance(result,
      label === nothing ? String(operation) : String(label),
      source === nothing ? String(operation) : String(source);
      id=identifier,
      status=status,
      warnings=warnings,
      errors=errors,
      fallbacks=fallbacks,
      notes=notes,
      parameters=_provenance_parameters(parameters),
      parent_ids=parent_ids,
      timestamp=timestamp,
      provenance_hash=provenance_hash,
      _register_graph=false,
    )
  end
  return result
end

Base.summary(provenance::ResultProvenance) = provenance_summary(provenance)
Base.show(io::IO, provenance::ResultProvenance) = print(io, provenance_summary(provenance))
Base.show(io::IO, ::MIME"text/plain", provenance::ResultProvenance) = print(io, provenance_summary(provenance))

"""
    analysis_result_fields(::Type{T}) where T

Return the fields that should appear in the shared analysis-result summary and
fallback tabular view.
"""
analysis_result_fields(::Type{T}) where {T} = Tuple(name for name in fieldnames(T) if name != :provenance)

@inline function _analysis_result_value_preview(value)
  return (value isa AbstractArray || value isa AbstractDict) ? summary(value) : repr(value)
end

@inline function _analysis_result_field_names(result)
  return hasproperty(result, :fieldnames) ? propertynames(result) : fieldnames(typeof(result))
end

@inline function _analysis_result_safe_getproperty(result, field_name)
  return hasproperty(result, field_name) ? getproperty(result, field_name) : nothing
end

analysis_result_type(result) = analysis_result_type(typeof(result))
analysis_result_type(::Type) = :generic

function _analysis_result_label_parts(result)
  kind = analysis_result_type(result)
  kind === :generic ? String[] : [string(kind)]
end

analysis_result_provenance_label_parts(::Any) = String[]
analysis_result_provenance_parameters(::Any) = NamedTuple()

function analysis_result_provenance(result)
  explicit = _analysis_result_safe_getproperty(result, :provenance)
  if explicit isa ResultProvenance
    return explicit
  elseif explicit isa Union{AbstractDict,NamedTuple}
    return provenance_record(explicit)
  end

  result_type = typeof(result)
  type_name = string(nameof(result_type))
  source_name = string(parentmodule(result_type))
  warnings = String[]
  errors = String[]
  fallbacks = String[]
  notes = String[]
  parameters = NamedTuple()
  label_parts = _analysis_result_label_parts(result)

  method = _analysis_result_safe_getproperty(result, :method)
  if method !== nothing && !isempty(string(method))
    push!(label_parts, string(method))
    parameters = merge(parameters, (; method=method))
  end

  reference_name = _analysis_result_safe_getproperty(result, :reference_name)
  if reference_name !== nothing && !isempty(string(reference_name))
    push!(label_parts, string(reference_name))
    parameters = merge(parameters, (; reference_name=reference_name))
  end

  phenotype_name = _analysis_result_safe_getproperty(result, :phenotype_name)
  if phenotype_name !== nothing && !isempty(string(phenotype_name))
    push!(label_parts, string(phenotype_name))
    parameters = merge(parameters, (; phenotype_name=phenotype_name))
  end

  system = _analysis_result_safe_getproperty(result, :system)
  if system !== nothing && !isempty(string(system))
    push!(label_parts, string(system))
    parameters = merge(parameters, (; system=system))
  end

  target_gene = _analysis_result_safe_getproperty(result, :target_gene)
  if target_gene !== nothing && !isempty(string(target_gene))
    push!(label_parts, string("target=", target_gene))
    parameters = merge(parameters, (; target_gene=target_gene))
  end

  term_id = _analysis_result_safe_getproperty(result, :term_id)
  if term_id !== nothing && !isempty(string(term_id))
    push!(label_parts, string("term=", term_id))
    parameters = merge(parameters, (; term_id=term_id))
  end

  db = _analysis_result_safe_getproperty(result, :db)
  if db !== nothing && !isempty(string(db))
    push!(label_parts, string("db=", db))
    parameters = merge(parameters, (; db=db))
  end

  modalities = _analysis_result_safe_getproperty(result, :modalities)
  if modalities !== nothing
    push!(label_parts, "multimodal")
    parameters = merge(parameters, (; modalities=modalities))
  end

  extra_label_parts = analysis_result_provenance_label_parts(result)
  !isempty(extra_label_parts) && append!(label_parts, string.(extra_label_parts))
  parameters = merge(parameters, analysis_result_provenance_parameters(result))

  source = isempty(label_parts) ? source_name : string(source_name, "/", join(label_parts, "/"))
  label = isempty(label_parts) ? type_name : string(type_name, " [", join(label_parts, ", "), "]")

  converged = _analysis_result_safe_getproperty(result, :converged)
  if converged isa Bool && !converged
    push!(warnings, "fit did not converge")
    push!(fallbacks, "best available estimate was retained")
  end

  zero_inflated = _analysis_result_safe_getproperty(result, :zero_inflated)
  if zero_inflated isa Bool && zero_inflated
    push!(notes, "zero-inflated model path used")
  end

  passes_filters = _analysis_result_safe_getproperty(result, :passes_filters)
  if passes_filters isa Bool && !passes_filters
    push!(warnings, "result did not pass all filters")
  end

  status = isempty(errors) ? (isempty(warnings) ? :ok : :warn) : :error
  return provenance_record(label, source; status=status, warnings=warnings, errors=errors, fallbacks=fallbacks, notes=notes, parameters=parameters)
end

function analysis_result_summary(result)
  provenance = analysis_result_provenance(result)
  result_type = typeof(result)
  field_names = analysis_result_fields(result_type)
  parts = String[]
  for field_name in field_names
    value = getfield(result, field_name)
    push!(parts, string(field_name, "=", _analysis_result_value_preview(value)))
  end

  payload = string(provenance.label, " | id=", isempty(provenance.id) ? "untracked" : provenance.id, " | source=", provenance.source, " | status=", provenance.status, " | ", join(parts, ", "))
  if !isempty(provenance.warnings)
    payload = string(payload, " | warnings=", join(provenance.warnings, "; "))
  end
  if !isempty(provenance.fallbacks)
    payload = string(payload, " | fallbacks=", join(provenance.fallbacks, "; "))
  end
  if !isempty(provenance.notes)
    payload = string(payload, " | notes=", join(provenance.notes, "; "))
  end
  if length(provenance.parameters) > 0
    payload = string(payload, " | parameters=", repr(provenance.parameters))
  end
  return payload
end

function Base.summary(result::AbstractAnalysisResult)
  return analysis_result_summary(result)
end

function Base.show(io::IO, result::AbstractAnalysisResult)
  print(io, analysis_result_summary(result))
end

function Base.show(io::IO, ::MIME"text/plain", result::AbstractAnalysisResult)
  print(io, analysis_result_summary(result))
  provenance = analysis_result_provenance(result)
  if !isempty(provenance.errors) || !isempty(provenance.warnings) || !isempty(provenance.fallbacks)
    print(io, "\nprovenance:")
    if !isempty(provenance.warnings)
      print(io, "\n  warnings: ", join(provenance.warnings, "; "))
    end
    if !isempty(provenance.fallbacks)
      print(io, "\n  fallbacks: ", join(provenance.fallbacks, "; "))
    end
    if !isempty(provenance.errors)
      print(io, "\n  errors: ", join(provenance.errors, "; "))
    end
  end
end

"""
    DataFrames.DataFrame(result::AbstractAnalysisResult)

Convert an analysis result struct to a single-row `DataFrame`.
Provenance is attached via `DataAPI.metadata!` (`:note` style), **not** as
repeated columns. This avoids O(n×columns) memory waste for large result
tables in genomics workloads.

To explicitly export provenance as columns (e.g., for CSV export), call
`provenance_lineage_table(ctx)` or flatten with `hcat(df, provenance_df)`.
"""
function DataFrames.DataFrame(result::AbstractAnalysisResult)
  result_type = typeof(result)
  field_names = analysis_result_fields(result_type)
  columns = ntuple(index -> [getfield(result, field_names[index])], length(field_names))
  table = DataFrames.DataFrame(NamedTuple{field_names}(columns))
  # Attach provenance as metadata only — zero extra columns, zero RAM waste
  provenance = analysis_result_provenance(result)
  DataAPI.metadata!(table, PROVENANCE_METADATA_KEY, provenance; style=:note)
  return table
end

# ==============================================================================
# Extended Provenance Engine — high-level convenience APIs
# All added utilities are zero-cost when prov_ctx === nothing.
# ==============================================================================


"""
    provenance_chain(ctx, node_id) → Vector{ProvenanceNode}

Retrieve the full ancestor chain (breadth-first, oldest-first) of `node_id`
within `ctx`. Returns an empty vector if `node_id` is not found.
"""
function provenance_chain(ctx::ProvenanceContext, node_id::AbstractString)::Vector{ProvenanceNode}
  chain = ProvenanceNode[]
  visited = Set{String}()
  queue = String[String(node_id)]
  while !isempty(queue)
    id = popfirst!(queue)
    id ∈ visited && continue
    push!(visited, id)
    haskey(ctx.nodes, id) || continue
    node = ctx.nodes[id]
    push!(chain, node)
    append!(queue, node.parent_ids)
  end
  return reverse!(chain)
end

"""
    provenance_lineage_table(ctx) → DataFrame

Return a `DataFrame` summary of every node in `ctx`, one row per operation.
Columns: `id`, `operation`, `timestamp`, `n_parents`, `parents`, `parameters`.
"""
function provenance_lineage_table(ctx::ProvenanceContext)::DataFrames.DataFrame
  rows = NamedTuple[]
  for node in values(ctx.nodes)
    push!(rows, (
      id=node.id,
      operation=node.operation,
      timestamp=node.timestamp,
      n_parents=length(node.parent_ids),
      parents=join(node.parent_ids, ","),
      parameters=isempty(node.parameters) ? "{}" :
                 JSON.json(Dict(string(k) => v for (k, v) in pairs(node.parameters))),
    ))
  end
  if isempty(rows)
    return DataFrames.DataFrame(
      id=String[], operation=String[], timestamp=String[],
      n_parents=Int[], parents=String[], parameters=String[])
  end
  return DataFrames.DataFrame(rows)
end

"""
    provenance_diff(ctx_a, ctx_b)

Compare two `ProvenanceContext` objects. Returns a NamedTuple with:
- `only_in_a`  — nodes present only in `ctx_a`
- `only_in_b`  — nodes present only in `ctx_b`
- `common`     — nodes present in both (by id)
"""
function provenance_diff(ctx_a::ProvenanceContext, ctx_b::ProvenanceContext)
  ids_a = Set(keys(ctx_a.nodes))
  ids_b = Set(keys(ctx_b.nodes))
  only_a = [ctx_a.nodes[id] for id in sort!(collect(setdiff(ids_a, ids_b)))]
  only_b = [ctx_b.nodes[id] for id in sort!(collect(setdiff(ids_b, ids_a)))]
  common = [ctx_a.nodes[id] for id in sort!(collect(intersect(ids_a, ids_b)))]
  return (only_in_a=only_a, only_in_b=only_b, common=common)
end

function merge_provenance_contexts(contexts::ProvenanceContext...; max_nodes::Union{Int,Nothing}=nothing, max_depth::Union{Int,Nothing}=nothing)
  merged = ProvenanceContext(; max_nodes=max_nodes, max_depth=max_depth)
  for ctx in contexts
    for node in values(ctx.nodes)
      register_provenance!(merged, node)
    end
  end
  return merged
end

"""
    DistributedProvenanceContext

A serializable collection of worker-local provenance DAGs. A worker may create
nodes independently; merging is deterministic by `(worker_id, local_clock,
timestamp, node_id)` after respecting any explicit parent edges.
"""
mutable struct DistributedProvenanceContext
  worker_ids::Vector{Int}
  local_contexts::Dict{Int,ProvenanceContext}
  merge_strategy::Symbol
end

function DistributedProvenanceContext(worker_ids::AbstractVector{<:Integer}; merge_strategy::Symbol=:topological)
  merge_strategy == :topological || throw(ArgumentError("only merge_strategy=:topological is supported"))
  ids = sort!(unique(Int.(worker_ids)))
  isempty(ids) && throw(ArgumentError("DistributedProvenanceContext requires at least one worker ID"))
  return DistributedProvenanceContext(ids, Dict(id => ProvenanceContext() for id in ids), merge_strategy)
end

function worker_provenance_context(context::DistributedProvenanceContext, worker_id::Integer)
  id = Int(worker_id)
  haskey(context.local_contexts, id) || throw(ArgumentError("worker $(id) is not registered in this DistributedProvenanceContext"))
  return context.local_contexts[id]
end

function _node_equivalent(left::ProvenanceNode, right::ProvenanceNode)
  return left.id == right.id && left.operation == right.operation &&
         left.parameters == right.parameters && left.parent_ids == right.parent_ids &&
         left.timestamp == right.timestamp
end

function _distributed_merge_order(context::DistributedProvenanceContext)
  nodes = Dict{String,Tuple{ProvenanceNode,Int,Int}}()
  for worker_id in context.worker_ids
    local_context = context.local_contexts[worker_id]
    for (clock, node) in enumerate(values(local_context.nodes))
      existing = get(nodes, node.id, nothing)
      if existing === nothing
        nodes[node.id] = (node, worker_id, clock)
      elseif !_node_equivalent(existing[1], node)
        throw(ArgumentError("conflicting provenance node ID $(node.id) from workers $(existing[2]) and $(worker_id)"))
      end
    end
  end

  indegree = Dict{String,Int}(node_id => 0 for node_id in keys(nodes))
  children = Dict{String,Vector{String}}()
  for (node_id, (node, _, _)) in nodes
    for parent_id in node.parent_ids
      haskey(nodes, parent_id) || continue
      indegree[node_id] += 1
      push!(get!(children, parent_id, String[]), node_id)
    end
  end
  order_key = node_id -> begin
    node, worker_id, clock = nodes[node_id]
    (worker_id, clock, node.timestamp, node_id)
  end
  ready = sort!([node_id for (node_id, degree) in indegree if degree == 0]; by=order_key)
  ordered = ProvenanceNode[]
  while !isempty(ready)
    node_id = popfirst!(ready)
    push!(ordered, nodes[node_id][1])
    for child_id in get(children, node_id, String[])
      indegree[child_id] -= 1
      if indegree[child_id] == 0
        # Binary insertion into the sorted ready list instead of re-sorting
        # the entire list on every iteration — O(n log n) total vs O(n² log n).
        key = order_key(child_id)
        lo, hi = 1, length(ready) + 1
        while lo < hi
          mid = (lo + hi) >> 1
          if order_key(ready[mid]) < key
            lo = mid + 1
          else
            hi = mid
          end
        end
        insert!(ready, lo, child_id)
      end
    end
  end
  length(ordered) == length(nodes) || throw(ArgumentError("distributed provenance graph contains a cycle"))
  return ordered
end

"""
    merge_distributed_provenance!(master_ctx, distributed_ctx)

Merge worker-local provenance into `master_ctx`. Identical duplicate nodes are
idempotent; conflicting duplicate identifiers and cycles are errors. The method
never uses wall-clock ordering as causal evidence: timestamps only break ties
between independent worker events after local clocks and dependency edges.
"""
function merge_distributed_provenance!(master_ctx::ProvenanceContext, distributed_ctx::DistributedProvenanceContext)
  for worker_id in distributed_ctx.worker_ids
    worker_provenance_context(distributed_ctx, worker_id).schema_version == master_ctx.schema_version ||
      throw(ArgumentError("worker $(worker_id) uses a different provenance schema version"))
  end
  for node in _distributed_merge_order(distributed_ctx)
    existing = get(master_ctx.nodes, node.id, nothing)
    if existing === nothing
      register_provenance!(master_ctx, node)
    elseif !_node_equivalent(existing, node)
      throw(ArgumentError("master provenance context already contains conflicting node ID $(node.id)"))
    end
  end
  return master_ctx
end

function merge_distributed_provenance!(master_ctx::ThreadSafeProvenanceContext, distributed_ctx::DistributedProvenanceContext)
  return Base.lock(master_ctx.lock) do
    merge_distributed_provenance!(master_ctx.ctx, distributed_ctx)
  end
end

"""
    @dist_provenance distributed_ctx worker_id operation parameters expr

Register an expression result in the selected worker-local context. This macro
is serializable for use inside `Distributed.@everywhere` definitions because it
does not capture process-local state beyond the supplied context object.
"""
macro dist_provenance(context_expr, worker_expr, operation_expr, parameters_expr, body_expr)
  context_symbol = gensym("distributed_context")
  worker_symbol = gensym("worker_id")
  result_symbol = gensym("distributed_result")
  quote
    local $(context_symbol) = $(esc(context_expr))
    local $(worker_symbol) = $(esc(worker_expr))
    local $(result_symbol) = $(esc(body_expr))
    provenance_result!(
      worker_provenance_context($(context_symbol), $(worker_symbol)),
      $(result_symbol),
      String($(esc(operation_expr)));
      parameters=($(esc(parameters_expr))))
  end
end

function _provenance_graph_label(value::AbstractString)
  return replace(value, "\\" => "\\\\", "\"" => "\\\"", "\n" => " ")
end

function provenance_to_dot(ctx::ProvenanceContext)::String
  io = IOBuffer()
  println(io, "digraph provenance {")
  println(io, "  rankdir=LR;")
  seen_edges = Set{Tuple{String,String}}()
  for node in values(ctx.nodes)
    println(io, "  \"", node.id, "\" [label=\"", _provenance_graph_label(node.operation), "\"];")
    for parent_id in node.parent_ids
      edge = (parent_id, node.id)
      if !(edge in seen_edges)
        println(io, "  \"", parent_id, "\" -> \"", node.id, "\";")
        push!(seen_edges, edge)
      end
    end
  end
  println(io, "}")
  return String(take!(io))
end

function _provenance_mermaid_id(id::AbstractString)
  token = bytes2hex(sha256(codeunits(String(id))))[1:12]
  return string("n", token)
end

function provenance_to_mermaid(ctx::ProvenanceContext)::String
  io = IOBuffer()
  println(io, "flowchart LR")
  seen = Set{String}()
  for node in values(ctx.nodes)
    node_ref = _provenance_mermaid_id(node.id)
    if !(node.id in seen)
      println(io, "  ", node_ref, "[\"", _provenance_graph_label(node.operation), "\"]")
      push!(seen, node.id)
    end
    for parent_id in node.parent_ids
      parent_ref = _provenance_mermaid_id(parent_id)
      if !(parent_id in seen)
        # Use operation label when parent node exists; fall back to
        # full ID substring only for dangling references
        parent_label = if haskey(ctx.nodes, parent_id)
          _provenance_graph_label(ctx.nodes[parent_id].operation)
        else
          # `end` is not valid outside array indexing; use length()
          first(parent_id, min(length(parent_id), 16)) * "…"
        end
        println(io, "  ", parent_ref, "[\"", parent_label, "\"]")
        push!(seen, parent_id)
      end
      println(io, "  ", parent_ref, " --> ", node_ref)
    end
  end
  return String(take!(io))
end

"""
    provenance_structural_hash(data) → String

Fast **structural** hash — encodes only array dimensions, element type, and
scalar values. Does **not** inspect element values for arrays/matrices.

Use this for:
- Checking whether matrix *shapes* changed between pipeline steps.
- Annotating provenance nodes in hot paths where hashing bytes would be too slow.
- Any dataset larger than a few MB.

!!! warning "Not a content hash"
    Two arrays with identical shape but different values will return the same
    hash. For value-level integrity checking use `provenance_content_hash`.
"""
function provenance_structural_hash(data)::String
  token = if data isa AbstractString
    # Strings: full content is small, hash it for real
    String(data)
  elseif data isa AbstractArray
    # Shape + eltype only — O(1), safe for any array size
    string(size(data), ':', eltype(data), ':', length(data))
  elseif data isa Number
    string(data)
  elseif data isa AbstractDict
    string(length(data), ':', keytype(data), ':', valtype(data))
  elseif data isa DataFrames.AbstractDataFrame
    string(size(data), ':', names(data))
  else
    string(typeof(data), ':', objectid(data))
  end
  return bytes2hex(SHA.sha256(codeunits(token)))
end

"""
    provenance_content_hash(data) → String

True **SHA-256 content hash** — streams through element values so that two
datasets with the same shape but different values produce different hashes.

This is the correct function to use for **data integrity verification**.

!!! warning "Performance"
    This function iterates over all elements of arrays and DataFrames.
    For large genomics datasets (e.g., a 5 GB count matrix or VCF) it will
    block the calling thread for seconds to minutes. Profile before use in
    production pipelines. For large data, prefer `file_provenance_hash` on
    the source file, or `provenance_structural_hash` as a cheap shape check.

# Examples
```julia
# Small result — safe to content-hash
hash = provenance_content_hash(de_results)

# Large matrix — use structural instead
hash = provenance_structural_hash(count_matrix)

# Source file — stream without reading into RAM
hash = file_provenance_hash("/data/sample.vcf.gz")
```
"""
function provenance_content_hash(data)::String
  ctx_sha = SHA.SHA2_256_CTX()
  if data isa AbstractString
    SHA.update!(ctx_sha, codeunits(String(data)))
  elseif data isa AbstractVector{UInt8}
    SHA.update!(ctx_sha, data)
  elseif data isa AbstractArray{T} where T<:Union{Float16,Float32,Float64,
    Int8,Int16,Int32,Int64,Int128,
    UInt8,UInt16,UInt32,UInt64,UInt128}
    # Hash shape first so arrays with same values but different dimensions
    # (e.g. 2×3 vs 6-element vector) produce distinct hashes.
    SHA.update!(ctx_sha, codeunits(string(size(data))))
    # Restrict pointer fast path to concrete Array to avoid MethodError on
    # SubArray and other AbstractArray subtypes where sizeof/pointer may not
    # be defined.
    if data isa Array && IndexStyle(data) === IndexLinear()
      GC.@preserve data begin
        p = convert(Ptr{UInt8}, pointer(data))
        sz = sizeof(data)
        if p != C_NULL && sz > 0
          SHA.update!(ctx_sha, unsafe_wrap(Array, p, sz))
        end
      end
    else
      for val in data
        SHA.update!(ctx_sha, codeunits(repr(val)))
      end
    end
  elseif data isa AbstractArray
    # Non-numeric: hash shape + repr() which is at least stable for a given Julia version
    SHA.update!(ctx_sha, codeunits(string(size(data))))
    for val in data
      SHA.update!(ctx_sha, codeunits(repr(val)))
    end
  elseif data isa DataFrames.AbstractDataFrame
    SHA.update!(ctx_sha, codeunits(string(size(data), ':', names(data))))
    for col in eachcol(data)
      if isbitstype(eltype(col)) && eltype(col) <: Union{Float16,Float32,Float64,
        Int8,Int16,Int32,Int64,Int128,
        UInt8,UInt16,UInt32,UInt64,UInt128}
        SHA.update!(ctx_sha, reinterpret(UInt8, col))
      else
        for val in col
          SHA.update!(ctx_sha, codeunits(repr(val)))
        end
      end
    end
  elseif data isa AbstractDict
    for (k, v) in sort!(collect(pairs(data)); by=x->string(x[1]))
      SHA.update!(ctx_sha, codeunits(repr(k)))
      SHA.update!(ctx_sha, codeunits(repr(v)))
    end
  elseif data isa Union{Number,Bool}
    SHA.update!(ctx_sha, codeunits(repr(data)))
  else
    SHA.update!(ctx_sha, codeunits(repr(data)))
  end
  return bytes2hex(SHA.digest!(ctx_sha))
end

"""
    file_provenance_hash(path) → String

Compute a **streaming SHA-256 hash** of the file at `path` without reading the
entire file into memory. Safe to use on multi-GB FASTA, VCF, BAM, or HDF5 files.

This is the recommended way to fingerprint raw input files for provenance:
```julia
ctx = ProvenanceContext()
register_provenance!(ctx, "load_vcf";
    parameters=(path=path, sha256=file_provenance_hash(path),))
vcf = read_vcf(path)
```

!!! note
    I/O speed is typically the bottleneck, not CPU. On NVMe storage this
    processes ~2–4 GB/s; on network mounts it is bounded by network throughput.
"""
function file_provenance_hash(path::AbstractString)::String
  isfile(path) || throw(ArgumentError("file_provenance_hash: path does not exist: $path"))
  ctx_sha = SHA.SHA2_256_CTX()
  buf = Vector{UInt8}(undef, 65536)  # 64 KiB read buffer
  open(path, "r") do io
    while !eof(io)
      n = readbytes!(io, buf)
      SHA.update!(ctx_sha, view(buf, 1:n))
    end
  end
  return bytes2hex(SHA.digest!(ctx_sha))
end


"""
    lazy_provenance_id!(ctx, result, operation; parents, parameters) → result

Like `provenance_result!` but skips the SHA call and dict insert entirely when
`ctx === nothing`. Use this in very-high-frequency hot paths where even the
method dispatch of `provenance_result!(nothing, ...)` must be avoided.
"""
@inline function lazy_provenance_id!(
  ctx::Nothing, result, ::AbstractString;
  parents=String[], parameters=NamedTuple())
  return result
end

function lazy_provenance_id!(
  ctx::ProvenanceContext, result, operation::AbstractString;
  parents::AbstractVector{<:AbstractString}=String[],
  parameters=NamedTuple())
  return provenance_result!(ctx, result, operation; parents=parents, parameters=parameters)
end

function lazy_provenance_id!(
  ctx::ThreadSafeProvenanceContext, result, operation::AbstractString;
  parents::AbstractVector{<:AbstractString}=String[],
  parameters=NamedTuple())
  return provenance_result!(ctx, result, operation; parents=parents, parameters=parameters)
end

"""
    @provenance ctx operation parameters expr

Zero-cost macro wrapper that stamps provenance into `ctx` for the result of
`expr` when `ctx !== nothing`. When `ctx === nothing` (default), expands to a
no-op — no allocation, no hash, no timestamp. Julia will constant-fold the
branch away entirely when the context is provably `nothing`.

# Example
```julia
result = @provenance prov_ctx "normalize" (method=:pqn,) begin
    normalize_counts(data)
end
```
"""
function _prov_parent_stack()
  get!(task_local_storage(), :_prov_parent_stack) do
    String[]
  end::Vector{String}
end

function _push_parent!(id::String)
  push!(_prov_parent_stack(), id)
end
function _pop_parent!()
  pop!(_prov_parent_stack())
end
function _current_parents()
  copy(_prov_parent_stack())
end

"""
    @provenance ctx operation parents parameters expr

...
!!! warning "Task-local scope"
    Auto-parent linking uses `task_local_storage` and works correctly with
    `@threads` but **not** with `@spawn`. For `@spawn`-based parallelism,
    pass `parents=` explicitly.
"""
macro provenance(ctx_expr, operation_expr, parameters_expr, body_expr)
  ctx_sym = gensym("_prov_ctx")
  res_sym = gensym("_result")
  id_sym = gensym("_prov_id")
  quote
    local $(ctx_sym) = $(esc(ctx_expr))
    if $(ctx_sym) === nothing
      $(esc(body_expr))
    else
      ($(ctx_sym) isa ProvenanceContext || $(ctx_sym) isa ThreadSafeProvenanceContext) ||
        throw(ArgumentError("@provenance ctx must be a ProvenanceContext, ThreadSafeProvenanceContext, or nothing"))
      local $(id_sym) = new_provenance_id()
      _push_parent!($(id_sym))
      try
        local $(res_sym) = $(esc(body_expr))
        provenance_result!($(ctx_sym), $(res_sym), String($(esc(operation_expr)));
          id=($(id_sym)),
          parents=_current_parents()[1:(end-1)],
          parameters=($(esc(parameters_expr))))
      finally
        _pop_parent!()
      end
    end
  end
end

"""
    @provenance ctx operation parameters expr

...
!!! warning "Task-local scope"
    Auto-parent linking uses `task_local_storage` and works correctly with
    `@threads` but **not** with `@spawn`. For `@spawn`-based parallelism,
    pass `parents=` explicitly.
"""
macro provenance(ctx_expr, operation_expr, parents_expr, parameters_expr, body_expr)
  ctx_sym = gensym("_prov_ctx")
  res_sym = gensym("_result")
  id_sym = gensym("_prov_id")
  quote
    local $(ctx_sym) = $(esc(ctx_expr))
    if $(ctx_sym) === nothing
      $(esc(body_expr))
    else
      ($(ctx_sym) isa ProvenanceContext || $(ctx_sym) isa ThreadSafeProvenanceContext) ||
        throw(ArgumentError("@provenance ctx must be a ProvenanceContext, ThreadSafeProvenanceContext, or nothing"))
      local $(id_sym) = new_provenance_id()
      _push_parent!($(id_sym))
      try
        local $(res_sym) = $(esc(body_expr))
        provenance_result!($(ctx_sym), $(res_sym), String($(esc(operation_expr)));
          id=($(id_sym)),
          parents=String.($(esc(parents_expr))),
          parameters=($(esc(parameters_expr))))
      finally
        _pop_parent!()
      end
    end
  end
end

# Core types — exported so users write ProvenanceContext() not BioToolkit.ProvenanceContext()
export AbstractAnalysisResult
export AbstractProvenanceBackend, InMemoryProvenanceBackend, SQLiteProvenanceBackend
export register_provenance_backend!, provenance_backend_nodes, WhyProvenanceNode, LineageTracer, trace_filter, why_in_result
export ProvenanceMigrator, DEFAULT_PROVENANCE_MIGRATOR, register_migration!, migrate_provenance_payload
export UncertaintyPropagationResult, propagate_uncertainty!
export ModelFitDiagnostics, ConvergenceReport, ValidationBenchmarkResult, CalibrationCurveResult, BackendConfig
export ProvenanceContext, ProvenanceNode, ResultProvenance, EnvironmentSnapshot
export ThreadSafeProvenanceContext
export provenance_chain, provenance_lineage_table, provenance_diff
export import_provenance_json, export_provenance_json, capture_environment_snapshot
export merge_provenance_contexts, provenance_to_dot, provenance_to_mermaid, provenance_to_html, export_provenance_html
export DistributedProvenanceContext, worker_provenance_context, merge_distributed_provenance!, @dist_provenance
export analysis_result_type, analysis_result_fields, analysis_result_summary, analysis_result_provenance
export provenance_structural_hash, provenance_content_hash, file_provenance_hash
export lazy_provenance_id!, new_provenance_id
export @provenance
# Provenance stamping / querying — used across all BioToolkit modules
export stamp_provenance!, update_provenance!, with_provenance, provenance_result!
export register_provenance!, register_container_provenance!
export active_provenance_context
export container_provenance_id, container_provenance_summary, ensure_provenance_id!
export metadata_provenance, provenance_record, provenance_summary, provenance_parent_ids

# ------------------------------------------------------------------------------
# 1. ProvenanceParams
# ------------------------------------------------------------------------------
"""
    ProvenanceParams

Union type accepted by all provenance APIs for the `parameters` argument.

- `NamedTuple`       — preferred for hot paths; struct-like, compiler-optimised,
                       zero allocation on the fast path.
- `Dict{String,Any}` — preferred for dynamic pipelines where parameter sets are
                       not known at compile time (e.g., loaded from a config file).
"""
const ProvenanceParams = Union{NamedTuple,Dict{String,Any}}
export ProvenanceParams

# ------------------------------------------------------------------------------
# 2. Type-safe dispatch: has_provenance / requires_provenance
# ------------------------------------------------------------------------------
"""
    has_provenance(result) → Bool

Return `true` if `result` carries embedded provenance metadata (either as a
`ResultProvenance` field or via `DataAPI.metadata!`). Use as a lightweight
compile-friendly guard in visualisation or validation functions.
"""
function has_provenance(result)::Bool
  if hasproperty(result, :provenance) && getproperty(result, :provenance) isa ResultProvenance
    return true
  end
  if result isa DataFrames.AbstractDataFrame
    return metadata_provenance(result) !== nothing
  end
  d = container_provenance_dict(result)
  return d !== nothing && metadata_provenance(d) !== nothing
end

"""
    requires_provenance(result) → result

Assert that `result` carries provenance, throwing `ArgumentError` otherwise.
Pipeline guard for functions that require full lineage traceability.
"""
function requires_provenance(result)
  has_provenance(result) && return result
  T = typeof(result)
  throw(ArgumentError(
    "Result of type $T does not carry provenance metadata. " *
    "Pass a ProvenanceContext via prov_ctx= when calling the upstream function."
  ))
end
export has_provenance, requires_provenance

# ------------------------------------------------------------------------------
# 3. DAG traversal: provenance_ancestors / provenance_descendants
# ------------------------------------------------------------------------------
"""
    provenance_ancestors(ctx, node_id) → Vector{ProvenanceNode}

All ancestor nodes of `node_id` in `ctx` (BFS, oldest-first, excluding self).
Complements `provenance_chain` for inverse / upstream traversal.
"""
function provenance_ancestors(ctx::ProvenanceContext, node_id::AbstractString)::Vector{ProvenanceNode}
  visited = Set{String}()
  queue = String[node_id]
  result = ProvenanceNode[]
  while !isempty(queue)
    cid = popfirst!(queue)
    cid in visited && continue
    push!(visited, cid)
    node = get(ctx.nodes, cid, nothing)
    node === nothing && continue
    cid != node_id && push!(result, node)
    for pid in node.parent_ids
      pid in visited || push!(queue, pid)
    end
  end
  return result
end

"""
    provenance_descendants(ctx, node_id) → Vector{ProvenanceNode}

All descendant nodes of `node_id` in `ctx` (BFS, shallowest-first).
Useful for forward impact analysis: "what downstream results depended on this step?"
"""
function provenance_descendants(ctx::ProvenanceContext, node_id::AbstractString)::Vector{ProvenanceNode}
  children = Dict{String,Vector{String}}()
  for (id, node) in ctx.nodes
    for pid in node.parent_ids
      push!(get!(children, pid, String[]), id)
    end
  end
  visited = Set{String}()
  queue = String[node_id]
  result = ProvenanceNode[]
  while !isempty(queue)
    cid = popfirst!(queue)
    cid in visited && continue
    push!(visited, cid)
    for child_id in get(children, cid, String[])
      child_id in visited && continue
      node = get(ctx.nodes, child_id, nothing)
      node === nothing && continue
      push!(result, node)
      push!(queue, child_id)
    end
  end
  return result
end
export provenance_ancestors, provenance_descendants

# ------------------------------------------------------------------------------
# 4. @pure_provenance — deterministic zero-cost macro for pure functions
# ------------------------------------------------------------------------------
"""
    @pure_provenance ctx operation content_key parameters expr

Like `@provenance` but generates a **deterministic** content-addressable
provenance ID by hashing `content_key`. Two calls with the same `content_key`
always produce the same node ID — enabling reproducible audit trails across
independent runs. The macro stamps the result container when possible and
records the corresponding provenance node in `ctx`.

```julia
result = @pure_provenance prov_ctx "normalize" string(method, size(mat)) (method=:zscore,) begin
    normalize_counts(mat)
end
```
"""
macro pure_provenance(ctx_expr, operation_expr, content_key_expr, parameters_expr, body_expr)
  quote
    local _prov_ctx_v = $(esc(ctx_expr))
    if _prov_ctx_v === nothing
      $(esc(body_expr))
    else
      (_prov_ctx_v isa ProvenanceContext || _prov_ctx_v isa ThreadSafeProvenanceContext) ||
        throw(ArgumentError("@pure_provenance ctx must be a ProvenanceContext, ThreadSafeProvenanceContext, or nothing"))
      local _det_id = new_provenance_id(String($(esc(content_key_expr))))
      _push_parent!(_det_id)
      try
        local _result_v = $(esc(body_expr))
        provenance_result!(_prov_ctx_v, _result_v, String($(esc(operation_expr)));
          id=_det_id,
          parents=_current_parents()[1:(end-1)],
          parameters=($(esc(parameters_expr))))
      finally
        _pop_parent!()
      end
    end
  end
end
export @pure_provenance

# ------------------------------------------------------------------------------
# 5. stamp_provenance_with_hash! — content-hash reference instead of data copy
# ------------------------------------------------------------------------------
"""
    stamp_provenance_with_hash!(container, data; kwargs...) → container

Stamp `container` with provenance **and** attach a SHA-256 content-hash of
`data` as `provenance_hash`. This avoids memory-doubling large sequences,
matrices, or DataFrames by storing only a compact fingerprint instead of the
raw object.

```julia
stamp_provenance_with_hash!(result_df, input_matrix;
    label="normalize_output", source="Metabolomics/normalize")
```
"""
function stamp_provenance_with_hash!(container, data;
  label::AbstractString="data",
  source::AbstractString="unknown",
  status::Symbol=:ok,
  warnings::AbstractVector{<:AbstractString}=String[],
  errors::AbstractVector{<:AbstractString}=String[],
  fallbacks::AbstractVector{<:AbstractString}=String[],
  notes::AbstractVector{<:AbstractString}=String[],
  parameters::ProvenanceParams=NamedTuple(),
  parent_ids::AbstractVector{<:AbstractString}=String[],
  id::Union{Nothing,AbstractString}=nothing,
  hash_mode::Symbol=:structural)

  # Choose hashing strategy based on hash_mode:
  #   :structural (default) — O(1), hashes shape/type only. Safe for any size.
  #   :content              — O(n), true SHA-256 over element values.
  #                           Only use for small-to-medium data where correctness
  #                           matters more than speed.
  content_hash = if hash_mode === :content
    provenance_content_hash(data)
  else
    # :structural or any unrecognised symbol → structural (fast, safe)
    hash_mode === :structural || @warn "stamp_provenance_with_hash!: unknown hash_mode=$hash_mode, defaulting to :structural"
    provenance_structural_hash(data)
  end

  _set_container_provenance_hash!(container, content_hash)
  params_nt = parameters isa NamedTuple ? parameters :
              NamedTuple(Symbol(k) => v for (k, v) in parameters)
  return stamp_provenance!(container;
    label=label, source=source, status=status,
    warnings=warnings, errors=errors, fallbacks=fallbacks,
    notes=notes, parameters=params_nt,
    parent_ids=parent_ids, id=id)
end
export stamp_provenance_with_hash!

# ==============================================================================
# 6. generate_methods_section — provenance → academic paper draft
# ==============================================================================

# Internal: render a parameter dict as a compact, human-readable string.
# e.g. {"method"=>"tmm", "n_genes"=>500} → "method: tmm, n_genes: 500"
function _methods_params_string(params::Dict{String,Any})::String
  isempty(params) && return ""
  parts = String[]
  for k in sort!(collect(keys(params)))
    v = params[k]
    push!(parts, string(k, ": ", v isa AbstractDict ? "{…}" : string(v)))
  end
  return join(parts, ", ")
end

# Internal: convert a snake_case / CamelCase operation name into readable prose.
# "normalize_counts"       → "normalize counts"
# "pairwise_align_codons"  → "pairwise align codons"
# "SCTransform"            → "SCTransform"
function _methods_operation_label(op::AbstractString)::String
  # Replace underscores with spaces; keep existing capitalisation
  label = replace(String(op), '_' => ' ')
  # Strip internal module prefix (e.g. "BioToolkit/normalize" → "normalize")
  if '/' in label
    label = split(label, '/')[end]
  end
  return strip(label)
end

# Internal: topological sort of DAG nodes (Kahn's algorithm, root-to-leaf).
# Returns nodes ordered from earliest (no parents) to latest.
function _methods_topological_order(ctx::ProvenanceContext)::Vector{ProvenanceNode}
  in_degree = Dict{String,Int}()
  children_of = Dict{String,Vector{String}}()
  for (id, node) in ctx.nodes
    in_degree[id] = count(pid -> haskey(ctx.nodes, pid), node.parent_ids)
    for pid in node.parent_ids
      push!(get!(children_of, pid, String[]), id)
    end
  end

  # Sort roots by timestamp so the Methods section is deterministic across runs
  roots = [id for (id, deg) in in_degree if deg == 0]
  sort!(roots; by=id -> ctx.nodes[id].timestamp)
  queue = roots
  ordered = ProvenanceNode[]
  while !isempty(queue)
    id = popfirst!(queue)
    haskey(ctx.nodes, id) || continue
    push!(ordered, ctx.nodes[id])
    children = get(children_of, id, String[])
    # Sort children by timestamp for determinism within each level.
    sort!(children; by=cid -> begin
      n = get(ctx.nodes, cid, nothing)
      n !== nothing ? n.timestamp : ""
    end)
    for child_id in children
      if haskey(in_degree, child_id)
        in_degree[child_id] -= 1
        in_degree[child_id] == 0 && push!(queue, child_id)
      end
    end
  end
  visited_ids = Set{String}(n.id for n in ordered)
  for (id, node) in ctx.nodes
    id in visited_ids || push!(ordered, node)
  end
  return ordered
end

# Internal: render one ProvenanceNode as a Methods sentence.
# Picks the lead verb intelligently:
#   - First step with a loading operation → "Data were obtained using"
#   - First step with computation         → "Data were processed using"
#   - All subsequent steps                → "This was followed by"
const _LOAD_VERB_RE = r"load|read|import|parse|fetch|stream|open"i
function _methods_sentence(node::ProvenanceNode, index::Int)::String
  label = _methods_operation_label(node.operation)
  params = _methods_params_string(node.parameters)
  lead = if index == 1 && occursin(_LOAD_VERB_RE, label)
    "Data were obtained using"
  elseif index == 1
    "Data were processed using"
  else
    "This was followed by"
  end
  if isempty(params)
    return string(lead, " **", label, "**.")
  else
    return string(lead, " **", label, "** (", params, ").")
  end
end

# Internal: build the closing citation / software sentence from the environment.
function _methods_citation_sentence(env::EnvironmentSnapshot)::String
  parts = String["BioToolkit.jl"]
  env.julia_version != "" && push!(parts, string("Julia v", env.julia_version))
  if env.manifest_hash !== nothing
    push!(parts, string("Manifest SHA-256: ", env.manifest_hash[1:12], "…"))
  end
  if env.git_commit !== nothing
    push!(parts, string("git commit: ", env.git_commit[1:8], "…"))
  end
  return string("All analyses were performed using ",
    join(parts, "; "), ".")
end

"""
    generate_methods_section(ctx::ProvenanceContext; kwargs...) → String
    generate_methods_section(result;               kwargs...) → String

Generate a **draft Methods-section paragraph** from a provenance graph,
suitable for pasting directly into an academic manuscript.

The function walks the provenance DAG in chronological order (roots first,
using a topological sort) and renders each recorded operation as a
natural-language sentence.  The closing sentence cites the BioToolkit version,
Julia version, Manifest hash, and git commit captured at the time of analysis.

# Arguments
- `ctx`     — a `ProvenanceContext` built up during the analysis pipeline.
- `result`  — any object that carries an embedded `ProvenanceContext` via
              `DataAPI.metadata` or a `:provenance` field.  The context is
              extracted automatically.
- `title`   — optional section heading line (default: `"## Methods"`).
              Pass `title=""` to suppress the heading.
- `max_ops` — cap the number of rendered operations (default: `nothing` = all).

# Example
```julia
ctx = ProvenanceContext()
counts = normalize_counts(raw; method=:tmm, prov_ctx=ctx)
pcs    = run_pca(counts; n_components=50, prov_ctx=ctx)
result = find_clusters(pcs; resolution=0.5, prov_ctx=ctx)

println(generate_methods_section(ctx))
```

Produces output similar to:
```
## Methods

Data were processed using **normalize counts** (method: tmm).
This was followed by **run pca** (n_components: 50).
This was followed by **find clusters** (resolution: 0.5).
All analyses were performed using BioToolkit.jl; Julia v1.12.5;
Manifest SHA-256: 3a9f2c11b4d0…; git commit: 4a1bcdef….
```
"""
function generate_methods_section(
  ctx::ProvenanceContext;
  title::AbstractString="## Methods",
  max_ops::Union{Nothing,Int}=nothing,
)::String
  io = IOBuffer()

  # Heading
  if !isempty(title)
    println(io, title)
    println(io)
  end

  if isempty(ctx.nodes)
    println(io, "_No provenance operations were recorded._")
    println(io)
    println(io, _methods_citation_sentence(ctx.environment))
    return String(take!(io))
  end

  ordered = _methods_topological_order(ctx)
  n_show = max_ops === nothing ? length(ordered) : min(max_ops, length(ordered))

  for (i, node) in enumerate(ordered[1:n_show])
    println(io, _methods_sentence(node, i))
  end

  if max_ops !== nothing && length(ordered) > max_ops
    omitted = length(ordered) - max_ops
    println(io, "_… (", omitted, " additional step",
      omitted == 1 ? "" : "s", " omitted; see full provenance graph)._")
  end

  println(io)
  println(io, _methods_citation_sentence(ctx.environment))

  return String(take!(io))
end

# Convenience overload for result containers (DataFrame, AbstractAnalysisResult, etc.)
function generate_methods_section(
  result;
  title::AbstractString="## Methods",
  max_ops::Union{Nothing,Int}=nothing,
)::String
  # Try extracting context from DataAPI metadata
  if result isa DataFrames.AbstractDataFrame
    for key in ("_provenance_context", "provenance_context")
      ctx_candidate = try
        DataAPI.metadata(result, key)
      catch
        ;
        nothing
      end
      if ctx_candidate isa ProvenanceContext
        return generate_methods_section(ctx_candidate; title=title, max_ops=max_ops)
      end
    end
  end
  # Try :provenance_context field
  if hasproperty(result, :provenance_context) &&
     getproperty(result, :provenance_context) isa ProvenanceContext
    return generate_methods_section(
      getproperty(result, :provenance_context); title=title, max_ops=max_ops)
  end
  # Fallback: build a stub section from the embedded ResultProvenance if present
  prov = analysis_result_provenance(result)
  env = capture_environment_snapshot()
  io = IOBuffer()
  !isempty(title) && (println(io, title); println(io))
  println(io, _methods_sentence(
    ProvenanceNode(prov.id, prov.label, _provenance_parameters_dict(prov.parameters),
      prov.parent_ids, prov.timestamp), 1))
  println(io)
  println(io, _methods_citation_sentence(env))
  return String(take!(io))
end

function provenance_chain(ts::ThreadSafeProvenanceContext, node_id::AbstractString)::Vector{ProvenanceNode}
  Base.lock(ts.lock) do
    return provenance_chain(ts.ctx, node_id)
  end
end

function provenance_lineage_table(ts::ThreadSafeProvenanceContext)::DataFrames.DataFrame
  snapshot = Base.lock(ts.lock) do
    copy(ts.ctx)
  end
  return provenance_lineage_table(snapshot)
end

function provenance_diff(ctx_a::ThreadSafeProvenanceContext, ctx_b::ProvenanceContext)
  Base.lock(ctx_a.lock) do
    return provenance_diff(ctx_a.ctx, ctx_b)
  end
end

function provenance_diff(ctx_a::ProvenanceContext, ctx_b::ThreadSafeProvenanceContext)
  Base.lock(ctx_b.lock) do
    return provenance_diff(ctx_a, ctx_b.ctx)
  end
end

function provenance_diff(ctx_a::ThreadSafeProvenanceContext, ctx_b::ThreadSafeProvenanceContext)
  first, second = objectid(ctx_a) <= objectid(ctx_b) ? (ctx_a, ctx_b) : (ctx_b, ctx_a)
  Base.lock(first.lock) do
    Base.lock(second.lock) do
      return provenance_diff(ctx_a.ctx, ctx_b.ctx)
    end
  end
end

function merge_provenance_contexts(contexts::Union{ProvenanceContext,ThreadSafeProvenanceContext}...; max_nodes::Union{Int,Nothing}=nothing, max_depth::Union{Int,Nothing}=nothing)
  merged = ProvenanceContext(; max_nodes=max_nodes, max_depth=max_depth)
  for ctx in contexts
    if ctx isa ThreadSafeProvenanceContext
      Base.lock(ctx.lock) do
        for node in values(ctx.ctx.nodes)
          register_provenance!(merged, node)
        end
      end
    else
      for node in values(ctx.nodes)
        register_provenance!(merged, node)
      end
    end
  end
  return merged
end

function provenance_to_dot(ts::ThreadSafeProvenanceContext)::String
  Base.lock(ts.lock) do
    return provenance_to_dot(ts.ctx)
  end
end

function provenance_to_mermaid(ts::ThreadSafeProvenanceContext)::String
  Base.lock(ts.lock) do
    return provenance_to_mermaid(ts.ctx)
  end
end

function provenance_ancestors(ts::ThreadSafeProvenanceContext, node_id::AbstractString)::Vector{ProvenanceNode}
  Base.lock(ts.lock) do
    return provenance_ancestors(ts.ctx, node_id)
  end
end

function provenance_descendants(ts::ThreadSafeProvenanceContext, node_id::AbstractString)::Vector{ProvenanceNode}
  Base.lock(ts.lock) do
    return provenance_descendants(ts.ctx, node_id)
  end
end

function generate_methods_section(
  ts::ThreadSafeProvenanceContext;
  title::AbstractString="## Methods",
  max_ops::Union{Nothing,Int}=nothing,
)::String
  snapshot = Base.lock(ts.lock) do
    return copy(ts.ctx)
  end
  return generate_methods_section(snapshot; title=title, max_ops=max_ops)
end

export generate_methods_section

# Snapshot Context
Base.copy(ctx::ProvenanceContext) = ProvenanceContext(
  OrderedDict(
    k => ProvenanceNode(
      v.id,
      v.operation,
      deepcopy(v.parameters),
      copy(v.parent_ids),
      v.timestamp,
    ) for (k, v) in ctx.nodes
  ),
  EnvironmentSnapshot(
    ctx.environment.julia_version,
    ctx.environment.manifest_hash,
    copy(ctx.environment.package_versions),
    ctx.environment.cpu_threads,
    ctx.environment.git_commit,
    ctx.environment.random_state_hash,
    ctx.environment.timestamp,
  ),
  ctx.max_nodes,
  ctx.max_depth,
  ctx.schema_version,
  copy(ctx.placeholder_ids),
  copy(ctx.node_depths),
)

# Flatten Provenance for CSV Export
export flatten_provenance!
function flatten_provenance!(df::DataFrames.AbstractDataFrame)
  pid = container_provenance_id(df)
  phash = _container_provenance_hash(df)
  if pid !== nothing
    df[!, :provenance_id] .= pid
  end
  if phash !== nothing
    df[!, :provenance_hash] .= phash
  end
  return df
end

Base.iterate(ctx::ProvenanceContext) = iterate(ctx.nodes)
Base.iterate(ctx::ProvenanceContext, state) = iterate(ctx.nodes, state)
Base.length(ctx::ProvenanceContext) = length(ctx.nodes)

# ==============================================================================
# Interactive HTML5 Provenance Visualizer
# ==============================================================================

# HTML/script escaping helpers — same pattern already applied in align.jl and
# annotation.jl visualizers, now consistently applied here.
function _prov_escape_html(s::AbstractString)
  s = replace(s, '&' => "&amp;")
  s = replace(s, '<' => "&lt;")
  s = replace(s, '>' => "&gt;")
  s = replace(s, '"' => "&quot;")
  s = replace(s, '\'' => "&#39;")
  return s
end

@inline function _prov_escape_script_json(s::AbstractString)
  return replace(s, "</" => "<\\/")
end

"""
    provenance_to_html(ctx; title="BioToolkit Provenance Visualizer", width="100%", height="650px", dark_mode=true, standalone=true) → String

Generate an interactive HTML5 visualizer for `ctx` featuring DAG graph visualization,
search/filtering, node details sidebar, PROV-JSON view, Mermaid view, and Methods text.
"""
function provenance_to_html(ctx::ProvenanceContext;
  title::AbstractString="BioToolkit Provenance Visualizer",
  width::AbstractString="100%",
  height::AbstractString="650px",
  dark_mode::Bool=true,
  standalone::Bool=true)

  inst_id = "bioprov_" * bytes2hex(rand(UInt8, 4))

  nodes_json = JSON.json([
    Dict(
      "id" => node.id,
      "label" => node.operation,
      "operation" => node.operation,
      "timestamp" => node.timestamp,
      "parent_ids" => node.parent_ids,
      "parameters" => node.parameters
    ) for node in values(ctx.nodes)
  ])

  edges_json = JSON.json([
    Dict(
      "from" => pid,
      "to" => node.id
    ) for node in values(ctx.nodes) for pid in node.parent_ids
  ])

  prov_json_val = JSON.json(export_provenance_json(ctx))
  mermaid_val = JSON.json(provenance_to_mermaid(ctx))
  methods_val = JSON.json(generate_methods_section(ctx; title=""))

  bg_color = dark_mode ? "#0f172a" : "#f8fafc"
  card_bg = dark_mode ? "#1e293b" : "#ffffff"
  text_color = dark_mode ? "#f1f5f9" : "#0f172a"
  border_color = dark_mode ? "#334155" : "#e2e8f0"
  accent_color = dark_mode ? "#38bdf8" : "#0284c7"
  sub_text = dark_mode ? "#94a3b8" : "#64748b"
  node_bg = dark_mode ? "#3b82f6" : "#2563eb"

  safe_title = _prov_escape_html(title)
  # Escape JSON payloads for safe embedding in <script> blocks — prevent
  # early tag termination if the serialised data happens to contain "</script>".
  safe_nodes_json = _prov_escape_script_json(nodes_json)
  safe_edges_json = _prov_escape_script_json(edges_json)
  safe_prov_json_val = _prov_escape_script_json(prov_json_val)
  safe_mermaid_val = _prov_escape_script_json(mermaid_val)
  safe_methods_val = _prov_escape_script_json(methods_val)

  html = """
  $(standalone ? "<!DOCTYPE html><html><head><meta charset='utf-8'><title>$(safe_title)</title>" : "")
  <style>
    #container-$(inst_id) {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
      background: $(bg_color);
      color: $(text_color);
      border: 1px solid $(border_color);
      border-radius: 12px;
      overflow: hidden;
      width: $(width);
      box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.3);
      margin: 10px 0;
    }
    #container-$(inst_id) .bioprov-header {
      padding: 16px 20px;
      background: $(card_bg);
      border-bottom: 1px solid $(border_color);
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 12px;
    }
    #container-$(inst_id) .bioprov-title-area h2 {
      margin: 0;
      font-size: 1.15rem;
      font-weight: 700;
      color: $(accent_color);
      display: flex;
      align-items: center;
      gap: 8px;
    }
    #container-$(inst_id) .bioprov-meta-badges {
      display: flex;
      gap: 8px;
      font-size: 0.78rem;
    }
    #container-$(inst_id) .bioprov-badge {
      background: $(border_color);
      color: $(text_color);
      padding: 4px 10px;
      border-radius: 9999px;
      font-weight: 500;
    }
    #container-$(inst_id) .bioprov-toolbar {
      display: flex;
      gap: 10px;
      align-items: center;
    }
    #container-$(inst_id) .bioprov-search {
      padding: 6px 12px;
      border-radius: 6px;
      border: 1px solid $(border_color);
      background: $(bg_color);
      color: $(text_color);
      font-size: 0.85rem;
      outline: none;
      width: 180px;
    }
    #container-$(inst_id) .bioprov-tab-btn {
      background: transparent;
      border: 1px solid $(border_color);
      color: $(text_color);
      padding: 6px 12px;
      border-radius: 6px;
      cursor: pointer;
      font-size: 0.85rem;
      transition: all 0.2s;
    }
    #container-$(inst_id) .bioprov-tab-btn.active {
      background: $(accent_color);
      color: #ffffff;
      border-color: $(accent_color);
    }
    #container-$(inst_id) .bioprov-body {
      display: flex;
      height: $(height);
      position: relative;
    }
    #container-$(inst_id) .bioprov-graph-view {
      flex: 1;
      height: 100%;
      background: $(bg_color);
      position: relative;
    }
    #container-$(inst_id) .bioprov-sidebar {
      width: 320px;
      background: $(card_bg);
      border-left: 1px solid $(border_color);
      padding: 16px;
      overflow-y: auto;
      font-size: 0.85rem;
    }
    #container-$(inst_id) .bioprov-sidebar h3 {
      margin-top: 0;
      font-size: 1rem;
      color: $(accent_color);
      border-bottom: 1px solid $(border_color);
      padding-bottom: 8px;
    }
    #container-$(inst_id) .bioprov-param-table {
      width: 100%;
      border-collapse: collapse;
      margin-top: 10px;
    }
    #container-$(inst_id) .bioprov-param-table th, #container-$(inst_id) .bioprov-param-table td {
      border: 1px solid $(border_color);
      padding: 6px 8px;
      text-align: left;
      font-size: 0.8rem;
    }
    #container-$(inst_id) .bioprov-text-pane {
      display: none;
      width: 100%;
      height: 100%;
      padding: 20px;
      box-sizing: border-box;
      background: $(card_bg);
      color: $(text_color);
      font-family: 'Fira Code', 'Courier New', monospace;
      font-size: 0.85rem;
      white-space: pre-wrap;
      overflow: auto;
    }
  </style>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.js"></script>
  <link href="https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.css" rel="stylesheet" type="text/css" />

  <div class="bioprov-container" id="container-$(inst_id)">
    <div class="bioprov-header">
      <div class="bioprov-title-area">
        <h2>🧬 $(safe_title)</h2>
        <div class="bioprov-meta-badges">
          <span class="bioprov-badge">Nodes: $(length(ctx.nodes))</span>
          <span class="bioprov-badge">Julia v$(ctx.environment.julia_version)</span>
          <span class="bioprov-badge">Threads: $(ctx.environment.cpu_threads)</span>
        </div>
      </div>
      <div class="bioprov-toolbar">
        <input type="text" class="bioprov-search" id="search-$(inst_id)" placeholder="Search operations..." oninput="window.bioprovFilterNodes_$(inst_id)(this.value)">
        <button class="bioprov-tab-btn active" onclick="window.bioprovSetTab_$(inst_id)('graph', this)">DAG Graph</button>
        <button class="bioprov-tab-btn" onclick="window.bioprovSetTab_$(inst_id)('json', this)">PROV-JSON</button>
        <button class="bioprov-tab-btn" onclick="window.bioprovSetTab_$(inst_id)('mermaid', this)">Mermaid</button>
        <button class="bioprov-tab-btn" onclick="window.bioprovSetTab_$(inst_id)('methods', this)">Methods</button>
      </div>
    </div>
    <div class="bioprov-body">
      <div id="graph-$(inst_id)" class="bioprov-graph-view"></div>
      <div id="sidebar-$(inst_id)" class="bioprov-sidebar">
        <h3>Node Inspector</h3>
        <p style="color:$(sub_text)">Click on any node in the graph to inspect its parameters, parent lineage, and timestamp.</p>
        <div id="details-$(inst_id)"></div>
      </div>
      <div id="pane-json-$(inst_id)" class="bioprov-text-pane"></div>
      <div id="pane-mermaid-$(inst_id)" class="bioprov-text-pane"></div>
      <div id="pane-methods-$(inst_id)" class="bioprov-text-pane"></div>
    </div>
  </div>

  <script>
    (function() {
      function escapeHtml(str) {
        return String(str)
          .replace(/&/g, "&amp;")
          .replace(/</g, "&lt;")
          .replace(/>/g, "&gt;")
          .replace(/"/g, "&quot;")
          .replace(/'/g, "&#039;");
      }

      const nodesData = $(safe_nodes_json);
      const edgesData = $(safe_edges_json);
      const provJsonData = $(safe_prov_json_val);
      const mermaidData = $(safe_mermaid_val);
      const methodsData = $(safe_methods_val);

      const rootElem = document.getElementById('container-$(inst_id)');
      if (!rootElem) return;

      rootElem.querySelector('#pane-json-$(inst_id)').innerText = provJsonData;
      rootElem.querySelector('#pane-mermaid-$(inst_id)').innerText = mermaidData;
      rootElem.querySelector('#pane-methods-$(inst_id)').innerText = methodsData;

      let network = null;

      function initNetwork() {
        if (!window.vis) return;
        const visNodes = new vis.DataSet(nodesData.map(n => ({
          id: n.id,
          label: n.operation,
          shape: 'box',
          color: {
            background: '$(node_bg)',
            border: '$(accent_color)',
            highlight: { background: '#8b5cf6', border: '#a78bfa' }
          },
          font: { color: '#ffffff', size: 14 }
        })));

        const visEdges = new vis.DataSet(edgesData.map(e => ({
          from: e.from,
          to: e.to,
          arrows: 'to',
          color: { color: '$(sub_text)', highlight: '#8b5cf6' }
        })));

        const graphContainer = rootElem.querySelector('#graph-$(inst_id)');
        network = new vis.Network(graphContainer, { nodes: visNodes, edges: visEdges }, {
          layout: { hierarchical: { direction: 'LR', sortMethod: 'directed' } },
          physics: { enabled: false },
          interaction: { hover: true, dragNodes: true, zoomView: true }
        });

        network.on("selectNode", function(params) {
          const nodeId = params.nodes[0];
          const node = nodesData.find(n => n.id === nodeId);
          if (!node) return;

          let detailsHtml = '<strong>Operation:</strong> <span style="color:$(accent_color)">' + escapeHtml(node.operation) + '</span><br>';
          detailsHtml += '<strong>ID:</strong> <code style="font-size:0.75rem">' + escapeHtml(node.id) + '</code><br>';
          detailsHtml += '<strong>Timestamp:</strong> ' + escapeHtml(node.timestamp) + '<br>';
          detailsHtml += '<strong>Parents:</strong> ' + (node.parent_ids.length > 0 ? node.parent_ids.map(p => '<code>' + escapeHtml(p.substring(0, 8)) + '...</code>').join(', ') : '<em>None (Root Node)</em>') + '<br><br>';
          
          detailsHtml += '<strong>Parameters:</strong>';
          if (Object.keys(node.parameters).length === 0) {
            detailsHtml += '<br><em>None recorded</em>';
          } else {
            detailsHtml += '<table class="bioprov-param-table"><thead><tr><th>Key</th><th>Value</th></tr></thead><tbody>';
            for (const [k, v] of Object.entries(node.parameters)) {
              const valStr = typeof v === 'object' ? JSON.stringify(v) : String(v);
              detailsHtml += '<tr><td>' + escapeHtml(k) + '</td><td>' + escapeHtml(valStr) + '</td></tr>';
            }
            detailsHtml += '</tbody></table>';
          }
          rootElem.querySelector('#details-$(inst_id)').innerHTML = detailsHtml;
        });
      }

      if (typeof vis !== 'undefined') {
        initNetwork();
      } else {
        window.addEventListener('load', initNetwork);
        setTimeout(initNetwork, 500);
      }

      window.bioprovFilterNodes_$(inst_id) = function(query) {
        if (!network) return;
        if (!query.trim()) {
          network.unselectAll();
          return;
        }
        const matches = nodesData.filter(n => n.operation.toLowerCase().includes(query.toLowerCase())).map(n => n.id);
        network.selectNodes(matches);
      };

      window.bioprovSetTab_$(inst_id) = function(tab, btnElem) {
        const btns = rootElem.querySelectorAll('.bioprov-tab-btn');
        btns.forEach(b => b.classList.remove('active'));
        if (btnElem) btnElem.classList.add('active');

        rootElem.querySelector('#graph-$(inst_id)').style.display = tab === 'graph' ? 'block' : 'none';
        rootElem.querySelector('#sidebar-$(inst_id)').style.display = tab === 'graph' ? 'block' : 'none';
        rootElem.querySelector('#pane-json-$(inst_id)').style.display = tab === 'json' ? 'block' : 'none';
        rootElem.querySelector('#pane-mermaid-$(inst_id)').style.display = tab === 'mermaid' ? 'block' : 'none';
        rootElem.querySelector('#pane-methods-$(inst_id)').style.display = tab === 'methods' ? 'block' : 'none';
      };
    })();
  </script>
  $(standalone ? "</body></html>" : "")
  """

  return html
end

function provenance_to_html(ts::ThreadSafeProvenanceContext; kwargs...)
  snapshot = Base.lock(ts.lock) do
    copy(ts.ctx)
  end
  return provenance_to_html(snapshot; kwargs...)
end

function export_provenance_html(ctx::ProvenanceContext, filepath::AbstractString; title::AbstractString="BioToolkit Provenance Report", kwargs...)
  html_content = provenance_to_html(ctx; title=title, standalone=true, kwargs...)
  write(filepath, html_content)
  return String(filepath)
end
function export_provenance_html(ts::ThreadSafeProvenanceContext, filepath::AbstractString; kwargs...)
  snapshot = Base.lock(ts.lock) do
    copy(ts.ctx)
  end
  return export_provenance_html(snapshot, filepath; kwargs...)
end

to_html(ctx::ProvenanceContext) = provenance_to_html(ctx; standalone=true)
to_html(ts::ThreadSafeProvenanceContext) = provenance_to_html(ts; standalone=true)

function Base.show(io::IO, ::MIME"text/html", ctx::ProvenanceContext)
  print(io, provenance_to_html(ctx; standalone=false))
end
function Base.show(io::IO, ::MIME"text/html", ts::ThreadSafeProvenanceContext)
  snapshot = Base.lock(ts.lock) do
    copy(ts.ctx)
  end
  print(io, provenance_to_html(snapshot; standalone=false))
end

export provenance_to_html, export_provenance_html

