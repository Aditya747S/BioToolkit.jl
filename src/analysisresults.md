# BioToolkit.jl - `analysisresults.jl` API Reference

> **Module:** `BioToolkit` · **File:** `analysisresults.jl`

---

## Table of Contents

1. [Purpose and Scope](#purpose-and-scope)
2. [Design Summary](#design-summary)
3. [Analysis Result Protocol](#analysis-result-protocol)
4. [Core Provenance Types](#core-provenance-types)
5. [Global and Scoped Tracking](#global-and-scoped-tracking)
6. [Manual Registration and Result Stamping](#manual-registration-and-result-stamping)
7. [Container Metadata API](#container-metadata-api)
8. [Graph Inspection and Export](#graph-inspection-and-export)
9. [Hashing and Data Integrity](#hashing-and-data-integrity)
10. [Provenance Macros](#provenance-macros)
11. [Methods Section Generation](#methods-section-generation)
12. [Thread Safety](#thread-safety)
13. [Implementation Notes](#implementation-notes)
14. [Quick Reference Card](#quick-reference-card)

---

## Purpose and Scope

`analysisresults.jl` provides two shared facilities used across BioToolkit:

- A common analysis-result protocol for compact display, provenance-aware summaries, and one-row `DataFrame` conversion.
- A provenance engine for recording analysis DAGs, stamping result containers, exporting PROV-JSON, visualizing lineage, computing integrity hashes, and generating draft manuscript methods text.

The file is intentionally broad because nearly every domain module benefits from the same conventions. Sequence analysis, omics workflows, clinical utilities, plotting helpers, and result objects can all participate in the same audit trail without each module reimplementing provenance storage.

---

## Design Summary

| Decision | Rationale |
|---|---|
| **Instrumentation is opt-in** | Provenance is inactive unless the caller supplies or enables a context. Functions can resolve `nothing` and return with minimal overhead. |
| **One context resolver** | Domain functions call `active_provenance_context(...)` rather than inspecting globals or task-local storage directly. |
| **Embedded metadata plus graph nodes** | Results can carry compact `ResultProvenance` metadata while a `ProvenanceContext` stores the operation DAG. |
| **Container-oriented stamping** | DataFrames, dictionaries, and objects with `metadata` or `annotations` dictionaries can be stamped uniformly. |
| **No data copying for large objects** | Provenance parameters are JSON-normalized with truncation, and large data are represented by structural or content hashes. |
| **DataFrame provenance stays in metadata** | `DataFrame(result::AbstractAnalysisResult)` attaches provenance through `DataAPI.metadata!`, avoiding repeated columns in large tables. |
| **Thread safety is explicit** | `ProvenanceContext` is lightweight for single-threaded work; `ThreadSafeProvenanceContext` wraps mutations with a lock. |
| **Limit violations throw** | `max_nodes` and `max_depth` prevent silent lineage loss. When limits are hit, registration raises an `ArgumentError`. |

---

## Analysis Result Protocol

### `AbstractAnalysisResult`

```julia
abstract type AbstractAnalysisResult end
```

Common supertype for domain-specific result objects. Subtyping it gives a result:

- `summary(result)` and `show(result)` based on `analysis_result_summary`.
- `DataFrame(result)` conversion into a one-row table.
- Automatic provenance inference via `analysis_result_provenance`.

Domain modules normally define concrete result structs like this:

```julia
struct MyResult <: AbstractAnalysisResult
    statistic::Float64
    pvalue::Float64
    method::Symbol
    provenance::ResultProvenance
end
```

If a result has a field named `provenance` containing a `ResultProvenance`, that record is used directly. Otherwise BioToolkit builds a compact inferred record from recognized fields such as `method`, `reference_name`, `phenotype_name`, `system`, `target_gene`, `term_id`, `db`, `modalities`, `converged`, `zero_inflated`, and `passes_filters`.

### `analysis_result_type`

```julia
analysis_result_type(result)
analysis_result_type(::Type) = :generic
```

Returns a symbolic category used in summary and provenance labels. The default is `:generic`.

Override this when a result type should display a domain-specific category:

```julia
BioToolkit.analysis_result_type(::Type{MyResult}) = :differential_expression
```

### `analysis_result_fields`

```julia
analysis_result_fields(::Type{T}) where T
```

Returns the fields included in summaries and fallback `DataFrame` conversion. By default, every field except `:provenance` is included.

Override this when a result contains large arrays or implementation details:

```julia
BioToolkit.analysis_result_fields(::Type{MyResult}) = (:statistic, :pvalue, :method)
```

### `analysis_result_provenance`

```julia
analysis_result_provenance(result) -> ResultProvenance
```

Returns the result's provenance record. The function first checks for an explicit `provenance` field. If absent, it synthesizes a `ResultProvenance` from the result type, module, recognized fields, warning flags, and optional extension hooks.

Extension hooks:

```julia
analysis_result_provenance_label_parts(result) = String[]
analysis_result_provenance_parameters(result) = NamedTuple()
```

Define these for a result type when standard fields are not enough to build a meaningful provenance label.

### `analysis_result_summary`

```julia
analysis_result_summary(result) -> String
```

Builds a compact one-line summary containing the inferred label, provenance id, source, status, selected fields, warnings, fallbacks, notes, and parameters.

```julia
summary(result)
show(result)
```

For `AbstractAnalysisResult`, both methods delegate to `analysis_result_summary`.

### `DataFrame(result::AbstractAnalysisResult)`

```julia
DataFrames.DataFrame(result::AbstractAnalysisResult) -> DataFrame
```

Converts the result to a one-row table using `analysis_result_fields`. Provenance is attached as DataFrame metadata under the `"provenance"` key rather than repeated as columns.

```julia
using BioToolkit, DataFrames, DataAPI

df = DataFrame(result)
prov = DataAPI.metadata(df, "provenance")
```

Use `flatten_provenance!(df)` when provenance identifiers must become ordinary columns for CSV-style export.

---

## Core Provenance Types

### `ResultProvenance`

```julia
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
```

Compact provenance metadata for one result or table. It is designed for display and metadata storage, not as the full graph. The graph lives in `ProvenanceContext`.

Construct records with `provenance_record` rather than calling the constructor directly.

### `ProvenanceNode`

```julia
struct ProvenanceNode
    id::String
    operation::String
    parameters::Dict{String,Any}
    parent_ids::Vector{String}
    timestamp::String
end
```

One node in the provenance DAG. A node records an operation, JSON-compatible parameter payload, parent node identifiers, and timestamp.

### `EnvironmentSnapshot`

```julia
struct EnvironmentSnapshot
    julia_version::String
    manifest_hash::Union{String,Nothing}
    package_versions::Dict{String,String}
    cpu_threads::Int
    git_commit::Union{String,Nothing}
    random_state_hash::String
    timestamp::String
end
```

Captures the computational environment when a context is created or reset. The snapshot includes Julia version, active manifest hash when available, package versions, thread count, git commit when available, default RNG state hash, and timestamp.

### `ProvenanceContext`

```julia
mutable struct ProvenanceContext
    nodes::OrderedDict{String,ProvenanceNode}
    environment::EnvironmentSnapshot
    max_nodes::Union{Int,Nothing}
    max_depth::Union{Int,Nothing}
end
```

Stores the provenance DAG as insertion-ordered nodes.

Constructors:

```julia
ProvenanceContext(; max_nodes=nothing, max_depth=nothing)
ProvenanceContext(nodes::OrderedDict{String,ProvenanceNode})
```

The optional limits are audit safeguards:

- `max_nodes`: maximum number of nodes allowed.
- `max_depth`: maximum DAG depth allowed.

When a new node would exceed either limit, registration throws rather than dropping lineage.

### `ThreadSafeProvenanceContext`

```julia
mutable struct ThreadSafeProvenanceContext
    ctx::ProvenanceContext
    lock::ReentrantLock
end

ThreadSafeProvenanceContext()
```

Thread-safe wrapper around `ProvenanceContext`. It serializes mutation and graph inspection through a `ReentrantLock`. Use this when multiple Julia tasks or threads write to a shared context.

---

## Global and Scoped Tracking

### `enable_provenance!`

```julia
enable_provenance!() -> ProvenanceContext
enable_provenance!(ctx) -> ctx
```

Enables global automatic provenance tracking. When enabled, BioToolkit functions that call `active_provenance_context()` will record into the global context.

```julia
using BioToolkit

ctx = enable_provenance!()

# BioToolkit functions that support provenance will record into ctx.
# result = some_biotoolkit_function(input)

json = export_provenance_json(ctx)
disable_provenance!()
```

Pass an existing context to resume or share tracking:

```julia
ctx = ProvenanceContext()
enable_provenance!(ctx)
```

For threaded workloads, pass a `ThreadSafeProvenanceContext`.

### `disable_provenance!`

```julia
disable_provenance!() -> nothing
```

Disables global provenance tracking and clears the global context reference.

### `get_provenance_context`

```julia
get_provenance_context() -> Union{Nothing, ProvenanceContext, ThreadSafeProvenanceContext}
```

Returns the active global context, or `nothing` if provenance is disabled.

### `with_global_provenance`

```julia
with_global_provenance(f)
with_global_provenance(f, ctx)
```

Runs `f()` with a temporary global context and restores the previous global state afterward.

```julia
result = with_global_provenance() do
    # pipeline steps
end
```

### `with_provenance` for scoped tracking

```julia
with_provenance(f::Function, ctx)
with_provenance(f::Function)
```

Pushes a task-local scoped context for the duration of `f`. This is the preferred per-analysis API because it avoids global state while still allowing nested BioToolkit calls to find a context.

```julia
ctx = ProvenanceContext()

result = with_provenance(ctx) do
    # nested BioToolkit functions record into ctx
end
```

Calling `with_provenance(f)` creates a fresh context internally.

### `active_provenance_context`

```julia
active_provenance_context([explicit])
```

Resolves the context that a public BioToolkit function should use:

1. Return `explicit` when it is not `nothing`.
2. Otherwise use the task-local scoped context.
3. Otherwise use the enabled global context.
4. Otherwise return `nothing`.

BioToolkit implementers should call this once near the public boundary:

```julia
function my_operation(x; prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    result = compute_result(x)
    return provenance_result!(_ctx, result, "my_operation";
        parents=provenance_parent_ids(x),
        parameters=(n=length(x),))
end
```

---

## Manual Registration and Result Stamping

### `new_provenance_id`

```julia
new_provenance_id() -> String
new_provenance_id(content::AbstractString) -> String
new_provenance_id(content::AbstractVector{UInt8}) -> String
```

Generates SHA-256-style provenance identifiers.

- Without content, the id is non-deterministic and suitable for runtime operations.
- With content, the id is deterministic and suitable for pure functions.

### `register_provenance!`

```julia
register_provenance!(ctx, node::ProvenanceNode)
register_provenance!(ctx, operation; parents=String[], parameters=NamedTuple())
register_provenance!(ctx, node_id, operation; parents=String[], parameters=NamedTuple())
```

Adds a node to a `ProvenanceContext` or `ThreadSafeProvenanceContext`.

```julia
ctx = ProvenanceContext()

node = register_provenance!(ctx, "load_counts";
    parameters=(path="counts.csv",))

register_provenance!(ctx, "normalize_counts";
    parents=[node.id],
    parameters=(method=:median_ratio,))
```

### `provenance_result!`

```julia
provenance_result!(ctx, result, operation; kwargs...) -> result
provenance_result!(nothing, result, operation; kwargs...) -> result
```

Registers a graph node and stamps the result container when possible. Passing `nothing` is a no-op.

Important keywords:

| Keyword | Meaning |
|---|---|
| `parents` | Parent provenance node ids. |
| `parameters` | Operation parameters, accepted as `NamedTuple`, `Dict`, pairs, or compatible payloads. |
| `id` | Optional explicit node/container id. |
| `label` | Human-facing result label. |
| `source` | Source operation/module label. |
| `status` | `:ok`, `:warn`, `:error`, or other symbolic status. |
| `warnings`, `errors`, `fallbacks`, `notes` | Audit messages stored in `ResultProvenance`. |
| `provenance_hash` | Optional attached hash string. |

```julia
df = DataFrame(gene=["A", "B"], score=[1.0, 2.0])
ctx = ProvenanceContext()

provenance_result!(ctx, df, "score_genes";
    parameters=(method=:zscore,),
    label="gene scores",
    source="BioToolkit/score_genes")
```

### `with_provenance` for existing results

```julia
with_provenance(result, label, source; kwargs...) -> result
```

Stamps an existing result, table, dictionary, or supported container with `ResultProvenance`. If a scoped or global context is active, it can also register a graph node.

```julia
metadata = Dict{Symbol,Any}()

with_provenance(metadata, "sample metadata", "read_samples";
    parameters=(n_samples=12,))
```

### `provenance_record`

```julia
provenance_record(label, source; kwargs...) -> ResultProvenance
provenance_record(provenance::ResultProvenance; kwargs...) -> ResultProvenance
provenance_record(payload::Union{AbstractDict,NamedTuple}) -> ResultProvenance
```

Normalizes user-facing provenance payloads into `ResultProvenance`.

```julia
prov = provenance_record("PCA result", "run_pca";
    status=:ok,
    parameters=(n_components=50,))
```

### `provenance_summary`

```julia
provenance_summary(provenance) -> String
```

Returns a compact display string for a `ResultProvenance`, dictionary payload, named tuple payload, or `nothing`.

### `ProvenanceParams`

```julia
const ProvenanceParams = Union{NamedTuple, Dict{String,Any}}
```

Type alias used by provenance APIs for parameter payloads. Prefer `NamedTuple` in hot paths and `Dict{String,Any}` for dynamic pipelines.

---

## Container Metadata API

### Supported containers

Provenance can be stored in:

- `DataFrames.AbstractDataFrame` metadata via `DataAPI.metadata!`.
- Dictionaries with compatible key/value types.
- Objects with a `metadata` dictionary property.
- Objects with an `annotations` dictionary property.
- Vectors of supported containers, where each supported element can be stamped individually.

### `metadata_provenance`

```julia
metadata_provenance(metadata::AbstractDict)
metadata_provenance(table::DataFrames.AbstractDataFrame)
```

Returns embedded `ResultProvenance` if present, otherwise `nothing`.

### `stamp_provenance!`

```julia
stamp_provenance!(metadata::AbstractDict, provenance::ResultProvenance)
stamp_provenance!(table::DataFrame, provenance::ResultProvenance)
stamp_provenance!(metadata_or_table; kwargs...)
```

Stores provenance in metadata. If the record has no id, an id is generated and recorded.

### `update_provenance!`

```julia
update_provenance!(metadata_or_table; kwargs...)
```

Updates existing embedded provenance, appending warnings, errors, fallbacks, notes, and parent ids while preserving current fields when replacement keywords are omitted.

### `container_provenance_id`

```julia
container_provenance_id(container) -> Union{String,Nothing}
```

Returns the stored provenance id for a supported container.

### `ensure_provenance_id!`

```julia
ensure_provenance_id!(metadata::AbstractDict) -> String
ensure_provenance_id!(table::DataFrame) -> String
```

Returns an existing id or creates one.

### `provenance_parent_ids`

```julia
provenance_parent_ids(values...) -> Vector{String}
```

Extracts unique provenance ids from containers, vectors, tuples, and sets. This is the preferred way to link outputs to input containers.

```julia
parents = provenance_parent_ids(raw_counts, sample_metadata)
```

### `register_container_provenance!`

```julia
register_container_provenance!(ctx, container, operation; parents=String[], parameters=NamedTuple(), provenance_hash=nothing)
```

Ensures a container id, optionally stores a hash, and registers a graph node with the same id.

### `container_provenance_summary`

```julia
container_provenance_summary(container) -> String
```

Returns a compact `provenance=id=..., hash=...` summary or `provenance=untracked`.

### `container_provenance_hash`

```julia
container_provenance_hash(container) -> Union{String,Nothing}
```

Returns a stored provenance hash if present.

### `has_provenance` and `requires_provenance`

```julia
has_provenance(result) -> Bool
requires_provenance(result) -> result
```

`has_provenance` checks whether a result carries embedded provenance metadata. `requires_provenance` returns the result when provenance exists and throws `ArgumentError` otherwise.

### `stamp_provenance_with_hash!`

```julia
stamp_provenance_with_hash!(container, data; hash_mode=:structural, kwargs...) -> container
```

Stamps provenance and stores a hash of `data` under the provenance hash metadata key.

Hash modes:

- `:structural`: O(1) shape/type hash for arrays and DataFrames. Safe for large data but not value-level integrity.
- `:content`: O(n) true SHA-256 over content. Use for small to medium data when value integrity matters.

```julia
stamp_provenance_with_hash!(result_df, input_matrix;
    label="normalized counts",
    source="normalize_counts",
    hash_mode=:structural)
```

### `flatten_provenance!`

```julia
flatten_provenance!(df::DataFrame) -> df
```

Copies DataFrame-level `provenance_id` and `provenance_hash` metadata into ordinary columns. This is useful before writing CSV or any format that does not preserve DataAPI metadata.

---

## Graph Inspection and Export

### `find_nodes`

```julia
find_nodes(ctx, pattern) -> Vector{ProvenanceNode}
```

Finds nodes whose operation matches a regex-compatible pattern or exact string.

```julia
find_nodes(ctx, "normalize")
```

### `provenance_chain`

```julia
provenance_chain(ctx, node_id) -> Vector{ProvenanceNode}
```

Returns the ancestor chain for `node_id`, oldest first. The queried node is included when present.

### `provenance_ancestors`

```julia
provenance_ancestors(ctx, node_id) -> Vector{ProvenanceNode}
```

Returns all upstream nodes excluding the queried node.

### `provenance_descendants`

```julia
provenance_descendants(ctx, node_id) -> Vector{ProvenanceNode}
```

Returns downstream nodes that depend on `node_id`.

### `provenance_lineage_table`

```julia
provenance_lineage_table(ctx) -> DataFrame
```

Builds a one-row-per-node table with columns:

- `id`
- `operation`
- `timestamp`
- `n_parents`
- `parents`
- `parameters`

### `provenance_diff`

```julia
provenance_diff(ctx_a, ctx_b)
```

Compares two contexts by node id and returns:

```julia
(only_in_a=..., only_in_b=..., common=...)
```

### `merge_provenance_contexts`

```julia
merge_provenance_contexts(contexts...; max_nodes=nothing, max_depth=nothing)
```

Combines nodes from multiple contexts into a fresh `ProvenanceContext`. Limit checks apply during merge.

### `export_provenance_json`

```julia
export_provenance_json(ctx) -> String
```

Serializes the context to a PROV-JSON-like payload with:

- schema version
- PROV namespace prefix
- environment snapshot
- entities
- activities
- generation relationships
- derivation relationships

### `import_provenance_json`

```julia
import_provenance_json(json_text) -> ProvenanceContext
import_provenance_json(ThreadSafeProvenanceContext, json_text) -> ThreadSafeProvenanceContext
```

Restores a context from exported JSON. The schema version must match the current provenance schema.

### `provenance_to_dot`

```julia
provenance_to_dot(ctx) -> String
```

Renders the provenance DAG as Graphviz DOT.

### `provenance_to_mermaid`

```julia
provenance_to_mermaid(ctx) -> String
```

Renders the provenance DAG as a Mermaid `flowchart LR`.

### `Base.empty!`

```julia
empty!(ctx::ProvenanceContext) -> ctx
empty!(ctx::ThreadSafeProvenanceContext) -> ctx
```

Clears all nodes and refreshes the environment snapshot.

### `Base.copy`

```julia
copy(ctx::ProvenanceContext) -> ProvenanceContext
```

Creates a snapshot copy of a context, including copied nodes, parameters, parent ids, environment metadata, and limits.

### Iteration and length

```julia
iterate(ctx::ProvenanceContext)
length(ctx::ProvenanceContext)
```

Iterating a context delegates to its ordered node dictionary.

---

## Hashing and Data Integrity

### `provenance_structural_hash`

```julia
provenance_structural_hash(data) -> String
```

Computes a fast hash for structural identity.

Behavior:

- Strings: hashes full string content.
- Arrays: hashes shape, element type, and length only.
- DataFrames: hashes table size and column names.
- Dictionaries: hashes length, key type, and value type.
- Numbers and other objects: hashes a compact representation.

Use for large data where a cheap audit fingerprint is enough.

### `provenance_content_hash`

```julia
provenance_content_hash(data) -> String
```

Computes a true SHA-256 content hash by iterating through values. This is appropriate for integrity verification, but it can be expensive for large arrays and DataFrames.

### `file_provenance_hash`

```julia
file_provenance_hash(path::AbstractString) -> String
```

Streams a file in 64 KiB chunks and computes SHA-256 without reading the full file into memory. This is the preferred way to fingerprint raw FASTA, VCF, BAM, HDF5, or count-matrix input files.

### `set_provenance_limits!`

```julia
set_provenance_limits!(; max_array_len::Int, max_string_len::Int)
```

Configures JSON-normalization truncation limits for arrays and strings stored inside provenance parameter payloads.

---

## Provenance Macros

### `@provenance`

```julia
@provenance ctx operation parameters expr
@provenance ctx operation parents parameters expr
```

Runs `expr`, stamps the result, and records a graph node when `ctx` is a `ProvenanceContext` or `ThreadSafeProvenanceContext`. When `ctx === nothing`, the macro expands to a no-op path that simply evaluates the expression.

```julia
result = @provenance ctx "normalize_counts" (method=:median_ratio,) begin
    normalize_counts(counts)
end
```

Use the five-argument form to pass explicit parent ids:

```julia
result = @provenance ctx "normalize_counts" parent_ids (method=:median_ratio,) begin
    normalize_counts(counts)
end
```

Auto-parent linking uses task-local storage. It works for ordinary nested calls and threaded loops, but `@spawn` workflows should pass parents explicitly.

### `@pure_provenance`

```julia
@pure_provenance ctx operation content_key parameters expr
```

Like `@provenance`, but uses a deterministic id from `new_provenance_id(String(content_key))`. Use this for deterministic pure functions where identical inputs should produce identical provenance node ids.

```julia
result = @pure_provenance ctx "scale_matrix" string(size(mat), ":zscore") (method=:zscore,) begin
    scale_matrix(mat)
end
```

### `lazy_provenance_id!`

```julia
lazy_provenance_id!(ctx, result, operation; parents=String[], parameters=NamedTuple())
```

Like `provenance_result!`, but explicitly skips work when `ctx === nothing`. This is intended for very hot paths.

---

## Methods Section Generation

### `generate_methods_section`

```julia
generate_methods_section(ctx::ProvenanceContext; title="## Methods", max_ops=nothing) -> String
generate_methods_section(ctx::ThreadSafeProvenanceContext; title="## Methods", max_ops=nothing) -> String
generate_methods_section(result; title="## Methods", max_ops=nothing) -> String
```

Generates draft methods prose from a provenance graph.

The function:

1. Topologically orders nodes from roots to leaves.
2. Converts operation names such as `normalize_counts` into readable labels.
3. Renders parameter dictionaries as compact parenthetical text.
4. Adds a closing software/environment sentence from `EnvironmentSnapshot`.

```julia
ctx = ProvenanceContext()

raw = register_provenance!(ctx, "load_counts"; parameters=(path="counts.csv",))
norm = register_provenance!(ctx, "normalize_counts";
    parents=[raw.id],
    parameters=(method=:median_ratio,))

println(generate_methods_section(ctx))
```

Output shape:

```markdown
## Methods

Data were obtained using **load counts** (path: counts.csv).
This was followed by **normalize counts** (method: median_ratio).

All analyses were performed using BioToolkit.jl; Julia v...
```

The generated text is a draft. Review it before including it in a manuscript.

---

## Thread Safety

`ProvenanceContext` is not thread-safe for concurrent writes. It is the right default for REPL, scripts, notebooks, and single-threaded pipelines.

Use `ThreadSafeProvenanceContext` when:

- multiple threads register nodes into one shared context;
- a global provenance context is enabled during threaded work;
- tasks share one context and may mutate it concurrently.

Thread-safe overloads exist for registration, result stamping, lineage queries, diffing, merging, export, graph rendering, ancestor/descendant traversal, and methods generation.

When locking two thread-safe contexts for `provenance_diff`, the implementation orders locks by `objectid` to avoid deadlock.

---

## Implementation Notes

### Parameter serialization

Parameter payloads are normalized to JSON-compatible values:

- dictionaries become string-keyed dictionaries;
- named tuples become dictionaries;
- arrays are recursively serialized until length limits are exceeded;
- large arrays are replaced by a compact truncated summary and hash;
- long strings are truncated with a short hash suffix;
- symbols become strings;
- unrecognized objects fall back to `string(value)`.

### Placeholder pruning

The registration layer prunes generic placeholder self-registrations when a more specific node for the same operation appears. This keeps graphs readable when early placeholder nodes are later replaced by result-specific provenance.

### Metadata keys

Core metadata keys:

| Constant | Stored key |
|---|---|
| `PROVENANCE_METADATA_KEY` | `"provenance"` |
| `PROVENANCE_ID_KEY` | `:provenance_id` or `"provenance_id"` |
| `PROVENANCE_HASH_KEY` | `:provenance_hash` or `"provenance_hash"` |

The implementation adapts to symbol-keyed and string-keyed metadata dictionaries.

### Scientific software cautions

- Provenance must not change numerical results.
- Resolve the context once at the public boundary and pass it inward.
- Prefer compact parameter summaries over storing full datasets.
- Use `file_provenance_hash` or `provenance_content_hash` for integrity claims.
- Treat `provenance_structural_hash` as a shape/type check, not a value check.
- Use `ThreadSafeProvenanceContext` for concurrent writers.

---

## Quick Reference Card

### Result protocol

| API | Purpose |
|---|---|
| `AbstractAnalysisResult` | Supertype for result objects with shared display/export behavior. |
| `analysis_result_type` | Return symbolic result category. |
| `analysis_result_fields` | Select fields shown/exported by fallback methods. |
| `analysis_result_provenance` | Return explicit or inferred `ResultProvenance`. |
| `analysis_result_summary` | Build compact display string. |
| `DataFrame(result)` | One-row table with provenance in metadata. |

### Context lifecycle

| API | Purpose |
|---|---|
| `ProvenanceContext()` | Create a provenance DAG context. |
| `ThreadSafeProvenanceContext()` | Create a locked context for concurrent writers. |
| `enable_provenance!()` | Enable global automatic provenance. |
| `disable_provenance!()` | Disable global provenance. |
| `get_provenance_context()` | Return the current global context. |
| `with_global_provenance` | Temporarily replace global context. |
| `with_provenance(f, ctx)` | Use task-local scoped provenance. |
| `active_provenance_context` | Resolve explicit, scoped, or global context. |

### Stamping and registration

| API | Purpose |
|---|---|
| `new_provenance_id` | Create random or deterministic ids. |
| `register_provenance!` | Add a node to a context. |
| `provenance_result!` | Stamp a result and register a node. |
| `register_container_provenance!` | Ensure container id and register node. |
| `stamp_provenance!` | Store embedded provenance. |
| `update_provenance!` | Merge provenance updates. |
| `with_provenance(result, label, source)` | Stamp an existing result. |
| `stamp_provenance_with_hash!` | Stamp provenance plus hash. |
| `flatten_provenance!` | Copy DataFrame provenance metadata into columns. |

### Metadata queries

| API | Purpose |
|---|---|
| `metadata_provenance` | Read embedded `ResultProvenance`. |
| `container_provenance_id` | Read container id. |
| `container_provenance_hash` | Read container hash. |
| `container_provenance_summary` | Display provenance id/hash summary. |
| `ensure_provenance_id!` | Create or return metadata id. |
| `provenance_parent_ids` | Extract parent ids from inputs. |
| `has_provenance` | Check embedded provenance. |
| `requires_provenance` | Assert embedded provenance exists. |

### Graph operations

| API | Purpose |
|---|---|
| `find_nodes` | Search by operation. |
| `provenance_chain` | Return ancestor chain including queried node. |
| `provenance_ancestors` | Return upstream nodes excluding queried node. |
| `provenance_descendants` | Return downstream nodes. |
| `provenance_lineage_table` | Tabular graph summary. |
| `provenance_diff` | Compare contexts. |
| `merge_provenance_contexts` | Merge contexts. |
| `export_provenance_json` | Serialize graph. |
| `import_provenance_json` | Restore graph. |
| `provenance_to_dot` | Graphviz output. |
| `provenance_to_mermaid` | Mermaid output. |
| `generate_methods_section` | Draft methods text. |

### Hashing

| API | Purpose |
|---|---|
| `provenance_structural_hash` | Fast shape/type hash. |
| `provenance_content_hash` | Full value-level SHA-256. |
| `file_provenance_hash` | Streaming file SHA-256. |
| `set_provenance_limits!` | Configure parameter truncation limits. |

### Macros

| API | Purpose |
|---|---|
| `@provenance` | Runtime provenance wrapper with random id. |
| `@pure_provenance` | Deterministic provenance wrapper with content-derived id. |
| `lazy_provenance_id!` | Hot-path result stamping helper. |

---

## Detailed API Reference

This section records the exported API in a per-symbol format: kind, signature, purpose, behavior, and examples. The earlier sections give the workflow view; this section is intended for source-adjacent reference documentation.

---

## 15. Analysis Result Protocol - Detailed Reference

### `AbstractAnalysisResult`

```julia
abstract type AbstractAnalysisResult end
```

**Kind:** Abstract supertype

**Description:** Common supertype for BioToolkit domain result objects. Concrete subtypes receive shared summary, display, provenance inference, and fallback `DataFrame` conversion behavior.

**Provided behavior:**

| Method | Behavior |
|---|---|
| `summary(result)` | Delegates to `analysis_result_summary(result)`. |
| `show(io, result)` | Prints the compact summary. |
| `show(io, MIME"text/plain"(), result)` | Prints summary plus provenance warnings/fallbacks/errors when present. |
| `DataFrame(result)` | Converts selected fields to a one-row table and stores provenance in metadata. |

**Example:**

```julia
struct MyResult <: AbstractAnalysisResult
    statistic::Float64
    pvalue::Float64
    method::Symbol
    provenance::ResultProvenance
end
```

---

### `analysis_result_fields(::Type{T})`

```julia
analysis_result_fields(::Type{T}) where {T}
```

**Kind:** Generic function

**Description:** Selects result fields included in summaries and fallback tabular export. The default includes every field except `:provenance`.

**Returns:** `Tuple` of field names.

**Default behavior:**

```julia
Tuple(name for name in fieldnames(T) if name != :provenance)
```

**Override example:**

```julia
BioToolkit.analysis_result_fields(::Type{MyLargeResult}) =
    (:method, :n_features, :score)
```

---

### `analysis_result_type(result)` / `analysis_result_type(::Type)`

```julia
analysis_result_type(result)
analysis_result_type(::Type) = :generic
```

**Kind:** Generic function

**Description:** Returns a symbolic result category used in inferred labels and sources.

**Returns:** `Symbol`

**Override example:**

```julia
BioToolkit.analysis_result_type(::Type{DifferentialExpressionResult}) =
    :differential_expression
```

---

### `analysis_result_provenance(result)`

```julia
analysis_result_provenance(result) -> ResultProvenance
```

**Kind:** Generic function

**Description:** Returns a result provenance record. Explicit `result.provenance::ResultProvenance` is used directly. Dictionary or named-tuple provenance payloads are normalized with `provenance_record`. If no explicit provenance exists, the function infers a record from result type, parent module, recognized fields, and warning/status flags.

**Recognized inferred fields:**

| Field | Effect |
|---|---|
| `method` | Added to label and parameters. |
| `reference_name` | Added to label and parameters. |
| `phenotype_name` | Added to label and parameters. |
| `system` | Added to label and parameters. |
| `target_gene` | Added as `target=...`. |
| `term_id` | Added as `term=...`. |
| `db` | Added as `db=...`. |
| `modalities` | Adds `multimodal` label and parameter. |
| `converged=false` | Adds warning and fallback. |
| `zero_inflated=true` | Adds note. |
| `passes_filters=false` | Adds warning. |

**Extension hooks:**

```julia
analysis_result_provenance_label_parts(result) = String[]
analysis_result_provenance_parameters(result) = NamedTuple()
```

---

### `analysis_result_summary(result)`

```julia
analysis_result_summary(result) -> String
```

**Kind:** Generic function

**Description:** Builds a compact one-line result summary from selected fields and provenance. Arrays and dictionaries are represented by `summary(value)` instead of full contents.

**Returns:** `String`

**Includes:** label, id, source, status, selected fields, warnings, fallbacks, notes, and parameters.

---

### `DataFrames.DataFrame(result::AbstractAnalysisResult)`

```julia
DataFrames.DataFrame(result::AbstractAnalysisResult) -> DataFrames.DataFrame
```

**Kind:** Conversion method

**Description:** Converts a result object into a one-row `DataFrame` using `analysis_result_fields`. Provenance is attached via `DataAPI.metadata!(table, "provenance", provenance; style=:note)` rather than repeated as columns.

**Why metadata instead of columns:** Large genomics result tables can have millions of rows. Repeating provenance values as columns wastes memory and makes provenance look like per-row data when it is table-level metadata.

---

## 16. Provenance Type Reference

### `ResultProvenance`

```julia
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
```

**Kind:** Immutable concrete struct

**Description:** Compact provenance record intended for embedding in result objects or metadata dictionaries. It summarizes a single result's origin and status. It is not the full lineage graph; full graph lineage is stored in `ProvenanceContext`.

**Fields:**

| Field | Type | Description |
|---|---|---|
| `id` | `String` | Provenance id. Empty string means untracked until stamped. |
| `label` | `String` | Human-facing result label. |
| `source` | `String` | Source operation, module, or data origin. |
| `status` | `Symbol` | Typical values are `:ok`, `:warn`, or `:error`. |
| `warnings` | `Vector{String}` | Non-fatal warnings. |
| `errors` | `Vector{String}` | Error messages associated with the result. |
| `fallbacks` | `Vector{String}` | Fallback algorithms or degraded paths used. |
| `notes` | `Vector{String}` | Additional audit notes. |
| `parameters` | `NamedTuple` | Compact operation parameters. |
| `parent_ids` | `Vector{String}` | Upstream provenance ids. |
| `timestamp` | `String` | UTC timestamp in ISO-like format. |

**Constructor recommendation:** Use `provenance_record(...)` rather than calling this constructor directly.

---

### `ProvenanceNode`

```julia
struct ProvenanceNode
    id::String
    operation::String
    parameters::Dict{String,Any}
    parent_ids::Vector{String}
    timestamp::String
end
```

**Kind:** Immutable concrete struct

**Description:** One operation node in a provenance DAG. Nodes are keyed by `id` inside a `ProvenanceContext`.

**Fields:**

| Field | Type | Description |
|---|---|---|
| `id` | `String` | Node id. |
| `operation` | `String` | Operation label such as `read_fasta` or `normalize_counts`. |
| `parameters` | `Dict{String,Any}` | JSON-compatible parameter summary. |
| `parent_ids` | `Vector{String}` | Upstream nodes used to create this node. |
| `timestamp` | `String` | Registration time. |

**Display:** `summary(node)` prints `ProvenanceNode(operation, parents=n)`.

---

### `EnvironmentSnapshot`

```julia
struct EnvironmentSnapshot
    julia_version::String
    manifest_hash::Union{String,Nothing}
    package_versions::Dict{String,String}
    cpu_threads::Int
    git_commit::Union{String,Nothing}
    random_state_hash::String
    timestamp::String
end
```

**Kind:** Immutable concrete struct

**Description:** Captures computational environment metadata at context creation/reset time.

**Fields:**

| Field | Description |
|---|---|
| `julia_version` | Active Julia version. |
| `manifest_hash` | SHA-256 hash of active `Manifest.toml`, if available. |
| `package_versions` | Installed dependency versions from `Pkg.dependencies()`. |
| `cpu_threads` | `Threads.nthreads()`. |
| `git_commit` | Current repository commit if discoverable. |
| `random_state_hash` | Hash of the default RNG representation. |
| `timestamp` | UTC snapshot timestamp. |

---

### `ProvenanceContext`

```julia
mutable struct ProvenanceContext
    nodes::OrderedDict{String,ProvenanceNode}
    environment::EnvironmentSnapshot
    max_nodes::Union{Int,Nothing}
    max_depth::Union{Int,Nothing}
end
```

**Kind:** Mutable concrete struct

**Description:** Stores a provenance DAG as ordered nodes plus environment metadata. This is the main graph object for single-threaded scripts, notebooks, and pipelines.

**Constructors:**

```julia
ProvenanceContext(; max_nodes=nothing, max_depth=nothing)
ProvenanceContext(nodes::OrderedDict{String,ProvenanceNode})
```

**Limit semantics:**

| Limit | Behavior |
|---|---|
| `max_nodes` | Registration throws when adding a new node would exceed the limit. |
| `max_depth` | Registration throws when the new node lineage would exceed depth. |

**Display:** `ProvenanceContext(n nodes)`.

---

### `ThreadSafeProvenanceContext`

```julia
mutable struct ThreadSafeProvenanceContext
    ctx::ProvenanceContext
    lock::ReentrantLock
end

ThreadSafeProvenanceContext()
```

**Kind:** Mutable wrapper type

**Description:** Lock-protected wrapper for shared provenance contexts in threaded code. Mutating operations acquire `lock` and delegate to the inner `ctx`.

**Use when:** multiple threads or tasks write to the same context.

**Do not assume:** `ProvenanceContext` itself is safe for concurrent writes.

---

## 17. Context Lifecycle API

### `enable_provenance!`

```julia
enable_provenance!() -> ProvenanceContext
enable_provenance!(ctx) -> ctx
```

**Kind:** Mutating global-state function

**Description:** Enables global provenance tracking. Functions that resolve context with `active_provenance_context()` will record into the global context when no explicit or scoped context is present.

**Returns:** The enabled context.

**Example:**

```julia
ctx = enable_provenance!()
# run BioToolkit workflow
graph_json = export_provenance_json(ctx)
disable_provenance!()
```

**Threading note:** Use `ThreadSafeProvenanceContext()` for global tracking in concurrent workloads.

---

### `disable_provenance!`

```julia
disable_provenance!() -> nothing
```

**Kind:** Mutating global-state function

**Description:** Clears the global context and disables global provenance recording.

---

### `get_provenance_context`

```julia
get_provenance_context()
```

**Kind:** Global-state query

**Description:** Returns the currently configured global provenance context, or `nothing` when disabled.

---

### `with_global_provenance`

```julia
with_global_provenance(f)
with_global_provenance(f, ctx)
```

**Kind:** Scoped global-state helper

**Description:** Temporarily replaces the global provenance context for the dynamic extent of `f()`, then restores the previous state in a `finally` block.

**Example:**

```julia
result = with_global_provenance() do
    # pipeline using temporary global context
end
```

---

### `with_provenance(f::Function, ctx)`

```julia
with_provenance(f::Function, ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext})
with_provenance(f::Function)
```

**Kind:** Task-local scoped helper

**Description:** Pushes a task-local provenance context while `f()` runs, then restores the previous stack state. This is the preferred scoped API for per-analysis tracking.

**Example:**

```julia
ctx = ProvenanceContext()
result = with_provenance(ctx) do
    # nested BioToolkit calls record into ctx
end
```

---

### `active_provenance_context`

```julia
active_provenance_context([explicit])
```

**Kind:** Context resolver

**Description:** Resolves which context a BioToolkit public function should use.

**Resolution order:**

1. Explicit non-`nothing` argument.
2. Task-local context from `with_provenance(ctx) do ... end`.
3. Enabled global context.
4. `nothing`.

**Implementer pattern:**

```julia
function public_api(x; prov_ctx=nothing)
    _ctx = active_provenance_context(prov_ctx)
    result = compute(x)
    return provenance_result!(_ctx, result, "public_api";
        parents=provenance_parent_ids(x),
        parameters=(n=length(x),))
end
```

---

## 18. Registration and Stamping API

### `new_provenance_id`

```julia
new_provenance_id() -> String
new_provenance_id(content::AbstractString) -> String
new_provenance_id(content::AbstractVector{UInt8}) -> String
```

**Kind:** Identifier generator

**Description:** Generates SHA-256-style provenance ids. The zero-argument form includes UUID/time/randomness and is non-deterministic. Content forms are deterministic.

**Use cases:**

| Form | Use case |
|---|---|
| `new_provenance_id()` | Runtime I/O, random models, non-deterministic operations. |
| `new_provenance_id(content)` | Pure deterministic operations and content-addressable lineage. |

---

### `provenance_record`

```julia
provenance_record(label, source; kwargs...) -> ResultProvenance
provenance_record(provenance::ResultProvenance; kwargs...) -> ResultProvenance
provenance_record(payload::Union{AbstractDict,NamedTuple}) -> ResultProvenance
```

**Kind:** Constructor/normalizer

**Description:** Creates or normalizes a compact `ResultProvenance` record.

**Common keywords:** `id`, `status`, `warnings`, `errors`, `fallbacks`, `notes`, `parameters`, `parent_ids`, `timestamp`.

**Example:**

```julia
prov = provenance_record("PCA result", "run_pca";
    parameters=(n_components=50,), status=:ok)
```

---

### `register_provenance!`

```julia
register_provenance!(ctx, node::ProvenanceNode)
register_provenance!(ctx, operation; parents=String[], parameters=NamedTuple())
register_provenance!(ctx, node_id, operation; parents=String[], parameters=NamedTuple())
```

**Kind:** Mutating graph registration function

**Description:** Adds a node to a provenance context, enforcing context limits. Thread-safe overloads lock the wrapper context before mutating.

**Returns:** `ProvenanceNode`

**Errors:** Throws `ArgumentError` when `max_nodes` or `max_depth` would be exceeded.

---

### `provenance_result!`

```julia
provenance_result!(ctx, result, operation; kwargs...) -> result
provenance_result!(nothing, result, operation; kwargs...) -> result
```

**Kind:** Result stamping and graph registration function

**Description:** Registers an operation node and stamps supported result containers. Passing `nothing` is a no-op.

**Important keywords:**

| Keyword | Description |
|---|---|
| `parents` | Parent node ids. |
| `parameters` | Operation parameters. |
| `id` | Explicit node/container id. |
| `label` | Human-facing result label. |
| `source` | Source label stored in embedded provenance. |
| `status` | Status symbol. |
| `warnings`, `errors`, `fallbacks`, `notes` | Audit messages. |
| `provenance_hash` | Optional data hash reference. |

---

### `stamp_provenance!` and `update_provenance!`

```julia
stamp_provenance!(metadata::AbstractDict, provenance::ResultProvenance)
stamp_provenance!(table::DataFrames.AbstractDataFrame, provenance::ResultProvenance)
stamp_provenance!(container; kwargs...)

update_provenance!(container; kwargs...)
```

**Kind:** Metadata mutation functions

**Description:** `stamp_provenance!` stores a provenance record in supported metadata. `update_provenance!` merges updates into an existing record, appending warning/error/fallback/note vectors and merging parameters.

**Supported containers:** dictionaries, DataFrames, and objects exposing `metadata` or `annotations` dictionaries.

---

### `with_provenance(result, label, source; kwargs...)`

```julia
with_provenance(result, label::AbstractString, source::AbstractString; kwargs...) -> result
```

**Kind:** Existing-result stamping helper

**Description:** Attaches provenance metadata to an existing result/table/container. If a scoped or global context is active, the helper can also register a graph node unless `_register_graph=false` is used internally.

---

## 19. Metadata Query API

| Function | Signature | Description |
|---|---|---|
| `metadata_provenance` | `(metadata_or_table)` | Returns embedded `ResultProvenance` or `nothing`. |
| `container_provenance_id` | `(container)` | Returns stored provenance id or `nothing`. |
| `ensure_provenance_id!` | `(metadata_or_table)` | Returns existing id or creates one. |
| `container_provenance_hash` | `(container)` | Returns stored hash or `nothing`. |
| `container_provenance_summary` | `(container)` | Returns compact id/hash summary. |
| `provenance_parent_ids` | `(values...)` | Extracts unique provenance ids from nested containers. |
| `has_provenance` | `(result)` | Checks whether embedded provenance exists. |
| `requires_provenance` | `(result)` | Returns result or throws `ArgumentError`. |
| `flatten_provenance!` | `(df::DataFrame)` | Copies DataFrame metadata id/hash into ordinary columns. |

---

## 20. Graph Inspection API

| Function | Returns | Description |
|---|---|---|
| `find_nodes(ctx, pattern)` | `Vector{ProvenanceNode}` | Finds operations matching a regex-compatible pattern. |
| `provenance_chain(ctx, node_id)` | `Vector{ProvenanceNode}` | Ancestor chain, oldest first, including queried node when present. |
| `provenance_ancestors(ctx, node_id)` | `Vector{ProvenanceNode}` | Upstream nodes excluding the queried node. |
| `provenance_descendants(ctx, node_id)` | `Vector{ProvenanceNode}` | Downstream dependent nodes. |
| `provenance_lineage_table(ctx)` | `DataFrame` | One row per node. |
| `provenance_diff(ctx_a, ctx_b)` | `NamedTuple` | Nodes only in A, only in B, and common by id. |
| `merge_provenance_contexts(contexts...)` | `ProvenanceContext` | Merges node sets into a new context. |
| `provenance_to_dot(ctx)` | `String` | Graphviz DOT graph. |
| `provenance_to_mermaid(ctx)` | `String` | Mermaid `flowchart LR`. |

---

## 21. Serialization and Hashing API

### `export_provenance_json` / `import_provenance_json`

```julia
export_provenance_json(ctx) -> String
import_provenance_json(json_text) -> ProvenanceContext
import_provenance_json(ThreadSafeProvenanceContext, json_text) -> ThreadSafeProvenanceContext
```

**Description:** Serializes/restores provenance graphs using the package's PROV-JSON-style schema. Import validates the `schema_version` field.

---

### Hashing functions

| Function | Description | Cost |
|---|---|---|
| `provenance_structural_hash(data)` | Hashes shape/type-style structure for arrays/tables. Does not inspect array values. | Usually O(1) for arrays. |
| `provenance_content_hash(data)` | Streams through values and computes true SHA-256. | O(n). |
| `file_provenance_hash(path)` | Streams a file in chunks and computes SHA-256. | O(file size). |
| `stamp_provenance_with_hash!(container, data; hash_mode=:structural)` | Stamps metadata and stores a structural or content hash. | Depends on hash mode. |
| `set_provenance_limits!(; max_array_len, max_string_len)` | Sets truncation limits for JSON-normalized parameters. | O(1). |

---

## 22. Macro API

### `@provenance`

```julia
@provenance ctx operation parameters expr
@provenance ctx operation parents parameters expr
```

**Kind:** Macro

**Description:** Runs `expr`, stamps the result, and registers an operation when `ctx` is a provenance context. If `ctx === nothing`, it simply evaluates `expr`.

**Parent behavior:** The four-argument form uses task-local auto-parent tracking. The five-argument form uses explicit parent ids.

---

### `@pure_provenance`

```julia
@pure_provenance ctx operation content_key parameters expr
```

**Kind:** Macro

**Description:** Like `@provenance`, but uses `new_provenance_id(String(content_key))` so the node id is deterministic for the same content key.

---

### `lazy_provenance_id!`

```julia
lazy_provenance_id!(ctx, result, operation; parents=String[], parameters=NamedTuple())
```

**Kind:** Hot-path helper

**Description:** Avoids provenance work when `ctx === nothing`; otherwise delegates to `provenance_result!`.

---

## 23. Methods Text API

### `generate_methods_section`

```julia
generate_methods_section(ctx::ProvenanceContext; title="## Methods", max_ops=nothing) -> String
generate_methods_section(ctx::ThreadSafeProvenanceContext; title="## Methods", max_ops=nothing) -> String
generate_methods_section(result; title="## Methods", max_ops=nothing) -> String
```

**Kind:** Report-generation function

**Description:** Converts a provenance DAG into draft manuscript methods text. Nodes are topologically ordered from roots to leaves; operation names are converted to readable labels; parameters are rendered compactly; and a closing software/environment sentence is appended.

**Parameters:**

| Parameter | Description |
|---|---|
| `title` | Markdown heading. Pass `""` to suppress. |
| `max_ops` | Optional maximum number of operations to render. |

**Caution:** The output is a draft. Review and edit before publication.

---

## 24. Complete Usage Example

```julia
using BioToolkit
using DataFrames

ctx = ProvenanceContext()

raw = DataFrame(gene=["A", "B", "C"], count=[10, 5, 0])
stamp_provenance!(raw;
    label="raw counts",
    source="read_counts",
    parameters=(path="counts.csv",))

norm = DataFrame(gene=raw.gene, value=raw.count ./ maximum(raw.count))
provenance_result!(ctx, norm, "normalize_counts";
    parents=provenance_parent_ids(raw),
    parameters=(method=:max_scale,),
    label="normalized counts",
    source="BioToolkit/normalize_counts")

id = container_provenance_id(norm)
lineage = provenance_lineage_table(ctx)
dot = provenance_to_dot(ctx)
json = export_provenance_json(ctx)
methods = generate_methods_section(ctx)

requires_provenance(norm)
flatten_provenance!(norm)
```
