# `pipeline.jl` - Workflow Pipeline Utilities

## Overview

`pipeline.jl` defines file artifacts, pipeline nodes, dependency graphs, validation, execution levels, simple read/alignment/count/DE templates, distributed execution helpers, retry wrappers, containerized nodes, and Terra/AnVIL/SLURM planning utilities.

### Purpose

This page is a hand-authored reference for `pipeline.jl`, grouped around the exported user workflows. Internal helper functions are omitted unless the module exposes them as part of the public API.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Source-matched APIs** | Entries use names implemented in the corresponding `.jl` file. |
| **Workflow sections** | Related functions are documented together so the analysis path is clear. |
| **Concrete result descriptions** | Structs and result types are described by their role in downstream workflows. |
| **Julia-first data flow** | APIs compose through normal Julia arrays, tables, graphs, and BioToolkit objects. |
| **Export-oriented coverage** | The page covers the public functions users are likely to import from BioToolkit. |

---

## 1. Core Graph Types

Pipeline graphs track file inputs/outputs, commands/functions, and dependencies.

| API | Description |
|---|---|
| `FileArtifact` | Path plus metadata for a file consumed or produced by a workflow. |
| `PipelineNode` | One workflow step with inputs, outputs, command/function, resources, and metadata. |
| `PipelineGraph` | Collection of nodes and dependency edges. |
| `file_artifact` | Constructs a `FileArtifact`. |
| `pipeline_node` | Constructs a `PipelineNode`. |
| `add_node!` | Adds a node to a graph. |
| `add_dependency!` | Adds an edge between pipeline nodes. |
| `validate_pipeline` | Checks missing nodes, cycles, and artifact consistency. |
| `execution_levels` | Computes topological execution batches. |

## 2. Execution and Templates

Execution helpers run or describe common bioinformatics steps.

| API | Description |
|---|---|
| `execute_pipeline` | Runs a validated pipeline graph by dependency level. |
| `align_reads` | Template/action for read alignment steps. |
| `count_features` | Template/action for feature counting steps. |
| `template_read_fasta_node` | Creates a FASTA-reading pipeline node. |
| `template_align_reads_node` | Creates an alignment pipeline node. |
| `template_count_features_node` | Creates a counting pipeline node. |
| `template_differential_expression_node` | Creates a differential-expression node. |

## 3. Platforms and Reliability

Planning helpers prepare workflows for local distributed execution and external platforms.

| API | Description |
|---|---|
| `terra_workspace_plan` | Builds a Terra workspace planning payload. |
| `anvil_workspace_manifest` | Builds an AnVIL workspace manifest. |
| `distributed_pipeline_map` | Maps pipeline tasks across workers/processes. |
| `retry_execute_pipeline` | Runs a pipeline with retry handling for failed nodes. |
| `containerized_node` | Wraps a node in container execution metadata. |
| `slurm_array_plan` | Creates SLURM array-job planning records. |

---

## Complete Usage Example

```julia
using BioToolkit

g = PipelineGraph()
reads = file_artifact("reads.fq.gz")
aln = template_align_reads_node(reads, "ref.fa", "aligned.bam")
add_node!(g, aln)
validate_pipeline(g)
levels = execution_levels(g)
```

