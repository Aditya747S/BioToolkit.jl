# `BioToolkitBody.jl` - Module Assembly and Public Namespace

## Overview

`BioToolkitBody.jl` is the package assembly file. It is included from `BioToolkit.jl` inside the top-level `BioToolkit` module and is responsible for loading dependencies, including source files, importing nested-module APIs, and exporting public bindings.

### Purpose

BioToolkit is organized as one broad user-facing namespace backed by many domain files and several nested modules. `BioToolkitBody.jl` defines the load order that makes this possible. Foundational files such as `analysisresults.jl`, `biotypes.jl`, `schema.jl`, and `io.jl` are included first; higher-level biology modules are included later.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Centralized include order** | Files that define shared infrastructure are loaded before files that depend on them. |
| **One public namespace** | Selected APIs from nested modules are imported into `BioToolkit`, giving users a cohesive toolkit rather than many separate module imports. |
| **Automatic export sweep** | `_export_public_bindings!` exports functions and types defined in `BioToolkit` and nested BioToolkit modules, reducing manual export drift. |
| **Explicit submodule imports** | Large nested modules such as `DifferentialExpression`, `SingleCell`, `GWAS`, and `Clinical` expose curated bindings via `using .Module: ...`. |
| **Browser loaded last** | `browser.jl` is included near the end after the APIs it may browse or inspect have been loaded. |

---

## 1. Dependency Imports

```julia
using Arrow
using Tables
using LinearAlgebra
using Distributions
using Statistics
using Random
using SpecialFunctions
using Printf
using BGZFStreams
using BigWig
```

**Kind:** Package dependency imports

**Description:** Loads shared dependencies used by multiple BioToolkit source files.

**Major roles:**

| Dependency | Typical role |
|---|---|
| `Arrow` | Columnar table storage and Arrow interoperability. |
| `Tables` | Table schema and table interface support. |
| `LinearAlgebra` | Matrix operations, decompositions, distances, model fitting. |
| `Distributions` | Statistical models and tests. |
| `Statistics` | Means, variances, quantiles, correlations. |
| `Random` | Simulations, initialization, reproducibility. |
| `SpecialFunctions` | Gamma/beta/log-special functions used in statistics. |
| `Printf` | Formatting. |
| `BGZFStreams` | Block-gzip genomics file streams. |
| `BigWig` | BigWig coverage track I/O. |

---

## 2. Source Include Order

The body file includes source modules in a dependency-aware order.

### Foundation layer

```julia
include("analysisresults.jl")
include("biotypes.jl")
include("schema.jl")
include("io.jl")
include("lazy_gpu.jl")
include("query.jl")
```

**Description:** Defines provenance, typed biological sequences, shared schemas, I/O utilities, GPU dispatch helpers, and query utilities.

### Core genomics and sequence layer

```julia
include("genomicranges.jl")
include("bam.jl")
include("record.jl")
include("quality.jl")
include("sequence.jl")
include("longread.jl")
include("protein.jl")
include("annotation.jl")
include("align.jl")
include("msa.jl")
include("motif.jl")
include("search.jl")
include("databases.jl")
```

**Description:** Loads sequence records, sequence analysis, alignment, MSA, motifs, database access, and genomic interval infrastructure.

### Modeling and domain modules

```julia
include("hmm.jl")
include("phylo.jl")
include("structure.jl")
include("coevolution.jl")
include("gene_prediction.jl")
include("popgen.jl")
include("differentialexpression.jl")
include("enrichment.jl")
include("singlecell.jl")
include("immunology.jl")
include("spatial.jl")
include("deeplearning.jl")
include("epigenetics.jl")
include("clinical.jl")
include("microbiome.jl")
include("proteomics.jl")
include("metabolomics.jl")
include("systemsbio.jl")
include("somatic.jl")
include("crispr.jl")
include("gwas.jl")
include("bioplotting.jl")
include("bioconductor_compat.jl")
include("pipeline.jl")
```

**Description:** Loads higher-level analysis modules spanning omics, population genetics, clinical genomics, systems biology, plotting, and workflow execution.

---

## 3. Nested Module Imports

After includes, `BioToolkitBody.jl` imports curated names from nested modules into the top-level namespace.

Examples:

```julia
using .Restriction: RestrictionEnzyme, restriction_sites, digest_sequence
using .DifferentialExpression: CountMatrix, DEResult, DESeq, results
using .SingleCell: SingleCellExperiment, normalize_counts, run_pca
using .GWAS: GenotypeMatrix, GWASResult, gwas_linear_scan
```

**Kind:** Namespace re-export preparation

**Description:** These imports make nested-module APIs available as top-level bindings. The later export sweep then exports public functions and types.

---

## 4. Browser Include

```julia
include("browser.jl")
```

**Description:** Loads the browser module after the broader API surface is assembled. This allows browser functionality to reference package types and workflows already loaded.

---

## 5. `_export_public_bindings!`

```julia
function _export_public_bindings!(mod::Module)
    for name in names(mod; all=true, imported=true)
        startswith(string(name), "_") && continue
        isdefined(mod, name) || continue
        value = getfield(mod, name)
        value isa Union{Function,Type} || continue
        parent = try
            parentmodule(value)
        catch
            nothing
        end
        parent === mod || continue
        @eval export $(name)
    end
end
```

**Kind:** Internal export helper

**Description:** Exports public functions and types whose parent module is the module being scanned. Names beginning with `_` are skipped. Non-function and non-type values are skipped.

**Important behavior:**

- avoids exporting private underscore-prefixed helpers;
- avoids exporting bindings merely imported from unrelated modules;
- exports functions and types defined in BioToolkit-owned modules.

---

## 6. Top-Level Export Sweep

```julia
_export_public_bindings!(@__MODULE__)
for name in names(@__MODULE__; all=true, imported=true)
    # find nested BioToolkit modules and export their public bindings
end
```

**Description:** First exports public functions/types defined directly in `BioToolkit`, then scans nested `BioToolkit.*` modules and exports their public functions/types.

**Effect for users:**

```julia
using BioToolkit

DNASeq("ACGT")
run_pca(...)
gwas_linear_scan(...)
```

The user does not normally need `using BioToolkit.SingleCell` or `using BioToolkit.GWAS`.

---

## Implementation Notes

### Include order matters

Moving includes can break downstream modules. For example, `record.jl`, `sequence.jl`, and `align.jl` rely on types from `biotypes.jl` and provenance helpers from `analysisresults.jl`.

### Export sweep is broad

The export sweep is intentionally broad and can expose new public functions automatically. New helper functions should use a leading underscore if they are internal.

### Nested module APIs are curated before export

The explicit `using .Module: ...` statements define the intended top-level API for many nested modules. When adding a new function to a nested module, add it there only if it should become part of the top-level public namespace.

---

## Quick Reference

| Section | Purpose |
|---|---|
| Dependency imports | Load shared package dependencies. |
| `include(...)` list | Load source files in dependency order. |
| `using .Module: ...` | Bring nested module APIs into `BioToolkit`. |
| `include("browser.jl")` | Load browser support after core APIs. |
| `_export_public_bindings!` | Export public functions and types. |

---

## Complete Usage Example

```julia
using BioToolkit

seq = DNASeq("ATGGCC")
protein = translate_dna(seq)
gc = gc_content(seq)
```

Those names become available because `BioToolkitBody.jl` loaded the relevant source files and exported their public bindings.
