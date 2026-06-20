# `substitution_matrices_data.jl` - Amino-Acid Substitution Matrix Data

## Overview

`substitution_matrices_data.jl` stores static amino-acid substitution matrices used by alignment, protein comparison, and phylogenetic distance routines. It intentionally contains data definitions rather than public analysis functions.

### Purpose

Keeping the matrices in a separate source file avoids mixing large lookup tables into algorithm modules. Alignment and protein utilities can include or import the data without duplicating matrix constants.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Data-only module file** | The file is a lookup-table source, not a workflow API. |
| **Centralized matrices** | Shared substitution scores live in one place for consistency. |
| **Algorithm-neutral storage** | Alignment, search, and phylogeny code can choose how to consume the tables. |
| **Stable amino-acid ordering** | Matrices rely on predictable residue ordering for fast score lookup. |

---

## 1. Matrix Data

| Data | Description |
|---|---|
| Substitution matrix constants | Amino-acid pair scoring tables used by alignment and protein distance functions. |
| Residue ordering metadata | Ordering information used to map residues to matrix rows/columns. |

---

## Complete Usage Example

```julia
using BioToolkit

# Users normally access substitution scores through alignment/search APIs,
# rather than reading this data file directly.
```

