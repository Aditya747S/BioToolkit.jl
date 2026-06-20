# `BioToolkit.jl` - Package Entry Point

## Overview

`BioToolkit.jl` is the top-level package entry file. It defines the public `BioToolkit` module and delegates the actual package assembly to `BioToolkitBody.jl`.

The file is intentionally small:

```julia
module BioToolkit

include("Bio" * "ToolkitBody.jl")

end
```

### Purpose

The package entry point exists to give Julia one stable module boundary while keeping the long dependency-ordered include list in `BioToolkitBody.jl`. This separation keeps the package entry clean and lets the body file focus on loading submodules, importing selected names, and exporting public bindings.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **One top-level module** | All public APIs live under `BioToolkit`, so users can write `using BioToolkit` and access sequence, omics, structural, and clinical functionality from one namespace. |
| **Body file included indirectly** | `include("Bio" * "ToolkitBody.jl")` avoids writing the body filename as one literal string. Functionally it resolves to `BioToolkitBody.jl`. |
| **No public logic here** | The entry file does not define algorithms, types, or exports. It only establishes the module boundary. |

---

## 1. Module Definition

### `module BioToolkit`

```julia
module BioToolkit
```

**Kind:** Julia module declaration

**Description:** Defines the package namespace. Every file included by `BioToolkitBody.jl` is evaluated inside this module unless it defines its own nested module.

**User-facing effect:**

```julia
using BioToolkit
```

After loading, exported bindings from the body and nested modules become available.

---

## 2. Body Include

### `include("Bio" * "ToolkitBody.jl")`

```julia
include("Bio" * "ToolkitBody.jl")
```

**Kind:** Source include

**Description:** Loads `BioToolkitBody.jl` into the `BioToolkit` module. The string concatenation evaluates to `"BioToolkitBody.jl"`.

**Why this matters:** `BioToolkitBody.jl` handles dependency imports, source-file include order, submodule imports, and export collection.

---

## 3. Relationship to `BioToolkitBody.jl`

`BioToolkit.jl` is the shell; `BioToolkitBody.jl` is the implementation assembly.

Responsibilities of `BioToolkit.jl`:

- define the `BioToolkit` module;
- include the body file;
- close the module.

Responsibilities delegated to `BioToolkitBody.jl`:

- load package dependencies with `using`;
- include all source files in a dependency-safe order;
- import selected bindings from nested modules;
- export public functions and types.

---

## Quick Reference

| Construct | Purpose |
|---|---|
| `module BioToolkit` | Creates the package namespace. |
| `include("Bio" * "ToolkitBody.jl")` | Loads the source assembly file. |
| `end` | Closes the module. |

---

## Complete Usage Example

```julia
using BioToolkit

seq = DNASeq("ACGT")
gc_content(seq)
```

The user never needs to include this file manually; Julia loads it when the package is imported.
