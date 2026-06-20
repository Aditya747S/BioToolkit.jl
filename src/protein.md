# `protein.jl` - Protein Property Calculations

## Overview

`protein.jl` implements high-performance ProtParam-style protein property calculations using byte-indexed lookup tables.

### Purpose

Protein analysis often needs quick estimates of mass, extinction coefficient, instability, isoelectric point, hydropathicity, and aliphatic index. This file provides those calculations for `BioSequence{AminoAcidAlphabet}` and string inputs converted to `AASeq`.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Byte-indexed lookup tables** | 256-element arrays provide O(1) residue lookup with good cache locality. |
| **Typed protein inputs** | Primary APIs accept `BioSequence{AminoAcidAlphabet}`. |
| **String compatibility retained** | String overloads convert to `AASeq` for convenience. |
| **Separate internal calculators** | Public methods handle typing/provenance; `_protein_*` helpers operate on bytes. |
| **ProtParam conventions** | Calculations follow common ExPASy/ProtParam-style formulas and tables. |

---

## 1. Mass

### `protein_mass`

```julia
protein_mass(sequence; type="monoisotopic")
```

**Description:** Computes protein molecular mass from amino-acid residue masses.

**Parameters:**

| Parameter | Default | Description |
|---|---|---|
| `type` | `"monoisotopic"` | `"monoisotopic"` or average mass mode. |

---

## 2. Extinction Coefficient

### `extinction_coefficient`

```julia
extinction_coefficient(sequence)
```

**Description:** Estimates extinction coefficient from aromatic/cysteine content using Pace-style constants.

---

## 3. Instability

### `instability_index`

```julia
instability_index(sequence)
```

**Description:** Computes the Guruprasad dipeptide instability index from adjacent residue pairs.

---

## 4. Isoelectric Point

### `isoelectric_point`

```julia
isoelectric_point(sequence; precision=0.01)
```

**Description:** Estimates protein pI by searching pH values using pK tables.

---

## 5. Hydropathicity

### `gravy`

```julia
gravy(sequence)
```

**Description:** Computes GRAVY, the average Kyte-Doolittle hydropathicity.

---

## 6. Aliphatic Index

### `aliphatic_index`

```julia
aliphatic_index(sequence)
```

**Description:** Computes Ikai-style aliphatic index from alanine, valine, isoleucine, and leucine content.

---

## 7. Combined ProtParam Summary

### `protparam`

```julia
protparam(sequence)
```

**Description:** Returns a combined property summary containing the core protein metrics.

**Typical contents:**

- mass;
- extinction coefficient;
- instability index;
- isoelectric point;
- GRAVY;
- aliphatic index.

---

## Quick Reference

| API | Purpose |
|---|---|
| `protein_mass` | Monoisotopic or average protein mass. |
| `extinction_coefficient` | Aromatic/cysteine extinction estimate. |
| `instability_index` | Dipeptide instability score. |
| `isoelectric_point` | Estimated pI. |
| `gravy` | Average hydropathicity. |
| `aliphatic_index` | Thermostability-related aliphatic index. |
| `protparam` | Combined protein property summary. |

---

## Complete Usage Example

```julia
protein = AASeq("MVLSPADKTNVKAAW")

protein_mass(protein)
extinction_coefficient(protein)
instability_index(protein)
isoelectric_point(protein)
gravy(protein)
aliphatic_index(protein)

summary = protparam(protein)
```
