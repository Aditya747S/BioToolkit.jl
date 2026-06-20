# `structure.jl` - Macromolecular Structure Utilities

## Overview

`structure.jl` provides atom/residue/chain/model/structure containers, PDB/mmCIF parsing and writing, coordinate extraction, Kabsch superposition, RMSD, ensemble statistics, torsion angles, KD-tree radius search, geometry payloads, point clouds, and wrappers for DSSP/PDB2PQR-style tools.

### Purpose

This documentation covers the public API implemented in `structure.jl`, with concrete descriptions for the result types and workflow functions exposed by BioToolkit.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Concrete workflow coverage** | Each section maps to a real analysis task implemented by the module. |
| **Typed results** | Important model outputs and summaries are represented with explicit structs. |
| **Downstream compatibility** | Results are shaped for plotting, tabulation, or reuse in other BioToolkit modules. |
| **Source-aligned API names** | Entries use implemented/exported names rather than speculative aliases. |
| **Readable defaults** | Examples show the minimal flow without hiding required biological inputs. |

---

## 1. Structure Types and I/O

Structure objects preserve hierarchy from atoms up to models and files.

| API | Description |
|---|---|
| `Atom` | Atomic coordinate, element, occupancy/B-factor, charge, and identifiers. |
| `Residue` | Residue with name, sequence number, insertion code, and atoms. |
| `Chain` | Chain containing residues. |
| `Model` | One structural model containing chains. |
| `Structure` | Full structure with id, models, metadata, and provenance. |
| `read_pdb` | Parses PDB text or file data. |
| `write_pdb` | Writes PDB text. |
| `read_mmcif` | Parses mmCIF text or file data. |
| `write_mmcif` | Writes mmCIF text. |

## 2. Coordinates and Alignment

Coordinate helpers support superposition, RMSD, and ensemble analysis.

| API | Description |
|---|---|
| `structure_atoms` | Returns all atoms from a structure/model/chain. |
| `atom_coordinates` | Returns `(x, y, z)` for an atom. |
| `coordinate_matrix` | Builds an `N x 3` coordinate matrix. |
| `kabsch` | Computes optimal rotation/translation between point sets. |
| `rmsd` | Computes RMSD with optional superposition. |
| `superpose` | Returns transformed mobile structure/coordinates. |
| `superpose!` | Applies superposition in place. |
| `superpose_models!` | Aligns models in an ensemble to a reference model. |
| `ensemble_rmsd_matrix` | Computes pairwise model RMSDs. |
| `trajectory_statistics` | Summarizes per-model RMSD and trajectory-like variation. |

## 3. Geometry and External Tools

Geometry functions derive bonds, torsions, point clouds, and nearby atoms.

| API | Description |
|---|---|
| `torsion_angle` | Computes dihedral angle from four points. |
| `phi_psi` | Computes backbone phi/psi for a residue in a chain. |
| `backbone_torsions` | Computes backbone torsions for a chain. |
| `KDTreeNode` | KD-tree node for atom spatial indexing. |
| `AtomKDTree` | KD-tree wrapper for atom searches. |
| `build_atom_kdtree` | Builds atom KD-tree. |
| `atoms_within_radius` | Finds atoms near a point or atom. |
| `structure_geometry` | Builds plot-ready atom/bond geometry. |
| `structure_pointcloud` | Builds point-cloud coordinates and annotations. |
| `run_dssp` | Runs DSSP/mkdssp on a path or temporary PDB from a structure. |
| `run_pdb2pqr` | Runs PDB2PQR-style conversion. |

---

## Complete Usage Example

```julia
using BioToolkit

structure = read_pdb("model.pdb")
atoms = structure_atoms(structure)
coords = coordinate_matrix(atoms)
nearby = atoms_within_radius(structure, atom_coordinates(first(atoms)); radius=5.0)
```

