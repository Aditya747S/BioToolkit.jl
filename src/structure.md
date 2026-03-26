# structure.jl

## Purpose
This file handles protein structure representation and geometry-heavy analysis. It combines atom-level records, chain/model hierarchy, distance calculations, structural superposition, contact analysis, residue properties, and multiple visualization helpers.

## Main types
- `Atom`, `Residue`, `Chain`, `Model`, and `Structure` model the structural hierarchy.
- `AtomSelectionPolicy` controls altloc and occupancy selection behavior.
- `SuperpositionResult` stores a rotation, translation, and RMSD value.
- `AtomKDTree` and its internal node type support spatial lookup.
- `DSSPEntry` stores secondary-structure annotations.
- `HydrogenBond` stores donor/acceptor geometry.

## Public functions
- Reading and writing: `read_pdb`, `read_mmcif`, `write_pdb`, and `write_mmcif`.
- Hierarchy and coordinate access: `structure_models`, `structure_chains`, `structure_residues`, `structure_atoms`, `atom_coordinates`, and `coordinate_matrix`.
- Geometry and superposition: `torsion_angle`, `phi_psi`, `backbone_torsions`, `kabsch`, `rmsd`, `superpose`, and `superpose!`.
- Spatial search: `build_atom_kdtree` and `atoms_within_radius`.
- Sequence and residue helpers: `residue_one_letter`, `sequence_from_structure`, `structure_sequences`, `residue_bfactor`, `residue_bfactors`, `flexible_residues`, `residue_property`, `select_residues`, `select_atoms`, and `collapse_altlocs`.
- Contacts and surface metrics: `residue_contacts`, `contact_map`, `interface_residues`, `residue_free_sasa`, `buried_surface_area`, `interface_profile`, `calculate_interface_residues`, and `residues_within_radius`.
- Annotation and external tools: `read_dssp`, `annotate_dssp!`, `run_dssp`, and `run_pdb2pqr`.
- Hydrogen bonds and rotamers: `hydrogen_bonds`, `ramachandran_region`, `ramachandran_profile`, `chi_angles`, `rotamer_label`, `rotamer_state`, and `rotamer_statistics`.
- Graphs, plots, and summaries: `structure_contact_graph`, `structure_contact_mermaid`, `write_structure_mermaid`, `contact_map_svg`, `write_contact_map_svg`, `plot_contact_graph`, `plot_contact_graph!`, `plot_contact_map`, `plot_contact_map!`, `plot_structure_atoms`, `plot_structure_atoms!`, `plot_backbone_trace`, `plot_backbone_trace!`, `plot_chain_ribbon`, `plot_chain_ribbon!`, `residue_pick_hooks`, `connect_residue_picking!`, `plot_structure_viewer`, `atomic_mass`, `center_of_mass`, `bounding_box`, `radius_of_gyration`, `atom_distance_matrix`, `residue_distance_matrix`, `chain_contact_matrix`, `structure_summary`, `chain_summary`, `superpose_models!`, `ensemble_rmsd_matrix`, and `trajectory_statistics`.

## Threading notes
- `atom_distance_matrix` and `residue_distance_matrix` now default to threaded CPU execution and still accept a `threaded` keyword for compatibility.
- The pairwise matrix fill is shared so the atom-level and residue-level distance helpers use the same symmetric threaded pattern.

## How it is used
The standard workflow is to read a structure, inspect its hierarchy, extract sequences or coordinates, and then use the geometry helpers to compare conformations or locate contacts. The plotting and Mermaid helpers make it easier to inspect the result without hand-building visualization code.

`AtomSelectionPolicy` exists to manage alternate locations and occupancy decisions consistently, while `AtomKDTree` supports repeated neighborhood searches over large structures.

## Implementation notes
- The hierarchy is mutable, which makes it easier to attach annotations and build up parsed structures incrementally.
- Residues and atoms carry enough metadata to support both sequence-centric and geometry-centric workflows.
- The file exposes a wide surface area because structural biology often requires both analysis and visualization utilities in the same module.

## Why it matters
Structure analysis only works well when the coordinate geometry and the biological hierarchy are both first-class. This file provides that combined representation and the common geometric operations that downstream structural workflows depend on.
