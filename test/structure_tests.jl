using Test
using BioToolkit

@testset "Structural biology" begin
    pdb_text = """
HEADER    TEST STRUCTURE
TITLE     METADATA ROUND TRIP
REMARK    1 EXAMPLE METADATA LINE
ATOM      1  N   ALA A   1      11.104  13.207   9.104  1.00 20.00           N  
ATOM      2  CA  ALA A   1      12.208  12.297   9.504  1.00 20.00           C  
ATOM      3  C   ALA A   1      13.410  12.980  10.104  1.00 20.00           C  
ATOM      4  N   GLY A   2      14.510  12.120  10.604  1.00 20.00           N  
ATOM      5  CA  GLY A   2      15.702  12.740  11.204  1.00 20.00           C  
ATOM      6  C   GLY A   2      16.900  11.900  11.804  1.00 20.00           C  
END
"""

    mktempdir() do dir
        pdb_path = joinpath(dir, "example.pdb")
        write(pdb_path, pdb_text)

        structure = BioToolkit.read_pdb(pdb_path)
        @test length(BioToolkit.structure_models(structure)) == 1
        @test length(BioToolkit.structure_chains(structure.models[1])) == 1
        @test length(BioToolkit.structure_residues(structure.models[1].chains[1])) == 2
        @test length(BioToolkit.structure_atoms(structure)) == 6
        @test occursin("TEST STRUCTURE", structure.metadata["pdb:HEADER"])
        @test occursin("METADATA ROUND TRIP", structure.metadata["pdb:TITLE"])
        @test size(BioToolkit.coordinate_matrix(structure)) == (6, 3)

        roundtrip_path = joinpath(dir, "roundtrip.pdb")
        BioToolkit.write_pdb(roundtrip_path, structure)
        roundtrip = BioToolkit.read_pdb(roundtrip_path)
        @test length(BioToolkit.structure_atoms(roundtrip)) == 6
        @test BioToolkit.structure_residues(roundtrip.models[1].chains[1])[1].name == "ALA"
        roundtrip_text = read(roundtrip_path, String)
        @test occursin("TITLE", roundtrip_text)
        @test occursin("REMARK", roundtrip_text)

        shifted = deepcopy(structure)
        for atom in BioToolkit.structure_atoms(shifted)
            atom.x += 5.0
            atom.y -= 2.0
        end

        @test BioToolkit.rmsd(structure, shifted; superpose=false) > 0.0
        @test BioToolkit.rmsd(structure, shifted) < 1e-6

        tree = BioToolkit.build_atom_kdtree(structure)
        neighbors = BioToolkit.atoms_within_radius(tree, BioToolkit.structure_atoms(structure)[1]; radius=2.5)
        @test !isempty(neighbors)

        geometry = BioToolkit.structure_geometry(structure)
        @test size(geometry.positions) == (6, 3)
        @test length(geometry.colors) == 6
        @test !isempty(geometry.bonds)
    end

    mmcif_text = """
    data_example
    loop_
    _entity.id
    _entity.type
    1 polymer
    #
    loop_
    _atom_site.group_PDB
    _atom_site.id
    _atom_site.type_symbol
    _atom_site.label_atom_id
    _atom_site.label_comp_id
    _atom_site.label_asym_id
    _atom_site.label_seq_id
    _atom_site.Cartn_x
    _atom_site.Cartn_y
    _atom_site.Cartn_z
    _atom_site.occupancy
    _atom_site.B_iso_or_equiv
    ATOM 1 N 'N' 'ALA' A 1 1.0 2.0 3.0 1.0 10.0
    ATOM 2 C 'CA' 'ALA' A 1 2.0 3.0 4.0 1.0 10.0
    ATOM 3 C C ALA A 1 3.0 4.0 5.0 1.0 10.0
    #
    """

    mktempdir() do dir
        cif_path = joinpath(dir, "example.cif")
        write(cif_path, mmcif_text)

        structure = BioToolkit.read_mmcif(cif_path)
        @test length(BioToolkit.structure_models(structure)) == 1
        @test length(BioToolkit.structure_atoms(structure)) == 3
        @test BioToolkit.structure_residues(structure.models[1].chains[1])[1].name == "ALA"
        @test occursin("_entity.id", structure.metadata["mmcif:raw_blocks"])
        @test occursin("_entity.id", structure.metadata["mmcif:category:_entity"])

        mmcif_roundtrip_path = joinpath(dir, "roundtrip.cif")
        BioToolkit.write_mmcif(mmcif_roundtrip_path, structure)
        mmcif_roundtrip_text = read(mmcif_roundtrip_path, String)
        @test occursin("_entity.id", mmcif_roundtrip_text)
        @test occursin("_atom_site.pdbx_PDB_model_num", mmcif_roundtrip_text)
    end

    chain = BioToolkit.Chain("A")
    push!(chain.residues, BioToolkit.Residue("ALA", 1, ' ' , [
        BioToolkit.Atom(1, "N", 0.0, 0.0, 0.0),
        BioToolkit.Atom(2, "CA", 1.0, 0.0, 0.0),
        BioToolkit.Atom(3, "C", 1.0, 1.0, 0.0),
    ]))
    push!(chain.residues, BioToolkit.Residue("GLY", 2, ' ' , [
        BioToolkit.Atom(4, "N", 2.0, 1.0, 0.0),
        BioToolkit.Atom(5, "CA", 2.0, 2.0, 0.0),
        BioToolkit.Atom(6, "C", 3.0, 2.0, 1.0),
    ]))
    push!(chain.residues, BioToolkit.Residue("SER", 3, ' ' , [
        BioToolkit.Atom(7, "N", 4.0, 2.0, 1.0),
        BioToolkit.Atom(8, "CA", 4.0, 3.0, 1.0),
        BioToolkit.Atom(9, "C", 5.0, 3.0, 2.0),
    ]))

    torsions = BioToolkit.phi_psi(chain, 2)
    @test !ismissing(torsions.phi)
    @test !ismissing(torsions.psi)
    @test isfinite(BioToolkit.torsion_angle((1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 1.0, 1.0)))

    chain_b = BioToolkit.Chain("B")
    push!(chain_b.residues, BioToolkit.Residue("LYS", 1, ' ', [
        BioToolkit.Atom(10, "N", 6.0, 3.0, 2.0),
        BioToolkit.Atom(11, "CA", 6.5, 3.5, 2.5),
        BioToolkit.Atom(12, "C", 7.0, 4.0, 3.0),
    ]))

    for atom in chain.residues[3].atoms
        atom.bfactor = 40.0
    end

    model = BioToolkit.Model(1)
    push!(model.chains, chain)
    push!(model.chains, chain_b)
    structure = BioToolkit.Structure("example")
    push!(structure.models, model)

    @test BioToolkit.sequence_from_structure(chain) == "AGS"
    @test BioToolkit.sequence_from_structure(structure; model_index=1) == "AGSK"
    @test BioToolkit.structure_sequences(structure)["1:A"] == "AGS"
    @test BioToolkit.residue_property(chain.residues[1], :hydrophobic)
    @test BioToolkit.residue_property(chain.residues[3], :polar)

    selection_chain = BioToolkit.Chain("Z")
    push!(selection_chain.residues, BioToolkit.Residue("GLY", 10, 'A', [
        BioToolkit.Atom(13, "CA", 0.0, 0.0, 0.0),
    ]))
    push!(selection_chain.residues, BioToolkit.Residue("GLY", 11, 'A', [
        BioToolkit.Atom(16, "CA", 1.0, 0.0, 0.0),
    ]))
    selection_model = BioToolkit.Model(1)
    push!(selection_model.chains, selection_chain)
    selection_structure = BioToolkit.Structure("selection")
    push!(selection_structure.models, selection_model)
    @test length(BioToolkit.select_residues(selection_structure; selector="chain Z and residue 10 and icode A")) == 1
    @test length(BioToolkit.select_residues(selection_structure; selector="chain Z and name GLY")) == 2
    @test length(BioToolkit.select_residues(selection_structure; selector="chain Z and (residue 10 or residue 11)")) == 2

    atom_chain = BioToolkit.Chain("Q")
    push!(atom_chain.residues, BioToolkit.Residue("ALA", 1, ' ', [
        BioToolkit.Atom(17, "N", 0.0, 0.0, 0.0),
        BioToolkit.Atom(18, "CA", 1.0, 0.0, 0.0),
    ]))
    atom_model = BioToolkit.Model(1)
    push!(atom_model.chains, atom_chain)
    atom_structure = BioToolkit.Structure("atom-selection")
    push!(atom_structure.models, atom_model)
    @test length(BioToolkit.select_atoms(atom_structure; selector="chain Q and atom name CA")) == 1
    @test BioToolkit.select_atoms(atom_structure; selector="chain Q and atom name CA")[1].name == "CA"

    altloc_chain = BioToolkit.Chain("Y")
    push!(altloc_chain.residues, BioToolkit.Residue("SER", 1, ' ', [
        BioToolkit.Atom(14, "CA", 1.0, 0.0, 0.0; altloc='A', occupancy=0.60),
        BioToolkit.Atom(15, "CA", 1.1, 0.0, 0.0; altloc='B', occupancy=0.40),
    ]))
    altloc_collapsed = BioToolkit.collapse_altlocs(altloc_chain)
    @test length(altloc_collapsed.residues[1].atoms) == 1
    @test altloc_collapsed.residues[1].atoms[1].altloc == 'A'
    altloc_model = BioToolkit.Model(1)
    push!(altloc_model.chains, altloc_chain)
    altloc_structure = BioToolkit.Structure("altloc")
    push!(altloc_structure.models, altloc_model)
    @test length(BioToolkit.select_atoms(altloc_structure; selector="name SER", policy=BioToolkit.AtomSelectionPolicy(altloc=:all, occupancy=:all))) == 2

    @test BioToolkit.residue_bfactor(chain.residues[3]) == 40.0
    bfactors = BioToolkit.residue_bfactors(structure)
    @test length(bfactors) == 4
    @test BioToolkit.flexible_residues(structure; percentile=90)[1].residue == "SER"

    contacts = BioToolkit.residue_contacts(chain; cutoff=5.0)
    @test !isempty(contacts)
    @test BioToolkit.contact_map(chain; cutoff=5.0)[1, 2]

    whole_contacts = BioToolkit.contact_map(structure; cutoff=8.0)
    @test size(whole_contacts.matrix) == (4, 4)
    @test any(startswith(label, "A:ALA1") for label in whole_contacts.labels)

    interfaces = BioToolkit.interface_residues(structure; cutoff=8.0)
    @test any(residue -> residue.name == "ALA", interfaces)
    @test any(residue -> residue.name == "LYS", interfaces)
    interface_profile = BioToolkit.interface_profile(structure; cutoff=8.0, samples=24)
    @test length(interface_profile) == 4
    @test any(entry -> entry.is_interface, interface_profile)
    @test BioToolkit.buried_surface_area(structure; samples=24) >= 0

    @test BioToolkit.ramachandran_region(-60.0, -45.0) == :alpha_helix
    @test BioToolkit.ramachandran_region(-135.0, 135.0) == :beta_sheet
    @test BioToolkit.ramachandran_region(60.0, 30.0) == :left_handed_helix
    @test length(BioToolkit.ramachandran_profile(chain)) == 3

    lys = BioToolkit.Residue("LYS", 10, ' ', [
        BioToolkit.Atom(20, "N", 0.0, 0.0, 0.0),
        BioToolkit.Atom(21, "CA", 1.0, 0.0, 0.0),
        BioToolkit.Atom(22, "CB", 1.0, 1.0, 0.0),
        BioToolkit.Atom(23, "CG", 2.0, 1.0, 0.0),
        BioToolkit.Atom(24, "CD", 2.0, 2.0, 0.0),
        BioToolkit.Atom(25, "CE", 3.0, 2.0, 0.0),
        BioToolkit.Atom(26, "NZ", 3.0, 3.0, 0.0),
    ])
    chis = BioToolkit.chi_angles(lys)
    @test !ismissing(chis.chi1)
    @test BioToolkit.rotamer_label(60.0) == :gauche_plus
    @test BioToolkit.rotamer_label(-60.0) == :gauche_minus
    @test BioToolkit.rotamer_label(180.0) == :trans
    @test !ismissing(BioToolkit.rotamer_state(lys).chi1_label)

    hb_chain = BioToolkit.Chain("H")
    push!(hb_chain.residues, BioToolkit.Residue("ALA", 1, ' ', [
        BioToolkit.Atom(30, "N", 0.0, 0.0, 0.0),
        BioToolkit.Atom(31, "CA", -1.0, 0.0, 0.0),
        BioToolkit.Atom(32, "C", -1.0, 0.0, 0.0),
        BioToolkit.Atom(33, "O", -0.5, 1.0, 0.0),
    ]))
    push!(hb_chain.residues, BioToolkit.Residue("GLY", 2, ' ', [
        BioToolkit.Atom(34, "N", 3.5, 1.0, 0.0),
        BioToolkit.Atom(35, "CA", 4.0, 1.0, 0.0),
        BioToolkit.Atom(36, "C", 4.0, 0.0, 0.0),
        BioToolkit.Atom(37, "O", 3.0, 0.0, 0.0),
    ]))

    hb_model = BioToolkit.Model(1)
    push!(hb_model.chains, hb_chain)
    hb_structure = BioToolkit.Structure("hb")
    push!(hb_structure.models, hb_model)

    hbonds = BioToolkit.hydrogen_bonds(hb_structure)
    @test length(hbonds) == 1
    @test hbonds[1].donor_residue == "ALA"
    @test hbonds[1].acceptor_residue == "GLY"

    rot_chain = BioToolkit.Chain("R")
    push!(rot_chain.residues, lys)

    rotamer_stats_chain = BioToolkit.rotamer_statistics(rot_chain)
    @test rotamer_stats_chain.total == 1
    @test rotamer_stats_chain.with_chi1 == 1
    rot_structure = BioToolkit.Structure("rot")
    rot_model = BioToolkit.Model(1)
    push!(rot_model.chains, chain)
    push!(rot_model.chains, rot_chain)
    push!(rot_structure.models, rot_model)

    rotamer_stats_structure = BioToolkit.rotamer_statistics(rot_structure)
    @test rotamer_stats_structure.total == 4
    @test rotamer_stats_structure.with_chi1 == 1
    @test haskey(rotamer_stats_structure.chains, "R")

    atom_matrix = BioToolkit.atom_distance_matrix(hb_structure)
    @test size(atom_matrix) == (8, 8)
    @test atom_matrix[1, 1] == 0.0
    @test atom_matrix[1, 2] ≈ atom_matrix[2, 1]

    residue_matrix_chain = BioToolkit.residue_distance_matrix(hb_chain)
    @test size(residue_matrix_chain.matrix) == (2, 2)
    @test length(residue_matrix_chain.labels) == 2

    residue_matrix_structure = BioToolkit.residue_distance_matrix(hb_structure)
    @test size(residue_matrix_structure.matrix) == (2, 2)
    @test startswith(residue_matrix_structure.labels[1], "H:")

    chain_matrix = BioToolkit.chain_contact_matrix(structure; cutoff=8.0)
    @test size(chain_matrix.matrix) == (2, 2)
    @test chain_matrix.matrix[1, 2] > 0

    com = BioToolkit.center_of_mass(hb_structure)
    @test com.x isa Float64
    @test com.y isa Float64
    @test com.z isa Float64

    bbox = BioToolkit.bounding_box(hb_structure)
    @test bbox.min[1] <= bbox.max[1]
    @test bbox.min[2] <= bbox.max[2]
    @test bbox.min[3] <= bbox.max[3]

    rg = BioToolkit.radius_of_gyration(hb_structure)
    @test rg > 0

    summary = BioToolkit.structure_summary(hb_structure)
    @test summary.chain_count == 1
    @test summary.residue_count == 2
    @test summary.atom_count == 8
    @test summary.center_of_mass !== nothing
    @test summary.bounding_box !== nothing

    summary_chain = BioToolkit.chain_summary(hb_chain; structure=hb_structure)
    @test summary_chain.chain_id == "H"
    @test summary_chain.residue_count == 2
    @test summary_chain.atom_count == 8

    hetero_chain = BioToolkit.Chain("X")
    push!(hetero_chain.residues, BioToolkit.Residue("HOH", 1, ' ', [
        BioToolkit.Atom(40, "O", 0.0, 0.0, 0.0; hetatm=true, element="O"),
    ]))
    push!(hetero_chain.residues, BioToolkit.Residue("LIG", 2, ' ', [
        BioToolkit.Atom(41, "C1", 2.0, 0.0, 0.0; hetatm=true, element="C"),
        BioToolkit.Atom(42, "C2", 2.5, 0.5, 0.0; hetatm=true, element="C"),
    ]))

    hetero_model = BioToolkit.Model(1)
    push!(hetero_model.chains, hetero_chain)
    hetero_structure = BioToolkit.Structure("hetero")
    push!(hetero_structure.models, hetero_model)

    @test BioToolkit.is_water_residue(hetero_chain.residues[1])
    @test BioToolkit.is_ligand_residue(hetero_chain.residues[2])
    @test length(BioToolkit.select_residues(hetero_structure; property=:water)) == 1
    @test length(BioToolkit.select_residues(hetero_structure; property=:ligand)) == 1
    @test length(BioToolkit.select_residues(hetero_structure; property=:hetero)) == 2

    sasa_total = BioToolkit.structure_sasa(hetero_structure; samples=24)
    @test sasa_total > 0
    sasa_profile = BioToolkit.sasa_profile(hetero_structure; samples=24)
    @test length(sasa_profile) == 2
    @test all(entry.kind in (:water, :ligand, :protein) for entry in sasa_profile)

    disulfide_chain = BioToolkit.Chain("D")
    push!(disulfide_chain.residues, BioToolkit.Residue("CYS", 1, ' ', [
        BioToolkit.Atom(50, "SG", 0.0, 0.0, 0.0; element="S"),
    ]))
    push!(disulfide_chain.residues, BioToolkit.Residue("CYS", 2, ' ', [
        BioToolkit.Atom(51, "SG", 2.0, 0.0, 0.0; element="S"),
    ]))
    disulfide_model = BioToolkit.Model(1)
    push!(disulfide_model.chains, disulfide_chain)
    disulfide_structure = BioToolkit.Structure("ss")
    push!(disulfide_structure.models, disulfide_model)

    ss_bonds = BioToolkit.disulfide_bonds(disulfide_structure)
    @test length(ss_bonds) == 1
    ss_summary = BioToolkit.structure_summary(disulfide_structure)
    @test ss_summary.disulfide_count == 1
    @test ss_summary.protein_count == 2

    disulfide_chain_summary = BioToolkit.chain_summary(disulfide_chain; structure=disulfide_structure)
    @test disulfide_chain_summary.disulfide_count == 1

    ensemble_structure = deepcopy(structure)
    push!(ensemble_structure.models, deepcopy(structure.models[1]))
    for residue in ensemble_structure.models[2].chains[1].residues
        for atom in residue.atoms
            atom.x += 2.0
        end
    end
    ensemble_matrix = BioToolkit.ensemble_rmsd_matrix(ensemble_structure; superpose=false)
    @test size(ensemble_matrix.matrix) == (2, 2)
    @test ensemble_matrix.matrix[1, 2] > 0
    trajectory_stats = BioToolkit.trajectory_statistics(ensemble_structure)
    @test trajectory_stats.model_count == 2
    @test length(trajectory_stats.per_model) == 2
    superposed = deepcopy(ensemble_structure)
    superpose_results = BioToolkit.superpose_models!(superposed)
    @test length(superpose_results) == 1
    superposed_matrix = BioToolkit.ensemble_rmsd_matrix(superposed; superpose=false)
    @test superposed_matrix.matrix[1, 2] < ensemble_matrix.matrix[1, 2]

    contact_graph = BioToolkit.structure_contact_graph(hb_structure; cutoff=5.0, hbond_cutoff=3.5)
    @test length(contact_graph.nodes) == 2
    @test length(contact_graph.generic_edges) == 1
    @test length(contact_graph.hydrogen_bond_edges) == 1
    @test length(contact_graph.layout) == 2

    svg_map = BioToolkit.contact_map_svg(hb_structure; cutoff=5.0, hbond_cutoff=3.5)
    @test startswith(svg_map, "<svg")
    @test occursin("Residue contact map", svg_map)

    if Base.find_package("Makie") !== nothing
        @eval using Makie
        graph_fig = BioToolkit.plot_contact_graph(hb_structure; cutoff=5.0, hbond_cutoff=3.5)
        map_fig = BioToolkit.plot_contact_map(hb_structure; cutoff=5.0, hbond_cutoff=3.5)
        atoms_fig = BioToolkit.plot_structure_atoms(hb_structure; atom_size=10)
        trace_fig = BioToolkit.plot_backbone_trace(hb_structure; line_width=3)
        ribbon_fig = BioToolkit.plot_chain_ribbon(hb_structure; ribbon_width=4)
        viewer_fig = BioToolkit.plot_structure_viewer(hb_structure; atom_size=8, ribbon_width=3)
        @test graph_fig isa Makie.Figure
        @test map_fig isa Makie.Figure
        @test atoms_fig isa Makie.Figure
        @test trace_fig isa Makie.Figure
        @test ribbon_fig isa Makie.Figure
        @test viewer_fig isa Makie.Figure

        hooks = BioToolkit.residue_pick_hooks(hb_structure)
        @test haskey(hooks.lookup, "r1")
        picked = hooks.resolve("r1")
        @test picked !== nothing
        @test picked.residue.name == "ALA"

        mktempdir() do dir
            graph_path = joinpath(dir, "graph.png")
            map_path = joinpath(dir, "map.png")
            atoms_path = joinpath(dir, "atoms.png")
            trace_path = joinpath(dir, "trace.png")
            ribbon_path = joinpath(dir, "ribbon.png")
            viewer_path = joinpath(dir, "viewer.png")
            Makie.save(graph_path, graph_fig)
            Makie.save(map_path, map_fig)
            Makie.save(atoms_path, atoms_fig)
            Makie.save(trace_path, trace_fig)
            Makie.save(ribbon_path, ribbon_fig)
            Makie.save(viewer_path, viewer_fig)
            @test isfile(graph_path)
            @test isfile(map_path)
            @test isfile(atoms_path)
            @test isfile(trace_path)
            @test isfile(ribbon_path)
            @test isfile(viewer_path)
        end
    else
        @test hasmethod(BioToolkit.plot_contact_graph, Tuple{BioToolkit.Structure})
        @test hasmethod(BioToolkit.plot_contact_map, Tuple{BioToolkit.Structure})
        @test hasmethod(BioToolkit.plot_structure_atoms, Tuple{BioToolkit.Structure})
        @test hasmethod(BioToolkit.plot_backbone_trace, Tuple{BioToolkit.Structure})
        @test hasmethod(BioToolkit.plot_chain_ribbon, Tuple{BioToolkit.Structure})
        @test hasmethod(BioToolkit.plot_structure_viewer, Tuple{BioToolkit.Structure})
        @test hasmethod(BioToolkit.residue_pick_hooks, Tuple{BioToolkit.Structure})
    end

    mktempdir() do dir
        svg_path = joinpath(dir, "contact-map.svg")
        BioToolkit.write_contact_map_svg(svg_path, hb_structure; cutoff=5.0, hbond_cutoff=3.5)
        svg_text = read(svg_path, String)
        @test startswith(svg_text, "<svg")
        @test occursin("#ef4444", svg_text)
    end

    dssp_entries = [
        BioToolkit.DSSPEntry("H", 1, ' ', 'A', 'H', 12.0),
        BioToolkit.DSSPEntry("H", 2, ' ', 'G', 'E', 34.0),
    ]
    BioToolkit.annotate_dssp!(hb_structure, dssp_entries)
    @test hb_structure.models[1].chains[1].residues[1].secondary_structure == 'H'
    @test hb_structure.models[1].chains[1].residues[2].secondary_structure == 'E'

    diagram = BioToolkit.structure_contact_mermaid(hb_structure; cutoff=4.0, hbond_cutoff=3.5)
    @test occursin("flowchart LR", diagram)
    @test occursin("-.->", diagram)
    @test occursin("<br/>H", diagram)

    mktempdir() do dir
        mermaid_path = joinpath(dir, "structure.mmd")
        BioToolkit.write_structure_mermaid(mermaid_path, hb_structure; cutoff=4.0, hbond_cutoff=3.5)
        mermaid_text = read(mermaid_path, String)
        @test startswith(mermaid_text, "flowchart LR")
        @test occursin("-.->", mermaid_text)
    end

    dssp_line = function (seqnum, chain_id, amino_acid, secondary_structure, accessibility; insertion=' ')
        chars = fill(' ', 40)
        seq_text = collect(lpad(string(seqnum), 5))
        chars[6:10] = seq_text
        chars[11] = insertion
        chars[12] = chain_id
        chars[14] = amino_acid
        chars[17] = secondary_structure
        acc_text = collect(lpad(string(accessibility), 3))
        chars[36:38] = acc_text
        return String(chars)
    end

    dssp_text = "\n" * join([
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC",
        dssp_line(1, 'A', 'A', 'H', 12),
        dssp_line(2, 'A', 'G', 'E', 34),
    ], "\n") * "\n"

    mktempdir() do dir
        dssp_path = joinpath(dir, "example.dssp")
        write(dssp_path, dssp_text)

        entries = BioToolkit.read_dssp(dssp_path)
        @test length(entries) == 2

        annotated = deepcopy(structure)
        BioToolkit.annotate_dssp!(annotated, entries)
        @test annotated.models[1].chains[1].residues[1].secondary_structure == 'H'
        @test annotated.models[1].chains[1].residues[2].secondary_structure == 'E'
        @test annotated.models[1].chains[1].residues[1].accessibility == 12.0
    end
end
