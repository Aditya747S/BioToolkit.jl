using Printf
using Statistics

export Atom, Residue, Chain, Model, Structure, SuperpositionResult, AtomKDTree, DSSPEntry
export read_pdb, read_mmcif, write_pdb, write_mmcif
export structure_models, structure_chains, structure_residues, structure_atoms
export atom_coordinates, coordinate_matrix
export torsion_angle, phi_psi, backbone_torsions
export kabsch, rmsd, superpose, superpose!
export build_atom_kdtree, atoms_within_radius
export structure_geometry, structure_pointcloud
export residue_one_letter, sequence_from_structure, structure_sequences
export residue_bfactor, residue_bfactors, flexible_residues
export residue_property, select_residues, residue_contacts, contact_map, interface_residues
export select_atoms, collapse_altlocs, AtomSelectionPolicy
export read_dssp, annotate_dssp!, run_dssp, run_pdb2pqr
export HydrogenBond, hydrogen_bonds
export ramachandran_region, ramachandran_profile
export chi_angles, rotamer_label, rotamer_state
export rotamer_statistics
export structure_contact_graph, structure_contact_mermaid, write_structure_mermaid
export contact_map_svg, write_contact_map_svg
export plot_contact_graph!, plot_contact_graph, plot_contact_map!, plot_contact_map
export plot_structure_atoms!, plot_structure_atoms
export plot_backbone_trace!, plot_backbone_trace
export plot_chain_ribbon!, plot_chain_ribbon
export residue_pick_hooks, connect_residue_picking!, plot_structure_viewer
export atomic_mass, center_of_mass, bounding_box, radius_of_gyration
export atom_distance_matrix, residue_distance_matrix, chain_contact_matrix, structure_summary, chain_summary
export plot_structure_atoms!, plot_structure_atoms
export plot_backbone_trace!, plot_backbone_trace
export plot_chain_ribbon!, plot_chain_ribbon
export residue_pick_hooks, connect_residue_picking!, plot_structure_viewer
export residue_free_sasa, buried_surface_area, interface_profile, calculate_interface_residues, residues_within_radius
export superpose_models!, ensemble_rmsd_matrix, trajectory_statistics

mutable struct Atom
    serial::Int
    name::String
    altloc::Char
    element::String
    x::Float64
    y::Float64
    z::Float64
    occupancy::Float64
    bfactor::Float64
    charge::String
    hetatm::Bool
end

mutable struct Residue
    name::String
    seqnum::Int
    insertion_code::Char
    secondary_structure::Union{Missing,Char}
    accessibility::Union{Missing,Float64}
    atoms::Vector{Atom}
end

mutable struct Chain
    id::String
    residues::Vector{Residue}
end

mutable struct Model
    id::Int
    chains::Vector{Chain}
end

mutable struct Structure
    id::String
    models::Vector{Model}
    metadata::Dict{String,String}
end

struct AtomSelectionPolicy
    altloc::Symbol
    occupancy::Symbol
    insertion_code::Symbol
end

AtomSelectionPolicy(; altloc::Symbol=:highest_occupancy, occupancy::Symbol=:highest, insertion_code::Symbol=:keep) = AtomSelectionPolicy(altloc, occupancy, insertion_code)

struct SuperpositionResult
    rotation::Matrix{Float64}
    translation::Vector{Float64}
    rmsd::Float64
end

mutable struct KDTreeNode{T}
    point::NTuple{3,Float64}
    payload::T
    axis::Int
    left::Union{Nothing,KDTreeNode{T}}
    right::Union{Nothing,KDTreeNode{T}}
end

struct AtomKDTree
    root::Union{Nothing,KDTreeNode{Atom}}
end

struct DSSPEntry
    chain::String
    seqnum::Int
    insertion_code::Char
    amino_acid::Char
    secondary_structure::Char
    accessibility::Float64
end

struct HydrogenBond
    donor_chain::String
    donor_residue::String
    donor_seqnum::Int
    donor_insertion_code::Char
    donor_atom::String
    acceptor_chain::String
    acceptor_residue::String
    acceptor_seqnum::Int
    acceptor_insertion_code::Char
    acceptor_atom::String
    distance::Float64
    angle::Float64
end

Structure(id::AbstractString; metadata::AbstractDict=Dict{String,String}()) = Structure(String(id), Model[], Dict{String,String}(string(key) => string(value) for (key, value) in metadata))
Model(id::Integer) = Model(Int(id), Chain[])
Chain(id::AbstractString) = Chain(String(id), Residue[])
Residue(name::AbstractString, seqnum::Integer, insertion_code::Char=' '; secondary_structure::Union{Missing,Char}=missing, accessibility::Union{Missing,Real}=missing, atoms::AbstractVector{Atom}=Atom[]) = Residue(String(name), Int(seqnum), insertion_code, secondary_structure, accessibility === missing ? missing : Float64(accessibility), Vector{Atom}(atoms))
Residue(name::AbstractString, seqnum::Integer, insertion_code::Char, atoms::AbstractVector{Atom}) = Residue(String(name), Int(seqnum), insertion_code, missing, missing, Vector{Atom}(atoms))
Atom(serial::Integer, name::AbstractString, x::Real, y::Real, z::Real; altloc::Char=' ', element::AbstractString="", occupancy::Real=1.0, bfactor::Real=0.0, charge::AbstractString="", hetatm::Bool=false) = Atom(Int(serial), String(name), altloc, String(element), Float64(x), Float64(y), Float64(z), Float64(occupancy), Float64(bfactor), String(charge), Bool(hetatm))

function _metadata_append!(metadata::Dict{String,String}, key::AbstractString, value::AbstractString)
    existing = get(metadata, String(key), "")
    metadata[String(key)] = isempty(existing) ? String(value) : string(existing, "\n", value)
    return metadata
end

function _metadata_push_line!(metadata::Dict{String,String}, key::AbstractString, value::AbstractString)
    _metadata_append!(metadata, key, strip(String(value)))
end

function _metadata_block!(metadata::Dict{String,String}, key::AbstractString, block::AbstractString)
    metadata[String(key)] = String(block)
    return metadata
end

Base.show(io::IO, atom::Atom) = print(io, "Atom(", atom.name, ", ", round(atom.x, digits=3), ", ", round(atom.y, digits=3), ", ", round(atom.z, digits=3), ")")
Base.show(io::IO, residue::Residue) = print(io, "Residue(", residue.name, residue.seqnum, ", ss=", residue.secondary_structure, ", atoms=", length(residue.atoms), ")")
Base.show(io::IO, chain::Chain) = print(io, "Chain(", chain.id, ", residues=", length(chain.residues), ")")
Base.show(io::IO, model::Model) = print(io, "Model(", model.id, ", chains=", length(model.chains), ")")
Base.show(io::IO, structure::Structure) = print(io, "Structure(", structure.id, ", models=", length(structure.models), ")")
Base.show(io::IO, result::SuperpositionResult) = print(io, "SuperpositionResult(rmsd=", round(result.rmsd, digits=4), ")")
Base.show(io::IO, tree::AtomKDTree) = print(io, "AtomKDTree(", tree.root === nothing ? 0 : 1, ")")

structure_models(structure::Structure) = structure.models
structure_chains(model::Model) = model.chains
structure_residues(chain::Chain) = chain.residues
structure_atoms(residue::Residue) = residue.atoms

function structure_atoms(structure::Structure)
    atoms = Atom[]
    for model in structure.models
        for chain in model.chains
            for residue in chain.residues
                append!(atoms, residue.atoms)
            end
        end
    end
    return atoms
end

const _AA3_TO1 = Dict{String,Char}(
    "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D', "CYS" => 'C',
    "GLN" => 'Q', "GLU" => 'E', "GLY" => 'G', "HIS" => 'H', "ILE" => 'I',
    "LEU" => 'L', "LYS" => 'K', "MET" => 'M', "PHE" => 'F', "PRO" => 'P',
    "SER" => 'S', "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
    "ASX" => 'B', "GLX" => 'Z', "SEC" => 'U', "PYL" => 'O', "UNK" => 'X'
)

const _HYDROPHOBIC_RESIDUES = Set(["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "PRO", "CYS"])
const _CHARGED_RESIDUES = Set(["ASP", "GLU", "LYS", "ARG", "HIS"])
const _AROMATIC_RESIDUES = Set(["PHE", "TYR", "TRP", "HIS"])
const _POLAR_RESIDUES = Set(["ASN", "GLN", "SER", "THR", "CYS", "TYR", "HIS"])
const _WATER_RESIDUES = Set(["HOH", "H2O", "WAT", "DOD", "TIP", "SOL", "OH2"])

function residue_one_letter(residue::Residue)
    return get(_AA3_TO1, uppercase(residue.name), 'X')
end

function sequence_from_structure(chain::Chain)
    io = IOBuffer()
    for residue in chain.residues
        print(io, residue_one_letter(residue))
    end
    return String(take!(io))
end

function sequence_from_structure(structure::Structure; model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    return join(sequence_from_structure(chain) for chain in structure.models[model_index].chains)
end

function structure_sequences(structure::Structure)
    sequences = Dict{String,String}()
    for model in structure.models
        for chain in model.chains
            sequences[string(model.id, ":", chain.id)] = sequence_from_structure(chain)
        end
    end
    return sequences
end

const _ATOMIC_MASS = Dict{String,Float64}(
    "H" => 1.008,
    "C" => 12.011,
    "N" => 14.007,
    "O" => 15.999,
    "P" => 30.974,
    "S" => 32.06,
    "SE" => 78.971,
    "F" => 18.998,
    "CL" => 35.45,
    "BR" => 79.904,
    "I" => 126.904,
)

function _atom_element_symbol(atom::Atom)
    element = uppercase(strip(atom.element))
    isempty(element) && return ""
    return element
end

function atomic_mass(atom::Atom)
    element = _atom_element_symbol(atom)
    if isempty(element)
        cleaned = filter(isletter, atom.name)
        isempty(cleaned) && return 12.0
        element = uppercase(string(cleaned[1]))
    end
    return get(_ATOMIC_MASS, element, 12.0)
end

function _entity_atoms(entity)
    entity isa Atom && return Atom[entity]
    entity isa Residue && return entity.atoms
    entity isa Chain && return structure_atoms(entity)
    entity isa Structure && return structure_atoms(entity)
    entity isa AbstractVector{Atom} && return Vector{Atom}(entity)
    throw(ArgumentError("unsupported entity for atom operations"))
end

function center_of_mass(entity)
    atoms = _entity_atoms(entity)
    isempty(atoms) && throw(ArgumentError("center of mass requires at least one atom"))
    total_mass = 0.0
    x = 0.0
    y = 0.0
    z = 0.0
    for atom in atoms
        mass = atomic_mass(atom)
        total_mass += mass
        x += atom.x * mass
        y += atom.y * mass
        z += atom.z * mass
    end
    total_mass == 0 && return (x=0.0, y=0.0, z=0.0)
    return (x=x / total_mass, y=y / total_mass, z=z / total_mass)
end

function bounding_box(entity)
    atoms = _entity_atoms(entity)
    isempty(atoms) && throw(ArgumentError("bounding box requires at least one atom"))
    min_x = atoms[1].x
    max_x = atoms[1].x
    min_y = atoms[1].y
    max_y = atoms[1].y
    min_z = atoms[1].z
    max_z = atoms[1].z
    for atom in atoms[2:end]
        min_x = min(min_x, atom.x)
        max_x = max(max_x, atom.x)
        min_y = min(min_y, atom.y)
        max_y = max(max_y, atom.y)
        min_z = min(min_z, atom.z)
        max_z = max(max_z, atom.z)
    end
    return (min=(min_x, min_y, min_z), max=(max_x, max_y, max_z))
end

function radius_of_gyration(entity)
    atoms = _entity_atoms(entity)
    isempty(atoms) && throw(ArgumentError("radius of gyration requires at least one atom"))
    com = center_of_mass(atoms)
    total_mass = 0.0
    weighted_sum = 0.0
    for atom in atoms
        mass = atomic_mass(atom)
        total_mass += mass
        dx = atom.x - com.x
        dy = atom.y - com.y
        dz = atom.z - com.z
        weighted_sum += mass * (dx * dx + dy * dy + dz * dz)
    end
    total_mass == 0 && return 0.0
    return sqrt(weighted_sum / total_mass)
end

function _pairwise_distance_matrix(points::AbstractVector{<:Any}; threaded::Bool=true)
    matrix = Matrix{Float64}(undef, length(points), length(points))
    if threaded && length(points) > 1 && Threads.nthreads() > 1
        Threads.@threads for i in eachindex(points)
            matrix[i, i] = 0.0
            @inbounds for j in (i + 1):lastindex(points)
                dx = points[i][1] - points[j][1]
                dy = points[i][2] - points[j][2]
                dz = points[i][3] - points[j][3]
                distance = sqrt(dx * dx + dy * dy + dz * dz)
                matrix[i, j] = distance
                matrix[j, i] = distance
            end
        end
    else
        for i in eachindex(points)
            matrix[i, i] = 0.0
            @inbounds for j in (i + 1):lastindex(points)
                dx = points[i][1] - points[j][1]
                dy = points[i][2] - points[j][2]
                dz = points[i][3] - points[j][3]
                distance = sqrt(dx * dx + dy * dy + dz * dz)
                matrix[i, j] = distance
                matrix[j, i] = distance
            end
        end
    end
    return matrix
end

function atom_distance_matrix(atoms::AbstractVector{Atom}; threaded::Bool=true)
    points = [(atom.x, atom.y, atom.z) for atom in atoms]
    return _pairwise_distance_matrix(points; threaded=threaded)
end

atom_distance_matrix(structure::Structure; threaded::Bool=true) = atom_distance_matrix(structure_atoms(structure); threaded=threaded)

function _residue_reference_point(residue::Residue; atom_name::AbstractString="CA")
    atom = _residue_atom(residue, atom_name)
    atom !== nothing && return atom_coordinates(atom)
    isempty(residue.atoms) && return nothing
    return _residue_center(residue)
end

function residue_distance_matrix(chain::Chain; atom_name::AbstractString="CA", threaded::Bool=true)
    labels = String[]
    points = NTuple{3,Float64}[]
    for residue in chain.residues
        point = _residue_reference_point(residue; atom_name=atom_name)
        point === nothing && continue
        push!(labels, string(residue.name, residue.seqnum, residue.insertion_code))
        push!(points, point)
    end
    matrix = _pairwise_distance_matrix(points; threaded=threaded)
    return (matrix=matrix, labels=labels)
end

function residue_distance_matrix(structure::Structure; model_index::Int=1, atom_name::AbstractString="CA", threaded::Bool=true)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    labels = String[]
    points = NTuple{3,Float64}[]
    for chain in structure.models[model_index].chains
        for residue in chain.residues
            point = _residue_reference_point(residue; atom_name=atom_name)
            point === nothing && continue
            push!(labels, string(chain.id, ":", residue.name, residue.seqnum, residue.insertion_code))
            push!(points, point)
        end
    end
    matrix = _pairwise_distance_matrix(points; threaded=threaded)
    return (matrix=matrix, labels=labels)
end

function chain_contact_matrix(structure::Structure; model_index::Int=1, cutoff::Real=8.0)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    model = structure.models[model_index]
    labels = [chain.id for chain in model.chains]
    matrix = zeros(Int, length(labels), length(labels))
    lookup = Dict{String,Int}(chain.id => index for (index, chain) in enumerate(model.chains))
    for contact in residue_contacts(structure; cutoff=cutoff, model_index=model_index)
        left = lookup[contact.left_chain]
        right = lookup[contact.right_chain]
        matrix[left, right] += 1
        matrix[right, left] += 1
    end
    return (matrix=matrix, labels=labels)
end

function structure_summary(structure::Structure; model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    model = structure.models[model_index]
    atoms = _model_atoms(structure, model_index)
    residues = [residue for chain in model.chains for residue in chain.residues]
    bbox = isempty(atoms) ? nothing : bounding_box(atoms)
    com = isempty(atoms) ? nothing : center_of_mass(atoms)
    rog = isempty(atoms) ? nothing : radius_of_gyration(atoms)
    residue_counts = Dict(chain.id => length(chain.residues) for chain in model.chains)
    protein_count = count(is_protein_residue(residue) for residue in residues)
    hetero_count = count(is_hetero_residue(residue) for residue in residues)
    water_count = count(is_water_residue(residue) for residue in residues)
    ligand_count = count(is_ligand_residue(residue) for residue in residues)
    interface_entries = interface_profile(structure; model_index=model_index)
    return (
        structure_id=structure.id,
        model_id=model.id,
        chain_count=length(model.chains),
        residue_count=length(residues),
        atom_count=length(atoms),
        protein_count=protein_count,
        hetero_count=hetero_count,
        water_count=water_count,
        ligand_count=ligand_count,
        interface_count=count(entry -> entry.is_interface, interface_entries),
        disulfide_count=length(disulfide_bonds(structure; model_index=model_index)),
        buried_surface_area=buried_surface_area(structure; model_index=model_index),
        chain_ids=[chain.id for chain in model.chains],
        residue_counts=residue_counts,
        center_of_mass=com,
        radius_of_gyration=rog,
        bounding_box=bbox,
    )
end

function chain_summary(chain::Chain; structure::Union{Nothing,Structure}=nothing, model_index::Int=1)
    residues = chain.residues
    atoms = Atom[]
    for residue in residues
        append!(atoms, residue.atoms)
    end
    interface_count = 0
    buried_area = nothing
    if structure !== nothing
        interface_count = count(entry -> entry.chain == chain.id && entry.is_interface, interface_profile(structure; model_index=model_index))
        buried_area = buried_surface_area(chain, structure; model_index=model_index)
    end
    return (
        chain_id=chain.id,
        residue_count=length(residues),
        atom_count=length(atoms),
        protein_count=count(is_protein_residue(residue) for residue in residues),
        hetero_count=count(is_hetero_residue(residue) for residue in residues),
        water_count=count(is_water_residue(residue) for residue in residues),
        ligand_count=count(is_ligand_residue(residue) for residue in residues),
        interface_count=interface_count,
        disulfide_count=structure === nothing ? 0 : count(bond -> bond.left_chain == chain.id || bond.right_chain == chain.id, disulfide_bonds(structure; model_index=model_index)),
        buried_surface_area=buried_area,
        center_of_mass=isempty(atoms) ? nothing : center_of_mass(atoms),
        radius_of_gyration=isempty(atoms) ? nothing : radius_of_gyration(atoms),
    )
end

function residue_bfactor(residue::Residue)
    isempty(residue.atoms) && return missing
    return mean(atom.bfactor for atom in residue.atoms)
end

function residue_bfactors(structure::Structure; model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    values = NamedTuple{(:chain, :residue, :seqnum, :bfactor, :secondary_structure, :accessibility), Tuple{String,String,Int,Union{Missing,Float64},Union{Missing,Char},Union{Missing,Float64}}}[]
    for chain in structure.models[model_index].chains
        for residue in chain.residues
            push!(values, (chain=chain.id, residue=residue.name, seqnum=residue.seqnum, bfactor=residue_bfactor(residue), secondary_structure=residue.secondary_structure, accessibility=residue.accessibility))
        end
    end
    return values
end

function flexible_residues(structure::Structure; percentile::Real=90, model_index::Int=1)
    entries = residue_bfactors(structure; model_index=model_index)
    numeric = Float64[entry.bfactor for entry in entries if entry.bfactor !== missing]
    isempty(numeric) && return typeof(entries)([])
    threshold = sort(numeric)[clamp(ceil(Int, length(numeric) * percentile / 100), 1, length(numeric))]
    return [entry for entry in entries if entry.bfactor !== missing && entry.bfactor >= threshold]
end

function residue_property(residue::Residue, property::Symbol)
    key = uppercase(residue.name)
    property === :hydrophobic && return key in _HYDROPHOBIC_RESIDUES
    property === :charged && return key in _CHARGED_RESIDUES
    property === :aromatic && return key in _AROMATIC_RESIDUES
    property === :polar && return key in _POLAR_RESIDUES
    property === :unknown && return !(key in keys(_AA3_TO1))
    throw(ArgumentError("unsupported residue property: $(property)"))
end

function is_protein_residue(residue::Residue)
    return uppercase(residue.name) in keys(_AA3_TO1)
end

function is_water_residue(residue::Residue)
    return uppercase(residue.name) in _WATER_RESIDUES
end

function is_hetero_residue(residue::Residue)
    return !is_protein_residue(residue)
end

function is_ligand_residue(residue::Residue)
    return is_hetero_residue(residue) && !is_water_residue(residue)
end

function _selection_compare(left::Real, operator::AbstractString, right::Real)
    operator == "<" && return left < right
    operator == "<=" && return left <= right
    operator == ">" && return left > right
    operator == ">=" && return left >= right
    operator == "==" && return left == right
    operator == "!=" && return left != right
    throw(ArgumentError("unsupported comparison operator: $(operator)"))
end

function _selection_tokens(selector::AbstractString)
    tokens = String[]
    buffer = IOBuffer()
    function flush_buffer!()
        if position(buffer) > 0
            push!(tokens, String(take!(buffer)))
        end
        return nothing
    end
    for char in selector
        if isspace(char)
            flush_buffer!()
        elseif char == '(' || char == ')'
            flush_buffer!()
            push!(tokens, string(char))
        else
            write(buffer, char)
        end
    end
    flush_buffer!()
    return tokens
end

function _selection_leaf_matches(chain::Chain, residue::Residue, clause::AbstractString, atom::Union{Nothing,Atom}=nothing)
    cleaned = lowercase(strip(String(clause)))
    isempty(cleaned) && return true
    cleaned in ("protein", "water", "ligand", "hetero") && return (
        cleaned == "protein" ? is_protein_residue(residue) :
        cleaned == "water" ? is_water_residue(residue) :
        cleaned == "ligand" ? is_ligand_residue(residue) :
        is_hetero_residue(residue)
    )

    atom_match = match(r"^(?:atom|atomname|atom_name|atom_id|atomid|serial)\s*[:= ]\s*(.+)$", cleaned)
    if atom_match !== nothing
        atom_clause = strip(atom_match.captures[1])
        atom_kind = if startswith(cleaned, "serial") || startswith(cleaned, "atom_id") || startswith(cleaned, "atomid")
            :serial
        elseif startswith(atom_clause, "element ")
            :element
        else
            :name
        end
        if atom_kind == :element
            atom_clause = strip(atom_clause[9:end])
        elseif startswith(atom_clause, "name ")
            atom_clause = strip(atom_clause[6:end])
        elseif startswith(atom_clause, "serial ")
            atom_clause = strip(atom_clause[8:end])
        end
        atom_predicate = candidate_atom -> begin
            atom_kind === :serial && return candidate_atom.serial in _selection_int_set(atom_clause)
            atom_kind === :element && return uppercase(candidate_atom.element) in _selection_string_set(atom_clause)
            return strip(candidate_atom.name) in _selection_string_set(atom_clause)
        end
        atom !== nothing && return atom_predicate(atom)
        return any(atom_predicate(candidate_atom) for candidate_atom in residue.atoms)
    end

    chain_match = match(r"^chain\s*[:= ]\s*(.+)$", cleaned)
    if chain_match !== nothing
        return uppercase(strip(chain.id)) in _selection_string_set(chain_match.captures[1])
    end

    residue_match = match(r"^(?:residue|seqnum|res)\s*[:= ]\s*(.+)$", cleaned)
    if residue_match !== nothing
        return residue.seqnum in _selection_int_set(residue_match.captures[1])
    end

    name_match = match(r"^name\s*[:= ]\s*(.+)$", cleaned)
    if name_match !== nothing
        return uppercase(residue.name) in _selection_string_set(name_match.captures[1])
    end

    icode_match = match(r"^(?:icode|insertion|insertion_code)\s*[:= ]\s*(.+)$", cleaned)
    if icode_match !== nothing
        wanted = uppercase(first(strip(icode_match.captures[1])))
        return residue.insertion_code == wanted
    end

    altloc_match = match(r"^altloc\s*[:= ]\s*(.+)$", cleaned)
    if altloc_match !== nothing
        wanted = _selection_string_set(altloc_match.captures[1])
        return any(atom.altloc == ' ' ? "" in wanted : uppercase(string(atom.altloc)) in wanted for atom in residue.atoms)
    end

    occupancy_match = match(r"^occupancy\s*(<=|>=|==|!=|<|>)\s*([0-9.]+)$", cleaned)
    if occupancy_match !== nothing
        threshold = something(tryparse(Float64, occupancy_match.captures[2]), 0.0)
        best = isempty(residue.atoms) ? 0.0 : maximum(atom.occupancy for atom in residue.atoms)
        return _selection_compare(best, occupancy_match.captures[1], threshold)
    end

    shorthand = match(r"^([A-Za-z0-9])\s*:\s*([0-9,\-]+)$", cleaned)
    if shorthand !== nothing
        return uppercase(chain.id) == uppercase(shorthand.captures[1]) && residue.seqnum in _selection_int_set(shorthand.captures[2])
    end

    return residue_property(residue, Symbol(cleaned))
end

function _selection_parse_primary(tokens::Vector{String}, index::Int, chain::Chain, residue::Residue, atom::Union{Nothing,Atom})
    index > length(tokens) && return true, index
    token = tokens[index]
    if token == "("
        value, next_index = _selection_parse_or(tokens, index + 1, chain, residue, atom)
        next_index <= length(tokens) || throw(ArgumentError("unmatched '(' in selection expression"))
        tokens[next_index] == ")" || throw(ArgumentError("unmatched '(' in selection expression"))
        return value, next_index + 1
    end

    parts = String[]
    current_index = index
    while current_index <= length(tokens)
        token = tokens[current_index]
        lower = lowercase(token)
        (token == ")" || lower == "and" || lower == "or") && break
        push!(parts, token)
        current_index += 1
    end
    isempty(parts) && return true, current_index
    return _selection_leaf_matches(chain, residue, join(parts, " "), atom), current_index
end

function _selection_parse_not(tokens::Vector{String}, index::Int, chain::Chain, residue::Residue, atom::Union{Nothing,Atom})
    index > length(tokens) && return true, index
    if lowercase(tokens[index]) == "not"
        value, next_index = _selection_parse_not(tokens, index + 1, chain, residue, atom)
        return !value, next_index
    end
    return _selection_parse_primary(tokens, index, chain, residue, atom)
end

function _selection_parse_and(tokens::Vector{String}, index::Int, chain::Chain, residue::Residue, atom::Union{Nothing,Atom})
    value, next_index = _selection_parse_not(tokens, index, chain, residue, atom)
    while next_index <= length(tokens) && lowercase(tokens[next_index]) == "and"
        rhs, rhs_index = _selection_parse_not(tokens, next_index + 1, chain, residue, atom)
        value = value && rhs
        next_index = rhs_index
    end
    return value, next_index
end

function _selection_parse_or(tokens::Vector{String}, index::Int, chain::Chain, residue::Residue, atom::Union{Nothing,Atom})
    value, next_index = _selection_parse_and(tokens, index, chain, residue, atom)
    while next_index <= length(tokens) && lowercase(tokens[next_index]) == "or"
        rhs, rhs_index = _selection_parse_and(tokens, next_index + 1, chain, residue, atom)
        value = value || rhs
        next_index = rhs_index
    end
    return value, next_index
end

function _selection_matches(chain::Chain, residue::Residue, selector::AbstractString)
    return _selection_matches(chain, residue, selector, nothing)
end

function _selection_matches(chain::Chain, residue::Residue, selector::AbstractString, atom::Union{Nothing,Atom})
    tokens = _selection_tokens(selector)
    value, next_index = _selection_parse_or(tokens, 1, chain, residue, atom)
    next_index <= length(tokens) && throw(ArgumentError("unexpected trailing tokens in selection expression"))
    return value
end

function _selection_int_set(text::AbstractString)
    values = Set{Int}()
    for raw_item in Base.split(String(text), ',')
        item = strip(raw_item)
        isempty(item) && continue
        if occursin('-', item)
            bounds = Base.split(item, '-', limit=2)
            length(bounds) == 2 || continue
            start_value = tryparse(Int, strip(bounds[1]))
            stop_value = tryparse(Int, strip(bounds[2]))
            start_value === nothing && continue
            stop_value === nothing && continue
            for value in start_value:stop_value
                push!(values, value)
            end
        else
            parsed = tryparse(Int, item)
            parsed === nothing || push!(values, parsed)
        end
    end
    return values
end

function _selection_string_set(text::AbstractString)
    values = Set{String}()
    for raw_item in Base.split(String(text), ',')
        item = uppercase(strip(raw_item))
        isempty(item) && continue
        push!(values, item)
    end
    return values
end

function _residue_atom_selection(atoms::AbstractVector{Atom}, policy::AtomSelectionPolicy)
    isempty(atoms) && return Atom[]
    if policy.occupancy === :all || policy.altloc === :all
        return Vector{Atom}(atoms)
    end
    filtered = policy.occupancy === :nonzero ? [atom for atom in atoms if atom.occupancy > 0] : Vector{Atom}(atoms)
    isempty(filtered) && return Atom[]
    if policy.altloc === :first
        return [first(filtered)]
    end
    if policy.altloc === :blank
        for atom in filtered
            atom.altloc == ' ' && return [atom]
        end
    end
    chosen = filtered[1]
    chosen_score = (chosen.occupancy, chosen.altloc == ' ' ? 1 : 0, -Int(chosen.altloc), -chosen.serial)
    for atom in filtered[2:end]
        score = (atom.occupancy, atom.altloc == ' ' ? 1 : 0, -Int(atom.altloc), -atom.serial)
        if score > chosen_score
            chosen = atom
            chosen_score = score
        end
    end
    return [chosen]
end

function collapse_altlocs(residue::Residue; policy::AtomSelectionPolicy=AtomSelectionPolicy())
    collapsed = deepcopy(residue)
    grouped = Dict{String,Vector{Atom}}()
    order = String[]
    for atom in collapsed.atoms
        key = strip(atom.name)
        if !haskey(grouped, key)
            grouped[key] = Atom[]
            push!(order, key)
        end
        push!(grouped[key], atom)
    end
    atoms = Atom[]
    for key in order
        append!(atoms, _residue_atom_selection(grouped[key], policy))
    end
    collapsed.atoms = atoms
    return collapsed
end

function collapse_altlocs(chain::Chain; policy::AtomSelectionPolicy=AtomSelectionPolicy())
    collapsed = deepcopy(chain)
    collapsed.residues = [collapse_altlocs(residue; policy=policy) for residue in collapsed.residues]
    return collapsed
end

function collapse_altlocs(structure::Structure; policy::AtomSelectionPolicy=AtomSelectionPolicy())
    collapsed = deepcopy(structure)
    for model in collapsed.models
        for chain in model.chains
            chain.residues = [collapse_altlocs(residue; policy=policy) for residue in chain.residues]
        end
    end
    return collapsed
end


function _residue_selector_matches(chain::Chain, residue::Residue, selector::AbstractString)
    return _selection_matches(chain, residue, selector, nothing)
end

function select_residues(chain::Chain; selector::Union{Nothing,AbstractString}=nothing, property::Symbol=:hydrophobic)
    selected = Residue[]
    for residue in chain.residues
        if selector === nothing
            if property === :protein
                is_protein_residue(residue) && push!(selected, residue)
            elseif property === :ligand
                is_ligand_residue(residue) && push!(selected, residue)
            elseif property === :water
                is_water_residue(residue) && push!(selected, residue)
            elseif property === :hetero
                is_hetero_residue(residue) && push!(selected, residue)
            else
                residue_property(residue, property) && push!(selected, residue)
            end
        else
            _residue_selector_matches(chain, residue, selector) && push!(selected, residue)
        end
    end
    return selected
end

function select_residues(structure::Structure; selector::Union{Nothing,AbstractString}=nothing, property::Symbol=:hydrophobic, model_index::Int=1, chain_id::Union{Nothing,AbstractString}=nothing)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    selected = Residue[]
    for chain in structure.models[model_index].chains
        chain_id !== nothing && chain.id != String(chain_id) && continue
        append!(selected, select_residues(chain; selector=selector, property=property))
    end
    return selected
end

function select_atoms(structure::Structure; selector::Union{Nothing,AbstractString}=nothing, model_index::Int=1, policy::AtomSelectionPolicy=AtomSelectionPolicy(), chain_id::Union{Nothing,AbstractString}=nothing)
    atoms = Atom[]
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    for chain in structure.models[model_index].chains
        chain_id !== nothing && chain.id != String(chain_id) && continue
        for residue in chain.residues
            selector !== nothing && !_selection_matches(chain, residue, selector, nothing) && continue
            if selector === nothing
                append!(atoms, _residue_atom_selection(residue.atoms, policy))
            else
                for atom in residue.atoms
                    _selection_matches(chain, residue, selector, atom) || continue
                    push!(atoms, atom)
                end
            end
        end
    end
    return atoms
end

function _mmcif_split_row(line::AbstractString)
    tokens = String[]
    current = IOBuffer()
    in_quote = false
    quote_char = '\0'
    for char in line
        if in_quote
            if char == quote_char
                push!(tokens, String(take!(current)))
                in_quote = false
                quote_char = '\0'
            else
                write(current, char)
            end
        elseif isspace(char)
            if position(current) > 0
                push!(tokens, String(take!(current)))
            end
        elseif char == '"' || char == '\''
            if position(current) > 0
                push!(tokens, String(take!(current)))
            end
            in_quote = true
            quote_char = char
        else
            write(current, char)
        end
    end
    if position(current) > 0
        push!(tokens, String(take!(current)))
    end
    return tokens
end

const _VDW_RADIUS = Dict{String,Float64}(
    "H" => 1.20,
    "C" => 1.70,
    "N" => 1.55,
    "O" => 1.52,
    "F" => 1.47,
    "P" => 1.80,
    "S" => 1.80,
    "CL" => 1.75,
    "BR" => 1.85,
    "SE" => 1.90,
    "I" => 1.98,
    "NA" => 2.27,
    "MG" => 1.73,
    "CA" => 2.31,
    "ZN" => 1.39,
    "FE" => 1.56,
    "CU" => 1.40,
)
const _MAX_VDW_RADIUS = maximum(values(_VDW_RADIUS))
const _FIBONACCI_CACHE = Dict{Int,Vector{NTuple{3,Float64}}}()

function _atom_vdw_radius(atom::Atom)
    element = _atom_element_symbol(atom)
    if isempty(element)
        cleaned = filter(isletter, atom.name)
        isempty(cleaned) && return 1.70
        element = uppercase(string(cleaned[1]))
    end
    return get(_VDW_RADIUS, element, 1.70)
end

function _fibonacci_sphere_points(samples::Int)
    samples > 0 || throw(ArgumentError("sample count must be positive"))
    return get!(_FIBONACCI_CACHE, samples) do
        points = Vector{NTuple{3,Float64}}(undef, samples)
        golden_angle = pi * (3 - sqrt(5.0))
        for index in 1:samples
            y = 1 - 2 * (index - 0.5) / samples
            radius = sqrt(max(0.0, 1 - y * y))
            theta = golden_angle * index
            points[index] = (cos(theta) * radius, y, sin(theta) * radius)
        end
        points
    end
end

function _model_atoms(structure::Structure, model_index::Int)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    atoms = Atom[]
    for chain in structure.models[model_index].chains
        for residue in chain.residues
            append!(atoms, residue.atoms)
        end
    end
    return atoms
end

function _atom_sasa(atom::Atom, atoms::AbstractVector{Atom}, points::AbstractVector{NTuple{3,Float64}}; probe::Real=1.4)
    radius = _atom_vdw_radius(atom) + Float64(probe)
    candidate_radii = Float64[_atom_vdw_radius(other) + Float64(probe) for other in atoms]
    accessible = 0
    for point in points
        px = atom.x + point[1] * radius
        py = atom.y + point[2] * radius
        pz = atom.z + point[3] * radius
        blocked = false
        for (other_index, other) in pairs(atoms)
            other === atom && continue
            other_radius = candidate_radii[other_index]
            dx = px - other.x
            dy = py - other.y
            dz = pz - other.z
            if dx * dx + dy * dy + dz * dz <= other_radius * other_radius
                blocked = true
                break
            end
        end
        accessible += blocked ? 0 : 1
    end
    return 4 * pi * radius * radius * (accessible / length(points))
end

function _atom_sasa(atom::Atom, atoms::AbstractVector{Atom}; probe::Real=1.4, samples::Int=96)
    return _atom_sasa(atom, atoms, _fibonacci_sphere_points(samples); probe=probe)
end

function _atom_sasa(atom::Atom, atoms::AbstractVector{Atom}, tree::AtomKDTree, points::AbstractVector{NTuple{3,Float64}}; probe::Real=1.4)
    search_radius = _atom_vdw_radius(atom) + Float64(probe) + Float64(probe) + _MAX_VDW_RADIUS
    candidates = atoms_within_radius(tree, atom; radius=search_radius)
    return _atom_sasa(atom, candidates, points; probe=probe)
end

function atom_sasa(atom::Atom, structure::Structure; probe::Real=1.4, samples::Int=96, model_index::Int=1)
    atoms = _model_atoms(structure, model_index)
    isempty(atoms) && return 0.0
    tree = build_atom_kdtree(atoms)
    points = _fibonacci_sphere_points(samples)
    return _atom_sasa(atom, atoms, tree, points; probe=probe)
end

function residue_sasa(residue::Residue, structure::Structure; probe::Real=1.4, samples::Int=96, model_index::Int=1)
    isempty(residue.atoms) && return 0.0
    atoms = _model_atoms(structure, model_index)
    tree = build_atom_kdtree(atoms)
    points = _fibonacci_sphere_points(samples)
    return sum(_atom_sasa(atom, atoms, tree, points; probe=probe) for atom in residue.atoms)
end

function chain_sasa(chain::Chain, structure::Structure; probe::Real=1.4, samples::Int=96, model_index::Int=1)
    atoms = _model_atoms(structure, model_index)
    tree = build_atom_kdtree(atoms)
    points = _fibonacci_sphere_points(samples)
    return sum(sum(_atom_sasa(atom, atoms, tree, points; probe=probe) for atom in residue.atoms) for residue in chain.residues)
end

function structure_sasa(structure::Structure; probe::Real=1.4, samples::Int=96, model_index::Int=1)
    atoms = _model_atoms(structure, model_index)
    isempty(atoms) && return 0.0
    tree = build_atom_kdtree(atoms)
    points = _fibonacci_sphere_points(samples)
    return sum(_atom_sasa(atom, atoms, tree, points; probe=probe) for atom in atoms)
end

function sasa_profile(structure::Structure; probe::Real=1.4, samples::Int=96, model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    atoms = _model_atoms(structure, model_index)
    tree = build_atom_kdtree(atoms)
    points = _fibonacci_sphere_points(samples)
    entries = NamedTuple{(:chain, :residue, :seqnum, :sasa, :relative_accessibility, :kind), Tuple{String, String, Int, Float64, Float64, Symbol}}[]
    for chain in structure.models[model_index].chains
        for residue in chain.residues
            sasa = sum(_atom_sasa(atom, atoms, tree, points; probe=probe) for atom in residue.atoms)
            max_sasa = sum(4 * pi * (_atom_vdw_radius(atom) + Float64(probe))^2 for atom in residue.atoms)
            kind = is_water_residue(residue) ? :water : (is_ligand_residue(residue) ? :ligand : :protein)
            push!(entries, (
                chain=chain.id,
                residue=residue.name,
                seqnum=residue.seqnum,
                sasa=sasa,
                relative_accessibility=max_sasa == 0 ? 0.0 : sasa / max_sasa,
                kind=kind,
            ))
        end
    end
    return entries
end

function _isolated_chain_structure(chain::Chain)
    isolated = Structure("chain")
    model = Model(1)
    push!(model.chains, deepcopy(chain))
    push!(isolated.models, model)
    return isolated
end

function residue_free_sasa(residue::Residue, chain::Chain; probe::Real=1.4, samples::Int=96)
    atoms = Atom[]
    for chain_residue in chain.residues
        append!(atoms, chain_residue.atoms)
    end
    tree = build_atom_kdtree(atoms)
    points = _fibonacci_sphere_points(samples)
    return sum(_atom_sasa(atom, atoms, tree, points; probe=probe) for atom in residue.atoms)
end

function buried_surface_area(residue::Residue, chain::Chain, structure::Structure; probe::Real=1.4, samples::Int=96, model_index::Int=1)
    free_sasa = residue_free_sasa(residue, chain; probe=probe, samples=samples)
    complex_sasa = residue_sasa(residue, structure; probe=probe, samples=samples, model_index=model_index)
    return max(0.0, free_sasa - complex_sasa)
end

function buried_surface_area(chain::Chain, structure::Structure; probe::Real=1.4, samples::Int=96, model_index::Int=1)
    return sum(buried_surface_area(residue, chain, structure; probe=probe, samples=samples, model_index=model_index) for residue in chain.residues)
end

function buried_surface_area(structure::Structure; probe::Real=1.4, samples::Int=96, model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    total = 0.0
    for chain in structure.models[model_index].chains
        total += buried_surface_area(chain, structure; probe=probe, samples=samples, model_index=model_index)
    end
    return total
end

function _residue_contact_counts(structure::Structure; model_index::Int=1, cutoff::Real=8.0)
    counts = Dict{Tuple{String,Int,Char},Int}()
    for contact in residue_contacts(structure; cutoff=cutoff, model_index=model_index)
        left_key = (contact.left_chain, contact.left_seqnum, contact.left_insertion_code)
        right_key = (contact.right_chain, contact.right_seqnum, contact.right_insertion_code)
        counts[left_key] = get(counts, left_key, 0) + 1
        counts[right_key] = get(counts, right_key, 0) + 1
    end
    return counts
end

function interface_profile(structure::Structure; probe::Real=1.4, samples::Int=96, model_index::Int=1, cutoff::Real=8.0, minimum_buried_sasa::Real=1.0)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    model = structure.models[model_index]
    contact_counts = _residue_contact_counts(structure; model_index=model_index, cutoff=cutoff)
    complex_atoms = _model_atoms(structure, model_index)
    complex_tree = build_atom_kdtree(complex_atoms)
    points = _fibonacci_sphere_points(samples)
    profile = NamedTuple{(:chain, :residue, :seqnum, :insertion_code, :complex_sasa, :free_sasa, :buried_sasa, :contact_count, :is_interface), Tuple{String, String, Int, Char, Float64, Float64, Float64, Int, Bool}}[]
    for chain in model.chains
        isolated_atoms = Atom[]
        for chain_residue in chain.residues
            append!(isolated_atoms, chain_residue.atoms)
        end
        isolated_tree = build_atom_kdtree(isolated_atoms)
        for residue in chain.residues
            free_sasa = sum(_atom_sasa(atom, isolated_atoms, isolated_tree, points; probe=probe) for atom in residue.atoms)
            complex_sasa = sum(_atom_sasa(atom, complex_atoms, complex_tree, points; probe=probe) for atom in residue.atoms)
            buried_sasa = max(0.0, free_sasa - complex_sasa)
            contact_count = get(contact_counts, (chain.id, residue.seqnum, residue.insertion_code), 0)
            is_interface = contact_count > 0 || buried_sasa > minimum_buried_sasa
            push!(profile, (
                chain=chain.id,
                residue=residue.name,
                seqnum=residue.seqnum,
                insertion_code=residue.insertion_code,
                complex_sasa=complex_sasa,
                free_sasa=free_sasa,
                buried_sasa=buried_sasa,
                contact_count=contact_count,
                is_interface=is_interface,
            ))
        end
    end
    return profile
end

function disulfide_bonds(structure::Structure; model_index::Int=1, cutoff::Real=2.2)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    model = structure.models[model_index]
    cysteines = NamedTuple{(:chain, :residue), Tuple{Chain, Residue}}[]
    for chain in model.chains
        for residue in chain.residues
            uppercase(residue.name) == "CYS" || continue
            any(atom.name == "SG" for atom in residue.atoms) || continue
            push!(cysteines, (chain=chain, residue=residue))
        end
    end
    pairs = NamedTuple{(:left_chain, :left_residue, :left_seqnum, :left_insertion_code, :right_chain, :right_residue, :right_seqnum, :right_insertion_code, :distance), Tuple{String, String, Int, Char, String, String, Int, Char, Float64}}[]
    cutoff2 = Float64(cutoff)^2
    for left_index in firstindex(cysteines):(lastindex(cysteines) - 1)
        left_item = cysteines[left_index]
        left_atom = first(atom for atom in left_item.residue.atoms if atom.name == "SG")
        for right_index in (left_index + 1):lastindex(cysteines)
            right_item = cysteines[right_index]
            right_atom = first(atom for atom in right_item.residue.atoms if atom.name == "SG")
            distance2 = _distance2((left_atom.x, left_atom.y, left_atom.z), (right_atom.x, right_atom.y, right_atom.z))
            distance2 <= cutoff2 || continue
            push!(pairs, (
                left_chain=left_item.chain.id,
                left_residue=left_item.residue.name,
                left_seqnum=left_item.residue.seqnum,
                left_insertion_code=left_item.residue.insertion_code,
                right_chain=right_item.chain.id,
                right_residue=right_item.residue.name,
                right_seqnum=right_item.residue.seqnum,
                right_insertion_code=right_item.residue.insertion_code,
                distance=sqrt(distance2),
            ))
        end
    end
    return pairs
end

function _residue_center(residue::Residue)
    isempty(residue.atoms) && return nothing
    x = 0.0
    y = 0.0
    z = 0.0
    for atom in residue.atoms
        x += atom.x
        y += atom.y
        z += atom.z
    end
    invn = 1 / length(residue.atoms)
    return (x * invn, y * invn, z * invn)
end

function _distance2(point_a, point_b)
    dx = point_a[1] - point_b[1]
    dy = point_a[2] - point_b[2]
    dz = point_a[3] - point_b[3]
    return dx * dx + dy * dy + dz * dz
end

function _residue_contact_pairs(residues::AbstractVector{Residue}; cutoff::Real=8.0)
    centers = [_residue_center(residue) for residue in residues]
    bucket_size = Float64(cutoff)
    bucket_index(center) = (floor(Int, center[1] / bucket_size), floor(Int, center[2] / bucket_size), floor(Int, center[3] / bucket_size))
    buckets = Dict{NTuple{3,Int},Vector{Int}}()
    for (index, center) in enumerate(centers)
        center === nothing && continue
        push!(get!(buckets, bucket_index(center), Int[]), index)
    end

    pairs = NamedTuple{(:left, :right, :distance), Tuple{Int, Int, Float64}}[]
    seen = Set{Tuple{Int,Int}}()
    cutoff2 = bucket_size^2
    for (left_index, left_center) in enumerate(centers)
        left_center === nothing && continue
        left_bucket = bucket_index(left_center)
        for dx in -1:1, dy in -1:1, dz in -1:1
            neighbor_indices = get(buckets, (left_bucket[1] + dx, left_bucket[2] + dy, left_bucket[3] + dz), nothing)
            neighbor_indices === nothing && continue
            for right_index in neighbor_indices
                right_index <= left_index && continue
                pair = (left_index, right_index)
                pair in seen && continue
                right_center = centers[right_index]
                right_center === nothing && continue
                distance2 = _distance2(left_center, right_center)
                distance2 <= cutoff2 || continue
                push!(pairs, (left=left_index, right=right_index, distance=sqrt(distance2)))
                push!(seen, pair)
            end
        end
    end
    return pairs
end

function residue_contacts(chain::Chain; cutoff::Real=8.0)
    return _residue_contact_pairs(chain.residues; cutoff=cutoff)
end

function residue_contacts(structure::Structure; cutoff::Real=8.0, model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    model = structure.models[model_index]
    residues = NamedTuple{(:chain, :residue), Tuple{Chain, Residue}}[]
    for chain in model.chains
        for residue in chain.residues
            push!(residues, (chain=chain, residue=residue))
        end
    end

    local_pairs = _residue_contact_pairs([item.residue for item in residues]; cutoff=cutoff)
    pairs = NamedTuple{(:left_chain, :left_residue, :left_seqnum, :left_insertion_code, :right_chain, :right_residue, :right_seqnum, :right_insertion_code, :distance), Tuple{String, String, Int, Char, String, String, Int, Char, Float64}}[]
    for pair in local_pairs
        left_item = residues[pair.left]
        right_item = residues[pair.right]
        push!(pairs, (
            left_chain=left_item.chain.id,
            left_residue=left_item.residue.name,
            left_seqnum=left_item.residue.seqnum,
            left_insertion_code=left_item.residue.insertion_code,
            right_chain=right_item.chain.id,
            right_residue=right_item.residue.name,
            right_seqnum=right_item.residue.seqnum,
            right_insertion_code=right_item.residue.insertion_code,
            distance=pair.distance,
        ))
    end
    return pairs
end

function contact_map(chain::Chain; cutoff::Real=8.0)
    matrix = falses(length(chain.residues), length(chain.residues))
    for pair in residue_contacts(chain; cutoff=cutoff)
        matrix[pair.left, pair.right] = true
        matrix[pair.right, pair.left] = true
    end
    return matrix
end

function contact_map(structure::Structure; cutoff::Real=8.0, model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    model = structure.models[model_index]
    labels = String[]
    for chain in model.chains
        for residue in chain.residues
            push!(labels, string(chain.id, ":", residue.name, residue.seqnum, residue.insertion_code))
        end
    end

    matrix = falses(length(labels), length(labels))
    residues = residue_contacts(structure; cutoff=cutoff, model_index=model_index)
    lookup = Dict{Tuple{String,String,Int,Char},Int}()
    index = 1
    for chain in model.chains
        for residue in chain.residues
            lookup[(chain.id, residue.name, residue.seqnum, residue.insertion_code)] = index
            index += 1
        end
    end

    for contact in residues
        left_index = lookup[(contact.left_chain, contact.left_residue, contact.left_seqnum, contact.left_insertion_code)]
        right_index = lookup[(contact.right_chain, contact.right_residue, contact.right_seqnum, contact.right_insertion_code)]
        matrix[left_index, right_index] = true
        matrix[right_index, left_index] = true
    end

    return (matrix=matrix, labels=labels)
end

function interface_residues(structure::Structure; cutoff::Real=8.0, model_index::Int=1, minimum_buried_sasa::Real=1.0)
    profile = interface_profile(structure; cutoff=cutoff, model_index=model_index, minimum_buried_sasa=minimum_buried_sasa)
    marked = Set{Tuple{String,Int,Char}}()
    residues = Residue[]
    model = structure.models[model_index]
    for entry in profile
        entry.is_interface || continue
        push!(marked, (entry.chain, entry.seqnum, entry.insertion_code))
    end
    for chain in model.chains
        for residue in chain.residues
            (chain.id, residue.seqnum, residue.insertion_code) in marked && push!(residues, residue)
        end
    end
    return residues
end

function calculate_interface_residues(chain_a::Chain, chain_b::Chain; probe::Real=1.4, samples::Int=96, cutoff::Real=8.0, minimum_buried_sasa::Real=1.0)
    structure = Structure("interface", [Model(1, [chain_a, chain_b])], Dict{String,String}())
    profile = interface_profile(structure; probe=probe, samples=samples, model_index=1, cutoff=cutoff, minimum_buried_sasa=minimum_buried_sasa)
    return [entry for entry in profile if entry.is_interface]
end

function residues_within_radius(chain::Chain, ligand::Residue; radius::Real=5.0)
    ligand_atoms = ligand.atoms
    isempty(ligand_atoms) && return Residue[]
    radius2 = Float64(radius)^2
    hits = Residue[]

    for residue in chain.residues
        residue === ligand && continue
        for atom in residue.atoms
            any(ligand_atom -> _distance2((atom.x, atom.y, atom.z), (ligand_atom.x, ligand_atom.y, ligand_atom.z)) <= radius2, ligand_atoms) && (push!(hits, residue); break)
        end
    end

    return hits
end

function residues_within_radius(structure::Structure, ligand::Residue; radius::Real=5.0, model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    hits = Residue[]
    for chain in structure.models[model_index].chains
        append!(hits, residues_within_radius(chain, ligand; radius=radius))
    end
    return hits
end

function _dssp_char(value::AbstractString, default::Char=' ')
    isempty(value) && return default
    return first(value)
end

function read_dssp(filepath::AbstractString)
    open(filepath, "r") do io
        return read_dssp(io)
    end
end

function read_dssp(io::IO)
    entries = DSSPEntry[]
    in_table = false
    for raw_line in eachline(io)
        line = rstrip(raw_line)
        if !in_table
            occursin("RESIDUE", line) && occursin("STRUCTURE", line) && (in_table = true)
            continue
        end
        isempty(strip(line)) && continue
        startswith(strip(line), "#") && continue

        seqnum = tryparse(Int, _pdb_field(line, 6, 10))
        seqnum === nothing && continue
        chain_id = _pdb_field(line, 12, 12)
        amino_acid = _dssp_char(_pdb_field(line, 14, 14), 'X')
        secondary_structure = _dssp_char(_pdb_field(line, 17, 17), ' ')
        accessibility_text = _pdb_field(line, 36, 38)
        accessibility = something(tryparse(Float64, accessibility_text), 0.0)
        insertion_code = _dssp_char(_pdb_field(line, 11, 11), ' ')
        push!(entries, DSSPEntry(chain_id, seqnum, insertion_code, amino_acid, secondary_structure, accessibility))
    end
    return entries
end

function annotate_dssp!(structure::Structure, entries::AbstractVector{DSSPEntry}; model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    residues = Dict{Tuple{String,Int,Char},Residue}()
    for chain in structure.models[model_index].chains
        for residue in chain.residues
            residues[(chain.id, residue.seqnum, residue.insertion_code)] = residue
        end
    end
    for entry in entries
        residue = get(residues, (entry.chain, entry.seqnum, entry.insertion_code), nothing)
        residue === nothing && continue
        residue.secondary_structure = entry.secondary_structure == ' ' ? missing : entry.secondary_structure
        residue.accessibility = entry.accessibility
    end
    return structure
end

function _atom_name(atom::Atom)
    return uppercase(strip(atom.name))
end

function _residue_atom(residue::Residue, names::Vararg{AbstractString})
    target_names = Set(uppercase(strip(String(name))) for name in names)
    for atom in residue.atoms
        _atom_name(atom) in target_names && return atom
    end
    return nothing
end

function _residue_center_atom(residue::Residue, names::Vararg{AbstractString})
    atom = _residue_atom(residue, names...)
    atom === nothing && return nothing
    return atom
end

function _vector_between(atom_a::Atom, atom_b::Atom)
    return (atom_b.x - atom_a.x, atom_b.y - atom_a.y, atom_b.z - atom_a.z)
end

function _vector_norm(vector)
    return sqrt(sum(value^2 for value in vector))
end

function _angle_degrees(atom_a::Atom, atom_b::Atom, atom_c::Atom)
    vector_ab = _vector_between(atom_b, atom_a)
    vector_cb = _vector_between(atom_b, atom_c)
    norm_ab = _vector_norm(vector_ab)
    norm_cb = _vector_norm(vector_cb)
    norm_ab == 0 && return NaN
    norm_cb == 0 && return NaN
    cosine = clamp((vector_ab[1] * vector_cb[1] + vector_ab[2] * vector_cb[2] + vector_ab[3] * vector_cb[3]) / (norm_ab * norm_cb), -1.0, 1.0)
    return acos(cosine) * 180 / π
end

function hydrogen_bonds(structure::Structure; model_index::Int=1, cutoff::Real=3.5, angle_cutoff::Real=120)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    model = structure.models[model_index]
    residues = Vector{NamedTuple{(:chain, :residue), Tuple{Chain, Residue}}}()
    for chain in model.chains
        for residue in chain.residues
            push!(residues, (chain=chain, residue=residue))
        end
    end

    bonds = HydrogenBond[]
    cutoff2 = Float64(cutoff)^2
    for donor_index in eachindex(residues)
        donor_item = residues[donor_index]
        donor_n = _residue_atom(donor_item.residue, "N")
        donor_c = _residue_atom(donor_item.residue, "C")
        donor_n === nothing && continue
        donor_c === nothing && continue
        for acceptor_index in eachindex(residues)
            donor_index == acceptor_index && continue
            acceptor_item = residues[acceptor_index]
            acceptor_o = _residue_atom(acceptor_item.residue, "O", "OXT")
            acceptor_c = _residue_atom(acceptor_item.residue, "C")
            acceptor_o === nothing && continue
            acceptor_c === nothing && continue
            distance2 = _distance2(atom_coordinates(donor_n), atom_coordinates(acceptor_o))
            distance2 <= cutoff2 || continue
            donor_angle = _angle_degrees(donor_c, donor_n, acceptor_o)
            acceptor_angle = _angle_degrees(donor_n, acceptor_o, acceptor_c)
            if donor_angle >= angle_cutoff && acceptor_angle >= angle_cutoff
                push!(bonds, HydrogenBond(
                    donor_item.chain.id,
                    donor_item.residue.name,
                    donor_item.residue.seqnum,
                    donor_item.residue.insertion_code,
                    donor_n.name,
                    acceptor_item.chain.id,
                    acceptor_item.residue.name,
                    acceptor_item.residue.seqnum,
                    acceptor_item.residue.insertion_code,
                    acceptor_o.name,
                    sqrt(distance2),
                    min(donor_angle, acceptor_angle),
                ))
            end
        end
    end
    return bonds
end

function hydrogen_bonds(chain::Chain; cutoff::Real=3.5, angle_cutoff::Real=120)
    structure = Structure("chain")
    model = Model(1)
    push!(model.chains, chain)
    push!(structure.models, model)
    return hydrogen_bonds(structure; model_index=1, cutoff=cutoff, angle_cutoff=angle_cutoff)
end

function ramachandran_region(phi::Real, psi::Real; residue_name::AbstractString="ALA")
    isfinite(phi) && isfinite(psi) || return :unknown
    residue = uppercase(String(residue_name))
    if residue == "PRO"
        if phi >= -105 && phi <= -35 && psi >= -80 && psi <= 180
            return :proline_preferred
        end
    elseif residue == "GLY"
        if phi >= -180 && phi <= 180 && psi >= -180 && psi <= 180
            if phi >= 0 && psi >= -40 && psi <= 180
                return :glycine_allowed
            end
        end
    end

    if phi >= -160 && phi <= -30 && psi >= -90 && psi <= 45
        return :alpha_helix
    end
    if phi >= -180 && phi <= -40 && psi >= 90 && psi <= 180
        return :beta_sheet
    end
    if phi >= 30 && phi <= 100 && psi >= -90 && psi <= 90
        return :left_handed_helix
    end
    if residue == "GLY"
        return :glycine_allowed
    end
    return :loop
end

function ramachandran_profile(chain::Chain)
    profile = NamedTuple{(:index, :residue, :phi, :psi, :region), Tuple{Int, String, Union{Missing,Float64}, Union{Missing,Float64}, Symbol}}[]
    for index in eachindex(chain.residues)
        residue = chain.residues[index]
        torsions = phi_psi(chain, index)
        region = ismissing(torsions.phi) || ismissing(torsions.psi) ? :unknown : ramachandran_region(torsions.phi, torsions.psi; residue_name=residue.name)
        push!(profile, (index=index, residue=residue.name, phi=torsions.phi, psi=torsions.psi, region=region))
    end
    return profile
end

function ramachandran_profile(structure::Structure; model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    profile = NamedTuple{(:chain, :index, :residue, :phi, :psi, :region), Tuple{String, Int, String, Union{Missing,Float64}, Union{Missing,Float64}, Symbol}}[]
    for chain in structure.models[model_index].chains
        for entry in ramachandran_profile(chain)
            push!(profile, (chain=chain.id, index=entry.index, residue=entry.residue, phi=entry.phi, psi=entry.psi, region=entry.region))
        end
    end
    return profile
end

const _CHI_DEFINITIONS = Dict{String,Vector{NTuple{4,String}}}(
    "ARG" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "NE"), ("CG", "CD", "NE", "CZ")],
    "ASN" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "ASP" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "CYS" => [("N", "CA", "CB", "SG")],
    "GLN" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
    "GLU" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
    "HIS" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")],
    "ILE" => [("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")],
    "LEU" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "LYS" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "CE"), ("CG", "CD", "CE", "NZ")],
    "MET" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "SD"), ("CB", "CG", "SD", "CE")],
    "PHE" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "PRO" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")],
    "SER" => [("N", "CA", "CB", "OG")],
    "THR" => [("N", "CA", "CB", "OG1")],
    "TRP" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "TYR" => [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "VAL" => [("N", "CA", "CB", "CG1")],
)

function _chi_angle(residue::Residue, atoms::NTuple{4,String})
    atom1 = _residue_atom(residue, atoms[1])
    atom2 = _residue_atom(residue, atoms[2])
    atom3 = _residue_atom(residue, atoms[3])
    atom4 = _residue_atom(residue, atoms[4])
    if atom1 === nothing || atom2 === nothing || atom3 === nothing || atom4 === nothing
        return missing
    end
    return torsion_angle(atom_coordinates(atom1), atom_coordinates(atom2), atom_coordinates(atom3), atom_coordinates(atom4))
end

function chi_angles(residue::Residue)
    definitions = get(_CHI_DEFINITIONS, uppercase(residue.name), NTuple{4,String}[])
    chi_values = Union{Missing,Float64}[missing, missing, missing, missing]
    for (index, atoms) in enumerate(definitions)
        index > 4 && break
        chi_values[index] = _chi_angle(residue, atoms)
    end
    return (chi1=chi_values[1], chi2=chi_values[2], chi3=chi_values[3], chi4=chi_values[4])
end

function rotamer_label(angle::Real)
    isfinite(angle) || return :unknown
    normalized = mod(angle + 180, 360) - 180
    abs(abs(normalized) - 180) <= 60 && return :trans
    normalized >= 0 && normalized <= 120 && return :gauche_plus
    normalized < 0 && normalized >= -120 && return :gauche_minus
    return :other
end

function rotamer_state(residue::Residue)
    chis = chi_angles(residue)
    return (
        chi1=chis.chi1,
        chi2=chis.chi2,
        chi3=chis.chi3,
        chi4=chis.chi4,
        chi1_label=ismissing(chis.chi1) ? missing : rotamer_label(chis.chi1),
        chi2_label=ismissing(chis.chi2) ? missing : rotamer_label(chis.chi2),
        chi3_label=ismissing(chis.chi3) ? missing : rotamer_label(chis.chi3),
        chi4_label=ismissing(chis.chi4) ? missing : rotamer_label(chis.chi4),
    )
end

function rotamer_statistics(chain::Chain)
    counts = Dict{Symbol,Int}(
        :trans => 0,
        :gauche_plus => 0,
        :gauche_minus => 0,
        :other => 0,
        :unknown => 0,
    )
    total = 0
    with_chi1 = 0
    for residue in chain.residues
        state = rotamer_state(residue)
        ismissing(state.chi1_label) && begin
            counts[:unknown] += 1
            total += 1
            continue
        end
        with_chi1 += 1
        counts[state.chi1_label] = get(counts, state.chi1_label, 0) + 1
        total += 1
    end
    fractions = Dict(symbol => (total == 0 ? 0.0 : count / total) for (symbol, count) in counts)
    return (total=total, with_chi1=with_chi1, counts=counts, fractions=fractions)
end

function rotamer_statistics(structure::Structure; model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    chain_stats = Dict{String,NamedTuple}()
    overall_counts = Dict{Symbol,Int}(
        :trans => 0,
        :gauche_plus => 0,
        :gauche_minus => 0,
        :other => 0,
        :unknown => 0,
    )
    total = 0
    with_chi1 = 0
    for chain in structure.models[model_index].chains
        stats = rotamer_statistics(chain)
        chain_stats[chain.id] = stats
        total += stats.total
        with_chi1 += stats.with_chi1
        for (symbol, count) in stats.counts
            overall_counts[symbol] = get(overall_counts, symbol, 0) + count
        end
    end
    fractions = Dict(symbol => (total == 0 ? 0.0 : count / total) for (symbol, count) in overall_counts)
    return (total=total, with_chi1=with_chi1, counts=overall_counts, fractions=fractions, chains=chain_stats)
end

function _residue_node_label(chain::String, residue::Residue)
    ss = ismissing(residue.secondary_structure) ? "?" : string(residue.secondary_structure)
    acc = ismissing(residue.accessibility) ? "" : @sprintf(" %.1f", residue.accessibility)
    return string(chain, ":", residue.name, residue.seqnum, residue.insertion_code, "<br/>", ss, acc)
end

function _structure_node_label(chain::String, residue::Residue)
    ss = ismissing(residue.secondary_structure) ? "?" : string(residue.secondary_structure)
    return string(chain, ":", residue.name, residue.seqnum, residue.insertion_code, "\n", ss)
end

function _residue_graph_nodes(structure::Structure; model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    model = structure.models[model_index]
    nodes = Vector{NamedTuple{(:id, :label, :chain, :residue, :seqnum, :insertion_code, :secondary_structure), Tuple{String, String, String, String, Int, Char, Union{Missing,Char}}}}()
    lookup = Dict{Tuple{String,Int,Char},String}()
    counter = 1
    for chain in model.chains
        for residue in chain.residues
            node_id = string("r", counter)
            push!(nodes, (id=node_id, label=_structure_node_label(chain.id, residue), chain=chain.id, residue=residue.name, seqnum=residue.seqnum, insertion_code=residue.insertion_code, secondary_structure=residue.secondary_structure))
            lookup[(chain.id, residue.seqnum, residue.insertion_code)] = node_id
            counter += 1
        end
    end
    return nodes, lookup, model
end

function _graph_layout(nnodes::Int)
    if nnodes == 0
        return Vector{Tuple{Float64,Float64}}()
    end
    if nnodes == 1
        return [(0.0, 0.0)]
    end
    layout = Vector{Tuple{Float64,Float64}}(undef, nnodes)
    radius = 1.0 + nnodes / 6
    for index in 1:nnodes
        angle = 2π * (index - 1) / nnodes
        layout[index] = (radius * cos(angle), radius * sin(angle))
    end
    return layout
end

function structure_contact_graph(structure::Structure; model_index::Int=1, cutoff::Real=8.0, hbond_cutoff::Real=3.5)
    nodes, lookup, model = _residue_graph_nodes(structure; model_index=model_index)
    layout = _graph_layout(length(nodes))
    generic_edges = NamedTuple{(:source, :target, :distance), Tuple{String, String, Float64}}[]
    hbond_edges = NamedTuple{(:source, :target, :distance, :angle), Tuple{String, String, Float64, Float64}}[]

    for contact in residue_contacts(structure; cutoff=cutoff, model_index=model_index)
        source = lookup[(contact.left_chain, contact.left_seqnum, contact.left_insertion_code)]
        target = lookup[(contact.right_chain, contact.right_seqnum, contact.right_insertion_code)]
        push!(generic_edges, (source=source, target=target, distance=contact.distance))
    end

    for bond in hydrogen_bonds(structure; model_index=model_index, cutoff=hbond_cutoff)
        source = lookup[(bond.donor_chain, bond.donor_seqnum, bond.donor_insertion_code)]
        target = lookup[(bond.acceptor_chain, bond.acceptor_seqnum, bond.acceptor_insertion_code)]
        push!(hbond_edges, (source=source, target=target, distance=bond.distance, angle=bond.angle))
    end

    return (nodes=nodes, layout=layout, generic_edges=generic_edges, hydrogen_bond_edges=hbond_edges, model=model)
end

function _structure_secondary_structure_color(value::Union{Missing,Char})
    ismissing(value) && return "#64748b"
    value == 'H' && return "#ef4444"
    value == 'E' && return "#2563eb"
    value == 'G' && return "#f97316"
    value == 'I' && return "#8b5cf6"
    value == 'B' && return "#06b6d4"
    value == 'T' && return "#10b981"
    value == 'S' && return "#84cc16"
    return "#64748b"
end

function plot_contact_graph!(axis, structure::Structure; model_index::Int=1, cutoff::Real=8.0, hbond_cutoff::Real=3.5, node_size::Real=22)
    graph = structure_contact_graph(structure; model_index=model_index, cutoff=cutoff, hbond_cutoff=hbond_cutoff)
    ids = [node.id for node in graph.nodes]
    labels = [node.label for node in graph.nodes]
    positions = Dict{String,Tuple{Float64,Float64}}(id => graph.layout[index] for (index, id) in enumerate(ids))
    xs = [positions[id][1] for id in ids]
    ys = [positions[id][2] for id in ids]
    node_colors = [_structure_secondary_structure_color(node.secondary_structure) for node in graph.nodes]

    for edge in graph.generic_edges
        p1 = positions[edge.source]
        p2 = positions[edge.target]
        Makie.lines!(axis, [p1[1], p2[1]], [p1[2], p2[2]]; color=Makie.RGBAf(0.62, 0.62, 0.62, 0.65), linewidth=2)
    end

    for edge in graph.hydrogen_bond_edges
        p1 = positions[edge.source]
        p2 = positions[edge.target]
        Makie.lines!(axis, [p1[1], p2[1]], [p1[2], p2[2]]; color=Makie.RGBAf(0.93, 0.27, 0.27, 0.9), linewidth=2)
    end

    Makie.scatter!(axis, xs, ys; markersize=node_size, color=node_colors, strokecolor=:white, strokewidth=1.5)
    Makie.text!(axis, xs, ys; text=labels, align=(:center, :center), fontsize=10, color=:black)
    axis.aspect = Makie.DataAspect()
    return graph
end

function plot_contact_graph(structure::Structure; kwargs...)
    fig = Makie.Figure()
    axis = Makie.Axis(fig[1, 1]; title="Residue contact graph", xlabel="graph x", ylabel="graph y")
    plot_contact_graph!(axis, structure; kwargs...)
    return fig
end

function plot_contact_map!(axis, structure::Structure; model_index::Int=1, cutoff::Real=8.0, hbond_cutoff::Real=3.5)
    graph = structure_contact_graph(structure; model_index=model_index, cutoff=cutoff, hbond_cutoff=hbond_cutoff)
    n = length(graph.nodes)
    values = zeros(Int, n, n)
    ids = [node.id for node in graph.nodes]
    lookup = Dict{String,Int}(id => index for (index, id) in enumerate(ids))

    for edge in graph.generic_edges
        left = lookup[edge.source]
        right = lookup[edge.target]
        values[left, right] = max(values[left, right], 1)
        values[right, left] = max(values[right, left], 1)
    end

    for edge in graph.hydrogen_bond_edges
        left = lookup[edge.source]
        right = lookup[edge.target]
        values[left, right] = 2
        values[right, left] = 2
    end

    labels = [node.label for node in graph.nodes]
    if n == 0
        values = zeros(Int, 1, 1)
        labels = ["empty"]
    end

    Makie.heatmap!(axis, values; colormap=["#fffdf8", "#9ca3af", "#ef4444"], colorrange=(0, 2), interpolate=false)
    axis.xticks = (1:length(labels), labels)
    axis.yticks = (1:length(labels), labels)
    axis.xticklabelrotation = π / 4
    axis.xlabel = "residues"
    axis.ylabel = "residues"
    return graph
end

function plot_contact_map(structure::Structure; kwargs...)
    fig = Makie.Figure()
    axis = Makie.Axis(fig[1, 1]; title="Residue contact map")
    plot_contact_map!(axis, structure; kwargs...)
    return fig
end

function _residue_color(residue::Residue)
    region = ismissing(residue.secondary_structure) ? :unknown : residue.secondary_structure
    region == 'H' && return "#ef4444"
    region == 'E' && return "#2563eb"
    region == 'G' && return "#f97316"
    region == 'I' && return "#8b5cf6"
    region == 'B' && return "#06b6d4"
    region == 'T' && return "#10b981"
    region == 'S' && return "#84cc16"
    residue.name == "GLY" && return "#7c3aed"
    residue.name == "PRO" && return "#f59e0b"
    return "#64748b"
end

function _residue_label(chain::Chain, residue::Residue)
    return string(chain.id, ":", residue.name, residue.seqnum, residue.insertion_code)
end

function _residue_id(chain::Chain, residue::Residue)
    return (chain.id, residue.seqnum, residue.insertion_code)
end

function _protein_backbone_points(chain::Chain)
    points = Float64[]
    labels = String[]
    residues = Residue[]
    for residue in chain.residues
        ca = _residue_atom(residue, "CA")
        ca === nothing && continue
        push!(points, ca.x, ca.y, ca.z)
        push!(labels, _residue_label(chain, residue))
        push!(residues, residue)
    end
    if isempty(residues)
        return Matrix{Float64}(undef, 0, 3), String[], Residue[]
    end
    return reshape(points, 3, :)', labels, residues
end

function _residue_midpoint(chain::Chain, residue_index::Int)
    1 <= residue_index <= length(chain.residues) || return nothing
    residue = chain.residues[residue_index]
    ca = _residue_atom(residue, "CA")
    ca !== nothing && return (ca.x, ca.y, ca.z)
    atom_count = length(residue.atoms)
    atom_count == 0 && return nothing
    x = 0.0
    y = 0.0
    z = 0.0
    for atom in residue.atoms
        x += atom.x
        y += atom.y
        z += atom.z
    end
    invn = 1 / atom_count
    return (x * invn, y * invn, z * invn)
end

function _chain_backbone_points(chain::Chain)
    points = NamedTuple{(:x, :y, :z), Tuple{Float64, Float64, Float64}}[]
    labels = String[]
    residue_refs = Residue[]
    for residue in chain.residues
        midpoint = _residue_midpoint(chain, length(residue_refs) + 1)
        midpoint === nothing && continue
        push!(points, (x=midpoint[1], y=midpoint[2], z=midpoint[3]))
        push!(labels, _residue_label(chain, residue))
        push!(residue_refs, residue)
    end
    return points, labels, residue_refs
end

function _chain_backbone_coords(chain::Chain)
    coords = Float64[]
    labels = String[]
    residues = Residue[]
    for residue in chain.residues
        ca = _residue_atom(residue, "CA")
        ca === nothing && continue
        push!(coords, ca.x, ca.y, ca.z)
        push!(labels, _residue_label(chain, residue))
        push!(residues, residue)
    end
    if isempty(residues)
        return Matrix{Float64}(undef, 0, 3), String[], Residue[]
    end
    return reshape(coords, 3, :)', labels, residues
end

function _chain_backbone_pairs(coords::AbstractMatrix)
    if size(coords, 1) < 2
        return Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    end
    x1 = Float64[]
    y1 = Float64[]
    z1 = Float64[]
    x2 = Float64[]
    y2 = Float64[]
    z2 = Float64[]
    for index in 1:(size(coords, 1) - 1)
        push!(x1, coords[index, 1]); push!(y1, coords[index, 2]); push!(z1, coords[index, 3])
        push!(x2, coords[index + 1, 1]); push!(y2, coords[index + 1, 2]); push!(z2, coords[index + 1, 3])
    end
    return x1, y1, z1, x2, y2, z2
end

function _binding_axis3(parent)
    if isdefined(Main, :Makie) && parent isa Makie.Axis3
        return parent
    end
    return parent
end

function plot_structure_atoms!(axis, structure::Structure; model_index::Int=1, atom_size::Real=12)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    colors = Any[]
    for chain in structure.models[model_index].chains
        for residue in chain.residues
            color = _residue_color(residue)
            for atom in residue.atoms
                push!(xs, atom.x)
                push!(ys, atom.y)
                push!(zs, atom.z)
                push!(colors, color)
            end
        end
    end
    Makie.meshscatter!(axis, xs, ys, zs; markersize=atom_size, color=colors)
    return axis
end

function plot_structure_atoms(structure::Structure; kwargs...)
    fig = Makie.Figure()
    axis = Makie.Axis3(fig[1, 1]; title="Atom cloud")
    plot_structure_atoms!(axis, structure; kwargs...)
    return fig
end

function plot_backbone_trace!(axis, structure::Structure; model_index::Int=1, line_width::Real=4)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    model = structure.models[model_index]
    for chain in model.chains
        coords, labels, residues = _chain_backbone_coords(chain)
        size(coords, 1) < 2 && continue
        xs = coords[:, 1]
        ys = coords[:, 2]
        zs = coords[:, 3]
        colors = [_residue_color(residue) for residue in residues]
        Makie.lines!(axis, xs, ys, zs; color=first(colors), linewidth=line_width)
        Makie.scatter!(axis, xs, ys, zs; color=colors, markersize=line_width + 2)
    end
    return axis
end

function plot_backbone_trace(structure::Structure; kwargs...)
    fig = Makie.Figure()
    axis = Makie.Axis3(fig[1, 1]; title="Backbone trace")
    plot_backbone_trace!(axis, structure; kwargs...)
    return fig
end

function _ribbon_segment_points(chain::Chain)
    coords, labels, residues = _chain_backbone_coords(chain)
    segments = NamedTuple{(:start, :stop, :midpoint, :label, :color), Tuple{NTuple{3,Float64}, NTuple{3,Float64}, NTuple{3,Float64}, String, String}}[]
    if size(coords, 1) < 2
        return segments
    end
    for index in 1:(size(coords, 1) - 1)
        start_point = (coords[index, 1], coords[index, 2], coords[index, 3])
        stop_point = (coords[index + 1, 1], coords[index + 1, 2], coords[index + 1, 3])
        midpoint = ((start_point[1] + stop_point[1]) / 2, (start_point[2] + stop_point[2]) / 2, (start_point[3] + stop_point[3]) / 2)
        region = ismissing(residues[index].secondary_structure) ? :unknown : residues[index].secondary_structure
        color = _residue_color(residues[index])
        push!(segments, (start=start_point, stop=stop_point, midpoint=midpoint, label=labels[index], color=color))
    end
    return segments
end

function plot_chain_ribbon!(axis, chain::Chain; ribbon_width::Real=10, point_size::Real=16)
    coords, labels, residues = _chain_backbone_coords(chain)
    size(coords, 1) == 0 && return axis
    colors = [_residue_color(residue) for residue in residues]
    if size(coords, 1) > 1
        for index in 1:(size(coords, 1) - 1)
            Makie.lines!(axis, coords[index:index+1, 1], coords[index:index+1, 2], coords[index:index+1, 3]; color=colors[index], linewidth=ribbon_width)
        end
    end
    Makie.scatter!(axis, coords[:, 1], coords[:, 2], coords[:, 3]; color=colors, markersize=point_size)
    return axis
end

function plot_chain_ribbon(structure::Structure; model_index::Int=1, chain_id::Union{Nothing,AbstractString}=nothing, kwargs...)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    fig = Makie.Figure()
    axis = Makie.Axis3(fig[1, 1]; title="Chain ribbon")
    model = structure.models[model_index]
    for chain in model.chains
        chain_id !== nothing && chain.id != String(chain_id) && continue
        plot_chain_ribbon!(axis, chain; kwargs...)
    end
    return fig
end

function residue_pick_hooks(structure::Structure; model_index::Int=1)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    graph = structure_contact_graph(structure; model_index=model_index)
    residue_lookup = Dict{String,Tuple{String,Int,Char,Residue}}()
    for chain in structure.models[model_index].chains
        for residue in chain.residues
            node_id = nothing
            for node in graph.nodes
                if node.chain == chain.id && node.seqnum == residue.seqnum && node.insertion_code == residue.insertion_code
                    node_id = node.id
                    break
                end
            end
            node_id === nothing && continue
            residue_lookup[node_id] = (chain.id, residue.seqnum, residue.insertion_code, residue)
        end
    end

    function resolve_node(node_id::AbstractString)
        entry = get(residue_lookup, String(node_id), nothing)
        entry === nothing && return nothing
        return (chain=entry[1], seqnum=entry[2], insertion_code=entry[3], residue=entry[4])
    end

    function pick_callback(node_id::AbstractString)
        return resolve_node(node_id)
    end

    return (lookup=residue_lookup, resolve=resolve_node, on_pick=pick_callback)
end

function connect_residue_picking!(axis, structure::Structure; model_index::Int=1, on_pick=nothing, cutoff::Real=8.0, hbond_cutoff::Real=3.5)
    hooks = residue_pick_hooks(structure; model_index=model_index)
    if on_pick === nothing
        on_pick = residue -> nothing
    end
    try
        scene = axis.scene
        click_state = Ref(false)
        Makie.on(events(scene).mousebutton) do event
            if event.button == Makie.Mouse.left && event.action == Makie.Mouse.press
                click_state[] = true
            elseif event.button == Makie.Mouse.left && event.action == Makie.Mouse.release && click_state[]
                click_state[] = false
                try
                    selection = Makie.pick(scene, events(scene).mouseposition[])
                    if selection !== nothing && hasproperty(selection, :plot) && hasproperty(selection.plot, :id)
                        picked = hooks.resolve(string(selection.plot.id))
                        picked !== nothing && on_pick(picked)
                    end
                catch
                end
            end
        end
    catch
    end
    return hooks
end

function plot_structure_viewer(structure::Structure; model_index::Int=1, atom_size::Real=8, ribbon_width::Real=5, on_pick=nothing)
    fig = Makie.Figure()
    axis = Makie.Axis3(fig[1, 1]; title="Structure viewer")
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    for chain in structure.models[model_index].chains
        plot_chain_ribbon!(axis, chain; ribbon_width=ribbon_width)
    end
    plot_structure_atoms!(axis, structure; model_index=model_index, atom_size=atom_size)
    connect_residue_picking!(axis, structure; model_index=model_index, on_pick=on_pick)
    return fig
end

function _svg_escape(text::AbstractString)
    escaped = replace(String(text), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", '"' => "&quot;")
    return escaped
end

function contact_map_svg(structure::Structure; model_index::Int=1, cutoff::Real=8.0, hbond_cutoff::Real=3.5, cell_size::Int=14, margin::Int=36)
    graph = structure_contact_graph(structure; model_index=model_index, cutoff=cutoff, hbond_cutoff=hbond_cutoff)
    n = length(graph.nodes)
    width = margin * 2 + max(n, 1) * cell_size
    height = margin * 2 + max(n, 1) * cell_size
    labels = [node.label for node in graph.nodes]
    node_ids = [node.id for node in graph.nodes]
    generic_pairs = Set{Tuple{String,String}}()
    for edge in graph.generic_edges
        source = edge.source
        target = edge.target
        source < target ? push!(generic_pairs, (source, target)) : push!(generic_pairs, (target, source))
    end
    hbond_pairs = Set{Tuple{String,String}}()
    for edge in graph.hydrogen_bond_edges
        source = edge.source
        target = edge.target
        source < target ? push!(hbond_pairs, (source, target)) : push!(hbond_pairs, (target, source))
    end

    io = IOBuffer()
    println(io, @sprintf("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 %d %d\" width=\"%d\" height=\"%d\">", width, height, width, height))
    println(io, "<rect width=\"100%%\" height=\"100%%\" fill=\"#fffdf8\"/>")
    println(io, @sprintf("<text x=\"%d\" y=\"22\" font-family=\"monospace\" font-size=\"12\" fill=\"#2a2a2a\">Residue contact map</text>", margin))

    for index in 1:n
        x = margin + (index - 1) * cell_size
        y = margin + (index - 1) * cell_size
        println(io, @sprintf("<text x=\"%d\" y=\"%d\" font-family=\"monospace\" font-size=\"8\" fill=\"#4a4a4a\" transform=\"rotate(-45 %d %d)\">%s</text>", x + 2, margin - 6, x + 2, margin - 6, _svg_escape(labels[index])))
        println(io, @sprintf("<text x=\"%d\" y=\"%d\" font-family=\"monospace\" font-size=\"8\" fill=\"#4a4a4a\">%s</text>", 6, y + 9, _svg_escape(labels[index])))
    end

    for i in 1:n
        for j in 1:n
            x = margin + (j - 1) * cell_size
            y = margin + (i - 1) * cell_size
            pair = i < j ? (node_ids[i], node_ids[j]) : (node_ids[j], node_ids[i])
            if pair in generic_pairs
                fill = pair in hbond_pairs ? "#ef4444" : "#9ca3af"
                println(io, @sprintf("<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"%s\" opacity=\"0.85\"/>", x, y, cell_size - 1, cell_size - 1, fill))
            end
        end
    end

    println(io, "</svg>")
    return String(take!(io))
end

function write_contact_map_svg(filepath::AbstractString, structure::Structure; kwargs...)
    open(filepath, "w") do io
        write(io, contact_map_svg(structure; kwargs...))
    end
    return filepath
end

function structure_contact_mermaid(structure::Structure; model_index::Int=1, cutoff::Real=8.0, hbond_cutoff::Real=3.5)
    1 <= model_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    model = structure.models[model_index]
    nodes = Vector{NamedTuple{(:id, :label), Tuple{String, String}}}()
    lookup = Dict{Tuple{String,Int,Char},String}()
    counter = 1
    for chain in model.chains
        for residue in chain.residues
            node_id = string("r", counter)
            push!(nodes, (id=node_id, label=_residue_node_label(chain.id, residue)))
            lookup[(chain.id, residue.seqnum, residue.insertion_code)] = node_id
            counter += 1
        end
    end

    lines = String["flowchart LR"]
    for node in nodes
        push!(lines, string("    ", node.id, "[\"", replace(node.label, "\"" => "'"), "\"]"))
    end

    for contact in residue_contacts(structure; cutoff=cutoff, model_index=model_index)
        left = lookup[(contact.left_chain, contact.left_seqnum, contact.left_insertion_code)]
        right = lookup[(contact.right_chain, contact.right_seqnum, contact.right_insertion_code)]
        push!(lines, string("    ", left, " --- ", right))
    end

    for bond in hydrogen_bonds(structure; model_index=model_index, cutoff=hbond_cutoff)
        left = lookup[(bond.donor_chain, bond.donor_seqnum, bond.donor_insertion_code)]
        right = lookup[(bond.acceptor_chain, bond.acceptor_seqnum, bond.acceptor_insertion_code)]
        push!(lines, string("    ", left, " -.-> ", right))
    end

    return join(lines, "\n")
end

function write_structure_mermaid(filepath::AbstractString, structure::Structure; kwargs...)
    open(filepath, "w") do io
        write(io, structure_contact_mermaid(structure; kwargs...))
    end
    return filepath
end

function _pdb_field(line::AbstractString, start_index::Int, stop_index::Int)
    start_index > lastindex(line) && return ""
    stop_index = min(stop_index, lastindex(line))
    return String(strip(String(SubString(line, start_index, stop_index))))
end

function _pdb_char(line::AbstractString, index::Int, default::Char=' ')
    index > lastindex(line) && return default
    return line[index]
end

function _ensure_model!(structure::Structure, model_map::Dict{Int,Model}, model_id::Int)
    model = get(model_map, model_id, nothing)
    model !== nothing && return model
    model = Model(model_id)
    push!(structure.models, model)
    model_map[model_id] = model
    return model
end

function _ensure_chain!(model::Model, chain_map::Dict{String,Chain}, chain_id::AbstractString)
    key = String(chain_id)
    chain = get(chain_map, key, nothing)
    chain !== nothing && return chain
    chain = Chain(key)
    push!(model.chains, chain)
    chain_map[key] = chain
    return chain
end

function _ensure_residue!(chain::Chain, residue_map::Dict{Tuple{String,Int,Char,String},Residue}, residue_name::AbstractString, seqnum::Int, insertion_code::Char)
    key = (chain.id, seqnum, insertion_code, String(residue_name))
    residue = get(residue_map, key, nothing)
    residue !== nothing && return residue
    residue = Residue(String(residue_name), seqnum, insertion_code)
    push!(chain.residues, residue)
    residue_map[key] = residue
    return residue
end

function _atom_element(atom_name::AbstractString, raw_element::AbstractString)
    isempty(raw_element) || return uppercase(raw_element)
    cleaned = filter(isletter, atom_name)
    isempty(cleaned) && return ""
    return uppercase(String(cleaned[1:min(length(cleaned), 2)]))
end

function _pdb_atom_from_line(line::AbstractString)
    record_name = _pdb_field(line, 1, 6)
    record_name in ("ATOM", "HETATM") || return nothing

    serial = tryparse(Int, _pdb_field(line, 7, 11))
    serial === nothing && return nothing

    atom_name = _pdb_field(line, 13, 16)
    altloc = _pdb_char(line, 17)
    residue_name = _pdb_field(line, 18, 20)
    chain_id = _pdb_field(line, 22, 22)
    chain_id = isempty(chain_id) ? " " : chain_id
    seqnum = tryparse(Int, _pdb_field(line, 23, 26))
    seqnum === nothing && return nothing
    insertion_code = _pdb_char(line, 27)
    x = tryparse(Float64, _pdb_field(line, 31, 38))
    y = tryparse(Float64, _pdb_field(line, 39, 46))
    z = tryparse(Float64, _pdb_field(line, 47, 54))
    x === nothing && return nothing
    y === nothing && return nothing
    z === nothing && return nothing
    occupancy = something(tryparse(Float64, _pdb_field(line, 55, 60)), 1.0)
    bfactor = something(tryparse(Float64, _pdb_field(line, 61, 66)), 0.0)
    raw_element = _pdb_field(line, 77, 78)
    charge = _pdb_field(line, 79, 80)
    hetatm = record_name == "HETATM"

    atom = Atom(serial, atom_name, altloc, _atom_element(atom_name, raw_element), x, y, z, occupancy, bfactor, charge, hetatm)
    return (atom=atom, residue_name=residue_name, chain_id=chain_id, seqnum=seqnum, insertion_code=insertion_code)
end

const _PDB_METADATA_RECORDS = Set(["HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "REMARK", "JRNL"])

function _capture_pdb_metadata!(metadata::Dict{String,String}, line::AbstractString)
    record = uppercase(_pdb_field(line, 1, 6))
    if record in _PDB_METADATA_RECORDS
        payload = strip(_pdb_field(line, 11, lastindex(line)))
        _metadata_append!(metadata, string("pdb:", record), payload)
        return true
    end
    return false
end

function read_pdb(filepath::AbstractString)
    open(filepath, "r") do io
        return read_pdb(io; structure_id=splitext(basename(filepath))[1])
    end
end

function read_pdb(io::IO; structure_id::AbstractString="structure")
    structure = Structure(structure_id)
    model_map = Dict{Int,Model}()
    chain_maps = Dict{Int,Dict{String,Chain}}()
    residue_maps = Dict{Int,Dict{Tuple{String,Int,Char,String},Residue}}()
    current_model_id = 1
    current_model = _ensure_model!(structure, model_map, current_model_id)
    chain_maps[current_model_id] = Dict{String,Chain}()
    residue_maps[current_model_id] = Dict{Tuple{String,Int,Char,String},Residue}()

    for raw_line in eachline(io)
        line = rpad(raw_line, 80)
        _capture_pdb_metadata!(structure.metadata, line) && continue
        startswith(line, "MODEL") && begin
            parsed = tryparse(Int, strip(_pdb_field(line, 11, 14)))
            current_model_id = parsed === nothing ? length(structure.models) + 1 : parsed
            current_model = _ensure_model!(structure, model_map, current_model_id)
            chain_maps[current_model_id] = get(chain_maps, current_model_id, Dict{String,Chain}())
            residue_maps[current_model_id] = get(residue_maps, current_model_id, Dict{Tuple{String,Int,Char,String},Residue}())
            continue
        end
        startswith(line, "ENDMDL") && begin
            current_model_id = length(structure.models) + 1
            continue
        end
        startswith(line, "ATOM") || startswith(line, "HETATM") || continue
        parsed_atom = _pdb_atom_from_line(line)
        parsed_atom === nothing && continue
        atom = parsed_atom.atom
        chain_id = parsed_atom.chain_id
        residue_name = parsed_atom.residue_name
        seqnum = parsed_atom.seqnum
        insertion_code = parsed_atom.insertion_code

        chain_map = get!(chain_maps, current_model_id, Dict{String,Chain}())
        residue_map = get!(residue_maps, current_model_id, Dict{Tuple{String,Int,Char,String},Residue}())
        chain = _ensure_chain!(current_model, chain_map, chain_id)
        residue = _ensure_residue!(chain, residue_map, residue_name, seqnum, insertion_code)
        push!(residue.atoms, atom)
    end

    isempty(structure.models) && push!(structure.models, Model(1))
    return structure
end

function _write_pdb_metadata!(io::IO, structure::Structure)
    metadata = structure.metadata
    if haskey(metadata, "pdb:HEADER")
        println(io, @sprintf("HEADER    %s", metadata["pdb:HEADER"]))
    end
    for record in ("TITLE", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "JRNL", "REMARK")
        key = string("pdb:", record)
        haskey(metadata, key) || continue
        for line in Base.split(metadata[key], "\n")
            isempty(strip(line)) && continue
            println(io, @sprintf("%-6s    %s", record, line))
        end
    end
    return nothing
end

function _mmcif_atom_value(value)
    text = String(value)
    isempty(text) && return "."
    return replace(text, ' ' => '_')
end

function _write_mmcif_metadata!(io::IO, structure::Structure)
    metadata = structure.metadata
    category_keys = sort([key for key in keys(metadata) if startswith(key, "mmcif:category:") || startswith(key, "mmcif:item:")])
    if !isempty(category_keys)
        for key in category_keys
            block = metadata[key]
            isempty(strip(block)) && continue
            println(io, block)
            endswith(strip(block), "#") || println(io, "#")
        end
        return nothing
    end
    haskey(metadata, "mmcif:raw_blocks") || return nothing
    for block in Base.split(metadata["mmcif:raw_blocks"], "\n\n")
        isempty(strip(block)) && continue
        println(io, block)
        endswith(strip(block), "#") || println(io, "#")
    end
    return nothing
end

function _mmcif_atom_fields(field_names::Vector{String})
    lookup = Dict{String,Int}()
    for (index, name) in pairs(field_names)
        lookup[name] = index
    end
    return lookup
end

function _mmcif_first_index(lookup::Dict{String,Int}, keys::Vector{String})
    for key in keys
        index = get(lookup, key, 0)
        index > 0 && return index
    end
    return 0
end

function _mmcif_value(row::AbstractVector{<:AbstractString}, lookup::Dict{String,Int}, keys::Vector{String}, default::String="")
    for key in keys
        index = get(lookup, key, 0)
        index > 0 && index <= length(row) && return row[index]
    end
    return default
end

function _mmcif_value(row::AbstractVector{<:AbstractString}, index::Int, default::String="")
    index > 0 && index <= length(row) && return row[index]
    return default
end

function read_mmcif(filepath::AbstractString)
    open(filepath, "r") do io
        return read_mmcif(io; structure_id=splitext(basename(filepath))[1])
    end
end

function read_mmcif(io::IO; structure_id::AbstractString="structure")
    lines = collect(eachline(io))
    structure = Structure(structure_id)
    model_map = Dict{Int,Model}()
    chain_maps = Dict{Int,Dict{String,Chain}}()
    residue_maps = Dict{Int,Dict{Tuple{String,Int,Char,String},Residue}}()
    raw_blocks = String[]
    i = 1

    while i <= length(lines)
        line = strip(lines[i])
        if startswith(line, "data_")
            structure.id = isempty(structure.id) ? line[6:end] : structure.id
            i += 1
            continue
        end

        if isempty(line) || line == "#"
            i += 1
            continue
        end

        if line == "loop_"
            block_start = i
            i += 1
            field_names = String[]
            while i <= length(lines) && startswith(strip(lines[i]), "_")
                push!(field_names, strip(lines[i]))
                i += 1
            end
            if !isempty(field_names) && any(startswith(field, "_atom_site.") for field in field_names)
                atom_fields = _mmcif_atom_fields(field_names)
                model_id_index = _mmcif_first_index(atom_fields, ["_atom_site.pdbx_PDB_model_num", "_atom_site.pdbx_model_num", "_atom_site.model_id"])
                atom_name_index = _mmcif_first_index(atom_fields, ["_atom_site.auth_atom_id", "_atom_site.label_atom_id"])
                residue_name_index = _mmcif_first_index(atom_fields, ["_atom_site.auth_comp_id", "_atom_site.label_comp_id"])
                chain_id_index = _mmcif_first_index(atom_fields, ["_atom_site.auth_asym_id", "_atom_site.label_asym_id"])
                seqnum_index = _mmcif_first_index(atom_fields, ["_atom_site.auth_seq_id", "_atom_site.label_seq_id"])
                insertion_code_index = _mmcif_first_index(atom_fields, ["_atom_site.pdbx_PDB_ins_code"])
                altloc_index = _mmcif_first_index(atom_fields, ["_atom_site.pdbx_PDB_alt_id", "_atom_site.label_alt_id"])
                x_index = _mmcif_first_index(atom_fields, ["_atom_site.Cartn_x"])
                y_index = _mmcif_first_index(atom_fields, ["_atom_site.Cartn_y"])
                z_index = _mmcif_first_index(atom_fields, ["_atom_site.Cartn_z"])
                occupancy_index = _mmcif_first_index(atom_fields, ["_atom_site.occupancy"])
                bfactor_index = _mmcif_first_index(atom_fields, ["_atom_site.B_iso_or_equiv"])
                element_index = _mmcif_first_index(atom_fields, ["_atom_site.type_symbol"])
                group_index = _mmcif_first_index(atom_fields, ["_atom_site.group_PDB"])
                charge_index = _mmcif_first_index(atom_fields, ["_atom_site.pdbx_formal_charge"])
                serial_index = _mmcif_first_index(atom_fields, ["_atom_site.id"])
                current_model_id = 1
                current_model = _ensure_model!(structure, model_map, current_model_id)
                chain_maps[current_model_id] = get(chain_maps, current_model_id, Dict{String,Chain}())
                residue_maps[current_model_id] = get(residue_maps, current_model_id, Dict{Tuple{String,Int,Char,String},Residue}())

                while i <= length(lines)
                    row_line = strip(lines[i])
                    isempty(row_line) && break
                    startswith(row_line, "loop_") && break
                    startswith(row_line, "_") && break
                    startswith(row_line, "#") && break

                    row = _mmcif_split_row(row_line)
                    if length(row) < length(field_names)
                        i += 1
                        continue
                    end

                    model_id = tryparse(Int, _mmcif_value(row, model_id_index, "1"))
                    current_model_id = model_id === nothing ? 1 : model_id
                    current_model = _ensure_model!(structure, model_map, current_model_id)
                    chain_map = get!(chain_maps, current_model_id, Dict{String,Chain}())
                    residue_map = get!(residue_maps, current_model_id, Dict{Tuple{String,Int,Char,String},Residue}())

                    atom_name = String(_mmcif_value(row, atom_name_index, ""))
                    residue_name = String(_mmcif_value(row, residue_name_index, ""))
                    chain_id = String(_mmcif_value(row, chain_id_index, " "))
                    seqnum_value = String(_mmcif_value(row, seqnum_index, "0"))
                    seqnum = something(tryparse(Int, seqnum_value), 0)
                    insertion_code_value = String(_mmcif_value(row, insertion_code_index, " "))
                    insertion_code = isempty(insertion_code_value) ? ' ' : first(insertion_code_value)
                    altloc_value = String(_mmcif_value(row, altloc_index, "."))
                    altloc = isempty(altloc_value) || altloc_value in (".", "?") ? ' ' : first(altloc_value)
                    x = something(tryparse(Float64, _mmcif_value(row, x_index, "0.0")), 0.0)
                    y = something(tryparse(Float64, _mmcif_value(row, y_index, "0.0")), 0.0)
                    z = something(tryparse(Float64, _mmcif_value(row, z_index, "0.0")), 0.0)
                    occupancy = something(tryparse(Float64, _mmcif_value(row, occupancy_index, "1.0")), 1.0)
                    bfactor = something(tryparse(Float64, _mmcif_value(row, bfactor_index, "0.0")), 0.0)
                    element = String(_mmcif_value(row, element_index, ""))
                    hetatm = uppercase(String(_mmcif_value(row, group_index, "ATOM"))) == "HETATM"
                    charge = String(_mmcif_value(row, charge_index, ""))
                    serial = something(tryparse(Int, _mmcif_value(row, serial_index, "0")), 0)
                    atom = Atom(serial, atom_name, altloc, element, x, y, z, occupancy, bfactor, charge, hetatm)
                    chain = _ensure_chain!(current_model, chain_map, chain_id)
                    residue = _ensure_residue!(chain, residue_map, residue_name, seqnum, insertion_code)
                    push!(residue.atoms, atom)
                    i += 1
                end
                continue
            else
                block = join(lines[block_start:i-1], "\n")
                push!(raw_blocks, block)
                if !isempty(field_names)
                    category_name = first(Base.split(first(field_names), "."; limit=2))
                    if !startswith(category_name, "_atom_site.") && category_name != "_atom_site"
                        _metadata_append!(structure.metadata, string("mmcif:category:", category_name), block)
                    end
                end
                continue
            end
        end

        if startswith(line, "_")
            push!(raw_blocks, lines[i])
            _metadata_append!(structure.metadata, string("mmcif:item:", first(Base.split(lines[i], " "; limit=2))), lines[i])
            i += 1
            continue
        end

        i += 1
    end

    isempty(raw_blocks) || _metadata_block!(structure.metadata, "mmcif:raw_blocks", join(raw_blocks, "\n"))

    isempty(structure.models) && push!(structure.models, Model(1))
    return structure
end

function write_pdb(filepath::AbstractString, structure::Structure)
    open(filepath, "w") do io
        write_pdb(io, structure)
    end
    return filepath
end

function write_pdb(io::IO, structure::Structure)
    _write_pdb_metadata!(io, structure)
    for model in structure.models
        if length(structure.models) > 1
            println(io, @sprintf("MODEL     %4d", model.id))
        end
        for chain in model.chains
            for residue in chain.residues
                for atom in residue.atoms
                    record_name = atom.hetatm ? "HETATM" : "ATOM"
                    atom_name = isempty(atom.name) ? "    " : rpad(atom.name[1:min(length(atom.name), 4)], 4)
                    element = rpad(atom.element, 2)
                    charge = rpad(atom.charge, 2)
                    println(io, @sprintf("%-6s%5d %-4s%c%3s %1s%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %-2s%2s",
                        record_name,
                        atom.serial,
                        atom_name,
                        atom.altloc,
                        residue.name,
                        isempty(chain.id) ? " " : chain.id[1],
                        residue.seqnum,
                        residue.insertion_code,
                        atom.x,
                        atom.y,
                        atom.z,
                        atom.occupancy,
                        atom.bfactor,
                        element,
                        charge))
                end
            end
            println(io, "TER")
        end
        if length(structure.models) > 1
            println(io, "ENDMDL")
        end
    end
    println(io, "END")
    return nothing
end

function write_mmcif(filepath::AbstractString, structure::Structure)
    open(filepath, "w") do io
        write_mmcif(io, structure)
    end
    return filepath
end

function write_mmcif(io::IO, structure::Structure)
    println(io, string("data_", isempty(structure.id) ? "structure" : structure.id))
    _write_mmcif_metadata!(io, structure)
    println(io, "loop_")
    atom_fields = [
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.auth_atom_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id",
        "_atom_site.auth_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.pdbx_PDB_model_num",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.pdbx_formal_charge",
    ]
    for field in atom_fields
        println(io, field)
    end
    for model in structure.models
        for chain in model.chains
            for residue in chain.residues
                for atom in residue.atoms
                    group_pdb = atom.hetatm ? "HETATM" : "ATOM"
                    altloc = atom.altloc == ' ' ? "." : string(atom.altloc)
                    insertion_code = residue.insertion_code == ' ' ? "?" : string(residue.insertion_code)
                    println(io, join([
                        group_pdb,
                        string(atom.serial),
                        _mmcif_atom_value(atom.element),
                        _mmcif_atom_value(atom.name),
                        _mmcif_atom_value(altloc),
                        _mmcif_atom_value(residue.name),
                        _mmcif_atom_value(chain.id),
                        string(residue.seqnum),
                        _mmcif_atom_value(atom.name),
                        _mmcif_atom_value(residue.name),
                        _mmcif_atom_value(chain.id),
                        string(residue.seqnum),
                        insertion_code,
                        string(model.id),
                        @sprintf("%.3f", atom.x),
                        @sprintf("%.3f", atom.y),
                        @sprintf("%.3f", atom.z),
                        @sprintf("%.2f", atom.occupancy),
                        @sprintf("%.2f", atom.bfactor),
                        isempty(atom.charge) ? "?" : atom.charge,
                    ], " "))
                end
            end
        end
    end
    println(io, "#")
    return nothing
end

atom_coordinates(atom::Atom) = (atom.x, atom.y, atom.z)

function coordinate_matrix(atoms::AbstractVector{Atom})
    points = Matrix{Float64}(undef, length(atoms), 3)
    for (index, atom) in enumerate(atoms)
        points[index, 1] = atom.x
        points[index, 2] = atom.y
        points[index, 3] = atom.z
    end
    return points
end

coordinate_matrix(structure::Structure) = coordinate_matrix(structure_atoms(structure))

function _centroid(points::AbstractMatrix)
    count = size(points, 1)
    count > 0 || throw(ArgumentError("at least one point is required"))
    centroid = zeros(Float64, 3)
    for row in 1:count
        centroid[1] += points[row, 1]
        centroid[2] += points[row, 2]
        centroid[3] += points[row, 3]
    end
    centroid ./= count
    return centroid
end

function kabsch(reference::AbstractMatrix, mobile::AbstractMatrix)
    size(reference, 2) == 3 && size(mobile, 2) == 3 || throw(ArgumentError("coordinates must be N×3 matrices"))
    size(reference, 1) == size(mobile, 1) || throw(ArgumentError("point sets must have the same length"))
    npoints = size(reference, 1)
    npoints > 0 || throw(ArgumentError("point sets must not be empty"))

    reference_centroid = _centroid(reference)
    mobile_centroid = _centroid(mobile)

    reference_centered = reference .- reference_centroid'
    mobile_centered = mobile .- mobile_centroid'

    covariance = mobile_centered' * reference_centered
    decomposition = svd(covariance)
    correction = Diagonal([1.0, 1.0, sign(det(decomposition.U * decomposition.Vt))])
    rotation = decomposition.U * correction * decomposition.Vt
    translation = reference_centroid .- vec(mobile_centroid' * rotation)
    aligned = mobile * rotation .+ translation'
    rmsd_value = sqrt(sum((reference .- aligned) .^ 2) / npoints)
    return SuperpositionResult(rotation, translation, rmsd_value), aligned
end

function _coordinates(entity)
    entity isa Structure && return coordinate_matrix(entity)
    entity isa AbstractVector{Atom} && return coordinate_matrix(entity)
    entity isa AbstractMatrix && return Matrix{Float64}(entity)
    throw(ArgumentError("unsupported coordinate container"))
end

function rmsd(reference, mobile; superpose::Bool=true)
    reference_points = _coordinates(reference)
    mobile_points = _coordinates(mobile)
    if superpose
        result, _ = kabsch(reference_points, mobile_points)
        return result.rmsd
    end
    size(reference_points) == size(mobile_points) || throw(ArgumentError("point sets must have the same shape"))
    count = size(reference_points, 1)
    count > 0 || return 0.0
    total = 0.0
    for row in 1:count
        dx = reference_points[row, 1] - mobile_points[row, 1]
        dy = reference_points[row, 2] - mobile_points[row, 2]
        dz = reference_points[row, 3] - mobile_points[row, 3]
        total += dx * dx + dy * dy + dz * dz
    end
    return sqrt(total / count)
end

function _apply_transform!(atoms::AbstractVector{Atom}, rotation::Matrix{Float64}, translation::Vector{Float64})
    for atom in atoms
        x = atom.x * rotation[1, 1] + atom.y * rotation[2, 1] + atom.z * rotation[3, 1] + translation[1]
        y = atom.x * rotation[1, 2] + atom.y * rotation[2, 2] + atom.z * rotation[3, 2] + translation[2]
        z = atom.x * rotation[1, 3] + atom.y * rotation[2, 3] + atom.z * rotation[3, 3] + translation[3]
        atom.x = x
        atom.y = y
        atom.z = z
    end
    return atoms
end

function superpose(reference, mobile)
    result, _ = kabsch(_coordinates(reference), _coordinates(mobile))
    return result
end

function superpose!(reference, mobile)
    reference_points = _coordinates(reference)
    mobile_points = _coordinates(mobile)
    result, _ = kabsch(reference_points, mobile_points)
    mobile isa Structure && _apply_transform!(structure_atoms(mobile), result.rotation, result.translation)
    mobile isa AbstractVector{Atom} && _apply_transform!(mobile, result.rotation, result.translation)
    return result
end

function _model_atoms(model::Model)
    atoms = Atom[]
    for chain in model.chains
        for residue in chain.residues
            append!(atoms, residue.atoms)
        end
    end
    return atoms
end

function _selected_model_atoms(model::Model; atom_selector::Union{Symbol,AbstractString}=Symbol("CA"))
    atoms = Atom[]
    normalized = atom_selector isa Symbol ? atom_selector : Symbol(uppercase(String(atom_selector)))
    for chain in model.chains
        for residue in chain.residues
            if normalized === :all
                append!(atoms, residue.atoms)
            elseif normalized === :backbone
                for atom_name in ("N", "CA", "C", "O")
                    atom = _atom_lookup(residue, atom_name)
                    atom === nothing || push!(atoms, atom)
                end
            else
                atom = _atom_lookup(residue, String(normalized))
                atom === nothing || push!(atoms, atom)
            end
        end
    end
    return atoms
end

function superpose_models!(structure::Structure; reference_index::Int=1, atom_selector::Union{Symbol,AbstractString}=Symbol("CA"))
    1 <= reference_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    reference_model = structure.models[reference_index]
    reference_points = coordinate_matrix(_selected_model_atoms(reference_model; atom_selector=atom_selector))
    results = NamedTuple{(:model_id, :rmsd), Tuple{Int, Float64}}[]
    for (index, model) in enumerate(structure.models)
        index == reference_index && continue
        mobile_points = coordinate_matrix(_selected_model_atoms(model; atom_selector=atom_selector))
        common = min(size(reference_points, 1), size(mobile_points, 1))
        common == 0 && continue
        result, _ = kabsch(reference_points[1:common, :], mobile_points[1:common, :])
        _apply_transform!(_model_atoms(model), result.rotation, result.translation)
        push!(results, (model_id=model.id, rmsd=result.rmsd))
    end
    return results
end

function ensemble_rmsd_matrix(structure::Structure; atom_selector::Union{Symbol,AbstractString}=Symbol("CA"), superpose::Bool=true)
    count = length(structure.models)
    labels = [string("model ", model.id) for model in structure.models]
    matrix = zeros(Float64, count, count)
    selected = [coordinate_matrix(_selected_model_atoms(model; atom_selector=atom_selector)) for model in structure.models]
    for i in 1:count
        for j in i+1:count
            common = min(size(selected[i], 1), size(selected[j], 1))
            if common == 0
                matrix[i, j] = NaN
                matrix[j, i] = NaN
                continue
            end
            matrix[i, j] = superpose ? kabsch(selected[i][1:common, :], selected[j][1:common, :])[1].rmsd : rmsd(selected[i][1:common, :], selected[j][1:common, :]; superpose=false)
            matrix[j, i] = matrix[i, j]
        end
    end
    return (matrix=matrix, labels=labels, atom_selector=atom_selector)
end

function trajectory_statistics(structure::Structure; reference_index::Int=1, atom_selector::Union{Symbol,AbstractString}=Symbol("CA"))
    1 <= reference_index <= length(structure.models) || throw(ArgumentError("model index out of bounds"))
    reference = structure.models[reference_index]
    reference_points = coordinate_matrix(_selected_model_atoms(reference; atom_selector=atom_selector))
    per_model = NamedTuple{(:model_id, :rmsd), Tuple{Int, Float64}}[]
    for (index, model) in enumerate(structure.models)
        index == reference_index && push!(per_model, (model_id=model.id, rmsd=0.0))
        index == reference_index && continue
        selected = coordinate_matrix(_selected_model_atoms(model; atom_selector=atom_selector))
        common = min(size(reference_points, 1), size(selected, 1))
        common == 0 && push!(per_model, (model_id=model.id, rmsd=NaN))
        common == 0 && continue
        push!(per_model, (model_id=model.id, rmsd=rmsd(reference_points[1:common, :], selected[1:common, :]; superpose=true)))
    end
    rmsd_values = [entry.rmsd for entry in per_model if isfinite(entry.rmsd)]
    return (
        reference_model=reference.id,
        atom_selector=atom_selector,
        model_count=length(structure.models),
        per_model=per_model,
        mean_rmsd=isempty(rmsd_values) ? NaN : mean(rmsd_values),
        max_rmsd=isempty(rmsd_values) ? NaN : maximum(rmsd_values),
    )
end

function _atom_lookup(residue::Residue, name::AbstractString)
    for atom in residue.atoms
        strip(atom.name) == strip(name) && return atom
    end
    return nothing
end

function torsion_angle(p1, p2, p3, p4)
    v1 = collect(p2) .- collect(p1)
    v2 = collect(p3) .- collect(p2)
    v3 = collect(p4) .- collect(p3)
    n1 = cross(v1, v2)
    n2 = cross(v2, v3)
    norm_v2 = norm(v2)
    norm_v2 == 0 && return NaN
    m1 = cross(n1, v2 ./ norm_v2)
    return atan(dot(m1, n2), dot(n1, n2)) * 180 / π
end

function phi_psi(chain::Chain, residue_index::Int)
    1 <= residue_index <= length(chain.residues) || throw(ArgumentError("residue index out of bounds"))
    residue = chain.residues[residue_index]
    residue_n = _atom_lookup(residue, "N")
    residue_ca = _atom_lookup(residue, "CA")
    residue_c = _atom_lookup(residue, "C")
    if residue_n === nothing || residue_ca === nothing || residue_c === nothing
        return (phi=missing, psi=missing)
    end

    phi = missing
    if residue_index > 1
        previous_c = _atom_lookup(chain.residues[residue_index - 1], "C")
        previous_c !== nothing && (phi = torsion_angle(atom_coordinates(previous_c), atom_coordinates(residue_n), atom_coordinates(residue_ca), atom_coordinates(residue_c)))
    end

    psi = missing
    if residue_index < length(chain.residues)
        next_n = _atom_lookup(chain.residues[residue_index + 1], "N")
        next_n !== nothing && (psi = torsion_angle(atom_coordinates(residue_n), atom_coordinates(residue_ca), atom_coordinates(residue_c), atom_coordinates(next_n)))
    end

    return (phi=phi, psi=psi)
end

function backbone_torsions(chain::Chain)
    values = Vector{NamedTuple{(:index, :phi, :psi), Tuple{Int, Any, Any}}}(undef, length(chain.residues))
    for index in eachindex(chain.residues)
        torsions = phi_psi(chain, index)
        values[index] = (index=index, phi=torsions.phi, psi=torsions.psi)
    end
    return values
end

function _kd_build(points::Vector{Atom}, depth::Int)
    isempty(points) && return nothing
    axis = mod1(depth, 3)
    sorted_points = copy(points)
    sort!(sorted_points; by = atom -> atom_coordinates(atom)[axis])
    median_index = cld(length(sorted_points), 2)
    payload = sorted_points[median_index]
    left = _kd_build(sorted_points[1:median_index-1], depth + 1)
    right = _kd_build(sorted_points[median_index+1:end], depth + 1)
    return KDTreeNode(atom_coordinates(payload), payload, axis, left, right)
end

function build_atom_kdtree(atoms::AbstractVector{Atom})
    return AtomKDTree(_kd_build(Vector{Atom}(atoms), 1))
end

build_atom_kdtree(structure::Structure) = build_atom_kdtree(structure_atoms(structure))

function _radius_search!(matches::Vector{Atom}, node::Union{Nothing,KDTreeNode{Atom}}, point::NTuple{3,Float64}, radius2::Float64)
    node === nothing && return matches
    dx = point[1] - node.point[1]
    dy = point[2] - node.point[2]
    dz = point[3] - node.point[3]
    if dx * dx + dy * dy + dz * dz <= radius2
        push!(matches, node.payload)
    end
    axis_delta = point[node.axis] - node.point[node.axis]
    if axis_delta <= sqrt(radius2)
        _radius_search!(matches, node.left, point, radius2)
    end
    if axis_delta >= -sqrt(radius2)
        _radius_search!(matches, node.right, point, radius2)
    end
    return matches
end

function atoms_within_radius(tree::AtomKDTree, point::NTuple{3,Real}; radius::Real=5.0)
    matches = Atom[]
    _radius_search!(matches, tree.root, (Float64(point[1]), Float64(point[2]), Float64(point[3])), Float64(radius)^2)
    return matches
end

atoms_within_radius(structure::Structure, point::NTuple{3,Real}; radius::Real=5.0) = atoms_within_radius(build_atom_kdtree(structure), point; radius=radius)
atoms_within_radius(tree::AtomKDTree, atom::Atom; radius::Real=5.0) = atoms_within_radius(tree, atom_coordinates(atom); radius=radius)
atoms_within_radius(structure::Structure, atom::Atom; radius::Real=5.0) = atoms_within_radius(structure, atom_coordinates(atom); radius=radius)

function _element_color(element::String)
    upper = uppercase(element)
    upper == "C" && return "#4c78a8"
    upper == "N" && return "#1f77b4"
    upper == "O" && return "#d62728"
    upper == "S" && return "#ff7f0e"
    upper == "P" && return "#9467bd"
    return "#7f7f7f"
end

function structure_geometry(structure::Structure; bond_cutoff::Real=1.9)
    atoms = structure_atoms(structure)
    positions = coordinate_matrix(atoms)
    colors = [_element_color(atom.element) for atom in atoms]
    chains = String[]
    residues = String[]
    for model in structure.models
        for chain in model.chains
            for residue in chain.residues
                for _ in residue.atoms
                    push!(chains, chain.id)
                    push!(residues, string(residue.name, residue.seqnum, residue.insertion_code))
                end
            end
        end
    end

    atom_offsets = Dict{Atom,Int}()
    for (atom_index, atom) in enumerate(atoms)
        atom_offsets[atom] = atom_index
    end

    bonds = Tuple{Int,Int}[]
    index = 1
    for model in structure.models
        for chain in model.chains
            for residue_index in eachindex(chain.residues)
                residue = chain.residues[residue_index]
                atom_start = index
                atom_stop = index + length(residue.atoms) - 1
                for current in atom_start:(atom_stop - 1)
                    next_index = current + 1
                    dx = positions[current, 1] - positions[next_index, 1]
                    dy = positions[current, 2] - positions[next_index, 2]
                    dz = positions[current, 3] - positions[next_index, 3]
                    if dx * dx + dy * dy + dz * dz <= bond_cutoff^2
                        push!(bonds, (current, next_index))
                    end
                end
                if residue_index < length(chain.residues)
                    c_atom = _atom_lookup(residue, "C")
                    n_atom = _atom_lookup(chain.residues[residue_index + 1], "N")
                    if c_atom !== nothing && n_atom !== nothing && haskey(atom_offsets, c_atom) && haskey(atom_offsets, n_atom)
                        c_index = atom_offsets[c_atom]
                        next_index = atom_offsets[n_atom]
                        dx = positions[c_index, 1] - positions[next_index, 1]
                        dy = positions[c_index, 2] - positions[next_index, 2]
                        dz = positions[c_index, 3] - positions[next_index, 3]
                        if dx * dx + dy * dy + dz * dz <= bond_cutoff^2
                            push!(bonds, (c_index, next_index))
                        end
                    end
                end
                index = atom_stop + 1
            end
        end
    end

    return (positions=positions, colors=colors, bonds=bonds, chains=chains, residues=residues)
end

structure_pointcloud(structure::Structure) = (points=coordinate_matrix(structure), colors=[_element_color(atom.element) for atom in structure_atoms(structure)])

function _structure_to_temp_pdb(input::Structure)
    temp_path, io = mktemp(suffix=".pdb")
    try
        write_pdb(io, input)
    finally
        close(io)
    end
    return temp_path
end

function _structure_input_path(input::AbstractString)
    return input
end

function _structure_input_path(input::Structure)
    return _structure_to_temp_pdb(input)
end

function run_dssp(input::Union{AbstractString,Structure}; command::AbstractString="mkdssp", args::AbstractVector{<:AbstractString}=String[])
    input_path = _structure_input_path(input)
    try
        cmd = Cmd(vcat(String(command), String.(collect(args)), String(input_path)))
        return read(cmd, String)
    finally
        input isa Structure && isfile(input_path) && rm(input_path, force=true)
    end
end

function run_pdb2pqr(input::Union{AbstractString,Structure}, output_path::AbstractString; command::AbstractString="pdb2pqr30", args::AbstractVector{<:AbstractString}=String[])
    input_path = _structure_input_path(input)
    try
        cmd = Cmd(vcat(String(command), String.(collect(args)), String(input_path), String(output_path)))
        run(cmd)
        return output_path
    finally
        input isa Structure && isfile(input_path) && rm(input_path, force=true)
    end
end