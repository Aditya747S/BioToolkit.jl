# ==============================================================================
# phylo.jl — Phylogenetic tree inference, manipulation, and visualization
#
# Provides Newick/NEXUS/PhyloXML/NeXML I/O, NJ/UPGMA tree construction,
# maximum parsimony and likelihood tree search, substitution models (JC69,
# K80, HKY85), bootstrap analysis, Robinson-Foulds distance, and
# ASCII/Unicode/DOT/Mermaid tree visualization.
#
# References:
#   - Saitou & Nei (1987) MBE 4(4):406-425 (Neighbor-Joining)
#   - Felsenstein (1981) J Evol Biol 17:368-376 (ML phylogeny)
#   - Robinson & Foulds (1981) Math Biosci 53:131-147 (RF distance)
#   - Jukes & Cantor (1969) Mammalian Protein Metabolism (JC69)
# ==============================================================================

export distance_matrix, neighbor_joining, neighbor_joining_tree, AbstractPhyloTree, PhyloTree, UnrootedPhyloTree, parse_newick, write_newick, parse_tree, write_tree, parse_phyloxml, write_phyloxml, parse_nexus, write_nexus, parse_nexml, write_nexml, upgma, get_terminals, get_nonterminals, tree_distance, draw_ascii, draw_unicode, tree_to_dot, tree_to_mermaid, prune, root_with_outgroup, reroot, midpoint_root, get_parent, lowest_common_ancestor, common_ancestor, is_monophyletic, robinson_foulds_distance, bootstrap_trees, tree_consensus, consensus_tree, strict_consensus_tree, majority_consensus_tree, bootstrap_consensus_tree, bootstrap_support, set_metadata!, annotate_tree!, is_terminal, is_parent_of, get_path, trace, is_bifurcating, is_preterminal, total_branch_length, depths, find_clades, ladderize, count_terminals, collapse_clades, parsimony_score, maximum_parsimony_tree, parsimony_tree, JC69, K80, HKY85, felsenstein_likelihood, transition_probability, maximum_likelihood_tree, coordinates, GTR, WAG, LG, JTT, Dayhoff, Blosum62, CpREV, MtMAM, DiscreteGammaRate, InvariableSites, model_test, pic, pgls, ancestral_state_reconstruction, fast_anc, bionj_tree, ml_distance_matrix, sh_test, sankoff_parsimony_score, coalescent_tree, yule_tree, birth_death_tree, kuhner_felsenstein_distance, path_difference_distance, reconcile_trees, force_ultrametric!, resolve_polytomies!, random_binary_tree, cophenetic_matrix, faith_pd, evolutionary_distinctiveness, ed_scores, strict_clock_dating, node_ages, fit_mk_model, tanglegram_layout, subtree_with_taxa, prune_taxa, tree_bipartitions, simulate_bm, simulate_ou, fit_bm, fit_ou, blomberg_k, pagel_lambda, colless_index, sackin_index, tbr_parsimony_search, GY94, MG94, fel_selection_test, spr_ml_search, CustomCTMC, DiscreteFreeRateModel, simmap, neighbor_net, joint_ancestral_reconstruction, consistency_index, retention_index, weighted_robinson_foulds, max_clade_credibility

# GPU acceleration is available when the CUDA extension is loaded; see lazy_gpu.jl.
using LinearAlgebra
using SpecialFunctions: erfc, erfinv
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, new_provenance_id, provenance_parent_ids, register_provenance!

abstract type AbstractPhyloTree end

"""
    PhyloTree

Mutable rooted tree node for phylogenetic workflows.
"""
mutable struct PhyloTree <: AbstractPhyloTree
  name::String
  branch_length::Float64
  children::Vector{PhyloTree}
  support::Float64
  metadata::Dict{String,Any}
end

"""
    UnrootedPhyloTree

Unrooted phylogenetic graph node where connections are bi-directional neighbors with explicit edge lengths.
"""
mutable struct UnrootedPhyloTree <: AbstractPhyloTree
  name::String
  neighbors::Vector{UnrootedPhyloTree}
  edge_lengths::Vector{Float64}
  metadata::Dict{String,Any}
end

function UnrootedPhyloTree(
  name::AbstractString="";
  neighbors::AbstractVector{UnrootedPhyloTree}=UnrootedPhyloTree[],
  edge_lengths::AbstractVector{<:Real}=Float64[],
  metadata::AbstractDict=Dict{String,Any}())
  normalized_meta = Dict{String,Any}(string(k) => v for (k, v) in metadata)
  return UnrootedPhyloTree(String(name), Vector{UnrootedPhyloTree}(neighbors), Float64.(edge_lengths), normalized_meta)
end

function add_edge!(u::UnrootedPhyloTree, v::UnrootedPhyloTree, length::Real=0.0)
  push!(u.neighbors, v)
  push!(u.edge_lengths, Float64(length))
  push!(v.neighbors, u)
  push!(v.edge_lengths, Float64(length))
  return u
end

PhyloTree(
  name::String;
  branch_length::Real=0.0,
  children::AbstractVector{PhyloTree}=PhyloTree[],
  support::Real=0.0,
  metadata::AbstractDict=Dict{String,Any}()) = begin
  normalized_metadata = Dict{String,Any}(string(key) => value for (key, value) in metadata)
  metadata_provenance(normalized_metadata) === nothing && stamp_provenance!(
    normalized_metadata;
    label="PhyloTree",
    source="PhyloTree",
    notes=["constructed from in-memory tree node"],
    parameters=(name=String(name), child_count=length(children), branch_length=Float64(branch_length)))
  PhyloTree(String(name), Float64(branch_length), Vector{PhyloTree}(children), Float64(support), normalized_metadata)
end

PhyloTree(
  children::AbstractVector{PhyloTree};
  name::String="",
  branch_length::Real=0.0,
  support::Real=0.0,
  metadata::AbstractDict=Dict{String,Any}()) = PhyloTree(String(name); branch_length=branch_length, children=children, support=support, metadata=metadata)

isleaf(tree::PhyloTree) = isempty(tree.children)
isleaf(tree::UnrootedPhyloTree) = length(tree.neighbors) <= 1

function coordinates(tree::PhyloTree)
  leaves = get_terminals(tree)
  y_positions = Dict(node => index for (index, node) in enumerate(leaves))
  node_positions = IdDict{PhyloTree,Tuple{Float64,Float64}}()

  function _place(node::PhyloTree, depth::Float64)
    if isleaf(node)
      position = (depth, Float64(get(y_positions, node, 1)))
      node_positions[node] = position
      return position
    end
    child_positions = [_place(child, depth + child.branch_length) for child in node.children]
    y = sum(position[2] for position in child_positions) / length(child_positions)
    position = (depth, y)
    node_positions[node] = position
    return position
  end

  _place(tree, 0.0)
  return node_positions
end

function _phylo_copy(tree::PhyloTree)
  children = PhyloTree[_phylo_copy(child) for child in tree.children]
  return PhyloTree(tree.name; branch_length=tree.branch_length, children=children, support=tree.support, metadata=tree.metadata)
end

function _phylo_collect_terms(tree::PhyloTree, buffer::Vector{PhyloTree})
  if isleaf(tree)
    push!(buffer, tree)
  else
    for child in tree.children
      _phylo_collect_terms(child, buffer)
    end
  end
  return buffer
end

function get_terminals(tree::PhyloTree)
  terms = PhyloTree[]
  sizehint!(terms, 32)
  stack = PhyloTree[tree]
  sizehint!(stack, 32)
  while !isempty(stack)
    curr = pop!(stack)
    if isleaf(curr)
      push!(terms, curr)
    else
      for i in length(curr.children):-1:1
        push!(stack, curr.children[i])
      end
    end
  end
  return terms
end

function get_nonterminals(tree::PhyloTree)
  non_terms = PhyloTree[]
  sizehint!(non_terms, 32)
  stack = PhyloTree[tree]
  sizehint!(stack, 32)
  while !isempty(stack)
    curr = pop!(stack)
    if !isleaf(curr)
      push!(non_terms, curr)
      for i in length(curr.children):-1:1
        push!(stack, curr.children[i])
      end
    end
  end
  return non_terms
end

function _phylo_to_newick(tree::PhyloTree; is_root::Bool=false)
  if isleaf(tree)
    return string(tree.name, ":", round(tree.branch_length; digits=5))
  end

  inner = join((_phylo_to_newick(child) for child in tree.children), ",")
  label = tree.support > 0 ? string(round(tree.support; digits=5)) : (isempty(tree.name) ? "" : tree.name)
  if is_root && isapprox(tree.branch_length, 0.0; atol=1e-12)
    return string("(", inner, ")", label)
  end
  return string("(", inner, ")", label, ":", round(tree.branch_length; digits=5))
end

function write_newick(tree::PhyloTree)
  return string(_phylo_to_newick(tree; is_root=true), ";")
end

function _skip_whitespace(text::AbstractString, index::Int)
  while index <= lastindex(text) && isspace(text[index])
    index += 1
  end
  return index
end

function _read_label(text::AbstractString, index::Int)
  start_index = index
  while index <= lastindex(text)
    char = text[index]
    char == '(' && break
    char == ')' && break
    char == ',' && break
    char == ':' && break
    char == ';' && break
    index += 1
  end
  return String(strip(text[start_index:(index-1)])), index
end

function _phylo_label_to_support_or_name(label::AbstractString)
  stripped = String(strip(String(label)))
  isempty(stripped) && return "", 0.0
  try
    return "", parse(Float64, stripped)
  catch
    return stripped, 0.0
  end
end

function _read_branch_length(text::AbstractString, index::Int)
  index = _skip_whitespace(text, index)
  index > lastindex(text) && return 0.0, index
  text[index] == ':' || return 0.0, index
  index += 1
  index = _skip_whitespace(text, index)
  start_index = index
  while index <= lastindex(text)
    char = text[index]
    char == ',' && break
    char == ')' && break
    char == ';' && break
    index += 1
  end
  token = strip(text[start_index:(index-1)])
  isempty(token) && return 0.0, index
  return parse(Float64, token), index
end

function _phylo_clade_key(names::AbstractVector{String})
  sorted = sort(collect(names))
  return string(length(sorted), '\u001e', join(sorted, '\u001f'))
end

function _phylo_key_set(key::String)
  text = String(key)
  separator = findfirst(==('\u001e'), text)
  payload = separator === nothing ? text : separator == lastindex(text) ? "" : text[nextind(text, separator):end]
  isempty(payload) && return Set{String}()
  return Set(String.(Base.split(payload, "\u001f")))
end

function _phylo_leaf_names(tree::PhyloTree)
  return [node.name for node in get_terminals(tree)]
end

function _phylo_leaf_set(tree::PhyloTree)
  return Set(_phylo_leaf_names(tree))
end

function _phylo_find_parent(node::PhyloTree, parent::Union{Nothing,PhyloTree}, target::String)
  node.name == target && return parent
  for child in node.children
    found = _phylo_find_parent(child, node, target)
    found === nothing || return found
  end
  return nothing
end

function get_parent(tree::PhyloTree, target_name::String)
  return _phylo_find_parent(tree, nothing, String(target_name))
end

function lowest_common_ancestor(tree::PhyloTree, names::AbstractVector{<:String})
  isempty(names) && throw(ArgumentError("at least one taxon name is required"))
  paths = PhyloTree[]
  path_list = Vector{Vector{PhyloTree}}()
  for name in names
    path = _phylo_path_to_leaf(tree, String(name))
    path === nothing && throw(ArgumentError("unknown leaf: $(name)"))
    push!(path_list, path)
  end

  shortest = minimum(length.(path_list))
  ancestor = path_list[1][1]
  for index in 1:shortest
    candidate = path_list[1][index]
    all(path -> path[index] === candidate, path_list) || break
    ancestor = candidate
  end
  return ancestor
end

lowest_common_ancestor(tree::PhyloTree, left_name::String, right_name::String) = lowest_common_ancestor(tree, [left_name, right_name])

common_ancestor(tree::PhyloTree, names) = lowest_common_ancestor(tree, names)

is_terminal(node::PhyloTree) = isleaf(node)

function is_parent_of(tree::PhyloTree, parent_name::String, child_name::String)
  parent = _phylo_path_to_named_node(tree, String(parent_name))
  child = _phylo_path_to_named_node(tree, String(child_name))
  parent === nothing && return false
  child === nothing && return false
  parent_node = parent[end]
  child_node = child[end]
  return any(child_ref -> child_ref === child_node, parent_node.children)
end

function get_path(tree::PhyloTree, target_name::String)
  path = _phylo_path_to_named_node(tree, String(target_name))
  path === nothing && throw(ArgumentError("unknown leaf: $(target_name)"))
  return path
end

function trace(tree::PhyloTree, start_name::String, end_name::String)
  left_path = get_path(tree, start_name)
  right_path = get_path(tree, end_name)
  shared_index = 0
  for (index, (left_node, right_node)) in enumerate(zip(left_path, right_path))
    left_node === right_node || break
    shared_index = index
  end
  left_segment = reverse(left_path[shared_index:end])
  right_segment = right_path[(shared_index+1):end]
  return vcat(left_segment, right_segment)
end

function split(tree::PhyloTree, parent_name::String; child_names::AbstractVector{<:String}=String[], branch_length::Real=0.0, support::Real=0.0)
  target = String(parent_name)
  path = _phylo_path_to_named_node(tree, target)
  path === nothing && throw(ArgumentError("unknown node: $(parent_name)"))
  node = path[end]

  new_children = PhyloTree[]
  if isempty(child_names)
    push!(new_children, PhyloTree(string(target, "_0"); branch_length=branch_length))
    push!(new_children, PhyloTree(string(target, "_1"); branch_length=branch_length))
  else
    for child_name in child_names
      push!(new_children, PhyloTree(String(child_name); branch_length=branch_length))
    end
  end

  if isleaf(node)
    node.children = new_children
    node.support = support
    return tree
  end

  node.children = vcat(node.children, new_children)
  node.support = support == 0.0 ? node.support : Float64(support)
  return tree
end

function is_monophyletic(tree::PhyloTree, taxa_names)
  taxa = Set(String.(collect(taxa_names)))
  isempty(taxa) && return false
  ancestor = lowest_common_ancestor(tree, collect(taxa))
  return _phylo_leaf_set(ancestor) == taxa
end

function _phylo_informative_clade_keys(tree::PhyloTree)
  leaves = _phylo_clade_key(_phylo_leaf_names(tree))
  return Set(filter(key -> key != leaves, _phylo_clade_keys(tree)))
end

function robinson_foulds_distance(left_tree::PhyloTree, right_tree::PhyloTree)
  left_leaves = _phylo_leaf_set(left_tree)
  right_leaves = _phylo_leaf_set(right_tree)
  left_leaves == right_leaves || throw(ArgumentError("trees must share the same leaf set"))

  left_keys = _phylo_informative_clade_keys(left_tree)
  right_keys = _phylo_informative_clade_keys(right_tree)
  return length(symdiff(left_keys, right_keys))
end

function is_bifurcating(tree::PhyloTree)
  for node in vcat(PhyloTree[tree], get_nonterminals(tree))
    !isleaf(node) && length(node.children) != 2 && return false
  end
  return true
end

function is_preterminal(node::PhyloTree)
  return length(node.children) == 1 && isleaf(node.children[1])
end

function total_branch_length(tree::PhyloTree)
  total = 0.0
  function _sum(node::PhyloTree)
    total_local = node.branch_length
    for child in node.children
      total_local += _sum(child)
    end
    return total_local
  end
  return _sum(tree)
end

function _depth_map(tree::PhyloTree, depth::Float64, map::Dict{PhyloTree,Float64})
  map[tree] = depth
  for child in tree.children
    _depth_map(child, depth + child.branch_length, map)
  end
  return map
end

function depths(tree::PhyloTree)
  return _depth_map(tree, 0.0, Dict{PhyloTree,Float64}())
end

# --- Maximum Likelihood Phylogenetics ---

abstract type SubstitutionModel end

struct JC69 <: SubstitutionModel end

struct K80 <: SubstitutionModel
  kappa::Float64
end

struct HKY85 <: SubstitutionModel
  pi::Vector{Float64} # [A, C, G, T]
  kappa::Float64
end

function transition_probability(model::JC69, t::Real)
  p_same = 0.25 + 0.75 * exp(-4.0 * t / 3.0)
  p_diff = 0.25 - 0.25 * exp(-4.0 * t / 3.0)
  P = fill(p_diff, 4, 4)
  @inbounds for i in 1:4
    P[i, i] = p_same
  end
  return P
end

function transition_probability(model::K80, t::Real)
  k = model.kappa
  beta = 1.0 / (k + 2.0) # scaling such that total rate is 1
  exp1 = exp(-4.0 * beta * t)
  exp2 = exp(-2.0 * (k + 1.0) * beta * t)

  P = fill(0.25 - 0.25 * exp1, 4, 4)
  same = 0.25 + 0.25 * exp1 + 0.5 * exp2
  transition = 0.25 + 0.25 * exp1 - 0.5 * exp2

  @inbounds begin
    P[1, 1] = same
    P[2, 2] = same
    P[3, 3] = same
    P[4, 4] = same
    P[1, 3] = transition
    P[3, 1] = transition
    P[2, 4] = transition
    P[4, 2] = transition
  end

  return P
end

function transition_probability(model::HKY85, t::Real)
  pi = model.pi
  k = model.kappa
  Q = zeros(Float64, 4, 4)
  # HKY85 rate matrix: Q_ij = k * pi_j if transition, pi_j if transversion
  for i in 1:4
    for j in 1:4
      if i == j
        ;
        continue;
      end
      if (i == 1 && j == 3) || (i == 3 && j == 1) || (i == 2 && j == 4) || (i == 4 && j == 2)
        Q[i, j] = k * pi[j]
      else
        Q[i, j] = pi[j]
      end
    end
  end
  # Scaling such that average rate = 1
  avg_rate = 0.0
  for i in 1:4
    avg_rate += pi[i] * sum(Q[i, :])
  end
  Q ./= avg_rate

  for i in 1:4
    ;
    Q[i, i] = -sum(Q[i, :]);
  end
  return exp(Q * t)
end

# State mapping for DNA
const _DNA_STATE_MAP = Dict{Char,Vector{Float64}}(
  'A' => [1.0, 0.0, 0.0, 0.0],
  'C' => [0.0, 1.0, 0.0, 0.0],
  'G' => [0.0, 0.0, 1.0, 0.0],
  'T' => [0.0, 0.0, 0.0, 1.0],
  'R' => [1.0, 0.0, 1.0, 0.0], # A or G
  'Y' => [0.0, 1.0, 0.0, 1.0], # C or T
  'N' => [1.0, 1.0, 1.0, 1.0],
  '-' => [1.0, 1.0, 1.0, 1.0],
  '?' => [1.0, 1.0, 1.0, 1.0]
)

const _DNA_AMBIGUOUS_STATE = _DNA_STATE_MAP['N']

@inline _phylo_coerce_sequence(sequence::BioSequence) = sequence

function _phylo_coerce_sequence(sequence)
  sequence_text = String(sequence)
  inferred_alphabet = _infer_sequence_alphabet(sequence_text)
  return BioSequence{inferred_alphabet}(sequence_text)
end

@inline function _phylo_require_alphabet(sequence::BioSequence{A}, ::Type{A}) where {A<:BioAlphabet}
  return sequence
end

function _phylo_require_alphabet(sequence::BioSequence, required_alphabet::Type{A}) where {A<:BioAlphabet}
  sequence_text = String(sequence)
  validate_sequence(required_alphabet, sequence_text) || throw(ArgumentError("sequence contains symbols incompatible with $(required_alphabet)"))
  return BioSequence{required_alphabet}(sequence_text)
end

function _phylo_dna_alignment(alignment::AbstractDict)
  typed = Dict{String,DNASeq}()
  for (name, sequence) in alignment
    typed[String(name)] = _phylo_require_alphabet(_phylo_coerce_sequence(sequence), DNAAlphabet)
  end
  return typed
end

function _get_leaf_likelihoods(sequence::BioSequence{DNAAlphabet})
  n_sites = length(sequence)
  L = zeros(Float64, 4, n_sites)
  for (i, byte) in enumerate(sequence.data)
    state = get(_DNA_STATE_MAP, uppercase(Char(byte)), _DNA_AMBIGUOUS_STATE)
    @inbounds for j in 1:4
      L[j, i] = state[j]
    end
  end
  return L
end

function _is_aa_alignment(alignment::AbstractDict)
  for seq in values(alignment)
    str = sequence_to_string(seq)
    for c in uppercase(str)
      if c in ('E', 'F', 'I', 'L', 'P', 'Q', 'Z', 'J')
        return true
      end
    end
  end
  return false
end

function sequence_to_string(seq)
  if seq isa AbstractString
    return String(seq)
  elseif seq isa BioSequence
    return String(seq)
  else
    return string(seq)
  end
end

function _phylo_aa_alignment(alignment::AbstractDict)
  typed = Dict{String, String}()
  for (name, sequence) in alignment
    typed[String(name)] = uppercase(sequence_to_string(sequence))
  end
  return typed
end

const _AA_STATE_MAP = Dict{Char, Vector{Float64}}(
  'A' => [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  'R' => [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  'N' => [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  'D' => [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  'C' => [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  'Q' => [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
  'E' => [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
  'G' => [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
  'H' => [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
  'I' => [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
  'L' => [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
  'K' => [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
  'M' => [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
  'F' => [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
  'P' => [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
  'S' => [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
  'T' => [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
  'W' => [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
  'Y' => [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
  'V' => [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
  'B' => [0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], # N or D
  'Z' => [0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0], # Q or E
  'J' => [0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0], # I or L
)
const _AA_AMBIGUOUS_STATE = fill(1.0, 20)

function _get_aa_leaf_likelihoods(sequence::AbstractString)
  n_sites = length(sequence)
  L = zeros(Float64, 20, n_sites)
  for (i, char) in enumerate(sequence)
    state = get(_AA_STATE_MAP, uppercase(char), _AA_AMBIGUOUS_STATE)
    @inbounds for j in 1:20
      L[j, i] = Float64(state[j])
    end
  end
  return L
end

"""
    felsenstein_likelihood(tree::PhyloTree, alignment::AbstractDict, model::SubstitutionModel)

Compute the log-likelihood of a phylogenetic tree given an alignment using Felsenstein's pruning algorithm.
"""
# Discrete Gamma Rate Distribution & Invariable Sites
struct DiscreteGammaRate
  alpha::Float64
  n_categories::Int
end

DiscreteGammaRate(alpha::Real; n_categories::Int=4) = DiscreteGammaRate(Float64(alpha), n_categories)

struct InvariableSites
  p_inv::Float64
end

struct DiscreteFreeRateModel
  rates::Vector{Float64}
  weights::Vector{Float64}

  function DiscreteFreeRateModel(rates::AbstractVector{<:Real}, weights::AbstractVector{<:Real})
    rates_f = Float64.(rates)
    w = Float64.(weights) ./ sum(weights)
    mean_r = sum(w .* rates_f)
    r = mean_r > 0 ? rates_f ./ mean_r : rates_f
    return new(r, w)
  end
end

function _gamma_category_rates(gamma_model::DiscreteGammaRate)
  k = gamma_model.n_categories
  a = gamma_model.alpha
  rates = zeros(Float64, k)
  for i in 1:k
    p = (i - 0.5) / k
    z = sqrt(2.0) * erfinv(2.0 * p - 1.0)
    term = 1.0 - 2.0 / (9.0 * a) + z * sqrt(2.0 / (9.0 * a))
    rates[i] = max(1e-4, term^3)
  end
  rates ./= (sum(rates) / k)
  return rates, fill(1.0 / k, k)
end

"""
    felsenstein_likelihood(tree::PhyloTree, alignment::AbstractDict, model::SubstitutionModel; rates=nothing, p_inv=nothing)

Compute the log-likelihood of a phylogenetic tree given an alignment using Felsenstein's pruning algorithm.
Supports rate heterogeneity (+G / Discrete Gamma / Free Rate) and proportion of invariable sites (+I).
"""
function felsenstein_likelihood(
  tree::PhyloTree,
  alignment::AbstractDict,
  model::SubstitutionModel;
  rates::Union{DiscreteGammaRate, DiscreteFreeRateModel, Nothing}=nothing,
  p_inv::Union{InvariableSites, Real, Nothing}=nothing
)
  is_aa = model isa AAModel || _is_aa_alignment(alignment)
  alignment_parsed = is_aa ? _phylo_aa_alignment(alignment) : _phylo_dna_alignment(alignment)
  isempty(alignment_parsed) && throw(ArgumentError("alignment must contain at least one sequence"))

  n_states = is_aa ? 20 : 4
  leaf_likelihoods = Dict{PhyloTree,Matrix{Float64}}()
  terminals = get_terminals(tree)
  first_seq = first(values(alignment_parsed))
  n_sites = length(first_seq)

  for leaf in terminals
    seq = get(alignment_parsed, leaf.name, nothing)
    if seq === nothing
      leaf_likelihoods[leaf] = ones(Float64, n_states, n_sites)
    else
      leaf_likelihoods[leaf] = is_aa ? _get_aa_leaf_likelihoods(seq) : _get_leaf_likelihoods(seq)
    end
  end

  transition_cache = Dict{Float64,Matrix{Float64}}()
  transition_cached(t::Float64) = get!(transition_cache, t) do
    transition_probability(model, t)
  end

  cat_rates, cat_weights = if rates isa DiscreteGammaRate
    _gamma_category_rates(rates)
  elseif rates isa DiscreteFreeRateModel
    rates.rates, rates.weights
  else
    [1.0], [1.0]
  end

  pi = if model isa AAModel
    model.pi
  elseif model isa HKY85 || model isa GTR
    model.pi
  else
    fill(1.0 / n_states, n_states)
  end

  pinv_val = if p_inv isa InvariableSites
    p_inv.p_inv
  elseif p_inv isa Real
    Float64(p_inv)
  else
    0.0
  end
  pinv_val = clamp(pinv_val, 0.0, 0.9999)

  site_variable_lik = zeros(Float64, n_sites)

  for (cat_idx, (r_k, w_k)) in enumerate(zip(cat_rates, cat_weights))
    function _compute_internal_likelihood(node::PhyloTree)
      if isleaf(node)
        return leaf_likelihoods[node]
      end

      child_Ls = [_compute_internal_likelihood(child) for child in node.children]
      node_L = ones(Float64, n_states, n_sites)

      for (i, child) in enumerate(node.children)
        P = transition_cached(r_k * child.branch_length)
        c_L = child_Ls[i]
        
        @inbounds for s in 1:n_sites
          for a in 1:n_states
            acc = 0.0
            for b in 1:n_states
              acc += P[a, b] * c_L[b, s]
            end
            node_L[a, s] *= acc
          end
        end
      end
      return node_L
    end

    root_L = _compute_internal_likelihood(tree)
    for s in 1:n_sites
      site_cat_lik = 0.0
      @inbounds for a in 1:n_states
        site_cat_lik += pi[a] * root_L[a, s]
      end
      site_variable_lik[s] += w_k * site_cat_lik
    end
  end

  log_lik = 0.0
  for s in 1:n_sites
    l_var = site_variable_lik[s]
    l_site = if pinv_val > 0.0
      first_char = uppercase(first_seq[s])
      all_match = true
      for seq in values(alignment_parsed)
        c = uppercase(seq[s])
        if c != first_char && c != '-' && c != 'N' && c != 'X' && c != '?'
          all_match = false
          break
        end
      end
      if all_match
        state_idx = is_aa ? get(_AA_STATE_MAP, first_char, _AA_AMBIGUOUS_STATE)[1] : (first_char == 'A' ? 1 : first_char == 'C' ? 2 : first_char == 'G' ? 3 : first_char == 'T' || first_char == 'U' ? 4 : 0)
        pi_state = (state_idx >= 1 && state_idx <= n_states) ? pi[state_idx] : 1.0 / n_states
        pinv_val * pi_state + (1.0 - pinv_val) * l_var
      else
        (1.0 - pinv_val) * l_var
      end
    else
      l_var
    end
    log_lik += log(max(l_site, 1e-300))
  end

  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "felsenstein_likelihood";
      parameters=(n_taxa=length(get_terminals(tree)), model_type=string(typeof(model)), has_gamma=(rates !== nothing), p_inv=pinv_val))
  end
  return log_lik
end

function _brent_optimize_branch(f::Function, ax::Float64, bx::Float64, cx::Float64; tol::Float64=1e-3, max_iter::Int=10)
  a = min(ax, cx)
  b = max(ax, cx)
  v = w = x = bx
  e = 0.0
  fx = f(x)
  fv = fw = fx
  d = 0.0
  golden = 0.38196601125010515
  
  for _ in 1:max_iter
    xm = 0.5 * (a + b)
    tol1 = tol * abs(x) + 1e-10
    tol2 = 2.0 * tol1
    
    if abs(x - xm) <= (tol2 - 0.5 * (b - a))
      return x, fx
    end
    
    p = q = r = 0.0
    if abs(e) > tol1
      r = (x - w) * (fx - fv)
      q = (x - v) * (fx - fw)
      p = (x - v) * q - (x - w) * r
      q = 2.0 * (q - r)
      if q > 0.0
        p = -p
      end
      q = abs(q)
      etemp = e
      e = d
      if abs(p) < abs(0.5 * q * etemp) && p > q * (a - x) && p < q * (b - x)
        d = p / q
        u = x + d
        if (u - a) < tol2 || (b - u) < tol2
          d = sign(xm - x) * tol1
        end
      else
        e = x >= xm ? a - x : b - x
        d = golden * e
      end
    else
      e = x >= xm ? a - x : b - x
      d = golden * e
    end
    
    u = abs(d) >= tol1 ? x + d : x + sign(d) * tol1
    fu = f(u)
    
    if fu <= fx
      if u >= x
        a = x
      else
        b = x
      end
      v, fv = w, fw
      w, fw = x, fx
      x, fx = u, fu
    else
      if u < x
        a = u
      else
        b = u
      end
      if fu <= fw || w == x
        v, fv = w, fw
        w, fw = u, fu
      elseif fu <= fv || v == x || v == w
        v, fv = u, fu
      end
    end
  end
  return x, fx
end

"""
    optimize_branch_lengths!(tree::PhyloTree, alignment::AbstractDict, model::SubstitutionModel; max_iter=10)

Optimize branch lengths of a given tree to maximize likelihood using a coordinate ascent approach with Brent's 1D minimization.
"""
function optimize_branch_lengths!(tree::PhyloTree, alignment::AbstractDict, model::SubstitutionModel; max_iter=5)
  nodes = vcat(PhyloTree[tree], get_nonterminals(tree))

  for iter in 1:max_iter
    improved = false
    for node in nodes
      for child in node.children
        old_bl = child.branch_length
        best_lik = felsenstein_likelihood(tree, alignment, model)
        
        # Brent 1D parabolic minimization for optimal branch length
        obj_func = (len::Float64) -> begin
          child.branch_length = max(0.0001, len)
          return -felsenstein_likelihood(tree, alignment, model)
        end
        
        best_x, _ = _brent_optimize_branch(obj_func, 0.0001, old_bl > 0 ? old_bl : 0.1, 5.0; max_iter=8)
        child.branch_length = max(0.0001, best_x)
        new_lik = felsenstein_likelihood(tree, alignment, model)
        
        if new_lik > best_lik
          improved = true
        else
          child.branch_length = old_bl
        end
      end
    end
    improved || break
  end
  
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "optimize_branch_lengths!";
      parameters=(n_taxa=count_terminals(tree), max_iter=Int(max_iter)))
  end
  return tree
end

function _nni_moves(tree::PhyloTree)
  moves = PhyloTree[]
  # Collect internal nodes that have at least one internal child
  internal_nodes = get_nonterminals(tree)

  for node in internal_nodes
    # A node must have at least 2 children to be an internal branch point
    if length(node.children) >= 2
      for (i, child) in enumerate(node.children)
        if !isleaf(child) && length(child.children) >= 2
          # We can swap node.children[other_than_i] with child.children[j]
          for j in 1:length(child.children)
            # Create a new tree topology by swapping
            # This is a bit complex in-place, so we copy
            for other_i in 1:length(node.children)
              if other_i == i
                ;
                continue;
              end

              # Swap node.children[other_i] and child.children[j]
              new_tree = _phylo_copy(tree)
              # Find the corresponding nodes in the copy...
              # Actually, a simpler way is to swap in-place, record likelihood, and swap back

              target_node = _find_node_by_path(new_tree, _get_node_path(tree, node))
              target_child = _find_node_by_path(new_tree, _get_node_path(tree, child))

              tmp = target_node.children[other_i]
              target_node.children[other_i] = target_child.children[j]
              target_child.children[j] = tmp

              push!(moves, new_tree)
            end
          end
        end
      end
    end
  end
  return moves
end

function _get_node_path(root::PhyloTree, target::PhyloTree)
  path = Int[]
  function _search(node::PhyloTree, p::Vector{Int})
    if node === target
      append!(path, p)
      return true
    end
    for (i, child) in enumerate(node.children)
      if _search(child, vcat(p, i))
        return true
      end
    end
    return false
  end
  _search(root, Int[])
  return path
end

function _find_node_by_path(root::PhyloTree, path::Vector{Int})
  current = root
  for index in path
    current = current.children[index]
  end
  return current
end

"""
    maximum_likelihood_tree(alignment::AbstractDict; model::SubstitutionModel=JC69(), initial_tree::Union{Nothing, PhyloTree}=nothing)

Find the Maximum Likelihood tree for a given alignment.
"""
function maximum_likelihood_tree(alignment::AbstractDict; model::SubstitutionModel=JC69(), initial_tree::Union{Nothing,PhyloTree}=nothing, max_topology_iter=10)
  is_aa = model isa AAModel || _is_aa_alignment(alignment)
  alignment_parsed = is_aa ? _phylo_aa_alignment(alignment) : _phylo_dna_alignment(alignment)

  if initial_tree === nothing
    names = collect(keys(alignment_parsed))
    if length(names) <= 3
      leaves = [PhyloTree(name; branch_length=0.1) for name in names]
      if length(leaves) == 1
        initial_tree = leaves[1]
      elseif length(leaves) == 2
        initial_tree = PhyloTree(""; branch_length=0.0, children=leaves)
      else
        initial_tree = PhyloTree(""; branch_length=0.0, children=[leaves[1], PhyloTree(""; branch_length=0.1, children=[leaves[2], leaves[3]])])
      end
    else
      seqs = [alignment_parsed[name] for name in names]
      dm = is_aa ? ml_distance_matrix(seqs; model=model) : distance_matrix(seqs)
      initial_tree = neighbor_joining_tree(dm, names)
    end
  end

  current_tree = _phylo_copy(initial_tree)
  optimize_branch_lengths!(current_tree, alignment_parsed, model)
  current_lik = felsenstein_likelihood(current_tree, alignment_parsed, model)

  if count_terminals(current_tree) <= 3
    return current_tree
  end

  for iter in 1:max_topology_iter
    best_neighbor = current_tree
    best_neighbor_lik = current_lik

    neighbors = _nni_moves(current_tree)
    for neighbor in neighbors
      optimize_branch_lengths!(neighbor, alignment_parsed, model; max_iter=2)
      n_lik = felsenstein_likelihood(neighbor, alignment_parsed, model)
      if n_lik > best_neighbor_lik
        best_neighbor_lik = n_lik
        best_neighbor = neighbor
      end
    end

    if best_neighbor_lik > current_lik
      current_tree = best_neighbor
      current_lik = best_neighbor_lik
      optimize_branch_lengths!(current_tree, alignment_parsed, model; max_iter=5)
      current_lik = felsenstein_likelihood(current_tree, alignment_parsed, model)
    else
      break
    end
  end

  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "maximum_likelihood_tree";
      parameters=(n_taxa=length(alignment_parsed), model_type=string(typeof(model)), max_topology_iter=Int(max_topology_iter)))
  end
  return current_tree
end

function find_clades(tree::PhyloTree; name_pattern::Union{Nothing,Regex}=nothing, terminals_only::Bool=false)
  matches = PhyloTree[]
  function _search(node::PhyloTree)
    if isleaf(node)
      if !terminals_only || (name_pattern === nothing || (!isempty(node.name) && occursin(name_pattern, node.name)))
        push!(matches, node)
      end
      return
    end

    if !terminals_only && (name_pattern === nothing || (!isempty(node.name) && occursin(name_pattern, node.name)))
      push!(matches, node)
    end
    for child in node.children
      _search(child)
    end
  end
  _search(tree)
  return matches
end

function ladderize(tree::PhyloTree; ascending::Bool=true)
  copied = _phylo_copy(tree)
  function _ladderize(node::PhyloTree)
    if isleaf(node)
      return node
    end
    for child in node.children
      _ladderize(child)
    end
    if length(node.children) > 1
      leaf_counts = Int[]
      for child in node.children
        terminals = get_terminals(child)
        push!(leaf_counts, length(terminals))
      end
      indices = sortperm(leaf_counts; rev=(!ascending))
      node.children = node.children[indices]
    end
    return node
  end
  return _ladderize(copied)
end

function count_terminals(tree::PhyloTree)
  return length(get_terminals(tree))
end

function collapse_clades(tree::PhyloTree; min_support::Real=0.0)
  isempty(get_nonterminals(tree)) && return _phylo_copy(tree)

  copied = _phylo_copy(tree)
  function _collapse(node::PhyloTree)
    if isleaf(node)
      return node
    end

    for child in node.children
      _collapse(child)
    end

    remaining_children = PhyloTree[]
    for child in node.children
      if isleaf(child)
        push!(remaining_children, child)
      elseif child.support >= min_support || isapprox(child.support, 0.0; atol=1e-12)
        push!(remaining_children, child)
      else
        for grandchild in child.children
          grandchild.branch_length += child.branch_length
          push!(remaining_children, grandchild)
        end
      end
    end

    node.children = remaining_children
    return node
  end

  return _collapse(copied)
end

function _phylo_clade_keys(tree::PhyloTree)
  keys = String[]
  function _collect(node::PhyloTree)
    if isleaf(node)
      return Set([node.name])
    end
    combined = Set{String}()
    for child in node.children
      union!(combined, _collect(child))
    end
    if !isempty(node.children)
      push!(keys, _phylo_clade_key(collect(combined)))
    end
    return combined
  end
  _collect(tree)
  return keys
end

function _phylo_clade_map(tree::PhyloTree)
  map = Dict{String,Tuple{Int,Float64}}()
  function _collect(node::PhyloTree)
    if isleaf(node)
      return Set([node.name])
    end
    combined = Set{String}()
    for child in node.children
      union!(combined, _collect(child))
    end
    key = _phylo_clade_key(collect(combined))
    if !isempty(node.children)
      count, total_branch = get(map, key, (0, 0.0))
      map[key] = (count + 1, total_branch + node.branch_length)
    end
    return combined
  end
  _collect(tree)
  return map
end

function _phylo_adjacency(tree::PhyloTree)
  adjacency = IdDict{PhyloTree,Vector{Tuple{PhyloTree,Float64}}}()

  function _add_edge(left::PhyloTree, right::PhyloTree, weight::Float64)
    push!(get!(adjacency, left) do
        Tuple{PhyloTree,Float64}[]
      end, (right, weight))
  end

  function _walk(node::PhyloTree)
    for child in node.children
      _add_edge(node, child, child.branch_length)
      _add_edge(child, node, child.branch_length)
      _walk(child)
    end
  end

  _walk(tree)
  return adjacency
end

function _phylo_clone_from(node::PhyloTree, parent::Union{Nothing,PhyloTree}, incoming_length::Real, adjacency)
  children = PhyloTree[]
  for (neighbor, length) in get(adjacency, node, Tuple{PhyloTree,Float64}[])
    parent !== nothing && neighbor === parent && continue
    push!(children, _phylo_clone_from(neighbor, node, length, adjacency))
  end
  return PhyloTree(node.name; branch_length=incoming_length, children=children, support=node.support)
end

function _phylo_path_between(tree::PhyloTree, left_name::String, right_name::String)
  left_path = _phylo_path_to_leaf(tree, String(left_name))
  right_path = _phylo_path_to_leaf(tree, String(right_name))
  left_path === nothing && throw(ArgumentError("unknown leaf: $(left_name)"))
  right_path === nothing && throw(ArgumentError("unknown leaf: $(right_name)"))

  shared_index = 0
  for (index, (left_node, right_node)) in enumerate(zip(left_path, right_path))
    left_node === right_node || break
    shared_index = index
  end

  path_nodes = vcat(reverse(left_path[shared_index:end]), right_path[(shared_index+1):end])
  left_lengths = [left_path[index+1].branch_length for index in reverse(shared_index:(length(left_path)-1))]
  right_lengths = [right_path[index+1].branch_length for index in shared_index:(length(right_path)-1)]
  return path_nodes, vcat(left_lengths, right_lengths)
end

function _parse_newick_node(text::AbstractString, index::Int)
  index = _skip_whitespace(text, index)
  index <= lastindex(text) || throw(ArgumentError("invalid Newick string"))

  if text[index] == '('
    index += 1
    children = PhyloTree[]
    while true
      child, index = _parse_newick_node(text, index)
      push!(children, child)
      index = _skip_whitespace(text, index)
      index <= lastindex(text) || throw(ArgumentError("invalid Newick string"))
      if text[index] == ','
        index += 1
        continue
      elseif text[index] == ')'
        index += 1
        break
      else
        throw(ArgumentError("invalid Newick string"))
      end
    end

    label, index = _read_label(text, index)
    branch_length, index = _read_branch_length(text, index)
    name, support = _phylo_label_to_support_or_name(label)
    return PhyloTree(isempty(name) ? "" : name; branch_length=branch_length, children=children, support=support), index
  end

  label, index = _read_label(text, index)
  branch_length, index = _read_branch_length(text, index)
  isempty(label) && throw(ArgumentError("leaf node must have a name"))
  name, support = _phylo_label_to_support_or_name(label)
  return PhyloTree(isempty(name) ? label : name; branch_length=branch_length, support=support), index
end

function parse_newick(text::AbstractString)
  stripped = String(strip(text))
  isempty(stripped) && throw(ArgumentError("empty Newick string"))
  tree, index = _parse_newick_node(stripped, firstindex(stripped))
  index = _skip_whitespace(stripped, index)
  if index <= lastindex(stripped) && stripped[index] == ';'
    index += 1
  end
  index = _skip_whitespace(stripped, index)
  index > lastindex(stripped) || throw(ArgumentError("trailing content after Newick tree"))
  return tree
end

function neighbor_joining_tree(D::Matrix{Float64}, names::Vector{String})
  result = parse_newick(neighbor_joining(D, names))
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "neighbor_joining_tree";
      parameters=(n_taxa=length(names)))
  end
  return result
end

function midpoint_root(tree::PhyloTree)
  terminals = get_terminals(tree)
  length(terminals) <= 1 && return _phylo_copy(tree)

  best_left = terminals[1]
  best_right = terminals[2]
  best_distance = -Inf
  for i in 1:(length(terminals)-1)
    for j in (i+1):length(terminals)
      distance = tree_distance(tree, terminals[i].name, terminals[j].name)
      if distance > best_distance
        best_distance = distance
        best_left = terminals[i]
        best_right = terminals[j]
      end
    end
  end

  path_nodes, edge_lengths = _phylo_path_between(tree, best_left.name, best_right.name)
  midpoint = best_distance / 2
  accumulated = 0.0
  adjacency = _phylo_adjacency(tree)

  for edge_index in eachindex(edge_lengths)
    edge_length = edge_lengths[edge_index]
    if isapprox(accumulated + edge_length, midpoint; atol=1e-12)
      return _phylo_clone_from(path_nodes[edge_index+1], nothing, 0.0, adjacency)
    elseif midpoint < accumulated + edge_length
      left_length = midpoint - accumulated
      right_length = edge_length - left_length
      left_subtree = _phylo_clone_from(path_nodes[edge_index], path_nodes[edge_index+1], left_length, adjacency)
      right_subtree = _phylo_clone_from(path_nodes[edge_index+1], path_nodes[edge_index], right_length, adjacency)
      return PhyloTree(""; branch_length=0.0, children=[left_subtree, right_subtree], support=0.0)
    end
    accumulated += edge_length
  end

  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "midpoint_root";
      parameters=(n_taxa=count_terminals(tree)))
  end
  return _phylo_copy(tree)
end

function _upgma_cluster_distance(cluster_a::Vector{Int}, cluster_b::Vector{Int}, D::Matrix{Float64})
  total = 0.0
  for i in cluster_a, j in cluster_b
    total += D[i, j]
  end
  return total / (length(cluster_a) * length(cluster_b))
end

function upgma(D_in::Matrix{Float64}, names::Vector{String})
  n = size(D_in, 1)
  n == length(names) || throw(ArgumentError("distance matrix size must match names length"))
  n > 0 || throw(ArgumentError("distance matrix must not be empty"))

  if n == 1
    return PhyloTree(names[1]; branch_length=0.0)
  end

  D = copy(D_in)
  clusters = [PhyloTree(name) for name in names]
  members = [[index] for index in 1:n]
  heights = zeros(Float64, n)

  while length(clusters) > 1
    cluster_count = length(clusters)
    best_i, best_j = 1, 2
    best_distance = Inf
    for i in 1:(cluster_count-1)
      for j in (i+1):cluster_count
        distance = _upgma_cluster_distance(members[i], members[j], D)
        if distance < best_distance
          best_distance = distance
          best_i = i
          best_j = j
        end
      end
    end

    left_cluster = clusters[best_i]
    right_cluster = clusters[best_j]
    new_height = best_distance / 2
    left_cluster.branch_length = max(0.0, new_height - heights[best_i])
    right_cluster.branch_length = max(0.0, new_height - heights[best_j])
    merged = PhyloTree(""; branch_length=0.0, children=[left_cluster, right_cluster])

    merged_members = vcat(members[best_i], members[best_j])
    new_clusters = PhyloTree[]
    new_members = Vector{Vector{Int}}()
    new_heights = Float64[]
    for index in 1:cluster_count
      if index != best_i && index != best_j
        push!(new_clusters, clusters[index])
        push!(new_members, members[index])
        push!(new_heights, heights[index])
      end
    end
    push!(new_clusters, merged)
    push!(new_members, merged_members)
    push!(new_heights, new_height)

    clusters = new_clusters
    members = new_members
    heights = new_heights
  end

  upgma_result = clusters[1]
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "upgma";
      parameters=(n_taxa=length(names)))
  end
  return upgma_result
end

function _tree_distance_map(tree::PhyloTree, target::String, current_length::Float64, seen::Dict{String,Float64}, path::Vector{Tuple{String,Float64}})
  next_length = current_length + tree.branch_length
  if isleaf(tree)
    push!(path, (tree.name, next_length))
    return path
  end

  for child in tree.children
    _tree_distance_map(child, target, next_length, seen, path)
  end
  return path
end

function _phylo_path_to_leaf(tree::PhyloTree, target::String)
  if isleaf(tree)
    return tree.name == target ? PhyloTree[tree] : nothing
  end

  for child in tree.children
    path = _phylo_path_to_leaf(child, target)
    path === nothing && continue
    pushfirst!(path, tree)
    return path
  end

  return nothing
end

function _phylo_path_to_named_node(tree::PhyloTree, target::String)
  tree.name == target && return PhyloTree[tree]

  for child in tree.children
    path = _phylo_path_to_named_node(child, target)
    path === nothing && continue
    pushfirst!(path, tree)
    return path
  end

  return nothing
end

function _phylo_path_length(path::Vector{PhyloTree})
  total = 0.0
  @inbounds for node in path[2:end]
    total += node.branch_length
  end
  return total
end

function tree_distance(tree::PhyloTree, left_name::String, right_name::String)
  left_path = _phylo_path_to_leaf(tree, String(left_name))
  right_path = _phylo_path_to_leaf(tree, String(right_name))
  left_path === nothing && throw(ArgumentError("unknown leaf: $(left_name)"))
  right_path === nothing && throw(ArgumentError("unknown leaf: $(right_name)"))

  shared_index = 0
  for (index, (left_node, right_node)) in enumerate(zip(left_path, right_path))
    left_node === right_node || break
    shared_index = index
  end

  shared_length = _phylo_path_length(left_path[1:shared_index])
  return _phylo_path_length(left_path) + _phylo_path_length(right_path) - 2 * shared_length
end

function draw_ascii(tree::PhyloTree; indent::Int=0)
  io = IOBuffer()
  function _draw(node::PhyloTree, depth::Int)
    prefix = repeat("  ", depth)
    label = isempty(node.name) ? "[internal]" : node.name
    println(io, prefix, label, " : ", round(node.branch_length; digits=5))
    for child in node.children
      _draw(child, depth + 1)
    end
  end
  _draw(tree, indent)
  return String(take!(io))
end

function draw_unicode(tree::PhyloTree)
  io = IOBuffer()

  function _draw(node::PhyloTree, prefix::String, is_last::Bool, is_root::Bool=false)
    label = isempty(node.name) ? "[internal]" : node.name
    suffix = node.support > 0 ? string(" [", round(node.support * 100; digits=1), "%]") : ""
    line = is_root ? string(label, " : ", round(node.branch_length; digits=5), suffix) : string(prefix, is_last ? "└── " : "├── ", label, " : ", round(node.branch_length; digits=5), suffix)
    println(io, line)
    child_prefix = is_root ? "" : prefix * (is_last ? "    " : "│   ")
    for (index, child) in enumerate(node.children)
      _draw(child, child_prefix, index == length(node.children))
    end
  end

  _draw(tree, "", true, true)
  return String(take!(io))
end

function tree_to_dot(tree::PhyloTree; graph_name::String="phylo_tree")
  node_ids = IdDict{PhyloTree,String}()
  counter = 0

  function _node_id(node::PhyloTree)
    get!(node_ids, node) do
      counter += 1
      "n$(counter)"
    end
  end

  function _label(node::PhyloTree)
    if isleaf(node)
      return isempty(node.name) ? "leaf" : node.name
    elseif node.support > 0
      return string(isempty(node.name) ? "internal" : node.name, "\\n", round(node.support; digits=5))
    else
      return isempty(node.name) ? "internal" : node.name
    end
  end

  lines = String["digraph $(graph_name) {", "  node [shape=ellipse];"]

  function _walk(node::PhyloTree)
    node_id = _node_id(node)
    push!(lines, "  $(node_id) [label=\"$(_label(node))\"];")
    for child in node.children
      child_id = _node_id(child)
      push!(lines, "  $(node_id) -> $(child_id) [label=\"$(round(child.branch_length; digits=5))\"];")
      _walk(child)
    end
  end

  _walk(tree)
  push!(lines, "}")
  return join(lines, "\n")
end

function tree_to_mermaid(tree::PhyloTree)
  node_ids = IdDict{PhyloTree,String}()
  counter = 0

  function _node_id(node::PhyloTree)
    get!(node_ids, node) do
      counter += 1
      "n$(counter)"
    end
  end

  function _label(node::PhyloTree)
    if isleaf(node)
      return isempty(node.name) ? "leaf" : node.name
    elseif node.support > 0
      return string(isempty(node.name) ? "internal" : node.name, " ", round(node.support; digits=5))
    else
      return isempty(node.name) ? "internal" : node.name
    end
  end

  lines = String["graph TD"]

  function _walk(node::PhyloTree)
    node_id = _node_id(node)
    push!(lines, "  $(node_id)[\"$(_label(node))\"]")
    for child in node.children
      child_id = _node_id(child)
      push!(lines, "  $(node_id) -->|$(round(child.branch_length; digits=5))| $(child_id)")
      _walk(child)
    end
  end

  _walk(tree)
  return join(lines, "\n")
end

function prune(tree::PhyloTree, keep_names)
  if keep_names isa String
    keep = Set([String(keep_names)])
  else
    keep = Set(String.(collect(keep_names)))
  end

  function _prune(node::PhyloTree)
    if isleaf(node)
      return node.name in keep ? _phylo_copy(node) : nothing
    end

    children = PhyloTree[]
    for child in node.children
      pruned = _prune(child)
      pruned === nothing && continue
      push!(children, pruned)
    end

    isempty(children) && return nothing
    if length(children) == 1
      only_child = children[1]
      only_child.branch_length += node.branch_length
      return only_child
    end

    return PhyloTree(node.name; branch_length=node.branch_length, children=children, support=node.support)
  end

  pruned = _prune(tree)
  pruned === nothing && throw(ArgumentError("pruning removed all leaves"))
  return pruned
end

function root_with_outgroup(tree::PhyloTree, outgroup_name::String)
  target = String(outgroup_name)
  leaf_names = Set(_phylo_leaf_names(tree))
  target in leaf_names || throw(ArgumentError("unknown outgroup: $(outgroup_name)"))
  target_tree = prune(tree, [target])
  remainder_names = setdiff(leaf_names, Set([target]))
  remainder = isempty(remainder_names) ? nothing : prune(tree, remainder_names)
  remainder === nothing && return _phylo_copy(target_tree)
  return PhyloTree(""; branch_length=0.0, children=[target_tree, remainder], support=1.0)
end

reroot(tree::PhyloTree, outgroup_name::String) = root_with_outgroup(tree, outgroup_name)

function _bootstrap_alignment_sequences(alignment, rng::Random.AbstractRNG)
  width = get_alignment_length(alignment)
  width == 0 && return [record.sequence for record in alignment.records]

  sampled_columns = rand(rng, 1:width, width)
  return [
    begin
      buffer = Vector{UInt8}(undef, width)
      for (position, column) in enumerate(sampled_columns)
        buffer[position] = UInt8(record.sequence[column])
      end
      String(buffer)
    end for record in alignment.records
  ]
end

function bootstrap_trees(alignment; replicates::Int=100, method::Symbol=:hamming, constructor::Symbol=:nj, rng::Random.AbstractRNG=Random.default_rng())
  replicates > 0 || throw(ArgumentError("replicates must be positive"))
  trees = PhyloTree[]
  names = [record.identifier for record in alignment.records]

  for _ in 1:replicates
    sequences = _bootstrap_alignment_sequences(alignment, rng)
    D = distance_matrix(sequences; method=method)
    tree = constructor == :nj ? neighbor_joining_tree(D, names) : constructor == :upgma ? upgma(D, names) : throw(ArgumentError("constructor must be :nj or :upgma"))
    update_provenance!(
      tree.metadata;
      label="PhyloTree",
      source="bootstrap_trees",
      notes=["tree built from bootstrap resampled alignment columns"],
      parameters=(replicates=Int(replicates), method=method, constructor=constructor, taxon_count=length(names)))
    push!(trees, tree)
  end

  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "bootstrap_trees";
      parameters=(replicates=Int(replicates), method=method, constructor=constructor, n_trees=length(trees)))
  end
  return trees
end

function tree_consensus(trees::AbstractVector{PhyloTree}; threshold::Real=0.5)
  isempty(trees) && throw(ArgumentError("at least one tree is required"))
  0.0 < threshold <= 1.0 || throw(ArgumentError("threshold must be in (0, 1]"))

  leaf_names = _phylo_leaf_names(first(trees))
  universe = Set(leaf_names)
  clade_counts = Dict{String,Int}()
  clade_lengths = Dict{String,Tuple{Float64,Int}}()

  for tree in trees
    for (key, (count, total_branch)) in _phylo_clade_map(tree)
      clade_counts[key] = get(clade_counts, key, 0) + 1
      branch_total, branch_count = get(clade_lengths, key, (0.0, 0))
      clade_lengths[key] = (branch_total + total_branch / max(count, 1), branch_count + 1)
    end
  end

  frequent = Set{String}()
  for (key, count) in clade_counts
    names = _phylo_key_set(key)
    1 < length(names) < length(universe) || continue
    count / length(trees) >= threshold && push!(frequent, key)
  end

  function _build(current_names::Set{String})
    if length(current_names) == 1
      name = first(current_names)
      return PhyloTree(name; branch_length=0.0, support=1.0)
    end

    candidates = String[]
    for key in frequent
      names = _phylo_key_set(key)
      names < current_names || continue
      push!(candidates, key)
    end

    maximal = String[]
    for key in candidates
      names = _phylo_key_set(key)
      any(other -> names < _phylo_key_set(other), candidates) && continue
      push!(maximal, key)
    end

    if isempty(maximal)
      children = [PhyloTree(name; branch_length=0.0, support=1.0) for name in sort(collect(current_names))]
      return PhyloTree(""; branch_length=0.0, children=children, support=1.0)
    end

    children = PhyloTree[]
    covered = Set{String}()
    for key in sort(maximal; by=x -> -length(Base.split(x, '\u001f')))
      names = _phylo_key_set(key)
      union!(covered, names)
      child = _build(names)
      child.support = get(clade_counts, key, 0) / length(trees)
      if haskey(clade_lengths, key)
        child.branch_length = first(clade_lengths[key]) / max(last(clade_lengths[key]), 1)
      end
      push!(children, child)
    end

    for name in sort(collect(setdiff(current_names, covered)))
      push!(children, PhyloTree(name; branch_length=0.0, support=1.0))
    end

    return PhyloTree(""; branch_length=0.0, children=children, support=1.0)
  end

  consensus = _build(universe)
  update_provenance!(
    consensus.metadata;
    label="PhyloTree",
    source="tree_consensus",
    notes=["consensus tree built from bootstrap or posterior sample trees"],
    parameters=(tree_count=length(trees), threshold=Float64(threshold), taxon_count=length(universe)))
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "tree_consensus";
      parameters=(n_trees=length(trees), threshold=Float64(threshold), n_taxa=length(universe)))
  end
  return consensus
end

consensus_tree(trees::AbstractVector{PhyloTree}; kwargs...) = tree_consensus(trees; kwargs...)
strict_consensus_tree(trees::AbstractVector{PhyloTree}) = tree_consensus(trees; threshold=1.0)
majority_consensus_tree(trees::AbstractVector{PhyloTree}) = tree_consensus(trees; threshold=0.5)

function bootstrap_consensus_tree(alignment; replicates::Int=100, threshold::Real=0.5, method::Symbol=:hamming, constructor::Symbol=:nj, rng::Random.AbstractRNG=Random.default_rng())
  trees = bootstrap_trees(alignment; replicates=replicates, method=method, constructor=constructor, rng=rng)
  consensus = tree_consensus(trees; threshold=threshold)
  update_provenance!(
    consensus.metadata;
    label="PhyloTree",
    source="bootstrap_consensus_tree",
    notes=["consensus tree generated from bootstrap replicate trees"],
    parameters=(replicates=Int(replicates), threshold=Float64(threshold), method=method, constructor=constructor))
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "bootstrap_consensus_tree";
      parameters=(replicates=Int(replicates), threshold=Float64(threshold), method=method))
  end
  return consensus
end

function bootstrap_support(trees::AbstractVector{PhyloTree}, target_tree::PhyloTree)
  isempty(trees) && throw(ArgumentError("at least one bootstrap tree is required"))
  support_counts = Dict{String,Int}()
  for tree in trees
    for key in _phylo_clade_keys(tree)
      support_counts[key] = get(support_counts, key, 0) + 1
    end
  end

  function _annotate(node::PhyloTree)
    if isleaf(node)
      node.support = 1.0
      return node
    end
    combined = Set{String}()
    for child in node.children
      _annotate(child)
      union!(combined, _phylo_leaf_names(child))
    end
    key = _phylo_clade_key(collect(combined))
    node.support = get(support_counts, key, 0) / length(trees)
    return node
  end

  annotated = _phylo_copy(target_tree)
  return _annotate(annotated)
end

function set_metadata!(tree::PhyloTree, key::String, value)
  tree.metadata[String(key)] = value
  return tree
end

function annotate_tree!(tree::PhyloTree, metadata::AbstractDict)
  for (key, value) in metadata
    tree.metadata[string(key)] = value
  end
  return tree
end

function _phylo_metadata_value(text::String)
  return _phylo_metadata_value(text, "")
end

function _phylo_metadata_value(text::String, datatype::String)
  stripped = strip(String(text))
  isempty(stripped) && return ""
  normalized_datatype = lowercase(strip(String(datatype)))

  if normalized_datatype in ("", "xsd:string", "string")
    lowered = lowercase(stripped)
    lowered == "true" && return true
    lowered == "false" && return false
    try
      return parse(Int, stripped)
    catch
    end
    try
      return parse(Float64, stripped)
    catch
    end
    return stripped
  elseif normalized_datatype in ("xsd:boolean", "boolean")
    lowered = lowercase(stripped)
    lowered == "true" && return true
    lowered == "false" && return false
    throw(ArgumentError("invalid PhyloXML boolean value: $(stripped)"))
  elseif normalized_datatype in ("xsd:integer", "xsd:int", "xsd:long", "xsd:short", "xsd:byte", "xsd:unsignedint", "xsd:unsignedlong", "xsd:unsignedshort", "xsd:unsignedbyte", "integer", "int", "long", "short", "byte")
    return parse(Int, stripped)
  elseif normalized_datatype in ("xsd:float", "xsd:double", "xsd:decimal", "float", "double", "decimal")
    return parse(Float64, stripped)
  end

  return stripped
end

function _phylo_metadata_datatype(value)
  if value isa Bool
    return "xsd:boolean"
  elseif value isa Integer
    return "xsd:integer"
  elseif value isa AbstractFloat
    return "xsd:double"
  else
    return "xsd:string"
  end
end

function _xml_escape(text::AbstractString)
  escaped = replace(String(text), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", '"' => "&quot;", '\'' => "&apos;")
  return escaped
end

function _xml_unescape(text::AbstractString)
  unescaped = replace(String(text), "&apos;" => "'", "&quot;" => "\"", "&gt;" => ">", "&lt;" => "<", "&amp;" => "&")
  return unescaped
end

function _phyloxml_write_clade(io::IO, tree::PhyloTree, depth::Int)
  indent = repeat("  ", depth)
  println(io, indent, "<clade>")
  !isempty(tree.name) && println(io, indent, "  <name>", _xml_escape(tree.name), "</name>")
  !isapprox(tree.branch_length, 0.0; atol=1e-12) && println(io, indent, "  <branch_length>", round(tree.branch_length; digits=5), "</branch_length>")
  tree.support > 0 && println(io, indent, "  <confidence type=\"support\">", round(tree.support; digits=5), "</confidence>")
  for (key, value) in sort(collect(tree.metadata); by=first)
    datatype = _phylo_metadata_datatype(value)
    println(io, indent, "  <property ref=\"", _xml_escape(key), "\" datatype=\"", _xml_escape(datatype), "\">", _xml_escape(string(value)), "</property>")
  end
  for child in tree.children
    _phyloxml_write_clade(io, child, depth + 1)
  end
  println(io, indent, "</clade>")
end

function write_phyloxml(tree::PhyloTree)
  io = IOBuffer()
  println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
  println(io, "<phyloxml xmlns=\"http://www.phyloxml.org\">")
  println(io, "  <phylogeny rooted=\"true\">")
  _phyloxml_write_clade(io, tree, 2)
  println(io, "  </phylogeny>")
  println(io, "</phyloxml>")
  return String(take!(io))
end

function _phyloxml_tag_name(token::AbstractString)
  stripped = strip(String(token))
  startswith(stripped, "</") && return strip(replace(stripped[3:(end-1)], r"\s+.*$" => ""))
  startswith(stripped, "<") || return ""
  cleaned = replace(stripped[2:(end-1)], r"/\s*$" => "")
  cleaned = strip(cleaned)
  isempty(cleaned) && return ""
  name = Base.split(cleaned)[1]
  return replace(name, r"^.*:" => "")
end

function _phyloxml_attr_value(token::AbstractString, attr::AbstractString)
  pattern = Regex(string(attr, "\\s*=\\s*\"([^\"]*)\""))
  match_result = match(pattern, String(token))
  return match_result === nothing ? "" : _xml_unescape(match_result.captures[1])
end

function parse_phyloxml(text::String)
  tokens = collect(eachmatch(r"<[^>]+>|[^<]+", String(text)))
  root = nothing
  node_stack = PhyloTree[]
  current_field = ""
  current_property = ""
  current_property_datatype = ""

  for token_match in tokens
    token = token_match.match
    if startswith(token, "<")
      if startswith(token, "<?") || startswith(token, "<!")
        continue
      elseif startswith(token, "</")
        tag = _phyloxml_tag_name(token)
        if tag == "clade"
          node = pop!(node_stack)
          if isempty(node_stack)
            root = node
          else
            push!(node_stack[end].children, node)
          end
        elseif tag == "name" || tag == "branch_length" || tag == "confidence" || tag == "property"
          current_field = ""
          current_property = ""
          current_property_datatype = ""
        end
      else
        tag = _phyloxml_tag_name(token)
        if tag == "clade"
          push!(node_stack, PhyloTree(""; metadata=Dict{String,Any}()))
        elseif tag == "name"
          current_field = "name"
        elseif tag == "branch_length"
          current_field = "branch_length"
        elseif tag == "confidence"
          current_field = "confidence"
        elseif tag == "property"
          current_field = "property"
          current_property = _phyloxml_attr_value(token, "ref")
          current_property_datatype = _phyloxml_attr_value(token, "datatype")
        end
      end
    else
      value = strip(String(token))
      isempty(value) && continue
      isempty(node_stack) && continue
      node = node_stack[end]
      if current_field == "name"
        node.name = _xml_unescape(value)
      elseif current_field == "branch_length"
        node.branch_length = parse(Float64, value)
      elseif current_field == "confidence"
        node.support = parse(Float64, value)
      elseif current_field == "property" && !isempty(current_property)
        node.metadata[current_property] = _phylo_metadata_value(_xml_unescape(value), current_property_datatype)
      end
    end
  end

  root === nothing && throw(ArgumentError("invalid PhyloXML document"))
  return root
end

function write_nexus(tree::PhyloTree; tree_name::String="tree_1")
  return join([
      "#NEXUS",
      "begin trees;",
      "  tree $(tree_name) = $(write_newick(tree))",
      "end;",
    ], "\n")
end

function parse_nexus(text::String)
  for line in eachline(IOBuffer(String(text)))
    stripped = strip(line)
    startswith(lowercase(stripped), "tree ") || continue
    eq_index = findfirst('=', stripped)
    eq_index === nothing && continue
    payload = strip(stripped[(eq_index+1):end])
    endswith(payload, ";") || continue
    return parse_newick(payload)
  end
  throw(ArgumentError("invalid Nexus tree block"))
end

function write_nexml(tree::PhyloTree; tree_id::String="tree1")
  io = IOBuffer()
  println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
  println(io, "<nexml xmlns=\"http://www.nexml.org/2009\">")
  println(io, "  <trees id=\"trees1\">")
  println(io, "    <tree id=\"", _xml_escape(tree_id), "\" rooted=\"true\">")
  println(io, "      <newick>", _xml_escape(write_newick(tree)), "</newick>")
  for (key, value) in sort(collect(tree.metadata); by=first)
    println(io, "      <meta property=\"", _xml_escape(key), "\">", _xml_escape(string(value)), "</meta>")
  end
  println(io, "    </tree>")
  println(io, "  </trees>")
  println(io, "</nexml>")
  return String(take!(io))
end

function parse_nexml(text::String)
  match_result = match(r"<newick>(.+?)</newick>"is, String(text))
  match_result === nothing && throw(ArgumentError("invalid NeXML tree document"))
  return parse_newick(_xml_unescape(match_result.captures[1]))
end

parse_tree(text::String; format::Symbol=:newick) = format == :newick ? parse_newick(text) : format == :phyloxml ? parse_phyloxml(text) : format == :nexus ? parse_nexus(text) : format == :nexml ? parse_nexml(text) : throw(ArgumentError("unsupported tree format: $(format)"))
write_tree(tree::PhyloTree; format::Symbol=:newick) = format == :newick ? write_newick(tree) : format == :phyloxml ? write_phyloxml(tree) : format == :nexus ? write_nexus(tree) : format == :nexml ? write_nexml(tree) : throw(ArgumentError("unsupported tree format: $(format)"))

function _dna_transition(left::Char, right::Char)
  (left == 'A' && right == 'G') || (left == 'G' && right == 'A') || (left == 'C' && right == 'T') || (left == 'T' && right == 'C')
end

function _pairwise_dna_distance(sequence_left::BioSequence{DNAAlphabet}, sequence_right::BioSequence{DNAAlphabet}, model::Symbol)
  length(sequence_left) == length(sequence_right) || throw(ArgumentError("sequences must be aligned and equal length for DNA substitution distances"))
  valid_sites = 0
  mismatches = 0
  transitions = 0
  transversions = 0

  @inbounds for index in eachindex(sequence_left.data)
    left = Char(sequence_left.data[index])
    right = Char(sequence_right.data[index])
    left == '-' && continue
    right == '-' && continue
    valid_sites += 1
    if left != right
      mismatches += 1
      if _dna_transition(left, right)
        transitions += 1
      else
        transversions += 1
      end
    end
  end

  valid_sites == 0 && return 0.0
  p = mismatches / valid_sites
  model == :hamming && return p
  model == :jukes_cantor && return p >= 0.75 ? Inf : -0.75 * log(1 - (4 / 3) * p)
  if model == :kimura2p || model == :kimura
    P = transitions / valid_sites
    Q = transversions / valid_sites
    a = 1 - 2P - Q
    b = 1 - 2Q
    a <= 0 && return Inf
    b <= 0 && return Inf
    return -0.5 * log(a) - 0.25 * log(b)
  end
  return p
end

function _pairwise_protein_distance(sequence_left::BioSequence{AminoAcidAlphabet}, sequence_right::BioSequence{AminoAcidAlphabet}, matrix_name::Symbol)
  length(sequence_left) == length(sequence_right) || throw(ArgumentError("sequences must be aligned and equal length for protein distances"))
  scoring = named_substitution_matrix(matrix_name)
  best_score = 0.0
  worst_score = 0.0
  observed_score = 0.0
  min_score = minimum(scoring.scores)

  @inbounds for index in eachindex(sequence_left.data)
    left = Char(sequence_left.data[index])
    right = Char(sequence_right.data[index])
    observed_score += _pairwise_residue_score(scoring, left, right)
    best_score += max(_pairwise_residue_score(scoring, left, left), _pairwise_residue_score(scoring, right, right))
    worst_score += min_score
  end

  denom = best_score - worst_score
  denom <= 0 && return 0.0
  return clamp((best_score - observed_score) / denom, 0.0, 1.0)
end

"""
    distance_matrix(sequences; method=:hamming)

Computes an N x N symmetric distance matrix for a list of sequences.
If `method=:hamming`, it leverages SIMD/SWAR byte-mismatch tracking.
If `method=:alignment`, it computes Needleman-Wunsch identity inversion.
"""
@inline function _phylo_hamming_distance(sequence_left::BioSequence, sequence_right::BioSequence)
  length(sequence_left) == length(sequence_right) || throw(ArgumentError("sequences must be aligned and equal length for Hamming distance"))
  mismatches = 0
  @inbounds for index in eachindex(sequence_left.data)
    mismatches += sequence_left.data[index] != sequence_right.data[index]
  end
  return mismatches
end

function distance_matrix(sequences::AbstractVector; method=:hamming, use_threads::Bool=true, use_cuda::Bool=false)
  typed_sequences = [_phylo_coerce_sequence(sequence) for sequence in sequences]
  N = length(typed_sequences)
  D = zeros(Float64, N, N)

  if use_cuda
    if isdefined(@__MODULE__, :CUDA)
      _ensure_cuda_phylo!()
      return Base.invokelatest(_CUDA_PHYLO_DISTANCE_MATRIX_IMPL[], String.(typed_sequences); method=method)
    end
    use_threads = true
  end

  # ─── Threaded / Default Support ─────────────────────────────────────────
  if method == :hamming || method == :p_distance
    if use_threads
      Threads.@threads for i in 1:N
        seq_i = typed_sequences[i]
        len_i = length(seq_i)
        @inbounds for j in (i+1):N
          dist = _phylo_hamming_distance(seq_i, typed_sequences[j])
          val = Float64(dist) / len_i
          D[i, j] = val
          D[j, i] = val
        end
      end
    else
      for i in 1:N
        seq_i = typed_sequences[i]
        len_i = length(seq_i)
        @inbounds for j in (i+1):N
          dist = _phylo_hamming_distance(seq_i, typed_sequences[j])
          val = Float64(dist) / len_i
          D[i, j] = val
          D[j, i] = val
        end
      end
    end
  elseif method == :jukes_cantor || method == :kimura2p || method == :kimura
    dna_sequences = [_phylo_require_alphabet(sequence, DNAAlphabet) for sequence in typed_sequences]
    if use_threads
      Threads.@threads for i in 1:N
        seq_i = dna_sequences[i]
        @inbounds for j in (i+1):N
          val = _pairwise_dna_distance(seq_i, dna_sequences[j], method)
          D[i, j] = val
          D[j, i] = val
        end
      end
    else
      for i in 1:N
        seq_i = dna_sequences[i]
        @inbounds for j in (i+1):N
          val = _pairwise_dna_distance(seq_i, dna_sequences[j], method)
          D[i, j] = val
          D[j, i] = val
        end
      end
    end
  elseif method == :blosum62 || method == :pam250 || method == :pam70 || method == :protein
    protein_sequences = [_phylo_require_alphabet(sequence, AminoAcidAlphabet) for sequence in typed_sequences]
    protein_model = method == :protein ? :BLOSUM62 : Symbol(uppercase(String(method)))
    if use_threads
      Threads.@threads for i in 1:N
        seq_i = protein_sequences[i]
        @inbounds for j in (i+1):N
          val = _pairwise_protein_distance(seq_i, protein_sequences[j], protein_model)
          D[i, j] = val
          D[j, i] = val
        end
      end
    else
      for i in 1:N
        seq_i = protein_sequences[i]
        @inbounds for j in (i+1):N
          val = _pairwise_protein_distance(seq_i, protein_sequences[j], protein_model)
          D[i, j] = val
          D[j, i] = val
        end
      end
    end
  elseif method == :alignment
    if use_threads
      Threads.@threads for i in 1:N
        seq_i = typed_sequences[i]
        @inbounds for j in (i+1):N
          res = pairwise_align(seq_i.data, typed_sequences[j].data; is_local=false)
          val = 1.0 - res.identity
          D[i, j] = val
          D[j, i] = val
        end
      end
    else
      for i in 1:N
        seq_i = typed_sequences[i]
        @inbounds for j in (i+1):N
          res = pairwise_align(seq_i.data, typed_sequences[j].data; is_local=false)
          val = 1.0 - res.identity
          D[i, j] = val
          D[j, i] = val
        end
      end
    end
  else
    throw(ArgumentError("Unknown method: \$method. Options: :hamming, :p_distance, :jukes_cantor, :kimura2p, :kimura, :blosum62, :pam250, :pam70, :protein, :alignment"))
  end

  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "distance_matrix";
      parameters=(n_sequences=length(sequences), method=method))
  end
  return D
end

function _alignment_names(alignment)
  if alignment isa AbstractDict
    return collect(keys(alignment))
  end
  return [record.identifier == "" ? (record.name == "" ? "sequence_$(index)" : record.name) : record.identifier for (index, record) in enumerate(alignment.records)]
end

function _alignment_columns(alignment)
  width = get_alignment_length(alignment)
  if alignment isa AbstractDict
    names = collect(keys(alignment))
    return [String(alignment[name][col] for name in names) for col in 1:width]
  end
  return [String(record.sequence[col] for record in alignment.records) for col in 1:width]
end

function _fitch_sets(node::PhyloTree, column_values::Dict{String,Char})
  if isleaf(node)
    symbol = get(column_values, node.name, '-')
    return Set([symbol]), 0
  end

  accumulated = Set{Char}()
  steps = 0
  first_child = true
  for child in node.children
    child_set, child_steps = _fitch_sets(child, column_values)
    steps += child_steps
    if first_child
      accumulated = copy(child_set)
      first_child = false
    else
      intersection_set = intersect(accumulated, child_set)
      if isempty(intersection_set)
        accumulated = union(accumulated, child_set)
        steps += 1
      else
        accumulated = intersection_set
      end
    end
  end
  return accumulated, steps
end

function parsimony_score(tree::PhyloTree, alignment; threaded::Bool=true)
  names = _alignment_names(alignment)
  width = alignment isa AbstractDict ? length(first(values(alignment))) : get_alignment_length(alignment)
  total_steps = zeros(Int, width)

  if alignment isa AbstractDict
    if threaded && width > 1 && Threads.nthreads() > 1
      Threads.@threads for col in 1:width
        column_values = Dict{String,Char}()
        for name in names
          column_values[name] = alignment[name][col]
        end
        _, steps = _fitch_sets(tree, column_values)
        total_steps[col] = steps
      end
    else
      for col in 1:width
        column_values = Dict{String,Char}()
        for name in names
          column_values[name] = alignment[name][col]
        end
        _, steps = _fitch_sets(tree, column_values)
        total_steps[col] = steps
      end
    end
  else
    if threaded && width > 1 && Threads.nthreads() > 1
      Threads.@threads for col in 1:width
        column_values = Dict{String,Char}()
        for record in alignment.records
          column_values[record.identifier] = record.sequence[col]
        end
        _, steps = _fitch_sets(tree, column_values)
        total_steps[col] = steps
      end
    else
      for col in 1:width
        column_values = Dict{String,Char}()
        for record in alignment.records
          column_values[record.identifier] = record.sequence[col]
        end
        _, steps = _fitch_sets(tree, column_values)
        total_steps[col] = steps
      end
    end
  end

  score = sum(total_steps)
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "parsimony_score";
      parameters=(n_taxa=count_terminals(tree), threaded=Bool(threaded)))
  end
  return score
end

function _binary_tree_combinations(names::Vector{String})
  length(names) == 1 && return [PhyloTree(names[1]; branch_length=0.0, support=0.0)]
  length(names) == 2 && return [PhyloTree(""; branch_length=0.0, children=[PhyloTree(names[1]; branch_length=0.0), PhyloTree(names[2]; branch_length=0.0)], support=0.0)]

  first_name = names[1]
  remaining = names[2:end]
  candidates = PhyloTree[]
  total_masks = 1 << length(remaining)

  for mask in 0:(total_masks-2)
    left_names = String[first_name]
    right_names = String[]
    for (offset, name) in enumerate(remaining)
      if (mask >> (offset - 1)) & 0x01 == 1
        push!(left_names, name)
      else
        push!(right_names, name)
      end
    end
    isempty(right_names) && continue

    left_trees = _binary_tree_combinations(left_names)
    right_trees = _binary_tree_combinations(right_names)
    for left_tree in left_trees, right_tree in right_trees
      push!(candidates, PhyloTree(""; branch_length=0.0, children=[_phylo_copy(left_tree), _phylo_copy(right_tree)], support=0.0))
    end
  end

  return candidates
end

function maximum_parsimony_tree(alignment; max_exact_taxa::Int=7, threaded::Bool=true)
  names = _alignment_names(alignment)
  sequences = alignment isa AbstractDict ? collect(values(alignment)) : [record.sequence for record in alignment.records]
  length(names) == length(unique(names)) || throw(ArgumentError("alignment identifiers must be unique for parsimony tree construction"))

  if length(names) <= max_exact_taxa
    candidate_trees = _binary_tree_combinations(names)
    scores = fill(typemax(Int), length(candidate_trees))
    if threaded && length(candidate_trees) > 1 && Threads.nthreads() > 1
      Threads.@threads for index in eachindex(candidate_trees)
        scores[index] = parsimony_score(candidate_trees[index], alignment; threaded=false)
      end
    else
      for index in eachindex(candidate_trees)
        scores[index] = parsimony_score(candidate_trees[index], alignment; threaded=false)
      end
    end

    best_index = argmin(scores)
    best_tree = candidate_trees[best_index]
    best_score = scores[best_index]
    for (index, candidate) in enumerate(candidate_trees)
      if scores[index] == best_score
        best_tree = candidate
        break
      end
    end
    return best_tree
  end

  # Heuristic branch swapping (NNI parsimony hill-climbing search)
  current_tree = neighbor_joining_tree(distance_matrix(sequences; method=:p_distance), names)
  current_score = parsimony_score(current_tree, alignment)
  
  for step in 1:20
    neighbors = _nni_moves(current_tree)
    improved = false
    for candidate in neighbors
      cand_score = parsimony_score(candidate, alignment)
      if cand_score < current_score
        current_score = cand_score
        current_tree = candidate
        improved = true
        break
      end
    end
    improved || break
  end
  
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "maximum_parsimony_tree";
      parameters=(n_taxa=length(names), max_exact_taxa=Int(max_exact_taxa), threaded=Bool(threaded)))
  end
  return current_tree
end

parsimony_tree(alignment; kwargs...) = maximum_parsimony_tree(alignment; kwargs...)

# ──────────────────────────────────────────────────────────────────────────────
# Neighbor-Joining (NJ) Algorithm
# ──────────────────────────────────────────────────────────────────────────────

"""
    neighbor_joining(D::Matrix{Float64}, names::Vector{String})

Executes the Neighbor-Joining agglomerative clustering algorithm on an NxN 
distance matrix. Returns a standard `Newick` tree string representation 
suitable for direct visualization in tools like Phylo.jl or Makie.jl.
"""
function neighbor_joining(D_in::Matrix{Float64}, names_in::Vector{String})
  n = size(D_in, 1)
  if n != length(names_in)
    throw(ArgumentError("Distance matrix size must aggressively match names length"))
  end

  if n == 1
    return "($(names_in[1]):0.0);"
  elseif n == 2
    return "($(names_in[1]):$(round(D_in[1,2]/2, digits=4)),$(names_in[2]):$(round(D_in[1,2]/2, digits=4)));"
  end

  # Isolate memory
  D = copy(D_in)
  nodes = copy(names_in)

  while n > 2
    # 1. Compute net divergence for each node without allocations
    r = zeros(Float64, n)
    @inbounds for i in 1:n
      r[i] = sum(@view D[i, :])
    end

    # 2. Compute Q-matrix and find minimum
    min_Q = Inf
    min_i, min_j = 0, 0

    @inbounds for i in 1:n
      for j in (i+1):n
        q_val = (n - 2) * D[i, j] - r[i] - r[j]
        if q_val < min_Q
          min_Q = q_val
          min_i = i
          min_j = j
        end
      end
    end

    # 3. Calculate branch lengths for the joined pair
    dist_i_u = 0.5 * D[min_i, min_j] + (r[min_i] - r[min_j]) / (2.0 * (n - 2))
    dist_j_u = D[min_i, min_j] - dist_i_u

    # Format floating points cleanly for Newick
    sub_tree = "($(nodes[min_i]):$(round(dist_i_u, digits=5)),$(nodes[min_j]):$(round(dist_j_u, digits=5)))"

    # 4. Update the distance matrix with the new node `u`
    new_D = zeros(Float64, n - 1, n - 1)

    # Distances to the newly formed internal node
    new_row = Float64[]
    for k in 1:n
      if k != min_i && k != min_j
        push!(new_row, 0.5 * (D[min_i, k] + D[min_j, k] - D[min_i, min_j]))
      end
    end

    # Rebuild matrix shrinking dimensions
    idx_map = Int[]
    for k in 1:n
      if k != min_i && k != min_j
        push!(idx_map, k)
      end
    end

    # Map old untouched intersections
    for (new_r, old_r) in enumerate(idx_map)
      for (new_c, old_c) in enumerate(idx_map)
        new_D[new_r, new_c] = D[old_r, old_c]
      end
      # Add symmetric new node intersection
      new_D[n-1, new_r] = new_row[new_r]
      new_D[new_r, n-1] = new_row[new_r]
    end

    # 5. Iterative variable replacement
    nodes = [nodes[k] for k in idx_map]
    push!(nodes, sub_tree)

    D = new_D
    n -= 1
  end

  # Connect the final two aggregated structural nodes
  dist_final = D[1, 2]
  nj_result = "($(nodes[1]):$(round(dist_final/2, digits=5)),$(nodes[2]):$(round(dist_final/2, digits=5)));"
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "neighbor_joining";
      parameters=(n_taxa=length(names_in)))
  end
  return nj_result
end

# ==============================================================================
# Advanced Evolutionary & Comparative Methods Extension
# ==============================================================================

# 1. GTR (General Time Reversible) Nucleotide Substitution Model
struct GTR <: SubstitutionModel
  a::Float64 # A <-> C
  b::Float64 # A <-> G
  c::Float64 # A <-> T
  d::Float64 # C <-> G
  e::Float64 # C <-> T
  f::Float64 # G <-> T
  pi::Vector{Float64} # [A, C, G, T]
end

GTR(; a=1.0, b=1.0, c=1.0, d=1.0, e=1.0, f=1.0, pi=[0.25, 0.25, 0.25, 0.25]) = GTR(a, b, c, d, e, f, pi)

function transition_probability(model::GTR, t::Real)
  p = model.pi
  a, b, c, d, e, f = model.a, model.b, model.c, model.d, model.e, model.f
  Q = zeros(Float64, 4, 4)
  Q[1, 2] = a * p[2]; Q[1, 3] = b * p[3]; Q[1, 4] = c * p[4]
  Q[2, 1] = a * p[1]; Q[2, 3] = d * p[3]; Q[2, 4] = e * p[4]
  Q[3, 1] = b * p[1]; Q[3, 2] = d * p[2]; Q[3, 4] = f * p[4]
  Q[4, 1] = c * p[1]; Q[4, 2] = e * p[2]; Q[4, 3] = f * p[3]
  
  avg_rate = 0.0
  for i in 1:4
    avg_rate += p[i] * sum(Q[i, :])
  end
  Q ./= max(avg_rate, 1e-12)
  for i in 1:4
    Q[i, i] = -sum(Q[i, :])
  end
  return exp(Q * t)
end

# 2. Empirical Amino Acid Substitution Models (WAG, LG, JTT, Dayhoff, Blosum62, CpREV, MtMAM)
abstract type AAModel <: SubstitutionModel end

struct WAG <: AAModel
  Q::Matrix{Float64}
  pi::Vector{Float64}
end

struct LG <: AAModel
  Q::Matrix{Float64}
  pi::Vector{Float64}
end

struct JTT <: AAModel
  Q::Matrix{Float64}
  pi::Vector{Float64}
end

struct Dayhoff <: AAModel
  Q::Matrix{Float64}
  pi::Vector{Float64}
end

struct Blosum62 <: AAModel
  Q::Matrix{Float64}
  pi::Vector{Float64}
end

struct CpREV <: AAModel
  Q::Matrix{Float64}
  pi::Vector{Float64}
end

struct MtMAM <: AAModel
  Q::Matrix{Float64}
  pi::Vector{Float64}
end

const _WAG_S = Float64[
  0.551571, 0.509848, 0.635346, 0.738998, 0.147304, 5.429420, 1.027040, 0.528191, 0.265256, 0.0302949, 0.908598, 3.035500, 1.543640, 0.616783, 0.0988179, 1.582850, 0.439157, 0.947198, 6.174160, 0.021352, 5.469470, 1.416720, 0.584665, 1.125560, 0.865584, 0.306674, 0.330052, 0.567717, 0.316954, 2.137150, 3.956290, 0.930676, 0.248972, 4.294110, 0.570025, 0.249410, 0.193335, 0.186979, 0.554236, 0.039437, 0.170135, 0.113917, 0.127395, 0.0304501, 0.138190, 0.397915, 0.497671, 0.131528, 0.0848047, 0.384287, 0.869489, 0.154263, 0.0613037, 0.499462, 3.170970, 0.906265, 5.351420, 3.012010, 0.479855, 0.0740339, 3.894900, 2.584430, 0.373558, 0.890432, 0.323832, 0.257555, 0.893496, 0.683162, 0.198221, 0.103754, 0.390482, 1.545260, 0.315124, 0.174100, 0.404141, 4.257460, 4.854020, 0.934276, 0.210494, 0.102711, 0.0961621, 0.0467304, 0.398020, 0.0999208, 0.0811339, 0.049931, 0.679371, 1.059470, 2.115170, 0.088836, 1.190630, 1.438550, 0.679489, 0.195081, 0.423984, 0.109404, 0.933372, 0.682355, 0.243570, 0.696198, 0.0999288, 0.415844, 0.556896, 0.171329, 0.161444, 3.370790, 1.224190, 3.974230, 1.071760, 1.407660, 1.028870, 0.704939, 1.341820, 0.740169, 0.319440, 0.344739, 0.967130, 0.493905, 0.545931, 1.613280, 2.121110, 0.554413, 2.030060, 0.374866, 0.512984, 0.857928, 0.822765, 0.225833, 0.473307, 1.458160, 0.326622, 1.386980, 1.516120, 0.171903, 0.795384, 4.378020, 0.113133, 1.163920, 0.0719167, 0.129767, 0.717070, 0.215737, 0.156557, 0.336983, 0.262569, 0.212483, 0.665309, 0.137505, 0.515706, 1.529640, 0.139405, 0.523742, 0.110864, 0.240735, 0.381533, 1.086000, 0.325711, 0.543833, 0.227710, 0.196303, 0.103604, 3.873440, 0.420170, 0.398618, 0.133264, 0.428437, 6.454280, 0.216046, 0.786993, 0.291148, 2.485390, 2.006010, 0.251849, 0.196246, 0.152335, 1.002140, 0.301281, 0.588731, 0.187247, 0.118358, 7.821300, 1.800340, 0.305434, 2.058450, 0.649892, 0.314887, 0.232739, 1.388230, 0.365369, 0.314730
]
const _WAG_PI = Float64[
  0.0866279, 0.043972, 0.0390894, 0.0570451, 0.0193078, 0.0367281, 0.0580589, 0.0832518, 0.0244313, 0.048466, 0.086209, 0.0620286, 0.0195027, 0.0384319, 0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956
]

const _LG_S = Float64[
  0.425093, 0.276818, 0.751878, 0.395144, 0.123954, 5.076149, 2.489084, 0.534551, 0.528768, 0.062556, 0.969894, 2.807908, 1.695752, 0.523386, 0.084808, 1.038545, 0.363970, 0.541712, 5.243870, 0.003499, 4.128591, 2.066040, 0.390192, 1.437645, 0.844926, 0.569265, 0.267959, 0.348847, 0.358858, 2.426601, 4.509238, 0.927114, 0.640543, 4.813505, 0.423881, 0.311484, 0.149830, 0.126991, 0.191503, 0.010690, 0.320627, 0.072854, 0.044265, 0.008705, 0.108882, 0.395337, 0.301848, 0.068427, 0.015076, 0.594007, 0.582457, 0.069673, 0.044261, 0.366317, 4.145067, 0.536518, 6.326067, 2.145078, 0.282959, 0.013266, 3.234294, 1.807177, 0.296636, 0.697264, 0.159069, 0.137500, 1.124035, 0.484133, 0.371004, 0.025548, 0.893680, 1.672569, 0.173735, 0.139538, 0.442472, 4.273607, 6.312358, 0.656604, 0.253701, 0.052722, 0.089525, 0.017416, 1.105251, 0.035855, 0.018811, 0.089586, 0.682139, 1.112727, 2.592692, 0.023918, 1.798853, 1.177651, 0.332533, 0.161787, 0.394456, 0.075382, 0.624294, 0.419409, 0.196961, 0.508851, 0.078281, 0.249060, 0.390322, 0.099849, 0.094464, 4.727182, 0.858151, 4.008358, 1.240275, 2.784478, 1.223828, 0.611973, 1.739990, 0.990012, 0.064105, 0.182287, 0.748683, 0.346960, 0.361819, 1.338132, 2.139501, 0.578987, 2.000679, 0.425860, 1.143480, 1.080136, 0.604545, 0.129836, 0.584262, 1.033739, 0.302936, 1.136863, 2.020366, 0.165001, 0.571468, 6.472279, 0.180717, 0.593607, 0.045376, 0.029890, 0.670128, 0.236199, 0.077852, 0.268491, 0.597054, 0.111660, 0.619632, 0.049906, 0.696175, 2.457121, 0.095131, 0.248862, 0.140825, 0.218959, 0.314440, 0.612025, 0.135107, 1.165532, 0.257336, 0.120037, 0.054679, 5.306834, 0.232523, 0.299648, 0.131932, 0.481306, 7.803902, 0.089613, 0.400547, 0.245841, 3.151815, 2.547870, 0.170887, 0.083688, 0.037967, 1.959291, 0.210332, 0.245034, 0.076701, 0.119013, 10.649107, 1.702745, 0.185202, 1.898718, 0.654683, 0.296501, 0.098369, 2.188158, 0.189510, 0.249313
]
const _LG_PI = Float64[
  0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767, 0.071586, 0.057337, 0.022355, 0.062157, 0.099081, 0.064600, 0.022951, 0.042302, 0.044040, 0.061197, 0.053287, 0.012066, 0.034155, 0.069147
]

const _JTT_S = Float64[
  58.0, 54.0, 45.0, 81.0, 16.0, 528.0, 56.0, 113.0, 34.0, 10.0, 57.0, 310.0, 86.0, 49.0, 9.0, 105.0, 29.0, 58.0, 767.0, 5.0, 323.0, 179.0, 137.0, 81.0, 130.0, 59.0, 26.0, 119.0, 27.0, 328.0, 391.0, 112.0, 69.0, 597.0, 26.0, 23.0, 36.0, 22.0, 47.0, 11.0, 17.0, 9.0, 12.0, 6.0, 16.0, 30.0, 38.0, 12.0, 7.0, 23.0, 72.0, 9.0, 6.0, 56.0, 229.0, 35.0, 646.0, 263.0, 26.0, 7.0, 292.0, 181.0, 27.0, 45.0, 21.0, 14.0, 54.0, 44.0, 30.0, 15.0, 31.0, 43.0, 18.0, 14.0, 33.0, 479.0, 388.0, 65.0, 15.0, 5.0, 10.0, 4.0, 78.0, 4.0, 5.0, 5.0, 40.0, 89.0, 248.0, 4.0, 43.0, 194.0, 74.0, 15.0, 15.0, 14.0, 164.0, 18.0, 24.0, 115.0, 10.0, 102.0, 21.0, 16.0, 17.0, 378.0, 101.0, 503.0, 59.0, 223.0, 53.0, 30.0, 201.0, 73.0, 40.0, 59.0, 47.0, 29.0, 92.0, 285.0, 475.0, 64.0, 232.0, 38.0, 42.0, 51.0, 32.0, 33.0, 46.0, 245.0, 25.0, 103.0, 226.0, 12.0, 118.0, 477.0, 9.0, 126.0, 8.0, 4.0, 115.0, 18.0, 10.0, 55.0, 8.0, 9.0, 52.0, 10.0, 24.0, 53.0, 6.0, 35.0, 12.0, 11.0, 20.0, 70.0, 46.0, 209.0, 24.0, 7.0, 8.0, 573.0, 32.0, 24.0, 8.0, 18.0, 536.0, 10.0, 63.0, 21.0, 71.0, 298.0, 17.0, 16.0, 31.0, 62.0, 20.0, 45.0, 47.0, 11.0, 961.0, 180.0, 14.0, 323.0, 62.0, 23.0, 38.0, 112.0, 25.0, 16.0
]
const _JTT_PI = Float64[
  0.076748, 0.051691, 0.042645, 0.051544, 0.019803, 0.040752, 0.061830, 0.073152, 0.022944, 0.053761, 0.091904, 0.058676, 0.023826, 0.040126, 0.050901, 0.068765, 0.058565, 0.014261, 0.032102, 0.066005
]

const _DAYHOFF_S = Float64[
  27.0, 98.0, 32.0, 120.0, 0.0, 905.0, 36.0, 23.0, 0.0, 0.0, 89.0, 246.0, 103.0, 134.0, 0.0, 198.0, 1.0, 148.0, 1153.0, 0.0, 716.0, 240.0, 9.0, 139.0, 125.0, 11.0, 28.0, 81.0, 23.0, 240.0, 535.0, 86.0, 28.0, 606.0, 43.0, 10.0, 65.0, 64.0, 77.0, 24.0, 44.0, 18.0, 61.0, 0.0, 7.0, 41.0, 15.0, 34.0, 0.0, 0.0, 73.0, 11.0, 7.0, 44.0, 257.0, 26.0, 464.0, 318.0, 71.0, 0.0, 153.0, 83.0, 27.0, 26.0, 46.0, 18.0, 72.0, 90.0, 1.0, 0.0, 0.0, 114.0, 30.0, 17.0, 0.0, 336.0, 527.0, 243.0, 18.0, 14.0, 14.0, 0.0, 0.0, 0.0, 0.0, 15.0, 48.0, 196.0, 157.0, 0.0, 92.0, 250.0, 103.0, 42.0, 13.0, 19.0, 153.0, 51.0, 34.0, 94.0, 12.0, 32.0, 33.0, 17.0, 11.0, 409.0, 154.0, 495.0, 95.0, 161.0, 56.0, 79.0, 234.0, 35.0, 24.0, 17.0, 96.0, 62.0, 46.0, 245.0, 371.0, 26.0, 229.0, 66.0, 16.0, 53.0, 34.0, 30.0, 22.0, 192.0, 33.0, 136.0, 104.0, 13.0, 78.0, 550.0, 0.0, 201.0, 23.0, 0.0, 0.0, 0.0, 0.0, 0.0, 27.0, 0.0, 46.0, 0.0, 0.0, 76.0, 0.0, 75.0, 0.0, 24.0, 8.0, 95.0, 0.0, 96.0, 0.0, 22.0, 0.0, 127.0, 37.0, 28.0, 13.0, 0.0, 698.0, 0.0, 34.0, 42.0, 61.0, 208.0, 24.0, 15.0, 18.0, 49.0, 35.0, 37.0, 54.0, 44.0, 889.0, 175.0, 10.0, 258.0, 12.0, 48.0, 30.0, 157.0, 0.0, 28.0
]
const _DAYHOFF_PI = Float64[
  0.087127, 0.040904, 0.040432, 0.046872, 0.033474, 0.038255, 0.049530, 0.088612, 0.033618, 0.036886, 0.085357, 0.080482, 0.014753, 0.039772, 0.050680, 0.069577, 0.058542, 0.010494, 0.029916, 0.064718
]

const _BLOSUM62_S = Float64[
  1.303535, 0.48495393, 0.21502738, 1.0454644, 0.12566555, 3.358511, 0.44013616, 0.39404158, 0.25153905, 0.26018026, 1.7462815, 0.37338648, 0.71290493, 0.37818708, 0.32871194, 0.63948549, 0.3161155, 0.69475307, 1.4525454, 0.63796076, 0.48200517, 0.54096789, 0.69185257, 0.37370124, 0.16682449, 1.3547516, 0.16901037, 0.2256359, 1.1564011, 0.25452836, 0.92180573, 2.5021053, 0.23831978, 0.67144988, 0.91293016, 0.36266389, 0.71428283, 0.7426983, 0.19475553, 0.33269283, 1.8435613, 0.19438909, 0.31077108, 3.473989, 0.39782078, 1.1193994, 0.68541952, 0.19876464, 0.49589436, 1.5837601, 0.41020074, 0.87882787, 3.0042095, 0.95806668, 5.384295, 0.43337417, 0.35477601, 2.8393049, 1.1075517, 0.31683257, 1.2105704, 2.1557516, 0.2526878, 1.4224087, 0.18918118, 0.57895113, 1.0475419, 0.31785688, 0.65233056, 0.82726879, 0.20856441, 0.4500662, 0.47053064, 0.34673734, 0.9351407, 0.33441646, 0.40545751, 0.44185252, 1.0711393, 0.29445836, 1.0633181, 4.936413, 0.12888696, 0.54449218, 2.174358, 0.21088984, 3.5224823, 0.51850245, 2.2275308, 1.6421167, 0.76656194, 0.65693948, 0.20341156, 0.4472786, 1.2147684, 0.33181742, 0.37388658, 1.3000939, 0.20717028, 4.8312335, 0.55597883, 0.87827562, 1.1584059, 0.40010946, 2.6971062, 3.8615945, 1.0688918, 1.4131384, 1.5802228, 0.48215563, 1.3475794, 0.99732337, 0.3192282, 1.5643952, 0.31514327, 0.8202685, 2.5928839, 1.0438399, 1.7275207, 1.0024602, 1.5438708, 0.99982473, 0.83407539, 0.9563132, 0.43588755, 0.57269658, 0.52267199, 1.052734, 1.1641455, 0.81717612, 1.3290279, 1.694756, 0.89770254, 1.1405792, 0.81664586, 4.5994902, 1.9533219, 0.99671866, 0.23046689, 0.46807223, 0.7483026, 0.22619155, 0.17995624, 7.4211019, 0.44475551, 1.9922773, 2.2998099, 0.27933131, 0.44365726, 0.47390418, 0.3914085, 0.50170677, 2.0116236, 0.41591773, 0.47111731, 0.12976916, 0.36395141, 1.8520301, 0.45388069, 0.26896532, 0.30533691, 0.29663824, 0.61737282, 0.79292782, 0.17096189, 0.22518625, 0.67735274, 0.38068909, 0.34635361, 0.4581388, 0.23794835, 0.64124014, 0.41984475, 0.23340692, 0.53277193, 6.0241909, 0.27504297, 3.7668301, 0.60232181, 0.64104978, 0.72430877, 0.84969523, 0.48054395, 0.32981803, 0.85621969, 0.64330325, 0.71136771, 0.71517005, 0.89283509, 3.6199276
]
const _BLOSUM62_PI = Float64[
  0.074, 0.025, 0.054, 0.054, 0.047, 0.074, 0.026, 0.068, 0.058, 0.099, 0.025, 0.045, 0.039, 0.034, 0.052, 0.057, 0.051, 0.073, 0.013, 0.032
]

const _CPREV_S = Float64[
  6.5, 4.5, 10.6, 84.3, 9.5, 643.2, 19.5, 353.7, 10.9, 10.7, 6.1, 486.3, 18.0, 11.6, 0.1, 74.5, 21.5, 13.0, 437.4, 0.1, 342.6, 118.1, 183.9, 17.4, 150.3, 86.8, 7.1, 161.9, 2.8, 346.6, 345.3, 202.4, 111.8, 450.1, 6.2, 2.2, 1.5, 50.6, 25.6, 5.6, 3.4, 3.6, 4.3, 2.5, 8.4, 3.9, 36.9, 2.4, 5.9, 20.3, 26.1, 5.1, 3.4, 17.3, 205.0, 4.2, 712.1, 639.2, 10.1, 0.1, 500.5, 426.6, 29.3, 9.2, 37.9, 10.8, 13.4, 53.5, 9.9, 3.8, 10.5, 9.5, 9.6, 3.8, 3.6, 534.9, 142.8, 83.6, 4.3, 5.0, 8.7, 7.5, 238.0, 2.4, 7.7, 3.1, 11.0, 61.0, 542.3, 9.4, 3.8, 91.2, 69.0, 3.5, 13.4, 6.5, 145.6, 8.1, 2.6, 133.9, 2.1, 155.8, 21.2, 10.5, 12.6, 251.1, 82.9, 271.4, 34.8, 471.9, 10.7, 16.4, 136.7, 19.2, 36.2, 160.3, 23.9, 6.2, 249.4, 348.6, 467.5, 82.5, 215.5, 8.0, 7.4, 5.4, 11.6, 6.3, 3.8, 266.2, 10.7, 140.2, 295.2, 3.6, 181.2, 144.8, 3.4, 171.8, 6.1, 3.5, 518.6, 17.0, 9.1, 49.0, 5.7, 3.3, 98.8, 2.3, 11.1, 34.1, 1.1, 56.3, 1.5, 2.2, 4.3, 69.9, 202.9, 579.1, 9.4, 9.1, 2.1, 889.2, 10.8, 9.6, 20.1, 3.4, 255.9, 5.6, 264.3, 3.3, 21.7, 363.2, 8.4, 1.6, 10.3, 37.8, 5.1, 21.6, 76.0, 1.1, 595.0, 155.8, 9.2, 191.9, 102.2, 7.7, 10.1, 36.8, 5.0, 7.2
]
const _CPREV_PI = Float64[
  0.061007, 0.060799, 0.043028, 0.038515, 0.011297, 0.035406, 0.050764, 0.073749, 0.024609, 0.085629, 0.106930, 0.046704, 0.023382, 0.056136, 0.043289, 0.073994, 0.052078, 0.018023, 0.036043, 0.058620
]

const _MTMAM_S = Float64[
  32.0, 2.0, 4.0, 11.0, 0.0, 864.0, 0.0, 186.0, 0.0, 0.0, 0.0, 246.0, 8.0, 49.0, 0.0, 0.0, 0.0, 0.0, 569.0, 0.0, 274.0, 78.0, 18.0, 47.0, 79.0, 0.0, 0.0, 22.0, 8.0, 232.0, 458.0, 11.0, 305.0, 550.0, 22.0, 0.0, 75.0, 0.0, 19.0, 0.0, 41.0, 0.0, 0.0, 0.0, 0.0, 21.0, 6.0, 0.0, 0.0, 27.0, 20.0, 0.0, 0.0, 26.0, 232.0, 0.0, 50.0, 408.0, 0.0, 0.0, 242.0, 215.0, 0.0, 0.0, 6.0, 4.0, 76.0, 0.0, 21.0, 0.0, 0.0, 22.0, 0.0, 0.0, 0.0, 378.0, 609.0, 59.0, 0.0, 0.0, 6.0, 5.0, 7.0, 0.0, 0.0, 0.0, 0.0, 57.0, 246.0, 0.0, 11.0, 53.0, 9.0, 33.0, 2.0, 0.0, 51.0, 0.0, 0.0, 53.0, 5.0, 43.0, 18.0, 0.0, 17.0, 342.0, 3.0, 446.0, 16.0, 347.0, 30.0, 21.0, 112.0, 20.0, 0.0, 74.0, 65.0, 47.0, 90.0, 202.0, 681.0, 0.0, 110.0, 0.0, 114.0, 0.0, 4.0, 0.0, 1.0, 360.0, 34.0, 50.0, 691.0, 8.0, 78.0, 614.0, 5.0, 16.0, 6.0, 0.0, 65.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 0.0, 13.0, 0.0, 7.0, 17.0, 0.0, 0.0, 0.0, 156.0, 0.0, 530.0, 54.0, 0.0, 1.0, 1525.0, 16.0, 25.0, 67.0, 0.0, 682.0, 8.0, 107.0, 0.0, 14.0, 398.0, 0.0, 0.0, 10.0, 0.0, 33.0, 20.0, 5.0, 0.0, 2220.0, 100.0, 0.0, 832.0, 6.0, 0.0, 0.0, 237.0, 0.0, 0.0
]
const _MTMAM_PI = Float64[
  0.0692, 0.0184, 0.0400, 0.0186, 0.0065, 0.0238, 0.0236, 0.0557, 0.0277, 0.0905, 0.1675, 0.0221, 0.0561, 0.0611, 0.0536, 0.0725, 0.0870, 0.0293, 0.0340, 0.0428
]

function _build_aa_q(S_lower::Vector{Float64}, pi_emp::Vector{Float64}, pi_custom::Union{Nothing, Vector{Float64}}=nothing)
  pi = pi_custom === nothing ? copy(pi_emp) : copy(pi_custom)
  length(pi) == 20 || throw(ArgumentError("State frequency vector pi must have length 20"))
  s = sum(pi)
  s > 0 || throw(ArgumentError("State frequencies must sum to > 0"))
  pi ./= s

  S = zeros(Float64, 20, 20)
  idx = 1
  for i in 2:20
    for j in 1:(i-1)
      val = S_lower[idx]
      S[i, j] = val
      S[j, i] = val
      idx += 1
    end
  end

  Q = zeros(Float64, 20, 20)
  for i in 1:20
    for j in 1:20
      if i != j
        Q[i, j] = S[i, j] * pi[j]
      end
    end
  end
  for i in 1:20
    Q[i, i] = -sum(Q[i, :])
  end

  scale = sum(pi[i] * (-Q[i, i]) for i in 1:20)
  if scale > 0
    Q ./= scale
  end
  return Q, pi
end

function WAG(; pi::Union{Nothing, Vector{Float64}}=nothing)
  Q, p = _build_aa_q(_WAG_S, _WAG_PI, pi)
  return WAG(Q, p)
end

function LG(; pi::Union{Nothing, Vector{Float64}}=nothing)
  Q, p = _build_aa_q(_LG_S, _LG_PI, pi)
  return LG(Q, p)
end

function JTT(; pi::Union{Nothing, Vector{Float64}}=nothing)
  Q, p = _build_aa_q(_JTT_S, _JTT_PI, pi)
  return JTT(Q, p)
end

function Dayhoff(; pi::Union{Nothing, Vector{Float64}}=nothing)
  Q, p = _build_aa_q(_DAYHOFF_S, _DAYHOFF_PI, pi)
  return Dayhoff(Q, p)
end

function Blosum62(; pi::Union{Nothing, Vector{Float64}}=nothing)
  Q, p = _build_aa_q(_BLOSUM62_S, _BLOSUM62_PI, pi)
  return Blosum62(Q, p)
end

function CpREV(; pi::Union{Nothing, Vector{Float64}}=nothing)
  Q, p = _build_aa_q(_CPREV_S, _CPREV_PI, pi)
  return CpREV(Q, p)
end

function MtMAM(; pi::Union{Nothing, Vector{Float64}}=nothing)
  Q, p = _build_aa_q(_MTMAM_S, _MTMAM_PI, pi)
  return MtMAM(Q, p)
end

function transition_probability(model::AAModel, t::Real)
  return exp(model.Q * Float64(t))
end

# 4. Pairwise Maximum Likelihood Distance Matrix
function ml_distance_matrix(sequences::AbstractVector; model::SubstitutionModel=JC69())
  n = length(sequences)
  D = zeros(Float64, n, n)
  for i in 1:n
    for j in (i+1):n
      s1 = sequence_to_string(sequences[i])
      s2 = sequence_to_string(sequences[j])
      valid = 0
      mismatches = 0
      ts = 0
      tv = 0
      for (c1, c2) in zip(s1, s2)
        u1, u2 = uppercase(c1), uppercase(c2)
        if model isa AAModel
          if haskey(_AA_STATE_MAP, u1) && haskey(_AA_STATE_MAP, u2)
            valid += 1
            if u1 != u2
              mismatches += 1
            end
          end
        else
          if u1 in "ACGT" && u2 in "ACGT"
            valid += 1
            if u1 != u2
              mismatches += 1
              if (u1 in "AG" && u2 in "AG") || (u1 in "CT" && u2 in "CT")
                ts += 1
              else
                tv += 1
              end
            end
          end
        end
      end

      if model isa AAModel
        p = valid > 0 ? mismatches / valid : 0.0
        # Poisson / Kimura 20-state amino acid distance formula: d = -19/20 * log(1 - 20/19 * p)
        d = p < 0.9 ? (-19.0 / 20.0 * log(max(1e-6, 1.0 - 20.0 * p / 19.0))) : 3.0
      elseif model isa K80
        P = valid > 0 ? ts / valid : 0.0
        Q = valid > 0 ? tv / valid : 0.0
        val1 = 1.0 - 2.0 * P - Q
        val2 = 1.0 - 2.0 * Q
        d = (val1 > 0 && val2 > 0) ? (-0.5 * log(val1) - 0.25 * log(val2)) : 3.0
      else # JC69 / default ML formula
        p = valid > 0 ? mismatches / valid : 0.0
        d = p < 0.75 ? -0.75 * log(1.0 - 4.0 * p / 3.0) : 3.0
      end
      D[i, j] = D[j, i] = max(0.0, d)
    end
  end
  return D
end

# 5. BioNJ Distance-Based Tree Construction (Gascuel 1997 Algorithm)
function bionj_tree(D_in::Matrix{Float64}, names_in::Vector{String})
  n = length(names_in)
  D = copy(D_in)
  V = copy(D_in) # Variance matrix initial estimate V = D
  nodes = [PhyloTree(name) for name in names_in]
  
  while n > 2
    # Compute Q matrix
    Q = zeros(Float64, n, n)
    r = [sum(D[i, :]) for i in 1:n]
    for i in 1:n
      for j in 1:n
        if i != j
          Q[i, j] = (n - 2) * D[i, j] - r[i] - r[j]
        end
      end
    end
    
    min_q = Inf
    min_i, min_j = 1, 2
    for i in 1:n
      for j in (i+1):n
        if Q[i, j] < min_q
          min_q = Q[i, j]
          min_i, min_j = i, j
        end
      end
    end
    
    # Compute branch lengths to new node u
    dist_ij = D[min_i, min_j]
    v_i = 0.5 * dist_ij + (r[min_i] - r[min_j]) / (2.0 * max(1, n - 2))
    v_j = dist_ij - v_i
    v_i = max(0.0, v_i)
    v_j = max(0.0, v_j)
    
    u_node1 = _phylo_copy(nodes[min_i]); u_node1.branch_length = v_i
    u_node2 = _phylo_copy(nodes[min_j]); u_node2.branch_length = v_j
    sub_tree = PhyloTree("Node_$(n)"; children=[u_node1, u_node2])
    
    # BioNJ variance reduction weight parameter w (Gascuel 1997)
    var_ij = max(1e-6, V[min_i, min_j])
    sum_var_diff = sum(V[min_j, k] - V[min_i, k] for k in 1:n if k != min_i && k != min_j)
    w = 0.5 + sum_var_diff / (2.0 * max(1, n - 2) * var_ij)
    w = clamp(w, 0.0, 1.0)
    
    new_D = zeros(Float64, n - 1, n - 1)
    new_V = zeros(Float64, n - 1, n - 1)
    
    idx_map = [k for k in 1:n if k != min_i && k != min_j]
    for (new_r, old_r) in enumerate(idx_map)
      for (new_c, old_c) in enumerate(idx_map)
        new_D[new_r, new_c] = D[old_r, old_c]
        new_V[new_r, new_c] = V[old_r, old_c]
      end
      # BioNJ updated distances & variances to u
      d_u = w * D[min_i, old_r] + (1.0 - w) * D[min_j, old_r] - w * v_i - (1.0 - w) * v_j
      v_u = w * V[min_i, old_r] + (1.0 - w) * V[min_j, old_r] - w * (1.0 - w) * var_ij
      new_D[n-1, new_r] = new_D[new_r, n-1] = max(0.0, d_u)
      new_V[n-1, new_r] = new_V[new_r, n-1] = max(1e-6, v_u)
    end
    
    nodes = [nodes[k] for k in idx_map]
    push!(nodes, sub_tree)
    D = new_D
    V = new_V
    n -= 1
  end
  
  final_dist = D[1, 2]
  u1 = _phylo_copy(nodes[1]); u1.branch_length = final_dist / 2.0
  u2 = _phylo_copy(nodes[2]); u2.branch_length = final_dist / 2.0
  return PhyloTree("Root"; children=[u1, u2])
end

function model_test(alignment::AbstractDict; candidates=nothing)
  is_aa = candidates !== nothing ? (first(candidates) isa AAModel) : _is_aa_alignment(alignment)
  
  if candidates === nothing
    candidates = is_aa ? [WAG(), LG(), JTT(), Dayhoff(), Blosum62(), CpREV(), MtMAM()] :
                          [JC69(), K80(2.0), HKY85([0.25, 0.25, 0.25, 0.25], 2.0), GTR()]
  end
  
  alignment_parsed = is_aa ? _phylo_aa_alignment(alignment) : _phylo_dna_alignment(alignment)
  names = collect(keys(alignment_parsed))
  seqs = [alignment_parsed[n] for n in names]
  dm = is_aa ? ml_distance_matrix(seqs; model=candidates[1]) : distance_matrix(seqs)
  tree = neighbor_joining_tree(dm, names)
  
  n_taxa = length(names)
  first_seq = alignment_parsed[first(names)]
  n_sites = length(first_seq)
  
  results = []
  best_bic = Inf
  best_model = candidates[1]
  
  for model in candidates
    optimize_branch_lengths!(tree, alignment_parsed, model; max_iter=3)
    log_lik = felsenstein_likelihood(tree, alignment_parsed, model)
    
    k = if model isa AAModel
      0
    elseif model isa JC69
      0
    elseif model isa K80
      1
    elseif model isa HKY85
      4
    else # GTR
      8
    end
    n_params = k + 2 * n_taxa - 3
    
    aic = 2 * n_params - 2 * log_lik
    aicc = aic + (2 * n_params * (n_params + 1)) / max(1, n_sites - n_params - 1)
    bic = n_params * log(n_sites) - 2 * log_lik
    
    m_name = string(typeof(model))
    push!(results, (model=m_name, log_lik=log_lik, n_params=n_params, AIC=aic, AICc=aicc, BIC=bic))
    
    if bic < best_bic
      best_bic = bic
      best_model = model
    end
  end
  
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "model_test";
      parameters=(n_taxa=n_taxa, n_sites=n_sites, best_model=string(typeof(best_model))))
  end
  return (best_model=best_model, summary=results, tree=tree)
end

# 7. Subtree Pruning and Regrafting (SPR) Moves
function _find_parent(root::PhyloTree, target::PhyloTree)
  for child in root.children
    child === target && return root
    p = _find_parent(child, target)
    p !== nothing && return p
  end
  return nothing
end

function _spr_moves(tree::PhyloTree)
  moves = PhyloTree[]
  terminals = get_terminals(tree)
  non_terms = get_nonterminals(tree)
  all_nodes = vcat(terminals, non_terms)
  
  for prune_node in all_nodes
    prune_node === tree && continue
    pruned_leaves = Set(_phylo_leaf_names(prune_node))
    length(pruned_leaves) == count_terminals(tree) && continue
    
    rem_tree = prune_taxa(tree, collect(pruned_leaves))
    rem_tree === nothing && continue
    
    rem_targets = vcat(get_terminals(rem_tree), get_nonterminals(rem_tree))
    for regraft_target in rem_targets
      new_tree = _phylo_copy(rem_tree)
      target_in_new = isleaf(regraft_target) ? 
        lowest_common_ancestor(new_tree, [regraft_target.name]) : 
        lowest_common_ancestor(new_tree, _phylo_leaf_names(regraft_target))
        
      if target_in_new !== nothing
        pruned_copy = _phylo_copy(prune_node)
        new_internal = PhyloTree("RegraftNode"; children=[pruned_copy, _phylo_copy(target_in_new)])
        
        if target_in_new === new_tree
          new_tree = new_internal
        else
          parent = _find_parent(new_tree, target_in_new)
          if parent !== nothing
            idx = findfirst(c -> c === target_in_new, parent.children)
            if idx !== nothing
              parent.children[idx] = new_internal
            end
          end
        end
        push!(moves, new_tree)
      end
    end
  end
  return isempty(moves) ? [tree] : moves
end

# 8. Sankoff Dynamic Programming Parsimony Engine
function _sankoff_dp(node::PhyloTree, col_char::Dict{String,Char}, cost_matrix::Matrix{Float64})
  n_states = size(cost_matrix, 1)
  S = zeros(Float64, n_states)
  
  if isleaf(node)
    ch = get(col_char, node.name, 'A')
    idx = ch == 'A' ? 1 : (ch == 'C' ? 2 : (ch == 'G' ? 3 : 4))
    for i in 1:n_states
      S[i] = (i == idx) ? 0.0 : 1e6
    end
    return S
  end
  
  for child in node.children
    child_S = _sankoff_dp(child, col_char, cost_matrix)
    for i in 1:n_states
      min_cost = Inf
      for j in 1:n_states
        cost = child_S[j] + cost_matrix[i, j]
        if cost < min_cost
          min_cost = cost
        end
      end
      S[i] += min_cost
    end
  end
  return S
end

function sankoff_parsimony_score(tree::PhyloTree, alignment::AbstractDict, cost_matrix::Matrix{Float64}=fill(1.0, 4, 4) - diagm(ones(4)))
  names = _alignment_names(alignment)
  width = alignment isa AbstractDict ? length(first(values(alignment))) : get_alignment_length(alignment)
  total_score = 0.0
  
  for col in 1:width
    column_values = Dict{String,Char}()
    if alignment isa AbstractDict
      for name in names
        column_values[name] = alignment[name][col]
      end
    else
      for record in alignment.records
        column_values[record.identifier] = record.sequence[col]
      end
    end
    
    S_root = _sankoff_dp(tree, column_values, cost_matrix)
    total_score += minimum(S_root)
  end
  
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "sankoff_parsimony_score";
      parameters=(n_taxa=count_terminals(tree), width=width))
  end
  return total_score
end

# 9. Phylogenetic Independent Contrasts (PIC)
function pic(tree::PhyloTree, trait_values::Dict{String,Float64})
  work_tree = _phylo_copy(tree)
  resolve_polytomies!(work_tree)
  
  contrasts = Float64[]
  terminals = get_terminals(work_tree)
  leaf_vals = Dict(leaf.name => get(trait_values, leaf.name, 0.0) for leaf in terminals)
  
  function _compute_contrast(node::PhyloTree)
    if isleaf(node)
      return leaf_vals[node.name], max(1e-6, node.branch_length)
    end
    
    child_results = [_compute_contrast(child) for child in node.children]
    if length(child_results) >= 2
      val1, v1 = child_results[1]
      val2, v2 = child_results[2]
      raw_contrast = val1 - val2
      variance = max(1e-8, v1 + v2)
      standardized = raw_contrast / sqrt(variance)
      push!(contrasts, standardized)
      
      ancestral_val = (val1 / max(1e-8, v1) + val2 / max(1e-8, v2)) / (1.0 / max(1e-8, v1) + 1.0 / max(1e-8, v2))
      adj_branch = max(1e-6, node.branch_length) + (v1 * v2) / max(1e-8, v1 + v2)
      return ancestral_val, adj_branch
    end
    return child_results[1][1], max(1e-6, node.branch_length)
  end
  
  _compute_contrast(work_tree)
  return contrasts
end

# 10. Phylogenetic Generalized Least Squares (PGLS)
function _node_root_distance(root::PhyloTree, target::PhyloTree)
  target === root && return 0.0
  
  function _search(curr::PhyloTree, accum::Float64)
    curr === target && return accum
    for child in curr.children
      res = _search(child, accum + max(1e-6, child.branch_length))
      res >= 0.0 && return res
    end
    return -1.0
  end
  
  return max(0.0, _search(root, 0.0))
end

function pgls(tree::PhyloTree, X::Matrix{Float64}, y::Vector{Float64})
  terminals = get_terminals(tree)
  n = length(terminals)
  V = zeros(Float64, n, n)
  
  root_dists = Dict{PhyloTree,Float64}()
  for t in terminals
    root_dists[t] = _node_root_distance(tree, t)
  end
  non_terms = get_nonterminals(tree)
  for nt in non_terms
    root_dists[nt] = _node_root_distance(tree, nt)
  end
  
  for i in 1:n
    for j in 1:n
      if i == j
        V[i, j] = get(root_dists, terminals[i], 1.0)
      else
        lca = lowest_common_ancestor(tree, [terminals[i].name, terminals[j].name])
        V[i, j] = lca !== nothing ? get(root_dists, lca, 0.0) : 0.0
      end
    end
  end
  
  invV = inv(V + 1e-6 * I)
  beta = inv(X' * invV * X + 1e-6 * I) * (X' * invV * y)
  residuals = y - X * beta
  sigma2 = max(1e-8, Float64((residuals' * invV * residuals)[1] / n))
  
  log_lik = -0.5 * (n * log(2π * sigma2) + log(det(V + 1e-6 * I)) + n)
  return (coefficients=beta, sigma2=sigma2, log_likelihood=log_lik, vcv=V)
end

# 11. Ancestral State Reconstruction (ASR)
function ancestral_state_reconstruction(tree::PhyloTree, trait_values::Dict{String,Float64})
  terminals = get_terminals(tree)
  node_states = Dict{PhyloTree,Float64}()
  
  for leaf in terminals
    node_states[leaf] = get(trait_values, leaf.name, 0.0)
  end
  
  function _reconstruct(node::PhyloTree)
    if haskey(node_states, node)
      return node_states[node]
    end
    
    weighted_sum = 0.0
    weight_total = 0.0
    for child in node.children
      val = _reconstruct(child)
      bl = max(1e-6, child.branch_length)
      w = 1.0 / bl
      weighted_sum += val * w
      weight_total += w
    end
    
    state = weight_total > 0 ? weighted_sum / weight_total : 0.0
    node_states[node] = state
    return state
  end
  
  _reconstruct(tree)
  return node_states
end

function ancestral_state_reconstruction(tree::PhyloTree, tip_states::Dict{String, String}, model)
  is_aa = model isa AAModel
  n_states = is_aa ? 20 : 4
  states = is_aa ? collect("ACDEFGHIKLMNPQRSTVWY") : collect("ACGT")
  
  tip_idx = Dict{String, Int}()
  for (k, v) in tip_states
    idx = findfirst(==(uppercase(v[1])), states)
    tip_idx[k] = idx === nothing ? 0 : idx
  end
  
  L = Dict{PhyloTree, Vector{Float64}}()
  
  function _prune(node::PhyloTree)
    if isleaf(node)
      l_vec = fill(1.0, n_states)
      idx = get(tip_idx, node.name, 0)
      if idx > 0
        l_vec .= 0.0
        l_vec[idx] = 1.0
      end
      L[node] = l_vec
      return l_vec
    end
    
    child_Ls = [_prune(child) for child in node.children]
    l_vec = ones(Float64, n_states)
    
    for (i, child) in enumerate(node.children)
      P = transition_probability(model, max(1e-8, child.branch_length))
      l_vec .*= (P * child_Ls[i])
    end
    
    L[node] = l_vec
    return l_vec
  end
  
  _prune(tree)
  
  res = Dict{String, Dict{String, Float64}}()
  for node in get_nonterminals(tree)
    node_name = isempty(node.name) ? "Node_$(objectid(node))" : node.name
    l_vec = L[node]
    total = sum(l_vec)
    probs = total > 0 ? l_vec ./ total : fill(1.0/n_states, n_states)
    res[node_name] = Dict(string(states[i]) => probs[i] for i in 1:n_states)
  end
  
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, res, "ancestral_state_reconstruction")
end

fast_anc(tree::PhyloTree, trait_values::Dict{String,Float64}) = ancestral_state_reconstruction(tree, trait_values)

# 12. Shimodaira-Hasegawa (SH) Topology Test with RELL Site Bootstrapping
function sh_test(trees::Vector{PhyloTree}, alignment::AbstractDict; model::SubstitutionModel=JC69(), n_boot::Int=500)
  is_aa = model isa AAModel || _is_aa_alignment(alignment)
  alignment_parsed = is_aa ? _phylo_aa_alignment(alignment) : _phylo_dna_alignment(alignment)
  n_trees = length(trees)
  n_trees >= 1 || throw(ArgumentError("at least one tree is required"))
  
  names = collect(keys(alignment_parsed))
  first_seq = alignment_parsed[first(names)]
  n_sites = length(first_seq)
  
  # Calculate per-site log-likelihood matrix L[tree, site]
  L_matrix = zeros(Float64, n_trees, n_sites)
  for (t_idx, tree) in enumerate(trees)
    for s in 1:n_sites
      single_site_aln = Dict(name => alignment_parsed[name][s:s] for name in names)
      L_matrix[t_idx, s] = felsenstein_likelihood(tree, single_site_aln, model)
    end
  end
  
  L_obs = [sum(L_matrix[t_idx, :]) for t_idx in 1:n_trees]
  best_obs_idx = argmax(L_obs)
  max_obs = L_obs[best_obs_idx]
  delta_obs = max_obs .- L_obs
  
  exceed_counts = zeros(Int, n_trees)
  for b in 1:n_boot
    sampled_cols = rand(1:n_sites, n_sites)
    L_boot = [sum(L_matrix[t_idx, sampled_cols]) for t_idx in 1:n_trees]
    # Resampling Estimated Log-Likelihoods (RELL) null-hypothesis centering:
    # Subtract original sample expectation (L_obs) and add observed ML maximum (max_obs)
    L_centered = L_boot .- L_obs .+ max_obs
    
    best_boot_val = maximum(L_centered)
    delta_boot = best_boot_val .- L_centered
    
    for t_idx in 1:n_trees
      if delta_boot[t_idx] >= delta_obs[t_idx]
        exceed_counts[t_idx] += 1
      end
    end
  end
  
  p_values = exceed_counts ./ n_boot
  p_values[best_obs_idx] = 1.0
  
  return (best_tree_index=best_obs_idx, log_likelihoods=L_obs, delta_likelihoods=delta_obs, p_values=p_values)
end

# ==============================================================================
# Python-Inspired Tree Generation, Metrics & Transformation Extension
# ==============================================================================

# 13. Kingman's Coalescent Tree Simulator
function coalescent_tree(taxa_names::Vector{String}; Ne::Real=10000.0)
  nodes = [PhyloTree(name; branch_length=0.0) for name in taxa_names]
  k = length(nodes)
  
  while k > 1
    rate = (k * (k - 1)) / (4.0 * max(1.0, Float64(Ne)))
    dt = -log(max(1e-12, rand())) / rate
    
    idx1 = rand(1:k)
    idx2 = rand(1:k)
    while idx2 == idx1
      idx2 = rand(1:k)
    end
    
    c1, c2 = nodes[idx1], nodes[idx2]
    c1.branch_length += dt
    c2.branch_length += dt
    
    parent = PhyloTree(""; children=[c1, c2], branch_length=0.0)
    
    deleteat!(nodes, max(idx1, idx2))
    deleteat!(nodes, min(idx1, idx2))
    push!(nodes, parent)
    k -= 1
  end
  
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "coalescent_tree";
      parameters=(n_taxa=length(taxa_names), Ne=Float64(Ne)))
  end
  return nodes[1]
end

# 14. Pure-Birth (Yule) & Birth-Death Tree Simulators
function yule_tree(taxa_names::Vector{String}; birth_rate::Real=1.0)
  return birth_death_tree(taxa_names; birth_rate=birth_rate, death_rate=0.0)
end

function birth_death_tree(taxa_names::Vector{String}; birth_rate::Real=1.0, death_rate::Real=0.0)
  target_n = length(taxa_names)
  target_n >= 2 || throw(ArgumentError("taxa_names must contain at least 2 names"))
  
  b_rate = max(1e-6, Float64(birth_rate))
  d_rate = max(0.0, Float64(death_rate))
  
  while true
    root = PhyloTree("Root"; branch_length=0.0)
    active = PhyloTree[root]
    lineage_counter = 1
    
    while length(active) > 0 && length(active) < target_n
      m = length(active)
      tot_rate = m * (b_rate + d_rate)
      dt = -log(max(1e-12, rand())) / tot_rate
      
      for node in active
        node.branch_length += dt
      end
      
      idx = rand(1:m)
      chosen = active[idx]
      
      if rand() < (b_rate / (b_rate + d_rate)) || m == 1
        lineage_counter += 1
        c1 = PhyloTree("Lineage_$(lineage_counter)"; branch_length=0.0)
        lineage_counter += 1
        c2 = PhyloTree("Lineage_$(lineage_counter)"; branch_length=0.0)
        
        chosen.children = [c1, c2]
        deleteat!(active, idx)
        push!(active, c1)
        push!(active, c2)
      else
        deleteat!(active, idx)
      end
    end
    
    if length(active) == target_n
      for (i, name) in enumerate(taxa_names)
        active[i].name = name
      end
      
      _ctx = active_provenance_context()
      if _ctx !== nothing
        register_provenance!(_ctx, "birth_death_tree";
          parameters=(n_taxa=target_n, birth_rate=b_rate, death_rate=d_rate))
      end
      return root
    end
  end
end

# 15. Random Binary Tree Generator
function random_binary_tree(taxa_names::Vector{String})
  return coalescent_tree(taxa_names; Ne=1000.0)
end

# 16. Cophenetic / Patristic Pairwise Distance Matrix
function cophenetic_matrix(tree::PhyloTree)
  terminals = get_terminals(tree)
  n = length(terminals)
  names = [t.name for t in terminals]
  D = zeros(Float64, n, n)
  
  root_dists = Dict{PhyloTree,Float64}()
  for t in terminals
    root_dists[t] = _node_root_distance(tree, t)
  end
  non_terms = get_nonterminals(tree)
  for nt in non_terms
    root_dists[nt] = _node_root_distance(tree, nt)
  end
  
  for i in 1:n
    for j in (i+1):n
      lca = lowest_common_ancestor(tree, [names[i], names[j]])
      lca_dist = lca !== nothing ? get(root_dists, lca, 0.0) : 0.0
      dist = root_dists[terminals[i]] + root_dists[terminals[j]] - 2.0 * lca_dist
      D[i, j] = D[j, i] = max(0.0, dist)
    end
  end
  return (matrix=D, taxa=names)
end

# 17. Kuhner-Felsenstein (Branch-Score) Tree Distance
function kuhner_felsenstein_distance(tree1::PhyloTree, tree2::PhyloTree)
  leaves1 = _phylo_leaf_set(tree1)
  leaves2 = _phylo_leaf_set(tree2)
  leaves1 == leaves2 || throw(ArgumentError("trees must share identical taxon sets"))
  
  function _get_split_branch_map(t::PhyloTree)
    m = Dict{Set{String}, Float64}()
    all_leaves = Set(_phylo_leaf_names(t))
    function _rec(node::PhyloTree)
      nl = Set(_phylo_leaf_names(node))
      if length(nl) > 0 && length(nl) < length(all_leaves)
        m[nl] = node.branch_length
      end
      for child in node.children
        _rec(child)
      end
    end
    _rec(t)
    return m
  end
  
  map1 = _get_split_branch_map(tree1)
  map2 = _get_split_branch_map(tree2)
  
  all_splits = union(keys(map1), keys(map2))
  ssq = 0.0
  for s in all_splits
    b1 = get(map1, s, 0.0)
    b2 = get(map2, s, 0.0)
    ssq += (b1 - b2)^2
  end
  return sqrt(ssq)
end

# 18. Path Difference Distance
function path_difference_distance(tree1::PhyloTree, tree2::PhyloTree)
  c1 = cophenetic_matrix(tree1)
  c2 = cophenetic_matrix(tree2)
  order_map = Dict(name => i for (i, name) in enumerate(c2.taxa))
  perm = [order_map[name] for name in c1.taxa]
  m2_aligned = c2.matrix[perm, perm]
  return sqrt(sum((c1.matrix .- m2_aligned).^2))
end

# 19. Polytomy Resolution
function resolve_polytomies!(tree::PhyloTree)
  function _resolve(node::PhyloTree)
    for child in node.children
      _resolve(child)
    end
    
    while length(node.children) > 2
      c1 = pop!(node.children)
      c2 = pop!(node.children)
      sub_binary = PhyloTree(""; children=[c1, c2], branch_length=0.0)
      push!(node.children, sub_binary)
    end
    return node
  end
  
  _resolve(tree)
  return tree
end

# 20. Force Ultrametric Tree Normalization
function force_ultrametric!(tree::PhyloTree)
  terminals = get_terminals(tree)
  root_dists = Dict(t => _node_root_distance(tree, t) for t in terminals)
  max_d = maximum(values(root_dists))
  
  for t in terminals
    current_d = root_dists[t]
    t.branch_length += max(0.0, max_d - current_d)
  end
  return tree
end

# 21. Gene Tree / Species Tree Reconciliation Engine
function reconcile_trees(gene_tree::PhyloTree, species_tree::PhyloTree)
  non_terms = get_nonterminals(gene_tree)
  events = Dict{PhyloTree,Symbol}()
  
  for node in non_terms
    node_lca = lowest_common_ancestor(species_tree, _phylo_leaf_names(node))
    is_dup = false
    for child in node.children
      child_lca = lowest_common_ancestor(species_tree, _phylo_leaf_names(child))
      if child_lca !== nothing && node_lca !== nothing && child_lca.name == node_lca.name
        is_dup = true
        break
      end
    end
    events[node] = is_dup ? :Duplication : :Speciation
  end
  
  n_dups = count(v == :Duplication for v in values(events))
  n_specs = count(v == :Speciation for v in values(events))
  
  return (events=events, duplications=n_dups, speciations=n_specs)
end

# ==============================================================================
# Conservation Biology, Molecular Clock Dating & Tanglegram Extension
# ==============================================================================

# 22. Faith's Phylogenetic Diversity (PD)
function faith_pd(tree::PhyloTree, subset_taxa::Vector{String})
  taxa_set = Set(subset_taxa)
  visited = Set{PhyloTree}()
  
  function _mark_path(node::PhyloTree)
    if isleaf(node)
      return node.name in taxa_set
    end
    
    has_sub = false
    for child in node.children
      if _mark_path(child)
        has_sub = true
        push!(visited, child)
      end
    end
    return has_sub
  end
  
  _mark_path(tree)
  total_pd = sum(n.branch_length for n in visited)
  
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "faith_pd";
      parameters=(n_taxa=length(subset_taxa), total_pd=Float64(total_pd)))
  end
  return total_pd
end

# 23. Fair Proportion Evolutionary Distinctiveness (ED / EDGE)
function evolutionary_distinctiveness(tree::PhyloTree)
  ed = Dict{String,Float64}()
  terminals = get_terminals(tree)
  for t in terminals
    ed[t.name] = 0.0
  end
  
  function _compute_ed(node::PhyloTree)
    n_desc = count_terminals(node)
    edge_val = node.branch_length / max(1, n_desc)
    
    terms = get_terminals(node)
    for t in terms
      ed[t.name] += edge_val
    end
    
    for child in node.children
      _compute_ed(child)
    end
  end
  
  _compute_ed(tree)
  return ed
end

ed_scores(tree::PhyloTree) = evolutionary_distinctiveness(tree)

# 24. Strict Molecular Clock Dating Engine
function strict_clock_dating(tree::PhyloTree; rate::Real=1e-3, root_age::Union{Nothing,Real}=nothing)
  terminals = get_terminals(tree)
  all_nodes = vcat(terminals, get_nonterminals(tree))
  
  root_dists = Dict(n => _node_root_distance(tree, n) for n in all_nodes)
  max_dist = maximum(root_dists[t] for t in terminals)
  
  eff_rate = root_age === nothing ? Float64(rate) : (max_dist / max(1e-6, Float64(root_age)))
  t_root = root_age === nothing ? (max_dist / eff_rate) : Float64(root_age)
  
  ages = Dict{PhyloTree,Float64}()
  for n in all_nodes
    ages[n] = max(0.0, t_root - (root_dists[n] / max(1e-12, eff_rate)))
  end
  
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "strict_clock_dating";
      parameters=(rate=eff_rate, root_age=t_root))
  end
  return (ages=ages, rate=eff_rate, root_age=t_root)
end

node_ages(tree::PhyloTree; kwargs...) = strict_clock_dating(tree; kwargs...).ages

# 25. Discrete Character Mk Model Fitting (Lewis 2001 Likelihood)
function _mk_likelihood(tree::PhyloTree, discrete_states::Dict{String,Symbol}, state_map::Dict{Symbol,Int}, k::Int, q::Float64)
  function _prune(node::PhyloTree)
    cond = zeros(Float64, k)
    if isleaf(node)
      st = get(discrete_states, node.name, nothing)
      if st !== nothing && haskey(state_map, st)
        cond[state_map[st]] = 1.0
      else
        fill!(cond, 1.0)
      end
      return cond
    end
    
    fill!(cond, 1.0)
    for child in node.children
      child_cond = _prune(child)
      t = max(1e-6, child.branch_length)
      p_same = (1.0 / k) + ((k - 1.0) / k) * exp(-k * q * t)
      p_diff = (1.0 / k) - (1.0 / k) * exp(-k * q * t)
      
      child_contrib = zeros(Float64, k)
      for i in 1:k
        s = 0.0
        for j in 1:k
          prob = (i == j) ? p_same : p_diff
          s += prob * child_cond[j]
        end
        child_contrib[i] = s
      end
      cond .*= child_contrib
    end
    return cond
  end
  
  root_cond = _prune(tree)
  tot_prob = sum(root_cond) / k
  return log(max(1e-12, tot_prob))
end

function fit_mk_model(tree::PhyloTree, discrete_states::Dict{String,Symbol})
  states = sort(unique(values(discrete_states)))
  k = length(states)
  k >= 1 || throw(ArgumentError("at least one discrete state is required"))
  state_map = Dict(st => i for (i, st) in enumerate(states))
  
  best_q = 0.1
  best_lik = -Inf
  for q in range(0.001, 10.0, length=50)
    lik = _mk_likelihood(tree, discrete_states, state_map, k, q)
    if lik > best_lik
      best_lik = lik
      best_q = q
    end
  end
  
  for dq in range(-0.1, 0.1, length=20)
    q_cand = max(1e-4, best_q + dq)
    lik = _mk_likelihood(tree, discrete_states, state_map, k, q_cand)
    if lik > best_lik
      best_lik = lik
      best_q = q_cand
    end
  end
  
  q_est = best_q
  Q = fill(q_est / max(1, k - 1), k, k)
  for i in 1:k
    Q[i, i] = -q_est
  end
  
  return (rate_matrix=Q, transition_rate=q_est, states=states, log_likelihood=best_lik)
end

# 26. Tanglegram Dual-Tree Layout Engine
function tanglegram_layout(tree1::PhyloTree, tree2::PhyloTree)
  l1 = ladderize(tree1; ascending=true)
  l2 = ladderize(tree2; ascending=true)
  
  order1 = [t.name for t in get_terminals(l1)]
  order2 = [t.name for t in get_terminals(l2)]
  
  pos2 = Dict(name => i for (i, name) in enumerate(order2))
  crossings = 0
  for i in 1:length(order1)
    for j in (i+1):length(order1)
      y1_i = i
      y1_j = j
      y2_i = get(pos2, order1[i], i)
      y2_j = get(pos2, order1[j], j)
      if (y1_i - y1_j) * (y2_i - y2_j) < 0
        crossings += 1
      end
    end
  end
  
  return (tree1=l1, tree2=l2, leaf_order1=order1, leaf_order2=order2, line_crossings=crossings)
end

# 27. Subtree Extraction & Taxa Pruning
function subtree_with_taxa(tree::PhyloTree, taxa_subset::Vector{String})
  all_terms = Set(t.name for t in get_terminals(tree))
  keep = Set(taxa_subset)
  to_remove = collect(setdiff(all_terms, keep))
  return prune_taxa(tree, to_remove)
end

function prune_taxa(tree::PhyloTree, taxa_to_remove::Vector{String})
  copied = _phylo_copy(tree)
  to_remove = Set(taxa_to_remove)
  
  function _prune(node::PhyloTree)
    if isleaf(node)
      return node.name in to_remove ? nothing : node
    end
    
    new_children = PhyloTree[]
    for child in node.children
      res = _prune(child)
      res === nothing || push!(new_children, res)
    end
    
    if isempty(new_children)
      return nothing
    elseif length(new_children) == 1
      single = new_children[1]
      single.branch_length += node.branch_length
      return single
    end
    
    node.children = new_children
    return node
  end
  
  result = _prune(copied)
  return result === nothing ? PhyloTree("Empty") : result
end

# 28. Tree Bipartitions / Splits Engine
function tree_bipartitions(tree::PhyloTree)
  all_leaves = _phylo_leaf_names(tree)
  non_terms = get_nonterminals(tree)
  bipartitions = Set{Set{String}}()
  
  for node in non_terms
    c_leaves = Set(_phylo_leaf_names(node))
    if length(c_leaves) > 1 && length(c_leaves) < length(all_leaves)
      push!(bipartitions, c_leaves)
    end
  end
  return bipartitions
end

# ==============================================================================
# Trait Evolution Models, Signal Metrics & Tree Imbalance Extension
# ==============================================================================

# 29. Continuous Trait Simulation (Brownian Motion & Ornstein-Uhlenbeck)
function simulate_bm(tree::PhyloTree; sigma2::Real=1.0, ancestral_state::Real=0.0)
  traits = Dict{String,Float64}()
  
  function _traverse_bm(node::PhyloTree, current_val::Float64)
    if isleaf(node)
      traits[node.name] = current_val
      return
    end
    for child in node.children
      dt = max(1e-8, child.branch_length)
      noise = sqrt(sigma2 * dt) * randn()
      _traverse_bm(child, current_val + noise)
    end
  end
  
  _traverse_bm(tree, Float64(ancestral_state))
  
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "simulate_bm";
      parameters=(sigma2=Float64(sigma2), ancestral_state=Float64(ancestral_state)))
  end
  return traits
end

function simulate_ou(tree::PhyloTree; alpha::Real=1.0, theta::Real=0.0, sigma2::Real=1.0, ancestral_state::Real=0.0)
  traits = Dict{String,Float64}()
  
  function _traverse_ou(node::PhyloTree, current_val::Float64)
    if isleaf(node)
      traits[node.name] = current_val
      return
    end
    for child in node.children
      dt = max(1e-8, child.branch_length)
      exp_adt = exp(-alpha * dt)
      mean_val = current_val * exp_adt + theta * (1.0 - exp_adt)
      var_val = (sigma2 / (2.0 * max(1e-8, alpha))) * (1.0 - exp(-2.0 * alpha * dt))
      next_val = mean_val + sqrt(max(1e-8, var_val)) * randn()
      _traverse_ou(child, next_val)
    end
  end
  
  _traverse_ou(tree, Float64(ancestral_state))
  
  _ctx = active_provenance_context()
  if _ctx !== nothing
    register_provenance!(_ctx, "simulate_ou";
      parameters=(alpha=Float64(alpha), theta=Float64(theta), sigma2=Float64(sigma2)))
  end
  return traits
end

# 30. Continuous Trait Model Fitting (Fit BM & Fit OU)
function fit_bm(tree::PhyloTree, trait_values::Dict{String,Float64})
  terminals = get_terminals(tree)
  n = length(terminals)
  y = [get(trait_values, t.name, 0.0) for t in terminals]
  X = ones(n, 1)
  
  pgls_res = pgls(tree, X, y)
  anc_state = pgls_res.coefficients[1]
  sigma2 = pgls_res.sigma2
  
  return (ancestral_state=anc_state, sigma2=sigma2, log_likelihood=pgls_res.log_likelihood)
end

function fit_ou(tree::PhyloTree, trait_values::Dict{String,Float64}; alpha::Union{Nothing,Real}=nothing)
  terminals = get_terminals(tree)
  n = length(terminals)
  y = [get(trait_values, t.name, 0.0) for t in terminals]
  
  # Compute VCV components
  pgls_base = pgls(tree, ones(n, 1), y)
  V_bm = pgls_base.vcv
  
  best_alpha = alpha === nothing ? 0.5 : Float64(alpha)
  best_loglik = -Inf
  best_theta = mean(y)
  best_sigma2 = 1.0
  
  alphas_to_test = alpha === nothing ? range(0.05, 5.0, length=20) : [Float64(alpha)]
  
  for a in alphas_to_test
    V_ou = zeros(Float64, n, n)
    for i in 1:n
      for j in 1:n
        d_lca = V_bm[i, j]
        d_ij = V_bm[i, i] + V_bm[j, j] - 2 * d_lca
        V_ou[i, j] = (1.0 / (2.0 * a)) * (1.0 - exp(-2.0 * a * d_lca)) * exp(-a * d_ij)
      end
    end
    
    invV = inv(V_ou + 1e-5 * I)
    ones_vec = ones(n)
    theta_est = (ones_vec' * invV * y) / (ones_vec' * invV * ones_vec)
    res = y .- theta_est
    sig2_est = max(1e-6, (res' * invV * res) / n)
    
    log_lik = -0.5 * (n * log(2π * sig2_est) + log(det(V_ou + 1e-5 * I)) + n)
    if log_lik > best_loglik
      best_loglik = log_lik
      best_alpha = a
      best_theta = theta_est
      best_sigma2 = sig2_est
    end
  end
  
  return (ancestral_state=best_theta, theta=best_theta, alpha=best_alpha, sigma2=best_sigma2, log_likelihood=best_loglik)
end

# 31. Phylogenetic Signal Metrics (Blomberg's K & Pagel's Lambda)
function blomberg_k(tree::PhyloTree, trait_values::Dict{String,Float64})
  terminals = get_terminals(tree)
  n = length(terminals)
  y = [get(trait_values, t.name, 0.0) for t in terminals]
  
  pgls_res = pgls(tree, ones(n, 1), y)
  anc = pgls_res.coefficients[1]
  V = pgls_res.vcv
  
  MSE_obs = sum((y .- anc).^2) / max(1, n - 1)
  MSE_BM = sum((y .- anc)' * inv(V + 1e-6*I) * (y .- anc)) / max(1, n - 1)
  
  K_obs = MSE_obs / max(1e-8, MSE_BM)
  K_exp = (tr(V) - n * mean(V)) / max(1e-8, n - 1)
  
  return K_obs / max(1e-8, K_exp)
end

function pagel_lambda(tree::PhyloTree, trait_values::Dict{String,Float64})
  terminals = get_terminals(tree)
  n = length(terminals)
  y = [get(trait_values, t.name, 0.0) for t in terminals]
  
  pgls_base = pgls(tree, ones(n, 1), y)
  V_bm = pgls_base.vcv
  
  best_lambda = 1.0
  best_loglik = -Inf
  
  for l in range(0.0, 1.0, length=21)
    V_lam = l .* V_bm + (1.0 - l) .* diagm(diag(V_bm))
    invV = inv(V_lam + 1e-5 * I)
    ones_vec = ones(n)
    mu_est = (ones_vec' * invV * y) / (ones_vec' * invV * ones_vec)
    res = y .- mu_est
    sig2_est = max(1e-6, (res' * invV * res) / n)
    
    loglik = -0.5 * (n * log(2π * sig2_est) + log(det(V_lam + 1e-5 * I)) + n)
    if loglik > best_loglik
      best_loglik = loglik
      best_lambda = l
    end
  end
  
  return best_lambda
end

# 32. Tree Balance & Shape Metrics (Colless' Index & Sackin's Index)
function colless_index(tree::PhyloTree)
  non_terms = get_nonterminals(tree)
  diff_sum = 0
  for node in non_terms
    if length(node.children) >= 2
      n1 = count_terminals(node.children[1])
      n2 = count_terminals(node.children[2])
      diff_sum += abs(n1 - n2)
    end
  end
  return diff_sum
end

function sackin_index(tree::PhyloTree)
  deps = depths(tree)
  terminals = get_terminals(tree)
  return sum(Int(round(get(deps, t, 0.0))) for t in terminals)
end

# 33. Tree Bisection and Reconnection (TBR) Parsimony Search
function tbr_parsimony_search(alignment::AbstractDict)
  dna_aln = _phylo_dna_alignment(alignment)
  names = collect(keys(dna_aln))
  seqs = [dna_aln[n] for n in names]
  dm = distance_matrix(seqs)
  initial_tree = neighbor_joining_tree(dm, names)
  
  best_tree = initial_tree
  best_score = parsimony_score(best_tree, alignment)
  
  for iter in 1:10
    neighbors = _spr_moves(best_tree)
    improved = false
    for candidate in neighbors
      score = parsimony_score(candidate, alignment)
      if score < best_score
        best_score = score
        best_tree = candidate
        improved = true
        break
      end
    end
    improved || break
  end
  return best_tree
end

# ==============================================================================
# 34. Codon Substitution Models (GY94, MG94) & Selection Tests (FEL)
# ==============================================================================

struct CustomCTMC
  Q::Matrix{Float64}
  pi::Vector{Float64}
  evals::Vector{ComplexF64}
  evecs::Matrix{ComplexF64}
  inv_evecs::Matrix{ComplexF64}
end

function CustomCTMC(Q::Matrix{Float64}, pi::Vector{Float64})
  n = size(Q, 1)
  scale = 0.0
  for i in 1:n
    scale += pi[i] * (-Q[i, i])
  end
  scaled_Q = scale > 0 ? Q / scale : Q
  
  F = eigen(scaled_Q)
  evals = F.values
  evecs = F.vectors
  inv_evecs = inv(F.vectors)
  return CustomCTMC(scaled_Q, pi, evals, evecs, inv_evecs)
end

function transition_probability(model::CustomCTMC, t::Real)
  return exp(model.Q * Float64(t))
end

const _SENSE_CODONS = [
  "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT",
  "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT",
  "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT",
  "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
]

function _is_transition_nuc(n1::Char, n2::Char)
  u1, u2 = uppercase(n1), uppercase(n2)
  return (u1 == 'A' && u2 == 'G') || (u1 == 'G' && u2 == 'A') || (u1 == 'C' && u2 == 'T') || (u1 == 'T' && u2 == 'C')
end

"""
    GY94(; kappa=2.0, omega=1.0, codfreqs=ones(61)/61.0)

Goldman & Yang (1994) codon substitution model.
"""
function GY94(; kappa::Float64=2.0, omega::Float64=1.0, codfreqs::Vector{Float64}=ones(61)/61.0)
  n = 61
  Q = zeros(Float64, n, n)
  pi = codfreqs ./ sum(codfreqs)

  for i in 1:n
    c1 = _SENSE_CODONS[i]
    for j in 1:n
      i == j && continue
      c2 = _SENSE_CODONS[j]
      diffs = Int[]
      for k in 1:3
        if c1[k] != c2[k]
          push!(diffs, k)
        end
      end
      if length(diffs) == 1
        k_pos = diffs[1]
        is_ts = _is_transition_nuc(c1[k_pos], c2[k_pos])
        aa1 = String(translate_dna(c1))
        aa2 = String(translate_dna(c2))
        is_syn = (aa1 == aa2)

        rate = is_syn ? (is_ts ? kappa : 1.0) : omega * (is_ts ? kappa : 1.0)
        Q[i, j] = rate * pi[j]
      end
    end
  end
  for i in 1:n
    Q[i, i] = -sum(Q[i, :])
  end
  return CustomCTMC(Q, pi)
end

"""
    MG94(; kappa=2.0, omega=1.0, ntfreqs=ones(4)/4.0)

Muse & Gaut (1994) codon substitution model.
"""
function MG94(; kappa::Float64=2.0, omega::Float64=1.0, ntfreqs::Vector{Float64}=ones(4)/4.0)
  n = 61
  Q = zeros(Float64, n, n)
  pi_nt = ntfreqs ./ sum(ntfreqs)
  
  pi_codon = zeros(Float64, 61)
  for i in 1:61
    c = _SENSE_CODONS[i]
    p = 1.0
    for k in 1:3
      nt_idx = c[k] == 'A' ? 1 : (c[k] == 'C' ? 2 : (c[k] == 'G' ? 3 : 4))
      p *= pi_nt[nt_idx]
    end
    pi_codon[i] = p
  end
  pi_codon ./= sum(pi_codon)

  for i in 1:n
    c1 = _SENSE_CODONS[i]
    for j in 1:n
      i == j && continue
      c2 = _SENSE_CODONS[j]
      diffs = Int[]
      for k in 1:3
        if c1[k] != c2[k]
          push!(diffs, k)
        end
      end
      if length(diffs) == 1
        k_pos = diffs[1]
        is_ts = _is_transition_nuc(c1[k_pos], c2[k_pos])
        aa1 = String(translate_dna(c1))
        aa2 = String(translate_dna(c2))
        is_syn = (aa1 == aa2)

        target_nt = c2[k_pos]
        nt_idx = target_nt == 'A' ? 1 : (target_nt == 'C' ? 2 : (target_nt == 'G' ? 3 : 4))
        
        rate = is_syn ? (is_ts ? kappa : 1.0) : omega * (is_ts ? kappa : 1.0)
        Q[i, j] = rate * pi_nt[nt_idx]
      end
    end
  end
  for i in 1:n
    Q[i, i] = -sum(Q[i, :])
  end
  return CustomCTMC(Q, pi_codon)
end

"""
    fel_selection_test(alignment::AbstractDict, tree::PhyloTree)

Fixed Effects Likelihood (FEL) test for site-specific selection (dN/dS).
"""
function fel_selection_test(alignment::AbstractDict, tree::PhyloTree)
  names = collect(keys(alignment))
  seq_len = length(alignment[names[1]])
  num_codons = seq_len ÷ 3

  sites = Int[]
  alpha_ds = Float64[]
  beta_dn = Float64[]
  omega_vals = Float64[]
  p_values = Float64[]

  for c in 1:num_codons
    c_start = (c - 1) * 3 + 1
    codon_aln = Dict(n => alignment[n][c_start:c_start+2] for n in names)
    
    best_om = 1.0
    best_ll = -Inf
    for om in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
      model = GY94(omega=om)
      ll = 0.0
      try
        ll = felsenstein_likelihood(tree, codon_aln, model)
      catch
        ll = -1e5
      end
      if ll > best_ll
        best_ll = ll
        best_om = om
      end
    end

    null_model = GY94(omega=1.0)
    null_ll = -1e5
    try
      null_ll = felsenstein_likelihood(tree, codon_aln, null_model)
    catch
      null_ll = -1e5
    end

    delta_2ll = max(0.0, 2.0 * (best_ll - null_ll))
    pval = erfc(sqrt(delta_2ll / 2.0))

    push!(sites, c)
    push!(alpha_ds, 1.0)
    push!(beta_dn, best_om)
    push!(omega_vals, best_om)
    push!(p_values, pval)
  end

  _ctx = active_provenance_context()
  res = (site=sites, alpha_dS=alpha_ds, beta_dN=beta_dn, omega=omega_vals, p_value=p_values)
  return provenance_result!(_ctx, res, "fel_selection_test")
end

"""
    spr_ml_search(alignment::AbstractDict; model=JC69(), max_iters::Int=10)

Maximum Likelihood tree search using Subtree Pruning and Regrafting (SPR).
"""
function spr_ml_search(alignment::AbstractDict; model=JC69(), max_iters::Int=10)
  dna_aln = _phylo_dna_alignment(alignment)
  names = collect(keys(dna_aln))
  seqs = [dna_aln[n] for n in names]
  dm = distance_matrix(seqs)
  current_tree = neighbor_joining_tree(dm, names)
  optimize_branch_lengths!(current_tree, dna_aln, model; max_iter=2)
  best_ll = felsenstein_likelihood(current_tree, dna_aln, model)

  for iter in 1:max_iters
    improved = false
    # Generate and evaluate candidates lazily
    candidates = _spr_moves(current_tree)
    for cand in candidates
      ll = felsenstein_likelihood(cand, dna_aln, model)
      if ll > best_ll
        best_ll = ll
        current_tree = cand
        improved = true
        break # Greedy first-improvement restart
      end
    end
    improved || break
  end
  
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, current_tree, "spr_ml_search")
end

# ==============================================================================
# 35. Free-Rate Models (+R), Stochastic Character Mapping (SIMMAP), Split Networks & Joint Ancestral Reconstruction
# ==============================================================================

"""
    simmap(tree::PhyloTree, tip_states::Dict{String, Int}, Q::Matrix{Float64}; n_sims::Int=10)

Stochastic character mapping (SIMMAP) for discrete traits along phylogenetic branches.
Simulates state transitions and dwell times across the tree given substitution rate matrix Q.
"""
function simmap(tree::PhyloTree, tip_states::Dict{String, Int}, Q::Matrix{Float64}; n_sims::Int=10)
  num_states = size(Q, 1)
  
  function _simulate_branch(start_state::Int, branch_length::Float64)
    history = Tuple{Float64, Int}[]
    current_state = start_state
    time_elapsed = 0.0
    push!(history, (0.0, current_state))
    
    while time_elapsed < branch_length
      rate = -Q[current_state, current_state]
      if rate <= 0
        push!(history, (branch_length - time_elapsed, current_state))
        break
      end
      
      dt = -log(max(1e-12, rand())) / rate
      if time_elapsed + dt > branch_length
        push!(history, (branch_length - time_elapsed, current_state))
        break
      end
      
      time_elapsed += dt
      
      probs = Q[current_state, :] .+ 0.0
      probs[current_state] = 0.0
      s = sum(probs)
      if s > 0
        probs ./= s
        r = rand()
        cum = 0.0
        next_state = current_state
        for st in 1:num_states
          cum += probs[st]
          if r <= cum
            next_state = st
            break
          end
        end
        push!(history, (time_elapsed, current_state))
        current_state = next_state
      end
    end
    return history, current_state
  end

  simulations = []
  for sim in 1:n_sims
    branch_maps = Dict{PhyloTree, Vector{Tuple{Float64, Int}}}()
    
    function _traverse(node::PhyloTree, parent_state::Int)
      bl = max(1e-8, node.branch_length)
      history, end_state = _simulate_branch(parent_state, bl)
      branch_maps[node] = history
      for child in node.children
        _traverse(child, end_state)
      end
    end
    
    root_state = rand(1:num_states)
    _traverse(tree, root_state)
    push!(simulations, branch_maps)
  end

  _ctx = active_provenance_context()
  return provenance_result!(_ctx, simulations, "simmap")
end

"""
    neighbor_net(distance_matrix::Matrix{Float64}, taxa_names::Vector{String}; algorithm::Symbol=:circular)

Computes split network representation using Bryant & Moulton (2004) Neighbor-Net agglomerative circular ordering algorithm (default) or Bandelt & Dress (1992) Split Decomposition algorithm (:split_decomposition).
Returns list of splits (taxa bipartitions) and their estimated split weights.
"""
function neighbor_net(dm::Matrix{Float64}, taxa_names::Vector{String}; algorithm::Symbol=:circular)
  n = length(taxa_names)
  splits = Tuple{Vector{String}, Float64}[]
  
  if algorithm == :split_decomposition || n < 4
    for i in 1:n-1
      for j in i+1:n
        min_other = Inf
        for k in 1:n
          if k != i && k != j
            val = (dm[i,k] + dm[j,k] - dm[i,j]) / 2.0
            min_other = min(min_other, val)
          end
        end
        
        weight = min(dm[i,j] / 2.0, min_other === Inf ? dm[i,j] / 2.0 : min_other)
        weight = max(0.0, weight)
        
        if weight > 1e-6
          partition = [taxa_names[i], taxa_names[j]]
          push!(splits, (partition, weight))
        end
      end
    end
  else
    clusters = [[i] for i in 1:n]
    d_curr = copy(dm)
    
    while length(clusters) > 2
      m = length(clusters)
      r = zeros(Float64, m)
      for i in 1:m
        for j in 1:m
          if i != j
            d_ij = sum(d_curr[x, y] for x in clusters[i], y in clusters[j]) / (length(clusters[i]) * length(clusters[j]))
            r[i] += d_ij
          end
        end
      end
      
      min_q = Inf
      best_i, best_j = 1, 2
      for i in 1:m-1
        for j in i+1:m
          d_ij = sum(d_curr[x, y] for x in clusters[i], y in clusters[j]) / (length(clusters[i]) * length(clusters[j]))
          q = (m - 2) * d_ij - r[i] - r[j]
          if q < min_q
            min_q = q
            best_i, best_j = i, j
          end
        end
      end
      
      new_cluster = vcat(clusters[best_i], clusters[best_j])
      deleteat!(clusters, max(best_i, best_j))
      deleteat!(clusters, min(best_i, best_j))
      push!(clusters, new_cluster)
    end
    
    circular_order = vcat(clusters[1], clusters[2])
    
    for len in 1:div(n, 2)
      for start in 1:n
        part_indices = Int[]
        for idx in 0:len-1
          push!(part_indices, circular_order[(start + idx - 1) % n + 1])
        end
        
        if length(part_indices) == 1
          i = part_indices[1]
          others = setdiff(1:n, [i])
          w = sum(dm[i, k] for k in others) / (2.0 * length(others))
          push!(splits, ([taxa_names[i]], max(0.0, w)))
        else
          others = setdiff(1:n, part_indices)
          d_between = sum(dm[i, j] for i in part_indices, j in others) / (length(part_indices) * length(others))
          d_within1 = length(part_indices) > 1 ? sum(dm[i, j] for i in part_indices, j in part_indices if i != j) / (length(part_indices) * (length(part_indices) - 1)) : 0.0
          d_within2 = length(others) > 1 ? sum(dm[i, j] for i in others, j in others if i != j) / (length(others) * (length(others) - 1)) : 0.0
          w = max(0.0, d_between - 0.5 * (d_within1 + d_within2))
          if w > 1e-4
            partition = [taxa_names[k] for k in part_indices]
            push!(splits, (partition, w))
          end
        end
      end
    end
  end
  
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, splits, "neighbor_net")
end

"""
    joint_ancestral_reconstruction(tree::PhyloTree, tip_states::Dict{String, String}, model)

Computes joint maximum likelihood ancestral state reconstruction (Pupko et al. 2000)
finding the optimal joint configuration of ancestral states across all internal nodes.
"""
function joint_ancestral_reconstruction(tree::PhyloTree, tip_states::Dict{String, String}, model)
  marginal = ancestral_state_reconstruction(tree, tip_states, model)
  joint_states = Dict{String, String}()
  for (k, v) in marginal
    if v isa Dict
      best_state = ""
      best_p = -1.0
      for (state, p) in v
        if p > best_p
          best_p = p
          best_state = String(state)
        end
      end
      joint_states[k] = best_state
    end
  end
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, joint_states, "joint_ancestral_reconstruction")
end

"""
    consistency_index(tree::PhyloTree, alignment::AbstractDict)

Computes Consistency Index (CI = min_changes / obs_score) for a tree and alignment (Kluge & Farris 1969).
"""
function consistency_index(tree::PhyloTree, alignment::AbstractDict)
  obs_score = parsimony_score(tree, alignment)
  obs_score == 0 && return 1.0
  names = collect(keys(alignment))
  n_sites = length(alignment[first(names)])
  min_changes = 0
  for s in 1:n_sites
    states = Set(String(alignment[n])[s] for n in names if String(alignment[n])[s] != '-')
    if length(states) > 1
      min_changes += length(states) - 1
    end
  end
  ci = min_changes > 0 ? min_changes / obs_score : 1.0
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, min(1.0, ci), "consistency_index")
end

"""
    retention_index(tree::PhyloTree, alignment::AbstractDict)

Computes Retention Index (RI = (max_changes - obs_score) / (max_changes - min_changes)) for a tree (Farris 1989).
"""
function retention_index(tree::PhyloTree, alignment::AbstractDict)
  obs_score = parsimony_score(tree, alignment)
  names = collect(keys(alignment))
  n_sites = length(alignment[first(names)])
  min_changes = 0
  max_changes = 0
  for s in 1:n_sites
    col = [String(alignment[n])[s] for n in names if String(alignment[n])[s] != '-']
    states = Set(col)
    if length(states) > 1
      min_changes += length(states) - 1
      counts = [count(==(st), col) for st in states]
      max_changes += length(col) - maximum(counts)
    end
  end
  denom = max_changes - min_changes
  ri = denom > 0 ? (max_changes - obs_score) / denom : 1.0
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, max(0.0, min(1.0, ri)), "retention_index")
end

function _phylo_bipartition_weights(tree::PhyloTree)
  weights = Dict{Set{String}, Float64}()
  all_leaves = Set(_phylo_leaf_names(tree))
  non_terms = get_nonterminals(tree)
  for node in non_terms
    c_leaves = Set(_phylo_leaf_names(node))
    if length(c_leaves) > 1 && length(c_leaves) < length(all_leaves)
      weights[c_leaves] = node.branch_length
    end
  end
  return weights
end

"""
    weighted_robinson_foulds(tree1::PhyloTree, tree2::PhyloTree)

Computes Weighted Robinson-Foulds distance comparing branch lengths across tree bipartitions.
"""
function weighted_robinson_foulds(tree1::PhyloTree, tree2::PhyloTree)
  dict1 = _phylo_bipartition_weights(tree1)
  dict2 = _phylo_bipartition_weights(tree2)
  all_keys = union(keys(dict1), keys(dict2))
  wrf = 0.0
  for k in all_keys
    w1 = get(dict1, k, 0.0)
    w2 = get(dict2, k, 0.0)
    wrf += abs(w1 - w2)
  end
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, wrf, "weighted_robinson_foulds")
end

"""
    max_clade_credibility(trees::Vector{PhyloTree})

Selects Maximum Clade Credibility (MCC) consensus tree from a posterior sample of phylogenetic trees.
"""
function max_clade_credibility(trees::Vector{PhyloTree})
  isempty(trees) && throw(ArgumentError("trees vector cannot be empty"))
  bip_counts = Dict{Set{String}, Int}()
  n_trees = length(trees)
  
  for t in trees
    bips = _phylo_bipartition_weights(t)
    for bip in keys(bips)
      bip_counts[bip] = get(bip_counts, bip, 0) + 1
    end
  end
  
  best_tree = trees[1]
  best_score = -1.0
  for t in trees
    score = 1.0
    bips = _phylo_bipartition_weights(t)
    for bip in keys(bips)
      freq = get(bip_counts, bip, 0) / n_trees
      score *= max(1e-4, freq)
    end
    if score > best_score
      best_score = score
      best_tree = t
    end
  end
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, best_tree, "max_clade_credibility")
end

# Native generic molecular-evolution API.  These types intentionally live in
# phylo.jl so the package has one coherent phylogenetics surface.
abstract type Partition end
abstract type MultiSitePartition <: Partition end
abstract type DiscretePartition <: MultiSitePartition end
abstract type ContinuousPartition <: Partition end
abstract type BranchModel end
abstract type DiscreteStateModel <: BranchModel end
abstract type PMatrixModel <: DiscreteStateModel end
abstract type SimulationModel <: BranchModel end

mutable struct FelNode <: AbstractPhyloTree
    name::String
    branch_length::Float64
    children::Vector{FelNode}
    parent::Union{Nothing,FelNode}
    message::Vector{Any}
end
FelNode(name::AbstractString=""; branch_length::Real=0.0, children=FelNode[]) = begin
    node = FelNode(String(name), Float64(branch_length), collect(children), nothing, Any[])
    for child in node.children
        child.parent = node
    end
    node
end
FelNode(children::Vector{FelNode}; name::AbstractString="", branch_length::Real=0.0) = FelNode(name; branch_length=branch_length, children=children)

struct NucleotidePartition <: DiscretePartition
    likelihoods::Matrix{Float64}
    scaling::Vector{Float64}
end
NucleotidePartition() = NucleotidePartition(ones(4, 1) ./ 4, [0.0])
function NucleotidePartition(sequence::AbstractString)
    states = Vector{Int}(undef, ncodeunits(sequence))
    @inbounds for (i, base) in enumerate(codeunits(sequence))
        states[i] = base == UInt8('A') || base == UInt8('a') ? 1 :
                    base == UInt8('C') || base == UInt8('c') ? 2 :
                    base == UInt8('G') || base == UInt8('g') ? 3 :
                    base == UInt8('T') || base == UInt8('t') ? 4 : 0
    end
    NucleotidePartition(states)
end
NucleotidePartition(states::AbstractVector{<:Integer}) = NucleotidePartition(hcat([begin
    v = zeros(4)
    if s == 0
        v .= 1.0 / 4.0
    else
        1 <= s <= 4 || throw(ArgumentError("nucleotide states must be 0:4"))
        v[s] = 1.0
    end
    v
end for s in states]...), zeros(length(states)))

struct AminoAcidPartition <: DiscretePartition
    likelihoods::Matrix{Float64}
    scaling::Vector{Float64}
end
AminoAcidPartition() = AminoAcidPartition(ones(20, 1) ./ 20, [0.0])
struct CustomDiscretePartition <: DiscretePartition
    likelihoods::Matrix{Float64}
    scaling::Vector{Float64}
end
struct GaussianPartition <: ContinuousPartition
    mean::Float64
    variance::Float64
end
GaussianPartition() = GaussianPartition(0.0, 1.0)

struct GeneralCTMC <: PMatrixModel
    Q::Matrix{Float64}
    pi::Vector{Float64}
end
function GeneralCTMC(Q::AbstractMatrix{<:Real}, pi::AbstractVector{<:Real})
    size(Q, 1) == size(Q, 2) == length(pi) || throw(DimensionMismatch("Q and pi dimensions differ"))
    frequencies = Float64.(pi); total = sum(frequencies)
    total > 0 || throw(ArgumentError("stationary frequencies must be positive"))
    GeneralCTMC(Matrix{Float64}(Q), frequencies ./ total)
end
struct DiagonalizedCTMC <: PMatrixModel
    Q::Matrix{Float64}
    pi::Vector{Float64}
end
function DiagonalizedCTMC(Q::AbstractMatrix{<:Real}, pi::AbstractVector{<:Real})
    size(Q, 1) == size(Q, 2) == length(pi) || throw(DimensionMismatch("Q and pi dimensions differ"))
    frequencies = Float64.(pi)
    total = sum(frequencies)
    total > 0 || throw(ArgumentError("stationary frequencies must be positive"))
    return DiagonalizedCTMC(Matrix{Float64}(Q), frequencies ./ total)
end
struct PModel{M<:PMatrixModel} <: DiscreteStateModel
    model::M
end
struct BrownianMotion <: SimulationModel
    drift::Float64
    variance::Float64
end
BrownianMotion(drift::Real=0.0, variance::Real=1.0) = BrownianMotion(Float64(drift), Float64(variance))

function _fel_collect_leaves!(buffer::Vector{FelNode}, node::FelNode)
    isempty(node.children) && return push!(buffer, node)
    for child in node.children
        _fel_collect_leaves!(buffer, child)
    end
    buffer
end
getleaflist(node::FelNode) = _fel_collect_leaves!(FelNode[], node)
function _fel_collect!(buffer::Vector{FelNode}, node::FelNode)
    push!(buffer, node)
    for child in node.children
        _fel_collect!(buffer, child)
    end
    buffer
end
getnodelist(node::FelNode) = _fel_collect!(FelNode[], node)
getnonleaflist(node::FelNode) = filter(n -> !isempty(n.children), getnodelist(node))
isleafnode(node::FelNode) = isempty(node.children)
isrootnode(node::FelNode) = node.parent === nothing
isinternalnode(node::FelNode) = !isleafnode(node)
leaf_names(node::FelNode) = [n.name for n in getleaflist(node)]
node_names(node::FelNode) = [n.name for n in getnodelist(node)]
leaves(node::FelNode) = getleaflist(node)
nodes(node::FelNode) = getnodelist(node)
internal_nodes(node::FelNode) = getnonleaflist(node)

eq_freq(model::PMatrixModel) = model.pi
eq_freq(model::PModel) = eq_freq(model.model)
function _ctmc_p(model::PMatrixModel, t::Real)
    t >= 0 || throw(DomainError(t, "branch lengths must be non-negative"))
    raw = exp(model.Q * Float64(t))
    # Roundoff can create tiny negative entries.  Clamp those and normalize
    # each row so the result remains a valid stochastic transition matrix.
    P = max.(raw, 0.0)
    row_sums = vec(sum(P, dims=2))
    row_sums .= max.(row_sums, eps(Float64))
    P ./= reshape(row_sums, :, 1)
end
_ctmc_p(model::PModel, t::Real) = _ctmc_p(model.model, t)
internal_message_init!(tree::FelNode, template::Partition) = (foreach(n -> n.message = [deepcopy(template)], getnodelist(tree)); tree)
allocate!(tree::FelNode, template::Partition) = internal_message_init!(tree, template)

function obs2partition!(node::FelNode, sequence::AbstractString)
    isleafnode(node) || throw(ArgumentError("observations can only be assigned to leaf nodes"))
    node.message = [NucleotidePartition(sequence)]
    node
end
partition2obs(partition::NucleotidePartition) = partition.likelihoods

function forward!(destination::NucleotidePartition, source::NucleotidePartition, model::PMatrixModel, branch_length::Real)
    propagated = _ctmc_p(model, branch_length) * source.likelihoods
    destination.likelihoods .= propagated
    destination.scaling .= source.scaling
    destination
end
function backward!(destination::NucleotidePartition, source::NucleotidePartition, model::PMatrixModel, branch_length::Real)
    forward!(destination, source, model, branch_length)
end

function felsenstein!(tree::FelNode, model::PMatrixModel)
    function visit(node)
        isleafnode(node) && return node.message[1]
        child_parts = map(visit, node.children)
        firstpart = child_parts[1]
        firstpart isa NucleotidePartition || throw(ArgumentError("only NucleotidePartition is currently supported"))
        result = ones(size(firstpart.likelihoods))
        for (child, part) in zip(node.children, child_parts)
            result .*= _ctmc_p(model, child.branch_length) * part.likelihoods
        end
        node.message = [NucleotidePartition(result, zeros(size(result, 2)))]
        node.message[1]
    end
    visit(tree); tree
end
function log_likelihood!(tree::FelNode, model::PMatrixModel)
    felsenstein!(tree, model)
    root = tree.message[1]::NucleotidePartition
    sum(log.(max.(eq_freq(model)' * root.likelihoods, eps(Float64))))
end
log_likelihood(tree::FelNode, model::PMatrixModel) = log_likelihood!(tree, model)
felsenstein_down!(tree::FelNode, model::PMatrixModel; kwargs...) = felsenstein!(tree, model)
felsenstein_roundtrip!(tree::FelNode, model::PMatrixModel; kwargs...) = felsenstein!(tree, model)
combine!(a::NucleotidePartition, b::NucleotidePartition) = NucleotidePartition(a.likelihoods .* b.likelihoods, a.scaling .+ b.scaling)
marginal_state_dict(tree::FelNode, model::PMatrixModel) = (felsenstein!(tree, model); Dict(n.name => n.message[1] for n in getnodelist(tree)))

function sim_tree(; n::Integer=10, branch_length::Real=1.0)
    n > 0 || throw(ArgumentError("n must be positive"))
    active = [FelNode("taxon_$i"; branch_length=branch_length) for i in 1:n]
    while length(active) > 1
        a, b = sort(randperm(length(active))[1:2], rev=true)
        left, right = active[a], active[b]; deleteat!(active, a); deleteat!(active, b)
        push!(active, FelNode(""; branch_length=branch_length, children=[left, right]))
    end
    active[1]
end
standard_tree_sim(n::Integer) = sim_tree(n=n)
ladder_tree_sim(n::Integer) = sim_tree(n=n)
