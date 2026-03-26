export distance_matrix, neighbor_joining, neighbor_joining_tree, PhyloTree, parse_newick, write_newick, parse_tree, write_tree, parse_phyloxml, write_phyloxml, parse_nexus, write_nexus, parse_nexml, write_nexml, upgma, get_terminals, get_nonterminals, tree_distance, draw_ascii, draw_unicode, tree_to_dot, tree_to_mermaid, prune, root_with_outgroup, reroot, midpoint_root, get_parent, lowest_common_ancestor, common_ancestor, is_monophyletic, robinson_foulds_distance, bootstrap_trees, tree_consensus, consensus_tree, strict_consensus_tree, majority_consensus_tree, bootstrap_consensus_tree, bootstrap_support, set_metadata!, annotate_tree!, is_terminal, is_parent_of, get_path, trace, is_bifurcating, is_preterminal, total_branch_length, depths, find_clades, ladderize, count_terminals, collapse_clades, parsimony_score, maximum_parsimony_tree, parsimony_tree, JC69, K80, HKY85, felsenstein_likelihood, transition_probability, maximum_likelihood_tree, coordinates

const _HAS_CUDA_PHYLO = false
using LinearAlgebra

mutable struct PhyloTree
    name::String
    branch_length::Float64
    children::Vector{PhyloTree}
    support::Float64
    metadata::Dict{String,String}
end

PhyloTree(
    name::AbstractString;
    branch_length::Real=0.0,
    children::AbstractVector{PhyloTree}=PhyloTree[],
    support::Real=0.0,
    metadata::AbstractDict=Dict{String,String}(),
) = PhyloTree(String(name), Float64(branch_length), Vector{PhyloTree}(children), Float64(support), Dict{String,String}(string.(keys(metadata)) .=> string.(values(metadata))))

PhyloTree(
    children::AbstractVector{PhyloTree};
    name::AbstractString="",
    branch_length::Real=0.0,
    support::Real=0.0,
    metadata::AbstractDict=Dict{String,String}(),
) = PhyloTree(String(name), Float64(branch_length), Vector{PhyloTree}(children), Float64(support), Dict{String,String}(string.(keys(metadata)) .=> string.(values(metadata))))

isleaf(tree::PhyloTree) = isempty(tree.children)

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
    return _phylo_collect_terms(tree, PhyloTree[])
end

function _phylo_collect_internal(tree::PhyloTree, buffer::Vector{PhyloTree})
    isleaf(tree) || push!(buffer, tree)
    for child in tree.children
        _phylo_collect_internal(child, buffer)
    end
    return buffer
end

function get_nonterminals(tree::PhyloTree)
    return _phylo_collect_internal(tree, PhyloTree[])
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
    return strip(text[start_index:index-1]), index
end

function _phylo_label_to_support_or_name(label::AbstractString)
    stripped = strip(String(label))
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
    token = strip(text[start_index:index-1])
    isempty(token) && return 0.0, index
    return parse(Float64, token), index
end

function _phylo_clade_key(names::AbstractVector{String})
    return join(sort(collect(names)), '\u001f')
end

function _phylo_key_set(key::AbstractString)
    parts = Base.split(String(key), '\u001f')
    return Set(String.(parts))
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

function get_parent(tree::PhyloTree, target_name::AbstractString)
    return _phylo_find_parent(tree, nothing, String(target_name))
end

function lowest_common_ancestor(tree::PhyloTree, names::AbstractVector{<:AbstractString})
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

lowest_common_ancestor(tree::PhyloTree, left_name::AbstractString, right_name::AbstractString) = lowest_common_ancestor(tree, [left_name, right_name])

common_ancestor(tree::PhyloTree, names) = lowest_common_ancestor(tree, names)

is_terminal(node::PhyloTree) = isleaf(node)

function is_parent_of(tree::PhyloTree, parent_name::AbstractString, child_name::AbstractString)
    parent = _phylo_path_to_named_node(tree, String(parent_name))
    child = _phylo_path_to_named_node(tree, String(child_name))
    parent === nothing && return false
    child === nothing && return false
    parent_node = parent[end]
    child_node = child[end]
    return any(child_ref -> child_ref === child_node, parent_node.children)
end

function get_path(tree::PhyloTree, target_name::AbstractString)
    path = _phylo_path_to_named_node(tree, String(target_name))
    path === nothing && throw(ArgumentError("unknown leaf: $(target_name)"))
    return path
end

function trace(tree::PhyloTree, start_name::AbstractString, end_name::AbstractString)
    left_path = get_path(tree, start_name)
    right_path = get_path(tree, end_name)
    shared_index = 0
    for (index, (left_node, right_node)) in enumerate(zip(left_path, right_path))
        left_node === right_node || break
        shared_index = index
    end
    left_segment = reverse(left_path[shared_index:end])
    right_segment = right_path[shared_index+1:end]
    return vcat(left_segment, right_segment)
end

function split(tree::PhyloTree, parent_name::AbstractString; child_names::AbstractVector{<:AbstractString}=String[], branch_length::Real=0.0, support::Real=0.0)
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
            if i == j; continue; end
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
    
    for i in 1:4; Q[i, i] = -sum(Q[i, :]); end
    return exp(Q * t)
end

# State mapping for DNA
const _DNA_STATE_MAP = Dict{Char, Vector{Float64}}(
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

function _get_leaf_likelihoods(sequence::AbstractString)
    n_sites = ncodeunits(sequence)
    L = zeros(Float64, 4, n_sites)
    for (i, byte) in enumerate(codeunits(sequence))
        state = get(_DNA_STATE_MAP, uppercase(Char(byte)), _DNA_AMBIGUOUS_STATE)
        @inbounds for j in 1:4
            L[j, i] = state[j]
        end
    end
    return L
end

"""
    felsenstein_likelihood(tree::PhyloTree, alignment::Dict{String, String}, model::SubstitutionModel)

Compute the log-likelihood of a phylogenetic tree given an alignment using Felsenstein's pruning algorithm.
"""
function felsenstein_likelihood(tree::PhyloTree, alignment::Dict{String, String}, model::SubstitutionModel)
    leaf_likelihoods = Dict{PhyloTree, Matrix{Float64}}()
    terminals = get_terminals(tree)
    for leaf in terminals
        seq = get(alignment, leaf.name, "")
        if isempty(seq)
            # Handle missing sequences as all ambiguities
            n_sites = length(first(values(alignment)))
            leaf_likelihoods[leaf] = ones(Float64, 4, n_sites)
        else
            leaf_likelihoods[leaf] = _get_leaf_likelihoods(seq)
        end
    end

    n_sites = size(first(values(leaf_likelihoods)), 2)
    transition_cache = Dict{Float64, Matrix{Float64}}()

    transition_cached(t::Real) = get!(transition_cache, Float64(t)) do
        transition_probability(model, t)
    end
    
    function _compute_internal_likelihood(node::PhyloTree)
        if isleaf(node)
            return leaf_likelihoods[node]
        end

        # Post-order: compute children first
        child_Ls = [_compute_internal_likelihood(child) for child in node.children]
        
        node_L = ones(Float64, 4, n_sites)
        
        for (i, child) in enumerate(node.children)
            P = transition_cached(child.branch_length)
            # Sum over possible states at child node
            # L_p(i) = sum_j P_ij(t) * L_c(j)
            node_L .*= (P * child_Ls[i])
        end
        return node_L
    end

    root_L = _compute_internal_likelihood(tree)
    
    # Final likelihood: sum over root states (assuming equal frequencies 0.25 for JC69/K80)
    # For HKY85, use the model's pi
    pi = model isa HKY85 ? model.pi : fill(0.25, 4)
    
    log_lik = 0.0
    for s in 1:n_sites
        site_lik = sum(pi .* root_L[:, s])
        log_lik += log(max(site_lik, 1e-300))
    end
    
    return log_lik
end

"""
    optimize_branch_lengths!(tree::PhyloTree, alignment::Dict{String, String}, model::SubstitutionModel; max_iter=10)

Optimize branch lengths of a given tree to maximize likelihood using a simple coordinate ascent approach.
"""
function optimize_branch_lengths!(tree::PhyloTree, alignment::Dict{String, String}, model::SubstitutionModel; max_iter=5)
    nodes = vcat(PhyloTree[tree], get_nonterminals(tree))
    
    for iter in 1:max_iter
        improved = false
        for node in nodes
            for child in node.children
                # Simple grid search or golden section could be used
                # Here we use a very simple approach for demonstration
                best_bl = child.branch_length
                best_lik = felsenstein_likelihood(tree, alignment, model)
                
                for scale in [0.8, 0.9, 1.1, 1.2]
                    old_bl = child.branch_length
                    child.branch_length = max(0.0001, old_bl * scale)
                    new_lik = felsenstein_likelihood(tree, alignment, model)
                    if new_lik > best_lik
                        best_lik = new_lik
                        best_bl = child.branch_length
                        improved = true
                    else
                        child.branch_length = old_bl
                    end
                end
            end
        end
        if !improved; break; end
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
                            if other_i == i; continue; end
                            
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
    maximum_likelihood_tree(alignment::Dict{String, String}; model::SubstitutionModel=JC69(), initial_tree::Union{Nothing, PhyloTree}=nothing)

Find the Maximum Likelihood tree for a given alignment.
"""
function maximum_likelihood_tree(alignment::Dict{String, String}; model::SubstitutionModel=JC69(), initial_tree::Union{Nothing, PhyloTree}=nothing, max_topology_iter=10)
    if initial_tree === nothing
        names = collect(keys(alignment))
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
            # Use NJ as starting point for larger alignments
            seqs = String[alignment[name] for name in names]
            dm = distance_matrix(seqs)
            initial_tree = neighbor_joining_tree(dm, names)
        end
    end
    
    current_tree = _phylo_copy(initial_tree)
    optimize_branch_lengths!(current_tree, alignment, model)
    current_lik = felsenstein_likelihood(current_tree, alignment, model)

        if count_terminals(current_tree) <= 3
            return current_tree
        end
    
    for iter in 1:max_topology_iter
        best_neighbor = current_tree
        best_neighbor_lik = current_lik
        
        # Explore NNI space
        neighbors = _nni_moves(current_tree)
        for neighbor in neighbors
            optimize_branch_lengths!(neighbor, alignment, model; max_iter=2)
            n_lik = felsenstein_likelihood(neighbor, alignment, model)
            if n_lik > best_neighbor_lik
                best_neighbor_lik = n_lik
                best_neighbor = neighbor
            end
        end
        
        if best_neighbor_lik > current_lik
            current_tree = best_neighbor
            current_lik = best_neighbor_lik
            optimize_branch_lengths!(current_tree, alignment, model; max_iter=5)
            current_lik = felsenstein_likelihood(current_tree, alignment, model)
        else
            break
        end
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
            indices = sortperm(leaf_counts; rev=!ascending)
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

function _phylo_path_between(tree::PhyloTree, left_name::AbstractString, right_name::AbstractString)
    left_path = _phylo_path_to_leaf(tree, String(left_name))
    right_path = _phylo_path_to_leaf(tree, String(right_name))
    left_path === nothing && throw(ArgumentError("unknown leaf: $(left_name)"))
    right_path === nothing && throw(ArgumentError("unknown leaf: $(right_name)"))

    shared_index = 0
    for (index, (left_node, right_node)) in enumerate(zip(left_path, right_path))
        left_node === right_node || break
        shared_index = index
    end

    path_nodes = vcat(reverse(left_path[shared_index:end]), right_path[shared_index+1:end])
    left_lengths = [left_path[index + 1].branch_length for index in reverse(shared_index:(length(left_path) - 1))]
    right_lengths = [right_path[index + 1].branch_length for index in shared_index:(length(right_path) - 1)]
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
    stripped = strip(text)
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
    return parse_newick(neighbor_joining(D, names))
end

function midpoint_root(tree::PhyloTree)
    terminals = get_terminals(tree)
    length(terminals) <= 1 && return _phylo_copy(tree)

    best_left = terminals[1]
    best_right = terminals[2]
    best_distance = -Inf
    for i in 1:length(terminals)-1
        for j in i+1:length(terminals)
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
            return _phylo_clone_from(path_nodes[edge_index + 1], nothing, 0.0, adjacency)
        elseif midpoint < accumulated + edge_length
            left_length = midpoint - accumulated
            right_length = edge_length - left_length
            left_subtree = _phylo_clone_from(path_nodes[edge_index], path_nodes[edge_index + 1], left_length, adjacency)
            right_subtree = _phylo_clone_from(path_nodes[edge_index + 1], path_nodes[edge_index], right_length, adjacency)
            return PhyloTree(""; branch_length=0.0, children=[left_subtree, right_subtree], support=0.0)
        end
        accumulated += edge_length
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
        for i in 1:cluster_count-1
            for j in i+1:cluster_count
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

    return clusters[1]
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

function tree_distance(tree::PhyloTree, left_name::AbstractString, right_name::AbstractString)
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

function tree_to_dot(tree::PhyloTree; graph_name::AbstractString="phylo_tree")
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
    if keep_names isa AbstractString
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

function root_with_outgroup(tree::PhyloTree, outgroup_name::AbstractString)
    target = String(outgroup_name)
    leaf_names = Set(_phylo_leaf_names(tree))
    target in leaf_names || throw(ArgumentError("unknown outgroup: $(outgroup_name)"))
    target_tree = prune(tree, [target])
    remainder_names = setdiff(leaf_names, Set([target]))
    remainder = isempty(remainder_names) ? nothing : prune(tree, remainder_names)
    remainder === nothing && return _phylo_copy(target_tree)
    return PhyloTree(""; branch_length=0.0, children=[target_tree, remainder], support=1.0)
end

reroot(tree::PhyloTree, outgroup_name::AbstractString) = root_with_outgroup(tree, outgroup_name)

function _bootstrap_alignment_sequences(alignment)
    width = get_alignment_length(alignment)
    width == 0 && return [record.sequence for record in alignment.records]

    sampled_columns = rand(1:width, width)
    return [begin
        buffer = Vector{UInt8}(undef, width)
        for (position, column) in enumerate(sampled_columns)
            buffer[position] = UInt8(record.sequence[column])
        end
        String(buffer)
    end for record in alignment.records]
end

function bootstrap_trees(alignment; replicates::Int=100, method::Symbol=:hamming, constructor::Symbol=:nj)
    replicates > 0 || throw(ArgumentError("replicates must be positive"))
    trees = PhyloTree[]
    names = [record.identifier for record in alignment.records]

    for _ in 1:replicates
        sequences = _bootstrap_alignment_sequences(alignment)
        D = distance_matrix(sequences; method=method)
        tree = constructor == :nj ? neighbor_joining_tree(D, names) : constructor == :upgma ? upgma(D, names) : throw(ArgumentError("constructor must be :nj or :upgma"))
        push!(trees, tree)
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
        for key in sort(maximal; by = x -> -length(Base.split(x, '\u001f')))
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

    return _build(universe)
end

consensus_tree(trees::AbstractVector{PhyloTree}; kwargs...) = tree_consensus(trees; kwargs...)
strict_consensus_tree(trees::AbstractVector{PhyloTree}) = tree_consensus(trees; threshold=1.0)
majority_consensus_tree(trees::AbstractVector{PhyloTree}) = tree_consensus(trees; threshold=0.5)

function bootstrap_consensus_tree(alignment; replicates::Int=100, threshold::Real=0.5, method::Symbol=:hamming, constructor::Symbol=:nj)
    trees = bootstrap_trees(alignment; replicates=replicates, method=method, constructor=constructor)
    return tree_consensus(trees; threshold=threshold)
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

function set_metadata!(tree::PhyloTree, key::AbstractString, value::AbstractString)
    tree.metadata[String(key)] = String(value)
    return tree
end

function annotate_tree!(tree::PhyloTree, metadata::AbstractDict)
    for (key, value) in metadata
        tree.metadata[string(key)] = string(value)
    end
    return tree
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
        println(io, indent, "  <property ref=\"", _xml_escape(key), "\">", _xml_escape(value), "</property>")
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
    startswith(stripped, "</") && return strip(replace(stripped[3:end-1], r"\s+.*$" => ""))
    startswith(stripped, "<") || return ""
    cleaned = replace(stripped[2:end-1], r"/\s*$" => "")
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

function parse_phyloxml(text::AbstractString)
    tokens = collect(eachmatch(r"<[^>]+>|[^<]+", String(text)))
    root = nothing
    node_stack = PhyloTree[]
    current_field = ""
    current_property = ""

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
                end
            else
                tag = _phyloxml_tag_name(token)
                if tag == "clade"
                    push!(node_stack, PhyloTree(""; metadata=Dict{String,String}()))
                elseif tag == "name"
                    current_field = "name"
                elseif tag == "branch_length"
                    current_field = "branch_length"
                elseif tag == "confidence"
                    current_field = "confidence"
                elseif tag == "property"
                    current_field = "property"
                    current_property = _phyloxml_attr_value(token, "ref")
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
                node.metadata[current_property] = _xml_unescape(value)
            end
        end
    end

    root === nothing && throw(ArgumentError("invalid PhyloXML document"))
    return root
end

function write_nexus(tree::PhyloTree; tree_name::AbstractString="tree_1")
    return join([
        "#NEXUS",
        "begin trees;",
        "  tree $(tree_name) = $(write_newick(tree))",
        "end;",
    ], "\n")
end

function parse_nexus(text::AbstractString)
    for line in eachline(IOBuffer(String(text)))
        stripped = strip(line)
        startswith(lowercase(stripped), "tree ") || continue
        eq_index = findfirst('=', stripped)
        eq_index === nothing && continue
        payload = strip(stripped[eq_index + 1:end])
        endswith(payload, ";") || continue
        return parse_newick(payload)
    end
    throw(ArgumentError("invalid Nexus tree block"))
end

function write_nexml(tree::PhyloTree; tree_id::AbstractString="tree1")
    io = IOBuffer()
    println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
    println(io, "<nexml xmlns=\"http://www.nexml.org/2009\">")
    println(io, "  <trees id=\"trees1\">")
    println(io, "    <tree id=\"", _xml_escape(tree_id), "\" rooted=\"true\">")
    println(io, "      <newick>", _xml_escape(write_newick(tree)), "</newick>")
    for (key, value) in sort(collect(tree.metadata); by=first)
        println(io, "      <meta property=\"", _xml_escape(key), "\">", _xml_escape(value), "</meta>")
    end
    println(io, "    </tree>")
    println(io, "  </trees>")
    println(io, "</nexml>")
    return String(take!(io))
end

function parse_nexml(text::AbstractString)
    match_result = match(r"<newick>(.+?)</newick>"is, String(text))
    match_result === nothing && throw(ArgumentError("invalid NeXML tree document"))
    return parse_newick(_xml_unescape(match_result.captures[1]))
end

parse_tree(text::AbstractString; format::Symbol=:newick) = format == :newick ? parse_newick(text) : format == :phyloxml ? parse_phyloxml(text) : format == :nexus ? parse_nexus(text) : format == :nexml ? parse_nexml(text) : throw(ArgumentError("unsupported tree format: $(format)"))
write_tree(tree::PhyloTree; format::Symbol=:newick) = format == :newick ? write_newick(tree) : format == :phyloxml ? write_phyloxml(tree) : format == :nexus ? write_nexus(tree) : format == :nexml ? write_nexml(tree) : throw(ArgumentError("unsupported tree format: $(format)"))

function _dna_transition(left::Char, right::Char)
    (left == 'A' && right == 'G') || (left == 'G' && right == 'A') || (left == 'C' && right == 'T') || (left == 'T' && right == 'C')
end

function _pairwise_dna_distance(sequence_left::AbstractString, sequence_right::AbstractString, model::Symbol)
    length(sequence_left) == length(sequence_right) || throw(ArgumentError("sequences must be aligned and equal length for DNA substitution distances"))
    valid_sites = 0
    mismatches = 0
    transitions = 0
    transversions = 0

    @inbounds for index in eachindex(sequence_left)
        left = uppercase(sequence_left[index])
        right = uppercase(sequence_right[index])
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

function _pairwise_protein_distance(sequence_left::AbstractString, sequence_right::AbstractString, matrix_name::Symbol)
    length(sequence_left) == length(sequence_right) || throw(ArgumentError("sequences must be aligned and equal length for protein distances"))
    scoring = named_substitution_matrix(matrix_name)
    best_score = 0.0
    worst_score = 0.0
    observed_score = 0.0
    min_score = minimum(scoring.scores)

    @inbounds for index in eachindex(sequence_left)
        left = sequence_left[index]
        right = sequence_right[index]
        observed_score += _pairwise_residue_score(scoring, left, right)
        best_score += max(_pairwise_residue_score(scoring, left, left), _pairwise_residue_score(scoring, right, right))
        worst_score += min_score
    end

    denom = best_score - worst_score
    denom <= 0 && return 0.0
    return clamp((best_score - observed_score) / denom, 0.0, 1.0)
end

"""
    distance_matrix(sequences::Vector{<:AbstractString}; method=:hamming)

Computes an N x N symmetric distance matrix for a list of sequences.
If `method=:hamming`, it leverages SIMD/SWAR byte-mismatch tracking.
If `method=:alignment`, it computes Needleman-Wunsch identity inversion.
"""
function distance_matrix(sequences::Vector{<:AbstractString}; method=:hamming, use_threads::Bool=true, use_cuda::Bool=false)
    N = length(sequences)
    D = zeros(Float64, N, N)

    if use_cuda
        if isdefined(@__MODULE__, :CUDA)
            _ensure_cuda_phylo!()
            return Base.invokelatest(_CUDA_PHYLO_DISTANCE_MATRIX_IMPL[], sequences; method=method)
        end
        use_threads = true
    end

    # ─── Threaded / Default Support ─────────────────────────────────────────
    if method == :hamming || method == :p_distance
        if use_threads
            Threads.@threads for i in 1:N
                seq_i = sequences[i]
                len_i = length(seq_i)
                @inbounds for j in (i+1):N
                    dist = hamming_distance(seq_i, sequences[j])
                    val = Float64(dist) / len_i
                    D[i, j] = val
                    D[j, i] = val
                end
            end
        else
            for i in 1:N
                seq_i = sequences[i]
                len_i = length(seq_i)
                @inbounds for j in (i+1):N
                    dist = hamming_distance(seq_i, sequences[j])
                    val = Float64(dist) / len_i
                    D[i, j] = val
                    D[j, i] = val
                end
            end
        end
    elseif method == :jukes_cantor || method == :kimura2p || method == :kimura
        if use_threads
            Threads.@threads for i in 1:N
                seq_i = sequences[i]
                @inbounds for j in (i+1):N
                    val = _pairwise_dna_distance(seq_i, sequences[j], method)
                    D[i, j] = val
                    D[j, i] = val
                end
            end
        else
            for i in 1:N
                seq_i = sequences[i]
                @inbounds for j in (i+1):N
                    val = _pairwise_dna_distance(seq_i, sequences[j], method)
                    D[i, j] = val
                    D[j, i] = val
                end
            end
        end
    elseif method == :blosum62 || method == :pam250 || method == :pam70 || method == :protein
        protein_model = method == :protein ? :BLOSUM62 : Symbol(uppercase(String(method)))
        if use_threads
            Threads.@threads for i in 1:N
                seq_i = sequences[i]
                @inbounds for j in (i+1):N
                    val = _pairwise_protein_distance(seq_i, sequences[j], protein_model)
                    D[i, j] = val
                    D[j, i] = val
                end
            end
        else
            for i in 1:N
                seq_i = sequences[i]
                @inbounds for j in (i+1):N
                    val = _pairwise_protein_distance(seq_i, sequences[j], protein_model)
                    D[i, j] = val
                    D[j, i] = val
                end
            end
        end
    elseif method == :alignment
        if use_threads
            Threads.@threads for i in 1:N
                seq_i = sequences[i]
                @inbounds for j in (i+1):N
                    res = pairwise_align(seq_i, sequences[j]; is_local=false)
                    val = 1.0 - res.identity
                    D[i, j] = val
                    D[j, i] = val
                end
            end
        else
            for i in 1:N
                seq_i = sequences[i]
                @inbounds for j in (i+1):N
                    res = pairwise_align(seq_i, sequences[j]; is_local=false)
                    val = 1.0 - res.identity
                    D[i, j] = val
                    D[j, i] = val
                end
            end
        end
    else
        throw(ArgumentError("Unknown method: \$method. Options: :hamming, :p_distance, :jukes_cantor, :kimura2p, :kimura, :blosum62, :pam250, :pam70, :protein, :alignment"))
    end
    
    return D
end

function _alignment_names(alignment)
    return [record.identifier == "" ? (record.name == "" ? "sequence_$(index)" : record.name) : record.identifier for (index, record) in enumerate(alignment.records)]
end

function _alignment_columns(alignment)
    width = get_alignment_length(alignment)
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
    width = get_alignment_length(alignment)
    total_steps = zeros(Int, width)

    if threaded && width > 1 && Threads.nthreads() > 1
        Threads.@threads for col in 1:width
            column_values = Dict{String,Char}()
            for (name, record) in zip(names, alignment.records)
                column_values[name] = record.sequence[col]
            end
            _, steps = _fitch_sets(tree, column_values)
            total_steps[col] = steps
        end
    else
        for col in 1:width
            column_values = Dict{String,Char}()
            for (name, record) in zip(names, alignment.records)
                column_values[name] = record.sequence[col]
            end
            _, steps = _fitch_sets(tree, column_values)
            total_steps[col] = steps
        end
    end

    return sum(total_steps)
end

function _binary_tree_combinations(names::Vector{String})
    length(names) == 1 && return [PhyloTree(names[1]; branch_length=0.0, support=0.0)]
    length(names) == 2 && return [PhyloTree(""; branch_length=0.0, children=[PhyloTree(names[1]; branch_length=0.0), PhyloTree(names[2]; branch_length=0.0)], support=0.0)]

    first_name = names[1]
    remaining = names[2:end]
    candidates = PhyloTree[]
    total_masks = 1 << length(remaining)

    for mask in 0:(total_masks - 2)
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
    sequences = [record.sequence for record in alignment.records]
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

    seed = neighbor_joining_tree(distance_matrix(sequences; method=:p_distance), names)
    return seed
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
    return "($(nodes[1]):$(round(dist_final/2, digits=5)),$(nodes[2]):$(round(dist_final/2, digits=5)));"
end
