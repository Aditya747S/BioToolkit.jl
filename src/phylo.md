# phylo.jl

## Purpose
This file implements the phylogenetic tree layer of BioToolkit. It defines the recursive tree model, tree-building algorithms, tree import and export functions, comparison metrics, and coordinate helpers for plotting.

## Main struct
- PhyloTree: a recursive node with a name, branch length, children, support value, and metadata dictionary.

## Public functions
- get_terminals and get_nonterminals: collect leaf and internal nodes.
- neighbor_joining, neighbor_joining_tree, upgma: build trees from distance information.
- parse_newick, write_newick, parse_phyloxml, write_phyloxml, parse_nexus, write_nexus, parse_nexml, write_nexml: tree I/O for common phylogenetic formats.
- prune, root_with_outgroup, reroot, midpoint_root: manipulate tree topology and rooting.
- get_parent, lowest_common_ancestor, common_ancestor: relationship queries.
- is_monophyletic, is_bifurcating, is_preterminal, is_terminal, is_parent_of: tree property checks.
- tree_distance and robinson_foulds_distance: compare trees.
- bootstrap_trees, tree_consensus, consensus_tree, strict_consensus_tree, majority_consensus_tree, bootstrap_consensus_tree, bootstrap_support: bootstrap and consensus helpers.
- set_metadata!, annotate_tree!: add labels or metadata to trees.
- trace, get_path, find_clades, ladderize, collapse_clades, count_terminals, total_branch_length, depths: traversal and structural helpers.
- parsimony_score, maximum_parsimony_tree, parsimony_tree: parsimony-related routines.
- JC69, K80, HKY85, felsenstein_likelihood, transition_probability, maximum_likelihood_tree: substitution-model and likelihood routines.
- coordinates: compute 2D coordinates for plotting.

## What the module does
The module can build trees from distances, load trees from standard file formats, and manipulate the resulting tree topology. It also gives BioToolkit the comparison and inference routines needed for evolutionary analysis, from simple distance methods to likelihood and parsimony workflows.

## How the struct is used
PhyloTree is the central data structure. Every tree function either builds it, transforms it, traverses it, or compares it. The coordinates helper turns a tree into a drawable layout, which is why the plotting extension can render trees directly from the same structure.

## Typical usage
1. Parse or construct a PhyloTree.
2. Use get_terminals, prune, reroot, or ladderize to inspect or reshape the tree.
3. Use comparison functions such as robinson_foulds_distance when comparing trees.
4. Build consensus trees or bootstrap support summaries when analyzing replicate trees.
5. Call coordinates or the plotting extension when you need a visual representation.

## Important implementation details
- The tree is mutable, which makes rerooting and metadata annotation practical.
- The module supports multiple common phylogenetic file formats so trees can be exchanged with other tools.
- Likelihood and parsimony methods are included so the package is not limited to simple topology operations.

## Threading notes
- `distance_matrix(...; use_threads=true)` already defaults to threaded pairwise distance evaluation for the heavy matrix cases.
- `parsimony_score()` now defaults to threaded column-wise scoring.
- `maximum_parsimony_tree()` now defaults to threaded scoring of exact candidate trees when the exact-search branch is used.

## Why this file matters
Phylogenetic analysis is a major biological workflow, and this module captures both inference and representation in one place. It is the structural basis for tree parsing, comparison, and plotting in BioToolkit.
# Phylogenetic Networking (phylo.jl)

BioToolkit provides native structural Neighbor-Joining phylogeny derived from lightning-fast DNA pairwise distance matrices.

### Distance Matrices
The `distance_matrix(sequences::Vector{String}; method=:hamming)` dynamically maps input sequence arrays into an `N x N` Float64 mathematical matrix.
- `method=:hamming` invokes our `hamming_distance` kernel, processing aligned reads via SWAR-vectorization. Highly optimized for identical-length MSA columns.
- `method=:alignment` runs `pairwise_align` on every combination and computes total length identity ratios mapping to a physical distance bound.

### Neighbor-Joining (NJ)
The `neighbor_joining(D::Matrix{Float64}, names::Vector{String})` function consumes the raw symmetric matrix map to construct topological dependencies un-rooted phylogenetic structures.

### Tree Objects and Utilities
- `PhyloTree` stores rooted tree structure with names, branch lengths, and children.
- `parse_newick` and `write_newick` round-trip simple Newick trees.
- `parse_tree` and `write_tree` dispatch across Newick, PhyloXML, Nexus, and NeXML.
- `parse_phyloxml`, `write_phyloxml`, `parse_nexus`, `write_nexus`, `parse_nexml`, and `write_nexml` provide format-specific tree I/O.
- `neighbor_joining_tree` converts the NJ output into a tree object.
- `upgma` builds a rooted ultrametric tree from a distance matrix.
- `set_metadata!` and `annotate_tree!` attach arbitrary metadata to tree nodes.
- `get_terminals` and `get_nonterminals` expose traversal helpers.
- `count_terminals` returns the number of leaf nodes in a tree.
- `tree_distance` computes the path length between two named tips.
- `distance_matrix` supports Hamming, p-distance, Jukes-Cantor, Kimura2P, and protein substitution-model distances such as BLOSUM62 and PAM250.
- `parsimony_score` and `maximum_parsimony_tree` provide a practical maximum-parsimony workflow for aligned sequences.
- `draw_ascii` emits a quick text rendering for REPL inspection.
- `draw_unicode` emits a richer box-drawing tree rendering.
- `tree_to_dot` exports tree to Graphviz DOT format.
- `tree_to_mermaid` exports tree to Mermaid diagram format.
- `get_parent` returns the parent node of a named taxon or internal clade.
- `lowest_common_ancestor` finds the deepest shared ancestor for one or more taxa.
- `common_ancestor` is an alias for `lowest_common_ancestor`.
- `is_terminal` checks whether a node is a leaf.
- `is_parent_of` checks whether one named node is the direct parent of another.
- `get_path` returns the path from the root to a named taxon or internal clade.
- `trace` returns the nodes connecting two named taxa.
- `split` adds child nodes under a named clade.
- `is_monophyletic` checks whether a taxon set forms a clade.
- `robinson_foulds_distance` compares two trees by informative clade splits.
- `is_bifurcating` checks whether all internal nodes have exactly two children.
- `is_preterminal` checks whether a node has a single leaf child.
- `total_branch_length` sums all branch lengths in the tree.
- `depths` returns a dictionary mapping each node to its distance from the root.
- `find_clades` searches for nodes by name pattern or filters terminals.
- `ladderize` sorts branches for aesthetic display (commonly used before publication).
- `collapse_clades` removes clades below a support threshold, promoting weakly-supported internal nodes to polytomies.
- `prune` keeps only selected leaves and collapses single-child branches.
- `root_with_outgroup` and `reroot` rebuild the tree around a chosen outgroup.
- `midpoint_root` roots the tree at the midpoint of the longest path.
- `bootstrap_trees` resamples alignment columns and reconstructs replicate trees.
- `tree_consensus`, `consensus_tree`, and `bootstrap_consensus_tree` summarize tree sets and attach support values.
- `bootstrap_support` annotates a target tree with clade frequencies from bootstrap replicates.

**Usage:**
```julia
using BioToolkit

seqs = ["ACGT", "ACCT", "ACCA", "AAAA"]
labels = ["Seq1", "Seq2", "Seq3", "Seq4"]

# Map all sequences to a purely continuous relative distance matrix
matrix = distance_matrix(seqs, method=:hamming)

# Iteratively collapse nearest taxonomic branches into a Newick string
newick_tree = neighbor_joining(matrix, labels)

println(newick_tree)
# Output: (((Seq3:0.25,Seq2:0.125):0.125,Seq4:0.625),Seq1:0.0);

tree = parse_newick(newick_tree)
println(draw_ascii(tree))
```

### Tree Queries
Use `lowest_common_ancestor(tree, ["A", "B"])` to recover shared ancestry, `is_monophyletic(tree, ["A", "B"])` to validate clades, and `robinson_foulds_distance(tree1, tree2)` to compare tree topologies.

### Tree Statistics and Validation
Check tree structure with `is_bifurcating(tree)` to confirm all internal nodes have exactly two children, get leaf counts with `count_terminals(tree)`, compute total evolutionary distance with `total_branch_length(tree)`, and retrieve node depths with `depths(tree)`. Use `is_preterminal(node)` to identify nodes with a single leaf child.

### Tree Search and Organization
Find nodes matching patterns with `find_clades(tree; name_pattern=r"pattern")`, search for leaf nodes only with `find_clades(tree; terminals_only=true)`, and improve readability with `ladderize(tree)` before publishing phylograms. Ladderize sorts child branches by the number of descendant leaves (by default, ascending; use `ascending=false` to reverse).

### Visualizing
Because `neighbor_joining` returns a strictly compliant **Newick** string standard string, you can effortlessly pipe the output into `Phylo.jl`:
```julia
using Phylo

tree = parsenewick(newick_tree)
plot(tree)
```
