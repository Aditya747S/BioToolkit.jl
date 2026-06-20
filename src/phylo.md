# `phylo.jl` - Phylogenetics

## Overview

`phylo.jl` implements tree containers, Newick/PhyloXML/Nexus/NeXML parsing and writing, tree traversal, rooting/pruning, distance matrices, neighbor joining, UPGMA, maximum parsimony/likelihood helpers, consensus trees, bootstrap support, and text graph exports.

### Purpose

This page is a hand-authored reference for `phylo.jl`, grouped around the exported user workflows. Internal helper functions are omitted unless the module exposes them as part of the public API.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Source-matched APIs** | Entries use names implemented in the corresponding `.jl` file. |
| **Workflow sections** | Related functions are documented together so the analysis path is clear. |
| **Concrete result descriptions** | Structs and result types are described by their role in downstream workflows. |
| **Julia-first data flow** | APIs compose through normal Julia arrays, tables, graphs, and BioToolkit objects. |
| **Export-oriented coverage** | The page covers the public functions users are likely to import from BioToolkit. |

---

## 1. Tree Types and Traversal

Tree helpers inspect topology, labels, branch lengths, and ancestry relationships.

| API | Description |
|---|---|
| `PhyloTree` | Tree node with name, branch length, support, metadata, and children. |
| `isleaf` | Tests whether a tree node has no children. |
| `coordinates` | Computes plotting/layout coordinates for a tree. |
| `get_terminals` | Returns terminal leaf nodes. |
| `get_nonterminals` | Returns internal nodes. |
| `get_parent` | Finds the parent of a named node. |
| `lowest_common_ancestor` | Finds the LCA of two or more named leaves. |
| `common_ancestor` | Alias for common ancestor lookup. |
| `is_terminal` | Tests terminal status. |
| `is_parent_of` | Tests whether one named node is ancestral to another. |
| `get_path` | Returns path from root to a named node. |
| `trace` | Returns path between two named nodes. |
| `find_clades` | Finds clades by name pattern or terminal status. |
| `count_terminals` | Counts terminal leaves. |
| `depths` | Computes root-to-node branch-length depths. |
| `total_branch_length` | Sums branch lengths. |
| `tree_distance` | Computes path length between two leaves. |

## 2. Models and Tree Building

Distance, parsimony, and likelihood helpers construct or score trees from alignments.

| API | Description |
|---|---|
| `JC69` | Jukes-Cantor substitution model. |
| `K80` | Kimura two-parameter substitution model. |
| `HKY85` | HKY85 substitution model. |
| `transition_probability` | Returns transition matrix/probabilities for a model and branch length. |
| `distance_matrix` | Computes sequence distance matrices. |
| `neighbor_joining` | Neighbor-joining tree builder. |
| `neighbor_joining_tree` | Neighbor-joining constructor from distance matrix and names. |
| `upgma` | UPGMA clustering tree builder. |
| `felsenstein_likelihood` | Computes pruning-algorithm likelihood for a tree/alignment/model. |
| `optimize_branch_lengths!` | Updates branch lengths to improve likelihood. |
| `maximum_likelihood_tree` | Builds/refines a maximum-likelihood tree. |
| `parsimony_score` | Computes Fitch parsimony score. |
| `maximum_parsimony_tree` | Finds a maximum-parsimony tree for small alignments. |
| `parsimony_tree` | Alias for maximum parsimony tree construction. |

## 3. Tree Editing and Consensus

Editing functions change topology or summarize replicate trees.

| API | Description |
|---|---|
| `midpoint_root` | Roots a tree at the midpoint of the longest leaf-to-leaf path. |
| `root_with_outgroup` | Roots a tree using a named outgroup. |
| `reroot` | Alias for outgroup rooting. |
| `prune` | Keeps selected leaves and removes others. |
| `split` | Adds child branches below a named parent. |
| `is_monophyletic` | Tests whether a set of taxa forms a clade. |
| `robinson_foulds_distance` | Computes RF distance between two trees. |
| `is_bifurcating` | Tests whether every internal node has two children. |
| `is_preterminal` | Tests whether a node directly parents only terminal nodes. |
| `ladderize` | Orders child clades by size. |
| `collapse_clades` | Collapses low-support clades. |
| `bootstrap_trees` | Builds bootstrap replicate trees. |
| `tree_consensus` | Creates a consensus tree from many trees. |
| `consensus_tree` | Consensus tree alias. |
| `strict_consensus_tree` | Consensus requiring all replicates. |
| `majority_consensus_tree` | Majority-rule consensus. |
| `bootstrap_consensus_tree` | Runs bootstrap and consensus construction. |
| `bootstrap_support` | Maps replicate support onto a target tree. |
| `set_metadata!` | Sets one metadata field. |
| `annotate_tree!` | Adds metadata fields to a tree. |

## 4. I/O and Text Export

Multiple tree text formats and simple diagram exports are supported.

| API | Description |
|---|---|
| `parse_newick` | Parses Newick text. |
| `write_newick` | Writes Newick text. |
| `parse_tree` | Parses Newick, PhyloXML, Nexus, or NeXML by format keyword. |
| `write_tree` | Writes supported tree formats by format keyword. |
| `write_phyloxml` | Writes PhyloXML text. |
| `parse_phyloxml` | Parses PhyloXML text. |
| `write_nexus` | Writes Nexus tree text. |
| `parse_nexus` | Parses Nexus tree text. |
| `write_nexml` | Writes NeXML text. |
| `parse_nexml` | Parses NeXML text. |
| `draw_ascii` | Creates ASCII tree drawing. |
| `draw_unicode` | Creates Unicode tree drawing. |
| `tree_to_dot` | Exports Graphviz DOT. |
| `tree_to_mermaid` | Exports Mermaid graph syntax. |

---

## Complete Usage Example

```julia
using BioToolkit

tree = parse_newick("((A:0.1,B:0.1):0.2,C:0.3);")
leaves = get_terminals(tree)
rooted = midpoint_root(tree)
println(write_newick(rooted))
println(tree_to_mermaid(rooted))
```

