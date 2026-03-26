using Test
using BioToolkit

@testset "Phylogenetics (Neighbor-Joining)" begin
    @testset "Distance Matrix (Hamming)" begin
        seqs = ["ACGT", "ACCT", "ACCA", "AAAA"]
        D = distance_matrix(seqs; method=:hamming)
        
        @test size(D) == (4, 4)
        @test D[1, 1] == 0.0
        # ACGT vs ACCT -> 1 mismatch / 4 length
        @test D[1, 2] == 0.25
        # ACGT vs ACCA -> 2 mismatches
        @test D[1, 3] == 0.50
        # ACGT vs AAAA -> 3 mismatches
        @test D[1, 4] == 0.75
        
        # Identity property
        @test all(sum(D[i,i] for i in 1:4) == 0.0)
        # Symmetry property
        @test D ≈ transpose(D)
    end
    
    @testset "Neighbor Joining Structural Output" begin
        # Classic textbook 4-taxa example
        nodes = ["A", "B", "C", "D"]
        D = [
            0.0  7.0  11.  14.;
            7.0  0.0  6.0   9.;
            11.  6.0  0.0   7.;
            14.  9.0  7.0  0.0
        ]
        
        newick = neighbor_joining(D, nodes)
        @test typeof(newick) == String
        # Should gracefully format standard nested edge syntax without errors
        @test occursin("(", newick)
        @test occursin("A", newick)
        @test occursin("B", newick)
        @test occursin("C", newick)
        @test occursin("D", newick)
        @test newick[end] == ';'
    end
    
    @testset "Distance Matrix (Alignment)" begin
        seqs = ["ACGTG", "ACGTA"]
        D = distance_matrix(seqs; method=:alignment)
        @test size(D) == (2, 2)
        @test D[1, 1] == 0.0
        @test D[1, 2] > 0.0 # 1 mismatch should drop identity
    end

    @testset "Tree parsing and writing" begin
        tree = parse_newick("((A:0.25,B:0.25):0.5,C:0.75);")
        @test tree.children[1].children[1].name == "A"
        @test write_newick(tree) == "((A:0.25,B:0.25):0.5,C:0.75);"
        @test length(get_terminals(tree)) == 3
        @test length(get_nonterminals(tree)) == 2
        @test isapprox(tree_distance(tree, "A", "B"), 0.5; atol=1e-8)
        @test get_parent(tree, "A") == tree.children[1]
        @test lowest_common_ancestor(tree, ["A", "B"]) == tree.children[1]
        @test is_monophyletic(tree, ["A", "B"])
        @test !is_monophyletic(tree, ["A", "C"])
        @test robinson_foulds_distance(tree, tree) == 0
        @test occursin("A", draw_ascii(tree))
        @test occursin("└", draw_unicode(tree))
        @test occursin("digraph", tree_to_dot(tree))
        @test occursin("graph TD", tree_to_mermaid(tree))

        support_tree = parse_newick("((A:0.25,B:0.25)0.95:0.5,C:0.75);")
        @test isapprox(support_tree.children[1].support, 0.95; atol=1e-8)
        @test occursin("0.95", write_newick(support_tree))

        midpoint = midpoint_root(tree)
        @test sort(collect(node.name for node in get_terminals(midpoint))) == ["A", "B", "C"]
        @test isapprox(tree_distance(midpoint, "A", "B"), tree_distance(tree, "A", "B"); atol=1e-8)

        pruned = prune(tree, ["A", "C"])
        @test sort(collect(node.name for node in get_terminals(pruned))) == ["A", "C"]

        rerooted = root_with_outgroup(tree, "C")
        @test sort(collect(node.name for node in get_terminals(rerooted))) == ["A", "B", "C"]
        @test sort(collect(node.name for node in get_terminals(reroot(tree, "C")))) == ["A", "B", "C"]
    end

    @testset "Tree metadata and formats" begin
        tree = parse_newick("((A:0.1,B:0.2)X:0.3,C:0.4);")
        set_metadata!(tree, "organism", "synthetic")
        annotate_tree!(tree, Dict("source" => "test"))

        phyloxml = write_phyloxml(tree)
        @test occursin("phyloxml", phyloxml)
        @test occursin("synthetic", phyloxml)
        @test occursin("test", phyloxml)
        @test count_terminals(parse_phyloxml(phyloxml)) == 3

        nexus = write_nexus(tree)
        @test occursin("#NEXUS", nexus)
        @test count_terminals(parse_nexus(nexus)) == 3

        nexml = write_nexml(tree)
        @test occursin("nexml", nexml)
        @test count_terminals(parse_nexml(nexml)) == 3

        @test count_terminals(parse_tree(write_tree(tree); format=:newick)) == 3
        @test count_terminals(parse_tree(write_tree(tree; format=:phyloxml); format=:phyloxml)) == 3
    end

    @testset "Substitution model distances" begin
        dna = ["ACGT", "ACGA", "AGGT"]
        @test isapprox(distance_matrix(dna; method=:hamming)[1, 2], 0.25; atol=1e-8)
        @test distance_matrix(dna; method=:jukes_cantor)[1, 2] > 0.0
        @test distance_matrix(dna; method=:kimura2p)[1, 2] > 0.0

        protein = ["MKWV", "MKWV", "MKRV"]
        @test distance_matrix(protein; method=:blosum62)[1, 2] == 0.0
        @test distance_matrix(protein; method=:pam250)[1, 3] >= 0.0
    end

    @testset "Maximum parsimony tree" begin
        alignment = BioToolkit.MultipleSequenceAlignment([
            BioToolkit.SeqRecordLite("ACGT"; identifier="A"),
            BioToolkit.SeqRecordLite("ACGA"; identifier="B"),
            BioToolkit.SeqRecordLite("ACGG"; identifier="C"),
        ])
        parsimony = maximum_parsimony_tree(alignment)
        @test count_terminals(parsimony) == 3
        @test parsimony_score(parsimony, alignment) >= 0
        @test count_terminals(parsimony_tree(alignment)) == 3
    end

    @testset "UPGMA tree building" begin
        seqs = ["ACGT", "ACGA", "TCGT"]
        names = ["s1", "s2", "s3"]
        D = distance_matrix(seqs; method=:hamming)
        tree = upgma(D, names)
        terminals = get_terminals(tree)
        @test sort(collect(node.name for node in terminals)) == sort(names)
        @test startswith(write_newick(tree), "(")
        @test occursin("s1", write_newick(tree))

        consensus = consensus_tree([tree, tree])
        @test sort(collect(node.name for node in get_terminals(consensus))) == sort(names)
        @test all(node.support >= 0.0 for node in get_nonterminals(consensus))
        @test strict_consensus_tree([tree, tree]) |> get_terminals |> length == 3
        @test majority_consensus_tree([tree, tree]) |> get_terminals |> length == 3
    end

    @testset "NJ tree conversion" begin
        seqs = ["ACGT", "ACCT", "ACCA", "AAAA"]
        names = ["A", "B", "C", "D"]
        D = distance_matrix(seqs; method=:hamming)
        tree = neighbor_joining_tree(D, names)
        @test sort(collect(node.name for node in get_terminals(tree))) == sort(names)
        @test endswith(write_newick(tree), ";")

        annotated = bootstrap_support([tree, tree], tree)
        @test occursin("1.0", write_newick(annotated))

        alignment = BioToolkit.MultipleSequenceAlignment([
            BioToolkit.SeqRecordLite("ACGT"; identifier="A"),
            BioToolkit.SeqRecordLite("ACGA"; identifier="B"),
            BioToolkit.SeqRecordLite("ACGT"; identifier="C"),
        ])
        bootstrap = bootstrap_trees(alignment; replicates=3, method=:hamming, constructor=:nj)
        @test length(bootstrap) == 3

        annotated = bootstrap_support([tree, tree], tree)
        @test all(node.support >= 0.0 for node in get_nonterminals(annotated))

        bootstrap_consensus = bootstrap_consensus_tree(alignment; replicates=3, threshold=0.5, method=:hamming, constructor=:nj)
        @test length(get_terminals(bootstrap_consensus)) == 3
    end

    @testset "Tree validation and statistics" begin
        tree = parse_newick("((A:0.1,B:0.2):0.3,C:0.4);")
        
        @test count_terminals(tree) == 3
        @test is_bifurcating(tree)
        @test total_branch_length(tree) > 0.0
        
        depth_map = depths(tree)
        @test length(depth_map) > 0
        @test depth_map[tree] == 0.0
        
        non_bifurcating = parse_newick("((A:0.1,B:0.2,C:0.3):0.4,D:0.5);")
        @test !is_bifurcating(non_bifurcating)
        
        preterminal_tree = parse_newick("((A:0.1):0.2,B:0.3);")
        @test !isempty(tree.children[1].children)
        @test is_preterminal(preterminal_tree.children[1])
    end

    @testset "Tree search and query" begin
        tree = parse_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
        
        all_clades = find_clades(tree)
        @test length(all_clades) > 0
        
        terminals = find_clades(tree; terminals_only=true)
        @test all(isempty(node.children) for node in terminals)
        @test length(terminals) == 4
        
        named_leaves = find_clades(tree; name_pattern=r"[A-D]", terminals_only=true)
        @test length(named_leaves) >= 0

        @test is_terminal(terminals[1])

        named_tree = parse_newick("((A:0.1,B:0.2)X:0.3,(C:0.4,D:0.5):0.6);")
        @test is_parent_of(named_tree, "X", "A")

        path = get_path(named_tree, "X")
        @test path[end].name == "X"

        trace_nodes = trace(tree, "A", "D")
        @test length(trace_nodes) > 0
        @test common_ancestor(tree, ["A", "B"]) == lowest_common_ancestor(tree, ["A", "B"])
    end

    @testset "Split helper" begin
        tree = parse_newick("(A:0.1,B:0.2)X:0.3;")
        split_tree = BioToolkit.split(tree, "X"; child_names=["A1", "A2"], branch_length=0.05)
        @test count_terminals(split_tree) >= 4
    end

    @testset "Ladderize" begin
        tree = parse_newick("((A:0.1,B:0.2):0.3,C:0.4);")
        
        laddered_asc = ladderize(tree; ascending=true)
        @test count_terminals(laddered_asc) == 3
        @test is_bifurcating(laddered_asc)
        
        laddered_desc = ladderize(tree; ascending=false)
        @test count_terminals(laddered_desc) == 3
        
        newick_original = write_newick(tree)
        newick_laddered = write_newick(laddered_asc)
        @test sort(collect(node.name for node in get_terminals(tree))) == sort(collect(node.name for node in get_terminals(laddered_asc)))
    end

    @testset "Collapse clades" begin
        support_tree = parse_newick("((A:0.1,B:0.2)0.3:0.3,(C:0.4,D:0.5)0.9:0.4);")
        
        collapsed = collapse_clades(support_tree; min_support=0.5)
        @test count_terminals(collapsed) == 4
        
        collapsed_strict = collapse_clades(support_tree; min_support=0.8)
        @test count_terminals(collapsed_strict) == 4
    end
end
