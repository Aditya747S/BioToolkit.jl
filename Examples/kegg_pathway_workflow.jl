using BioToolkit

pathway_text = """
ENTRY       map00010                      Pathway
NAME        Glycolysis / Gluconeogenesis
ENZYME      1.1.1.1 2.7.1.1
COMPOUND    C00031 C00022
GENE        b0001  gapA
///
"""

fixture_path = joinpath(@__DIR__, "fixtures", "kegg_pathway_expected.mmd")
mktempdir() do dir
    input_path = joinpath(dir, "pathway.keg")
    output_path = joinpath(dir, "pathway.mmd")

    write(input_path, pathway_text)
    record = BioToolkit.read_kegg_pathway(input_path)
    graph = BioToolkit.read_pathway_graph(record)

    BioToolkit.write_kegg_pathway_mermaid(output_path, graph; title="KEGG Pathway Example")

    expected = read(fixture_path, String)
    actual = read(output_path, String)
    @assert actual == expected

    println("KEGG pathway workflow")
    println("  entry: ", record.entry)
    println("  nodes: ", length(BioToolkit.pathway_nodes(graph)))
    println("  edges: ", length(BioToolkit.pathway_edges(graph)))
    println("  mermaid file: ", output_path)
end
