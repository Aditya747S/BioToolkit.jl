using Test
using BioToolkit
using BioToolkit.Restriction
using BioToolkit.Entrez
using BioToolkit.Medline
using BioToolkit.Compass
using BioToolkit.KEGG
using BioToolkit.Pathway
using BioToolkit.SCOP
using BioToolkit.CATH

@testset "Database & Restriction Engine Hardening Tests" begin

    @testset "Restriction Module" begin
        # Test Catalog & Queries
        eco = restriction_enzyme("EcoRI")
        @test eco.name == "EcoRI"
        @test eco.overhang == 4
        @test !eco.is_blunt

        sma = restriction_enzyme("SmaI")
        @test sma.is_blunt
        @test sma.overhang == 0

        blunts = blunt_enzymes()
        @test any(e -> e.name == "SmaI", blunts)
        @test !any(e -> e.name == "EcoRI", blunts)

        stickies = sticky_enzymes()
        @test any(e -> e.name == "EcoRI", stickies)

        # Isoschizomers
        # XmaI and SmaI both recognize CCCGGG
        iso = isoschizomers("SmaI")
        @test any(e -> e.name == "XmaI", iso)

        # Linear vs Circular Restriction Scanning
        # Sequence: GAATTC (EcoRI site)
        seq_linear = "CCGAATTCGG"
        sites = find_restriction_sites(seq_linear, "EcoRI"; circular=false)
        @test length(sites) == 1
        @test sites[1].position == 3
        @test sites[1].cut_position == 3

        # Circular Plasmid site spanning origin
        # Origin boundary: ATTC...GA (EcoRI site GAATTC splits across origin: last 2 bp GA at end, first 4 bp ATTC at start)
        seq_circ = "ATTCGGGGGGGA"
        sites_circ = find_restriction_sites(seq_circ, "EcoRI"; circular=true)
        @test length(sites_circ) == 1
        @test sites_circ[1].position == 11 # Site starts at position 11

        # Unique Cutters
        seq_unique = "CCGAATTCGGTT"
        uniques = find_unique_cutters(seq_unique, [eco, sma]; circular=false)
        @test length(uniques) == 1
        @test uniques[1].name == "EcoRI"

        # Digestion: Linear
        # EcoRI cuts GAATTC (cut offset 1) at pos 3 -> cut_pos = 3.
        # Sequence: CCG|AATTCGG (len 10) -> fragments: CCG (3 bp) and AATTCGG (7 bp)
        frags_linear = digest_sequence(seq_linear, ["EcoRI"]; circular=false)
        @test length(frags_linear) == 2
        @test String(frags_linear[1]) == "CCG"
        @test String(frags_linear[2]) == "AATTCGG"

        # Digestion: Circular Plasmid (1 cut gives 1 full-length linear fragment)
        frags_circ = digest_sequence(seq_linear, ["EcoRI"]; circular=true)
        @test length(frags_circ) == 1
        @test length(frags_circ[1]) == 10

        # HTML Plasmid Map Visualizer Generation
        html_map = restriction_map_html(seq_linear, sites; circular=true, title="pUC19 Plasmid Map")
        @test occursin("pUC19 Plasmid Map", html_map)
        @test occursin("EcoRI", html_map)
        @test occursin("<canvas", html_map)

        tmp_html = tempname() * ".html"
        try
            export_restriction_map_html(seq_linear, tmp_html; circular=true)
            @test isfile(tmp_html)
            @test filesize(tmp_html) > 0
        finally
            rm(tmp_html; force=true)
        end
    end

    @testset "Entrez Module" begin
        mock_search_json = """
        {
          "header": { "type": "esearch", "version": "0.8" },
          "esearchresult": {
            "count": "2",
            "retmax": "2",
            "retstart": "0",
            "idlist": ["123456", "789012"],
            "querykey": "1",
            "webenv": "NCID_1_2_3"
          }
        }
        """
        parsed_search = parse_entrez_search_response(mock_search_json)
        @test parsed_search.count == 2
        @test parsed_search.ids == ["123456", "789012"]
        @test parsed_search.webenv == "NCID_1_2_3"

        mock_post_json = """
        {
          "epostresult": {
            "querykey": "2",
            "webenv": "NCID_POST_456"
          }
        }
        """
        parsed_post = parse_entrez_post_response(mock_post_json)
        @test parsed_post.query_key == 2
        @test parsed_post.webenv == "NCID_POST_456"
    end

    @testset "Medline Module" begin
        sample_medline = """
        PMID- 35000001
        OWN - NLM
        STAT- MEDLINE
        TI  - Comprehensive Genomic Profiling in Clinical Oncology.
        AB  - Next-generation sequencing allows precise tumor profiling.
        FAU - Smith, John
        AU  - Smith J
        FAU - Doe, Jane
        AU  - Doe J
        MH  - Neoplasms/genetics
        MH  - High-Throughput Nucleotide Sequencing
        DP  - 2024 Jan 15
        AID - 10.1038/s41588-023-00001-x [doi]

        PMID- 35000002
        TI  - BioToolkit: Next-Generation Bioinformatics in Julia.
        AB  - High performance sequence analysis and database integration.
        DP  - 2025
        """
        records = parse_medline_text(sample_medline)
        @test length(records) == 2
        @test records[1].pmid == "35000001"
        @test records[1].title == "Comprehensive Genomic Profiling in Clinical Oncology."
        @test records[1].year == 2024
        @test records[1].doi == "10.1038/s41588-023-00001-x"
        @test length(records[1].authors) >= 2
        @test length(records[1].mesh_terms) == 2

        @test records[2].pmid == "35000002"
        @test records[2].year == 2025

        xml_sample = """
        <PubmedArticleSet>
          <PubmedArticle>
            <MedlineCitation>
              <PMID>999999</PMID>
            </MedlineCitation>
            <Article>
              <ArticleTitle>XML Parsing Test Title</ArticleTitle>
              <AbstractText>Testing Medline XML parser fidelity.</AbstractText>
              <Journal><Title>Journal of Julia Bio</Title></Journal>
              <JournalIssue><PubDate><Year>2026</Year></PubDate></JournalIssue>
              <AuthorList>
                <Author><LastName>Turing</LastName><ForeName>Alan</ForeName></Author>
              </AuthorList>
            </Article>
          </PubmedArticle>
        </PubmedArticleSet>
        """
        xml_recs = parse_medline_xml(xml_sample)
        @test length(xml_recs) == 1
        @test xml_recs[1].pmid == "999999"
        @test xml_recs[1].title == "XML Parsing Test Title"
        @test xml_recs[1].authors == ["Alan Turing"]
    end

    @testset "Compass Native Fallback" begin
        s1 = BioSequence{DNAAlphabet}("ATGCGATCG")
        s2 = BioSequence{DNAAlphabet}("ATGCGATAG")
        
        align_needle = run_needle(s1, s2)
        @test occursin("Alignment", align_needle)

        align_water = run_water(s1, s2)
        @test occursin("Alignment", align_water)
    end

    @testset "KEGG & Pathway Graph & KGML" begin
        kegg_text = """
        ENTRY       map00010                    Pathway
        NAME        Glycolysis / Gluconeogenesis
        ENZYME      1.1.1.1  2.7.1.1  5.3.1.9
        COMPOUND    C00031  C00022  C00024
        GENE        1234  Alpha-enolase
                    5678  Beta-enolase
                    ENO1  Gamma-enolase
        ///
        """
        rec = read_kegg_pathway(IOBuffer(kegg_text))
        @test rec.entry == "map00010"
        @test rec.title == "Glycolysis / Gluconeogenesis"
        @test length(rec.enzymes) == 3
        @test length(rec.compounds) == 3

        graph = read_pathway_graph(rec)
        @test length(graph.nodes) > 5
        @test length(pathway_genes(graph)) == 3
        @test length(pathway_enzymes(graph)) == 3
        @test length(pathway_compounds(graph)) == 3

        sub = pathway_subgraph(graph, ["1234", "C00031"])
        @test length(sub.nodes) == 2

        mermaid_str = kegg_pathway_mermaid(graph)
        @test startswith(mermaid_str, "flowchart LR")
        @test occursin("ENO1", mermaid_str)

        # KGML XML Parsing
        kgml_xml = """
        <pathway name="path:map00010" org="map" number="00010" title="Glycolysis">
          <entry id="1" name="path:map00010" type="map"/>
          <entry id="2" name="ec:1.1.1.1" type="enzyme"/>
          <entry id="3" name="cpd:C00031" type="compound"/>
          <relation entry1="2" entry2="3" type="EC-rel"/>
        </pathway>
        """
        kgml_graph = read_kgml(kgml_xml)
        @test length(kgml_graph.nodes) == 3
        @test length(kgml_graph.edges) == 1
        @test kgml_graph.edges[1].relation == "EC-rel"

        # HTML Pathway Visualizer
        html_pw = pathway_to_html(kgml_graph; title="Glycolysis Pathway")
        @test occursin("Glycolysis Pathway", html_pw)
        @test occursin("vis.Network", html_pw)
    end

    @testset "SCOP & CATH Module" begin
        scop_line = "d1dlwa_ 1dlw A:1-112 a.1.1.1 SCOP Domain Description Text"
        scop_rec = parse_scop_record(scop_line)
        @test scop_rec.sid == "d1dlwa_"
        @test scop_rec.pdb_id == "1dlw"
        @test scop_rec.class_id == "a.1.1.1"
        @test scop_rec.hierarchy == ["a", "1", "1", "1"]

        found_scop = find_scop_by_pdb([scop_rec], "1dlw")
        @test length(found_scop) == 1

        cath_line = "1oaiA00 1oai A 1 10 8 10 CATH Domain Description"
        cath_rec = parse_cath_record(cath_line)
        @test cath_rec.domain == "1oaiA00"
        @test cath_rec.pdb_id == "1oai"
        @test cath_rec.hierarchy == ["1", "10", "8", "10"]

        found_cath = find_cath_by_pdb([cath_rec], "1OAI")
        @test length(found_cath) == 1
    end

end
