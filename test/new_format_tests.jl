using BioToolkit
using Test

@testset "EMBL & SwissProt Parsing" begin
    embl_data = """ID   XYZ123; SV 1; linear; DNA; ; ; 30 BP.
XX
AC   AC123456;
XX
DE   Test EMBL record.
XX
FT   source          1..30
FT                   /organism="Test organism"
FT   gene            10..20
FT                   /gene="test_gene"
XX
SQ   Sequence 30 BP; 10 A; 10 C; 10 G; 0 T; 0 OTHER;
     aaaaaaaaaa cccccccccc gggggggggg                                30
//
"""
    tmp_embl, io = mktemp()
    write(io, embl_data)
    close(io)
    
    try
        records = BioToolkit.read_embl(tmp_embl)
        @test length(records) == 1
        @test records[1].identifier == "XYZ123;"
        @test records[1].sequence == "aaaaaaaaaaccccccccccgggggggggg"
        @test records[1].annotations[:accession] == "AC123456;"
        @test length(records[1].features) == 2
        @test records[1].features[1].feature_type == "source"
        @test feature_start(records[1].features[1]) == 1
        @test feature_stop(records[1].features[1]) == 30
        @test records[1].features[1].qualifiers["organism"] == ["Test organism"]
        
        # SwissProt uses the same parser engine
        records_sp = BioToolkit.read_swissprot(tmp_embl)
        @test length(records_sp) == 1
    finally
        rm(tmp_embl)
    end
end

@testset "Database-style modules" begin
    @test isdefined(BioToolkit, :Restriction)
    @test isdefined(BioToolkit, :Entrez)
    @test isdefined(BioToolkit, :KEGG)
    @test isdefined(BioToolkit, :Pathway)
    @test isdefined(BioToolkit, :SCOP)
    @test isdefined(BioToolkit, :CATH)
    @test isdefined(BioToolkit, :entrez_link)
    @test isdefined(BioToolkit, :entrez_link_ids)
    @test isdefined(BioToolkit, :entrez_post)
    @test isdefined(BioToolkit, :entrez_elink)
    @test isdefined(BioToolkit, :entrez_genome_search)
    @test isdefined(BioToolkit, :entrez_genome_fetch)

    @test length(BioToolkit.restriction_enzyme_names()) >= 60
    enzyme = BioToolkit.restriction_enzyme("EcoRI")
    hits = BioToolkit.find_restriction_sites("TTTGAATTCTTT", enzyme)
    @test length(hits) == 1
    @test hits[1].cut_position == 5
    @test BioToolkit.restriction_sites("TTTGAATTCTTT", enzyme) == hits
    @test BioToolkit.digest_sequence("TTTGAATTCTTT", "EcoRI") == ["TTTG", "AATTCTTT"]
    @test occursin("EcoRI", sprint(show, enzyme))
    @test length(BioToolkit.restriction_digest_map("TTTGAATTCGAATTC", ["EcoRI"])["EcoRI"]) == 2
    @test occursin("RestrictionSite", sprint(show, hits[1]))

    entrez_payload = """{"header":{"type":"esearch","version":"0.3"},"esearchresult":{"count":"2","retmax":"2","retstart":"0","idlist":["123","456"],"translationset":[],"translationstack":[],"querytranslation":"BRCA1"}}"""
    entrez_result = BioToolkit.parse_entrez_search_response(entrez_payload)
    @test entrez_result.count == 2
    @test entrez_result.ids == ["123", "456"]
    @test occursin("count=2", sprint(show, entrez_result))
    @test BioToolkit.entrez_search_count isa Function
    @test BioToolkit.entrez_fetch_fasta isa Function

    entrez_post_payload = """{"header":{"type":"epost","version":"0.3"},"epostresult":{"querykey":"7","webenv":"NCBI-WEBENV-123","idlist":["321","654"]}}"""
    entrez_post_result = BioToolkit.parse_entrez_post_response(entrez_post_payload)
    @test entrez_post_result.query_key == 7
    @test entrez_post_result.webenv == "NCBI-WEBENV-123"
    @test BioToolkit.entrez_post_ids(entrez_post_result) == ["321", "654"]
    @test occursin("EntrezPostResult", sprint(show, entrez_post_result))

    pathway_text = """
ENTRY       map00010                      Pathway
NAME        Glycolysis / Gluconeogenesis
ENZYME      1.1.1.1 2.7.1.1
COMPOUND    C00031 C00022
GENE        b0001  gene1
///
"""
    pathway_file, pathway_io = mktemp()
    write(pathway_io, pathway_text)
    close(pathway_io)
    try
        pathway_record = BioToolkit.read_kegg_pathway(pathway_file)
        @test pathway_record.entry == "map00010"
        @test occursin("Glycolysis", pathway_record.title)
        @test "1.1.1.1" in pathway_record.enzymes
        pathway_graph = BioToolkit.read_pathway_graph(pathway_record)
        @test length(BioToolkit.pathway_nodes(pathway_graph)) >= 2
        @test length(BioToolkit.pathway_edges(pathway_graph)) >= length(pathway_record.enzymes)
        @test occursin("PathwayGraph", sprint(show, pathway_graph))
        mermaid = BioToolkit.kegg_pathway_mermaid(pathway_graph)
        @test occursin("flowchart LR", mermaid)
        @test occursin("map00010", mermaid)
        buffer = IOBuffer()
        BioToolkit.write_kegg_pathway_mermaid(buffer, pathway_graph)
        @test occursin("flowchart LR", String(take!(buffer)))
        fixture = read(joinpath(@__DIR__, "..", "Examples", "fixtures", "kegg_pathway_expected.mmd"), String)
        @test BioToolkit.kegg_pathway_mermaid(pathway_graph; title="KEGG Pathway Example") == fixture
    finally
        rm(pathway_file)
    end

    enzyme_text = """
ENTRY       EC 1.1.1.1
NAME        alcohol dehydrogenase; ADH
REACTION    ethanol + NAD+ = acetaldehyde + NADH
PATHWAY     map00010  Glycolysis / Gluconeogenesis
///
"""
    enzyme_file, enzyme_io = mktemp()
    write(enzyme_io, enzyme_text)
    close(enzyme_io)
    try
        enzyme_record = BioToolkit.read_kegg_enzyme(enzyme_file)
        @test startswith(enzyme_record.entry, "EC")
        @test occursin("alcohol dehydrogenase", enzyme_record.names[1])
        @test !isempty(enzyme_record.pathways)
    finally
        rm(enzyme_file)
    end

    scop_record = BioToolkit.parse_scop_record("d1tq3a_ 1tq3 A:1-100 a.1.1.1 Protein description")
    @test scop_record.sid == "d1tq3a_"
    @test scop_record.class_id == "a.1.1.1"
    @test scop_record.hierarchy == ["a", "1", "1", "1"]
    @test occursin("SCOPRecord", sprint(show, scop_record))

    cath_record = BioToolkit.parse_cath_record("1abcA00 1abc A 1 10 20 30 Example domain")
    @test cath_record.domain == "1abcA00"
    @test cath_record.homologous_superfamily == "30"
    @test cath_record.hierarchy == ["1", "10", "20", "30"]
    @test occursin("CATHRecord", sprint(show, cath_record))
end

@testset "MAF Alignment Parsing" begin
    maf_data = """##maf version=1
a score=100
s target.chr1 10 5 + 100 ACGTC
s query.chr2 20 5 - 200 ACG-C
"""
    tmp_maf, io = mktemp()
    write(io, maf_data)
    close(io)
    
    try
        msa = BioToolkit.read_alignment(tmp_maf, "maf")
        @test length(msa) == 2
        @test msa[1].identifier == "target.chr1"
        @test msa[1].annotations[:start] == 10
        @test msa[2].annotations[:strand] == "-"
        @test msa[2].sequence == "ACG-C"
    finally
        rm(tmp_maf)
    end
end

@testset "GCG Alignment Parsing" begin
    gcg_data = """!!NA_SEQUENCE 1.0
Test GCG sequence.
  Check: 1234  Length: 20 ..
1  ACGT ACGT GTGT GTGT
"""
    tmp_gcg, io = mktemp()
    write(io, gcg_data)
    close(io)
    
    try
        msa = BioToolkit.read_alignment(tmp_gcg, "gcg")
        @test length(msa) == 1
        @test msa[1].identifier == "Check: 1234  Length: 20"
        @test msa[1].sequence == "ACGTACGTGTGTGTGT"
    finally
        rm(tmp_gcg)
    end
end
