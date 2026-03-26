# databases.jl

## Purpose
This file is a collection of sequence-oriented database helpers. It combines restriction-enzyme lookup and digestion logic with Entrez, MEDLINE, COMPASS, KEGG, SCOP, and CATH record parsers so the package can work with external biological databases from one place.

## Module layout
- `Restriction` provides restriction-enzyme catalog and digestion support.
- `Entrez` wraps NCBI E-utilities for search, fetch, post, and link operations.
- `Medline` parses PubMed-style text and XML records into structured objects.
- `Compass` wraps common EMBOSS command-line tools such as Needle, Water, Transeq, Revseq, Compseq, and Seqret.
- `KEGG` parses KEGG flat-file records and pathway/enzyme records.
- `Pathway` converts KEGG pathway data into a graph and Mermaid diagram output.
- `SCOP` and `CATH` parse structural classification records.

## Main types
- `RestrictionEnzyme` and `RestrictionSite` model recognition patterns and cut positions.
- `EntrezSearchResult` and `EntrezPostResult` model NCBI search/post responses.
- `MedlineRecord` stores article metadata, abstract text, author lists, MeSH terms, and raw fields.
- `KEGGRecord`, `KEGGPathwayRecord`, and `KEGGEnzymeRecord` store parsed KEGG flat-file data.
- `PathwayNode`, `PathwayEdge`, and `PathwayGraph` model a pathway graph derived from KEGG data.
- `SCOPRecord` and `CATHRecord` store structural classification rows.

## Public functions
- Restriction helpers: `restriction_enzymes`, `restriction_enzyme`, `restriction_enzyme_names`, `restriction_catalog`, `restriction_sites`, `restriction_digest_map`, `find_restriction_sites`, and `digest_sequence`.
- Entrez helpers: `entrez_search`, `entrez_search_ids`, `entrez_search_count`, `entrez_fetch`, `entrez_fetch_fasta`, `entrez_fetch_genbank`, `entrez_summary`, `entrez_post`, `entrez_post_ids`, `entrez_post_webenv`, `entrez_post_query_key`, `entrez_link`, `entrez_elink`, `entrez_link_ids`, `entrez_link_linksets`, `entrez_link_records`, `entrez_pubmed_search`, `entrez_nuccore_fetch`, `entrez_nucleotide_search`, `entrez_protein_search`, `entrez_gene_search`, `entrez_taxonomy_search`, `entrez_genome_search`, `entrez_genome_fetch`, `entrez_pubmed_fetch`.
- MEDLINE helpers: `parse_medline`, `parse_medline_xml`, `parse_medline_text`, and `read_medline`.
- KEGG helpers: `read_kegg_record`, `read_kegg_pathway`, `read_kegg_enzyme`, `kegg_field`, `kegg_entries`, and `kegg_entry_id`.
- Pathway helpers: `read_pathway_graph`, `pathway_nodes`, `pathway_edges`, `kegg_pathway_mermaid`, and `write_kegg_pathway_mermaid`.
- Structural classification helpers: `parse_scop_record`, `read_scop_records`, `parse_cath_record`, and `read_cath_records`.

## How it is used
For restriction analysis, users can fetch the catalog with `restriction_catalog`, look up one enzyme with `restriction_enzyme`, and scan a DNA string with `find_restriction_sites` or `digest_sequence`.

For literature and database workflows, `entrez_search` and `entrez_fetch` provide the basic NCBI access layer, `parse_medline` converts PubMed records into structured metadata, and the KEGG/Pathway helpers turn flat files into graph-like structures that are easy to render.

## Implementation notes
- Restriction patterns support IUPAC ambiguity codes by compiling regexes from the recognition sequence.
- Entrez helpers build URLs explicitly and download responses as text or JSON.
- MEDLINE parsing supports both text and XML representations.
- Pathway rendering sanitizes identifiers before emitting Mermaid syntax.

## Why it matters
This file is the package's database integration layer. It gives BioToolkit practical access to sequence databases, literature metadata, and pathway catalogs without forcing callers to write their own parsers or HTTP wrappers.
