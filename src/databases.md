# `databases.jl` - Database and Pathway Parsers

## Overview

`databases.jl` collects parsers and lightweight clients for NCBI Entrez, MEDLINE, KEGG, pathway graphs, SCOP/CATH classification records, and restriction-enzyme catalogs.

### Purpose

This page documents the user-facing types and workflows implemented in `databases.jl`. The emphasis is on public APIs exported through BioToolkit and the biological analysis step each one supports.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| **Biological workflow names** | APIs are named after common analysis tasks so code reads like a methods section. |
| **Typed summaries where useful** | Reusable outputs have named structs instead of anonymous dictionaries. |
| **Matrix/table interoperability** | Functions accept standard Julia matrices, vectors, and table-like records for easy integration. |
| **Deterministic defaults** | Randomized or approximate methods expose seeds, thresholds, or iteration controls where relevant. |
| **Downstream-ready results** | Outputs can feed plotting, enrichment, annotation, or reporting modules. |

---

## 1. NCBI Entrez

Entrez helpers wrap search, fetch, summary, post/history, and link workflows.

| API | Description |
|---|---|
| `EntrezSearchResult` | Parsed ESearch result with ids, count, and query metadata. |
| `EntrezPostResult` | Parsed EPost result containing WebEnv/query-key information. |
| `entrez_search` | Runs an Entrez search request. |
| `entrez_search_ids` | Returns only ids from an Entrez search. |
| `entrez_search_count` | Returns only hit count. |
| `entrez_fetch` | Fetches raw records from an Entrez database. |
| `entrez_fetch_fasta` | Fetches FASTA records. |
| `entrez_fetch_genbank` | Fetches GenBank records. |
| `entrez_summary` | Retrieves ESummary metadata. |
| `entrez_post` | Posts ids to Entrez history. |
| `entrez_link` | Finds linked records between Entrez databases. |
| `entrez_pubmed_search` | Convenience PubMed search. |
| `entrez_pubmed_fetch` | Convenience PubMed fetch. |
| `entrez_nuccore_fetch` | Convenience nucleotide fetch. |
| `entrez_gene_search` | Convenience gene search. |
| `entrez_taxonomy_search` | Convenience taxonomy search. |

## 2. MEDLINE and KEGG

Flat-file parsers turn biomedical records and KEGG entries into typed records.

| API | Description |
|---|---|
| `MedlineRecord` | Parsed MEDLINE citation fields. |
| `parse_medline` | Parses MEDLINE text records. |
| `parse_medline_xml` | Parses MEDLINE XML records. |
| `parse_medline_text` | Parses plain MEDLINE text. |
| `read_medline` | Reads MEDLINE records from path or stream. |
| `KEGGRecord` | Generic KEGG flat-file record. |
| `KEGGPathwayRecord` | KEGG pathway-specific record. |
| `KEGGEnzymeRecord` | KEGG enzyme-specific record. |
| `read_kegg_record` | Reads a generic KEGG record. |
| `read_kegg_pathway` | Reads a KEGG pathway record. |
| `read_kegg_enzyme` | Reads a KEGG enzyme record. |
| `kegg_field` | Retrieves a named KEGG field. |
| `kegg_entries` | Splits multi-entry KEGG text. |
| `kegg_entry_id` | Extracts an entry identifier. |

## 3. Pathways, Domains, and Restriction Enzymes

These helpers build graph exports and parse classification/catalog records.

| API | Description |
|---|---|
| `PathwayNode` | Node in a pathway graph. |
| `PathwayEdge` | Edge in a pathway graph. |
| `PathwayGraph` | Pathway graph container. |
| `read_pathway_graph` | Parses pathway graph data. |
| `pathway_nodes` | Returns graph nodes. |
| `pathway_edges` | Returns graph edges. |
| `kegg_pathway_mermaid` | Exports KEGG pathway graph as Mermaid. |
| `write_kegg_pathway_mermaid` | Writes Mermaid pathway text. |
| `SCOPRecord` | SCOP classification record. |
| `read_scop_records` | Reads SCOP records. |
| `parse_scop_record` | Parses one SCOP line/record. |
| `CATHRecord` | CATH classification record. |
| `read_cath_records` | Reads CATH records. |
| `parse_cath_record` | Parses one CATH line/record. |
| `RestrictionEnzyme` | Restriction enzyme definition. |
| `RestrictionSite` | Restriction-site hit. |
| `restriction_enzymes` | Returns catalog enzymes. |
| `restriction_enzyme` | Looks up one enzyme. |
| `restriction_sites` | Finds enzyme cut sites. |
| `digest_sequence` | Simulates restriction digest fragments. |

---

## Complete Usage Example

```julia
using BioToolkit

ids = entrez_search_ids("nuccore", "BRCA1[Gene] AND human[orgn]")
fasta = entrez_fetch_fasta("nuccore", ids[1:1])
record = read_kegg_pathway(kegg_text)
mermaid = kegg_pathway_mermaid(read_pathway_graph(record))
```

