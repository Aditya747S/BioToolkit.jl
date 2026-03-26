# enrichment.jl

## Purpose
This file implements functional enrichment analysis for BioToolkit. It provides the data model for annotated terms, the mapping layer for gene identifiers, the overrepresentation statistics, and helper functions for GO and KEGG style analyses.

## Main structs
- IDMapper: maps input identifiers to canonical IDs and stores reverse lookup information.
- EnrichmentTerm: one annotation term with its ID, name, namespace, member genes, and parent terms.
- EnrichmentDatabase: a collection of enrichment terms plus the identifier mapper used by the database.
- EnrichmentResult: one enriched term with overlap statistics, p-values, adjusted p-values, odds ratio, and the matched genes.

## Public functions
- IDMapper(): create an empty mapper.
- build_annotation_database(terms; mapper): construct an enrichment database from a list of terms.
- builtin_annotation_terms and builtin_annotation_database: provide a built-in GO/KEGG-style example database.
- load_annotation_database(path) and save_annotation_database(path, database): serialize and restore annotation databases as JSON.
- map_id and map_ids: normalize one identifier or many identifiers through the mapper.
- enrichment_test(query_genes, database; background, namespace, elim, threaded, min_overlap): perform general overrepresentation testing.
- go_enrichment and kegg_enrichment: namespace-specific wrappers around enrichment_test.
- dotplot: plotting helper for displaying enriched results.

## What the module does
The module takes a query gene list and compares it against term gene sets stored in the database. It filters the query into the background universe, computes a right-tail Fisher-style enrichment p-value for each term, adjusts p-values with Benjamini-Hochberg correction, and returns a sorted list of EnrichmentResult objects.

## How the structs work together
EnrichmentTerm holds the biological annotation itself. EnrichmentDatabase stores the full term map and the IDMapper used to translate gene symbols or other identifiers. EnrichmentResult is the final statistical output that a user inspects, plots, or filters for biological interpretation.

## Typical usage
1. Build or load an EnrichmentDatabase.
2. Optionally map a gene list through map_ids so identifiers are consistent with the database.
3. Call enrichment_test, go_enrichment, or kegg_enrichment.
4. Sort or filter the returned EnrichmentResult objects by padj or overlap.
5. Pass the result list to dotplot if you want a visualization of the top terms.

## Important implementation details
- _background_set collects the universe of genes represented in the database.
- _effective_gene_sets supports elim-style term pruning when descendants or significant terms should be removed from parents.
- _fisher_right_tail computes the right-tail hypergeometric enrichment probability.
- _odds_ratio computes a smoothed odds ratio for the contingency table.
- _benjamini_hochberg provides multiple-testing correction for the final term list.

## Why this file matters
This module is the pathway and ontology interpretation layer of BioToolkit. It turns a gene list into ranked, statistically supported biological terms that can be reported or plotted.
