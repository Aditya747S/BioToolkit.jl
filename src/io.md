# io.jl

## Purpose
This file implements text-format I/O for common genomics records. It covers VCF, BED, GFF, and GenBank-style structures, plus the parsing and rendering helpers needed to move between raw text and BioToolkit's internal record types.

## Main structs
- VariantTextRecord: a compact text-level VCF record with chromosome, position, ID, ref, alt, and quality.
- BedRecord: a BED interval with chromosome, start, and stop.
- GffRecord: a full GFF record with source, feature type, score, strand, phase, attributes, and parsed attribute map.
- GenBankFeature: a GenBank feature with a key, a textual location, qualifier map, and parsed location placeholder.
- GenBankRecord: a structured GenBank entry with locus information, sequence, and a list of features.
- GenBankArrowRecord: a compact Arrow-friendly representation of a GenBank record.

## Public functions
- parse_vcf_record(line): parse one VCF record line into a VariantTextRecord.
- read_vcf(input or io): read all VCF records from a file or IO stream.
- write_vcf(output or io, records): write VariantTextRecord objects back to VCF text.
- parse_bed_record(line): parse one BED record.
- read_bed(input or io): read a BED file or stream.
- GffRecord(...): construct a typed GFF record from field values.
- GenBankFeature(...): construct a GenBank feature and optionally attach parsed locations.
- GenBankRecord(...): construct a full GenBank record.
- parse_gff_attributes(attributes): decode GFF attribute strings into a dictionary.
- _render_gff_attributes(attributes, attribute_map): turn parsed attributes back into text.

## What the module does
The module converts between plain text file formats and typed Julia records. It validates the minimum field counts, normalizes coordinate fields, and keeps the parsed attribute maps around so other modules can inspect or rewrite annotations without re-parsing the original strings.

## How the structs work together
VariantTextRecord is the intermediate VCF object used before compact encoding or downstream analysis. BedRecord and GffRecord represent standard interval formats. GenBankFeature and GenBankRecord store richer annotation data for GenBank-style sequence files. GenBankArrowRecord provides a compact export format for tabular or Arrow-based workflows.

## Typical usage
1. Read a text file with read_vcf or read_bed.
2. Convert parsed records into genomic interval or annotation objects in other modules.
3. Use parse_gff_attributes to inspect annotation tags.
4. Use write_vcf when you need to export processed variant records.
5. Convert GenBank or GFF objects into BioToolkit annotation types for downstream slicing or feature querying.

## Important implementation details
- parse_vcf_record and parse_bed_record are defensive: malformed lines become explicit errors instead of silent failures.
- GffRecord stores both the original text and the parsed attribute map so the record can be round-tripped or inspected.
- GenBankFeature keeps a parsed_location placeholder so higher-level annotation code can attach structured locations later.
- The module is designed as a format layer, not a biological analysis layer; it feeds richer modules with clean input objects.

## Why this file matters
This module is the package's text-format boundary. It makes it possible to load real genomics files into typed Julia structures and then hand those structures to the rest of BioToolkit.
