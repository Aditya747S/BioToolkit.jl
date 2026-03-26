# schema.jl

## Purpose
This file defines the compact variant schema used for efficient storage and transport of genomic variant events. It is designed for workflows that need to convert text-based VCF-like records into a smaller, Arrow-friendly representation.

## Main structs
- VariantEvent: a compact encoded variant record containing chromosome code, position, hashed ID, encoded ref/alt bases, quality, and a quality-present flag.
- VariantTextRecord: the text-level variant representation used as the intermediate conversion format.

## Public functions
- arrow_schema(::Type{VariantEvent}): expose the Arrow/Tables schema for VariantEvent.
- encode_chromosome(chrom): map chromosome labels such as chr1, X, Y, and MT into compact integer codes.
- encode_base(base): encode nucleotide bases into small integer codes.
- encode_identifier(identifier): hash a string identifier into a compact UInt32.
- compact_variant_event(record): convert a VariantTextRecord into a compact VariantEvent.

## What the module does
The module compresses a text variant record into a small fixed-width structure that is easier to store in memory and easier to pass through Arrow-compatible pipelines. It keeps the important coordinates and identity information, but replaces the textual fields with numeric codes where possible.

## How the structs work together
VariantTextRecord is the human-readable or parser-facing representation. VariantEvent is the compact representation used after encoding. The conversion function bridges the two, while the encoding helpers keep chromosome labels, bases, and identifiers standardized.

## Typical usage
1. Parse or construct VariantTextRecord objects from VCF-like input.
2. Call compact_variant_event to convert them to VariantEvent.
3. Use the VariantEvent array in Arrow-compatible storage or high-throughput processing.
4. Decode or reformat as needed through the text-based pathway in io.jl.

## Important implementation details
- encode_chromosome handles standard autosomes plus X, Y, and mitochondrial labels.
- encode_base maps bases to compact codes and falls back to 0 for unknown values.
- encode_identifier uses an FNV-style hash to create stable compact IDs.
- compact_variant_event stores both the numeric quality and a Boolean indicating whether quality was present in the source.

## Why this file matters
This module is the space-efficient variant representation layer of BioToolkit. It makes large-scale variant data easier to store, move, and analyze without losing the ability to round-trip from text records.
