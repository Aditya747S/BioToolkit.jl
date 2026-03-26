pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
using BioSequences
using BioAlignments
using BED
using GFF3

include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function main()
    dna_left = repeat("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 500)
    dna_right_chars = collect(dna_left)
    for position in (1, 101, 201, 301, 401)
        dna_right_chars[position] = 'T'
    end
    dna_right = String(dna_right_chars)

    bed_line = "chrX\t151080532\t151081699"
    gff3_line = "CCDS1.1\tCCDS\tgene\t801943\t802434\t.\t-\t.\tNAME=LINC00115"

    biosequence_left = LongDNA{2}(dna_left)
    biosequence_right = LongDNA{2}(dna_right)

    reverse_ms, reverse_result = repeat_elapsed_ms(() -> BioToolkit.reverse_complement(dna_left), 400)
    biosequences_reverse_ms, biosequences_reverse_result = repeat_elapsed_ms(() -> BioSequences.reverse_complement(copy(biosequence_left)), 400)


    hamming_ms, hamming_result = repeat_elapsed_ms(() -> BioToolkit.hamming_distance(dna_left, dna_right), 1_000)
    biosequences_hamming_ms, biosequences_hamming_result = repeat_elapsed_ms(() -> BioSequences.mismatches(biosequence_left, biosequence_right), 1_000)
    bioalignments_hamming_ms, bioalignments_hamming_result = repeat_elapsed_ms(() -> distance(pairalign(HammingDistance(), biosequence_left, biosequence_right, distance_only=true)), 1_000)

    bed_parse_ms, bed_record = repeat_elapsed_ms(() -> BioToolkit.parse_bed_record(bed_line), 50_000)
    bed_package_ms, bed_package_record = repeat_elapsed_ms(() -> BED.Record(bed_line), 50_000)

    gff3_parse_ms, gff3_record = repeat_elapsed_ms(() -> BioToolkit.parse_gff_record(gff3_line), 50_000)
    gff3_package_ms, gff3_package_record = repeat_elapsed_ms(() -> GFF3.Record(gff3_line), 50_000)

    @assert reverse_result == String(biosequences_reverse_result)
    @assert hamming_result == biosequences_hamming_result
    @assert bed_record.chrom == BED.chrom(bed_package_record)
    @assert bed_record.start + 1 == Int32(BED.chromstart(bed_package_record))
    @assert bed_record.stop == Int32(BED.chromend(bed_package_record))
    @assert gff3_record.chrom == GFF3.seqid(gff3_package_record)
    @assert gff3_record.source == GFF3.source(gff3_package_record)
    @assert gff3_record.feature == GFF3.featuretype(gff3_package_record)
    @assert gff3_record.start == Int32(GFF3.seqstart(gff3_package_record))
    @assert gff3_record.stop == Int32(GFF3.seqend(gff3_package_record))
    @assert gff3_record.score === missing
    @assert gff3_record.phase === missing

    println("Julia package comparison benchmark")
    println("  reverse_ms=", round(reverse_ms, digits=4))
    println("  biosequences_reverse_ms=", round(biosequences_reverse_ms, digits=4))
    println("  hamming_ms=", round(hamming_ms, digits=4))
    println("  biosequences_hamming_ms=", round(biosequences_hamming_ms, digits=4))
    println("  bioalignments_hamming_ms=", round(bioalignments_hamming_ms, digits=4))
    println("  bed_parse_ms=", round(bed_parse_ms, digits=4))
    println("  bed_package_ms=", round(bed_package_ms, digits=4))
    println("  gff3_parse_ms=", round(gff3_parse_ms, digits=4))
    println("  gff3_package_ms=", round(gff3_package_ms, digits=4))
    println("  bed_summary=", bed_record.chrom, ":", bed_record.start, "-", bed_record.stop)
    println("  gff3_summary=", gff3_record.chrom, ":", gff3_record.start, "-", gff3_record.stop)
end

main()
