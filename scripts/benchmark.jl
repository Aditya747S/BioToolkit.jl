pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..")))

using BioToolkit
include(joinpath(@__DIR__, "benchmark_helpers.jl"))

function write_sample_vcf(path::AbstractString, n::Integer)
    open(path, "w") do io
        write(io, "##fileformat=VCFv4.2\n")
        write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")
        for i in 1:n
            chrom = isodd(i) ? "chr1" : "chr2"
            pos = i
            id = "rs$i"
            ref = isodd(i) ? "A" : "C"
            alt = isodd(i) ? "G" : "T"
            qual = 50.0 + Float64(mod(i, 50))
            write(io, "$(chrom)\t$(pos)\t$(id)\t$(ref)\t$(alt)\t$(qual)\n")
        end
    end
end

function write_sample_bed(path::AbstractString, n::Integer)
    open(path, "w") do io
        for i in 1:n
            chrom = isodd(i) ? "chr1" : "chr2"
            start = 10 * (i - 1)
            stop = start + 10
            write(io, "$(chrom)\t$(start)\t$(stop)\n")
        end
    end
end

function elapsed_ms(f)
    GC.gc()
    start_time = time_ns()
    result = f()
    return (time_ns() - start_time) / 1e6, result
end

mktempdir() do dir
    vcf_path = joinpath(dir, "sample.vcf")
    bed_path = joinpath(dir, "sample.bed")
    vcf_arrow = joinpath(dir, "sample-vcf.arrow")
    bed_arrow = joinpath(dir, "sample-bed.arrow")

    write_sample_vcf(vcf_path, 20_000)
    write_sample_bed(bed_path, 20_000)

    ingest_ms, _ = elapsed_ms() do
        BioToolkit.ingest_vcf(vcf_path, vcf_arrow; chunk_size=2_000)
    end

    bed_ms, _ = elapsed_ms() do
        BioToolkit.ingest_bed(bed_path, bed_arrow; chunk_size=2_000)
    end

    table = BioToolkit.load_arrow_table(vcf_arrow)

    filter_ms, subset = elapsed_ms() do
        BioToolkit.filter_region(table, "chr1", 5_000, 15_000; sorted=false)
    end

    threaded_hist_ms, _ = elapsed_ms() do
        BioToolkit.bin_positions(subset.pos, 100; threaded=true)
    end

    serial_hist_ms, _ = elapsed_ms() do
        BioToolkit.bin_positions(subset.pos, 100; threaded=false)
    end

    println("BioToolkit benchmark")
    println("  ingest_vcf: $(round(ingest_ms, digits=2)) ms")
    println("  ingest_bed: $(round(bed_ms, digits=2)) ms")
    println("  filter_region: $(round(filter_ms, digits=2)) ms")
    println("  bin_positions threaded: $(round(threaded_hist_ms, digits=2)) ms")
    println("  bin_positions serial: $(round(serial_hist_ms, digits=2)) ms")

    restriction_sequence = repeat("GATTACAGAATTCCTGA", 5_000)
    restriction_ms, restriction_sites = repeat_elapsed_ms(() -> BioToolkit.find_restriction_sites(restriction_sequence, "EcoRI"), 20)

    motif_counts = BioToolkit.MotifCounts(['A', 'C', 'G', 'T'], [10 2 1 1; 2 10 1 1; 1 1 10 2; 1 1 2 10])
    logo_ms, _ = repeat_elapsed_ms(() -> BioToolkit.sequence_logo_svg(motif_counts; width=480, height=180, title="BioToolkit"), 20)

    entrez_payload = """
{"header":{"type":"esearch","version":"0.3"},"esearchresult":{"count":"2","retmax":"2","retstart":"0","idlist":["123","456"],"translationset":[],"translationstack":[],"querytranslation":"BRCA1"}}
"""
    entrez_ms, _ = repeat_elapsed_ms(() -> BioToolkit.parse_entrez_search_response(entrez_payload), 1_000)

    pathway_text = """
ENTRY       map00010                      Pathway
NAME        Glycolysis / Gluconeogenesis
ENZYME      1.1.1.1 2.7.1.1
COMPOUND    C00031 C00022
GENE        b0001  gene1
///
"""
    pathway_ms, _ = repeat_elapsed_ms(() -> BioToolkit.read_kegg_pathway(IOBuffer(pathway_text)), 1_000)

    scop_ms, _ = repeat_elapsed_ms(() -> BioToolkit.parse_scop_record("d1tq3a_ 1tq3 A:1-100 a.1.1.1 Protein description"), 10_000)
    cath_ms, _ = repeat_elapsed_ms(() -> BioToolkit.parse_cath_record("1abcA00 1abc A 1 10 20 30 Example domain"), 10_000)

    println("  restriction_sites (20 reps): $(round(restriction_ms / 20, digits=4)) ms per call; hits=$(length(restriction_sites))")
    println("  sequence_logo_svg (20 reps): $(round(logo_ms / 20, digits=4)) ms per call")
    println("  parse_entrez_search_response (1000 reps): $(round(entrez_ms / 1000, digits=6)) ms per call")
    println("  read_kegg_pathway (1000 reps): $(round(pathway_ms / 1000, digits=6)) ms per call")
    println("  parse_scop_record (10000 reps): $(round(scop_ms / 10000, digits=6)) ms per call")
    println("  parse_cath_record (10000 reps): $(round(cath_ms / 10000, digits=6)) ms per call")
end
