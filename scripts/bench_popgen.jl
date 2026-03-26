using BioToolkit
using BenchmarkTools

# 1. Allele Frequencies
function bench_allele_freq()
    l = Locus(('A', 'B'))
    inds = [PopGenIndividual("ind$i", [l]) for i in 1:10000]
    pop = Population("BenchPop", inds)
    @btime allele_frequencies($pop, 1)
end

# 2. Tajima's D
function bench_tajimas_d()
    # Create a mock MSA with 1k sequences of 1k bp
    seq = "A"^1000
    records = [SeqRecordLite(seq, identifier="s$i") for i in 1:1000]
    # Add some mutations
    for i in 1:100
        records[i].sequence = "T" * records[i].sequence[2:end]
    end
    msa = MultipleSequenceAlignment(records)
    @btime tajimas_d($msa)
end

println("Allele Frequency (10k individuals):")
bench_allele_freq()

println("\nTajima's D (1k sequences, 1k bp):")
bench_tajimas_d()
