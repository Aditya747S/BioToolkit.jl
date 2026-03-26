using BioToolkit

sequence = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

println("DNA workflow example")
println("  valid DNA: ", BioToolkit.validate_dna(sequence))
println("  nucleotide counts: ", BioToolkit.count_nucleotides(sequence))
println("  transcribed RNA: ", BioToolkit.transcribe_dna(sequence))
println("  reverse complement: ", BioToolkit.reverse_complement(sequence))
println("  GC content: ", round(BioToolkit.gc_content(sequence), digits=3))
println("  translated protein: ", BioToolkit.translate_dna(sequence; stop_at_stop=true))
println("  ORFs found: ", length(BioToolkit.find_orfs(sequence; min_aa=5)))
println("  3-mer frequency: ", BioToolkit.kmer_frequency(sequence, 3)["ATG"])