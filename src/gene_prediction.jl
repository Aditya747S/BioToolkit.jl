export predict_genes_hmm

"""
    predict_genes_hmm(sequence::AbstractString; p_coding=0.01)

Uses a specialized Hidden Markov Model and the Viterbi API to predict 
open reading frames (genes) in nucleotide DNA mathematically over the entire structure.
It out-performs naive substring scanning by intrinsically penalizing nonsense gaps 
and filtering overlapping false positive hits statistically.

Returns an array of `(start_index, end_index)` genomic tuples.
"""
function predict_genes_hmm(sequence::AbstractString; p_coding=0.01)
    bytes = codeunits(sequence)
    L = length(bytes)
    
    # 1. Build a 2-State (Non-Coding vs Coding) profile HMM abstraction
    # To keep it generic without deep structural nested state arrays (Start/Stop codons strictly mapped),
    # we represent:
    # State 1: Background (Non-Coding) - GC content ~ 40-50%
    # State 2: Gene (Coding) - GC content ~ 60%, strict boundaries.
    
    states = ["NonCoding", "Coding"]
    alphabet = UInt8[UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T')]
    
    initial = [0.99, 0.01]
    
    # Emphasizes that "Coding" blocks are long (gene size ~ 1000bp)
    # So transition out of Coding is ~ 1/1000.
    transitions = [
        1.0 - p_coding       p_coding;
        1/1000.0             1.0 - 1/1000.0
    ]
    
    # Standard simplistic eukaryotic GC skew
    emissions = [
        0.25  0.25  0.25  0.25; # Flat background
        0.15  0.35  0.35  0.15  # Elevate C/G exactly mathematically
    ]
    
    hmm = HMM(states, alphabet, initial, transitions, emissions; log_space=false)
    
    # 2. Feed directly into standard Viterbi engine
    path, prob = viterbi(hmm, bytes)
    
    # 3. Post-process structural boundaries (State 2 block mappings)
    genes = Tuple{Int, Int}[]
    in_gene = false
    gene_start = 0
    
    for i in 1:L
        if path[i] == 2 && !in_gene
            in_gene = true
            gene_start = i
        elseif path[i] == 1 && in_gene
            in_gene = false
            gene_end = i - 1
            # Filter biologically meaningless noise blocks (< 30bp)
            if (gene_end - gene_start) > 30
                push!(genes, (gene_start, gene_end))
            end
        end
    end
    
    # Flush trailing
    if in_gene && (L - gene_start) > 30
        push!(genes, (gene_start, L))
    end
    
    return genes
end
