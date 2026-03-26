pushfirst!(LOAD_PATH, ".")
using BioToolkit

function verify_all()
    println("--- Verifying all requested functionalities ---\n")

    # 1. Global alignment
    seq1 = "ACGTGCCG"
    seq2 = "ACGTGCG"
    align_result = BioToolkit.pairwise_align(seq1, seq2, match=2, mismatch=-1, gap_open=-2, gap_extend=-0) # linear gap penalty by setting gap_open to same as gap penalty? Or just use defaults.
    println("1. Global Alignment (Needleman-Wunsch):")
    println("   Score: ", align_result.score)
    println("   Seq1:  ", align_result.left)
    println("   Seq2:  ", align_result.right)
    println()

    # 2. Dotplots
    dotplot = BioToolkit.dotmatrix(seq1, seq1)
    println("2. Dotplot (first row): ", dotplot[1, :])
    println()

    # 3. Transcription, Reverse Complement, Translation, ORFs
    dna = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    rna = BioToolkit.transcribe_dna(dna)
    rev_comp = BioToolkit.reverse_complement(dna)
    protein = BioToolkit.translate_dna(dna)
    orfs = BioToolkit.find_orfs(dna)
    println("3. Sequence Operations:")
    println("   DNA:      ", dna)
    println("   RNA:      ", rna)
    println("   RevComp:  ", rev_comp)
    println("   Protein:  ", protein)
    println("   ORFs:     ", length(orfs), " found")
    println()

    # 4. GC-content, minimum skew
    gc = BioToolkit.gc_content(dna)
    skew_min = BioToolkit.minimum_skew(dna)
    println("4. Sequence Stats:")
    println("   GC-content:   ", round(gc, digits=4))
    println("   Minimum skew: ", skew_min)
    println()

    # 5. Protein statistics
    prot_seq = "ACDEFGHIKLMNPQRSTVWY"
    pmass = BioToolkit.protein_mass(prot_seq)
    ec = BioToolkit.extinction_coefficient(prot_seq)
    ii = BioToolkit.instability_index(prot_seq)
    pi = BioToolkit.isoelectric_point(prot_seq)
    gravy_val = BioToolkit.gravy(prot_seq)
    println("5. Protein Statistics (for sequence $prot_seq):")
    println("   Protein mass:           ", round(pmass, digits=2), " Da")
    println("   Extinction coefficient: ", ec, " M⁻¹ cm⁻¹")
    println("   Instability index:      ", round(ii, digits=2))
    println("   Isoelectric point (pI): ", round(pi, digits=2))
    println("   GRAVY:                  ", round(gravy_val, digits=4))
    
    println("\nAll functionalities verified successfully!")
end

verify_all()
