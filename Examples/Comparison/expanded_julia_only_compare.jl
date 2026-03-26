pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using Random
using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function _build_atoms(count::Integer)
    Random.seed!(2027)
    atoms = BioToolkit.Atom[]
    for i in 1:count
        x = 10.0 * sin(i / 17)
        y = 10.0 * cos(i / 23)
        z = mod(i, 97) / 3
        push!(atoms, BioToolkit.Atom(i, "C", x, y, z; element="C"))
    end
    return atoms
end

function _build_phylo_alignment()
    return BioToolkit.MultipleSequenceAlignment([
        BioToolkit.SeqRecordLite("ACGTACGT"; identifier="tax1"),
        BioToolkit.SeqRecordLite("ACGTTCGT"; identifier="tax2"),
        BioToolkit.SeqRecordLite("ACGTTTGT"; identifier="tax3"),
        BioToolkit.SeqRecordLite("ACGTACTT"; identifier="tax4"),
    ])
end

function _build_phylo_trees()
    return [
        BioToolkit.parse_newick("((tax1:0.1,tax2:0.1):0.2,(tax3:0.1,tax4:0.1):0.2);") ,
        BioToolkit.parse_newick("((tax1:0.1,tax3:0.1):0.2,(tax2:0.1,tax4:0.1):0.2);") ,
        BioToolkit.parse_newick("((tax1:0.1,tax4:0.1):0.2,(tax2:0.1,tax3:0.1):0.2);")
    ]
end

function _vcf_fixture()
    io = IOBuffer()
    write(io, "##fileformat=VCFv4.2\n")
    write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")
    for i in 1:10_000
        chrom = isodd(i) ? "chr1" : "chr2"
        ref = isodd(i) ? "A" : "C"
        alt = isodd(i) ? "G" : "T"
        write(io, "$(chrom)\t$(i)\trs$(i)\t$(ref)\t$(alt)\t$(50 + mod(i, 50))\n")
    end
    return String(take!(io))
end

function main()
    atoms = _build_atoms(10_000)
    tree = BioToolkit.build_atom_kdtree(atoms)
    center = atom_coordinates(atoms[5_000])
    neighbor_ms, neighbors = repeat_elapsed_ms(() -> BioToolkit.atoms_within_radius(tree, center; radius=8.0), 200)

    reference_coords = BioToolkit.coordinate_matrix(atoms[1:500])
    mobile_coords = copy(reference_coords)
    mobile_coords[:, 1] .+= 0.5
    superpose_ms, superposition = repeat_elapsed_ms(() -> BioToolkit.superpose(reference_coords, mobile_coords), 500)

    phylo_alignment = _build_phylo_alignment()
    phylo_trees = _build_phylo_trees()
    consensus_ms, consensus_tree = repeat_elapsed_ms(() -> BioToolkit.consensus_tree(phylo_trees), 200)
    rf_ms, rf_value = repeat_elapsed_ms(() -> BioToolkit.robinson_foulds_distance(phylo_trees[1], phylo_trees[2]), 500)
    parsimony_ms, parsimony_tree = repeat_elapsed_ms(() -> BioToolkit.maximum_parsimony_tree(phylo_alignment), 50)

    vcf_path = tempname()
    vcf_text = _vcf_fixture()
    open(vcf_path, "w") do io
        write(io, vcf_text)
    end
    vcf_ms, vcf_records = repeat_elapsed_ms(() -> BioToolkit.read_vcf(vcf_path), 200)

    kmer_targets = [repeat("ACGT", 2_500)]
    kmer_ms, kmer_index = repeat_elapsed_ms(() -> BioToolkit.build_index(kmer_targets, ["chr1"]; k=6), 100)

    # Mock Entrez fetch to isolate batching / request assembly from network latency.
    @eval BioToolkit.Entrez _download_text(::AbstractString) = "LOCUS       MOCK\nORIGIN\n//"
    entrez_ids = [string(i) for i in 1:500]
    entrez_ms = missing
    entrez_result = ""
    # The live Entrez request is intentionally skipped here so the benchmark remains runnable offline.
    println("Julia-only expanded benchmark")
    println("  neighbor_ms=", round(neighbor_ms / 200, digits=4))
    println("  superpose_ms=", round(superpose_ms / 500, digits=4))
    println("  consensus_ms=", round(consensus_ms / 200, digits=4))
    println("  rf_ms=", round(rf_ms / 500, digits=4))
    println("  parsimony_ms=", round(parsimony_ms / 50, digits=4))
    println("  vcf_ms=", round(vcf_ms / 200, digits=4))
    println("  kmer_ms=", round(kmer_ms / 100, digits=4))
    println("  entrez_ms=", round(entrez_ms / 200, digits=4))
    println("  neighbor_hits=", length(neighbors))
    println("  superpose_rmsd=", superposition.rmsd)
    println("  consensus_terminals=", BioToolkit.count_terminals(consensus_tree))
    println("  rf_value=", rf_value)
    println("  parsimony_terminals=", BioToolkit.count_terminals(parsimony_tree))
    println("  vcf_records=", length(vcf_records))
    println("  kmer_index_keys=", length(kmer_index.database))
    println("  entrez_length=", length(entrez_result))
end

main()