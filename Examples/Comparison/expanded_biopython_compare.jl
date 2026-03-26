pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using Random
using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function _write_python_sequence_fixture(path::AbstractString, sequence::AbstractString, blast_xml::AbstractString, blast_tabular::AbstractString, vcf_text::AbstractString)
    open(path, "w") do io
        write(io, "from time import perf_counter\n")
        write(io, "from pathlib import Path\n")
        write(io, "from Bio.Seq import Seq\n")
        write(io, "from Bio.SeqUtils import MeltingTemp\n")
        write(io, "sequence = Seq('$(sequence)')\n")
        write(io, "blast_xml = $(repr(blast_xml))\n")
        write(io, "blast_tabular = $(repr(blast_tabular))\n")
        write(io, "vcf_text = $(repr(vcf_text))\n")
        write(io, "def time_call(func, repetitions=100):\n")
        write(io, "    start = perf_counter()\n")
        write(io, "    result = None\n")
        write(io, "    for _ in range(repetitions):\n")
        write(io, "        result = func()\n")
        write(io, "    return (perf_counter() - start) * 1000.0, result\n")
        write(io, "def safe_import(module_name):\n")
        write(io, "    try:\n")
        write(io, "        module = __import__(module_name)\n")
        write(io, "        return module\n")
        write(io, "    except Exception:\n")
        write(io, "        return None\n")
        write(io, "melting_ms, melting = time_call(lambda: MeltingTemp.Tm_Wallace(str(sequence)))\n")
        write(io, "slice_ms, slice_value = time_call(lambda: str(sequence[99:5000]))\n")
        write(io, "with open('blast.xml', 'w') as handle:\n")
        write(io, "    handle.write(blast_xml)\n")
        write(io, "with open('blast.tab', 'w') as handle:\n")
        write(io, "    handle.write(blast_tabular)\n")
        write(io, "from Bio import SearchIO\n")
        write(io, "searchio_xml_ms, searchio_xml = time_call(lambda: list(SearchIO.parse('blast.xml', 'blast-xml')), 20)\n")
        write(io, "searchio_tab_ms, searchio_tab = time_call(lambda: list(SearchIO.parse('blast.tab', 'blast-tab')), 200)\n")
        write(io, "print(f'python_melting_ms={melting_ms:.4f}')\n")
        write(io, "print(f'python_slice_ms={slice_ms:.4f}')\n")
        write(io, "print(f'python_searchio_xml_ms={searchio_xml_ms:.4f}')\n")
        write(io, "print(f'python_searchio_tab_ms={searchio_tab_ms:.4f}')\n")
        write(io, "print(f'python_slice_len={len(slice_value)}')\n")
        write(io, "print(f'python_searchio_xml_hits={len(searchio_xml[0]) if searchio_xml else 0}')\n")
        write(io, "print(f'python_searchio_tab_hits={len(searchio_tab[0]) if searchio_tab else 0}')\n")
        write(io, "pvcf = safe_import('vcf')\n")
        write(io, "if pvcf is None:\n")
        write(io, "    print('python_vcf_ms=skipped')\n")
        write(io, "else:\n")
        write(io, "    with open('sample.vcf', 'w') as handle:\n")
        write(io, "        handle.write(vcf_text)\n")
        write(io, "    vcf_ms, vcf_records = time_call(lambda: list(pvcf.Reader(filename='sample.vcf')), 10)\n")
        write(io, "    print(f'python_vcf_ms={vcf_ms:.4f}')\n")
        write(io, "    print(f'python_vcf_records={len(vcf_records)}')\n")
    end
end

function _blast_xml_fixture()
    return """<?xml version=\"1.0\"?>
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.15.0+</BlastOutput_version>
  <BlastOutput_db>synthetic_db</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>synthetic query</BlastOutput_query-def>
  <BlastOutput_query-len>12</BlastOutput_query-len>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>synthetic query</Iteration_query-def>
      <Iteration_query-len>12</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_id>subject1</Hit_id>
          <Hit_def>subject one</Hit_def>
          <Hit_len>12</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>42.0</Hsp_bit-score>
              <Hsp_evalue>1e-10</Hsp_evalue>
              <Hsp_query-from>2</Hsp_query-from>
              <Hsp_query-to>9</Hsp_query-to>
              <Hsp_hit-from>3</Hsp_hit-from>
              <Hsp_hit-to>10</Hsp_hit-to>
              <Hsp_identity>8</Hsp_identity>
              <Hsp_positive>8</Hsp_positive>
              <Hsp_align-len>8</Hsp_align-len>
              <Hsp_qseq>ACGTACGT</Hsp_qseq>
              <Hsp_hseq>ACGTACGT</Hsp_hseq>
              <Hsp_midline>||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""
end

function _blast_tabular_fixture()
    return """query1\tsubject1\t100.0\t8\t0\t0\t2\t9\t3\t10\t1e-10\t42.0
query1\tsubject2\t87.5\t8\t1\t0\t5\t12\t1\t8\t2e-5\t31.0
query2\tsubject3\t75.0\t4\t1\t0\t1\t4\t4\t7\t0.01\t20.0
"""
end

function _vcf_fixture()
    io = IOBuffer()
    write(io, "##fileformat=VCFv4.2\n")
    write(io, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\n")
    for i in 1:1000
        chrom = isodd(i) ? "chr1" : "chr2"
        ref = isodd(i) ? "A" : "C"
        alt = isodd(i) ? "G" : "T"
        write(io, "$(chrom)\t$(i)\trs$(i)\t$(ref)\t$(alt)\t$(50 + mod(i, 50))\n")
    end
    return String(take!(io))
end

function _build_atoms(count::Integer)
    Random.seed!(2026)
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
    sequences = [
        BioToolkit.SeqRecordLite("ACGTACGT"; identifier="tax1"),
        BioToolkit.SeqRecordLite("ACGTTCGT"; identifier="tax2"),
        BioToolkit.SeqRecordLite("ACGTTTGT"; identifier="tax3"),
        BioToolkit.SeqRecordLite("ACGTACTT"; identifier="tax4"),
    ]
    return BioToolkit.MultipleSequenceAlignment(sequences)
end

function _build_phylo_trees()
    trees = [
        BioToolkit.parse_newick("((tax1:0.1,tax2:0.1):0.2,(tax3:0.1,tax4:0.1):0.2);")
        BioToolkit.parse_newick("((tax1:0.1,tax3:0.1):0.2,(tax2:0.1,tax4:0.1):0.2);")
        BioToolkit.parse_newick("((tax1:0.1,tax4:0.1):0.2,(tax2:0.1,tax3:0.1):0.2);")
    ]
    return trees
end

function main()
    sequence = repeat("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 300)
    sequence_slice = sequence[100:5000]
    blast_xml = _blast_xml_fixture()
    blast_tabular = _blast_tabular_fixture()
    vcf_text = _vcf_fixture()

    slice_ms, slice_value = repeat_elapsed_ms(() -> sequence[100:5000], 200)
    melting_ms, melting_value = repeat_elapsed_ms(() -> BioToolkit.melting_temp(sequence), 200)
    blast_xml_ms, blast_xml_records = repeat_elapsed_ms(() -> BioToolkit.parse_blast_xml(IOBuffer(blast_xml)), 500)
    blast_tab_ms, blast_tab_records = repeat_elapsed_ms(() -> BioToolkit.parse_blast_tabular(IOBuffer(blast_tabular)), 500)

    atoms = _build_atoms(10_000)
    tree = BioToolkit.build_atom_kdtree(atoms)
    center = atom_coordinates(atoms[5_000])
    neighbor_ms, neighbors = repeat_elapsed_ms(() -> BioToolkit.atoms_within_radius(tree, center; radius=8.0), 100)

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
    open(vcf_path, "w") do io
        write(io, vcf_text)
    end
    vcf_ms, vcf_records = repeat_elapsed_ms(() -> BioToolkit.read_vcf(vcf_path), 200)

    kmer_targets = [repeat("ACGT", 2500)]
    kmer_ms, kmer_index = repeat_elapsed_ms(() -> BioToolkit.build_index(kmer_targets, ["chr1"]; k=6), 100)

    python_script = mktempdir() do dir
        script_path = joinpath(dir, "expanded_biopython_compare.py")
        _write_python_sequence_fixture(script_path, sequence, blast_xml, blast_tabular, vcf_text)
        read(`conda run -n general --no-capture-output python $script_path`, String)
    end

    println("Expanded BioToolkit vs Biopython benchmark")
    println("  slice_ms=", round(slice_ms / 200, digits=4))
    println("  melting_ms=", round(melting_ms / 200, digits=4))
    println("  blast_xml_ms=", round(blast_xml_ms / 500, digits=4))
    println("  blast_tab_ms=", round(blast_tab_ms / 500, digits=4))
    println("  neighbor_ms=", round(neighbor_ms / 100, digits=4))
    println("  superpose_ms=", round(superpose_ms / 500, digits=4))
    println("  consensus_ms=", round(consensus_ms / 200, digits=4))
    println("  rf_ms=", round(rf_ms / 500, digits=4))
    println("  parsimony_ms=", round(parsimony_ms / 50, digits=4))
    println("  vcf_ms=", round(vcf_ms / 200, digits=4))
    println("  kmer_ms=", round(kmer_ms / 100, digits=4))
    println("  slice_len=", length(slice_value))
    println("  melting=", round(melting_value, digits=4))
    println("  blast_xml_hits=", length(blast_xml_records[1].hits))
    println("  blast_tab_hits=", length(blast_tab_records[1].hits))
    println("  neighbor_hits=", length(neighbors))
    println("  superpose_rmsd=", superposition.rmsd)
    println("  consensus_terminals=", BioToolkit.count_terminals(consensus_tree))
    println("  rf_value=", rf_value)
    println("  parsimony_terminals=", BioToolkit.count_terminals(parsimony_tree))
    println("  vcf_records=", length(vcf_records))
    println("  kmer_index_keys=", length(kmer_index.database))
    print(python_script)
end

main()