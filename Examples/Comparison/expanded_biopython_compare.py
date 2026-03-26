from io import StringIO
from math import cos, sin
from pathlib import Path
from time import perf_counter

from Bio import Phylo, SearchIO
from Bio.Align import MultipleSeqAlignment
from Bio.PDB import Atom, NeighborSearch, Superimposer
from Bio.Phylo.Consensus import majority_consensus
from Bio.Phylo.TreeConstruction import ParsimonyScorer, NNITreeSearcher, ParsimonyTreeConstructor
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
from Bio.SeqRecord import SeqRecord


def time_call(func, repetitions=1):
    start = perf_counter()
    result = None
    for _ in range(repetitions):
        result = func()
    return (perf_counter() - start) * 1000.0, result


sequence = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" * 300)
script_dir = Path(__file__).resolve().parent
blast_xml = script_dir / "blast.xml"
blast_tab = script_dir / "blast.tab"
vcf_path = script_dir / "sample.vcf"

melting_ms, melting = time_call(lambda: MeltingTemp.Tm_Wallace(str(sequence)), 200)
slice_ms, slice_value = time_call(lambda: str(sequence[99:5000]), 200)

searchio_xml_ms, searchio_xml = time_call(lambda: list(SearchIO.parse(str(blast_xml), "blast-xml")), 20)
searchio_tab_ms, searchio_tab = time_call(lambda: list(SearchIO.parse(str(blast_tab), "blast-tab")), 200)

atoms = [Atom.Atom(f"C{i}", [10.0 * __import__("math").sin(i / 17), 10.0 * __import__("math").cos(i / 23), (i % 97) / 3], 1.0, 1.0, " ", f"C{i}", i) for i in range(1, 10_001)]
neighbor_search = NeighborSearch(atoms)
neighbor_ms, neighbor_hits = time_call(lambda: neighbor_search.search([10.0, 10.0, 10.0], 8.0, level="A"), 100)

reference_atoms = atoms[:500]
mobile_atoms = [Atom.Atom(atom.get_name(), [coord + 0.5 if index == 0 else coord for index, coord in enumerate(atom.get_coord())], 1.0, 1.0, " ", atom.get_fullname(), atom.get_serial_number()) for atom in reference_atoms]
superimposer = Superimposer()
superpose_ms, superpose_result = time_call(lambda: superimposer.set_atoms(reference_atoms, mobile_atoms), 500)

tree_a = Phylo.read(StringIO("((tax1:0.1,tax2:0.1):0.2,(tax3:0.1,tax4:0.1):0.2);"), "newick")
tree_b = Phylo.read(StringIO("((tax1:0.1,tax3:0.1):0.2,(tax2:0.1,tax4:0.1):0.2);"), "newick")
tree_c = Phylo.read(StringIO("((tax1:0.1,tax4:0.1):0.2,(tax2:0.1,tax3:0.1):0.2);"), "newick")
consensus_ms, consensus_tree = time_call(lambda: majority_consensus([tree_a, tree_b, tree_c], cutoff=0.5), 200)

try:
    import dendropy

    def tree_distance():
        left = dendropy.Tree.get(data="((tax1,tax2),(tax3,tax4));", schema="newick")
        right = dendropy.Tree.get(data="((tax1,tax3),(tax2,tax4));", schema="newick")
        return dendropy.calculate.treecompare.symmetric_difference(left, right)

    rf_ms, rf_value = time_call(tree_distance, 500)
except Exception:
    rf_ms, rf_value = float("nan"), "skipped"

scorer = ParsimonyScorer()
searcher = NNITreeSearcher(scorer)
parsimony_constructor = ParsimonyTreeConstructor(searcher, starting_tree=tree_a)
msa = MultipleSeqAlignment([
    SeqRecord(Seq("ACGTACGT"), id="tax1"),
    SeqRecord(Seq("ACGTTCGT"), id="tax2"),
    SeqRecord(Seq("ACGTTTGT"), id="tax3"),
    SeqRecord(Seq("ACGTACTT"), id="tax4"),
])
parsimony_ms, parsimony_tree = time_call(lambda: parsimony_constructor.build_tree(msa), 10)

print(f"python_melting_ms={melting_ms:.4f}")
print(f"python_slice_ms={slice_ms:.4f}")
print(f"python_searchio_xml_ms={searchio_xml_ms:.4f}")
print(f"python_searchio_tab_ms={searchio_tab_ms:.4f}")
print(f"python_neighbor_ms={neighbor_ms:.4f}")
print(f"python_superpose_ms={superpose_ms:.4f}")
print(f"python_slice_len={len(slice_value)}")
print(f"python_searchio_xml_hits={len(searchio_xml[0]) if searchio_xml else 0}")
print(f"python_searchio_tab_hits={len(searchio_tab[0]) if searchio_tab else 0}")
print(f"python_neighbor_hits={len(neighbor_hits)}")
print(f"python_superpose_rmsd={getattr(superpose_result, 'rms', '0.0')}")
print(f"python_consensus_terminals={len(consensus_tree.get_terminals())}")
print(f"python_rf_ms={rf_ms if isinstance(rf_ms, float) else float('nan'):.4f}")
print(f"python_rf_value={rf_value}")
print(f"python_parsimony_terminals={len(parsimony_tree.get_terminals())}")

if vcf_path.exists():
    try:
        import vcf
        vcf_ms, vcf_records = time_call(lambda: list(vcf.Reader(filename=str(vcf_path))), 10)
        print(f"python_vcf_ms={vcf_ms:.4f}")
        print(f"python_vcf_records={len(vcf_records)}")
    except Exception:
        print("python_vcf_ms=skipped")
