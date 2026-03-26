import sys
import os
import time

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.SeqUtils import GC_skew
except ImportError:
    print("Biopython not found.")
    sys.exit(1)

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
try:
    from benchmark_helpers import format_ms, repeat_elapsed_ms
except ImportError:
    # Inline repeat_elapsed_ms if not found
    def repeat_elapsed_ms(func, reps):
        import time
        start = time.perf_counter()
        for _ in range(reps):
            func()
        return ((time.perf_counter() - start) * 1000) / reps


    def aliphatic_index(sequence):
        length = len(sequence)
        counts = {"A": 0, "V": 0, "I": 0, "L": 0}
        for residue in sequence.upper():
            if residue in counts:
                counts[residue] += 1
        x_a = 100.0 * counts["A"] / length
        x_v = 100.0 * counts["V"] / length
        x_i = 100.0 * counts["I"] / length
        x_l = 100.0 * counts["L"] / length
        return x_a + 2.9 * x_v + 3.9 * (x_i + x_l)

    def format_ms(val):
        return f"{val:.4f}"

protein_seq = "ACDEFGHIKLMNPQRSTVWY" * 1000 # 20k residues
dna_seq = "ACGTGCGCAGCT" * 5000 # 60k bp

# Heat up
analysis = ProteinAnalysis(protein_seq)
weight = analysis.molecular_weight()
extinction = analysis.molar_extinction_coefficient()
instability = analysis.instability_index()
isoelectric = analysis.isoelectric_point()
gravy = analysis.gravy()
protparam = {
    "molecular_weight": analysis.molecular_weight(),
    "extinction_coefficient": analysis.molar_extinction_coefficient()[0],
    "instability_index": analysis.instability_index(),
    "isoelectric_point": analysis.isoelectric_point(),
    "aliphatic_index": aliphatic_index(protein_seq),
    "gravy": analysis.gravy(),
}
skew = GC_skew(dna_seq, window=10)

def calc_protein():
    a = ProteinAnalysis(protein_seq)
    _ = a.molecular_weight()
    _ = a.molar_extinction_coefficient()
    _ = a.instability_index()
    _ = a.isoelectric_point()
    _ = a.gravy()


def calc_protein_bundle():
    a = ProteinAnalysis(protein_seq)
    _ = {
        "molecular_weight": a.molecular_weight(),
        "extinction_coefficient": a.molar_extinction_coefficient()[0],
        "instability_index": a.instability_index(),
        "isoelectric_point": a.isoelectric_point(),
        "aliphatic_index": aliphatic_index(protein_seq),
        "gravy": a.gravy(),
    }

def calc_skew():
    _ = GC_skew(dna_seq, window=10)

reps = 100

protein_ms = repeat_elapsed_ms(calc_protein, reps)
bundle_ms = repeat_elapsed_ms(calc_protein_bundle, reps)
skew_ms = repeat_elapsed_ms(calc_skew, reps)

print("Python Protein Stats benchmark")
print(f"  reps: {reps}")
print(f"  python_protein_ms={format_ms(protein_ms)}")
print(f"  python_protein_bundle_ms={format_ms(bundle_ms)}")
print(f"  python_skew_ms={format_ms(skew_ms)}")
