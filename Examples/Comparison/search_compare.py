import sys
import os

try:
    from Bio import pairwise2
except ImportError:
    print("Biopython not found")
    sys.exit(1)

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "scripts"))
try:
    from benchmark_helpers import format_ms, repeat_elapsed_ms
except ImportError:
    def repeat_elapsed_ms(func, reps):
        import time
        start = time.perf_counter()
        for _ in range(reps):
            func()
        return ((time.perf_counter() - start) * 1000)

print("Generating sequences...")
# Target: 100,000 bp
target_base = "ACGT" * 25000
target = target_base[:50000] + "GGCCGATCGATCGATCGATCGA" + target_base[50022:]

# Query: 1,000 bp
query = ("T" * 500) + "GGCCGATCGATCGATCGATCGA" + ("A" * 478)

def do_python_local():
    # Match=2, mismatch=-1, open_penalty=-5, extend_penalty=-1 (Since biopython requires negative values)
    # Using pairwise2.align.localms
    _ = pairwise2.align.localms(query, target, 2, -1, -5, -1)

reps = 1

print(f"Running Biopython pairwise2 localms (target=100k, query=1k) ({reps} reps)...")
ms = repeat_elapsed_ms(do_python_local, reps)
ms_per_rep = ms / reps

print("\nPython Local Search (pairwise2) benchmark")
print(f"  reps={reps}")
print(f"  python_search_ms={ms_per_rep:.4f}")
