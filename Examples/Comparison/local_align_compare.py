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

seq1 = "ACGTGCCG" * 125 # 1000 bp
seq2 = "ACGTGCG" * 100  # 700 bp

def do_local():
    _ = pairwise2.align.localms(seq1, seq2, 2, -1, -5, -1)

reps = 10
do_local() # heat

ms = repeat_elapsed_ms(do_local, reps)
ms_per_rep = ms / reps
print("Python Local Align benchmark")
print(f"  reps={reps}")
print(f"  python_local_align_ms={ms_per_rep:.4f}")
