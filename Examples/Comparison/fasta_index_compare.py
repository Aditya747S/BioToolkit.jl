from hashlib import sha1
from sys import argv
from time import perf_counter

from Bio import SeqIO


def write_sample_fasta(path):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(">seq1\n")
        handle.write("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" * 1000)
        handle.write("\n")


def elapsed_ms(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


fasta_path = argv[1] if len(argv) > 1 else "sample.fasta"
write_sample_fasta(fasta_path)

index_ms, fasta = elapsed_ms(lambda: SeqIO.index(fasta_path, "fasta"))
fetch_ms, subsequence = elapsed_ms(lambda: str(fasta["seq1"].seq[9999:11000]))
subsequence_digest = sha1(subsequence.encode("utf-8")).hexdigest()

print("Python FASTA benchmark")
print(f"  index_ms={index_ms:.4f}")
print(f"  fetch_ms={fetch_ms:.4f}")
print(f"  subsequence_len={len(subsequence)}")
print(f"  subsequence_sha1={subsequence_digest}")