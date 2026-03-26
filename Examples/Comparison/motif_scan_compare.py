from time import perf_counter

from Bio import motifs
from Bio.Seq import Seq


def build_motif_inputs():
    instances = [Seq("TTTACGTTTT"), Seq("GGGACGTGGG"), Seq("CCCACGTCCC"), Seq("TTTACGTTTT")]
    sequence = "TTTACGTTTTGGGACGTGGGCCCACGTCCC" * 1000
    return instances, sequence


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def main():
    instances, sequence = build_motif_inputs()
    motif = motifs.create(instances)
    motif.pseudocounts = 0.5
    pssm = motif.pssm

    scan_ms, hits = time_call(lambda: list(pssm.search(Seq(sequence), threshold=0.0, both=True)))

    print(f"python_motif_scan_ms={scan_ms:.4f}")
    print(f"python_motif_hit_count={len(hits)}")
    print(f"python_motif_first_score={hits[0][1]:.6f}" if hits else "python_motif_first_score=0.000000")
    print(f"python_motif_scan_ok={str(True).lower()}")


if __name__ == "__main__":
    main()
