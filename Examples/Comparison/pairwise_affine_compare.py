from time import perf_counter

from Bio.Align import PairwiseAligner


def build_sequences():
    left = "".join("ACGT"[(index * 7 + 3) % 4] for index in range(1, 1201))
    right = left[:500] + left[570:]
    return left, right


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def main():
    left, right = build_sequences()

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    list(aligner.align(left, right))
    py_ms, alignments = time_call(lambda: list(aligner.align(left, right)))
    best = alignments[0]
    aligned_left = str(best[0])
    aligned_right = str(best[1])
    score = best.score
    matches = sum(1 for left_byte, right_byte in zip(aligned_left, aligned_right) if left_byte == right_byte)
    identity = matches / len(aligned_left) if aligned_left else 0.0

    print(f"python_pairwise_affine_ms={py_ms:.4f}")
    print(f"python_pairwise_affine_score={score}")
    print(f"python_pairwise_affine_matches={matches}")
    print(f"python_pairwise_affine_identity={identity:.6f}")
    print(f"python_pairwise_affine_length={len(aligned_left)}")


if __name__ == "__main__":
    main()
