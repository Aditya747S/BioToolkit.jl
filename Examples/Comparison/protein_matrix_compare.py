from time import perf_counter

from Bio.Align import PairwiseAligner, substitution_matrices


def build_sequences():
    left = "MTEITAAMVKELRESTGAGMMDCKNALSETQHEWAY" * 12
    right_chars = list(left)
    for index in range(0, len(right_chars), 3):
        if right_chars[index] == "A":
            right_chars[index] = "G"
        elif right_chars[index] == "G":
            right_chars[index] = "A"
        elif right_chars[index] == "M":
            right_chars[index] = "L"
    right = "".join(right_chars)
    return left, right


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def align(sequence_left, sequence_right, matrix_name, mode):
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load(matrix_name)
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    alignments = list(aligner.align(sequence_left, sequence_right))
    best = alignments[0]
    aligned_left = str(best[0])
    aligned_right = str(best[1])
    matches = sum(1 for left_byte, right_byte in zip(aligned_left, aligned_right) if left_byte == right_byte)
    identity = matches / len(aligned_left) if aligned_left else 0.0
    return best.score, matches, identity, len(aligned_left)


def main():
    sequence_left, sequence_right = build_sequences()

    names = list(substitution_matrices.load())
    alignable_names = [name for name in names if name not in {"GENETIC", "SCHNEIDER", "TRANS"}]
    load_ms, _ = time_call(lambda: [substitution_matrices.load(name) for name in alignable_names])

    blosum62_global_ms, blosum62_global = time_call(lambda: align(sequence_left, sequence_right, "BLOSUM62", "global"))
    blosum62_local_ms, blosum62_local = time_call(lambda: align(sequence_left, sequence_right, "BLOSUM62", "local"))
    pam250_global_ms, pam250_global = time_call(lambda: align(sequence_left, sequence_right, "PAM250", "global"))
    pam250_local_ms, pam250_local = time_call(lambda: align(sequence_left, sequence_right, "PAM250", "local"))

    print("python_named_matrix_inventory=" + ",".join(names))
    print(f"python_alignable_matrix_count={len(alignable_names)}")
    print(f"python_named_matrix_load_ms={load_ms:.4f}")
    print(f"python_blosum62_global_ms={blosum62_global_ms:.4f}")
    print(f"python_blosum62_global_score={blosum62_global[0]}")
    print(f"python_blosum62_global_identity={blosum62_global[2]:.6f}")
    print(f"python_blosum62_global_length={blosum62_global[3]}")
    print(f"python_blosum62_local_ms={blosum62_local_ms:.4f}")
    print(f"python_blosum62_local_score={blosum62_local[0]}")
    print(f"python_blosum62_local_identity={blosum62_local[2]:.6f}")
    print(f"python_blosum62_local_length={blosum62_local[3]}")
    print(f"python_pam250_global_ms={pam250_global_ms:.4f}")
    print(f"python_pam250_global_score={pam250_global[0]}")
    print(f"python_pam250_global_identity={pam250_global[2]:.6f}")
    print(f"python_pam250_global_length={pam250_global[3]}")
    print(f"python_pam250_local_ms={pam250_local_ms:.4f}")
    print(f"python_pam250_local_score={pam250_local[0]}")
    print(f"python_pam250_local_identity={pam250_local[2]:.6f}")
    print(f"python_pam250_local_length={pam250_local[3]}")


if __name__ == "__main__":
    main()
