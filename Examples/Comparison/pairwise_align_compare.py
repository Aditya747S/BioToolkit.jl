from time import perf_counter

from Bio.Align import PairwiseAligner


def build_sequences():
    left = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" * 30
    right = left.replace("A", "T", 180)
    return left, right


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def needleman_wunsch(left, right, match=2, mismatch=-1, gap=-2):
    left_length = len(left)
    right_length = len(right)
    scores = [[0] * (right_length + 1) for _ in range(left_length + 1)]
    trace = [[0] * (right_length + 1) for _ in range(left_length + 1)]

    for i in range(1, left_length + 1):
        scores[i][0] = scores[i - 1][0] + gap
        trace[i][0] = 2
    for j in range(1, right_length + 1):
        scores[0][j] = scores[0][j - 1] + gap
        trace[0][j] = 3

    for i in range(1, left_length + 1):
        left_byte = left[i - 1]
        for j in range(1, right_length + 1):
            right_byte = right[j - 1]
            diag = scores[i - 1][j - 1] + (match if left_byte == right_byte else mismatch)
            up = scores[i - 1][j] + gap
            left_score = scores[i][j - 1] + gap

            best = diag
            direction = 1
            if up > best:
                best = up
                direction = 2
            if left_score > best:
                best = left_score
                direction = 3

            scores[i][j] = best
            trace[i][j] = direction

    aligned_left = []
    aligned_right = []
    i = left_length
    j = right_length

    while i > 0 or j > 0:
        direction = trace[i][j]
        if direction == 1:
            aligned_left.append(left[i - 1])
            aligned_right.append(right[j - 1])
            i -= 1
            j -= 1
        elif direction == 2:
            aligned_left.append(left[i - 1])
            aligned_right.append("-")
            i -= 1
        else:
            aligned_left.append("-")
            aligned_right.append(right[j - 1])
            j -= 1

    aligned_left.reverse()
    aligned_right.reverse()
    aligned_left = "".join(aligned_left)
    aligned_right = "".join(aligned_right)
    score = scores[left_length][right_length]
    matches = sum(1 for left_byte, right_byte in zip(aligned_left, aligned_right) if left_byte == right_byte)
    identity = matches / len(aligned_left) if aligned_left else 0.0
    return aligned_left, aligned_right, score, matches, identity


def main():
    left, right = build_sequences()

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -2

    list(aligner.align(left, right))
    py_ms, alignments = time_call(lambda: list(aligner.align(left, right)))
    best = alignments[0]
    aligned_left = str(best[0])
    aligned_right = str(best[1])
    score = best.score
    matches = sum(1 for left_byte, right_byte in zip(aligned_left, aligned_right) if left_byte == right_byte)
    identity = matches / len(aligned_left) if aligned_left else 0.0

    print(f"python_pairwise_ms={py_ms:.4f}")
    print(f"python_pairwise_score={score}")
    print(f"python_pairwise_matches={matches}")
    print(f"python_pairwise_identity={identity:.6f}")
    print(f"python_pairwise_length={len(aligned_left)}")


if __name__ == "__main__":
    main()
