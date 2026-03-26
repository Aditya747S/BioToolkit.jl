from copy import deepcopy
from io import StringIO
from time import perf_counter

from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def build_sequences(sequence_count=60, width=160):
    template = ("ACGT" * ((width + 3) // 4))[:width]
    bases = "ACGT"
    records = []
    for sequence_index in range(1, sequence_count + 1):
        chars = list(template)
        for position in range(1, width + 1):
            if (sequence_index + position) % 13 == 0:
                chars[position - 1] = bases[(sequence_index + position - 2) % len(bases)]
            elif (sequence_index * position) % 29 == 0:
                chars[position - 1] = bases[(sequence_index + 2 * position - 3) % len(bases)]
        records.append(SeqRecord(Seq("".join(chars)), id=f"seq{sequence_index}", description="synthetic"))
    return records


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def build_tree():
    alignment = MultipleSeqAlignment(build_sequences())
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor()
    distance_matrix = calculator.get_distance(alignment)
    return constructor.nj(distance_matrix)


def main():
    tree = build_tree()

    def midpoint_root():
        rooted = deepcopy(tree)
        rooted.root_at_midpoint()
        return rooted

    midpoint_ms, rooted = time_call(midpoint_root)

    assert len(rooted.get_terminals()) == len(tree.get_terminals())

    print(f"python_midpoint_root_ms={midpoint_ms:.4f}")
    print(f"python_midpoint_root_ok={str(True).lower()}")


if __name__ == "__main__":
    main()
