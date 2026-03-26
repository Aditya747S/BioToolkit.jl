from copy import deepcopy
from time import perf_counter

from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def build_records(sequence_count=80, width=140):
    template = ("ACGT" * ((width + 3) // 4))[:width]
    bases = "ACGT"
    records = []
    for sequence_index in range(1, sequence_count + 1):
        chars = list(template)
        for position in range(1, width + 1):
            if (sequence_index + position) % 11 == 0:
                chars[position - 1] = bases[(sequence_index + position - 2) % len(bases)]
            elif (sequence_index * position) % 17 == 0:
                chars[position - 1] = bases[(sequence_index + 2 * position - 3) % len(bases)]
        records.append(SeqRecord(Seq("".join(chars)), id=f"seq{sequence_index}", description="prune_input"))
    return records


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def build_tree():
    alignment = MultipleSeqAlignment(build_records())
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator)
    distance_matrix = calculator.get_distance(alignment)
    return constructor.nj(distance_matrix)


def prune_tree(tree, keep_names):
    rooted = deepcopy(tree)
    keep = set(keep_names)
    for terminal in list(rooted.get_terminals()):
        if terminal.name not in keep:
            rooted.prune(terminal.name)
    return rooted


def reroot_tree(tree):
    rooted = deepcopy(tree)
    rooted.root_with_outgroup("seq1")
    return rooted


def main():
    tree = build_tree()
    keep_names = [f"seq{index}" for index in range(1, 41, 2)]

    prune_ms, pruned = time_call(lambda: prune_tree(tree, keep_names))
    reroot_ms, rerooted = time_call(lambda: reroot_tree(tree))

    assert len(pruned.get_terminals()) == len(keep_names)
    assert any(terminal.name == "seq1" for terminal in rerooted.get_terminals())

    print(f"python_prune_ms={prune_ms:.4f}")
    print(f"python_reroot_ms={reroot_ms:.4f}")
    print(f"python_prune_ok={str(True).lower()}")
    print(f"python_reroot_ok={str(True).lower()}")


if __name__ == "__main__":
    main()