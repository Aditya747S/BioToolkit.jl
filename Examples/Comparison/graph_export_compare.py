from contextlib import redirect_stdout
from io import StringIO
from time import perf_counter

from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from copy import deepcopy


def build_records(sequence_count=40, width=100):
    template = ("ACGT" * ((width + 3) // 4))[:width]
    bases = "ACGT"
    records = []
    for sequence_index in range(1, sequence_count + 1):
        chars = list(template)
        for position in range(1, width + 1):
            if (sequence_index + position) % 13 == 0:
                chars[position - 1] = bases[(sequence_index + position - 2) % len(bases)]
            elif (sequence_index * position) % 23 == 0:
                chars[position - 1] = bases[(sequence_index + 2 * position - 3) % len(bases)]
        records.append(SeqRecord(Seq("".join(chars)), id=f"seq{sequence_index}", description="export_input"))
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


def export_ascii(tree):
    sio = StringIO()
    with redirect_stdout(sio):
        Phylo.draw_ascii(deepcopy(tree))
    return sio.getvalue()


def export_newick(tree):
    sio = StringIO()
    Phylo.write([deepcopy(tree)], sio, "newick")
    return sio.getvalue()


def main():
    tree = build_tree()

    ascii_ms, ascii_output = time_call(lambda: export_ascii(tree))
    newick_ms, newick_output = time_call(lambda: export_newick(tree))

    assert "|" in ascii_output or "+" in ascii_output
    assert newick_output.strip().endswith(";")

    print(f"python_tree_ascii_ms={ascii_ms:.4f}")
    print(f"python_tree_newick_ms={newick_ms:.4f}")
    print(f"python_tree_export_ok={str(True).lower()}")


if __name__ == "__main__":
    main()
