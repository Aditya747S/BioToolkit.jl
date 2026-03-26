import random
from io import StringIO
from time import perf_counter

from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo import Consensus
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def build_records(record_count=120, width=80):
    template = ("ACGT" * ((width + 3) // 4))[:width]
    bases = "ACGT"
    records = []
    for record_index in range(1, record_count + 1):
        chars = list(template)
        for position in range(1, width + 1):
            if (record_index + position) % 17 == 0:
                chars[position - 1] = "-"
            elif (record_index * position) % 19 == 0:
                chars[position - 1] = bases[(record_index + position - 2) % len(bases)]
        records.append(SeqRecord(Seq("".join(chars)), id=f"seq{record_index}", description="bootstrap_input"))
    return records


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def bootstrap_consensus(alignment):
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator)
    trees = list(Consensus.bootstrap_trees(alignment, 25, constructor))
    consensus = Consensus.majority_consensus(trees, cutoff=0.5)
    output = StringIO()
    Phylo.write([consensus], output, "newick")
    return output.getvalue()


def main():
    random.seed(2024)
    alignment = MultipleSeqAlignment(build_records())
    consensus_ms, newick = time_call(lambda: bootstrap_consensus(alignment))

    assert "(" in newick
    assert len(alignment) == 120

    print(f"python_consensus_support_newick_ms={consensus_ms:.4f}")
    print(f"python_consensus_support_newick_ok={str(True).lower()}")


if __name__ == "__main__":
    main()
