from time import perf_counter
from io import StringIO

from Bio import motifs
from Bio.Seq import Seq


def build_sequences(count, motif_length):
    template = list("ACGTACGTACGT")
    bases = "ACGT"
    if motif_length != len(template):
        raise ValueError(f"motif_length must equal {len(template)}")
    instances = []
    for index in range(1, count + 1):
        sequence = template.copy()
        if index % 7 == 0:
            sequence[3] = bases[index % len(bases)]
        if index % 11 == 0:
            sequence[8] = bases[(index + 1) % len(bases)]
        instances.append(Seq("".join(sequence)))
    return instances


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def main():
    record_count = 20000
    motif_length = 12
    instances = build_sequences(record_count, motif_length)
    jaspar_data = ">MA0001.1 Example motif\nA [ 3 1 0 0 ]\nC [ 0 2 4 0 ]\nG [ 1 1 0 5 ]\nT [ 0 0 0 0 ]\n"

    motif = motifs.create(instances)
    motif.counts
    motif.pwm
    list(motifs.parse(StringIO(jaspar_data), "jaspar"))

    counts_ms, motif_object = time_call(lambda: motifs.create(instances))
    pwm_ms, pwm = time_call(lambda: motif_object.pwm)
    jaspar_ms, jaspar_motif_list = time_call(lambda: list(motifs.parse(StringIO(jaspar_data), "jaspar")))
    jaspar_motif = jaspar_motif_list[0]

    print(f"python_counts_ms={counts_ms:.4f}")
    print(f"python_pwm_ms={pwm_ms:.4f}")
    print(f"python_jaspar_ms={jaspar_ms:.4f}")
    print(f"python_consensus={motif_object.consensus}")
    print(f"python_pwm_first={pwm['A'][0]:.6f}")
    print(f"python_jaspar_name={jaspar_motif.name}")
    print(f"python_jaspar_consensus={jaspar_motif.consensus}")
    print(f"python_jaspar_pwm_first={jaspar_motif.pwm['A'][0]:.6f}")


if __name__ == "__main__":
    main()
