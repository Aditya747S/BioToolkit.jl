from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from time import perf_counter
import os
import sys

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def build_sequences(record_count=8, width=180):
    template = ("ACGT" * ((width + 3) // 4))[:width]
    bases = "ACGT"
    records = []
    for record_index in range(1, record_count + 1):
        chars = list(template)
        for position in range(1, width + 1):
            if (record_index + position) % 19 == 0:
                chars[position - 1] = bases[(record_index + position - 2) % len(bases)]
            elif (record_index + position) % 13 == 0:
                chars[position - 1] = bases[(record_index + position) % len(bases)]
        records.append(SeqRecord(Seq("".join(chars)), id=f"seq{record_index}", description="synthetic_sequence"))
    return records


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def resolve_binary(candidate):
    candidates = [candidate]
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.append(str(Path(conda_prefix) / "bin" / candidate))
    home = Path.home()
    for root in [home / "miniconda3" / "envs" / "general", home / ".conda" / "envs" / "general", home / "anaconda3" / "envs" / "general"]:
        candidates.append(str(root / "bin" / candidate))
    for path in candidates:
        if Path(path).is_file():
            return path
    raise FileNotFoundError(candidate)


def write_fasta(records, path):
    AlignIO.write(MultipleSeqAlignment(records), str(path), "fasta")


def run_clustalo(binary, input_path, output_path):
    run([binary, "--infile", str(input_path), "--outfile", str(output_path), "--outfmt", "fasta", "--force"], check=True, capture_output=True)


def run_muscle(binary, input_path, output_path):
    attempts = [
        [binary, "-align", str(input_path), "-output", str(output_path)],
        [binary, "-in", str(input_path), "-out", str(output_path)],
    ]
    last_error = None
    for command in attempts:
        try:
            run(command, check=True, capture_output=True)
            return
        except Exception as exc:
            last_error = exc
    raise last_error


def roundtrip(binary, runner, records, suffix, temp_dir):
    input_path = Path(temp_dir) / f"input{suffix}.fasta"
    output_path = Path(temp_dir) / f"output{suffix}.fasta"
    write_fasta(records, input_path)
    runner(binary, input_path, output_path)
    return AlignIO.read(str(output_path), "fasta")


def assert_alignment(alignment, expected_count):
    assert len(alignment) == expected_count
    assert alignment.get_alignment_length() > 0
    assert all(str(record.seq) for record in alignment)


def main():
    records = build_sequences()
    clustalo = resolve_binary("clustalo")
    muscle = resolve_binary("muscle")

    with TemporaryDirectory() as temp_dir:
        clustal_ms, clustal_alignment = time_call(lambda: roundtrip(clustalo, run_clustalo, records, "clustal", temp_dir))
        muscle_ms, muscle_alignment = time_call(lambda: roundtrip(muscle, run_muscle, records, "muscle", temp_dir))

        assert_alignment(clustal_alignment, len(records))
        assert_alignment(muscle_alignment, len(records))

        print(f"python_clustal_msa_ms={clustal_ms:.4f}")
        print(f"python_muscle_msa_ms={muscle_ms:.4f}")
        print(f"python_external_msa_records={len(clustal_alignment)}")
        print(f"python_external_msa_columns={clustal_alignment.get_alignment_length()}")
        print("python_clustal_msa_ok=1")
        print("python_muscle_msa_ok=1")


if __name__ == "__main__":
    main()
