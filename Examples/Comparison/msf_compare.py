from pathlib import Path
from sys import argv
from tempfile import TemporaryDirectory
from time import perf_counter

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def build_aligned_records(record_count=500, width=120):
    template = ("ACGT" * ((width + 3) // 4))[:width]
    bases = "ACGT"
    records = []
    for record_index in range(1, record_count + 1):
        chars = list(template)
        for position in range(1, width + 1):
            if (record_index + position) % 17 == 0:
                chars[position - 1] = "-"
            elif (record_index + position) % 11 == 0:
                chars[position - 1] = bases[(record_index + position - 2) % len(bases)]
        records.append(SeqRecord(Seq("".join(chars)), id=f"seq{record_index}", description="aligned_sequence"))
    return records


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def main():
    if len(argv) != 2:
        raise SystemExit("usage: msf_compare.py <msf-path>")

    input_path = Path(argv[1])
    records = build_aligned_records()

    with TemporaryDirectory() as temp_dir:
        alignment = MultipleSeqAlignment(records)
        roundtrip_ms, roundtrip_alignment_result = time_call(lambda: AlignIO.read(str(input_path), "msf"))

        assert len(alignment) == 500
        assert alignment.get_alignment_length() == 120
        assert len(roundtrip_alignment_result) == 500
        assert str(roundtrip_alignment_result[0].seq) == str(alignment[0].seq)
        assert all(expected.id == actual.id for expected, actual in zip(alignment, roundtrip_alignment_result))
        assert all(str(expected.seq) == str(actual.seq) for expected, actual in zip(alignment, roundtrip_alignment_result))

        print(f"python_msf_roundtrip_ms={roundtrip_ms:.4f}")
        print(f"python_msf_roundtrip_ok={str(True).lower()}")


if __name__ == "__main__":
    main()
