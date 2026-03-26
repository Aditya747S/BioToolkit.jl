from collections import Counter
from pathlib import Path
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
        records.append(
            SeqRecord(Seq("".join(chars)), id=f"seq{record_index}", description="aligned_sequence")
        )
    return records


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def roundtrip_alignment(alignment, format_name, suffix, temp_dir):
    path = Path(temp_dir) / f"alignment{suffix}"
    AlignIO.write(alignment, str(path), format_name)
    return AlignIO.read(str(path), format_name)


def assert_sequences(expected_records, actual_records):
    assert len(expected_records) == len(actual_records)
    for expected_record, actual_record in zip(expected_records, actual_records):
        assert str(expected_record.seq) == str(actual_record.seq)


def main():
    records = build_aligned_records()

    with TemporaryDirectory() as temp_dir:
        construct_ms, alignment = time_call(
            lambda: MultipleSeqAlignment(
                records,
                annotations={"source": "synthetic"},
                column_annotations={"consensus": "A" * 120, "clustal_consensus": "*" * 120, "emboss_consensus": "|" * 120},
            )
        )

        column_ms, column = time_call(lambda: alignment[:, 25])
        row_ms, row = time_call(lambda: alignment[249])

        appended = MultipleSeqAlignment(
            records,
            annotations={"source": "synthetic"},
            column_annotations={"consensus": "A" * 120, "clustal_consensus": "*" * 120, "emboss_consensus": "|" * 120},
        )
        append_ms, _ = time_call(
            lambda: appended.append(SeqRecord(Seq("ACGT" * 30), id="seq501", description="aligned sequence"))
        )
        frequency_ms, frequency = time_call(lambda: {base: count / len(alignment) for base, count in Counter(str(alignment[:, 25])).items()})

        fasta_roundtrip_ms, fasta_roundtrip_alignment = time_call(
            lambda: roundtrip_alignment(alignment, "fasta", ".fasta", temp_dir)
        )
        clustal_roundtrip_ms, clustal_roundtrip_alignment = time_call(
            lambda: roundtrip_alignment(alignment, "clustal", ".aln", temp_dir)
        )
        stockholm_roundtrip_ms, stockholm_roundtrip_alignment = time_call(
            lambda: roundtrip_alignment(alignment, "stockholm", ".sto", temp_dir)
        )

        assert len(alignment) == 500
        assert alignment.get_alignment_length() == 120
        assert len(column) == 500
        assert row.id == "seq250"
        assert len(appended) == 501
        assert frequency["A"] > 0.0
        assert_sequences(alignment, fasta_roundtrip_alignment)
        assert_sequences(alignment, clustal_roundtrip_alignment)
        assert_sequences(alignment, stockholm_roundtrip_alignment)
        assert all(expected.id == actual.id for expected, actual in zip(alignment, fasta_roundtrip_alignment))
        assert all(expected.id == actual.id for expected, actual in zip(alignment, clustal_roundtrip_alignment))
        assert all(expected.id == actual.id for expected, actual in zip(alignment, stockholm_roundtrip_alignment))
        assert fasta_roundtrip_alignment.annotations == {}
        assert clustal_roundtrip_alignment.annotations == {}
        assert stockholm_roundtrip_alignment.annotations == {}
        assert fasta_roundtrip_alignment.column_annotations == {}
        assert clustal_roundtrip_alignment.column_annotations == {"clustal_consensus": "*" * 120}
        assert stockholm_roundtrip_alignment.column_annotations == {}

        print(f"python_msa_construct_ms={construct_ms:.4f}")
        print(f"python_msa_column_ms={column_ms:.4f}")
        print(f"python_msa_row_ms={row_ms:.4f}")
        print(f"python_msa_append_ms={append_ms:.4f}")
        print(f"python_msa_frequency_ms={frequency_ms:.4f}")
        print(f"python_msa_fasta_roundtrip_ms={fasta_roundtrip_ms:.4f}")
        print(f"python_msa_clustal_roundtrip_ms={clustal_roundtrip_ms:.4f}")
        print(f"python_msa_stockholm_roundtrip_ms={stockholm_roundtrip_ms:.4f}")
        print(f"python_msa_rows={len(alignment)}")
        print(f"python_msa_columns={alignment.get_alignment_length()}")
        print("python_msa_fasta_roundtrip_ok=1")
        print("python_msa_clustal_roundtrip_ok=1")
        print("python_msa_stockholm_roundtrip_ok=1")


if __name__ == "__main__":
    main()
