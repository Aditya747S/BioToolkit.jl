from sys import argv
from tempfile import NamedTemporaryFile
from time import perf_counter

from Bio import SeqIO


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def main():
    if len(argv) != 2:
        raise SystemExit("usage: genbank_write_compare.py <genbank-path>")

    input_path = argv[1]
    records = list(SeqIO.parse(input_path, "genbank"))

    with NamedTemporaryFile(mode="w+", encoding="utf-8", suffix=".gb", delete=True) as handle:
        py_ms, _ = time_call(lambda: SeqIO.write(records, handle, "genbank"))
        handle.flush()
        handle.seek(0)
        output_text = handle.read()
        handle.seek(0)
        roundtrip = list(SeqIO.parse(handle, "genbank"))

    print(f"python_genbank_write_ms={py_ms:.4f}")
    print(f"python_genbank_roundtrip_ok={str(len(roundtrip) == len(records) and roundtrip[0].name == records[0].name and roundtrip[0].annotations['accessions'][0] == records[0].annotations['accessions'][0]).lower()}")
    print(f"python_genbank_output_keywords={str('KEYWORDS' in output_text).lower()}")


if __name__ == "__main__":
    main()
