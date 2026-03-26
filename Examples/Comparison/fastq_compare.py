from sys import argv
from time import perf_counter

from Bio import SeqIO


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def parse_fastq(fastq_path):
    with open(fastq_path) as handle:
        return list(SeqIO.parse(handle, "fastq"))


def write_fastq(records, output_path):
    with open(output_path, "w") as handle:
        return SeqIO.write(records, handle, "fastq")


def main():
    if len(argv) != 2:
        raise SystemExit("usage: fastq_compare.py <fastq-path>")
    fastq_path = argv[1]
    parse_fastq(fastq_path)
    py_ms, py_records = time_call(lambda: parse_fastq(fastq_path))
    output_path = fastq_path + ".written"
    write_fastq(py_records, output_path)
    py_write_ms, py_written = time_call(lambda: write_fastq(py_records, output_path))
    print(f"python_fastq_ms={py_ms:.4f}")
    print(f"python_fastq_records={len(py_records)}")
    print(f"python_fastq_bases={sum(len(record.seq) for record in py_records)}")
    print(f"python_fastq_write_ms={py_write_ms:.4f}")
    print(f"python_fastq_write_records={py_written}")


if __name__ == "__main__":
    main()
