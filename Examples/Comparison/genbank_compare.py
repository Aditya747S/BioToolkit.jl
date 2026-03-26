from sys import argv
from time import perf_counter

from Bio import SeqIO


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def main():
    if len(argv) != 2:
        raise SystemExit("usage: genbank_compare.py <genbank-path>")

    input_path = argv[1]
    list(SeqIO.parse(input_path, "genbank"))
    py_ms, records = time_call(lambda: list(SeqIO.parse(input_path, "genbank")))
    record = records[0]
    gene_features = [feature for feature in record.features if feature.type == "gene"]

    print(f"python_genbank_ms={py_ms:.4f}")
    print(f"python_genbank_locus={record.name}")
    print(f"python_genbank_accession={record.annotations['accessions'][0]}")
    print(f"python_genbank_version={record.id}")
    print(f"python_genbank_feature_count={len(record.features)}")
    print(f"python_genbank_sequence_length={len(record.seq)}")
    print(f"python_genbank_gene={gene_features[0].qualifiers['gene'][0]}")


if __name__ == "__main__":
    main()
