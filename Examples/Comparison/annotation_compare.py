from sys import argv
from time import perf_counter

from Bio import SeqIO


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def main():
    if len(argv) != 2:
        raise SystemExit("usage: annotation_compare.py <genbank-path>")

    input_path = argv[1]
    list(SeqIO.parse(input_path, "genbank"))
    py_ms, records = time_call(lambda: list(SeqIO.parse(input_path, "genbank")))
    record = records[0]
    source_feature = record.features[0]
    feature1 = record.features[1]
    feature2 = record.features[2]

    print(f"python_annotation_ms={py_ms:.4f}")
    print(f"python_annotation_locus={record.name}")
    print(f"python_annotation_feature_count={len(record.features)}")
    print(f"python_annotation_source_type={source_feature.type}")
    print(f"python_annotation_source_extract={source_feature.extract(record.seq)}")
    print(f"python_annotation_feature1_type={feature1.type}")
    print(f"python_annotation_feature1_strand={feature1.location.strand}")
    print(f"python_annotation_feature1_extract={feature1.extract(record.seq)}")
    print(f"python_annotation_feature2_type={feature2.type}")
    print(f"python_annotation_feature2_strand={feature2.location.strand}")
    print(f"python_annotation_feature2_extract={feature2.extract(record.seq)}")
    print(f"python_annotation_gene={feature1.qualifiers['gene'][0]}")


if __name__ == "__main__":
    main()
