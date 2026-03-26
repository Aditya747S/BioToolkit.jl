from copy import deepcopy
import os
from io import StringIO
from pathlib import Path
from tempfile import mkstemp
from time import perf_counter

from Bio.PopGen import GenePop
from Bio.PopGen.GenePop import FileParser


def build_genepop_text(population_count=3, individual_count=240, locus_count=24):
    lines = ["BioInfomatics popgen benchmark\n"]
    lines.extend(f"Locus_{index}\n" for index in range(1, locus_count + 1))
    for population_index in range(1, population_count + 1):
        lines.append("Pop\n")
        for individual_index in range(1, individual_count + 1):
            markers = []
            for locus_index in range(1, locus_count + 1):
                base_value = 1 if (population_index + individual_index + locus_index) % 2 else 2
                alternate_value = 2 if base_value == 1 else 1
                if (population_index + individual_index + locus_index) % 7 == 0:
                    genotype = (base_value, alternate_value)
                elif (population_index + locus_index) % 5 == 0:
                    genotype = (alternate_value, alternate_value)
                else:
                    genotype = (base_value, base_value)
                markers.append(f"{genotype[0]:03d}{genotype[1]:03d}")
            lines.append(f"P{population_index}_I{individual_index}, {' '.join(markers)}\n")
    return "".join(lines)


def measure(func, repetitions):
    start = perf_counter()
    result = None
    for _ in range(repetitions):
        result = func()
    return ((perf_counter() - start) * 1000.0) / repetitions, result


def measure_file_operation(func, repetitions):
    start = perf_counter()
    for _ in range(repetitions):
        fd, raw_path = mkstemp(suffix=".gen")
        os.close(fd)
        output_path = Path(raw_path)
        try:
            func(output_path)
        finally:
            output_path.unlink(missing_ok=True)
    return ((perf_counter() - start) * 1000.0) / repetitions


def main():
    population_count = 3
    individual_count = 240
    locus_count = 24
    genepop_text = build_genepop_text(population_count, individual_count, locus_count)

    fd, raw_path = mkstemp(suffix=".gen")
    os.close(fd)
    temp_path = Path(raw_path)
    try:
        temp_path.write_text(genepop_text)

        with temp_path.open() as handle:
            record = GenePop.read(handle)
        file_record = FileParser.read(str(temp_path))

        parse_ms, _ = measure(lambda: GenePop.read(StringIO(genepop_text)), 10)
        file_parse_ms, _ = measure(lambda: FileParser.read(str(temp_path)), 10)
        stringify_ms, _ = measure(lambda: str(deepcopy(record)), 10)
        split_pops_ms, _ = measure(lambda: record.split_in_pops([f"pop_{index}" for index in range(len(record.populations))]), 25)
        split_loci_ms, _ = measure(lambda: record.split_in_loci(record), 25)
        remove_pop_ms, _ = measure(lambda: deepcopy(record).remove_population(0), 25)
        remove_locus_ms, _ = measure(lambda: deepcopy(record).remove_locus_by_position(0), 25)

        file_stringify_ms, _ = measure(lambda: str(FileParser.read(str(temp_path))), 10)
        file_skip_ms, _ = measure(lambda: FileParser.read(str(temp_path)).skip_population(), 25)
        file_remove_pop_ms = measure_file_operation(lambda output_path: file_record.remove_population(0, str(output_path)), 10)
        file_remove_locus_ms = measure_file_operation(lambda output_path: file_record.remove_locus_by_position(0, str(output_path)), 10)

        print("Python GenePop benchmark")
        print(f"  populations={population_count}")
        print(f"  individuals_per_pop={individual_count}")
        print(f"  loci={locus_count}")
        print(f"  python_genepop_read_ms={parse_ms:.4f}")
        print(f"  python_genepop_file_read_ms={file_parse_ms:.4f}")
        print(f"  python_genepop_stringify_ms={stringify_ms:.4f}")
        print(f"  python_genepop_split_pops_ms={split_pops_ms:.4f}")
        print(f"  python_genepop_split_loci_ms={split_loci_ms:.4f}")
        print(f"  python_genepop_remove_pop_ms={remove_pop_ms:.4f}")
        print(f"  python_genepop_remove_locus_ms={remove_locus_ms:.4f}")
        print(f"  python_genepop_file_stringify_ms={file_stringify_ms:.4f}")
        print(f"  python_genepop_skip_population_ms={file_skip_ms:.4f}")
        print(f"  python_genepop_file_remove_pop_ms={file_remove_pop_ms:.4f}")
        print(f"  python_genepop_file_remove_locus_ms={file_remove_locus_ms:.4f}")
    finally:
        temp_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()