from io import StringIO
from time import perf_counter

from Bio import SeqIO


def build_records(record_count, read_length):
    bases = "ACGT"
    high_quality = "I" * (read_length - 12)
    low_quality = "!" * 12
    records = []
    for index in range(1, record_count + 1):
        sequence = "".join(bases[(index + position - 2) % len(bases)] for position in range(1, read_length + 1))
        quality = high_quality + low_quality if index % 5 == 0 else "I" * read_length
        records.append((f"read{index}", sequence, quality))
    return records


def records_to_fastq_text(records):
    lines = []
    for identifier, sequence, quality in records:
        lines.append(f"@{identifier}\n")
        lines.append(f"{sequence}\n")
        lines.append("+\n")
        lines.append(f"{quality}\n")
    return "".join(lines)


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def phred_scores(quality, offset=33):
    if quality and isinstance(quality[0], int):
        return list(quality)
    return [ord(character) - offset for character in quality]


def mean_quality(quality, offset=33):
    scores = phred_scores(quality, offset=offset)
    return sum(scores) / len(scores) if scores else 0.0


def phred_string(scores, offset=33):
    return "".join(chr(score + offset) for score in scores)


def trim_low_quality(record, window=8, threshold=25, offset=33):
    identifier, sequence, quality = record
    scores = phred_scores(quality, offset=offset)
    if len(scores) < window:
        return None

    first_keep = None
    last_keep = None
    threshold_sum = threshold * window
    window_sum = sum(scores[:window])

    for start in range(0, len(scores) - window + 1):
        if start > 0:
            window_sum += scores[start + window - 1] - scores[start - 1]
        if window_sum >= threshold_sum:
            if first_keep is None:
                first_keep = start
            last_keep = start + window - 1

    if first_keep is None:
        return None

    return identifier, sequence[first_keep:last_keep + 1], quality[first_keep:last_keep + 1]


def parse_biopython_records(records):
    fastq_text = records_to_fastq_text(records)
    return list(SeqIO.parse(StringIO(fastq_text), "fastq"))


def main():
    record_count = 20000
    read_length = 100
    records = build_records(record_count, read_length)
    biopython_records = parse_biopython_records(records)
    score_vector = list(range(40))

    phred_scores(biopython_records[0].letter_annotations["phred_quality"])
    mean_quality(records[0][2])
    phred_string(score_vector)
    trim_low_quality(records[0])

    decode_ms, decode_total = time_call(lambda: sum(len(phred_scores(quality)) for _, _, quality in records))
    mean_ms, mean_total = time_call(lambda: sum(mean_quality(quality) for _, _, quality in records))
    string_ms, string_value = time_call(lambda: phred_string(score_vector))
    trim_ms, trim_total = time_call(lambda: sum(len(trimmed[1]) if trimmed is not None else 0 for trimmed in (trim_low_quality(record) for record in records)))

    bio_decode_total = sum(len(record.letter_annotations["phred_quality"]) for record in biopython_records)
    bio_mean_total = sum(sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"]) for record in biopython_records)
    bio_trim_total = 0
    for record in biopython_records:
        quality_string = "".join(chr(score + 33) for score in record.letter_annotations["phred_quality"])
        trimmed = trim_low_quality((record.id, str(record.seq), quality_string))
        bio_trim_total += len(trimmed[1]) if trimmed is not None else 0

    print("Quality toolkit benchmark")
    print(f"  records={record_count}")
    print(f"  read_length={read_length}")
    print(f"  python_decode_ms={decode_ms:.4f}")
    print(f"  python_decode_total={decode_total}")
    print(f"  python_mean_ms={mean_ms:.4f}")
    print(f"  python_mean_total={mean_total:.4f}")
    print(f"  python_phred_string_ms={string_ms:.4f}")
    print(f"  python_phred_string_len={len(string_value)}")
    print(f"  python_trim_ms={trim_ms:.4f}")
    print(f"  python_trim_total={trim_total}")
    print(f"  biopython_decode_total={bio_decode_total}")
    print(f"  biopython_mean_total={bio_mean_total:.4f}")
    print(f"  biopython_trim_total={bio_trim_total}")


if __name__ == "__main__":
    main()
