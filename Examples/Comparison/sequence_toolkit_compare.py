from collections import Counter
from time import perf_counter

from Bio.Seq import Seq
from Bio.SeqUtils import CodonAdaptationIndex


def count_nucleotides(sequence):
    counts = Counter(sequence.upper())
    return {base: counts.get(base, 0) for base in ("A", "C", "G", "T", "N")}


def gc_content(sequence):
    sequence = sequence.upper()
    total = sum(1 for base in sequence if base in "ACGTN")
    if total == 0:
        return 0.0
    return sum(1 for base in sequence if base in "GC") / total


def transcribe_dna(sequence):
    return str(Seq(sequence).transcribe())


def hamming_distance(left, right):
    if len(left) != len(right):
        raise ValueError("sequences must have the same length")
    return sum(1 for a, b in zip(left, right) if a != b)


def kmer_frequency(sequence, k):
    sequence = sequence.upper()
    return Counter(sequence[index : index + k] for index in range(0, len(sequence) - k + 1))


def codon_usage(sequence):
    sequence = sequence.upper()
    return Counter(sequence[index : index + 3] for index in range(0, len(sequence) - 2, 3))


def codon_usage_table(sequence):
    counts = codon_usage(sequence)
    total = sum(counts.values())
    if total == 0:
        return {}
    return {codon: count / total for codon, count in counts.items()}


def elapsed_ms(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def repeat_elapsed_ms(func, repetitions):
    start = perf_counter()
    result = None
    for _ in range(repetitions):
        result = func()
    return (perf_counter() - start) * 1000.0, result


sequence = "ATGGCCCTG" * 2000
other_sequence = sequence.replace("A", "T", 25)
repetitions = 200
seq = Seq(sequence)
codon_sequence = "ATGGCCATG" * 2000
codon_reference = "ATGGCCATG" * 500
codon_reference_cai = CodonAdaptationIndex([Seq(codon_reference)])

count_ms, counts = repeat_elapsed_ms(lambda: count_nucleotides(sequence), repetitions)
gc_ms, gc = repeat_elapsed_ms(lambda: gc_content(sequence), repetitions)
transcribe_ms, transcript = repeat_elapsed_ms(lambda: transcribe_dna(sequence), repetitions)
reverse_ms, reverse = repeat_elapsed_ms(lambda: str(seq.reverse_complement()), repetitions)
translate_ms, protein = repeat_elapsed_ms(lambda: str(seq.translate(to_stop=True)), repetitions)
hamming_ms, distance = repeat_elapsed_ms(lambda: hamming_distance(sequence, other_sequence), repetitions)
kmer_ms, kmers = repeat_elapsed_ms(lambda: kmer_frequency(sequence, 3), repetitions)
codon_usage_ms, codon_counts = repeat_elapsed_ms(lambda: codon_usage(codon_sequence), repetitions)
codon_usage_table_ms, codon_table = repeat_elapsed_ms(lambda: codon_usage_table(codon_sequence), repetitions)
codon_adaptiveness_ms, codon_adaptiveness = repeat_elapsed_ms(lambda: CodonAdaptationIndex([Seq(codon_reference)]), repetitions)
cai_ms, cai_value = repeat_elapsed_ms(lambda: codon_reference_cai.calculate(codon_sequence), repetitions)

print("Python sequence toolkit benchmark")
print(f"  repetitions: {repetitions}")
print(f"  count_ms={count_ms:.4f}")
print(f"  gc_ms={gc_ms:.4f}")
print(f"  transcribe_ms={transcribe_ms:.4f}")
print(f"  reverse_ms={reverse_ms:.4f}")
print(f"  translate_ms={translate_ms:.4f}")
print(f"  hamming_ms={hamming_ms:.4f}")
print(f"  kmer_ms={kmer_ms:.4f}")
print(f"  codon_usage_ms={codon_usage_ms:.4f}")
print(f"  codon_usage_table_ms={codon_usage_table_ms:.4f}")
print(f"  codon_adaptiveness_ms={codon_adaptiveness_ms:.4f}")
print(f"  cai_ms={cai_ms:.4f}")
print(f"  count_A={counts['A']}")
print(f"  count_C={counts['C']}")
print(f"  count_G={counts['G']}")
print(f"  count_T={counts['T']}")
print(f"  count_N={counts['N']}")
print(f"  gc={gc:.6f}")
print(f"  transcribe_len={len(transcript)}")
print(f"  reverse_len={len(reverse)}")
print(f"  protein_len={len(protein)}")
print(f"  hamming_distance={distance}")
print(f"  kmer_unique={len(kmers)}")
print(f"codon_usage_ATG={codon_counts['ATG']}")
print(f"codon_usage_GCC={codon_counts['GCC']}")
print(f"codon_usage_table_ATG={codon_table['ATG']:.6f}")
print(f"codon_adaptiveness_ms={codon_adaptiveness_ms:.4f}")
print(f"cai={cai_value:.6f}")