from time import perf_counter


def hamming_distance(left, right):
    if len(left) != len(right):
        raise ValueError("sequences must have the same length")
    return sum(1 for a, b in zip(left, right) if a != b)


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


left = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG" * 500
right = left.replace("A", "T", 5)

repetitions = 1000

distance_ms, distance = repeat_elapsed_ms(lambda: hamming_distance(left, right), repetitions)

print("Python Hamming benchmark")
print(f"  repetitions: {repetitions}")
print(f"  distance_ms={distance_ms:.4f}")
print(f"  distance={distance}")