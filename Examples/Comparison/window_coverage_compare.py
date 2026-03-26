from collections import defaultdict
from time import perf_counter


def build_intervals(n, window_size):
    starts = []
    stops = []
    for index in range(1, n + 1):
        start = ((index * 37) % 100000) + 1
        width = ((index * 11) % (window_size * 2)) + 1
        starts.append(start)
        stops.append(start + width)
    return starts, stops


def time_call(func):
    start = perf_counter()
    result = func()
    return (perf_counter() - start) * 1000.0, result


def window_coverage(starts, stops, window_size):
    coverage = defaultdict(int)
    for start, stop in zip(starts, stops):
        if stop <= start:
            continue
        first_window = start // window_size
        last_window = (stop - 1) // window_size
        for window_index in range(first_window, last_window + 1):
            window_start = window_index * window_size
            window_stop = window_start + window_size
            overlap_start = max(start, window_start)
            overlap_stop = min(stop, window_stop)
            overlap = overlap_stop - overlap_start
            if overlap > 0:
                coverage[window_index] += overlap
    return dict(coverage)


def main():
    n = 25000
    window_size = 100
    starts, stops = build_intervals(n, window_size)

    py_ms, py_cov = time_call(lambda: window_coverage(starts, stops, window_size))

    print("Window coverage benchmark")
    print(f"  intervals={n}")
    print(f"  window_size={window_size}")
    print(f"  python_window_ms={py_ms:.4f}")
    print(f"  python_window_total={sum(py_cov.values())}")


if __name__ == "__main__":
    main()
