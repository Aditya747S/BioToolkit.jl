import sys
import time
import random
import os

try:
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
except ImportError:
    print("python_nj_ok=false")
    print("python_nj_ms=0.0")
    sys.exit(0)

def build_distance_matrix(n=200):
    random.seed(2024)
    names = [f"Taxon_{i+1}" for i in range(n)]
    matrix = []
    
    # Biopython requires lower triangular matrix for DistanceMatrix
    for i in range(n):
        row = []
        for j in range(i+1):
            if i == j:
                row.append(0.0)
            else:
                row.append(random.random())
        matrix.append(row)
        
    return DistanceMatrix(names=names, matrix=matrix)

def repeat_elapsed_ms(func, reps=10):
    start_time = time.perf_counter()
    for _ in range(reps):
        func()
    end_time = time.perf_counter()
    return ((end_time - start_time) / reps) * 1000

def main():
    dm = build_distance_matrix(200)
    
    constructor = DistanceTreeConstructor()
    
    # Benchmark
    nj_ms = repeat_elapsed_ms(lambda: constructor.nj(dm), 10)
    
    print(f"python_nj_ms={nj_ms:.4f}")
    print("python_nj_ok=true")

if __name__ == "__main__":
    main()
