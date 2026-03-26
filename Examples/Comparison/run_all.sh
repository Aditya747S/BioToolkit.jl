#!/bin/bash
cd "/home/aditya-sharma/Coding/Computation Biology/Examples/Comparison"
echo "--- Julia Benchmarks ---" > all_results.txt
for f in clustal_muscle_compare.jl consensus_support_newick_compare.jl graph_export_compare.jl midpoint_root_compare.jl motif_scan_compare.jl msf_compare.jl phylip_compare.jl popgen_compare.jl protein_matrix_compare.jl prune_reroot_compare.jl cuda_histogram_compare.jl; do
    echo "==== $f ====" >> all_results.txt
    julia --project=../../ $f >> all_results.txt 2>&1
done
echo "--- Python Benchmarks ---" >> all_results.txt
for f in clustal_muscle_compare.py consensus_support_newick_compare.py graph_export_compare.py midpoint_root_compare.py motif_scan_compare.py msf_compare.py phylip_compare.py protein_matrix_compare.py prune_reroot_compare.py; do
    echo "==== $f ====" >> all_results.txt
    python3 $f >> all_results.txt 2>&1
done
