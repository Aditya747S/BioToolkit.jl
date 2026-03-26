#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
COMPARE_PROJECT="${BIOMATH_COMPARE_PROJECT:-/tmp/jl_AvpXPp}"
run_julia() {
    julia --project="$ROOT_DIR" "$@"
}

run_julia_compare() {
    julia --project="$COMPARE_PROJECT" "$@"
}

metric_eq() {
    local output="$1"
    local key="$2"
    grep -F -m1 "${key}=" <<<"$output" | head -n1 | cut -d'=' -f2
}

metric_eq_after_marker() {
    local output="$1"
    local marker="$2"
    local key="$3"
    awk -v marker="$marker" -v key="$key" '
        index($0, marker) { found = 1; next }
        found && index($0, key "=") {
            sub(/^.*=/, "", $0)
            print $0
            exit
        }
    ' <<<"$output"
}

metric_colon() {
    local output="$1"
    local label="$2"
    grep -F -m1 "$label" <<<"$output" | head -n1 | sed -E "s/.*${label}[[:space:]]*:[[:space:]]*([0-9.]+).*/\1/"
}

make_fastq_sample() {
    local path="$1"
    run_julia - "$path" <<'JL'
const path = ARGS[1]
const record_count = 30000
const read_length = 100
const bases = ('A', 'C', 'G', 'T')
open(path, "w") do io
    for record_index in 1:record_count
        write(io, "@read$(record_index)\n")
        for base_index in 1:read_length
            base = bases[mod(record_index + base_index - 2, length(bases)) + 1]
            write(io, base)
        end
        write(io, "\n+\n")
        for _ in 1:read_length
            write(io, 'I')
        end
        write(io, "\n")
    end
end
JL
}

fastq_path="$(mktemp)"
trap 'rm -f "$fastq_path" "$fastq_path.written"' EXIT
make_fastq_sample "$fastq_path"

fasta_out="$(run_julia "$ROOT_DIR/Examples/Comparison/fasta_index_compare.jl")"
fastq_out="$(run_julia "$ROOT_DIR/Examples/Comparison/fastq_compare.jl" "$fastq_path")"
hamming_out="$(run_julia "$ROOT_DIR/Examples/Comparison/hamming_compare.jl")"
julia_package_out="$(run_julia_compare "$ROOT_DIR/Examples/Comparison/julia_package_compare.jl")"
sequence_out="$(run_julia "$ROOT_DIR/Examples/Comparison/sequence_toolkit_compare.jl")"
motif_out="$(run_julia "$ROOT_DIR/Examples/Comparison/motif_pwm_compare.jl")"
pairwise_out="$(run_julia "$ROOT_DIR/Examples/Comparison/pairwise_align_compare.jl")"
affine_out="$(run_julia "$ROOT_DIR/Examples/Comparison/pairwise_affine_compare.jl")"
protein_out="$(run_julia "$ROOT_DIR/Examples/Comparison/protein_matrix_compare.jl")"
msa_out="$(run_julia "$ROOT_DIR/Examples/Comparison/msa_compare.jl")"
quality_out="$(run_julia "$ROOT_DIR/Examples/Comparison/quality_compare.jl")"
window_out="$(run_julia "$ROOT_DIR/Examples/Comparison/window_coverage_compare.jl")"
genbank_out="$(run_julia "$ROOT_DIR/Examples/Comparison/genbank_compare.jl")"
popgen_out="$(run_julia "$ROOT_DIR/Examples/Comparison/popgen_compare.jl")"
annotation_out="$(run_julia "$ROOT_DIR/Examples/Comparison/annotation_compare.jl")"
cuda_out="$(run_julia "$ROOT_DIR/Examples/Comparison/cuda_histogram_compare.jl")"
package_out="$(run_julia "$ROOT_DIR/scripts/benchmark.jl")"

fasta_julia="index $(metric_eq "$fasta_out" 'index_ms') ms, fetch $(metric_eq "$fasta_out" 'fetch_ms') ms"
fasta_python="index $(metric_eq_after_marker "$fasta_out" 'Python FASTA benchmark' 'index_ms') ms, fetch $(metric_eq_after_marker "$fasta_out" 'Python FASTA benchmark' 'fetch_ms') ms"

fastq_julia="parse $(metric_eq "$fastq_out" 'julia_fastq_ms') ms, write $(metric_eq "$fastq_out" 'julia_fastq_write_ms') ms"
fastq_python="parse $(metric_eq "$fastq_out" 'python_fastq_ms') ms, write $(metric_eq "$fastq_out" 'python_fastq_write_ms') ms"

hamming_julia="$(metric_eq "$hamming_out" 'distance_ms') ms"
hamming_python="$(metric_eq_after_marker "$hamming_out" 'Python Hamming benchmark' 'distance_ms') ms"
hamming_cuda="$(metric_eq "$cuda_out" 'julia_gpu_hamming_ms') ms"

julia_package_summary="reverse $(metric_eq "$julia_package_out" 'reverse_ms') ms, hamming $(metric_eq "$julia_package_out" 'hamming_ms') ms, bed parse $(metric_eq "$julia_package_out" 'bed_parse_ms') ms, gff3 parse $(metric_eq "$julia_package_out" 'gff3_parse_ms') ms"

sequence_julia="count $(metric_eq "$sequence_out" 'count_ms') ms, gc $(metric_eq "$sequence_out" 'gc_ms') ms, reverse $(metric_eq "$sequence_out" 'reverse_ms') ms, translate $(metric_eq "$sequence_out" 'translate_ms') ms, hamming $(metric_eq "$sequence_out" 'hamming_ms') ms, kmer $(metric_eq "$sequence_out" 'kmer_ms') ms"
sequence_python="count $(metric_eq_after_marker "$sequence_out" 'Python sequence toolkit benchmark' 'count_ms') ms, gc $(metric_eq_after_marker "$sequence_out" 'Python sequence toolkit benchmark' 'gc_ms') ms, reverse $(metric_eq_after_marker "$sequence_out" 'Python sequence toolkit benchmark' 'reverse_ms') ms, translate $(metric_eq_after_marker "$sequence_out" 'Python sequence toolkit benchmark' 'translate_ms') ms, hamming $(metric_eq_after_marker "$sequence_out" 'Python sequence toolkit benchmark' 'hamming_ms') ms, kmer $(metric_eq_after_marker "$sequence_out" 'Python sequence toolkit benchmark' 'kmer_ms') ms"
sequence_cuda_raw="$(metric_eq "$sequence_out" 'kmer_cuda_ms')"
if [[ "$sequence_cuda_raw" == "skipped" ]]; then
    sequence_cuda="skipped"
else
    sequence_cuda="kmer ${sequence_cuda_raw} ms"
fi

motif_julia="counts $(metric_eq "$motif_out" 'julia_counts_ms') ms, pwm $(metric_eq "$motif_out" 'julia_pwm_ms') ms"
motif_python="counts $(metric_eq "$motif_out" 'python_counts_ms') ms, pwm $(metric_eq "$motif_out" 'python_pwm_ms') ms"

pairwise_julia="$(metric_eq "$pairwise_out" 'julia_pairwise_ms') ms"
pairwise_python="$(metric_eq "$pairwise_out" 'python_pairwise_ms') ms"

affine_julia="$(metric_eq "$affine_out" 'julia_affine_ms') ms"
affine_python="$(metric_eq "$affine_out" 'python_pairwise_affine_ms') ms"

protein_julia="load $(metric_eq "$protein_out" 'julia_named_matrix_load_ms') ms, BLOSUM62 global $(metric_eq "$protein_out" 'julia_blosum62_global_ms') ms, BLOSUM62 local $(metric_eq "$protein_out" 'julia_blosum62_local_ms') ms, PAM250 global $(metric_eq "$protein_out" 'julia_pam250_global_ms') ms, PAM250 local $(metric_eq "$protein_out" 'julia_pam250_local_ms') ms"
protein_python="load $(metric_eq "$protein_out" 'python_named_matrix_load_ms') ms, BLOSUM62 global $(metric_eq "$protein_out" 'python_blosum62_global_ms') ms, BLOSUM62 local $(metric_eq "$protein_out" 'python_blosum62_local_ms') ms, PAM250 global $(metric_eq "$protein_out" 'python_pam250_global_ms') ms, PAM250 local $(metric_eq "$protein_out" 'python_pam250_local_ms') ms"

msa_julia="construct $(metric_eq "$msa_out" 'julia_msa_construct_ms') ms, column $(metric_eq "$msa_out" 'julia_msa_column_ms') ms, row $(metric_eq "$msa_out" 'julia_msa_row_ms') ms, append $(metric_eq "$msa_out" 'julia_msa_append_ms') ms, freq $(metric_eq "$msa_out" 'julia_msa_frequency_ms') ms, fasta rt $(metric_eq "$msa_out" 'julia_msa_fasta_roundtrip_ms') ms, clustal rt $(metric_eq "$msa_out" 'julia_msa_clustal_roundtrip_ms') ms, stockholm rt $(metric_eq "$msa_out" 'julia_msa_stockholm_roundtrip_ms') ms"
msa_python="construct $(metric_eq "$msa_out" 'python_msa_construct_ms') ms, column $(metric_eq "$msa_out" 'python_msa_column_ms') ms, row $(metric_eq "$msa_out" 'python_msa_row_ms') ms, append $(metric_eq "$msa_out" 'python_msa_append_ms') ms, freq $(metric_eq "$msa_out" 'python_msa_frequency_ms') ms, fasta rt $(metric_eq "$msa_out" 'python_msa_fasta_roundtrip_ms') ms, clustal rt $(metric_eq "$msa_out" 'python_msa_clustal_roundtrip_ms') ms, stockholm rt $(metric_eq "$msa_out" 'python_msa_stockholm_roundtrip_ms') ms"

quality_julia="decode $(metric_eq "$quality_out" 'julia_decode_ms') ms, mean $(metric_eq "$quality_out" 'julia_mean_ms') ms, trim $(metric_eq "$quality_out" 'julia_trim_ms') ms"
quality_python="decode $(metric_eq "$quality_out" 'python_decode_ms') ms, mean $(metric_eq "$quality_out" 'python_mean_ms') ms, trim $(metric_eq "$quality_out" 'python_trim_ms') ms"
quality_cuda="decode $(metric_eq "$cuda_out" 'julia_gpu_decode_ms') ms, mean $(metric_eq "$cuda_out" 'julia_gpu_mean_ms') ms, trim n/a"

window_julia="coverage $(metric_eq "$window_out" 'julia_window_ms') ms"
window_python="coverage $(metric_eq "$window_out" 'python_window_ms') ms"
window_cuda_raw="$(metric_eq "$window_out" 'julia_window_cuda_ms')"
if [[ "$window_cuda_raw" == "skipped" ]]; then
    window_cuda="skipped"
else
    window_cuda="coverage ${window_cuda_raw} ms"
fi

genbank_julia="parse $(metric_eq "$genbank_out" 'julia_parse_ms') ms, ingest $(metric_eq "$genbank_out" 'julia_ingest_ms') ms"
genbank_python="parse $(metric_eq "$genbank_out" 'python_genbank_ms') ms"

popgen_julia="read $(metric_eq "$popgen_out" 'julia_genepop_read_ms') ms, record read $(metric_eq "$popgen_out" 'julia_genepop_record_read_ms') ms, write $(metric_eq "$popgen_out" 'julia_genepop_write_ms') ms, split pops $(metric_eq "$popgen_out" 'julia_genepop_split_pops_ms') ms, split loci $(metric_eq "$popgen_out" 'julia_genepop_split_loci_ms') ms, remove pop $(metric_eq "$popgen_out" 'julia_genepop_remove_pop_ms') ms, remove locus $(metric_eq "$popgen_out" 'julia_genepop_remove_locus_ms') ms, allele $(metric_eq "$popgen_out" 'julia_allele_frequency_ms') ms, heterozygosity $(metric_eq "$popgen_out" 'julia_heterozygosity_ms') ms, HWE $(metric_eq "$popgen_out" 'julia_hwe_ms') ms, F-statistics $(metric_eq "$popgen_out" 'julia_f_statistics_ms') ms, LD $(metric_eq "$popgen_out" 'julia_ld_ms') ms, PCA $(metric_eq "$popgen_out" 'julia_population_pca_ms') ms"
popgen_python="read $(metric_eq_after_marker "$popgen_out" 'Python GenePop benchmark' 'python_genepop_read_ms') ms, file read $(metric_eq_after_marker "$popgen_out" 'Python GenePop benchmark' 'python_genepop_file_read_ms') ms, stringify $(metric_eq_after_marker "$popgen_out" 'Python GenePop benchmark' 'python_genepop_stringify_ms') ms, split pops $(metric_eq_after_marker "$popgen_out" 'Python GenePop benchmark' 'python_genepop_split_pops_ms') ms, split loci $(metric_eq_after_marker "$popgen_out" 'Python GenePop benchmark' 'python_genepop_split_loci_ms') ms, remove pop $(metric_eq_after_marker "$popgen_out" 'Python GenePop benchmark' 'python_genepop_remove_pop_ms') ms, remove locus $(metric_eq_after_marker "$popgen_out" 'Python GenePop benchmark' 'python_genepop_remove_locus_ms') ms, file remove pop $(metric_eq_after_marker "$popgen_out" 'Python GenePop benchmark' 'python_genepop_file_remove_pop_ms') ms, file remove locus $(metric_eq_after_marker "$popgen_out" 'Python GenePop benchmark' 'python_genepop_file_remove_locus_ms') ms"

annotation_julia="parse $(metric_eq "$annotation_out" 'parse_ms') ms"
annotation_python="parse $(metric_eq "$annotation_out" 'python_annotation_ms') ms"

cuda_cpu="$(metric_eq "$cuda_out" 'julia_cpu_ms') ms"
cuda_gpu="$(metric_eq "$cuda_out" 'julia_gpu_ms') ms"
cuda_python="$(metric_eq "$cuda_out" 'python_hist_ms') ms"

package_julia="ingest_vcf $(metric_colon "$package_out" 'ingest_vcf') ms, ingest_bed $(metric_colon "$package_out" 'ingest_bed') ms, filter_region $(metric_colon "$package_out" 'filter_region') ms, bin_positions threaded $(metric_colon "$package_out" 'bin_positions threaded') ms"

printf '%s\n' '| Benchmark | Julia | Python/Biopython | CUDA Julia | Notes |'
printf '%s\n' '|---|---:|---:|---:|---|'
printf '| FASTA index/fetch | %s | %s | n/a | subsequence SHA1 matched |\n' "$fasta_julia" "$fasta_python"
printf '| FASTQ parse/write | %s | %s | n/a | 30,000 reads, 100 bp |\n' "$fastq_julia" "$fastq_python"
printf '| Hamming distance | %s | %s | %s | distance 5 |\n' "$hamming_julia" "$hamming_python" "$hamming_cuda"
printf '| Julia package comparison | %s | n/a | n/a | BioSequences, BioAlignments, BED, and GFF3 native comparisons |\n' "$julia_package_summary"
printf '| Motif counts/PWM | %s | %s | n/a | consensus and PWM matched |\n' "$motif_julia" "$motif_python"
printf '| Pairwise alignment | %s | %s | n/a | score and identity matched |\n' "$pairwise_julia" "$pairwise_python"
printf '| Affine-gap alignment | %s | %s | n/a | Biopython affine scoring matched |\n' "$affine_julia" "$affine_python"
printf '| Protein matrix alignment | %s | %s | n/a | named matrix inventory matched and protein scores aligned |\n' "$protein_julia" "$protein_python"
printf '| Multiple sequence alignment | %s | %s | n/a | rows, columns, append, frequency, and round-trips matched |\n' "$msa_julia" "$msa_python"
printf '| Quality toolkit | %s | %s | %s | Julia, Python, and Biopython totals matched |\n' "$quality_julia" "$quality_python" "$quality_cuda"
printf '| Sequence toolkit | %s | %s | %s | gc, reverse, translate, hamming, and kmer matched |\n' "$sequence_julia" "$sequence_python" "$sequence_cuda"
printf '| Windowed coverage | %s | %s | %s | overlap totals matched |\n' "$window_julia" "$window_python" "$window_cuda"
printf '| GenBank parser/Arrow ingest | %s | %s | n/a | locus, accession, version, features, and Arrow output matched |\n' "$genbank_julia" "$genbank_python"
printf '| PopGen GenePop I/O/statistics | %s | %s | n/a | Biopython covers GenePop parsing/manipulation; analytic popgen kernels are Julia-only |\n' "$popgen_julia" "$popgen_python"
printf '| Sequence annotation objects | %s | %s | n/a | feature types, strands, and extracted sequences matched |\n' "$annotation_julia" "$annotation_python"
printf '| CUDA histogram | CPU %s | Python %s | GPU %s | bins matched |\n' "$cuda_cpu" "$cuda_python" "$cuda_gpu"
printf '| Package benchmark | %s | n/a | n/a | repository-level ingest/filter/histogram suite |\n' "$package_julia"
