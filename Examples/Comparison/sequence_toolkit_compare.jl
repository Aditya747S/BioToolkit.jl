pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

const has_cuda_package = Base.find_package("CUDA") !== nothing
if has_cuda_package
    @eval using CUDA
end

sequence = repeat("ATGGCCCTG", 2_000)
other_sequence = replace(sequence, 'A' => 'T'; count=25)
repetitions = 200
sequence_bytes = Vector{UInt8}(codeunits(sequence))
translate_buffer = Vector{UInt8}(undef, length(sequence_bytes) ÷ 3)
codon_sequence = repeat("ATGGCCATG", 2_000)
codon_reference = repeat("ATGGCCATG", 500)

function main()
    counts = BioToolkit.count_nucleotides(sequence)
    gc_value = BioToolkit.gc_content(sequence)
    transcript = BioToolkit.transcribe_dna(sequence)
    reverse = BioToolkit.reverse_complement(sequence)
    BioToolkit.translate_dna!(translate_buffer, sequence_bytes; stop_at_stop=true)
    distance = BioToolkit.hamming_distance(sequence, other_sequence)
    kmers = BioToolkit.kmer_frequency(sequence, 3)
    codon_usage = BioToolkit.codon_usage(codon_sequence)
    codon_usage_table = BioToolkit.codon_usage_table(codon_sequence)
    codon_adaptiveness = BioToolkit.relative_codon_adaptiveness(codon_reference)
    cai_value = BioToolkit.codon_adaptation_index(codon_sequence; reference=codon_reference)

    count_ms, counts = repeat_elapsed_ms(() -> BioToolkit.count_nucleotides(sequence), repetitions)
    gc_ms, gc = repeat_elapsed_ms(() -> BioToolkit.gc_content(sequence), repetitions)
    transcribe_ms, transcript = repeat_elapsed_ms(() -> BioToolkit.transcribe_dna(sequence), repetitions)
    reverse_ms, reverse = repeat_elapsed_ms(() -> BioToolkit.reverse_complement(sequence), repetitions)
    translate_ms, protein = repeat_elapsed_ms(() -> begin
        protein_length = BioToolkit.translate_dna!(translate_buffer, sequence_bytes; stop_at_stop=true)
        unsafe_string(pointer(translate_buffer), protein_length)
    end, repetitions)
    hamming_ms, distance = repeat_elapsed_ms(() -> BioToolkit.hamming_distance(sequence, other_sequence), repetitions)
    kmer_ms, kmers = repeat_elapsed_ms(() -> BioToolkit.kmer_frequency(sequence, 3), repetitions)
    codon_usage_ms, codon_usage = repeat_elapsed_ms(() -> BioToolkit.codon_usage(codon_sequence), repetitions)
    codon_usage_table_ms, codon_usage_table = repeat_elapsed_ms(() -> BioToolkit.codon_usage_table(codon_sequence), repetitions)
    codon_adaptiveness_ms, codon_adaptiveness = repeat_elapsed_ms(() -> BioToolkit.relative_codon_adaptiveness(codon_reference), repetitions)
    cai_ms, cai_value = repeat_elapsed_ms(() -> BioToolkit.codon_adaptation_index(codon_sequence; reference=codon_reference), repetitions)

    kmer_cuda_ms = missing
    kmer_cuda_unique = missing
    if has_cuda_package && CUDA.functional()
        gpu_sequence = CuArray(sequence_bytes)
        CUDA.@sync BioToolkit.kmer_frequency(gpu_sequence, 3)
        kmer_cuda_ms, gpu_kmers = elapsed_ms() do
            CUDA.@sync BioToolkit.kmer_frequency(gpu_sequence, 3)
        end
        @assert gpu_kmers == kmers
        kmer_cuda_unique = length(gpu_kmers)
    end

    python_output = read(`conda run -n general --no-capture-output python $(joinpath(@__DIR__, "sequence_toolkit_compare.py"))`, String)

    python_count_a = parse(Int, match(r"count_A=(\d+)", python_output).captures[1])
    python_count_c = parse(Int, match(r"count_C=(\d+)", python_output).captures[1])
    python_count_g = parse(Int, match(r"count_G=(\d+)", python_output).captures[1])
    python_count_t = parse(Int, match(r"count_T=(\d+)", python_output).captures[1])
    python_count_n = parse(Int, match(r"count_N=(\d+)", python_output).captures[1])
    python_gc = parse(Float64, match(r"gc=([0-9.]+)", python_output).captures[1])
    python_transcribe_len = parse(Int, match(r"transcribe_len=(\d+)", python_output).captures[1])
    python_reverse_len = parse(Int, match(r"reverse_len=(\d+)", python_output).captures[1])
    python_protein_len = parse(Int, match(r"protein_len=(\d+)", python_output).captures[1])
    python_hamming_distance = parse(Int, match(r"hamming_distance=(\d+)", python_output).captures[1])
    python_kmer_unique = parse(Int, match(r"kmer_unique=(\d+)", python_output).captures[1])
    python_codon_usage_atg = parse(Int, match(r"codon_usage_ATG=(\d+)", python_output).captures[1])
    python_codon_usage_gcc = parse(Int, match(r"codon_usage_GCC=(\d+)", python_output).captures[1])
    python_codon_usage_table_atg = parse(Float64, match(r"codon_usage_table_ATG=([0-9.]+)", python_output).captures[1])
    python_codon_adaptiveness_ms = parse(Float64, match(r"codon_adaptiveness_ms=([0-9.]+)", python_output).captures[1])
    python_cai_ms = parse(Float64, match(r"cai_ms=([0-9.]+)", python_output).captures[1])
    python_cai = parse(Float64, match(r"cai=([0-9.]+)", python_output).captures[1])

    @assert counts == (A=python_count_a, C=python_count_c, G=python_count_g, T=python_count_t, N=python_count_n)
    @assert round(gc_value, digits=6) == python_gc
    @assert length(transcript) == python_transcribe_len
    @assert length(reverse) == python_reverse_len
    @assert length(protein) == python_protein_len
    @assert distance == python_hamming_distance
    @assert length(kmers) == python_kmer_unique
    @assert codon_usage["ATG"] == python_codon_usage_atg
    @assert codon_usage["GCC"] == python_codon_usage_gcc
    @assert isapprox(codon_usage_table["ATG"], python_codon_usage_table_atg; atol=1e-6)
    @assert isapprox(cai_value, python_cai; atol=1e-12)

    println("Julia sequence toolkit benchmark")
    println("  repetitions: ", repetitions)
    println("  count_ms=", round(count_ms, digits=4))
    println("  gc_ms=", round(gc_ms, digits=4))
    println("  transcribe_ms=", round(transcribe_ms, digits=4))
    println("  reverse_ms=", round(reverse_ms, digits=4))
    println("  translate_ms=", round(translate_ms, digits=4))
    println("  hamming_ms=", round(hamming_ms, digits=4))
    println("  kmer_ms=", round(kmer_ms, digits=4))
    println("  codon_usage_ms=", round(codon_usage_ms, digits=4))
    println("  codon_usage_table_ms=", round(codon_usage_table_ms, digits=4))
    println("  codon_adaptiveness_ms=", round(codon_adaptiveness_ms, digits=4))
    println("  cai_ms=", round(cai_ms, digits=4))
    if kmer_cuda_ms === missing
        println("  kmer_cuda_ms=skipped")
        println("  kmer_cuda_unique=skipped")
    else
        println("  kmer_cuda_ms=", round(kmer_cuda_ms, digits=4))
        println("  kmer_cuda_unique=", kmer_cuda_unique)
    end
    println("  gc=", round(gc_value, digits=6))
    println("  transcribe_len=", length(transcript))
    println("  reverse_len=", length(reverse))
    println("  protein_len=", length(protein))
    println("  hamming_distance=", distance)
    println("  kmer_unique=", length(kmers))
    println("  codon_usage_ATG=", codon_usage["ATG"])
    println("  codon_usage_GCC=", codon_usage["GCC"])
    println("  codon_usage_table_ATG=", round(codon_usage_table["ATG"], digits=6))
    println("  codon_adaptiveness_ms=", round(codon_adaptiveness_ms, digits=4))
    println("  cai=", round(cai_value, digits=6))
    print(python_output)
end

main()