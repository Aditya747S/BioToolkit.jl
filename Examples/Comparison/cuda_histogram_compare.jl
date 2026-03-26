pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit

include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

const has_cuda_package = Base.find_package("CUDA") !== nothing
if has_cuda_package
    @eval using CUDA
end

function build_positions(n::Integer, max_position::Integer)
    positions = Vector{Int32}(undef, n)
    @inbounds for index in 1:n
        positions[index] = Int32(mod(index * 97, max_position) + 1)
    end
    return positions
end

function cpu_histogram(positions::Vector{Int32}, bin_size::Integer, bin_count::Integer)
    bins = zeros(Int32, bin_count)
    @inbounds for position in positions
        bin_index = fld(Int(position) - 1, bin_size) + 1
        bins[bin_index] += 1
    end
    return bins
end

function write_python_benchmark_script(path::AbstractString, n::Integer, max_position::Integer, bin_size::Integer)
    open(path, "w") do io
        write(io, "from collections import Counter\n")
        write(io, "from time import perf_counter\n")
        write(io, "n = $(n)\n")
        write(io, "max_position = $(max_position)\n")
        write(io, "bin_size = $(bin_size)\n")
        write(io, "positions = [((index * 97) % max_position) + 1 for index in range(1, n + 1)]\n")
        write(io, "def time_call(func):\n")
        write(io, "    start = perf_counter()\n")
        write(io, "    result = func()\n")
        write(io, "    return (perf_counter() - start) * 1000.0, result\n")
        write(io, "def histogram(values, bin_size):\n")
        write(io, "    counts = Counter()\n")
        write(io, "    for position in values:\n")
        write(io, "        counts[(position - 1) // bin_size] += 1\n")
        write(io, "    return counts\n")
        write(io, "py_ms, py_hist = time_call(lambda: histogram(positions, bin_size))\n")
        write(io, "print(f'python_hist_ms={py_ms:.4f}')\n")
        write(io, "print(f'python_bins={len(py_hist)}')\n")
    end
end

function gpu_histogram_kernel!(bins, positions, bin_size::Int32)
    index = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
    if index <= length(positions)
        bin_index = Int(fld(Int(positions[index]) - 1, Int(bin_size))) + 1
        CUDA.@atomic bins[bin_index] += Int32(1)
    end
    return
end

function main()
    n = 2_000_000
    max_position = 2_000_000
    bin_size = 128
    bin_count = cld(max_position, bin_size)
    positions = build_positions(n, max_position)
    quality_bytes = fill(UInt8('I'), 2_000_000)
    left = repeat("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 500)
    right = replace(left, 'A' => 'T'; count=5)

    cpu_ms, cpu_hist = repeat_elapsed_ms(() -> cpu_histogram(positions, bin_size, bin_count), 1)

    gpu_ms = missing
    gpu_hist = nothing
    gpu_hamming_ms = missing
    gpu_hamming = missing
    gpu_decode_ms = missing
    gpu_decode_total = missing
    gpu_decoded_scores = nothing
    gpu_mean_ms = missing
    gpu_mean = missing
    if has_cuda_package && CUDA.functional()
        gpu_positions = CuArray(positions)
        gpu_bins = CUDA.zeros(Int32, bin_count)
        gpu_left = CuArray(collect(codeunits(left)))
        gpu_right = CuArray(collect(codeunits(right)))
        gpu_quality = CuArray(quality_bytes)
        threads = 256
        blocks = cld(length(gpu_positions), threads)

        cpu_hist = cpu_histogram(positions, bin_size, bin_count)
        cpu_hamming = BioToolkit.hamming_distance(left, right)
        cpu_quality_scores = BioToolkit.phred_scores(quality_bytes)
        cpu_quality_mean = BioToolkit.mean_quality(quality_bytes)

        # warm up the GPU kernel
        CUDA.fill!(gpu_bins, 0)
        CUDA.@sync @cuda threads=threads blocks=blocks gpu_histogram_kernel!(gpu_bins, gpu_positions, Int32(bin_size))
        CUDA.@sync BioToolkit.hamming_distance(gpu_left, gpu_right)
        CUDA.@sync BioToolkit.phred_scores(gpu_quality)
        CUDA.@sync BioToolkit.mean_quality(gpu_quality)

        CUDA.fill!(gpu_bins, 0)
        gpu_ms, _ = elapsed_ms() do
            CUDA.@sync @cuda threads=threads blocks=blocks gpu_histogram_kernel!(gpu_bins, gpu_positions, Int32(bin_size))
        end

        gpu_hamming_ms, gpu_hamming = elapsed_ms() do
            CUDA.@sync BioToolkit.hamming_distance(gpu_left, gpu_right)
        end
        @assert gpu_hamming == cpu_hamming

        gpu_decode_ms, gpu_decoded_scores = elapsed_ms() do
            CUDA.@sync BioToolkit.phred_scores(gpu_quality)
        end
        gpu_decode_total = length(gpu_decoded_scores)
        @assert Array(gpu_decoded_scores) == cpu_quality_scores

        gpu_mean_ms, gpu_mean = elapsed_ms() do
            CUDA.@sync BioToolkit.mean_quality(gpu_quality)
        end
        @assert isapprox(gpu_mean, cpu_quality_mean; atol=1e-12)

        gpu_hist = Array(gpu_bins)
        @assert gpu_hist == cpu_hist
    end

    python_output = mktempdir() do dir
        python_script = joinpath(dir, "cuda_histogram_compare.py")
        write_python_benchmark_script(python_script, n, max_position, bin_size)
        read(`conda run -n general --no-capture-output python $python_script`, String)
    end

    println("CUDA histogram benchmark")
    println("  items: ", n)
    println("  bin_size: ", bin_size)
    println("  julia_cpu_ms=", round(cpu_ms, digits=4))
    println("  julia_cpu_bins=", length(cpu_hist))
    println("  julia_cpu_total=", sum(cpu_hist))
    if gpu_ms === missing
        println("  julia_gpu_ms=skipped")
        println("  julia_gpu_bins=skipped")
        println("  julia_gpu_hamming_ms=skipped")
        println("  julia_gpu_hamming=skipped")
        println("  julia_gpu_decode_ms=skipped")
        println("  julia_gpu_decode_total=skipped")
        println("  julia_gpu_mean_ms=skipped")
        println("  julia_gpu_mean=skipped")
    else
        println("  julia_gpu_ms=", round(gpu_ms, digits=4))
        println("  julia_gpu_bins=", length(gpu_hist))
        println("  julia_gpu_total=", sum(gpu_hist))
        println("  julia_gpu_hamming_ms=", round(gpu_hamming_ms, digits=4))
        println("  julia_gpu_hamming=", gpu_hamming)
        println("  julia_gpu_decode_ms=", round(gpu_decode_ms, digits=4))
        println("  julia_gpu_decode_total=", gpu_decode_total)
        println("  julia_gpu_mean_ms=", round(gpu_mean_ms, digits=4))
        println("  julia_gpu_mean=", round(gpu_mean, digits=4))
    end
    print(python_output)
end

main()