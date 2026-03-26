pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

const has_cuda_package = Base.find_package("CUDA") !== nothing
if has_cuda_package
    @eval using CUDA
end

function build_intervals(n::Integer, window_size::Integer)
    starts = Vector{Int32}(undef, n)
    stops = Vector{Int32}(undef, n)
    for index in 1:n
        start = Int32(mod(index * 37, 100_000) + 1)
        width = Int32(mod(index * 11, window_size * 2) + 1)
        starts[index] = start
        stops[index] = start + width
    end
    return starts, stops
end

function python_list_literal(values)
    return "[" * join(values, ", ") * "]"
end

function write_python_benchmark_script(path::AbstractString, starts, stops, window_size::Integer)
    open(path, "w") do io
        write(io, "from collections import defaultdict\n")
        write(io, "from time import perf_counter\n")
        write(io, "starts = $(python_list_literal(collect(starts)))\n")
        write(io, "stops = $(python_list_literal(collect(stops)))\n")
        write(io, "window_size = $(window_size)\n")
        write(io, "def time_call(func):\n")
        write(io, "    start = perf_counter()\n")
        write(io, "    result = func()\n")
        write(io, "    return (perf_counter() - start) * 1000.0, result\n")
        write(io, "def window_coverage(starts, stops, window_size):\n")
        write(io, "    coverage = defaultdict(int)\n")
        write(io, "    for start, stop in zip(starts, stops):\n")
        write(io, "        if stop <= start:\n")
        write(io, "            continue\n")
        write(io, "        first_window = start // window_size\n")
        write(io, "        last_window = (stop - 1) // window_size\n")
        write(io, "        for window_index in range(first_window, last_window + 1):\n")
        write(io, "            window_start = window_index * window_size\n")
        write(io, "            window_stop = window_start + window_size\n")
        write(io, "            overlap_start = max(start, window_start)\n")
        write(io, "            overlap_stop = min(stop, window_stop)\n")
        write(io, "            overlap = overlap_stop - overlap_start\n")
        write(io, "            if overlap > 0:\n")
        write(io, "                coverage[window_index] += overlap\n")
        write(io, "    return dict(coverage)\n")
        write(io, "py_ms, py_cov = time_call(lambda: window_coverage(starts, stops, window_size))\n")
        write(io, "print(f'python_window_ms={py_ms:.4f}')\n")
        write(io, "print(f'python_window_total={sum(py_cov.values())}')\n")
    end
end

function main()
    n = 25_000
    window_size = 100
    starts, stops = build_intervals(n, window_size)

    cpu_ms, cpu_cov = repeat_elapsed_ms(() -> BioToolkit.window_coverage(starts, stops, window_size), 1)

    gpu_ms = missing
    gpu_cov = nothing
    if has_cuda_package && CUDA.functional()
        try
            gpu_starts = CuArray(starts)
            gpu_stops = CuArray(stops)
            CUDA.@sync BioToolkit.window_coverage(gpu_starts, gpu_stops, window_size)
            gpu_ms, gpu_cov = elapsed_ms() do
                CUDA.@sync BioToolkit.window_coverage(gpu_starts, gpu_stops, window_size)
            end
            @assert gpu_cov == cpu_cov
        catch err
            @warn "Skipping CUDA window coverage benchmark" exception=(err, catch_backtrace())
            gpu_ms = missing
            gpu_cov = nothing
        end
    end

    python_output = mktempdir() do dir
        python_script = joinpath(dir, "window_coverage_compare.py")
        write_python_benchmark_script(python_script, starts, stops, window_size)
        read(`conda run -n general --no-capture-output python $python_script`, String)
    end

    python_total = parse(Int, match(r"python_window_total=(\d+)", python_output).captures[1])
    @assert python_total == sum(values(cpu_cov))

    println("Window coverage benchmark")
    println("  intervals=", n)
    println("  window_size=", window_size)
    println("  julia_window_ms=", round(cpu_ms, digits=4))
    println("  julia_window_total=", sum(values(cpu_cov)))
    if gpu_ms === missing
        println("  julia_window_cuda_ms=skipped")
        println("  julia_window_cuda_total=skipped")
    else
        println("  julia_window_cuda_ms=", round(gpu_ms, digits=4))
        println("  julia_window_cuda_total=", sum(values(gpu_cov)))
    end
    print(python_output)
end

main()
