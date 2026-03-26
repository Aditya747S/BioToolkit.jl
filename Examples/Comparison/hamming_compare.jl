pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit

function elapsed_ms(f)
    GC.gc()
    start_time = time_ns()
    result = f()
    return (time_ns() - start_time) / 1e6, result
end

function repeat_elapsed_ms(f, repetitions::Integer)
    GC.gc()
    start_time = time_ns()
    result = nothing

    for _ in 1:repetitions
        result = f()
    end

    return (time_ns() - start_time) / 1e6, result
end

left = repeat("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 500)
right = replace(left, 'A' => 'T'; count=5)
python_script = joinpath(@__DIR__, "hamming_compare.py")

BioToolkit.hamming_distance(left, right)

repetitions = 1_000

distance_ms, distance = repeat_elapsed_ms(() -> BioToolkit.hamming_distance(left, right), repetitions)

python_output = read(`conda run -n general --no-capture-output python $python_script`, String)

println("Julia Hamming benchmark")
println("  repetitions=", repetitions)
println("  distance_ms=", round(distance_ms, digits=4))
println("  distance=", distance)
print(python_output)