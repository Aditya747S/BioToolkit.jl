pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..")))

using BioToolkit

include(joinpath(@__DIR__, "benchmark_helpers.jl"))

function write_python_benchmark_script(path::AbstractString, sequence::AbstractString)
    open(path, "w") do io
        write(io, "from time import perf_counter\n")
        write(io, "from Bio.Seq import Seq\n")
        write(io, "repetitions = 200\n")
        write(io, "sequence = Seq('$(sequence)')\n")
        write(io, "seq = sequence\n")
        write(io, "def time_call(func):\n")
        write(io, "    start = perf_counter()\n")
        write(io, "    result = None\n")
        write(io, "    for _ in range(repetitions):\n")
        write(io, "        result = func()\n")
        write(io, "    return (perf_counter() - start) * 1000.0, result\n")
        write(io, "gc_ms, gc_value = time_call(lambda: float(sequence.count('G') + sequence.count('C')) / len(sequence))\n")
        write(io, "rc_ms, rc_value = time_call(lambda: str(seq.reverse_complement()))\n")
        write(io, "tr_ms, tr_value = time_call(lambda: str(seq.translate(to_stop=True)))\n")
        write(io, "print(f'repetitions={repetitions}')\n")
        write(io, "print(f'biopython_gc_ms={gc_ms:.4f}')\n")
        write(io, "print(f'biopython_rc_ms={rc_ms:.4f}')\n")
        write(io, "print(f'biopython_translate_ms={tr_ms:.4f}')\n")
        write(io, "print(f'biopython_gc={gc_value:.6f}')\n")
        write(io, "print(f'biopython_rc_len={len(rc_value)}')\n")
        write(io, "print(f'biopython_translate_len={len(tr_value)}')\n")
    end
end

sequence = repeat("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 2_000)
repetitions = 200
sequence_bytes = Vector{UInt8}(codeunits(sequence))
translate_buffer = Vector{UInt8}(undef, length(sequence_bytes) ÷ 3)

mktempdir() do dir
    python_script = joinpath(dir, "sequence_benchmark.py")
    write_python_benchmark_script(python_script, sequence)

    BioToolkit.gc_content(sequence)
    BioToolkit.reverse_complement(sequence)
    BioToolkit.translate_dna!(translate_buffer, sequence_bytes; stop_at_stop=true)

    julia_gc_ms, julia_gc = repeat_elapsed_ms(() -> BioToolkit.gc_content(sequence), repetitions)

    julia_rc_ms, julia_rc = repeat_elapsed_ms(() -> BioToolkit.reverse_complement(sequence), repetitions)

    julia_translate_ms, julia_translate = repeat_elapsed_ms(() -> begin
        protein_length = BioToolkit.translate_dna!(translate_buffer, sequence_bytes; stop_at_stop=true)
        unsafe_string(pointer(translate_buffer), protein_length)
    end, repetitions)

    python_output = read(`conda run -n general --no-capture-output python $python_script`, String)

    println("Sequence benchmark")
    println("  repetitions: ", repetitions)
    println("  julia_gc_ms=", round(julia_gc_ms, digits=4))
    println("  julia_rc_ms=", round(julia_rc_ms, digits=4))
    println("  julia_translate_ms=", round(julia_translate_ms, digits=4))
    println("  julia_gc=", round(julia_gc, digits=6))
    println("  julia_rc_len=", length(julia_rc))
    println("  julia_translate_len=", length(julia_translate))
    print(python_output)
end