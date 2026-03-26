pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function write_mock_embl(path::AbstractString)
    embl_data = """ID   XYZ123; SV 1; linear; DNA; ; ; 30 BP.
XX
AC   AC123456;
XX
DE   Test EMBL record.
XX
FT   source          1..30
FT                   /organism=\"Test organism\"
FT   gene            10..20
FT                   /gene=\"test_gene\"
XX
SQ   Sequence 30 BP; 10 A; 10 C; 10 G; 0 T; 0 OTHER;
     aaaaaaaaaa cccccccccc gggggggggg                                30
//
"""
    open(path, "w") do io
        write(io, embl_data)
    end
end

function write_python_script(path::AbstractString, embl_path::AbstractString, repetitions::Integer)
    open(path, "w") do io
        write(io, "from time import perf_counter\n")
        write(io, "from Bio import SeqIO\n")
        write(io, "embl_path = r'$(embl_path)'\n")
        write(io, "repetitions = $(repetitions)\n")
        write(io, "def repeat_elapsed_ms(func, repetitions):\n")
        write(io, "    start = perf_counter()\n")
        write(io, "    result = None\n")
        write(io, "    for _ in range(repetitions):\n")
        write(io, "        result = func()\n")
        write(io, "    return (perf_counter() - start) * 1000.0, result\n")
        write(io, "embl_ms, embl_record = repeat_elapsed_ms(lambda: SeqIO.read(embl_path, 'embl'), repetitions)\n")
        write(io, "print('Python flatfile benchmark')\n")
        write(io, "print(f'  repetitions={repetitions}')\n")
        write(io, "print(f'  python_embl_ms={embl_ms:.4f}')\n")
        write(io, "print(f'  python_embl_seq_len={len(embl_record.seq)}')\n")
        write(io, "print(f'  python_embl_feat_count={len(embl_record.features)}')\n")
        write(io, "print(f'  python_embl_id={embl_record.id}')\n")
    end
end

function main()
    repetitions = 200

    mktempdir() do dir
        embl_path = joinpath(dir, "mock.embl")
        python_script = joinpath(dir, "flatfile_compare.py")

        write_mock_embl(embl_path)
        write_python_script(python_script, embl_path, repetitions)

        embl_ms, embl_record = repeat_elapsed_ms(() -> first(read_embl(embl_path)), repetitions)

        python_output = read(`conda run -n general --no-capture-output python $python_script`, String)

        python_embl_ms = parse(Float64, match(r"python_embl_ms=([0-9.]+)", python_output).captures[1])
        python_embl_seq_len = parse(Int, match(r"python_embl_seq_len=(\d+)", python_output).captures[1])
        python_embl_id = match(r"python_embl_id=(.*)", python_output).captures[1]

        @assert length(embl_record.sequence) == python_embl_seq_len

        println("Julia flatfile benchmark")
        println("  repetitions=", repetitions)
        println("  julia_embl_ms=", round(embl_ms, digits=4))
        println("  julia_embl_seq_len=", length(embl_record.sequence))
        println("  julia_embl_feat_count=", length(embl_record.features))
        println("  julia_embl_id=", embl_record.identifier)
        print(python_output)
        println("  python_embl_ms=", round(python_embl_ms, digits=4))
    end
end

main()
