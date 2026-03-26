pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
using Mmap
using SHA

function write_sample_fasta(path::AbstractString)
    open(path, "w") do io
        write(io, ">seq1\n")
        write(io, repeat("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 1000))
        write(io, "\n")
    end
end

function elapsed_ms(f)
    GC.gc()
    start_time = time_ns()
    result = f()
    return (time_ns() - start_time) / 1e6, result
end

mktempdir() do dir
    fasta_path = joinpath(dir, "sample.fasta")
    python_script = joinpath(@__DIR__, "fasta_index_compare.py")
    write_sample_fasta(fasta_path)

    index_ms, index = elapsed_ms() do
        BioToolkit.fasta_index(fasta_path)
    end

    record = index["seq1"]
    open(fasta_path, "r") do io
        filebytes = Mmap.mmap(io)

        fetch_ms, subsequence = elapsed_ms() do
            BioToolkit.fetch_fasta_sequence(filebytes, record, 10_000, 11_000)
        end

        python_output = read(`conda run -n general --no-capture-output python $python_script $fasta_path`, String)
        python_subsequence_len = parse(Int, match(r"subsequence_len=(\d+)", python_output).captures[1])
        python_subsequence_sha1 = match(r"subsequence_sha1=([0-9a-f]+)", python_output).captures[1]
        julia_subsequence_sha1 = bytes2hex(sha1(subsequence))

        @assert length(subsequence) == python_subsequence_len
        @assert julia_subsequence_sha1 == python_subsequence_sha1

        println("Julia FASTA benchmark")
        println("  index_ms=", round(index_ms, digits=4))
        println("  fetch_ms=", round(fetch_ms, digits=4))
        println("  subsequence_len=", length(subsequence))
        print(python_output)
    end
end