pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

function write_mock_abif_v3(path::AbstractString)
    open(path, "w") do io
        write(io, "ABIF")
        write(io, hton(UInt16(100)))
        write(io, "tdir")
        write(io, hton(UInt32(1)))
        write(io, hton(UInt16(1023)))
        write(io, hton(UInt16(28)))
        write(io, hton(UInt32(2)))
        write(io, hton(UInt32(56)))
        write(io, hton(UInt32(128)))

        curr = position(io)
        for _ in 1:(128 - curr)
            write(io, UInt8(0))
        end

        # PBAS2: sequence data
        write(io, "PBAS")
        write(io, hton(UInt32(2)))
        write(io, hton(UInt16(2)))
        write(io, hton(UInt16(1)))
        write(io, hton(UInt32(4)))
        write(io, hton(UInt32(4)))
        write(io, "TAGC")

        # PCON2: quality values
        write(io, "PCON")
        write(io, hton(UInt32(2)))
        write(io, hton(UInt16(1)))
        write(io, hton(UInt16(1)))
        write(io, hton(UInt32(4)))
        write(io, hton(UInt32(4)))
        write(io, UInt8(10))
        write(io, UInt8(20))
        write(io, UInt8(30))
        write(io, UInt8(40))
    end
end

function main()
    repetitions = 500

    mktempdir() do dir
        abif_path = joinpath(dir, "mock.ab1")

        write_mock_abif_v3(abif_path)

        julia_ms, trace = repeat_elapsed_ms(() -> read_abif(abif_path), repetitions)

        println("Julia ABIF benchmark")
        println("  repetitions=", repetitions)
        println("  julia_abif_ms=", round(julia_ms, digits=4))
        println("  julia_abif_seq_len=", length(trace.sequence))
        println("  julia_abif_qual_len=", length(trace.qualities))
    end
end

main()
