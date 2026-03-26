using BioToolkit

function write_mock_abif_v3(filename)
    open(filename, "w") do io
        # Header (128 bytes)
        write(io, "ABIF")
        write(io, hton(UInt16(100))) # version
        
        # Directory entry in header (Starts at offset 6, 28 bytes)
        write(io, "tdir")
        write(io, hton(UInt32(1)))
        write(io, hton(UInt16(1023))) # type
        write(io, hton(UInt16(28)))   # size
        write(io, hton(UInt32(2)))    # 2 elements
        write(io, hton(UInt32(2 * 28))) # data size
        write(io, hton(UInt32(128)))  # data offset (points to directory)
        write(io, hton(UInt32(0)))    # spare/handle
        
        # Padding header to 128
        curr = position(io)
        for _ in 1:(128 - curr)
            write(io, UInt8(0))
        end
        
        # Directory at 128
        # Entry 1: PBAS1
        write(io, "PBAS")
        write(io, hton(UInt32(1)))
        write(io, hton(UInt16(2))) # Type 2 (Char)
        write(io, hton(UInt16(1))) # Element size
        write(io, hton(UInt32(4))) # 4 elements
        write(io, hton(UInt32(4))) # Total size
        write(io, "TAGC")           # DATA in offset field
        
        # Entry 2: FWO_1
        write(io, "FWO_")
        write(io, hton(UInt32(1)))
        write(io, hton(UInt16(2))) # Type 2 (Char)
        write(io, hton(UInt16(1)))
        write(io, hton(UInt32(4)))
        write(io, hton(UInt32(4)))
        write(io, "GATC")
    end
end

filename = "test_mock_v3.ab1"
write_mock_abif_v3(filename)

try
    println("Reading mock ABIF V3...")
    trace = BioToolkit.read_abif(filename)
    println("ABIF Sequence: ", trace.sequence)
    println("ABIF FWO: ", trace.annotations["fwo"])
    if trace.sequence == "TAGC" && trace.annotations["fwo"] == "GATC"
        println("ABIF Mock Test V3 Success!")
    else
        println("ABIF Mock Test V3 Failed: Got sequence '", trace.sequence, "'")
    end
finally
    rm(filename, force=true)
end
