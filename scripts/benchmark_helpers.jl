using BioToolkit: SeqRecordLite, FastqRecord, MultipleSequenceAlignment, read_abif, read_embl, read_genepop, write_genepop, linkage_disequilibrium, population_pca, feature_sequence, get_alignment_length, motif_consensus

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