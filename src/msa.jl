abstract type AbstractMultipleSequenceAlignment end

mutable struct MultipleSequenceAlignment <: AbstractMultipleSequenceAlignment
    records::Vector{SeqRecordLite}
    annotations::Dict{Symbol,Any}
    column_annotations::Dict{Symbol,Any}
end

function _copy_seqrecord(record::SeqRecordLite, sequence::AbstractString=record.sequence)
    return SeqRecordLite(
        sequence;
        identifier=record.identifier,
        name=record.name,
        description=record.description,
        annotations=record.annotations,
        letter_annotations=record.letter_annotations,
    )
end

function _coerce_msa_record(record::SeqRecordLite)
    return _copy_seqrecord(record)
end

function _coerce_msa_record(record::AbstractString, index::Integer)
    return SeqRecordLite(String(record); identifier="sequence_$(index)", name="sequence_$(index)", description="")
end

function _coerce_msa_records(records)
    aligned_records = SeqRecordLite[]
    for (index, record) in enumerate(records)
        if record isa SeqRecordLite
            push!(aligned_records, _coerce_msa_record(record))
        elseif record isa AbstractString
            push!(aligned_records, _coerce_msa_record(record, index))
        else
            throw(ArgumentError("MultipleSequenceAlignment expects SeqRecordLite or string records"))
        end
    end
    return aligned_records
end

function _column_annotation_length(value)
    value isa AbstractString && return lastindex(value)
    value isa AbstractVector && return length(value)
    return nothing
end

function _validate_column_annotations(column_annotations::AbstractDict, alignment_length::Integer)
    for (key, value) in column_annotations
        annotation_length = _column_annotation_length(value)
        annotation_length === nothing && continue
        annotation_length == alignment_length || throw(ArgumentError("column annotation $(key) must match the alignment length"))
    end
    return column_annotations
end

function _slice_column_annotation(value, cols)
    value isa AbstractString && return value[cols]
    value isa AbstractVector && return value[cols]
    return value
end

function _column_annotation_text(value)
    value isa AbstractString && return value
    value isa AbstractVector && return join(string.(value))
    return nothing
end

function _subset_column_annotations(column_annotations::AbstractDict, cols)
    subset = Dict{Symbol,Any}()
    for (key, value) in column_annotations
        subset[key] = _slice_column_annotation(value, cols)
    end
    return subset
end

function _concat_column_annotation(left_value, right_value)
    left_value isa AbstractString && right_value isa AbstractString && return string(left_value, right_value)
    left_value isa AbstractVector && right_value isa AbstractVector && return vcat(left_value, right_value)
    return string(left_value, right_value)
end

function _column_symbol_counts(alignment::MultipleSequenceAlignment, col::Integer; gap::Char='-')
    counts = Dict{Char,Int}()
    @inbounds for record in alignment.records
        symbol = record.sequence[col]
        symbol == gap && continue
        counts[symbol] = get(counts, symbol, 0) + 1
    end
    return counts
end

function _alignment_format_symbol(identity::Bool, strong::Bool, weak::Bool, style::Symbol)
    if identity
        return style === :emboss ? '|' : '*'
    elseif strong
        return ':'
    elseif weak
        return '.'
    else
        return ' '
    end
end

function _pairwise_residue_score(scoring::SubstitutionMatrix, left::Char, right::Char)
    left_byte = UInt8(left)
    right_byte = UInt8(right)
    left_index = get(scoring.lookup, left_byte, 0)
    right_index = get(scoring.lookup, right_byte, 0)
    return left_index == 0 || right_index == 0 ? scoring.default : scoring.scores[left_index, right_index]
end

function _column_consensus_symbol(
    alignment::MultipleSequenceAlignment,
    col::Integer;
    gap::Char='-',
    scoring::Union{Nothing,SubstitutionMatrix}=nothing,
    strong_threshold::Int=1,
    weak_threshold::Int=0,
    style::Symbol=:clustal,
)
    counts = _column_symbol_counts(alignment, col; gap=gap)
    isempty(counts) && return ' '
    length(counts) == 1 && return _alignment_format_symbol(true, false, false, style)
    scoring === nothing && return ' '

    residues = collect(keys(counts))
    all_scores = Int[]
    for i in 1:length(residues)
        for j in i + 1:length(residues)
            push!(all_scores, _pairwise_residue_score(scoring, residues[i], residues[j]))
        end
    end

    isempty(all_scores) && return _alignment_format_symbol(true, false, false, style)
    strong = all(score >= strong_threshold for score in all_scores)
    weak = any(score >= weak_threshold for score in all_scores)
    return _alignment_format_symbol(false, strong, weak, style)
end

function clustal_consensus(
    alignment::MultipleSequenceAlignment;
    gap::Char='-',
    scoring::Union{Nothing,SubstitutionMatrix}=nothing,
    strong_threshold::Int=1,
    weak_threshold::Int=0,
)
    return alignment_symbol_line(
        alignment;
        gap=gap,
        scoring=scoring,
        strong_threshold=strong_threshold,
        weak_threshold=weak_threshold,
        style=:clustal,
    )
end

function alignment_symbol_line(
    alignment::MultipleSequenceAlignment;
    gap::Char='-',
    scoring::Union{Nothing,SubstitutionMatrix}=nothing,
    strong_threshold::Int=1,
    weak_threshold::Int=0,
    style::Symbol=:clustal,
    threaded::Bool=true,
)
    if haskey(alignment.column_annotations, :clustal_consensus) && style === :clustal
        value = alignment.column_annotations[:clustal_consensus]
        value isa AbstractString && return value
        value isa AbstractVector && return join(value)
    end

    width = get_alignment_length(alignment)
    symbols = Vector{Char}(undef, width)
    if threaded && width > 1 && Threads.nthreads() > 1
        Threads.@threads for column in 1:width
            symbols[column] = _column_consensus_symbol(
                alignment,
                column;
                gap=gap,
                scoring=scoring,
                strong_threshold=strong_threshold,
                weak_threshold=weak_threshold,
                style=style,
            )
        end
    else
        for column in 1:width
            symbols[column] = _column_consensus_symbol(
                alignment,
                column;
                gap=gap,
                scoring=scoring,
                strong_threshold=strong_threshold,
                weak_threshold=weak_threshold,
                style=style,
            )
        end
    end
    return String(symbols)
end

function emboss_consensus(
    alignment::MultipleSequenceAlignment;
    gap::Char='-',
    scoring::Union{Nothing,SubstitutionMatrix}=nothing,
    strong_threshold::Int=1,
    weak_threshold::Int=0,
)
    return alignment_symbol_line(
        alignment;
        gap=gap,
        scoring=scoring,
        strong_threshold=strong_threshold,
        weak_threshold=weak_threshold,
        style=:emboss,
    )
end

function alignment_column_counts(alignment::MultipleSequenceAlignment, col::Integer; gap::Char='-')
    return _column_symbol_counts(alignment, col; gap=gap)
end

function alignment_column_frequencies(alignment::MultipleSequenceAlignment, col::Integer; gap::Char='-')
    counts = alignment_column_counts(alignment, col; gap=gap)
    total = sum(values(counts))
    total == 0 && return Dict{Char,Float64}()

    frequencies = Dict{Char,Float64}()
    for (symbol, count) in counts
        frequencies[symbol] = count / total
    end
    return frequencies
end

function alignment_column_symbol(
    alignment::MultipleSequenceAlignment,
    col::Integer;
    gap::Char='-',
    scoring::Union{Nothing,SubstitutionMatrix}=nothing,
    strong_threshold::Int=1,
    weak_threshold::Int=0,
    style::Symbol=:clustal,
)
    return _column_consensus_symbol(
        alignment,
        col;
        gap=gap,
        scoring=scoring,
        strong_threshold=strong_threshold,
        weak_threshold=weak_threshold,
        style=style,
    )
end

function _read_alignment_fasta(lines::AbstractVector{<:AbstractString})
    records = SeqRecordLite[]
    current_identifier = ""
    current_description = ""
    sequence = IOBuffer()

    function flush_record!()
        isempty(current_identifier) && return nothing
        push!(records, SeqRecordLite(String(take!(sequence)); identifier=current_identifier, description=current_description))
        return nothing
    end

    for raw_line in lines
        line = strip(raw_line)
        isempty(line) && continue
        if startswith(line, '>')
            flush_record!()
            current_description = strip(line[2:end])
            separator = findfirst(isspace, current_description)
            current_identifier = separator === nothing ? current_description : current_description[1:prevind(current_description, separator)]
        else
            write(sequence, line)
        end
    end

    flush_record!()
    return MultipleSequenceAlignment(records)
end

function _read_alignment_clustal(lines::AbstractVector{<:AbstractString})
    fragments = Dict{String,IOBuffer}()
    order = String[]
    consensus = IOBuffer()

    for raw_line in lines
        line = rstrip(raw_line)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "CLUSTAL") && continue
        startswith(stripped, "MUSCLE") && continue
        startswith(stripped, "#") && continue

        if all(ch -> ch in ('*', ':', '.', ' ', '|'), stripped)
            write(consensus, stripped)
            continue
        end

        separator = findfirst(isspace, stripped)
        separator === nothing && continue
        identifier = stripped[firstindex(stripped):prevind(stripped, separator)]
        fragment_start = findnext(!isspace, stripped, separator)
        fragment_start === nothing && continue
        fragment_stop = findnext(isspace, stripped, fragment_start)
        fragment = fragment_stop === nothing ? stripped[fragment_start:end] : stripped[fragment_start:prevind(stripped, fragment_stop)]

        buffer = get!(fragments, identifier) do
            push!(order, identifier)
            IOBuffer()
        end
        write(buffer, fragment)
    end

    records = SeqRecordLite[]
    for identifier in order
        sequence = String(take!(fragments[identifier]))
        push!(records, SeqRecordLite(sequence; identifier=identifier, description=identifier))
    end

    clustal_consensus = String(take!(consensus))
    column_annotations = isempty(clustal_consensus) ? Dict{Symbol,Any}() : Dict(:clustal_consensus => clustal_consensus)
    return MultipleSequenceAlignment(records; column_annotations=column_annotations)
end

function _read_alignment_stockholm(lines::AbstractVector{<:AbstractString})
    fragments = Dict{String,IOBuffer}()
    order = String[]
    column_annotation_buffers = Dict{Symbol,IOBuffer}()

    for raw_line in lines
        line = strip(raw_line)
        isempty(line) && continue
        line == "//" && break
        line == "# STOCKHOLM 1.0" && continue

        if startswith(line, "#=GC")
            fields = Base.split(line, limit=3)
            length(fields) == 3 || continue
            key = Symbol(fields[2])
            value = fields[3]
            buffer = get!(column_annotation_buffers, key) do
                IOBuffer()
            end
            write(buffer, value)
            continue
        end

        startswith(line, "#") && continue

        fields = Base.split(line)
        length(fields) < 2 && continue
        identifier = fields[1]
        fragment = fields[2]
        buffer = get!(fragments, identifier) do
            push!(order, identifier)
            IOBuffer()
        end
        write(buffer, fragment)
    end

    records = SeqRecordLite[]
    for identifier in order
        sequence = String(take!(fragments[identifier]))
        push!(records, SeqRecordLite(sequence; identifier=identifier, description=identifier))
    end

    column_annotations = Dict{Symbol,Any}(key => String(take!(buffer)) for (key, buffer) in column_annotation_buffers)
    return MultipleSequenceAlignment(records; column_annotations=column_annotations)
end

function read_alignment(source::AbstractString, format::AbstractString="auto")
    open(source, "r") do io
        return read_alignment(io, format)
    end
end

function read_alignment(io::IO, format::AbstractString="auto")
    lines = readlines(io)
    normalized = lowercase(strip(String(format)))
    if normalized == "auto"
        normalized = isempty(lines) ? "fasta" : startswith(uppercase(strip(first(lines))), "#NEXUS") ? "nexus" : startswith(strip(first(lines)), "!!") ? "msf" : startswith(strip(first(lines)), "!!") && occursin("Check:", first(lines)) ? "gcg" : startswith(strip(first(lines)), "##maf") ? "maf" : startswith(strip(first(lines)), ">P1;") || startswith(strip(first(lines)), ">F1;") ? "pir" : startswith(strip(first(lines)), "# STOCKHOLM") ? "stockholm" : startswith(strip(first(lines)), "CLUSTAL") ? "clustal" : occursin(r"^\s*\d+\s+\d+\s*$", strip(first(lines))) ? "phylip" : "fasta"
    end
    normalized in ("fasta", "fa") && return _read_alignment_fasta(lines)
    normalized == "clustal" && return _read_alignment_clustal(lines)
    normalized in ("stockholm", "sto") && return _read_alignment_stockholm(lines)
    normalized == "pir" && return _read_alignment_pir(lines)
    normalized == "nexus" && return _read_alignment_nexus(lines)
    normalized == "msf" && return _read_alignment_msf(lines)
    normalized == "phylip" && return _read_alignment_phylip(lines)
    normalized == "phylip-relaxed" && return _read_alignment_phylip(lines)
    normalized == "maf" && return _read_alignment_maf(lines)
    normalized == "gcg" && return _read_alignment_gcg(lines)
    throw(ArgumentError("unsupported alignment format: $(format)"))
end

function _validate_msa_records(records::Vector{SeqRecordLite})
    isempty(records) && return 0
    expected_length = length(records[1].sequence)
    @inbounds for record in records
        length(record.sequence) == expected_length || throw(ArgumentError("sequences must all have the same length"))
    end
    return expected_length
end

function MultipleSequenceAlignment(
    records::AbstractVector;
    annotations::AbstractDict=Dict{Symbol,Any}(),
    column_annotations::AbstractDict=Dict{Symbol,Any}(),
)
    aligned_records = _coerce_msa_records(records)
    alignment_length = _validate_msa_records(aligned_records)
    _validate_column_annotations(column_annotations, alignment_length)
    return MultipleSequenceAlignment(
        aligned_records,
        Dict{Symbol,Any}(annotations),
        Dict{Symbol,Any}(column_annotations),
    )
end

function Base.length(alignment::MultipleSequenceAlignment)
    return length(alignment.records)
end

function get_alignment_length(alignment::MultipleSequenceAlignment)
    isempty(alignment.records) && return 0
    return length(alignment.records[1].sequence)
end

Base.size(alignment::MultipleSequenceAlignment) = (length(alignment), get_alignment_length(alignment))

function Base.iterate(alignment::MultipleSequenceAlignment, state::Int=1)
    state > length(alignment.records) && return nothing
    return alignment.records[state], state + 1
end

function Base.getindex(alignment::MultipleSequenceAlignment, index::Integer)
    return alignment.records[index]
end

function Base.getindex(alignment::MultipleSequenceAlignment, ::Colon)
    return MultipleSequenceAlignment(
        alignment.records;
        annotations=alignment.annotations,
        column_annotations=alignment.column_annotations,
    )
end

function Base.copy(alignment::MultipleSequenceAlignment)
    return MultipleSequenceAlignment(
        alignment.records;
        annotations=alignment.annotations,
        column_annotations=alignment.column_annotations,
    )
end

function _slice_record(record::SeqRecordLite, cols)
    return _copy_seqrecord(record, record.sequence[cols])
end

function _normalize_rows(rows)
    rows isa Integer && return [rows]
    rows isa AbstractVector && return collect(rows)
    rows isa AbstractRange && return collect(rows)
    rows isa Colon && return nothing
    throw(ArgumentError("invalid row index"))
end

function _column_string(alignment::MultipleSequenceAlignment, row_indices, col::Integer)
    buffer = Vector{UInt8}(undef, length(row_indices))
    for (i, row_index) in enumerate(row_indices)
        buffer[i] = UInt8(alignment.records[row_index].sequence[col])
    end
    return String(buffer)
end

function Base.getindex(alignment::MultipleSequenceAlignment, rows::Colon, col::Integer)
    return _column_string(alignment, eachindex(alignment.records), col)
end

function Base.getindex(alignment::MultipleSequenceAlignment, row::Integer, col::Integer)
    return alignment.records[row].sequence[col]
end

function Base.getindex(alignment::MultipleSequenceAlignment, row::Integer, ::Colon)
    return alignment.records[row]
end

function Base.getindex(alignment::MultipleSequenceAlignment, row::Integer, cols)
    return _slice_record(alignment.records[row], cols)
end

function Base.getindex(alignment::MultipleSequenceAlignment, rows::Colon, cols)
    if cols isa Integer
        return _column_string(alignment, eachindex(alignment.records), cols)
    end

    selected = SeqRecordLite[]
    for record in alignment.records
        push!(selected, _slice_record(record, cols))
    end

    return MultipleSequenceAlignment(
        selected;
        annotations=alignment.annotations,
        column_annotations=_subset_column_annotations(alignment.column_annotations, cols),
    )
end

function Base.getindex(alignment::MultipleSequenceAlignment, rows, cols::Colon)
    row_indices = _normalize_rows(rows)
    row_indices === nothing && return MultipleSequenceAlignment(alignment.records; annotations=alignment.annotations, column_annotations=alignment.column_annotations)
    return MultipleSequenceAlignment(alignment.records[row_indices]; annotations=alignment.annotations, column_annotations=alignment.column_annotations)
end

function Base.getindex(alignment::MultipleSequenceAlignment, rows::AbstractVector)
    return MultipleSequenceAlignment(
        alignment.records[rows];
        annotations=alignment.annotations,
        column_annotations=_subset_column_annotations(alignment.column_annotations, :),
    )
end

function Base.getindex(alignment::MultipleSequenceAlignment, rows::AbstractRange)
    return MultipleSequenceAlignment(
        alignment.records[collect(rows)];
        annotations=alignment.annotations,
        column_annotations=_subset_column_annotations(alignment.column_annotations, :),
    )
end

function Base.getindex(alignment::MultipleSequenceAlignment, rows, cols)
    if rows isa Integer && cols isa Integer
        return alignment.records[rows].sequence[cols]
    elseif rows isa Integer
        return _slice_record(alignment.records[rows], cols)
    else
        row_indices = _normalize_rows(rows)
        row_indices === nothing && return MultipleSequenceAlignment(alignment.records; annotations=alignment.annotations, column_annotations=alignment.column_annotations)
        if cols isa Integer
            return _column_string(alignment, row_indices, cols)
        end
        selected = SeqRecordLite[]
        for row_index in row_indices
            push!(selected, _slice_record(alignment.records[row_index], cols))
        end
        return MultipleSequenceAlignment(
            selected;
            annotations=alignment.annotations,
            column_annotations=_subset_column_annotations(alignment.column_annotations, cols),
        )
    end
end

function Base.push!(alignment::MultipleSequenceAlignment, record::SeqRecordLite)
    if isempty(alignment.records)
        push!(alignment.records, _copy_seqrecord(record))
    else
        expected_length = get_alignment_length(alignment)
        length(record.sequence) == expected_length || throw(ArgumentError("sequence length must match the alignment"))
        push!(alignment.records, _copy_seqrecord(record))
    end
    return alignment
end

function Base.append!(alignment::MultipleSequenceAlignment, record::SeqRecordLite)
    return push!(alignment, record)
end

function Base.deleteat!(alignment::MultipleSequenceAlignment, index)
    deleteat!(alignment.records, index)
    return alignment
end

function Base.sort!(alignment::MultipleSequenceAlignment; key=nothing, reverse::Bool=false)
    if key === nothing
        sort!(alignment.records; by = record -> record.identifier, rev=reverse)
    else
        sort!(alignment.records; by = key, rev=reverse)
    end
    return alignment
end

function Base.:+(left::MultipleSequenceAlignment, right::MultipleSequenceAlignment)
    length(left) == length(right) || throw(ArgumentError("alignments must have the same number of rows"))
    merged = SeqRecordLite[]
    for (left_record, right_record) in zip(left.records, right.records)
        merged_sequence = string(left_record.sequence, right_record.sequence)
        push!(merged, _copy_seqrecord(left_record, merged_sequence))
    end
    annotations = Dict{Symbol,Any}()
    for (key, value) in left.annotations
        haskey(right.annotations, key) && right.annotations[key] == value && (annotations[key] = value)
    end
    column_annotations = Dict{Symbol,Any}()
    for (key, value) in left.column_annotations
        if haskey(right.column_annotations, key)
            column_annotations[key] = _concat_column_annotation(value, right.column_annotations[key])
        end
    end
    return MultipleSequenceAlignment(merged; annotations=annotations, column_annotations=column_annotations)
end

function _alignment_block_lines(alignment::MultipleSequenceAlignment, start_col::Integer, stop_col::Integer)
    lines = String[]
    for (key, value) in alignment.column_annotations
        block = _column_annotation_text(value)
        block === nothing && continue
        push!(lines, string(key, " ", block[start_col:stop_col]))
    end

    for record in alignment.records
        push!(lines, string(record.identifier, " ", record.sequence[start_col:stop_col]))
    end

    return lines
end

function _alignment_text(alignment::MultipleSequenceAlignment; width::Integer=60)
    rows = length(alignment)
    cols = get_alignment_length(alignment)
    io = IOBuffer()
    print(io, "Alignment with ", rows, " rows and ", cols, " columns")
    rows == 0 && return String(take!(io))

    width = max(1, width)
    for start_col in 1:width:cols
        stop_col = min(start_col + width - 1, cols)
        print(io, "\n")
        for line in _alignment_block_lines(alignment, start_col, stop_col)
            print(io, line, "\n")
        end
    end

    return String(take!(io))
end

function _wrap_sequence(sequence::AbstractString, width::Integer)
    width <= 0 && return [sequence]
    pieces = String[]
    start_index = firstindex(sequence)
    while start_index <= lastindex(sequence)
        stop_index = min(start_index + width - 1, lastindex(sequence))
        push!(pieces, sequence[start_index:stop_index])
        start_index = stop_index + 1
    end
    return pieces
end

function _write_fasta(io::IO, alignment::MultipleSequenceAlignment; width::Integer=60)
    isempty(alignment.records) && throw(ArgumentError("cannot format an empty multiple sequence alignment"))
    for record in alignment.records
        # Bug fix: only append description when it is non-empty AND differs from
        # the identifier; otherwise FASTA headers come out as ">seq1 seq1".
        header = (isempty(record.description) || record.description == record.identifier) ?
                 record.identifier : string(record.identifier, " ", record.description)
        print(io, ">", header, "\n")
        sequence = record.sequence
        if width <= 0
            print(io, sequence, "\n")
            continue
        end

        start_index = firstindex(sequence)
        while start_index <= lastindex(sequence)
            stop_index = min(start_index + width - 1, lastindex(sequence))
            print(io, SubString(sequence, start_index, stop_index), "\n")
            start_index = stop_index + 1
        end
    end
end

function _write_stockholm(io::IO, alignment::MultipleSequenceAlignment; width::Integer=60)
    isempty(alignment.records) && throw(ArgumentError("cannot format an empty multiple sequence alignment"))
    alignment_width = get_alignment_length(alignment)
    max_identifier = maximum(length(record.identifier) for record in alignment.records)
    block_width = max(1, width)
    padded_identifiers = [rpad(record.identifier, max_identifier + 2) for record in alignment.records]
    annotation_blocks = Pair{Symbol,String}[]
    for (key, value) in alignment.column_annotations
        block = _column_annotation_text(value)
        block === nothing && continue
        push!(annotation_blocks, key => block)
    end

    print(io, "# STOCKHOLM 1.0\n\n")
    for start_col in 1:block_width:alignment_width
        stop_col = min(start_col + block_width - 1, alignment_width)
        for (row_index, record) in enumerate(alignment.records)
            print(io, padded_identifiers[row_index])
            print(io, SubString(record.sequence, start_col, stop_col), "\n")
        end
        for (key, block) in annotation_blocks
            print(io, "#=GC ", key, " ", SubString(block, start_col, stop_col), "\n")
        end
        print(io, "\n")
    end
    print(io, "//")
end

function _msf_chunk(sequence::AbstractString, start_index::Int, width::Int)
    stop_index = min(start_index + width - 1, lastindex(sequence))
    chunk = sequence[start_index:stop_index]
    groups = String[]
    for group_start in 1:10:length(chunk)
        push!(groups, chunk[group_start:min(group_start + 9, lastindex(chunk))])
    end
    return join(groups, " ")
end

function _write_msf(io::IO, alignment::MultipleSequenceAlignment)
    isempty(alignment.records) && throw(ArgumentError("cannot format an empty multiple sequence alignment"))
    alignment_width = get_alignment_length(alignment)
    sequence_type = all(record -> all(base -> base in ('A', 'C', 'G', 'T', '-', 'N', 'a', 'c', 'g', 't', 'n'), record.sequence), alignment.records) ? "N" : "P"

    println(io, "!!$(sequence_type)A_MULTIPLE_ALIGNMENT 1.0")
    println(io)
    println(io, " MSF: ", alignment_width, " Type: ", sequence_type, " Check: 0 ..")
    println(io)
    for record in alignment.records
        println(io, " Name: ", record.identifier, " Len: ", alignment_width, " Check: 0 Weight: 1.00")
    end
    println(io)
    println(io, "//")

    block_width = 50
    for start_index in 1:block_width:alignment_width
        println(io)
        for record in alignment.records
            print(io, rpad(record.identifier, 16), _msf_chunk(record.sequence, start_index, block_width), "\n")
        end
    end
end

function _pir_identifier(identifier::AbstractString)
    return replace(String(identifier), r"\s+" => "_")
end

function _write_pir(io::IO, alignment::MultipleSequenceAlignment; width::Integer=60)
    isempty(alignment.records) && throw(ArgumentError("cannot format an empty multiple sequence alignment"))
    block_width = max(1, width)

    for record in alignment.records
        println(io, ">P1;", _pir_identifier(record.identifier))
        description = isempty(record.description) ? record.identifier : record.description
        println(io, description)
        for start_index in 1:block_width:lastindex(record.sequence)
            stop_index = min(start_index + block_width - 1, lastindex(record.sequence))
            println(io, SubString(record.sequence, start_index, stop_index))
        end
        println(io, "*")
    end
end

function _read_alignment_pir(lines::AbstractVector{<:AbstractString})
    records = Dict{String,IOBuffer}()
    order = String[]
    current_identifier = ""
    reading_description = false

    for raw_line in lines
        line = strip(raw_line)
        isempty(line) && continue

        if startswith(line, ">")
            header = line[2:end]
            separator = findfirst(';', header)
            separator === nothing && throw(ArgumentError("invalid PIR header"))
            current_identifier = strip(header[separator + 1:end])
            buffer = get!(records, current_identifier) do
                push!(order, current_identifier)
                IOBuffer()
            end
            reading_description = true
            continue
        end

        current_identifier == "" && continue
        reading_description && (reading_description = false; continue)

        stripped_sequence = replace(line, r"\s+" => "")
        stripped_sequence == "*" && continue
        endswith(stripped_sequence, "*") && (stripped_sequence = stripped_sequence[1:end-1])
        stripped_sequence == "" && continue
        write(records[current_identifier], stripped_sequence)
    end

    parsed_records = SeqRecordLite[]
    alignment_length = nothing
    for identifier in order
        sequence = String(take!(records[identifier]))
        alignment_length === nothing && (alignment_length = length(sequence))
        length(sequence) == alignment_length || throw(ArgumentError("PIR sequence length mismatch"))
        push!(parsed_records, SeqRecordLite(sequence; identifier=identifier, description=identifier))
    end

    return MultipleSequenceAlignment(parsed_records)
end

function _read_alignment_msf(lines::AbstractVector{<:AbstractString})
    records = Dict{String,IOBuffer}()
    order = String[]
    in_matrix = false
    alignment_length = nothing

    for raw_line in lines
        line = rstrip(raw_line)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "//") && (in_matrix = true; continue)

        if !in_matrix
            if occursin(r"\bMSF:\s*(\d+)", uppercase(stripped))
                if (m = match(r"MSF:\s*(\d+)", uppercase(stripped))) !== nothing
                    alignment_length = parse(Int, m.captures[1])
                end
            end
            continue
        end

        startswith(stripped, "Name:") && continue
        startswith(stripped, "!!") && continue

        fields = Base.split(stripped)
        length(fields) < 2 && continue
        identifier = fields[1]
        fragment = join(fields[2:end], "")
        buffer = get!(records, identifier) do
            push!(order, identifier)
            IOBuffer()
        end
        write(buffer, fragment)
    end

    parsed_records = SeqRecordLite[]
    for identifier in order
        sequence = String(take!(records[identifier]))
        alignment_length !== nothing && length(sequence) == alignment_length || throw(ArgumentError("MSF sequence length mismatch"))
        push!(parsed_records, SeqRecordLite(sequence; identifier=identifier, description=identifier))
    end

    return MultipleSequenceAlignment(parsed_records)
end

function _normalize_alignment_format(format::AbstractString)
    normalized = lowercase(strip(String(format)))
    normalized in ("fasta", "fa") && return "fasta"
    normalized == "clustal" && return "clustal"
    normalized in ("stockholm", "sto") && return "stockholm"
    normalized in ("pir", "nbrf") && return "pir"
    normalized == "nexus" && return "nexus"
    normalized == "msf" && return "msf"
    normalized in ("phylip", "phylip-sequential", "phylip_sequential") && return "phylip"
    normalized in ("phylip-relaxed", "phylip_relaxed") && return "phylip-relaxed"
    throw(ArgumentError("unsupported alignment format: $(format)"))
end

function _nexus_identifier(identifier::AbstractString)
    return replace(String(identifier), r"\s+" => "_")
end

function _write_nexus(io::IO, alignment::MultipleSequenceAlignment)
    isempty(alignment.records) && throw(ArgumentError("cannot format an empty multiple sequence alignment"))
    alignment_width = get_alignment_length(alignment)

    println(io, "#NEXUS")
    println(io, "Begin data;")
    println(io, "Dimensions ntax=", length(alignment.records), " nchar=", alignment_width, ";")
    println(io, "Format datatype=standard missing=? gap=- interleave=no;")
    println(io, "Matrix")
    for record in alignment.records
        println(io, _nexus_identifier(record.identifier), " ", record.sequence)
    end
    println(io, ";")
    println(io, "End;")
end

function _read_alignment_nexus(lines::AbstractVector{<:AbstractString})
    records = Dict{String,IOBuffer}()
    order = String[]
    in_data = false
    in_matrix = false
    sequence_count = nothing
    alignment_length = nothing
    last_identifier = ""

    for raw_line in lines
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, "[") && continue

        upper_line = uppercase(line)
        if startswith(upper_line, "#NEXUS")
            continue
        elseif startswith(upper_line, "BEGIN DATA")
            in_data = true
            continue
        elseif in_data && startswith(upper_line, "DIMENSIONS")
            if (m = match(r"NTAX\s*=\s*(\d+)", upper_line)) !== nothing
                sequence_count = parse(Int, m.captures[1])
            end
            if (m = match(r"NCHAR\s*=\s*(\d+)", upper_line)) !== nothing
                alignment_length = parse(Int, m.captures[1])
            end
            continue
        elseif in_data && startswith(upper_line, "MATRIX")
            in_matrix = true
            continue
        elseif in_matrix && line == ";"
            break
        elseif in_data && startswith(upper_line, "END")
            break
        end

        if in_matrix
            if startswith(raw_line, ' ') || startswith(raw_line, '\t')
                last_identifier == "" && continue
                fragment = replace(line, r"\s+" => "")
                write(records[last_identifier], fragment)
                continue
            end

            fields = Base.split(line)
            length(fields) >= 2 || continue
            identifier = fields[1]
            fragment = join(fields[2:end], "")
            buffer = get!(records, identifier) do
                push!(order, identifier)
                IOBuffer()
            end
            write(buffer, fragment)
            last_identifier = identifier
        end
    end

    parsed_records = SeqRecordLite[]
    for identifier in order
        sequence = String(take!(records[identifier]))
        if alignment_length !== nothing && length(sequence) != alignment_length
            throw(ArgumentError("NEXUS sequence length mismatch"))
        end
        push!(parsed_records, SeqRecordLite(sequence; identifier=identifier, description=identifier))
    end

    if sequence_count !== nothing && length(parsed_records) != sequence_count
        throw(ArgumentError("NEXUS record count mismatch"))
    end
    return MultipleSequenceAlignment(parsed_records)
end

function _phylip_identifier(identifier::AbstractString; relaxed::Bool=false)
    relaxed && return String(identifier)
    stop_index = min(lastindex(identifier), firstindex(identifier) + 9)
    return String(identifier[firstindex(identifier):stop_index])
end

function _write_phylip(io::IO, alignment::MultipleSequenceAlignment; relaxed::Bool=false)
    isempty(alignment.records) && throw(ArgumentError("cannot format an empty multiple sequence alignment"))
    alignment_width = get_alignment_length(alignment)
    name_width = relaxed ? max(10, maximum(length(record.identifier) for record in alignment.records) + 2) : 10

    println(io, length(alignment.records), " ", alignment_width)
    for record in alignment.records
        identifier = _phylip_identifier(record.identifier; relaxed=relaxed)
        print(io, rpad(identifier, name_width), record.sequence, "\n")
    end
end

function _read_alignment_phylip(lines::AbstractVector{<:AbstractString})
    stripped_lines = [strip(line) for line in lines if !isempty(strip(line))]
    isempty(stripped_lines) && return MultipleSequenceAlignment(SeqRecordLite[])

    header = Base.split(stripped_lines[1])
    length(header) >= 2 || throw(ArgumentError("invalid PHYLIP header"))
    sequence_count = parse(Int, header[1])
    alignment_length = parse(Int, header[2])

    records = Dict{String,IOBuffer}()
    order = String[]
    continuation_index = 1

    for raw_line in @view lines[2:end]
        line = rstrip(raw_line)
        stripped = strip(line)
        isempty(stripped) && (continuation_index = 1; continue)

        if !isempty(line) && isspace(first(line))
            continuation_index > length(order) && continue
            fragment = replace(stripped, r"\s+" => "")
            identifier = order[continuation_index]
            write(records[identifier], fragment)
            continuation_index += 1
            continuation_index > length(order) && (continuation_index = 1)
            continue
        end

        fields = Base.split(stripped)
        length(fields) >= 2 || continue
        identifier = fields[1]
        fragment = join(fields[2:end], "")
        buffer = get!(records, identifier) do
            push!(order, identifier)
            IOBuffer()
        end
        write(buffer, fragment)
    end

    parsed_records = SeqRecordLite[]
    for identifier in order
        sequence = String(take!(records[identifier]))
        length(sequence) == alignment_length || throw(ArgumentError("PHYLIP sequence length mismatch"))
        push!(parsed_records, SeqRecordLite(sequence; identifier=identifier, description=identifier))
    end

    length(parsed_records) == sequence_count || throw(ArgumentError("PHYLIP record count mismatch"))
    return MultipleSequenceAlignment(parsed_records)
end

function _write_clustal(
    io::IO,
    alignment::MultipleSequenceAlignment;
    width::Integer=60,
    scoring::Union{Nothing,SubstitutionMatrix}=nothing,
    strong_threshold::Int=1,
    weak_threshold::Int=0,
    symbol_style::Symbol=:clustal,
)
    isempty(alignment.records) && throw(ArgumentError("cannot format an empty multiple sequence alignment"))
    alignment_width = get_alignment_length(alignment)
    consensus = alignment_symbol_line(alignment; scoring=scoring, strong_threshold=strong_threshold, weak_threshold=weak_threshold, style=symbol_style)
    max_identifier = maximum(length(record.identifier) for record in alignment.records)
    block_width = max(1, width)
    padded_identifiers = [rpad(record.identifier, max_identifier + 2) for record in alignment.records]

    print(io, "CLUSTAL multiple sequence alignment\n\n")
    for start_col in 1:block_width:alignment_width
        stop_col = min(start_col + block_width - 1, alignment_width)
        for (row_index, record) in enumerate(alignment.records)
            print(io, padded_identifiers[row_index])
            print(io, SubString(record.sequence, start_col, stop_col), "\n")
        end
        print(io, rpad("", max_identifier + 2))
        print(io, SubString(consensus, start_col, stop_col), "\n\n")
    end

end

function _format_fasta(alignment::MultipleSequenceAlignment; width::Integer=60)
    return sprint(io -> _write_fasta(io, alignment; width=width))
end

function _format_clustal(
    alignment::MultipleSequenceAlignment;
    width::Integer=60,
    scoring::Union{Nothing,SubstitutionMatrix}=nothing,
    strong_threshold::Int=1,
    weak_threshold::Int=0,
    symbol_style::Symbol=:clustal,
)
    return sprint(io -> _write_clustal(io, alignment; width=width, scoring=scoring, strong_threshold=strong_threshold, weak_threshold=weak_threshold, symbol_style=symbol_style))
end

function _format_stockholm(alignment::MultipleSequenceAlignment; width::Integer=60)
    return sprint(io -> _write_stockholm(io, alignment; width=width))
end

function _format_pir(alignment::MultipleSequenceAlignment; width::Integer=60)
    return sprint(io -> _write_pir(io, alignment; width=width))
end

function _format_msf(alignment::MultipleSequenceAlignment)
    return sprint(io -> _write_msf(io, alignment))
end

function _format_nexus(alignment::MultipleSequenceAlignment)
    return sprint(io -> _write_nexus(io, alignment))
end

function _format_phylip(alignment::MultipleSequenceAlignment; relaxed::Bool=false)
    return sprint(io -> _write_phylip(io, alignment; relaxed=relaxed))
end

function format_alignment(
    alignment::MultipleSequenceAlignment,
    format::AbstractString="clustal";
    width::Integer=60,
    scoring::Union{Nothing,SubstitutionMatrix}=nothing,
    strong_threshold::Int=1,
    weak_threshold::Int=0,
    symbol_style::Symbol=:clustal,
)
    normalized = _normalize_alignment_format(format)
    if normalized in ("fasta", "fa")
        return _format_fasta(alignment; width=width)
    elseif normalized == "clustal"
        return _format_clustal(alignment; width=width, scoring=scoring, strong_threshold=strong_threshold, weak_threshold=weak_threshold, symbol_style=symbol_style)
    elseif normalized in ("stockholm", "sto")
        return _format_stockholm(alignment; width=width)
    elseif normalized == "pir"
        return _format_pir(alignment; width=width)
    elseif normalized == "nexus"
        return _format_nexus(alignment)
    elseif normalized == "msf"
        return _format_msf(alignment)
    elseif normalized == "phylip"
        return _format_phylip(alignment; relaxed=false)
    elseif normalized == "phylip-relaxed"
        return _format_phylip(alignment; relaxed=true)
    end
    throw(ArgumentError("unsupported alignment format: $(format)"))
end

function write_alignment(
    io::IO,
    alignment::MultipleSequenceAlignment,
    format::AbstractString="clustal";
    width::Integer=60,
    scoring::Union{Nothing,SubstitutionMatrix}=nothing,
    strong_threshold::Int=1,
    weak_threshold::Int=0,
    symbol_style::Symbol=:clustal,
)
    normalized = _normalize_alignment_format(format)
    if normalized in ("fasta", "fa")
        _write_fasta(io, alignment; width=width)
    elseif normalized == "clustal"
        _write_clustal(io, alignment; width=width, scoring=scoring, strong_threshold=strong_threshold, weak_threshold=weak_threshold, symbol_style=symbol_style)
    elseif normalized in ("stockholm", "sto")
        _write_stockholm(io, alignment; width=width)
    elseif normalized == "pir"
        _write_pir(io, alignment; width=width)
    elseif normalized == "nexus"
        _write_nexus(io, alignment)
    elseif normalized == "msf"
        _write_msf(io, alignment)
    elseif normalized == "phylip"
        _write_phylip(io, alignment; relaxed=false)
    elseif normalized == "phylip-relaxed"
        _write_phylip(io, alignment; relaxed=true)
    else
        throw(ArgumentError("unsupported alignment format: $(format)"))
    end
    return alignment
end

function write_alignment(
    path::AbstractString,
    alignment::MultipleSequenceAlignment,
    format::AbstractString="clustal";
    width::Integer=60,
    scoring::Union{Nothing,SubstitutionMatrix}=nothing,
    strong_threshold::Int=1,
    weak_threshold::Int=0,
    symbol_style::Symbol=:clustal,
)
    open(path, "w") do io
        write_alignment(io, alignment, format; width=width, scoring=scoring, strong_threshold=strong_threshold, weak_threshold=weak_threshold, symbol_style=symbol_style)
    end
    return alignment
end

function _normalize_external_alignment_format(format::AbstractString)
    normalized = _normalize_alignment_format(format)
    normalized in ("fasta", "fa") && return "fasta"
    normalized == "clustal" && return "clustal"
    normalized in ("stockholm", "sto") && return "stockholm"
    throw(ArgumentError("unsupported external alignment format: $(format)"))
end

function _alignment_output_extension(format::AbstractString)
    normalized = _normalize_external_alignment_format(format)
    normalized == "fasta" && return "fasta"
    normalized == "clustal" && return "aln"
    return "sto"
end

function _external_msa_records(records)
    records isa MultipleSequenceAlignment && return records.records
    records isa AbstractVector || throw(ArgumentError("external MSA wrappers expect a vector of sequences or a MultipleSequenceAlignment"))

    coerced = SeqRecordLite[]
    for (index, record) in enumerate(records)
        if record isa SeqRecordLite
            push!(coerced, _copy_seqrecord(record))
        elseif record isa AbstractString
            push!(coerced, SeqRecordLite(String(record); identifier="sequence_$(index)", name="sequence_$(index)", description=""))
        else
            throw(ArgumentError("external MSA wrappers expect SeqRecordLite or string records"))
        end
    end
    return coerced
end

function _write_fasta_records(io::IO, records::AbstractVector{SeqRecordLite})
    for record in records
        header = (isempty(record.description) || record.description == record.identifier) ? record.identifier : string(record.identifier, " ", record.description)
        print(io, ">", header, "\n", record.sequence, "\n")
    end
    return nothing
end

function _external_msa_search_paths(candidate::AbstractString)
    paths = String[]
    conda_prefix = get(ENV, "CONDA_PREFIX", "")
    if !isempty(conda_prefix)
        push!(paths, joinpath(conda_prefix, "bin", candidate))
    end

    home = homedir()
    for root in (
        joinpath(home, "miniconda3", "envs"),
        joinpath(home, ".conda", "envs"),
        joinpath(home, "anaconda3", "envs"),
    )
        isdir(root) || continue
        for env_name in readdir(root)
            push!(paths, joinpath(root, env_name, "bin", candidate))
        end
    end

    push!(paths, candidate)
    return paths
end

function _resolve_external_msa_executable(candidates::AbstractVector{<:AbstractString}, executable::Union{Nothing,AbstractString})
    if executable !== nothing
        executable_string = String(executable)
        if occursin('/', executable_string) || occursin('\\', executable_string)
            isfile(executable_string) && return executable_string
        else
            resolved = Sys.which(executable_string)
            resolved !== nothing && return resolved
        end
        throw(ArgumentError("unable to locate external MSA executable: $(executable)"))
    end

    for candidate in candidates
        for path in _external_msa_search_paths(candidate)
            isfile(path) && return path
        end
        resolved = Sys.which(candidate)
        resolved !== nothing && return resolved
    end

    throw(ArgumentError("no external MSA executable found on PATH"))
end

function _run_external_msa(
    tool::Symbol,
    records;
    executable::Union{Nothing,AbstractString}=nothing,
    output_format::AbstractString="clustal",
    extra_args::AbstractVector{<:AbstractString}=String[],
)
    input_records = _external_msa_records(records)
    normalized_output_format = _normalize_external_alignment_format(output_format)
    output_extension = _alignment_output_extension(normalized_output_format)

    mktempdir() do dir
        input_path = joinpath(dir, "input.fasta")
        output_path = joinpath(dir, "output.$(output_extension)")
        open(input_path, "w") do io
            _write_fasta_records(io, input_records)
        end

        executable_path = _resolve_external_msa_executable(
            tool === :clustal ? ["clustalo", "clustalw2", "clustalw"] : ["muscle"],
            executable,
        )

        if tool === :clustal
            command_variants = [
                Cmd(vcat([executable_path, "--infile", input_path, "--outfile", output_path, "--outfmt", normalized_output_format, "--force"], String.(extra_args))),
                Cmd(vcat([executable_path, "-infile=$(input_path)", "-outfile=$(output_path)", "-output=$(uppercase(normalized_output_format))"], String.(extra_args))),
            ]
            last_error = nothing
            for command in command_variants
                try
                    run(command)
                    return read_alignment(output_path, normalized_output_format)
                catch err
                    last_error = err
                end
            end
            throw(ArgumentError("clustal execution failed: $(last_error === nothing ? "unknown error" : sprint(showerror, last_error))"))
        elseif tool === :muscle
            normalized_output_format == "fasta" || throw(ArgumentError("muscle currently supports FASTA output only"))
            command_variants = [
                Cmd(vcat([executable_path, "-align", input_path, "-output", output_path], String.(extra_args))),
                Cmd(vcat([executable_path, "-in", input_path, "-out", output_path], String.(extra_args))),
            ]
            last_error = nothing
            for command in command_variants
                try
                    run(command)
                    return read_alignment(output_path, normalized_output_format)
                catch err
                    last_error = err
                end
            end
            throw(ArgumentError("muscle execution failed: $(last_error === nothing ? "unknown error" : sprint(showerror, last_error))"))
        else
            throw(ArgumentError("unsupported external MSA tool: $(tool)"))
        end
    end
end

clustal_available(; executable::Union{Nothing,AbstractString}=nothing) = try
    _resolve_external_msa_executable(["clustalo", "clustalw2", "clustalw"], executable)
    true
catch
    false
end

muscle_available(; executable::Union{Nothing,AbstractString}=nothing) = try
    _resolve_external_msa_executable(["muscle"], executable)
    true
catch
    false
end

function clustal_msa(records; executable::Union{Nothing,AbstractString}=nothing, output_format::AbstractString="clustal", extra_args::AbstractVector{<:AbstractString}=String[])
    return _run_external_msa(:clustal, records; executable=executable, output_format=output_format, extra_args=extra_args)
end

function muscle_msa(records; executable::Union{Nothing,AbstractString}=nothing, output_format::AbstractString="clustal", extra_args::AbstractVector{<:AbstractString}=String[])
    return _run_external_msa(:muscle, records; executable=executable, output_format=output_format, extra_args=extra_args)
end

function Base.show(io::IO, alignment::MultipleSequenceAlignment)
    print(io, _alignment_text(alignment))
end

function Base.show(io::IO, ::MIME"text/plain", alignment::MultipleSequenceAlignment)
    print(io, _alignment_text(alignment))
end

function Base.show(io::IO, ::MIME"text/markdown", alignment::MultipleSequenceAlignment)
    print(io, "```text\n", _alignment_text(alignment), "\n```")
end

function Base.show(io::IO, ::MIME"text/html", alignment::MultipleSequenceAlignment)
    print(io, "<pre>", _alignment_text(alignment), "</pre>")
end

function Base.show(io::IO, ::MIME"text/latex", alignment::MultipleSequenceAlignment)
    print(io, "\\begin{verbatim}\n", _alignment_text(alignment), "\n\\end{verbatim}")
end

function Base.show(io::IO, ::MIME"application/vnd.julia.vscode.stderr", alignment::MultipleSequenceAlignment)
    rows = length(alignment)
    cols = get_alignment_length(alignment)
    print(io, _alignment_text(alignment; width=min(cols, 60)))
end

function Base.repr(alignment::MultipleSequenceAlignment)
    return sprint(show, alignment)
end

function alignment_column(alignment::MultipleSequenceAlignment, col::Integer)
    return _column_string(alignment, collect(eachindex(alignment.records)), col)
end

function consensus_sequence(alignment::MultipleSequenceAlignment; gap::Char='-', ambiguous::Char='N', threaded::Bool=true)
    isempty(alignment.records) && return ""
    width    = get_alignment_length(alignment)
    gap_byte = UInt8(gap)
    amb_byte = UInt8(ambiguous)
    consensus = Vector{UInt8}(undef, width)

    if threaded && width > 1 && Threads.nthreads() > 1
        Threads.@threads for column in 1:width
            byte_counts = zeros(Int, 256)
            @inbounds for record in alignment.records
                byte = UInt8(record.sequence[column])
                byte == gap_byte && continue
                byte_counts[Int(byte) + 1] += 1
            end

            best_byte  = amb_byte
            best_count = -1
            for b in 0x00:0xff
                cnt = byte_counts[Int(b) + 1]
                cnt == 0 && continue
                if cnt > best_count || (cnt == best_count && b < best_byte)
                    best_byte  = b
                    best_count = cnt
                end
            end
            consensus[column] = best_count < 0 ? amb_byte : best_byte
        end
    else
        @inbounds for column in 1:width
            byte_counts = zeros(Int, 256)
            for record in alignment.records
                byte = UInt8(record.sequence[column])
                byte == gap_byte && continue
                byte_counts[Int(byte) + 1] += 1
            end

            best_byte  = amb_byte
            best_count = -1
            for b in 0x00:0xff
                cnt = byte_counts[Int(b) + 1]
                cnt == 0 && continue
                if cnt > best_count || (cnt == best_count && b < best_byte)
                    best_byte  = b
                    best_count = cnt
                end
            end
            consensus[column] = best_count < 0 ? amb_byte : best_byte
        end
    end

    return String(consensus)
end

progressive_multiple_sequence_alignment(records::AbstractVector; kwargs...) = MultipleSequenceAlignment(records; kwargs...)

function _read_alignment_maf(lines::Vector{String})
    records = SeqRecordLite[]
    for line in lines
        stripped = strip(line)
        if startswith(stripped, "s ")
            parts = Base.split(stripped)
            # s src start size strand srcSize sequence
            length(parts) < 7 && continue
            src = parts[2]
            start = parse(Int, parts[3])
            size = parse(Int, parts[4])
            strand = parts[5]
            src_size = parse(Int, parts[6])
            sequence = parts[7]
            
            annotations = Dict{Symbol, Any}(
                :src => src,
                :start => start,
                :size => size,
                :strand => strand,
                :srcSize => src_size
            )
            push!(records, SeqRecordLite(sequence, identifier=src, name=src, annotations=annotations))
        end
    end
    return MultipleSequenceAlignment(records)
end

function _read_alignment_gcg(lines::Vector{String})
    records = SeqRecordLite[]
    sequence_data = IOBuffer()
    id = "gcg_sequence"
    in_sequence = false
    
    for line in lines
        if occursin("..", line)
            id = strip(Base.split(line, "..")[1])
            in_sequence = true
            continue
        end
        
        if in_sequence
            # Strip numbers and spaces
            for c in line
                if isletter(c) || c == '-' || c == '.'
                    print(sequence_data, c == '.' ? '-' : c)
                end
            end
        end
    end
    
    push!(records, SeqRecordLite(String(take!(sequence_data)), identifier=id, name=id))
    return MultipleSequenceAlignment(records)
end
multiple_sequence_alignment(records::AbstractVector; kwargs...) = MultipleSequenceAlignment(records; kwargs...)