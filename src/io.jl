"""
    parse_vcf_record(line)

Parse a single VCF text line into a structured variant record.
"""
function parse_vcf_record(line::AbstractString)
    fields = Base.split(line, '\t'; limit=7)
    length(fields) < 6 && return nothing

    try
        chrom = fields[1]
        pos = Int32(parse(Int, fields[2]))
        pos > 0 || throw(ArgumentError("VCF position must be positive"))
        id = fields[3]
        ref = fields[4]
        alt = fields[5]
        qual = fields[6] == "." ? missing : Float32(parse(Float64, fields[6]))

        return VariantTextRecord(chrom, pos, id, ref, alt, qual)
    catch err
        err isa ArgumentError && rethrow()
        throw(ArgumentError("malformed VCF record: $(line)"))
    end
end

"""
    read_vcf(input)

Read VCF records from a file path.
"""
function read_vcf(input::AbstractString)
    open(input, "r") do io
        return read_vcf(io)
    end
end

"""
    read_vcf(io)

Read VCF records from an IO stream.
"""
function read_vcf(io::IO)
    records = VariantTextRecord[]

    for raw_line in eachline(io)
        startswith(raw_line, '#') && continue
        record = parse_vcf_record(raw_line)
        record === nothing && continue
        push!(records, record)
    end

    return records
end

"""
    write_vcf(output, records)

Write VCF records to a file path.
"""
function write_vcf(output::AbstractString, records)
    open(output, "w") do io
        write_vcf(io, records)
    end

    return output
end

"""
    write_vcf(io, records)

Write VCF records to an IO stream.
"""
function write_vcf(io::IO, records)
    for record in records
        qual = record.qual === missing ? "." : string(record.qual)
        # Performance: single print per record avoids many small syscalls
        print(io, record.chrom, '\t', record.pos, '\t', record.id, '\t',
                  record.ref, '\t', record.alt, '\t', qual, '\n')
    end
    return nothing
end

"""
    BedRecord

Compact BED interval record.
"""
struct BedRecord
    chrom::String
    start::Int32
    stop::Int32
end

"""
    GffRecord

Structured GFF3 record with parsed attributes.
"""
struct GffRecord
    chrom::String
    source::String
    feature::String
    start::Int32
    stop::Int32
    score::Union{Missing,Float32}
    strand::String
    phase::Union{Missing,Int8}
    attributes::String
    attribute_map::Dict{String,Vector{String}}
end

"""
    Base.:(==)(left, right)

Test two GFF records for value equality.
"""
function Base.:(==)(left::GffRecord, right::GffRecord)
    return isequal(left.chrom, right.chrom) && isequal(left.source, right.source) && isequal(left.feature, right.feature) && isequal(left.start, right.start) && isequal(left.stop, right.stop) && isequal(left.score, right.score) && isequal(left.strand, right.strand) && isequal(left.phase, right.phase) && isequal(left.attributes, right.attributes) && isequal(left.attribute_map, right.attribute_map)
end

"""
    _parse_gff_attributes(attributes)

Parse the attribute column of a GFF record into a dictionary.
"""
function _parse_gff_attributes(attributes::AbstractString)
    parsed = Dict{String,Vector{String}}()
    stripped = strip(attributes)
    isempty(stripped) && return parsed

    for field in Base.split(stripped, ';')
        pair = strip(field)
        isempty(pair) && continue
        if occursin("=", pair)
            key, value = Base.split(pair, "=", limit=2)
            parsed[String(strip(key))] = isempty(value) ? String[""] : String[strip(entry) for entry in Base.split(value, ',')]
        else
            parsed[String(pair)] = [""]
        end
    end

    return parsed
end

"""
    _render_gff_attributes(attributes, attribute_map)

Render GFF attributes back into a text column.
"""
function _render_gff_attributes(attributes::String, attribute_map::Dict{String,Vector{String}})
    isempty(strip(attributes)) || return attributes
    isempty(attribute_map) && return "."

    parts = String[]
    for (key, values) in attribute_map
        isempty(values) ? push!(parts, key) : push!(parts, string(key, "=", join(values, ',')))
    end

    return join(parts, ';')
end

"""
    GffRecord(chrom, source, feature, start, stop, score, strand, phase, attributes, attribute_map=...)

Construct a normalized GFF record from text or parsed components.
"""
function GffRecord(
    chrom::AbstractString,
    source::AbstractString,
    feature::AbstractString,
    start::Integer,
    stop::Integer,
    score::Union{Missing,Float32},
    strand::AbstractString,
    phase::Union{Missing,Int8},
    attributes::AbstractString,
    attribute_map::AbstractDict=Dict{String,Vector{String}}(),
)
    return GffRecord(
        String(chrom),
        String(source),
        String(feature),
        Int32(start),
        Int32(stop),
        score,
        String(strand),
        phase,
        String(attributes),
        Dict{String,Vector{String}}(attribute_map),
    )
end

"""
    GenBankFeature

Parsed GenBank feature with raw and parsed location information.
"""
struct GenBankFeature
    key::String
    location::String
    qualifiers::Dict{String,Vector{String}}
    parsed_location::Any
end

"""
    Base.:(==)(left, right)

Test two GenBank features for value equality.
"""
function Base.:(==)(left::GenBankFeature, right::GenBankFeature)
    return isequal(left.key, right.key) && isequal(left.location, right.location) && isequal(left.qualifiers, right.qualifiers) && isequal(left.parsed_location, right.parsed_location)
end

GenBankFeature(key::AbstractString, location::AbstractString, qualifiers::AbstractDict=Dict{String,Vector{String}}()) = GenBankFeature(
    String(key),
    String(location),
    Dict{String,Vector{String}}(qualifiers),
    nothing,
)

GenBankFeature(key::AbstractString, location::AbstractString, qualifiers::AbstractDict, parsed_location) = GenBankFeature(
    String(key),
    String(location),
    Dict{String,Vector{String}}(qualifiers),
    parsed_location,
)

"""
    GenBankRecord

Structured GenBank record containing metadata, sequence, and features.
"""
struct GenBankRecord
    locus::String
    locus_line::String
    definition::String
    accession::String
    version::String
    keywords::String
    source::String
    organism::String
    comment::String
    sequence::String
    features::Vector{GenBankFeature}
end

"""
    Base.:(==)(left, right)

Test two GenBank records for value equality.
"""
function Base.:(==)(left::GenBankRecord, right::GenBankRecord)
    return isequal(left.locus, right.locus) && isequal(left.locus_line, right.locus_line) && isequal(left.definition, right.definition) && isequal(left.accession, right.accession) && isequal(left.version, right.version) && isequal(left.keywords, right.keywords) && isequal(left.source, right.source) && isequal(left.organism, right.organism) && isequal(left.comment, right.comment) && isequal(left.sequence, right.sequence) && isequal(left.features, right.features)
end

"""
    GenBankRecord(...)

Construct a normalized GenBank record from parsed components.
"""
function GenBankRecord(
    locus::AbstractString,
    definition::AbstractString,
    accession::AbstractString,
    version::AbstractString,
    source::AbstractString,
    organism::AbstractString,
    sequence::AbstractString,
    features::AbstractVector{GenBankFeature},
    locus_line::AbstractString="",
    keywords::AbstractString="",
    comment::AbstractString="",
)
    return GenBankRecord(
        String(locus),
        String(locus_line),
        String(definition),
        String(accession),
        String(version),
        String(keywords),
        String(source),
        String(organism),
        String(comment),
        String(sequence),
        features isa Vector{GenBankFeature} ? features : GenBankFeature[feature for feature in features],
    )
end

"""
    GenBankArrowRecord

Flattened GenBank record used for Arrow ingestion.
"""
struct GenBankArrowRecord
    locus::String
    locus_line::String
    accession::String
    version::String
    definition::String
    keywords::String
    source::String
    organism::String
    comment::String
    sequence::String
    feature_count::Int32
    feature_keys::String
    feature_locations::String
end

const _GENBANK_FEATURE_KEY_PREFIX = "     "
const _GENBANK_FEATURE_QUALIFIER_PREFIX = "                     "

"""
    parse_bed_record(line)

Parse a single BED text line into a BED record.
"""
function parse_bed_record(line::AbstractString)
    fields = Base.split(strip(line), '\t'; limit=4)
    length(fields) < 3 && return nothing

    try
        chrom = fields[1]
        start = parse(Int, fields[2])
        stop = parse(Int, fields[3])
        start >= 0 || throw(ArgumentError("BED start must be nonnegative"))
        stop > start || throw(ArgumentError("BED stop must be greater than start"))
        return BedRecord(chrom, Int32(start), Int32(stop))
    catch err
        err isa ArgumentError && rethrow()
        throw(ArgumentError("malformed BED record: $(line)"))
    end
end

"""
    read_bed(input)

Read BED records from a file path.
"""
function read_bed(input::AbstractString)
    open(input, "r") do io
        return read_bed(io)
    end
end

"""
    read_bed(io)

Read BED records from an IO stream.
"""
function read_bed(io::IO)
    records = BedRecord[]

    for (line_number, raw_line) in enumerate(eachline(io))
        startswith(raw_line, '#') && continue
        record = try
            parse_bed_record(raw_line)
        catch err
            throw(ArgumentError("malformed BED record on line $(line_number): $(err.msg)"))
        end
        record === nothing && continue
        push!(records, record)
    end

    return records
end

"""
    write_bed(output, records)

Write BED records to a file path.
"""
function write_bed(output::AbstractString, records)
    open(output, "w") do io
        write_bed(io, records)
    end

    return output
end

"""
    write_bed(io, records)

Write BED records to an IO stream.
"""
function write_bed(io::IO, records)
    for record in records
        print(io, record.chrom, '\t', record.start, '\t', record.stop, '\n')
    end
    return nothing
end

"""
    parse_gff_record(line)

Parse a single GFF3 text line into a structured record.
"""
function parse_gff_record(line::AbstractString)
    first_index = firstindex(line)
    tab1 = findnext('\t', line, first_index)
    tab1 === nothing && return nothing
    tab2 = findnext('\t', line, nextind(line, tab1))
    tab2 === nothing && return nothing
    tab3 = findnext('\t', line, nextind(line, tab2))
    tab3 === nothing && return nothing
    tab4 = findnext('\t', line, nextind(line, tab3))
    tab4 === nothing && return nothing
    tab5 = findnext('\t', line, nextind(line, tab4))
    tab5 === nothing && return nothing
    tab6 = findnext('\t', line, nextind(line, tab5))
    tab6 === nothing && return nothing
    tab7 = findnext('\t', line, nextind(line, tab6))
    tab7 === nothing && return nothing
    tab8 = findnext('\t', line, nextind(line, tab7))
    tab8 === nothing && return nothing

    try
        chrom = SubString(line, first_index, prevind(line, tab1))
        source = SubString(line, nextind(line, tab1), prevind(line, tab2))
        feature = SubString(line, nextind(line, tab2), prevind(line, tab3))
        start = parse(Int, SubString(line, nextind(line, tab3), prevind(line, tab4)))
        stop = parse(Int, SubString(line, nextind(line, tab4), prevind(line, tab5)))
        start > 0 || throw(ArgumentError("GFF start must be positive"))
        stop >= start || throw(ArgumentError("GFF stop must be >= start"))

        score_field = SubString(line, nextind(line, tab5), prevind(line, tab6))
        score = score_field == "." ? missing : Float32(parse(Float64, score_field))

        strand = SubString(line, nextind(line, tab6), prevind(line, tab7))
        strand in ("+", "-", ".") || throw(ArgumentError("GFF strand must be +, -, or ."))

        phase_field = SubString(line, nextind(line, tab7), prevind(line, tab8))
        phase = phase_field == "." ? missing : Int8(parse(Int, phase_field))

        attributes = SubString(line, nextind(line, tab8), lastindex(line))
        attribute_map = _parse_gff_attributes(attributes)

        return GffRecord(chrom, source, feature, Int32(start), Int32(stop), score, strand, phase, attributes, attribute_map)
    catch err
        err isa ArgumentError && rethrow()
        throw(ArgumentError("malformed GFF record: $(line)"))
    end
end

"""
    read_gff(input)

Read GFF records from a file path.
"""
function read_gff(input::AbstractString)
    open(input, "r") do io
        return read_gff(io)
    end
end

"""
    read_gff(io)

Read GFF records from an IO stream.
"""
function read_gff(io::IO)
    records = GffRecord[]

    for (line_number, raw_line) in enumerate(eachline(io))
        startswith(raw_line, '#') && continue
        record = try
            parse_gff_record(raw_line)
        catch err
            throw(ArgumentError("malformed GFF record on line $(line_number): $(err.msg)"))
        end
        record === nothing && continue
        push!(records, record)
    end

    return records
end

"""
    write_gff(output, records)

Write GFF records to a file path.
"""
function write_gff(output::AbstractString, records)
    open(output, "w") do io
        write_gff(io, records)
    end

    return output
end

"""
    write_gff(io, records)

Write GFF records to an IO stream.
"""
function write_gff(io::IO, records)
    for record in records
        score = record.score === missing ? "." : string(record.score)
        phase = record.phase === missing ? "." : string(record.phase)
        attributes = _render_gff_attributes(record.attributes, record.attribute_map)
        print(io, record.chrom, '\t', record.source, '\t', record.feature, '\t',
                  record.start, '\t', record.stop, '\t', score, '\t',
                  record.strand, '\t', phase, '\t', attributes, '\n')
    end
    return nothing
end

"""
    _genbank_push_qualifier!(qualifiers, key, value)

Append a qualifier value to a GenBank qualifier dictionary.
"""
function _genbank_push_qualifier!(qualifiers::Dict{String,Vector{String}}, key::String, value::String)
    push!(get!(qualifiers, key, String[]), value)
    return nothing
end

"""
    _genbank_parse_qualifier(text)

Parse a single GenBank qualifier line into a key-value pair.
"""
function _genbank_parse_qualifier(text::AbstractString)
    stripped = strip(text)
    startswith(stripped, "/") || return nothing

    payload = stripped[2:end]
    if occursin("=", payload)
        key, value = Base.split(payload, "=", limit=2)
        value = strip(value)
        if startswith(value, "\"") && endswith(value, "\"") && length(value) >= 2
            value = value[2:end-1]
        end
        return strip(key), value
    end

    return strip(payload), ""
end

"""
    parse_genbank_record(lines)

Parse a GenBank flat-file record from its raw text lines.
"""
function parse_genbank_record(lines::AbstractVector{<:AbstractString})
    isempty(lines) && return nothing

    locus = ""
    locus_line = ""
    definition = IOBuffer()
    accession = ""
    version = ""
    keywords = IOBuffer()
    source = IOBuffer()
    organism = IOBuffer()
    comment = IOBuffer()
    sequence = IOBuffer()
    features = GenBankFeature[]

    current_field = :none
    in_features = false
    in_origin = false
    current_feature_key = ""
    current_feature_location = IOBuffer()
    current_feature_qualifiers = Dict{String,Vector{String}}()
    current_qualifier_key = nothing

    function flush_feature!()
        if !isempty(current_feature_key)
            qualifiers = Dict{String,Vector{String}}(
                key => copy(values) for (key, values) in current_feature_qualifiers
            )
            location_text = strip(String(take!(current_feature_location)))
            parsed_location = isempty(location_text) ? nothing : parse_feature_location(location_text)
            push!(features, GenBankFeature(current_feature_key, location_text, qualifiers, parsed_location))
        end
        empty!(current_feature_qualifiers)
        current_feature_location = IOBuffer()
        current_feature_key = ""
        current_qualifier_key = nothing
        return nothing
    end

    for raw_line in lines
        line = rstrip(raw_line)
        line == "//" && break

        if startswith(line, "LOCUS")
            fields = Base.split(strip(line))
            length(fields) >= 2 && (locus = fields[2])
            locus_line = strip(line)
            current_field = :none
            continue
        elseif startswith(line, "DEFINITION")
            print(definition, strip(line[11:end]))
            current_field = :definition
            continue
        elseif startswith(line, "ACCESSION")
            accession = strip(line[10:end])
            current_field = :accession
            continue
        elseif startswith(line, "VERSION")
            version = strip(line[8:end])
            current_field = :version
            continue
        elseif startswith(line, "KEYWORDS")
            print(keywords, strip(line[9:end]))
            current_field = :keywords
            continue
        elseif startswith(line, "SOURCE")
            print(source, strip(line[7:end]))
            current_field = :source
            continue
        elseif startswith(line, "  ORGANISM")
            print(organism, strip(line[11:end]))
            current_field = :organism
            continue
        elseif startswith(line, "COMMENT")
            print(comment, strip(line[8:end]))
            current_field = :comment
            continue
        elseif startswith(line, "FEATURES")
            in_features = true
            in_origin = false
            current_field = :features
            continue
        elseif startswith(line, "ORIGIN")
            flush_feature!()
            in_features = false
            in_origin = true
            current_field = :origin
            continue
        end

        if in_origin
            for byte in codeunits(line)
                if (UInt8('a') <= byte <= UInt8('z')) || (UInt8('A') <= byte <= UInt8('Z'))
                    write(sequence, uppercase(Char(byte)))
                end
            end
            continue
        end

        if in_features
            if startswith(line, _GENBANK_FEATURE_QUALIFIER_PREFIX)
                qualifier_line = line[length(_GENBANK_FEATURE_QUALIFIER_PREFIX) + 1:end]
                parsed = _genbank_parse_qualifier(qualifier_line)
                if parsed !== nothing
                    key, value = parsed
                    _genbank_push_qualifier!(current_feature_qualifiers, String(key), String(value))
                    current_qualifier_key = String(key)
                end
                continue
            end

            if startswith(line, _GENBANK_FEATURE_KEY_PREFIX)
                parts = Base.split(strip(line), limit=2)
                if length(parts) >= 2
                    flush_feature!()
                    current_feature_key = parts[1]
                    print(current_feature_location, parts[2])
                    current_qualifier_key = nothing
                    continue
                end
            end

            if current_qualifier_key !== nothing
                values = current_feature_qualifiers[current_qualifier_key]
                stripped = strip(line)
                isempty(stripped) || (values[end] = string(values[end], " ", stripped))
                continue
            end

            stripped = strip(line)
            if !isempty(stripped)
                write(current_feature_location, ' ')
                write(current_feature_location, stripped)
            end
            continue
        end

        if current_field == :definition
            print(definition, " ", strip(line))
        elseif current_field == :keywords
            print(keywords, " ", strip(line))
        elseif current_field == :source
            print(source, " ", strip(line))
        elseif current_field == :organism
            print(organism, " ", strip(line))
        elseif current_field == :comment
            print(comment, " ", strip(line))
        end
    end

    flush_feature!()

    keywords_text = strip(String(take!(keywords)))
    if endswith(keywords_text, ".")
        keywords_text = strip(keywords_text[1:end-1])
    end

    return GenBankRecord(
        locus,
        locus_line,
        String(take!(definition)),
        accession,
        version,
        keywords_text,
        String(take!(source)),
        String(take!(organism)),
        String(take!(comment)),
        String(take!(sequence)),
        features,
    )
end

"""
    _wrap_genbank_text(prefix, text; continuation_prefix=..., width=79, empty_text="")

Wrap a GenBank text field to the expected line width.
"""
function _wrap_genbank_text(prefix::AbstractString, text::AbstractString; continuation_prefix::AbstractString=repeat(" ", ncodeunits(prefix)), width::Int=79, empty_text::AbstractString="")
    payload = strip(text)
    isempty(payload) && return isempty(empty_text) ? String[] : [string(prefix, empty_text)]

    words = Base.split(payload)
    lines = String[]
    current_prefix = String(prefix)
    current = String(prefix)
    current_length = ncodeunits(current_prefix)

    for word in words
        needs_space = current_length > ncodeunits(current_prefix)
        projected = current_length + (needs_space ? 1 : 0) + ncodeunits(word)
        if projected > width && current != current_prefix
            push!(lines, current)
            current_prefix = String(continuation_prefix)
            current = string(current_prefix, word)
            current_length = ncodeunits(current)
        else
            if needs_space
                current = string(current, " ")
                current_length += 1
            end
            current = string(current, word)
            current_length += ncodeunits(word)
        end
    end

    push!(lines, current)
    return lines
end

"""
    _format_genbank_locus_line(record)

Render the LOCUS line for a GenBank record.
"""
function _format_genbank_locus_line(record::GenBankRecord)
    isempty(strip(record.locus_line)) && return string("LOCUS       ", rpad(record.locus, 16), lpad(string(ncodeunits(record.sequence)), 11), " bp    DNA     linear   UNK")
    return record.locus_line
end

"""
    _write_genbank_wrapped_section(io, prefix, text; continuation_prefix=..., width=79, empty_text="")

Write a wrapped GenBank text section to an IO stream.
"""
function _write_genbank_wrapped_section(io::IO, prefix::AbstractString, text::AbstractString; continuation_prefix::AbstractString=repeat(" ", ncodeunits(prefix)), width::Int=79, empty_text::AbstractString="")
    for line in _wrap_genbank_text(prefix, text; continuation_prefix=continuation_prefix, width=width, empty_text=empty_text)
        isempty(line) && continue
        println(io, line)
    end
end

"""
    _render_genbank_qualifier(key, values)

Render a GenBank qualifier and its values as text.
"""
function _render_genbank_qualifier(key::String, values::Vector{String})
    isempty(values) && return "                     /$(key)"
    if length(values) == 1
        value = values[1]
        isempty(value) && return "                     /$(key)"
        return "                     /$(key)=\"$(value)\""
    end
    return join(("                     /$(key)=\"$(value)\"" for value in values), '\n')
end

"""
    _render_genbank_sequence(sequence)

Render a GenBank ORIGIN sequence block.
"""
function _render_genbank_sequence(sequence::AbstractString)
    buffer = IOBuffer()
    for index in 1:60:length(sequence)
        chunk = sequence[index:min(index + 59, lastindex(sequence))]
        print(buffer, lpad(string(div(index - 1, 60) + 1), 9), " ")
        for chunk_index in 1:10:length(chunk)
            print(buffer, lowercase(chunk[chunk_index:min(chunk_index + 9, lastindex(chunk))]), ' ')
        end
        print(buffer, '\n')
    end
    return String(take!(buffer))
end

"""
    write_genbank(output_path, records)

Write GenBank records to a file path.
"""
function write_genbank(output_path::AbstractString, records)
    open(output_path, "w") do io
        write_genbank(io, records)
    end

    return output_path
end

"""
    write_genbank(io, records)

Write GenBank records to an IO stream.
"""
function write_genbank(io::IO, records)
    for record in records
        println(io, _format_genbank_locus_line(record))
        _write_genbank_wrapped_section(io, "DEFINITION  ", record.definition; continuation_prefix="            ")
        println(io, "ACCESSION   ", record.accession)
        println(io, "VERSION     ", record.version)
        keyword_text = strip(record.keywords)
        isempty(keyword_text) ? println(io, "KEYWORDS    .") : _write_genbank_wrapped_section(io, "KEYWORDS    ", endswith(keyword_text, ".") ? keyword_text : string(keyword_text, "."); continuation_prefix="            ")
        _write_genbank_wrapped_section(io, "SOURCE      ", record.source; continuation_prefix="            ")
        _write_genbank_wrapped_section(io, "  ORGANISM  ", record.organism; continuation_prefix="            ")
        isempty(strip(record.comment)) || _write_genbank_wrapped_section(io, "COMMENT     ", record.comment; continuation_prefix="            ")
        println(io, "FEATURES             Location/Qualifiers")
        for feature in record.features
            println(io, "     ", rpad(feature.key, 15), feature.location)
            for (key, values) in feature.qualifiers
                for rendered in Base.split(_render_genbank_qualifier(key, values), '\n')
                    println(io, rendered)
                end
            end
        end
        println(io, "ORIGIN")
        print(io, _render_genbank_sequence(record.sequence))
        println(io, "//")
    end
    return nothing
end

"""
    _genbank_flatten_record(record)

Flatten a GenBank record into an Arrow-ready storage layout.
"""
function _genbank_flatten_record(record::GenBankRecord)
    feature_keys = join((feature.key for feature in record.features), ";")
    feature_locations = join((feature.location for feature in record.features), ";")
    return GenBankArrowRecord(
        record.locus,
        record.locus_line,
        record.accession,
        record.version,
        record.definition,
        record.keywords,
        record.source,
        record.organism,
        record.comment,
        record.sequence,
        Int32(length(record.features)),
        feature_keys,
        feature_locations,
    )
end

"""
    read_genbank(input_path)

Read GenBank records from a file path.
"""
function read_genbank(input_path::AbstractString)
    records = GenBankRecord[]
    current_lines = String[]

    open(input_path, "r") do io
        for raw_line in eachline(io)
            push!(current_lines, raw_line)
            if strip(raw_line) == "//"
                record = parse_genbank_record(current_lines)
                record === nothing || push!(records, record)
                empty!(current_lines)
            end
        end
    end

    isempty(current_lines) || push!(records, parse_genbank_record(current_lines))
    return records
end

"""
    _write_genbank_chunk!(writer, loci, locus_lines, accessions, versions, definitions, keywords, sources, organisms, comments, sequences, feature_counts, feature_keys, feature_locations)

Write a chunk of flattened GenBank records to an Arrow writer.
"""
function _write_genbank_chunk!(writer, loci, locus_lines, accessions, versions, definitions, keywords, sources, organisms, comments, sequences, feature_counts, feature_keys, feature_locations)
    isempty(loci) && return nothing

    Arrow.write(
        writer,
        (
            locus = copy(loci),
            locus_line = copy(locus_lines),
            accession = copy(accessions),
            version = copy(versions),
            definition = copy(definitions),
            keywords = copy(keywords),
            source = copy(sources),
            organism = copy(organisms),
            comment = copy(comments),
            sequence = copy(sequences),
            feature_count = copy(feature_counts),
            feature_keys = copy(feature_keys),
            feature_locations = copy(feature_locations),
        ),
    )

    empty!(loci)
    empty!(locus_lines)
    empty!(accessions)
    empty!(versions)
    empty!(definitions)
    empty!(keywords)
    empty!(sources)
    empty!(organisms)
    empty!(comments)
    empty!(sequences)
    empty!(feature_counts)
    empty!(feature_keys)
    empty!(feature_locations)

    return nothing
end

"""
    ingest_genbank(input_path, output_path; chunk_size=100)

Convert GenBank records into a chunked Arrow table.
"""
function ingest_genbank(
    input_path::AbstractString,
    output_path::AbstractString;
    chunk_size::Integer=100,
)
    loci = String[]
    locus_lines = String[]
    accessions = String[]
    versions = String[]
    definitions = String[]
    keywords = String[]
    sources = String[]
    organisms = String[]
    comments = String[]
    sequences = String[]
    feature_counts = Int32[]
    feature_keys = String[]
    feature_locations = String[]

    current_lines = String[]

    open(Arrow.Writer, output_path; file=true) do writer
        open(input_path, "r") do io
            for raw_line in eachline(io)
                push!(current_lines, raw_line)
                if strip(raw_line) == "//"
                    record = parse_genbank_record(current_lines)
                    if record !== nothing
                        flattened = _genbank_flatten_record(record)
                        push!(loci, flattened.locus)
                        push!(locus_lines, flattened.locus_line)
                        push!(accessions, flattened.accession)
                        push!(versions, flattened.version)
                        push!(definitions, flattened.definition)
                        push!(keywords, flattened.keywords)
                        push!(sources, flattened.source)
                        push!(organisms, flattened.organism)
                        push!(comments, flattened.comment)
                        push!(sequences, flattened.sequence)
                        push!(feature_counts, flattened.feature_count)
                        push!(feature_keys, flattened.feature_keys)
                        push!(feature_locations, flattened.feature_locations)

                        if length(loci) >= chunk_size
                            _write_genbank_chunk!(writer, loci, locus_lines, accessions, versions, definitions, keywords, sources, organisms, comments, sequences, feature_counts, feature_keys, feature_locations)
                        end
                    end
                    empty!(current_lines)
                end
            end
        end

        if !isempty(current_lines)
            record = parse_genbank_record(current_lines)
            if record !== nothing
                flattened = _genbank_flatten_record(record)
                push!(loci, flattened.locus)
                push!(locus_lines, flattened.locus_line)
                push!(accessions, flattened.accession)
                push!(versions, flattened.version)
                push!(definitions, flattened.definition)
                push!(keywords, flattened.keywords)
                push!(sources, flattened.source)
                push!(organisms, flattened.organism)
                push!(comments, flattened.comment)
                push!(sequences, flattened.sequence)
                push!(feature_counts, flattened.feature_count)
                push!(feature_keys, flattened.feature_keys)
                push!(feature_locations, flattened.feature_locations)
            end
        end

        _write_genbank_chunk!(writer, loci, locus_lines, accessions, versions, definitions, keywords, sources, organisms, comments, sequences, feature_counts, feature_keys, feature_locations)
    end

    return output_path
end

"""
    _write_vcf_chunk!(writer, chroms, positions, ids, refs, alts, quals)

Write a chunk of flattened VCF records to an Arrow writer.
"""
function _write_vcf_chunk!(writer, chroms, positions, ids, refs, alts, quals)
    isempty(chroms) && return nothing

    chunk = (
        chrom = copy(chroms),
        pos = copy(positions),
        id = copy(ids),
        ref = copy(refs),
        alt = copy(alts),
        qual = copy(quals),
    )

    Arrow.write(
        writer,
        chunk,
    )

    empty!(chroms)
    empty!(positions)
    empty!(ids)
    empty!(refs)
    empty!(alts)
    empty!(quals)

    return nothing
end

"""
    _write_bed_chunk!(writer, chroms, starts, stops)

Write a chunk of flattened BED records to an Arrow writer.
"""
function _write_bed_chunk!(writer, chroms, starts, stops)
    isempty(chroms) && return nothing

    Arrow.write(
        writer,
        (
            chrom = copy(chroms),
            start = copy(starts),
            stop = copy(stops),
        ),
    )

    empty!(chroms)
    empty!(starts)
    empty!(stops)

    return nothing
end

"""
    ingest_vcf(input_path, output_path; chunk_size=10_000)

Convert VCF records into a chunked Arrow table.
"""
function ingest_vcf(
    input_path::AbstractString,
    output_path::AbstractString;
    chunk_size::Integer=10_000,
)
    chroms = String[]
    positions = Int32[]
    ids = String[]
    refs = String[]
    alts = String[]
    quals = Union{Missing,Float32}[]

    open(Arrow.Writer, output_path; file=true) do writer
        open(input_path, "r") do io
            for raw_line in eachline(io)
                startswith(raw_line, '#') && continue
                record = parse_vcf_record(raw_line)
                record === nothing && continue

                push!(chroms, record.chrom)
                push!(positions, record.pos)
                push!(ids, record.id)
                push!(refs, record.ref)
                push!(alts, record.alt)
                push!(quals, record.qual)

                if length(chroms) >= chunk_size
                    _write_vcf_chunk!(writer, chroms, positions, ids, refs, alts, quals)
                end
            end
        end

        _write_vcf_chunk!(writer, chroms, positions, ids, refs, alts, quals)
    end

    return output_path
end

"""
    ingest_bed(input_path, output_path; chunk_size=10_000)

Convert BED records into a chunked Arrow table.
"""
function ingest_bed(
    input_path::AbstractString,
    output_path::AbstractString;
    chunk_size::Integer=10_000,
)
    chroms = String[]
    starts = Int32[]
    stops = Int32[]

    open(Arrow.Writer, output_path; file=true) do writer
        open(input_path, "r") do io
            for raw_line in eachline(io)
                startswith(raw_line, '#') && continue
                record = parse_bed_record(raw_line)
                record === nothing && continue

                push!(chroms, record.chrom)
                push!(starts, record.start)
                push!(stops, record.stop)

                if length(chroms) >= chunk_size
                    _write_bed_chunk!(writer, chroms, starts, stops)
                end
            end
        end

        _write_bed_chunk!(writer, chroms, starts, stops)
    end

    return output_path
end

"""
    _write_gff_chunk!(writer, chroms, sources, features, starts, stops, scores, strands, phases, attributes)

Write a chunk of flattened GFF records to an Arrow writer.
"""
function _write_gff_chunk!(writer, chroms, sources, features, starts, stops, scores, strands, phases, attributes)
    isempty(chroms) && return nothing

    Arrow.write(
        writer,
        (
            chrom = copy(chroms),
            source = copy(sources),
            feature = copy(features),
            start = copy(starts),
            stop = copy(stops),
            score = copy(scores),
            strand = copy(strands),
            phase = copy(phases),
            attributes = copy(attributes),
        ),
    )

    empty!(chroms)
    empty!(sources)
    empty!(features)
    empty!(starts)
    empty!(stops)
    empty!(scores)
    empty!(strands)
    empty!(phases)
    empty!(attributes)

    return nothing
end

"""
    ingest_gff(input_path, output_path; chunk_size=10_000)

Convert GFF records into a chunked Arrow table.
"""
function ingest_gff(
    input_path::AbstractString,
    output_path::AbstractString;
    chunk_size::Integer=10_000,
)
    chroms = String[]
    sources = String[]
    features = String[]
    starts = Int32[]
    stops = Int32[]
    scores = Union{Missing,Float32}[]
    strands = String[]
    phases = Union{Missing,Int8}[]
    attributes = String[]

    open(Arrow.Writer, output_path; file=true) do writer
        open(input_path, "r") do io
            for raw_line in eachline(io)
                startswith(raw_line, '#') && continue
                record = parse_gff_record(raw_line)
                record === nothing && continue

                push!(chroms, record.chrom)
                push!(sources, record.source)
                push!(features, record.feature)
                push!(starts, record.start)
                push!(stops, record.stop)
                push!(scores, record.score)
                push!(strands, record.strand)
                push!(phases, record.phase)
                push!(attributes, record.attributes)

                if length(chroms) >= chunk_size
                    _write_gff_chunk!(writer, chroms, sources, features, starts, stops, scores, strands, phases, attributes)
                end
            end
        end

        _write_gff_chunk!(writer, chroms, sources, features, starts, stops, scores, strands, phases, attributes)
    end

    return output_path
end

"""
    load_arrow_table(path)

Load an Arrow table from disk.
"""
function load_arrow_table(path::AbstractString)
    return Arrow.Table(path)
end

"""
    write_arrow_table(output_path, table)

Write a table to an Arrow file.
"""
function write_arrow_table(output_path::AbstractString, table)
    Arrow.write(output_path, table)
    return output_path
end

# 10. EMBL & SwissProt Support
# ------------------------------------------------------------------------------

"""
    read_embl(filepath::String)

Parses an EMBL flat file into a vector of AnnotatedSeqRecord objects.
"""
function read_embl(filepath::String)
    open(filepath, "r") do io
        return read_embl(io)
    end
end

"""
    read_embl(io)

Read EMBL records from an IO stream.
"""
function read_embl(io::IO)
    lines = readlines(io)
    records = AnnotatedSeqRecord[]
    sizehint!(records, 1)

    record_start = 1

    for (index, line) in pairs(lines)
        if startswith(line, "//")
            if record_start <= index - 1
                push!(records, _parse_embl_record(@view lines[record_start:index - 1]))
            end
            record_start = index + 1
        end
    end

    if record_start <= length(lines)
        push!(records, _parse_embl_record(@view lines[record_start:end]))
    end

    return records
end

"""
    _embl_first_token(text)

Return the first whitespace-delimited token from an EMBL field.
"""
@inline function _embl_first_token(text::AbstractString)
    token_end = findfirst(isspace, text)
    token_end === nothing && return String(text)
    return String(text[firstindex(text):prevind(text, token_end)])
end

"""
    _embl_compact_tail(text)

Strip spaces from an EMBL continuation field.
"""
@inline function _embl_compact_tail(text::AbstractString)
    isempty(text) && return ""
    buffer = IOBuffer()
    for character in text
        character == ' ' && continue
        print(buffer, character)
    end
    return String(take!(buffer))
end

"""
    _parse_embl_record(lines)

Parse a single EMBL record from raw text lines.
"""
function _parse_embl_record(lines::AbstractVector{<:AbstractString})
    id = ""
    accession = ""
    description = IOBuffer()
    sequence = IOBuffer()
    features = SeqFeatureLite[]
    annotations = Dict{Symbol, Any}()
    
    current_feature_type = ""
    current_feature_loc = ""
    current_feature_qualifiers = Dict{String, Vector{String}}()
    
    in_sequence = false
    
    for line in lines
        if length(line) < 2; continue; end
        code = line[1:2]
        content = length(line) > 5 ? strip(line[6:end]) : ""
        
        if code == "ID"
            id = _embl_first_token(content)
        elseif code == "AC"
            accession = isempty(accession) ? content : accession * " " * content
        elseif code == "DE"
            print(description, content, " ")
        elseif code == "FT"
            # Feature Table
            # EMBL format: 
            # FT   feature_name    location
            # FT                   /qualifier="value"
            # FT                   location_continuation
            
            # Check if this is a qualifier line
            line_content = line[6:end]
            trimmed_content = strip(line_content)
            
            if startswith(trimmed_content, "/")
                # Qualifier
                qual_part = trimmed_content[2:end] # remove /
                eq_index = findfirst(==( '=' ), qual_part)
                if eq_index !== nothing
                    q_key = String(qual_part[firstindex(qual_part):prevind(qual_part, eq_index)])
                    q_val = strip(String(qual_part[nextind(qual_part, eq_index):end]), '"')
                    push!(get!(current_feature_qualifiers, q_key, String[]), q_val)
                else
                    push!(get!(current_feature_qualifiers, String(qual_part), String[]), "")
                end
            elseif !isempty(trimmed_content) && isspace(line[6]) && isspace(line[7]) && isspace(line[8])
                current_feature_loc *= trimmed_content
            elseif !isempty(trimmed_content)
                space_index = findfirst(isspace, trimmed_content)
                if space_index !== nothing
                    if !isempty(current_feature_type)
                        loc_obj = parse_feature_location(current_feature_loc)
                        push!(features, SeqFeatureLite(current_feature_type, loc_obj, qualifiers=current_feature_qualifiers))
                    end
                    current_feature_type = String(trimmed_content[firstindex(trimmed_content):prevind(trimmed_content, space_index)])
                    current_feature_loc = _embl_compact_tail(trimmed_content[nextind(trimmed_content, space_index):end])
                    current_feature_qualifiers = Dict{String, Vector{String}}()
                end
            end
        elseif code == "SQ"
            in_sequence = true
        elseif in_sequence && code != "//"
            # Sequence data: strip spaces and numbers
            for c in content
                if isletter(c)
                    print(sequence, c)
                end
            end
        end
    end
    
    # Flush last feature
    if !isempty(current_feature_type)
        loc_obj = parse_feature_location(current_feature_loc)
        push!(features, SeqFeatureLite(current_feature_type, loc_obj, qualifiers=current_feature_qualifiers))
    end
    
    if !isempty(accession)
        annotations[:accession] = accession
    end

    return AnnotatedSeqRecord(
        String(take!(sequence));
        identifier=id,
        name=id,
        description=strip(String(take!(description))),
        annotations=annotations,
        features=features,
    )
end

"""
    read_swissprot(filepath::String)

Parses a SwissProt/UniProt .dat file into a vector of AnnotatedSeqRecord objects.
"""
function read_swissprot(filepath::String)
    return read_embl(filepath) # Formats are extremely similar in structure
end

# --- ABIF (.ab1) Sanger Trace Parser ---

"""
    read_abif(filepath)

Read an ABIF chromatogram file from a path.
"""
function read_abif(filepath::AbstractString)
    open(filepath, "r") do io
        return read_abif(io)
    end
end

"""
    read_abif(io)

Read an ABIF chromatogram file from an IO stream.
"""
function read_abif(io::IO)
    magic = String(read(io, 4))
    if magic != "ABIF"
        throw(ArgumentError("Not a valid ABIF file (magic: \$magic)"))
    end
    
    version = ntoh(read(io, UInt16))
    
    # Directory entry (28 bytes) begins at offset 6
    dir_name = String(read(io, 4))
    dir_number = ntoh(read(io, UInt32))
    dir_type = ntoh(read(io, UInt16))
    dir_elem_size = ntoh(read(io, UInt16))
    dir_elements = ntoh(read(io, UInt32))
    dir_data_size = ntoh(read(io, UInt32))
    dir_data_offset = ntoh(read(io, UInt32))
    
    # Header is 128 bytes total, seek past it
    seek(io, 128)
    
    # Then seek to the actual directory location
    seek(io, dir_data_offset)
    sequence_bytes = UInt8[]
    qualities = UInt8[]
    fwo = "GATC"
    trace_data = Dict{Int, Vector{UInt16}}()

    @inline function read_tag_bytes(size::UInt32, data_offset::UInt32, data_pos::Int)
        if size <= 4
            seek(io, data_pos - 4)
            data = read(io, 4)[1:size]
        else
            seek(io, data_offset)
            data = read(io, size)
        end
        seek(io, data_pos)
        return data
    end
    
    for i in 1:dir_elements
        tag_name = String(read(io, 4))
        tag_num = ntoh(read(io, UInt32))
        tag_type = ntoh(read(io, UInt16))
        tag_elem_size = ntoh(read(io, UInt16))
        tag_num_elems = ntoh(read(io, UInt32))
        tag_data_size = ntoh(read(io, UInt32))
        tag_data_offset = ntoh(read(io, UInt32))
        
        pos = position(io)

        if tag_name == "PBAS" && (tag_num == 1 || tag_num == 2)
            sequence_bytes = read_tag_bytes(tag_data_size, tag_data_offset, pos)
        elseif (tag_name == "PQC" || tag_name == "PCON") && (tag_num == 1 || tag_num == 2)
            qualities = read_tag_bytes(tag_data_size, tag_data_offset, pos)
        elseif tag_name == "FWO_" && tag_num == 1
            fwo = String(read_tag_bytes(tag_data_size, tag_data_offset, pos))
        elseif tag_name == "DATA"
            trace_data[Int(tag_num)] = [ntoh(val) for val in reinterpret(UInt16, read_tag_bytes(tag_data_size, tag_data_offset, pos))]
        else
            seek(io, pos)
        end
    end
    
    sequence = isempty(sequence_bytes) ? "" : String(sequence_bytes)
    qualities = qualities
    
    trace_map = Dict{Char, Vector{UInt16}}()
    for (i, char) in enumerate(fwo)
        trace_map[char] = get(trace_data, i, UInt16[])
    end
    
    return SangerTrace(
        sequence,
        qualities,
        get(trace_map, 'A', UInt16[]),
        get(trace_map, 'C', UInt16[]),
        get(trace_map, 'G', UInt16[]),
        get(trace_map, 'T', UInt16[]),
        Dict{String, Any}("fwo" => fwo, "version" => version)
    )
end
