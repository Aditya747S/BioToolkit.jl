using DataFrames

struct MotifCounts
    alphabet::Vector{Char}
    counts::Matrix{Int}
end

struct MotifFrequencyMatrix
    alphabet::Vector{Char}
    values::Matrix{Float64}
end

struct MotifPWM
    alphabet::Vector{Char}
    values::Matrix{Float64}
end

struct MotifHit
    start::Int
    strand::Int8
    score::Float64
    window::String
end

struct MotifSite
    sequence_index::Int
    start::Int
    strand::Int8
    mismatches::Int
    window::String
end

struct MotifDiscoveryResult
    seed::String
    alphabet::Vector{Char}
    counts::MotifCounts
    pwm::MotifPWM
    sites::Vector{MotifSite}
    support::Int
    information_content::Float64
end

const _MOTIF_CODE = fill(UInt8(255), 256)

for byte in (UInt8('A'), UInt8('a'))
    _MOTIF_CODE[Int(byte) + 1] = UInt8(1)
end
for byte in (UInt8('C'), UInt8('c'))
    _MOTIF_CODE[Int(byte) + 1] = UInt8(2)
end
for byte in (UInt8('G'), UInt8('g'))
    _MOTIF_CODE[Int(byte) + 1] = UInt8(3)
end
for byte in (UInt8('T'), UInt8('t'), UInt8('U'), UInt8('u'))
    _MOTIF_CODE[Int(byte) + 1] = UInt8(4)
end

function Base.show(io::IO, profile::MotifCounts)
    print(io, "MotifCounts(", join(profile.alphabet), ", ", size(profile.counts, 2), " positions)")
end

function Base.show(io::IO, profile::MotifFrequencyMatrix)
    print(io, "MotifFrequencyMatrix(", join(profile.alphabet), ", ", size(profile.values, 2), " positions)")
end

function Base.show(io::IO, pwm::MotifPWM)
    print(io, "MotifPWM(", join(pwm.alphabet), ", ", size(pwm.values, 2), " positions)")
end

function Base.show(io::IO, hit::MotifHit)
    print(io, "MotifHit(start=", hit.start, ", strand=", hit.strand, ", score=", round(hit.score, digits=3), ")")
end

function Base.show(io::IO, site::MotifSite)
    print(io, "MotifSite(seq=", site.sequence_index, ", start=", site.start, ", strand=", site.strand, ", mismatches=", site.mismatches, ")")
end

function Base.show(io::IO, result::MotifDiscoveryResult)
    print(io, "MotifDiscoveryResult(seed=", result.seed, ", support=", result.support, ", ic=", round(result.information_content, digits=3), ")")
end

function DataFrames.DataFrame(hits::AbstractVector{<:MotifHit})
    rows = collect(hits)
    return DataFrames.DataFrame(
        start = [hit.start for hit in rows],
        strand = [hit.strand for hit in rows],
        score = [hit.score for hit in rows],
        window = [hit.window for hit in rows],
    )
end

function DataFrames.DataFrame(sites::AbstractVector{<:MotifSite})
    rows = collect(sites)
    return DataFrames.DataFrame(
        sequence_index = [site.sequence_index for site in rows],
        start = [site.start for site in rows],
        strand = [site.strand for site in rows],
        mismatches = [site.mismatches for site in rows],
        window = [site.window for site in rows],
    )
end

function _motif_alphabet_index(alphabet::AbstractString)
    alphabet_chars = collect(uppercase(alphabet))
    index = Dict{Char,Int}()

    for (position, character) in pairs(alphabet_chars)
        index[character] = position
    end

    return alphabet_chars, index
end

function motif_counts(sequences::AbstractVector{<:AbstractString}; alphabet::AbstractString="ACGT")
    isempty(sequences) && throw(ArgumentError("at least one sequence is required"))

    alphabet_chars, alphabet_index = _motif_alphabet_index(alphabet)
    motif_length = ncodeunits(first(sequences))
    counts = zeros(Int, length(alphabet_chars), motif_length)
    row_lookup = zeros(Int, 4)

    for (position, character) in pairs(alphabet_chars)
        character == 'A' && (row_lookup[1] = position)
        character == 'C' && (row_lookup[2] = position)
        character == 'G' && (row_lookup[3] = position)
        character == 'T' && (row_lookup[4] = position)
        character == 'U' && (row_lookup[4] = position)
    end

    for sequence in sequences
        ncodeunits(sequence) == motif_length || throw(ArgumentError("all motif sequences must have the same length"))

        @inbounds for (position, byte) in enumerate(codeunits(sequence))
            code = _MOTIF_CODE[Int(byte) + 1]
            code == 0xff && throw(ArgumentError("unsupported motif character $(Char(byte))"))
            row = row_lookup[Int(code)]
            row == 0 && throw(ArgumentError("unsupported motif character $(Char(byte))"))
            counts[row, position] += 1
        end
    end

    return MotifCounts(alphabet_chars, counts)
end

function motif_counts(records::AbstractVector{SeqRecordLite}; alphabet::AbstractString="ACGT")
    return motif_counts([record.sequence for record in records]; alphabet=alphabet)
end

function _motif_frequency_matrix(profile::MotifCounts; pseudocount::Real=0.0)
    pseudocount >= 0 || throw(ArgumentError("pseudocount must be nonnegative"))
    nrows = size(profile.counts, 1)
    ncols = size(profile.counts, 2)
    values = Matrix{Float64}(undef, nrows, ncols)

    @inbounds for column in 1:ncols
        column_total = 0
        for row in 1:nrows
            column_total += profile.counts[row, column]
        end
        denominator = column_total + pseudocount * nrows
        denominator > 0 || throw(ArgumentError("motif column has no mass"))

        for row in 1:nrows
            values[row, column] = (profile.counts[row, column] + pseudocount) / denominator
        end
    end

    return MotifFrequencyMatrix(profile.alphabet, values)
end

function motif_frequency_matrix(profile::MotifCounts; pseudocount::Real=0.0)
    return _motif_frequency_matrix(profile; pseudocount=pseudocount)
end

function motif_frequency_matrix(sequences::AbstractVector{<:AbstractString}; alphabet::AbstractString="ACGT", pseudocount::Real=0.0)
    return motif_frequency_matrix(motif_counts(sequences; alphabet=alphabet); pseudocount=pseudocount)
end

function motif_frequency_matrix(records::AbstractVector{SeqRecordLite}; alphabet::AbstractString="ACGT", pseudocount::Real=0.0)
    return motif_frequency_matrix(motif_counts(records; alphabet=alphabet); pseudocount=pseudocount)
end

function motif_entropy(profile::MotifFrequencyMatrix)
    ncols = size(profile.values, 2)
    total = 0.0

    @inbounds for column in axes(profile.values, 2)
        column_entropy = 0.0
        for row in eachindex(profile.alphabet)
            probability = profile.values[row, column]
            probability <= 0 && continue
            column_entropy -= probability * log2(probability)
        end
        total += column_entropy
    end

    return ncols == 0 ? 0.0 : total / ncols
end

function motif_relative_entropy(profile::MotifFrequencyMatrix; background=nothing)
    background_vector = _motif_background_vector(profile.alphabet, background)
    ncols = size(profile.values, 2)
    total = 0.0

    @inbounds for column in axes(profile.values, 2)
        column_re = 0.0
        for row in eachindex(profile.alphabet)
            probability = profile.values[row, column]
            probability <= 0 && continue
            background_probability = background_vector[row]
            background_probability > 0 || throw(ArgumentError("background entries must be positive"))
            column_re += probability * log2(probability / background_probability)
        end
        total += column_re
    end

    return ncols == 0 ? 0.0 : total / ncols
end

@inline function _motif_hamming_distance(left::AbstractString, right::AbstractString)
    ncodeunits(left) == ncodeunits(right) || throw(ArgumentError("motif windows must have the same length"))
    mismatches = 0
    @inbounds for (left_byte, right_byte) in zip(codeunits(left), codeunits(right))
        left_byte == right_byte && continue
        mismatches += 1
    end
    return mismatches
end

function _motif_lookup_array(alphabet::Vector{Char})
    lookup = zeros(Int, 256)
    @inbounds for (position, character) in enumerate(alphabet)
        lookup[Int(UInt8(character)) + 1] = position
    end
    return lookup
end

function _motif_canonical_kmer(kmer::AbstractString; reverse_complements::Bool=true)
    reverse_complements || return String(kmer)
    reverse_kmer = reverse_complement(kmer)
    return reverse_kmer < kmer ? reverse_kmer : String(kmer)
end

function _motif_background_vector_from_pwm(alphabet::Vector{Char}, background)
    return _motif_background_vector(alphabet, background)
end

function _motif_seed_counts(sequences::AbstractVector{<:AbstractString}, k::Int; reverse_complements::Bool=true)
    seed_counts = Dict{String,Int}()
    for sequence in sequences
        sequence_length = ncodeunits(sequence)
        sequence_length >= k || throw(ArgumentError("motif length must be no greater than each sequence length"))
        for start_index in 1:(sequence_length - k + 1)
            window = String(sequence[start_index:start_index + k - 1])
            canonical = _motif_canonical_kmer(window; reverse_complements=reverse_complements)
            seed_counts[canonical] = get(seed_counts, canonical, 0) + 1
        end
    end
    return seed_counts
end

function _motif_best_site(sequence::AbstractString, seed::AbstractString, seed_rc::AbstractString, k::Int, max_mismatches::Int)
    sequence_length = ncodeunits(sequence)
    best_site = nothing
    best_mismatches = max_mismatches + 1
    best_strand = Int8(1)

    for start_index in 1:(sequence_length - k + 1)
        window = sequence[start_index:start_index + k - 1]
        mismatches = _motif_hamming_distance(window, seed)
        strand = Int8(1)
        oriented_window = String(window)

        if mismatches > max_mismatches
            rc_mismatches = _motif_hamming_distance(window, seed_rc)
            if rc_mismatches <= max_mismatches
                mismatches = rc_mismatches
                strand = Int8(-1)
                oriented_window = reverse_complement(window)
            else
                continue
            end
        end

        if best_site === nothing || mismatches < best_mismatches
            best_site = MotifSite(0, start_index, strand, mismatches, oriented_window)
            best_mismatches = mismatches
            best_strand = strand
        end
    end

    best_site === nothing && return nothing
    return MotifSite(0, best_site.start, best_strand, best_site.mismatches, best_site.window)
end

function _motif_information_content(profile::MotifPWM; background=nothing)
    background_vector = _motif_background_vector_from_pwm(profile.alphabet, background)
    total = 0.0
    ncols = size(profile.values, 2)

    for column in axes(profile.values, 2)
        column_ic = 0.0
        for row in eachindex(profile.alphabet)
            probability = profile.values[row, column]
            probability <= 0 && continue
            background_probability = background_vector[row]
            background_probability > 0 || throw(ArgumentError("background entries must be positive"))
            column_ic += probability * log2(probability / background_probability)
        end
        total += column_ic
    end

    return ncols == 0 ? 0.0 : total / ncols
end

function motif_information_content(profile::MotifPWM; background=nothing)
    return _motif_information_content(profile; background=background)
end

function motif_scan(sequence::AbstractString, pwm::MotifPWM; threshold::Real=0.0, background=nothing)
    threshold isa Real || throw(ArgumentError("threshold must be real"))
    motif_length = size(pwm.values, 2)
    sequence_bytes = codeunits(sequence)
    sequence_length = length(sequence_bytes)
    sequence_length >= motif_length || return MotifHit[]

    row_lookup = _motif_lookup_array(pwm.alphabet)
    background_vector = _motif_background_vector_from_pwm(pwm.alphabet, background)
    hits = MotifHit[]
    window_buffer = Vector{UInt8}(undef, motif_length)

    for start_index in 1:(sequence_length - motif_length + 1)
        score = 0.0
        valid = true

        @inbounds for column_index in 1:motif_length
            byte = sequence_bytes[start_index + column_index - 1]
            row_index = row_lookup[Int(byte) + 1]
            if row_index == 0
                valid = false
                break
            end

            probability = pwm.values[row_index, column_index]
            probability <= 0 && (valid = false; break)
            background_probability = background_vector[row_index]
            background_probability > 0 || throw(ArgumentError("background entries must be positive"))
            score += log2(probability / background_probability)
            window_buffer[column_index] = byte
        end

        if valid && score >= threshold
            push!(hits, MotifHit(start_index, Int8(1), score, String(copy(window_buffer))))
        end
    end

    return hits
end

function _motif_is_dna_alphabet(alphabet::Vector{Char})
    for character in alphabet
        character in ('A', 'C', 'G', 'T', 'U', 'N') || return false
    end
    return true
end

function _motif_scan_both_strands_dna(sequence::AbstractString, pwm::MotifPWM; threshold::Real=0.0, background=nothing)
    threshold isa Real || throw(ArgumentError("threshold must be real"))
    motif_length = size(pwm.values, 2)
    sequence_bytes = codeunits(sequence)
    sequence_length = length(sequence_bytes)
    sequence_length >= motif_length || return MotifHit[]

    row_lookup = _motif_lookup_array(pwm.alphabet)

    background_vector = _motif_background_vector_from_pwm(pwm.alphabet, background)
    complement = _DNA_COMPLEMENT
    hits = MotifHit[]
    window_buffer = Vector{UInt8}(undef, motif_length)

    @inbounds for start_index in 1:(sequence_length - motif_length + 1)
        forward_score = 0.0
        forward_valid = true
        for column_index in 1:motif_length
            byte = sequence_bytes[start_index + column_index - 1]
            row_index = row_lookup[Int(byte) + 1]
            if row_index == 0
                forward_valid = false
                break
            end

            probability = pwm.values[row_index, column_index]
            probability <= 0 && (forward_valid = false; break)
            background_probability = background_vector[row_index]
            background_probability > 0 || throw(ArgumentError("background entries must be positive"))
            forward_score += log2(probability / background_probability)
            window_buffer[column_index] = byte
        end

        if forward_valid && forward_score >= threshold
            push!(hits, MotifHit(start_index, Int8(1), forward_score, String(copy(window_buffer))))
        end

        reverse_score = 0.0
        reverse_valid = true
        for column_index in 1:motif_length
            byte = sequence_bytes[start_index + motif_length - column_index]
            complemented = complement[Int(byte) + 1]
            row_index = row_lookup[Int(complemented) + 1]
            if row_index == 0
                reverse_valid = false
                break
            end

            probability = pwm.values[row_index, column_index]
            probability <= 0 && (reverse_valid = false; break)
            background_probability = background_vector[row_index]
            background_probability > 0 || throw(ArgumentError("background entries must be positive"))
            reverse_score += log2(probability / background_probability)
            window_buffer[column_index] = complemented
        end

        if reverse_valid && reverse_score >= threshold
            push!(hits, MotifHit(start_index, Int8(-1), reverse_score, String(copy(window_buffer))))
        end
    end

    sort!(hits; by = hit -> (hit.start, -hit.strand, -hit.score))
    return hits
end

function motif_scan_both_strands(sequence::AbstractString, pwm::MotifPWM; threshold::Real=0.0, background=nothing)
    if _motif_is_dna_alphabet(pwm.alphabet)
        return _motif_scan_both_strands_dna(sequence, pwm; threshold=threshold, background=background)
    end

    forward_hits = motif_scan(sequence, pwm; threshold=threshold, background=background)
    motif_length = size(pwm.values, 2)
    sequence_length = ncodeunits(sequence)
    reverse_sequence = reverse_complement(sequence)
    reverse_hits = motif_scan(reverse_sequence, pwm; threshold=threshold, background=background)

    if isempty(reverse_hits)
        return forward_hits
    end

    hits = copy(forward_hits)
    for hit in reverse_hits
        start_index = sequence_length - hit.start - motif_length + 2
        push!(hits, MotifHit(start_index, Int8(-1), hit.score, hit.window))
    end

    sort!(hits; by = hit -> (hit.start, -hit.strand, -hit.score))
    return hits
end

function discover_motifs(
    sequences::AbstractVector{<:AbstractString};
    alphabet::AbstractString="ACGT",
    k::Int=6,
    top_n::Int=3,
    min_support::Int=2,
    max_mismatches::Int=1,
    pseudocount::Real=0.5,
    background=nothing,
    reverse_complements::Bool=true,
)
    isempty(sequences) && throw(ArgumentError("at least one sequence is required"))
    k > 0 || throw(ArgumentError("motif length must be positive"))
    top_n > 0 || throw(ArgumentError("top_n must be positive"))
    min_support > 0 || throw(ArgumentError("min_support must be positive"))
    max_mismatches >= 0 || throw(ArgumentError("max_mismatches must be nonnegative"))

    sequence_lengths = map(ncodeunits, sequences)
    minimum(sequence_lengths) >= k || throw(ArgumentError("motif length must be no greater than each sequence length"))

    seed_counts = _motif_seed_counts(sequences, k; reverse_complements=reverse_complements)
    ranked_seeds = collect(seed_counts)
    sort!(ranked_seeds; by = item -> (-item[2], item[1]))

    selected = MotifDiscoveryResult[]
    for (seed, seed_support) in ranked_seeds
        seed_support < min_support && break
        seed_rc = reverse_complement(seed)
        sites = MotifSite[]

        for (sequence_index, sequence) in enumerate(sequences)
            best_site = _motif_best_site(sequence, seed, seed_rc, k, max_mismatches)
            best_site === nothing && continue
            push!(sites, MotifSite(sequence_index, best_site.start, best_site.strand, best_site.mismatches, best_site.window))
        end

        support = length(sites)
        support < min_support && continue

        site_windows = [site.window for site in sites]
        counts = motif_counts(site_windows; alphabet=alphabet)
        pwm = motif_pwm(counts; pseudocount=pseudocount, background=background)
        information_content = _motif_information_content(pwm; background=background)
        push!(selected, MotifDiscoveryResult(seed, counts.alphabet, counts, pwm, sites, support, information_content))
        length(selected) >= top_n && break
    end

    return selected
end

function discover_motifs(records::AbstractVector{SeqRecordLite}; kwargs...)
    return discover_motifs([record.sequence for record in records]; kwargs...)
end

function _motif_background_vector(alphabet::Vector{Char}, background)
    if background === nothing
        return fill(1.0 / length(alphabet), length(alphabet))
    elseif background isa AbstractVector
        length(background) == length(alphabet) || throw(ArgumentError("background vector length must match the alphabet length"))
        total = sum(Float64.(background))
        total > 0 || throw(ArgumentError("background vector must have positive mass"))
        return Float64.(background) ./ total
    elseif background isa AbstractDict
        values = [Float64(get(background, character, 0.0)) for character in alphabet]
        total = sum(values)
        total > 0 || throw(ArgumentError("background dictionary must have positive mass"))
        return values ./ total
    else
        throw(ArgumentError("background must be nothing, a vector, or a dictionary"))
    end
end

function motif_pwm(
    profile::MotifCounts;
    pseudocount::Real=0.0,
    background=nothing,
)
    pseudocount >= 0 || throw(ArgumentError("pseudocount must be nonnegative"))
    background_vector = _motif_background_vector(profile.alphabet, background)
    nrows  = size(profile.counts, 1)
    ncols  = size(profile.counts, 2)
    values = Matrix{Float64}(undef, nrows, ncols)

    @inbounds for column in 1:ncols
        # Performance: manual sum avoids allocating a column-slice vector
        column_total = 0
        for row in 1:nrows
            column_total += profile.counts[row, column]
        end

        denominator = background === nothing ? column_total + pseudocount * nrows : column_total + pseudocount
        denominator > 0 || throw(ArgumentError("motif column has no mass"))

        for row in 1:nrows
            base_background = background_vector[row]
            base_background > 0 || throw(ArgumentError("background entries must be positive"))
            weight = background === nothing ? pseudocount : pseudocount * base_background
            values[row, column] = (profile.counts[row, column] + weight) / denominator
        end
    end

    return MotifPWM(profile.alphabet, values)
end

function motif_pwm(sequences::AbstractVector{<:AbstractString}; alphabet::AbstractString="ACGT", pseudocount::Real=0.0, background=nothing)
    return motif_pwm(motif_counts(sequences; alphabet=alphabet); pseudocount=pseudocount, background=background)
end

function motif_pwm(records::AbstractVector{SeqRecordLite}; alphabet::AbstractString="ACGT", pseudocount::Real=0.0, background=nothing)
    return motif_pwm(motif_counts(records; alphabet=alphabet); pseudocount=pseudocount, background=background)
end

function motif_consensus(profile::MotifCounts; threshold::Real=0.5)
    0 <= threshold <= 1 || throw(ArgumentError("threshold must be between 0 and 1"))
    consensus = Vector{Char}(undef, size(profile.counts, 2))
    total_sequences = sum(profile.counts[:, 1])

    for column in axes(profile.counts, 2)
        column_counts = profile.counts[:, column]
        best_count, best_index = findmax(column_counts)
        if total_sequences == 0 || best_count / total_sequences < threshold || count(==(best_count), column_counts) != 1
            consensus[column] = 'N'
        else
            consensus[column] = profile.alphabet[best_index]
        end
    end

    return String(consensus)
end

struct MotifOccurrence
    sequence::String
    start::Int
    strand::Int8
    score::Float64
    pvalue::Float64
    site::String
end

struct MotifProfile
    name::String
    alphabet::Vector{Char}
    counts::MotifCounts
    pwm::MotifPWM
    occurrences::Vector{MotifOccurrence}
    metadata::Dict{String,String}
end

function Base.show(io::IO, occurrence::MotifOccurrence)
    print(io, "MotifOccurrence(seq=", occurrence.sequence, ", start=", occurrence.start, ", strand=", occurrence.strand, ", site=", occurrence.site, ")")
end

function Base.show(io::IO, profile::MotifProfile)
    print(io, "MotifProfile(name=", profile.name, ", sites=", length(profile.occurrences), ", width=", size(profile.pwm.values, 2), ")")
end

motif_counts(profile::MotifProfile) = profile.counts
motif_pwm(profile::MotifProfile) = profile.pwm
motif_frequency_matrix(profile::MotifProfile; pseudocount::Real=0.0) = motif_frequency_matrix(profile.counts; pseudocount=pseudocount)
motif_consensus(profile::MotifProfile; threshold::Real=0.5) = motif_consensus(profile.counts; threshold=threshold)
motif_information_content(profile::MotifProfile; background=nothing) = motif_information_content(profile.pwm; background=background)
motif_scan(sequence::AbstractString, profile::MotifProfile; kwargs...) = motif_scan(sequence, profile.pwm; kwargs...)
motif_scan_both_strands(sequence::AbstractString, profile::MotifProfile; kwargs...) = motif_scan_both_strands(sequence, profile.pwm; kwargs...)

function _motif_profile_alphabet(alphabet::AbstractString, alength::Int)
    alphabet_chars = collect(uppercase(alphabet))
    length(alphabet_chars) >= alength || throw(ArgumentError("alphabet is shorter than the motif alphabet length"))
    return alphabet_chars[1:alength]
end

function _motif_parse_float_row(line::AbstractString)
    values = Float64[]
    for match in eachmatch(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
        push!(values, parse(Float64, match.match))
    end
    return values
end

function _motif_parse_int(header::AbstractString, pattern::Regex, default::Int)
    found = Base.match(pattern, header)
    found === nothing && return default
    return parse(Int, found.captures[1])
end

function _motif_parse_float(header::AbstractString, pattern::Regex, default::Float64)
    found = Base.match(pattern, header)
    found === nothing && return default
    return parse(Float64, found.captures[1])
end

function _motif_profile_from_pwm(name::AbstractString, alphabet::Vector{Char}, pwm::MotifPWM; counts=nothing, occurrences=MotifOccurrence[], metadata=Dict{String,String}())
    motif_counts_profile = counts === nothing ? MotifCounts(alphabet, round.(Int, pwm.values .* 1.0)) : counts
    return MotifProfile(String(name), alphabet, motif_counts_profile, pwm, Vector{MotifOccurrence}(occurrences), Dict{String,String}(metadata))
end

function _motif_logo_pwm(profile::Union{MotifCounts, MotifFrequencyMatrix, MotifPWM, MotifProfile})
    profile isa MotifPWM && return profile
    profile isa MotifProfile && return profile.pwm
    profile isa MotifCounts && return motif_pwm(profile; pseudocount=0.5)
    return MotifPWM(profile.alphabet, profile.values)
end

function _motif_column_information(profile::MotifPWM; background=nothing)
    background_vector = _motif_background_vector_from_pwm(profile.alphabet, background)
    ncols = size(profile.values, 2)
    info = zeros(Float64, ncols)

    for column in axes(profile.values, 2)
        column_info = 0.0
        for row in eachindex(profile.alphabet)
            probability = profile.values[row, column]
            probability <= 0 && continue
            background_probability = background_vector[row]
            background_probability > 0 || throw(ArgumentError("background entries must be positive"))
            column_info += probability * log2(probability / background_probability)
        end
        info[column] = column_info
    end

    return info
end

function _motif_color(letter::Char)
    letter == 'A' && return "#2ca02c"
    letter == 'C' && return "#1f77b4"
    letter == 'G' && return "#ff7f0e"
    letter == 'T' && return "#d62728"
    letter == 'U' && return "#d62728"
    return "#7f7f7f"
end

function _motif_profile_label(profile)
    profile isa MotifProfile && return profile.name
    return "Motif Logo"
end

function sequence_logo_svg(profile::Union{MotifCounts, MotifFrequencyMatrix, MotifPWM, MotifProfile}; width::Integer=720, height::Integer=220, background=nothing, title::Union{Nothing,AbstractString}=nothing)
    pwm = _motif_logo_pwm(profile)
    info = _motif_column_information(pwm; background=background)
    ncols = length(info)
    ncols > 0 || return "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$(width)\" height=\"$(height)\"></svg>"

    left_margin = 36.0
    right_margin = 12.0
    top_margin = title === nothing ? 12.0 : 32.0
    bottom_margin = 24.0
    plot_width = max(1.0, Float64(width) - left_margin - right_margin)
    plot_height = max(1.0, Float64(height) - top_margin - bottom_margin)
    column_width = plot_width / ncols
    max_information = max(maximum(info), 1e-9)
    scale = plot_height / max_information

    io = IOBuffer()
    print(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", width, "\" height=\"", height, "\" viewBox=\"0 0 ", width, " ", height, "\">")
    print(io, "<rect x=\"0\" y=\"0\" width=\"", width, "\" height=\"", height, "\" fill=\"white\"/>")
    if title !== nothing
        print(io, "<text x=\"", width / 2, "\" y=\"18\" text-anchor=\"middle\" font-family=\"DejaVu Sans, Arial, sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#111\">", title, "</text>")
    end
    print(io, "<line x1=\"", left_margin, "\" y1=\"", height - bottom_margin, "\" x2=\"", width - right_margin, "\" y2=\"", height - bottom_margin, "\" stroke=\"#444\" stroke-width=\"1\"/>")
    print(io, "<line x1=\"", left_margin, "\" y1=\"", top_margin, "\" x2=\"", left_margin, "\" y2=\"", height - bottom_margin, "\" stroke=\"#444\" stroke-width=\"1\"/>")

    sorted_rows = collect(1:length(pwm.alphabet))
    for column in 1:ncols
        x0 = left_margin + (column - 1) * column_width
        y_base = height - bottom_margin
        column_info = info[column]
        if column_info <= 0
            continue
        end
        sort!(sorted_rows; by = row -> pwm.values[row, column])
        for row in sorted_rows
            probability = pwm.values[row, column]
            probability <= 0 && continue
            block_height = probability * column_info * scale
            block_height <= 0 && continue
            y_base -= block_height
            letter = pwm.alphabet[row]
            print(io, "<rect x=\"", x0 + 1.0, "\" y=\"", y_base, "\" width=\"", max(1.0, column_width - 2.0), "\" height=\"", block_height, "\" rx=\"2\" ry=\"2\" fill=\"", _motif_color(letter), "\" opacity=\"0.92\"/>")
            font_size = min(24.0, max(8.0, block_height * 0.72))
            print(io, "<text x=\"", x0 + column_width / 2, "\" y=\"", y_base + block_height * 0.72, "\" text-anchor=\"middle\" font-family=\"DejaVu Sans, Arial, sans-serif\" font-size=\"", font_size, "\" font-weight=\"700\" fill=\"white\">", letter, "</text>")
        end
        print(io, "<text x=\"", x0 + column_width / 2, "\" y=\"", height - 6, "\" text-anchor=\"middle\" font-family=\"DejaVu Sans, Arial, sans-serif\" font-size=\"10\" fill=\"#444\">", column, "</text>")
    end

    print(io, "</svg>")
    return String(take!(io))
end

motif_logo_svg(profile::Union{MotifCounts, MotifFrequencyMatrix, MotifPWM, MotifProfile}; kwargs...) = sequence_logo_svg(profile; kwargs...)
sequence_logo(profile::Union{MotifCounts, MotifFrequencyMatrix, MotifPWM, MotifProfile}; kwargs...) = sequence_logo_svg(profile; kwargs...)

function _motif_finalize_profile(name::AbstractString, alphabet::Vector{Char}, counts_matrix::Matrix{Int}, pwm_values::Matrix{Float64}, occurrences::Vector{MotifOccurrence}, metadata::Dict{String,String})
    counts = MotifCounts(alphabet, counts_matrix)
    pwm = MotifPWM(alphabet, pwm_values)
    return MotifProfile(String(name), alphabet, counts, pwm, occurrences, metadata)
end

function read_meme(filepath::AbstractString; alphabet::AbstractString="ACGT")
    open(filepath, "r") do io
        return read_meme(io; alphabet=alphabet)
    end
end

function read_meme(io::IO; alphabet::AbstractString="ACGT")
    lines = collect(eachline(io))
    profiles = MotifProfile[]
    current_name = ""
    i = 1

    while i <= length(lines)
        line = strip(lines[i])
        if startswith(line, "MOTIF")
            parts = Base.split(line)
            current_name = length(parts) > 1 ? join(parts[2:end], " ") : "motif$(length(profiles) + 1)"
            i += 1
            continue
        end

        if occursin("letter-probability matrix:", line)
            alength = _motif_parse_int(line, r"alength\s*=\s*(\d+)", length(collect(uppercase(alphabet))))
            width = _motif_parse_int(line, r"w\s*=\s*(\d+)", 0)
            nsites = _motif_parse_float(line, r"nsites\s*=\s*([0-9eE+\-.]+)", Float64(width == 0 ? 1 : width))
            motif_alphabet = _motif_profile_alphabet(alphabet, alength)
            width > 0 || throw(ArgumentError("MEME matrix width could not be determined"))

            matrix = zeros(Float64, alength, width)
            for column in 1:width
                i += 1
                i <= length(lines) || throw(ArgumentError("unexpected end of MEME motif matrix"))
                row_values = _motif_parse_float_row(lines[i])
                length(row_values) >= alength || throw(ArgumentError("MEME matrix row has too few probabilities"))
                @inbounds for row in 1:alength
                    matrix[row, column] = row_values[row]
                end
            end

            counts_matrix = round.(Int, matrix .* max(nsites, 1.0))
            metadata = Dict{String,String}("source" => "MEME", "nsites" => string(nsites), "width" => string(width))
            occurrences = MotifOccurrence[]

            j = i + 1
            while j <= length(lines)
                site_line = strip(lines[j])
                isempty(site_line) && break
                startswith(site_line, "MOTIF") && break
                startswith(site_line, "letter-probability matrix:") && break
                tokens = Base.split(site_line)
                if length(tokens) >= 4 && tokens[1] != "Motif" && tryparse(Int, tokens[2]) !== nothing
                    start = parse(Int, tokens[2])
                    pvalue = tryparse(Float64, tokens[3]) === nothing ? 0.0 : parse(Float64, tokens[3])
                    site = uppercase(tokens[end])
                    if all(character -> character in ('A', 'C', 'G', 'T', 'U', 'N', '-') , site)
                        strand = any(token -> token == "-", tokens) ? Int8(-1) : Int8(1)
                        push!(occurrences, MotifOccurrence(tokens[1], start, strand, 0.0, pvalue, site))
                    end
                end
                j += 1
            end
            push!(profiles, _motif_finalize_profile(current_name == "" ? "MEME_$(length(profiles) + 1)" : current_name, motif_alphabet, counts_matrix, matrix, occurrences, metadata))
            current_name = ""
            i = j
            continue
        end

        i += 1
    end

    return profiles
end

function read_alignace(filepath::AbstractString; alphabet::AbstractString="ACGT")
    open(filepath, "r") do io
        return read_alignace(io; alphabet=alphabet)
    end
end

function read_alignace(io::IO; alphabet::AbstractString="ACGT")
    motif_profiles = MotifProfile[]
    current_name = Ref{String}("")
    current_sites = String[]
    current_occurrences = MotifOccurrence[]
    current_metadata = Dict{String,String}()

    function finalize_current()
        isempty(current_sites) && isempty(current_occurrences) && return nothing
        motif_name = isempty(current_name[]) ? "ALIGNACE_$(length(motif_profiles) + 1)" : current_name[]
        counts = motif_counts(current_sites; alphabet=alphabet)
        pwm = motif_pwm(counts; pseudocount=0.5)
        metadata = Dict{String,String}(current_metadata)
        metadata["source"] = "ALIGNACE"
        push!(motif_profiles, MotifProfile(motif_name, counts.alphabet, counts, pwm, Vector{MotifOccurrence}(current_occurrences), metadata))
        empty!(current_sites)
        empty!(current_occurrences)
        empty!(current_metadata)
        current_name[] = ""
        return nothing
    end

    for raw_line in eachline(io)
        line = strip(raw_line)
        isempty(line) && (finalize_current(); continue)
        startswith(line, "#") && continue

        upper = uppercase(line)
        if startswith(upper, "ALIGNACE")
            current_metadata["header"] = line
            continue
        end
        if startswith(upper, "MOTIF")
            finalize_current()
            parts = Base.split(line)
            current_name[] = length(parts) > 1 ? join(parts[2:end], " ") : ""
            continue
        end
        if occursin("CONSENSUS", upper)
            parts = Base.split(line, ['=', ':']; limit=2)
            if length(parts) == 2
                current_metadata["consensus"] = strip(parts[2])
            end
            continue
        end

        tokens = Base.split(line)
        if isempty(tokens)
            continue
        end

        site_candidate = uppercase(tokens[end])
        if all(character -> character in ('A', 'C', 'G', 'T', 'U', 'N', '-') , site_candidate)
            push!(current_sites, site_candidate)
            sequence_name = length(tokens) > 1 ? tokens[1] : "site$(length(current_sites))"
            start = length(tokens) > 2 && tryparse(Int, tokens[2]) !== nothing ? parse(Int, tokens[2]) : length(current_sites)
            strand = occursin('-', line) ? Int8(-1) : Int8(1)
            pvalue = length(tokens) > 3 && tryparse(Float64, tokens[3]) !== nothing ? parse(Float64, tokens[3]) : 0.0
            push!(current_occurrences, MotifOccurrence(sequence_name, start, strand, 0.0, pvalue, site_candidate))
        end
    end

    finalize_current()
    return motif_profiles
end

function _jaspar_parse_row(line::AbstractString)
    values = Float64[]
    for match in eachmatch(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", line)
        push!(values, parse(Float64, match.match))
    end
    return values
end

function read_jaspar(filepath::AbstractString; alphabet::AbstractString="ACGT")
    open(filepath, "r") do io
        return read_jaspar(io; alphabet=alphabet)
    end
end

function read_jaspar(io::IO; alphabet::AbstractString="ACGT")
    profiles = MotifProfile[]
    lines = collect(eachline(io))
    i = 1

    while i <= length(lines)
        line = strip(lines[i])
        isempty(line) && (i += 1; continue)
        startswith(line, "#") && (i += 1; continue)

        if startswith(line, ">")
            header = Base.split(String(line[2:end]))
            jaspar_id = isempty(header) ? "" : header[1]
            jaspar_name = length(header) > 1 ? join(header[2:end], " ") : jaspar_id
            i += 1

            row_map = Dict{Char,Vector{Float64}}()
            while i <= length(lines)
                row_line = strip(lines[i])
                isempty(row_line) && break
                startswith(row_line, ">") && break
                row_char = uppercase(first(row_line))
                if row_char in uppercase(alphabet)
                    row_map[row_char] = _jaspar_parse_row(row_line)
                end
                i += 1
            end

            isempty(row_map) && continue

            motif_alphabet = collect(uppercase(alphabet))
            width = maximum(length(values) for values in values(row_map))
            counts_matrix = zeros(Int, length(motif_alphabet), width)

            for (row_index, character) in pairs(motif_alphabet)
                row_values = get(row_map, character, Float64[])
                @inbounds for column in eachindex(row_values)
                    counts_matrix[row_index, column] = round(Int, row_values[column])
                end
            end

            counts = MotifCounts(motif_alphabet, counts_matrix)
            pwm = motif_pwm(counts; pseudocount=0.0)
            metadata = Dict{String,String}("source" => "JASPAR")
            !isempty(jaspar_id) && (metadata["jaspar_id"] = jaspar_id)
            !isempty(jaspar_name) && (metadata["name"] = jaspar_name)
            push!(profiles, MotifProfile(!isempty(jaspar_name) ? jaspar_name : jaspar_id, motif_alphabet, counts, pwm, MotifOccurrence[], metadata))
            continue
        end

        i += 1
    end

    return profiles
end
