# ==============================================================================
# motif.jl — Sequence motif discovery and scanning
#
# Provides position weight matrix (PWM) construction, JASPAR/TRANSFAC I/O,
# motif scanning with IUPAC consensus, information content computation,
# and de novo motif discovery via expectation-maximization.
#
# References:
#   - Stormo (2000) Bioinformatics 16(1):16-23 (weight matrices)
#   - Sandelin et al. (2004) NAR 32(D91-D94) (JASPAR database)
# ==============================================================================

using DataFrames
using SHA
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, register_provenance!

"""
    MotifCounts{A <: BioAlphabet}

Count matrix for symbols in alphabet `A` across a motif. 
The `alphabet` field stores the symbol order (rows of the `counts` matrix).
"""
struct MotifCounts{A <: BioAlphabet}
    alphabet::Vector{UInt8}
    counts::Matrix{Int}
end

"""
    MotifFrequencyMatrix{A <: BioAlphabet}

Probability matrix (PFM) for symbols in alphabet `A`.
"""
struct MotifFrequencyMatrix{A <: BioAlphabet}
    alphabet::Vector{UInt8}
    values::Matrix{Float64}
end

"""
    MotifPWM{A <: BioAlphabet}

Position weight matrix (PWM) for symbols in alphabet `A`.
Scores are usually in log2-odds format.
"""
struct MotifPWM{A <: BioAlphabet}
    alphabet::Vector{UInt8}
    values::Matrix{Float64}
end

function MotifCounts(alphabet::AbstractVector{<:Union{Char,UInt8}}, counts::AbstractMatrix{<:Integer})
    alphabet_bytes = Vector{UInt8}(codeunits(String(collect(alphabet))))
    return MotifCounts{_motif_infer_alphabet(String(alphabet_bytes))}(alphabet_bytes, Matrix{Int}(counts))
end

function MotifFrequencyMatrix(alphabet::AbstractVector{<:Union{Char,UInt8}}, values::AbstractMatrix{<:Real})
    alphabet_bytes = Vector{UInt8}(codeunits(String(collect(alphabet))))
    return MotifFrequencyMatrix{_motif_infer_alphabet(String(alphabet_bytes))}(alphabet_bytes, Matrix{Float64}(values))
end

function MotifPWM(alphabet::AbstractVector{<:Union{Char,UInt8}}, values::AbstractMatrix{<:Real})
    alphabet_bytes = Vector{UInt8}(codeunits(String(collect(alphabet))))
    return MotifPWM{_motif_infer_alphabet(String(alphabet_bytes))}(alphabet_bytes, Matrix{Float64}(values))
end

struct MotifHit{A <: BioAlphabet}
    start::Int
    strand::Int8
    score::Float64
    window::BioSequence{A}
end

struct MotifSite{A <: BioAlphabet}
    sequence_index::Int
    start::Int
    strand::Int8
    mismatches::Int
    window::BioSequence{A}
end

"""
    MotifDiscoveryResult{A <: BioAlphabet}

Result of de novo motif discovery, including the PWM and identified sites.
"""
struct MotifDiscoveryResult{A <: BioAlphabet} <: AbstractAnalysisResult
    seed::BioSequence{A}
    alphabet::Vector{UInt8}
    counts::MotifCounts{A}
    pwm::MotifPWM{A}
    sites::Vector{MotifSite{A}}
    support::Int
    information_content::Float64
end

function _motif_infer_alphabet(sequence)
    upper = uppercase(String(sequence))
    if all(character -> character in ('A', 'C', 'G', 'T', 'U', 'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', '-'), upper)
        return DNAAlphabet
    end
    return AminoAcidAlphabet
end

MotifHit(start::Integer, strand::Integer, score::Real, window::BioSequence{A}) where {A <: BioAlphabet} = MotifHit{A}(Int(start), Int8(strand), Float64(score), window)

function MotifHit(start::Integer, strand::Integer, score::Real, window)
    window_text = String(window)
    alphabet = _motif_infer_alphabet(window_text)
    return MotifHit(Int(start), Int8(strand), Float64(score), BioSequence{alphabet}(window_text))
end

MotifSite(sequence_index::Integer, start::Integer, strand::Integer, mismatches::Integer, window::BioSequence{A}) where {A <: BioAlphabet} = MotifSite{A}(Int(sequence_index), Int(start), Int8(strand), Int(mismatches), window)

MotifDiscoveryResult(seed::BioSequence{A}, alphabet::Vector{UInt8}, counts::MotifCounts{A}, pwm::MotifPWM{A}, sites::Vector{MotifSite{A}}, support::Integer, information_content::Real) where {A <: BioAlphabet} = MotifDiscoveryResult{A}(seed, alphabet, counts, pwm, sites, Int(support), Float64(information_content))

# Manual byte mappings are replaced by the parametric symbol system
# where possible, keeping compat for IUPAC letters.

function Base.show(io::IO, profile::MotifCounts{A}) where {A}
    print(io, "MotifCounts{", A, "}(", String(profile.alphabet), ", ", size(profile.counts, 2), " positions)")
end

function Base.show(io::IO, profile::MotifFrequencyMatrix{A}) where {A}
    print(io, "MotifFrequencyMatrix{", A, "}(", String(profile.alphabet), ", ", size(profile.values, 2), " positions)")
end

function Base.show(io::IO, pwm::MotifPWM{A}) where {A}
    print(io, "MotifPWM{", A, "}(", String(pwm.alphabet), ", ", size(pwm.values, 2), " positions)")
end

function Base.show(io::IO, hit::MotifHit)
    print(io, "MotifHit(start=", hit.start, ", strand=", hit.strand, ", score=", round(hit.score, digits=3), ")")
end

function Base.show(io::IO, site::MotifSite)
    print(io, "MotifSite(seq=", site.sequence_index, ", start=", site.start, ", strand=", site.strand, ", mismatches=", site.mismatches, ")")
end

function Base.show(io::IO, result::MotifDiscoveryResult)
    print(io, analysis_result_summary(result))
end

function DataFrames.DataFrame(hits::AbstractVector{<:MotifHit})
    rows = collect(hits)
    return DataFrames.DataFrame(
        start = [hit.start for hit in rows],
        strand = [hit.strand for hit in rows],
        score = [hit.score for hit in rows],
        window = [String(hit.window) for hit in rows])
end

function DataFrames.DataFrame(sites::AbstractVector{<:MotifSite})
    rows = collect(sites)
    return DataFrames.DataFrame(
        sequence_index = [site.sequence_index for site in rows],
        start = [site.start for site in rows],
        strand = [site.strand for site in rows],
        mismatches = [site.mismatches for site in rows],
        window = [String(site.window) for site in rows])
end

function _motif_alphabet_index(A::Type{<:BioAlphabet}, alphabet_str::String)
    alphabet_bytes = collect(codeunits(uppercase(alphabet_str)))
    index = Dict{UInt8,Int}()

    valid = symbols(A)
    for (position, byte) in enumerate(alphabet_bytes)
        byte in valid || throw(ArgumentError("symbol '$(Char(byte))' not in alphabet $(A)"))
        index[byte] = position
    end

    return alphabet_bytes, index
end

"""
    motif_counts(sequences; alphabet="ACGT")

Count symbol frequencies at each position across a set of aligned sequences.
The primary API accepts `AbstractVector{BioSequence{A}}`.
"""
function motif_counts(sequences::AbstractVector{BioSequence{A}}; alphabet::String="ACGT", prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx)) where {A <: BioAlphabet}
    isempty(sequences) && throw(ArgumentError("at least one sequence is required"))

    alphabet_bytes, alphabet_index = _motif_alphabet_index(A, alphabet)
    motif_length = length(first(sequences))
    counts = zeros(Int, length(alphabet_bytes), motif_length)
    
    # Pre-calculate mapping for speed
    lookup = fill(0, 256)
    for (byte, pos) in alphabet_index
        lookup[Int(byte) + 1] = pos
    end

    for sequence in sequences
        length(sequence) == motif_length || throw(ArgumentError("all motif sequences must have the same length"))
        bytes = sequence.data
        @inbounds for (position, byte) in enumerate(bytes)
            row = lookup[Int(byte) + 1]
            row == 0 && throw(ArgumentError("unsupported motif character $(Char(byte)) for alphabet $(A)"))
            counts[row, position] += 1
        end
    end

    return MotifCounts{A}(alphabet_bytes, counts)
end

# Removed motif_counts(::AbstractVector{<:AbstractString}) - use Vector{BioSequence} instead
function motif_counts(sequences::AbstractVector{<:AbstractString}; alphabet::String="ACGT")
    return motif_counts(DNASeq.(String.(sequences)); alphabet=alphabet)
end

function motif_counts(records::AbstractVector{<:SeqRecordLite}; alphabet::String="ACGT")
    return motif_counts([record.sequence for record in records]; alphabet=alphabet)
end

function _motif_frequency_matrix(profile::MotifCounts{A}; pseudocount::Real=0.0) where {A}
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

    return MotifFrequencyMatrix{A}(profile.alphabet, values)
end

function motif_frequency_matrix(profile::MotifCounts{A}; pseudocount::Real=0.0) where {A}
    return _motif_frequency_matrix(profile; pseudocount=pseudocount)
end

function motif_frequency_matrix(records::AbstractVector{<:SeqRecordLite}; alphabet::String="ACGT", pseudocount::Real=0.0)
    return motif_frequency_matrix(motif_counts(records; alphabet=alphabet); pseudocount=pseudocount)
end

function motif_entropy(profile::MotifFrequencyMatrix{A}) where {A}
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

function motif_relative_entropy(profile::MotifFrequencyMatrix{A}; background=nothing) where {A}
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

@inline function _motif_hamming_distance(left::BioSequence, right::BioSequence)
    length(left) == length(right) || throw(ArgumentError("motif windows must have the same length"))
    mismatches = 0
    @inbounds for (left_byte, right_byte) in zip(left.data, right.data)
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

function _motif_canonical_kmer(kmer::String; reverse_complements::Bool=true)
    reverse_complements || return String(kmer)
    reverse_kmer = reverse_complement(kmer)
    return reverse_kmer < kmer ? reverse_kmer : String(kmer)
end

function _motif_background_vector_from_pwm(alphabet::Vector{Char}, background)
    return _motif_background_vector(alphabet, background)
end

function _motif_background_vector_from_pwm(alphabet::Vector{UInt8}, background)
    return _motif_background_vector(alphabet, background)
end

function _motif_seed_counts(sequences::AbstractVector{<:String}, k::Int; reverse_complements::Bool=true)
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

function _motif_best_site(sequence::BioSequence{A}, seed::BioSequence{A}, seed_rc::BioSequence{A}, k::Int, max_mismatches::Int) where {A <: BioAlphabet}
    sequence_length = length(sequence)
    best_site = nothing
    best_mismatches = max_mismatches + 1
    best_strand = Int8(1)

    for start_index in 1:(sequence_length - k + 1)
        window = sequence[start_index:start_index + k - 1]
        mismatches = _motif_hamming_distance(window, seed)
        strand = Int8(1)
        oriented_window = window

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

function _motif_information_content(profile::MotifPWM{A}; background=nothing) where {A}
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

function motif_information_content(profile::MotifPWM{A}; background=nothing) where {A}
    return _motif_information_content(profile; background=background)
end

"""
    motif_scan(sequence, pwm; threshold=0.0, background=nothing)

Scan a biological sequence for occurrences of a motif PWM.
Primary API accepts `BioSequence{A}`.
"""
function motif_scan(sequence::BioSequence{A}, pwm::MotifPWM{A}; threshold::Real=0.0, background=nothing, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx)) where {A <: BioAlphabet}
    threshold isa Real || throw(ArgumentError("threshold must be real"))
    motif_length = size(pwm.values, 2)
    sequence_bytes = sequence.data
    sequence_length = length(sequence_bytes)
    sequence_length >= motif_length || return MotifHit{A}[]

    # Pre-build lookup for the PWM's symbol order
    row_lookup = fill(0, 256)
    for (i, byte) in enumerate(pwm.alphabet)
        row_lookup[Int(byte) + 1] = i
    end

    background_vector = _motif_background_vector_from_pwm(pwm.alphabet, background)
    hits = MotifHit{A}[]
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
            push!(hits, MotifHit(start_index, Int8(1), score, BioSequence{A}(copy(window_buffer); validate=false)))
        end
    end

    if _ctx !== nothing
        register_provenance!(_ctx, "motif_scan";
        parents=provenance_parent_ids(sequence, pwm),
        parameters=(threshold=threshold, hit_count=length(hits), alphabet=string(A)))
    end
    return hits
end

# Removed motif_scan(::AbstractString, ::MotifPWM) - use BioSequence instead

"""
    motif_scan_both_strands(sequence, pwm; threshold=0.0, background=nothing)

Scan both strands of a DNA/RNA sequence for a motif.
"""
function motif_scan_both_strands(sequence::BioSequence{A}, pwm::MotifPWM{A}; threshold::Real=0.0, background=nothing, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx)) where {A <: Union{DNAAlphabet, RNAAlphabet}}
    hits = _motif_scan_both_strands_nucleotide(sequence, pwm; threshold=threshold, background=background)
    _ctx !== nothing && register_provenance!(_ctx, "motif_scan_both_strands"; parents=provenance_parent_ids(sequence, pwm), parameters=(threshold=threshold, hit_count=length(hits), alphabet=string(A)))
    return hits
end

# Removed motif_scan_both_strands(::AbstractString, ::MotifPWM) - use BioSequence instead
function motif_scan_both_strands(sequence::AbstractString, pwm::MotifPWM{A}; threshold::Real=0.0, background=nothing) where {A <: Union{DNAAlphabet, RNAAlphabet}}
    typed_sequence = A === DNAAlphabet ? DNASeq(sequence; validate=false) : RNASeq(sequence; validate=false)
    _ctx = active_provenance_context()


    return motif_scan_both_strands(typed_sequence, pwm; threshold=threshold, background=background, _ctx=_ctx)
end

function _motif_scan_both_strands_nucleotide(sequence::BioSequence{A}, pwm::MotifPWM{A}; threshold::Real=0.0, background=nothing) where {A}
    threshold isa Real || throw(ArgumentError("threshold must be real"))
    motif_length = size(pwm.values, 2)
    sequence_bytes = sequence.data
    sequence_length = length(sequence_bytes)
    sequence_length >= motif_length || return MotifHit{A}[]

    # Pre-build lookup for the PWM's symbol order
    row_lookup = fill(0, 256)
    for (i, byte) in enumerate(pwm.alphabet)
        row_lookup[Int(byte) + 1] = i
    end

    background_vector = _motif_background_vector_from_pwm(pwm.alphabet, background)
    complement = _DNA_COMPLEMENT
    hits = MotifHit{A}[]
    window_buffer = Vector{UInt8}(undef, motif_length)

    @inbounds for start_index in 1:(sequence_length - motif_length + 1)
        # Forward strand
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
            push!(hits, MotifHit(start_index, Int8(1), forward_score, BioSequence{A}(copy(window_buffer); validate=false)))
        end

        # Reverse strand
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
            push!(hits, MotifHit(start_index, Int8(-1), reverse_score, BioSequence{A}(copy(window_buffer); validate=false)))
        end
    end

    sort!(hits; by = hit -> (hit.start, -hit.strand, -hit.score))
    return hits
end

"""
    discover_motifs(sequences; kwargs...)

Perform de novo motif discovery using a seed-and-extend approach.
"""
function discover_motifs(
    sequences::AbstractVector{BioSequence{A}};
    alphabet::String=A <: DNAAlphabet ? "ACGT" : "ACDEFGHIKLMNPQRSTVWY",
    k::Int=6,
    top_n::Int=3,
    min_support::Int=2,
    max_mismatches::Int=1,
    pseudocount::Real=0.5,
    background=nothing,
    reverse_complements::Bool=true) where {A <: BioAlphabet}
    isempty(sequences) && throw(ArgumentError("at least one sequence is required"))
    k > 0 || throw(ArgumentError("motif length must be positive"))
    top_n > 0 || throw(ArgumentError("top_n must be positive"))
    min_support > 0 || throw(ArgumentError("min_support must be positive"))
    max_mismatches >= 0 || throw(ArgumentError("max_mismatches must be nonnegative"))

    sequence_lengths = map(length, sequences)
    minimum(sequence_lengths) >= k || throw(ArgumentError("motif length must be no greater than each sequence length"))

    # Convert to String briefly for seed counting (dict keys)
    str_seqs = [String(s) for s in sequences]
    seed_counts = _motif_seed_counts(str_seqs, k; reverse_complements=reverse_complements)
    ranked_seeds = collect(seed_counts)
    sort!(ranked_seeds; by = item -> (-item[2], item[1]))

    selected = MotifDiscoveryResult{A}[]
    typed_seed_cache = Dict{String,BioSequence{A}}()
    for (seed, seed_support) in ranked_seeds
        seed_support < min_support && break
        seed_seq = get!(typed_seed_cache, seed) do
            BioSequence{A}(seed)
        end
        seed_rc = reverse_complement(seed_seq)
        sites = MotifSite{A}[]

        for (sequence_index, sequence) in enumerate(sequences)
            best_site = _motif_best_site(sequence, seed_seq, seed_rc, k, max_mismatches)
            best_site === nothing && continue
            push!(sites, MotifSite(sequence_index, best_site.start, best_site.strand, best_site.mismatches, best_site.window))
        end

        support = length(sites)
        support < min_support && continue

        site_seqs = [site.window for site in sites]
        counts = motif_counts(site_seqs; alphabet=alphabet)
        pwm = motif_pwm(counts; pseudocount=pseudocount, background=background)
        information_content = _motif_information_content(pwm; background=background)
        push!(selected, MotifDiscoveryResult(seed_seq, counts.alphabet, counts, pwm, sites, support, information_content))
        length(selected) >= top_n && break
    end

    _ctx = active_provenance_context()
    _ctx !== nothing && register_provenance!(_ctx, "discover_motifs"; parents=provenance_parent_ids(sequences), parameters=(alphabet=alphabet, k=k, top_n=top_n, min_support=min_support, max_mismatches=max_mismatches, motif_count=length(selected)))
    return selected
end

function discover_motifs(records::AbstractVector{<:SeqRecordLite}, kwargs...)
    _ctx = active_provenance_context()


    return discover_motifs([record.sequence for record in records]; _ctx=_ctx, kwargs...)
end

function _motif_background_vector(alphabet::Vector{UInt8}, background)
    if background === nothing
        return fill(1.0 / length(alphabet), length(alphabet))
    elseif background isa AbstractVector
        length(background) == length(alphabet) || throw(ArgumentError("background vector length must match the alphabet length"))
        total = sum(Float64.(background))
        total > 0 || throw(ArgumentError("background vector must have positive mass"))
        return Float64.(background) ./ total
    elseif background isa AbstractDict
        # Support both Char and UInt8 keys
        values = Float64[]
        for byte in alphabet
            val = get(background, byte, get(background, Char(byte), 0.0))
            push!(values, Float64(val))
        end
        total = sum(values)
        total > 0 || throw(ArgumentError("background dictionary must have positive mass"))
        return values ./ total
    else
        throw(ArgumentError("background must be nothing, a vector, or a dictionary"))
    end
end

function motif_pwm(
    profile::MotifCounts{A};
    pseudocount::Real=0.0,
    background=nothing,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx)) where {A}
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

    pwm = MotifPWM{A}(profile.alphabet, values)
    _ctx !== nothing && register_provenance!(_ctx, "motif_pwm"; parents=provenance_parent_ids(profile), parameters=(pseudocount=pseudocount, background=background === nothing ? "none" : "provided"))
    return pwm
end

function motif_pwm(sequences::AbstractVector{BioSequence{A}}; alphabet::String="ACGT", pseudocount::Real=0.0, background=nothing) where {A}
    _ctx = active_provenance_context()


    return motif_pwm(motif_counts(sequences; alphabet=alphabet); pseudocount=pseudocount, background=background, _ctx=_ctx)
end

function motif_pwm(records::AbstractVector{<:SeqRecordLite}; alphabet::String="ACGT", pseudocount::Real=0.0, background=nothing)
    _ctx = active_provenance_context()


    return motif_pwm(motif_counts(records; alphabet=alphabet); pseudocount=pseudocount, background=background, _ctx=_ctx)
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
            consensus[column] = Char(profile.alphabet[best_index])
        end
    end

    result = String(consensus)
    _ctx = active_provenance_context()
    _ctx !== nothing && register_provenance!(_ctx, "motif_consensus"; parents=provenance_parent_ids(profile), parameters=(threshold=threshold, consensus=result))
    return result
end

struct MotifOccurrence{A <: BioAlphabet}
    sequence_id::String
    start::Int
    strand::Int8
    score::Float64
    pvalue::Float64
    site::BioSequence{A}
end

MotifOccurrence(sequence_id::AbstractString, start::Integer, strand::Integer, score::Real, pvalue::Real, site::BioSequence{A}) where {A <: BioAlphabet} =
    MotifOccurrence{A}(String(sequence_id), Int(start), Int8(strand), Float64(score), Float64(pvalue), site)

struct MotifProfile
    name::String
    alphabet::Vector{Char}
    counts::MotifCounts
    pwm::MotifPWM
    occurrences::Vector{MotifOccurrence}
    metadata::Dict{String,String}
end

function Base.show(io::IO, occurrence::MotifOccurrence)
    print(io, "MotifOccurrence(seq=", occurrence.sequence_id, ", start=", occurrence.start, ", strand=", occurrence.strand, ", site=", String(occurrence.site), ")")
end

function Base.show(io::IO, profile::MotifProfile)
    print(io, "MotifProfile(name=", profile.name, ", sites=", length(profile.occurrences), ", width=", size(profile.pwm.values, 2), ")")
end

motif_counts(profile::MotifProfile) = profile.counts
motif_pwm(profile::MotifProfile) = profile.pwm
motif_frequency_matrix(profile::MotifProfile; pseudocount::Real=0.0) = motif_frequency_matrix(profile.counts; pseudocount=pseudocount)
motif_consensus(profile::MotifProfile; threshold::Real=0.5) = motif_consensus(profile.counts; threshold=threshold)
motif_information_content(profile::MotifProfile; background=nothing) = motif_information_content(profile.pwm; background=background)
motif_scan(sequence::BioSequence, profile::MotifProfile; kwargs...) = motif_scan(sequence, profile.pwm; kwargs...)
motif_scan_both_strands(sequence::BioSequence, profile::MotifProfile; kwargs...) = motif_scan_both_strands(sequence, profile.pwm; kwargs...)
# Removed motif_scan(::AbstractString, ::MotifProfile) - use BioSequence instead
# Removed motif_scan_both_strands(::AbstractString, ::MotifProfile) - use BioSequence instead

function _motif_profile_alphabet(alphabet::String, alength::Int)
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

function _motif_profile_from_pwm(name::String, alphabet::Vector{Char}, pwm::MotifPWM; counts=nothing, occurrences=MotifOccurrence[], metadata=Dict{String,String}())
    motif_counts_profile = counts === nothing ? MotifCounts{_motif_infer_alphabet(String(alphabet))}(Vector{UInt8}(codeunits(String(alphabet))), round.(Int, pwm.values .* 1.0)) : counts
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

_motif_color(letter::UInt8) = _motif_color(Char(letter))

function _motif_profile_label(profile)
    profile isa MotifProfile && return profile.name
    return "Motif Logo"
end

function sequence_logo_svg(profile::Union{MotifCounts, MotifFrequencyMatrix, MotifPWM, MotifProfile}; width::Integer=720, height::Integer=220, background=nothing, title::Union{Nothing,String}=nothing)
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
            letter = Char(pwm.alphabet[row])
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

function _motif_finalize_profile(name::String, alphabet::AbstractVector{<:Union{Char,UInt8}}, counts_matrix::Matrix{Int}, pwm_values::Matrix{Float64}, occurrences::Vector{MotifOccurrence}, metadata::Dict{String,String}; provenance_id::Union{Nothing,String}=nothing, provenance_hash::Union{Nothing,String}=nothing)
    alphabet_chars = alphabet isa Vector{Char} ? alphabet : Char.(alphabet)
    counts = MotifCounts{_motif_infer_alphabet(String(alphabet_chars))}(Vector{UInt8}(codeunits(String(alphabet_chars))), counts_matrix)
    pwm = MotifPWM(alphabet_chars, pwm_values)
    provenance_id === nothing || (metadata["provenance_id"] = provenance_id)
    provenance_hash === nothing || (metadata["provenance_hash"] = provenance_hash)
    return MotifProfile(String(name), alphabet_chars, counts, pwm, occurrences, metadata)
end

function read_meme(filepath::String; alphabet::String="ACGT", prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    if _ctx === nothing
        open(filepath, "r") do io
            return read_meme(io; alphabet=alphabet)
        end
    end
    raw_bytes = read(filepath)
    provenance_hash = bytes2hex(sha256(raw_bytes))


    return read_meme(IOBuffer(raw_bytes); alphabet=alphabet, _ctx=_ctx, provenance_hash=provenance_hash, provenance_source=filepath)
end

function read_meme(io::IO; alphabet::String="ACGT", provenance_hash::Union{Nothing,AbstractString}=nothing, provenance_source::AbstractString="read_meme", prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
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
                        push!(occurrences, MotifOccurrence(tokens[1], start, strand, 0.0, pvalue, BioSequence{DNAAlphabet}(site)))
                    end
                end
                j += 1
            end
            profile_name = current_name == "" ? "MEME_$(length(profiles) + 1)" : current_name
            profile_id = _ctx === nothing ? nothing : register_provenance!(_ctx, "read_meme"; parents=String[], parameters=(source=provenance_source, motif_name=profile_name, hash=provenance_hash)).id
            push!(profiles, _motif_finalize_profile(profile_name, motif_alphabet, counts_matrix, matrix, occurrences, metadata; provenance_id=profile_id, provenance_hash=provenance_hash === nothing ? nothing : String(provenance_hash)))
            current_name = ""
            i = j
            continue
        end

        i += 1
    end

    _ctx !== nothing && register_provenance!(_ctx, "read_meme"; parents=String[], parameters=(source=provenance_source, profile_count=length(profiles), hash=provenance_hash))
    return profiles
end

function read_alignace(filepath::String; alphabet::String="ACGT", prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    if _ctx === nothing
        open(filepath, "r") do io
            return read_alignace(io; alphabet=alphabet)
        end
    end
    raw_bytes = read(filepath)
    provenance_hash = bytes2hex(sha256(raw_bytes))


    return read_alignace(IOBuffer(raw_bytes); alphabet=alphabet, _ctx=_ctx, provenance_hash=provenance_hash, provenance_source=filepath)
end

function read_alignace(io::IO; alphabet::String="ACGT", provenance_hash::Union{Nothing,AbstractString}=nothing, provenance_source::AbstractString="read_alignace", prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
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
        profile_id = _ctx === nothing ? nothing : register_provenance!(_ctx, "read_alignace"; parents=String[], parameters=(source=provenance_source, motif_name=motif_name, hash=provenance_hash)).id
        push!(motif_profiles, _motif_finalize_profile(motif_name, counts.alphabet, Matrix{Int}(counts.counts), Matrix{Float64}(pwm.values), Vector{MotifOccurrence}(current_occurrences), metadata; provenance_id=profile_id, provenance_hash=provenance_hash === nothing ? nothing : String(provenance_hash)))
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
            push!(current_occurrences, MotifOccurrence(sequence_name, start, strand, 0.0, pvalue, BioSequence{DNAAlphabet}(site_candidate)))
        end
    end

    finalize_current()
    _ctx !== nothing && register_provenance!(_ctx, "read_alignace"; parents=String[], parameters=(source=provenance_source, profile_count=length(motif_profiles), hash=provenance_hash))
    return motif_profiles
end

function _jaspar_parse_row(line::AbstractString)
    values = Float64[]
    for match in eachmatch(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", line)
        push!(values, parse(Float64, match.match))
    end
    return values
end

function read_jaspar(filepath::String; alphabet::String="ACGT", prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    if _ctx === nothing
        open(filepath, "r") do io
            return read_jaspar(io; alphabet=alphabet)
        end
    end
    raw_bytes = read(filepath)
    provenance_hash = bytes2hex(sha256(raw_bytes))


    return read_jaspar(IOBuffer(raw_bytes); alphabet=alphabet, _ctx=_ctx, provenance_hash=provenance_hash, provenance_source=filepath)
end

function read_jaspar(io::IO; alphabet::String="ACGT", provenance_hash::Union{Nothing,AbstractString}=nothing, provenance_source::AbstractString="read_jaspar", prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
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

            counts = MotifCounts{_motif_infer_alphabet(String(motif_alphabet))}(Vector{UInt8}(codeunits(String(motif_alphabet))), counts_matrix)
            pwm = motif_pwm(counts; pseudocount=0.0)
            metadata = Dict{String,String}("source" => "JASPAR")
            !isempty(jaspar_id) && (metadata["jaspar_id"] = jaspar_id)
            !isempty(jaspar_name) && (metadata["name"] = jaspar_name)
            profile_name = !isempty(jaspar_name) ? jaspar_name : jaspar_id
            profile_id = _ctx === nothing ? nothing : register_provenance!(_ctx, "read_jaspar"; parents=String[], parameters=(source=provenance_source, motif_name=profile_name, hash=provenance_hash)).id
            push!(profiles, _motif_finalize_profile(profile_name, motif_alphabet, counts_matrix, pwm.values, MotifOccurrence[], metadata; provenance_id=profile_id, provenance_hash=provenance_hash === nothing ? nothing : String(provenance_hash)))
            continue
        end

        i += 1
    end

    _ctx !== nothing && register_provenance!(_ctx, "read_jaspar"; parents=String[], parameters=(source=provenance_source, profile_count=length(profiles), hash=provenance_hash))
    return profiles
end

# ---- Typed BioSequence{DNAAlphabet} dispatch ---------------------------------
# These overloads ensure that typed DNA sequences can be used directly with
# motif scanning functions without explicit String conversion.
