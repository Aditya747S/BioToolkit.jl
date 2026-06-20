# ==============================================================================
# coevolution.jl — Co-evolutionary contact inference, DCA, and MSA analytics
#
# References:
#   - Ekeberg et al. (2013) Phys Rev E 87:012707 (PLM-DCA)
#   - Morcos et al. (2011) PNAS 108:E1294-E1301 (Direct Coupling Analysis)
#   - Shannon (1948) Bell System Tech J 27:379-423 (entropy)
#   - Weigt et al. (2009) PNAS 106:67-72 (MI-based pairing)
#   - Dunn et al. (2007) Bioinformatics 23:1684-1691 (MIPIE / mutual info correction)
#   - Jones et al. (2012) Bioinformatics 28:184-190 (PSICOV)
#   - Ledoit & Wolf (2004) J Multivariate Anal 88:365-411 (covariance shrinkage)
# ==============================================================================

module Coevolution

using LinearAlgebra
using Random
using Statistics

using ..BioToolkit: Atom, Chain, Model, MultipleSequenceAlignment, Residue, SeqRecordLite, Structure
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!

@inline function _register_coevolution_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export ContactMap, PseudoLikelihoodModel
export filter_alignment_for_dca, sequence_reweighting
export fit_pseudolikelihood_model, compute_contact_scores, predict_contact_map, top_contact_pairs
export fold_from_contacts

export mutual_information_contacts
export direct_information_contacts
export column_conservation_scores
export sequence_logo_entropy
export evolutionary_coupling_network
export contact_enrichment_statistics
export shrinkage_precision_contacts
export phylogenetic_correction
export positional_covariation_matrix
export gap_analysis
export contact_precision_recall
export alignment_quality_report

# ---------------------------------------------------------------------------
# Core data structures
# ---------------------------------------------------------------------------

struct ContactMap
    scores::Matrix{Float64}
    residue_ids::Vector{Int}
end

struct PseudoLikelihoodModel
    fields::Matrix{Float64}
    couplings::Array{Float64,4}
    alphabet::Vector{Char}
    weights::Vector{Float64}
    effective_sequences::Float64
    raw_scores::Matrix{Float64}
    apc_scores::Matrix{Float64}
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

function _alignment_strings(alignment::MultipleSequenceAlignment)
    n = length(alignment)
    n > 0 || throw(ArgumentError("alignment must contain at least one sequence"))
    return [uppercase(String(alignment.records[index].sequence)) for index in 1:n]
end

function filter_alignment_for_dca(alignment::MultipleSequenceAlignment; max_gap_fraction::Real=0.5, min_sequence_coverage::Real=0.5)
    strings = _alignment_strings(alignment)
    l = ncodeunits(strings[1])
    all(ncodeunits(seq) == l for seq in strings) || throw(ArgumentError("all alignment sequences must have equal length"))

    keep_columns = Int[]
    for col in 1:l
        gap_fraction = mean(seq[col] == '-' for seq in strings)
        if gap_fraction <= max_gap_fraction
            push!(keep_columns, col)
        end
    end
    isempty(keep_columns) && throw(ArgumentError("all columns were removed by the gap filter"))

    filtered_records = SeqRecordLite[]
    for record in alignment.records
        seq = uppercase(String(record.sequence))
        filtered_seq = String([seq[col] for col in keep_columns])
        sequence_coverage = mean(char != '-' for char in filtered_seq)
        if sequence_coverage >= min_sequence_coverage
            push!(filtered_records, SeqRecordLite(filtered_seq; identifier=record.identifier, name=record.name, description=record.description))
        end
    end

    isempty(filtered_records) && throw(ArgumentError("all sequences were removed by coverage filtering"))

    return MultipleSequenceAlignment(filtered_records)
end

function _alphabet(strings::Vector{String})
    states = Set{Char}()
    for seq in strings
        foreach(ch -> push!(states, ch), seq)
    end
    push!(states, '-')

    chars = collect(states)
    sort!(chars)

    if '-' in chars
        deleteat!(chars, findfirst(==('-'), chars))
        pushfirst!(chars, '-')
    end
    return chars
end

function _encode_alignment(strings::Vector{String}, alphabet::Vector{Char})
    n = length(strings)
    l = ncodeunits(strings[1])
    lookup = Dict{Char,Int}(char => index for (index, char) in enumerate(alphabet))
    matrix = Matrix{Int}(undef, n, l)

    for i in 1:n
        seq = strings[i]
        for j in 1:l
            matrix[i, j] = get(lookup, seq[j], lookup['-'])
        end
    end

    return matrix
end

function sequence_reweighting(encoded_alignment::AbstractMatrix{<:Integer}; identity_threshold::Real=0.8)
    n = size(encoded_alignment, 1)
    l = size(encoded_alignment, 2)

    weights = ones(Float64, n)
    for i in 1:n
        similar_count = 0
        for j in 1:n
            identity = count(encoded_alignment[i, pos] == encoded_alignment[j, pos] for pos in 1:l) / l
            if identity >= identity_threshold
                similar_count += 1
            end
        end
        weights[i] = 1.0 / max(similar_count, 1)
    end

    return weights
end

function _one_hot(encoded_alignment::AbstractMatrix{<:Integer}, q::Int)
    n, l = size(encoded_alignment)
    one_hot = zeros(Float64, n, l * q)
    for i in 1:n
        for j in 1:l
            state = encoded_alignment[i, j]
            one_hot[i, (j - 1) * q + state] = 1.0
        end
    end
    return one_hot
end

function _apc_correct(scores::Matrix{Float64})
    row_mean = vec(mean(scores, dims=2))
    col_mean = vec(mean(scores, dims=1))
    overall = mean(scores)
    overall <= 0 && return copy(scores)

    corrected = copy(scores)
    for i in 1:size(scores, 1)
        for j in 1:size(scores, 2)
            corrected[i, j] = scores[i, j] - (row_mean[i] * col_mean[j]) / overall
        end
    end

    for i in axes(corrected, 1)
        corrected[i, i] = 0.0
    end

    return corrected
end

function fit_pseudolikelihood_model(alignment::MultipleSequenceAlignment; max_gap_fraction::Real=0.5, min_sequence_coverage::Real=0.5, identity_threshold::Real=0.8, regularization::Real=0.01, pseudocount::Real=0.5)
    filtered = filter_alignment_for_dca(alignment; max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage)
    strings = _alignment_strings(filtered)
    alphabet = _alphabet(strings)
    encoded = _encode_alignment(strings, alphabet)

    n, l = size(encoded)
    q = length(alphabet)

    weights = sequence_reweighting(encoded; identity_threshold=identity_threshold)
    effective_n = sum(weights)

    one_hot = _one_hot(encoded, q)
    weighted_mean = vec((weights' * one_hot) ./ effective_n)
    centered = one_hot .- reshape(weighted_mean, 1, :)
    weighted_centered = centered .* reshape(sqrt.(weights), :, 1)

    covariance = (weighted_centered' * weighted_centered) ./ effective_n
    covariance += regularization * I
    precision = inv(Symmetric(covariance))

    frequencies = zeros(Float64, l, q)
    for i in 1:n
        for j in 1:l
            frequencies[j, encoded[i, j]] += weights[i]
        end
    end
    frequencies .+= pseudocount
    frequencies ./= sum(frequencies, dims=2)
    fields = log.(frequencies)

    couplings = zeros(Float64, l, l, q, q)
    raw_scores = zeros(Float64, l, l)

    for i in 1:l-1
        i_range = (i - 1) * q + 1:i * q
        for j in i+1:l
            j_range = (j - 1) * q + 1:j * q
            block = -Matrix(precision[i_range, j_range])
            couplings[i, j, :, :] .= block
            couplings[j, i, :, :] .= block'

            score = norm(block)
            raw_scores[i, j] = score
            raw_scores[j, i] = score
        end
    end

    apc_scores = _apc_correct(raw_scores)

    return PseudoLikelihoodModel(fields, couplings, alphabet, weights, effective_n, raw_scores, apc_scores)
end

function compute_contact_scores(model::PseudoLikelihoodModel; apc::Bool=true, min_separation::Integer=5)
    scores = apc ? copy(model.apc_scores) : copy(model.raw_scores)
    l = size(scores, 1)

    for i in 1:l
        scores[i, i] = 0.0
        lo = max(1, i - min_separation)
        hi = min(l, i + min_separation)
        for j in lo:hi
            scores[i, j] = 0.0
            scores[j, i] = 0.0
        end
    end

    return scores
end

function _normalize_scores(scores::Matrix{Float64})
    upper_values = Float64[]
    l = size(scores, 1)
    for i in 1:l-1
        for j in i+1:l
            if scores[i, j] > 0
                push!(upper_values, scores[i, j])
            end
        end
    end

    isempty(upper_values) && return zeros(Float64, size(scores))
    min_val = minimum(upper_values)
    max_val = maximum(upper_values)
    scale = max(max_val - min_val, eps(Float64))

    normalized = max.((scores .- min_val) ./ scale, 0.0)
    for i in 1:l
        normalized[i, i] = 0.0
    end
    return normalized
end

function predict_contact_map(alignment::MultipleSequenceAlignment; top_l::Union{Nothing,Int}=nothing, min_separation::Integer=5, return_model::Bool=false, kwargs...)
    model = fit_pseudolikelihood_model(alignment; kwargs...)
    scores = compute_contact_scores(model; apc=true, min_separation=min_separation)

    if top_l !== nothing && top_l > 0
        l = size(scores, 1)
        pairs = Tuple{Int,Int,Float64}[]
        for i in 1:l-1
            for j in i+1:l
                if scores[i, j] > 0
                    push!(pairs, (i, j, scores[i, j]))
                end
            end
        end
        sort!(pairs; by=pair -> pair[3], rev=true)

        keep = Set{Tuple{Int,Int}}()
        for pair in pairs[1:min(top_l, length(pairs))]
            push!(keep, (pair[1], pair[2]))
        end

        filtered = zeros(Float64, l, l)
        for (i, j) in keep
            filtered[i, j] = scores[i, j]
            filtered[j, i] = scores[j, i]
        end
        scores = filtered
    end

    contact_map = ContactMap(_normalize_scores(scores), collect(1:size(scores, 1)))
    if return_model
        return contact_map, model
    end

    return contact_map
end

function top_contact_pairs(contact_map::ContactMap; top_n::Integer=10, min_separation::Integer=5)
    l = size(contact_map.scores, 1)
    pairs = Tuple{Int,Int,Float64}[]

    for i in 1:l-1
        for j in i+1:l
            if abs(i - j) >= min_separation && contact_map.scores[i, j] > 0
                push!(pairs, (i, j, contact_map.scores[i, j]))
            end
        end
    end

    sort!(pairs; by=pair -> pair[3], rev=true)

    return pairs[1:min(top_n, length(pairs))]
end

function fold_from_contacts(contact_map::ContactMap; top_n::Union{Nothing,Int}=nothing, contact_distance::Real=7.5, backbone_distance::Real=3.8, iterations::Integer=2000, learning_rate::Real=0.01, seed::Integer=1)
    l = size(contact_map.scores, 1)
    l >= 3 || throw(ArgumentError("contact map must contain at least three residues"))

    candidate_pairs = Tuple{Int,Int,Float64}[]
    for i in 1:l-1
        for j in i+1:l
            weight = contact_map.scores[i, j]
            if weight > 0
                push!(candidate_pairs, (i, j, weight))
            end
        end
    end
    sort!(candidate_pairs; by=pair -> pair[3], rev=true)

    if top_n !== nothing && top_n > 0
        candidate_pairs = candidate_pairs[1:min(top_n, length(candidate_pairs))]
    end

    rng = MersenneTwister(seed)
    coords = randn(rng, l, 3)
    for i in 1:l
        coords[i, 1] += 3.8 * (i - 1)
    end

    for _ in 1:iterations
        gradients = zeros(Float64, l, 3)

        for (i, j, weight) in candidate_pairs
            diff = coords[i, :] .- coords[j, :]
            distance = sqrt(sum(diff .^ 2) + 1e-10)
            error = distance - contact_distance
            g = (2.0 * weight * error / distance) .* diff
            gradients[i, :] .+= g
            gradients[j, :] .-= g
        end

        for i in 1:l-1
            diff = coords[i, :] .- coords[i + 1, :]
            distance = sqrt(sum(diff .^ 2) + 1e-10)
            error = distance - backbone_distance
            g = (0.3 * 2.0 * error / distance) .* diff
            gradients[i, :] .+= g
            gradients[i + 1, :] .-= g
        end

        coords .-= learning_rate .* gradients
        coords .-= mean(coords, dims=1)
    end

    model = Model(1)
    chain = Chain("A")

    for residue_idx in 1:l
        atom = Atom(
            residue_idx,
            "CA",
            coords[residue_idx, 1],
            coords[residue_idx, 2],
            coords[residue_idx, 3];
            element="C",
            occupancy=1.0,
            bfactor=10.0,
            hetatm=false)
        residue = Residue("GLY", residue_idx, ' ', Atom[atom])
        push!(chain.residues, residue)
    end

    push!(model.chains, chain)
    structure = Structure("PredictedContactFold")
    push!(structure.models, model)

    return structure
end

# ---------------------------------------------------------------------------
# Mutual Information Contacts
# ---------------------------------------------------------------------------

"""
    mutual_information_contacts(alignment; pseudocount=0.5, min_separation=5, apc=true)

Compute residue-residue mutual information (MI) scores from a multiple
sequence alignment, analogous to `PSICOV` MI and `DCA` raw MI outputs.

Optionally applies Average Product Correction (APC) to remove phylogenetic
background signal (Dunn et al. 2007).

Returns a `ContactMap` with MI-based scores.
"""
function mutual_information_contacts(
    alignment::MultipleSequenceAlignment;
    pseudocount::Real=0.5,
    min_separation::Int=5,
    apc::Bool=true,
    max_gap_fraction::Real=0.5,
    min_sequence_coverage::Real=0.5,
    identity_threshold::Real=0.8)
    filtered = filter_alignment_for_dca(alignment; max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage)
    strings = _alignment_strings(filtered)
    alphabet = _alphabet(strings)
    encoded = _encode_alignment(strings, alphabet)
    n, l = size(encoded)
    q = length(alphabet)

    weights = sequence_reweighting(encoded; identity_threshold=identity_threshold)
    eff_n = max(sum(weights), eps(Float64))

    # Marginal frequencies with pseudocount
    f1 = zeros(Float64, l, q)   # single-site
    f2 = zeros(Float64, l, l, q, q)  # pairwise

    for s in 1:n
        w = weights[s]
        for i in 1:l
            f1[i, encoded[s, i]] += w
        end
        for i in 1:l, j in i+1:l
            f2[i, j, encoded[s, i], encoded[s, j]] += w
            f2[j, i, encoded[s, j], encoded[s, i]] += w
        end
    end

    # Normalise + pseudocount
    pc = Float64(pseudocount)
    f1 = (f1 .+ pc / q) ./ (eff_n + pc)
    f2 = (f2 .+ pc / (q * q)) ./ (eff_n + pc)

    mi_scores = zeros(Float64, l, l)
    for i in 1:l-1
        for j in i+1:l
            mi = 0.0
            for a in 1:q, b in 1:q
                pij = f2[i, j, a, b]
                pi  = f1[i, a]
                pj  = f1[j, b]
                if pij > 0 && pi > 0 && pj > 0
                    mi += pij * log(pij / (pi * pj))
                end
            end
            mi_scores[i, j] = mi
            mi_scores[j, i] = mi
        end
    end

    apc_mi = apc ? _apc_correct(mi_scores) : mi_scores

    # Apply separation mask
    for i in 1:l
        for j in max(1, i - min_separation):min(l, i + min_separation)
            apc_mi[i, j] = 0.0
        end
    end

    result = ContactMap(_normalize_scores(apc_mi), collect(1:l))
    _ctx = active_provenance_context()


    return _register_coevolution_result!(_ctx, result, "mutual_information_contacts"; parents=String[], parameters=(n_seqs=n, n_cols=l, apc=apc, min_separation=min_separation))
end

# ---------------------------------------------------------------------------
# Direct Information (DI) Contacts — mean-field DCA
# ---------------------------------------------------------------------------

"""
    direct_information_contacts(alignment; min_separation=5, regularization=0.05)

Compute Direct Information (DI) scores using the mean-field approximation
to Direct Coupling Analysis (mfDCA), analogous to the `evfold` DI output.

Mean-field DCA inverts the connected correlation matrix C to obtain a direct
coupling matrix J, from which DI is computed via 2-site marginals.

Returns a `ContactMap` with DI scores.
"""
function direct_information_contacts(
    alignment::MultipleSequenceAlignment;
    min_separation::Int=5,
    regularization::Real=0.05,
    pseudocount::Real=0.5,
    max_gap_fraction::Real=0.5,
    min_sequence_coverage::Real=0.5,
    identity_threshold::Real=0.8)
    filtered = filter_alignment_for_dca(alignment; max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage)
    strings = _alignment_strings(filtered)
    alphabet = _alphabet(strings)
    encoded  = _encode_alignment(strings, alphabet)
    n, l     = size(encoded)
    q        = length(alphabet)

    weights  = sequence_reweighting(encoded; identity_threshold=identity_threshold)
    eff_n    = max(sum(weights), eps(Float64))

    # Single-site and connected correlation matrices
    one_hot = _one_hot(encoded, q)
    wmean   = vec((weights' * one_hot) ./ eff_n)
    centered = one_hot .- reshape(wmean, 1, :)
    wcentered = centered .* reshape(sqrt.(weights), :, 1)

    C = (wcentered' * wcentered) ./ eff_n
    C_reg = C + Float64(regularization) * I
    J = -inv(Symmetric(C_reg))   # direct coupling matrix (N*q × N*q)

    # Direct Information per pair
    di_scores = zeros(Float64, l, l)
    for i in 1:l-1
        i_range = (i - 1) * q + 1:i * q
        pi = wmean[i_range] .+ Float64(pseudocount) / q
        pi ./= sum(pi)
        for j in i+1:l
            j_range = (j - 1) * q + 1:j * q
            pj = wmean[j_range] .+ Float64(pseudocount) / q
            pj ./= sum(pj)

            Jij = J[i_range, j_range]
            # Compute 2-site marginals via message passing (1-step approximation)
            pij = zeros(Float64, q, q)
            for a in 1:q, b in 1:q
                pij[a, b] = pi[a] * pj[b] * exp(Jij[a, b])
            end
            pij ./= max(sum(pij), eps(Float64))

            di = 0.0
            for a in 1:q, b in 1:q
                if pij[a, b] > 0
                    di += pij[a, b] * log(pij[a, b] / max(pi[a] * pj[b], eps(Float64)))
                end
            end
            di_scores[i, j] = di
            di_scores[j, i] = di
        end
    end

    apc_di = _apc_correct(di_scores)
    for i in 1:l
        for j in max(1, i - min_separation):min(l, i + min_separation)
            apc_di[i, j] = 0.0
        end
    end

    result = ContactMap(_normalize_scores(apc_di), collect(1:l))
    _ctx = active_provenance_context()


    return _register_coevolution_result!(_ctx, result, "direct_information_contacts"; parents=String[], parameters=(n_seqs=n, n_cols=l, min_separation=min_separation, regularization=Float64(regularization)))
end

# ---------------------------------------------------------------------------
# Column Conservation Scores
# ---------------------------------------------------------------------------

"""
    column_conservation_scores(alignment; pseudocount=0.5, gap_penalise=true)

Compute per-column conservation (Shannon entropy–based) scores for a multiple
sequence alignment, analogous to `scorecons` / `al2co`.

Conservation = 1 - H / H_max, where H is per-column Shannon entropy.

Returns a `Vector{Float64}` of conservation scores (0=variable, 1=identical).
"""
function column_conservation_scores(
    alignment::MultipleSequenceAlignment;
    pseudocount::Real=0.5,
    gap_penalise::Bool=true)
    strings = _alignment_strings(alignment)
    l = ncodeunits(strings[1])
    all(ncodeunits(s) == l for s in strings) || throw(ArgumentError("all sequences must have equal length"))
    n = length(strings)

    scores = zeros(Float64, l)
    for col in 1:l
        char_counts = Dict{Char,Float64}()
        for seq in strings
            c = seq[col]
            (!gap_penalise && c == '-') && continue
            char_counts[c] = get(char_counts, c, 0.0) + 1.0
        end
        total = sum(values(char_counts)) + Float64(pseudocount) * length(char_counts)
        total <= 0 && continue
        H = 0.0
        for v in values(char_counts)
            p = (v + Float64(pseudocount) / length(char_counts)) / total
            p > 0 && (H -= p * log2(p))
        end
        k = length(char_counts)
        H_max = k > 1 ? log2(Float64(k)) : 1.0
        scores[col] = clamp(1.0 - H / max(H_max, eps(Float64)), 0.0, 1.0)
    end
    result = scores
    _ctx = active_provenance_context()


    return _register_coevolution_result!(_ctx, result, "column_conservation_scores"; parents=String[], parameters=(n_seqs=n, n_cols=l, gap_penalise=gap_penalise))
end

# ---------------------------------------------------------------------------
# Sequence Logo Entropy
# ---------------------------------------------------------------------------

"""
    sequence_logo_entropy(alignment; pseudocount=0.5, information_content=true)

Compute per-column Shannon entropy and information content for generating
sequence logos, analogous to `seqLogo` / `ggseqlogo` in Bioconductor.

Returns a `DataFrame` with columns: `position`, `entropy`, `information_content`,
plus one column per amino acid/nucleotide giving its relative frequency.
"""
function sequence_logo_entropy(
    alignment::MultipleSequenceAlignment;
    pseudocount::Real=0.5,
    information_content::Bool=true)
    using_df = @isdefined(DataFrames)
    strings = _alignment_strings(alignment)
    l = ncodeunits(strings[1])
    all(ncodeunits(s) == l for s in strings) || throw(ArgumentError("all sequences must have equal length"))
    n = length(strings)

    # Collect all characters (non-gap)
    all_chars = sort!(unique(collect(Iterators.flatten(strings))))
    filter!(c -> c != '-', all_chars)

    pos_vec  = Int[]
    ent_vec  = Float64[]
    ic_vec   = Float64[]
    freq_mat = zeros(Float64, l, length(all_chars))

    for col in 1:l
        counts = Dict{Char,Float64}()
        for seq in strings
            c = seq[col]
            counts[c] = get(counts, c, 0.0) + 1.0
        end
        total = n + Float64(pseudocount) * length(all_chars)
        H = 0.0
        for (k_idx, c) in enumerate(all_chars)
            cnt = get(counts, c, 0.0) + Float64(pseudocount)
            p   = cnt / total
            freq_mat[col, k_idx] = p
            p > 0 && (H -= p * log2(p))
        end
        H_max = log2(max(length(all_chars), 2.0))
        IC = H_max - H
        push!(pos_vec, col)
        push!(ent_vec, H)
        push!(ic_vec, IC)
    end

    # Build a simple NamedTuple table (avoid DataFrames dependency here)
    char_freqs = NamedTuple{Tuple(Symbol.(string.(all_chars)))}(Tuple(freq_mat[:, k] for k in 1:length(all_chars)))
    result = merge((position=pos_vec, entropy=ent_vec, information_content=ic_vec), char_freqs)
    _ctx = active_provenance_context()


    return _register_coevolution_result!(_ctx, result, "sequence_logo_entropy"; parents=String[], parameters=(n_seqs=n, n_cols=l, n_chars=length(all_chars)))
end

# ---------------------------------------------------------------------------
# Evolutionary Coupling Network
# ---------------------------------------------------------------------------

"""
    evolutionary_coupling_network(contact_map; score_threshold=0.4, min_separation=5)

Build an evolutionary coupling network from a `ContactMap`, where edges connect
residue pairs with scores above `score_threshold`.

Analogous to the network view in `EVcouplings` / `gremlin`.

Returns `(edges=NamedTuple, degree=Vector, hub_residues=Vector{Int})`.
"""
function evolutionary_coupling_network(
    contact_map::ContactMap;
    score_threshold::Real=0.4,
    min_separation::Int=5)
    scores = contact_map.scores
    l = size(scores, 1)

    node_i = Int[]
    node_j = Int[]
    weights = Float64[]

    for i in 1:l-1, j in (i+1):l
        abs(j - i) < min_separation && continue
        scores[i, j] > Float64(score_threshold) || continue
        push!(node_i, i)
        push!(node_j, j)
        push!(weights, scores[i, j])
    end

    degree = zeros(Int, l)
    for (i, j) in zip(node_i, node_j)
        degree[i] += 1
        degree[j] += 1
    end

    hub_threshold = quantile(Float64.(degree[degree .> 0]), 0.9)
    hub_residues = findall(d -> d >= hub_threshold, degree)

    edges = (node_i=node_i, node_j=node_j, weights=weights)
    result = (edges=edges, degree=degree, hub_residues=hub_residues, n_edges=length(node_i))
    _ctx = active_provenance_context()


    return _register_coevolution_result!(_ctx, result, "evolutionary_coupling_network"; parents=String[], parameters=(l=l, score_threshold=Float64(score_threshold), n_edges=length(node_i), n_hubs=length(hub_residues)))
end

# ---------------------------------------------------------------------------
# Contact Enrichment Statistics
# ---------------------------------------------------------------------------

"""
    contact_enrichment_statistics(contact_map, true_contacts; top_fractions=[0.5,1.0,2.0])

Evaluate predicted contacts against known structural contacts by computing
precision at different L/k contact fractions, analogous to `metapsicov` and
`EVcouplings` benchmarking output.

`true_contacts`: binary matrix (l × l) where 1 = true structural contact.
`top_fractions`: multiples of L (sequence length) to evaluate at.

Returns a `NamedTuple` with `precision_at_k`, `ppv_auc`, `mcc` per threshold.
"""
function contact_enrichment_statistics(contact_map::ContactMap, true_contacts::AbstractMatrix{<:Real};
    top_fractions=[0.5, 1.0, 2.0],
    min_separation::Int=5)
    l = size(contact_map.scores, 1)
    size(true_contacts) == (l, l) || throw(DimensionMismatch("true_contacts must match contact map size"))

    # Collect all predicted pairs
    pairs = Tuple{Int,Int,Float64}[]
    for i in 1:l-1, j in (i+1):l
        abs(j - i) < min_separation && continue
        push!(pairs, (i, j, contact_map.scores[i, j]))
    end
    sort!(pairs; by=p -> p[3], rev=true)

    # True contact pairs
    true_set = Set{Tuple{Int,Int}}()
    for i in 1:l-1, j in (i+1):l
        true_contacts[i, j] > 0.5 && push!(true_set, (i, j))
    end
    n_true = length(true_set)

    prec_at_k = Dict{Float64,Float64}()
    for frac in top_fractions
        k = max(1, round(Int, frac * l))
        k = min(k, length(pairs))
        tp = count(p -> (p[1], p[2]) in true_set, pairs[1:k])
        prec_at_k[frac] = k > 0 ? tp / k : 0.0
    end

    # PPV AUC (area under precision-recall for top L pairs)
    max_k = min(l, length(pairs))
    precisions = Float64[]
    for k in 1:max_k
        tp = count(p -> (p[1], p[2]) in true_set, pairs[1:k])
        push!(precisions, tp / k)
    end
    ppv_auc = max_k > 0 ? mean(precisions) : 0.0

    # MCC at L contacts
    k_L = min(l, length(pairs))
    tp_L = count(p -> (p[1], p[2]) in true_set, pairs[1:k_L])
    fp_L = k_L - tp_L
    fn_L = max(n_true - tp_L, 0)
    n_pairs = (l * (l - 1)) ÷ 2 - sum(abs(p[2]-p[1]) < min_separation ? 1 : 0 for p in pairs; init=0)
    tn_L = max(n_pairs - tp_L - fp_L - fn_L, 0)
    mcc_denom = sqrt(Float64((tp_L + fp_L) * (tp_L + fn_L) * (tn_L + fp_L) * (tn_L + fn_L)))
    mcc = mcc_denom > 0 ? (tp_L * tn_L - fp_L * fn_L) / mcc_denom : 0.0

    return (precision_at_k=prec_at_k, ppv_auc=ppv_auc, mcc=mcc, n_true_contacts=n_true)
end

# ---------------------------------------------------------------------------
# Ledoit-Wolf Shrinkage for Precision Contacts (PSICOV-style)
# ---------------------------------------------------------------------------

"""
    shrinkage_precision_contacts(alignment; shrinkage=:ledoit_wolf, min_separation=5)

Compute contact scores from a regularised (shrunk) precision matrix, analogous
to `PSICOV` (Jones et al. 2012). Uses Ledoit-Wolf optimal shrinkage to
improve covariance matrix conditioning.

Returns a `ContactMap` with shrinkage-precision scores.
"""
function shrinkage_precision_contacts(
    alignment::MultipleSequenceAlignment;
    shrinkage::Symbol=:ledoit_wolf,
    min_separation::Int=5,
    max_gap_fraction::Real=0.5,
    min_sequence_coverage::Real=0.5,
    identity_threshold::Real=0.8)
    filtered = filter_alignment_for_dca(alignment; max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage)
    strings = _alignment_strings(filtered)
    alphabet = _alphabet(strings)
    encoded = _encode_alignment(strings, alphabet)
    n, l = size(encoded)
    q = length(alphabet)

    weights = sequence_reweighting(encoded; identity_threshold=identity_threshold)
    eff_n = max(sum(weights), eps(Float64))

    one_hot = _one_hot(encoded, q)
    wmean = vec((weights' * one_hot) ./ eff_n)
    centered = one_hot .- reshape(wmean, 1, :)
    wcentered = centered .* reshape(sqrt.(weights), :, 1)

    S = (wcentered' * wcentered) ./ eff_n   # sample covariance
    p = size(S, 1)

    # Ledoit-Wolf shrinkage intensity
    if shrinkage == :ledoit_wolf
        μ = tr(S) / p
        # Simplified Oracle Approximating Shrinkage (OAS) estimate
        rho = clamp(((1 - 2/p) * tr(S * S) + tr(S)^2) / ((eff_n + 1 - 2/p) * max(tr(S * S) - tr(S)^2 / p, eps(Float64))), 0.0, 1.0)
        Σ_shrunk = (1 - rho) .* S .+ rho * μ * I
    else
        reg = 0.05
        Σ_shrunk = S + reg * I
    end

    precision = inv(Symmetric(Matrix(Σ_shrunk)))

    scores = zeros(Float64, l, l)
    for i in 1:l-1
        i_range = (i - 1) * q + 1:i * q
        for j in i+1:l
            j_range = (j - 1) * q + 1:j * q
            sc = norm(precision[i_range, j_range])
            scores[i, j] = sc
            scores[j, i] = sc
        end
    end

    apc = _apc_correct(scores)
    for i in 1:l
        for j in max(1, i - min_separation):min(l, i + min_separation)
            apc[i, j] = 0.0
        end
    end

    return ContactMap(_normalize_scores(apc), collect(1:l))
end

# ---------------------------------------------------------------------------
# Phylogenetic Correction (effective sequence weighting)
# ---------------------------------------------------------------------------

"""
    phylogenetic_correction(alignment; identity_threshold=0.8, method=:henikoff)

Compute phylogenetically corrected sequence weights using column-based
position-specific weighting (Henikoff & Henikoff 1994) or standard
maximum-identity-threshold clustering.

`method`: `:henikoff` (position specific) or `:clustering` (threshold-based).

Returns `(weights=Vector{Float64}, effective_sequences=Float64)`.
"""
function phylogenetic_correction(
    alignment::MultipleSequenceAlignment;
    identity_threshold::Real=0.8,
    method::Symbol=:henikoff)
    strings = _alignment_strings(alignment)
    l = ncodeunits(strings[1])
    all(ncodeunits(s) == l for s in strings) || throw(ArgumentError("all sequences must have equal length"))
    n = length(strings)

    if method == :henikoff
        # Position-specific weighting: w_i = sum_j 1/(r_ij * s_j)
        # where r_ij = count of residue type at position j, s_j = number of different residue types
        weights = zeros(Float64, n)
        for col in 1:l
            col_chars = [strings[i][col] for i in 1:n]
            char_counts = Dict{Char,Int}()
            for c in col_chars
                char_counts[c] = get(char_counts, c, 0) + 1
            end
            s_j = length(char_counts)
            s_j == 0 && continue
            for i in 1:n
                c = col_chars[i]
                r_ij = char_counts[c]
                weights[i] += 1.0 / (r_ij * s_j)
            end
        end
        # Normalise so mean weight = 1
        w_mean = mean(weights)
        w_mean > 0 && (weights ./= w_mean)
    else
        # Standard threshold-based
        alphabet = _alphabet(strings)
        encoded  = _encode_alignment(strings, alphabet)
        weights  = sequence_reweighting(encoded; identity_threshold=identity_threshold)
        weights .*= n / max(sum(weights), eps(Float64))   # rescale to n
    end

    eff_n = sum(weights)
    return (weights=weights, effective_sequences=eff_n)
end

# ---------------------------------------------------------------------------
# Positional Covariation Matrix
# ---------------------------------------------------------------------------

"""
    positional_covariation_matrix(alignment; pseudocount=0.5, metric=:pearson)

Compute a symmetric l × l covariation (correlation) matrix between alignment
columns, analogous to `covariation` in Rfam tools.

`metric`: `:pearson` (Pearson correlation of one-hot encodings),
          `:mutual_info` (mutual information).

Returns a named tuple `(covariation=Matrix, positions=Vector{Int})`.
"""
function positional_covariation_matrix(
    alignment::MultipleSequenceAlignment;
    pseudocount::Real=0.5,
    metric::Symbol=:pearson,
    max_gap_fraction::Real=0.5,
    min_sequence_coverage::Real=0.5,
    identity_threshold::Real=0.8)
    filtered = filter_alignment_for_dca(alignment; max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage)
    strings  = _alignment_strings(filtered)
    alphabet = _alphabet(strings)
    encoded  = _encode_alignment(strings, alphabet)
    n, l     = size(encoded)
    q        = length(alphabet)

    weights  = sequence_reweighting(encoded; identity_threshold=identity_threshold)
    eff_n    = max(sum(weights), eps(Float64))

    if metric == :pearson
        # Column-wise one-hot, compute weighted Pearson correlation between columns
        # Summarise each column as its dominant-state indicator (reduces to scalar per seq)
        col_vectors = zeros(Float64, n, l)
        for i in 1:n, j in 1:l
            col_vectors[i, j] = Float64(encoded[i, j])
        end
        col_means = vec((weights' * col_vectors) ./ eff_n)
        centered  = col_vectors .- col_means'
        weighted_c = centered .* sqrt.(weights)
        cov_mat = (weighted_c' * weighted_c) ./ eff_n

        # Pearson correlation
        std_vec = sqrt.(max.(diag(cov_mat), eps(Float64)))
        cor_mat = cov_mat ./ (std_vec * std_vec')
        cor_mat[diagind(cor_mat)] .= 1.0
        return (covariation=cor_mat, positions=collect(1:l))
    else
        # Reuse MI contact scores without APC
        mi_map = mutual_information_contacts(filtered; pseudocount=pseudocount, min_separation=0, apc=false, max_gap_fraction=1.0, min_sequence_coverage=0.0, identity_threshold=identity_threshold)
        return (covariation=mi_map.scores, positions=collect(1:l))
    end
end

# ---------------------------------------------------------------------------
# Gap Analysis
# ---------------------------------------------------------------------------

"""
    gap_analysis(alignment)

Analyse gap patterns in a multiple sequence alignment, analogous to
`ggmsa` gap visualisation and `al2co` gap statistics.

Returns a named tuple with:
- `column_gap_fraction`: per-column gap fraction
- `sequence_gap_fraction`: per-sequence gap fraction
- `gap_blocks`: runs of consecutive gapped columns per sequence
- `total_gap_fraction`: overall gap fraction
"""
function gap_analysis(alignment::MultipleSequenceAlignment)
    strings = _alignment_strings(alignment)
    n = length(strings)
    l = ncodeunits(strings[1])
    all(ncodeunits(s) == l for s in strings) || throw(ArgumentError("all sequences must have equal length"))

    col_gap = [mean(strings[i][col] == '-' for i in 1:n) for col in 1:l]
    seq_gap = [mean(strings[i][col] == '-' for col in 1:l) for i in 1:n]

    # Gap blocks: (start, end, length) per sequence
    gap_blocks = Vector{Vector{Tuple{Int,Int,Int}}}(undef, n)
    for i in 1:n
        blocks = Tuple{Int,Int,Int}[]
        in_gap = false
        start_col = 0
        for col in 1:l
            is_gap = strings[i][col] == '-'
            if is_gap && !in_gap
                in_gap = true
                start_col = col
            elseif !is_gap && in_gap
                push!(blocks, (start_col, col - 1, col - start_col))
                in_gap = false
            end
        end
        in_gap && push!(blocks, (start_col, l, l - start_col + 1))
        gap_blocks[i] = blocks
    end

    total_gap = mean(col_gap)

    return (
        column_gap_fraction  = col_gap,
        sequence_gap_fraction = seq_gap,
        gap_blocks           = gap_blocks,
        total_gap_fraction   = total_gap)
end

# ---------------------------------------------------------------------------
# Contact Precision-Recall
# ---------------------------------------------------------------------------

"""
    contact_precision_recall(contact_map, true_contacts; min_separation=5, n_points=50)

Compute a precision-recall curve for a predicted `ContactMap` against known
structural contacts, analogous to the benchmarking in `EVcouplings` and CASP.

Returns `(thresholds, precision, recall, auc_pr)`.
"""
function contact_precision_recall(
    contact_map::ContactMap,
    true_contacts::AbstractMatrix{<:Real};
    min_separation::Int=5,
    n_points::Int=50)
    l = size(contact_map.scores, 1)
    size(true_contacts) == (l, l) || throw(DimensionMismatch("true_contacts must match contact map size"))

    pairs = Tuple{Int,Int,Float64}[]
    for i in 1:l-1, j in (i+1):l
        abs(j - i) < min_separation && continue
        push!(pairs, (i, j, contact_map.scores[i, j]))
    end
    sort!(pairs; by=p -> p[3], rev=true)

    true_set = Set{Tuple{Int,Int}}()
    for i in 1:l-1, j in (i+1):l
        true_contacts[i, j] > 0.5 && push!(true_set, (i, j))
    end
    n_true = length(true_set)
    n_true == 0 && return (thresholds=Float64[], precision=Float64[], recall=Float64[], auc_pr=0.0)

    thresholds = range(0.0, 1.0, length=n_points)
    prec_vec = Float64[]
    rec_vec  = Float64[]

    for thr in thresholds
        pred = [(p[1], p[2]) for p in pairs if p[3] >= thr]
        tp = count(p -> p in true_set, pred)
        fp = length(pred) - tp
        push!(prec_vec, length(pred) > 0 ? tp / length(pred) : 1.0)
        push!(rec_vec, tp / n_true)
    end

    # AUC using trapezoid rule
    auc = sum(
        0.5 * (prec_vec[k] + prec_vec[k+1]) * abs(rec_vec[k+1] - rec_vec[k])
        for k in 1:(n_points - 1)
    )

    return (thresholds=collect(thresholds), precision=prec_vec, recall=rec_vec, auc_pr=auc)
end

# ---------------------------------------------------------------------------
# Alignment Quality Report
# ---------------------------------------------------------------------------

"""
    alignment_quality_report(alignment)

Generate a comprehensive quality summary for a multiple sequence alignment,
analogous to `trimal` quality metrics and `NCBI MSA viewer` statistics.

Returns a `NamedTuple` with:
- `n_sequences`, `alignment_length`
- `effective_sequences` (at 80% identity threshold)
- `mean_pairwise_identity`, `min_identity`, `max_identity`
- `column_gap_fraction` (per-column), `mean_gap_fraction`
- `conservation_scores` (per-column), `mean_conservation`
- `entropy_per_column`
"""
function alignment_quality_report(alignment::MultipleSequenceAlignment)
    strings = _alignment_strings(alignment)
    n = length(strings)
    l = ncodeunits(strings[1])
    all(ncodeunits(s) == l for s in strings) || throw(ArgumentError("sequences must have equal length"))

    alphabet = _alphabet(strings)
    encoded  = _encode_alignment(strings, alphabet)

    weights  = sequence_reweighting(encoded; identity_threshold=0.8)
    eff_n    = sum(weights)

    # Pairwise identity
    identities = Float64[]
    for i in 1:n-1, j in (i+1):n
        id = count(encoded[i, k] == encoded[j, k] && strings[i][k] != '-' for k in 1:l) /
             max(count(strings[i][k] != '-' || strings[j][k] != '-' for k in 1:l), 1)
        push!(identities, id)
    end
    mean_id = isempty(identities) ? 0.0 : mean(identities)
    min_id  = isempty(identities) ? 0.0 : minimum(identities)
    max_id  = isempty(identities) ? 0.0 : maximum(identities)

    # Per-column stats
    col_gap  = [mean(strings[i][col] == '-' for i in 1:n) for col in 1:l]
    cons     = column_conservation_scores(alignment)
    logo     = sequence_logo_entropy(alignment)

    return (
        n_sequences           = n,
        alignment_length      = l,
        effective_sequences   = eff_n,
        mean_pairwise_identity = mean_id,
        min_pairwise_identity  = min_id,
        max_pairwise_identity  = max_id,
        column_gap_fraction    = col_gap,
        mean_gap_fraction      = mean(col_gap),
        conservation_scores    = cons,
        mean_conservation      = mean(cons),
        entropy_per_column     = logo.entropy)
end

end
