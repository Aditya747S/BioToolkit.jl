# ==============================================================================
# coevolution.jl — Co-evolutionary contact inference, DCA, and MSA analytics
#
# References:
#   - Morcos et al. (2011) PNAS 108:E1294-E1301 (Direct Coupling Analysis)
#   - Ekeberg et al. (2013) Phys Rev E 87:012707 (PLM-DCA)
#   - Shannon (1948) Bell System Tech J 27:379-423 (entropy)
#   - Weigt et al. (2009) PNAS 106:67-72 (MI-based pairing)
#   - Dunn et al. (2007) Bioinformatics 23:1684-1691 (MIPIE / mutual info correction)
#   - Jones et al. (2012) Bioinformatics 28:184-190 (PSICOV)
#   - Ledoit & Wolf (2004) J Multivariate Anal 88:365-411 (covariance shrinkage)
#   - Henikoff & Henikoff (1994) J Mol Biol 243:574-578 (position-based sequence weights)
# ==============================================================================

module Coevolution

using LinearAlgebra
using Random
using Statistics

using ..BioToolkit: Atom, Chain, Model, MultipleSequenceAlignment, Residue, SeqRecordLite, Structure, BlenderIntegrator
using ..BlenderIntegrator: BlenderContactPayload, BlenderMaterial

using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!
import ..BioToolkit: to_html, export_html

@inline function _register_coevolution_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export ContactMap, PseudoLikelihoodModel
export filter_alignment_for_dca, sequence_reweighting
export fit_pseudolikelihood_model, compute_contact_scores, predict_contact_map, top_contact_pairs
export fold_from_contacts, structure_to_contact_matrix, filter_contacts_by_sequence_distance

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

export visualize_contact_map_html
export visualize_coevolution_network_html
export visualize_alignment_quality_html
export visualize_sequence_logo_html
export to_html, export_html

# ---------------------------------------------------------------------------
# Core data structures
# ---------------------------------------------------------------------------

struct ContactMap
    scores::Matrix{Float64}
    residue_ids::Vector{Int}
end

function BlenderIntegrator.to_blender_payload(cmap::ContactMap, coords::Matrix{Float64}; tube_radius::Float64=0.2, top_n::Int=50, name::String="ContactNetwork")
    pairs_info = top_contact_pairs(cmap; top_n=top_n)
    pairs = [(p[1], p[2]) for p in pairs_info]
    scores = [p[3] for p in pairs_info]
    mat = BlenderMaterial(name=name * "_mat", color=(1.0, 0.4, 0.1, 1.0), roughness=0.2)
    return BlenderContactPayload(name, coords, pairs, scores, tube_radius, mat)
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

"""
    filter_alignment_for_dca(alignment; max_gap_fraction=0.5, min_sequence_coverage=0.5)

Filter multiple sequence alignment columns and sequences based on gap thresholds.
"""
function filter_alignment_for_dca(alignment::MultipleSequenceAlignment;
    max_gap_fraction::Real=0.5,
    min_sequence_coverage::Real=0.5,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    strings = _alignment_strings(alignment)
    l = ncodeunits(strings[1])
    all(ncodeunits(seq) == l for seq in strings) || throw(ArgumentError("all alignment sequences must have equal length"))
    n = length(strings)

    keep_columns = Int[]
    for col in 1:l
        gap_count = 0
        @inbounds for seq in strings
            if seq[col] == '-'
                gap_count += 1
            end
        end
        if (gap_count / n) <= max_gap_fraction
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

    result = MultipleSequenceAlignment(filtered_records)
    return _register_coevolution_result!(_ctx, result, "filter_alignment_for_dca"; parents=provenance_parent_ids(alignment), parameters=(max_gap_fraction=Float64(max_gap_fraction), min_sequence_coverage=Float64(min_sequence_coverage), retained_columns=length(keep_columns), retained_sequences=length(filtered_records)))
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
    lookup_table = fill(0, 256)
    for (index, char) in enumerate(alphabet)
        c_code = Int(char)
        if 1 <= c_code <= 256
            lookup_table[c_code] = index
        end
    end

    matrix = Matrix{Int}(undef, n, l)
    for i in 1:n
        seq = strings[i]
        for j in 1:l
            c_code = Int(seq[j])
            if 1 <= c_code <= 256 && lookup_table[c_code] > 0
                matrix[i, j] = lookup_table[c_code]
            else
                throw(ArgumentError("Invalid character '$(seq[j])' at sequence $i, position $j (code $c_code) not in alignment alphabet"))
            end
        end
    end

    return matrix
end

"""
    sequence_reweighting(encoded_alignment; identity_threshold=0.8)

Compute phylogenetically corrected sequence weights (1 / N_similar) using ultra-fast SIMD-friendly comparison with early exit.
"""
function sequence_reweighting(encoded_alignment::AbstractMatrix{<:Integer};
    identity_threshold::Real=0.8,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    identity_threshold > 1.0 && throw(ArgumentError("identity_threshold must be in (0, 1], got $identity_threshold"))
    identity_threshold <= 0 && throw(ArgumentError("identity_threshold must be in (0, 1], got $identity_threshold"))

    n, l = size(encoded_alignment)
    weights = ones(Float64, n)
    max_mismatches = floor(Int, (1.0 - Float64(identity_threshold)) * l + 1e-9)

    counts = ones(Int, n)
    @inbounds for i in 1:n-1
        for j in i+1:n
            mismatches = 0
            for pos in 1:l
                if encoded_alignment[i, pos] != encoded_alignment[j, pos]
                    mismatches += 1
                    if mismatches > max_mismatches
                        break
                    end
                end
            end
            if mismatches <= max_mismatches
                counts[i] += 1
                counts[j] += 1
            end
        end
    end

    for i in 1:n
        weights[i] = 1.0 / counts[i]
    end

    return _register_coevolution_result!(_ctx, weights, "sequence_reweighting"; parents=provenance_parent_ids(encoded_alignment), parameters=(n_seqs=n, n_cols=l, identity_threshold=Float64(identity_threshold)))
end

function _one_hot(encoded_alignment::AbstractMatrix{<:Integer}, q::Int)
    n, l = size(encoded_alignment)
    one_hot = zeros(Float64, n, l * q)
    @inbounds for i in 1:n
        for j in 1:l
            state = encoded_alignment[i, j]
            one_hot[i, (j - 1) * q + state] = 1.0
        end
    end
    return one_hot
end

"""
Zero-Sum Gauge transformation (Ising / Frobenius gauge) for coupling block J_ij.
Enforces sum_b J_ij(a, b) = 0 and sum_a J_ij(a, b) = 0.
"""
function _zero_sum_gauge!(J_block::AbstractMatrix{Float64})
    q1, q2 = size(J_block)
    row_means = vec(mean(J_block, dims=2))
    col_means = vec(mean(J_block, dims=1))
    total_mean = mean(J_block)
    for a in 1:q1
        for b in 1:q2
            J_block[a, b] -= row_means[a] + col_means[b] - total_mean
        end
    end
    return J_block
end

"""
Average Product Correction (APC) according to Dunn et al. (2007).
Calculates row and overall averages strictly excluding self-coupling diagonal terms.
"""
function _apc_correct(scores::Matrix{Float64})
    l = size(scores, 1)
    l <= 1 && return copy(scores)

    row_sum = zeros(Float64, l)
    for i in 1:l
        for j in 1:l
            if i != j
                row_sum[i] += scores[i, j]
            end
        end
    end
    total_sum = sum(row_sum)
    mean_val = (l * (l - 1)) > 0 ? total_sum / (l * (l - 1)) : 0.0
    mean_val <= 0 && return copy(scores)

    row_mean = row_sum ./ (l - 1)
    corrected = copy(scores)
    for i in 1:l
        for j in 1:l
            if i != j
                corrected[i, j] = scores[i, j] - (row_mean[i] * row_mean[j]) / mean_val
            else
                corrected[i, i] = 0.0
            end
        end
    end

    return corrected
end

"""
Normalize contact scores across unmasked pairs (|i - j| >= min_separation) into [0, 1].
"""
function _normalize_scores(scores::Matrix{Float64}; min_separation::Int=1)
    l = size(scores, 1)
    unmasked = Float64[]
    for i in 1:l-1
        for j in i+1:l
            if abs(i - j) >= min_separation
                push!(unmasked, scores[i, j])
            end
        end
    end

    isempty(unmasked) && return zeros(Float64, size(scores))
    min_val = minimum(unmasked)
    max_val = maximum(unmasked)
    scale = max(max_val - min_val, eps(Float64))

    normalized = zeros(Float64, l, l)
    for i in 1:l-1
        for j in i+1:l
            if abs(i - j) >= min_separation
                val = max(0.0, (scores[i, j] - min_val) / scale)
                normalized[i, j] = val
                normalized[j, i] = val
            end
        end
    end
    return normalized
end

# ---------------------------------------------------------------------------
# Pseudolikelihood & Direct Coupling Model Fitting
# ---------------------------------------------------------------------------

"""
    fit_pseudolikelihood_model(alignment; max_gap_fraction=0.5, min_sequence_coverage=0.5, identity_threshold=0.8, regularization=0.01, pseudocount=0.5, algorithm=:mean_field, weights=nothing)

Fit a Co-evolutionary Direct Coupling Analysis (DCA) model to a multiple sequence alignment.

Algorithms supported:
- `:plm` (or `:pseudolikelihood`): Pseudolikelihood Maximization (PLM-DCA; Ekeberg et al. 2013). Per-column conditional log-likelihood gradient optimization over couplings and fields.
- `:mean_field` (or `:mfdca`): Classical Mean-Field DCA (mfDCA; Morcos et al. 2011). Global inverse covariance matrix precision estimation with zero-sum gauge.
- `:direct_correlation` (or `:local_correlation`): Pairwise connected covariance Frobenius norm estimation.
"""
function fit_pseudolikelihood_model(alignment::MultipleSequenceAlignment;
    max_gap_fraction::Real=0.5,
    min_sequence_coverage::Real=0.5,
    identity_threshold::Real=0.8,
    regularization::Real=0.01,
    pseudocount::Real=0.5,
    algorithm::Symbol=:mean_field,
    weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    filtered = filter_alignment_for_dca(alignment; max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage, _ctx=_ctx)
    strings = _alignment_strings(filtered)
    alphabet = _alphabet(strings)
    encoded = _encode_alignment(strings, alphabet)

    n, l = size(encoded)
    q = length(alphabet)

    weights_vec = weights !== nothing ? Float64.(weights) : sequence_reweighting(encoded; identity_threshold=identity_threshold, _ctx=_ctx)
    length(weights_vec) == n || throw(ArgumentError("weights length ($(length(weights_vec))) must match the filtered alignment sequence count ($n), not the original alignment"))
    effective_n = max(sum(weights_vec), eps(Float64))

    frequencies = zeros(Float64, l, q)
    for i in 1:n
        w = weights_vec[i]
        for j in 1:l
            frequencies[j, encoded[i, j]] += w
        end
    end
    frequencies .+= Float64(pseudocount)
    frequencies ./= sum(frequencies, dims=2)
    fields = log.(frequencies)

    couplings = zeros(Float64, l, l, q, q)
    raw_scores = zeros(Float64, l, l)

    if algorithm in (:plm, :pseudolikelihood)
        lr = 0.05
        l2_reg = Float64(regularization)
        n_iters = 60

        Threads.@threads for r in 1:l
            h_r = copy(fields[r, :])
            J_r = zeros(Float64, l, q, q)

            for _ in 1:n_iters
                grad_h = l2_reg .* h_r
                grad_J = l2_reg .* J_r

                for k in 1:n
                    w_k = weights_vec[k]
                    yr = encoded[k, r]

                    eta = copy(h_r)
                    for s in 1:l
                        if s != r
                            ys = encoded[k, s]
                            for a in 1:q
                                eta[a] += J_r[s, a, ys]
                            end
                        end
                    end

                    max_eta = maximum(eta)
                    exp_eta = exp.(eta .- max_eta)
                    probs = exp_eta ./ sum(exp_eta)

                    for a in 1:q
                        delta = w_k * (probs[a] - (a == yr ? 1.0 : 0.0)) / effective_n
                        grad_h[a] += delta
                        for s in 1:l
                            if s != r
                                ys = encoded[k, s]
                                grad_J[s, a, ys] += delta
                            end
                        end
                    end
                end

                h_r .-= lr .* grad_h
                J_r .-= lr .* grad_J
            end

            fields[r, :] .= h_r
            for s in 1:l
                if s != r
                    block = copy(J_r[s, :, :])
                    _zero_sum_gauge!(block)
                    couplings[r, s, :, :] .= block
                end
            end
        end

        for i in 1:l-1
            for j in i+1:l
                sym_block = 0.5 .* (couplings[i, j, :, :] .+ couplings[j, i, :, :]')
                _zero_sum_gauge!(sym_block)
                couplings[i, j, :, :] .= sym_block
                couplings[j, i, :, :] .= sym_block'
                score = norm(sym_block)
                raw_scores[i, j] = score
                raw_scores[j, i] = score
            end
        end

    elseif algorithm in (:direct_correlation, :local_correlation)
        one_hot = _one_hot(encoded, q)
        for i in 1:l-1
            for j in i+1:l
                c_block = zeros(Float64, q, q)
                @inbounds for s in 1:n
                    w = weights_vec[s]
                    ai = encoded[s, i]
                    bj = encoded[s, j]
                    c_block[ai, bj] += w
                end
                c_block ./= effective_n
                c_block .-= (frequencies[i, :] * frequencies[j, :]')
                _zero_sum_gauge!(c_block)

                couplings[i, j, :, :] .= c_block
                couplings[j, i, :, :] .= c_block'
                score = norm(c_block)
                raw_scores[i, j] = score
                raw_scores[j, i] = score
            end
        end
    else
        one_hot = _one_hot(encoded, q)
        weighted_mean = vec((weights_vec' * one_hot) ./ effective_n)
        centered = one_hot .- reshape(weighted_mean, 1, :)
        weighted_centered = centered .* reshape(sqrt.(weights_vec), :, 1)

        covariance = (weighted_centered' * weighted_centered) ./ effective_n
        covariance += Float64(regularization) * I
        precision = inv(Symmetric(covariance))

        for i in 1:l-1
            i_range = (i - 1) * q + 1:i * q
            for j in i+1:l
                j_range = (j - 1) * q + 1:j * q
                block = -Matrix(precision[i_range, j_range])
                _zero_sum_gauge!(block)
                couplings[i, j, :, :] .= block
                couplings[j, i, :, :] .= block'
                score = norm(block)
                raw_scores[i, j] = score
                raw_scores[j, i] = score
            end
        end
    end

    apc_scores = _apc_correct(raw_scores)
    model = PseudoLikelihoodModel(fields, couplings, alphabet, weights_vec, effective_n, raw_scores, apc_scores)

    return _register_coevolution_result!(_ctx, model, "fit_pseudolikelihood_model"; parents=provenance_parent_ids(alignment), parameters=(n_seqs=n, n_cols=l, effective_sequences=effective_n, algorithm=algorithm))
end

"""
    compute_contact_scores(model; apc=true, min_separation=5)

Extract residue-residue contact score matrix from a fitted `PseudoLikelihoodModel`.
Applies sequence separation mask (|i - j| < min_separation) and optional Average Product Correction.
"""
function compute_contact_scores(model::PseudoLikelihoodModel;
    apc::Bool=true,
    min_separation::Integer=5,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    scores = apc ? copy(model.apc_scores) : copy(model.raw_scores)
    l = size(scores, 1)

    for i in 1:l
        scores[i, i] = 0.0
        lo = max(1, i - min_separation + 1)
        hi = min(l, i + min_separation - 1)
        for j in lo:hi
            scores[i, j] = 0.0
            scores[j, i] = 0.0
        end
    end

    result = scores
    return _register_coevolution_result!(_ctx, result, "compute_contact_scores"; parents=provenance_parent_ids(model), parameters=(apc=apc, min_separation=Int(min_separation)))
end

"""
    predict_contact_map(alignment; top_l=nothing, min_separation=5, return_model=false, kwargs...)

Predict residue-residue contact map from an alignment using direct coupling analysis.
Returns normalized contact map scores in [0, 1] range.
"""
function predict_contact_map(alignment::MultipleSequenceAlignment;
    top_l::Union{Nothing,Int}=nothing,
    min_separation::Integer=5,
    return_model::Bool=false,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx),
    kwargs...)

    model = fit_pseudolikelihood_model(alignment; _ctx=_ctx, kwargs...)
    scores = compute_contact_scores(model; apc=true, min_separation=min_separation, _ctx=_ctx)

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

    cmap_scores = _normalize_scores(scores; min_separation=Int(min_separation))
    contact_map = ContactMap(cmap_scores, collect(1:size(scores, 1)))

    if return_model
        return contact_map, model
    end

    return contact_map
end

"""
    top_contact_pairs(contact_map; top_n=10, min_separation=5)

Retrieve top N contact pairs sorted by co-evolutionary score.
"""
function top_contact_pairs(contact_map::ContactMap;
    top_n::Integer=10,
    min_separation::Integer=5,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

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
    result = pairs[1:min(top_n, length(pairs))]
    return _register_coevolution_result!(_ctx, result, "top_contact_pairs"; parents=provenance_parent_ids(contact_map), parameters=(top_n=Int(top_n), min_separation=Int(min_separation)))
end

"""
    fold_from_contacts(contact_map; top_n=nothing, contact_distance=7.5, backbone_distance=3.8, iterations=2000, learning_rate=0.01, seed=1)

Perform 3D structure generation from predicted residue contacts using distance geometry optimization.
Incorporates contact distance restraints, Cα-Cα backbone connectivity, and steric clash avoidance.
"""
function fold_from_contacts(contact_map::ContactMap;
    top_n::Union{Nothing,Int}=nothing,
    contact_distance::Real=7.5,
    backbone_distance::Real=3.8,
    iterations::Integer=2000,
    learning_rate::Real=0.01,
    seed::Integer=1,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

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

    steric_min = 3.2
    for _ in 1:iterations
        gradients = zeros(Float64, l, 3)

        # Contact restraints
        for (i, j, weight) in candidate_pairs
            diff = coords[i, :] .- coords[j, :]
            distance = sqrt(sum(diff .^ 2) + 1e-10)
            error = distance - contact_distance
            g = (2.0 * weight * error / distance) .* diff
            gradients[i, :] .+= g
            gradients[j, :] .-= g
        end

        # Backbone connectivity
        for i in 1:l-1
            diff = coords[i, :] .- coords[i + 1, :]
            distance = sqrt(sum(diff .^ 2) + 1e-10)
            error = distance - backbone_distance
            g = (0.6 * 2.0 * error / distance) .* diff
            gradients[i, :] .+= g
            gradients[i + 1, :] .-= g
        end

        # Steric clash repulsion
        steric_sq = steric_min * steric_min
        for i in 1:l-2
            xi, yi, zi = coords[i, 1], coords[i, 2], coords[i, 3]
            for j in i+2:l
                dx = xi - coords[j, 1]
                dy = yi - coords[j, 2]
                dz = zi - coords[j, 3]
                dist_sq = dx*dx + dy*dy + dz*dz
                if dist_sq < steric_sq
                    dist = sqrt(dist_sq + 1e-10)
                    overlap = steric_min - dist
                    inv_dist = 1.0 / dist
                    gx = (-overlap * inv_dist) * dx
                    gy = (-overlap * inv_dist) * dy
                    gz = (-overlap * inv_dist) * dz
                    gradients[i, 1] += gx
                    gradients[i, 2] += gy
                    gradients[i, 3] += gz
                    gradients[j, 1] -= gx
                    gradients[j, 2] -= gy
                    gradients[j, 3] -= gz
                end
            end
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

    return _register_coevolution_result!(_ctx, structure, "fold_from_contacts"; parents=provenance_parent_ids(contact_map), parameters=(residues=l, restraints=length(candidate_pairs), iterations=Int(iterations)))
end

# ---------------------------------------------------------------------------
# Mutual Information Contacts (Streaming O(q^2) Memory)
# ---------------------------------------------------------------------------

"""
    mutual_information_contacts(alignment; pseudocount=0.5, min_separation=5, apc=true, weights=nothing)

Compute residue-residue mutual information (MI) contact scores with Average Product Correction (APC).
Uses streaming pairwise contingency tables to maintain O(q^2) memory footprint (Kilobytes vs Gigabytes).
"""
function mutual_information_contacts(alignment::MultipleSequenceAlignment;
    pseudocount::Real=0.5,
    min_separation::Int=5,
    apc::Bool=true,
    max_gap_fraction::Real=0.5,
    min_sequence_coverage::Real=0.5,
    identity_threshold::Real=0.8,
    weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    filtered = filter_alignment_for_dca(alignment; max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage, _ctx=_ctx)
    strings = _alignment_strings(filtered)
    alphabet = _alphabet(strings)
    encoded = _encode_alignment(strings, alphabet)
    n, l = size(encoded)
    q = length(alphabet)

    weights_vec = weights !== nothing ? Float64.(weights) : sequence_reweighting(encoded; identity_threshold=identity_threshold, _ctx=_ctx)
    length(weights_vec) == n || throw(ArgumentError("weights length ($(length(weights_vec))) must match the filtered alignment sequence count ($n), not the original alignment"))
    eff_n = max(sum(weights_vec), eps(Float64))

    # Single-site frequencies f1
    f1 = zeros(Float64, l, q)
    for s in 1:n
        w = weights_vec[s]
        for i in 1:l
            f1[i, encoded[s, i]] += w
        end
    end
    pc = Float64(pseudocount)
    f1 = (f1 .+ pc / q) ./ (eff_n + pc)

    mi_scores = zeros(Float64, l, l)
    c_ij = zeros(Float64, q, q)

    for i in 1:l-1
        for j in i+1:l
            fill!(c_ij, 0.0)
            @inbounds for s in 1:n
                c_ij[encoded[s, i], encoded[s, j]] += weights_vec[s]
            end

            pij = (c_ij .+ pc / (q * q)) ./ (eff_n + pc)

            mi = 0.0
            for a in 1:q, b in 1:q
                p_ab = pij[a, b]
                p_a = f1[i, a]
                p_b = f1[j, b]
                if p_ab > 0 && p_a > 0 && p_b > 0
                    mi += p_ab * log(p_ab / (p_a * p_b))
                end
            end
            mi_scores[i, j] = mi
            mi_scores[j, i] = mi
        end
    end

    apc_mi = apc ? _apc_correct(mi_scores) : mi_scores

    norm_scores = _normalize_scores(apc_mi; min_separation=min_separation)
    result = ContactMap(norm_scores, collect(1:l))

    return _register_coevolution_result!(_ctx, result, "mutual_information_contacts"; parents=provenance_parent_ids(alignment), parameters=(n_seqs=n, n_cols=l, apc=apc, min_separation=min_separation))
end

# ---------------------------------------------------------------------------
# Direct Information (DI) Contacts — mean-field DCA
# ---------------------------------------------------------------------------

"""
    direct_information_contacts(alignment; min_separation=5, regularization=0.05, pseudocount=0.5, weights=nothing)

Compute Direct Information (DI) contact scores using mean-field Direct Coupling Analysis (mfDCA) with RAS / IPFP 2-site marginal consistency.
Covariance matrix C and marginals are built from pseudocount-regularized joint and marginal frequencies according to Morcos et al. (2011).
"""
function direct_information_contacts(alignment::MultipleSequenceAlignment;
    min_separation::Int=5,
    regularization::Real=0.05,
    pseudocount::Real=0.5,
    max_gap_fraction::Real=0.5,
    min_sequence_coverage::Real=0.5,
    identity_threshold::Real=0.8,
    weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    filtered = filter_alignment_for_dca(alignment; max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage, _ctx=_ctx)
    strings = _alignment_strings(filtered)
    alphabet = _alphabet(strings)
    encoded  = _encode_alignment(strings, alphabet)
    n, l     = size(encoded)
    q        = length(alphabet)

    weights_vec = weights !== nothing ? Float64.(weights) : sequence_reweighting(encoded; identity_threshold=identity_threshold, _ctx=_ctx)
    length(weights_vec) == n || throw(ArgumentError("weights length ($(length(weights_vec))) must match the filtered alignment sequence count ($n), not the original alignment"))
    eff_n    = max(sum(weights_vec), eps(Float64))
    pc       = Float64(pseudocount)

    # Pseudocount-regularized single-site frequencies Pi(a)
    f1 = zeros(Float64, l, q)
    @inbounds for s in 1:n
        w = weights_vec[s]
        for i in 1:l
            f1[i, encoded[s, i]] += w
        end
    end
    pi_mat = (f1 .+ pc / q) ./ (eff_n + pc) # size (l, q)

    # Pseudocount-regularized joint covariance C
    C = zeros(Float64, l * q, l * q)

    for i in 1:l
        i_range = (i - 1) * q + 1:i * q
        for a in 1:q
            C[i_range[a], i_range[a]] = pi_mat[i, a] * (1.0 - pi_mat[i, a])
            for a2 in 1:q
                if a != a2
                    C[i_range[a], i_range[a2]] = -pi_mat[i, a] * pi_mat[i, a2]
                end
            end
        end
    end

    c_ij = zeros(Float64, q, q)
    for i in 1:l-1
        i_range = (i - 1) * q + 1:i * q
        for j in i+1:l
            j_range = (j - 1) * q + 1:j * q
            fill!(c_ij, 0.0)
            @inbounds for s in 1:n
                c_ij[encoded[s, i], encoded[s, j]] += weights_vec[s]
            end
            pij = (c_ij .+ pc / (q * q)) ./ (eff_n + pc)
            for a in 1:q, b in 1:q
                cov_val = pij[a, b] - pi_mat[i, a] * pi_mat[j, b]
                C[i_range[a], j_range[b]] = cov_val
                C[j_range[b], i_range[a]] = cov_val
            end
        end
    end

    C_reg = C + Float64(regularization) * I
    J = -inv(Symmetric(C_reg))

    di_scores = zeros(Float64, l, l)
    for i in 1:l-1
        i_range = (i - 1) * q + 1:i * q
        pi = pi_mat[i, :]
        for j in i+1:l
            j_range = (j - 1) * q + 1:j * q
            pj = pi_mat[j, :]

            Jij = copy(J[i_range, j_range])
            _zero_sum_gauge!(Jij)

            # RAS / IPFP 2-site marginal self-consistency iterations with tolerance exit
            pij = zeros(Float64, q, q)
            for a in 1:q, b in 1:q
                pij[a, b] = pi[a] * pj[b] * exp(Jij[a, b])
            end
            pij ./= max(sum(pij), eps(Float64))

            for iter in 1:100
                r_sums = vec(sum(pij, dims=2))
                max_diff = maximum(abs.(r_sums .- pi))
                if max_diff < 1e-6 && iter > 1
                    break
                end
                for a in 1:q
                    scale = r_sums[a] > 0 ? pi[a] / r_sums[a] : 0.0
                    pij[a, :] .*= scale
                end
                c_sums = vec(sum(pij, dims=1))
                for b in 1:q
                    scale = c_sums[b] > 0 ? pj[b] / c_sums[b] : 0.0
                    pij[:, b] .*= scale
                end
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
        lo = max(1, i - min_separation + 1)
        hi = min(l, i + min_separation - 1)
        for j in lo:hi
            apc_di[i, j] = 0.0
        end
    end

    norm_scores = _normalize_scores(apc_di; min_separation=min_separation)
    result = ContactMap(norm_scores, collect(1:l))

    return _register_coevolution_result!(_ctx, result, "direct_information_contacts"; parents=provenance_parent_ids(alignment), parameters=(n_seqs=n, n_cols=l, min_separation=min_separation, regularization=Float64(regularization)))
end

# ---------------------------------------------------------------------------
# Column Conservation Scores
# ---------------------------------------------------------------------------

"""
    column_conservation_scores(alignment; pseudocount=0.5, gap_penalise=true, method=:shannon)

Compute per-column conservation scores for a multiple sequence alignment using normalized Shannon entropy or Valdar score.
"""
function column_conservation_scores(alignment::MultipleSequenceAlignment;
    pseudocount::Real=0.5,
    gap_penalise::Bool=true,
    method::Symbol=:shannon,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    strings = _alignment_strings(alignment)
    l = ncodeunits(strings[1])
    all(ncodeunits(s) == l for s in strings) || throw(ArgumentError("all sequences must have equal length"))
    n = length(strings)

    alphabet = _alphabet(strings)
    non_gap_alphabet = filter(!=('-'), alphabet)
    q_non_gap = max(length(non_gap_alphabet), 2)
    H_max = log2(Float64(q_non_gap))

    scores = zeros(Float64, l)
    for col in 1:l
        counts = Dict{Char,Float64}()
        for seq in strings
            c = seq[col]
            if !gap_penalise && c == '-'
                continue
            end
            counts[c] = get(counts, c, 0.0) + 1.0
        end

        total = sum(values(counts)) + Float64(pseudocount) * length(counts)
        total <= 0 && continue

        H = 0.0
        for (c, cnt) in counts
            p = (cnt + Float64(pseudocount) / length(counts)) / total
            if p > 0
                H -= p * log2(p)
            end
        end

        score = clamp(1.0 - (H / H_max), 0.0, 1.0)
        scores[col] = score
    end

    result = scores
    return _register_coevolution_result!(_ctx, result, "column_conservation_scores"; parents=provenance_parent_ids(alignment), parameters=(n_seqs=n, n_cols=l, gap_penalise=gap_penalise, method=method))
end

# ---------------------------------------------------------------------------
# Sequence Logo Entropy
# ---------------------------------------------------------------------------

"""
    sequence_logo_entropy(alignment; pseudocount=0.5, information_content=true)

Compute per-column Shannon entropy and Information Content (height of stack) with small-sample error correction (Miller-Madow).
Retains all sequence characters (including gaps '-') in alphabet so column probability mass correctly sums to 1.0.
"""
function sequence_logo_entropy(alignment::MultipleSequenceAlignment;
    pseudocount::Real=0.5,
    information_content::Bool=true,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    strings = _alignment_strings(alignment)
    l = ncodeunits(strings[1])
    all(ncodeunits(s) == l for s in strings) || throw(ArgumentError("all sequences must have equal length"))
    n = length(strings)

    all_chars = sort!(unique(collect(Iterators.flatten(strings))))
    q = max(length(all_chars), 2)

    pos_vec  = Int[]
    ent_vec  = Float64[]
    ic_vec   = Float64[]
    freq_mat = zeros(Float64, l, length(all_chars))

    small_sample_err = (q - 1) / (2 * log(2) * n)

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
            if p > 0
                H -= p * log2(p)
            end
        end
        H_max = log2(Float64(q))
        IC = max(0.0, H_max - (H + small_sample_err))
        push!(pos_vec, col)
        push!(ent_vec, H)
        push!(ic_vec, IC)
    end

    char_freqs = NamedTuple{Tuple(Symbol.(string.(all_chars)))}(Tuple(freq_mat[:, k] for k in 1:length(all_chars)))
    result = merge((position=pos_vec, entropy=ent_vec, information_content=ic_vec), char_freqs)

    return _register_coevolution_result!(_ctx, result, "sequence_logo_entropy"; parents=provenance_parent_ids(alignment), parameters=(n_seqs=n, n_cols=l, n_chars=length(all_chars)))
end

# ---------------------------------------------------------------------------
# Evolutionary Coupling Network
# ---------------------------------------------------------------------------

"""
    evolutionary_coupling_network(contact_map; score_threshold=0.4, min_separation=5)

Build an evolutionary coupling network connecting residue pairs with contact scores above `score_threshold`.
Identifies high-degree hub residues.
"""
function evolutionary_coupling_network(contact_map::ContactMap;
    score_threshold::Real=0.4,
    min_separation::Int=5,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

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

    active_degrees = degree[degree .> 0]
    hub_threshold = !isempty(active_degrees) ? quantile(Float64.(active_degrees), 0.9) : 0.0
    hub_residues = !isempty(active_degrees) ? findall(d -> d >= hub_threshold && d > 0, degree) : Int[]

    edges = (node_i=node_i, node_j=node_j, weights=weights)
    result = (edges=edges, degree=degree, hub_residues=hub_residues, n_edges=length(node_i))

    return _register_coevolution_result!(_ctx, result, "evolutionary_coupling_network"; parents=provenance_parent_ids(contact_map), parameters=(l=l, score_threshold=Float64(score_threshold), n_edges=length(node_i), n_hubs=length(hub_residues)))
end

# ---------------------------------------------------------------------------
# Contact Enrichment Statistics
# ---------------------------------------------------------------------------

"""
    contact_enrichment_statistics(contact_map, true_contacts; top_fractions=[0.5, 1.0, 2.0], min_separation=5)

Benchmarking predicted contacts against true structural contacts (Precision at L/k, PPV AUC, Matthews Correlation Coefficient).
"""
function contact_enrichment_statistics(contact_map::ContactMap, true_contacts::AbstractMatrix{<:Real};
    top_fractions=[0.5, 1.0, 2.0],
    min_separation::Int=5,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

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
        abs(j - i) >= min_separation && true_contacts[i, j] > 0.5 && push!(true_set, (i, j))
    end
    n_true = length(true_set)

    prec_at_k = Dict{Float64,Float64}()
    for frac in top_fractions
        k = max(1, round(Int, frac * l))
        k = min(k, length(pairs))
        tp = count(p -> (p[1], p[2]) in true_set, pairs[1:k])
        prec_at_k[frac] = k > 0 ? tp / k : 0.0
    end

    max_k = min(l, length(pairs))
    precisions = Float64[]
    for k in 1:max_k
        tp = count(p -> (p[1], p[2]) in true_set, pairs[1:k])
        push!(precisions, tp / k)
    end
    ppv_auc = max_k > 0 ? mean(precisions) : 0.0

    k_L = min(l, length(pairs))
    tp_L = count(p -> (p[1], p[2]) in true_set, pairs[1:k_L])
    fp_L = k_L - tp_L
    fn_L = max(n_true - tp_L, 0)
    n_candidate_pairs = length(pairs)
    tn_L = max(n_candidate_pairs - tp_L - fp_L - fn_L, 0)

    mcc_denom = sqrt(Float64((tp_L + fp_L) * (tp_L + fn_L) * (tn_L + fp_L) * (tn_L + fn_L)))
    mcc = mcc_denom > 0 ? (tp_L * tn_L - fp_L * fn_L) / mcc_denom : 0.0

    result = (precision_at_k=prec_at_k, ppv_auc=ppv_auc, mcc=mcc, n_true_contacts=n_true)
    return _register_coevolution_result!(_ctx, result, "contact_enrichment_statistics"; parents=provenance_parent_ids(contact_map), parameters=(l=l, min_separation=min_separation, n_true=n_true))
end

# ---------------------------------------------------------------------------
# Ledoit-Wolf Shrinkage for Precision Contacts (PSICOV)
# ---------------------------------------------------------------------------

"""
    shrinkage_precision_contacts(alignment; shrinkage=:ledoit_wolf, min_separation=5, weights=nothing)

Compute contact scores using Ledoit-Wolf optimal covariance shrinkage for matrix conditioning.
Supports shrinkage methods:
- `:ledoit_wolf`: Fourth-moment Ledoit-Wolf estimator (Ledoit & Wolf 2004).
- `:oas`: Oracle Approximating Shrinkage estimator (Chen et al. 2010).
- `:ridge`: Fixed L2 ridge regularization.
"""
function shrinkage_precision_contacts(alignment::MultipleSequenceAlignment;
    shrinkage::Symbol=:ledoit_wolf,
    min_separation::Int=5,
    max_gap_fraction::Real=0.5,
    min_sequence_coverage::Real=0.5,
    identity_threshold::Real=0.8,
    weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    filtered = filter_alignment_for_dca(alignment; max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage, _ctx=_ctx)
    strings = _alignment_strings(filtered)
    alphabet = _alphabet(strings)
    encoded = _encode_alignment(strings, alphabet)
    n, l = size(encoded)
    q = length(alphabet)

    weights_vec = weights !== nothing ? Float64.(weights) : sequence_reweighting(encoded; identity_threshold=identity_threshold, _ctx=_ctx)
    length(weights_vec) == n || throw(ArgumentError("weights length ($(length(weights_vec))) must match the filtered alignment sequence count ($n), not the original alignment"))
    eff_n = max(sum(weights_vec), eps(Float64))

    one_hot = _one_hot(encoded, q)
    wmean = vec((weights_vec' * one_hot) ./ eff_n)
    centered = one_hot .- reshape(wmean, 1, :)
    wcentered = centered .* reshape(sqrt.(weights_vec), :, 1)

    S = (wcentered' * wcentered) ./ eff_n
    p = size(S, 1)

    if shrinkage == :ledoit_wolf
        μ = tr(S) / p
        d2 = max(sum(abs2, S) - (tr(S)^2 / p), eps(Float64))
        Y = wcentered * S
        b2_bar = 0.0
        for k in 1:n
            xk = @view wcentered[k, :]
            yk = @view Y[k, :]
            norm_sq = sum(abs2, xk)
            xSx = dot(xk, yk)
            b2_bar += (norm_sq^2 - 2.0 * xSx + sum(abs2, S))
        end
        b2_bar /= (eff_n^2)
        b2 = min(d2, b2_bar)
        rho = clamp(b2 / d2, 0.0, 1.0)
        Σ_shrunk = (1.0 - rho) .* S + (rho * μ) * I
    elseif shrinkage == :oas
        μ = tr(S) / p
        tr_S2 = sum(abs2, S)
        tr_S = tr(S)
        denom = max((eff_n + 1 - 2/p) * (tr_S2 - (tr_S^2 / p)), eps(Float64))
        num = (1 - 2/p) * tr_S2 + (tr_S^2)
        rho = clamp(num / denom, 0.0, 1.0)
        Σ_shrunk = (1.0 - rho) .* S + (rho * μ) * I
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
            block = copy(precision[i_range, j_range])
            _zero_sum_gauge!(block)
            sc = norm(block)
            scores[i, j] = sc
            scores[j, i] = sc
        end
    end

    apc = _apc_correct(scores)
    for i in 1:l
        lo = max(1, i - min_separation + 1)
        hi = min(l, i + min_separation - 1)
        for j in lo:hi
            apc[i, j] = 0.0
        end
    end

    norm_scores = _normalize_scores(apc; min_separation=min_separation)
    result = ContactMap(norm_scores, collect(1:l))
    return _register_coevolution_result!(_ctx, result, "shrinkage_precision_contacts"; parents=provenance_parent_ids(alignment), parameters=(shrinkage=shrinkage, min_separation=min_separation))
end

# ---------------------------------------------------------------------------
# Phylogenetic Correction
# ---------------------------------------------------------------------------

"""
    phylogenetic_correction(alignment; identity_threshold=0.8, method=:henikoff)

Compute phylogenetically corrected sequence weights using Henikoff position-based weighting or threshold clustering.
"""
function phylogenetic_correction(alignment::MultipleSequenceAlignment;
    identity_threshold::Real=0.8,
    method::Symbol=:henikoff,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    strings = _alignment_strings(alignment)
    l = ncodeunits(strings[1])
    all(ncodeunits(s) == l for s in strings) || throw(ArgumentError("all sequences must have equal length"))
    n = length(strings)

    if method == :henikoff
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
        w_mean = mean(weights)
        w_mean > 0 && (weights ./= w_mean)
    else
        alphabet = _alphabet(strings)
        encoded  = _encode_alignment(strings, alphabet)
        weights  = sequence_reweighting(encoded; identity_threshold=identity_threshold, _ctx=_ctx)
        weights .*= n / max(sum(weights), eps(Float64))
    end

    eff_n = sum(weights)
    result = (weights=weights, effective_sequences=eff_n)
    return _register_coevolution_result!(_ctx, result, "phylogenetic_correction"; parents=provenance_parent_ids(alignment), parameters=(method=method, identity_threshold=Float64(identity_threshold)))
end

# ---------------------------------------------------------------------------
# Positional Covariation Matrix
# ---------------------------------------------------------------------------

"""
    positional_covariation_matrix(alignment; pseudocount=0.5, metric=:frobenius, weights=nothing)

Compute a symmetric l x l covariation matrix between alignment columns.
Metrics supported:
- `:frobenius` (or `:frobenius_correlation`): Categorical one-hot matrix Frobenius correlation norm ||C_ij||_F / sqrt(||C_ii||_F * ||C_jj||_F) in [0, 1].
- `:mi` (or `:mutual_information`): Direct mutual information matrix.
"""
function positional_covariation_matrix(alignment::MultipleSequenceAlignment;
    pseudocount::Real=0.5,
    metric::Symbol=:frobenius,
    max_gap_fraction::Real=0.5,
    min_sequence_coverage::Real=0.5,
    identity_threshold::Real=0.8,
    weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    filtered = filter_alignment_for_dca(alignment; max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage, _ctx=_ctx)
    strings  = _alignment_strings(filtered)
    alphabet = _alphabet(strings)
    encoded  = _encode_alignment(strings, alphabet)
    n, l     = size(encoded)
    q        = length(alphabet)

    weights_vec  = weights !== nothing ? Float64.(weights) : sequence_reweighting(encoded; identity_threshold=identity_threshold, _ctx=_ctx)
    length(weights_vec) == n || throw(ArgumentError("weights length ($(length(weights_vec))) must match the filtered alignment sequence count ($n), not the original alignment"))
    eff_n    = max(sum(weights_vec), eps(Float64))

    if metric in (:frobenius, :frobenius_correlation)
        one_hot = _one_hot(encoded, q)
        wmean   = vec((weights_vec' * one_hot) ./ eff_n)
        centered = one_hot .- reshape(wmean, 1, :)
        wcentered = centered .* reshape(sqrt.(weights_vec), :, 1)

        cov_mat = (wcentered' * wcentered) ./ eff_n
        cor_mat = zeros(Float64, l, l)

        for i in 1:l
            i_range = (i - 1) * q + 1:i * q
            cov_ii = norm(cov_mat[i_range, i_range])
            for j in 1:l
                j_range = (j - 1) * q + 1:j * q
                cov_jj = norm(cov_mat[j_range, j_range])
                cov_ij = norm(cov_mat[i_range, j_range])
                denom = sqrt(cov_ii * cov_jj)
                cor_mat[i, j] = denom > 0 ? clamp(cov_ij / denom, 0.0, 1.0) : 0.0
            end
            cor_mat[i, i] = 1.0
        end

        result = (covariation=cor_mat, positions=collect(1:l))
        return _register_coevolution_result!(_ctx, result, "positional_covariation_matrix"; parents=provenance_parent_ids(alignment), parameters=(metric=metric, n_cols=l))
    elseif metric in (:mi, :mutual_information)
        mi_map = mutual_information_contacts(filtered; pseudocount=pseudocount, min_separation=min_separation, apc=apc, max_gap_fraction=max_gap_fraction, min_sequence_coverage=min_sequence_coverage, identity_threshold=identity_threshold, weights=weights_vec, _ctx=_ctx)
        result = (covariation=mi_map.scores, positions=collect(1:l))
        return _register_coevolution_result!(_ctx, result, "positional_covariation_matrix"; parents=provenance_parent_ids(alignment), parameters=(metric=metric, n_cols=l))
    else
        throw(ArgumentError("Unsupported metric: $metric. Supported metrics are :frobenius, :frobenius_correlation, :mi, :mutual_information"))
    end
end

# ---------------------------------------------------------------------------
# Gap Analysis
# ---------------------------------------------------------------------------

"""
    gap_analysis(alignment)

Analyse gap distribution per column and sequence in a multiple sequence alignment.
"""
function gap_analysis(alignment::MultipleSequenceAlignment;
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    strings = _alignment_strings(alignment)
    n = length(strings)
    l = ncodeunits(strings[1])
    all(ncodeunits(s) == l for s in strings) || throw(ArgumentError("all sequences must have equal length"))

    col_gap = [mean(strings[i][col] == '-' for i in 1:n) for col in 1:l]
    seq_gap = [mean(strings[i][col] == '-' for col in 1:l) for i in 1:n]

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
    result = (
        column_gap_fraction   = col_gap,
        sequence_gap_fraction = seq_gap,
        gap_blocks            = gap_blocks,
        total_gap_fraction    = total_gap
    )

    return _register_coevolution_result!(_ctx, result, "gap_analysis"; parents=provenance_parent_ids(alignment), parameters=(n_seqs=n, n_cols=l, total_gap=total_gap))
end

# ---------------------------------------------------------------------------
# Contact Precision-Recall
# ---------------------------------------------------------------------------

"""
    contact_precision_recall(contact_map, true_contacts; min_separation=5, n_points=50)

Compute Precision-Recall curve statistics and Area Under PR Curve (AUC-PR).
"""
function contact_precision_recall(contact_map::ContactMap, true_contacts::AbstractMatrix{<:Real};
    min_separation::Int=5,
    n_points::Int=50,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

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
        abs(j - i) >= min_separation && true_contacts[i, j] > 0.5 && push!(true_set, (i, j))
    end
    n_true = length(true_set)
    n_true == 0 && return (thresholds=Float64[], precision=Float64[], recall=Float64[], auc_pr=0.0)

    max_s = isempty(pairs) ? 1.0 : maximum(p -> p[3], pairs)
    min_s = isempty(pairs) ? 0.0 : minimum(p -> p[3], pairs)
    if min_s >= max_s
        min_s = max_s - 1.0
    end

    thresholds = range(max_s, min_s, length=n_points)
    prec_vec = Float64[]
    rec_vec  = Float64[]

    push!(prec_vec, 1.0)
    push!(rec_vec, 0.0)

    for thr in thresholds
        pred = [(p[1], p[2]) for p in pairs if p[3] >= thr]
        tp = count(p -> p in true_set, pred)
        prec = !isempty(pred) ? tp / length(pred) : 1.0
        rec = tp / n_true
        push!(prec_vec, prec)
        push!(rec_vec, rec)
    end

    auc = 0.0
    for k in 1:(length(rec_vec) - 1)
        dr = rec_vec[k+1] - rec_vec[k]
        if dr > 0
            auc += 0.5 * (prec_vec[k] + prec_vec[k+1]) * dr
        end
    end

    result = (thresholds=collect(thresholds), precision=prec_vec, recall=rec_vec, auc_pr=auc)
    return _register_coevolution_result!(_ctx, result, "contact_precision_recall"; parents=provenance_parent_ids(contact_map), parameters=(n_points=n_points, n_true=n_true, auc_pr=auc))
end

# ---------------------------------------------------------------------------
# Alignment Quality Report
# ---------------------------------------------------------------------------

"""
    alignment_quality_report(alignment)

Generate a comprehensive quality report for a multiple sequence alignment.
"""
function alignment_quality_report(alignment::MultipleSequenceAlignment;
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    strings = _alignment_strings(alignment)
    n = length(strings)
    l = ncodeunits(strings[1])
    all(ncodeunits(s) == l for s in strings) || throw(ArgumentError("sequences must have equal length"))

    alphabet = _alphabet(strings)
    encoded  = _encode_alignment(strings, alphabet)

    weights  = sequence_reweighting(encoded; identity_threshold=0.8, _ctx=_ctx)
    eff_n    = sum(weights)

    identities = Float64[]
    for i in 1:n-1, j in (i+1):n
        denom = count(strings[i][k] != '-' || strings[j][k] != '-' for k in 1:l)
        if denom > 0
            id = count(encoded[i, k] == encoded[j, k] && strings[i][k] != '-' for k in 1:l) / denom
            push!(identities, id)
        end
    end
    mean_id = isempty(identities) ? 0.0 : mean(identities)
    min_id  = isempty(identities) ? 0.0 : minimum(identities)
    max_id  = isempty(identities) ? 0.0 : maximum(identities)

    col_gap  = [mean(strings[i][col] == '-' for i in 1:n) for col in 1:l]
    cons     = column_conservation_scores(alignment; _ctx=_ctx)
    logo     = sequence_logo_entropy(alignment; _ctx=_ctx)

    result = (
        n_sequences            = n,
        alignment_length       = l,
        effective_sequences    = eff_n,
        mean_pairwise_identity = mean_id,
        min_pairwise_identity  = min_id,
        max_pairwise_identity  = max_id,
        column_gap_fraction     = col_gap,
        mean_gap_fraction       = mean(col_gap),
        conservation_scores     = cons,
        mean_conservation       = mean(cons),
        entropy_per_column      = logo.entropy
    )

    return _register_coevolution_result!(_ctx, result, "alignment_quality_report"; parents=provenance_parent_ids(alignment), parameters=(n_seqs=n, n_cols=l, effective_sequences=eff_n))
end

# ---------------------------------------------------------------------------
# Interactive HTML Visualizations
# ---------------------------------------------------------------------------

"""
    to_html(cmap::ContactMap) -> String

Generate a standalone interactive HTML Canvas report for a ContactMap.
Provides dark mode styling, score threshold sliders, matrix heatmaps, hover tooltips, and top contact pair inspection.
"""
function to_html(cmap::ContactMap)
    l = size(cmap.scores, 1)
    matrix_json = "[" * join(["[" * join([string(round(cmap.scores[i, j]; digits=4)) for j in 1:l], ",") * "]" for i in 1:l], ",") * "]"

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BioToolkit — Co-evolutionary Contact Map</title>
    <style>
        :root { --bg: #0f172a; --panel: #1e293b; --text: #f8fafc; --accent: #38bdf8; --border: #334155; --glow: #818cf8; }
        body { margin: 0; font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); padding: 20px; }
        .header { display: flex; align-items: center; justify-content: space-between; background: var(--panel); padding: 16px 24px; border-radius: 12px; border: 1px solid var(--border); margin-bottom: 20px; box-shadow: 0 4px 20px rgba(0,0,0,0.3); }
        .title { font-size: 1.4rem; font-weight: 700; color: var(--accent); display: flex; align-items: center; gap: 10px; }
        .controls { display: flex; gap: 16px; align-items: center; }
        label { font-size: 0.9rem; color: #94a3b8; font-weight: 500; }
        input[type=range] { accent-color: var(--accent); cursor: pointer; }
        .main-grid { display: grid; grid-template-columns: 1fr 340px; gap: 20px; }
        .card { background: var(--panel); border-radius: 12px; border: 1px solid var(--border); padding: 20px; box-shadow: 0 4px 20px rgba(0,0,0,0.3); position: relative; }
        canvas { width: 100%; display: block; border-radius: 8px; cursor: crosshair; }
        .tooltip { position: absolute; background: rgba(15, 23, 42, 0.95); border: 1px solid var(--accent); padding: 8px 12px; border-radius: 6px; font-size: 12px; pointer-events: none; display: none; z-index: 100; box-shadow: 0 4px 12px rgba(0,0,0,0.5); }
        .contact-list { max-height: 540px; overflow-y: auto; font-family: monospace; font-size: 13px; }
        .contact-row { display: flex; justify-content: space-between; padding: 8px 12px; border-bottom: 1px solid var(--border); border-radius: 4px; }
        .contact-row:hover { background: #334155; }
        .badge { background: #1e1b4b; color: var(--glow); padding: 2px 8px; border-radius: 12px; font-weight: 600; }
    </style>
</head>
<body>
    <div class="header">
        <div class="title">🧬 Residue Contact Map Visualizer <span style="font-size: 0.8rem; color: #94a3b8;">($(l) × $(l) Residues)</span></div>
        <div class="controls">
            <label for="threshold">Score Cutoff: <span id="thresh-val" style="color: var(--accent); font-weight: bold;">0.20</span></label>
            <input type="range" id="threshold" min="0.0" max="1.0" step="0.01" value="0.20" oninput="updatePlot()">
        </div>
    </div>
    <div class="main-grid">
        <div class="card">
            <canvas id="contactCanvas"></canvas>
            <div id="tooltip" class="tooltip"></div>
        </div>
        <div class="card">
            <h3 style="margin-top: 0; color: var(--accent);">Top Co-evolving Pairs</h3>
            <div id="contactList" class="contact-list"></div>
        </div>
    </div>
    <script>
        const scores = $(matrix_json);
        const L = $(l);
        const canvas = document.getElementById('contactCanvas');
        const ctx = canvas.getContext('2d');
        const tooltip = document.getElementById('tooltip');

        function resizeCanvas() {
            const size = Math.min(canvas.parentElement.clientWidth - 40, 600);
            canvas.width = size * window.devicePixelRatio;
            canvas.height = size * window.devicePixelRatio;
            canvas.style.width = size + 'px';
            canvas.style.height = size + 'px';
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
            render();
        }
        window.addEventListener('resize', resizeCanvas);

        function getColor(val, cutoff) {
            if (val < cutoff) return '#1e293b';
            const norm = (val - cutoff) / (1.0 - cutoff + 1e-6);
            const r = Math.round(56 + norm * 199);
            const g = Math.round(189 - norm * 50);
            const b = Math.round(248 - norm * 100);
            return `rgb(\${r},\${g},\${b})`;
        }

        function render() {
            const cutoff = parseFloat(document.getElementById('threshold').value);
            document.getElementById('thresh-val').textContent = cutoff.toFixed(2);
            const displaySize = parseFloat(canvas.style.width);
            const cell = displaySize / L;

            ctx.clearRect(0, 0, displaySize, displaySize);

            for (let i = 0; i < L; i++) {
                for (let j = 0; j < L; j++) {
                    const val = scores[i][j];
                    ctx.fillStyle = getColor(val, cutoff);
                    ctx.fillRect(j * cell, i * cell, cell, cell);
                }
            }
            updateContactList(cutoff);
        }

        function updateContactList(cutoff) {
            const pairs = [];
            for (let i = 0; i < L - 1; i++) {
                for (let j = i + 1; j < L; j++) {
                    if (scores[i][j] >= cutoff) {
                        pairs.push({ i: i + 1, j: j + 1, score: scores[i][j] });
                    }
                }
            }
            pairs.sort((a, b) => b.score - a.score);

            const listEl = document.getElementById('contactList');
            listEl.innerHTML = pairs.slice(0, 30).map(p => `
                <div class="contact-row">
                    <span>Residue \${p.i} — \${p.j}</span>
                    <span class="badge">\${p.score.toFixed(3)}</span>
                </div>
            `).join('');
        }

        canvas.addEventListener('mousemove', (e) => {
            const rect = canvas.getBoundingClientRect();
            const displaySize = parseFloat(canvas.style.width);
            const cell = displaySize / L;
            const x = e.clientX - rect.left;
            const y = e.clientY - rect.top;
            const col = Math.floor(x / cell);
            const row = Math.floor(y / cell);

            if (row >= 0 && row < L && col >= 0 && col < L) {
                tooltip.style.display = 'block';
                tooltip.style.left = (x + 15) + 'px';
                tooltip.style.top = (y + 15) + 'px';
                tooltip.innerHTML = `<strong>Residues (\${row + 1}, \${col + 1})</strong><br>Score: \${scores[row][col].toFixed(4)}`;
            }
        });

        canvas.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });
        function updatePlot() { render(); }
        setTimeout(resizeCanvas, 50);
    </script>
</body>
</html>
"""
end

"""
    to_html(model::PseudoLikelihoodModel) -> String

Generate an interactive HTML inspection report for a PseudoLikelihoodModel object.
"""
function to_html(model::PseudoLikelihoodModel)
    return to_html(ContactMap(model.apc_scores, collect(1:size(model.apc_scores, 1))))
end

"""
    visualize_contact_map_html(cmap::ContactMap) -> String

Alias function for interactive HTML ContactMap report generation.
"""
function visualize_contact_map_html(cmap::ContactMap)
    return to_html(cmap)
end

"""
    visualize_coevolution_network_html(network) -> String

Generate an interactive HTML5 network visualizer for an evolutionary coupling network.
"""
function visualize_coevolution_network_html(network)
    edges_json = "[" * join(["{\"source\":$(network.edges.node_i[k]),\"target\":$(network.edges.node_j[k]),\"weight\":$(round(network.edges.weights[k]; digits=4))}" for k in 1:length(network.edges.node_i)], ",") * "]"
    degree_json = "[" * join([string(d) for d in network.degree], ",") * "]"
    hubs_json = "[" * join([string(h) for h in network.hub_residues], ",") * "]"

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BioToolkit — Evolutionary Coupling Network</title>
    <style>
        body { margin: 0; font-family: system-ui, sans-serif; background: #0f172a; color: #f8fafc; padding: 20px; }
        .card { background: #1e293b; border-radius: 12px; border: 1px solid #334155; padding: 20px; }
        h2 { color: #38bdf8; margin-top: 0; }
        canvas { width: 100%; height: 500px; display: block; border-radius: 8px; }
    </style>
</head>
<body>
    <div class="card">
        <h2>🌐 Evolutionary Coupling Network</h2>
        <canvas id="netCanvas"></canvas>
    </div>
    <script>
        const edges = $(edges_json);
        const degrees = $(degree_json);
        const hubs = new Set($(hubs_json));
        const L = degrees.length;
        const canvas = document.getElementById('netCanvas');
        const ctx = canvas.getContext('2d');

        function draw() {
            canvas.width = canvas.clientWidth * window.devicePixelRatio;
            canvas.height = 500 * window.devicePixelRatio;
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
            const width = canvas.clientWidth;
            const height = 500;

            ctx.clearRect(0, 0, width, height);
            const cx = width / 2;
            const cy = height / 2;
            const radius = Math.min(width, height) * 0.38;

            const nodes = [];
            for (let i = 0; i < L; i++) {
                const angle = (2 * Math.PI * i) / L - Math.PI / 2;
                nodes.push({ x: cx + radius * Math.cos(angle), y: cy + radius * Math.sin(angle), id: i + 1 });
            }

            // Draw edges
            edges.forEach(e => {
                const n1 = nodes[e.source - 1];
                const n2 = nodes[e.target - 1];
                ctx.beginPath();
                ctx.moveTo(n1.x, n1.y);
                ctx.lineTo(n2.x, n2.y);
                ctx.strokeStyle = `rgba(56, 189, 248, \${Math.min(1.0, e.weight * 1.5)})`;
                ctx.lineWidth = Math.max(1, e.weight * 3);
                ctx.stroke();
            });

            // Draw nodes
            nodes.forEach((n, i) => {
                const isHub = hubs.has(n.id);
                ctx.beginPath();
                ctx.arc(n.x, n.y, isHub ? 8 : 5, 0, 2 * Math.PI);
                ctx.fillStyle = isHub ? '#818cf8' : '#38bdf8';
                ctx.fill();
                ctx.strokeStyle = '#0f172a';
                ctx.lineWidth = 1.5;
                ctx.stroke();
            });
        }
        setTimeout(draw, 50);
    </script>
</body>
</html>
"""
end

"""
    visualize_alignment_quality_html(report) -> String

Generate an interactive HTML quality report for a multiple sequence alignment.
"""
function visualize_alignment_quality_html(report)
    cons_json = "[" * join([string(round(c; digits=4)) for c in report.conservation_scores], ",") * "]"
    gaps_json = "[" * join([string(round(g; digits=4)) for g in report.column_gap_fraction], ",") * "]"

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>BioToolkit — Alignment Quality Report</title>
    <style>
        body { font-family: system-ui, sans-serif; background: #0f172a; color: #f8fafc; padding: 20px; }
        .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 16px; margin-bottom: 20px; }
        .metric { background: #1e293b; border: 1px solid #334155; padding: 16px; border-radius: 8px; text-align: center; }
        .val { font-size: 1.6rem; font-weight: bold; color: #38bdf8; }
        .label { font-size: 0.85rem; color: #94a3b8; }
    </style>
</head>
<body>
    <h2>📊 MSA Quality Dashboard</h2>
    <div class="grid">
        <div class="metric"><div class="val">$(report.n_sequences)</div><div class="label">Sequences</div></div>
        <div class="metric"><div class="val">$(report.alignment_length)</div><div class="label">Alignment Length</div></div>
        <div class="metric"><div class="val">$(round(report.effective_sequences; digits=1))</div><div class="label">Effective Seqs (Neff)</div></div>
        <div class="metric"><div class="val">$(round(report.mean_pairwise_identity * 100; digits=1))%</div><div class="label">Mean Identity</div></div>
    </div>
</body>
</html>
"""
end

"""
    visualize_sequence_logo_html(logo) -> String

Generate an interactive HTML sequence logo viewer.
"""
function visualize_sequence_logo_html(logo)
    return """
<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><title>BioToolkit — Sequence Logo</title></head>
<body style="background: #0f172a; color: #f8fafc; font-family: system-ui, sans-serif; padding: 20px;">
    <h2>🔤 Sequence Logo Entropy View</h2>
    <p>Alignment Position Count: $(length(logo.position))</p>
</body>
</html>
"""
end

# ---------------------------------------------------------------------------
# Helper Functions: PDB Contact Matrix & Distance Filtering
# ---------------------------------------------------------------------------

"""
    structure_to_contact_matrix(structure::Structure; cutoff::Real=8.0, atom_type::String="CA")

Extract heavy atom (default "CA") 3D coordinates per residue from a PDB `Structure` object and build an L x L symmetric binary contact matrix where `M[i,j] = 1.0` if inter-residue 3D distance <= `cutoff`, else `0.0`.
"""
function structure_to_contact_matrix(structure::Structure;
    cutoff::Real=8.0,
    atom_type::String="CA")

    residues = Residue[]
    for model in structure.models
        for chain in model.chains
            for res in chain.residues
                push!(residues, res)
            end
        end
    end

    l = length(residues)
    coords = Vector{Union{Nothing, Vector{Float64}}}(nothing, l)

    for (i, res) in enumerate(residues)
        target_atom = nothing
        for atom in res.atoms
            if atom.name == atom_type || (atom_type == "CA" && target_atom === nothing)
                target_atom = atom
                if atom.name == atom_type
                    break
                end
            end
        end
        if target_atom !== nothing
            coords[i] = [target_atom.x, target_atom.y, target_atom.z]
        end
    end

    matrix = zeros(Float64, l, l)
    cutoff_val = Float64(cutoff)

    for i in 1:l-1
        c_i = coords[i]
        c_i === nothing && continue
        for j in i+1:l
            c_j = coords[j]
            c_j === nothing && continue
            dx = c_i[1] - c_j[1]
            dy = c_i[2] - c_j[2]
            dz = c_i[3] - c_j[3]
            dist = sqrt(dx*dx + dy*dy + dz*dz)
            if dist <= cutoff_val
                matrix[i, j] = 1.0
                matrix[j, i] = 1.0
            end
        end
    end

    return matrix
end

"""
    filter_contacts_by_sequence_distance(contact_map::ContactMap; min_separation::Int=5)

Filter contact map scores to retain only long-range coevolutionary signals where |i - j| >= min_separation.
Returns a new ContactMap with short-range contacts zeroed out.
"""
function filter_contacts_by_sequence_distance(contact_map::ContactMap; min_separation::Int=5)
    scores = copy(contact_map.scores)
    l = size(scores, 1)
    for i in 1:l
        lo = max(1, i - min_separation + 1)
        hi = min(l, i + min_separation - 1)
        for j in lo:hi
            scores[i, j] = 0.0
        end
    end
    return ContactMap(scores, copy(contact_map.residue_ids))
end

end
