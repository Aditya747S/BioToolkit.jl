# ==============================================================================
# hmm.jl — Hidden Markov Model algorithms
#
# Provides: Viterbi, Forward, Backward (original), plus:
#   - Baum-Welch EM training
#   - Posterior state decoding
#   - Profile HMMs (HMMER-style: match/insert/delete states)
#   - Pair-HMMs (for pairwise alignment)
#   - HMM-based sequence segmentation
#   - Scaled forward/backward for long sequences
#   - Viterbi training (hard-assignment EM)
#   - HMM I/O (save/load)
#
# References:
#   - Rabiner (1989) Proc IEEE 77(2):257-286
#   - Durbin et al. (1998) Biological Sequence Analysis (Cambridge UP)
#   - Eddy (1998) Bioinformatics 14(9):755-763 (profile HMMs)
#   - Eddy (2011) PLoS Comput Biol 7(10):e1002195 (HMMER3)
# ==============================================================================

export HMM, viterbi, forward, backward
export baum_welch!, baum_welch_train
export posterior_decode, posterior_state_probabilities
export HMMSegment, segment_sequence
export ProfileHMM, ProfileHMMNode, build_profile_hmm, score_profile_hmm
export PairHMM, pair_hmm_align, pair_hmm_score
export viterbi_train!, hmm_log_likelihood
export save_hmm, load_hmm
export baum_welch_from_se  # SummarizedExperiment convenience

using Statistics
using LinearAlgebra
using Random

# Biotype imports for type-safe dispatch
using .BioToolkit: BioSequence, DNAAlphabet, RNAAlphabet, AminoAcidAlphabet,
                   DNASeq, RNASeq, AASeq, SummarizedExperiment, assay
using .BioToolkit: ProvenanceContext, provenance_result!, provenance_parent_ids, register_provenance!, ThreadSafeProvenanceContext, new_provenance_id, ProvenanceParams, active_provenance_context

@inline function _register_hmm_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

# ===========================================================================
# Core HMM type (unchanged from original)
# ===========================================================================

"""
    HMM{T<:AbstractFloat}

A Hidden Markov Model operating in log-space for numerical stability.
"""
struct HMM{T<:AbstractFloat}
    states::Vector{String}
    alphabet::Vector{UInt8}
    initial::Vector{T}
    transitions::Matrix{T}     # [from, to] — log-space
    emissions::Matrix{T}       # [state, alphabet_idx] — log-space
    _emission_lookup::Vector{Int}
    _unknown_emission::T

    function HMM{T}(states, alphabet, initial, transitions, emissions, unknown_prob=T(-Inf)) where {T}
        length(states) == length(initial) || throw(ArgumentError("initial must match states"))
        size(transitions) == (length(states), length(states)) || throw(ArgumentError("transitions must be N x N"))
        size(emissions) == (length(states), length(alphabet)) || throw(ArgumentError("emissions must be N x V"))
        lookup = fill(0, 256)
        for (i, b) in enumerate(alphabet); lookup[Int(b)+1] = i; end
        new{T}(states, alphabet, initial, transitions, emissions, lookup, unknown_prob)
    end
end

function HMM(states::Vector{String}, alphabet::Vector{UInt8}, initial::AbstractVector,
             transitions::AbstractMatrix, emissions::AbstractMatrix; log_space::Bool=false)
    T = Float64
    ini = T.(log_space ? initial : log.(max.(initial, eps(T))))
    trn = Matrix{T}(log_space ? transitions : log.(max.(transitions, eps(T))))
    ems = Matrix{T}(log_space ? emissions  : log.(max.(emissions,  eps(T))))
    return HMM{T}(states, alphabet, ini, trn, ems)
end

@inline function _get_emission_logprob(hmm::HMM{T}, state_idx::Int, byte::UInt8) where {T}
    col = hmm._emission_lookup[Int(byte)+1]
    col == 0 && return hmm._unknown_emission
    return hmm.emissions[state_idx, col]
end

@inline function logaddexp(x::T, y::T) where {T<:AbstractFloat}
    x == -Inf && return y; y == -Inf && return x
    m = max(x, y); return m + log1p(exp(min(x,y) - m))
end

# ===========================================================================
# Viterbi (original; unchanged)
# ===========================================================================

"""
    viterbi(hmm, sequence) → (path, log_prob)

Compute the most likely hidden state path using the Viterbi algorithm.
"""
function viterbi(hmm::HMM{T}, sequence::AbstractVector{UInt8}) where {T}
    L = length(sequence); N = length(hmm.states)
    _ctx = active_provenance_context()
    L == 0 && return _register_hmm_result!(_ctx, (Int[], T(-Inf)), "viterbi"; parameters=(seq_length=0, n_states=N))
    v = Matrix{T}(undef, N, L)
    tb = Matrix{Int}(undef, N, L)
    @inbounds for i in 1:N
        v[i,1] = hmm.initial[i] + _get_emission_logprob(hmm, i, sequence[1])
        tb[i,1] = 0
    end
    @inbounds for t in 2:L
        b = sequence[t]
        for j in 1:N
            best_p, best_i = T(-Inf), 0
            for i in 1:N
                p = v[i,t-1] + hmm.transitions[i,j]
                p > best_p && (best_p = p; best_i = i)
            end
            v[j,t] = best_p + _get_emission_logprob(hmm, j, b)
            tb[j,t] = best_i
        end
    end
    best_p, best_last = T(-Inf), 0
    @inbounds for i in 1:N; v[i,L] > best_p && (best_p = v[i,L]; best_last = i); end
    path = Vector{Int}(undef, L)
    if best_last != 0
        curr = best_last
        @inbounds for t in L:-1:1; path[t] = curr; curr = tb[curr,t]; end
    end
    result = (path, best_p)
    return _register_hmm_result!(_ctx, result, "viterbi"; parameters=(seq_length=L, n_states=N, log_prob=Float64(best_p)))
end

viterbi(hmm::HMM, seq::BioSequence; _ctx=nothing)     = viterbi(hmm, seq.data)
viterbi(hmm::HMM, seq::String; _ctx=nothing)          = viterbi(hmm, collect(codeunits(String(seq))))

# ===========================================================================
# Forward (original; log-space, rolling buffer)
# ===========================================================================

"""
    forward(hmm, sequence) → log_probability
"""
function forward(hmm::HMM{T}, sequence::AbstractVector{UInt8}) where {T}
    L = length(sequence); N = length(hmm.states)
    _ctx = active_provenance_context()
    L == 0 && return _register_hmm_result!(_ctx, T(-Inf), "forward"; parameters=(seq_length=0, n_states=N))
    alpha = Matrix{T}(undef, N, 2); cc, pc = 1, 2
    @inbounds for i in 1:N
        alpha[i,cc] = hmm.initial[i] + _get_emission_logprob(hmm, i, sequence[1])
    end
    @inbounds for t in 2:L
        cc, pc = pc, cc; b = sequence[t]
        for j in 1:N
            s = T(-Inf)
            for i in 1:N; s = logaddexp(s, alpha[i,pc] + hmm.transitions[i,j]); end
            alpha[j,cc] = s + _get_emission_logprob(hmm, j, b)
        end
    end
    total = T(-Inf)
    @inbounds for i in 1:N; total = logaddexp(total, alpha[i,cc]); end
    return _register_hmm_result!(_ctx, total, "forward"; parameters=(seq_length=L, n_states=N, log_likelihood=Float64(total)))
end

forward(hmm::HMM, seq::BioSequence; _ctx=nothing)    = forward(hmm, seq.data)
forward(hmm::HMM, seq::String; _ctx=nothing)         = forward(hmm, collect(codeunits(String(seq))))

# ===========================================================================
# Backward (original)
# ===========================================================================

"""
    backward(hmm, sequence) → N × L log-probability matrix
"""
function backward(hmm::HMM{T}, sequence::AbstractVector{UInt8}) where {T}
    L = length(sequence); N = length(hmm.states)
    beta = Matrix{T}(undef, N, L)
    L == 0 && return beta
    @inbounds for i in 1:N; beta[i,L] = T(0); end
    @inbounds for t in L-1:-1:1
        nb = sequence[t+1]
        for i in 1:N
            s = T(-Inf)
            for j in 1:N
                s = logaddexp(s, hmm.transitions[i,j] + _get_emission_logprob(hmm, j, nb) + beta[j,t+1])
            end
            beta[i,t] = s
        end
    end
    return beta
end

backward(hmm::HMM, seq::BioSequence)    = backward(hmm, seq.data)
backward(hmm::HMM, seq::String)        = backward(hmm, collect(codeunits(String(seq))))
# Removed backward(::HMM, ::AbstractString) - use BioSequence instead

# ===========================================================================
# Full forward table (N × L) — needed for Baum-Welch
# ===========================================================================

function _forward_table(hmm::HMM{T}, sequence::AbstractVector{UInt8}) where {T}
    L = length(sequence); N = length(hmm.states)
    alpha = fill(T(-Inf), N, L)
    @inbounds for i in 1:N
        alpha[i,1] = hmm.initial[i] + _get_emission_logprob(hmm, i, sequence[1])
    end
    @inbounds for t in 2:L, j in 1:N
        s = T(-Inf)
        for i in 1:N; s = logaddexp(s, alpha[i,t-1] + hmm.transitions[i,j]); end
        alpha[j,t] = s + _get_emission_logprob(hmm, j, sequence[t])
    end
    return alpha
end

# ===========================================================================
# NEW: HMM log-likelihood
# ===========================================================================

"""
    hmm_log_likelihood(hmm, sequences) → Float64

Compute the total log-likelihood of a set of observed sequences given `hmm`.
"""
function hmm_log_likelihood(hmm::HMM, sequences)
    total = 0.0
    for seq in sequences
        bytes = seq isa AbstractVector{UInt8} ? seq :
                seq isa BioSequence ? seq.data : collect(codeunits(String(seq)))
        total += forward(hmm, bytes)
    end

    return total
end

# ===========================================================================
# NEW: Posterior state probabilities (smoothed decoding)
# ===========================================================================

"""
    posterior_state_probabilities(hmm, sequence) → N × L Matrix{Float64}

Compute γ_t(i) = P(state=i at t | sequence, hmm) using forward-backward.
Equivalent to the "E-step" responsibilities without accumulating.
"""
function posterior_state_probabilities(hmm::HMM{T}, sequence::AbstractVector{UInt8}) where {T}
    L = length(sequence); N = length(hmm.states)
    alpha = _forward_table(hmm, sequence)
    beta  = backward(hmm, sequence)
    log_px = forward(hmm, sequence)
    gamma = Matrix{Float64}(undef, N, L)
    @inbounds for t in 1:L
        s = -Inf
        for i in 1:N; s = logaddexp(s, alpha[i,t] + beta[i,t]); end
        for i in 1:N
            gamma[i,t] = exp(alpha[i,t] + beta[i,t] - s)
        end
    end
    return gamma
end

# Removed posterior_state_probabilities(::HMM, ::AbstractString) - use BioSequence instead

# ===========================================================================
# NEW: Posterior decode (MAP assignment per position)
# ===========================================================================

"""
    posterior_decode(hmm, sequence) → (assignments::Vector{Int}, gamma::Matrix{Float64})

Assign each position to the highest-posterior state. Softer than Viterbi;
does not guarantee globally-consistent paths. Analogous to `hmmlearn.decode`.
"""
function posterior_decode(hmm::HMM, sequence)
    bytes = sequence isa AbstractVector{UInt8} ? sequence :
            sequence isa BioSequence           ? sequence.data : collect(codeunits(String(sequence)))
    gamma = posterior_state_probabilities(hmm, bytes)
    assignments = [argmax(gamma[:,t]) for t in axes(gamma,2)]

    return assignments, gamma
end

# ===========================================================================
# NEW: Baum-Welch EM training
# ===========================================================================

"""
    baum_welch!(hmm, sequences; max_iter=100, tol=1e-4, verbose=false) → log_likelihoods

Train `hmm` in-place using the Baum-Welch (EM) algorithm on `sequences`.
Returns the per-iteration log-likelihood history.

Analogous to `pomegranate.HiddenMarkovModel.fit` / `hmmlearn.GaussianHMM.fit`.
"""
function baum_welch!(hmm::HMM{T}, sequences; max_iter::Int=100, tol::Real=1e-4, verbose::Bool=false) where {T}
    N  = length(hmm.states)
    V  = length(hmm.alphabet)
    seqs = [s isa AbstractVector{UInt8} ? s :
            s isa BioSequence ? s.data : collect(codeunits(String(s))) for s in sequences]

    log_likelihoods = Float64[]

    # Working copies (exp space for accumulation, then log for storage)
    ini_acc  = zeros(Float64, N)
    trans_acc = zeros(Float64, N, N)
    emis_acc  = zeros(Float64, N, V)

    for iter in 1:max_iter
        fill!(ini_acc,   0.0)
        fill!(trans_acc, 0.0)
        fill!(emis_acc,  0.0)
        total_ll = 0.0

        for s in seqs
            L   = length(s)
            L < 1 && continue
            alpha = _forward_table(hmm, s)
            beta  = backward(hmm, s)

            # Log-probability of this sequence
            log_px = -Inf
            for i in 1:N; log_px = logaddexp(log_px, alpha[i,L]); end
            total_ll += log_px

            # γ_t(i) accumulation
            for t in 1:L
                denom = -Inf
                for i in 1:N; denom = logaddexp(denom, alpha[i,t]+beta[i,t]); end
                for i in 1:N
                    γ = exp(alpha[i,t] + beta[i,t] - denom)
                    t == 1 && (ini_acc[i] += γ)
                    # Emission accumulation
                    col = hmm._emission_lookup[Int(s[t])+1]
                    col > 0 && (emis_acc[i,col] += γ)
                end
            end

            # ξ_t(i,j) transition accumulation
            for t in 1:(L-1)
                denom = -Inf
                for i in 1:N, j in 1:N
                    denom = logaddexp(denom,
                        alpha[i,t] + hmm.transitions[i,j] +
                        _get_emission_logprob(hmm, j, s[t+1]) + beta[j,t+1])
                end
                for i in 1:N, j in 1:N
                    ξ = exp(alpha[i,t] + hmm.transitions[i,j] +
                            _get_emission_logprob(hmm, j, s[t+1]) + beta[j,t+1] - denom)
                    trans_acc[i,j] += ξ
                end
            end
        end

        push!(log_likelihoods, total_ll)
        verbose && @info "Baum-Welch iter $iter  LL = $(round(total_ll; digits=4))"

        # Convergence check
        if length(log_likelihoods) >= 2 &&
           abs(log_likelihoods[end] - log_likelihoods[end-1]) < Float64(tol)
            break
        end

        # M-step: update parameters in-place using HMM field mutation
        ini_sum  = max(sum(ini_acc),  eps(Float64))
        new_ini  = log.(ini_acc  ./ ini_sum)
        trans_row_sums = max.(vec(sum(trans_acc, dims=2)), eps(Float64))
        new_trans = log.(trans_acc ./ trans_row_sums)
        emis_row_sums  = max.(vec(sum(emis_acc,  dims=2)), eps(Float64))
        new_emis  = log.(emis_acc  ./ emis_row_sums)

        # Mutate the HMM fields (they're Vector/Matrix so in-place assignment works)
        hmm.initial    .= new_ini
        hmm.transitions .= new_trans
        hmm.emissions   .= new_emis
    end
    return log_likelihoods
end

"""
    baum_welch_train(n_states, alphabet, sequences; kwargs...) → (hmm, log_likelihoods)

Randomly initialise and Baum-Welch train a new HMM from scratch.
"""
function baum_welch_train(n_states::Int, alphabet::Vector{UInt8}, sequences;
                          max_iter::Int=100, tol::Real=1e-4, seed::Int=1, verbose::Bool=false)
    rng   = Random.MersenneTwister(seed)
    N, V  = n_states, length(alphabet)
    ini   = normalize(rand(rng, Float64, N), 1)
    trans = Matrix{Float64}(undef, N, N)
    for i in 1:N; trans[i,:] = normalize(rand(rng, Float64, N), 1); end
    emis  = Matrix{Float64}(undef, N, V)
    for i in 1:N; emis[i,:]  = normalize(rand(rng, Float64, V), 1); end
    hmm   = HMM(["state_$i" for i in 1:n_states], alphabet, ini, trans, emis)
    lls   = baum_welch!(hmm, sequences; max_iter=max_iter, tol=tol, verbose=verbose)
    _ctx = active_provenance_context()
    _register_hmm_result!(_ctx, hmm, "baum_welch_train"; parameters=(n_states=n_states, n_symbols=V, max_iter=max_iter, n_iterations=length(lls)))


    return hmm, lls
end

# ===========================================================================
# NEW: Viterbi training (hard EM — faster approximation of Baum-Welch)
# ===========================================================================

"""
    viterbi_train!(hmm, sequences; max_iter=50, tol=1e-4)

Hard-assignment EM: decode the Viterbi path, then re-estimate parameters
from assignments. Faster but less accurate than Baum-Welch.
"""
function viterbi_train!(hmm::HMM{T}, sequences; max_iter::Int=50, tol::Real=1e-4) where {T}
    N = length(hmm.states); V = length(hmm.alphabet)
    seqs = [s isa AbstractVector{UInt8} ? s :
            s isa BioSequence ? s.data : collect(codeunits(String(s))) for s in sequences]

    prev_ll = -Inf
    for iter in 1:max_iter
        ini_acc   = zeros(Float64, N)
        trans_acc = zeros(Float64, N, N)
        emis_acc  = zeros(Float64, N, V)

        for s in seqs
            path, _ = viterbi(hmm, s)
            isempty(path) && continue
            ini_acc[path[1]] += 1.0
            for t in 1:length(s)
                col = hmm._emission_lookup[Int(s[t])+1]
                col > 0 && (emis_acc[path[t], col] += 1.0)
                t < length(path) && (trans_acc[path[t], path[t+1]] += 1.0)
            end
        end

        ini_sum  = max(sum(ini_acc),  eps())
        hmm.initial    .= log.(ini_acc  ./ ini_sum)
        for i in 1:N
            rs = max(sum(trans_acc[i,:]), eps())
            hmm.transitions[i,:] .= log.(trans_acc[i,:] ./ rs)
        end
        for i in 1:N
            rs = max(sum(emis_acc[i,:]), eps())
            hmm.emissions[i,:]   .= log.(emis_acc[i,:]  ./ rs)
        end

        ll = hmm_log_likelihood(hmm, seqs)
        abs(ll - prev_ll) < Float64(tol) && break
        prev_ll = ll
    end
    return hmm
end

# ===========================================================================
# NEW: HMM sequence segmentation
# ===========================================================================

"""
    HMMSegment

A contiguous segment of sequence assigned to a single HMM state.
"""
struct HMMSegment
    start::Int
    stop::Int
    state_index::Int
    state_name::String
    log_probability::Float64
end
# ===========================================================================
# NEW: Profile HMMs (HMMER-style)
# ===========================================================================

"""
    ProfileHMMNode

One column of a profile HMM: match, insert, and delete state parameters.
"""
struct ProfileHMMNode
    position::Int
    match_emissions::Vector{Float64}    # log-prob over alphabet; length = V
    insert_emissions::Vector{Float64}
    match_to_match::Float64             # log-prob transitions
    match_to_insert::Float64
    match_to_delete::Float64
    insert_to_match::Float64
    insert_to_insert::Float64
    delete_to_match::Float64
    delete_to_delete::Float64
end

"""
    ProfileHMM

HMMER-style profile HMM for a conserved sequence family.
"""
struct ProfileHMM
    name::String
    length::Int               # number of match columns
    alphabet::Vector{UInt8}
    nodes::Vector{ProfileHMMNode}
    null_log_odds::Float64    # log-odds null model correction
end

# ---------------------------------------------------------------------------
# SummarizedExperiment convenience wrapper
# ---------------------------------------------------------------------------

"""
    baum_welch_from_se(se, n_states; assay_name="counts", kwargs...) → (HMM, log_likelihoods)

Train an HMM on sequences extracted from a `SummarizedExperiment`.
Expects the assay matrix to encode discretised observations (one column = one sequence).
This provides a biotype-native entry point for HMM training on genomic data.
"""
function baum_welch_from_se(
    se::SummarizedExperiment,
    n_states::Int;
    assay_name::String="counts",
    kwargs...)
    X = assay(se, assay_name)
    # Treat each column as a discrete observation sequence (quantise if needed)
    max_val = max(maximum(X), 1)
    alphabet = UInt8.(1:min(ceil(Int, max_val), 255))
    sequences = [UInt8.(clamp.(round.(Int, X[:, j]), 1, length(alphabet))) for j in axes(X, 2)]
    return baum_welch_train(n_states, alphabet, sequences; kwargs...)
end

# ---------------------------------------------------------------------------
# AASeq / RNASeq dispatch for profile HMM scoring
# ---------------------------------------------------------------------------

"""
    score_profile_hmm(phmm, sequence::AASeq) → Float64

Type-safe overload for scoring a protein sequence against a profile HMM.
"""
score_profile_hmm(phmm::ProfileHMM, seq::BioSequence{AminoAcidAlphabet}) =
    score_profile_hmm(phmm, String(seq))

"""
    score_profile_hmm(phmm, sequence::DNASeq) → Float64

Type-safe overload for scoring a DNA sequence against a profile HMM.
"""
score_profile_hmm(phmm::ProfileHMM, seq::BioSequence{DNAAlphabet}) =
    score_profile_hmm(phmm, String(seq))

# Segment sequence with BioSequence input
segment_sequence(hmm::HMM, seq::BioSequence; kwargs...) =
    segment_sequence(hmm, seq.data; kwargs...)

posterior_decode(hmm::HMM, seq::BioSequence; kwargs...) =
    posterior_decode(hmm, seq.data; kwargs...)

posterior_state_probabilities(hmm::HMM, seq::BioSequence; kwargs...) =
    posterior_state_probabilities(hmm, seq.data; kwargs...)

"""
    segment_sequence(hmm, sequence; method=:viterbi) → Vector{HMMSegment}

Segment a sequence into contiguous runs of the same hidden state.
`method` can be `:viterbi` or `:posterior`.

Analogous to `hmmlearn` `predict` + run-length encoding.
"""
function segment_sequence(hmm::HMM, sequence; method::Symbol=:viterbi)
    bytes = sequence isa AbstractVector{UInt8} ? sequence :
            sequence isa BioSequence           ? sequence.data : collect(codeunits(String(sequence)))

    if method == :posterior
        assignments, gamma = posterior_decode(hmm, bytes)
        probs = [gamma[assignments[t], t] for t in eachindex(assignments)]
    else
        assignments, log_p = viterbi(hmm, bytes)
        probs = fill(exp(log_p / max(length(bytes), 1)), length(bytes))
    end

    segments = HMMSegment[]
    isempty(assignments) && return segments

    seg_start = 1; prev_state = assignments[1]
    for t in 2:length(assignments)
        if assignments[t] != prev_state
            push!(segments, HMMSegment(seg_start, t-1, prev_state,
                hmm.states[prev_state], log(max(mean(probs[seg_start:t-1]), eps()))))
            seg_start = t; prev_state = assignments[t]
        end
    end
    push!(segments, HMMSegment(seg_start, length(assignments), prev_state,
        hmm.states[prev_state], log(max(mean(probs[seg_start:end]), eps()))))

    return segments
end

"""
    build_profile_hmm(msa; alphabet=nothing, pseudocount=1.0, name="profile") → ProfileHMM

Build a position-specific profile HMM from a multiple sequence alignment.
Uses Dirichlet pseudocounts for probability estimation.

Analogous to `hmmbuild` / `Bio.HMM.Profile`.
"""
function build_profile_hmm(msa; alphabet::Union{Nothing,Vector{UInt8}}=nothing,
                            pseudocount::Real=1.0, name::String="profile")
    seqs = msa isa Vector{String} ? msa :
           [String(s) for s in (hasproperty(msa, :sequences) ? msa.sequences : msa)]
    isempty(seqs) && throw(ArgumentError("MSA must contain at least one sequence"))
    L   = length(seqs[1])
    all(s -> length(s) == L, seqs) ||
        throw(ArgumentError("All MSA sequences must be the same length"))

    if alphabet === nothing
        # Amino acid alphabet
        alphabet = UInt8.(collect(b"ACDEFGHIKLMNPQRSTVWY"))
    end
    V  = length(alphabet)
    al = Dict(b => i for (i,b) in enumerate(alphabet))
    pc = Float64(pseudocount)

    nodes = ProfileHMMNode[]
    for pos in 1:L
        col = [uppercase(s[pos]) for s in seqs]
        gap_frac = count(c -> c == '-', col) / length(col)

        match_cnts  = fill(pc, V)
        insert_cnts = fill(pc, V)
        for c in col
            b = UInt8(c)
            get(al, b, 0) > 0 && (match_cnts[al[b]] += 1.0)
        end
        ms = sum(match_cnts)
        log_match  = log.(match_cnts ./ ms)
        log_insert = log.(insert_cnts ./ sum(insert_cnts))

        # Simple transition estimates (gap-adjusted)
        m2m = log(max(1.0 - gap_frac - 0.03, 0.01))
        m2i = log(max(0.03, 0.001))
        m2d = log(max(gap_frac, 0.01))
        i2m = log(0.7); i2i = log(0.3)
        d2m = log(0.8); d2d = log(0.2)

        push!(nodes, ProfileHMMNode(pos, log_match, log_insert,
            m2m, m2i, m2d, i2m, i2i, d2m, d2d))
    end

    # Null model: uniform log-odds in length units
    null_lo = -L * log(V)
    return ProfileHMM(name, L, alphabet, nodes, null_lo)
end

"""
    score_profile_hmm(phmm, sequence; return_per_position=false) → Float64

Score a sequence against a profile HMM using the Viterbi algorithm over
match/insert/delete states. Returns log-odds score (bits if divided by log(2)).

Analogous to `hmmsearch` / `hmmscan` score.
"""
function score_profile_hmm(phmm::ProfileHMM, sequence::BioSequence; return_per_position::Bool=false)
    seq = String(sequence) |> uppercase
    n   = length(seq)
    L   = phmm.length
    al  = Dict(b => i for (i,b) in enumerate(phmm.alphabet))

    # States: M (match), I (insert), D (delete) for each position + flanking
    # DP table: 3 states × (L+1) positions
    INF = -Inf
    M   = fill(INF, n+1, L+2)
    I   = fill(INF, n+1, L+2)
    D   = fill(INF, n+1, L+2)
    M[1,1] = 0.0     # start in begin state

    al_lookup = fill(0, 256)
    for (b,i) in al; al_lookup[Int(b)+1] = i; end

    for j in 1:L
        nd = phmm.nodes[j]
        for i in 1:n
            ch  = UInt8(seq[i])
            emi = al_lookup[Int(ch)+1]
            e_m = emi > 0 ? nd.match_emissions[emi]  : INF
            e_i = emi > 0 ? nd.insert_emissions[emi] : INF

            # Match(i,j): came from M/I/D at (i-1, j-1)
            prev = max(M[i,j], I[i,j], D[i,j])
            M[i+1,j+1] = prev + nd.match_to_match + e_m

            # Insert(i,j): consume seq char; stay in column j
            M_to_I = M[i,j]  + nd.match_to_insert + e_i
            I_to_I = I[i,j]  + nd.insert_to_insert + e_i
            I[i+1,j+1] = max(M_to_I, I_to_I)

            # Delete(i,j): skip seq char in column j
            M_to_D = M[i,j]  + nd.match_to_delete
            D_to_D = D[i,j]  + nd.delete_to_delete
            D[i+1,j+1] = max(M_to_D, D_to_D)
        end
    end

    best = max(M[n+1,L+1], I[n+1,L+1], D[n+1,L+1])

    return best - phmm.null_log_odds
end

# ===========================================================================
# NEW: Pair-HMM (pairwise sequence alignment)
# ===========================================================================

"""
    PairHMM

A three-state pair-HMM for aligning two sequences: Match (M), Insert (I), Delete (D).
"""
struct PairHMM
    match_score::Matrix{Float64}      # V × V substitution score in log-space
    gap_open::Float64                 # log-probability of opening a gap
    gap_extend::Float64               # log-probability of extending a gap
    alphabet::Vector{UInt8}
end

"""
    pair_hmm_align(phmm, seq1, seq2) → (score, alignment1, alignment2)

Align two sequences using a pair-HMM with affine gap penalties.
Equivalent to Gotoh's algorithm (Smith-Waterman with affine gaps).
"""
function pair_hmm_align(phmm::PairHMM, seq1::BioSequence, seq2::BioSequence)
    s1 = collect(uppercase(String(seq1)))
    s2 = collect(uppercase(String(seq2)))
    m, n = length(s1), length(s2)
    al = Dict(b => i for (i,b) in enumerate(Char.(phmm.alphabet)))

    INF = -1e15
    # M[i,j] = best score aligning s1[1:i] to s2[1:j] ending in match
    M = fill(INF, m+1, n+1)
    Ix = fill(INF, m+1, n+1)  # gap in s2 (delete)
    Iy = fill(INF, m+1, n+1)  # gap in s1 (insert)
    M[1,1] = 0.0

    for i in 1:m+1; Ix[i,1] = (i == 1 ? 0.0 : phmm.gap_open + (i-1)*phmm.gap_extend); end
    for j in 1:n+1; Iy[1,j] = (j == 1 ? 0.0 : phmm.gap_open + (j-1)*phmm.gap_extend); end

    for i in 1:m, j in 1:n
        ai = get(al, s1[i], 0); aj = get(al, s2[j], 0)
        sub = (ai > 0 && aj > 0) ? phmm.match_score[ai, aj] : INF
        M[i+1,j+1]  = sub + max(M[i,j], Ix[i,j], Iy[i,j])
        Ix[i+1,j+1] = max(M[i+1,j] + phmm.gap_open, Ix[i+1,j] + phmm.gap_extend)
        Iy[i+1,j+1] = max(M[i,j+1] + phmm.gap_open, Iy[i,j+1] + phmm.gap_extend)
    end

    score = max(M[m+1,n+1], Ix[m+1,n+1], Iy[m+1,n+1])

    # Traceback
    aln1 = Char[]; aln2 = Char[]
    i, j = m, n
    state = :M
    while i > 0 || j > 0
        if state == :M
            i == 0 || j == 0 && break
            push!(aln1, s1[i]); push!(aln2, s2[j]); i -= 1; j -= 1
        elseif state == :Ix
            push!(aln1, s1[i]); push!(aln2, '-'); i -= 1
            Ix[i+1,j+1] ≈ Ix[i+2,j+1] + phmm.gap_extend && (state = :Ix) || (state = :M)
        else
            push!(aln1, '-'); push!(aln2, s2[j]); j -= 1
            Iy[i+1,j+1] ≈ Iy[i+1,j+2] + phmm.gap_extend && (state = :Iy) || (state = :M)
        end
    end

    return score, String(reverse(aln1)), String(reverse(aln2))
end

"""
    pair_hmm_score(phmm, seq1, seq2) → Float64

Compute the pair-HMM alignment score without returning the alignment.
"""
function pair_hmm_score(phmm::PairHMM, seq1, seq2)
    score, _, _ = pair_hmm_align(phmm, seq1, seq2)

    return score
end

# ===========================================================================
# NEW: HMM I/O (JSON-based save/load)
# ===========================================================================

"""
    save_hmm(path, hmm)

Save an HMM to a JSON file.
"""
function save_hmm(path::AbstractString, hmm::HMM)
    # Lightweight serialisation without JSON dependency
    open(String(path), "w") do io
        println(io, "# BioToolkit HMM v1")
        println(io, "states=", join(hmm.states, ","))
        println(io, "alphabet=", join(string.(Int.(hmm.alphabet)), ","))
        println(io, "initial=", join(hmm.initial, ","))
        for i in axes(hmm.transitions, 1)
            println(io, "trans_$(i)=", join(hmm.transitions[i,:], ","))
        end
        for i in axes(hmm.emissions, 1)
            println(io, "emis_$(i)=", join(hmm.emissions[i,:], ","))
        end
    end

    return path
end

"""
    load_hmm(path) → HMM

Load an HMM previously saved by `save_hmm`.
"""
function load_hmm(path::AbstractString)
    lines = readlines(String(path))
    d = Dict{String,String}()
    for l in lines
        startswith(l, "#") && continue
        kv = split(l, "=", limit=2)
        length(kv) == 2 && (d[strip(kv[1])] = strip(kv[2]))
    end
    states  = split(d["states"], ",")
    alphabet = UInt8.(parse.(Int, split(d["alphabet"], ",")))
    initial = parse.(Float64, split(d["initial"], ","))
    N = length(states)
    V = length(alphabet)
    transitions = Matrix{Float64}(undef, N, N)
    emissions   = Matrix{Float64}(undef, N, V)
    for i in 1:N
        transitions[i,:] = parse.(Float64, split(d["trans_$i"], ","))
        emissions[i,:]   = parse.(Float64, split(d["emis_$i"],  ","))
    end

    return HMM(String.(states), alphabet, initial, transitions, emissions; log_space=true)
end
