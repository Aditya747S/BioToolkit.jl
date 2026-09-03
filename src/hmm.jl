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

export AbstractHMM, HMM, GaussianHMM, ChromHMM, GenericEmissionHMM
export log_emission, n_states, initial_distribution, transition_matrix, states
export viterbi, forward, backward, forward_backward, posteriors, loglikelihood, joint_loglikelihood
export baum_welch!, baum_welch_train
export posterior_decode, posterior_state_probabilities
export HMMSegment, segment_sequence, segment_chromatin
export ProfileHMM, ProfileHMMNode, build_profile_hmm, score_profile_hmm
export PairHMM, pair_hmm_align, pair_hmm_score
export aic, bic, state_entropy, cross_validate_hmm
export fit!, logdensityof
export viterbi_train!, hmm_log_likelihood
export save_hmm, load_hmm, to_blender_payload
export baum_welch_from_se  # SummarizedExperiment convenience
export statdists, nparams

using Statistics
using LinearAlgebra
using Random: Random, MersenneTwister, AbstractRNG, default_rng, randperm

@inline function _register_hmm_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

# ===========================================================================
# Core HMM type (unchanged from original)
# ===========================================================================

@inline function _hmm_logsumexp(values)
    m = maximum(values)
    m == -Inf && return m
    total = 0.0
    @inbounds for value in values
        total += exp(value - m)
    end
    return m + log(total)
end

function _hmm_log_probabilities(values::AbstractVector, name::AbstractString)
    any(x -> !isfinite(x) || x < 0, values) && throw(ArgumentError("$name probabilities must be finite and non-negative"))
    total = sum(values)
    total > 0 || throw(ArgumentError("$name probabilities must have positive mass"))
    result = Vector{Float64}(undef, length(values))
    @inbounds for i in eachindex(values)
        result[i] = values[i] == 0 ? -Inf : log(Float64(values[i]) / total)
    end
    return result
end

function _hmm_normalize_log_vector(values::AbstractVector, name::AbstractString)
    any(x -> isnan(x) || x == Inf, values) && throw(ArgumentError("$name log probabilities cannot contain NaN or +Inf"))
    z = _hmm_logsumexp(values)
    isfinite(z) || throw(ArgumentError("$name probabilities must have positive mass"))
    return Float64.(values) .- z
end

function _hmm_log_matrix(matrix::AbstractMatrix, name::AbstractString)
    result = Matrix{Float64}(undef, size(matrix))
    for i in axes(matrix, 1)
        result[i, :] .= _hmm_log_probabilities(view(matrix, i, :), "$name row $i")
    end
    return result
end

function _hmm_log_matrix_from_logs(matrix::AbstractMatrix, name::AbstractString)
    result = Matrix{Float64}(undef, size(matrix))
    for i in axes(matrix, 1)
        result[i, :] .= _hmm_normalize_log_vector(view(matrix, i, :), "$name row $i")
    end
    return result
end

function _hmm_update_log_vector!(destination::AbstractVector, counts::AbstractVector)
    total = sum(counts)
    total <= 0 && return destination
    @inbounds for i in eachindex(destination, counts)
        destination[i] = counts[i] == 0 ? -Inf : log(counts[i] / total)
    end
    return destination
end

function _hmm_update_log_rows!(destination::AbstractMatrix, counts::AbstractMatrix)
    for i in axes(destination, 1)
        _hmm_update_log_vector!(view(destination, i, :), view(counts, i, :))
    end
    return destination
end

"""
    AbstractHMM{T<:AbstractFloat}

Abstract supertype for all Hidden Markov Model implementations in BioToolkit.
Subtypes must implement:
  - `states(hmm)` → Vector{String}
  - `initial_distribution(hmm)` → Vector{T} (log-space)
  - `transition_matrix(hmm)` → Matrix{T} (log-space N x N)
  - `log_emission(hmm, state_idx, observation)` → T (log-probability)
"""
abstract type AbstractHMM{T<:AbstractFloat} end

n_states(hmm::AbstractHMM) = length(states(hmm))
states(hmm::AbstractHMM) = hmm.states
initial_distribution(hmm::AbstractHMM) = hmm.initial
transition_matrix(hmm::AbstractHMM) = hmm.transitions

"""
    GenericEmissionHMM{T, E}

An HMM supporting arbitrary user-provided emission functions or distributions per state.
`emissions` is a vector of callables where `emissions[i](obs)` returns the log-probability of `obs`.
"""
struct GenericEmissionHMM{T<:AbstractFloat, E} <: AbstractHMM{T}
    states::Vector{String}
    initial::Vector{T}        # Log-space
    transitions::Matrix{T}    # Log-space N x N
    emissions::Vector{E}      # State emission functions / distributions
end

function GenericEmissionHMM(states::AbstractVector{<:AbstractString}, initial::AbstractVector,
                            transitions::AbstractMatrix, emissions::AbstractVector; log_space::Bool=false)
    T = Float64
    ini = log_space ? _hmm_normalize_log_vector(initial, "initial") : _hmm_log_probabilities(initial, "initial")
    trn = log_space ? _hmm_log_matrix_from_logs(transitions, "transition") : _hmm_log_matrix(transitions, "transition")
    return GenericEmissionHMM{T, eltype(emissions)}(collect(states), ini, trn, collect(emissions))
end

@inline log_emission(gem::GenericEmissionHMM, state_idx::Int, obs) = gem.emissions[state_idx](obs)

"""
    HMM{T<:AbstractFloat}

A discrete Hidden Markov Model operating in log-space for numerical stability.
"""
struct HMM{T<:AbstractFloat} <: AbstractHMM{T}
    states::Vector{String}
    alphabet::Vector{UInt8}
    initial::Vector{T}
    transitions::Matrix{T}     # [from, to] — log-space
    emissions::Matrix{T}       # [state, alphabet_idx] — log-space
    _emission_lookup::Vector{Int}
    _unknown_emission::T

    function HMM{T}(states, alphabet, initial, transitions, emissions, unknown_prob=T(-Inf)) where {T}
        isempty(states) && throw(ArgumentError("HMM must contain at least one state"))
        isempty(alphabet) && throw(ArgumentError("HMM alphabet must not be empty"))
        length(unique(states)) == length(states) || throw(ArgumentError("state names must be unique"))
        length(unique(alphabet)) == length(alphabet) || throw(ArgumentError("alphabet symbols must be unique"))
        length(states) == length(initial) || throw(ArgumentError("initial must match states"))
        size(transitions) == (length(states), length(states)) || throw(ArgumentError("transitions must be N x N"))
        size(emissions) == (length(states), length(alphabet)) || throw(ArgumentError("emissions must be N x V"))
        (any(isnan, initial) || any(isnan, transitions) || any(isnan, emissions)) && throw(ArgumentError("HMM log probabilities cannot be NaN"))
        lookup = fill(0, 256)
        for (i, b) in enumerate(alphabet); lookup[Int(b)+1] = i; end
        new{T}(states, alphabet, initial, transitions, emissions, lookup, unknown_prob)
    end
end

function HMM(states::AbstractVector{<:AbstractString}, alphabet::AbstractVector{UInt8}, initial::AbstractVector,
             transitions::AbstractMatrix, emissions::AbstractMatrix; log_space::Bool=false)
    T = Float64
    if log_space
        ini = _hmm_normalize_log_vector(initial, "initial")
        trn = _hmm_log_matrix_from_logs(transitions, "transition")
        ems = _hmm_log_matrix_from_logs(emissions, "emission")
    else
        ini = _hmm_log_probabilities(initial, "initial")
        trn = _hmm_log_matrix(transitions, "transition")
        ems = _hmm_log_matrix(emissions, "emission")
    end
    return HMM{T}(String.(states), collect(alphabet), ini, trn, ems)
end

Base.size(hmm::HMM) = (length(hmm.states), length(hmm.alphabet))
Base.size(hmm::HMM, dim::Integer) = size(hmm)[dim]

Base.copy(hmm::HMM{T}) where {T} = HMM{T}(
    copy(hmm.states),
    copy(hmm.alphabet),
    copy(hmm.initial),
    copy(hmm.transitions),
    copy(hmm.emissions),
    hmm._unknown_emission
)

"""
    statdists(hmm) → Vector{Vector{Float64}}

Return the stationary state distribution(s) of `hmm`.
"""
function statdists(hmm::HMM)
    A_prob = exp.(hmm.transitions)
    eig = eigen(collect(transpose(A_prob)))
    dists = Vector{Float64}[]
    for (i, val) in enumerate(eig.values)
        if isapprox(real(val), 1.0; atol=1e-5) && isapprox(imag(val), 0.0; atol=1e-5)
            vec_r = real.(eig.vectors[:, i])
            total = sum(vec_r)
            total > 0 && push!(dists, vec_r ./ total)
        end
    end
    return dists
end

"""
    nparams(hmm) → Int

Return the number of free parameters in `hmm`.
"""
function nparams(hmm::HMM)
    n_states = length(hmm.states)
    n_alpha  = length(hmm.alphabet)
    return (n_states - 1) + n_states * (n_states - 1) + n_states * (n_alpha - 1)
end

function _hmm_sample_categorical(rng::AbstractRNG, log_probs::AbstractVector{T}) where {T}
    m = maximum(log_probs)
    probs = exp.(log_probs .- m)
    total = sum(probs)
    total > 0 || throw(ArgumentError("cannot sample from zero probability mass distribution"))
    u = rand(rng, T) * total
    acc = zero(T)
    for (i, p) in enumerate(probs)
        acc += p
        u <= acc && return i
    end
    return length(log_probs)
end

"""
    rand([rng=default_rng()], hmm::HMM, T::Integer; seq::Bool=false) → Vector{UInt8} | (Vector{Int}, Vector{UInt8})

Simulate a trajectory of `T` timesteps from `hmm`.
"""
function Random.rand(rng::AbstractRNG, hmm::HMM, T::Integer; seq::Bool=false)
    T >= 0 || throw(ArgumentError("sequence length must be non-negative"))
    z = Vector{Int}(undef, T)
    obs = Vector{UInt8}(undef, T)
    T == 0 && return seq ? (z, obs) : obs

    z[1] = _hmm_sample_categorical(rng, hmm.initial)
    col1 = _hmm_sample_categorical(rng, view(hmm.emissions, z[1], :))
    obs[1] = hmm.alphabet[col1]

    @inbounds for t in 2:T
        z[t] = _hmm_sample_categorical(rng, view(hmm.transitions, z[t-1], :))
        col_t = _hmm_sample_categorical(rng, view(hmm.emissions, z[t], :))
        obs[t] = hmm.alphabet[col_t]
    end

    return seq ? (z, obs) : obs
end

Random.rand(hmm::HMM, T::Integer; kwargs...) = Random.rand(Random.default_rng(), hmm, T; kwargs...)

"""
    joint_loglikelihood(hmm, sequence, states) → Float64

Compute the joint log-likelihood log P(sequence, states | HMM).
"""
function joint_loglikelihood(hmm::HMM{T}, sequence::AbstractVector{UInt8}, states::AbstractVector{<:Integer}) where {T}
    L = length(sequence)
    length(states) == L || throw(ArgumentError("sequence and states must have matching length"))
    L == 0 && return zero(T)

    s1 = states[1]
    (1 <= s1 <= length(hmm.states)) || throw(BoundsError(hmm.states, s1))
    ll = hmm.initial[s1] + _get_emission_logprob(hmm, s1, sequence[1])

    @inbounds for t in 2:L
        prev_s = states[t-1]
        curr_s = states[t]
        (1 <= curr_s <= length(hmm.states)) || throw(BoundsError(hmm.states, curr_s))
        ll += hmm.transitions[prev_s, curr_s] + _get_emission_logprob(hmm, curr_s, sequence[t])
    end

    return ll
end

joint_loglikelihood(hmm::HMM, seq::BioSequence, states) = joint_loglikelihood(hmm, _hmm_bytes(seq), states)
joint_loglikelihood(hmm::HMM, seq::AbstractString, states) = joint_loglikelihood(hmm, _hmm_bytes(seq), states)

@inline function _get_emission_logprob(hmm::HMM{T}, state_idx::Int, byte::UInt8) where {T}
    col = hmm._emission_lookup[Int(byte)+1]
    col == 0 && return hmm._unknown_emission
    return hmm.emissions[state_idx, col]
end

@inline log_emission(hmm::HMM, state_idx::Int, byte::UInt8) = _get_emission_logprob(hmm, state_idx, byte)

@inline _hmm_bytes(sequence::AbstractVector{UInt8}) = sequence
@inline _hmm_bytes(sequence::BioSequence) = sequence.data
@inline _hmm_bytes(sequence::AbstractString) = collect(codeunits(String(sequence)))
@inline _hmm_bytes(sequence) = collect(codeunits(String(sequence)))

@inline function logaddexp(x::T, y::T) where {T<:AbstractFloat}
    (isnan(x) || isnan(y)) && throw(DomainError((x, y), "logaddexp received NaN"))
    x == -Inf && return y; y == -Inf && return x
    m = max(x, y); return m + log1p(exp(min(x,y) - m))
end

# ===========================================================================
# Viterbi Algorithm (Generic for AbstractHMM and specialized for discrete HMM)
# ===========================================================================

"""
    viterbi(hmm, sequence) → (path, log_prob)

Compute the most likely hidden state path using the Viterbi algorithm.
"""
function viterbi(hmm::AbstractHMM{T}, sequence::AbstractVector) where {T}
    L = length(sequence); N = n_states(hmm)
    _ctx = active_provenance_context()
    L == 0 && return _register_hmm_result!(_ctx, (Int[], T(-Inf)), "viterbi"; parameters=(seq_length=0, n_states=N))
    v = Matrix{T}(undef, N, L)
    tb = Matrix{Int}(undef, N, L)
    ini = initial_distribution(hmm)
    trn = transition_matrix(hmm)
    @inbounds for i in 1:N
        v[i,1] = ini[i] + log_emission(hmm, i, sequence[1])
        tb[i,1] = 0
    end
    @inbounds for t in 2:L
        obs_t = sequence[t]
        for j in 1:N
            best_p, best_i = T(-Inf), 0
            for i in 1:N
                p = v[i,t-1] + trn[i,j]
                p > best_p && (best_p = p; best_i = i)
            end
            v[j,t] = best_p + log_emission(hmm, j, obs_t)
            tb[j,t] = best_i
        end
    end
    best_p, best_last = T(-Inf), 0
    @inbounds for i in 1:N; v[i,L] > best_p && (best_p = v[i,L]; best_last = i); end
    best_last == 0 && throw(ArgumentError("observation sequence has zero probability under this HMM"))
    path = Vector{Int}(undef, L)
    curr = best_last
    @inbounds for t in L:-1:1; path[t] = curr; curr = tb[curr,t]; end
    result = (path, best_p)
    return _register_hmm_result!(_ctx, result, "viterbi"; parameters=(seq_length=L, n_states=N, log_prob=Float64(best_p)))
end

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
    best_last == 0 && throw(ArgumentError("observation sequence has zero probability under this HMM"))
    path = Vector{Int}(undef, L)
    curr = best_last
    @inbounds for t in L:-1:1; path[t] = curr; curr = tb[curr,t]; end
    result = (path, best_p)
    return _register_hmm_result!(_ctx, result, "viterbi"; parameters=(seq_length=L, n_states=N, log_prob=Float64(best_p)))
end

viterbi(hmm::AbstractHMM, seq::BioSequence; _ctx=nothing)     = viterbi(hmm, _hmm_bytes(seq))
viterbi(hmm::AbstractHMM, seq::AbstractString; _ctx=nothing) = viterbi(hmm, _hmm_bytes(seq))

# ===========================================================================
# Forward Algorithm (Generic AbstractHMM & Specialized HMM)
# ===========================================================================

"""
    forward(hmm, sequence) → log_probability
"""
function forward(hmm::AbstractHMM{T}, sequence::AbstractVector) where {T}
    L = length(sequence); N = n_states(hmm)
    _ctx = active_provenance_context()
    L == 0 && return _register_hmm_result!(_ctx, T(-Inf), "forward"; parameters=(seq_length=0, n_states=N))
    alpha = Matrix{T}(undef, N, 2); cc, pc = 1, 2
    ini = initial_distribution(hmm)
    trn = transition_matrix(hmm)
    @inbounds for i in 1:N
        alpha[i,cc] = ini[i] + log_emission(hmm, i, sequence[1])
    end
    @inbounds for t in 2:L
        cc, pc = pc, cc
        obs_t = sequence[t]
        for j in 1:N
            s = T(-Inf)
            for i in 1:N; s = logaddexp(s, alpha[i,pc] + trn[i,j]); end
            alpha[j,cc] = s + log_emission(hmm, j, obs_t)
        end
    end
    total = T(-Inf)
    @inbounds for i in 1:N; total = logaddexp(total, alpha[i,cc]); end
    return _register_hmm_result!(_ctx, total, "forward"; parameters=(seq_length=L, n_states=N, log_likelihood=Float64(total)))
end

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

forward(hmm::AbstractHMM, seq::BioSequence; _ctx=nothing)    = forward(hmm, _hmm_bytes(seq))
forward(hmm::AbstractHMM, seq::AbstractString; _ctx=nothing) = forward(hmm, _hmm_bytes(seq))

# ===========================================================================
# Backward Algorithm
# ===========================================================================

"""
    backward(hmm, sequence) → N × L log-probability matrix
"""
function backward(hmm::AbstractHMM{T}, sequence::AbstractVector) where {T}
    L = length(sequence); N = n_states(hmm)
    beta = Matrix{T}(undef, N, L)
    L == 0 && return beta
    trn = transition_matrix(hmm)
    @inbounds for i in 1:N; beta[i,L] = T(0); end
    @inbounds for t in L-1:-1:1
        obs_next = sequence[t+1]
        for i in 1:N
            s = T(-Inf)
            for j in 1:N
                s = logaddexp(s, trn[i,j] + log_emission(hmm, j, obs_next) + beta[j,t+1])
            end
            beta[i,t] = s
        end
    end
    return beta
end

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

backward(hmm::AbstractHMM, seq::BioSequence)    = backward(hmm, _hmm_bytes(seq))
backward(hmm::AbstractHMM, seq::AbstractString) = backward(hmm, _hmm_bytes(seq))

# ===========================================================================
# Full forward table (N × L) — needed for Baum-Welch & Posterior Decoding
# ===========================================================================

function _forward_table(hmm::AbstractHMM{T}, sequence::AbstractVector) where {T}
    L = length(sequence); N = n_states(hmm)
    alpha = fill(T(-Inf), N, L)
    ini = initial_distribution(hmm)
    trn = transition_matrix(hmm)
    @inbounds for i in 1:N
        alpha[i,1] = ini[i] + log_emission(hmm, i, sequence[1])
    end
    @inbounds for t in 2:L, j in 1:N
        s = T(-Inf)
        for i in 1:N; s = logaddexp(s, alpha[i,t-1] + trn[i,j]); end
        alpha[j,t] = s + log_emission(hmm, j, sequence[t])
    end
    return alpha
end

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
# HMM Log-Likelihood
# ===========================================================================

"""
    hmm_log_likelihood(hmm, sequences) → Float64

Compute the total log-likelihood of a set of observed sequences given `hmm`.
"""
function hmm_log_likelihood(hmm::AbstractHMM, sequences)
    total = 0.0
    for seq in sequences
        total += forward(hmm, seq)
    end
    return total
end

# ===========================================================================
# Posterior State Probabilities (Smoothed Decoding)
# ===========================================================================

"""
    posterior_state_probabilities(hmm, sequence) → N × L Matrix{Float64}

Compute γ_t(i) = P(state=i at t | sequence, hmm) using forward-backward.
Equivalent to the "E-step" responsibilities without accumulating.
"""
function posterior_state_probabilities(hmm::AbstractHMM{T}, sequence::AbstractVector) where {T}
    L = length(sequence); N = n_states(hmm)
    L == 0 && return Matrix{Float64}(undef, N, 0)
    alpha = _forward_table(hmm, sequence)
    beta  = backward(hmm, sequence)
    log_px = forward(hmm, sequence)
    isfinite(log_px) || throw(ArgumentError("observation sequence has zero probability under this HMM"))
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

posterior_state_probabilities(hmm::AbstractHMM, seq::AbstractString; _ctx=nothing) =
    posterior_state_probabilities(hmm, _hmm_bytes(seq))

"""Return smoothed state posteriors and the sequence log likelihood."""
function forward_backward(hmm::AbstractHMM, sequence)
    return posterior_state_probabilities(hmm, sequence), forward(hmm, sequence)
end

posteriors(hmm::AbstractHMM, sequence) = first(forward_backward(hmm, sequence))
loglikelihood(hmm::AbstractHMM, sequence) = forward(hmm, sequence)

# ===========================================================================
# Posterior Decode (MAP assignment per position)
# ===========================================================================

"""
    posterior_decode(hmm, sequence) → (assignments::Vector{Int}, gamma::Matrix{Float64})

Assign each position to the highest-posterior state. Softer than Viterbi;
does not guarantee globally-consistent paths. Analogous to `hmmlearn.decode`.
"""
function posterior_decode(hmm::AbstractHMM, sequence)
    gamma = posterior_state_probabilities(hmm, sequence)
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
            isfinite(log_px) || throw(ArgumentError("observation sequence has zero probability under this HMM"))
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
        # Preserve distributions for states/rows with zero expected count;
        # manufacturing eps-sized probabilities makes EM silently invalid.
        _hmm_update_log_vector!(hmm.initial, ini_acc)
        _hmm_update_log_rows!(hmm.transitions, trans_acc)
        _hmm_update_log_rows!(hmm.emissions, emis_acc)
    end
    return log_likelihoods
end

"""
    baum_welch_train(n_states, alphabet, sequences; kwargs...) → (hmm, log_likelihoods)

Randomly initialise and Baum-Welch train a new HMM from scratch.
"""
function baum_welch_train(n_states::Int, alphabet::Vector{UInt8}, sequences;
                          max_iter::Int=100, tol::Real=1e-4, seed::Int=1, verbose::Bool=false)
    rng   = MersenneTwister(seed)
    N, V  = n_states, length(alphabet)
    ini_values = rand(rng, Float64, N)
    ini   = ini_values ./ sum(ini_values)
    trans = Matrix{Float64}(undef, N, N)
    for i in 1:N
        values = rand(rng, Float64, N)
        trans[i, :] = values ./ sum(values)
    end
    emis  = Matrix{Float64}(undef, N, V)
    for i in 1:N
        values = rand(rng, Float64, V)
        emis[i, :] = values ./ sum(values)
    end
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
    for _ in 1:max_iter
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

        _hmm_update_log_vector!(hmm.initial, ini_acc)
        _hmm_update_log_rows!(hmm.transitions, trans_acc)
        _hmm_update_log_rows!(hmm.emissions, emis_acc)

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
    bytes = _hmm_bytes(sequence)

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
                hmm.states[prev_state], log(max(sum(probs[seg_start:t-1]) / (t - seg_start), eps()))))
            seg_start = t; prev_state = assignments[t]
        end
    end
    push!(segments, HMMSegment(seg_start, length(assignments), prev_state,
        hmm.states[prev_state], log(max(sum(probs[seg_start:end]) / (length(probs) - seg_start + 1), eps()))))

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
    score = best - phmm.null_log_odds

    # Prefix scores are useful for locating the strongest local hit while
    # retaining the scalar return value used by the ordinary API.
    if return_per_position
        per_position = [max(M[i+1,L+1], I[i+1,L+1], D[i+1,L+1]) - phmm.null_log_odds
                        for i in 1:n]
        return (score=score, per_position=per_position)
    end
    return score
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

# ===========================================================================
# Gaussian HMM (Continuous Observation Signals: CNV, Intensity, Mass Spec)
# ===========================================================================

"""
    GaussianHMM{T<:AbstractFloat}

A Hidden Markov Model with 1D Gaussian emissions per state for continuous signals.
"""
struct GaussianHMM{T<:AbstractFloat} <: AbstractHMM{T}
    states::Vector{String}
    initial::Vector{T}        # Log-space
    transitions::Matrix{T}    # Log-space N x N
    means::Vector{T}          # State means μ_k
    vars::Vector{T}           # State variances σ_k^2 (>0)

    function GaussianHMM{T}(states, initial, transitions, means, vars) where {T}
        N = length(states)
        length(initial) == N || throw(ArgumentError("initial must match number of states"))
        size(transitions) == (N, N) || throw(ArgumentError("transitions matrix must be N x N"))
        length(means) == N || throw(ArgumentError("means vector must match number of states"))
        length(vars) == N || throw(ArgumentError("vars vector must match number of states"))
        any(v -> v <= 0 || !isfinite(v), vars) && throw(ArgumentError("variances must be strictly positive and finite"))
        new{T}(String.(states), T.(initial), T.(transitions), T.(means), T.(vars))
    end
end

@inline log_emission(ghmm::GaussianHMM, state_idx::Int, x::Real) = _gaussian_logpdf(ghmm.means[state_idx], ghmm.vars[state_idx], x)

function GaussianHMM(states::AbstractVector{<:AbstractString}, initial::AbstractVector,
                     transitions::AbstractMatrix, means::AbstractVector, vars::AbstractVector;
                     log_space::Bool=false)
    T = Float64
    ini = log_space ? _hmm_normalize_log_vector(initial, "initial") : _hmm_log_probabilities(initial, "initial")
    trn = log_space ? _hmm_log_matrix_from_logs(transitions, "transition") : _hmm_log_matrix(transitions, "transition")
    return GaussianHMM{T}(collect(states), ini, trn, collect(means), collect(vars))
end

@inline function _gaussian_logpdf(mean::T, var::T, x::Real) where {T}
    dx = Float64(x) - mean
    return -0.5 * (log(2π * var) + (dx * dx) / var)
end

function forward(ghmm::GaussianHMM{T}, sequence::AbstractVector{<:Real}) where {T}
    L = length(sequence); N = length(ghmm.states)
    L == 0 && return T(-Inf)
    alpha = Matrix{T}(undef, N, 2); cc, pc = 1, 2
    @inbounds for i in 1:N
        alpha[i,cc] = ghmm.initial[i] + _gaussian_logpdf(ghmm.means[i], ghmm.vars[i], sequence[1])
    end
    @inbounds for t in 2:L
        cc, pc = pc, cc
        x_t = sequence[t]
        for j in 1:N
            s = T(-Inf)
            for i in 1:N; s = logaddexp(s, alpha[i,pc] + ghmm.transitions[i,j]); end
            alpha[j,cc] = s + _gaussian_logpdf(ghmm.means[j], ghmm.vars[j], x_t)
        end
    end
    total = T(-Inf)
    @inbounds for i in 1:N; total = logaddexp(total, alpha[i,cc]); end
    return total
end

function backward(ghmm::GaussianHMM{T}, sequence::AbstractVector{<:Real}) where {T}
    L = length(sequence); N = length(ghmm.states)
    beta = Matrix{T}(undef, N, L)
    L == 0 && return beta
    @inbounds for i in 1:N; beta[i,L] = T(0); end
    @inbounds for t in L-1:-1:1
        x_next = sequence[t+1]
        for i in 1:N
            s = T(-Inf)
            for j in 1:N
                s = logaddexp(s, ghmm.transitions[i,j] + _gaussian_logpdf(ghmm.means[j], ghmm.vars[j], x_next) + beta[j,t+1])
            end
            beta[i,t] = s
        end
    end
    return beta
end

function viterbi(ghmm::GaussianHMM{T}, sequence::AbstractVector{<:Real}) where {T}
    L = length(sequence); N = length(ghmm.states)
    L == 0 && return (Int[], T(-Inf))
    v = Matrix{T}(undef, N, L)
    tb = Matrix{Int}(undef, N, L)
    @inbounds for i in 1:N
        v[i,1] = ghmm.initial[i] + _gaussian_logpdf(ghmm.means[i], ghmm.vars[i], sequence[1])
        tb[i,1] = 0
    end
    @inbounds for t in 2:L
        x_t = sequence[t]
        for j in 1:N
            best_p, best_i = T(-Inf), 0
            for i in 1:N
                p = v[i,t-1] + ghmm.transitions[i,j]
                p > best_p && (best_p = p; best_i = i)
            end
            v[j,t] = best_p + _gaussian_logpdf(ghmm.means[j], ghmm.vars[j], x_t)
            tb[j,t] = best_i
        end
    end
    best_p, best_last = T(-Inf), 0
    @inbounds for i in 1:N; v[i,L] > best_p && (best_p = v[i,L]; best_last = i); end
    best_last == 0 && throw(ArgumentError("observation sequence has zero probability under Gaussian HMM"))
    path = Vector{Int}(undef, L)
    curr = best_last
    @inbounds for t in L:-1:1; path[t] = curr; curr = tb[curr,t]; end
    return path, best_p
end

function posterior_state_probabilities(ghmm::GaussianHMM{T}, sequence::AbstractVector{<:Real}) where {T}
    L = length(sequence); N = length(ghmm.states)
    L == 0 && return Matrix{Float64}(undef, N, 0)
    alpha = fill(T(-Inf), N, L)
    @inbounds for i in 1:N
        alpha[i,1] = ghmm.initial[i] + _gaussian_logpdf(ghmm.means[i], ghmm.vars[i], sequence[1])
    end
    @inbounds for t in 2:L, j in 1:N
        s = T(-Inf)
        for i in 1:N; s = logaddexp(s, alpha[i,t-1] + ghmm.transitions[i,j]); end
        alpha[j,t] = s + _gaussian_logpdf(ghmm.means[j], ghmm.vars[j], sequence[t])
    end
    beta = backward(ghmm, sequence)
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

function baum_welch!(ghmm::GaussianHMM{T}, sequences; max_iter::Int=100, tol::Real=1e-4) where {T}
    N = length(ghmm.states)
    seqs = [collect(Float64.(s)) for s in sequences]
    log_likelihoods = Float64[]
    ini_acc  = zeros(Float64, N)
    trans_acc = zeros(Float64, N, N)

    for iter in 1:max_iter
        fill!(ini_acc, 0.0)
        fill!(trans_acc, 0.0)
        mean_num = zeros(Float64, N)
        mean_den = zeros(Float64, N)
        total_ll = 0.0

        for s in seqs
            L = length(s)
            L < 1 && continue
            gamma = posterior_state_probabilities(ghmm, s)
            log_px = forward(ghmm, s)
            total_ll += log_px

            for t in 1:L
                for i in 1:N
                    γ = gamma[i,t]
                    t == 1 && (ini_acc[i] += γ)
                    mean_num[i] += γ * s[t]
                    mean_den[i] += γ
                end
            end

            alpha = fill(T(-Inf), N, L)
            @inbounds for i in 1:N; alpha[i,1] = ghmm.initial[i] + _gaussian_logpdf(ghmm.means[i], ghmm.vars[i], s[1]); end
            @inbounds for t in 2:L, j in 1:N
                st = T(-Inf)
                for i in 1:N; st = logaddexp(st, alpha[i,t-1] + ghmm.transitions[i,j]); end
                alpha[j,t] = st + _gaussian_logpdf(ghmm.means[j], ghmm.vars[j], s[t])
            end
            beta = backward(ghmm, s)

            for t in 1:(L-1)
                denom = -Inf
                for i in 1:N, j in 1:N
                    denom = logaddexp(denom, alpha[i,t] + ghmm.transitions[i,j] + _gaussian_logpdf(ghmm.means[j], ghmm.vars[j], s[t+1]) + beta[j,t+1])
                end
                for i in 1:N, j in 1:N
                    ξ = exp(alpha[i,t] + ghmm.transitions[i,j] + _gaussian_logpdf(ghmm.means[j], ghmm.vars[j], s[t+1]) + beta[j,t+1] - denom)
                    trans_acc[i,j] += ξ
                end
            end
        end

        push!(log_likelihoods, total_ll)

        _hmm_update_log_vector!(ghmm.initial, ini_acc)
        _hmm_update_log_rows!(ghmm.transitions, trans_acc)
        for i in 1:N
            if mean_den[i] > 0
                new_m = mean_num[i] / mean_den[i]
                ghmm.means[i] = new_m
                v_num = 0.0
                for s in seqs
                    gamma_s = posterior_state_probabilities(ghmm, s)
                    for t in 1:length(s)
                        dx = s[t] - new_m
                        v_num += gamma_s[i,t] * dx * dx
                    end
                end
                ghmm.vars[i] = max(v_num / mean_den[i], 1e-4)
            end
        end

        if length(log_likelihoods) >= 2 && abs(log_likelihoods[end] - log_likelihoods[end-1]) < Float64(tol)
            break
        end
    end
    return log_likelihoods
end

# ===========================================================================
# ChromHMM (Multivariate Epigenomic Binary Signal HMM)
# ===========================================================================

"""
    ChromHMM{T<:AbstractFloat}

Multivariate Bernoulli HMM for chromatin state learning across multi-track histone mark signals (ChromHMM architecture).
"""
struct ChromHMM{T<:AbstractFloat} <: AbstractHMM{T}
    states::Vector{String}
    mark_names::Vector{String}
    initial::Vector{T}              # Log-space
    transitions::Matrix{T}          # Log-space N x N
    emission_probs::Matrix{Float64} # N states x M marks (in [0,1])

    function ChromHMM{T}(states, mark_names, initial, transitions, emission_probs) where {T}
        N = length(states); M = length(mark_names)
        length(initial) == N || throw(ArgumentError("initial must match number of states"))
        size(transitions) == (N, N) || throw(ArgumentError("transitions matrix must be N x N"))
        size(emission_probs) == (N, M) || throw(ArgumentError("emission_probs matrix must be N x M"))
        any(p -> p < 0 || p > 1, emission_probs) && throw(ArgumentError("emission probabilities must be in [0,1]"))
        new{T}(String.(states), String.(mark_names), T.(initial), T.(transitions), Float64.(emission_probs))
    end
end

function ChromHMM(states::AbstractVector{<:AbstractString}, mark_names::AbstractVector{<:AbstractString},
                  initial::AbstractVector, transitions::AbstractMatrix, emission_probs::AbstractMatrix;
                  log_space::Bool=false)
    T = Float64
    ini = log_space ? _hmm_normalize_log_vector(initial, "initial") : _hmm_log_probabilities(initial, "initial")
    trn = log_space ? _hmm_log_matrix_from_logs(transitions, "transition") : _hmm_log_matrix(transitions, "transition")
    return ChromHMM{T}(collect(states), collect(mark_names), ini, trn, collect(emission_probs))
end

@inline function _chrom_emission_logprob(chmm::ChromHMM, state_idx::Int, mark_vector::AbstractVector{<:Real})
    logp = 0.0
    @inbounds for m in 1:length(mark_vector)
        p = clamp(chmm.emission_probs[state_idx, m], 1e-6, 1.0 - 1e-6)
        logp += mark_vector[m] > 0 ? log(p) : log(1.0 - p)
    end
    return logp
end

@inline log_emission(chmm::ChromHMM, state_idx::Int, mark_vector::AbstractVector{<:Real}) = _chrom_emission_logprob(chmm, state_idx, mark_vector)

function viterbi(chmm::ChromHMM{T}, mark_matrix::AbstractMatrix{<:Real}) where {T}
    M, L = size(mark_matrix); N = length(chmm.states)
    L == 0 && return (Int[], T(-Inf))
    v = Matrix{T}(undef, N, L)
    tb = Matrix{Int}(undef, N, L)
    @inbounds for i in 1:N
        v[i,1] = chmm.initial[i] + _chrom_emission_logprob(chmm, i, view(mark_matrix, :, 1))
        tb[i,1] = 0
    end
    @inbounds for t in 2:L
        m_t = view(mark_matrix, :, t)
        for j in 1:N
            best_p, best_i = T(-Inf), 0
            for i in 1:N
                p = v[i,t-1] + chmm.transitions[i,j]
                p > best_p && (best_p = p; best_i = i)
            end
            v[j,t] = best_p + _chrom_emission_logprob(chmm, j, m_t)
            tb[j,t] = best_i
        end
    end
    best_p, best_last = T(-Inf), 0
    @inbounds for i in 1:N; v[i,L] > best_p && (best_p = v[i,L]; best_last = i); end
    path = Vector{Int}(undef, L)
    curr = best_last
    @inbounds for t in L:-1:1; path[t] = curr; curr = tb[curr,t]; end
    return path, best_p
end

"""
    segment_chromatin(chmm, mark_matrix) → Vector{HMMSegment}

Segment genome into chromatin states based on multi-track histone mark binary signals.
"""
function segment_chromatin(chmm::ChromHMM, mark_matrix::AbstractMatrix{<:Real})
    path, log_p = viterbi(chmm, mark_matrix)
    L = length(path)
    segments = HMMSegment[]
    L == 0 && return segments
    seg_start = 1; prev_state = path[1]
    for t in 2:L
        if path[t] != prev_state
            push!(segments, HMMSegment(seg_start, t-1, prev_state, chmm.states[prev_state], log_p / L))
            seg_start = t; prev_state = path[t]
        end
    end
    push!(segments, HMMSegment(seg_start, L, prev_state, chmm.states[prev_state], log_p / L))
    return segments
end

# ===========================================================================
# Model Selection, Entropy & StatsAPI Wrappers
# ===========================================================================

"""
    aic(hmm, sequences) → Float64

Akaike Information Criterion: AIC = 2*k - 2*logL.
"""
function aic(hmm::HMM, sequences)
    k = nparams(hmm)
    ll = hmm_log_likelihood(hmm, sequences)
    return 2 * k - 2 * ll
end

"""
    bic(hmm, sequences) → Float64

Bayesian Information Criterion: BIC = k*ln(N_obs) - 2*logL.
"""
function bic(hmm::HMM, sequences)
    k = nparams(hmm)
    total_obs = sum(length(s) for s in sequences)
    ll = hmm_log_likelihood(hmm, sequences)
    return k * log(max(total_obs, 1)) - 2 * ll
end

"""
    state_entropy(gamma::AbstractMatrix{Float64}) → Vector{Float64}

Compute Shannon state assignment entropy per position from posterior probabilities γ (N x L).
"""
function state_entropy(gamma::AbstractMatrix{Float64})
    L = size(gamma, 2)
    entropies = Vector{Float64}(undef, L)
    for t in 1:L
        h = 0.0
        for i in 1:size(gamma, 1)
            p = gamma[i,t]
            p > 0 && (h -= p * log2(p))
        end
        entropies[t] = h
    end
    return entropies
end

"""
    cross_validate_hmm(sequences, n_states_range; k_folds=5, seed=1) → Dict{Int,Float64}

Perform K-fold cross-validation over `n_states_range` to select the optimal number of HMM states.
"""
function cross_validate_hmm(sequences, n_states_range; k_folds::Int=5, seed::Int=1)
    rng = MersenneTwister(seed)
    n_seqs = length(sequences)
    n_seqs >= k_folds || throw(ArgumentError("Number of sequences ($n_seqs) must be >= k_folds ($k_folds)"))
    shuffled_idx = randperm(rng, n_seqs)
    fold_size = ceil(Int, n_seqs / k_folds)

    results = Dict{Int, Float64}()
    for k in n_states_range
        val_lls = Float64[]
        for f in 1:k_folds
            val_indices = shuffled_idx[((f-1)*fold_size + 1):min(f*fold_size, n_seqs)]
            train_indices = setdiff(shuffled_idx, val_indices)
            train_seqs = sequences[train_indices]
            val_seqs = sequences[val_indices]

            alphabet = UInt8.(collect(b"ACGT"))
            if !isempty(sequences) && sequences[1] isa AbstractVector{UInt8}
                alphabet = unique(vcat(sequences...))
            end

            trained_hmm, _ = baum_welch_train(k, alphabet, train_seqs; seed=seed+f, max_iter=30)
            push!(val_lls, hmm_log_likelihood(trained_hmm, val_seqs))
        end
        results[k] = mean(val_lls)
    end
    return results
end

"""
    fit!(hmm::HMM, sequences; kwargs...)

In-place Baum-Welch training wrapper for standard StatsAPI interface.
"""
fit!(hmm::HMM, sequences; kwargs...) = baum_welch!(hmm, sequences; kwargs...)

"""
    logdensityof(hmm::HMM, sequence) → Float64

Alias for forward sequence log-likelihood.
"""
logdensityof(hmm::HMM, sequence) = forward(hmm, sequence)

# ===========================================================================
# Blender 3D Integrator Payload Overload
# ===========================================================================

import .BlenderIntegrator: BlenderMaterial, BlenderHMMPayload, to_blender_payload

"""
    to_blender_payload(hmm::AbstractHMM, sequence; name="HMM_Landscape") → BlenderHMMPayload

Convert any HMM model and sequence into a 3D BlenderHMMPayload visualization node.
"""
function to_blender_payload(hmm::AbstractHMM, sequence; name::String="HMM_Landscape")
    gamma = posterior_state_probabilities(hmm, sequence)
    st = states(hmm)
    trans = exp.(transition_matrix(hmm))
    mat = BlenderMaterial(name="HMMMaterial", color=(0.2, 0.6, 0.9, 1.0))
    return BlenderHMMPayload(name, st, trans, gamma, mat)
end


