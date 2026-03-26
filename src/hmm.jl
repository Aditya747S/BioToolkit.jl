export HMM, viterbi, forward, backward

"""
    HMM{T<:AbstractFloat}

A Hidden Markov Model operating entirely in log-space for numerical stability.
"""
struct HMM{T<:AbstractFloat}
    states::Vector{String}
    alphabet::Vector{UInt8}
    
    # Probabilities strictly in log-space to prevent underflow
    initial::Vector{T}
    transitions::Matrix{T}  # [from_state, to_state]
    emissions::Matrix{T}    # [state, alphabet_index]
    
    # Fast O(1) lookup table mapping an arbitrary UInt8 byte to its alphabet column index
    _emission_lookup::Vector{Int}
    _unknown_emission::T
    
    function HMM{T}(states::Vector{String}, alphabet::Vector{UInt8}, initial::Vector{T}, transitions::Matrix{T}, emissions::Matrix{T}, unknown_prob::T=T(-Inf)) where {T}
        length(states) == length(initial) || throw(ArgumentError("initial must match states"))
        size(transitions) == (length(states), length(states)) || throw(ArgumentError("transitions must be N x N"))
        size(emissions) == (length(states), length(alphabet)) || throw(ArgumentError("emissions must be N x V"))
        
        lookup = fill(0, 256)
        for (i, byte) in enumerate(alphabet)
            lookup[Int(byte) + 1] = i
        end
        
        new{T}(states, alphabet, initial, transitions, emissions, lookup, unknown_prob)
    end
end

function HMM(states::Vector{String}, alphabet::Vector{UInt8}, initial::AbstractVector, transitions::AbstractMatrix, emissions::AbstractMatrix; log_space::Bool=false)
    T = Float64
    ini = T.(log_space ? initial : log.(initial))
    trn = Matrix{T}(log_space ? transitions : log.(transitions))
    ems = Matrix{T}(log_space ? emissions : log.(emissions))
    return HMM{T}(states, alphabet, ini, trn, ems)
end

@inline function _get_emission_logprob(hmm::HMM{T}, state_idx::Int, byte::UInt8) where {T}
    col = hmm._emission_lookup[Int(byte) + 1]
    col == 0 && return hmm._unknown_emission
    return hmm.emissions[state_idx, col]
end

# Log-Sum-Exp trick for numerical stability
@inline function logaddexp(x::T, y::T) where {T<:AbstractFloat}
    x == -Inf && return y
    y == -Inf && return x
    max_val = max(x, y)
    return max_val + log1p(exp(min(x, y) - max_val))
end

"""
    viterbi(hmm::HMM, sequence::AbstractVector{UInt8})

Computes the most likely sequence of hidden states (Viterbi path) and its log-probability.
Returns `(path::Vector{Int}, log_prob::Float64)`.
"""
function viterbi(hmm::HMM{T}, sequence::AbstractVector{UInt8}) where {T}
    L = length(sequence)
    N = length(hmm.states)
    
    if L == 0
        return Int[], T(-Inf)
    end
    
    v_scores = Matrix{T}(undef, N, L)
    v_trace = Matrix{Int}(undef, N, L)
    
    # Initialization
    first_byte = sequence[1]
    @inbounds for i in 1:N
        v_scores[i, 1] = hmm.initial[i] + _get_emission_logprob(hmm, i, first_byte)
        v_trace[i, 1] = 0
    end
    
    # Iteration
    @inbounds for t in 2:L
        byte = sequence[t]
        for j in 1:N # Current state
            max_p = T(-Inf)
            max_i = 0
            
            for i in 1:N # Previous state
                p = v_scores[i, t - 1] + hmm.transitions[i, j]
                if p > max_p
                    max_p = p
                    max_i = i
                end
            end
            
            v_scores[j, t] = max_p + _get_emission_logprob(hmm, j, byte)
            v_trace[j, t] = max_i
        end
    end
    
    # Termination & Traceback
    best_p = T(-Inf)
    best_last_state = 0
    @inbounds for i in 1:N
        if v_scores[i, L] > best_p
            best_p = v_scores[i, L]
            best_last_state = i
        end
    end
    
    path = Vector{Int}(undef, L)
    if best_last_state != 0
        curr_state = best_last_state
        @inbounds for t in L:-1:1
            path[t] = curr_state
            curr_state = v_trace[curr_state, t]
        end
    end
    
    return path, best_p
end

viterbi(hmm::HMM, sequence::AbstractString) = viterbi(hmm, codeunits(sequence))

"""
    forward(hmm::HMM, sequence::AbstractVector{UInt8})

Computes the total log-probability of the sequence given the model over all valid paths.
"""
function forward(hmm::HMM{T}, sequence::AbstractVector{UInt8}) where {T}
    L = length(sequence)
    N = length(hmm.states)
    
    if L == 0
        return T(-Inf)
    end
    
    alpha = Matrix{T}(undef, N, 2)
    curr_col = 1
    prev_col = 2
    
    # Initialization
    first_byte = sequence[1]
    @inbounds for i in 1:N
        alpha[i, curr_col] = hmm.initial[i] + _get_emission_logprob(hmm, i, first_byte)
    end
    
    # Iteration
    @inbounds for t in 2:L
        curr_col, prev_col = prev_col, curr_col
        byte = sequence[t]
        for j in 1:N
            sum_p = T(-Inf)
            for i in 1:N
                p = alpha[i, prev_col] + hmm.transitions[i, j]
                sum_p = logaddexp(sum_p, p)
            end
            alpha[j, curr_col] = sum_p + _get_emission_logprob(hmm, j, byte)
        end
    end
    
    # Termination
    total_p = T(-Inf)
    @inbounds for i in 1:N
        total_p = logaddexp(total_p, alpha[i, curr_col])
    end
    
    return total_p
end

forward(hmm::HMM, sequence::AbstractString) = forward(hmm, codeunits(sequence))

"""
    backward(hmm::HMM, sequence::AbstractVector{UInt8})

Computes the backward probability (log p(X_{t+1:T} | z_t = i)) array.
Returns an NxL matrix.
"""
function backward(hmm::HMM{T}, sequence::AbstractVector{UInt8}) where {T}
    L = length(sequence)
    N = length(hmm.states)
    
    beta = Matrix{T}(undef, N, L)
    if L == 0
        return beta
    end
    
    # Init
    @inbounds for i in 1:N
        beta[i, L] = T(0.0) # log(1)
    end
    
    # Iteration
    @inbounds for t in L-1:-1:1
        next_byte = sequence[t+1]
        for i in 1:N
            sum_p = T(-Inf)
            for j in 1:N
                p = hmm.transitions[i, j] + _get_emission_logprob(hmm, j, next_byte) + beta[j, t+1]
                sum_p = logaddexp(sum_p, p)
            end
            beta[i, t] = sum_p
        end
    end
    
    return beta
end

backward(hmm::HMM, sequence::AbstractString) = backward(hmm, codeunits(sequence))

