# ==============================================================================
# align.jl — Pairwise and codon-aware sequence alignment
#
# Implements Needleman-Wunsch (global) and Smith-Waterman (local) alignment
# with both linear and affine gap penalty models.  The affine gap model uses
# the Gotoh (1982) three-matrix formulation with separate gap-open and
# gap-extend scores.
# Functions accept Vector{UInt8} and BioSequence typed inputs via multiple
# dispatch.
#
# References:
#   - Needleman & Wunsch (1970) JMB 48(3):443-453
#   - Smith & Waterman (1981) JMB 147(1):195-197
#   - Gotoh (1982) JMB 162(3):705-708 (affine gap extension)
#   - Henikoff & Henikoff (1992) PNAS 89(22):10915-10919 (BLOSUM matrices)
# ==============================================================================

"""
    PairwiseAlignmentResult{A <: BioAlphabet}

Result object returned by the pairwise alignment routines. It stores the aligned
left and right sequences as `BioSequence{A}`, the final score, the number of
exact matches, and the identity fraction.
"""
struct PairwiseAlignmentResult{A<:BioAlphabet} <: AbstractAnalysisResult
  left::BioSequence{A}
  right::BioSequence{A}
  score::Int
  matches::Int
  identity::Float64
  metadata::Dict{Symbol,Any}
end

function PairwiseAlignmentResult(left::BioSequence{A}, right::BioSequence{A}, score::Integer, matches::Integer, identity::Real; metadata::AbstractDict=Dict{Symbol,Any}()) where {A<:BioAlphabet}
  metadata_copy = Dict{Symbol,Any}(metadata)
  ensure_provenance_id!(metadata_copy)
  return PairwiseAlignmentResult{A}(left, right, Int(score), Int(matches), Float64(identity), metadata_copy)
end

analysis_result_fields(::Type{PairwiseAlignmentResult{A}}) where {A} = (:left, :right, :score, :matches, :identity)

"""
    AbstractPairwiseScoring

Abstract parent type for all pairwise alignment scoring models.
"""
abstract type AbstractPairwiseScoring end


"""
    PairHMMScoring{S}

Probability model for global pair-HMM alignment. `match_prob` and
`mismatch_prob` are emission probabilities for the match state. `gap_open_prob`
controls transitions from match to insertion/deletion states and
`gap_extend_prob` controls self-transitions in gap states.
"""
struct PairHMMScoring{S<:AbstractFloat} <: AbstractPairwiseScoring
  match_prob::S
  mismatch_prob::S
  gap_open_prob::S
  gap_extend_prob::S
end

function PairHMMScoring(; match_prob::Real=0.97, mismatch_prob::Real=0.03, gap_open_prob::Real=0.02, gap_extend_prob::Real=0.70)
  return PairHMMScoring(Float64(match_prob), Float64(mismatch_prob), Float64(gap_open_prob), Float64(gap_extend_prob))
end

"""
    PosteriorAlignmentResult{A}

Posterior-decoded pair-HMM alignment summary. `posterior_matrix[i,j,k]` stores
posterior probability for state `k` at DP cell `(i,j)`, where `k=1` is match,
`k=2` is a gap in the right sequence, and `k=3` is a gap in the left sequence.
"""
struct PosteriorAlignmentResult{A<:BioAlphabet} <: AbstractAnalysisResult
  posterior_matrix::Array{Float32,3}
  expected_accuracy::Float64
  consensus_alignment::PairwiseAlignmentResult{A}
  log_likelihood::Float64
  state_order::Vector{Symbol}
  metadata::Dict{Symbol,Any}
end

function PosteriorAlignmentResult(posterior_matrix::Array{Float32,3}, expected_accuracy::Real, consensus_alignment::PairwiseAlignmentResult{A}, log_likelihood::Real; metadata::AbstractDict=Dict{Symbol,Any}()) where {A<:BioAlphabet}
  metadata_copy = Dict{Symbol,Any}(metadata)
  ensure_provenance_id!(metadata_copy)
  return PosteriorAlignmentResult{A}(posterior_matrix, Float64(expected_accuracy), consensus_alignment, Float64(log_likelihood), [:match, :gap_right, :gap_left], metadata_copy)
end

analysis_result_fields(::Type{PosteriorAlignmentResult{A}}) where {A} = (:expected_accuracy, :log_likelihood, :state_order)

"""
    DifferentiableScoring{F}

Smooth alignment scoring parameters for differentiable dynamic programming.
`weights` are `[match, mismatch, gap]` by default.
"""
struct DifferentiableScoring{F<:AbstractFloat} <: AbstractPairwiseScoring
  weights::Vector{F}
  temperature::F
end

function DifferentiableScoring(weights::AbstractVector{<:Real}; temperature::Real=1.0)
  length(weights) >= 3 || throw(ArgumentError("DifferentiableScoring requires at least [match, mismatch, gap] weights"))
  temperature > 0 || throw(ArgumentError("temperature must be positive"))
  return DifferentiableScoring(Float64.(weights), Float64(temperature))
end

@inline function _logsumexp3(a::Float64, b::Float64, c::Float64)
  m = max(a, b, c)
  isfinite(m) || return m
  return m + log(exp(a - m) + exp(b - m) + exp(c - m))
end

@inline function _logaddexp(a::Float64, b::Float64)
  m = max(a, b)
  isfinite(m) || return m
  return m + log(exp(a - m) + exp(b - m))
end

function _validate_pairhmm_scoring(scoring::PairHMMScoring)
  0 < scoring.match_prob <= 1 || throw(ArgumentError("match_prob must be in (0, 1]"))
  0 < scoring.mismatch_prob <= 1 || throw(ArgumentError("mismatch_prob must be in (0, 1]"))
  0 < scoring.gap_open_prob < 0.5 || throw(ArgumentError("gap_open_prob must be in (0, 0.5)"))
  0 < scoring.gap_extend_prob < 1 || throw(ArgumentError("gap_extend_prob must be in (0, 1)"))
  return nothing
end

@inline function _pairhmm_emit_log(left_byte::UInt8, right_byte::UInt8, scoring::PairHMMScoring)
  return log(left_byte == right_byte ? Float64(scoring.match_prob) : Float64(scoring.mismatch_prob))
end

function _pairhmm_transition_logs(scoring::PairHMMScoring)
  go = Float64(scoring.gap_open_prob)
  ge = Float64(scoring.gap_extend_prob)
  return (
    mm=log1p(-2go), mix=log(go), miy=log(go),
    ixm=log1p(-ge), ixx=log(ge),
    iym=log1p(-ge), iyy=log(ge),
  )
end

function _pairhmm_forward(left::AbstractVector{UInt8}, right::AbstractVector{UInt8}, scoring::PairHMMScoring)
  _validate_pairhmm_scoring(scoring)
  n = length(left)
  m = length(right)
  neg = -Inf
  M = fill(neg, n + 1, m + 1)
  Ix = fill(neg, n + 1, m + 1)
  Iy = fill(neg, n + 1, m + 1)
  T = _pairhmm_transition_logs(scoring)
  M[1, 1] = 0.0
  @inbounds for i in 0:n
    for j in 0:m
      ii = i + 1
      jj = j + 1
      if i > 0 && j > 0
        emit = _pairhmm_emit_log(left[i], right[j], scoring)
        M[ii, jj] = emit + _logsumexp3(M[ii-1, jj-1] + T.mm, Ix[ii-1, jj-1] + T.ixm, Iy[ii-1, jj-1] + T.iym)
      end
      if i > 0
        Ix[ii, jj] = _logaddexp(M[ii-1, jj] + T.mix, Ix[ii-1, jj] + T.ixx)
      end
      if j > 0
        Iy[ii, jj] = _logaddexp(M[ii, jj-1] + T.miy, Iy[ii, jj-1] + T.iyy)
      end
    end
  end
  return M, Ix, Iy
end

function _pairhmm_backward(left::AbstractVector{UInt8}, right::AbstractVector{UInt8}, scoring::PairHMMScoring)
  _validate_pairhmm_scoring(scoring)
  n = length(left)
  m = length(right)
  neg = -Inf
  BM = fill(neg, n + 1, m + 1)
  BIx = fill(neg, n + 1, m + 1)
  BIy = fill(neg, n + 1, m + 1)
  T = _pairhmm_transition_logs(scoring)
  BM[n+1, m+1] = 0.0
  BIx[n+1, m+1] = 0.0
  BIy[n+1, m+1] = 0.0
  @inbounds for i in n:-1:0
    for j in m:-1:0
      ii = i + 1
      jj = j + 1
      if i == n && j == m
        continue
      end
      vals_m = Float64[]
      vals_ix = Float64[]
      vals_iy = Float64[]
      if i < n && j < m
        emit = _pairhmm_emit_log(left[i+1], right[j+1], scoring)
        push!(vals_m, T.mm + emit + BM[ii+1, jj+1])
        push!(vals_ix, T.ixm + emit + BM[ii+1, jj+1])
        push!(vals_iy, T.iym + emit + BM[ii+1, jj+1])
      end
      if i < n
        push!(vals_m, T.mix + BIx[ii+1, jj])
        push!(vals_ix, T.ixx + BIx[ii+1, jj])
      end
      if j < m
        push!(vals_m, T.miy + BIy[ii, jj+1])
        push!(vals_iy, T.iyy + BIy[ii, jj+1])
      end
      !isempty(vals_m) && (BM[ii, jj] = reduce(_logaddexp, vals_m))
      !isempty(vals_ix) && (BIx[ii, jj] = reduce(_logaddexp, vals_ix))
      !isempty(vals_iy) && (BIy[ii, jj] = reduce(_logaddexp, vals_iy))
    end
  end
  return BM, BIx, BIy
end

function _posterior_consensus_alignment(::Type{A}, left::AbstractVector{UInt8}, right::AbstractVector{UInt8}, posterior::Array{Float32,3}) where {A<:BioAlphabet}
  n = length(left)
  m = length(right)
  
  # Optimal Accuracy DP (Holmes 1998, Do et al. 2005)
  dp = zeros(Float64, n + 1, m + 1)
  trace = zeros(UInt8, n + 1, m + 1) # 1=Match, 2=Up(GapRight), 3=Left(GapLeft)
  
  @inbounds for i in 1:n
      for j in 1:m
          p_match = posterior[i+1, j+1, 1]
          
          match_val = dp[i, j] + p_match
          up_val = dp[i, j+1]
          left_val = dp[i+1, j]
          
          if match_val >= up_val && match_val >= left_val
              dp[i+1, j+1] = match_val
              trace[i+1, j+1] = 0x01
          elseif up_val >= left_val
              dp[i+1, j+1] = up_val
              trace[i+1, j+1] = 0x02
          else
              dp[i+1, j+1] = left_val
              trace[i+1, j+1] = 0x03
          end
      end
  end
  
  aligned_left = UInt8[]
  aligned_right = UInt8[]
  matches = 0
  
  i, j = n + 1, m + 1
  while i > 1 || j > 1
      if i > 1 && j > 1 && trace[i, j] == 0x01
          push!(aligned_left, left[i-1])
          push!(aligned_right, right[j-1])
          matches += left[i-1] == right[j-1] ? 1 : 0
          i -= 1; j -= 1
      elseif i > 1 && (j == 1 || trace[i, j] == 0x02)
          push!(aligned_left, left[i-1])
          push!(aligned_right, UInt8('-'))
          i -= 1
      else
          push!(aligned_left, UInt8('-'))
          push!(aligned_right, right[j-1])
          j -= 1
      end
  end
  
  reverse!(aligned_left)
  reverse!(aligned_right)
  
  aln_len = length(aligned_left)
  identity = aln_len == 0 ? 0.0 : matches / aln_len
  score = round(Int, dp[n+1, m+1] * 100)
  
  return PairwiseAlignmentResult(BioSequence{A}(aligned_left; validate=false), BioSequence{A}(aligned_right; validate=false), score, matches, identity)
end

function _pairhmm_posterior(::Type{A}, left::AbstractVector{UInt8}, right::AbstractVector{UInt8}, scoring::PairHMMScoring) where {A<:BioAlphabet}
  F = _pairhmm_forward(left, right, scoring)
  B = _pairhmm_backward(left, right, scoring)
  n = length(left)
  m = length(right)
  logz = _logsumexp3(F[1][n+1, m+1], F[2][n+1, m+1], F[3][n+1, m+1])
  posterior = zeros(Float32, n + 1, m + 1, 3)
  @inbounds for k in 1:3, i in 1:(n+1), j in 1:(m+1)
    value = F[k][i, j] + B[k][i, j] - logz
    posterior[i, j, k] = isfinite(value) ? Float32(clamp(exp(value), 0.0, 1.0)) : 0.0f0
  end
  expected_accuracy = mean(Float64.(posterior[:, :, 1]))
  consensus = _posterior_consensus_alignment(A, left, right, posterior)
  return PosteriorAlignmentResult(posterior, expected_accuracy, consensus, logz; metadata=Dict{Symbol,Any}(:method => :pairhmm_posterior))
end

"""
    pairhmm_align(left, right; scoring=PairHMMScoring())

Run global pair-HMM forward/backward alignment and return posterior state
probabilities plus a posterior-decoded consensus alignment.
"""
function pairhmm_align(left::BioSequence{A}, right::BioSequence{A}; scoring::PairHMMScoring=PairHMMScoring()) where {A<:BioAlphabet}
  result = _pairhmm_posterior(A, left.data, right.data, scoring)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, result, "pairhmm_align"; parents=provenance_parent_ids(left, right), parameters=(match_prob=Float64(scoring.match_prob), mismatch_prob=Float64(scoring.mismatch_prob), gap_open_prob=Float64(scoring.gap_open_prob), gap_extend_prob=Float64(scoring.gap_extend_prob)))
end

pairhmm_align(left::AbstractString, right::AbstractString; kwargs...) = pairhmm_align(DNASeq(left; validate=false), DNASeq(right; validate=false); kwargs...)

"""
    soft_alignment_score(left, right, scoring)

Differentiable global alignment objective using log-sum-exp dynamic programming.
The return value approaches the hard Needleman-Wunsch score as temperature goes
to zero.
"""
function soft_alignment_score(left::BioSequence{A}, right::BioSequence{A}, scoring::DifferentiableScoring) where {A<:BioAlphabet}
  return soft_alignment_score(left.data, right.data, scoring)
end

function soft_alignment_score(left::AbstractVector{UInt8}, right::AbstractVector{UInt8}, scoring::DifferentiableScoring)
  n = length(left)
  m = length(right)
  match_score, mismatch_score, gap_score = scoring.weights[1], scoring.weights[2], scoring.weights[3]
  τ = scoring.temperature
  dp = fill(-Inf, n + 1, m + 1)
  dp[1, 1] = 0.0
  for i in 1:n
    dp[i+1, 1] = dp[i, 1] + gap_score
  end
  for j in 1:m
    dp[1, j+1] = dp[1, j] + gap_score
  end
  @inbounds for i in 1:n, j in 1:m
    sub = left[i] == right[j] ? match_score : mismatch_score
    a = (dp[i, j] + sub) / τ
    b = (dp[i, j+1] + gap_score) / τ
    c = (dp[i+1, j] + gap_score) / τ
    dp[i+1, j+1] = τ * _logsumexp3(a, b, c)
  end
  return dp[n+1, m+1]
end

soft_alignment_score(left::AbstractString, right::AbstractString, scoring::DifferentiableScoring) = soft_alignment_score(codeunits(String(left)), codeunits(String(right)), scoring)




"""
    SequenceGraph{A}

Directed acyclic sequence graph for pangenome-style alignment. Nodes carry
sequence segments and edges encode allowed traversal between segments.
"""
struct SequenceGraph{A<:BioAlphabet} <: AbstractAnalysisResult
  nodes::Vector{BioSequence{A}}
  edges::Vector{Tuple{Int,Int}}
  metadata::Dict{Symbol,Any}
end

function SequenceGraph(nodes::AbstractVector{<:BioSequence{A}}, edges::AbstractVector{<:Tuple{Int,Int}}; metadata::AbstractDict=Dict{Symbol,Any}()) where {A<:BioAlphabet}
  n = length(nodes)
  for (u, v) in edges
    1 <= u <= n || throw(ArgumentError("SequenceGraph edge source $u is out of bounds"))
    1 <= v <= n || throw(ArgumentError("SequenceGraph edge target $v is out of bounds"))
    u == v && throw(ArgumentError("SequenceGraph must be acyclic; self edge at node $u"))
  end
  metadata_copy = Dict{Symbol,Any}(metadata)
  ensure_provenance_id!(metadata_copy)
  graph = SequenceGraph{A}(collect(nodes), collect(edges), metadata_copy)
  _sequence_graph_toposort(graph)  # validates DAG
  return graph
end

struct GraphAlignmentResult{A<:BioAlphabet} <: AbstractAnalysisResult
  graph_path::Vector{Int}
  graph_sequence::BioSequence{A}
  query_alignment::PairwiseAlignmentResult{A}
  node_scores::Vector{Float64}
  metadata::Dict{Symbol,Any}
end

analysis_result_fields(::Type{GraphAlignmentResult{A}}) where {A} = (:graph_path, :query_alignment, :node_scores)

function _sequence_graph_toposort(graph::SequenceGraph)
  n = length(graph.nodes)
  indeg = zeros(Int, n)
  outgoing = [Int[] for _ in 1:n]
  for (u, v) in graph.edges
    push!(outgoing[u], v)
    indeg[v] += 1
  end
  queue = [i for i in 1:n if indeg[i] == 0]
  order = Int[]
  head = 1
  while head <= length(queue)
    u = queue[head]
    head += 1
    push!(order, u)
    for v in outgoing[u]
      indeg[v] -= 1
      indeg[v] == 0 && push!(queue, v)
    end
  end
  length(order) == n || throw(ArgumentError("SequenceGraph edges contain a cycle"))
  return order, outgoing
end

function _best_graph_path(graph::SequenceGraph{A}; match::Int=1, mismatch::Int=-1, edge_penalty::Real=0.0) where {A<:BioAlphabet}
  order, outgoing = _sequence_graph_toposort(graph)
  n = length(graph.nodes)
  preds = [Int[] for _ in 1:n]
  for (u, v) in graph.edges
    push!(preds[v], u)
  end
  node_scores = [sum(b == UInt8('-') ? 0.0 : match for b in node.data) for node in graph.nodes]
  best = fill(-Inf, n)
  parent = fill(0, n)
  for u in order
    local_score = node_scores[u]
    if isempty(preds[u])
      best[u] = local_score
    else
      pred_scores = [best[pred] + Float64(edge_penalty) for pred in preds[u]]
      idx = argmax(pred_scores)
      best[u] = pred_scores[idx] + local_score
      parent[u] = preds[u][idx]
    end
  end
  finish = argmax(best)
  path = Int[]
  cur = finish
  while cur != 0
    pushfirst!(path, cur)
    cur = parent[cur]
  end
  bytes = UInt8[]
  for node_id in path
    append!(bytes, graph.nodes[node_id].data)
  end
  return path, BioSequence{A}(bytes; validate=false), best
end

"""
    align_to_graph(query, graph; kwargs...)

Align a query sequence to the best scoring path through a DAG sequence graph,
then run the existing pairwise aligner against that path. The result keeps the
selected graph path and per-node path scores for diagnostics.
"""
function align_to_graph(query::BioSequence{A}, graph::SequenceGraph{A}; match::Int=1, mismatch::Int=-1, gap::Int=-1, edge_penalty::Real=0.0, kwargs...) where {A<:BioAlphabet}
  path, graph_sequence, node_scores = _best_graph_path(graph; match=match, mismatch=mismatch, edge_penalty=edge_penalty)
  aln = pairwise_align(query, graph_sequence; match=match, mismatch=mismatch, gap=gap, kwargs...)
  metadata = Dict{Symbol,Any}(:method => :best_path_dag_dp, :edge_penalty => Float64(edge_penalty))
  ensure_provenance_id!(metadata)
  result = GraphAlignmentResult{A}(path, graph_sequence, aln, node_scores, metadata)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, result, "align_to_graph"; parents=provenance_parent_ids(query, graph), parameters=(node_count=length(graph.nodes), edge_count=length(graph.edges), path_length=length(path), match=match, mismatch=mismatch, gap=gap, edge_penalty=Float64(edge_penalty)))
end

align_to_graph(query::AbstractString, graph::SequenceGraph{DNAAlphabet}; kwargs...) = align_to_graph(DNASeq(query; validate=false), graph; kwargs...)

"""
    AlignmentProfileHMM{A}

Profile representation for profile-profile alignment. Emission matrices are
positions by alphabet symbols; transition probabilities are retained for
provenance and future richer HMM alignment.
"""
struct AlignmentProfileHMM{A<:BioAlphabet} <: AbstractAnalysisResult
  alphabet::Vector{UInt8}
  match_emissions::Matrix{Float32}
  insert_emissions::Matrix{Float32}
  transition_probs::Array{Float32,3}
  metadata::Dict{Symbol,Any}
end

function _profile_hmm_construct(::Type{A}, alphabet::AbstractVector{UInt8}, match_emissions::AbstractMatrix{<:Real}; insert_emissions=nothing, transition_probs=nothing, metadata::AbstractDict=Dict{Symbol,Any}()) where {A<:BioAlphabet}
  alpha = collect(UInt8.(alphabet))
  match = Float32.(match_emissions)
  size(match, 2) == length(alpha) || throw(DimensionMismatch("match_emissions columns must match alphabet length"))
  insert = insert_emissions === nothing ? fill(Float32(1 / length(alpha)), size(match)) : Float32.(insert_emissions)
  size(insert) == size(match) || throw(DimensionMismatch("insert_emissions must match match_emissions dimensions"))
  trans = transition_probs === nothing ? fill(Float32(1/3), size(match, 1), 3, 3) : Float32.(transition_probs)
  metadata_copy = Dict{Symbol,Any}(metadata)
  ensure_provenance_id!(metadata_copy)
  return AlignmentProfileHMM{A}(alpha, match, insert, trans, metadata_copy)
end

AlignmentProfileHMM(alphabet::AbstractVector{UInt8}, match_emissions::AbstractMatrix{<:Real}; bioalphabet::Type{A}=DNAAlphabet, kwargs...) where {A<:BioAlphabet} =
  _profile_hmm_construct(A, alphabet, match_emissions; kwargs...)

struct ProfileAlignmentResult <: AbstractAnalysisResult
  path::Vector{Tuple{Int,Int}}
  score::Float64
  posterior_confidence::Vector{Float64}
  metadata::Dict{Symbol,Any}
end

analysis_result_fields(::Type{ProfileAlignmentResult}) = (:path, :score, :posterior_confidence)

@inline function _profile_column_score(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
  return sum(Float64(a[i]) * Float64(b[i]) for i in eachindex(a))
end

"""
    align_profiles(profile_a, profile_b; gap=-0.1)

Align two profile HMMs using match-emission dot products and a global DP over
profile positions.
"""
function align_profiles(profile_a::AlignmentProfileHMM, profile_b::AlignmentProfileHMM; gap::Real=-0.1)
  profile_a.alphabet == profile_b.alphabet || throw(ArgumentError("profile alphabets must match"))
  n = size(profile_a.match_emissions, 1)
  m = size(profile_b.match_emissions, 1)
  dp = zeros(Float64, n + 1, m + 1)
  trace = zeros(UInt8, n + 1, m + 1)
  for i in 1:n
    dp[i+1, 1] = dp[i, 1] + Float64(gap)
    trace[i+1, 1] = 0x02
  end
  for j in 1:m
    dp[1, j+1] = dp[1, j] + Float64(gap)
    trace[1, j+1] = 0x03
  end
  for i in 1:n, j in 1:m
    diag = dp[i, j] + _profile_column_score(view(profile_a.match_emissions, i, :), view(profile_b.match_emissions, j, :))
    up = dp[i, j+1] + Float64(gap)
    left = dp[i+1, j] + Float64(gap)
    if diag >= up && diag >= left
      dp[i+1, j+1] = diag;
      trace[i+1, j+1] = 0x01
    elseif up >= left
      dp[i+1, j+1] = up;
      trace[i+1, j+1] = 0x02
    else
      dp[i+1, j+1] = left;
      trace[i+1, j+1] = 0x03
    end
  end
  path = Tuple{Int,Int}[]
  conf = Float64[]
  i = n + 1;
  j = m + 1
  while i > 1 || j > 1
    t = trace[i, j]
    if t == 0x01
      pushfirst!(path, (i - 1, j - 1))
      pushfirst!(conf, clamp(_profile_column_score(view(profile_a.match_emissions, i - 1, :), view(profile_b.match_emissions, j - 1, :)), 0.0, 1.0))
      i -= 1;
      j -= 1
    elseif t == 0x02
      pushfirst!(path, (i - 1, 0));
      pushfirst!(conf, 0.0);
      i -= 1
    else
      pushfirst!(path, (0, j - 1));
      pushfirst!(conf, 0.0);
      j -= 1
    end
  end
  metadata = Dict{Symbol,Any}(:method => :profile_profile_dp, :gap => Float64(gap))
  ensure_provenance_id!(metadata)
  result = ProfileAlignmentResult(path, dp[n+1, m+1], conf, metadata)
  _ctx = active_provenance_context()
  return provenance_result!(_ctx, result, "align_profiles"; parents=provenance_parent_ids(profile_a, profile_b), parameters=(length_a=n, length_b=m, gap=Float64(gap)))
end

"""
    LinearPairwiseScoring

Simple match/mismatch scoring model for pairwise alignment.
"""
struct LinearPairwiseScoring <: AbstractPairwiseScoring
  match::Int
  mismatch::Int
end

"""
    SubstitutionMatrix

General substitution-matrix container with alphabet lookup tables, score
storage. Unknown symbols are rejected during alignment instead of being scored
silently.
"""
struct SubstitutionMatrix
  alphabet::Vector{UInt8}
  scores::Matrix{Int}
  lookup::Dict{UInt8,Int}
  lookup_table::Vector{Int}
  default::Int
end

"""
    MatrixPairwiseScoring

Pairwise scoring wrapper around a `SubstitutionMatrix`.
"""
struct MatrixPairwiseScoring <: AbstractPairwiseScoring
  matrix::SubstitutionMatrix
end

"""
    CodonSubstitutionMatrix

Substitution matrix specialized for codon tokens encoded as packed byte values.
Unknown codon tokens are rejected during alignment instead of being scored
silently.
"""
struct CodonSubstitutionMatrix
  alphabet::Vector{UInt8}
  scores::Matrix{Int}
  lookup::Dict{UInt8,Int}
  lookup_table::Vector{Int}
  default::Int
end

"""
    CodonMatrixPairwiseScoring

Pairwise scoring wrapper around a `CodonSubstitutionMatrix`.
"""
struct CodonMatrixPairwiseScoring <: AbstractPairwiseScoring
  matrix::CodonSubstitutionMatrix
end

const _STANDARD_SUBSTITUTION_MATRIX_TEXT = Dict{String,String}(
  "BLOSUM62" => raw"""
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
  A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
  R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
  N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
  D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
  C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
  Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
  E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
  G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
  H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
  I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
  L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
  K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
  M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
  F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
  P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
  S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
  T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
  W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
  Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
  V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
  B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
  Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
  X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
  * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
  """,
  "BLOSUM80" => raw"""
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
  A  7 -3 -3 -3 -1 -2 -2  0 -3 -3 -3 -1 -2 -4 -1  2  0 -5 -4 -1 -3 -2 -1 -8
  R -3  9 -1 -3 -6  1 -1 -4  0 -5 -4  3 -3 -5 -3 -2 -2 -5 -4 -4 -2  0 -2 -8
  N -3 -1  9  2 -5  0 -1 -1  1 -6 -6  0 -4 -6 -4  1  0 -7 -4 -5  5 -1 -2 -8
  D -3 -3  2 10 -7 -1  2 -3 -2 -7 -7 -2 -6 -6 -3 -1 -2 -8 -6 -6  6  1 -3 -8
  C -1 -6 -5 -7 13 -5 -7 -6 -7 -2 -3 -6 -3 -4 -6 -2 -2 -5 -5 -2 -6 -7 -4 -8
  Q -2  1  0 -1 -5  9  3 -4  1 -5 -4  2 -1 -5 -3 -1 -1 -4 -3 -4 -1  5 -2 -8
  E -2 -1 -1  2 -7  3  8 -4  0 -6 -6  1 -4 -6 -2 -1 -2 -6 -5 -4  1  6 -2 -8
  G  0 -4 -1 -3 -6 -4 -4  9 -4 -7 -7 -3 -5 -6 -5 -1 -3 -6 -6 -6 -2 -4 -3 -8
  H -3  0  1 -2 -7  1  0 -4 12 -6 -5 -1 -4 -2 -4 -2 -3 -4  3 -5 -1  0 -2 -8
  I -3 -5 -6 -7 -2 -5 -6 -7 -6  7  2 -5  2 -1 -5 -4 -2 -5 -3  4 -6 -6 -2 -8
  L -3 -4 -6 -7 -3 -4 -6 -7 -5  2  6 -4  3  0 -5 -4 -3 -4 -2  1 -7 -5 -2 -8
  K -1  3  0 -2 -6  2  1 -3 -1 -5 -4  8 -3 -5 -2 -1 -1 -6 -4 -4 -1  1 -2 -8
  M -2 -3 -4 -6 -3 -1 -4 -5 -4  2  3 -3  9  0 -4 -3 -1 -3 -3  1 -5 -3 -2 -8
  F -4 -5 -6 -6 -4 -5 -6 -6 -2 -1  0 -5  0 10 -6 -4 -4  0  4 -2 -6 -6 -3 -8
  P -1 -3 -4 -3 -6 -3 -2 -5 -4 -5 -5 -2 -4 -6 12 -2 -3 -7 -6 -4 -4 -2 -3 -8
  S  2 -2  1 -1 -2 -1 -1 -1 -2 -4 -4 -1 -3 -4 -2  7  2 -6 -3 -3  0 -1 -1 -8
  T  0 -2  0 -2 -2 -1 -2 -3 -3 -2 -3 -1 -1 -4 -3  2  8 -5 -3  0 -1 -2 -1 -8
  W -5 -5 -7 -8 -5 -4 -6 -6 -4 -5 -4 -6 -3  0 -7 -6 -5 16  3 -5 -8 -5 -5 -8
  Y -4 -4 -4 -6 -5 -3 -5 -6  3 -3 -2 -4 -3  4 -6 -3 -3  3 11 -3 -5 -4 -3 -8
  V -1 -4 -5 -6 -2 -4 -4 -6 -5  4  1 -4  1 -2 -4 -3  0 -5 -3  7 -6 -4 -2 -8
  B -3 -2  5  6 -6 -1  1 -2 -1 -6 -7 -1 -5 -6 -4  0 -1 -8 -5 -6  6  0 -3 -8
  Z -2  0 -1  1 -7  5  6 -4  0 -6 -5  1 -3 -6 -2 -1 -2 -5 -4 -4  0  6 -1 -8
  X -1 -2 -2 -3 -4 -2 -2 -3 -2 -2 -2 -2 -2 -3 -3 -1 -1 -5 -3 -2 -3 -1 -2 -8
  * -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8  1
  """,
  "PAM250" => raw"""
      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
  A  2 -2  0  0 -2  0  0  1 -1 -1 -2 -1 -1 -3  1  1  1 -6 -3  0  0  0  0 -8
  R -2  6  0 -1 -4  1 -1 -3  2 -2 -3  3  0 -4  0  0 -1  2 -4 -2 -1  0 -1 -8
  N  0  0  2  2 -4  1  1  0  2 -2 -3  1 -2 -3  0  1  0 -4 -2 -2  2  1  0 -8
  D  0 -1  2  4 -5  2  3  1  1 -2 -4  0 -3 -6 -1  0  0 -7 -4 -2  3  3 -1 -8
  C -2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3  0 -2 -8  0 -2 -4 -5 -3 -8
  Q  0  1  1  2 -5  4  2 -1  3 -2 -2  1 -1 -5  0 -1 -1 -5 -4 -2  1  3 -1 -8
  E  0 -1  1  3 -5  2  4  0  1 -2 -3  0 -2 -5 -1  0  0 -7 -4 -2  3  3 -1 -8
  G  1 -3  0  1 -3 -1  0  5 -2 -3 -4 -2 -3 -5  0  1  0 -7 -5 -1  0  0 -1 -8
  H -1  2  2  1 -3  3  1 -2  6 -2 -2  0 -2 -2  0 -1 -1 -3  0 -2  1  2 -1 -8
  I -1 -2 -2 -2 -2 -2 -2 -3 -2  5  2 -2  2  1 -2 -1  0 -5 -1  4 -2 -2 -1 -8
  L -2 -3 -3 -4 -6 -2 -3 -4 -2  2  6 -3  4  2 -3 -3 -2 -2 -1  2 -3 -3 -1 -8
  K -1  3  1  0 -5  1  0 -2  0 -2 -3  5  0 -5 -1  0  0 -3 -4 -2  1  0 -1 -8
  M -1  0 -2 -3 -5 -1 -2 -3 -2  2  4  0  6  0 -2 -2 -1 -4 -2  2 -2 -2 -1 -8
  F -3 -4 -3 -6 -4 -5 -5 -5 -2  1  2 -5  0  9 -5 -3 -3  0  7 -1 -4 -5 -2 -8
  P  1  0  0 -1 -3  0 -1  0  0 -2 -3 -1 -2 -5  6  1  0 -6 -5 -1 -1  0 -1 -8
  S  1  0  1  0  0 -1  0  1 -1 -1 -3  0 -2 -3  1  2  1 -2 -3 -1  0  0  0 -8
  T  1 -1  0  0 -2 -1  0  0 -1  0 -2  0 -1 -3  0  1  3 -5 -3  0  0 -1  0 -8
  W -6  2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4  0 -6 -2 -5 17  0 -6 -5 -6 -4 -8
  Y -3 -4 -2 -4  0 -4 -4 -5  0 -1 -1 -4 -2  7 -5 -3 -3  0 10 -2 -3 -4 -2 -8
  V  0 -2 -2 -2 -2 -2 -2 -1 -2  4  2 -2  2 -1 -1 -1  0 -6 -2  4 -2 -2 -1 -8
  B  0 -1  2  3 -4  1  3  0  1 -2 -3  1 -2 -4 -1  0  0 -5 -3 -2  3  2 -1 -8
  Z  0  0  1  3 -5  3  3  0  2 -2 -3  0 -2 -5  0  0 -1 -6 -4 -2  2  3 -1 -8
  X  0 -1  0 -1 -3 -1 -1 -1 -1 -1 -1 -1 -1 -2 -1  0  0 -4 -2 -1 -1 -1 -1 -8
  * -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8  1
  """
)

include("substitution_matrices_data.jl")
merge!(_STANDARD_SUBSTITUTION_MATRIX_TEXT, _BIOPYTHON_NAMED_SUBSTITUTION_MATRIX_TEXT)

"""
    _normalize_substitution_matrix_name(name)

Normalize a matrix name by uppercasing it and removing spacing and separator
characters.
"""
@inline function _normalize_substitution_matrix_name(name::String)
  return replace(uppercase(strip(name)), r"[\s._-]+" => "")
end

"""
    _parse_standard_substitution_matrix(spec; threaded=true)

Parse a standard amino-acid substitution matrix from its textual representation.
"""
function _parse_standard_substitution_matrix(spec::String; threaded::Bool=true)
  lines = filter(!isempty, Base.split(chomp(spec), '\n'))
  length(lines) >= 2 || throw(ArgumentError("substitution matrix spec must include a header and at least one row"))

  header = Base.split(strip(lines[1]))
  alphabet = join(header)
  row_count = length(header)
  scores = Matrix{Int}(undef, row_count, row_count)

  if threaded && row_count > 1 && Threads.nthreads() > 1
    Threads.@threads for row_index in 1:row_count
      fields = Base.split(strip(lines[row_index+1]))
      length(fields) == row_count + 1 || throw(ArgumentError("substitution matrix row has the wrong width"))
      fields[1] == header[row_index] || throw(ArgumentError("substitution matrix row order does not match its header"))
      @inbounds for col_index in 1:row_count
        scores[row_index, col_index] = parse(Int, fields[col_index+1])
      end
    end
  else
    for (row_index, line) in enumerate(@view lines[2:end])
      fields = Base.split(strip(line))
      length(fields) == row_count + 1 || throw(ArgumentError("substitution matrix row has the wrong width"))
      fields[1] == header[row_index] || throw(ArgumentError("substitution matrix row order does not match its header"))
      @inbounds for col_index in 1:row_count
        scores[row_index, col_index] = parse(Int, fields[col_index+1])
      end
    end
  end

  alphabet_bytes = collect(codeunits(alphabet))
  lookup_table = fill(0, 256)
  lookup = Dict{UInt8,Int}()

  for (index, byte) in enumerate(alphabet_bytes)
    _register_casefold_lookup!(lookup_table, lookup, byte, index)
  end

  return SubstitutionMatrix(alphabet_bytes, scores, lookup, lookup_table, minimum(scores))
end

@inline function _register_casefold_lookup!(lookup_table::Vector{Int}, lookup::Dict{UInt8,Int}, byte::UInt8, index::Int)
  lookup_table[Int(byte)+1] = index
  lookup[byte] = index
  if 0x41 <= byte <= 0x5a
    lower = byte + 0x20
    lookup_table[Int(lower)+1] = index
    lookup[lower] = index
  elseif 0x61 <= byte <= 0x7a
    upper = byte - 0x20
    lookup_table[Int(upper)+1] = index
    lookup[upper] = index
  end
end

const _STANDARD_SUBSTITUTION_MATRIX_CACHE = let cache = Dict{String,SubstitutionMatrix}()
  for (matrix_name, matrix_spec) in _STANDARD_SUBSTITUTION_MATRIX_TEXT
    normalized_name = replace(uppercase(strip(matrix_name)), r"[\s._-]+" => "")
    try
      cache[normalized_name] = _parse_standard_substitution_matrix(matrix_spec; threaded=false)
    catch err
      err isa ArgumentError || rethrow()
    end
  end
  cache
end

const _CODON_DECODE_TABLE = let table = Vector{String}(undef, 64)
  for first in 0:3
    for second in 0:3
      for third in 0:3
        code = (first << 4) | (second << 2) | third
        table[code+1] = String(UInt8[_kmer_base_from_code(first), _kmer_base_from_code(second), _kmer_base_from_code(third)])
      end
    end
  end
  table
end

"""
    _codon_code(byte)

Map a nucleotide byte to its compact DNA code.
"""
@inline function _codon_code(byte::UInt8)
  return _DNA_CODE[Int(byte)+1]
end

"""
    _pack_codon_token(first, second, third)

Pack three nucleotide bytes into a single codon token.
"""
@inline function _pack_codon_token(first::UInt8, second::UInt8, third::UInt8)
  code1 = _codon_code(first)
  code2 = _codon_code(second)
  code3 = _codon_code(third)
  return code1 > 3 || code2 > 3 || code3 > 3 ? UInt8(255) : UInt8((Int(code1) << 4) | (Int(code2) << 2) | Int(code3))
end

"""
    _encode_codon_sequence(bytes)

Convert a nucleotide byte sequence into packed codon tokens.
"""
function _encode_codon_sequence(bytes::AbstractVector{UInt8})
  length(bytes) % 3 == 0 || throw(ArgumentError("codon sequence length must be divisible by 3"))
  token_count = length(bytes) ÷ 3
  tokens = Vector{UInt8}(undef, token_count)

  @inbounds for token_index in 1:token_count
    byte_index = (token_index - 1) * 3 + 1
    tokens[token_index] = _pack_codon_token(bytes[byte_index], bytes[byte_index+1], bytes[byte_index+2])
  end

  return tokens
end

"""
    _decode_codon_token(token)

Convert a packed codon token back into a three-character codon string.
"""
function _decode_codon_token(token::UInt8)
  return token > 63 ? "NNN" : _CODON_DECODE_TABLE[Int(token)+1]
end

"""
    _parse_codon_substitution_matrix(spec; threaded=true)

Parse a codon substitution matrix from its textual table form.
"""
function _parse_codon_substitution_matrix(spec::String; threaded::Bool=true)
  lines = filter(!isempty, Base.split(chomp(spec), '\n'))
  length(lines) >= 2 || throw(ArgumentError("substitution matrix spec must include a header and at least one row"))

  header = Base.split(strip(lines[1]))
  length(header) % 3 == 0 || throw(ArgumentError("codon substitution matrix header must contain codon triplets"))
  row_count = length(header) ÷ 3
  row_count == 64 || throw(ArgumentError("codon substitution matrix must contain exactly 64 codons"))

  alphabet = Vector{UInt8}(undef, row_count)
  lookup_table = fill(0, 64)
  lookup = Dict{UInt8,Int}()

  for index in 1:row_count
    token_start = (index - 1) * 3 + 1
    token1 = header[token_start]
    token2 = header[token_start+1]
    token3 = header[token_start+2]
    token = join((token1, token2, token3), " ")
    code = _pack_codon_token(UInt8(token1[1]), UInt8(token2[1]), UInt8(token3[1]))
    code == 255 && throw(ArgumentError("invalid codon token in matrix header: $token"))
    alphabet[index] = code
    lookup[code] = index
    lookup_table[Int(code)+1] = index
  end

  scores = Matrix{Int}(undef, row_count, row_count)

  if threaded && row_count > 1 && Threads.nthreads() > 1
    Threads.@threads for row_index in 1:row_count
      fields = Base.split(strip(lines[row_index+1]))
      length(fields) == row_count + 1 || throw(ArgumentError("substitution matrix row has the wrong width"))
      @inbounds for col_index in 1:row_count
        scores[row_index, col_index] = parse(Int, fields[col_index+1])
      end
    end
  else
    for (row_index, line) in enumerate(@view lines[2:end])
      fields = Base.split(strip(line))
      length(fields) == row_count + 1 || throw(ArgumentError("substitution matrix row has the wrong width"))
      @inbounds for col_index in 1:row_count
        scores[row_index, col_index] = parse(Int, fields[col_index+1])
      end
    end
  end

  return CodonSubstitutionMatrix(alphabet, scores, lookup, lookup_table, minimum(scores))
end

const _STANDARD_CODON_SUBSTITUTION_MATRIX_CACHE = Dict{String,CodonSubstitutionMatrix}(
  "SCHNEIDER" => _parse_codon_substitution_matrix(_STANDARD_SUBSTITUTION_MATRIX_TEXT["SCHNEIDER"]; threaded=false),
)

const _PAIRWISE_STATE_MATCH = UInt8(0x01)
const _PAIRWISE_STATE_GAP_LEFT = UInt8(0x02)
const _PAIRWISE_STATE_GAP_RIGHT = UInt8(0x03)
const _PAIRWISE_TRACE_NONE = UInt8(0x00)
const _PAIRWISE_AFFINE_NEGATIVE_INFINITY = typemin(Int) ÷ 4

@inline function _pairwise_affine_transition(score::Int, penalty::Int)
  return score == _PAIRWISE_AFFINE_NEGATIVE_INFINITY ? _PAIRWISE_AFFINE_NEGATIVE_INFINITY : score + penalty
end

"""
    Base.show(io, result)

Render a compact summary of a `PairwiseAlignmentResult`.
"""
function Base.show(io::IO, result::PairwiseAlignmentResult)
  print(io, analysis_result_summary(result), " | ", container_provenance_summary(result))
end

function Base.show(io::IO, ::MIME"text/plain", result::PairwiseAlignmentResult)
  print(io, analysis_result_summary(result), " | ", container_provenance_summary(result))
end

"""
    SubstitutionMatrix(alphabet; match=1, mismatch=-1, default=mismatch, threaded=true)

Build a simple square substitution matrix from an alphabet string.

ASCII letters are registered in both cases, so lowercase sequence input does not
silently fall through to the fallback score.
"""
function SubstitutionMatrix(alphabet::String; match::Int=1, mismatch::Int=-1, default::Int=mismatch, threaded::Bool=true)
  alphabet_bytes = collect(codeunits(alphabet))
  scores = Matrix{Int}(undef, length(alphabet_bytes), length(alphabet_bytes))
  lookup_table = fill(0, 256)

  if threaded && length(alphabet_bytes) > 1 && Threads.nthreads() > 1
    Threads.@threads for i in eachindex(alphabet_bytes)
      lookup_table[Int(alphabet_bytes[i])+1] = i
      @inbounds for j in eachindex(alphabet_bytes)
        scores[i, j] = alphabet_bytes[i] == alphabet_bytes[j] ? match : mismatch
      end
    end
  else
    @inbounds for i in eachindex(alphabet_bytes)
      lookup_table[Int(alphabet_bytes[i])+1] = i
      for j in eachindex(alphabet_bytes)
        scores[i, j] = alphabet_bytes[i] == alphabet_bytes[j] ? match : mismatch
      end
    end
  end

  lookup = Dict{UInt8,Int}()
  @inbounds for (index, byte) in pairs(alphabet_bytes)
    _register_casefold_lookup!(lookup_table, lookup, byte, index)
  end
  return SubstitutionMatrix(alphabet_bytes, scores, lookup, lookup_table, default)
end

"""
    SubstitutionMatrix(alphabet, scores; default=0, threaded=true)

Wrap an explicit square score matrix for the given alphabet.
"""
function SubstitutionMatrix(alphabet::String, scores::AbstractMatrix{<:Integer}; default::Int=0, threaded::Bool=true)
  alphabet_bytes = collect(codeunits(alphabet))
  size(scores, 1) == length(alphabet_bytes) || throw(ArgumentError("score matrix row count must match alphabet length"))
  size(scores, 2) == length(alphabet_bytes) || throw(ArgumentError("score matrix column count must match alphabet length"))
  lookup_table = fill(0, 256)
  if threaded && length(alphabet_bytes) > 1 && Threads.nthreads() > 1
    Threads.@threads for index in eachindex(alphabet_bytes)
      byte = alphabet_bytes[index]
      lookup_table[Int(byte)+1] = index
    end
  else
    @inbounds for (index, byte) in pairs(alphabet_bytes)
      lookup_table[Int(byte)+1] = index
    end
  end
  lookup = Dict{UInt8,Int}()
  @inbounds for (index, byte) in pairs(alphabet_bytes)
    _register_casefold_lookup!(lookup_table, lookup, byte, index)
  end
  return SubstitutionMatrix(alphabet_bytes, Matrix{Int}(scores), lookup, lookup_table, default)
end

"""
    substitution_matrix(name_or_alphabet; kwargs...)

Load a named substitution matrix (e.g. "BLOSUM62", "FOLDSEEK3DI") or build a matrix for an alphabet string.
"""
function substitution_matrix(name_or_alphabet::String; kwargs...)
    norm = _normalize_substitution_matrix_name(name_or_alphabet)
    if haskey(_STANDARD_SUBSTITUTION_MATRIX_CACHE, norm) || haskey(_STANDARD_SUBSTITUTION_MATRIX_TEXT, name_or_alphabet) || haskey(_STANDARD_SUBSTITUTION_MATRIX_TEXT, norm) || norm in ("DNA", "NUC", "NUCLEOTIDE", "IDENTITY")
        return named_substitution_matrix(name_or_alphabet)
    end
    return SubstitutionMatrix(name_or_alphabet; kwargs...)
end

"""
    substitution_matrix(alphabet, scores; kwargs...)

Convenience alias for `SubstitutionMatrix(alphabet, scores; kwargs...)`.
"""
substitution_matrix(alphabet::String, scores::AbstractMatrix{<:Integer}; kwargs...) = SubstitutionMatrix(alphabet, scores; kwargs...)

"""
    named_substitution_matrix(name)

Load one of the built-in named substitution matrices or a DNA identity matrix.
"""
function named_substitution_matrix(name::String)
  normalized = _normalize_substitution_matrix_name(name)
  if normalized == "DNA" || normalized == "NUC" || normalized == "NUCLEOTIDE" || normalized == "IDENTITY"
    return SubstitutionMatrix("ACGT"; match=1, mismatch=-1, default=-1)
  end

  matrix = get(_STANDARD_SUBSTITUTION_MATRIX_CACHE, normalized, nothing)
  if matrix !== nothing
    return matrix
  end

  spec = get(_STANDARD_SUBSTITUTION_MATRIX_TEXT, String(name), nothing)
  spec === nothing && (spec = get(_STANDARD_SUBSTITUTION_MATRIX_TEXT, normalized, nothing))
  spec === nothing && throw(ArgumentError("unknown named substitution matrix: $name"))

  parsed = _parse_standard_substitution_matrix(spec)
  _STANDARD_SUBSTITUTION_MATRIX_CACHE[normalized] = parsed
  _ctx = active_provenance_context()
  _ctx !== nothing && register_provenance!(_ctx, "named_substitution_matrix";
    parameters=(name=name,))
  return parsed
end

"""
    named_substitution_matrix(name::Symbol)

Symbol-based wrapper for `named_substitution_matrix`.
"""
named_substitution_matrix(name::Symbol) = named_substitution_matrix(String(name))

"""
    substitution_matrix(name::Symbol)

Symbol-based alias for `named_substitution_matrix`.
"""
substitution_matrix(name::Symbol) = named_substitution_matrix(name)

"""
    available_named_substitution_matrices()

List the built-in named substitution matrices available in this module.
"""
available_named_substitution_matrices() = sort!(collect(keys(_STANDARD_SUBSTITUTION_MATRIX_TEXT)))

"""
    named_codon_substitution_matrix(name)

Load a built-in codon substitution matrix by name.
"""
function named_codon_substitution_matrix(name::String)
  normalized = _normalize_substitution_matrix_name(name)
  normalized == "SCHNEIDER" || throw(ArgumentError("unknown named codon substitution matrix: $name"))

  matrix = get(_STANDARD_CODON_SUBSTITUTION_MATRIX_CACHE, normalized, nothing)
  matrix === nothing && throw(ArgumentError("unknown named codon substitution matrix: $name"))
  _ctx = active_provenance_context()
  _ctx !== nothing && register_provenance!(_ctx, "named_codon_substitution_matrix";
    parameters=(name=name,))
  return matrix
end

"""
    named_codon_substitution_matrix(name::Symbol)

Symbol-based wrapper for `named_codon_substitution_matrix`.
"""
named_codon_substitution_matrix(name::Symbol) = named_codon_substitution_matrix(String(name))

"""
    codon_substitution_matrix(alphabet, scores; default=0, threaded=true)

Construct a codon-aware substitution matrix from an explicit codon alphabet and
score matrix.
"""
function codon_substitution_matrix(alphabet::String, scores::AbstractMatrix{<:Integer}; default::Int=0, threaded::Bool=true)
  codon_tokens = Base.split(strip(alphabet))
  length(codon_tokens) == 64 || throw(ArgumentError("codon alphabet must contain exactly 64 codons"))

  alphabet_codes = Vector{UInt8}(undef, 64)
  lookup_table = fill(0, 64)

  if threaded && length(codon_tokens) > 1 && Threads.nthreads() > 1
    Threads.@threads for index in eachindex(codon_tokens)
      token = codon_tokens[index]
      length(token) == 3 || throw(ArgumentError("codon alphabet entries must be triplets"))
      code = _pack_codon_token(UInt8(token[1]), UInt8(token[2]), UInt8(token[3]))
      code == 255 && throw(ArgumentError("invalid codon token in alphabet: $token"))
      alphabet_codes[index] = code
      lookup_table[Int(code)+1] = index
    end
  else
    for (index, token) in enumerate(codon_tokens)
      length(token) == 3 || throw(ArgumentError("codon alphabet entries must be triplets"))
      code = _pack_codon_token(UInt8(token[1]), UInt8(token[2]), UInt8(token[3]))
      code == 255 && throw(ArgumentError("invalid codon token in alphabet: $token"))
      alphabet_codes[index] = code
      lookup_table[Int(code)+1] = index
    end
  end

  lookup = Dict{UInt8,Int}(code => index for (index, code) in pairs(alphabet_codes))

  size(scores, 1) == 64 || throw(ArgumentError("codon score matrix row count must be 64"))
  size(scores, 2) == 64 || throw(ArgumentError("codon score matrix column count must be 64"))
  _ctx = active_provenance_context()
  _ctx !== nothing && register_provenance!(_ctx, "codon_substitution_matrix";
    parameters=(n_codons=64, default=default))
  return CodonSubstitutionMatrix(alphabet_codes, Matrix{Int}(scores), lookup, lookup_table, default)
end

"""
    _pairwise_score(scoring, left_byte, right_byte)

Score a nucleotide pair under linear match/mismatch scoring.
"""
@inline function _pairwise_score(scoring::LinearPairwiseScoring, left_byte::UInt8, right_byte::UInt8)
  return left_byte == right_byte ? scoring.match : scoring.mismatch
end

"""
    _pairwise_score(scoring, left_byte, right_byte)

Score a nucleotide pair using a substitution matrix lookup.
"""
@inline function _pairwise_score(scoring::MatrixPairwiseScoring, left_byte::UInt8, right_byte::UInt8)
  left_index = scoring.matrix.lookup_table[Int(left_byte)+1]
  right_index = scoring.matrix.lookup_table[Int(right_byte)+1]
  left_index == 0 && throw(ArgumentError("unknown symbol in substitution-matrix alignment: '$(Char(left_byte))' (code $(Int(left_byte)))"))
  right_index == 0 && throw(ArgumentError("unknown symbol in substitution-matrix alignment: '$(Char(right_byte))' (code $(Int(right_byte)))"))
  return scoring.matrix.scores[left_index, right_index]
end

"""
    _pairwise_score(scoring, left_token, right_token)

Score a codon pair using a codon substitution matrix lookup.
"""
@inline function _pairwise_score(scoring::CodonMatrixPairwiseScoring, left_token::UInt8, right_token::UInt8)
  left_index = left_token > 63 ? 0 : scoring.matrix.lookup_table[Int(left_token)+1]
  right_index = right_token > 63 ? 0 : scoring.matrix.lookup_table[Int(right_token)+1]
  left_index == 0 && throw(ArgumentError("unknown codon token in substitution-matrix alignment: $(Int(left_token))"))
  right_index == 0 && throw(ArgumentError("unknown codon token in substitution-matrix alignment: $(Int(right_token))"))
  return scoring.matrix.scores[left_index, right_index]
end

"""
    _pairwise_gap_parameters(gap, gap_open, gap_extend)

Resolve linear and affine gap penalties into explicit open and extend values.
"""
@inline function _pairwise_gap_parameters(gap::Int, gap_open::Union{Nothing,Int}, gap_extend::Union{Nothing,Int})
  if gap_open === nothing && gap_extend === nothing
    gap < 0 || throw(ArgumentError("gap penalty must be negative"))
    return gap, gap
  end

  open_score = gap_open === nothing ? something(gap_extend, gap) : gap_open
  extend_score = gap_extend === nothing ? something(gap_open, gap) : gap_extend
  (open_score < 0 && extend_score < 0) || throw(ArgumentError("gap penalties must be negative"))
  return open_score, extend_score
end

"""
    pairwise_align(left_seq, right_seq; kwargs...)

Align two `BioSequence{A}` objects using linear or affine gap scoring.
Returns a `PairwiseAlignmentResult{A}`.
"""
function pairwise_align(
  left::BioSequence{A},
  right::BioSequence{A};
  kwargs...
) where {A<:BioAlphabet}
  result = _pairwise_align_dispatch(A, left.data, right.data; kwargs...)
  return _pairwise_align_finalize(result;
    parents=provenance_parent_ids(left, right), parameters=kwargs)
end

"""
    pairwise_align(left_bytes, right_bytes; kwargs...)

Align two byte sequences using linear or alphabet-less scoring.
Returns a `PairwiseAlignmentResult{DNAAlphabet}` (defaulting to DNA for raw bytes).
"""
function pairwise_align(
  left_bytes::AbstractVector{UInt8},
  right_bytes::AbstractVector{UInt8};
  kwargs...
)
  result = _pairwise_align_dispatch(DNAAlphabet, left_bytes, right_bytes; kwargs...)
  return _pairwise_align_finalize(result;
    parents=provenance_parent_ids(left_bytes, right_bytes), parameters=kwargs)
end

# Removed pairwise_align(::AbstractString, ::AbstractString) - use BioSequence types instead
function pairwise_align(left::AbstractString, right::AbstractString; kwargs...)
  return pairwise_align(DNASeq(left; validate=false), DNASeq(right; validate=false); kwargs...)
end

@inline function _pairwise_align_finalize(result::PairwiseAlignmentResult;
  parents::AbstractVector{<:AbstractString},
  parameters)
  _ctx = active_provenance_context()
  _ctx !== nothing && register_provenance!(_ctx, "pairwise_align";
    parents=parents,
    parameters=(score=result.score, identity=round(result.identity; digits=4),
      n_left=length(result.left), n_right=length(result.right)))
  return provenance_result!(_ctx, result, "pairwise_align";
    parents=parents, parameters=parameters)
end

function _pairwise_align_dispatch(
  ::Type{A},
  left_bytes::AbstractVector{UInt8},
  right_bytes::AbstractVector{UInt8};
  match::Int=1,
  mismatch::Int=-1,
  substitution_matrix::Union{Nothing,SubstitutionMatrix}=nothing,
  gap::Int=-1,
  gap_open::Union{Nothing,Int}=nothing,
  gap_extend::Union{Nothing,Int}=nothing,
  is_local::Bool=false,
) where {A<:BioAlphabet}
  scoring = substitution_matrix === nothing ? LinearPairwiseScoring(match, mismatch) : MatrixPairwiseScoring(substitution_matrix)
  if gap_open === nothing && gap_extend === nothing
    gap < 0 || throw(ArgumentError("gap penalty must be negative"))
    return is_local ? _pairwise_align_local(A, left_bytes, right_bytes, scoring, gap) : _pairwise_align_global(A, left_bytes, right_bytes, scoring, gap)
  end

  affine_gap_open, affine_gap_extend = _pairwise_gap_parameters(gap, gap_open, gap_extend)
  return is_local ? _pairwise_align_affine_local(A, left_bytes, right_bytes, scoring, affine_gap_open, affine_gap_extend) : _pairwise_align_affine_global(A, left_bytes, right_bytes, scoring, affine_gap_open, affine_gap_extend)
end

function pairwise_align(left::SeqRecordLite, right::SeqRecordLite; kwargs...)
  result = pairwise_align(left.sequence, right.sequence; kwargs...)
  return _pairwise_align_finalize(result;
    parents=provenance_parent_ids(left, right), parameters=kwargs)
end

function pairwise_align(left::FastqRecord, right::FastqRecord; kwargs...)
  result = pairwise_align(left.sequence, right.sequence; kwargs...)
  return _pairwise_align_finalize(result;
    parents=provenance_parent_ids(left, right), parameters=kwargs)
end

"""
    _write_codon_token!(buffer, write_index, codon)

Write a three-character codon into a preallocated byte buffer from the back.
"""
@inline function _write_codon_token!(buffer::Vector{UInt8}, write_index::Int, codon::String)
  codon_bytes = codeunits(codon)
  buffer[write_index-2] = codon_bytes[1]
  buffer[write_index-1] = codon_bytes[2]
  buffer[write_index] = codon_bytes[3]
  return write_index - 3
end

"""
    _pairwise_traceback_codon(...)

Reconstruct a codon alignment from a standard traceback matrix.
"""
function _pairwise_traceback_codon(
  left_tokens::AbstractVector{UInt8},
  right_tokens::AbstractVector{UInt8},
  scores::Union{Nothing,Matrix{Int}},
  trace::Matrix{UInt8},
  i::Int,
  j::Int,
  best_score::Int,
)
  aligned_left = Vector{UInt8}(undef, 3 * (length(left_tokens) + length(right_tokens)))
  aligned_right = Vector{UInt8}(undef, 3 * (length(left_tokens) + length(right_tokens)))
  left_write_index = length(aligned_left)
  right_write_index = length(aligned_right)
  aligned_length = 0
  matches = 0

  current_i = i
  current_j = j

  while current_i > 1 || current_j > 1
    if scores !== nothing && scores[current_i, current_j] == 0 && trace[current_i, current_j] == 0x00
      break
    end

    direction = trace[current_i, current_j]
    direction == 0x00 && break

    if direction == 0x01
      left_token = left_tokens[current_i-1]
      right_token = right_tokens[current_j-1]
      left_write_index = _write_codon_token!(aligned_left, left_write_index, _decode_codon_token(left_token))
      right_write_index = _write_codon_token!(aligned_right, right_write_index, _decode_codon_token(right_token))
      matches += left_token == right_token ? 1 : 0
      current_i -= 1
      current_j -= 1
    elseif direction == 0x02
      left_write_index = _write_codon_token!(aligned_left, left_write_index, _decode_codon_token(left_tokens[current_i-1]))
      right_write_index = _write_codon_token!(aligned_right, right_write_index, "---")
      current_i -= 1
    else
      left_write_index = _write_codon_token!(aligned_left, left_write_index, "---")
      right_write_index = _write_codon_token!(aligned_right, right_write_index, _decode_codon_token(right_tokens[current_j-1]))
      current_j -= 1
    end

    aligned_length += 1
  end

  left_start_index = left_write_index + 1
  right_start_index = right_write_index + 1
  final_left = BioSequence{DNAAlphabet}(aligned_left[left_start_index:end]; validate=false)
  final_right = BioSequence{DNAAlphabet}(aligned_right[right_start_index:end]; validate=false)
  identity = aligned_length == 0 ? 0.0 : matches / aligned_length

  return PairwiseAlignmentResult(final_left, final_right, best_score, matches, identity)
end

"""
    _pairwise_traceback_codon_affine(...)

Reconstruct a codon alignment from affine-gap traceback matrices.
"""
function _pairwise_traceback_codon_affine(
  left_tokens::AbstractVector{UInt8},
  right_tokens::AbstractVector{UInt8},
  match_trace::Matrix{UInt8},
  gap_left_trace::Matrix{UInt8},
  gap_right_trace::Matrix{UInt8},
  i::Int,
  j::Int,
  start_state::UInt8,
  best_score::Int,
)
  aligned_left = Vector{UInt8}(undef, 3 * (length(left_tokens) + length(right_tokens)))
  aligned_right = Vector{UInt8}(undef, 3 * (length(left_tokens) + length(right_tokens)))
  left_write_index = length(aligned_left)
  right_write_index = length(aligned_right)
  aligned_length = 0
  matches = 0

  current_i = i
  current_j = j
  current_state = start_state

  while current_state != _PAIRWISE_TRACE_NONE && (current_i > 1 || current_j > 1)
    if current_state == _PAIRWISE_STATE_MATCH
      left_token = left_tokens[current_i-1]
      right_token = right_tokens[current_j-1]
      left_write_index = _write_codon_token!(aligned_left, left_write_index, _decode_codon_token(left_token))
      right_write_index = _write_codon_token!(aligned_right, right_write_index, _decode_codon_token(right_token))
      matches += left_token == right_token ? 1 : 0
      current_state = match_trace[current_i, current_j]
      current_i -= 1
      current_j -= 1
    elseif current_state == _PAIRWISE_STATE_GAP_LEFT
      left_write_index = _write_codon_token!(aligned_left, left_write_index, _decode_codon_token(left_tokens[current_i-1]))
      right_write_index = _write_codon_token!(aligned_right, right_write_index, "---")
      current_state = gap_left_trace[current_i, current_j]
      current_i -= 1
    else
      left_write_index = _write_codon_token!(aligned_left, left_write_index, "---")
      right_write_index = _write_codon_token!(aligned_right, right_write_index, _decode_codon_token(right_tokens[current_j-1]))
      current_state = gap_right_trace[current_i, current_j]
      current_j -= 1
    end

    aligned_length += 1
  end

  left_start_index = left_write_index + 1
  right_start_index = right_write_index + 1
  final_left = BioSequence{DNAAlphabet}(aligned_left[left_start_index:end]; validate=false)
  final_right = BioSequence{DNAAlphabet}(aligned_right[right_start_index:end]; validate=false)
  identity = aligned_length == 0 ? 0.0 : matches / aligned_length

  return PairwiseAlignmentResult(final_left, final_right, best_score, matches, identity)
end

"""
    _pairwise_align_codon_global(left_tokens, right_tokens, scoring, gap)

Run global codon alignment with a linear gap model.
"""
function _pairwise_align_codon_global(left_tokens::AbstractVector{UInt8}, right_tokens::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap::Int)
  left_length = length(left_tokens)
  right_length = length(right_tokens)

  previous_scores = Vector{Int}(undef, left_length + 1)
  current_scores = similar(previous_scores)
  trace = Matrix{UInt8}(undef, left_length + 1, right_length + 1)

  @inbounds begin
    previous_scores[1] = 0
    trace[1, 1] = 0x00

    for i in 2:(left_length+1)
      previous_scores[i] = previous_scores[i-1] + gap
      trace[i, 1] = 0x02
    end

    for j in 2:(right_length+1)
      current_scores[1] = previous_scores[1] + gap
      trace[1, j] = 0x03
      right_token = right_tokens[j-1]
      for i in 2:(left_length+1)
        left_token = left_tokens[i-1]

        diag_score = previous_scores[i-1] + _pairwise_score(scoring, left_token, right_token)
        vertical_score = current_scores[i-1] + gap
        horizontal_score = previous_scores[i] + gap

        cell_score = diag_score
        cell_trace = 0x01

        if vertical_score > cell_score
          cell_score = vertical_score
          cell_trace = 0x02
        end
        if horizontal_score > cell_score
          cell_score = horizontal_score
          cell_trace = 0x03
        end

        current_scores[i] = cell_score
        trace[i, j] = cell_trace
      end

      previous_scores, current_scores = current_scores, previous_scores
    end
  end

  return _pairwise_traceback_codon(left_tokens, right_tokens, nothing, trace, left_length + 1, right_length + 1, previous_scores[left_length+1])
end

"""
    _pairwise_align_codon_local(left_tokens, right_tokens, scoring, gap)

Run local codon alignment with a linear gap model.
"""
function _pairwise_align_codon_local(left_tokens::AbstractVector{UInt8}, right_tokens::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap::Int)
  left_length = length(left_tokens)
  right_length = length(right_tokens)

  scores = Matrix{Int}(undef, left_length + 1, right_length + 1)
  trace = Matrix{UInt8}(undef, left_length + 1, right_length + 1)

  best_score = 0
  best_i = 1
  best_j = 1

  @inbounds begin
    scores[1, 1] = 0
    trace[1, 1] = 0x00

    for i in 2:(left_length+1)
      scores[i, 1] = 0
      trace[i, 1] = 0x00
    end

    for j in 2:(right_length+1)
      scores[1, j] = 0
      trace[1, j] = 0x00
    end

    for j in 2:(right_length+1)
      right_token = right_tokens[j-1]
      for i in 2:(left_length+1)
        left_token = left_tokens[i-1]

        diag_score = scores[i-1, j-1] + _pairwise_score(scoring, left_token, right_token)
        up_score = scores[i-1, j] + gap
        left_score = scores[i, j-1] + gap

        cell_score = 0
        cell_trace = 0x00

        if diag_score > cell_score
          cell_score = diag_score
          cell_trace = 0x01
        end
        if up_score > cell_score
          cell_score = up_score
          cell_trace = 0x02
        end
        if left_score > cell_score
          cell_score = left_score
          cell_trace = 0x03
        end

        scores[i, j] = cell_score
        trace[i, j] = cell_trace

        if cell_score > best_score
          best_score = cell_score
          best_i = i
          best_j = j
        end
      end
    end
  end

  return _pairwise_traceback_codon(left_tokens, right_tokens, scores, trace, best_i, best_j, best_score)
end

"""
    _pairwise_align_codon_affine_global(left_tokens, right_tokens, scoring, gap_open, gap_extend)

Run global codon alignment with an affine gap model.
"""
function _pairwise_align_codon_affine_global(left_tokens::AbstractVector{UInt8}, right_tokens::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap_open::Int, gap_extend::Int)
  left_length = length(left_tokens)
  right_length = length(right_tokens)
  negative_infinity = typemin(Int) ÷ 4

  previous_match_scores = fill(negative_infinity, left_length + 1)
  previous_gap_left_scores = fill(negative_infinity, left_length + 1)
  previous_gap_right_scores = fill(negative_infinity, left_length + 1)
  current_match_scores = similar(previous_match_scores)
  current_gap_left_scores = similar(previous_gap_left_scores)
  current_gap_right_scores = similar(previous_gap_right_scores)

  match_trace = zeros(UInt8, left_length + 1, right_length + 1)
  gap_left_trace = zeros(UInt8, left_length + 1, right_length + 1)
  gap_right_trace = zeros(UInt8, left_length + 1, right_length + 1)

  previous_match_scores[1] = 0
  previous_gap_left_scores[1] = _PAIRWISE_AFFINE_NEGATIVE_INFINITY
  previous_gap_right_scores[1] = _PAIRWISE_AFFINE_NEGATIVE_INFINITY

  @inbounds begin
    for i in 2:(left_length+1)
      from_match = _pairwise_affine_transition(previous_match_scores[i-1], gap_open)
      from_gap = _pairwise_affine_transition(previous_gap_left_scores[i-1], gap_extend)
      if from_match >= from_gap
        previous_gap_left_scores[i] = from_match
        gap_left_trace[i, 1] = _PAIRWISE_STATE_MATCH
      else
        previous_gap_left_scores[i] = from_gap
        gap_left_trace[i, 1] = _PAIRWISE_STATE_GAP_LEFT
      end
    end

    for j in 2:(right_length+1)
      current_match_scores[1] = _PAIRWISE_AFFINE_NEGATIVE_INFINITY
      current_gap_left_scores[1] = _PAIRWISE_AFFINE_NEGATIVE_INFINITY

      from_match = _pairwise_affine_transition(previous_match_scores[1], gap_open)
      from_gap = _pairwise_affine_transition(previous_gap_right_scores[1], gap_extend)
      if from_match >= from_gap
        current_gap_right_scores[1] = from_match
        gap_right_trace[1, j] = _PAIRWISE_STATE_MATCH
      else
        current_gap_right_scores[1] = from_gap
        gap_right_trace[1, j] = _PAIRWISE_STATE_GAP_RIGHT
      end

      right_token = right_tokens[j-1]
      for i in 2:(left_length+1)
        left_token = left_tokens[i-1]

        best_prev = previous_match_scores[i-1]
        best_state = _PAIRWISE_STATE_MATCH
        if previous_gap_left_scores[i-1] > best_prev
          best_prev = previous_gap_left_scores[i-1]
          best_state = _PAIRWISE_STATE_GAP_LEFT
        end
        if previous_gap_right_scores[i-1] > best_prev
          best_prev = previous_gap_right_scores[i-1]
          best_state = _PAIRWISE_STATE_GAP_RIGHT
        end
        current_match_scores[i] = best_prev + _pairwise_score(scoring, left_token, right_token)
        match_trace[i, j] = best_state

        from_match = _pairwise_affine_transition(current_match_scores[i-1], gap_open)
        from_gap = _pairwise_affine_transition(current_gap_left_scores[i-1], gap_extend)
        if from_match >= from_gap
          current_gap_left_scores[i] = from_match
          gap_left_trace[i, j] = _PAIRWISE_STATE_MATCH
        else
          current_gap_left_scores[i] = from_gap
          gap_left_trace[i, j] = _PAIRWISE_STATE_GAP_LEFT
        end

        from_match = _pairwise_affine_transition(previous_match_scores[i], gap_open)
        from_gap = _pairwise_affine_transition(previous_gap_right_scores[i], gap_extend)
        if from_match >= from_gap
          current_gap_right_scores[i] = from_match
          gap_right_trace[i, j] = _PAIRWISE_STATE_MATCH
        else
          current_gap_right_scores[i] = from_gap
          gap_right_trace[i, j] = _PAIRWISE_STATE_GAP_RIGHT
        end
      end

      previous_match_scores, current_match_scores = current_match_scores, previous_match_scores
      previous_gap_left_scores, current_gap_left_scores = current_gap_left_scores, previous_gap_left_scores
      previous_gap_right_scores, current_gap_right_scores = current_gap_right_scores, previous_gap_right_scores
    end
  end

  best_score = previous_match_scores[left_length+1]
  best_state = _PAIRWISE_STATE_MATCH
  if previous_gap_left_scores[left_length+1] > best_score
    best_score = previous_gap_left_scores[left_length+1]
    best_state = _PAIRWISE_STATE_GAP_LEFT
  end
  if previous_gap_right_scores[left_length+1] > best_score
    best_score = previous_gap_right_scores[left_length+1]
    best_state = _PAIRWISE_STATE_GAP_RIGHT
  end

  return _pairwise_traceback_codon_affine(
    left_tokens,
    right_tokens,
    match_trace,
    gap_left_trace,
    gap_right_trace,
    left_length + 1,
    right_length + 1,
    best_state,
    best_score,
  )
end

"""
    _pairwise_align_codon_affine_local(left_tokens, right_tokens, scoring, gap_open, gap_extend)

Run local codon alignment with an affine gap model.
"""
function _pairwise_align_codon_affine_local(left_tokens::AbstractVector{UInt8}, right_tokens::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap_open::Int, gap_extend::Int)
  left_length = length(left_tokens)
  right_length = length(right_tokens)

  previous_match_scores = zeros(Int, left_length + 1)
  previous_gap_left_scores = zeros(Int, left_length + 1)
  previous_gap_right_scores = zeros(Int, left_length + 1)
  current_match_scores = zeros(Int, left_length + 1)
  current_gap_left_scores = zeros(Int, left_length + 1)
  current_gap_right_scores = zeros(Int, left_length + 1)

  match_trace = zeros(UInt8, left_length + 1, right_length + 1)
  gap_left_trace = zeros(UInt8, left_length + 1, right_length + 1)
  gap_right_trace = zeros(UInt8, left_length + 1, right_length + 1)

  best_score = 0
  best_i = 1
  best_j = 1
  best_state = _PAIRWISE_TRACE_NONE

  @inbounds for j in 2:(right_length+1)
    right_token = right_tokens[j-1]
    current_match_scores[1] = 0
    current_gap_left_scores[1] = 0
    current_gap_right_scores[1] = 0

    for i in 2:(left_length+1)
      left_token = left_tokens[i-1]

      best_prev = 0
      prev_state = _PAIRWISE_TRACE_NONE

      candidate = previous_match_scores[i-1]
      if candidate > best_prev
        best_prev = candidate
        prev_state = _PAIRWISE_STATE_MATCH
      end
      candidate = previous_gap_left_scores[i-1]
      if candidate > best_prev
        best_prev = candidate
        prev_state = _PAIRWISE_STATE_GAP_LEFT
      end
      candidate = previous_gap_right_scores[i-1]
      if candidate > best_prev
        best_prev = candidate
        prev_state = _PAIRWISE_STATE_GAP_RIGHT
      end

      match_score = best_prev + _pairwise_score(scoring, left_token, right_token)
      if match_score > 0
        current_match_scores[i] = match_score
        match_trace[i, j] = prev_state
        if match_score > best_score
          best_score = match_score
          best_i = i
          best_j = j
          best_state = _PAIRWISE_STATE_MATCH
        end
      else
        current_match_scores[i] = 0
      end

      from_match = _pairwise_affine_transition(current_match_scores[i-1], gap_open)
      from_gap = _pairwise_affine_transition(current_gap_left_scores[i-1], gap_extend)
      gap_score = 0
      gap_state = _PAIRWISE_TRACE_NONE
      if from_match > gap_score
        gap_score = from_match
        gap_state = _PAIRWISE_STATE_MATCH
      end
      if from_gap > gap_score
        gap_score = from_gap
        gap_state = _PAIRWISE_STATE_GAP_LEFT
      end
      if gap_score > 0
        current_gap_left_scores[i] = gap_score
        gap_left_trace[i, j] = gap_state
        if gap_score > best_score
          best_score = gap_score
          best_i = i
          best_j = j
          best_state = _PAIRWISE_STATE_GAP_LEFT
        end
      else
        current_gap_left_scores[i] = 0
      end

      from_match = _pairwise_affine_transition(previous_match_scores[i], gap_open)
      from_gap = _pairwise_affine_transition(previous_gap_right_scores[i], gap_extend)
      gap_score = 0
      gap_state = _PAIRWISE_TRACE_NONE
      if from_match > gap_score
        gap_score = from_match
        gap_state = _PAIRWISE_STATE_MATCH
      end
      if from_gap > gap_score
        gap_score = from_gap
        gap_state = _PAIRWISE_STATE_GAP_RIGHT
      end
      if gap_score > 0
        current_gap_right_scores[i] = gap_score
        gap_right_trace[i, j] = gap_state
        if gap_score > best_score
          best_score = gap_score
          best_i = i
          best_j = j
          best_state = _PAIRWISE_STATE_GAP_RIGHT
        end
      else
        current_gap_right_scores[i] = 0
      end
    end

    previous_match_scores, current_match_scores = current_match_scores, previous_match_scores
    previous_gap_left_scores, current_gap_left_scores = current_gap_left_scores, previous_gap_left_scores
    previous_gap_right_scores, current_gap_right_scores = current_gap_right_scores, previous_gap_right_scores
  end

  return _pairwise_traceback_codon_affine(
    left_tokens,
    right_tokens,
    match_trace,
    gap_left_trace,
    gap_right_trace,
    best_i,
    best_j,
    best_state,
    best_score,
  )
end

"""
    pairwise_align_codons(left, right; kwargs...)

Align two nucleotide strings in codon space, preserving triplet boundaries.

The wrapper uppercases nucleotide input before codon packing for the same
reason as `pairwise_align`.
"""
function pairwise_align_codons(
  left::BioSequence{DNAAlphabet},
  right::BioSequence{DNAAlphabet};
  match::Int=1,
  mismatch::Int=-1,
  substitution_matrix::Union{Nothing,CodonSubstitutionMatrix}=nothing,
  gap::Int=-1,
  gap_open::Union{Nothing,Int}=nothing,
  gap_extend::Union{Nothing,Int}=nothing,
  is_local::Bool=false,
)
  left_tokens = _encode_codon_sequence(left.data)
  right_tokens = _encode_codon_sequence(right.data)
  scoring = substitution_matrix === nothing ?
            LinearPairwiseScoring(match, mismatch) : CodonMatrixPairwiseScoring(substitution_matrix)

  result = if gap_open === nothing && gap_extend === nothing
    is_local ? _pairwise_align_codon_local(left_tokens, right_tokens, scoring, gap) :
    _pairwise_align_codon_global(left_tokens, right_tokens, scoring, gap)
  else
    affine_gap_open, affine_gap_extend = _pairwise_gap_parameters(gap, gap_open, gap_extend)
    is_local ? _pairwise_align_codon_affine_local(left_tokens, right_tokens, scoring, affine_gap_open, affine_gap_extend) :
    _pairwise_align_codon_affine_global(left_tokens, right_tokens, scoring, affine_gap_open, affine_gap_extend)
  end

  _ctx = active_provenance_context()
  _ctx !== nothing && register_provenance!(_ctx, "pairwise_align_codons";
    parents=provenance_parent_ids(left, right),
    parameters=(match=match, mismatch=mismatch, gap=gap,
      affine=gap_open !== nothing || gap_extend !== nothing,
      is_local=is_local))
  return result
end

# Removed pairwise_align_codons(::AbstractString, ::AbstractString) - use DNASeq instead
function pairwise_align_codons(left::AbstractString, right::AbstractString; kwargs...)
  return pairwise_align_codons(DNASeq(left; validate=false), DNASeq(right; validate=false); kwargs...)
end

"""
    needleman_wunsch_codons(left, right; kwargs...)

Wrapper for global codon alignment.
"""
needleman_wunsch_codons(left, right; kwargs...) = pairwise_align_codons(left, right; is_local=false, kwargs...)

"""
    smith_waterman_codons(left, right; kwargs...)

Wrapper for local codon alignment.
"""
smith_waterman_codons(left, right; kwargs...) = pairwise_align_codons(left, right; is_local=true, kwargs...)

"""
    local_align_codons(left, right; kwargs...)

Backward-compatible alias for `smith_waterman_codons`.
"""
local_align_codons(left, right; kwargs...) = smith_waterman_codons(left, right; kwargs...)

"""
    needleman_wunsch(left, right; kwargs...)

Explicit wrapper for the global pairwise alignment path.
"""
needleman_wunsch(left, right; kwargs...) = pairwise_align(left, right; is_local=false, kwargs...)

"""
    smith_waterman(left, right; kwargs...)

Explicit wrapper for the local pairwise alignment path.
"""
smith_waterman(left, right; kwargs...) = pairwise_align(left, right; is_local=true, kwargs...)

"""
    local_align(left, right; kwargs...)

Backward-compatible alias for `smith_waterman`.
"""
local_align(left, right; kwargs...) = smith_waterman(left, right; kwargs...)

"""
    _pairwise_align_global(::Type{A}, left_bytes, right_bytes, scoring, gap)

Run global pairwise alignment with a linear gap model.
"""
function _pairwise_align_global(::Type{A}, left_bytes::AbstractVector{UInt8}, right_bytes::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap::Int) where {A<:BioAlphabet}
  left_length = length(left_bytes)
  right_length = length(right_bytes)

  previous_scores = Vector{Int}(undef, left_length + 1)
  current_scores = similar(previous_scores)
  trace = Matrix{UInt8}(undef, left_length + 1, right_length + 1)

  @inbounds begin
    previous_scores[1] = 0
    trace[1, 1] = 0x00

    for i in 2:(left_length+1)
      previous_scores[i] = previous_scores[i-1] + gap
      trace[i, 1] = 0x02
    end

    for j in 2:(right_length+1)
      current_scores[1] = previous_scores[1] + gap
      trace[1, j] = 0x03
      right_byte = right_bytes[j-1]
      for i in 2:(left_length+1)
        left_byte = left_bytes[i-1]

        diag_score = previous_scores[i-1] + _pairwise_score(scoring, left_byte, right_byte)
        # current_scores[i-1] == S[row=i-2, col=j-1]: same column, row above -> VERTICAL move (consumes left_byte only)
        vertical_score = current_scores[i-1] + gap
        # previous_scores[i]   == S[row=i-1, col=j-2]: same row, column to the left -> HORIZONTAL move (consumes right_byte only)
        horizontal_score = previous_scores[i] + gap

        cell_score = diag_score
        cell_trace = 0x01

        if vertical_score > cell_score
          cell_score = vertical_score
          cell_trace = 0x02
        end
        if horizontal_score > cell_score
          cell_score = horizontal_score
          cell_trace = 0x03
        end

        current_scores[i] = cell_score
        trace[i, j] = cell_trace
      end

      previous_scores, current_scores = current_scores, previous_scores
    end
  end

  return _pairwise_traceback(A, left_bytes, right_bytes, nothing, trace, left_length + 1, right_length + 1, previous_scores[left_length+1])
end


"""
    _pairwise_align_local(::Type{A}, left_bytes, right_bytes, scoring, gap)

Run local pairwise alignment with a linear gap model.
"""
function _pairwise_align_local(::Type{A}, left_bytes::AbstractVector{UInt8}, right_bytes::AbstractVector{UInt8}, scoring::AbstractPairwiseScoring, gap::Int) where {A<:BioAlphabet}
  left_length = length(left_bytes)
  right_length = length(right_bytes)

  scores = Matrix{Int}(undef, left_length + 1, right_length + 1)
  trace = Matrix{UInt8}(undef, left_length + 1, right_length + 1)

  best_score = 0
  best_i = 1
  best_j = 1

  @inbounds begin
    scores[1, 1] = 0
    trace[1, 1] = 0x00

    for i in 2:(left_length+1)
      scores[i, 1] = 0
      trace[i, 1] = 0x00
    end

    for j in 2:(right_length+1)
      scores[1, j] = 0
      trace[1, j] = 0x00
    end

    for j in 2:(right_length+1)
      right_byte = right_bytes[j-1]
      for i in 2:(left_length+1)
        left_byte = left_bytes[i-1]

        diag_score = scores[i-1, j-1] + _pairwise_score(scoring, left_byte, right_byte)
        up_score = scores[i-1, j] + gap
        left_score = scores[i, j-1] + gap

        cell_score = 0
        cell_trace = 0x00

        if diag_score > cell_score
          cell_score = diag_score
          cell_trace = 0x01
        end
        if up_score > cell_score
          cell_score = up_score
          cell_trace = 0x02
        end
        if left_score > cell_score
          cell_score = left_score
          cell_trace = 0x03
        end

        scores[i, j] = cell_score
        trace[i, j] = cell_trace

        if cell_score > best_score
          best_score = cell_score
          best_i = i
          best_j = j
        end
      end
    end
  end

  return _pairwise_traceback(A, left_bytes, right_bytes, scores, trace, best_i, best_j, best_score)
end

"""
    _pairwise_align_affine_global(::Type{A}, left_bytes, right_bytes, scoring, gap_open, gap_extend)

Run global pairwise alignment with an affine gap model.
"""
function _pairwise_align_affine_global(
  ::Type{A},
  left_bytes::AbstractVector{UInt8},
  right_bytes::AbstractVector{UInt8},
  scoring::AbstractPairwiseScoring,
  gap_open::Int,
  gap_extend::Int,
) where {A<:BioAlphabet}
  left_length = length(left_bytes)
  right_length = length(right_bytes)

  previous_match_scores = fill(_PAIRWISE_AFFINE_NEGATIVE_INFINITY, left_length + 1)
  previous_gap_left_scores = fill(_PAIRWISE_AFFINE_NEGATIVE_INFINITY, left_length + 1)
  previous_gap_right_scores = fill(_PAIRWISE_AFFINE_NEGATIVE_INFINITY, left_length + 1)
  current_match_scores = similar(previous_match_scores)
  current_gap_left_scores = similar(previous_gap_left_scores)
  current_gap_right_scores = similar(previous_gap_right_scores)

  match_trace = zeros(UInt8, left_length + 1, right_length + 1)
  gap_left_trace = zeros(UInt8, left_length + 1, right_length + 1)
  gap_right_trace = zeros(UInt8, left_length + 1, right_length + 1)

  previous_match_scores[1] = 0
  previous_gap_left_scores[1] = _PAIRWISE_AFFINE_NEGATIVE_INFINITY
  previous_gap_right_scores[1] = _PAIRWISE_AFFINE_NEGATIVE_INFINITY

  @inbounds begin
    for i in 2:(left_length+1)
      from_match = _pairwise_affine_transition(previous_match_scores[i-1], gap_open)
      from_gap = _pairwise_affine_transition(previous_gap_left_scores[i-1], gap_extend)
      if from_match >= from_gap
        previous_gap_left_scores[i] = from_match
        gap_left_trace[i, 1] = _PAIRWISE_STATE_MATCH
      else
        previous_gap_left_scores[i] = from_gap
        gap_left_trace[i, 1] = _PAIRWISE_STATE_GAP_LEFT
      end
    end

    for j in 2:(right_length+1)
      current_match_scores[1] = _PAIRWISE_AFFINE_NEGATIVE_INFINITY
      current_gap_left_scores[1] = _PAIRWISE_AFFINE_NEGATIVE_INFINITY

      from_match = _pairwise_affine_transition(previous_match_scores[1], gap_open)
      from_gap = _pairwise_affine_transition(previous_gap_right_scores[1], gap_extend)
      if from_match >= from_gap
        current_gap_right_scores[1] = from_match
        gap_right_trace[1, j] = _PAIRWISE_STATE_MATCH
      else
        current_gap_right_scores[1] = from_gap
        gap_right_trace[1, j] = _PAIRWISE_STATE_GAP_RIGHT
      end

      right_byte = right_bytes[j-1]
      for i in 2:(left_length+1)
        left_byte = left_bytes[i-1]

        best_prev = previous_match_scores[i-1]
        best_state = _PAIRWISE_STATE_MATCH
        if previous_gap_left_scores[i-1] > best_prev
          best_prev = previous_gap_left_scores[i-1]
          best_state = _PAIRWISE_STATE_GAP_LEFT
        end
        if previous_gap_right_scores[i-1] > best_prev
          best_prev = previous_gap_right_scores[i-1]
          best_state = _PAIRWISE_STATE_GAP_RIGHT
        end
        current_match_scores[i] = best_prev + _pairwise_score(scoring, left_byte, right_byte)
        match_trace[i, j] = best_state

        from_match = _pairwise_affine_transition(current_match_scores[i-1], gap_open)
        from_gap = _pairwise_affine_transition(current_gap_left_scores[i-1], gap_extend)
        if from_match >= from_gap
          current_gap_left_scores[i] = from_match
          gap_left_trace[i, j] = _PAIRWISE_STATE_MATCH
        else
          current_gap_left_scores[i] = from_gap
          gap_left_trace[i, j] = _PAIRWISE_STATE_GAP_LEFT
        end

        from_match = _pairwise_affine_transition(previous_match_scores[i], gap_open)
        from_gap = _pairwise_affine_transition(previous_gap_right_scores[i], gap_extend)
        if from_match >= from_gap
          current_gap_right_scores[i] = from_match
          gap_right_trace[i, j] = _PAIRWISE_STATE_MATCH
        else
          current_gap_right_scores[i] = from_gap
          gap_right_trace[i, j] = _PAIRWISE_STATE_GAP_RIGHT
        end
      end

      previous_match_scores, current_match_scores = current_match_scores, previous_match_scores
      previous_gap_left_scores, current_gap_left_scores = current_gap_left_scores, previous_gap_left_scores
      previous_gap_right_scores, current_gap_right_scores = current_gap_right_scores, previous_gap_right_scores
    end
  end

  best_score = previous_match_scores[left_length+1]
  best_state = _PAIRWISE_STATE_MATCH
  if previous_gap_left_scores[left_length+1] > best_score
    best_score = previous_gap_left_scores[left_length+1]
    best_state = _PAIRWISE_STATE_GAP_LEFT
  end
  if previous_gap_right_scores[left_length+1] > best_score
    best_score = previous_gap_right_scores[left_length+1]
    best_state = _PAIRWISE_STATE_GAP_RIGHT
  end

  return _pairwise_traceback_affine(
    A,
    left_bytes,
    right_bytes,
    match_trace,
    gap_left_trace,
    gap_right_trace,
    left_length + 1,
    right_length + 1,
    best_state,
    best_score,
  )
end

"""
    _pairwise_align_affine_local(::Type{A}, left_bytes, right_bytes, scoring, gap_open, gap_extend)

Run local pairwise alignment with an affine gap model.
"""
function _pairwise_align_affine_local(
  ::Type{A},
  left_bytes::AbstractVector{UInt8},
  right_bytes::AbstractVector{UInt8},
  scoring::AbstractPairwiseScoring,
  gap_open::Int,
  gap_extend::Int,
) where {A<:BioAlphabet}
  left_length = length(left_bytes)
  right_length = length(right_bytes)

  # O(N) memory buffers instead of O(N * M)
  previous_match_scores = zeros(Int, left_length + 1)
  previous_gap_left_scores = zeros(Int, left_length + 1)
  previous_gap_right_scores = zeros(Int, left_length + 1)

  current_match_scores = zeros(Int, left_length + 1)
  current_gap_left_scores = zeros(Int, left_length + 1)
  current_gap_right_scores = zeros(Int, left_length + 1)

  # Trace matrices still require O(N * M) but use only 1 byte per cell
  match_trace = zeros(UInt8, left_length + 1, right_length + 1)
  gap_left_trace = zeros(UInt8, left_length + 1, right_length + 1)
  gap_right_trace = zeros(UInt8, left_length + 1, right_length + 1)

  best_score = 0
  best_i = 1
  best_j = 1
  best_state = _PAIRWISE_TRACE_NONE

  @inbounds for j in 2:(right_length+1)
    right_byte = right_bytes[j-1]

    # Reset the first cell of the current column to 0
    current_match_scores[1] = 0
    current_gap_left_scores[1] = 0
    current_gap_right_scores[1] = 0

    for i in 2:(left_length+1)
      left_byte = left_bytes[i-1]

      best_prev = 0
      prev_state = _PAIRWISE_TRACE_NONE

      # Diagonal: previous_scores[i - 1] means scores[i-1, j-1]
      candidate = previous_match_scores[i-1]
      if candidate > best_prev
        best_prev = candidate
        prev_state = _PAIRWISE_STATE_MATCH
      end
      candidate = previous_gap_left_scores[i-1]
      if candidate > best_prev
        best_prev = candidate
        prev_state = _PAIRWISE_STATE_GAP_LEFT
      end
      candidate = previous_gap_right_scores[i-1]
      if candidate > best_prev
        best_prev = candidate
        prev_state = _PAIRWISE_STATE_GAP_RIGHT
      end

      match_score = best_prev + _pairwise_score(scoring, left_byte, right_byte)
      if match_score > 0
        current_match_scores[i] = match_score
        match_trace[i, j] = prev_state
        if match_score > best_score
          best_score = match_score
          best_i = i
          best_j = j
          best_state = _PAIRWISE_STATE_MATCH
        end
      else
        current_match_scores[i] = 0
      end

      # Up: current_scores[i-1] means scores[i-1, j]
      from_match = _pairwise_affine_transition(current_match_scores[i-1], gap_open)
      from_gap = _pairwise_affine_transition(current_gap_left_scores[i-1], gap_extend)
      gap_score = 0
      gap_state = _PAIRWISE_TRACE_NONE
      if from_match > gap_score
        gap_score = from_match
        gap_state = _PAIRWISE_STATE_MATCH
      end
      if from_gap > gap_score
        gap_score = from_gap
        gap_state = _PAIRWISE_STATE_GAP_LEFT
      end
      if gap_score > 0
        current_gap_left_scores[i] = gap_score
        gap_left_trace[i, j] = gap_state
        if gap_score > best_score
          best_score = gap_score
          best_i = i
          best_j = j
          best_state = _PAIRWISE_STATE_GAP_LEFT
        end
      else
        current_gap_left_scores[i] = 0
      end

      # Left: previous_scores[i] means scores[i, j-1]
      from_match = previous_match_scores[i] + gap_open
      from_gap = previous_gap_right_scores[i] + gap_extend
      gap_score = 0
      gap_state = _PAIRWISE_TRACE_NONE
      if from_match > gap_score
        gap_score = from_match
        gap_state = _PAIRWISE_STATE_MATCH
      end
      if from_gap > gap_score
        gap_score = from_gap
        gap_state = _PAIRWISE_STATE_GAP_RIGHT
      end
      if gap_score > 0
        current_gap_right_scores[i] = gap_score
        gap_right_trace[i, j] = gap_state
        if gap_score > best_score
          best_score = gap_score
          best_i = i
          best_j = j
          best_state = _PAIRWISE_STATE_GAP_RIGHT
        end
      else
        current_gap_right_scores[i] = 0
      end
    end

    # Swap current memory to previous
    previous_match_scores, current_match_scores = current_match_scores, previous_match_scores
    previous_gap_left_scores, current_gap_left_scores = current_gap_left_scores, previous_gap_left_scores
    previous_gap_right_scores, current_gap_right_scores = current_gap_right_scores, previous_gap_right_scores
  end

  return _pairwise_traceback_affine(
    A,
    left_bytes,
    right_bytes,
    match_trace,
    gap_left_trace,
    gap_right_trace,
    best_i,
    best_j,
    best_state,
    best_score,
  )
end

"""
    _pairwise_traceback(::Type{A}, left_bytes, right_bytes, scores, trace, i, j, best_score)

Reconstruct a nucleotide alignment from standard traceback data.
"""
function _pairwise_traceback(
  ::Type{A},
  left_bytes::AbstractVector{UInt8},
  right_bytes::AbstractVector{UInt8},
  scores::Union{Nothing,Matrix{Int}},
  trace::Matrix{UInt8},
  i::Int,
  j::Int,
  best_score::Int,
) where {A<:BioAlphabet}
  aligned_left = Vector{UInt8}(undef, length(left_bytes) + length(right_bytes))
  aligned_right = Vector{UInt8}(undef, length(left_bytes) + length(right_bytes))
  write_index = length(aligned_left)
  aligned_length = 0
  matches = 0

  current_i = i
  current_j = j

  while current_i > 1 || current_j > 1
    if scores !== nothing && scores[current_i, current_j] == 0 && trace[current_i, current_j] == 0x00
      break
    end

    direction = trace[current_i, current_j]
    direction == 0x00 && break

    if direction == 0x01
      left_byte = left_bytes[current_i-1]
      right_byte = right_bytes[current_j-1]
      aligned_left[write_index] = left_byte
      aligned_right[write_index] = right_byte
      matches += left_byte == right_byte ? 1 : 0
      current_i -= 1
      current_j -= 1
    elseif direction == 0x02
      aligned_left[write_index] = left_bytes[current_i-1]
      aligned_right[write_index] = UInt8('-')
      current_i -= 1
    else
      aligned_left[write_index] = UInt8('-')
      aligned_right[write_index] = right_bytes[current_j-1]
      current_j -= 1
    end

    aligned_length += 1
    write_index -= 1
  end

  start_index = write_index + 1
  final_left = BioSequence{A}(aligned_left[start_index:end]; validate=false)
  final_right = BioSequence{A}(aligned_right[start_index:end]; validate=false)
  identity = aligned_length == 0 ? 0.0 : matches / aligned_length

  return PairwiseAlignmentResult(final_left, final_right, best_score, matches, identity)
end

"""
    _pairwise_traceback_affine(::Type{A}, left_bytes, right_bytes, match_trace, gap_left_trace, gap_right_trace, i, j, start_state, best_score)

Reconstruct a nucleotide alignment from affine-gap traceback data.
"""
function _pairwise_traceback_affine(
  ::Type{A},
  left_bytes::AbstractVector{UInt8},
  right_bytes::AbstractVector{UInt8},
  match_trace::Matrix{UInt8},
  gap_left_trace::Matrix{UInt8},
  gap_right_trace::Matrix{UInt8},
  i::Int,
  j::Int,
  start_state::UInt8,
  best_score::Int,
) where {A<:BioAlphabet}
  aligned_left = Vector{UInt8}(undef, length(left_bytes) + length(right_bytes))
  aligned_right = Vector{UInt8}(undef, length(left_bytes) + length(right_bytes))
  write_index = length(aligned_left)
  aligned_length = 0
  matches = 0

  current_i = i
  current_j = j
  current_state = start_state

  while current_state != _PAIRWISE_TRACE_NONE && (current_i > 1 || current_j > 1)
    if current_state == _PAIRWISE_STATE_MATCH
      left_byte = left_bytes[current_i-1]
      right_byte = right_bytes[current_j-1]
      aligned_left[write_index] = left_byte
      aligned_right[write_index] = right_byte
      matches += left_byte == right_byte ? 1 : 0
      current_state = match_trace[current_i, current_j]
      current_i -= 1
      current_j -= 1
    elseif current_state == _PAIRWISE_STATE_GAP_LEFT
      aligned_left[write_index] = left_bytes[current_i-1]
      aligned_right[write_index] = UInt8('-')
      current_state = gap_left_trace[current_i, current_j]
      current_i -= 1
    else
      aligned_left[write_index] = UInt8('-')
      aligned_right[write_index] = right_bytes[current_j-1]
      current_state = gap_right_trace[current_i, current_j]
      current_j -= 1
    end

    aligned_length += 1
    write_index -= 1
  end

  start_index = write_index + 1
  final_left = BioSequence{A}(aligned_left[start_index:end]; validate=false)
  final_right = BioSequence{A}(aligned_right[start_index:end]; validate=false)
  identity = aligned_length == 0 ? 0.0 : matches / aligned_length

  return PairwiseAlignmentResult(final_left, final_right, best_score, matches, identity)
end

# ==============================================================================
# Interactive HTML5/Canvas Visualizers for Alignment Module
# ==============================================================================

"""
    _escape_html(s)

Escape a string for safe interpolation into HTML markup.
"""
function _escape_html(s::AbstractString)
    s = replace(s, '&' => "&amp;")
    s = replace(s, '<' => "&lt;")
    s = replace(s, '>' => "&gt;")
    s = replace(s, '"' => "&quot;")
    s = replace(s, '\'' => "&#39;")
    return s
end

"""
    _json_escape(s)

Escape a string for inclusion inside a JSON literal embedded in a `<script>` tag.
"""
@inline function _json_escape(s::AbstractString)
    s = escape_string(s)
    s = replace(s, "</" => "<\\/")
    return s
end

function _alignment_to_json(res::PairwiseAlignmentResult)
    left_str = _json_escape(String(res.left))
    right_str = _json_escape(String(res.right))
    return "{\"left\":\"$left_str\",\"right\":\"$right_str\",\"score\":$(res.score),\"matches\":$(res.matches),\"identity\":$(res.identity)}"
end

"""
    visualize_alignment_html(res::PairwiseAlignmentResult; title="Pairwise Sequence Alignment") -> String

Generate an interactive HTML5/Canvas visualization for pairwise sequence alignment results.
Displays aligned sequences with residue color-coding, consensus match marks, hover tooltips, and linear sequence viewport controls.
"""
function visualize_alignment_html(res::PairwiseAlignmentResult; title::String="Pairwise Sequence Alignment")
    aln_len = length(res.left)
    ident_pct = round(res.identity * 100, digits=1)
    ident_color = ident_pct >= 80.0 ? "#10b981" : ident_pct >= 50.0 ? "#f59e0b" : "#ef4444"
    aln_json = _alignment_to_json(res)

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(_escape_html(title))</title>
    <style>
        :root {
            --bg: #0f172a; --panel: #1e293b; --border: #334155;
            --text: #f8fafc; --muted: #94a3b8; --accent: #38bdf8;
            --a: #10b981; --c: #3b82f6; --g: #f59e0b; --t: #ef4444; --gap: #64748b;
        }
        body { font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 24px; }
        .card { background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 20px; margin-bottom: 20px; box-shadow: 0 4px 12px rgba(0,0,0,0.3); }
        .title { font-size: 1.5rem; font-weight: 700; color: var(--accent); margin: 0 0 12px 0; }
        .meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 14px; margin-bottom: 20px; }
        .meta-item { background: rgba(15,23,42,0.6); padding: 10px 14px; border-radius: 8px; border: 1px solid var(--border); }
        .meta-label { font-size: 0.75rem; text-transform: uppercase; color: var(--muted); letter-spacing: 0.5px; }
        .meta-val { font-size: 1.1rem; font-weight: 600; color: var(--text); margin-top: 2px; }
        canvas { width: 100%; height: 280px; display: block; border-radius: 8px; background: #0b1329; }
        .controls { display: flex; gap: 14px; align-items: center; margin-bottom: 14px; }
        .btn { background: #3b82f6; color: #fff; border: none; padding: 8px 16px; border-radius: 6px; font-weight: 600; cursor: pointer; transition: 0.2s; }
        .btn:hover { background: #2563eb; }
        .legend { display: flex; gap: 16px; font-weight: 600; font-size: 0.85rem; margin-top: 12px; }
        .leg-item { display: flex; align-items: center; gap: 6px; }
        .dot { width: 10px; height: 10px; border-radius: 50%; display: inline-block; }
    </style>
</head>
<body>
    <div class="card">
        <div class="title">⚔️ $(_escape_html(title))</div>
        <div class="meta-grid">
            <div class="meta-item"><div class="meta-label">Alignment Score</div><div class="meta-val">$(res.score)</div></div>
            <div class="meta-item"><div class="meta-label">Identity</div><div class="meta-val" style="color:$(ident_color)">$(ident_pct)%</div></div>
            <div class="meta-item"><div class="meta-label">Matches</div><div class="meta-val">$(res.matches) / $(aln_len)</div></div>
            <div class="meta-item"><div class="meta-label">Aligned Length</div><div class="meta-val">$(aln_len) bp</div></div>
        </div>

        <div class="controls">
            <label style="font-size:0.88rem; color:var(--muted)">Scroll Position:</label>
            <input type="range" id="posSlider" min="0" max="100" value="0" style="flex:1;" oninput="drawAlignment()">
            <label style="font-size:0.88rem; color:var(--muted)">Window:</label>
            <select id="winSelect" onchange="drawAlignment()" style="background:#0f172a; color:#fff; border:1px solid var(--border); padding:6px 10px; border-radius:6px;">
                <option value="40">40 bp</option>
                <option value="80" selected>80 bp</option>
                <option value="120">120 bp</option>
            </select>
        </div>

        <canvas id="alnCanvas"></canvas>

        <div class="legend">
            <div class="leg-item"><span class="dot" style="background:var(--a)"></span> A</div>
            <div class="leg-item"><span class="dot" style="background:var(--c)"></span> C</div>
            <div class="leg-item"><span class="dot" style="background:var(--g)"></span> G</div>
            <div class="leg-item"><span class="dot" style="background:var(--t)"></span> T / U</div>
            <div class="leg-item"><span class="dot" style="background:var(--gap)"></span> Gap (-)</div>
        </div>
    </div>

    <script>
        const aln = $(aln_json);
        const canvas = document.getElementById('alnCanvas');
        const ctx = canvas.getContext('2d');

        function getColor(b) {
            const u = b.toUpperCase();
            if (u === 'A') return '#10b981';
            if (u === 'C') return '#3b82f6';
            if (u === 'G') return '#f59e0b';
            if (u === 'T' || u === 'U') return '#ef4444';
            if (u === '-') return '#64748b';
            return '#38bdf8';
        }

        function drawAlignment() {
            canvas.width = canvas.parentElement.clientWidth * window.devicePixelRatio;
            canvas.height = 280 * window.devicePixelRatio;
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

            const w = canvas.parentElement.clientWidth;
            const h = 280;
            ctx.clearRect(0, 0, w, h);

            const alnLen = aln.left.length;
            if (alnLen === 0) return;

            const winSize = parseInt(document.getElementById('winSelect').value);
            const sliderVal = parseInt(document.getElementById('posSlider').value);
            const startIdx = Math.floor((sliderVal / 100) * Math.max(0, alnLen - winSize));
            const endIdx = Math.min(startIdx + winSize, alnLen);

            const step = (w - 100) / (endIdx - startIdx);

            ctx.font = 'bold 12px monospace';
            ctx.textAlign = 'center';

            // Labels
            ctx.fillStyle = '#94a3b8';
            ctx.textAlign = 'left';
            ctx.fillText('Seq 1:', 20, 65);
            ctx.fillText('Match:', 20, 115);
            ctx.fillText('Seq 2:', 20, 165);

            for (let i = startIdx; i < endIdx; i++) {
                const x = 90 + (i - startIdx) * step + step / 2;
                const b1 = aln.left[i];
                const b2 = aln.right[i];

                // Seq 1 tile
                ctx.fillStyle = getColor(b1);
                ctx.fillRect(x - step/2 + 2, 45, step - 4, 30);
                ctx.fillStyle = '#ffffff';
                ctx.textAlign = 'center';
                ctx.fillText(b1, x, 65);

                // Match symbol
                if (b1 === b2 && b1 !== '-') {
                    ctx.fillStyle = '#10b981';
                    ctx.fillText('|', x, 115);
                } else if (b1 !== '-' && b2 !== '-') {
                    ctx.fillStyle = '#ef4444';
                    ctx.fillText('•', x, 115);
                } else {
                    ctx.fillStyle = '#64748b';
                    ctx.fillText(' ', x, 115);
                }

                // Seq 2 tile
                ctx.fillStyle = getColor(b2);
                ctx.fillRect(x - step/2 + 2, 145, step - 4, 30);
                ctx.fillStyle = '#ffffff';
                ctx.fillText(b2, x, 165);

                // Coordinate tick
                if ((i + 1) % 10 === 0 || i === startIdx) {
                    ctx.fillStyle = '#64748b';
                    ctx.font = '10px system-ui';
                    ctx.fillText((i + 1) + '', x, 210);
                    ctx.font = 'bold 12px monospace';
                }
            }
        }

        window.addEventListener('resize', drawAlignment);
        setTimeout(drawAlignment, 50);
    </script>
</body>
</html>
"""
end

function to_html(res::PairwiseAlignmentResult)
    return visualize_alignment_html(res)
end

function _posterior_matrix_to_json(res::PosteriorAlignmentResult)
    mat = res.posterior_matrix[:, :, 1]
    M, N = size(mat)
    step_m = max(1, cld(M, 150))
    step_n = max(1, cld(N, 150))
    grid = [round(Float64(mat[i, j]), digits=3) for i in 1:step_m:M, j in 1:step_n:N]
    rows = [ "[" * join(grid[i, :], ",") * "]" for i in 1:size(grid, 1) ]
    matrix_json = "[" * join(rows, ",") * "]"
    aln_json = _alignment_to_json(res.consensus_alignment)
    return "{\"matrix\":$matrix_json,\"expected_accuracy\":$(res.expected_accuracy),\"log_likelihood\":$(res.log_likelihood),\"consensus\":$aln_json}"
end

"""
    visualize_posterior_alignment_html(res::PosteriorAlignmentResult; title="Posterior Alignment Matrix") -> String

Generate an interactive HTML5/Canvas 2D posterior match probability heatmap matrix.
"""
function visualize_posterior_alignment_html(res::PosteriorAlignmentResult; title::String="Posterior Alignment Matrix")
    exp_acc = round(res.expected_accuracy, digits=3)
    log_lh = round(res.log_likelihood, digits=2)
    post_json = _posterior_matrix_to_json(res)

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(_escape_html(title))</title>
    <style>
        :root {
            --bg: #0f172a; --panel: #1e293b; --border: #334155;
            --text: #f8fafc; --muted: #94a3b8; --accent: #38bdf8;
        }
        body { font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 24px; }
        .card { background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 20px; margin-bottom: 20px; box-shadow: 0 4px 12px rgba(0,0,0,0.3); }
        .title { font-size: 1.5rem; font-weight: 700; color: var(--accent); margin: 0 0 12px 0; }
        .meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 14px; margin-bottom: 20px; }
        .meta-item { background: rgba(15,23,42,0.6); padding: 10px 14px; border-radius: 8px; border: 1px solid var(--border); }
        .meta-label { font-size: 0.75rem; text-transform: uppercase; color: var(--muted); }
        .meta-val { font-size: 1.1rem; font-weight: 600; color: var(--text); margin-top: 2px; }
        canvas { width: 100%; height: 380px; display: block; border-radius: 8px; background: #0b1329; }
    </style>
</head>
<body>
    <div class="card">
        <div class="title">🔥 $(_escape_html(title))</div>
        <div class="meta-grid">
            <div class="meta-item"><div class="meta-label">Expected Accuracy</div><div class="meta-val">$(exp_acc)</div></div>
            <div class="meta-item"><div class="meta-label">Log Likelihood</div><div class="meta-val">$(log_lh)</div></div>
            <div class="meta-item"><div class="meta-label">Consensus Score</div><div class="meta-val">$(res.consensus_alignment.score)</div></div>
            <div class="meta-item"><div class="meta-label">Identity</div><div class="meta-val">$(round(res.consensus_alignment.identity * 100, digits=1))%</div></div>
        </div>

        <canvas id="matrixCanvas"></canvas>
    </div>

    <script>
        const data = $(post_json);
        const canvas = document.getElementById('matrixCanvas');
        const ctx = canvas.getContext('2d');

        function drawMatrix() {
            canvas.width = canvas.parentElement.clientWidth * window.devicePixelRatio;
            canvas.height = 380 * window.devicePixelRatio;
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

            const w = canvas.parentElement.clientWidth;
            const h = 380;
            ctx.clearRect(0, 0, w, h);

            const grid = data.matrix;
            const rows = grid.length;
            if (rows === 0) return;
            const cols = grid[0].length;

            const cellW = (w - 80) / cols;
            const cellH = (h - 60) / rows;

            for (let r = 0; r < rows; r++) {
                for (let c = 0; c < cols; c++) {
                    const p = grid[r][c];
                    const red = Math.round(p * 245);
                    const green = Math.round(p * 158 + (1 - p) * 19);
                    const blue = Math.round((1 - p) * 41 + p * 248);

                    ctx.fillStyle = `rgb(\${red},\${green},\${blue})`;
                    ctx.fillRect(50 + c * cellW, 30 + r * cellH, cellW, cellH);
                }
            }

            ctx.fillStyle = '#94a3b8';
            ctx.font = '11px system-ui';
            ctx.fillText('Target Sequence Position →', w / 2 - 60, h - 8);
        }

        window.addEventListener('resize', drawMatrix);
        setTimeout(drawMatrix, 50);
    </script>
</body>
</html>
"""
end

function to_html(res::PosteriorAlignmentResult)
    return visualize_posterior_alignment_html(res)
end

function _sequence_graph_to_json(graph::SequenceGraph)
    node_items = [ "{\"id\":$i,\"seq\":\"$(_json_escape(String(n)))\"}" for (i, n) in enumerate(graph.nodes) ]
    edge_items = [ "{\"u\":$(e[1]),\"v\":$(e[2])}" for e in graph.edges ]
    return "{\"nodes\":[" * join(node_items, ",") * "],\"edges\":[" * join(edge_items, ",") * "]}"
end

"""
    visualize_graph_alignment_html(res::GraphAlignmentResult; title="Sequence Graph Alignment") -> String

Generate an interactive HTML5/Canvas DAG sequence graph alignment viewer.
"""
function visualize_graph_alignment_html(res::GraphAlignmentResult; title::String="Sequence Graph Alignment")
    node_count = length(res.node_scores)
    path_len = length(res.graph_path)
    score = res.query_alignment.score
    path_json = "[" * join(res.graph_path, ",") * "]"

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(_escape_html(title))</title>
    <style>
        :root {
            --bg: #0f172a; --panel: #1e293b; --border: #334155;
            --text: #f8fafc; --muted: #94a3b8; --accent: #38bdf8;
        }
        body { font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 24px; }
        .card { background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 20px; margin-bottom: 20px; box-shadow: 0 4px 12px rgba(0,0,0,0.3); }
        .title { font-size: 1.5rem; font-weight: 700; color: var(--accent); margin: 0 0 12px 0; }
        .meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 14px; margin-bottom: 20px; }
        .meta-item { background: rgba(15,23,42,0.6); padding: 10px 14px; border-radius: 8px; border: 1px solid var(--border); }
        .meta-label { font-size: 0.75rem; text-transform: uppercase; color: var(--muted); }
        .meta-val { font-size: 1.1rem; font-weight: 600; color: var(--text); margin-top: 2px; }
        canvas { width: 100%; height: 340px; display: block; border-radius: 8px; background: #0b1329; }
    </style>
</head>
<body>
    <div class="card">
        <div class="title">🕸️ $(_escape_html(title))</div>
        <div class="meta-grid">
            <div class="meta-item"><div class="meta-label">Path Score</div><div class="meta-val">$(score)</div></div>
            <div class="meta-item"><div class="meta-label">Path Nodes</div><div class="meta-val">$(path_len)</div></div>
            <div class="meta-item"><div class="meta-label">Graph Nodes</div><div class="meta-val">$(node_count)</div></div>
            <div class="meta-item"><div class="meta-label">Identity</div><div class="meta-val">$(round(res.query_alignment.identity * 100, digits=1))%</div></div>
        </div>

        <canvas id="graphCanvas"></canvas>
    </div>

    <script>
        const path = $(path_json);
        const canvas = document.getElementById('graphCanvas');
        const ctx = canvas.getContext('2d');

        function drawGraph() {
            canvas.width = canvas.parentElement.clientWidth * window.devicePixelRatio;
            canvas.height = 340 * window.devicePixelRatio;
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

            const w = canvas.parentElement.clientWidth;
            const h = 340;
            ctx.clearRect(0, 0, w, h);

            if (path.length === 0) return;

            const step = (w - 100) / path.length;
            ctx.strokeStyle = '#10b981';
            ctx.lineWidth = 3;
            ctx.beginPath();

            path.forEach((nodeId, idx) => {
                const x = 50 + idx * step + step/2;
                const y = h / 2;
                if (idx === 0) ctx.moveTo(x, y);
                else ctx.lineTo(x, y);
            });
            ctx.stroke();

            path.forEach((nodeId, idx) => {
                const x = 50 + idx * step + step/2;
                const y = h / 2;
                ctx.fillStyle = '#38bdf8';
                ctx.beginPath();
                ctx.arc(x, y, 12, 0, 2 * Math.PI);
                ctx.fill();

                ctx.fillStyle = '#0f172a';
                ctx.font = 'bold 11px system-ui';
                ctx.textAlign = 'center';
                ctx.fillText(nodeId + '', x, y + 4);
            });
        }

        window.addEventListener('resize', drawGraph);
        setTimeout(drawGraph, 50);
    </script>
</body>
</html>
"""
end

function to_html(res::GraphAlignmentResult)
    return visualize_graph_alignment_html(res)
end

"""
    visualize_sequence_graph_html(graph::SequenceGraph; title="Sequence Graph Structure") -> String

Generate an interactive HTML5/Canvas DAG sequence graph structure viewer.
"""
function visualize_sequence_graph_html(graph::SequenceGraph; title::String="Sequence Graph Structure")
    node_count = length(graph.nodes)
    edge_count = length(graph.edges)
    graph_json = _sequence_graph_to_json(graph)

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(_escape_html(title))</title>
    <style>
        :root {
            --bg: #0f172a; --panel: #1e293b; --border: #334155;
            --text: #f8fafc; --muted: #94a3b8; --accent: #38bdf8;
        }
        body { font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 24px; }
        .card { background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 20px; margin-bottom: 20px; box-shadow: 0 4px 12px rgba(0,0,0,0.3); }
        .title { font-size: 1.5rem; font-weight: 700; color: var(--accent); margin: 0 0 12px 0; }
        .meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 14px; margin-bottom: 20px; }
        .meta-item { background: rgba(15,23,42,0.6); padding: 10px 14px; border-radius: 8px; border: 1px solid var(--border); }
        .meta-label { font-size: 0.75rem; text-transform: uppercase; color: var(--muted); }
        .meta-val { font-size: 1.1rem; font-weight: 600; color: var(--text); margin-top: 2px; }
        canvas { width: 100%; height: 380px; display: block; border-radius: 8px; background: #0b1329; }
    </style>
</head>
<body>
    <div class="card">
        <div class="title">🧬 $(_escape_html(title))</div>
        <div class="meta-grid">
            <div class="meta-item"><div class="meta-label">Nodes</div><div class="meta-val">$(node_count)</div></div>
            <div class="meta-item"><div class="meta-label">Edges</div><div class="meta-val">$(edge_count)</div></div>
        </div>

        <canvas id="graphCanvas"></canvas>
    </div>

    <script>
        const graph = $(graph_json);
        const canvas = document.getElementById('graphCanvas');
        const ctx = canvas.getContext('2d');

        function drawGraph() {
            canvas.width = canvas.parentElement.clientWidth * window.devicePixelRatio;
            canvas.height = 380 * window.devicePixelRatio;
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

            const w = canvas.parentElement.clientWidth;
            const h = 380;
            ctx.clearRect(0, 0, w, h);

            const nodes = graph.nodes || [];
            const edges = graph.edges || [];
            if (nodes.length === 0) return;

            const step = (w - 120) / Math.max(1, nodes.length - 1);

            ctx.strokeStyle = '#64748b';
            ctx.lineWidth = 2;
            edges.forEach(e => {
                const uIdx = e.u - 1;
                const vIdx = e.v - 1;
                if (uIdx >= 0 && uIdx < nodes.length && vIdx >= 0 && vIdx < nodes.length) {
                    const x1 = 60 + uIdx * step;
                    const y1 = h / 2;
                    const x2 = 60 + vIdx * step;
                    const y2 = h / 2;
                    ctx.beginPath();
                    ctx.moveTo(x1, y1);
                    ctx.lineTo(x2, y2);
                    ctx.stroke();
                }
            });

            nodes.forEach((n, idx) => {
                const x = 60 + idx * step;
                const y = h / 2;

                ctx.fillStyle = '#38bdf8';
                ctx.beginPath();
                ctx.arc(x, y, 16, 0, 2 * Math.PI);
                ctx.fill();

                ctx.fillStyle = '#0f172a';
                ctx.font = 'bold 11px system-ui';
                ctx.textAlign = 'center';
                ctx.fillText(n.id + '', x, y + 4);

                if (n.seq) {
                    ctx.fillStyle = '#94a3b8';
                    ctx.font = '10px monospace';
                    ctx.fillText(n.seq.length > 8 ? n.seq.substring(0, 6) + '..' : n.seq, x, y + 30);
                }
            });
        }

        window.addEventListener('resize', drawGraph);
        setTimeout(drawGraph, 50);
    </script>
</body>
</html>
"""
end

function to_html(graph::SequenceGraph)
    return visualize_sequence_graph_html(graph)
end

function _profile_result_to_json(res::ProfileAlignmentResult)
    path_items = [ "{\"u\":$(p[1]),\"v\":$(p[2])}" for p in res.path ]
    conf_items = [ round(c, digits=3) for c in res.posterior_confidence ]
    return "{\"path\":[" * join(path_items, ",") * "],\"score\":$(res.score),\"confidence\":[" * join(conf_items, ",") * "]}"
end

"""
    visualize_profile_alignment_html(res::ProfileAlignmentResult; title="Profile HMM Alignment Path") -> String

Generate an interactive HTML5/Canvas visualization for Profile HMM alignment path and confidence scores.
"""
function visualize_profile_alignment_html(res::ProfileAlignmentResult; title::String="Profile HMM Alignment Path")
    mean_conf = isempty(res.posterior_confidence) ? 0.0 : round(sum(res.posterior_confidence) / length(res.posterior_confidence), digits=3)
    prof_json = _profile_result_to_json(res)

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(_escape_html(title))</title>
    <style>
        :root {
            --bg: #0f172a; --panel: #1e293b; --border: #334155;
            --text: #f8fafc; --muted: #94a3b8; --accent: #38bdf8;
        }
        body { font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 24px; }
        .card { background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 20px; margin-bottom: 20px; box-shadow: 0 4px 12px rgba(0,0,0,0.3); }
        .title { font-size: 1.5rem; font-weight: 700; color: var(--accent); margin: 0 0 12px 0; }
        .meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 14px; margin-bottom: 20px; }
        .meta-item { background: rgba(15,23,42,0.6); padding: 10px 14px; border-radius: 8px; border: 1px solid var(--border); }
        .meta-label { font-size: 0.75rem; text-transform: uppercase; color: var(--muted); }
        .meta-val { font-size: 1.1rem; font-weight: 600; color: var(--text); margin-top: 2px; }
        canvas { width: 100%; height: 320px; display: block; border-radius: 8px; background: #0b1329; }
    </style>
</head>
<body>
    <div class="card">
        <div class="title">📊 $(_escape_html(title))</div>
        <div class="meta-grid">
            <div class="meta-item"><div class="meta-label">Profile Score</div><div class="meta-val">$(round(res.score, digits=2))</div></div>
            <div class="meta-item"><div class="meta-label">Mean Confidence</div><div class="meta-val">$(mean_conf)</div></div>
            <div class="meta-item"><div class="meta-label">Path Length</div><div class="meta-val">$(length(res.path))</div></div>
        </div>

        <canvas id="profCanvas"></canvas>
    </div>

    <script>
        const prof = $(prof_json);
        const canvas = document.getElementById('profCanvas');
        const ctx = canvas.getContext('2d');

        function drawProfile() {
            canvas.width = canvas.parentElement.clientWidth * window.devicePixelRatio;
            canvas.height = 320 * window.devicePixelRatio;
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

            const w = canvas.parentElement.clientWidth;
            const h = 320;
            ctx.clearRect(0, 0, w, h);

            const conf = prof.confidence;
            if (conf.length === 0) return;

            const step = (w - 80) / conf.length;

            ctx.fillStyle = '#38bdf8';
            conf.forEach((c, i) => {
                const x = 40 + i * step;
                const barH = c * (h - 80);
                ctx.fillRect(x, h - 40 - barH, Math.max(step - 2, 2), barH);
            });

            ctx.fillStyle = '#94a3b8';
            ctx.font = '11px system-ui';
            ctx.fillText('Alignment Position →', w / 2 - 50, h - 10);
        }

        window.addEventListener('resize', drawProfile);
        setTimeout(drawProfile, 50);
    </script>
</body>
</html>
"""
end

function to_html(res::ProfileAlignmentResult)
    return visualize_profile_alignment_html(res)
end

function to_html(profile::AlignmentProfileHMM)
    pos_count = size(profile.match_emissions, 1)
    alpha_len = length(profile.alphabet)
    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Profile HMM Model</title>
    <style>
        :root { --bg: #0f172a; --panel: #1e293b; --border: #334155; --text: #f8fafc; --accent: #38bdf8; }
        body { font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); padding: 24px; }
        .card { background: var(--panel); border: 1px solid var(--border); border-radius: 12px; padding: 20px; }
        .title { font-size: 1.4rem; font-weight: 700; color: var(--accent); margin-bottom: 12px; }
    </style>
</head>
<body>
    <div class="card">
        <div class="title">📊 Alignment Profile HMM</div>
        <p>Positions: <strong>$(pos_count)</strong> | Alphabet size: <strong>$(alpha_len)</strong></p>
    </div>
</body>
</html>
"""
end

# ==============================================================================
# CIGAR Engine, Coordinate Mapping, Banded & Specialized Alignments
# ==============================================================================

"""
    cigar(left_seq, right_seq) -> String
    cigar(alignment::PairwiseAlignmentResult) -> String

Generate a SAM-compliant CIGAR string from an alignment (e.g. `"10M2I5M3D"`).
"""
function cigar(left_seq::BioSequence, right_seq::BioSequence)
    length(left_seq) == length(right_seq) || throw(ArgumentError("aligned sequence lengths must match"))
    l_bytes = left_seq.data
    r_bytes = right_seq.data
    len = length(l_bytes)
    
    len == 0 && return ""
    
    ops = Char[]
    lens = Int[]
    
    @inbounds for i in 1:len
        l = l_bytes[i]
        r = r_bytes[i]
        op = (l != UInt8('-') && r != UInt8('-')) ? 'M' : (l == UInt8('-') ? 'D' : 'I')
        if !isempty(ops) && ops[end] == op
            lens[end] += 1
        else
            push!(ops, op)
            push!(lens, 1)
        end
    end
    
    buf = IOBuffer()
    for (op, l) in zip(ops, lens)
        print(buf, l, op)
    end
    return String(take!(buf))
end

cigar(aln::PairwiseAlignmentResult) = cigar(aln.left, aln.right)

"""
    parse_cigar(cigar_str::AbstractString) -> Vector{Tuple{Char, Int}}

Parse a CIGAR string (e.g. `"10M2I5M"`) into a vector of `(operation, length)` tuples.
"""
function parse_cigar(cigar_str::AbstractString)
    ops = Tuple{Char, Int}[]
    for m in eachmatch(r"(\d+)([MIDNSHP=X])", String(cigar_str))
        len = parse(Int, m.captures[1])
        op = m.captures[2][1]
        push!(ops, (op, len))
    end
    return ops
end

"""
    count_matches(aln::PairwiseAlignmentResult) -> Int
"""
count_matches(aln::PairwiseAlignmentResult) = aln.matches

"""
    count_mismatches(aln::PairwiseAlignmentResult) -> Int
"""
function count_mismatches(aln::PairwiseAlignmentResult)
    l_bytes = aln.left.data
    r_bytes = aln.right.data
    count = 0
    @inbounds for i in 1:length(l_bytes)
        if l_bytes[i] != UInt8('-') && r_bytes[i] != UInt8('-') && l_bytes[i] != r_bytes[i]
            count += 1
        end
    end
    return count
end

"""
    count_insertions(aln::PairwiseAlignmentResult) -> Int
"""
function count_insertions(aln::PairwiseAlignmentResult)
    r_bytes = aln.right.data
    count = 0
    @inbounds for i in 1:length(r_bytes)
        if r_bytes[i] == UInt8('-')
            count += 1
        end
    end
    return count
end

"""
    count_deletions(aln::PairwiseAlignmentResult) -> Int
"""
function count_deletions(aln::PairwiseAlignmentResult)
    l_bytes = aln.left.data
    count = 0
    @inbounds for i in 1:length(l_bytes)
        if l_bytes[i] == UInt8('-')
            count += 1
        end
    end
    return count
end

"""
    count_aligned(aln::PairwiseAlignmentResult) -> Int
"""
count_aligned(aln::PairwiseAlignmentResult) = length(aln.left)

# --- Alignment Coordinate Mapping Infrastructure ---

"""
    seq2ref(aln::PairwiseAlignmentResult, seq_pos::Int) -> Int
Map a 1-based index in the query sequence to the corresponding index in the reference sequence.
"""
function seq2ref(aln::PairwiseAlignmentResult, seq_pos::Int)
    l_bytes = aln.left.data
    r_bytes = aln.right.data
    curr_seq = 0
    curr_ref = 0
    @inbounds for i in 1:length(l_bytes)
        if l_bytes[i] != UInt8('-')
            curr_seq += 1
        end
        if r_bytes[i] != UInt8('-')
            curr_ref += 1
        end
        if curr_seq == seq_pos
            return r_bytes[i] == UInt8('-') ? 0 : curr_ref
        end
    end
    return 0
end

"""
    ref2seq(aln::PairwiseAlignmentResult, ref_pos::Int) -> Int
Map a 1-based index in the reference sequence to the corresponding index in the query sequence.
"""
function ref2seq(aln::PairwiseAlignmentResult, ref_pos::Int)
    l_bytes = aln.left.data
    r_bytes = aln.right.data
    curr_seq = 0
    curr_ref = 0
    @inbounds for i in 1:length(r_bytes)
        if l_bytes[i] != UInt8('-')
            curr_seq += 1
        end
        if r_bytes[i] != UInt8('-')
            curr_ref += 1
        end
        if curr_ref == ref_pos
            return l_bytes[i] == UInt8('-') ? 0 : curr_seq
        end
    end
    return 0
end

"""
    seq2aln(aln::PairwiseAlignmentResult, seq_pos::Int) -> Int
"""
function seq2aln(aln::PairwiseAlignmentResult, seq_pos::Int)
    l_bytes = aln.left.data
    curr_seq = 0
    @inbounds for i in 1:length(l_bytes)
        if l_bytes[i] != UInt8('-')
            curr_seq += 1
            if curr_seq == seq_pos
                return i
            end
        end
    end
    return 0
end

"""
    ref2aln(aln::PairwiseAlignmentResult, ref_pos::Int) -> Int
"""
function ref2aln(aln::PairwiseAlignmentResult, ref_pos::Int)
    r_bytes = aln.right.data
    curr_ref = 0
    @inbounds for i in 1:length(r_bytes)
        if r_bytes[i] != UInt8('-')
            curr_ref += 1
            if curr_ref == ref_pos
                return i
            end
        end
    end
    return 0
end

"""
    aln2seq(aln::PairwiseAlignmentResult, aln_pos::Int) -> Int
"""
function aln2seq(aln::PairwiseAlignmentResult, aln_pos::Int)
    l_bytes = aln.left.data
    aln_pos <= length(l_bytes) || return 0
    return count(b -> b != UInt8('-'), l_bytes[1:aln_pos])
end

"""
    aln2ref(aln::PairwiseAlignmentResult, aln_pos::Int) -> Int
"""
function aln2ref(aln::PairwiseAlignmentResult, aln_pos::Int)
    r_bytes = aln.right.data
    aln_pos <= length(r_bytes) || return 0
    return count(b -> b != UInt8('-'), r_bytes[1:aln_pos])
end

# --- Banded Needleman-Wunsch Alignment ($O(k * N)$ complexity) ---

"""
    banded_needleman_wunsch(seq1, seq2; k=50, match=2, mismatch=-1, gap_open=-5, gap_extend=-1)

Compute global alignment using a band width parameter `k`, reducing computation from O(M*N) to O(k*N).
"""
function banded_needleman_wunsch(seq1::BioSequence{A}, seq2::BioSequence{A}; k::Int=50, match::Int=2, mismatch::Int=-1, gap_open::Int=-5, gap_extend::Int=-1, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx)) where {A <: BioAlphabet}
    s1 = seq1.data
    s2 = seq2.data
    m = length(s1)
    n = length(s2)
    
    if abs(m - n) > k
        return needleman_wunsch(seq1, seq2; match=match, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend)
    end
    
    INF = -1000000000
    M_mat = fill(INF, m + 1, n + 1)
    I_mat = fill(INF, m + 1, n + 1)
    D_mat = fill(INF, m + 1, n + 1)
    
    M_mat[1, 1] = 0
    for i in 1:min(m, k)
        I_mat[i+1, 1] = gap_open + (i - 1) * gap_extend
    end
    for j in 1:min(n, k)
        D_mat[1, j+1] = gap_open + (j - 1) * gap_extend
    end
    
    for i in 1:m
        min_j = max(1, i - k)
        max_j = min(n, i + k)
        for j in min_j:max_j
            s_match = (s1[i] == s2[j]) ? match : mismatch
            prev_best = max(M_mat[i, j], I_mat[i, j], D_mat[i, j])
            if prev_best != INF
                M_mat[i+1, j+1] = prev_best + s_match
            end
            
            from_m = M_mat[i, j+1] != INF ? M_mat[i, j+1] + gap_open : INF
            from_i = I_mat[i, j+1] != INF ? I_mat[i, j+1] + gap_extend : INF
            from_d = D_mat[i, j+1] != INF ? D_mat[i, j+1] + gap_open : INF
            I_mat[i+1, j+1] = max(from_m, from_i, from_d)
            
            from_m_d = M_mat[i+1, j] != INF ? M_mat[i+1, j] + gap_open : INF
            from_d_d = D_mat[i+1, j] != INF ? D_mat[i+1, j] + gap_extend : INF
            from_i_d = I_mat[i+1, j] != INF ? I_mat[i+1, j] + gap_open : INF
            D_mat[i+1, j+1] = max(from_m_d, from_d_d, from_i_d)
        end
    end
    
    val = max(M_mat[m+1, n+1], I_mat[m+1, n+1], D_mat[m+1, n+1])
    state = val == M_mat[m+1, n+1] ? :M : (val == I_mat[m+1, n+1] ? :I : :D)
    
    i, j = m, n
    al1 = UInt8[]
    al2 = UInt8[]
    matches = 0
    
    while i > 0 || j > 0
        if state == :M && i > 0 && j > 0
            push!(al1, s1[i])
            push!(al2, s2[j])
            s1[i] == s2[j] && (matches += 1)
            prev_best = max(M_mat[i, j], I_mat[i, j], D_mat[i, j])
            state = (prev_best == M_mat[i, j]) ? :M : ((prev_best == I_mat[i, j]) ? :I : :D)
            i -= 1; j -= 1
        elseif state == :I && i > 0
            push!(al1, s1[i])
            push!(al2, UInt8('-'))
            if I_mat[i+1, j+1] == (I_mat[i, j+1] != INF ? I_mat[i, j+1] + gap_extend : INF)
                state = :I
            elseif I_mat[i+1, j+1] == (M_mat[i, j+1] != INF ? M_mat[i, j+1] + gap_open : INF)
                state = :M
            else
                state = :D
            end
            i -= 1
        elseif state == :D && j > 0
            push!(al1, UInt8('-'))
            push!(al2, s2[j])
            if D_mat[i+1, j+1] == (D_mat[i+1, j] != INF ? D_mat[i+1, j] + gap_extend : INF)
                state = :D
            elseif D_mat[i+1, j+1] == (M_mat[i+1, j] != INF ? M_mat[i+1, j] + gap_open : INF)
                state = :M
            else
                state = :I
            end
            j -= 1
        else
            if i > 0
                push!(al1, s1[i]); push!(al2, UInt8('-')); i -= 1
            elseif j > 0
                push!(al1, UInt8('-')); push!(al2, s2[j]); j -= 1
            else
                break
            end
        end
    end
    
    reverse!(al1)
    reverse!(al2)
    final_score = val
    aln_len = length(al1)
    identity = aln_len > 0 ? matches / aln_len : 0.0
    
    res = PairwiseAlignmentResult(BioSequence{A}(al1), BioSequence{A}(al2), final_score, matches, identity)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, res, "banded_needleman_wunsch")
end

banded_needleman_wunsch(seq1::AbstractString, seq2::AbstractString; kwargs...) =
    banded_needleman_wunsch(DNASeq(seq1; validate=false), DNASeq(seq2; validate=false); kwargs...)

# --- Semi-Global and Overlap Alignment Modes ---

"""
    semi_global_align(seq1, seq2; match=2, mismatch=-1, gap_open=-5, gap_extend=-1)

Perform semi-global alignment allowing free end-gaps at sequence boundaries.
"""
function semi_global_align(seq1::BioSequence{A}, seq2::BioSequence{A}; match::Int=2, mismatch::Int=-1, gap_open::Int=-5, gap_extend::Int=-1, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx)) where {A <: BioAlphabet}
    s1 = seq1.data
    s2 = seq2.data
    m = length(s1)
    n = length(s2)
    
    INF = -1000000000
    M_mat = fill(INF, m + 1, n + 1)
    I_mat = fill(INF, m + 1, n + 1)
    D_mat = fill(INF, m + 1, n + 1)
    
    M_mat[1, 1] = 0
    for i in 1:m; I_mat[i+1, 1] = 0; end
    for j in 1:n; D_mat[1, j+1] = 0; end
    
    for i in 1:m
        for j in 1:n
            s_match = (s1[i] == s2[j]) ? match : mismatch
            prev_best = max(M_mat[i, j], I_mat[i, j], D_mat[i, j])
            if prev_best != INF
                M_mat[i+1, j+1] = prev_best + s_match
            end
            
            cost_open = (i == m) ? 0 : gap_open
            cost_ext  = (i == m) ? 0 : gap_extend
            from_m = M_mat[i, j+1] != INF ? M_mat[i, j+1] + cost_open : INF
            from_i = I_mat[i, j+1] != INF ? I_mat[i, j+1] + cost_ext  : INF
            from_d = D_mat[i, j+1] != INF ? D_mat[i, j+1] + cost_open : INF
            I_mat[i+1, j+1] = max(from_m, from_i, from_d)
            
            cost_open_d = (j == n) ? 0 : gap_open
            cost_ext_d  = (j == n) ? 0 : gap_extend
            from_m_d = M_mat[i+1, j] != INF ? M_mat[i+1, j] + cost_open_d : INF
            from_d_d = D_mat[i+1, j] != INF ? D_mat[i+1, j] + cost_ext_d  : INF
            from_i_d = I_mat[i+1, j] != INF ? I_mat[i+1, j] + cost_open_d : INF
            D_mat[i+1, j+1] = max(from_m_d, from_d_d, from_i_d)
        end
    end
    
    max_score = INF
    best_i, best_j = m, n
    state = :M
    for j in 0:n
        val = max(M_mat[m+1, j+1], I_mat[m+1, j+1], D_mat[m+1, j+1])
        if val > max_score
            max_score = val
            best_i, best_j = m, j
            state = val == M_mat[m+1, j+1] ? :M : (val == I_mat[m+1, j+1] ? :I : :D)
        end
    end
    for i in 0:m
        val = max(M_mat[i+1, n+1], I_mat[i+1, n+1], D_mat[i+1, n+1])
        if val > max_score
            max_score = val
            best_i, best_j = i, n
            state = val == M_mat[i+1, n+1] ? :M : (val == I_mat[i+1, n+1] ? :I : :D)
        end
    end
    
    al1 = UInt8[]
    al2 = UInt8[]
    matches = 0
    
    for i in m:-1:best_i+1
        push!(al1, s1[i]); push!(al2, UInt8('-'))
    end
    for j in n:-1:best_j+1
        push!(al1, UInt8('-')); push!(al2, s2[j])
    end
    
    i, j = best_i, best_j
    while i > 0 || j > 0
        if state == :M && i > 0 && j > 0
            push!(al1, s1[i]); push!(al2, s2[j])
            s1[i] == s2[j] && (matches += 1)
            prev_best = max(M_mat[i, j], I_mat[i, j], D_mat[i, j])
            state = (prev_best == M_mat[i, j]) ? :M : ((prev_best == I_mat[i, j]) ? :I : :D)
            i -= 1; j -= 1
        elseif state == :I && i > 0
            push!(al1, s1[i]); push!(al2, UInt8('-'))
            cost_open = (i == m) ? 0 : gap_open
            cost_ext  = (i == m) ? 0 : gap_extend
            if I_mat[i+1, j+1] == (I_mat[i, j+1] != INF ? I_mat[i, j+1] + cost_ext : INF)
                state = :I
            elseif I_mat[i+1, j+1] == (M_mat[i, j+1] != INF ? M_mat[i, j+1] + cost_open : INF)
                state = :M
            else
                state = :D
            end
            i -= 1
        elseif state == :D && j > 0
            push!(al1, UInt8('-')); push!(al2, s2[j])
            cost_open_d = (j == n) ? 0 : gap_open
            cost_ext_d  = (j == n) ? 0 : gap_extend
            if D_mat[i+1, j+1] == (D_mat[i+1, j] != INF ? D_mat[i+1, j] + cost_ext_d : INF)
                state = :D
            elseif D_mat[i+1, j+1] == (M_mat[i+1, j] != INF ? M_mat[i+1, j] + cost_open_d : INF)
                state = :M
            else
                state = :I
            end
            j -= 1
        else
            if i > 0
                push!(al1, s1[i]); push!(al2, UInt8('-')); i -= 1
            elseif j > 0
                push!(al1, UInt8('-')); push!(al2, s2[j]); j -= 1
            else
                break
            end
        end
    end
    
    reverse!(al1)
    reverse!(al2)
    aln_len = length(al1)
    identity = aln_len > 0 ? matches / aln_len : 0.0
    
    res = PairwiseAlignmentResult(BioSequence{A}(al1), BioSequence{A}(al2), max_score, matches, identity)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, res, "semi_global_align")
end

semi_global_align(seq1::AbstractString, seq2::AbstractString; kwargs...) =
    semi_global_align(DNASeq(seq1; validate=false), DNASeq(seq2; validate=false); kwargs...)

"""
    overlap_align(seq1, seq2; match=2, mismatch=-1, gap_open=-5, gap_extend=-1)

Perform overlap alignment for sequence read assembly.
"""
function overlap_align(seq1::BioSequence{A}, seq2::BioSequence{A}; match::Int=2, mismatch::Int=-1, gap_open::Int=-5, gap_extend::Int=-1, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx)) where {A <: BioAlphabet}
    s1 = seq1.data
    s2 = seq2.data
    m = length(s1)
    n = length(s2)
    
    INF = -1000000000
    M_mat = fill(INF, m + 1, n + 1)
    I_mat = fill(INF, m + 1, n + 1)
    D_mat = fill(INF, m + 1, n + 1)
    
    M_mat[1, 1] = 0
    for i in 1:m; I_mat[i+1, 1] = 0; end
    for j in 1:n; D_mat[1, j+1] = 0; end
    
    for i in 1:m
        for j in 1:n
            s_match = (s1[i] == s2[j]) ? match : mismatch
            prev_best = max(M_mat[i, j], I_mat[i, j], D_mat[i, j])
            if prev_best != INF
                M_mat[i+1, j+1] = prev_best + s_match
            end
            
            # Free end-gap logic for terminal overhangs
            cost_open = (i == m) ? 0 : gap_open
            cost_ext  = (i == m) ? 0 : gap_extend
            from_m = M_mat[i, j+1] != INF ? M_mat[i, j+1] + cost_open : INF
            from_i = I_mat[i, j+1] != INF ? I_mat[i, j+1] + cost_ext  : INF
            from_d = D_mat[i, j+1] != INF ? D_mat[i, j+1] + cost_open : INF
            I_mat[i+1, j+1] = max(from_m, from_i, from_d)
            
            cost_open_d = (j == n) ? 0 : gap_open
            cost_ext_d  = (j == n) ? 0 : gap_extend
            from_m_d = M_mat[i+1, j] != INF ? M_mat[i+1, j] + cost_open_d : INF
            from_d_d = D_mat[i+1, j] != INF ? D_mat[i+1, j] + cost_ext_d  : INF
            from_i_d = I_mat[i+1, j] != INF ? I_mat[i+1, j] + cost_open_d : INF
            D_mat[i+1, j+1] = max(from_m_d, from_d_d, from_i_d)
        end
    end
    
    max_score = INF
    best_i, best_j = m, n
    state = :M
    for j in 1:n
        val = max(M_mat[m+1, j+1], I_mat[m+1, j+1], D_mat[m+1, j+1])
        if val > max_score
            max_score = val
            best_i, best_j = m, j
            state = val == M_mat[m+1, j+1] ? :M : (val == I_mat[m+1, j+1] ? :I : :D)
        end
    end
    for i in 1:m
        val = max(M_mat[i+1, n+1], I_mat[i+1, n+1], D_mat[i+1, n+1])
        if val > max_score
            max_score = val
            best_i, best_j = i, n
            state = val == M_mat[i+1, n+1] ? :M : (val == I_mat[i+1, n+1] ? :I : :D)
        end
    end
    
    al1 = UInt8[]
    al2 = UInt8[]
    matches = 0
    
    for i in m:-1:best_i+1
        push!(al1, s1[i]); push!(al2, UInt8('-'))
    end
    for j in n:-1:best_j+1
        push!(al1, UInt8('-')); push!(al2, s2[j])
    end
    
    i, j = best_i, best_j
    while i > 0 || j > 0
        if state == :M && i > 0 && j > 0
            push!(al1, s1[i]); push!(al2, s2[j])
            s1[i] == s2[j] && (matches += 1)
            prev_best = max(M_mat[i, j], I_mat[i, j], D_mat[i, j])
            state = (prev_best == M_mat[i, j]) ? :M : ((prev_best == I_mat[i, j]) ? :I : :D)
            i -= 1; j -= 1
        elseif state == :I && i > 0
            push!(al1, s1[i]); push!(al2, UInt8('-'))
            c_open = (i == m) ? 0 : gap_open
            c_ext  = (i == m) ? 0 : gap_extend
            if I_mat[i+1, j+1] == (I_mat[i, j+1] != INF ? I_mat[i, j+1] + c_ext : INF)
                state = :I
            elseif I_mat[i+1, j+1] == (M_mat[i, j+1] != INF ? M_mat[i, j+1] + c_open : INF)
                state = :M
            else
                state = :D
            end
            i -= 1
        elseif state == :D && j > 0
            push!(al1, UInt8('-')); push!(al2, s2[j])
            c_open = (j == n) ? 0 : gap_open
            c_ext  = (j == n) ? 0 : gap_extend
            if D_mat[i+1, j+1] == (D_mat[i+1, j] != INF ? D_mat[i+1, j] + c_ext : INF)
                state = :D
            elseif D_mat[i+1, j+1] == (M_mat[i+1, j] != INF ? M_mat[i+1, j] + c_open : INF)
                state = :M
            else
                state = :I
            end
            j -= 1
        else
            if i > 0
                push!(al1, s1[i]); push!(al2, UInt8('-')); i -= 1
            elseif j > 0
                push!(al1, UInt8('-')); push!(al2, s2[j]); j -= 1
            else
                break
            end
        end
    end
    
    reverse!(al1)
    reverse!(al2)
    aln_len = length(al1)
    identity = aln_len > 0 ? matches / aln_len : 0.0
    
    res = PairwiseAlignmentResult(BioSequence{A}(al1), BioSequence{A}(al2), max_score, matches, identity)
    _ctx = active_provenance_context()
    return provenance_result!(_ctx, res, "overlap_align")
end

overlap_align(seq1::AbstractString, seq2::AbstractString; kwargs...) =
    overlap_align(DNASeq(seq1; validate=false), DNASeq(seq2; validate=false); kwargs...)

"""
    differentiable_align(seq1::BioSequence, seq2::BioSequence; kwargs...)

Differentiable dynamic programming alignment wrapper.
"""
function differentiable_align(seq1::BioSequence{A}, seq2::BioSequence{A}; scoring::DifferentiableScoring=DifferentiableScoring([1.0, -1.0, -1.0])) where {A<:BioAlphabet}
    score = soft_alignment_score(seq1, seq2, scoring)
    return PairwiseAlignmentResult(seq1, seq2, round(Int, score), 0, 0.0)
end
differentiable_align(seq1::AbstractString, seq2::AbstractString; kwargs...) =
    differentiable_align(DNASeq(seq1; validate=false), DNASeq(seq2; validate=false); kwargs...)

# --- Unified Pairwise Alignment Dispatcher ---

"""
    pairalign(mode::Symbol, seq1, seq2; kwargs...)

Unified dispatcher for pairwise alignment:
- `:global`: Needleman-Wunsch global alignment
- `:semi_global`: Semi-global alignment with free end gaps
- `:overlap`: Overlap alignment
- `:local`: Smith-Waterman local alignment
- `:banded`: Banded global alignment
- `:pair_hmm`: Pair-HMM probabilistic posterior alignment
- `:differentiable`: Differentiable dynamic programming alignment
- `:codon`: Codon-aware frame-preserving alignment
"""
function pairalign(mode::Symbol, seq1, seq2; kwargs...)
    if mode === :global || mode === :needleman_wunsch
        return needleman_wunsch(seq1, seq2; kwargs...)
    elseif mode === :semi_global
        return semi_global_align(seq1, seq2; kwargs...)
    elseif mode === :overlap
        return overlap_align(seq1, seq2; kwargs...)
    elseif mode === :local || mode === :smith_waterman
        return smith_waterman(seq1, seq2; kwargs...)
    elseif mode === :banded
        return banded_needleman_wunsch(seq1, seq2; kwargs...)
    elseif mode === :pair_hmm
        return pairhmm_align(seq1, seq2; kwargs...)
    elseif mode === :differentiable
        return differentiable_align(seq1, seq2; kwargs...)
    elseif mode === :codon
        return pairwise_align_codons(seq1, seq2; kwargs...)
    else
        throw(ArgumentError("Unsupported pairalign mode: $mode"))
    end
end

