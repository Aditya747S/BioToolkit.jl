# ==============================================================================
# flowcytometry.jl — Flow Cytometry and CyTOF analysis
# ==============================================================================

module FlowCytometry

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, active_provenance_context, provenance_result!
import ..BioToolkit: analysis_result_fields
using LinearAlgebra
using Statistics
using Random

export FlowExperiment, GateResult, FlowSOMResult, QuadrantGateResult
export read_fcs, mock_flow_experiment, compensate_fcs, parse_spillover, rectangle_gate, quadrant_gate, polygon_gate, ellipsoid_gate, apply_gate, flowsom_cluster, xshift_cluster
export arcsinh_transform, logicle_transform, inverse_logicle_transform, estimate_logicle_params, hlog_transform
export fcs_stats, fcs_correlation

struct FlowExperiment <: AbstractAnalysisResult
  events::Matrix{Float64}
  channels::Vector{String}
  metadata::Dict{String,Any}
  provenance::ResultProvenance
end

struct GateResult <: AbstractAnalysisResult
  gate_type::Symbol
  indices::Vector{Int}
  provenance::ResultProvenance
end

struct QuadrantGateResult <: AbstractAnalysisResult
  quadrants::NTuple{4, Vector{Int}}
  indices::Vector{Int}
  provenance::ResultProvenance
end

struct FlowSOMResult <: AbstractAnalysisResult
  metaclusters::Vector{Int}
  codes::Matrix{Float64}
  cluster_centers::Matrix{Float64}
  provenance::ResultProvenance
end

analysis_result_fields(::Type{FlowExperiment}) = (:events, :channels, :metadata)
analysis_result_fields(::Type{GateResult}) = (:gate_type, :indices)
analysis_result_fields(::Type{QuadrantGateResult}) = (:quadrants, :indices)
analysis_result_fields(::Type{FlowSOMResult}) = (:metaclusters, :codes, :cluster_centers)

@inline function _parse_fcs_offset(header::String, range)
  txt = strip(header[range])
  isempty(txt) && return 0
  return parse(Int, txt)
end

@inline function _parse_fcs_text(text::String)
  isempty(text) && throw(ArgumentError("empty FCS TEXT segment"))
  delim = text[1]
  parts = split(text[2:end], delim; keepempty=true)
  meta = Dict{String,Any}()
  i = 1
  while i < length(parts)
    key = String(parts[i])
    val = String(parts[i+1])
    !isempty(key) && (meta[key] = val)
    i += 2
  end
  return meta
end

@inline function _fcs_int(bytes::Vector{UInt8}, byteord::String)
  if startswith(byteord, "4,3,2,1") || startswith(byteord, "2,1")
    bytes = reverse(bytes)
  end
  value = UInt64(0)
  for (i, b) in enumerate(bytes)
    value |= UInt64(b) << (8 * (i - 1))
  end
  return Float64(value)
end

function _read_fcs_data(raw::Vector{UInt8}, meta::Dict{String,Any}, n_events::Int, n_channels::Int)
  datatype = uppercase(String(get(meta, "\$DATATYPE", "")))
  byteord = String(get(meta, "\$BYTEORD", "1,2,3,4"))
  events = Matrix{Float64}(undef, n_events, n_channels)
  cursor = 1
  if datatype == "A"
    rows = split(String(raw), ['\n', '\r']; keepempty=false)
    length(rows) >= n_events || throw(ArgumentError("FCS ASCII data has fewer rows than \$TOT"))
    for i in 1:n_events
      vals = split(strip(rows[i]))
      length(vals) >= n_channels || throw(ArgumentError("FCS ASCII data row $i has fewer columns than \$PAR"))
      for j in 1:n_channels
        events[i, j] = parse(Float64, vals[j])
      end
    end
  elseif datatype in ("F", "D")
    width = datatype == "F" ? 4 : 8
    expected = n_events * n_channels * width
    length(raw) >= expected || throw(ArgumentError("FCS DATA segment is shorter than expected for datatype $datatype"))
    raw_copy = copy(raw[1:expected])
    if startswith(byteord, "4,3,2,1") || startswith(byteord, "2,1")
      for i in 1:width:length(raw_copy)
        reverse!(@view raw_copy[i:min(i+width-1, end)])
      end
    end
    if datatype == "F"
      vec = reinterpret(Float32, raw_copy)
      # FCS binary list-mode data is event-major (Event1Ch1, Event1Ch2, ..., EventNChK).
      # Julia matrices are column-major. Reshape to (n_channels, n_events) and transpose
      # to obtain the correct (n_events, n_channels) matrix without data scrambling.
      events .= transpose(reshape(Float64.(vec), n_channels, n_events))
    else
      vec = reinterpret(Float64, raw_copy)
      events .= transpose(reshape(vec, n_channels, n_events))
    end
  elseif datatype == "I"
    widths = [parse(Int, String(get(meta, "\$P$(j)B", "16"))) ÷ 8 for j in 1:n_channels]
    for i in 1:n_events, j in 1:n_channels
      width = widths[j]
      width > 0 || throw(ArgumentError("invalid integer bit width for FCS parameter $j"))
      cursor + width - 1 <= length(raw) || throw(ArgumentError("FCS DATA segment ended while reading event $i channel $j"))
      events[i, j] = _fcs_int(raw[cursor:(cursor+width-1)], byteord)
      cursor += width
    end
  else
    throw(ArgumentError("unsupported or missing FCS \$DATATYPE '$datatype'"))
  end
  return events
end

function mock_flow_experiment(; n_events::Int=1000, channels::Vector{String}=["FSC-A", "SSC-A", "FITC-A", "PE-A", "APC-A"], seed::Int=1)
  _ctx = active_provenance_context()
  rng = MersenneTwister(seed)
  events = randn(rng, n_events, length(channels)) .* 2.0 .+ 5.0
  length(channels) >= 1 && (events[:, 1] .+= 10.0)
  length(channels) >= 2 && (events[:, 2] .+= 8.0)
  metadata = Dict{String,Any}("source" => "explicit_mock", "seed" => seed)
  prov_rec = provenance_record("FlowExperiment", "FlowCytometry/mock_flow_experiment"; parameters=(n_events=n_events, channel_count=length(channels), seed=seed))
  result = FlowExperiment(events, channels, metadata, prov_rec)
  return provenance_result!(_ctx, result, "mock_flow_experiment"; parents=String[], parameters=(n_events=n_events, channel_count=length(channels), seed=seed))
end

function read_fcs(path::String; on_error::Symbol=:throw, mock::Bool=false)
  _ctx = active_provenance_context()
  if mock || on_error === :mock
    return mock_flow_experiment()
  end
  on_error === :throw || throw(ArgumentError("unsupported read_fcs on_error policy: $on_error"))
  isfile(path) || throw(ArgumentError("FCS file does not exist: $path. Use mock_flow_experiment() or read_fcs(path; mock=true) for demo data."))
  open(path, "r") do io
    header = String(read(io, 58))
    startswith(header, "FCS") || throw(ArgumentError("not a valid FCS file: missing FCS version header"))
    version = strip(header[1:6])
    version in ("FCS2.0", "FCS3.0", "FCS3.1") || throw(ArgumentError("unsupported FCS version $version"))
    text_start = _parse_fcs_offset(header, 11:18)
    text_end = _parse_fcs_offset(header, 19:26)
    data_start = _parse_fcs_offset(header, 27:34)
    data_end = _parse_fcs_offset(header, 35:42)
    text_start > 0 && text_end >= text_start || throw(ArgumentError("FCS header has invalid TEXT offsets"))
    seek(io, text_start)
    meta = _parse_fcs_text(String(read(io, text_end - text_start + 1)))
    data_start == 0 && haskey(meta, "\$BEGINDATA") && (data_start = parse(Int, meta["\$BEGINDATA"]))
    data_end == 0 && haskey(meta, "\$ENDDATA") && (data_end = parse(Int, meta["\$ENDDATA"]))
    data_end >= data_start || throw(ArgumentError("FCS header has invalid DATA offsets"))
    n_events = parse(Int, String(get(meta, "\$TOT", "0")))
    n_channels = parse(Int, String(get(meta, "\$PAR", "0")))
    n_events > 0 && n_channels > 0 || throw(ArgumentError("FCS metadata must define positive \$TOT and \$PAR"))
    channels = [String(get(meta, "\$P$(i)N", "Ch$(i)")) for i in 1:n_channels]
    seek(io, data_start)
    raw = read(io, data_end - data_start + 1)
    events = _read_fcs_data(raw, meta, n_events, n_channels)
    meta["version"] = version

    # FCS 3.1 & FCS 2.0 transformation keyword parsing
    if haskey(meta, "\$TRANSFORMATION")
      trans_val = uppercase(String(meta["\$TRANSFORMATION"]))
      if trans_val in ("\$LOG", "\$LN", "LOG", "LN")
        meta["is_transformed"] = true
      end
    end
    for j in 1:n_channels
      key = "\$P$(j)E"
      if haskey(meta, key)
        val = strip(String(meta[key]))
        parts = split(val, [',', ' ']; keepempty=false)
        if length(parts) >= 1 && parse(Float64, parts[1]) > 0
          meta["is_transformed"] = true
          meta["\$P$(j)E_log"] = true
        end
      end
    end

    prov_rec = provenance_record("FlowExperiment", "FlowCytometry/read_fcs"; parameters=(path=path, version=version, n_events=n_events, n_channels=n_channels))
    result = FlowExperiment(events, channels, meta, prov_rec)
    return provenance_result!(_ctx, result, "read_fcs"; parents=String[], parameters=(path=path, version=version, n_events=n_events, n_channels=n_channels))
  end
end

"""
    parse_spillover(meta::Dict{String,Any})

Extract and parse the spillover/compensation matrix from FCS metadata keywords (`\$SPILLOVER`, `SPILL`, `\$SPILL`, `COMP`).
Returns `(spillover_matrix::Matrix{Float64}, channel_names::Vector{String})`.
"""
function parse_spillover(meta::Dict{String,Any})
  spill_str = ""
  for key in ("\$SPILLOVER", "SPILL", "\$SPILL", "COMP")
    if haskey(meta, key) && !isempty(String(meta[key]))
      spill_str = String(meta[key])
      break
    end
  end
  isempty(spill_str) && throw(ArgumentError("no spillover/compensation matrix found in metadata keywords"))

  parts = split(spill_str, [',', ';', '\t', '\n', '\r', ' ']; keepempty=false)
  length(parts) >= 2 || throw(ArgumentError("invalid spillover string format"))

  n_ch = parse(Int, parts[1])
  length(parts) >= 1 + n_ch + n_ch * n_ch || throw(ArgumentError("spillover string has insufficient elements for $n_ch channels"))

  channel_names = String.(parts[2:(1 + n_ch)])
  val_strings = parts[(2 + n_ch):(1 + n_ch + n_ch * n_ch)]

  spill_mat = zeros(Float64, n_ch, n_ch)
  idx = 1
  for i in 1:n_ch, j in 1:n_ch
    spill_mat[i, j] = parse(Float64, val_strings[idx])
    idx += 1
  end

  return spill_mat, channel_names
end

"""
    compensate_fcs(fcs::FlowExperiment)
    compensate_fcs(fcs::FlowExperiment, spillover_matrix::Matrix{Float64})

Apply spillover matrix compensation to flow cytometry channels.
If called without a matrix argument, automatically parses compensation matrix from FCS metadata keywords (`\$SPILLOVER`, `SPILL`, etc.).
Uses triangular backslash / LU system solver for maximum numerical precision and stability.
"""
function compensate_fcs(fcs::FlowExperiment)
  spillover_matrix, spill_channels = parse_spillover(fcs.metadata)
  if spill_channels == fcs.channels
    return compensate_fcs(fcs, spillover_matrix)
  else
    n_ch = length(fcs.channels)
    full_spill = Matrix{Float64}(I, n_ch, n_ch)
    for (i, c1) in enumerate(spill_channels)
      idx1 = findfirst(==(c1), fcs.channels)
      idx1 === nothing && continue
      for (j, c2) in enumerate(spill_channels)
        idx2 = findfirst(==(c2), fcs.channels)
        idx2 === nothing && continue
        full_spill[idx1, idx2] = spillover_matrix[i, j]
      end
    end
    return compensate_fcs(fcs, full_spill)
  end
end

@inline function compensate_fcs(fcs::FlowExperiment, spillover_matrix::Matrix{Float64})
  _ctx = active_provenance_context()
  n_ch = length(fcs.channels)
  size(spillover_matrix) == (n_ch, n_ch) || throw(ArgumentError("spillover_matrix dimensions must match the number of channels"))

  # Numerically stable solve (events * inv(spillover)) without explicit matrix inversion
  events = (spillover_matrix' \ fcs.events')'
  prov_rec = provenance_record("FlowExperiment", "FlowCytometry/compensate_fcs"; parameters=(n_channels=n_ch, condition_number=cond(spillover_matrix)))
  result = FlowExperiment(events, fcs.channels, copy(fcs.metadata), prov_rec)
  return provenance_result!(_ctx, result, "compensate_fcs"; parents=[fcs.provenance.id], parameters=(n_channels=n_ch,))
end

@inline function rectangle_gate(fcs::FlowExperiment, channel1::String, min1::Real, max1::Real, channel2::String, min2::Real, max2::Real)
  _ctx = active_provenance_context()
  idx1 = findfirst(==(channel1), fcs.channels)
  idx2 = findfirst(==(channel2), fcs.channels)
  idx1 === nothing && throw(ArgumentError("channel '$channel1' not found"))
  idx2 === nothing && throw(ArgumentError("channel '$channel2' not found"))

  col1 = @view fcs.events[:, idx1]
  col2 = @view fcs.events[:, idx2]
  mask = (min1 .<= col1 .<= max1) .& (min2 .<= col2 .<= max2)
  indices = findall(mask)

  prov_rec = provenance_record("GateResult", "FlowCytometry/rectangle_gate"; parameters=(channel1=channel1, min1=Float64(min1), max1=Float64(max1), channel2=channel2, min2=Float64(min2), max2=Float64(max2), n_selected=length(indices)))
  result = GateResult(:rectangle, indices, prov_rec)
  return provenance_result!(_ctx, result, "rectangle_gate"; parents=[fcs.provenance.id], parameters=(channel1=channel1, channel2=channel2, n_selected=length(indices)))
end

@inline function quadrant_gate(fcs::FlowExperiment, channel1::String, threshold1::Real, channel2::String, threshold2::Real)
  _ctx = active_provenance_context()
  idx1 = findfirst(==(channel1), fcs.channels)
  idx2 = findfirst(==(channel2), fcs.channels)
  idx1 === nothing && throw(ArgumentError("channel '$channel1' not found"))
  idx2 === nothing && throw(ArgumentError("channel '$channel2' not found"))

  col1 = @view fcs.events[:, idx1]
  col2 = @view fcs.events[:, idx2]

  q1_mask = (col1 .< threshold1) .& (col2 .< threshold2)  # LL
  q2_mask = (col1 .< threshold1) .& (col2 .>= threshold2) # UL
  q3_mask = (col1 .>= threshold1) .& (col2 .>= threshold2) # UR
  q4_mask = (col1 .>= threshold1) .& (col2 .< threshold2)  # LR

  q1_idx = findall(q1_mask)
  q2_idx = findall(q2_mask)
  q3_idx = findall(q3_mask)
  q4_idx = findall(q4_mask)
  all_indices = vcat(q1_idx, q2_idx, q3_idx, q4_idx)

  prov_rec = provenance_record("QuadrantGateResult", "FlowCytometry/quadrant_gate"; parameters=(channel1=channel1, threshold1=Float64(threshold1), channel2=channel2, threshold2=Float64(threshold2), n_q1=length(q1_idx), n_q2=length(q2_idx), n_q3=length(q3_idx), n_q4=length(q4_idx)))
  result = QuadrantGateResult((q1_idx, q2_idx, q3_idx, q4_idx), all_indices, prov_rec)
  return provenance_result!(_ctx, result, "quadrant_gate"; parents=[fcs.provenance.id], parameters=(channel1=channel1, channel2=channel2, n_q1=length(q1_idx), n_q2=length(q2_idx), n_q3=length(q3_idx), n_q4=length(q4_idx)))
end

@inline function polygon_gate(fcs::FlowExperiment, channel1::String, channel2::String, vertices::Vector{Tuple{Float64,Float64}})
  _ctx = active_provenance_context()
  idx1 = findfirst(==(channel1), fcs.channels)
  idx2 = findfirst(==(channel2), fcs.channels)
  idx1 === nothing && throw(ArgumentError("channel '$channel1' not found"))
  idx2 === nothing && throw(ArgumentError("channel '$channel2' not found"))
  length(vertices) >= 3 || throw(ArgumentError("polygon must have at least 3 vertices"))

  col1 = @view fcs.events[:, idx1]
  col2 = @view fcs.events[:, idx2]

  indices = Int[]
  @inbounds for i in 1:size(fcs.events, 1)
    x = col1[i]
    y = col2[i]
    inside = false
    for j in 1:length(vertices)
      j_next = j == length(vertices) ? 1 : j + 1
      y1 = vertices[j][2]
      y2 = vertices[j_next][2]
      if (y1 > y) != (y2 > y)
        x_intersect = vertices[j][1] + (y - y1) * (vertices[j_next][1] - vertices[j][1]) / (y2 - y1)
        if x < x_intersect
          inside = !inside
        end
      end
    end
    inside && push!(indices, i)
  end

  prov_rec = provenance_record("GateResult", "FlowCytometry/polygon_gate"; parameters=(channel1=channel1, channel2=channel2, n_vertices=length(vertices), n_selected=length(indices)))
  result = GateResult(:polygon, indices, prov_rec)
  return provenance_result!(_ctx, result, "polygon_gate"; parents=[fcs.provenance.id], parameters=(channel1=channel1, channel2=channel2, n_selected=length(indices)))
end

@inline function _chisq_quantile(df::Int, p::Float64)
  if df <= 0
    return 0.0
  end
  z = _norm_quantile(p)
  c = 1.0 / (9.0 * df)
  return df * (1.0 - c + z * sqrt(c))^3
end

@inline function _norm_quantile(p::Float64)
  if p <= 0.0 || p >= 1.0
    throw(ArgumentError("p must be in (0, 1)"))
  end
  a1 = -3.969683028665376e+01
  a2 = 2.209460984245205e+02
  a3 = -2.759285104469687e+02
  a4 = 1.383577518672690e+02
  a5 = -3.066479806614716e+01
  a6 = 2.506628277459239e+00
  b1 = -5.447609879822406e+01
  b2 = 1.615858368580409e+02
  b3 = -1.556989798598866e+02
  b4 = 6.680131188771972e+01
  b5 = -1.328068155288572e+01
  c1 = -7.784894002430293e-03
  c2 = -3.223964580411365e-01
  c3 = -2.400758277161838e+00
  c4 = -2.549732539343734e+00
  c5 = 4.374664141464968e+00
  c6 = 2.938163982698783e+00
  d1 = 7.784695709041462e-03
  d2 = 3.224671290700398e-01
  d3 = 2.445134137142996e+00
  d4 = 3.754408661907416e+00
  p_low = 0.02425
  p_high = 1.0 - p_low
  q = p
  if p < p_low
    q = sqrt(-2.0 * log(p))
    return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
           ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)
  elseif p > p_high
    q = sqrt(-2.0 * log(1.0 - p))
    return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)
  else
    q = p - 0.5
    r = q * q
    return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
           (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0)
  end
end

@inline function ellipsoid_gate(fcs::FlowExperiment, channels::Vector{String}, center::Vector{Float64}, covariance::Matrix{Float64}; confidence::Float64=0.95)
  _ctx = active_provenance_context()
  n_ch = length(channels)
  length(center) == n_ch || throw(ArgumentError("center length must match number of channels"))
  size(covariance) == (n_ch, n_ch) || throw(ArgumentError("covariance must be n_channels x n_channels"))

  idxs = [findfirst(==(c), fcs.channels) for c in channels]
  any(idx -> idx === nothing, idxs) && throw(ArgumentError("one or more channels not found"))

  # Add tiny diagonal regularization to ensure positive-definiteness under floating point noise
  cov_reg = Symmetric(covariance + 1e-10 * I)
  L = cholesky(cov_reg).L
  threshold = _chisq_quantile(n_ch, confidence)

  indices = Int[]
  diff_vec = zeros(Float64, n_ch)
  @inbounds for i in 1:size(fcs.events, 1)
    for j in 1:n_ch
      diff_vec[j] = fcs.events[i, idxs[j]] - center[j]
    end
    # Triangular backslash solver (L \ diff) avoids explicit matrix inversion for numerical stability
    y = L \ diff_vec
    mahal_sq = sum(abs2, y)
    mahal_sq <= threshold && push!(indices, i)
  end

  prov_rec = provenance_record("GateResult", "FlowCytometry/ellipsoid_gate"; parameters=(channels=channels, confidence=confidence, threshold=threshold, n_selected=length(indices)))
  result = GateResult(:ellipsoid, indices, prov_rec)
  return provenance_result!(_ctx, result, "ellipsoid_gate"; parents=[fcs.provenance.id], parameters=(channels=channels, confidence=confidence, n_selected=length(indices)))
end

@inline function apply_gate(fcs::FlowExperiment, gate::Union{GateResult, QuadrantGateResult})
  _ctx = active_provenance_context()
  sub_events = fcs.events[gate.indices, :]
  prov_rec = provenance_record("FlowExperiment", "FlowCytometry/apply_gate"; parameters=(n_events=size(sub_events, 1), gate_type=gate.gate_type))
  result = FlowExperiment(sub_events, fcs.channels, copy(fcs.metadata), prov_rec)
  return provenance_result!(_ctx, result, "apply_gate"; parents=[fcs.provenance.id, gate.provenance.id], parameters=(n_events=size(sub_events, 1),))
end

@inline function _find_bmu(event::AbstractVector{Float64}, codes::Matrix{Float64})
  n_codes = size(codes, 1)
  best_idx = 1
  min_dist = Inf
  @inbounds for c in 1:n_codes
    dist = 0.0
    for j in eachindex(event)
      d = codes[c, j] - event[j]
      dist += d * d
    end
    if dist < min_dist
      min_dist = dist
      best_idx = c
    end
  end
  return best_idx
end

@inline function _update_codes!(codes::Matrix{Float64}, event::AbstractVector{Float64}, bmu_idx::Int, grid_size::Int, lr::Float64, rad_sq::Float64)
  bmu_x = (bmu_idx - 1) % grid_size
  bmu_y = div(bmu_idx - 1, grid_size)
  n_codes = size(codes, 1)
  @inbounds for c in 1:n_codes
    cx = (c - 1) % grid_size
    cy = div(c - 1, grid_size)
    dist_sq = (cx - bmu_x)^2 + (cy - bmu_y)^2
    if dist_sq <= rad_sq
      influence = exp(-dist_sq / (2.0 * rad_sq))
      @inbounds for j in 1:size(codes, 2)
        codes[c, j] += lr * influence * (event[j] - codes[c, j])
      end
    end
  end
end

function flowsom_cluster(fcs::FlowExperiment; n_metaclusters::Int=10, grid_size::Int=10, n_epochs::Int=10, learning_rate::Float64=0.5, initial_radius::Float64=5.0, seed::Int=42)
  _ctx = active_provenance_context()

  n_codes = grid_size * grid_size
  n_features = length(fcs.channels)
  n_events = size(fcs.events, 1)

  # Normalize events
  events_norm = copy(fcs.events)
  means = mean(events_norm, dims=1)
  stds = std(events_norm, dims=1)
  for j in 1:n_features
    s = stds[j] > 0 ? stds[j] : 1.0
    @views events_norm[:, j] .-= means[j]
    @views events_norm[:, j] ./= s
  end

  # Initialize codebook
  rng = MersenneTwister(seed)
  codes = randn(rng, n_codes, n_features)

  # SOM training
  for epoch in 1:n_epochs
    lr = learning_rate * exp(-epoch / n_epochs)
    rad = initial_radius * exp(-epoch / n_epochs)
    rad_sq = rad^2

    shuffled = shuffle(rng, 1:n_events)
    for i in shuffled
      event = @view events_norm[i, :]
      bmu_idx = _find_bmu(event, codes)
      _update_codes!(codes, event, bmu_idx, grid_size, lr, rad_sq)
    end
  end

  # High-performance O(N_codes^2) Hierarchical Metaclustering with Nearest-Neighbor Tracking
  active_mask = trues(n_codes)
  cluster_members = [[i] for i in 1:n_codes]
  cluster_means = [copy(codes[i, :]) for i in 1:n_codes]

  # Distance matrix between cluster centroids
  D = zeros(Float64, n_codes, n_codes)
  @inbounds for i in 1:n_codes, j in 1:n_codes
    D[i, j] = i == j ? Inf : norm(cluster_means[i] .- cluster_means[j])
  end

  # Cache minimum distance and nearest active neighbor per cluster
  min_dist = fill(Inf, n_codes)
  nearest = zeros(Int, n_codes)
  @inbounds for i in 1:n_codes
    for j in 1:n_codes
      if j != i && D[i, j] < min_dist[i]
        min_dist[i] = D[i, j]
        nearest[i] = j
      end
    end
  end

  n_active = n_codes
  while n_active > n_metaclusters
    # Find active cluster with globally minimum distance to its nearest neighbor (O(N) search)
    best_i = 1
    best_d = Inf
    @inbounds for i in 1:n_codes
      if active_mask[i] && min_dist[i] < best_d
        best_d = min_dist[i]
        best_i = i
      end
    end

    c_keep = best_i
    c_drop = nearest[best_i]
    if c_drop == 0 || !active_mask[c_drop]
      min_d_tmp = Inf
      best_j = 0
      for j in 1:n_codes
        if j != c_keep && active_mask[j] && D[c_keep, j] < min_d_tmp
          min_d_tmp = D[c_keep, j]
          best_j = j
        end
      end
      c_drop = best_j
    end

    c_drop > 0 || break

    # Deactivate c_drop
    active_mask[c_drop] = false
    n_active -= 1

    # Merge c_drop into c_keep
    append!(cluster_members[c_keep], cluster_members[c_drop])
    n_k = length(cluster_members[c_keep])

    # Update centroid of c_keep
    fill!(cluster_means[c_keep], 0.0)
    @inbounds for m in cluster_members[c_keep]
      cluster_means[c_keep] .+= @view(codes[m, :])
    end
    cluster_means[c_keep] ./= n_k

    # Update distances between c_keep and remaining active clusters (O(N) step)
    min_dist[c_keep] = Inf
    nearest[c_keep] = 0
    @inbounds for j in 1:n_codes
      if active_mask[j] && j != c_keep
        d = norm(cluster_means[c_keep] .- cluster_means[j])
        D[c_keep, j] = d
        D[j, c_keep] = d

        if d < min_dist[c_keep]
          min_dist[c_keep] = d
          nearest[c_keep] = j
        end

        if nearest[j] == c_drop || nearest[j] == c_keep || d < min_dist[j]
          min_dist[j] = Inf
          nearest[j] = 0
          for k in 1:n_codes
            if active_mask[k] && k != j && D[j, k] < min_dist[j]
              min_dist[j] = D[j, k]
              nearest[j] = k
            end
          end
        end
      end
    end
  end

  active_clusters = findall(active_mask)
  metacluster_map = zeros(Int, n_codes)
  for (meta_id, c_idx) in enumerate(active_clusters)
    for m in cluster_members[c_idx]
      metacluster_map[m] = meta_id
    end
  end

  # Assign events to metaclusters
  event_clusters = Vector{Int}(undef, n_events)
  @inbounds for i in 1:n_events
    event = @view events_norm[i, :]
    bmu_idx = _find_bmu(event, codes)
    event_clusters[i] = metacluster_map[bmu_idx]
  end

  # Compute cluster centers
  n_final_clusters = length(active_clusters)
  cluster_centers = zeros(Float64, n_final_clusters, n_features)
  cluster_counts = zeros(Int, n_final_clusters)
  @inbounds for i in 1:n_events
    cl = event_clusters[i]
    cluster_counts[cl] += 1
    @views cluster_centers[cl, :] .+= fcs.events[i, :]
  end
  @inbounds for cl in 1:n_final_clusters
    if cluster_counts[cl] > 0
      cluster_centers[cl, :] ./= cluster_counts[cl]
    end
  end

  prov_rec = provenance_record("FlowSOMResult", "FlowCytometry/flowsom_cluster"; parameters=(n_metaclusters=n_metaclusters, grid_size=grid_size, n_epochs=n_epochs, n_events=n_events, n_features=n_features))
  result = FlowSOMResult(event_clusters, codes, cluster_centers, prov_rec)
  return provenance_result!(_ctx, result, "flowsom_cluster"; parents=[fcs.provenance.id], parameters=(n_metaclusters=n_metaclusters, grid_size=grid_size, n_epochs=n_epochs))
end

function _exact_knn(events::Matrix{Float64}, K::Int)
  n_points, n_dim = size(events)
  K = min(K, n_points - 1)
  knn_indices = [Vector{Int}(undef, K) for _ in 1:n_points]

  chunk_size = 500
  for start_i in 1:chunk_size:n_points
    end_i = min(start_i + chunk_size - 1, n_points)
    chunk_len = end_i - start_i + 1

    dists = zeros(Float64, chunk_len, n_points)
    for i_idx in 1:chunk_len
      i = start_i + i_idx - 1
      @inbounds for j in 1:n_points
        if i == j
          dists[i_idx, j] = Inf
        else
          d = 0.0
          @inbounds for k in 1:n_dim
            d += (events[i, k] - events[j, k])^2
          end
          dists[i_idx, j] = d
        end
      end
      perm = partialsortperm(@view(dists[i_idx, :]), 1:K)
      knn_indices[i] = perm
    end
  end
  return knn_indices
end

# Multi-tree Random Projection Tree (RP-Tree) with local graph refinement for high-dimensional kNN
function _approximate_knn(fcs::FlowExperiment, K::Int; seed::Int=42)
  n_points, n_dim = size(fcs.events)
  K = min(K, n_points - 1)
  rng = MersenneTwister(seed)

  if n_points <= 1000
    return _exact_knn(fcs.events, K)
  end

  n_trees = min(12, max(6, Int(ceil(log2(n_points)))))
  leaf_size = max(2 * K, 32)
  candidate_sets = [Set{Int}() for _ in 1:n_points]

  function build_rptree!(indices::Vector{Int}, depth::Int)
    if length(indices) <= leaf_size || depth >= 20
      @inbounds for i in indices, j in indices
        i != j && push!(candidate_sets[i], j)
      end
      return
    end

    proj = randn(rng, n_dim)
    proj ./= max(norm(proj), 1e-12)

    vals = zeros(Float64, length(indices))
    @inbounds for (k, idx) in enumerate(indices)
      v = 0.0
      for d in 1:n_dim
        v += fcs.events[idx, d] * proj[d]
      end
      vals[k] = v
    end
    med = median(vals)

    left = Int[]
    right = Int[]
    sizehint!(left, length(indices))
    sizehint!(right, length(indices))
    @inbounds for (k, idx) in enumerate(indices)
      if vals[k] <= med
        push!(left, idx)
      else
        push!(right, idx)
      end
    end

    if isempty(left) || isempty(right)
      half = length(indices) ÷ 2
      left = indices[1:half]
      right = indices[(half+1):end]
    end

    build_rptree!(left, depth + 1)
    build_rptree!(right, depth + 1)
  end

  for _ in 1:n_trees
    build_rptree!(collect(1:n_points), 0)
  end

  # Local-join / NN-descent 1-pass refinement step: connect 2-hop candidates for high recall
  refined_sets = [copy(candidate_sets[i]) for i in 1:n_points]
  @inbounds for i in 1:n_points
    for c in candidate_sets[i]
      for c_neighbor in candidate_sets[c]
        if c_neighbor != i
          push!(refined_sets[i], c_neighbor)
        end
      end
    end
  end

  knn_indices = [Vector{Int}(undef, K) for _ in 1:n_points]
  @inbounds for i in 1:n_points
    cands = collect(refined_sets[i])
    if length(cands) < K
      needed = K - length(cands)
      extra = filter(x -> x != i && !(x in refined_sets[i]), randperm(rng, n_points)[1:min(n_points, needed + 20)])
      append!(cands, extra)
    end
    dists = [sum(abs2, @view(fcs.events[i, :]) .- @view(fcs.events[j, :])) for j in cands]
    perm = partialsortperm(dists, 1:min(K, length(dists)))
    knn_indices[i] = cands[perm]
  end

  return knn_indices
end

"""
    xshift_cluster(fcs::FlowExperiment; K::Int=20, approximate::Bool=true, seed::Int=42)

Perform density-based XShift clustering (Samusik et al., 2016) on flow/mass cytometry data.
Supports both exact brute-force kNN (`approximate=false`) and multi-tree Random Projection Tree + local graph refinement approximate kNN (`approximate=true`).
"""
function xshift_cluster(fcs::FlowExperiment; K::Int=20, approximate::Bool=true, seed::Int=42)
  _ctx = active_provenance_context()

  n_points, n_dim = size(fcs.events)
  n_points > 1 || throw(ArgumentError("need at least 2 events for XShift clustering"))

  knn_indices = if approximate
    _approximate_knn(fcs, K; seed=seed)
  else
    _exact_knn(fcs.events, K)
  end

  # Compute density based on average distance to K nearest neighbors
  density = zeros(Float64, n_points)
  @inbounds for i in 1:n_points
    knn = knn_indices[i]
    avg_d = 0.0
    @inbounds for j in knn
      d = 0.0
      @inbounds for k in 1:n_dim
        d += (fcs.events[i, k] - fcs.events[j, k])^2
      end
      avg_d += sqrt(d)
    end
    avg_d /= length(knn)
    density[i] = 1.0 / (avg_d + 1e-6)
  end

  # Build parent tree (shift to neighbor with maximum density)
  parent = collect(1:n_points)
  @inbounds for i in 1:n_points
    neighbors = vcat(i, knn_indices[i])
    best_idx = neighbors[argmax(density[neighbors])]
    parent[i] = best_idx
  end

  # Find density peaks (roots)
  peaks = findall(i -> parent[i] == i, 1:n_points)

  # Assign clusters by following parent links
  cluster_assignment = zeros(Int, n_points)
  peak_map = Dict(p => idx for (idx, p) in enumerate(peaks))

  @inbounds for i in 1:n_points
    curr = i
    visited = Set{Int}()
    while parent[curr] != curr && !(curr in visited)
      push!(visited, curr)
      curr = parent[curr]
    end
    cluster_assignment[i] = get(peak_map, curr, 1)
  end

  prov_rec = provenance_record("XShiftResult", "FlowCytometry/xshift_cluster"; parameters=(K=K, approximate=approximate, n_clusters=length(peaks), n_events=n_points))
  return provenance_result!(_ctx, cluster_assignment, "xshift_cluster"; parents=[fcs.provenance.id], parameters=(K=K, approximate=approximate, n_clusters=length(peaks)))
end

# --- Transformation Functions ---

"""
    arcsinh_transform(fcs::FlowExperiment; cofactor::Float64=5.0, channels::Vector{String}=fcs.channels)

Apply inverse hyperbolic sine (arcsinh) transformation to specified channels.
Standard for CyTOF and mass cytometry data.
"""
function arcsinh_transform(fcs::FlowExperiment; cofactor::Float64=5.0, channels::Vector{String}=fcs.channels)
  _ctx = active_provenance_context()
  idxs = [findfirst(==(c), fcs.channels) for c in channels]
  any(idx -> idx === nothing, idxs) && throw(ArgumentError("one or more channels not found"))

  events = copy(fcs.events)
  @inbounds for idx in idxs
    @views events[:, idx] = asinh.(events[:, idx] ./ cofactor)
  end

  metadata = copy(fcs.metadata)
  metadata["transform"] = "arcsinh"
  metadata["cofactor"] = cofactor
  metadata["transformed_channels"] = channels

  prov_rec = provenance_record("FlowExperiment", "FlowCytometry/arcsinh_transform"; parameters=(cofactor=cofactor, channels=channels))
  result = FlowExperiment(events, fcs.channels, metadata, prov_rec)
  return provenance_result!(_ctx, result, "arcsinh_transform"; parents=[fcs.provenance.id], parameters=(cofactor=cofactor, channels=channels))
end

"""
    estimate_logicle_params(fcs::FlowExperiment, channel::String; q::Float64=0.05, M::Float64=4.5, A::Float64=0.0)
    estimate_logicle_params(fcs::FlowExperiment; q::Float64=0.05, M::Float64=4.5, A::Float64=0.0)

Estimate Logicle parameters `(T, W, M, A)` automatically from the dataset distribution for a given channel or across all channels.
"""
function estimate_logicle_params(fcs::FlowExperiment, channel::String; q::Float64=0.05, M::Float64=4.5, A::Float64=0.0)
  idx = findfirst(==(channel), fcs.channels)
  idx === nothing && throw(ArgumentError("channel '$channel' not found"))

  col = @view fcs.events[:, idx]
  T = max(maximum(col), 1000.0)

  r = quantile(col, q)
  W = if r < 0
    max(0.1, (M - log10(T / abs(r))) / 2.0)
  else
    0.5
  end
  W = min(W, M / 2.0)

  return (T=T, W=W, M=M, A=A)
end

function estimate_logicle_params(fcs::FlowExperiment; q::Float64=0.05, M::Float64=4.5, A::Float64=0.0)
  Dict(ch => estimate_logicle_params(fcs, ch; q=q, M=M, A=A) for ch in fcs.channels)
end

"""
    logicle_transform(fcs::FlowExperiment; channels::Vector{String}=fcs.channels, auto_estimate::Bool=false, T::Float64=262144.0, W::Float64=0.5, M::Float64=4.5, A::Float64=0.0)

Apply Logicle transformation (Parks et al., 2006) to specified channels.
Supports per-channel parameter auto-estimation (`auto_estimate=true`).
Standard for traditional flow cytometry data with negative values.
"""
function logicle_transform(fcs::FlowExperiment; channels::Vector{String}=fcs.channels, auto_estimate::Bool=false, T::Float64=262144.0, W::Float64=0.5, M::Float64=4.5, A::Float64=0.0)
  _ctx = active_provenance_context()
  idxs = [findfirst(==(c), fcs.channels) for c in channels]
  any(idx -> idx === nothing, idxs) && throw(ArgumentError("one or more channels not found"))

  events = copy(fcs.events)
  logicle_params_dict = Dict{String, Any}()

  @inbounds for (i, idx) in enumerate(idxs)
    ch = channels[i]
    p_T, p_W, p_M, p_A = if auto_estimate
      p = estimate_logicle_params(fcs, ch)
      (p.T, p.W, p.M, p.A)
    else
      (T, W, M, A)
    end
    logicle_params_dict[ch] = (T=p_T, W=p_W, M=p_M, A=p_A)

    x2 = p_W * p_M
    x1 = p_A + (x2 / (exp(p_M / p_W) - 1.0))
    b = p_T / (exp(p_M / p_W) - 1.0)

    for ev in 1:size(events, 1)
      val = events[ev, idx]
      events[ev, idx] = if val <= x1
        p_A + (val - x1) / b
      else
        p_A + p_W * log((val - x1) / b + 1.0)
      end
    end
  end

  metadata = copy(fcs.metadata)
  metadata["transform"] = "logicle"
  metadata["logicle_params"] = logicle_params_dict
  metadata["transformed_channels"] = channels

  prov_rec = provenance_record("FlowExperiment", "FlowCytometry/logicle_transform"; parameters=(auto_estimate=auto_estimate, T=T, W=W, M=M, A=A, channels=channels))
  result = FlowExperiment(events, fcs.channels, metadata, prov_rec)
  return provenance_result!(_ctx, result, "logicle_transform"; parents=[fcs.provenance.id], parameters=(auto_estimate=auto_estimate, T=T, W=W, M=M, A=A, channels=channels))
end

"""
    inverse_logicle_transform(fcs::FlowExperiment; channels::Vector{String}=fcs.channels, T::Float64=262144.0, W::Float64=0.5, M::Float64=4.5, A::Float64=0.0)

Apply inverse Logicle transformation to recover original scale.
"""
function inverse_logicle_transform(fcs::FlowExperiment; channels::Vector{String}=fcs.channels, T::Float64=262144.0, W::Float64=0.5, M::Float64=4.5, A::Float64=0.0)
  _ctx = active_provenance_context()
  idxs = [findfirst(==(c), fcs.channels) for c in channels]
  any(idx -> idx === nothing, idxs) && throw(ArgumentError("one or more channels not found"))

  x2 = W * M
  x1 = A + (x2 / (exp(M / W) - 1.0))
  b = T / (exp(M / W) - 1.0)

  function inverse_logicle(y::Float64)
    if y <= A
      return x1 + b * (y - A)
    else
      return x1 + b * (exp((y - A) / W) - 1.0)
    end
  end

  events = copy(fcs.events)
  @inbounds for idx in idxs
    @views events[:, idx] = inverse_logicle.(events[:, idx])
  end

  metadata = copy(fcs.metadata)
  metadata["transform"] = "inverse_logicle"

  prov_rec = provenance_record("FlowExperiment", "FlowCytometry/inverse_logicle_transform"; parameters=(T=T, W=W, M=M, A=A, channels=channels))
  result = FlowExperiment(events, fcs.channels, metadata, prov_rec)
  return provenance_result!(_ctx, result, "inverse_logicle_transform"; parents=[fcs.provenance.id], parameters=(T=T, W=W, M=M, A=A, channels=channels))
end

"""
    hlog_transform(fcs::FlowExperiment; channels::Vector{String}=fcs.channels, b::Float64=1.0)

Apply hyperbolic log (hlog) transformation to specified channels.
Smoothly handles negative values and zero.
"""
function hlog_transform(fcs::FlowExperiment; channels::Vector{String}=fcs.channels, b::Float64=1.0)
  _ctx = active_provenance_context()
  idxs = [findfirst(==(c), fcs.channels) for c in channels]
  any(idx -> idx === nothing, idxs) && throw(ArgumentError("one or more channels not found"))

  function hlog(x::Float64)
    return asinh(x / (2.0 * b)) + log(2.0 * b)
  end

  events = copy(fcs.events)
  @inbounds for idx in idxs
    @views events[:, idx] = hlog.(events[:, idx])
  end

  metadata = copy(fcs.metadata)
  metadata["transform"] = "hlog"
  metadata["hlog_b"] = b
  metadata["transformed_channels"] = channels

  prov_rec = provenance_record("FlowExperiment", "FlowCytometry/hlog_transform"; parameters=(b=b, channels=channels))
  result = FlowExperiment(events, fcs.channels, metadata, prov_rec)
  return provenance_result!(_ctx, result, "hlog_transform"; parents=[fcs.provenance.id], parameters=(b=b, channels=channels))
end

# --- Statistics and Utilities ---

"""
    fcs_stats(fcs::FlowExperiment; channels::Vector{String}=fcs.channels)

Compute summary statistics for each channel.
"""
function fcs_stats(fcs::FlowExperiment; channels::Vector{String}=fcs.channels)
  idxs = [findfirst(==(c), fcs.channels) for c in channels]
  any(idx -> idx === nothing, idxs) && throw(ArgumentError("one or more channels not found"))

  stats = Dict{String, Any}()
  for (name, idx) in zip(channels, idxs)
    col = @view fcs.events[:, idx]
    stats[name] = Dict(
      "n" => length(col),
      "mean" => mean(col),
      "std" => std(col),
      "median" => median(col),
      "min" => minimum(col),
      "max" => maximum(col),
      "q25" => quantile(col, 0.25),
      "q75" => quantile(col, 0.75),
      "cv" => std(col) / mean(col),
    )
  end
  return stats
end

"""
    fcs_correlation(fcs::FlowExperiment; channels::Vector{String}=fcs.channels)

Compute correlation matrix between channels.
"""
function fcs_correlation(fcs::FlowExperiment; channels::Vector{String}=fcs.channels)
  idxs = [findfirst(==(c), fcs.channels) for c in channels]
  any(idx -> idx === nothing, idxs) && throw(ArgumentError("one or more channels not found"))

  data = fcs.events[:, idxs]
  return cor(data)
end

end # module FlowCytometry
