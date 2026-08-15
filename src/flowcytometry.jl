# ==============================================================================
# flowcytometry.jl — Flow Cytometry and CyTOF analysis
# ==============================================================================

module FlowCytometry

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, active_provenance_context, provenance_result!
using LinearAlgebra
using Statistics
using Random

export FlowExperiment, GateResult, FlowSOMResult
export read_fcs, mock_flow_experiment, compensate_fcs, rectangle_gate, quadrant_gate, apply_gate, flowsom_cluster, xshift_cluster

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

struct FlowSOMResult <: AbstractAnalysisResult
  metaclusters::Vector{Int}
  codes::Matrix{Float64}
  provenance::ResultProvenance
end

function _parse_fcs_offset(header::String, range)
  txt = strip(header[range])
  isempty(txt) && return 0
  return parse(Int, txt)
end

function _parse_fcs_text(text::String)
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

function _fcs_int(bytes::Vector{UInt8}, byteord::String)
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
  events = zeros(Float64, n_events, n_channels)
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
    for i in 1:n_events, j in 1:n_channels
      chunk = raw[cursor:(cursor+width-1)]
      if startswith(byteord, "4,3,2,1")
        chunk = reverse(chunk)
      end
      events[i, j] = datatype == "F" ? Float64(reinterpret(Float32, chunk)[1]) : reinterpret(Float64, chunk)[1]
      cursor += width
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
  result = FlowExperiment(events, channels, Dict{String,Any}("source" => "explicit_mock", "seed" => seed), provenance_record("FlowExperiment", "FlowCytometry/mock_flow_experiment"; parameters=(n_events=n_events, channel_count=length(channels), seed=seed)))
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
    result = FlowExperiment(events, channels, meta, provenance_record("FlowExperiment", "FlowCytometry/read_fcs"; parameters=(path=path, version=version, n_events=n_events, n_channels=n_channels)))
    return provenance_result!(_ctx, result, "read_fcs"; parents=String[], parameters=(path=path, version=version, n_events=n_events, n_channels=n_channels))
  end
end

function compensate_fcs(fcs::FlowExperiment, spillover_matrix::Matrix{Float64})
  _ctx = active_provenance_context()
  n_ch = length(fcs.channels)
  if size(spillover_matrix) != (n_ch, n_ch)
    throw(ArgumentError("spillover_matrix dimensions must match the number of channels"))
  end

  events = fcs.events * inv(spillover_matrix)
  result = FlowExperiment(events, fcs.channels, copy(fcs.metadata), provenance_record("FlowExperiment", "flowcytometry"))
  return provenance_result!(_ctx, result, "compensate_fcs"; parents=[fcs.provenance.id])
end

function rectangle_gate(fcs::FlowExperiment, channel1::String, min1::Real, max1::Real, channel2::String, min2::Real, max2::Real)
  _ctx = active_provenance_context()
  idx1 = findfirst(==(channel1), fcs.channels)
  idx2 = findfirst(==(channel2), fcs.channels)
  if idx1 === nothing || idx2 === nothing
    throw(ArgumentError("channels not found"))
  end

  indices = Int[]
  for i in 1:size(fcs.events, 1)
    v1 = fcs.events[i, idx1]
    v2 = fcs.events[i, idx2]
    if min1 <= v1 <= max1 && min2 <= v2 <= max2
      push!(indices, i)
    end
  end

  result = GateResult(:rectangle, indices, provenance_record("GateResult", "flowcytometry"))
  return provenance_result!(_ctx, result, "rectangle_gate"; parents=[fcs.provenance.id])
end

function quadrant_gate(fcs::FlowExperiment, channel1::String, threshold1::Real, channel2::String, threshold2::Real)
  _ctx = active_provenance_context()
  idx1 = findfirst(==(channel1), fcs.channels)
  idx2 = findfirst(==(channel2), fcs.channels)
  if idx1 === nothing || idx2 === nothing
    throw(ArgumentError("channels not found"))
  end

  q1_idx = Int[] # LL
  q2_idx = Int[] # UL
  q3_idx = Int[] # UR
  q4_idx = Int[] # LR

  for i in 1:size(fcs.events, 1)
    v1 = fcs.events[i, idx1]
    v2 = fcs.events[i, idx2]
    if v1 < threshold1 && v2 < threshold2
      push!(q1_idx, i)
    elseif v1 < threshold1 && v2 >= threshold2
      push!(q2_idx, i)
    elseif v1 >= threshold1 && v2 >= threshold2
      push!(q3_idx, i)
    else
      push!(q4_idx, i)
    end
  end

  res = [
    GateResult(:Q1, q1_idx, provenance_record("GateResult", "flowcytometry")),
    GateResult(:Q2, q2_idx, provenance_record("GateResult", "flowcytometry")),
    GateResult(:Q3, q3_idx, provenance_record("GateResult", "flowcytometry")),
    GateResult(:Q4, q4_idx, provenance_record("GateResult", "flowcytometry"))
  ]

  return provenance_result!(_ctx, res, "quadrant_gate"; parents=[fcs.provenance.id])
end

function apply_gate(fcs::FlowExperiment, gate::GateResult)
  _ctx = active_provenance_context()
  sub_events = fcs.events[gate.indices, :]
  result = FlowExperiment(sub_events, fcs.channels, copy(fcs.metadata), provenance_record("FlowExperiment", "flowcytometry"))
  return provenance_result!(_ctx, result, "apply_gate"; parents=[fcs.provenance.id, gate.provenance.id])
end

function flowsom_cluster(fcs::FlowExperiment; n_metaclusters::Int=10)
  _ctx = active_provenance_context()

  grid_size = 10
  n_codes = grid_size * grid_size
  n_features = length(fcs.channels)

  events_norm = copy(fcs.events)
  means = mean(events_norm, dims=1)
  stds = std(events_norm, dims=1)
  for j in 1:n_features
    s = stds[j] > 0 ? stds[j] : 1.0
    events_norm[:, j] = (events_norm[:, j] .- means[j]) ./ s
  end

  rng = MersenneTwister(42)
  codes = randn(rng, n_codes, n_features)

  n_events = size(events_norm, 1)
  learning_rate = 0.5
  radius = 5.0

  for epoch in 1:10
    lr = learning_rate * exp(-epoch / 10.0)
    rad = radius * exp(-epoch / 10.0)
    rad_sq = rad^2

    shuffled = shuffle(rng, 1:n_events)
    for i in shuffled
      event = @view events_norm[i, :]

      best_idx = 1
      min_dist = Inf
      for c in 1:n_codes
        dist = sum((codes[c, :] .- event) .^ 2)
        if dist < min_dist
          min_dist = dist
          best_idx = c
        end
      end

      bmu_x = (best_idx - 1) % grid_size
      bmu_y = div(best_idx - 1, grid_size)

      for c in 1:n_codes
        cx = (c - 1) % grid_size
        cy = div(c - 1, grid_size)
        dist_sq = (cx - bmu_x)^2 + (cy - bmu_y)^2
        if dist_sq <= rad_sq
          influence = exp(-dist_sq / (2.0 * rad_sq))
          codes[c, :] .+= lr * influence * (event .- codes[c, :])
        end
      end
    end
  end

  dist_matrix = zeros(n_codes, n_codes)
  for i in 1:n_codes
    for j in i:n_codes
      dist_matrix[i, j] = norm(codes[i, :] .- codes[j, :])
      dist_matrix[j, i] = dist_matrix[i, j]
    end
  end

  cluster_labels = collect(1:n_codes)
  n_clusters = n_codes

  while n_clusters > n_metaclusters
    min_d = Inf
    best_i, best_j = 1, 2
    for i in 1:n_codes
      for j in (i+1):n_codes
        if cluster_labels[i] != cluster_labels[j]
          d = dist_matrix[i, j]
          if d < min_d
            min_d = d
            best_i = i
            best_j = j
          end
        end
      end
    end
    target = cluster_labels[best_j]
    dest = cluster_labels[best_i]
    for idx in 1:n_codes
      if cluster_labels[idx] == target
        cluster_labels[idx] = dest
      end
    end
    n_clusters = length(unique(cluster_labels))
  end

  unique_lbls = unique(cluster_labels)
  lbl_map = Dict(lbl => idx for (idx, lbl) in enumerate(unique_lbls))
  metacluster_map = [lbl_map[lbl] for lbl in cluster_labels]

  event_clusters = zeros(Int, n_events)
  for i in 1:n_events
    event = @view events_norm[i, :]
    best_idx = 1
    min_dist = Inf
    for c in 1:n_codes
      dist = sum((codes[c, :] .- event) .^ 2)
      if dist < min_dist
        min_dist = dist
        best_idx = c
      end
    end
    event_clusters[i] = metacluster_map[best_idx]
  end

  result = FlowSOMResult(event_clusters, codes, provenance_record("FlowSOMResult", "flowcytometry"))
  return provenance_result!(_ctx, result, "flowsom_cluster"; parents=[fcs.provenance.id])
end

function xshift_cluster(fcs::FlowExperiment; K::Int=20)
  _ctx = active_provenance_context()

  n_points, n_dim = size(fcs.events)

  dists = zeros(n_points, n_points)
  for i in 1:n_points
    for j in i:n_points
      d = norm(fcs.events[i, :] .- fcs.events[j, :])
      dists[i, j] = d
      dists[j, i] = d
    end
  end

  density = zeros(Float64, n_points)
  knn_indices = [Int[] for _ in 1:n_points]

  for i in 1:n_points
    sorted_indices = sortperm(dists[i, :])
    knn = sorted_indices[2:min(K+1, n_points)]
    knn_indices[i] = knn
    avg_d = mean(dists[i, knn])
    density[i] = 1.0 / (avg_d + 1e-6)
  end

  parent = collect(1:n_points)
  for i in 1:n_points
    neighbors = vcat(i, knn_indices[i])
    best_idx = neighbors[argmax(density[neighbors])]
    parent[i] = best_idx
  end

  peaks = findall(i -> parent[i] == i, 1:n_points)

  cluster_assignment = zeros(Int, n_points)
  peak_map = Dict(p => idx for (idx, p) in enumerate(peaks))

  for i in 1:n_points
    curr = i
    visited = Set{Int}()
    while parent[curr] != curr && !(curr in visited)
      push!(visited, curr)
      curr = parent[curr]
    end
    cluster_assignment[i] = get(peak_map, curr, 1)
  end

  return provenance_result!(_ctx, cluster_assignment, "xshift_cluster"; parents=[fcs.provenance.id])
end

end # module FlowCytometry
