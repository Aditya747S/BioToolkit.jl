# ==============================================================================
# deeplearning.jl — Lightweight deep-learning style genomics helpers
#
# References:
#   - Lopez et al. (2018) Nat Methods 15:1053-1058 (scVI)
#   - Theodoris et al. (2023) Nature 618:355-364 (Geneformer context)
#   - Bergen et al. (2020) Nat Biotechnol 38:1408-1414 (scVelo)
#   - Hao et al. (2021) Cell 184:3573-3587 (Seurat WNN)
#   - Lin et al. (2023) Science 379:eadc8743 (ESM-2 protein LM)
#   - van Dijk et al. (2018) Cell 174:716-729 (MAGIC)
# ==============================================================================

module DeepLearning

using DataFrames
using LinearAlgebra
using Statistics
using Random
using JSON
import UUIDs: uuid4
import ..BioToolkit
using ..BioToolkit: AminoAcidAlphabet, BioAlphabet, BioSequence, flux_available, maybe_to_device, maybe_to_host, resolve_backend, threaded_foreach
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_record, provenance_result!, register_provenance!, with_provenance, BackendConfig, ConvergenceReport, ModelFitDiagnostics

@inline function _register_dl_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
  return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export scvi_like_embedding, cellassign_like_mapping, geneformer_like_embedding, scgpt_like_embedding, scbert_like_embedding, attention_grn
export TrainingHistory, LatentModelResult, CellTypeModelResult, PerturbationModelResult
export fit_scvi_model, transform_scvi, predict_scvi_latent, fit_cellassign_model, predict_cell_types, fit_scgen_model, predict_perturbation_response, fit_sequence_transformer_embedding
export batch_corrected_latent, contrastive_cell_embedding
export flux_autoencoder_embedding, flux_mlp_classifier
export scgen_like_perturbation, graphsca_label_transfer

# Advanced single-cell & multi-modal modules
export cell_type_denoising
export sparse_autoencoder_features
export trajectory_neural_ode
export multimodal_wnn_embedding
export protein_sequence_embedding
export zero_shot_cell_annotation
export gene_regulatory_network_gnn
export self_supervised_pretraining
export cell_cycle_regression
export deep_factorization_embedding

# Interactive HTML Visualizers
export embedding_to_html, export_embedding_html
export grn_to_html, export_grn_html
export trajectory_to_html, export_trajectory_html
export cell_type_annotation_html, export_cell_type_annotation_html

struct TrainingHistory <: BioToolkit.AbstractAnalysisResult
  loss::Vector{Float64}
  epochs::Int
  converged::Bool
  seed::Int
  backend::Symbol
  provenance::BioToolkit.ResultProvenance
end

struct LatentModelResult <: BioToolkit.AbstractAnalysisResult
  latent::Matrix{Float64}
  loadings::Matrix{Float64}
  feature_means::Vector{Float64}
  diagnostics::ModelFitDiagnostics
  training_history::TrainingHistory
  model::Dict{Symbol,Any}
  provenance::BioToolkit.ResultProvenance
end

struct CellTypeModelResult <: BioToolkit.AbstractAnalysisResult
  labels::Vector{String}
  marker_sets::Dict{String,Vector{String}}
  gene_ids::Vector{String}
  weights::Matrix{Float64}
  diagnostics::ModelFitDiagnostics
  provenance::BioToolkit.ResultProvenance
end

struct PerturbationModelResult <: BioToolkit.AbstractAnalysisResult
  control_centroid::Vector{Float64}
  treated_centroid::Vector{Float64}
  delta::Vector{Float64}
  diagnostics::ModelFitDiagnostics
  provenance::BioToolkit.ResultProvenance
end

function _dl_diagnostics(method::Symbol, backend::Symbol, seed::Int, losses::Vector{Float64}, hyper::Dict{Symbol,Any}; warnings::Vector{String}=String[])
  final = isempty(losses) ? NaN : losses[end]
  first_loss = isempty(losses) ? final : losses[1]
  delta = isfinite(first_loss) && isfinite(final) ? first_loss - final : 0.0
  conv = ConvergenceReport(isempty(losses) ? true : final <= first_loss + sqrt(eps(Float64)), length(losses), final, delta, isempty(warnings) ? "completed" : join(warnings, "; "))
  cfg = BackendConfig(backend; device=backend === :cuda ? :cuda : :cpu, deterministic=true, seed=seed, parameters=copy(hyper))
  return ModelFitDiagnostics(method, cfg, conv; hyperparameters=copy(hyper), training_loss=losses, warnings=warnings)
end

function _pairwise_sqdistances(X::AbstractMatrix{<:Real})
  data = Matrix{Float64}(X)
  norms2 = sum(abs2, data, dims=2)
  D2 = norms2 .+ norms2' .- 2.0 .* (data * data')
  for i in axes(D2, 1)
    D2[i, i] = 0.0
  end
  return max.(D2, 0.0)
end

function fit_scvi_model(counts::AbstractMatrix{<:Real}; n_latent::Int=10, backend::Symbol=:auto, max_epochs::Int=1, seed::Int=1)
  _ctx = active_provenance_context()
  X = log1p.(Float64.(counts))
  cells_by_gene = permutedims(X)
  means = vec(mean(cells_by_gene, dims=1))
  centered = cells_by_gene .- permutedims(means)
  selected = resolve_backend(; backend=backend)
  work = maybe_to_device(centered; backend=selected)
  fac = svd(work; full=false)
  U = Matrix{Float64}(maybe_to_host(fac.U))
  S = Vector{Float64}(maybe_to_host(fac.S))
  V = Matrix{Float64}(maybe_to_host(fac.V))
  used = min(n_latent, size(U, 2))
  latent = U[:, 1:used] * Diagonal(S[1:used])
  loadings = V[:, 1:used]
  recon = latent * loadings' .+ permutedims(means)
  loss = mean(abs2, recon .- cells_by_gene)
  losses = fill(Float64(loss), max(max_epochs, 1))
  hyper = Dict{Symbol,Any}(:n_latent => used, :max_epochs => max_epochs, :method => :regularized_svd_gaussian_latent)
  diag = _dl_diagnostics(:fit_scvi_model, selected, seed, losses, hyper)
  hist = TrainingHistory(losses, length(losses), diag.convergence.converged, seed, selected, provenance_record("TrainingHistory", "DeepLearning/fit_scvi_model"; parameters=(epochs=length(losses), backend=selected, seed=seed)))
  model = Dict{Symbol,Any}(:method => :regularized_svd_gaussian_latent, :backend => selected, :n_latent => used)
  result = LatentModelResult(latent, loadings, means, diag, hist, model, provenance_record("LatentModelResult", "DeepLearning/fit_scvi_model"; parameters=(n_latent=used, backend=selected, seed=seed)))
  return _register_dl_result!(_ctx, result, "fit_scvi_model"; parents=provenance_parent_ids(counts), parameters=(n_latent=used, backend=selected, seed=seed))
end

function transform_scvi(model::LatentModelResult, counts::AbstractMatrix{<:Real})
  X = permutedims(log1p.(Float64.(counts)))
  size(X, 2) == length(model.feature_means) || throw(DimensionMismatch("counts gene dimension must match fitted scVI model"))
  return (X .- permutedims(model.feature_means)) * model.loadings
end

predict_scvi_latent(model::LatentModelResult, counts::AbstractMatrix{<:Real}) = transform_scvi(model, counts)

function fit_cellassign_model(expression::AbstractMatrix{<:Real}, gene_ids::AbstractVector{<:AbstractString}, marker_sets::AbstractDict; seed::Int=1)
  _ctx = active_provenance_context()
  X = Matrix{Float64}(expression)
  n_genes, n_cells = size(X)
  length(gene_ids) == n_genes || throw(DimensionMismatch("gene_ids must match expression rows"))
  genes = String.(gene_ids)
  labels = sort!(String.(collect(keys(marker_sets))))
  idx = Dict(g => i for (i, g) in enumerate(genes))
  weights = zeros(Float64, n_genes, length(labels))
  warnings = String[]
  for (j, label) in enumerate(labels)
    markers = [String(g) for g in marker_sets[label] if haskey(idx, String(g))]
    isempty(markers) && push!(warnings, "cell type $label had no matched markers")
    for g in markers
      weights[idx[g], j] = 1.0 / max(length(markers), 1)
    end
  end
  scores = (X' * weights)
  losses = [mean(abs2, scores .- mean(scores, dims=1))]
  hyper = Dict{Symbol,Any}(:label_count => length(labels), :gene_count => n_genes, :method => :marker_nb_linear_score)
  diag = _dl_diagnostics(:fit_cellassign_model, :cpu, seed, losses, hyper; warnings=warnings)
  result = CellTypeModelResult(labels, Dict(String(k) => String.(v) for (k, v) in marker_sets), genes, weights, diag, provenance_record("CellTypeModelResult", "DeepLearning/fit_cellassign_model"; notes=warnings, parameters=(label_count=length(labels), gene_count=n_genes, cell_count=n_cells)))
  return _register_dl_result!(_ctx, result, "fit_cellassign_model"; parents=provenance_parent_ids(expression), parameters=(label_count=length(labels), gene_count=n_genes, cell_count=n_cells))
end

function predict_cell_types(model::CellTypeModelResult, expression::AbstractMatrix{<:Real})
  X = Matrix{Float64}(expression)
  size(X, 1) == length(model.gene_ids) || throw(DimensionMismatch("expression row count must match fitted CellAssign model"))
  scores = X' * model.weights
  probs = _softmax_rows(scores)
  predicted = [model.labels[argmax(@view probs[i, :])] for i in axes(probs, 1)]
  return (predicted_label=predicted, probability=probs, label_order=model.labels, provenance=provenance_record("CellTypePrediction", "DeepLearning/predict_cell_types"; parameters=(cell_count=size(X, 2), label_count=length(model.labels))))
end

function fit_scgen_model(control_expression::AbstractMatrix{<:Real}, treated_expression::AbstractMatrix{<:Real}; seed::Int=1, backend::Symbol=:auto)
  _ctx = active_provenance_context()
  size(control_expression, 1) == size(treated_expression, 1) || throw(DimensionMismatch("control and treated expression must have matching genes"))
  control = vec(mean(log1p.(Float64.(control_expression)), dims=2))
  treated = vec(mean(log1p.(Float64.(treated_expression)), dims=2))
  delta = treated .- control
  loss = mean(abs2, (control .+ delta) .- treated)
  selected = resolve_backend(; backend=backend)
  hyper = Dict{Symbol,Any}(:gene_count => length(delta), :method => :centroid_delta)
  diag = _dl_diagnostics(:fit_scgen_model, selected, seed, [loss], hyper)
  result = PerturbationModelResult(control, treated, delta, diag, provenance_record("PerturbationModelResult", "DeepLearning/fit_scgen_model"; parameters=(gene_count=length(delta), backend=selected, seed=seed)))
  return _register_dl_result!(_ctx, result, "fit_scgen_model"; parents=provenance_parent_ids(control_expression, treated_expression), parameters=(gene_count=length(delta), backend=selected, seed=seed))
end

function predict_perturbation_response(model::PerturbationModelResult, query_control::AbstractMatrix{<:Real})
  X = log1p.(Float64.(query_control))
  size(X, 1) == length(model.delta) || throw(DimensionMismatch("query_control gene count must match fitted scGen model"))
  predicted = expm1.(X .+ model.delta)
  return (predicted_response=predicted, delta=model.delta, provenance=provenance_record("PerturbationPredictionResult", "DeepLearning/predict_perturbation_response"; parameters=(gene_count=size(X, 1), cell_count=size(X, 2))))
end

function fit_sequence_transformer_embedding(sequences::AbstractVector{<:AbstractString}; token_dim::Int=128, k::Int=3, seed::Int=1, backend::Symbol=:cpu, threaded::Bool=true)
  emb = scgpt_like_embedding(sequences; token_dim=token_dim, k=k, threaded=threaded)
  latent = Matrix{Float64}(emb)
  loadings = Matrix{Float64}(I, size(latent, 2), size(latent, 2))
  losses = [0.0]
  hyper = Dict{Symbol,Any}(:token_dim => token_dim, :k => k, :method => :hashed_context_embedding)
  diag = _dl_diagnostics(:fit_sequence_transformer_embedding, backend, seed, losses, hyper)
  hist = TrainingHistory(losses, 1, true, seed, backend, provenance_record("TrainingHistory", "DeepLearning/fit_sequence_transformer_embedding"; parameters=(backend=backend, seed=seed)))
  return LatentModelResult(latent, loadings, zeros(Float64, size(latent, 2)), diag, hist, Dict{Symbol,Any}(:method => :hashed_context_embedding), provenance_record("LatentModelResult", "DeepLearning/fit_sequence_transformer_embedding"; parameters=(sequence_count=length(sequences), token_dim=token_dim, k=k)))
end

fit_sequence_transformer_embedding(sequences::AbstractVector{<:BioSequence{A}}; kwargs...) where {A<:BioAlphabet} = fit_sequence_transformer_embedding(String.(sequences); kwargs...)

function _zscore_columns(x::AbstractMatrix{<:Real})
  data = Matrix{Float64}(x)
  for j in axes(data, 2)
    col = @view data[:, j]
    μ = mean(col)
    σ = std(col)
    if isfinite(σ) && σ > 0
      col .= (col .- μ) ./ σ
    else
      col .= 0.0
    end
  end
  return data
end

function _softmax_rows(x::AbstractMatrix{<:Real})
  out = copy(Matrix{Float64}(x))
  for i in axes(out, 1)
    row = @view out[i, :]
    row .-= maximum(row)
    row .= exp.(row)
    row ./= max(sum(row), eps(Float64))
  end
  return out
end

function scvi_like_embedding(counts::AbstractMatrix{<:Real}; n_latent::Int=10, backend::Symbol=:auto)
  X = log1p.(Float64.(counts))
  cell_by_gene = permutedims(X)
  cell_by_gene .-= mean(cell_by_gene, dims=1)
  selected = resolve_backend(; backend=backend)
  work = maybe_to_device(cell_by_gene; backend=selected)
  fac = svd(work; full=false)
  U = Matrix{Float64}(maybe_to_host(fac.U))
  S = Vector{Float64}(maybe_to_host(fac.S))
  V = Matrix{Float64}(maybe_to_host(fac.V))
  used = min(n_latent, size(U, 2))
  latent = U[:, 1:used] * Diagonal(S[1:used])
  result = (latent=latent, loadings=V[:, 1:used], backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/scvi_like_embedding"; parameters=(n_latent=used, backend=selected)))
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, result, "scvi_like_embedding"; parents=provenance_parent_ids(counts), parameters=(n_latent=used, backend=selected))
end

function cellassign_like_mapping(expression::AbstractMatrix{<:Real}, gene_ids::AbstractVector{<:AbstractString}, marker_sets::AbstractDict)
  X = Matrix{Float64}(expression)
  n_genes, n_cells = size(X)
  length(gene_ids) == n_genes || throw(DimensionMismatch("gene_ids must match expression rows"))
  index = Dict(String(g) => i for (i, g) in enumerate(gene_ids))

  labels = sort!(collect(String.(keys(marker_sets))))
  scores = zeros(Float64, n_cells, length(labels))
  for (j, label) in enumerate(labels)
    markers = String.(marker_sets[label])
    idx = [index[g] for g in markers if haskey(index, g)]
    isempty(idx) && continue
    scores[:, j] .= vec(mean(X[idx, :], dims=1))
  end

  predicted = [labels[argmax(@view(scores[i, :]))] for i in 1:n_cells]
  result = (predicted_label=predicted, score_matrix=scores, label_order=labels, provenance=provenance_record("DeepLearningResult", "DeepLearning/cellassign_like_mapping"; parameters=(gene_count=n_genes, cell_count=n_cells, label_count=length(labels))))
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, result, "cellassign_like_mapping"; parents=provenance_parent_ids(expression), parameters=(gene_count=n_genes, cell_count=n_cells, label_count=length(labels)))
end

function geneformer_like_embedding(sequences::AbstractVector{<:AbstractString}; k::Int=3, dim::Int=64, threaded::Bool=true)
  k >= 1 || throw(ArgumentError("k must be >= 1"))
  dim >= 4 || throw(ArgumentError("dim must be >= 4"))
  n_seqs = length(sequences)
  emb = zeros(Float64, n_seqs, dim)

  threaded_foreach(n_seqs, i -> begin
      s = uppercase(String(sequences[i]))
      cu = codeunits(s)
      n = length(cu)
      if n >= k
        for start in 1:(n-k+1)
          h = UInt(1469598103934665603)
          @inbounds for p in start:(start+k-1)
            h = (h ⊻ UInt(cu[p])) * UInt(1099511628211)
          end
          idx = Int(mod(h, UInt(dim))) + 1
          emb[i, idx] += 1.0
        end
        normv = norm(@view emb[i, :])
        normv > eps(Float64) && (emb[i, :] ./= normv)
      end
    end; threaded=threaded)
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, emb, "geneformer_like_embedding"; parents=provenance_parent_ids(sequences), parameters=(k=k, dim=dim, seq_count=n_seqs))
end

geneformer_like_embedding(sequences::AbstractVector{<:BioSequence{A}}; k::Int=3, dim::Int=64, threaded::Bool=true) where {A<:BioAlphabet} =
  geneformer_like_embedding(String.(sequences); k=k, dim=dim, threaded=threaded)

function scgpt_like_embedding(sequences::AbstractVector{<:AbstractString}; token_dim::Int=128, k::Int=3, context_window::Int=8, threaded::Bool=true)
  token_dim >= 8 || throw(ArgumentError("token_dim must be >= 8"))
  base = geneformer_like_embedding(sequences; k=k, dim=token_dim, threaded=threaded)
  n_seqs = length(sequences)
  pos = zeros(Float64, n_seqs, token_dim)
  win = max(Int(context_window), 1)

  threaded_foreach(n_seqs, i -> begin
      s = uppercase(String(sequences[i]))
      cu = codeunits(s)
      n = length(cu)
      n == 0 && return
      center = (n + 1) / 2.0
      for p in 1:n
        c = Int(cu[p])
        bucket = mod(c + p - 1, token_dim) + 1
        pos[i, bucket] += exp(-abs(p - center) / win)
      end
      nrm = norm(@view pos[i, :])
      nrm > eps(Float64) && (pos[i, :] ./= nrm)
    end; threaded=threaded)

  emb = 0.7 .* base .+ 0.3 .* pos
  for i in axes(emb, 1)
    nrm = norm(@view emb[i, :])
    nrm > eps(Float64) && (emb[i, :] ./= nrm)
  end
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, emb, "scgpt_like_embedding"; parents=provenance_parent_ids(sequences), parameters=(token_dim=token_dim, k=k, context_window=context_window, seq_count=n_seqs))
end

scgpt_like_embedding(sequences::AbstractVector{<:BioSequence{A}}; token_dim::Int=128, k::Int=3, context_window::Int=8, threaded::Bool=true) where {A<:BioAlphabet} =
  scgpt_like_embedding(String.(sequences); token_dim=token_dim, k=k, context_window=context_window, threaded=threaded)

function scbert_like_embedding(sequences::AbstractVector{<:AbstractString}; dim::Int=96, max_len::Int=512, threaded::Bool=true)
  dim >= 8 || throw(ArgumentError("dim must be >= 8"))
  max_len >= 1 || throw(ArgumentError("max_len must be >= 1"))
  n_seqs = length(sequences)
  emb = zeros(Float64, n_seqs, dim)

  threaded_foreach(n_seqs, i -> begin
      s = uppercase(String(sequences[i]))
      cu = codeunits(s)
      n = min(length(cu), max_len)
      n == 0 && return
      emb[i, 1] += 1.0
      for p in 1:n
        c = Int(cu[p])
        bucket = mod((c * 131 + p * 17), dim - 1) + 2
        emb[i, bucket] += 1.0
      end
      nrm = norm(@view emb[i, :])
      nrm > eps(Float64) && (emb[i, :] ./= nrm)
    end; threaded=threaded)

  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, emb, "scbert_like_embedding"; parents=provenance_parent_ids(sequences), parameters=(dim=dim, max_len=max_len, seq_count=n_seqs))
end

scbert_like_embedding(sequences::AbstractVector{<:BioSequence{A}}; dim::Int=96, max_len::Int=512, threaded::Bool=true) where {A<:BioAlphabet} =
  scbert_like_embedding(String.(sequences); dim=dim, max_len=max_len, threaded=threaded)

function attention_grn(expression::AbstractMatrix{<:Real}; gene_ids::AbstractVector{<:AbstractString}=String[], top_k::Int=100, temperature::Real=1.0)
  X = _zscore_columns(permutedims(expression))
  n_genes = size(expression, 1)

  score = (X' * X) ./ sqrt(size(X, 1)) ./ max(Float64(temperature), eps(Float64))
  score[diagind(score)] .= -Inf

  names_vec = length(gene_ids) == n_genes ? String.(gene_ids) : ["gene_$(i)" for i in 1:n_genes]

  edges = DataFrame(source=String[], target=String[], weight=Float64[])
  for i in 1:n_genes
    row = vec(@view score[i, :])
    order = sortperm(row, rev=true)
    for j in Iterators.take(order, min(top_k, n_genes - 1))
      isfinite(row[j]) || continue
      push!(edges, (names_vec[i], names_vec[j], row[j]))
    end
  end
  sort!(edges, :weight, rev=true)
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, edges, "attention_grn"; parents=provenance_parent_ids(expression), parameters=(top_k=top_k, temperature=Float64(temperature), n_genes=n_genes, edge_count=nrow(edges)))
end

function _latent_from_cells(cells_by_gene::AbstractMatrix{<:Real}, n_latent::Int; backend::Symbol=:auto)
  X = Matrix{Float64}(cells_by_gene)
  X .-= mean(X, dims=1)
  selected = resolve_backend(; backend=backend)
  work = maybe_to_device(X; backend=selected)
  fac = svd(work; full=false)
  U = Matrix{Float64}(maybe_to_host(fac.U))
  S = Vector{Float64}(maybe_to_host(fac.S))
  used = min(n_latent, size(U, 2))
  return U[:, 1:used] * Diagonal(S[1:used])
end

function batch_corrected_latent(counts::AbstractMatrix{<:Real}, batches; n_latent::Int=10, ridge::Real=1e-3, backend::Symbol=:auto)
  emb = scvi_like_embedding(counts; n_latent=n_latent, backend=backend)
  Z = Matrix{Float64}(emb.latent)
  n = size(Z, 1)
  length(batches) == n || throw(DimensionMismatch("batches length must equal number of cells"))

  labels = sort!(unique(String.(batches)))
  B = zeros(Float64, n, length(labels))
  idx = Dict(label => i for (i, label) in enumerate(labels))
  for i in 1:n
    B[i, idx[String(batches[i])]] = 1.0
  end

  coef = (B' * B + Float64(ridge) * I) \ (B' * Z)
  corrected = Z .- (B * coef)
  corrected .+= mean(Z, dims=1)
  result = (latent=corrected, raw_latent=Z, batch_coefficients=coef, batch_levels=labels, backend=emb.backend, provenance=provenance_record("DeepLearningResult", "DeepLearning/batch_corrected_latent"; parameters=(n_latent=n_latent, batch_count=length(labels), backend=emb.backend)))
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, result, "batch_corrected_latent"; parents=provenance_parent_ids(counts), parameters=(n_latent=n_latent, batch_count=length(labels), backend=emb.backend))
end

function contrastive_cell_embedding(counts::AbstractMatrix{<:Real}; n_latent::Int=16, dropout::Real=0.15, seed::Int=1, backend::Symbol=:auto)
  0.0 <= dropout < 1.0 || throw(ArgumentError("dropout must be in [0, 1)"))
  rng = MersenneTwister(seed)
  cells = permutedims(log1p.(Float64.(counts)))

  keep1 = rand(rng, size(cells)...) .>= Float64(dropout)
  keep2 = rand(rng, size(cells)...) .>= Float64(dropout)
  aug1 = cells .* keep1
  aug2 = cells .* keep2

  z1 = _latent_from_cells(aug1, n_latent; backend=backend)
  z2 = _latent_from_cells(aug2, n_latent; backend=backend)
  latent = (z1 .+ z2) ./ 2
  alignment = mean(sum(abs2, z1 .- z2, dims=2))
  selected = resolve_backend(; backend=backend)
  result = (latent=latent, view1=z1, view2=z2, alignment_loss=alignment, backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/contrastive_cell_embedding"; parameters=(n_latent=n_latent, dropout=Float64(dropout), seed=Int(seed), backend=selected)))
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, result, "contrastive_cell_embedding"; parents=provenance_parent_ids(counts), parameters=(n_latent=n_latent, dropout=Float64(dropout), seed=seed, backend=selected))
end

function _require_flux_module()
  flux_available() || throw(ArgumentError("Flux extension is not loaded. Install/load Flux to use flux_* APIs."))
  return BioToolkit._FLUX_MODULE[]
end

function flux_autoencoder_embedding(counts::AbstractMatrix{<:Real}; latent_dim::Int=16, hidden_dim::Int=64, epochs::Int=25, lr::Real=1e-3, backend::Symbol=:auto, seed::Int=1)
  F = _require_flux_module()
  Random.seed!(seed)

  X = Float32.(log1p.(Float64.(counts)))
  n_genes, n_cells = size(X)
  ldim = clamp(latent_dim, 2, max(2, min(n_genes, n_cells)))
  hdim = clamp(hidden_dim, ldim, max(ldim, n_genes))

  encoder = F.Chain(F.Dense(n_genes => hdim, F.relu), F.Dense(hdim => ldim))
  decoder = F.Chain(F.Dense(ldim => hdim, F.relu), F.Dense(hdim => n_genes))
  model = F.Chain(encoder, decoder)

  selected = resolve_backend(; backend=backend)
  xdev = X
  if selected == :gpu && isdefined(BioToolkit, :CUDA)
    model = F.gpu(model)
    xdev = F.gpu(X)
  end

  opt_state = F.setup(F.Adam(Float32(lr)), model)
  for _ in 1:max(1, epochs)
    loss, grads = F.withgradient(model) do m
      F.Losses.mse(m(xdev), xdev)
    end
    F.update!(opt_state, model, grads[1])
    isfinite(Float64(loss)) || break
  end

  latent = encoder(xdev)
  latent_host = permutedims(Array(F.cpu(latent)))
  result = (latent=latent_host, backend=selected, latent_dim=ldim, provenance=provenance_record("DeepLearningResult", "DeepLearning/flux_autoencoder_embedding"; parameters=(latent_dim=ldim, hidden_dim=Int(hidden_dim), epochs=Int(epochs), backend=selected)))
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, result, "flux_autoencoder_embedding"; parents=provenance_parent_ids(counts), parameters=(latent_dim=ldim, hidden_dim=Int(hidden_dim), epochs=Int(epochs), backend=selected))
end

function flux_mlp_classifier(features::AbstractMatrix{<:Real}, labels; hidden_dim::Int=64, epochs::Int=25, lr::Real=1e-3, backend::Symbol=:auto, seed::Int=1)
  F = _require_flux_module()
  Random.seed!(seed)

  X = Float32.(features)
  n_features, n_samples = size(X)
  y_labels = String.(labels)
  length(y_labels) == n_samples || throw(DimensionMismatch("labels must match number of samples (columns)"))
  classes = sort!(unique(y_labels))
  class_to_ix = Dict(c => i for (i, c) in enumerate(classes))
  y_ix = [class_to_ix[c] for c in y_labels]
  y_oh = F.onehotbatch(y_ix, 1:length(classes))

  hdim = clamp(hidden_dim, 4, max(4, n_features))
  model = F.Chain(F.Dense(n_features => hdim, F.relu), F.Dense(hdim => length(classes)))

  selected = resolve_backend(; backend=backend)
  xdev = X
  ydev = y_oh
  if selected == :gpu && isdefined(BioToolkit, :CUDA)
    model = F.gpu(model)
    xdev = F.gpu(X)
    ydev = F.gpu(y_oh)
  end

  opt_state = F.setup(F.Adam(Float32(lr)), model)
  for _ in 1:max(1, epochs)
    loss, grads = F.withgradient(model) do m
      logits = m(xdev)
      F.Losses.logitcrossentropy(logits, ydev)
    end
    F.update!(opt_state, model, grads[1])
    isfinite(Float64(loss)) || break
  end

  logits = model(xdev)
  probs = Array(F.cpu(F.softmax(logits; dims=1)))
  pred_ix = [argmax(@view probs[:, i]) for i in 1:size(probs, 2)]
  pred = [classes[i] for i in pred_ix]
  result = (classes=classes, probabilities=probs, predicted_label=pred, backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/flux_mlp_classifier"; parameters=(class_count=length(classes), hidden_dim=Int(hidden_dim), epochs=Int(epochs), backend=selected)))
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, result, "flux_mlp_classifier"; parents=provenance_parent_ids(features), parameters=(class_count=length(classes), hidden_dim=Int(hidden_dim), epochs=Int(epochs), backend=selected))
end

function scgen_like_perturbation(control_expression::AbstractMatrix{<:Real}, treated_expression::AbstractMatrix{<:Real}, query_control::AbstractMatrix{<:Real}; n_latent::Int=20, backend::Symbol=:auto)
  C = Matrix{Float64}(control_expression)
  T = Matrix{Float64}(treated_expression)
  Q = Matrix{Float64}(query_control)
  size(C, 1) == size(T, 1) == size(Q, 1) || throw(DimensionMismatch("all inputs must have the same gene dimension (rows)"))

  all_cells = hcat(C, T, Q)
  emb = scvi_like_embedding(all_cells; n_latent=n_latent, backend=backend)
  latent = emb.latent

  n_c = size(C, 2)
  n_t = size(T, 2)
  n_q = size(Q, 2)
  zc = latent[1:n_c, :]
  zt = latent[(n_c+1):(n_c+n_t), :]
  zq = latent[(n_c+n_t+1):(n_c+n_t+n_q), :]

  Δ = vec(mean(zt, dims=1) .- mean(zc, dims=1))
  zq_treated = zq .+ permutedims(Δ)

  X = log1p.(Float64.(all_cells))
  Xcell = permutedims(X)
  Xcell .-= mean(Xcell, dims=1)
  fac = svd(Xcell; full=false)
  used = min(size(zq_treated, 2), size(fac.V, 2))
  recon = zq_treated[:, 1:used] * fac.V[:, 1:used]'
  recon .+= mean(Xcell, dims=1)
  pred = exp.(recon) .- 1.0
  pred = max.(pred, 0.0)

  result = (predicted_treated=permutedims(pred), latent_shift=Δ, backend=emb.backend, provenance=provenance_record("DeepLearningResult", "DeepLearning/scgen_like_perturbation"; parameters=(n_latent=n_latent, backend=emb.backend)))
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, result, "scgen_like_perturbation"; parents=provenance_parent_ids(control_expression, treated_expression), parameters=(n_latent=n_latent, backend=emb.backend))
end

function graphsca_label_transfer(expression::AbstractMatrix{<:Real}, adjacency::AbstractMatrix{<:Real}, known_labels::AbstractVector{<:AbstractString}; alpha::Real=0.8, n_iter::Int=50)
  X = Matrix{Float64}(expression)
  A = Matrix{Float64}(adjacency)
  n_cells = size(X, 2)
  size(A, 1) == n_cells == size(A, 2) || throw(DimensionMismatch("adjacency must be square and match cell count"))
  length(known_labels) == n_cells || throw(DimensionMismatch("known_labels must match number of cells"))

  labels = String.(known_labels)
  classes = sort!(filter(!=(("unknown")), unique(labels)))
  isempty(classes) && throw(ArgumentError("at least one known label required"))
  class_to_ix = Dict(c => i for (i, c) in enumerate(classes))

  Y0 = zeros(Float64, n_cells, length(classes))
  for i in 1:n_cells
    lab = labels[i]
    haskey(class_to_ix, lab) || continue
    Y0[i, class_to_ix[lab]] = 1.0
  end

  deg = vec(sum(abs.(A), dims=2))
  W = copy(A)
  for i in 1:n_cells
    d = deg[i] > 0 ? deg[i] : 1.0
    W[i, :] ./= d
  end

  F = copy(Y0)
  α = clamp(Float64(alpha), 0.0, 1.0)
  for _ in 1:max(1, n_iter)
    F = α * (W * F) + (1 - α) * Y0
    row_sum = vec(sum(F, dims=2))
    for i in 1:n_cells
      s = row_sum[i]
      s > 0 && (F[i, :] ./= s)
    end
  end

  pred_ix = [argmax(@view F[i, :]) for i in 1:n_cells]
  pred = [classes[ix] for ix in pred_ix]
  confidence = [maximum(@view F[i, :]) for i in 1:n_cells]
  result = (predicted_label=pred, confidence=confidence, class_order=classes, score_matrix=F, provenance=provenance_record("DeepLearningResult", "DeepLearning/graphsca_label_transfer"; parameters=(class_count=length(classes), n_iter=Int(n_iter), alpha=Float64(alpha))))
  _ctx = active_provenance_context()
  return _register_dl_result!(_ctx, result, "graphsca_label_transfer"; parents=provenance_parent_ids(expression), parameters=(class_count=length(classes), n_iter=Int(n_iter), alpha=Float64(alpha)))
end

# ---------------------------------------------------------------------------
# Cell-type Denoising (MAGIC-like diffusion on kNN graph)
# ---------------------------------------------------------------------------

function cell_type_denoising(
  counts::AbstractMatrix{<:Real};
  k::Int=15,
  t::Int=3,
  n_pcs::Int=20,
  backend::Symbol=:auto)
  k >= 1 || throw(ArgumentError("k must be >= 1"))
  t >= 1 || throw(ArgumentError("t must be >= 1"))

  X = log1p.(Float64.(counts))   # genes × cells
  cells = permutedims(X)          # cells × genes
  cells .-= mean(cells, dims=1)

  selected = resolve_backend(; backend=backend)
  work = maybe_to_device(cells; backend=selected)
  fac = svd(work; full=false)
  U = Matrix{Float64}(maybe_to_host(fac.U))
  S = Vector{Float64}(maybe_to_host(fac.S))
  used = min(n_pcs, size(U, 2))
  pca = U[:, 1:used] * Diagonal(S[1:used])   # cells × PCs

  n_cells = size(pca, 1)
  kk = min(k, n_cells - 1)

  # Vectorized kNN adjacency (Gaussian kernel)
  D2 = _pairwise_sqdistances(pca)
  W = zeros(Float64, n_cells, n_cells)
  for i in 1:n_cells
    row = @view D2[i, :]
    order = sortperm(row)
    kth_dist = row[order[kk+1]]
    denom = max(2.0 * kth_dist, eps(Float64))
    for idx in 2:(kk+1)
      j = order[idx]
      W[i, j] = exp(-row[j] / denom)
    end
  end
  # Symmetrise + row-normalise → Markov matrix
  W = 0.5 .* (W .+ W')
  row_sums = vec(sum(W, dims=2))
  for i in 1:n_cells
    row_sums[i] > 0 && (W[i, :] ./= row_sums[i])
  end

  # Diffuse: M^t * X
  Mt = copy(W)
  for _ in 1:(t-1)
    Mt = Mt * W
  end

  denoised_cells = Mt * cells                       # cells × genes (smoothed)
  denoised = permutedims(denoised_cells)            # genes × cells
  result = (denoised=denoised, diffusion_operator=Mt, provenance=provenance_record("DeepLearningResult", "DeepLearning/cell_type_denoising"))
  _ctx = active_provenance_context()

  return _register_dl_result!(_ctx, result, "cell_type_denoising"; parents=provenance_parent_ids(counts), parameters=(k=k, t=t, n_pcs=n_pcs, backend=string(selected)))
end

# ---------------------------------------------------------------------------
# Sparse Autoencoder Features (K-sparse)
# ---------------------------------------------------------------------------

function sparse_autoencoder_features(
  counts::AbstractMatrix{<:Real};
  n_features::Int=64,
  sparsity_k::Int=8,
  n_iter::Int=300,
  lr::Real=5e-3,
  seed::Int=1,
  backend::Symbol=:auto)
  n_features >= 1 || throw(ArgumentError("n_features must be >= 1"))
  sparsity_k = clamp(sparsity_k, 1, n_features)

  X = log1p.(Float64.(counts))   # genes × cells
  n_genes, n_cells = size(X)
  X_host = permutedims(X)         # cells × genes

  selected = resolve_backend(; backend=backend)
  rng = MersenneTwister(seed)
  D = randn(rng, n_genes, n_features) ./ sqrt(n_genes)

  function normalise_dict!(D_mat)
    for j in 1:size(D_mat, 2)
      col = @view D_mat[:, j]
      nrm = norm(col)
      nrm > eps(Float64) && (col ./= nrm)
    end
    return D_mat
  end
  normalise_dict!(D)

  lr_f = Float64(lr)
  best_loss = Inf

  for iter in 1:n_iter
    A = X_host * D
    for i in 1:n_cells
      row = @view A[i, :]
      threshold = partialsort(abs.(row), sparsity_k, rev=true)
      row[abs.(row) .< threshold] .= 0.0
    end

    Xhat = A * D'
    residual = Xhat .- X_host
    loss = mean(abs2, residual)
    best_loss = min(best_loss, loss)

    grad = A' * residual ./ n_cells
    D .-= lr_f .* grad'
    normalise_dict!(D)

    if mod(iter, 100) == 0
      lr_f *= 0.5
    end
  end

  A_final = X_host * D
  for i in 1:n_cells
    row = @view A_final[i, :]
    threshold = partialsort(abs.(row), sparsity_k, rev=true)
    row[abs.(row) .< threshold] .= 0.0
  end

  recon = A_final * D'
  final_loss = mean(abs2, recon .- X_host)

  result = (features=D, feature_activations=A_final, reconstruction_loss=final_loss, backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/sparse_autoencoder_features"; parameters=(backend=selected,)))
  _ctx = active_provenance_context()

  return _register_dl_result!(_ctx, result, "sparse_autoencoder_features"; parents=provenance_parent_ids(counts), parameters=(n_features=n_features, sparsity_k=sparsity_k, n_iter=n_iter, backend=selected))
end

# ---------------------------------------------------------------------------
# Trajectory Neural ODE (Euler approximation, scVelo-NN style)
# ---------------------------------------------------------------------------

function trajectory_neural_ode(
  spliced_counts::AbstractMatrix{<:Real},
  unspliced_counts::AbstractMatrix{<:Real};
  n_latent::Int=10,
  n_steps::Int=20,
  dt::Real=0.1,
  backend::Symbol=:auto)
  size(spliced_counts) == size(unspliced_counts) ||
    throw(DimensionMismatch("spliced and unspliced count matrices must have the same shape"))

  emb_s = scvi_like_embedding(spliced_counts; n_latent=n_latent, backend=backend)
  emb_u = scvi_like_embedding(unspliced_counts; n_latent=n_latent, backend=backend)

  Zs = Matrix{Float64}(emb_s.latent)   # cells × latent
  Zu = Matrix{Float64}(emb_u.latent)

  velocity = Zu .- Zs   # cells × latent

  n_cells = size(Zs, 1)
  dim_lat = size(Zs, 2)

  trajectory = zeros(Float64, n_cells, dim_lat, n_steps + 1)
  trajectory[:, :, 1] = Zs

  for step in 1:n_steps
    z_current = trajectory[:, :, step]
    z_next = z_current .+ Float64(dt) .* velocity
    trajectory[:, :, step+1] = z_next
  end

  final_pos = trajectory[:, :, end]
  pseudotime = [norm(final_pos[i, :] .- Zs[i, :]) for i in 1:n_cells]
  pseudotime ./= max(maximum(pseudotime), eps(Float64))

  result = (pseudotime=pseudotime, latent_velocity=velocity, trajectory=trajectory, backend=emb_s.backend, provenance=provenance_record("DeepLearningResult", "DeepLearning/trajectory_neural_ode"; parameters=(backend=emb_s.backend,)))
  _ctx = active_provenance_context()

  return _register_dl_result!(_ctx, result, "trajectory_neural_ode"; parents=provenance_parent_ids(spliced_counts, unspliced_counts), parameters=(n_latent=n_latent, n_steps=n_steps, dt=Float64(dt), backend=emb_s.backend))
end

# ---------------------------------------------------------------------------
# Weighted Nearest Neighbour Multi-modal Embedding (Seurat WNN style)
# ---------------------------------------------------------------------------

function multimodal_wnn_embedding(
  modalities::AbstractVector;
  n_latent::Int=20,
  n_neighbors::Int=20,
  backend::Symbol=:auto,
  seed::Int=1)
  isempty(modalities) && throw(ArgumentError("provide at least one modality"))
  selected = resolve_backend(; backend=backend)

  n_cells = size(modalities[1], 2)
  for (k, m) in enumerate(modalities)
    size(m, 2) == n_cells || throw(DimensionMismatch("all modalities must have the same number of cells"))
  end

  M = length(modalities)
  latents = Matrix{Float64}[]
  for mod_mat in modalities
    X = log1p.(Float64.(mod_mat))
    cells = permutedims(X)
    cells .-= mean(cells, dims=1)
    work = maybe_to_device(cells; backend=selected)
    fac = svd(work; full=false)
    U = Matrix{Float64}(maybe_to_host(fac.U))
    S_vals = Vector{Float64}(maybe_to_host(fac.S))
    used = min(n_latent, size(U, 2))
    push!(latents, U[:, 1:used] * Diagonal(S_vals[1:used]))
  end

  kk = min(n_neighbors, n_cells - 1)
  precision = zeros(Float64, n_cells, M)

  for (m_idx, Z) in enumerate(latents)
    D2 = _pairwise_sqdistances(Z)
    for i in 1:n_cells
      row = @view D2[i, :]
      order = sortperm(row)
      knn_dists = row[order[2:(kk+1)]]
      precision[i, m_idx] = mean(1.0 ./ max.(knn_dists, eps(Float64)))
    end
  end

  log_prec = log.(max.(precision, eps(Float64)))
  weights = copy(log_prec)
  for i in 1:n_cells
    row = @view weights[i, :]
    row .-= maximum(row)
    row .= exp.(row)
    row ./= max(sum(row), eps(Float64))
  end

  combined_dim = minimum(size(l, 2) for l in latents)
  combined = zeros(Float64, n_cells, combined_dim)
  for (m_idx, Z) in enumerate(latents)
    used_dim = min(combined_dim, size(Z, 2))
    for i in 1:n_cells
      combined[i, 1:used_dim] .+= weights[i, m_idx] .* Z[i, 1:used_dim]
    end
  end

  nn_graph = zeros(Float64, n_cells, n_cells)
  D2_comb = _pairwise_sqdistances(combined)
  for i in 1:n_cells
    row = @view D2_comb[i, :]
    order = sortperm(row)
    for idx in 2:(kk+1)
      nn_graph[i, order[idx]] = 1.0
    end
  end

  result = (latent=combined, modality_weights=weights, nn_graph=nn_graph, backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/multimodal_wnn_embedding"; parameters=(backend=selected,)))
  _ctx = active_provenance_context()

  return _register_dl_result!(_ctx, result, "multimodal_wnn_embedding"; parents=String[], parameters=(n_modalities=length(modalities), n_latent=n_latent, n_neighbors=n_neighbors, backend=selected))
end

# ---------------------------------------------------------------------------
# Protein Sequence Embedding (ESM-like)
# ---------------------------------------------------------------------------

function protein_sequence_embedding(
  sequences::AbstractVector{<:AbstractString};
  dim::Int=128,
  k::Int=3,
  include_properties::Bool=true,
  threaded::Bool=true)
  dim >= 8 || throw(ArgumentError("dim must be >= 8"))

  kd = Dict('I'=>4.5, 'V'=>4.2, 'L'=>3.8, 'F'=>2.8, 'C'=>2.5, 'M'=>1.9, 'A'=>1.8,
    'G'=>-0.4, 'T'=>-0.7, 'S'=>-0.8, 'W'=>-0.9, 'Y'=>-1.3, 'P'=>-1.6,
    'H'=>-3.2, 'E'=>-3.5, 'Q'=>-3.5, 'D'=>-3.5, 'N'=>-3.5, 'K'=>-3.9, 'R'=>-4.5)
  charge_aa = Dict('R'=>1.0, 'K'=>1.0, 'H'=>0.1, 'D'=>-1.0, 'E'=>-1.0)
  helix_prop = Dict('A'=>1.42, 'L'=>1.21, 'M'=>1.45, 'E'=>1.51, 'K'=>1.16, 'R'=>0.98)
  sheet_prop = Dict('V'=>1.70, 'I'=>1.60, 'T'=>1.19, 'Y'=>1.47, 'W'=>1.35, 'F'=>1.38)

  base_emb = geneformer_like_embedding(sequences; k=k, dim=dim, threaded=threaded)

  if !include_properties
    return base_emb
  end

  prop_buckets = 8
  prop_emb = zeros(Float64, length(sequences), prop_buckets)

  threaded_foreach(length(sequences), idx -> begin
      seq = uppercase(String(sequences[idx]))
      aa = collect(seq)
      L = length(aa)
      L == 0 && return

      hydro = mean(get(kd, c, 0.0) for c in aa)
      charge = sum(get(charge_aa, c, 0.0) for c in aa) / max(L, 1)
      helix = mean(get(helix_prop, c, 1.0) for c in aa)
      sheet = mean(get(sheet_prop, c, 1.0) for c in aa)
      mw_proxy = L * 110.0 / 1000.0   # rough MW in kDa

      prop_emb[idx, 1] = hydro / 5.0
      prop_emb[idx, 2] = clamp(charge, -1.0, 1.0)
      prop_emb[idx, 3] = (helix - 1.0) / 0.5
      prop_emb[idx, 4] = (sheet - 1.0) / 0.5
      prop_emb[idx, 5] = log(mw_proxy + 1) / 5.0
      prop_emb[idx, 6] = count(c -> c == 'C', aa) / max(L, 1)
      prop_emb[idx, 7] = count(c -> c == 'P', aa) / max(L, 1)
      prop_emb[idx, 8] = count(c -> c in ('W', 'Y', 'F'), aa) / max(L, 1)
    end; threaded=threaded)

  prop_emb_full = hcat(base_emb[:, 1:(dim-prop_buckets)], prop_emb)

  for i in axes(prop_emb_full, 1)
    nrm = norm(@view prop_emb_full[i, :])
    nrm > eps(Float64) && (prop_emb_full[i, :] ./= nrm)
  end

  return prop_emb_full
end

protein_sequence_embedding(
  sequences::AbstractVector{<:BioSequence{AminoAcidAlphabet}};
  dim::Int=128,
  k::Int=3,
  include_properties::Bool=true,
  threaded::Bool=true) = protein_sequence_embedding(String.(sequences); dim=dim, k=k, include_properties=include_properties, threaded=threaded)

# ---------------------------------------------------------------------------
# Zero-Shot Cell Annotation using Archetype Vectors
# ---------------------------------------------------------------------------

function zero_shot_cell_annotation(
  counts::AbstractMatrix{<:Real},
  gene_ids::AbstractVector{<:AbstractString},
  archetype_profiles::AbstractDict;
  top_k::Int=1,
  log_transform::Bool=true,
  normalise::Bool=true)
  n_genes, n_cells = size(counts)
  length(gene_ids) == n_genes || throw(DimensionMismatch("gene_ids must match expression rows"))

  X = log_transform ? log1p.(Float64.(counts)) : Float64.(counts)
  if normalise
    col_totals = vec(sum(X, dims=1))
    scale = 1e4
    X = X ./ max.(col_totals', eps(Float64)) .* scale
  end

  gene_index = Dict(String(g) => i for (i, g) in enumerate(gene_ids))
  labels = sort!(collect(String.(keys(archetype_profiles))))
  n_types = length(labels)

  score_matrix = zeros(Float64, n_cells, n_types)
  for (j, label) in enumerate(labels)
    markers = String.(archetype_profiles[label])
    idx = [gene_index[g] for g in markers if haskey(gene_index, g)]
    isempty(idx) && continue
    score_matrix[:, j] = vec(mean(X[idx, :], dims=1))
  end

  score_norm = _softmax_rows(score_matrix)

  pred_labels = String[]
  pred_scores = Float64[]
  for i in 1:n_cells
    row = score_norm[i, :]
    best_idx = argmax(row)
    push!(pred_labels, labels[best_idx])
    push!(pred_scores, row[best_idx])
  end

  cell_df = DataFrame(
    cell_index=1:n_cells,
    predicted_label=pred_labels,
    score=pred_scores)
  return (cells=with_provenance(cell_df, "DeepLearningTable", "DeepLearning/zero_shot_cell_annotation"; parameters=(label_count=length(labels),)), score_matrix=score_norm, label_order=labels, provenance=provenance_record("DeepLearningResult", "DeepLearning/zero_shot_cell_annotation"; parameters=(label_count=length(labels),)))
end

# ---------------------------------------------------------------------------
# Gene Regulatory Network via GNN propagation
# ---------------------------------------------------------------------------

function gene_regulatory_network_gnn(
  expression::AbstractMatrix{<:Real};
  gene_ids::AbstractVector{<:AbstractString}=String[],
  k_neighbors::Int=10,
  n_propagation::Int=3,
  top_k_edges::Int=200,
  temperature::Real=1.0,
  n_pcs::Int=20,
  backend::Symbol=:auto)
  k_neighbors >= 1 || throw(ArgumentError("k_neighbors must be >= 1"))
  n_propagation >= 1 || throw(ArgumentError("n_propagation must be >= 1"))

  n_genes, n_cells = size(expression)
  selected = resolve_backend(; backend=backend)

  X = log1p.(Float64.(expression))   # genes × cells
  cells = permutedims(X)              # cells × genes

  cells_centered = cells .- mean(cells, dims=1)
  work = maybe_to_device(cells_centered; backend=selected)
  fac = svd(work; full=false)
  U = Matrix{Float64}(maybe_to_host(fac.U))
  S_vals = Vector{Float64}(maybe_to_host(fac.S))
  pca = U[:, 1:min(n_pcs, size(U, 2))] * Diagonal(S_vals[1:min(n_pcs, size(U, 2))])

  kk = min(k_neighbors, n_cells - 1)
  D2 = _pairwise_sqdistances(pca)
  W = zeros(Float64, n_cells, n_cells)
  for i in 1:n_cells
    row = @view D2[i, :]
    order = sortperm(row)
    for idx in 2:(kk+1)
      W[i, order[idx]] = 1.0
    end
  end
  row_s = vec(sum(W, dims=2))
  for i in 1:n_cells
    row_s[i] > 0 && (W[i, :] ./= row_s[i])
  end

  Xz = _zscore_columns(X)   # genes × cells
  raw_score = (Xz * Xz') ./ sqrt(n_cells) ./ max(Float64(temperature), eps(Float64))
  raw_score[diagind(raw_score)] .= -Inf

  X_prop = copy(X)
  for _ in 1:n_propagation
    X_prop = X_prop * W'   # genes × cells, smoothed
  end
  Xpz = _zscore_columns(X_prop)
  prop_score = (Xpz * Xpz') ./ sqrt(n_cells) ./ max(Float64(temperature), eps(Float64))
  prop_score[diagind(prop_score)] .= -Inf

  raw_vec = Float64[]
  prop_vec = Float64[]
  src_idx = Int[]
  tgt_idx = Int[]

  for i in 1:n_genes
    row = @view prop_score[i, :]
    order = sortperm(vec(row), rev=true)
    for j in order[1:min(top_k_edges÷n_genes+1, n_genes-1)]
      isfinite(row[j]) || continue
      push!(src_idx, i)
      push!(tgt_idx, j)
      push!(raw_vec, raw_score[i, j])
      push!(prop_vec, row[j])
    end
  end

  gene_names = length(gene_ids) == n_genes ? String.(gene_ids) : ["gene_$(i)" for i in 1:n_genes]

  edges = DataFrame(
    source_gene=[gene_names[i] for i in src_idx],
    target_gene=[gene_names[j] for j in tgt_idx],
    raw_weight=raw_vec,
    propagated_weight=prop_vec)
  sort!(edges, :propagated_weight, rev=true)
  edges = first(edges, min(top_k_edges, nrow(edges)))
  edges[!, :rank] = 1:nrow(edges)
  return (edges=with_provenance(edges, "DeepLearningTable", "DeepLearning/gene_regulatory_network_gnn"), backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/gene_regulatory_network_gnn"; parameters=(backend=selected,)))
end

# ---------------------------------------------------------------------------
# Self-Supervised Pre-training (Masked Gene Modelling)
# ---------------------------------------------------------------------------

function self_supervised_pretraining(
  counts::AbstractMatrix{<:Real};
  mask_fraction::Real=0.15,
  n_latent::Int=32,
  n_iter::Int=200,
  lr::Real=1e-3,
  seed::Int=1,
  backend::Symbol=:auto)
  0 < mask_fraction < 1 || throw(ArgumentError("mask_fraction must be in (0, 1)"))
  n_latent >= 2 || throw(ArgumentError("n_latent must be >= 2"))

  X = log1p.(Float64.(counts))   # genes × cells
  n_genes, n_cells = size(X)
  cells = permutedims(X)          # cells × genes

  rng = MersenneTwister(seed)
  selected = resolve_backend(; backend=backend)

  lr_f = Float64(lr)
  W_enc = randn(rng, n_genes, n_latent) ./ sqrt(n_genes)
  W_dec = randn(rng, n_latent, n_genes) ./ sqrt(n_latent)

  loss_history = Float64[]

  for iter in 1:n_iter
    mask = rand(rng, n_cells, n_genes) .< Float64(mask_fraction)
    X_masked = copy(cells)
    X_masked[mask] .= 0.0

    Z = X_masked * W_enc                        # cells × latent
    Xhat = Z * W_dec                             # cells × genes

    residual = (Xhat .- cells) .* mask
    loss = sum(abs2, residual) / max(sum(mask), 1)
    push!(loss_history, loss)

    dL_dXhat = 2 .* residual ./ max(sum(mask), 1)
    dL_dWdec = Z' * dL_dXhat                    # latent × genes
    dL_dZ = dL_dXhat * W_dec'                # cells × latent
    dL_dWenc = X_masked' * dL_dZ               # genes × latent

    W_dec .-= lr_f .* dL_dWdec
    W_enc .-= lr_f .* dL_dWenc

    if mod(iter, 50) == 0
      lr_f *= 0.9
    end
  end

  latent = cells * W_enc   # cells × latent
  return (latent=latent, reconstruction_loss_history=loss_history, backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/self_supervised_pretraining"; parameters=(backend=selected,)))
end

# ---------------------------------------------------------------------------
# Cell Cycle Regression
# ---------------------------------------------------------------------------

function cell_cycle_regression(
  counts::AbstractMatrix{<:Real},
  gene_ids::AbstractVector{<:AbstractString};
  s_genes::AbstractVector{<:AbstractString}=String[],
  g2m_genes::AbstractVector{<:AbstractString}=String[],
  ridge::Real=1e-3)
  n_genes, n_cells = size(counts)
  length(gene_ids) == n_genes || throw(DimensionMismatch("gene_ids must match expression rows"))

  gene_index = Dict(uppercase(String(g)) => i for (i, g) in enumerate(gene_ids))
  X = log1p.(Float64.(counts))  # genes × cells

  function score_module(gene_list)
    idx = [gene_index[uppercase(String(g))] for g in gene_list if haskey(gene_index, uppercase(String(g)))]
    isempty(idx) && return zeros(Float64, n_cells)
    vec(mean(X[idx, :], dims=1))
  end

  s_score = score_module(s_genes)
  g2m_score = score_module(g2m_genes)

  phase = [s > g2m && s > 0.1 ? "S" :
           g2m > s && g2m > 0.1 ? "G2M" : "G1"
           for (s, g2m) in zip(s_score, g2m_score)]

  X_t = permutedims(X)   # cells × genes
  B = hcat(ones(Float64, n_cells), s_score, g2m_score)  # cells × 3
  coef = (B' * B + Float64(ridge) * I) \ (B' * X_t)     # 3 × genes
  residuals = X_t .- B * coef                            # cells × genes
  corrected = permutedims(residuals)                     # genes × cells

  return (corrected_counts=corrected, s_score=s_score, g2m_score=g2m_score, phase=phase, provenance=provenance_record("DeepLearningResult", "DeepLearning/cell_cycle_regression"))
end

# ---------------------------------------------------------------------------
# Deep Factorization Embedding (MOFA+/JIVE-style multi-factor)
# ---------------------------------------------------------------------------

function deep_factorization_embedding(
  data_matrices::AbstractVector;
  n_factors::Int=10,
  n_iter::Int=300,
  lr::Real=1e-2,
  seed::Int=1,
  backend::Symbol=:auto)
  isempty(data_matrices) && throw(ArgumentError("provide at least one data matrix"))
  n_cells = size(data_matrices[1], 2)
  for (k, m) in enumerate(data_matrices)
    size(m, 2) == n_cells || throw(DimensionMismatch("all matrices must share the same number of cells"))
  end

  rng = MersenneTwister(seed)
  selected = resolve_backend(; backend=backend)
  K = min(n_factors, n_cells)

  views = [log1p.(Float64.(m)) for m in data_matrices]
  Z = randn(rng, n_cells, K) ./ sqrt(K)
  loadings = [randn(rng, size(v, 1), K) ./ sqrt(K) for v in views]

  lr_f = Float64(lr)
  for iter in 1:n_iter
    for (v_idx, v) in enumerate(views)
      Y = permutedims(v)   # cells × features
      ZtZ = Z' * Z + 1e-4 * I
      loadings[v_idx] = (ZtZ \ (Z' * Y))'   # features × K
    end

    dZ = zeros(Float64, n_cells, K)
    total_weight = 0.0
    for (v_idx, v) in enumerate(views)
      Y = permutedims(v)       # cells × features
      W = loadings[v_idx]      # features × K
      resid = Z * W' .- Y      # cells × features
      dZ .+= resid * W
      total_weight += size(v, 1)
    end
    dZ ./= max(total_weight, 1.0)
    Z .-= lr_f .* dZ

    if mod(iter, 100) == 0
      lr_f *= 0.7
    end
  end

  total_err = 0.0
  total_elem = 0
  for (v_idx, v) in enumerate(views)
    Y = permutedims(v)
    err = sum(abs2, Z * loadings[v_idx]' .- Y)
    total_err += err
    total_elem += length(Y)
  end
  recon_err = total_err / max(total_elem, 1)

  return (factors=Z, loadings=loadings, reconstruction_error=recon_err, n_factors=K, backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/deep_factorization_embedding"; parameters=(n_factors=K, backend=selected)))
end

# ---------------------------------------------------------------------------
# DOM-Scoped Interactive HTML Visualizers for Deep Learning Outputs
# ---------------------------------------------------------------------------

function _extract_latent_matrix(embedding)
  if typeof(embedding) <: AbstractMatrix
    return Matrix{Float64}(embedding)
  elseif hasproperty(embedding, :latent)
    return Matrix{Float64}(getproperty(embedding, :latent))
  elseif hasproperty(embedding, :factors)
    return Matrix{Float64}(getproperty(embedding, :factors))
  else
    throw(ArgumentError("embedding object does not contain a matrix or .latent property"))
  end
end

function _escape_html_text(s::AbstractString)
  return replace(String(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;", "'" => "&#39;")
end

function _wrap_html_document(content::String, title::String; standalone::Bool=true)
  if standalone
    esc_title = _escape_html_text(title)
    return """<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>$(esc_title)</title>
</head>
<body style="background:#090d16; margin:0; padding:20px;">
$(content)
</body>
</html>"""
  else
    return content
  end
end

"""
    embedding_to_html(embedding; labels=nothing, title="Latent Embedding Scatter Plot")

Generate an interactive HTML5/Plotly scatter plot for 2D/3D latent representations
(e.g., scVI, WNN, SAE, Contrastive, or Protein embeddings).
"""
function embedding_to_html(embedding; labels=nothing, title::String="Deep Learning Latent Embedding", standalone::Bool=true)
  mat = _extract_latent_matrix(embedding)
  n_cells, n_dims = size(mat)
  n_dims >= 2 || throw(ArgumentError("embedding must have at least 2 dimensions for visualization"))

  x_vals = mat[:, 1]
  y_vals = mat[:, 2]
  z_vals = n_dims >= 3 ? mat[:, 3] : Float64[]

  cell_labels = labels !== nothing ? String.(labels) : ["Cell $i" for i in 1:n_cells]
  unique_groups = sort!(unique(cell_labels))

  container_id = string("dl-emb-", string(uuid4())[1:8])
  data_payload = Dict(
    "x" => x_vals,
    "y" => y_vals,
    "z" => z_vals,
    "labels" => cell_labels,
    "groups" => unique_groups,
    "is_3d" => n_dims >= 3,
    "n_cells" => n_cells
  )

  json_data = JSON.json(data_payload)
  esc_title = _escape_html_text(title)

  snippet = """
  <div id="$(container_id)-wrapper" style="width:100%; max-width:960px; margin:20px auto; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif; background:#0f172a; color:#f8fafc; border-radius:12px; padding:20px; box-shadow:0 10px 25px rgba(0,0,0,0.5);">
    <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #334155; padding-bottom:12px; margin-bottom:16px;">
      <h3 style="margin:0; font-size:1.25rem; font-weight:700; color:#38bdf8;">$(esc_title)</h3>
      <span style="font-size:0.85rem; background:#1e293b; padding:4px 10px; border-radius:20px; border:1px solid #475569; color:#94a3b8;">N = $(n_cells) cells | $(n_dims)D Latent Space</span>
    </div>
    <div id="$(container_id)" style="width:100%; height:550px; background:#0f172a; border-radius:8px;"></div>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <script>
      (function() {
        const payload = $(json_data);
        const container = document.getElementById('$(container_id)');
        const groups = payload.groups;
        const traces = [];

        const colorPalette = ['#38bdf8', '#f43f5e', '#10b981', '#a855f7', '#f59e0b', '#06b6d4', '#ec4899', '#8b5cf6', '#84cc16'];

        groups.forEach((grp, idx) => {
          const indices = [];
          for (let i = 0; i < payload.labels.length; i++) {
            if (payload.labels[i] === grp) indices.push(i);
          }
          const traceColor = colorPalette[idx % colorPalette.length];
          if (payload.is_3d) {
            traces.push({
              x: indices.map(i => payload.x[i]),
              y: indices.map(i => payload.y[i]),
              z: indices.map(i => payload.z[i]),
              mode: 'markers',
              name: grp,
              type: 'scatter3d',
              marker: { size: 4, color: traceColor, opacity: 0.8 }
            });
          } else {
            traces.push({
              x: indices.map(i => payload.x[i]),
              y: indices.map(i => payload.y[i]),
              mode: 'markers',
              name: grp,
              type: 'scatter',
              marker: { size: 6, color: traceColor, opacity: 0.85 }
            });
          }
        });

        const layout = {
          margin: { l: 40, r: 40, b: 40, t: 40 },
          paper_bgcolor: '#0f172a',
          plot_bgcolor: '#0f172a',
          font: { color: '#cbd5e1' },
          legend: { orientation: 'h', y: -0.15 },
          xaxis: { title: 'Latent Dim 1', gridcolor: '#1e293b', zerolinecolor: '#334155' },
          yaxis: { title: 'Latent Dim 2', gridcolor: '#1e293b', zerolinecolor: '#334155' }
        };

        Plotly.newPlot(container, traces, layout, { responsive: true });
      })();
    </script>
  </div>
  """
  return _wrap_html_document(snippet, title; standalone=standalone)
end

function export_embedding_html(embedding, filepath::AbstractString; kwargs...)
  html = embedding_to_html(embedding; standalone=true, kwargs...)
  write(filepath, html)
  return String(filepath)
end

"""
    grn_to_html(edges; title="Gene Regulatory Network")

Generate an interactive HTML graph visualization for Gene Regulatory Networks.
"""
function grn_to_html(edges::DataFrame; title::String="Gene Regulatory Network", standalone::Bool=true)
  container_id = string("dl-grn-", string(uuid4())[1:8])

  src_col = "source" in names(edges) ? :source : (:source_gene in names(edges) ? :source_gene : names(edges)[1])
  tgt_col = "target" in names(edges) ? :target : (:target_gene in names(edges) ? :target_gene : names(edges)[2])
  w_col = "weight" in names(edges) ? :weight : (:propagated_weight in names(edges) ? :propagated_weight : names(edges)[3])

  nodes_set = Set{String}()
  for row in eachrow(edges)
    push!(nodes_set, String(row[src_col]))
    push!(nodes_set, String(row[tgt_col]))
  end

  nodes_list = [Dict("id" => g, "label" => g) for g in nodes_set]
  edges_list = [Dict("from" => String(r[src_col]), "to" => String(r[tgt_col]), "value" => Float64(r[w_col]), "title" => "Weight: $(round(Float64(r[w_col]), digits=4))") for r in eachrow(edges)]

  data_payload = Dict("nodes" => nodes_list, "edges" => edges_list)
  json_data = JSON.json(data_payload)
  esc_title = _escape_html_text(title)

  snippet = """
  <div id="$(container_id)-wrapper" style="width:100%; max-width:960px; margin:20px auto; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif; background:#0f172a; color:#f8fafc; border-radius:12px; padding:20px; box-shadow:0 10px 25px rgba(0,0,0,0.5);">
    <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #334155; padding-bottom:12px; margin-bottom:16px;">
      <h3 style="margin:0; font-size:1.25rem; font-weight:700; color:#38bdf8;">$(esc_title)</h3>
      <span style="font-size:0.85rem; background:#1e293b; padding:4px 10px; border-radius:20px; border:1px solid #475569; color:#94a3b8;">Nodes: $(length(nodes_set)) | Edges: $(nrow(edges))</span>
    </div>
    <div id="$(container_id)" style="width:100%; height:550px; background:#020617; border-radius:8px; border:1px solid #1e293b;"></div>
    <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
    <script>
      (function() {
        const payload = $(json_data);
        const container = document.getElementById('$(container_id)');
        const data = {
          nodes: new vis.DataSet(payload.nodes),
          edges: new vis.DataSet(payload.edges)
        };
        const options = {
          nodes: {
            shape: 'dot',
            size: 16,
            font: { color: '#f8fafc', size: 14 },
            color: { background: '#0284c7', border: '#38bdf8', highlight: { background: '#f43f5e', border: '#fb7185' } }
          },
          edges: {
            color: { color: '#475569', highlight: '#f43f5e' },
            arrows: { to: { enabled: true, scaleFactor: 0.5 } },
            smooth: { type: 'continuous' }
          },
          physics: {
            barnesHut: { gravitationalConstant: -3000, centralGravity: 0.3, springLength: 95 }
          }
        };
        new vis.Network(container, data, options);
      })();
    </script>
  </div>
  """
  return _wrap_html_document(snippet, title; standalone=standalone)
end

function export_grn_html(edges::DataFrame, filepath::AbstractString; kwargs...)
  html = grn_to_html(edges; standalone=true, kwargs...)
  write(filepath, html)
  return String(filepath)
end

"""
    trajectory_to_html(trajectory_result; title="Cell Pseudotime Trajectory")

Generate an interactive HTML plot visualizing pseudotime progression and trajectory curves.
"""
function trajectory_to_html(traj_res; title::String="Cell Pseudotime Neural ODE Trajectory", standalone::Bool=true)
  container_id = string("dl-traj-", string(uuid4())[1:8])
  pseudotime = Float64.(traj_res.pseudotime)
  n_cells = length(pseudotime)

  data_payload = Dict(
    "pseudotime" => pseudotime,
    "n_cells" => n_cells
  )
  json_data = JSON.json(data_payload)
  esc_title = _escape_html_text(title)

  snippet = """
  <div id="$(container_id)-wrapper" style="width:100%; max-width:960px; margin:20px auto; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif; background:#0f172a; color:#f8fafc; border-radius:12px; padding:20px; box-shadow:0 10px 25px rgba(0,0,0,0.5);">
    <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #334155; padding-bottom:12px; margin-bottom:16px;">
      <h3 style="margin:0; font-size:1.25rem; font-weight:700; color:#38bdf8;">$(esc_title)</h3>
      <span style="font-size:0.85rem; background:#1e293b; padding:4px 10px; border-radius:20px; border:1px solid #475569; color:#94a3b8;">N = $(n_cells) cells</span>
    </div>
    <div id="$(container_id)" style="width:100%; height:450px; background:#0f172a; border-radius:8px;"></div>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <script>
      (function() {
        const payload = $(json_data);
        const container = document.getElementById('$(container_id)');
        const sortedPt = [...payload.pseudotime].sort((a,b) => a - b);
        const trace = {
          x: Array.from({length: payload.n_cells}, (_, i) => i + 1),
          y: sortedPt,
          mode: 'lines+markers',
          type: 'scatter',
          line: { color: '#38bdf8', width: 3 },
          marker: { size: 6, color: sortedPt, colorscale: 'Viridis', showscale: true }
        };
        const layout = {
          margin: { l: 40, r: 40, b: 40, t: 40 },
          paper_bgcolor: '#0f172a',
          plot_bgcolor: '#0f172a',
          font: { color: '#cbd5e1' },
          xaxis: { title: 'Cell Rank (Trajectory Order)', gridcolor: '#1e293b' },
          yaxis: { title: 'Estimated Pseudotime [0, 1]', gridcolor: '#1e293b' }
        };
        Plotly.newPlot(container, [trace], layout, { responsive: true });
      })();
    </script>
  </div>
  """
  return _wrap_html_document(snippet, title; standalone=standalone)
end

function export_trajectory_html(traj_res, filepath::AbstractString; kwargs...)
  html = trajectory_to_html(traj_res; standalone=true, kwargs...)
  write(filepath, html)
  return String(filepath)
end

"""
    cell_type_annotation_html(annotation_result; title="Cell Type Annotation Report")

Generate an interactive HTML breakdown of predicted cell types and assignment probabilities.
"""
function cell_type_annotation_html(res; title::String="Cell Type Annotation Report", standalone::Bool=true)
  container_id = string("dl-annot-", string(uuid4())[1:8])

  labels = hasproperty(res, :predicted_label) ? res.predicted_label : (hasproperty(res, :cells) ? res.cells.predicted_label : String[])
  n_cells = length(labels)
  counts_dict = Dict{String,Int}()
  for l in labels
    counts_dict[l] = get(counts_dict, l, 0) + 1
  end

  cell_types = sort!(collect(keys(counts_dict)))
  type_counts = [counts_dict[ct] for ct in cell_types]

  data_payload = Dict("types" => cell_types, "counts" => type_counts, "total" => n_cells)
  json_data = JSON.json(data_payload)
  esc_title = _escape_html_text(title)

  snippet = """
  <div id="$(container_id)-wrapper" style="width:100%; max-width:960px; margin:20px auto; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif; background:#0f172a; color:#f8fafc; border-radius:12px; padding:20px; box-shadow:0 10px 25px rgba(0,0,0,0.5);">
    <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #334155; padding-bottom:12px; margin-bottom:16px;">
      <h3 style="margin:0; font-size:1.25rem; font-weight:700; color:#38bdf8;">$(esc_title)</h3>
      <span style="font-size:0.85rem; background:#1e293b; padding:4px 10px; border-radius:20px; border:1px solid #475569; color:#94a3b8;">Total Cells: $(n_cells)</span>
    </div>
    <div id="$(container_id)" style="width:100%; height:450px; background:#0f172a; border-radius:8px;"></div>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <script>
      (function() {
        const payload = $(json_data);
        const container = document.getElementById('$(container_id)');
        const trace = {
          labels: payload.types,
          values: payload.counts,
          type: 'pie',
          hole: 0.45,
          textinfo: 'label+percent',
          marker: { colors: ['#38bdf8', '#f43f5e', '#10b981', '#a855f7', '#f59e0b', '#06b6d4'] }
        };
        const layout = {
          margin: { l: 40, r: 40, b: 40, t: 40 },
          paper_bgcolor: '#0f172a',
          font: { color: '#cbd5e1' },
          showlegend: true
        };
        Plotly.newPlot(container, [trace], layout, { responsive: true });
      })();
    </script>
  </div>
  """
  return _wrap_html_document(snippet, title; standalone=standalone)
end

function export_cell_type_annotation_html(res, filepath::AbstractString; kwargs...)
  html = cell_type_annotation_html(res; standalone=true, kwargs...)
  write(filepath, html)
  return String(filepath)
end

import ..BlenderIntegrator: to_blender_payload, BlenderMaterial, BlenderSpatialPayload

"""
    to_blender_payload(lm::LatentModelResult; name="LatentEmbedding", glyph_scale=0.25)

Convert a neural network `LatentModelResult` into a 3D `BlenderSpatialPayload` for Blender rendering.
"""
function to_blender_payload(lm::LatentModelResult; name::String="LatentEmbedding", glyph_scale::Real=0.25)
  latent = lm.latent
  n = size(latent, 1)
  coords = size(latent, 2) >= 3 ? latent[:, 1:3] : (size(latent, 2) == 2 ? hcat(latent, zeros(n)) : hcat(latent, zeros(n, 2)))
  labels = fill("LatentCell", n)
  colors = [begin
    hue = (i - 1) / max(n, 1)
    r = abs(hue * 6 - 3) - 1
    g = 2 - abs(hue * 6 - 2)
    b = 2 - abs(hue * 6 - 4)
    (clamp(r, 0.0, 1.0), clamp(g, 0.0, 1.0), clamp(b, 0.0, 1.0))
  end for i in 1:n]
  mat = BlenderMaterial(name=name * "_mat")
  return BlenderSpatialPayload(name, Matrix{Float64}(coords), labels, colors, Float64(glyph_scale), mat)
end


end

