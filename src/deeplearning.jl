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
import ..BioToolkit
using ..BioToolkit: AminoAcidAlphabet, BioAlphabet, BioSequence, flux_available, maybe_to_device, maybe_to_host, resolve_backend, threaded_foreach
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_record, provenance_result!, register_provenance!, with_provenance

@inline function _register_dl_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export scvi_like_embedding, cellassign_like_mapping, geneformer_like_embedding, scgpt_like_embedding, scbert_like_embedding, attention_grn
export batch_corrected_latent, contrastive_cell_embedding
export flux_autoencoder_embedding, flux_mlp_classifier
export scgen_like_perturbation, graphsca_label_transfer

# New exports
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

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

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

    function hash_kmer(kmer::AbstractString)
        h = UInt(1469598103934665603)
        for c in codeunits(kmer)
            h = (h ⊻ UInt(c)) * UInt(1099511628211)
        end
        return Int(mod(h, UInt(dim))) + 1
    end

    emb = zeros(Float64, length(sequences), dim)
    threaded_foreach(length(sequences), i -> begin
        s = uppercase(String(sequences[i]))
        if ncodeunits(s) < k
            return
        end
        for start in 1:(ncodeunits(s) - k + 1)
            idx = hash_kmer(s[start:(start + k - 1)])
            emb[i, idx] += 1.0
        end
        normv = norm(@view emb[i, :])
        normv > eps(Float64) && (emb[i, :] ./= normv)
    end; threaded=threaded)
    _ctx = active_provenance_context()


    return _register_dl_result!(_ctx, emb, "geneformer_like_embedding"; parents=provenance_parent_ids(sequences), parameters=(k=k, dim=dim, seq_count=length(sequences)))
end

geneformer_like_embedding(sequences::AbstractVector{<:BioSequence{A}}; k::Int=3, dim::Int=64, threaded::Bool=true) where {A <: BioAlphabet} =
    geneformer_like_embedding(String.(sequences); k=k, dim=dim, threaded=threaded)

"""
    scgpt_like_embedding(sequences; token_dim=128)

scGPT-style sequence embedding using k-mer token histograms plus positional
context channels.
"""
function scgpt_like_embedding(sequences::AbstractVector{<:AbstractString}; token_dim::Int=128, k::Int=3, context_window::Int=8, threaded::Bool=true)
    token_dim >= 8 || throw(ArgumentError("token_dim must be >= 8"))
    base = geneformer_like_embedding(sequences; k=k, dim=token_dim, threaded=threaded)
    pos = zeros(Float64, length(sequences), token_dim)
    win = max(Int(context_window), 1)

    threaded_foreach(length(sequences), i -> begin
        s = uppercase(String(sequences[i]))
        n = ncodeunits(s)
        n == 0 && return
        for p in 1:n
            c = Int(codeunit(s, p))
            bucket = mod(c + p - 1, token_dim) + 1
            pos[i, bucket] += exp(-abs(p - (n + 1) / 2) / win)
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


    return _register_dl_result!(_ctx, emb, "scgpt_like_embedding"; parents=provenance_parent_ids(sequences), parameters=(token_dim=token_dim, k=k, context_window=context_window, seq_count=length(sequences)))
end

scgpt_like_embedding(sequences::AbstractVector{<:BioSequence{A}}; token_dim::Int=128, k::Int=3, context_window::Int=8, threaded::Bool=true) where {A <: BioAlphabet} =
    scgpt_like_embedding(String.(sequences); token_dim=token_dim, k=k, context_window=context_window, threaded=threaded)

"""
    scbert_like_embedding(sequences; dim=96)

scBERT-style character-token embedding with lightweight positional hashing.
"""
function scbert_like_embedding(sequences::AbstractVector{<:AbstractString}; dim::Int=96, max_len::Int=512, threaded::Bool=true)
    dim >= 8 || throw(ArgumentError("dim must be >= 8"))
    max_len >= 1 || throw(ArgumentError("max_len must be >= 1"))
    emb = zeros(Float64, length(sequences), dim)

    threaded_foreach(length(sequences), i -> begin
        s = uppercase(String(sequences[i]))
        n = min(ncodeunits(s), max_len)
        n == 0 && return
        emb[i, 1] += 1.0
        for p in 1:n
            c = Int(codeunit(s, p))
            bucket = mod((c * 131 + p * 17), dim - 1) + 2
            emb[i, bucket] += 1.0
        end
        nrm = norm(@view emb[i, :])
        nrm > eps(Float64) && (emb[i, :] ./= nrm)
    end; threaded=threaded)

    _ctx = active_provenance_context()


    return _register_dl_result!(_ctx, emb, "scbert_like_embedding"; parents=provenance_parent_ids(sequences), parameters=(dim=dim, max_len=max_len, seq_count=length(sequences)))
end

scbert_like_embedding(sequences::AbstractVector{<:BioSequence{A}}; dim::Int=96, max_len::Int=512, threaded::Bool=true) where {A <: BioAlphabet} =
    scbert_like_embedding(String.(sequences); dim=dim, max_len=max_len, threaded=threaded)

function attention_grn(expression::AbstractMatrix{<:Real}; top_k::Int=100, temperature::Real=1.0)
    X = _zscore_columns(permutedims(expression))
    score = (X' * X) ./ sqrt(size(X, 1)) ./ max(Float64(temperature), eps(Float64))
    score[diagind(score)] .= -Inf

    n = size(score, 1)
    edges = DataFrame(source=String[], target=String[], weight=Float64[])
    for i in 1:n
        row = vec(@view score[i, :])
        order = sortperm(row, rev=true)
        for j in Iterators.take(order, min(top_k, n - 1))
            isfinite(row[j]) || continue
            push!(edges, ("gene_$(i)", "gene_$(j)", row[j]))
        end
    end
    sort!(edges, :weight, rev=true)
    _ctx = active_provenance_context()


    return _register_dl_result!(_ctx, edges, "attention_grn"; parents=provenance_parent_ids(expression), parameters=(top_k=top_k, temperature=Float64(temperature), n_genes=size(expression, 1), edge_count=nrow(edges)))
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

"""
    batch_corrected_latent(counts, batches; n_latent=10)

scVI-style latent extraction followed by linear residualization of batch effects.
"""
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

"""
    contrastive_cell_embedding(counts; n_latent=16, dropout=0.15)

Contrastive proxy embedding using two dropout-augmented views of log counts.
"""
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
    return getfield(BioToolkit, :Flux)
end

"""
    flux_autoencoder_embedding(counts; latent_dim=16)

Train a lightweight Flux autoencoder and return latent embeddings.
"""
function flux_autoencoder_embedding(counts::AbstractMatrix{<:Real}; latent_dim::Int=16, hidden_dim::Int=64, epochs::Int=25, lr::Real=1e-3, backend::Symbol=:auto, seed::Int=1)
    F = _require_flux_module()
    Random.seed!(seed)

    X = Float32.(permutedims(log1p.(Float64.(counts))))
    n_cells, n_genes = size(X)
    ldim = clamp(latent_dim, 2, max(2, min(n_genes, n_cells)))
    hdim = clamp(hidden_dim, ldim, max(ldim, n_genes))

    encoder = F.Chain(F.Dense(n_genes => hdim, relu), F.Dense(hdim => ldim))
    decoder = F.Chain(F.Dense(ldim => hdim, relu), F.Dense(hdim => n_genes))
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
    latent_host = Array(F.cpu(latent))
    result = (latent=latent_host, backend=selected, latent_dim=ldim, provenance=provenance_record("DeepLearningResult", "DeepLearning/flux_autoencoder_embedding"; parameters=(latent_dim=ldim, hidden_dim=Int(hidden_dim), epochs=Int(epochs), backend=selected)))
    _ctx = active_provenance_context()


    return _register_dl_result!(_ctx, result, "flux_autoencoder_embedding"; parents=provenance_parent_ids(counts), parameters=(latent_dim=ldim, hidden_dim=Int(hidden_dim), epochs=Int(epochs), backend=selected))
end

"""
    flux_mlp_classifier(features, labels)

Train a small Flux MLP classifier and return class probabilities.
"""
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
    model = F.Chain(F.Dense(n_features => hdim, relu), F.Dense(hdim => length(classes)))

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

"""
    scgen_like_perturbation(control_expression, treated_expression, query_control)

scGen-inspired perturbation transfer using latent shift vectors.
All matrices are genes x cells.
"""
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
    zt = latent[(n_c + 1):(n_c + n_t), :]
    zq = latent[(n_c + n_t + 1):(n_c + n_t + n_q), :]

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

"""
    graphsca_label_transfer(expression, adjacency, known_labels)

GraphSCA-like semi-supervised label transfer via graph propagation.
Expression is genes x cells, adjacency is cells x cells.
"""
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

"""
    cell_type_denoising(counts; k=15, t=3, n_pcs=20, backend=:auto)

Impute and denoise single-cell expression using graph diffusion, analogous to
MAGIC (van Dijk et al. 2018) and the scVI decoder denoising mode.

Constructs a k-NN graph in PCA space, builds a Markov diffusion operator,
and applies t diffusion steps to smooth expression.

Returns `(denoised=Matrix{Float64}, diffusion_operator=Matrix{Float64})`.
"""
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

    # kNN adjacency (Gaussian kernel)
    sigma = 1.0
    W = zeros(Float64, n_cells, n_cells)
    for i in 1:n_cells
        dists = [sum(abs2, pca[i, :] .- pca[j, :]) for j in 1:n_cells]
        dists[i] = Inf
        order = sortperm(dists)
        kth_dist = dists[order[kk]]
        for j in order[1:kk]
            W[i, j] = exp(-dists[j] / max(2 * kth_dist, eps(Float64)))
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
    for _ in 1:(t - 1)
        Mt = Mt * W
    end

    denoised_cells = Mt * cells                       # cells × genes (smoothed)
    denoised = permutedims(denoised_cells)            # genes × cells
    result = (denoised=denoised, diffusion_operator=Mt, provenance=provenance_record("DeepLearningResult", "DeepLearning/cell_type_denoising"))
    _ctx = active_provenance_context()


    return _register_dl_result!(_ctx, result, "cell_type_denoising"; parents=provenance_parent_ids(counts), parameters=(k=k, t=t, n_pcs=n_pcs, backend=string(resolve_backend(; backend=backend))))
end

# ---------------------------------------------------------------------------
# Sparse Autoencoder Features (K-sparse)
# ---------------------------------------------------------------------------

"""
    sparse_autoencoder_features(counts; n_features=64, sparsity_k=8, n_iter=300, lr=5e-3)

Train a k-sparse linear autoencoder to extract interpretable gene features,
analogous to sparse autoencoders used in Geneformer mechanistic interpretability.

At each step only the top-k activations are kept (k-sparse constraint).

Returns `(features=Matrix, feature_activations=Matrix, reconstruction_loss)`.
`features` is genes × n_features (dictionary atoms).
`feature_activations` is cells × n_features.
"""
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

    selected = resolve_backend(; backend=backend)
    Xdev = maybe_to_device(permutedims(X); backend=selected)   # cells × genes

    rng = MersenneTwister(seed)
    D = maybe_to_device(randn(rng, n_genes, n_features) ./ sqrt(n_genes); backend=selected)

    # Normalise columns of D
    function normalise_dict!(D)
        for j in 1:size(D, 2)
            col = @view D[:, j]
            nrm = norm(col)
            nrm > eps(Float64) && (col ./= nrm)
        end
        return D
    end
    normalise_dict!(D)

    lr_f = Float64(lr)
    best_loss = Inf

    for iter in 1:n_iter
        # Encode: A = X * D  (cells × features), apply k-sparse mask
        A = Matrix{Float64}(maybe_to_host(Xdev)) * Matrix{Float64}(maybe_to_host(D))
        # k-sparse: keep top-k per cell
        for i in 1:n_cells
            row = @view A[i, :]
            threshold = partialsort(abs.(row), sparsity_k, rev=true)
            row[abs.(row) .< threshold] .= 0.0
        end

        Adev = maybe_to_device(A; backend=selected)
        # Reconstruct
        Xhat = Adev * D'      # cells × genes
        residual = Matrix{Float64}(maybe_to_host(Xhat .- Xdev))
        loss = mean(abs2, residual)
        best_loss = min(best_loss, loss)

        # Gradient w.r.t. D: (Xhat - X)' * A / n_cells
        grad = Matrix{Float64}(maybe_to_host(Adev))' * residual ./ n_cells  # features × genes → need transpose
        grad = grad'   # genes × features
        D_host = Matrix{Float64}(maybe_to_host(D))
        D_host .-= lr_f .* grad
        D = maybe_to_device(Matrix{Float64}(D_host); backend=selected)
        normalise_dict!(D)

        # Decay lr
        if mod(iter, 100) == 0
            lr_f *= 0.5
        end
    end

    D_final = Matrix{Float64}(maybe_to_host(D))   # genes × features
    X_host = permutedims(X)   # cells × genes
    A_final = X_host * D_final
    for i in 1:n_cells
        row = @view A_final[i, :]
        threshold = partialsort(abs.(row), sparsity_k, rev=true)
        row[abs.(row) .< threshold] .= 0.0
    end

    recon = A_final * D_final'
    final_loss = mean(abs2, recon .- X_host)

    result = (features=D_final, feature_activations=A_final, reconstruction_loss=final_loss, backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/sparse_autoencoder_features"; parameters=(backend=selected,)))
    _ctx = active_provenance_context()


    return _register_dl_result!(_ctx, result, "sparse_autoencoder_features"; parents=provenance_parent_ids(counts), parameters=(n_features=n_features, sparsity_k=sparsity_k, n_iter=n_iter, backend=selected))
end

# ---------------------------------------------------------------------------
# Trajectory Neural ODE (Euler approximation, scVelo-NN style)
# ---------------------------------------------------------------------------

"""
    trajectory_neural_ode(spliced_counts, unspliced_counts; n_latent=10, n_steps=20, dt=0.1)

Estimate RNA velocity-informed pseudotime by projecting spliced and unspliced
counts into a shared latent space and computing a Euler-integrated trajectory,
analogous to scVelo's neural ODE mode.

Returns `(pseudotime=Vector, latent_velocity=Matrix, trajectory=Array)`.
"""
function trajectory_neural_ode(
    spliced_counts::AbstractMatrix{<:Real},
    unspliced_counts::AbstractMatrix{<:Real};
    n_latent::Int=10,
    n_steps::Int=20,
    dt::Real=0.1,
    backend::Symbol=:auto)
    size(spliced_counts) == size(unspliced_counts) ||
        throw(DimensionMismatch("spliced and unspliced count matrices must have the same shape"))

    S = log1p.(Float64.(spliced_counts))    # genes × cells
    U = log1p.(Float64.(unspliced_counts))

    # Embed both modalities via SVD
    emb_s = scvi_like_embedding(spliced_counts; n_latent=n_latent, backend=backend)
    emb_u = scvi_like_embedding(unspliced_counts; n_latent=n_latent, backend=backend)

    Zs = Matrix{Float64}(emb_s.latent)   # cells × latent
    Zu = Matrix{Float64}(emb_u.latent)

    # Velocity in latent space: Δz = Zu - Zs (approximation)
    velocity = Zu .- Zs   # cells × latent

    n_cells = size(Zs, 1)
    dim_lat = size(Zs, 2)

    # Euler integration from each cell's starting position
    trajectory = zeros(Float64, n_cells, dim_lat, n_steps + 1)
    trajectory[:, :, 1] = Zs

    for step in 1:n_steps
        z_current = trajectory[:, :, step]
        # Velocity field: interpolate from nearest cell in current position
        z_next = z_current .+ Float64(dt) .* velocity
        trajectory[:, :, step + 1] = z_next
    end

    # Pseudotime: distance from initial position (norm of displacement)
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

"""
    multimodal_wnn_embedding(modalities; n_latent=20, n_neighbors=20, backend=:auto)

Integrate multiple single-cell modalities using a Weighted Nearest Neighbour
(WNN) approach analogous to Seurat v4 / Hao et al. 2021.

`modalities`: `AbstractVector` of matrices (genes/features × cells) — e.g. [rna, atac, protein].

Each modality is first reduced via SVD. WNN weights are computed per-cell based
on within-modality kNN precision.

Returns `(latent=Matrix, modality_weights=Matrix, nn_graph=Matrix)`.
"""
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
    # Per-modality SVD embedding
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

    # Compute per-modality kNN precision (within-modality)
    kk = min(n_neighbors, n_cells - 1)
    precision = zeros(Float64, n_cells, M)

    for (m_idx, Z) in enumerate(latents)
        for i in 1:n_cells
            dists = [sum(abs2, Z[i, :] .- Z[j, :]) for j in 1:n_cells]
            dists[i] = Inf
            order = sortperm(dists)
            knn_dists = dists[order[1:kk]]
            # Mean inverse distance (precision proxy)
            precision[i, m_idx] = mean(1.0 ./ max.(knn_dists, eps(Float64)))
        end
    end

    # WNN weights: softmax of precision across modalities per cell
    log_prec = log.(max.(precision, eps(Float64)))
    weights = copy(log_prec)
    for i in 1:n_cells
        row = @view weights[i, :]
        row .-= maximum(row)
        row .= exp.(row)
        row ./= max(sum(row), eps(Float64))
    end

    # Weighted concatenation → combined latent
    combined_dim = minimum(size(l, 2) for l in latents)
    combined = zeros(Float64, n_cells, combined_dim)
    for (m_idx, Z) in enumerate(latents)
        used_dim = min(combined_dim, size(Z, 2))
        for i in 1:n_cells
            combined[i, 1:used_dim] .+= weights[i, m_idx] .* Z[i, 1:used_dim]
        end
    end

    # Final kNN graph in the combined space
    nn_graph = zeros(Float64, n_cells, n_cells)
    for i in 1:n_cells
        dists = [sum(abs2, combined[i, :] .- combined[j, :]) for j in 1:n_cells]
        dists[i] = Inf
        order = sortperm(dists)
        for j in order[1:kk]
            nn_graph[i, j] = 1.0
        end
    end

    result = (latent=combined, modality_weights=weights, nn_graph=nn_graph, backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/multimodal_wnn_embedding"; parameters=(backend=selected,)))
    _ctx = active_provenance_context()


    return _register_dl_result!(_ctx, result, "multimodal_wnn_embedding"; parents=String[], parameters=(n_modalities=length(modalities), n_latent=n_latent, n_neighbors=n_neighbors, backend=selected))
end

# ---------------------------------------------------------------------------
# Protein Sequence Embedding (ESM-like)
# ---------------------------------------------------------------------------

"""
    protein_sequence_embedding(sequences; dim=128, k=3, include_properties=true)

Compute protein sequence embeddings using k-mer hashing plus per-residue
physicochemical features (charge, hydrophobicity, secondary-structure propensity),
analogous to ESM-2 style embeddings for downstream tasks.

Returns an `n_sequences × dim` embedding matrix.
"""
function protein_sequence_embedding(
    sequences::AbstractVector{<:AbstractString};
    dim::Int=128,
    k::Int=3,
    include_properties::Bool=true,
    threaded::Bool=true)
    dim >= 8 || throw(ArgumentError("dim must be >= 8"))

    # AA property look-up tables
    kd = Dict('I'=>4.5,'V'=>4.2,'L'=>3.8,'F'=>2.8,'C'=>2.5,'M'=>1.9,'A'=>1.8,
              'G'=>-0.4,'T'=>-0.7,'S'=>-0.8,'W'=>-0.9,'Y'=>-1.3,'P'=>-1.6,
              'H'=>-3.2,'E'=>-3.5,'Q'=>-3.5,'D'=>-3.5,'N'=>-3.5,'K'=>-3.9,'R'=>-4.5)
    charge_aa = Dict('R'=>1.0,'K'=>1.0,'H'=>0.1,'D'=>-1.0,'E'=>-1.0)
    helix_prop = Dict('A'=>1.42,'L'=>1.21,'M'=>1.45,'E'=>1.51,'K'=>1.16,'R'=>0.98)
    sheet_prop = Dict('V'=>1.70,'I'=>1.60,'T'=>1.19,'Y'=>1.47,'W'=>1.35,'F'=>1.38)

    base_emb = geneformer_like_embedding(sequences; k=k, dim=dim, threaded=threaded)

    if !include_properties
        return base_emb
    end

    prop_buckets = 8
    prop_emb = zeros(Float64, length(sequences), prop_buckets)

    threaded_foreach(length(sequences), idx -> begin
        seq = uppercase(String(sequences[idx]))
        aa  = collect(seq)
        L   = length(aa)
        L == 0 && return

        hydro   = mean(get(kd, c, 0.0) for c in aa)
        charge  = sum(get(charge_aa, c, 0.0) for c in aa) / max(L, 1)
        helix   = mean(get(helix_prop, c, 1.0) for c in aa)
        sheet   = mean(get(sheet_prop, c, 1.0) for c in aa)
        mw_proxy = L * 110.0 / 1000.0   # rough MW in kDa

        prop_emb[idx, 1] = hydro / 5.0           # normalised hydrophobicity
        prop_emb[idx, 2] = clamp(charge, -1.0, 1.0)
        prop_emb[idx, 3] = (helix - 1.0) / 0.5
        prop_emb[idx, 4] = (sheet - 1.0) / 0.5
        prop_emb[idx, 5] = log(mw_proxy + 1) / 5.0
        prop_emb[idx, 6] = count(c -> c in Set(['C']), aa) / max(L, 1)  # Cys fraction
        prop_emb[idx, 7] = count(c -> c == 'P', aa) / max(L, 1)         # Pro fraction
        prop_emb[idx, 8] = count(c -> c in Set(['W','Y','F']), aa) / max(L, 1)
    end; threaded=threaded)

    # Concatenate base (dim-8 dims) + 8 property features, then project back to dim
    prop_emb_full = hcat(base_emb[:, 1:(dim-prop_buckets)], prop_emb)

    # L2 normalise each row
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

"""
    zero_shot_cell_annotation(counts, gene_ids, archetype_profiles; top_k=1)

Annotate cells without labelled training data by comparing each cell's
expression profile to precomputed archetype vectors (e.g. from CellTypist,
Tabula Sapiens compressed profiles), analogous to zero-shot scType / CellTypist.

`archetype_profiles`: A `Dict{String, Vector{String}}` mapping cell-type label
to a vector of marker gene names (positive markers).

Returns a `DataFrame` with `cell_index`, `predicted_label`, `score`,
and a full score matrix.
"""
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

    # Softmax normalise scores per cell
    score_norm = _softmax_rows(score_matrix)

    # Top-k predictions
    pred_labels = String[]
    pred_scores = Float64[]
    for i in 1:n_cells
        row = score_norm[i, :]
        best_idx = argmax(row)
        push!(pred_labels, labels[best_idx])
        push!(pred_scores, row[best_idx])
    end

    cell_df = DataFrame(
        cell_index      = 1:n_cells,
        predicted_label = pred_labels,
        score           = pred_scores)
    return (cells=with_provenance(cell_df, "DeepLearningTable", "DeepLearning/zero_shot_cell_annotation"; parameters=(label_count=length(labels),)), score_matrix=score_norm, label_order=labels, provenance=provenance_record("DeepLearningResult", "DeepLearning/zero_shot_cell_annotation"; parameters=(label_count=length(labels),)))
end

# ---------------------------------------------------------------------------
# Gene Regulatory Network via GNN propagation
# ---------------------------------------------------------------------------

"""
    gene_regulatory_network_gnn(expression; k_neighbors=10, n_propagation=3, top_k_edges=200, temperature=1.0)

Infer a gene regulatory network (GRN) using multi-step graph neural network
propagation on a cell-cell similarity graph in PCA space, integrating
co-expression across neighbourhoods — analogous to SCENIC+ / GRNBoost2 GNN mode.

Returns an edge `DataFrame` with `source_gene`, `target_gene`, `raw_weight`,
`propagated_weight`, and `rank`.
"""
function gene_regulatory_network_gnn(
    expression::AbstractMatrix{<:Real};
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

    # PCA for cell similarity
    cells_centered = cells .- mean(cells, dims=1)
    work = maybe_to_device(cells_centered; backend=selected)
    fac = svd(work; full=false)
    U = Matrix{Float64}(maybe_to_host(fac.U))
    S_vals = Vector{Float64}(maybe_to_host(fac.S))
    pca = U[:, 1:min(n_pcs, size(U, 2))] * Diagonal(S_vals[1:min(n_pcs, size(U, 2))])

    # kNN cell graph → row-normalised propagation matrix
    kk = min(k_neighbors, n_cells - 1)
    W = zeros(Float64, n_cells, n_cells)
    for i in 1:n_cells
        dists = [sum(abs2, pca[i, :] .- pca[j, :]) for j in 1:n_cells]
        dists[i] = Inf
        order = sortperm(dists)
        for j in order[1:kk]
            W[i, j] = 1.0
        end
    end
    row_s = vec(sum(W, dims=2))
    for i in 1:n_cells
        row_s[i] > 0 && (W[i, :] ./= row_s[i])
    end

    # Raw gene-gene co-expression (attention-like score)
    Xz = _zscore_columns(X)   # genes × cells (z-scored per cell)
    raw_score = (Xz * Xz') ./ sqrt(n_cells) ./ max(Float64(temperature), eps(Float64))
    raw_score[diagind(raw_score)] .= -Inf

    # Propagate: convolve gene expression over the cell graph n_propagation times
    X_prop = copy(X)
    for _ in 1:n_propagation
        X_prop = X_prop * W'   # genes × cells, now smoothed
    end
    Xpz = _zscore_columns(X_prop)
    prop_score = (Xpz * Xpz') ./ sqrt(n_cells) ./ max(Float64(temperature), eps(Float64))
    prop_score[diagind(prop_score)] .= -Inf

    # Extract top edges
    raw_vec   = Float64[]
    prop_vec  = Float64[]
    src_idx   = Int[]
    tgt_idx   = Int[]

    for i in 1:n_genes
        row = @view prop_score[i, :]
        order = sortperm(vec(row), rev=true)
        for j in order[1:min(top_k_edges ÷ n_genes + 1, n_genes - 1)]
            isfinite(row[j]) || continue
            push!(src_idx, i)
            push!(tgt_idx, j)
            push!(raw_vec, raw_score[i, j])
            push!(prop_vec, row[j])
        end
    end

    edges = DataFrame(
        source_gene       = ["gene_$(i)" for i in src_idx],
        target_gene       = ["gene_$(j)" for j in tgt_idx],
        raw_weight        = raw_vec,
        propagated_weight = prop_vec)
    sort!(edges, :propagated_weight, rev=true)
    edges = first(edges, min(top_k_edges, nrow(edges)))
    edges[!, :rank] = 1:nrow(edges)
    return (edges=with_provenance(edges, "DeepLearningTable", "DeepLearning/gene_regulatory_network_gnn"), backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/gene_regulatory_network_gnn"; parameters=(backend=selected,)))
end

# ---------------------------------------------------------------------------
# Self-Supervised Pre-training (Masked Gene Modelling)
# ---------------------------------------------------------------------------

"""
    self_supervised_pretraining(counts; mask_fraction=0.15, n_latent=32, n_iter=200, lr=1e-3)

Self-supervised masked gene modelling pre-training, analogous to BERT-style
masking applied in scBERT and Geneformer.

Randomly masks `mask_fraction` of genes for each cell, trains a linear encoder
to reconstruct masked values, and returns the learned latent representations.

Returns `(latent=Matrix, reconstruction_loss_history=Vector)`.
"""
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
    # Encoder: cells × genes → cells × n_latent   (linear)
    W_enc = randn(rng, n_genes, n_latent) ./ sqrt(n_genes)
    # Decoder: cells × n_latent → cells × n_genes  (linear)
    W_dec = randn(rng, n_latent, n_genes) ./ sqrt(n_latent)

    loss_history = Float64[]

    for iter in 1:n_iter
        # Generate random mask
        mask = rand(rng, n_cells, n_genes) .< Float64(mask_fraction)
        X_masked = copy(cells)
        X_masked[mask] .= 0.0

        # Forward pass
        Z = X_masked * W_enc                        # cells × latent
        Xhat = Z * W_dec                             # cells × genes

        # Loss only on masked positions
        residual = (Xhat .- cells) .* mask
        loss = sum(abs2, residual) / max(sum(mask), 1)
        push!(loss_history, loss)

        # Gradients (simple SGD)
        dL_dXhat = 2 .* residual ./ max(sum(mask), 1)
        dL_dWdec = Z' * dL_dXhat                    # latent × genes
        dL_dZ    = dL_dXhat * W_dec'                # cells × latent
        dL_dWenc = X_masked' * dL_dZ               # genes × latent

        W_dec .-= lr_f .* dL_dWdec
        W_enc .-= lr_f .* dL_dWenc

        if mod(iter, 50) == 0
            lr_f *= 0.9
        end
    end

    # Final latent representations (no masking)
    latent = cells * W_enc   # cells × latent
    return (latent=latent, reconstruction_loss_history=loss_history, backend=selected, provenance=provenance_record("DeepLearningResult", "DeepLearning/self_supervised_pretraining"; parameters=(backend=selected,)))
end

# ---------------------------------------------------------------------------
# Cell Cycle Regression
# ---------------------------------------------------------------------------

"""
    cell_cycle_regression(counts, gene_ids; s_genes=String[], g2m_genes=String[])

Score cells for cell-cycle phase and regress out cell-cycle effects from
expression, analogous to `Seurat::CellCycleScoring` + `ScaleData(vars.to.regress)`.

`s_genes`   : S-phase marker genes (e.g. from Tirosh et al. 2016).
`g2m_genes` : G2M-phase marker genes.

Returns `(corrected_counts=Matrix, s_score=Vector, g2m_score=Vector, phase=Vector)`.
"""
function cell_cycle_regression(
    counts::AbstractMatrix{<:Real},
    gene_ids::AbstractVector{<:AbstractString};
    s_genes::AbstractVector{<:AbstractString}=String[],
    g2m_genes::AbstractVector{<:AbstractString}=String[],
    ridge::Real=1e-3)
    n_genes, n_cells = size(counts)
    length(gene_ids) == n_genes || throw(DimensionMismatch("gene_ids must match expression rows"))

    gene_index = Dict(String(g) => i for (i, g) in enumerate(gene_ids))
    X = log1p.(Float64.(counts))  # genes × cells

    function score_module(gene_list)
        idx = [gene_index[g] for g in String.(gene_list) if haskey(gene_index, g)]
        isempty(idx) && return zeros(Float64, n_cells)
        vec(mean(X[idx, :], dims=1))
    end

    s_score   = score_module(s_genes)
    g2m_score = score_module(g2m_genes)

    # Phase assignment
    phase = [s > g2m && s > 0.1 ? "S" :
             g2m > s && g2m > 0.1 ? "G2M" : "G1"
             for (s, g2m) in zip(s_score, g2m_score)]

    # Regress out S and G2M scores from each gene via ridge regression
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

"""
    deep_factorization_embedding(data_matrices; n_factors=10, n_iter=300, lr=1e-2, seed=1)

Multi-omic factorization analogous to MOFA+ (Argelaguet et al. 2020).
Factorises multiple data matrices sharing the same cells into shared latent
factors and view-specific loadings using alternating least squares.

`data_matrices`: `AbstractVector` of matrices, each `features_k × n_cells`.

Returns `(factors=Matrix, loadings=Vector{Matrix}, reconstruction_error=Float64)`.
"""
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

    # Normalize each view
    views = [log1p.(Float64.(m)) for m in data_matrices]

    # Initialise shared factor matrix Z (cells × K)
    Z = randn(rng, n_cells, K) ./ sqrt(K)

    # Initialise view-specific loadings W_v (features_v × K)
    loadings = [randn(rng, size(v, 1), K) ./ sqrt(K) for v in views]

    lr_f = Float64(lr)
    for iter in 1:n_iter
        # Update loadings for each view (least squares given Z)
        for (v_idx, v) in enumerate(views)
            Y = permutedims(v)   # cells × features
            # W_v = (Z'Z + ridge*I)^{-1} Z'Y transposed → features × K
            ZtZ = Z' * Z + 1e-4 * I
            loadings[v_idx] = (ZtZ \ (Z' * Y))'   # features × K
        end

        # Update Z (shared factors) given all loadings
        dZ = zeros(Float64, n_cells, K)
        total_weight = 0.0
        for (v_idx, v) in enumerate(views)
            Y = permutedims(v)       # cells × features
            W = loadings[v_idx]      # features × K
            # Gradient: (Z W' - Y) W
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

    # Compute reconstruction error
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

end
