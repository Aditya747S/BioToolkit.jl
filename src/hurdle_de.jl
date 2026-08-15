# ==============================================================================
# hurdle_de.jl — MAST-style hurdle model and pseudobulk workflows
# ==============================================================================

using Distributions: Chisq, ccdf
using LinearAlgebra
using SparseArrays
using DataFrames

struct HurdleDEResult <: AbstractAnalysisResult
    gene_id::String
    log2_fc_continuous::Float64
    log2_fc_detection::Float64
    pvalue_cont::Float64
    pvalue_det::Float64
    pvalue_hurdle::Float64
    padj::Float64
    provenance::ResultProvenance
end

function fit_logistic(X::Matrix{Float64}, y::Vector{Float64}; max_iter::Int=25, tol::Float64=1e-6)
    n, p = size(X)
    beta = zeros(p)
    
    for iter in 1:max_iter
        p_val = 1.0 ./ (1.0 .+ exp.(-X * beta))
        p_val = clamp.(p_val, 1e-15, 1.0 - 1e-15)
        
        W = p_val .* (1.0 .- p_val)
        grad = X' * (y .- p_val)
        
        H = (X' .* W') * X
        H += 1e-4 * I  # Firth-like ridge stabilization
        
        step = H \ grad
        beta += step
        
        if norm(step) < tol
            break
        end
    end
    
    p_val = 1.0 ./ (1.0 .+ exp.(-X * beta))
    p_val = clamp.(p_val, 1e-15, 1.0 - 1e-15)
    loglik = sum(y .* log.(p_val) .+ (1.0 .- y) .* log.(1.0 .- p_val))
    
    return beta, loglik
end

function fit_linear(X::Matrix{Float64}, y::Vector{Float64})
    nc, p = size(X)
    if nc <= p || rank(X) < p
        return zeros(p), -Inf
    end
    XtX = X' * X
    XtX += 1e-6 * I
    beta = XtX \ (X' * y)
    residuals = y .- X * beta
    rss = sum(residuals.^2)
    sigma2 = rss / nc
    if sigma2 <= 1e-12
        sigma2 = 1e-12
    end
    loglik = -0.5 * nc * (log(2 * pi * sigma2) + 1.0)
    return beta, loglik
end

function mast_hurdle_test(counts::AbstractMatrix{<:Real}, design::AbstractVector; gene_ids=nothing, min_cells::Int=5)
    _ctx = active_provenance_context()
    
    n_genes, n_cells = size(counts)
    genes = gene_ids === nothing ? ["gene_$(i)" for i in 1:n_genes] : String.(gene_ids)
    
    unique_conds = unique(design)
    if length(unique_conds) < 2
        throw(ArgumentError("design must contain at least 2 distinct conditions"))
    end
    
    X_full = ones(n_cells, length(unique_conds))
    for j in 2:length(unique_conds)
        cond = unique_conds[j]
        for i in 1:n_cells
            X_full[i, j] = (design[i] == cond) ? 1.0 : 0.0
        end
    end
    X_reduced = ones(n_cells, 1)
    
    df_df = length(unique_conds) - 1
    
    results = HurdleDEResult[]
    pvalues_hurdle = Float64[]
    
    counts_dense = counts isa SparseMatrixCSC ? Matrix(counts) : counts

    for g in 1:n_genes
        y = Float64.(counts_dense[g, :])
        d_i = Float64.(y .> 0)
        
        n_pos = sum(d_i)
        n_neg = n_cells - n_pos
        
        if n_pos < min_cells || n_neg < min_cells
            push!(results, HurdleDEResult(genes[g], 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, provenance_record("HurdleDEResult", "differentialexpression")))
            push!(pvalues_hurdle, 1.0)
            continue
        end
        
        # Discrete part
        beta_det_full, loglik_det_full = fit_logistic(X_full, d_i)
        _, loglik_det_red = fit_logistic(X_reduced, d_i)
        lrt_det = 2.0 * (loglik_det_full - loglik_det_red)
        pvalue_det = ccdf(Chisq(df_df), max(0.0, lrt_det))
        
        # Continuous part
        I_pos = findall(y .> 0)
        z = log2.(y .+ 1.0)
        
        beta_cont_full, loglik_cont_full = fit_linear(X_full[I_pos, :], z[I_pos])
        _, loglik_cont_red = fit_linear(X_reduced[I_pos, :], z[I_pos])
        lrt_cont = 2.0 * (loglik_cont_full - loglik_cont_red)
        pvalue_cont = ccdf(Chisq(df_df), max(0.0, lrt_cont))
        
        # Joint hurdle
        lrt_hurdle = max(0.0, lrt_det) + max(0.0, lrt_cont)
        pvalue_hurdle_val = ccdf(Chisq(2 * df_df), lrt_hurdle)
        
        log2_fc_det = length(beta_det_full) >= 2 ? beta_det_full[2] : 0.0
        log2_fc_cont = length(beta_cont_full) >= 2 ? beta_cont_full[2] : 0.0
        
        push!(results, HurdleDEResult(
            genes[g],
            log2_fc_cont,
            log2_fc_det,
            pvalue_cont,
            pvalue_det,
            pvalue_hurdle_val,
            1.0,
            provenance_record("HurdleDEResult", "differentialexpression")
        ))
        push!(pvalues_hurdle, pvalue_hurdle_val)
    end
    
    padj = benjamini_hochberg(pvalues_hurdle)
    for g in 1:n_genes
        results[g] = HurdleDEResult(
            results[g].gene_id,
            results[g].log2_fc_continuous,
            results[g].log2_fc_detection,
            results[g].pvalue_cont,
            results[g].pvalue_det,
            results[g].pvalue_hurdle,
            padj[g],
            results[g].provenance
        )
    end
    
    return provenance_result!(_ctx, results, "mast_hurdle_test"; parents=String[])
end

function mast_hurdle_test(cm::CountMatrix, design::AbstractVector; kwargs...)
    return mast_hurdle_test(cm.counts, design; gene_ids=cm.gene_ids, kwargs...)
end

function pseudobulk_by_donor_celltype(sce; donor_col::Symbol, celltype_col::Symbol)
    _ctx = active_provenance_context()
    donor_vec = sce.metadata[string(donor_col)]
    celltype_vec = sce.metadata[string(celltype_col)]
    
    n_cells = length(sce.cell_ids)
    n_genes = length(sce.gene_ids)
    
    groups = Tuple{String, String}[]
    for i in 1:n_cells
        d = string(donor_vec[i])
        c = string(celltype_vec[i])
        push!(groups, (d, c))
    end
    unique_groups = sort!(unique(groups))
    
    group_map = Dict{Tuple{String,String}, Vector{Int}}()
    for (i, g) in enumerate(groups)
        push!(get!(group_map, g, Int[]), i)
    end

    n_groups = length(unique_groups)
    aggregated = zeros(Int, n_genes, n_groups)
    sample_ids = String[]
    
    for (idx, (d, c)) in enumerate(unique_groups)
        push!(sample_ids, "$(d)_$(c)")
        cell_indices = get(group_map, (d, c), Int[])
        if !isempty(cell_indices)
            aggregated[:, idx] .= vec(sum(sce.counts[:, cell_indices], dims=2))
        end
    end
    
    result = CountMatrix(sparse(aggregated), copy(sce.gene_ids), sample_ids)
    return provenance_result!(_ctx, result, "pseudobulk_by_donor_celltype"; parents=String[])
end

function pseudobulk_de(sce; group_by=[:donor, :celltype], contrast::Symbol, method::Symbol=:deseq2, min_total::Real=0)
    _ctx = active_provenance_context()
    donor_col = group_by[1]
    celltype_col = group_by[2]
    
    cm = pseudobulk_by_donor_celltype(sce; donor_col=Symbol(donor_col), celltype_col=Symbol(celltype_col))
    
    n_samples = length(cm.sample_ids)
    design_vec = Vector{Symbol}(undef, n_samples)
    
    donor_vec = sce.metadata[string(donor_col)]
    celltype_vec = sce.metadata[string(celltype_col)]
    contrast_vec = sce.metadata[string(contrast)]
    
    for idx in 1:n_samples
        sample_id = cm.sample_ids[idx]
        cell_idx = findfirst(i -> "$(donor_vec[i])_$(celltype_vec[i])" == sample_id, 1:length(sce.cell_ids))
        if cell_idx !== nothing
            design_vec[idx] = Symbol(contrast_vec[cell_idx])
        else
            design_vec[idx] = :unknown
        end
    end
    
    if method == :deseq2
        de_results = differential_expression(cm, design_vec; min_total=min_total)
    elseif method == :hurdle
        de_results = mast_hurdle_test(cm, design_vec)
    else
        throw(ArgumentError("unknown method: $method"))
    end
    
    return provenance_result!(_ctx, de_results, "pseudobulk_de"; parents=String[])
end

function _mast_markers(experiment, labels::AbstractVector{<:Integer}, ident_1::Integer, ident_2::Union{Nothing,Integer}; min_cells::Int=5)
    if ident_2 === nothing
        group1 = [index for (index, label) in enumerate(labels) if label == ident_1]
        group2 = [index for (index, label) in enumerate(labels) if label != ident_1]
    else
        group1 = [index for (index, label) in enumerate(labels) if label == ident_1]
        group2 = [index for (index, label) in enumerate(labels) if label == ident_2]
    end
    
    isempty(group1) && return HurdleDEResult[]
    isempty(group2) && return HurdleDEResult[]
    
    subset = vcat(group1, group2)
    subset_matrix = experiment.counts[:, subset]
    
    design = Symbol[i <= length(group1) ? :ident_1 : :ident_2 for i in 1:length(subset)]
    
    return mast_hurdle_test(subset_matrix, design; gene_ids=experiment.gene_ids, min_cells=min_cells)
end
