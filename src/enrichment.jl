module Enrichment

using Statistics
using LinearAlgebra
using Distributions
using DataFrames
using JSON
using Random

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!

export IDMapper, EnrichmentTerm, EnrichmentDatabase, EnrichmentResult
export load_annotation_database, save_annotation_database, build_annotation_database
export builtin_annotation_database, builtin_annotation_terms
export map_id, map_ids, enrichment_test, go_enrichment, kegg_enrichment, dotplot
export GSEAResult, fgsea_like, gsea_preranked
export GSVAResult, gsva_score
export network_propagation, heat_diffusion_enrichment
export reactome_enrichment, wikipathways_enrichment, msigdb_enrichment
export gene_set_overlap_matrix, jaccard_similarity_matrix
export leading_edge_genes, enrichment_map
export competitive_gene_set_test, self_contained_gene_set_test
export enrichment_heatmap_data, bubble_chart_data
export gene_ontology_semantic_similarity, go_slim_mapping
export term_to_gene_matrix, rank_genes_by_set_membership

"""
    IDMapper

Bidirectional identifier mapper used to normalize query and annotation IDs.
"""
struct IDMapper
    forward::Dict{String,String}
    reverse::Dict{String,Vector{String}}
end

IDMapper() = IDMapper(Dict{String,String}(), Dict{String,Vector{String}}())

"""
    EnrichmentTerm

Annotation term with identifier, name, namespace, gene membership, and parents.
"""
struct EnrichmentTerm
    id::String
    name::String
    namespace::String
    genes::Vector{String}
    parents::Vector{String}
end

"""
    EnrichmentDatabase

Collection of annotation terms together with the identifier mapper used by them.
"""
mutable struct EnrichmentDatabase
    terms::Dict{String,EnrichmentTerm}
    mapper::IDMapper
end

"""
    EnrichmentResult

Single enrichment test result with overlap statistics and adjusted significance.
"""
struct EnrichmentResult <: AbstractAnalysisResult
    term_id::String
    term_name::String
    namespace::String
    overlap::Int
    term_size::Int
    query_size::Int
    background_size::Int
    pvalue::Float64
    padj::Float64
    odds_ratio::Float64
    genes::Vector{String}
    provenance::ResultProvenance
end

EnrichmentResult(term_id, term_name, namespace, overlap, term_size, query_size, background_size, pvalue, padj, odds_ratio, genes) =
    EnrichmentResult(term_id, term_name, namespace, overlap, term_size, query_size, background_size, pvalue, padj, odds_ratio, genes, provenance_record("EnrichmentResult", "enrichment"))

"""
    Base.show(io::IO, result::EnrichmentResult)

Print a compact one-line summary for an enrichment result.
"""
function Base.show(io::IO, result::EnrichmentResult)
    print(io, analysis_result_summary(result))
end

"""
    map_id(mapper::IDMapper, identifier::String)

Map a single identifier through the database mapper, falling back to the input.
"""
function map_id(mapper::IDMapper, identifier::String)

    return get(mapper.forward, String(identifier), String(identifier))
end

"""
    map_ids(mapper::IDMapper, identifiers)

Map a collection of identifiers and return unique normalized IDs.
"""
function map_ids(mapper::IDMapper, identifiers)
    mapped = String[map_id(mapper, identifier) for identifier in identifiers]

    return unique(mapped)
end

"""
    build_annotation_database(terms; mapper=IDMapper())

Build an enrichment database from a vector of annotation terms.
"""
function build_annotation_database(terms::AbstractVector{<:EnrichmentTerm}; mapper::IDMapper=IDMapper(), prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    result = EnrichmentDatabase(Dict(term.id => term for term in terms), mapper)
    return _register_enrichment_result!(_ctx, result, "build_annotation_database"; parents=provenance_parent_ids(terms), parameters=(term_count=length(terms), mapped_count=length(mapper.forward)))
end

"""
    builtin_annotation_terms()

Return the built-in example annotation term set bundled with the module.
"""
function builtin_annotation_terms()

    return EnrichmentTerm[
        EnrichmentTerm("GO:0008150", "biological_process", "GO", ["TP53", "CDKN1A", "MDM2", "BRCA1", "BRCA2", "RAD51", "CDK1", "CDK2", "CCNB1", "CDC20", "XRCC1", "PARP1", "POLB"], String[]),
        EnrichmentTerm("GO:0006974", "response to DNA damage stimulus", "GO", ["TP53", "BRCA1", "BRCA2", "RAD51", "XRCC1", "PARP1"], ["GO:0008150"]),
        EnrichmentTerm("GO:0006281", "DNA repair", "GO", ["BRCA1", "BRCA2", "RAD51", "XRCC1", "PARP1"], ["GO:0006974"]),
        EnrichmentTerm("GO:0006302", "double-strand break repair", "GO", ["BRCA1", "BRCA2", "RAD51"], ["GO:0006281"]),
        EnrichmentTerm("GO:0007049", "cell cycle", "GO", ["CDK1", "CDK2", "CCNB1", "CDC20", "CDKN1A", "MDM2"], ["GO:0008150"]),
        EnrichmentTerm("GO:0000082", "G1/S transition of mitotic cell cycle", "GO", ["CDK2", "CDKN1A", "MDM2"], ["GO:0007049"]),
        EnrichmentTerm("KEGG:hsa04110", "Cell cycle", "KEGG", ["CDK1", "CDK2", "CCNB1", "CDC20", "CDKN1A"], ["KEGG:hsa05200"]),
        EnrichmentTerm("KEGG:hsa03410", "Base excision repair", "KEGG", ["TP53", "XRCC1", "PARP1", "POLB"], ["KEGG:hsa05200"]),
        EnrichmentTerm("KEGG:hsa05200", "Pathways in cancer", "KEGG", ["TP53", "MDM2", "BRCA1", "CDKN1A", "CDK1", "CDK2", "CCNB1", "CDC20", "XRCC1", "PARP1", "POLB"], ["KEGG:hsa00001"]),
        EnrichmentTerm("KEGG:hsa00001", "reference pathways", "KEGG", ["TP53", "MDM2", "BRCA1", "CDKN1A", "CDK1", "CDK2", "XRCC1", "PARP1", "POLB", "CCNB1", "CDC20"], String[]),
    ]
end

"""
    builtin_annotation_database()

Build an enrichment database from the built-in annotation terms.
"""
function builtin_annotation_database()
    _ctx = active_provenance_context()


    return build_annotation_database(builtin_annotation_terms(); _ctx=_ctx)
end

"""
    save_annotation_database(path, database)

Write an enrichment database to a JSON file.
"""
function save_annotation_database(path::String, database::EnrichmentDatabase)
    payload = Dict{String,Any}(
        "terms" => [Dict(
            "id" => term.id,
            "name" => term.name,
            "namespace" => term.namespace,
            "genes" => term.genes,
            "parents" => term.parents) for term in values(database.terms)],
        "mapper" => Dict(
            "forward" => database.mapper.forward,
            "reverse" => database.mapper.reverse))
    open(path, "w") do io
        JSON.print(io, payload)
    end
    _ctx = active_provenance_context()


    return _register_enrichment_result!(_ctx, path, "save_annotation_database"; parents=provenance_parent_ids(database), parameters=(path=path, term_count=length(database.terms)))
end

"""
    load_annotation_database(path)

Load an enrichment database from a JSON file.
"""
function load_annotation_database(path::String)
    payload = JSON.parsefile(path)
    terms = EnrichmentTerm[]
    for item in payload["terms"]
        push!(terms, EnrichmentTerm(String(item["id"]), String(item["name"]), String(item["namespace"]), String.(item["genes"]), String.(item["parents"])))
    end
    mapper_payload = get(payload, "mapper", Dict{String,Any}())
    forward = Dict{String,String}(String(key) => String(value) for (key, value) in get(mapper_payload, "forward", Dict{String,Any}()))
    reverse = Dict{String,Vector{String}}(String(key) => String.(value) for (key, value) in get(mapper_payload, "reverse", Dict{String,Any}()))
    result = EnrichmentDatabase(Dict(term.id => term for term in terms), IDMapper(forward, reverse))
    _ctx = active_provenance_context()


    return _register_enrichment_result!(_ctx, result, "load_annotation_database"; parents=String[], parameters=(path=path, term_count=length(terms)))
end

function _term_genes(term::EnrichmentTerm, universe::Set{String})
    return Set(gene for gene in term.genes if gene in universe)
end

function _background_set(database::EnrichmentDatabase)
    background = Set{String}()
    for term in values(database.terms)
        union!(background, term.genes)
    end
    return background
end

function _all_descendants(database::EnrichmentDatabase, term_id::String, cache::Dict{String,Set{String}})
    haskey(cache, term_id) && return cache[term_id]
    descendants = Set{String}()
    for child in values(database.terms)
        term_id in child.parents || continue
        push!(descendants, child.id)
        union!(descendants, _all_descendants(database, child.id, cache))
    end
    cache[term_id] = descendants
    return descendants
end

function _effective_gene_sets(database::EnrichmentDatabase, universe::Set{String}; elim::Bool=false, significant_terms::Set{String}=Set{String}())
    gene_sets = Dict{String,Set{String}}()
    for (term_id, term) in database.terms
        gene_sets[term_id] = _term_genes(term, universe)
    end
    elim || return gene_sets
    descendants_cache = Dict{String,Set{String}}()
    for term_id in significant_terms
        term = get(database.terms, term_id, nothing)
        term === nothing && continue
        descendants = _all_descendants(database, term_id, descendants_cache)
        for parent in term.parents
            haskey(gene_sets, parent) || continue
            setdiff!(gene_sets[parent], gene_sets[term_id])
        end
        for descendant in descendants
            haskey(gene_sets, descendant) || continue
            setdiff!(gene_sets[descendant], gene_sets[term_id])
        end
    end
    return gene_sets
end

function _fisher_right_tail(overlap::Int, term_size::Int, query_size::Int, background_size::Int)
    (term_size <= 0 || query_size <= 0 || background_size <= 0) && return 1.0
    overlap <= 0 && return 1.0
    successes = clamp(term_size, 0, background_size)
    draws = clamp(query_size, 0, background_size)
    distribution = Hypergeometric(successes, background_size - successes, draws)
    return ccdf(distribution, overlap - 1)
end

function _odds_ratio(overlap::Int, term_size::Int, query_size::Int, background_size::Int)
    a = overlap + 0.5
    b = query_size - overlap + 0.5
    c = term_size - overlap + 0.5
    d = background_size - term_size - query_size + overlap + 0.5
    return (a * d) / (b * c)
end

function _benjamini_hochberg(pvalues::AbstractVector{<:Real})
    n = length(pvalues)
    n == 0 && return Float64[]
    sanitized = [clamp(isfinite(Float64(p)) ? Float64(p) : 1.0, 0.0, 1.0) for p in pvalues]
    order = sortperm(sanitized)
    sorted = sanitized[order]
    adjusted_sorted = similar(sorted)
    running_min = 1.0
    for index in n:-1:1
        running_min = min(running_min, sorted[index] * n / index)
        adjusted_sorted[index] = min(running_min, 1.0)
    end
    adjusted = similar(sanitized)
    for (rank, original_index) in enumerate(order)
        adjusted[original_index] = adjusted_sorted[rank]
    end
    return adjusted
end

function _result_sort_key(result::EnrichmentResult)
    return (result.padj, result.pvalue, -result.overlap, result.term_id)
end

"""
    enrichment_test(query_genes, database; background=nothing, namespace=nothing, elim=false, threaded=true, min_overlap=1)

Run an over-representation enrichment test for a set of query genes.
"""
function enrichment_test(query_genes, database::EnrichmentDatabase; background=nothing, namespace::Union{Nothing,String}=nothing, elim::Bool=false, threaded::Bool=true, min_overlap::Int=1, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    query = Set(String.(query_genes))
    universe = background === nothing ? _background_set(database) : Set(String.(background))
    query = intersect(query, universe)
    gene_sets = _effective_gene_sets(database, universe; elim=elim)

    term_ids = String[]
    for (term_id, term) in database.terms
        namespace !== nothing && term.namespace != String(namespace) && continue
        push!(term_ids, term_id)
    end

    results = Vector{Union{Nothing,EnrichmentResult}}(undef, length(term_ids))
    worker = index -> begin
        term_id = term_ids[index]
        term = database.terms[term_id]
        genes = gene_sets[term_id]
        overlap_genes = sort!(collect(intersect(query, genes)))
        overlap = length(overlap_genes)
        if overlap < min_overlap || isempty(genes)
            results[index] = nothing
            return
        end
        pvalue = _fisher_right_tail(overlap, length(genes), length(query), length(universe))
        odds_ratio = _odds_ratio(overlap, length(genes), length(query), length(universe))
        results[index] = EnrichmentResult(term.id, term.name, term.namespace, overlap, length(genes), length(query), length(universe), pvalue, 1.0, odds_ratio, overlap_genes)
    end

    if threaded && length(term_ids) > 1 && Base.Threads.nthreads() > 1
        Base.Threads.@threads for index in eachindex(term_ids)
            worker(index)
        end
    else
        for index in eachindex(term_ids)
            worker(index)
        end
    end

    filtered = EnrichmentResult[result for result in results if result !== nothing]
    padj = _benjamini_hochberg([result.pvalue for result in filtered])
    scored = EnrichmentResult[
        EnrichmentResult(result.term_id, result.term_name, result.namespace, result.overlap, result.term_size, result.query_size, result.background_size, result.pvalue, padj[index], result.odds_ratio, result.genes)
        for (index, result) in enumerate(filtered)
    ]
    sort!(scored; by=_result_sort_key)


    return _register_enrichment_result!(_ctx, scored, "enrichment_test"; parents=provenance_parent_ids(query_genes, database, background), parameters=(query_count=length(query), universe_count=length(universe), term_count=length(term_ids), namespace=namespace === nothing ? "all" : String(namespace), elim=elim, min_overlap=min_overlap))
end

"""
    go_enrichment(query_genes, database; kwargs...)

Convenience wrapper around [`enrichment_test`](@ref) restricted to the GO namespace.
"""
function go_enrichment(query_genes, database::EnrichmentDatabase, kwargs...)
    _ctx = active_provenance_context()


    return enrichment_test(query_genes, database; namespace="GO", _ctx=_ctx, kwargs...)
end

"""
    kegg_enrichment(query_genes, database; kwargs...)

Convenience wrapper around [`enrichment_test`](@ref) restricted to the KEGG namespace.
"""
function kegg_enrichment(query_genes, database::EnrichmentDatabase, kwargs...)
    _ctx = active_provenance_context()


    return enrichment_test(query_genes, database; namespace="KEGG", _ctx=_ctx, kwargs...)
end

function _plot_available()
    return isdefined(Main, :Makie)
end

function _with_makie(f::Function)
    _plot_available() || return nothing
    return f(getfield(Main, :Makie))
end

"""
    dotplot(results; top_n=20, save_path=nothing)

Create a compact enrichment dot plot when Makie is available.
"""
function dotplot(results::AbstractVector{<:EnrichmentResult}; top_n::Int=20, save_path::Union{Nothing,String}=nothing)
    selected = first(sort(collect(results); by=_result_sort_key), min(top_n, length(results)))
    figure = _with_makie() do Makie
        fig = Makie.Figure()
        axis = Makie.Axis(fig[1, 1]; title="Enrichment dot plot", xlabel="gene ratio", ylabel="term")
        y_positions = collect(1:length(selected))
        ratios = [result.overlap / max(result.term_size, 1) for result in selected]
        magnitudes = [-log10(max(result.padj, eps(Float64))) for result in selected]
        labels = [result.term_name for result in selected]
        Makie.scatter!(axis, ratios, y_positions; markersize=10 .+ 4 .* magnitudes, color=magnitudes, colormap=:viridis)
        axis.yticks = (y_positions, labels)
        axis.xlabel = "Gene ratio"
        axis.ylabel = "Term"
        return fig
    end
    figure !== nothing && save_path !== nothing && _with_makie(Makie -> Makie.save(save_path, figure))

    return (results=selected, figure=figure)
end

# =============================================================================
# NEW: GSEA / fgsea-like preranked gene set enrichment
# =============================================================================

"""
    GSEAResult

Result of a preranked GSEA analysis for a single gene set.
"""
struct GSEAResult <: AbstractAnalysisResult
    term_id::String
    term_name::String
    enrichment_score::Float64    # ES (signed)
    normalized_es::Float64       # NES
    pvalue::Float64
    padj::Float64
    leading_edge::Vector{String}
    n_genes::Int
    correlation_sum::Float64     # sum of absolute rankings in leading edge
    provenance::ResultProvenance
end

GSEAResult(term_id, term_name, enrichment_score, normalized_es, pvalue, padj, leading_edge, n_genes, correlation_sum) =
    GSEAResult(term_id, term_name, enrichment_score, normalized_es, pvalue, padj, leading_edge, n_genes, correlation_sum, provenance_record("GSEAResult", "enrichment"))

"""
    fgsea_like(ranked_genes, gene_sets; n_permutations=1000, min_size=15, max_size=500, seed=1) → Vector{GSEAResult}

Preranked gene set enrichment analogous to `fgsea` (Korotkevich et al. 2021).

`ranked_genes`: vector of `(gene_id, rank_score)` pairs, sorted by score descending.
`gene_sets`: `Dict{String, Vector{String}}` mapping set name → gene IDs.

Returns a `Vector{GSEAResult}` sorted by NES.
"""
function fgsea_like(
    ranked_genes::AbstractVector,   # [(gene, score), ...]
    gene_sets::Dict{String,<:AbstractVector{<:AbstractString}};
    n_permutations::Int=1000,
    min_size::Int=15,
    max_size::Int=500,
    seed::Int=1,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))
    rng = Random.MersenneTwister(seed)
    genes = [String(first(r)) for r in ranked_genes]
    scores = Float64[Float64(last(r)) for r in ranked_genes]
    N = length(genes)
    gene_rank = Dict(g => i for (i, g) in enumerate(genes))

    results = GSEAResult[]

    for (set_name, set_genes) in gene_sets
        set_g = filter(g -> haskey(gene_rank, String(g)), String.(set_genes))
        k = length(set_g)
        (k < min_size || k > max_size) && continue

        # ---- Observed ES ----
        es, le = _compute_es(scores, gene_rank, set_g, N)

        # ---- Permutation NES ----
        null_es = Float64[]
        for _ in 1:n_permutations
            perm_set = randperm(rng, N)[1:k]
            perm_genes = genes[perm_set]
            es_null, _ = _compute_es(scores, gene_rank, perm_genes, N)
            push!(null_es, es_null)
        end

        pos_null = filter(>(0), null_es)
        neg_null = filter(<(0), null_es)

        nes = if es >= 0
            pos_mean = isempty(pos_null) ? eps() : mean(pos_null)
            es / max(pos_mean, eps())
        else
            neg_mean = isempty(neg_null) ? -eps() : mean(neg_null)
            -es / max(-neg_mean, eps())
        end

        pval = if es >= 0
            count(x -> x >= es, null_es) / max(n_permutations, 1)
        else
            count(x -> x <= es, null_es) / max(n_permutations, 1)
        end
        pval = clamp(pval, 1.0 / n_permutations, 1.0)

        push!(results, GSEAResult(set_name, set_name, es, nes, pval, pval, le, k, sum(abs.(scores[gene_rank[g]] for g in le))))
    end

    # BH FDR on pvalues
    sort!(results, by=r -> r.pvalue)
    n = length(results)
    padj_vals = [min(results[i].pvalue * n / i, 1.0) for i in 1:n]
    for i in (n-1):-1:1
        padj_vals[i] = min(padj_vals[i], padj_vals[i+1])
    end

    final = [GSEAResult(r.term_id, r.term_name, r.enrichment_score, r.normalized_es,
        r.pvalue, padj_vals[i], r.leading_edge, r.n_genes, r.correlation_sum)
             for (i, r) in enumerate(results)]
    sort!(final, by=r -> -r.normalized_es)
    return _register_enrichment_result!(_ctx, final, "fgsea_like"; parents=provenance_parent_ids(ranked_genes, gene_sets), parameters=(permutations=n_permutations, min_size=min_size, max_size=max_size, result_count=length(final)))
end

"""
    gsea_preranked(ranked_genes, database; kwargs...) → Vector{GSEAResult}

Convenience wrapper that extracts gene sets from an `EnrichmentDatabase`.
"""
function gsea_preranked(ranked_genes, database::EnrichmentDatabase; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx), kwargs...)
    gene_sets = Dict(id => term.genes for (id, term) in database.terms)


    return fgsea_like(ranked_genes, gene_sets; _ctx=_ctx, kwargs...)
end

function _compute_es(scores::Vector{Float64}, gene_rank::Dict{String,Int}, set_genes::Vector{String}, N::Int)
    k = length(set_genes)
    k == 0 && return 0.0, String[]

    is_set = falses(N)
    for g in set_genes
        r = get(gene_rank, g, 0)
        r > 0 && (is_set[r] = true)
    end

    # Weights: |rank score| for set members, uniform for non-members
    abs_scores = abs.(scores)
    set_weight = sum(abs_scores[i] for i in 1:N if is_set[i]; init=eps())
    non_set_k = max(N - k, 1)

    cum_es = 0.0
    max_es = 0.0
    min_es = 0.0
    last_peak_i = 1

    for i in 1:N
        if is_set[i]
            cum_es += abs_scores[i] / set_weight
        else
            cum_es -= 1.0 / non_set_k
        end
        cum_es > max_es && (max_es = cum_es; last_peak_i = i)
        cum_es < min_es && (min_es = cum_es)
    end

    es = abs(max_es) > abs(min_es) ? max_es : min_es
    le = [g for g in set_genes if get(gene_rank, g, N + 1) <= last_peak_i]
    return es, le
end

# =============================================================================
# NEW: GSVA
# =============================================================================

"""
    GSVAResult

Gene Set Variation Analysis result for all samples and gene sets.
"""
struct GSVAResult <: AbstractAnalysisResult
    scores::Matrix{Float64}      # gene_sets × samples
    gene_set_names::Vector{String}
    sample_names::Vector{String}
    provenance::ResultProvenance
end

GSVAResult(scores, gene_set_names, sample_names) =
    GSVAResult(scores, gene_set_names, sample_names, provenance_record("GSVAResult", "enrichment"))

@inline function _register_enrichment_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

"""
    gsva_score(expression_matrix, gene_sets; method=:gsva, kcdf=:gaussian, sample_names=nothing) → GSVAResult

Compute per-sample gene set variation scores (Hänzelmann et al. 2013).

`expression_matrix`: genes × samples.
`gene_sets`: `Dict{String, Vector{String}}` (set name → gene IDs).
`gene_names`: names corresponding to rows of the matrix.

Methods: `:gsva` (default), `:ssgsea` (Barbie 2009), `:zscore`.
"""
function gsva_score(
    expression_matrix::AbstractMatrix{<:Real},
    gene_sets::Dict{String,<:AbstractVector{<:AbstractString}},
    gene_names::AbstractVector{<:AbstractString};
    method::Symbol=:gsva,
    sample_names::Union{Nothing,AbstractVector{<:AbstractString}}=nothing,
    min_size::Int=5,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))
    X = Matrix{Float64}(expression_matrix)
    n_genes, n_samples = size(X)
    length(gene_names) == n_genes || throw(DimensionMismatch("gene_names must match rows"))
    gene_idx = Dict(String(g) => i for (i, g) in enumerate(gene_names))

    sn = sample_names !== nothing ? String.(sample_names) : ["S$j" for j in 1:n_samples]
    set_names = String[]
    set_indices = Vector{Vector{Int}}()

    for (sn_s, sg) in gene_sets
        idx = filter(>(0), [get(gene_idx, String(g), 0) for g in sg])
        length(idx) < min_size && continue
        push!(set_names, sn_s)
        push!(set_indices, idx)
    end

    n_sets = length(set_names)
    scores = zeros(Float64, n_sets, n_samples)

    Threads.@threads for j in 1:n_samples
        col = X[:, j]
        # Rank-normalise column
        ranked = _rank_normalise(col)
        for (si, idx) in enumerate(set_indices)
            if method == :ssgsea
                scores[si, j] = _ssgsea_sample_score(ranked, idx, n_genes)
            elseif method == :zscore
                in_set = mean(col[idx])
                out_set = mean(col[setdiff(1:n_genes, idx)])
                σ = std(col)
                scores[si, j] = σ > 0 ? (in_set - out_set) / σ : 0.0
            else   # :gsva
                scores[si, j] = _gsva_sample_score(ranked, idx, n_genes)
            end
        end
    end

    result = GSVAResult(scores, set_names, sn)


    return _register_enrichment_result!(_ctx, result, "gsva_score"; parents=provenance_parent_ids(expression_matrix, gene_sets, gene_names), parameters=(method=method, set_count=n_sets, sample_count=n_samples, min_size=min_size))
end

function _rank_normalise(v::Vector{Float64})
    n = length(v)
    ord = sortperm(v)
    rnk = zeros(Float64, n)
    for (r, i) in enumerate(ord)
        rnk[i] = r / n
    end
    return rnk
end

function _ssgsea_sample_score(ranked::Vector{Float64}, set_idx::Vector{Int}, N::Int)
    is_set = falses(N)
    for i in set_idx
        is_set[i] = true
    end
    kappa = 0.25
    cum = 0.0
    for i in 1:N
        if is_set[i]
            cum += ranked[i]^kappa / sum(ranked[set_idx] .^ kappa)
        else
            cum -= 1.0 / (N - length(set_idx))
        end
    end
    return cum
end

function _gsva_sample_score(ranked::Vector{Float64}, set_idx::Vector{Int}, N::Int)
    k = length(set_idx)
    is_set = falses(N)
    for i in set_idx
        is_set[i] = true
    end
    max_pos = -Inf
    max_neg = -Inf
    cum = 0.0
    for i in 1:N
        if is_set[i]
            cum += ranked[i] / sum(ranked[set_idx])
        else
            cum -= 1.0 / (N - k)
        end
        cum > max_pos && (max_pos = cum)
        -cum > max_neg && (max_neg = -cum)
    end
    return max_pos - max_neg
end

# =============================================================================
# NEW: Network propagation
# =============================================================================

"""
    network_propagation(seed_genes, adjacency_matrix, gene_names; alpha=0.5, n_iter=30) → DataFrame

Propagate gene scores through a protein-protein interaction network using
Random Walk with Restart (RWR), analogous to HotNet2 / network propagation methods.

`seed_genes`: `Dict{String, Float64}` of seed gene → initial score.
`adjacency_matrix`: symmetric, non-negative.
`gene_names`: labels for rows/columns.

Returns a DataFrame with propagated scores and ranks.
"""
function network_propagation(
    seed_genes::Dict{String,<:Real},
    adjacency_matrix::AbstractMatrix{<:Real},
    gene_names::AbstractVector{<:AbstractString};
    alpha::Real=0.5,
    n_iter::Int=30,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))
    n = length(gene_names)
    size(adjacency_matrix) == (n, n) || throw(DimensionMismatch("adjacency must be n × n"))
    gene_idx = Dict(String(g) => i for (i, g) in enumerate(gene_names))

    # Row-normalise adjacency (transition matrix)
    A = Matrix{Float64}(adjacency_matrix)
    row_sums = max.(vec(sum(A, dims=2)), eps(Float64))
    W = A ./ row_sums     # column-stochastic when transposed

    # Initial score vector
    y = zeros(Float64, n)
    for (g, s) in seed_genes
        i = get(gene_idx, String(g), 0)
        i > 0 && (y[i] = Float64(s))
    end
    y ./= max(sum(y), eps(Float64))

    # RWR: f = α * W' * f + (1-α) * y
    f = copy(y)
    α = Float64(alpha)
    Wt = permutedims(W)
    for _ in 1:n_iter
        f_new = α .* (Wt * f) .+ (1 - α) .* y
        maximum(abs.(f_new .- f)) < 1e-8 && (f = f_new; break)
        f = f_new
    end

    ord = sortperm(f, rev=true)
    result = DataFrame(
        gene=gene_names[ord],
        propagated_score=f[ord],
        rank=1:n,
        is_seed=[haskey(seed_genes, gene_names[i]) for i in ord])


    return _register_enrichment_result!(_ctx, result, "network_propagation"; parents=provenance_parent_ids(seed_genes, adjacency_matrix, gene_names), parameters=(alpha=float(alpha), n_iter=n_iter, gene_count=n))
end

"""
    heat_diffusion_enrichment(seed_scores, adjacency, gene_names, gene_sets; t=0.5) → DataFrame

Run network propagation then test enrichment of gene sets in the top-ranked genes.
"""
function heat_diffusion_enrichment(
    seed_genes::Dict{String,<:Real},
    adjacency_matrix::AbstractMatrix{<:Real},
    gene_names::AbstractVector{<:AbstractString},
    gene_sets::Dict{String,<:AbstractVector{<:AbstractString}};
    t::Real=0.5, top_k::Int=200, kwargs...
)
    prop = network_propagation(seed_genes, adjacency_matrix, gene_names; alpha=t, _ctx=_ctx)
    top_genes = prop.gene[1:min(top_k, nrow(prop))]
    # Build a simple enrichment database from the provided gene sets
    terms = [EnrichmentTerm(id, id, "NetProp", collect(String.(genes)), String[])
             for (id, genes) in gene_sets]
    db = build_annotation_database(terms; _ctx=_ctx)
    _ctx = active_provenance_context()


    return enrichment_test(top_genes, db; _ctx=_ctx, kwargs...)
end

# =============================================================================
# NEW: Convenience enrichment wrappers (Reactome / WikiPathways / MSigDB)
# =============================================================================

"""
    reactome_enrichment(query_genes, database; kwargs...) → Vector{EnrichmentResult}

Enrichment restricted to Reactome-namespace terms.
"""
reactome_enrichment(query_genes, db::EnrichmentDatabase, kwargs...) =
    enrichment_test(query_genes, db; namespace="Reactome", _ctx=_ctx, kwargs...)

"""
    wikipathways_enrichment(query_genes, database; kwargs...) → Vector{EnrichmentResult}
"""
wikipathways_enrichment(query_genes, db::EnrichmentDatabase, kwargs...) =
    enrichment_test(query_genes, db; namespace="WikiPathways", _ctx=_ctx, kwargs...)

"""
    msigdb_enrichment(query_genes, database; collection=nothing, kwargs...) → Vector{EnrichmentResult}

MSigDB-style enrichment. Filter by namespace collection name (e.g. "H", "C2", "C5").
"""
function msigdb_enrichment(query_genes, db::EnrichmentDatabase; collection::Union{Nothing,String}=nothing, kwargs...)
    if collection !== nothing
        return enrichment_test(query_genes, db; namespace=collection, _ctx=_ctx, kwargs...)
    end
    _ctx = active_provenance_context()


    return enrichment_test(query_genes, db; _ctx=_ctx, kwargs...)
end

# =============================================================================
# NEW: Gene set similarity / overlap
# =============================================================================

"""
    gene_set_overlap_matrix(database; namespace=nothing) → (Matrix, labels)

Compute pairwise gene set overlap (Jaccard index) for all terms in a database.
Useful for building enrichment maps.
"""
function gene_set_overlap_matrix(database::EnrichmentDatabase; namespace::Union{Nothing,String}=nothing)
    terms = collect(values(database.terms))
    namespace !== nothing && filter!(t -> t.namespace == namespace, terms)
    n = length(terms)
    M = zeros(Float64, n, n)
    for i in 1:n, j in i:n
        ji = jaccard_similarity_matrix(terms[i].genes, terms[j].genes)
        M[i, j] = M[j, i] = ji
    end

    return M, [t.name for t in terms]
end

"""
    jaccard_similarity_matrix(a, b) → Float64

Jaccard index between two gene sets.
"""
function jaccard_similarity_matrix(a::AbstractVector{<:AbstractString}, b::AbstractVector{<:AbstractString})
    sa = Set(String.(a))
    sb = Set(String.(b))
    inter = length(intersect(sa, sb))
    union_n = length(union(sa, sb))

    return union_n == 0 ? 0.0 : inter / union_n
end

# =============================================================================
# NEW: Leading edge & enrichment map
# =============================================================================

"""
    leading_edge_genes(gsea_result; min_freq=2) → DataFrame

Extract the core leading edge genes from a GSEA result vector.
Genes appearing in multiple leading edge sets are more likely to be drivers.
"""
function leading_edge_genes(gsea_results::AbstractVector{GSEAResult}; min_freq::Int=2)
    gene_counts = Dict{String,Int}()
    for r in gsea_results
        for g in r.leading_edge
            gene_counts[g] = get(gene_counts, g, 0) + 1
        end
    end
    pairs_sorted = sort(collect(gene_counts), by=kv -> -kv[2])
    filt = filter(kv -> kv[2] >= min_freq, pairs_sorted)

    return DataFrame(gene=first.(filt), n_sets=last.(filt))
end

"""
    enrichment_map(gsea_results; similarity_threshold=0.3) → DataFrame

Build an enrichment map edge list (Merico et al. 2010).
Nodes = significant gene sets; edges = Jaccard similarity between leading edges.
"""
function enrichment_map(gsea_results::AbstractVector{GSEAResult}; similarity_threshold::Real=0.3)
    n = length(gsea_results)
    rows = NamedTuple[]
    for i in 1:n-1, j in (i+1):n
        ji = jaccard_similarity_matrix(gsea_results[i].leading_edge, gsea_results[j].leading_edge)
        ji >= Float64(similarity_threshold) && push!(rows, (
            set1=gsea_results[i].term_name, set2=gsea_results[j].term_name,
            jaccard=ji, nes1=gsea_results[i].normalized_es, nes2=gsea_results[j].normalized_es))
    end

    return DataFrame(rows)
end

"""
    competitive_gene_set_test(ranked_genes, gene_set, gene_names; n_perm=1000, seed=1) → NamedTuple

Competitive test: is the gene set enriched relative to all other genes?
Uses a permutation of gene labels (competitive null).

Returns `(mean_rank_in_set, mean_rank_outside, pvalue, direction)`.
"""
function competitive_gene_set_test(
    ranked_genes::AbstractVector,
    gene_set::AbstractVector{<:AbstractString},
    gene_names::AbstractVector{<:AbstractString};
    n_perm::Int=1000, seed::Int=1)
    rng = Random.MersenneTwister(seed)
    genes = String[String(first(r)) for r in ranked_genes]
    scores = Float64[Float64(last(r)) for r in ranked_genes]
    n = length(genes)
    gene_idx = Dict(g => i for (i, g) in enumerate(genes))
    set_g = filter(g -> haskey(gene_idx, g), String.(gene_set))
    k = length(set_g)
    k == 0 && return (mean_rank_in_set=0.0, mean_rank_outside=0.0, pvalue=1.0, direction="none")

    in_idx = [gene_idx[g] for g in set_g]
    out_idx = setdiff(1:n, in_idx)

    obs_mean_in = isempty(in_idx) ? 0.0 : mean(scores[in_idx])
    obs_mean_out = isempty(out_idx) ? 0.0 : mean(scores[out_idx])
    obs_diff = obs_mean_in - obs_mean_out

    null_diffs = Float64[]
    for _ in 1:n_perm
        perm_idx = randperm(rng, n)[1:k]
        perm_in = mean(scores[perm_idx])
        perm_out = mean(scores[setdiff(1:n, perm_idx)])
        push!(null_diffs, perm_in - perm_out)
    end

    pval = obs_diff > 0 ? count(x -> x >= obs_diff, null_diffs) / n_perm :
           count(x -> x <= obs_diff, null_diffs) / n_perm
    pval = clamp(pval, 1 / n_perm, 1.0)

    return (mean_rank_in_set=obs_mean_in, mean_rank_outside=obs_mean_out,
        pvalue=pval, direction=obs_diff > 0 ? "enriched" : "depleted")
end

"""
    self_contained_gene_set_test(expression_matrix, gene_set, gene_names, group1_idx, group2_idx; n_perm=1000) → NamedTuple

Self-contained test: does the gene set show differential expression between groups?
Permutes sample labels (self-contained null).
"""
function self_contained_gene_set_test(
    expression_matrix::AbstractMatrix{<:Real},
    gene_set::AbstractVector{<:AbstractString},
    gene_names::AbstractVector{<:AbstractString},
    group1_idx::AbstractVector{<:Integer},
    group2_idx::AbstractVector{<:Integer};
    n_perm::Int=1000, seed::Int=1)
    rng = Random.MersenneTwister(seed)
    X = Matrix{Float64}(expression_matrix)
    gidx = Dict(String(g) => i for (i, g) in enumerate(gene_names))
    set_idx = filter(>(0), [get(gidx, String(g), 0) for g in gene_set])
    isempty(set_idx) && return (test_statistic=0.0, pvalue=1.0, direction="none")

    Xset = X[set_idx, :]
    score(g1, g2) = mean(mean(Xset[:, g1], dims=2) .- mean(Xset[:, g2], dims=2))
    obs = score(group1_idx, group2_idx)

    all_idx = vcat(collect(group1_idx), collect(group2_idx))
    n1 = length(group1_idx)
    null_scores = Float64[]
    for _ in 1:n_perm
        perm = shuffle!(rng, copy(all_idx))
        push!(null_scores, score(perm[1:n1], perm[n1+1:end]))
    end

    pval = obs > 0 ? count(x -> x >= obs, null_scores) / n_perm :
           count(x -> x <= obs, null_scores) / n_perm
    return (test_statistic=obs, pvalue=clamp(pval, 1 / n_perm, 1.0),
        direction=obs > 0 ? "group1_higher" : "group2_higher")
end

"""
    enrichment_heatmap_data(gsva_result; top_n=30) → NamedTuple

Prepare GSVA score matrix for heatmap rendering.
Returns `(matrix, row_labels, col_labels)`.
"""
function enrichment_heatmap_data(result::GSVAResult; top_n::Int=30)
    n_sets = length(result.gene_set_names)
    selected = 1:min(top_n, n_sets)

    return (
        matrix=result.scores[selected, :],
        row_labels=result.gene_set_names[selected],
        col_labels=result.sample_names)
end

"""
    bubble_chart_data(results; top_n=20) → DataFrame

Prepare an enrichment bubble chart data table compatible with any plotting tool.
"""
function bubble_chart_data(results::AbstractVector{<:EnrichmentResult}; top_n::Int=20)
    sel = first(sort(collect(results); by=_result_sort_key), min(top_n, length(results)))

    return DataFrame(
        term_name=[r.term_name for r in sel],
        gene_ratio=[r.overlap / max(r.term_size, 1) for r in sel],
        neg_log10_padj=[-log10(max(r.padj, eps())) for r in sel],
        overlap_count=[r.overlap for r in sel],
        odds_ratio=[r.odds_ratio for r in sel],
        padj=[r.padj for r in sel],
        namespace=[r.namespace for r in sel])
end

# =============================================================================
# NEW: GO semantic similarity (Wang's graph-based method)
# =============================================================================

"""
    gene_ontology_semantic_similarity(term_id1, term_id2, database; method=:wang) → Float64

Compute semantic similarity between two GO terms using the
Wang et al. (2007) graph-based method.

Approximated here by depth-weighted information content from the database.
"""
function gene_ontology_semantic_similarity(
    term_id1::AbstractString,
    term_id2::AbstractString,
    database::EnrichmentDatabase;
    method::Symbol=:wang)
    t1 = get(database.terms, String(term_id1), nothing)
    t2 = get(database.terms, String(term_id2), nothing)
    (t1 === nothing || t2 === nothing) && return 0.0

    # Gene overlap-based proxy (Lin semantic similarity)
    g1 = Set(t1.genes)
    g2 = Set(t2.genes)
    inter = length(intersect(g1, g2))
    total_bg = sum(length(t.genes) for t in values(database.terms))
    ic1 = inter > 0 ? -log(length(g1) / max(total_bg, 1)) : 0.0
    ic2 = inter > 0 ? -log(length(g2) / max(total_bg, 1)) : 0.0
    ic_lcs = inter > 0 ? -log(inter / max(total_bg, 1)) : 0.0
    denom = ic1 + ic2
    return denom > 0 ? 2 * ic_lcs / denom : 0.0
end

"""
    go_slim_mapping(genes, database; slim_terms=nothing) → DataFrame

Map genes to GO Slim categories (a reduced GO annotation set).
Returns a DataFrame of (gene, slim_term_id, slim_term_name, namespace).
"""
function go_slim_mapping(
    genes::AbstractVector{<:AbstractString},
    database::EnrichmentDatabase;
    slim_terms::Union{Nothing,AbstractVector{<:AbstractString}}=nothing)
    gs = Set(String.(genes))
    terms = slim_terms !== nothing ?
            [database.terms[t] for t in String.(slim_terms) if haskey(database.terms, t)] :
            collect(values(database.terms))

    rows = NamedTuple[]
    for g in gs
        for t in terms
            g in t.genes && push!(rows, (gene=g, term_id=t.id, term_name=t.name, namespace=t.namespace))
        end
    end
    return DataFrame(rows)
end

# =============================================================================
# NEW: Utility — term-to-gene matrix and ranking
# =============================================================================

"""
    term_to_gene_matrix(database, genes; namespace=nothing) → (Matrix{Bool}, term_labels, gene_labels)

Build a binary term × gene membership matrix.
"""
function term_to_gene_matrix(
    database::EnrichmentDatabase,
    genes::AbstractVector{<:AbstractString};
    namespace::Union{Nothing,String}=nothing)
    gs = String.(genes)
    terms = collect(values(database.terms))
    namespace !== nothing && filter!(t -> t.namespace == namespace, terms)
    gene_idx = Dict(g => i for (i, g) in enumerate(gs))
    n_terms, n_genes = length(terms), length(gs)
    M = falses(n_terms, n_genes)
    for (ti, t) in enumerate(terms)
        for g in t.genes
            ji = get(gene_idx, g, 0)
            ji > 0 && (M[ti, ji] = true)
        end
    end
    return M, [t.id for t in terms], gs
end

"""
    rank_genes_by_set_membership(database, genes; min_sets=2) → DataFrame

Rank genes by how many gene sets they participate in — a proxy for
network hub / pathway coverage.
"""
function rank_genes_by_set_membership(
    database::EnrichmentDatabase,
    genes::AbstractVector{<:AbstractString};
    min_sets::Int=1)
    gs = Set(String.(genes))
    counts = Dict{String,Int}()
    for t in values(database.terms)
        for g in t.genes
            g in gs || continue
            counts[g] = get(counts, g, 0) + 1
        end
    end
    pairs_sorted = sort(collect(counts), by=kv -> -kv[2])
    filt = filter(kv -> kv[2] >= min_sets, pairs_sorted)
    return DataFrame(gene=first.(filt), n_sets=last.(filt), rank=1:length(filt))
end

end  # module Enrichment
