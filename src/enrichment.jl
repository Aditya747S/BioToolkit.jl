module Enrichment

using Statistics
using LinearAlgebra
using Distributions
using JSON

export IDMapper, EnrichmentTerm, EnrichmentDatabase, EnrichmentResult
export load_annotation_database, save_annotation_database, build_annotation_database
export builtin_annotation_database, builtin_annotation_terms
export map_id, map_ids, enrichment_test, go_enrichment, kegg_enrichment, dotplot

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
struct EnrichmentResult
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
end

"""
    Base.show(io::IO, result::EnrichmentResult)

Print a compact one-line summary for an enrichment result.
"""
function Base.show(io::IO, result::EnrichmentResult)
    print(io, "EnrichmentResult(", result.term_id, ", padj=", round(result.padj, digits=4), ", overlap=", result.overlap, ")")
end

"""
    map_id(mapper::IDMapper, identifier::AbstractString)

Map a single identifier through the database mapper, falling back to the input.
"""
function map_id(mapper::IDMapper, identifier::AbstractString)
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
function build_annotation_database(terms::AbstractVector{<:EnrichmentTerm}; mapper::IDMapper=IDMapper())
    return EnrichmentDatabase(Dict(term.id => term for term in terms), mapper)
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
    return build_annotation_database(builtin_annotation_terms())
end

"""
    save_annotation_database(path, database)

Write an enrichment database to a JSON file.
"""
function save_annotation_database(path::AbstractString, database::EnrichmentDatabase)
    payload = Dict{String,Any}(
        "terms" => [Dict(
            "id" => term.id,
            "name" => term.name,
            "namespace" => term.namespace,
            "genes" => term.genes,
            "parents" => term.parents,
        ) for term in values(database.terms)],
        "mapper" => Dict(
            "forward" => database.mapper.forward,
            "reverse" => database.mapper.reverse,
        ),
    )
    open(path, "w") do io
        JSON.print(io, payload)
    end
    return path
end

"""
    load_annotation_database(path)

Load an enrichment database from a JSON file.
"""
function load_annotation_database(path::AbstractString)
    payload = JSON.parsefile(path)
    terms = EnrichmentTerm[]
    for item in payload["terms"]
        push!(terms, EnrichmentTerm(String(item["id"]), String(item["name"]), String(item["namespace"]), String.(item["genes"]), String.(item["parents"])))
    end
    mapper_payload = get(payload, "mapper", Dict{String,Any}())
    forward = Dict{String,String}(String(key) => String(value) for (key, value) in get(mapper_payload, "forward", Dict{String,Any}()))
    reverse = Dict{String,Vector{String}}(String(key) => String.(value) for (key, value) in get(mapper_payload, "reverse", Dict{String,Any}()))
    return EnrichmentDatabase(Dict(term.id => term for term in terms), IDMapper(forward, reverse))
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
function enrichment_test(query_genes, database::EnrichmentDatabase; background=nothing, namespace::Union{Nothing,AbstractString}=nothing, elim::Bool=false, threaded::Bool=true, min_overlap::Int=1)
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
    sort!(scored; by = _result_sort_key)
    return scored
end

"""
    go_enrichment(query_genes, database; kwargs...)

Convenience wrapper around [`enrichment_test`](@ref) restricted to the GO namespace.
"""
function go_enrichment(query_genes, database::EnrichmentDatabase; kwargs...)
    return enrichment_test(query_genes, database; namespace="GO", kwargs...)
end

"""
    kegg_enrichment(query_genes, database; kwargs...)

Convenience wrapper around [`enrichment_test`](@ref) restricted to the KEGG namespace.
"""
function kegg_enrichment(query_genes, database::EnrichmentDatabase; kwargs...)
    return enrichment_test(query_genes, database; namespace="KEGG", kwargs...)
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
function dotplot(results::AbstractVector{<:EnrichmentResult}; top_n::Int=20, save_path::Union{Nothing,AbstractString}=nothing)
    selected = first(sort(collect(results); by = _result_sort_key), min(top_n, length(results)))
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

end