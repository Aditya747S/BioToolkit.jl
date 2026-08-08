module GeneAnnotation

using DataFrames
using Arrow
using JSON
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, active_provenance_context, register_provenance!, provenance_result!, ensembl, entrez, symbol, refseq, uniprot, GeneIdType, SummarizedExperiment, EnsemblID, EntrezID, SymbolID, RefSeqID, UniProtID
using ..Enrichment: IDMapper, EnrichmentDatabase, EnrichmentTerm, build_annotation_database

export OrganismAnnotationDb, IdConversionResult, PathwayMappingResult
export AnnotationDbCache, AnnotationRecord, OrthologMap, AnnotationDbBuildResult
export download_annotation_db, load_annotation_db, save_annotation_db, convert_ids, map_to_organism, pathway_mappings, annotate_summarized_experiment, build_enrichment_database
export build_annotation_db_from_files, annotation_db_from_gtf


struct AnnotationDbCache <: AbstractAnalysisResult
    cache_dir::String
    files::Dict{Symbol,String}
    checksums::Dict{Symbol,UInt64}
    provenance::ResultProvenance
end

struct AnnotationRecord <: AbstractAnalysisResult
    gene_id::String
    transcript_id::Union{Nothing,String}
    symbol::Union{Nothing,String}
    entrez_id::Union{Nothing,String}
    refseq_id::Union{Nothing,String}
    uniprot_id::Union{Nothing,String}
    source::Symbol
    attributes::Dict{String,String}
    provenance::ResultProvenance
end

struct OrthologMap <: AbstractAnalysisResult
    source_organism::Symbol
    target_organism::Symbol
    map::Dict{String,Vector{String}}
    method::Symbol
    provenance::ResultProvenance
end

struct AnnotationDbBuildResult <: AbstractAnalysisResult
    db::Any
    records::Vector{AnnotationRecord}
    cache::AnnotationDbCache
    warnings::Vector{String}
    provenance::ResultProvenance
end

struct OrganismAnnotationDb <: AbstractAnalysisResult
    organism::Symbol
    id_mapper::IDMapper
    term_databases::Dict{Symbol,EnrichmentDatabase}
    provenance::ResultProvenance
end

OrganismAnnotationDb(organism::Symbol, id_mapper::IDMapper, term_databases::Dict{Symbol,EnrichmentDatabase}) =
    OrganismAnnotationDb(organism, id_mapper, term_databases, provenance_record("OrganismAnnotationDb", "GeneAnnotation"))

struct IdConversionResult <: AbstractAnalysisResult
    mapped::Dict{String,Vector{String}}
    unmapped::Vector{String}
    from_type::GeneIdType
    to_type::GeneIdType
    provenance::ResultProvenance
end

IdConversionResult(mapped, unmapped, from_type, to_type) =
    IdConversionResult(mapped, unmapped, from_type, to_type, provenance_record("IdConversionResult", "GeneAnnotation"))

struct PathwayMappingResult <: AbstractAnalysisResult
    mapped_pathways::Dict{String,Vector{String}}
    unmapped_genes::Vector{String}
    provenance::ResultProvenance
end

PathwayMappingResult(mapped_pathways, unmapped_genes) =
    PathwayMappingResult(mapped_pathways, unmapped_genes, provenance_record("PathwayMappingResult", "GeneAnnotation"))


function _push_mapping!(forward::Dict{String,String}, reverse::Dict{String,Vector{String}}, from::Union{Nothing,String}, to::Union{Nothing,String})
    (from === nothing || to === nothing || isempty(from) || isempty(to)) && return nothing
    forward[from] = to
    push!(get!(reverse, to, String[]), from)
    return nothing
end

function _annotation_checksum(path::String)
    h = zero(UInt64)
    open(path, "r") do io
        for b in read(io)
            h = hash(b, h)
        end
    end
    return h
end

function _parse_gtf_attributes(attr::AbstractString)
    attrs = Dict{String,String}()
    for token in split(attr, ';')
        item = strip(token)
        isempty(item) && continue
        if occursin("=", item)
            key, val = split(item, '='; limit=2)
        else
            parts = split(item; limit=2)
            length(parts) == 2 || continue
            key, val = parts
        end
        attrs[String(strip(key))] = String(strip(strip(val), ['"', '\'']))
    end
    return attrs
end

function _read_gtf_records(path::String, source::Symbol)
    isfile(path) || throw(ArgumentError("annotation source file does not exist: $path"))
    records = AnnotationRecord[]
    open(path, "r") do io
        for line in eachline(io)
            isempty(line) && continue
            startswith(line, '#') && continue
            cols = split(line, '\t')
            length(cols) >= 9 || throw(ArgumentError("malformed GTF/GFF3 row in $path: expected at least 9 columns"))
            attrs = _parse_gtf_attributes(cols[9])
            gene_id = get(attrs, "gene_id", get(attrs, "ID", get(attrs, "gene", "")))
            isempty(gene_id) && continue
            tx_id = get(attrs, "transcript_id", get(attrs, "Parent", nothing))
            sym = get(attrs, "gene_name", get(attrs, "Name", get(attrs, "gene_symbol", nothing)))
            push!(records, AnnotationRecord(gene_id, tx_id, sym, nothing, nothing, nothing, source, attrs, provenance_record("AnnotationRecord", "GeneAnnotation/$(source)")))
        end
    end
    return records
end

function _read_gene_info_records(path::String)
    isfile(path) || throw(ArgumentError("NCBI gene_info source file does not exist: $path"))
    records = AnnotationRecord[]
    open(path, "r") do io
        header = split(readline(io), '\t')
        geneid_i = findfirst(==("GeneID"), header)
        symbol_i = findfirst(==("Symbol"), header)
        geneid_i === nothing && throw(ArgumentError("gene_info file is missing GeneID column"))
        for line in eachline(io)
            isempty(line) && continue
            cols = split(line, '\t')
            length(cols) < geneid_i && continue
            entrez_id = String(cols[geneid_i])
            sym = symbol_i === nothing || length(cols) < symbol_i ? nothing : String(cols[symbol_i])
            push!(records, AnnotationRecord(entrez_id, nothing, sym, entrez_id, nothing, nothing, :ncbi_gene_info, Dict{String,String}(), provenance_record("AnnotationRecord", "GeneAnnotation/gene_info")))
        end
    end
    return records
end

function _read_uniprot_mapping_records(path::String)
    isfile(path) || throw(ArgumentError("UniProt mapping file does not exist: $path"))
    records = AnnotationRecord[]
    open(path, "r") do io
        header = lowercase.(split(readline(io), '\t'))
        uni_i = findfirst(x -> occursin("uniprot", x) || x in ("entry", "from"), header)
        gene_i = findfirst(x -> occursin("ensembl", x) || x in ("to", "geneid"), header)
        sym_i = findfirst(x -> x in ("symbol", "gene names", "genes"), header)
        (uni_i === nothing || gene_i === nothing) && throw(ArgumentError("UniProt mapping file requires UniProt and gene identifier columns"))
        for line in eachline(io)
            isempty(line) && continue
            cols = split(line, '\t')
            length(cols) < max(uni_i, gene_i) && continue
            gene_id = String(cols[gene_i])
            uni = String(cols[uni_i])
            sym = sym_i === nothing || length(cols) < sym_i ? nothing : String(cols[sym_i])
            push!(records, AnnotationRecord(gene_id, nothing, sym, nothing, nothing, uni, :uniprot_mapping, Dict{String,String}(), provenance_record("AnnotationRecord", "GeneAnnotation/uniprot_mapping")))
        end
    end
    return records
end

function _build_db_from_records(organism::Symbol, records::Vector{AnnotationRecord}, cache::AnnotationDbCache, warnings::Vector{String})
    forward = Dict{String,String}()
    reverse = Dict{String,Vector{String}}()
    for rec in records
        canonical = rec.symbol === nothing ? rec.gene_id : rec.symbol
        _push_mapping!(forward, reverse, rec.gene_id, canonical)
        _push_mapping!(forward, reverse, rec.transcript_id, rec.gene_id)
        _push_mapping!(forward, reverse, rec.entrez_id, canonical)
        _push_mapping!(forward, reverse, rec.refseq_id, canonical)
        _push_mapping!(forward, reverse, rec.uniprot_id, canonical)
        _push_mapping!(forward, reverse, rec.symbol, rec.gene_id)
    end
    db = OrganismAnnotationDb(organism, IDMapper(forward, reverse), Dict{Symbol,EnrichmentDatabase}())
    prov = provenance_record("AnnotationDbBuildResult", "GeneAnnotation/build_annotation_db_from_files"; notes=warnings, parameters=(organism=organism, record_count=length(records), source_count=length(cache.files)))
    return AnnotationDbBuildResult(db, records, cache, warnings, prov)
end

function build_annotation_db_from_files(organism::Symbol; gtf::Union{Nothing,String}=nothing, gff3::Union{Nothing,String}=nothing, gene_info::Union{Nothing,String}=nothing, uniprot_mapping::Union{Nothing,String}=nothing, go_gaf::Union{Nothing,String}=nothing, reactome::Union{Nothing,String}=nothing, kegg::Union{Nothing,String}=nothing, msigdb::Union{Nothing,String}=nothing, cache_dir=nothing)
    _ctx = active_provenance_context()
    records = AnnotationRecord[]
    files = Dict{Symbol,String}()
    checksums = Dict{Symbol,UInt64}()
    warnings = String[]
    for (key, path) in ((:gtf, gtf), (:gff3, gff3), (:gene_info, gene_info), (:uniprot_mapping, uniprot_mapping), (:go_gaf, go_gaf), (:reactome, reactome), (:kegg, kegg), (:msigdb, msigdb))
        path === nothing && continue
        isfile(path) || throw(ArgumentError("annotation source file does not exist: $path"))
        files[key] = path
        checksums[key] = _annotation_checksum(path)
    end
    gtf !== nothing && append!(records, _read_gtf_records(gtf, :gtf))
    gff3 !== nothing && append!(records, _read_gtf_records(gff3, :gff3))
    gene_info !== nothing && append!(records, _read_gene_info_records(gene_info))
    uniprot_mapping !== nothing && append!(records, _read_uniprot_mapping_records(uniprot_mapping))
    isempty(records) && throw(ArgumentError("no annotation records were parsed from the provided files"))
    for key in (:go_gaf, :reactome, :kegg, :msigdb)
        haskey(files, key) && push!(warnings, "$(key) source was checksummed but term parsing is not enabled in this build path")
    end
    cache = AnnotationDbCache(cache_dir === nothing ? "" : String(cache_dir), files, checksums, provenance_record("AnnotationDbCache", "GeneAnnotation/build_annotation_db_from_files"))
    result = _build_db_from_records(organism, records, cache, warnings)
    return provenance_result!(_ctx, result, "build_annotation_db_from_files"; parents=String[], parameters=(organism=organism, record_count=length(records)))
end

annotation_db_from_gtf(path::String; organism::Symbol=:unknown, cache_dir=nothing) = build_annotation_db_from_files(organism; gtf=path, cache_dir=cache_dir)

function download_annotation_db(organism::Symbol; sources=[:ensembl,:entrez,:symbol,:refseq,:uniprot], cache_dir=nothing, gtf=nothing, gff3=nothing, gene_info=nothing, uniprot_mapping=nothing)
    _ctx = active_provenance_context()
    if any(!isnothing, (gtf, gff3, gene_info, uniprot_mapping))
        return build_annotation_db_from_files(organism; gtf=gtf, gff3=gff3, gene_info=gene_info, uniprot_mapping=uniprot_mapping, cache_dir=cache_dir).db
    end
    # Build cache path
    base_dir = cache_dir === nothing ? joinpath(homedir(), ".biotoolkit", "annotation") : cache_dir
    mkpath(base_dir)
    
    # Simulated mapping tables for offline validation
    forward = Dict{String,String}()
    reverse = Dict{String,Vector{String}}()
    
    mock_genes = [
        ("TP53", "ENSG00000141510", "7157", "NM_000546", "P04637"),
        ("BRCA1", "ENSG00000012048", "672", "NM_007294", "P38398"),
        ("GAPDH", "ENSG00000111640", "2597", "NM_002046", "P04406"),
        ("ACTB", "ENSG00000075624", "60", "NM_001101", "P60709")
    ]
    
    for (sym, ens, ent, ref, uni) in mock_genes
        # map ensembl -> symbol
        forward[ens] = sym
        push!(get!(reverse, sym, String[]), ens)
        
        # map entrez -> symbol
        forward[ent] = sym
        push!(get!(reverse, sym, String[]), ent)
        
        # map refseq -> symbol
        forward[ref] = sym
        push!(get!(reverse, sym, String[]), ref)
        
        # map uniprot -> symbol
        forward[uni] = sym
        push!(get!(reverse, sym, String[]), uni)
        
        # map symbol -> ensembl
        forward[sym] = ens
        push!(get!(reverse, ens, String[]), sym)
        
        # cross mapping
        forward[ens * "_entrez"] = ent
        forward[ent * "_ensembl"] = ens
    end
    
    id_mapper = IDMapper(forward, reverse)
    
    terms = EnrichmentTerm[
        EnrichmentTerm("GO:0008150", "biological_process", "GO", ["TP53", "BRCA1", "GAPDH", "ACTB"], String[]),
        EnrichmentTerm("GO:0006974", "response to DNA damage stimulus", "GO", ["TP53", "BRCA1"], ["GO:0008150"]),
        EnrichmentTerm("KEGG:hsa04110", "Cell cycle", "KEGG", ["TP53", "BRCA1"], String[])
    ]
    
    go_db = build_annotation_database(terms; mapper=id_mapper, _ctx=_ctx)
    term_databases = Dict{Symbol,EnrichmentDatabase}(:go => go_db)
    
    result = OrganismAnnotationDb(organism, id_mapper, term_databases, provenance_record("OrganismAnnotationDb", "GeneAnnotation/download_annotation_db"; status=:warning, fallbacks=["no local annotation sources supplied; returned built-in demo annotation database"]))
    return provenance_result!(_ctx, result, "download_annotation_db"; parents=String[], parameters=(organism=organism, sources=string.(sources), mode="built_in_demo"))
end

function convert_ids(db::OrganismAnnotationDb, ids, from::GeneIdType, to::GeneIdType)
    _ctx = active_provenance_context()
    mapped = Dict{String,Vector{String}}()
    unmapped = String[]
    
    for id in String.(ids)
        res = String[]
        if haskey(db.id_mapper.forward, id)
            push!(res, db.id_mapper.forward[id])
        end
        if haskey(db.id_mapper.reverse, id)
            append!(res, db.id_mapper.reverse[id])
        end
        
        filtered_res = String[]
        for r in res
            if to isa SymbolID
                if !startswith(r, "ENS") && !all(isdigit, r) && !startswith(r, "NM_")
                    push!(filtered_res, r)
                end
            elseif to isa EnsemblID
                if startswith(r, "ENS")
                    push!(filtered_res, r)
                end
            elseif to isa EntrezID
                if all(isdigit, r)
                    push!(filtered_res, r)
                end
            elseif to isa RefSeqID
                if startswith(r, "NM_")
                    push!(filtered_res, r)
                end
            elseif to isa UniProtID
                if length(r) == 6 && isletter(r[1])
                    push!(filtered_res, r)
                end
            else
                push!(filtered_res, r)
            end
        end
        
        if isempty(filtered_res)
            # Hop mapping (multi-hop search)
            for intermediate in res
                if haskey(db.id_mapper.forward, intermediate)
                    push!(filtered_res, db.id_mapper.forward[intermediate])
                end
                if haskey(db.id_mapper.reverse, intermediate)
                    append!(filtered_res, db.id_mapper.reverse[intermediate])
                end
            end
            unique!(filtered_res)
            final_res = String[]
            for r in filtered_res
                if to isa SymbolID
                    if !startswith(r, "ENS") && !all(isdigit, r) && !startswith(r, "NM_")
                        push!(final_res, r)
                    end
                elseif to isa EnsemblID
                    if startswith(r, "ENS")
                        push!(final_res, r)
                    end
                elseif to isa EntrezID
                    if all(isdigit, r)
                        push!(final_res, r)
                    end
                elseif to isa RefSeqID
                    if startswith(r, "NM_")
                        push!(final_res, r)
                    end
                elseif to isa UniProtID
                    if length(r) == 6 && isletter(r[1])
                        push!(final_res, r)
                    end
                else
                    push!(final_res, r)
                end
            end
            filtered_res = final_res
        end
        
        if isempty(filtered_res)
            push!(unmapped, id)
        else
            mapped[id] = unique!(filtered_res)
        end
    end
    
    result = IdConversionResult(mapped, unmapped, from, to)
    return provenance_result!(_ctx, result, "convert_ids"; parents=[db.provenance.id], parameters=(from=string(Symbol(from)), to=string(Symbol(to))))
end

convert_ids(db::OrganismAnnotationDb, ids, from::Symbol, to::Symbol) =
    convert_ids(db, ids, GeneIdType(from), GeneIdType(to))

function save_annotation_db(db::OrganismAnnotationDb, path::String)
    _ctx = active_provenance_context()
    payload = Dict{String,Any}(
        "organism" => string(db.organism),
        "id_mapper" => Dict(
            "forward" => db.id_mapper.forward,
            "reverse" => db.id_mapper.reverse
        ),
        "term_databases" => Dict(
            string(k) => [
                Dict(
                    "id" => term.id,
                    "name" => term.name,
                    "namespace" => term.namespace,
                    "genes" => term.genes,
                    "parents" => term.parents
                ) for term in values(v.terms)
            ] for (k, v) in db.term_databases
        )
    )
    
    if endswith(path, ".arrow")
        df = DataFrame(
            from_id = collect(keys(db.id_mapper.forward)),
            to_id = collect(values(db.id_mapper.forward))
        )
        Arrow.write(path, df)
    else
        open(path, "w") do io
            JSON.print(io, payload)
        end
    end
    
    return provenance_result!(_ctx, path, "save_annotation_db"; parents=[db.provenance.id], parameters=(path=path,))
end

function load_annotation_db(path::String)
    _ctx = active_provenance_context()
    if endswith(path, ".arrow")
        table = Arrow.Table(path)
        df = DataFrame(table)
        forward = Dict{String,String}()
        reverse = Dict{String,Vector{String}}()
        for row in eachrow(df)
            forward[String(row.from_id)] = String(row.to_id)
            push!(get!(reverse, String(row.to_id), String[]), String(row.from_id))
        end
        id_mapper = IDMapper(forward, reverse)
        result = OrganismAnnotationDb(:unknown, id_mapper, Dict{Symbol,EnrichmentDatabase}())
    else
        payload = JSON.parsefile(path)
        organism = Symbol(payload["organism"])
        mapper_payload = payload["id_mapper"]
        forward = Dict{String,String}(String(k) => String(v) for (k, v) in mapper_payload["forward"])
        reverse = Dict{String,Vector{String}}(String(k) => String.(v) for (k, v) in mapper_payload["reverse"])
        id_mapper = IDMapper(forward, reverse)
        
        term_databases = Dict{Symbol,EnrichmentDatabase}()
        for (k, v) in payload["term_databases"]
            terms = EnrichmentTerm[]
            for item in v
                push!(terms, EnrichmentTerm(String(item["id"]), String(item["name"]), String(item["namespace"]), String.(item["genes"]), String.(item["parents"])))
            end
            db = build_annotation_database(terms; mapper=id_mapper, _ctx=_ctx)
            term_databases[Symbol(k)] = db
        end
        result = OrganismAnnotationDb(organism, id_mapper, term_databases)
    end
    return provenance_result!(_ctx, result, "load_annotation_db"; parents=String[], parameters=(path=path,))
end

function map_to_organism(db::OrganismAnnotationDb, ids, target::Symbol)
    _ctx = active_provenance_context()
    mapped = Dict{String,Vector{String}}()
    unmapped = String[]
    for id in String.(ids)
        if target === :mouse
            val = titlecase(lowercase(id))
            mapped[id] = [val]
        elseif target === :human
            val = uppercase(id)
            mapped[id] = [val]
        else
            mapped[id] = [id]
        end
    end
    result = IdConversionResult(mapped, unmapped, symbol, symbol)
    return provenance_result!(_ctx, result, "map_to_organism"; parents=[db.provenance.id], parameters=(target=target,))
end

function pathway_mappings(db::OrganismAnnotationDb, genes; databases=[:go, :kegg])
    _ctx = active_provenance_context()
    mapped_pathways = Dict{String,Vector{String}}()
    unmapped_genes = String[]
    
    for gene in String.(genes)
        res = [gene]
        if haskey(db.id_mapper.forward, gene)
            push!(res, db.id_mapper.forward[gene])
        end
        found = false
        for db_name in databases
            haskey(db.term_databases, db_name) || continue
            term_db = db.term_databases[db_name]
            for term in values(term_db.terms)
                for r in res
                    if r in term.genes
                        push!(get!(mapped_pathways, gene, String[]), term.id)
                        found = true
                    end
                end
            end
        end
        if found
            unique!(mapped_pathways[gene])
        else
            push!(unmapped_genes, gene)
        end
    end
    
    result = PathwayMappingResult(mapped_pathways, unmapped_genes)
    return provenance_result!(_ctx, result, "pathway_mappings"; parents=[db.provenance.id], parameters=(databases=string.(databases),))
end

function annotate_summarized_experiment(se::SummarizedExperiment; id_col::Symbol=:gene_id, db::OrganismAnnotationDb, to::Symbol=:symbol)
    _ctx = active_provenance_context()
    haskey(se.rowData, id_col) || throw(ArgumentError("rowData does not contain column $id_col"))
    
    ids = se.rowData[id_col]
    conversion = convert_ids(db, ids, :ensembl, to)
    
    new_vals = String[]
    for id in ids
        str_id = String(id)
        if haskey(conversion.mapped, str_id)
            push!(new_vals, first(conversion.mapped[str_id]))
        else
            push!(new_vals, str_id)
        end
    end
    
    new_rowData = copy(se.rowData)
    new_rowData[to] = new_vals
    
    result = SummarizedExperiment(se.assays, new_rowData, copy(se.colData), copy(se.metadata))
    # We can retrieve provenance id if available
    parent_id = haskey(se.metadata, :provenance) ? se.metadata[:provenance].id : ""
    return provenance_result!(_ctx, result, "annotate_summarized_experiment"; parents=[parent_id, db.provenance.id])
end

function build_enrichment_database(db::OrganismAnnotationDb; db_type::Symbol=:go)
    _ctx = active_provenance_context()
    if haskey(db.term_databases, db_type)
        orig_db = db.term_databases[db_type]
        result = EnrichmentDatabase(orig_db.terms, db.id_mapper)
        return provenance_result!(_ctx, result, "build_enrichment_database"; parents=[db.provenance.id], parameters=(db_type=db_type,))
    else
        throw(ArgumentError("Database type $db_type not found in OrganismAnnotationDb"))
    end
end

end
