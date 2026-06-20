module Restriction

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, BioSequence, DNAAlphabet
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, analysis_result_summary, new_provenance_id, register_provenance!

export RestrictionEnzyme, RestrictionSite, restriction_enzymes, restriction_enzyme, restriction_enzyme_names, restriction_catalog, restriction_sites, restriction_digest_map, find_restriction_sites, digest_sequence

struct RestrictionEnzyme
    name::String
    recognition_site::BioSequence{DNAAlphabet}
    cut_offset::Int
    regex::Regex
end

struct RestrictionSite
    enzyme::RestrictionEnzyme
    position::Int
    cut_position::Int
end

Base.show(io::IO, enzyme::RestrictionEnzyme) = print(io, "RestrictionEnzyme($(enzyme.name), $(enzyme.recognition_site), cut=$(enzyme.cut_offset))")
Base.show(io::IO, site::RestrictionSite) = print(io, "RestrictionSite($(site.enzyme.name)@$(site.position)->$(site.cut_position))")

const _IUPAC_PATTERNS = Dict{Char,String}(
    'A' => "A",
    'C' => "C",
    'G' => "G",
    'T' => "T",
    'U' => "T",
    'R' => "[AG]",
    'Y' => "[CT]",
    'S' => "[GC]",
    'W' => "[AT]",
    'K' => "[GT]",
    'M' => "[AC]",
    'B' => "[CGT]",
    'D' => "[AGT]",
    'H' => "[ACT]",
    'V' => "[ACG]",
    'N' => "[ACGT]")

function _restriction_regex(site::BioSequence{DNAAlphabet})
    # Regex search is still the most robust way to find overlapping IUPAC matches
    # but we store the site as a typed sequence for scientific parity.
    site_str = String(site)
    pattern = IOBuffer()
    print(pattern, "(?i)")
    for character in site_str
        print(pattern, get(_IUPAC_PATTERNS, character, string(character)))
    end
    return Regex(String(take!(pattern)))
end

function _restriction_enzyme(name::String, recognition_site::String, cut_offset::Integer)
    site_seq = BioSequence{DNAAlphabet}(recognition_site)
    return RestrictionEnzyme(String(name), site_seq, Int(cut_offset), _restriction_regex(site_seq))
end

const _RESTRICTION_DATABASE = Dict{String,RestrictionEnzyme}(
    enzyme.name => enzyme for enzyme in (
        _restriction_enzyme("EcoRI", "GAATTC", 1),
        _restriction_enzyme("BamHI", "GGATCC", 1),
        _restriction_enzyme("HindIII", "AAGCTT", 1),
        _restriction_enzyme("NotI", "GCGGCCGC", 2),
        _restriction_enzyme("PstI", "CTGCAG", 5),
        _restriction_enzyme("SmaI", "CCCGGG", 3),
        _restriction_enzyme("XhoI", "CTCGAG", 1),
        _restriction_enzyme("SalI", "GTCGAC", 1),
        _restriction_enzyme("NdeI", "CATATG", 2),
        _restriction_enzyme("NcoI", "CCATGG", 1),
        _restriction_enzyme("KpnI", "GGTACC", 5),
        _restriction_enzyme("ClaI", "ATCGAT", 2),
        _restriction_enzyme("PvuII", "CAGCTG", 3),
        _restriction_enzyme("HaeIII", "GGCC", 2),
        _restriction_enzyme("AluI", "AGCT", 2),
        _restriction_enzyme("DpnI", "GATC", 2),
        _restriction_enzyme("MspI", "CCGG", 1),
        _restriction_enzyme("HinfI", "GANTC", 1),
        _restriction_enzyme("TaqI", "TCGA", 1),
        _restriction_enzyme("BstEII", "GGTNACC", 1),
        _restriction_enzyme("AatII", "GACGTC", 1),
        _restriction_enzyme("AbaSI", "CACNNNNGTG", 1),
        _restriction_enzyme("ApaI", "GGGCCC", 1),
        _restriction_enzyme("AseI", "ATTAAT", 1),
        _restriction_enzyme("AvrII", "CCTAGG", 1),
        _restriction_enzyme("BclI", "TGATCA", 1),
        _restriction_enzyme("BglII", "AGATCT", 1),
        _restriction_enzyme("BsaI", "GGTCTC", 1),
        _restriction_enzyme("BsaXI", "ACNNNNGTAYC", 1),
        _restriction_enzyme("BseRI", "GAGGAG", 10),
        _restriction_enzyme("BsmBI", "CGTCTC", 1),
        _restriction_enzyme("DraI", "TTTAAA", 3),
        _restriction_enzyme("EcoRV", "GATATC", 3),
        _restriction_enzyme("FokI", "GGATG", 9),
        _restriction_enzyme("HaeII", "RGCGCY", 2),
        _restriction_enzyme("HpaI", "GTTAAC", 3),
        _restriction_enzyme("MluI", "ACGCGT", 1),
        _restriction_enzyme("NheI", "GCTAGC", 1),
        _restriction_enzyme("PacI", "TTAATTAA", 5),
        _restriction_enzyme("SacI", "GAGCTC", 1),
        _restriction_enzyme("SacII", "CCGCGG", 1),
        _restriction_enzyme("SpeI", "ACTAGT", 1),
        _restriction_enzyme("StuI", "AGGCCT", 3),
        _restriction_enzyme("XbaI", "TCTAGA", 1),
        _restriction_enzyme("XmaI", "CCCGGG", 1),
        _restriction_enzyme("XmnI", "GAANNNNTTC", 5),
        _restriction_enzyme("BspHI", "TCATGA", 1),
        _restriction_enzyme("BspEI", "TCCGGA", 1),
        _restriction_enzyme("BspMI", "ACCTGC", 1),
        _restriction_enzyme("CspCI", "CAANNNNNGTGG", 1),
        _restriction_enzyme("MfeI", "CAATTG", 1),
        _restriction_enzyme("MluCI", "AATT", 1),
        _restriction_enzyme("PciI", "ACATGT", 1),
        _restriction_enzyme("PmlI", "CACGTG", 3),
        _restriction_enzyme("SbfI", "CCTGCAGG", 1),
        _restriction_enzyme("SphI", "GCATGC", 5),
        _restriction_enzyme("Sse8387I", "CCTGCAGG", 2),
        _restriction_enzyme("AgeI", "ACCGGT", 1),
        _restriction_enzyme("AflII", "CTTAAG", 1),
        _restriction_enzyme("AscI", "GGCGCGCC", 2),
        _restriction_enzyme("AsiSI", "GCGATCGC", 1),
        _restriction_enzyme("BsaI", "GGTCTC", 1),
        _restriction_enzyme("BsiWI", "CGTACG", 1),
        _restriction_enzyme("BsrGI", "TGTACA", 1),
        _restriction_enzyme("BstBI", "TTCGAA", 1),
        _restriction_enzyme("FspI", "TGCGCA", 3),
        _restriction_enzyme("NruI", "TCGCGA", 3),
        _restriction_enzyme("PmeI", "GTTTAAAC", 4),
        _restriction_enzyme("RsrII", "CGGWCCG", 2),
        _restriction_enzyme("SfiI", "GGCCNNNNNGGCC", 8),
        _restriction_enzyme("SfoI", "GGCGCC", 3),
        _restriction_enzyme("SrfI", "GCCCGGGC", 4),
        _restriction_enzyme("SwaI", "ATTTAAAT", 4),
        _restriction_enzyme("TspMI", "CCGG", 1),
        _restriction_enzyme("XcmI", "CCANNNNNNNNNTGG", 1),
        _restriction_enzyme("XmaIII", "CGGCCG", 1),
        _restriction_enzyme("ZraI", "GACGTC", 3))
)

restriction_enzymes() = _RESTRICTION_DATABASE
restriction_catalog() = collect(values(_RESTRICTION_DATABASE))
restriction_enzyme_names() = sort!(collect(keys(_RESTRICTION_DATABASE)))

function restriction_enzyme(name::String)
    enzyme = get(_RESTRICTION_DATABASE, String(name), nothing)
    enzyme === nothing && throw(ArgumentError("unknown restriction enzyme: $(name)"))
    return enzyme
end

restriction_sites(sequence::BioSequence{DNAAlphabet}, enzyme::String) = find_restriction_sites(sequence, enzyme)
restriction_sites(sequence::BioSequence{DNAAlphabet}, enzyme::RestrictionEnzyme) = find_restriction_sites(sequence, enzyme)
restriction_sites(sequence::BioSequence{DNAAlphabet}) = find_restriction_sites(sequence)

@inline _restriction_as_dna(sequence::AbstractString) = BioSequence{DNAAlphabet}(String(sequence); validate=false)

restriction_sites(sequence::AbstractString, enzyme::String) = find_restriction_sites(_restriction_as_dna(sequence), enzyme)
restriction_sites(sequence::AbstractString, enzyme::RestrictionEnzyme) = find_restriction_sites(_restriction_as_dna(sequence), enzyme)
restriction_sites(sequence::AbstractString) = find_restriction_sites(_restriction_as_dna(sequence))

function find_restriction_sites(sequence::BioSequence{DNAAlphabet}, enzyme::RestrictionEnzyme)
    sites = RestrictionSite[]
    sequence_str = String(sequence)
    for match in eachmatch(enzyme.regex, sequence_str; overlap=true)
        pos = match.offset
        # cut_offset is 1-indexed relative to the match
        cut_pos = pos + enzyme.cut_offset - 1
        push!(sites, RestrictionSite(enzyme, pos, cut_pos))
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "find_restriction_sites";
        parameters=(enzyme=enzyme.name, sequence_length=length(sequence), n_sites=length(sites)))
    end
    return sites
end
find_restriction_sites(sequence::BioSequence{DNAAlphabet}, name::String) = find_restriction_sites(sequence, restriction_enzyme(name))
find_restriction_sites(sequence::AbstractString, enzyme::RestrictionEnzyme) = find_restriction_sites(_restriction_as_dna(sequence), enzyme)
find_restriction_sites(sequence::AbstractString, name::String) = find_restriction_sites(_restriction_as_dna(sequence), name)

function find_restriction_sites(sequence::BioSequence{DNAAlphabet})
    sites = RestrictionSite[]
    for enzyme in values(_RESTRICTION_DATABASE)
        append!(sites, find_restriction_sites(sequence, enzyme))
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "find_restriction_sites_all";
        parameters=(n_enzymes=length(_RESTRICTION_DATABASE), sequence_length=length(sequence), n_sites=length(sites)))
    end
    return sites
end
find_restriction_sites(sequence::AbstractString) = find_restriction_sites(_restriction_as_dna(sequence))

function digest_sequence(sequence::BioSequence{DNAAlphabet}, enzymes::AbstractVector{<:RestrictionEnzyme})
    hits = RestrictionSite[]
    for enzyme in enzymes
        append!(hits, find_restriction_sites(sequence, enzyme))
    end
    sort!(hits; by = hit -> hit.cut_position)

    length(sequence) == 0 && return [sequence]
    
    # Unique cut positions
    cut_points = Int[0, length(sequence)]
    for hit in hits
        if 1 < hit.cut_position <= length(sequence)
            push!(cut_points, hit.cut_position)
        end
    end
    unique!(sort!(cut_points))

    fragments = BioSequence{DNAAlphabet}[]
    for index in 1:length(cut_points)-1
        start_pos = cut_points[index] + 1
        stop_pos = cut_points[index + 1]
        push!(fragments, sequence[start_pos:stop_pos])
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "digest_sequence";
        parameters=(n_enzymes=length(enzymes), sequence_length=length(sequence), n_fragments=length(fragments)))
    end
    return fragments
end

function digest_sequence(sequence::BioSequence{DNAAlphabet}, enzyme_names::AbstractVector{<:String})
    enzymes = [restriction_enzyme(name) for name in enzyme_names]
    return digest_sequence(sequence, enzymes)
end

# Removed digest_sequence(::AbstractString) wrappers - use BioSequence{DNAAlphabet} instead

digest_sequence(sequence::BioSequence{DNAAlphabet}, enzyme::RestrictionEnzyme) = digest_sequence(sequence, [enzyme])
digest_sequence(sequence::BioSequence{DNAAlphabet}, name::String) = digest_sequence(sequence, [name])
digest_sequence(sequence::BioSequence{DNAAlphabet}) = digest_sequence(sequence, collect(values(_RESTRICTION_DATABASE)))

digest_sequence(sequence::AbstractString, enzyme::RestrictionEnzyme) = digest_sequence(_restriction_as_dna(sequence), [enzyme])
digest_sequence(sequence::AbstractString, name::String) = digest_sequence(_restriction_as_dna(sequence), [name])
digest_sequence(sequence::AbstractString) = digest_sequence(_restriction_as_dna(sequence), collect(values(_RESTRICTION_DATABASE)))

function restriction_digest_map(sequence::BioSequence{DNAAlphabet}, enzymes=restriction_enzyme_names())
    digest = Dict{String,Vector{RestrictionSite}}()
    for enzyme_name in enzymes
        digest[String(enzyme_name)] = find_restriction_sites(sequence, enzyme_name)
    end
    return digest
end

restriction_digest_map(sequence::AbstractString, enzymes=restriction_enzyme_names()) = restriction_digest_map(_restriction_as_dna(sequence), enzymes)

end

module Entrez

using JSON
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, active_provenance_context, analysis_result_summary, register_provenance!

export EntrezSearchResult, EntrezPostResult, entrez_search, entrez_search_ids, entrez_search_count, entrez_fetch, entrez_fetch_fasta, entrez_fetch_genbank, entrez_summary, entrez_post, entrez_post_ids, entrez_post_webenv, entrez_post_query_key, entrez_link, entrez_elink, entrez_link_ids, entrez_link_linksets, entrez_link_records, entrez_pubmed_search, entrez_nuccore_fetch, entrez_nucleotide_search, entrez_protein_search, entrez_gene_search, entrez_taxonomy_search, entrez_genome_search, entrez_genome_fetch, entrez_pubmed_fetch, parse_entrez_search_response, parse_entrez_post_response

struct EntrezSearchResult <: AbstractAnalysisResult
    db::String
    term::String
    ids::Vector{String}
    count::Int
    webenv::Union{Nothing,String}
    query_key::Union{Nothing,Int}
    raw::Dict{String,Any}
    provenance::ResultProvenance
end

EntrezSearchResult(db, term, ids, count, webenv, query_key, raw) =
    EntrezSearchResult(String(db), String(term), String.(ids), Int(count), webenv === nothing ? nothing : String(webenv), query_key === nothing ? nothing : Int(query_key), Dict{String,Any}(String(key) => value for (key, value) in pairs(raw)), provenance_record("EntrezSearchResult", "databases"))

Base.show(io::IO, result::EntrezSearchResult) = print(io, analysis_result_summary(result))

struct EntrezPostResult <: AbstractAnalysisResult
    ids::Vector{String}
    webenv::Union{Nothing,String}
    query_key::Union{Nothing,Int}
    raw::Dict{String,Any}
    provenance::ResultProvenance
end

EntrezPostResult(ids, webenv, query_key, raw) =
    EntrezPostResult(String.(ids), webenv === nothing ? nothing : String(webenv), query_key === nothing ? nothing : Int(query_key), Dict{String,Any}(String(key) => value for (key, value) in pairs(raw)), provenance_record("EntrezPostResult", "databases"))

Base.show(io::IO, result::EntrezPostResult) = print(io, analysis_result_summary(result))

function _urlencode(value::String)
    io = IOBuffer()
    for byte in codeunits(String(value))
        if (byte >= 0x30 && byte <= 0x39) || (byte >= 0x41 && byte <= 0x5a) || (byte >= 0x61 && byte <= 0x7a) || byte in (0x2d, 0x2e, 0x5f, 0x7e)
            write(io, byte)
        elseif byte == 0x20
            write(io, UInt8('+'))
        else
            print(io, '%', uppercase(string(byte, base=16, pad=2)))
        end
    end
    return String(take!(io))
end

function _entrez_url(endpoint::String; parameters)
    query = join(("$(String(key))=$(_urlencode(string(value)))" for (key, value) in pairs(parameters)), "&")
    return "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/$(endpoint)?$(query)"
end

function _download_text(url::String)
    path = tempname()
    try
        download(url, path)
        return read(path, String)
    finally
        isfile(path) && rm(path; force=true)
    end
end

function _download_json(url::String)
    return JSON.parse(_download_text(url))
end

function parse_entrez_search_response(payload::String)
    data = JSON.parse(payload)
    search = data["esearchresult"]
    ids = String.(search["idlist"])
    count = parse(Int, string(search["count"]))
    query_key = haskey(search, "querykey") ? parse(Int, string(search["querykey"])) : nothing
    webenv = haskey(search, "webenv") ? String(search["webenv"]) : nothing
    return EntrezSearchResult("", "", ids, count, webenv, query_key, data)
end

function parse_entrez_post_response(payload::String)
    data = JSON.parse(payload)
    result = get(data, "epostresult", Dict{String,Any}())
    ids = String.(get(result, "idlist", String[]))
    query_key = haskey(result, "querykey") ? parse(Int, string(result["querykey"])) : nothing
    webenv = haskey(result, "webenv") ? String(result["webenv"]) : nothing
    return EntrezPostResult(ids, webenv, query_key, data)
end

function entrez_search(db::String, term::String, retmax::Integer=20, retstart::Integer=0, email::Union{Nothing,String}=nothing, tool::String="BioToolkit", api_key::Union{Nothing,String}=nothing)
    parameters = Dict{String,Any}(
        "db" => db,
        "term" => term,
        "retmode" => "json",
        "retmax" => retmax,
        "retstart" => retstart,
        "tool" => tool)
    email !== nothing && (parameters["email"] = email)
    api_key !== nothing && (parameters["api_key"] = api_key)
    payload = _download_text(_entrez_url("esearch.fcgi"; parameters=parameters))
    result = parse_entrez_search_response(payload)
    search_result = EntrezSearchResult(String(db), String(term), result.ids, result.count, result.webenv, result.query_key, result.raw)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "entrez_search";
        parameters=(db=db, term=term, retmax=Int(retmax), n_results=search_result.count))
    end
    return search_result
end

entrez_search_ids(db::String, term::String; kwargs...) = entrez_search(db, term; kwargs...).ids
entrez_search_count(db::String, term::String; kwargs...) = entrez_search(db, term; kwargs...).count

function entrez_summary(db::String, ids, email::Union{Nothing,String}=nothing, tool::String="BioToolkit", api_key::Union{Nothing,String}=nothing)
    id_string = join(String.(ids), ",")
    parameters = Dict{String,Any}("db" => db, "id" => id_string, "retmode" => "json", "tool" => tool)
    email !== nothing && (parameters["email"] = email)
    api_key !== nothing && (parameters["api_key"] = api_key)
    return _download_json(_entrez_url("esummary.fcgi"; parameters=parameters))
end

function entrez_post(db::String, ids; email::Union{Nothing,String}=nothing, tool::String="BioToolkit", api_key::Union{Nothing,String}=nothing)
    id_string = join(String.(ids), ",")
    parameters = Dict{String,Any}("db" => db, "id" => id_string, "retmode" => "json", "tool" => tool)
    email !== nothing && (parameters["email"] = email)
    api_key !== nothing && (parameters["api_key"] = api_key)
    return parse_entrez_post_response(_download_text(_entrez_url("epost.fcgi"; parameters=parameters)))
end

entrez_post_ids(db::String, ids; kwargs...) = entrez_post(db, ids; kwargs...).ids
entrez_post_webenv(db::String, ids; kwargs...) = entrez_post(db, ids; kwargs...).webenv
entrez_post_query_key(db::String, ids; kwargs...) = entrez_post(db, ids; kwargs...).query_key
entrez_post_ids(result::EntrezPostResult) = result.ids
entrez_post_webenv(result::EntrezPostResult) = result.webenv
entrez_post_query_key(result::EntrezPostResult) = result.query_key

function entrez_link(dbfrom::String, dbto::String, ids, email::Union{Nothing,String}=nothing, tool::String="BioToolkit", api_key::Union{Nothing,String}=nothing)
    id_string = join(String.(ids), ",")
    parameters = Dict{String,Any}("dbfrom" => dbfrom, "db" => dbto, "id" => id_string, "retmode" => "json", "tool" => tool)
    email !== nothing && (parameters["email"] = email)
    api_key !== nothing && (parameters["api_key"] = api_key)
    return _download_json(_entrez_url("elink.fcgi"; parameters=parameters))
end

entrez_elink(dbfrom::String, dbto::String, ids; kwargs...) = entrez_link(dbfrom, dbto, ids; kwargs...)

function entrez_link_linksets(result::AbstractDict)
    return get(result, "linksets", Any[])
end

function entrez_link_records(result::AbstractDict)
    records = Any[]
    for record in entrez_link_linksets(result)
        append!(records, get(record, "linksetdbs", Any[]))
    end
    return records
end

function entrez_link_ids(dbfrom::String, dbto::String, ids, kwargs...)
    data = entrez_link(dbfrom, dbto, ids; kwargs...)
    links = String[]
    for db_links in entrez_link_records(data)
            append!(links, String.(get(db_links, "links", String[])))
    end
    return unique!(links)
end

function entrez_link_ids(result::AbstractDict)
    links = String[]
    for db_links in entrez_link_records(result)
        append!(links, String.(get(db_links, "links", String[])))
    end
    return unique!(links)
end

function entrez_fetch(db::String, ids, rettype::String="fasta", retmode::String="text", email::Union{Nothing,String}=nothing, tool::String="BioToolkit", api_key::Union{Nothing,String}=nothing)
    id_string = join(String.(ids), ",")
    parameters = Dict{String,Any}("db" => db, "id" => id_string, "rettype" => rettype, "retmode" => retmode, "tool" => tool)
    email !== nothing && (parameters["email"] = email)
    api_key !== nothing && (parameters["api_key"] = api_key)
    fetch_result = _download_text(_entrez_url("efetch.fcgi"; parameters=parameters))
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "entrez_fetch";
        parameters=(db=db, n_ids=length(ids), rettype=rettype, retmode=retmode))
    end
    return fetch_result
end

entrez_fetch_fasta(db::String, ids; kwargs...) = entrez_fetch(db, ids; rettype="fasta", retmode="text", kwargs...)
entrez_fetch_genbank(db::String, ids; kwargs...) = entrez_fetch(db, ids; rettype="gb", retmode="text", kwargs...)
entrez_pubmed_fetch(ids; kwargs...) = entrez_fetch("pubmed", ids; rettype="abstract", retmode="text", kwargs...)

entrez_pubmed_search(term::String; kwargs...) = entrez_search("pubmed", term; kwargs...)
entrez_nucleotide_search(term::String; kwargs...) = entrez_search("nuccore", term; kwargs...)
entrez_protein_search(term::String; kwargs...) = entrez_search("protein", term; kwargs...)
entrez_gene_search(term::String; kwargs...) = entrez_search("gene", term; kwargs...)
entrez_taxonomy_search(term::String; kwargs...) = entrez_search("taxonomy", term; kwargs...)
entrez_genome_search(term::String; kwargs...) = entrez_search("genome", term; kwargs...)
entrez_genome_fetch(ids; kwargs...) = entrez_fetch("genome", ids; kwargs...)
entrez_nuccore_fetch(ids; kwargs...) = entrez_fetch("nuccore", ids; kwargs...)

end

module Medline

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, active_provenance_context, register_provenance!
export MedlineRecord, parse_medline, parse_medline_xml, parse_medline_text, read_medline

struct MedlineRecord
    pmid::String
    title::String
    abstract::String
    journal::String
    year::Union{Nothing,Int}
    authors::Vector{String}
    mesh_terms::Vector{String}
    doi::String
    raw::Dict{String,Vector{String}}
end

@inline function _medline_blocks(text::String)
    stripped = replace(String(text), "\r\n" => "\n", "\r" => "\n")
    blocks = split(strip(stripped), r"\n\s*\n")
    return filter(!isempty, blocks)
end

function _medline_parse_fields(lines::AbstractVector{<:AbstractString})
    fields = Dict{String,Vector{String}}()
    current_key = ""

    for raw_line in lines
        line = rstrip(String(raw_line))
        isempty(line) && continue
        if length(line) >= 6 && line[5:6] == "- "
            current_key = strip(String(line[1:4]))
            value = strip(String(line[7:end]))
            push!(get!(fields, current_key, String[]), value)
        elseif !isempty(current_key)
            push!(fields[current_key], strip(line))
        end
    end

    return fields
end

function _medline_year(text::String)
    match = Base.match(r"(\d{4})", text)
    match === nothing && return nothing
    return parse(Int, match.captures[1])
end

function _medline_join(fields::Dict{String,Vector{String}}, keys::AbstractVector{<:String})
    values = String[]
    for key in keys
        append!(values, get(fields, String(key), String[]))
    end
    isempty(values) && return ""
    return join(values, " ")
end

function _medline_record(fields::Dict{String,Vector{String}})
    pmid = _medline_join(fields, ["PMID"])
    title = _medline_join(fields, ["TI", "TITLE"])
    abstract = _medline_join(fields, ["AB"])
    journal = _medline_join(fields, ["JT", "TA"])
    year = _medline_year(_medline_join(fields, ["DP", "DA", "EDAT"]))
    authors = vcat(get(fields, "FAU", String[]), get(fields, "AU", String[]))
    mesh_terms = get(fields, "MH", String[])
    doi = ""

    for aid in get(fields, "AID", String[])
        if occursin("[doi]", lowercase(aid))
            doi = replace(aid, "[doi]" => "") |> strip
            break
        end
    end

    return MedlineRecord(pmid, title, abstract, journal, year, authors, mesh_terms, doi, fields)
end

function parse_medline_text(text::String)
    records = MedlineRecord[]
    for block in _medline_blocks(text)
        fields = _medline_parse_fields(split(block, '\n'))
        push!(records, _medline_record(fields))
    end
    return records
end

function _medline_xml_blocks(text::AbstractString, tag::String)
    pattern = Regex("<$(tag)(?:\\s[^>]*)?>(.*?)</$(tag)>", "s")
    return [String(match.captures[1]) for match in eachmatch(pattern, String(text))]
end

_medline_xml_first(text::AbstractString, tag::String) = begin
    found = match(Regex("<$(tag)(?:\\s[^>]*)?>(.*?)</$(tag)>", "s"), String(text))
    found === nothing && return ""
    return replace(String(found.captures[1]), "&lt;" => "<", "&gt;" => ">", "&amp;" => "&", "&quot;" => "\"", "&apos;" => "'")
end

function parse_medline_xml(xml::String)
    records = MedlineRecord[]
    xml_text = String(xml)

    for article in _medline_xml_blocks(xml_text, "PubmedArticle")
        medline = _medline_xml_first(article, "MedlineCitation")
        art = _medline_xml_first(article, "Article")

        pmid = _medline_xml_first(medline, "PMID")
        title = _medline_xml_first(art, "ArticleTitle")
        abstract_text = join(_medline_xml_blocks(art, "AbstractText"), " ")
        journal = _medline_xml_first(_medline_xml_first(art, "Journal"), "Title")
        year = nothing
        year_text = _medline_xml_first(article, "PubDate")
        year_match = match(r"(\d{4})", year_text)
        year_match !== nothing && (year = parse(Int, year_match.captures[1]))

        authors = String[]
        for author_block in _medline_xml_blocks(art, "Author")
            last = _medline_xml_first(author_block, "LastName")
            fore = _medline_xml_first(author_block, "ForeName")
            collective = _medline_xml_first(author_block, "CollectiveName")
            if !isempty(collective)
                push!(authors, collective)
            elseif !isempty(last) || !isempty(fore)
                push!(authors, isempty(fore) ? last : "$fore $last")
            end
        end

        mesh_terms = String[]
        for mesh_block in _medline_xml_blocks(medline, "MeshHeading")
            term = _medline_xml_first(mesh_block, "DescriptorName")
            isempty(term) || push!(mesh_terms, term)
        end

        doi = ""
        doi_match = match(Regex("<ArticleId[^>]*?IdType=\"doi\"[^>]*>(.*?)</ArticleId>", "s"), article)
        doi_match !== nothing && (doi = replace(String(doi_match.captures[1]), "&lt;" => "<", "&gt;" => ">", "&amp;" => "&", "&quot;" => "\"", "&apos;" => "'") )

        raw = Dict{String,Vector{String}}("PubmedArticle" => [article])
        push!(records, MedlineRecord(pmid, title, abstract_text, journal, year, authors, mesh_terms, doi, raw))
    end
    return records
end

function parse_medline(text::String)
    records = if startswith(lstrip(text), "<") && occursin("<PubmedArticle", text)
        parse_medline_xml(text)
    else
        parse_medline_text(text)
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "parse_medline";
        parameters=(n_chars=length(text), n_records=length(records)))
    end
    return records
end

function read_medline(path::String)
    records = open(path, "r") do io
        parse_medline(read(io, String))
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "read_medline";
        parameters=(path=path, n_records=length(records)))
    end
    return records
end

end

module Compass

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, BioSequence
export run_needle, run_water, run_transeq, run_revseq, run_compseq, run_seqret

function _compass_run(command::String, args::AbstractVector{<:String}=String[]; input::Union{Nothing,String}=nothing, output::Union{Nothing,String}=nothing)
    cmd_args = String[command]
    append!(cmd_args, String.(args))
    input !== nothing && push!(cmd_args, String(input))
    output !== nothing && push!(cmd_args, String(output))
    cmd = Cmd(cmd_args)
    if output === nothing
        return read(cmd, String)
    end
    run(cmd)
    return String(output)
end

function _temp_fasta(sequence::BioSequence; identifier::String="sequence")
    path, io = mktemp()
    try
        write(io, ">$identifier\n", String(sequence), "\n")
    finally
        close(io)
    end
    return path
end

function _run_pairwise(command::String, first_sequence::BioSequence, second_sequence::BioSequence; args::AbstractVector{<:String}=String[])
    first_path = _temp_fasta(first_sequence; identifier="seq1")
    second_path = _temp_fasta(second_sequence; identifier="seq2")
    output_path = tempname()
    try
        return _compass_run(command, vcat(String.(args), ["-asequence", first_path, "-bsequence", second_path, "-outfile", output_path]); output=output_path)
    finally
        isfile(first_path) && rm(first_path, force=true)
        isfile(second_path) && rm(second_path, force=true)
        isfile(output_path) && rm(output_path, force=true)
    end
end

# Top-level BioSequence dispatch is handled outside the Compass module
# to avoid circular dependencies with the core sequence types.

end

module KEGG

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, active_provenance_context, register_provenance!
export KEGGRecord, KEGGPathwayRecord, KEGGEnzymeRecord, read_kegg_record, read_kegg_pathway, read_kegg_enzyme, kegg_field, kegg_entries, kegg_entry_id
export kegg_pathway_mermaid, write_kegg_pathway_mermaid

struct KEGGRecord
    database::String
    fields::Dict{String,Vector{String}}
end

Base.show(io::IO, record::KEGGRecord) = print(io, "KEGGRecord($(record.database), fields=$(length(record.fields)))")

struct KEGGPathwayRecord
    entry::String
    title::String
    enzymes::Vector{String}
    compounds::Vector{String}
    genes::Vector{String}
    fields::Dict{String,Vector{String}}
end

Base.show(io::IO, record::KEGGPathwayRecord) = print(io, "KEGGPathwayRecord($(record.entry), enzymes=$(length(record.enzymes)), genes=$(length(record.genes)))")

struct KEGGEnzymeRecord
    entry::String
    names::Vector{String}
    reactions::Vector{String}
    pathways::Vector{String}
    fields::Dict{String,Vector{String}}
end

Base.show(io::IO, record::KEGGEnzymeRecord) = print(io, "KEGGEnzymeRecord($(record.entry), names=$(length(record.names)), reactions=$(length(record.reactions)))")

function _kegg_line_key(line::AbstractString)
    length(line) < 12 && return ""
    return strip(String(line[1:12]))
end

function _kegg_line_value(line::AbstractString)
    length(line) <= 12 && return ""
    return strip(String(line[13:end]))
end

function read_kegg_record(io::IO, database::String)
    fields = Dict{String,Vector{String}}()
    current_key = ""
    for raw_line in eachline(io)
        line = rstrip(raw_line)
        line == "///" && break
        key = _kegg_line_key(line)
        value = _kegg_line_value(line)
        if !isempty(key)
            current_key = key
            push!(get!(fields, current_key, String[]), value)
        elseif !isempty(current_key)
            push!(fields[current_key], value)
        end
    end
    kegg_rec = KEGGRecord(String(database), fields)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "read_kegg_record";
        parameters=(database=database, n_fields=length(fields)))
    end
    return kegg_rec
end

read_kegg_record(io::IO; database::String="") = read_kegg_record(io, database)

read_kegg_record(path::String; database::String="") = open(path, "r") do io
    read_kegg_record(io; database=database)
end

kegg_field(record::KEGGRecord, name::String) = get(record.fields, String(name), String[])
kegg_entries(record::KEGGRecord) = keys(record.fields)
kegg_entry_id(record::KEGGRecord) = begin
    entry_fields = kegg_field(record, "ENTRY")
    isempty(entry_fields) && return ""
    return first(split(entry_fields[1]))
end

function _kegg_joined_field(record::KEGGRecord, name::String)
    values = kegg_field(record, name)
    isempty(values) && return ""
    return join(values, " ")
end

function _kegg_first_field(record::KEGGRecord, name::String)
    values = kegg_field(record, name)
    isempty(values) && return ""
    return values[1]
end

function read_kegg_pathway(io::IO)
    record = read_kegg_record(io; database="pathway")
    entry = kegg_entry_id(record)
    title = _kegg_joined_field(record, "NAME")
    enzymes = split(replace(_kegg_joined_field(record, "ENZYME"), "  " => " "))
    compounds = split(replace(_kegg_joined_field(record, "COMPOUND"), "  " => " "))
    genes = kegg_field(record, "GENE")
    pathway_result = KEGGPathwayRecord(String(entry), String(title), filter(!isempty, enzymes), filter(!isempty, compounds), genes, record.fields)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "read_kegg_pathway";
        parameters=(entry=entry, n_enzymes=length(enzymes), n_compounds=length(compounds)))
    end
    return pathway_result
end

read_kegg_pathway(path::String) = open(path, "r") do io
    read_kegg_pathway(io)
end

function read_kegg_enzyme(io::IO)
    record = read_kegg_record(io; database="enzyme")
    entry = _kegg_first_field(record, "ENTRY")
    names = split(_kegg_joined_field(record, "NAME"), ";")
    reactions = kegg_field(record, "REACTION")
    pathways = kegg_field(record, "PATHWAY")
    enzyme_result = KEGGEnzymeRecord(String(entry), filter(!isempty, strip.(names)), reactions, pathways, record.fields)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        register_provenance!(_ctx, "read_kegg_enzyme";
        parameters=(entry=entry, n_reactions=length(reactions), n_pathways=length(pathways)))
    end
    return enzyme_result
end

read_kegg_enzyme(path::String) = open(path, "r") do io
    read_kegg_enzyme(io)
end

end

module Pathway

using ..KEGG: KEGGPathwayRecord
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, active_provenance_context, register_provenance!

export PathwayNode, PathwayEdge, PathwayGraph, read_pathway_graph, pathway_nodes, pathway_edges, kegg_pathway_mermaid, write_kegg_pathway_mermaid

struct PathwayNode
    identifier::String
    label::String
    kind::String
end

Base.show(io::IO, node::PathwayNode) = print(io, "PathwayNode($(node.identifier), kind=$(node.kind))")

struct PathwayEdge
    source::String
    target::String
    relation::String
end

Base.show(io::IO, edge::PathwayEdge) = print(io, "PathwayEdge($(edge.source)->$(edge.target), relation=$(edge.relation))")

struct PathwayGraph
    entry::String
    nodes::Vector{PathwayNode}
    edges::Vector{PathwayEdge}
    metadata::Dict{String,Vector{String}}
end

Base.show(io::IO, graph::PathwayGraph) = print(io, "PathwayGraph($(graph.entry), nodes=$(length(graph.nodes)), edges=$(length(graph.edges)))")

Base.length(graph::PathwayGraph) = length(graph.nodes)

function read_pathway_graph(record::KEGGPathwayRecord)
    nodes = PathwayNode[]
    edges = PathwayEdge[]
    push!(nodes, PathwayNode(record.entry, record.title, "pathway"))
    for gene_line in record.genes
        tokens = split(strip(gene_line), r"\s+", limit=2)
        isempty(tokens) && continue
        identifier = first(tokens)
        label = length(tokens) > 1 ? tokens[2] : identifier
        push!(nodes, PathwayNode(identifier, label, "gene"))
        push!(edges, PathwayEdge(record.entry, identifier, "gene"))
    end
    for enzyme in record.enzymes
        push!(nodes, PathwayNode(enzyme, enzyme, "enzyme"))
        push!(edges, PathwayEdge(record.entry, enzyme, "enzyme"))
    end
    for compound in record.compounds
        push!(nodes, PathwayNode(compound, compound, "compound"))
        push!(edges, PathwayEdge(record.entry, compound, "compound"))
    end
    return PathwayGraph(record.entry, nodes, edges, record.fields)
end

pathway_nodes(graph::PathwayGraph) = graph.nodes
pathway_edges(graph::PathwayGraph) = graph.edges

function _pathway_mermaid_id(label::String)
    cleaned = replace(String(label), r"[^A-Za-z0-9_]" => "_")
    isempty(cleaned) && return "node"
    return cleaned
end

function kegg_pathway_mermaid(graph::PathwayGraph; title::Union{Nothing,String}=nothing)
    io = IOBuffer()
    print(io, "flowchart LR\n")
    if title !== nothing
        print(io, "  %% ", title, "\n")
    end
    for node in graph.nodes
        node_id = _pathway_mermaid_id(node.identifier)
        if node.kind == "pathway"
            print(io, "  ", node_id, "[\"", node.label, "\"]\n")
        elseif node.kind == "enzyme"
            print(io, "  ", node_id, "((\"", node.label, "\"))\n")
        elseif node.kind == "compound"
            print(io, "  ", node_id, "{\"", node.label, "\"}\n")
        else
            print(io, "  ", node_id, "[\"", node.label, "\"]\n")
        end
    end
    for edge in graph.edges
        print(io, "  ", _pathway_mermaid_id(edge.source), " -->|", edge.relation, "| ", _pathway_mermaid_id(edge.target), "\n")
    end
    return String(take!(io))
end

kegg_pathway_mermaid(record::KEGGPathwayRecord; kwargs...) = kegg_pathway_mermaid(read_pathway_graph(record); kwargs...)

function write_kegg_pathway_mermaid(io::IO, graph::PathwayGraph, kwargs...)
    write(io, kegg_pathway_mermaid(graph; kwargs...))
end

write_kegg_pathway_mermaid(path::String, graph::PathwayGraph; kwargs...) = open(path, "w") do io
    write_kegg_pathway_mermaid(io, graph; kwargs...)
end

write_kegg_pathway_mermaid(io::IO, record::KEGGPathwayRecord; kwargs...) = write_kegg_pathway_mermaid(io, read_pathway_graph(record); kwargs...)
write_kegg_pathway_mermaid(path::String, record::KEGGPathwayRecord; kwargs...) = write_kegg_pathway_mermaid(path, read_pathway_graph(record); kwargs...)

end

module SCOP

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, active_provenance_context, register_provenance!
export SCOPRecord, read_scop_records, parse_scop_record

struct SCOPRecord
    sid::String
    pdb_id::String
    chain::String
    residues::String
    class_id::String
    hierarchy::Vector{String}
    description::String
end

Base.show(io::IO, record::SCOPRecord) = print(io, "SCOPRecord($(record.sid), $(record.class_id), levels=$(length(record.hierarchy)))")

function parse_scop_record(line::String)
    fields = split(strip(line), r"\s+"; limit=6)
    length(fields) >= 6 || throw(ArgumentError("SCOP record must have at least 6 fields"))
    class_id = String(fields[4])
    hierarchy = filter(!isempty, split(class_id, '.'))
    return SCOPRecord(String(fields[1]), String(fields[2]), String(fields[3]), String(fields[4]), class_id, hierarchy, join(String.(fields[5:end]), " "))
end

function read_scop_records(io::IO)
    records = SCOPRecord[]
    for raw_line in eachline(io)
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, "#") && continue
        push!(records, parse_scop_record(line))
    end
    return records
end

read_scop_records(path::String) = open(path, "r") do io
    read_scop_records(io)
end

end

module CATH

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, analysis_result_summary, active_provenance_context, register_provenance!
export CATHRecord, read_cath_records, parse_cath_record

struct CATHRecord
    domain::String
    pdb_id::String
    chain::String
    class_id::String
    architecture::String
    topology::String
    homologous_superfamily::String
    hierarchy::Vector{String}
    description::String
end

Base.show(io::IO, record::CATHRecord) = print(io, "CATHRecord($(record.domain), $(record.homologous_superfamily), levels=$(length(record.hierarchy)))")

function parse_cath_record(line::String)
    fields = split(strip(line), r"\s+"; limit=8)
    length(fields) >= 8 || throw(ArgumentError("CATH record must have at least 8 fields"))
    hierarchy = String[String(fields[4]), String(fields[5]), String(fields[6]), String(fields[7])]
    return CATHRecord(String(fields[1]), String(fields[2]), String(fields[3]), String(fields[4]), String(fields[5]), String(fields[6]), String(fields[7]), hierarchy, String(fields[8]))
end

function read_cath_records(io::IO)
    records = CATHRecord[]
    for raw_line in eachline(io)
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, "#") && continue
        push!(records, parse_cath_record(line))
    end
    return records
end

read_cath_records(path::String) = open(path, "r") do io
    read_cath_records(io)
end

end

# ==============================================================================
# Typed BioSequence dispatch for databases functions
#
# Typed DNASeq inputs are automatically converted to String for the
# Restriction module's sequence analysis functions. These are defined
# at the top-level scope since the nested Restriction module cannot
# reference the BioSequence type system.
# ==============================================================================

# Compass module typed dispatch returning BioSequence objects where applicable
export run_needle, run_water, run_transeq, run_revseq, run_compseq, run_seqret

function run_needle(s1::BioSequence, s2::BioSequence, kw...)
    # needle/water output is usually a results file, but if it returns a sequence
    # we would wrap it. For now, alignment results are just returned as String.
    return Compass._run_pairwise("needle", s1, s2; kw...)
end

function run_water(s1::BioSequence, s2::BioSequence, kw...)
    return Compass._run_pairwise("water", s1, s2; kw...)
end

function run_transeq(seq::BioSequence, kw...)
    res = Compass._compass_run("transeq", ["-sequence", Compass._temp_fasta(seq)])
    # transeq output is FASTA format
    records = read_fasta(IOBuffer(res); alphabet=AminoAcidAlphabet)
    return isempty(records) ? nothing : records[1].sequence
end

function run_revseq(seq::BioSequence{A}, kw...) where {A}
    res = Compass._compass_run("revseq", ["-sequence", Compass._temp_fasta(seq)])
    records = read_fasta(IOBuffer(res); alphabet=A)
    return isempty(records) ? nothing : records[1].sequence
end

function run_compseq(seq::BioSequence, kw...)
    return Compass._compass_run("compseq", ["-sequence", Compass._temp_fasta(seq)])
end

function run_seqret(seq::BioSequence{A}, kw...) where {A}
    res = Compass._compass_run("seqret", ["-sequence", Compass._temp_fasta(seq)])
    records = read_fasta(IOBuffer(res); alphabet=A)
    return isempty(records) ? nothing : records[1].sequence
end

"""
    DatabaseDownloadResult

Container for a downloaded database bundle together with its metadata.
"""
struct DatabaseDownloadResult <: AbstractAnalysisResult
    db::String
    source::String
    paths::Vector{String}
    metadata::Dict{Any,Any}
    provenance::ResultProvenance
end

DatabaseDownloadResult(db, source, paths, metadata) =
    DatabaseDownloadResult(db, source, paths, metadata, provenance_record("DatabaseDownloadResult", "databases"))

@inline function _database_download_result(db::AbstractString, source::AbstractString, paths::AbstractVector{<:AbstractString}; metadata=Dict{Any,Any}())
    return DatabaseDownloadResult(String(db), String(source), String.(paths), Dict{Any,Any}(metadata))
end

@inline function _database_download_hash(paths::AbstractVector{<:AbstractString})
    hashes = String[]
    for path in paths
        isfile(path) || continue
        push!(hashes, open(path, "r") do io
            bytes2hex(sha256(read(io)))
        end)
    end
    isempty(hashes) && return nothing
    return bytes2hex(sha256(codeunits(join(hashes, "\n"))))
end

function _finalize_database_download!(ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result::DatabaseDownloadResult, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    ctx === nothing && return result
    register_container_provenance!(ctx, result, operation; parents=parents, parameters=parameters, provenance_hash=_database_download_hash(result.paths))
    return result
end

"""
    download_database(source; download_dir=tempdir(), filename=nothing, fetcher=download, _ctx=nothing)

Download a database payload from a URL or file-like source into a local path.
"""
function download_database(source::AbstractString; download_dir::AbstractString=tempdir(), filename::Union{Nothing,AbstractString}=nothing, fetcher=download, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    mkpath(download_dir)
    cleaned_source = replace(String(source), r"[?#].*$" => "")
    target_name = filename === nothing ? basename(cleaned_source) : String(filename)
    isempty(target_name) && (target_name = "database.dat")
    output_path = joinpath(download_dir, target_name)
    fetcher(String(source), output_path)
    result = _database_download_result("url", source, [output_path]; metadata=Dict{Any,Any}(
        :download_dir => String(download_dir),
        :filename => target_name))
    return _finalize_database_download!(_ctx, result, "download_database"; parameters=(kind=:url, source=String(source), download_dir=String(download_dir), filename=target_name))
end

download_database(kind::Symbol, args...; kwargs...) = download_database(Val(kind), args...; kwargs...)

function download_database(::Val{:tcga}, hits::AbstractVector; base_url::String="https://api.gdc.cancer.gov", download_dir::AbstractString=tempdir(), fetcher=download, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    downloads = tcga_download_files(hits; base_url=base_url, download_dir=String(download_dir), fetcher=fetcher)
    result = _database_download_result("tcga", base_url, downloads.file_paths; metadata=Dict{Any,Any}(
        :base_url => String(base_url),
        :download_dir => String(download_dir),
        :file_count => length(downloads.file_paths),
        :sample_ids => String.(downloads.sample_ids)))
    return _finalize_database_download!(_ctx, result, "download_database"; parameters=(kind=:tcga, base_url=String(base_url), download_dir=String(download_dir), file_count=length(downloads.file_paths)))
end

download_database(::Val{:tcga}, query::NamedTuple; kwargs...) = download_database(:tcga, getproperty(query, :hits); kwargs...)

function download_database(::Val{:entrez}, result::Entrez.EntrezSearchResult; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx), kwargs...)
    return download_database(:entrez, result.db, result.ids; prov_ctx=_ctx, kwargs...)
end

function download_database(::Val{:entrez}, db::AbstractString, ids; rettype::String="fasta", retmode::String="text", download_dir::AbstractString=tempdir(), filename::Union{Nothing,AbstractString}=nothing, fetcher=entrez_fetch, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx), kwargs...)
    mkpath(download_dir)
    id_list = String.(collect(ids))
    target_name = filename === nothing ? string(String(db), "_", rettype, ".", retmode == "text" ? "txt" : retmode) : String(filename)
    payload = fetcher(String(db), id_list; rettype=rettype, retmode=retmode, kwargs...)
    output_path = joinpath(download_dir, target_name)
    open(output_path, "w") do io
        write(io, payload)
    end
    result = _database_download_result("entrez", db, [output_path]; metadata=Dict{Any,Any}(
        :db => String(db),
        :download_dir => String(download_dir),
        :filename => target_name,
        :id_count => length(id_list),
        :ids => id_list,
        :rettype => String(rettype),
        :retmode => String(retmode)))
    return _finalize_database_download!(_ctx, result, "download_database"; parameters=(kind=:entrez, db=String(db), rettype=String(rettype), retmode=String(retmode), download_dir=String(download_dir), id_count=length(id_list)))
end
