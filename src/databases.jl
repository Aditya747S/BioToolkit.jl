module Restriction

export RestrictionEnzyme, RestrictionSite, restriction_enzymes, restriction_enzyme, restriction_enzyme_names, restriction_catalog, restriction_sites, restriction_digest_map, find_restriction_sites, digest_sequence

struct RestrictionEnzyme
    name::String
    recognition_site::String
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
    'N' => "[ACGT]",
)

function _restriction_regex(site::AbstractString)
    pattern = IOBuffer()
    print(pattern, "(?i)")
    for character in uppercase(String(site))
        print(pattern, get(_IUPAC_PATTERNS, character, string(character)))
    end
    return Regex(String(take!(pattern)))
end

function _restriction_enzyme(name::AbstractString, recognition_site::AbstractString, cut_offset::Integer)
    return RestrictionEnzyme(String(name), uppercase(String(recognition_site)), Int(cut_offset), _restriction_regex(recognition_site))
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
        _restriction_enzyme("ZraI", "GACGTC", 3),
    )
)

restriction_enzymes() = _RESTRICTION_DATABASE
restriction_catalog() = collect(values(_RESTRICTION_DATABASE))
restriction_enzyme_names() = sort!(collect(keys(_RESTRICTION_DATABASE)))

function restriction_enzyme(name::AbstractString)
    enzyme = get(_RESTRICTION_DATABASE, String(name), nothing)
    enzyme === nothing && throw(ArgumentError("unknown restriction enzyme: $(name)"))
    return enzyme
end

restriction_sites(sequence::AbstractString, enzyme::AbstractString) = find_restriction_sites(sequence, enzyme)
restriction_sites(sequence::AbstractString, enzyme::RestrictionEnzyme) = find_restriction_sites(sequence, enzyme)
restriction_sites(sequence::AbstractString) = find_restriction_sites(sequence)

function find_restriction_sites(sequence::AbstractString, enzyme::RestrictionEnzyme)
    dna = uppercase(String(sequence))
    sites = RestrictionSite[]
    for match in eachmatch(enzyme.regex, dna)
        position = match.offset
        push!(sites, RestrictionSite(enzyme, position, position + enzyme.cut_offset))
    end
    return sites
end

find_restriction_sites(sequence::AbstractString, name::AbstractString) = find_restriction_sites(sequence, restriction_enzyme(name))
find_restriction_sites(sequence::AbstractString) = reduce(vcat, (find_restriction_sites(sequence, enzyme) for enzyme in values(_RESTRICTION_DATABASE)); init=RestrictionSite[])

function digest_sequence(sequence::AbstractString, enzymes::AbstractVector{<:AbstractString})
    hits = RestrictionSite[]
    for enzyme_name in enzymes
        append!(hits, find_restriction_sites(sequence, enzyme_name))
    end
    sort!(hits; by = hit -> hit.cut_position)

    dna = String(sequence)
    length(dna) == 0 && return String[]
    cut_points = Int[0, length(dna)]
    append!(cut_points, (hit.cut_position - 1 for hit in hits if 1 < hit.cut_position <= length(dna)))
    sort!(cut_points)

    fragments = String[]
    for index in 1:length(cut_points)-1
        start_pos = cut_points[index] + 1
        stop_pos = cut_points[index + 1]
        start_pos <= stop_pos || continue
        push!(fragments, dna[start_pos:stop_pos])
    end
    return fragments
end

digest_sequence(sequence::AbstractString, enzyme::AbstractString) = digest_sequence(sequence, [enzyme])
digest_sequence(sequence::AbstractString) = digest_sequence(sequence, collect(keys(_RESTRICTION_DATABASE)))

function restriction_digest_map(sequence::AbstractString, enzymes=restriction_enzyme_names())
    digest = Dict{String,Vector{RestrictionSite}}()
    for enzyme_name in enzymes
        digest[String(enzyme_name)] = find_restriction_sites(sequence, enzyme_name)
    end
    return digest
end

end

module Entrez

using JSON

export EntrezSearchResult, EntrezPostResult, entrez_search, entrez_search_ids, entrez_search_count, entrez_fetch, entrez_fetch_fasta, entrez_fetch_genbank, entrez_summary, entrez_post, entrez_post_ids, entrez_post_webenv, entrez_post_query_key, entrez_link, entrez_elink, entrez_link_ids, entrez_link_linksets, entrez_link_records, entrez_pubmed_search, entrez_nuccore_fetch, entrez_nucleotide_search, entrez_protein_search, entrez_gene_search, entrez_taxonomy_search, entrez_genome_search, entrez_genome_fetch, entrez_pubmed_fetch, parse_entrez_search_response, parse_entrez_post_response

struct EntrezSearchResult
    db::String
    term::String
    ids::Vector{String}
    count::Int
    webenv::Union{Nothing,String}
    query_key::Union{Nothing,Int}
    raw::Dict{String,Any}
end

Base.show(io::IO, result::EntrezSearchResult) = print(io, "EntrezSearchResult(db=$(result.db), term=$(result.term), count=$(result.count), ids=$(length(result.ids)))")

struct EntrezPostResult
    ids::Vector{String}
    webenv::Union{Nothing,String}
    query_key::Union{Nothing,Int}
    raw::Dict{String,Any}
end

Base.show(io::IO, result::EntrezPostResult) = print(io, "EntrezPostResult(ids=$(length(result.ids)), query_key=$(result.query_key))")

function _urlencode(value::AbstractString)
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

function _entrez_url(endpoint::AbstractString; parameters)
    query = join(("$(String(key))=$(_urlencode(string(value)))" for (key, value) in pairs(parameters)), "&")
    return "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/$(endpoint)?$(query)"
end

function _download_text(url::AbstractString)
    path = tempname()
    try
        download(url, path)
        return read(path, String)
    finally
        isfile(path) && rm(path; force=true)
    end
end

function _download_json(url::AbstractString)
    return JSON.parse(_download_text(url))
end

function parse_entrez_search_response(payload::AbstractString)
    data = JSON.parse(payload)
    search = data["esearchresult"]
    ids = String.(search["idlist"])
    count = parse(Int, string(search["count"]))
    query_key = haskey(search, "querykey") ? parse(Int, string(search["querykey"])) : nothing
    webenv = haskey(search, "webenv") ? String(search["webenv"]) : nothing
    return EntrezSearchResult("", "", ids, count, webenv, query_key, data)
end

function parse_entrez_post_response(payload::AbstractString)
    data = JSON.parse(payload)
    result = get(data, "epostresult", Dict{String,Any}())
    ids = String.(get(result, "idlist", String[]))
    query_key = haskey(result, "querykey") ? parse(Int, string(result["querykey"])) : nothing
    webenv = haskey(result, "webenv") ? String(result["webenv"]) : nothing
    return EntrezPostResult(ids, webenv, query_key, data)
end

function entrez_search(db::AbstractString, term::AbstractString; retmax::Integer=20, retstart::Integer=0, email::Union{Nothing,AbstractString}=nothing, tool::AbstractString="BioToolkit", api_key::Union{Nothing,AbstractString}=nothing)
    parameters = Dict{String,Any}(
        "db" => db,
        "term" => term,
        "retmode" => "json",
        "retmax" => retmax,
        "retstart" => retstart,
        "tool" => tool,
    )
    email !== nothing && (parameters["email"] = email)
    api_key !== nothing && (parameters["api_key"] = api_key)
    payload = _download_text(_entrez_url("esearch.fcgi"; parameters=parameters))
    result = parse_entrez_search_response(payload)
    return EntrezSearchResult(String(db), String(term), result.ids, result.count, result.webenv, result.query_key, result.raw)
end

entrez_search_ids(db::AbstractString, term::AbstractString; kwargs...) = entrez_search(db, term; kwargs...).ids
entrez_search_count(db::AbstractString, term::AbstractString; kwargs...) = entrez_search(db, term; kwargs...).count

function entrez_summary(db::AbstractString, ids; email::Union{Nothing,AbstractString}=nothing, tool::AbstractString="BioToolkit", api_key::Union{Nothing,AbstractString}=nothing)
    id_string = join(String.(ids), ",")
    parameters = Dict{String,Any}("db" => db, "id" => id_string, "retmode" => "json", "tool" => tool)
    email !== nothing && (parameters["email"] = email)
    api_key !== nothing && (parameters["api_key"] = api_key)
    return _download_json(_entrez_url("esummary.fcgi"; parameters=parameters))
end

function entrez_post(db::AbstractString, ids; email::Union{Nothing,AbstractString}=nothing, tool::AbstractString="BioToolkit", api_key::Union{Nothing,AbstractString}=nothing)
    id_string = join(String.(ids), ",")
    parameters = Dict{String,Any}("db" => db, "id" => id_string, "retmode" => "json", "tool" => tool)
    email !== nothing && (parameters["email"] = email)
    api_key !== nothing && (parameters["api_key"] = api_key)
    return parse_entrez_post_response(_download_text(_entrez_url("epost.fcgi"; parameters=parameters)))
end

entrez_post_ids(db::AbstractString, ids; kwargs...) = entrez_post(db, ids; kwargs...).ids
entrez_post_webenv(db::AbstractString, ids; kwargs...) = entrez_post(db, ids; kwargs...).webenv
entrez_post_query_key(db::AbstractString, ids; kwargs...) = entrez_post(db, ids; kwargs...).query_key
entrez_post_ids(result::EntrezPostResult) = result.ids
entrez_post_webenv(result::EntrezPostResult) = result.webenv
entrez_post_query_key(result::EntrezPostResult) = result.query_key

function entrez_link(dbfrom::AbstractString, dbto::AbstractString, ids; email::Union{Nothing,AbstractString}=nothing, tool::AbstractString="BioToolkit", api_key::Union{Nothing,AbstractString}=nothing)
    id_string = join(String.(ids), ",")
    parameters = Dict{String,Any}("dbfrom" => dbfrom, "db" => dbto, "id" => id_string, "retmode" => "json", "tool" => tool)
    email !== nothing && (parameters["email"] = email)
    api_key !== nothing && (parameters["api_key"] = api_key)
    return _download_json(_entrez_url("elink.fcgi"; parameters=parameters))
end

entrez_elink(dbfrom::AbstractString, dbto::AbstractString, ids; kwargs...) = entrez_link(dbfrom, dbto, ids; kwargs...)

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

function entrez_link_ids(dbfrom::AbstractString, dbto::AbstractString, ids; kwargs...)
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

function entrez_fetch(db::AbstractString, ids; rettype::AbstractString="fasta", retmode::AbstractString="text", email::Union{Nothing,AbstractString}=nothing, tool::AbstractString="BioToolkit", api_key::Union{Nothing,AbstractString}=nothing)
    id_string = join(String.(ids), ",")
    parameters = Dict{String,Any}("db" => db, "id" => id_string, "rettype" => rettype, "retmode" => retmode, "tool" => tool)
    email !== nothing && (parameters["email"] = email)
    api_key !== nothing && (parameters["api_key"] = api_key)
    return _download_text(_entrez_url("efetch.fcgi"; parameters=parameters))
end

entrez_fetch_fasta(db::AbstractString, ids; kwargs...) = entrez_fetch(db, ids; rettype="fasta", retmode="text", kwargs...)
entrez_fetch_genbank(db::AbstractString, ids; kwargs...) = entrez_fetch(db, ids; rettype="gb", retmode="text", kwargs...)
entrez_pubmed_fetch(ids; kwargs...) = entrez_fetch("pubmed", ids; rettype="abstract", retmode="text", kwargs...)

entrez_pubmed_search(term::AbstractString; kwargs...) = entrez_search("pubmed", term; kwargs...)
entrez_nucleotide_search(term::AbstractString; kwargs...) = entrez_search("nuccore", term; kwargs...)
entrez_protein_search(term::AbstractString; kwargs...) = entrez_search("protein", term; kwargs...)
entrez_gene_search(term::AbstractString; kwargs...) = entrez_search("gene", term; kwargs...)
entrez_taxonomy_search(term::AbstractString; kwargs...) = entrez_search("taxonomy", term; kwargs...)
entrez_genome_search(term::AbstractString; kwargs...) = entrez_search("genome", term; kwargs...)
entrez_genome_fetch(ids; kwargs...) = entrez_fetch("genome", ids; kwargs...)
entrez_nuccore_fetch(ids; kwargs...) = entrez_fetch("nuccore", ids; kwargs...)

end

module Medline

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

@inline function _medline_blocks(text::AbstractString)
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

function _medline_year(text::AbstractString)
    match = Base.match(r"(\d{4})", text)
    match === nothing && return nothing
    return parse(Int, match.captures[1])
end

function _medline_join(fields::Dict{String,Vector{String}}, keys::AbstractVector{<:AbstractString})
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

function parse_medline_text(text::AbstractString)
    records = MedlineRecord[]
    for block in _medline_blocks(text)
        fields = _medline_parse_fields(split(block, '\n'))
        push!(records, _medline_record(fields))
    end
    return records
end

function _medline_xml_blocks(text::AbstractString, tag::AbstractString)
    pattern = Regex("<$(tag)(?:\\s[^>]*)?>(.*?)</$(tag)>", "s")
    return [match.captures[1] for match in eachmatch(pattern, text)]
end

_medline_xml_first(text::AbstractString, tag::AbstractString) = begin
    found = match(Regex("<$(tag)(?:\\s[^>]*)?>(.*?)</$(tag)>", "s"), text)
    found === nothing && return ""
    return replace(String(found.captures[1]), "&lt;" => "<", "&gt;" => ">", "&amp;" => "&", "&quot;" => "\"", "&apos;" => "'")
end

function parse_medline_xml(xml::AbstractString)
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

function parse_medline(text::AbstractString)
    startswith(lstrip(text), "<") && occursin("<PubmedArticle", text) && return parse_medline_xml(text)
    return parse_medline_text(text)
end

function read_medline(path::AbstractString)
    open(path, "r") do io
        return parse_medline(read(io, String))
    end
end

end

module Compass

export run_needle, run_water, run_transeq, run_revseq, run_compseq, run_seqret

function _compass_run(command::AbstractString, args::AbstractVector{<:AbstractString}=String[]; input::Union{Nothing,AbstractString}=nothing, output::Union{Nothing,AbstractString}=nothing)
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

function _temp_fasta(sequence::AbstractString; identifier::AbstractString="sequence")
    path, io = mktemp()
    try
        write(io, ">$identifier\n", String(sequence), "\n")
    finally
        close(io)
    end
    return path
end

function _run_pairwise(command::AbstractString, first_sequence::AbstractString, second_sequence::AbstractString; args::AbstractVector{<:AbstractString}=String[])
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

run_needle(first_sequence::AbstractString, second_sequence::AbstractString; kwargs...) = _run_pairwise("needle", first_sequence, second_sequence)
run_water(first_sequence::AbstractString, second_sequence::AbstractString; kwargs...) = _run_pairwise("water", first_sequence, second_sequence)

function run_transeq(sequence::AbstractString; kwargs...)
    input_path = _temp_fasta(sequence; identifier="sequence")
    try
        return _compass_run("transeq", ["-sequence", input_path])
    finally
        isfile(input_path) && rm(input_path, force=true)
    end
end

function run_revseq(sequence::AbstractString; kwargs...)
    input_path = _temp_fasta(sequence; identifier="sequence")
    try
        return _compass_run("revseq", ["-sequence", input_path])
    finally
        isfile(input_path) && rm(input_path, force=true)
    end
end

function run_compseq(sequence::AbstractString; kwargs...)
    input_path = _temp_fasta(sequence; identifier="sequence")
    try
        return _compass_run("compseq", ["-sequence", input_path])
    finally
        isfile(input_path) && rm(input_path, force=true)
    end
end

function run_seqret(sequence::AbstractString; kwargs...)
    input_path = _temp_fasta(sequence; identifier="sequence")
    try
        return _compass_run("seqret", ["-sequence", input_path])
    finally
        isfile(input_path) && rm(input_path, force=true)
    end
end

end

module KEGG

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

function read_kegg_record(io::IO; database::AbstractString="")
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
    return KEGGRecord(String(database), fields)
end

read_kegg_record(path::AbstractString; database::AbstractString="") = open(path, "r") do io
    read_kegg_record(io; database=database)
end

kegg_field(record::KEGGRecord, name::AbstractString) = get(record.fields, String(name), String[])
kegg_entries(record::KEGGRecord) = keys(record.fields)
kegg_entry_id(record::KEGGRecord) = begin
    entry_fields = kegg_field(record, "ENTRY")
    isempty(entry_fields) && return ""
    return first(split(entry_fields[1]))
end

function _kegg_joined_field(record::KEGGRecord, name::AbstractString)
    values = kegg_field(record, name)
    isempty(values) && return ""
    return join(values, " ")
end

function _kegg_first_field(record::KEGGRecord, name::AbstractString)
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
    return KEGGPathwayRecord(String(entry), String(title), filter(!isempty, enzymes), filter(!isempty, compounds), genes, record.fields)
end

read_kegg_pathway(path::AbstractString) = open(path, "r") do io
    read_kegg_pathway(io)
end

function read_kegg_enzyme(io::IO)
    record = read_kegg_record(io; database="enzyme")
    entry = _kegg_first_field(record, "ENTRY")
    names = split(_kegg_joined_field(record, "NAME"), ";")
    reactions = kegg_field(record, "REACTION")
    pathways = kegg_field(record, "PATHWAY")
    return KEGGEnzymeRecord(String(entry), filter(!isempty, strip.(names)), reactions, pathways, record.fields)
end

read_kegg_enzyme(path::AbstractString) = open(path, "r") do io
    read_kegg_enzyme(io)
end

end

module Pathway

using ..KEGG: KEGGPathwayRecord

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

function _pathway_mermaid_id(label::AbstractString)
    cleaned = replace(String(label), r"[^A-Za-z0-9_]" => "_")
    isempty(cleaned) && return "node"
    return cleaned
end

function kegg_pathway_mermaid(graph::PathwayGraph; title::Union{Nothing,AbstractString}=nothing)
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

function write_kegg_pathway_mermaid(io::IO, graph::PathwayGraph; kwargs...)
    write(io, kegg_pathway_mermaid(graph; kwargs...))
end

write_kegg_pathway_mermaid(path::AbstractString, graph::PathwayGraph; kwargs...) = open(path, "w") do io
    write_kegg_pathway_mermaid(io, graph; kwargs...)
end

write_kegg_pathway_mermaid(io::IO, record::KEGGPathwayRecord; kwargs...) = write_kegg_pathway_mermaid(io, read_pathway_graph(record); kwargs...)
write_kegg_pathway_mermaid(path::AbstractString, record::KEGGPathwayRecord; kwargs...) = write_kegg_pathway_mermaid(path, read_pathway_graph(record); kwargs...)

end

module SCOP

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

function parse_scop_record(line::AbstractString)
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

read_scop_records(path::AbstractString) = open(path, "r") do io
    read_scop_records(io)
end

end

module CATH

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

function parse_cath_record(line::AbstractString)
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

read_cath_records(path::AbstractString) = open(path, "r") do io
    read_cath_records(io)
end

end
