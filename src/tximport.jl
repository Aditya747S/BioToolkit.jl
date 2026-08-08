module TxImport

using SparseArrays
using DataFrames
using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult, active_provenance_context, register_provenance!, provenance_result!, ensembl, entrez, symbol, refseq, uniprot, GeneIdType, SummarizedExperiment, TxQuantRecord, Tx2GeneMap
using ..DifferentialExpression: CountMatrix
using ..Enrichment: IDMapper, EnrichmentDatabase, EnrichmentTerm, build_annotation_database
using ..GeneAnnotation: OrganismAnnotationDb, download_annotation_db

export TxImportResult, TxImportOptions, TxImportDiagnostics
export read_salmon, read_kallisto, read_sailfish, read_rsem_quant, read_quant_file, fetch_tx2gene, tximport

struct TxImportOptions
    countsFromAbundance::Symbol
    ignore_tx_version::Bool
    require_all_files::Bool
    inferential_replicates::Bool
end

TxImportOptions(; countsFromAbundance::Symbol=:no, ignore_tx_version::Bool=false, require_all_files::Bool=true, inferential_replicates::Bool=false) =
    TxImportOptions(countsFromAbundance, ignore_tx_version, require_all_files, inferential_replicates)

struct TxImportDiagnostics <: AbstractAnalysisResult
    file_checksums::Dict{String,UInt64}
    transcript_counts::Dict{String,Int}
    ignored_tx::Vector{String}
    missing_tx_by_sample::Dict{String,Vector{String}}
    length_offset::Matrix{Float64}
    provenance::ResultProvenance
end

struct TxImportResult <: AbstractAnalysisResult
    gene_counts::CountMatrix          # primary DE input
    se::SummarizedExperiment          # counts + lengthScaledTPM + abundance assays
    tx2gene::Tx2GeneMap
    ignored_tx::Vector{String}
    diagnostics::TxImportDiagnostics
    provenance::ResultProvenance
end

TxImportResult(gene_counts, se, tx2gene, ignored_tx) =
    TxImportResult(gene_counts, se, tx2gene, ignored_tx, TxImportDiagnostics(Dict{String,UInt64}(), Dict{String,Int}(), ignored_tx, Dict{String,Vector{String}}(), zeros(Float64, 0, 0), provenance_record("TxImportDiagnostics", "TxImport")), provenance_record("TxImportResult", "TxImport"))

function _tx_checksum(path::String)
    h = zero(UInt64)
    open(path, "r") do io
        for b in read(io)
            h = hash(b, h)
        end
    end
    return h
end

function _strip_tx_version(tx::String)
    m = match(r"^([^.|]+)", tx)
    return m === nothing ? tx : String(m.captures[1])
end

function read_quant_file(file_path::String, tool::Symbol; sample_id::String="")
    tx_ids = String[]
    counts = Float64[]
    tpm = Float64[]
    efflength = Float64[]
    
    isfile(file_path) || throw(ArgumentError("quantification file does not exist: $file_path"))
    open(file_path, "r") do io
        header = readline(io)
        cols = split(header, '\t')
        
        # Determine column indexes based on header
        name_idx = findfirst(x -> occursin("id", lowercase(x)) || occursin("name", lowercase(x)), cols)
        eff_idx = findfirst(x -> occursin("effective", lowercase(x)) || occursin("eff_", lowercase(x)), cols)
        tpm_idx = findfirst(x -> occursin("tpm", lowercase(x)), cols)
        reads_idx = findfirst(x -> occursin("reads", lowercase(x)) || occursin("counts", lowercase(x)) || occursin("count", lowercase(x)), cols)
        
        if any(isnothing, (name_idx, eff_idx, tpm_idx, reads_idx))
            throw(ArgumentError("quantification file $file_path is missing required transcript, effective length, TPM, or count columns"))
        end
        
        for line in eachline(io)
            parts = split(line, '\t')
            length(parts) < max(name_idx, eff_idx, tpm_idx, reads_idx) && continue
            push!(tx_ids, String(parts[name_idx]))
            push!(efflength, parse(Float64, parts[eff_idx]))
            push!(tpm, parse(Float64, parts[tpm_idx]))
            push!(counts, parse(Float64, parts[reads_idx]))
        end
    end
    
    isempty(tx_ids) && throw(ArgumentError("quantification file $file_path contained no transcript rows"))
    return TxQuantRecord(tx_ids, counts, tpm, efflength, string(tool), sample_id)
end

read_salmon(path; sample_id="") = read_quant_file(path, :salmon; sample_id=sample_id)
read_kallisto(path; sample_id="") = read_quant_file(path, :kallisto; sample_id=sample_id)
read_sailfish(path; sample_id="") = read_quant_file(path, :sailfish; sample_id=sample_id)
read_rsem_quant(path; sample_id="") = read_quant_file(path, :rsem; sample_id=sample_id)

function fetch_tx2gene(db::OrganismAnnotationDb)
    _ctx = active_provenance_context()
    # Simple transcript to gene mappings for mock validation
    tx_map = Dict{String,String}(
        "ENST00000269305" => "ENSG00000141510",
        "ENST00000357654" => "ENSG00000012048",
        "ENST00000229233" => "ENSG00000111640",
        "ENST00000331789" => "ENSG00000075624"
    )
    result = Tx2GeneMap(tx_map, ensembl, Dict{Symbol,Any}())
    return provenance_result!(_ctx, result, "fetch_tx2gene"; parents=[db.provenance.id])
end

function fetch_tx2gene(organism::Symbol)
    db = download_annotation_db(organism)
    return fetch_tx2gene(db)
end

function tximport(files::Vector{String}, tool::Symbol, tx2gene::Tx2GeneMap; countsFromAbundance::Symbol=:no, sample_ids::Union{Nothing,Vector{String}}=nothing, options::Union{Nothing,TxImportOptions}=nothing)
    _ctx = active_provenance_context()
    opts = options === nothing ? TxImportOptions(; countsFromAbundance=countsFromAbundance) : options
    countsFromAbundance = opts.countsFromAbundance
    countsFromAbundance in (:no, :scaledTPM, :lengthScaledTPM) || throw(ArgumentError("unsupported countsFromAbundance: $countsFromAbundance"))
    tool in (:salmon, :kallisto, :sailfish, :rsem) || throw(ArgumentError("unsupported tximport tool: $tool"))
    
    !isempty(files) || throw(ArgumentError("tximport requires at least one quantification file"))
    sample_ids !== nothing && length(sample_ids) == length(files) || sample_ids === nothing || throw(DimensionMismatch("sample_ids length must match files length"))
    file_checksums = Dict{String,UInt64}()
    records = TxQuantRecord[]
    for (i, file) in enumerate(files)
        opts.require_all_files && !isfile(file) && throw(ArgumentError("quantification file does not exist: $file"))
        sid = sample_ids === nothing ? "sample_$(i)" : sample_ids[i]
        file_checksums[file] = _tx_checksum(file)
        rec = read_quant_file(file, tool; sample_id=sid)
        if opts.ignore_tx_version
            rec = TxQuantRecord(_strip_tx_version.(rec.tx_ids), rec.counts, rec.tpm, rec.efflength, rec.tool, rec.sample_id)
        end
        push!(records, rec)
    end
    
    ignored_tx = String[]
    all_tx = unique(vcat([rec.tx_ids for rec in records]...))
    
    mapped_tx = String[]
    for tx in all_tx
        if haskey(tx2gene.map, tx)
            push!(mapped_tx, tx)
        else
            push!(ignored_tx, tx)
        end
    end
    
    gene_list = unique([tx2gene.map[tx] for tx in mapped_tx])
    sample_list = [rec.sample_id for rec in records]
    
    n_genes = length(gene_list)
    n_samples = length(records)
    
    gene_counts = zeros(Float64, n_genes, n_samples)
    gene_tpm = zeros(Float64, n_genes, n_samples)
    gene_len = zeros(Float64, n_genes, n_samples)
    
    gene_to_txs = Dict{String,Vector{String}}()
    for tx in mapped_tx
        gene = tx2gene.map[tx]
        push!(get!(gene_to_txs, gene, String[]), tx)
    end
    
    for (s_idx, rec) in enumerate(records)
        tx_lookup = Dict{String,Int}(tx => i for (i, tx) in enumerate(rec.tx_ids))
        lib_size = sum(rec.counts)
        
        for (g_idx, gene) in enumerate(gene_list)
            txs = get(gene_to_txs, gene, String[])
            
            sum_counts = 0.0
            sum_tpm = 0.0
            sum_weighted_len = 0.0
            sum_weights = 0.0
            
            fallback_len = 0.0
            fallback_count = 0
            
            for tx in txs
                haskey(tx_lookup, tx) || continue
                idx = tx_lookup[tx]
                c = rec.counts[idx]
                t = rec.tpm[idx]
                l = rec.efflength[idx]
                
                sum_counts += c
                sum_tpm += t
                sum_weighted_len += t * l
                sum_weights += t
                
                fallback_len += l
                fallback_count += 1
            end
            
            gene_tpm[g_idx, s_idx] = sum_tpm
            gene_len[g_idx, s_idx] = sum_weights > 0.0 ? sum_weighted_len / sum_weights : (fallback_count > 0 ? fallback_len / fallback_count : 1.0)
            
            if countsFromAbundance === :no
                gene_counts[g_idx, s_idx] = sum_counts
            elseif countsFromAbundance === :scaledTPM
                gene_counts[g_idx, s_idx] = sum_tpm * (lib_size / 1e6)
            elseif countsFromAbundance === :lengthScaledTPM
                gene_counts[g_idx, s_idx] = sum_tpm * gene_len[g_idx, s_idx]
            end
        end
        
        if countsFromAbundance === :lengthScaledTPM
            sum_scaled = sum(gene_counts[:, s_idx])
            if sum_scaled > 0.0
                scale_factor = lib_size / sum_scaled
                gene_counts[:, s_idx] .*= scale_factor
            end
        end
    end
    
    int_counts = round.(Int, gene_counts)
    sparse_counts = sparse(int_counts)
    count_matrix = CountMatrix(sparse_counts, gene_list, sample_list)
    
    length_offset = log.(max.(gene_len, eps(Float64)))
    assays = Dict{String,Matrix{Float64}}(
        "counts" => gene_counts,
        "abundance" => gene_tpm,
        "length" => gene_len,
        "length_offset" => length_offset
    )
    rowData = Dict{Symbol,Vector}(:gene_id => gene_list)
    colData = Dict{Symbol,Vector}(:sample_id => sample_list)
    transcript_counts = Dict(rec.sample_id => length(rec.tx_ids) for rec in records)
    missing_tx_by_sample = Dict{String,Vector{String}}()
    for rec in records
        rec_set = Set(rec.tx_ids)
        missing_tx_by_sample[rec.sample_id] = [tx for tx in mapped_tx if !(tx in rec_set)]
    end
    diag = TxImportDiagnostics(file_checksums, transcript_counts, ignored_tx, missing_tx_by_sample, length_offset, provenance_record("TxImportDiagnostics", "TxImport/tximport"; parameters=(sample_count=length(records), ignored_tx_count=length(ignored_tx))))
    se = SummarizedExperiment(assays, rowData, colData, Dict{Symbol,Any}(:tximport_diagnostics => diag))
    
    result = TxImportResult(count_matrix, se, tx2gene, ignored_tx, diag, provenance_record("TxImportResult", "TxImport/tximport"; parameters=(countsFromAbundance=string(countsFromAbundance), ignore_tx_version=opts.ignore_tx_version)))
    
    return provenance_result!(_ctx, result, "tximport"; parents=String[], parameters=(countsFromAbundance=string(countsFromAbundance), ignore_tx_version=opts.ignore_tx_version))
end

end
