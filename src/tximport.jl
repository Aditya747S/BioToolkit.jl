module TxImport

using SparseArrays
using DataFrames
using Statistics
using SHA
using JSON
using HDF5
using Parquet
using CodecZlib
using Base.Threads

using ..BioToolkit: ResultProvenance, provenance_record, AbstractAnalysisResult,
    active_provenance_context, provenance_result!,
    SummarizedExperiment, TxQuantRecord, Tx2GeneMap
using ..DifferentialExpression: CountMatrix

export TxImportResult, TxImportOptions, TxImportDiagnostics
export read_salmon, read_kallisto, read_sailfish, read_rsem_quant, read_stringtie,
    read_piscem, read_oarfish, read_alevin, read_quant_file,
    fetch_tx2gene, tximport, summarize_to_gene, make_counts_from_abundance
# kept for backwards compatibility with existing callers — these are internal helpers
export _strip_tx_version, _strip_after_bar, _clean_tx_id, _tx_checksum, row_vars,
    _make_counts_from_abundance, _median_length_over_isoform, _replace_missing_length,
    convert_stringtie_counts

# Design invariants:
#  * gene-level and dense tx_out import return TxImportResult
#  * sparse import and type=:alevin return Dict{String,Any} (like R's list return)
#  * inferential replicates are ALWAYS laid out as Vector{Matrix} with one matrix
#    per replicate, each (features × samples). (Per-sample layout is internal only.)
#  * matrices returned by summarize_to_gene have rows ordered by gene_data["unique_genes"]

const CFA_METHODS     = (:no, :scaledTPM, :lengthScaledTPM, :dtuScaledTPM)
const INFREP_TYPES    = (:salmon, :sailfish, :kallisto, :piscem, :oarfish)
const SUPPORTED_TYPES = (:salmon, :sailfish, :kallisto, :rsem, :stringtie,
                         :alevin, :piscem, :oarfish, :none)

_default_sample_ids(n) = [String("sample_$i") for i in 1:n]

# ==============================================================================
# Core Types
# ==============================================================================

struct TxImportOptions
    countsFromAbundance::Symbol
    ignore_tx_version::Bool
    ignore_after_bar::Bool
    require_all_files::Bool
    inferential_replicates::Bool   # true => error (not warn) if replicates are missing
    var_reduce::Bool
    drop_inf_reps::Bool
    sparse::Bool
    sparse_threshold::Float64
    tx_in::Bool
    tx_out::Bool
end

function TxImportOptions(;
    countsFromAbundance::Symbol=:no, ignore_tx_version::Bool=false,
    ignore_after_bar::Bool=false, require_all_files::Bool=true,
    inferential_replicates::Bool=false, var_reduce::Bool=false,
    drop_inf_reps::Bool=false, sparse::Bool=false, sparse_threshold::Float64=1.0,
    tx_in::Bool=true, tx_out::Bool=false)
    return TxImportOptions(countsFromAbundance, ignore_tx_version, ignore_after_bar,
                           require_all_files, inferential_replicates, var_reduce,
                           drop_inf_reps, sparse, sparse_threshold, tx_in, tx_out)
end

struct TxImportDiagnostics <: AbstractAnalysisResult
    file_checksums::Dict{String,UInt64}
    transcript_counts::Dict{String,Int}
    ignored_tx::Vector{String}
    missing_tx_by_sample::Dict{String,Vector{String}}
    length_offset::Matrix{Float64}
    inferential_rep_info::Union{Nothing,Dict{String,Any}}
    provenance::ResultProvenance
end

struct TxImportResult <: AbstractAnalysisResult
    gene_counts::Union{CountMatrix,Nothing}
    transcript_counts::Union{Matrix{Float64},Nothing}
    abundance::Matrix{Float64}
    length::Matrix{Float64}
    counts_from_abundance::Symbol
    se::Union{SummarizedExperiment,Nothing}
    tx2gene::Union{Tx2GeneMap,Nothing}
    ignored_tx::Vector{String}
    inf_reps::Union{Vector{Matrix{Float64}},Nothing}   # one (features × samples) matrix per replicate
    variance::Union{Matrix{Float64},Nothing}
    diagnostics::TxImportDiagnostics
    provenance::ResultProvenance
end

struct InfRepData
    vars::Vector{Float64}
    reps::Matrix{Float64}
end

# ==============================================================================
# Utilities
# ==============================================================================

# SHA-256 truncated to a UInt64: stable across Julia sessions (unlike hash()).
function _tx_checksum(path::AbstractString)
    open(path, "r") do io
        return reinterpret(UInt64, sha2_256(io))[1]
    end
end

function _strip_tx_version(tx::String)
    m = match(r"^([^.|]+)", tx)
    return m === nothing ? tx : String(m.captures[1])
end

function _strip_after_bar(tx::String)
    m = match(r"^([^|]+)", tx)
    return m === nothing ? tx : String(m.captures[1])
end

function _clean_tx_id(tx::String, ignore_tx_version::Bool, ignore_after_bar::Bool)
    ignore_tx_version && (tx = _strip_tx_version(tx))
    ignore_after_bar  && (tx = _strip_after_bar(tx))
    return tx
end

function _clean_tx2gene_dict(tx2gene::Tx2GeneMap, ignore_tx_version::Bool, ignore_after_bar::Bool)
    d = Dict{String,String}()
    for (k, v) in tx2gene.map
        ck = _clean_tx_id(k, ignore_tx_version, ignore_after_bar)
        if haskey(d, ck) && d[ck] != v
            @warn "tx2gene: cleaned transcript id '$ck' maps to multiple genes ('$(d[ck])' and '$v'); keeping the first"
        else
            d[ck] = v
        end
    end
    return d
end

row_vars(x::AbstractMatrix{<:Real}) =
    size(x, 2) <= 1 ? zeros(Float64, size(x, 1)) : vec(var(x, dims=2))

_read_lines(path::AbstractString) = String.(filter(!isempty, strip.(readlines(path))))

_open_maybe_gz(f, body) = open(f, "r") do io
    body(endswith(f, ".gz") ? GzipDecompressorStream(io) : io)
end

# TPM from counts and effective lengths (guarding zero lengths), used by the
# kallisto h5 reader, the oarfish reader and the inf_rep_stat path.
function _tpm_from_counts(counts::AbstractVector{<:Real}, efflens::AbstractVector{<:Real})
    length(counts) == length(efflens) ||
        throw(DimensionMismatch("counts ($(length(counts))) and efflens ($(length(efflens))) length mismatch"))
    denom = 0.0
    for i in eachindex(counts)
        efflens[i] > 0 && (denom += counts[i] / efflens[i])
    end
    denom > 0 || return zeros(Float64, length(counts))
    out = Vector{Float64}(undef, length(counts))
    for i in eachindex(counts)
        out[i] = efflens[i] > 0 ? (counts[i] / efflens[i]) * (1e6 / denom) : 0.0
    end
    return out
end

# Rescale each column of `m` so it sums to the corresponding column sum of `ref`.
function _rescale_to_library!(m::AbstractMatrix, ref::AbstractMatrix)
    size(ref, 2) == size(m, 2) || throw(DimensionMismatch("column count mismatch"))
    for j in axes(m, 2)
        s = sum(@view m[:, j])
        s > 0 && (m[:, j] .*= sum(@view ref[:, j]) / s)
    end
    return m
end

# ==============================================================================
# Quantification File Readers
# ==============================================================================

function _detect_columns(header::String, tool::Symbol)
    cols = split(header, '\t')
    cols_lower = lowercase.(cols)

    name_idx  = findfirst(x -> occursin("id", x) || occursin("name", x), cols_lower)
    eff_idx   = findfirst(x -> occursin("effective", x) || occursin("eff_", x) || occursin("length", x), cols_lower)
    tpm_idx   = findfirst(x -> occursin("tpm", x) || occursin("fpkm", x), cols_lower)
    reads_idx = findfirst(x -> occursin("reads", x) || occursin("count", x), cols_lower)

    if tool == :kallisto
        name_idx  = findfirst(==("target_id"), cols)
        eff_idx   = findfirst(==("eff_length"), cols)
        tpm_idx   = findfirst(==("tpm"), cols)
        reads_idx = findfirst(==("est_counts"), cols)
    elseif tool == :rsem
        name_idx  = "transcript_id" in cols ? findfirst(==("transcript_id"), cols) :
                                            findfirst(==("gene_id"), cols)
        eff_idx   = findfirst(==("effective_length"), cols)
        tpm_idx   = findfirst(==("TPM"), cols)
        reads_idx = findfirst(==("expected_count"), cols)
    elseif tool == :stringtie
        name_idx  = findfirst(==("t_name"), cols)
        tpm_idx   = findfirst(==("FPKM"), cols)
        reads_idx = findfirst(==("cov"), cols)
        eff_idx   = findfirst(==("length"), cols)
    elseif tool == :piscem
        name_idx  = findfirst(==("target_name"), cols)
        eff_idx   = findfirst(==("eelen"), cols)
        tpm_idx   = findfirst(==("tpm"), cols)
        reads_idx = findfirst(==("ecount"), cols)
    elseif tool == :oarfish
        name_idx  = findfirst(==("tname"), cols)
        eff_idx   = findfirst(==("len"), cols)
        reads_idx = findfirst(==("num_reads"), cols)
        tpm_idx   = reads_idx   # oarfish has no abundance column; TPM is derived below
    elseif tool == :salmon || tool == :sailfish
        name_idx  = findfirst(==("Name"), cols)
        eff_idx   = findfirst(==("EffectiveLength"), cols)
        tpm_idx   = findfirst(==("TPM"), cols)
        reads_idx = findfirst(==("NumReads"), cols)
    end

    return name_idx, eff_idx, tpm_idx, reads_idx, cols
end

function _resolve_col(cols, spec, what, file)
    if spec isa Integer
        1 <= spec <= length(cols) ||
            throw(ArgumentError("$what column index $spec out of range in $file (has $(length(cols)) columns)"))
        return Int(spec)
    end
    i = findfirst(==(String(spec)), cols)
    i === nothing && throw(ArgumentError("$what column '$spec' not found in $file (columns: $(join(cols, ", ")))"))
    return i
end

function _parse_quant_body(io::IO, i_id::Int, i_len::Int, i_tpm::Int, i_cnt::Int,
                           file_path::String, sample_id::String, tool::String,
                           ignore_tx_version::Bool, ignore_after_bar::Bool)
    tx_ids = String[]; counts = Float64[]; tpm = Float64[]; efflength = Float64[]
    maxc = max(i_id, i_len, i_tpm, i_cnt)
    lineno = 1
    for line in eachline(io)
        lineno += 1
        line = rstrip(line, '\r')                     # CRLF safety
        (isempty(line) || first(line) == '#') && continue
        parts = split(line, '\t')
        length(parts) >= maxc ||
            throw(ArgumentError("$file_path: line $lineno has $(length(parts)) fields, expected ≥ $maxc"))
        push!(tx_ids, _clean_tx_id(String(parts[i_id]), ignore_tx_version, ignore_after_bar))
        try
            push!(efflength, parse(Float64, parts[i_len]))
            push!(tpm,       parse(Float64, parts[i_tpm]))
            push!(counts,    parse(Float64, parts[i_cnt]))
        catch
            throw(ArgumentError("$file_path: line $lineno contains a non-numeric field: $line"))
        end
    end
    isempty(tx_ids) && throw(ArgumentError("quantification file $file_path contained no transcript rows"))
    return TxQuantRecord(tx_ids, counts, tpm, efflength, tool, sample_id)
end

function read_quant_file(file_path::AbstractString, tool::Symbol;
                         sample_id::String="", ignore_tx_version::Bool=false,
                         ignore_after_bar::Bool=false)
    isfile(file_path) || throw(ArgumentError("quantification file does not exist: $file_path"))
    open(file_path, "r") do io
        header = replace(rstrip(readline(io), '\r'), "\ufeff" => "")   # strip CR + UTF-8 BOM
        i_id, i_len, i_tpm, i_cnt, cols = _detect_columns(header, tool)
        if any(isnothing, (i_id, i_len, i_tpm, i_cnt))
            throw(ArgumentError("quantification file $file_path is missing required columns; found: $(join(cols, ", "))"))
        end
        rec = _parse_quant_body(io, i_id, i_len, i_tpm, i_cnt, file_path,
                                sample_id, string(tool), ignore_tx_version, ignore_after_bar)
        if tool == :oarfish
            # oarfish reports read counts only; derive TPM from counts / effective length
            rec = TxQuantRecord(rec.tx_ids, rec.counts, _tpm_from_counts(rec.counts, rec.efflength),
                                rec.efflength, rec.tool, rec.sample_id)
        end
        return rec
    end
end

# Reader honoring explicit column names or 1-based indices (R-style colSpec).
function _read_quant_by_cols(file_path::AbstractString; id_col, tpm_col, cnt_col, len_col,
                             sample_id::String="", tool::String="custom",
                             ignore_tx_version::Bool=false, ignore_after_bar::Bool=false)
    open(file_path, "r") do io
        header = replace(rstrip(readline(io), '\r'), "\ufeff" => "")
        cols = split(header, '\t')
        i_id  = _resolve_col(cols, id_col,  "id",        file_path)
        i_tpm = _resolve_col(cols, tpm_col, "abundance", file_path)
        i_cnt = _resolve_col(cols, cnt_col, "counts",    file_path)
        i_len = _resolve_col(cols, len_col, "length",    file_path)
        return _parse_quant_body(io, i_id, i_len, i_tpm, i_cnt, file_path,
                                 sample_id, tool, ignore_tx_version, ignore_after_bar)
    end
end

function read_kallisto_h5(file_path::AbstractString; sample_id::String="",
                          ignore_tx_version::Bool=false, ignore_after_bar::Bool=false)
    isfile(file_path) || throw(ArgumentError("kallisto HDF5 file does not exist: $file_path"))
    h5open(file_path, "r") do fid
        counts  = Float64.(read(fid, "est_counts"))
        ids     = [String(id) for id in read(fid, "aux/ids")]
        efflens = Float64.(read(fid, "aux/eff_lengths"))
        tx_ids  = _clean_tx_id.(ids, ignore_tx_version, ignore_after_bar)
        return TxQuantRecord(tx_ids, counts, _tpm_from_counts(counts, efflens), efflens,
                             "kallisto", sample_id)
    end
end

_is_kallisto_h5(f) = endswith(f, ".h5") || basename(f) == "abundance.h5"

for (fn, tool) in (:read_salmon => :salmon, :read_sailfish => :sailfish,
                   :read_rsem_quant => :rsem, :read_stringtie => :stringtie,
                   :read_piscem => :piscem, :read_oarfish => :oarfish)
    @eval $fn(path; sample_id="", ignore_tx_version=false, ignore_after_bar=false) =
        read_quant_file(path, $(QuoteNode(tool)); sample_id=sample_id,
                        ignore_tx_version=ignore_tx_version, ignore_after_bar=ignore_after_bar)
end

read_kallisto(path; sample_id="", ignore_tx_version=false, ignore_after_bar=false) =
    _is_kallisto_h5(path) ?
        read_kallisto_h5(path; sample_id=sample_id, ignore_tx_version=ignore_tx_version,
                         ignore_after_bar=ignore_after_bar) :
        read_quant_file(path, :kallisto; sample_id=sample_id,
                        ignore_tx_version=ignore_tx_version, ignore_after_bar=ignore_after_bar)

# ==============================================================================
# Inferential Replicates Readers
# ==============================================================================

function read_inf_rep_fish(dir::AbstractString, tool::Symbol)
    tool in (:salmon, :sailfish) ||
        throw(ArgumentError("read_inf_rep_fish expects :salmon or :sailfish"))
    json_path = joinpath(dir, "cmd_info.json")
    isfile(json_path) || return nothing
    cmd_info  = JSON.parsefile(json_path)
    aux_path  = joinpath(dir, get(cmd_info, "auxDir", tool == :salmon ? "aux_info" : "aux"))
    isdir(aux_path) || return nothing
    meta_path = joinpath(aux_path, "meta_info.json")
    isfile(meta_path) || return nothing
    minfo = JSON.parsefile(meta_path)

    for (key, minv) in ("salmon_version" => v"0.8.0", "sailfish_version" => v"0.9.0")
        if haskey(minfo, key)
            v = tryparse(VersionNumber, string(minfo[key]))
            (v === nothing || v < minv) && return nothing
        end
    end

    num_boot    = Int(get(minfo, "num_bootstraps", 0))
    num_boot > 0 || return nothing
    num_targets = Int(get(minfo, "num_valid_targets", get(minfo, "num_targets", 0)))
    num_targets > 0 || return nothing

    boot_path = joinpath(aux_path, "bootstrap", "bootstraps.gz")
    isfile(boot_path) || return nothing

    n_expected = num_targets * num_boot
    boots = open(boot_path, "r") do io
        gz  = GzipDecompressorStream(io)
        buf = Vector{Float64}(undef, n_expected)
        nb  = readbytes!(gz, reinterpret(UInt8, buf), n_expected * sizeof(Float64))
        if nb == n_expected * sizeof(Float64)
            buf
        else
            # older salmon wrote Int32 bootstraps; restart the stream cleanly
            seekstart(io)
            gz   = GzipDecompressorStream(io)
            ints = Vector{Int32}(undef, n_expected)
            nb   = readbytes!(gz, reinterpret(UInt8, ints), n_expected * sizeof(Int32))
            nb == n_expected * sizeof(Int32) ||
                throw(ArgumentError("unexpected bootstrap file size in $boot_path (expected $n_expected values)"))
            Float64.(ints)
        end
    end

    # bootstraps are stored target-major: column b of the reshape is bootstrap b
    reps = reshape(boots, num_targets, num_boot)
    return InfRepData(row_vars(reps), reps)
end

function read_inf_rep_kallisto(dir::AbstractString)
    h5_path = joinpath(dir, "abundance.h5")
    isfile(h5_path) || return nothing
    h5open(h5_path, "r") do fid
        haskey(fid, "bootstrap") || return nothing
        g = fid["bootstrap"]
        # count the bs* datasets inside the group, sorted numerically
        # (lexical order would give bs0, bs1, bs10, bs11, bs2, ...)
        boot_names = sort!(collect(String.(keys(g)));
                           by = n -> (m = match(r"\d+$", n); m === nothing ? 0 : parse(Int, m.match)))
        isempty(boot_names) && return nothing
        tx_count = length(read(fid, "aux/ids"))
        reps = Matrix{Float64}(undef, tx_count, length(boot_names))
        for (i, n) in enumerate(boot_names)
            reps[:, i] = Float64.(read(g[n]))
        end
        return InfRepData(row_vars(reps), reps)
    end
end

function read_inf_rep_piscem(file_path::AbstractString)
    # inferential replicates live as infreps.pq next to the quant file.
    # (Derived by filename join, NOT by regex-replacing "quant" — the old regex
    #  never matched "quant.sf" and handed the quant file itself to the parser.)
    d = dirname(abspath(file_path))
    path = nothing
    for cand in ("infreps.pq", "infreps.parquet")
        p = joinpath(d, cand)
        isfile(p) && (path = p; break)
    end
    path === nothing && return nothing

    df = DataFrame(Parquet.read_parquet(path))
    boot_cols = [c for c in names(df) if occursin(r"^bootstrap", c)]
    isempty(boot_cols) && return nothing
    sort!(boot_cols; by = c -> (m = match(r"(\d+)$", c); m === nothing ? 0 : parse(Int, m[1])))

    reps = Matrix{Float64}(undef, nrow(df), length(boot_cols))
    for (i, c) in enumerate(boot_cols)
        reps[:, i] = Float64.(df[!, c])
    end
    return InfRepData(row_vars(reps), reps)
end

_read_inf_reps_for(type::Symbol, file::AbstractString) =
    if type in (:piscem, :oarfish)
        read_inf_rep_piscem(file)
    elseif type in (:salmon, :sailfish)
        read_inf_rep_fish(dirname(file), type)
    elseif type == :kallisto
        read_inf_rep_kallisto(dirname(file))
    else
        nothing
    end

# ==============================================================================
# Alevin (EDS format)
# ==============================================================================

"""
    read_alevin(dir; filter_barcodes=false, tier_import=false, drop_mean_var=false, drop_inf_reps=false)

Read an alevin run directory into sparse gene × cell matrices. All returned
matrices stay sparse (scRNA data must not be densified).

EDS layout per cell: `ceil(n_genes/8)` bytes of presence bits (LSB first),
followed by the counts of the expressed genes as Float32 (see fishpond::readEDS,
which reads `size=4`). Reading them as Float64 desynchronizes the whole stream.
"""
struct AlevinData
    counts::SparseMatrixCSC{Float64,Int}
    tier::Union{SparseMatrixCSC{Float64,Int},Nothing}
    mean::Union{SparseMatrixCSC{Float64,Int},Nothing}
    variance::Union{SparseMatrixCSC{Float64,Int},Nothing}
    inf_reps::Union{Vector{SparseMatrixCSC{Float64,Int}},Nothing}
    gene_names::Vector{String}
    cell_names::Vector{String}
end

function _read_eds_matrix(io::IO, n_genes::Int, cell_names::AbstractVector{String})
    n_cells = length(cell_names)
    nbytes  = cld(n_genes, 8)
    bytes   = Vector{UInt8}(undef, nbytes)
    I = Int[]; J = Int[]; V = Float64[]
    for j in 1:n_cells
        nread = readbytes!(io, bytes, nbytes)
        nread == nbytes ||
            throw(ArgumentError("EDS stream truncated in cell $j (got $nread of $nbytes bit bytes)"))
        nexp = sum(count_ones, bytes)
        vals = nexp > 0 ? read(io, Float32, nexp) : Float32[]
        k = 0
        for (bi, b) in enumerate(bytes), bit in 0:7
            gene = (bi - 1) * 8 + bit + 1
            gene > n_genes && continue           # padding bits in the final byte
            if (b >> bit) & 0x01 == 0x01
                k += 1
                push!(I, gene); push!(J, j); push!(V, Float64(vals[k]))
            end
        end
        k == nexp || throw(ArgumentError("EDS bit/counts mismatch in cell $j (got $k values for $nexp set bits)"))
    end
    return sparse(I, J, V, n_genes, n_cells)
end

function _read_eds_file(path::AbstractString, gene_names::AbstractVector{String},
                        cell_names::AbstractVector{String}; matrices::Int=1)
    open(path, "r") do io
        gz = GzipDecompressorStream(io)
        n_genes = length(gene_names)
        return [_read_eds_matrix(gz, n_genes, cell_names) for _ in 1:matrices]
    end
end

function read_alevin(dir::AbstractString; filter_barcodes::Bool=false,
                     tier_import::Bool=false, drop_mean_var::Bool=false,
                     drop_inf_reps::Bool=false)
    barcode_file = joinpath(dir, "alevin/quants_mat_rows.txt")
    gene_file    = joinpath(dir, "alevin/quants_mat_cols.txt")
    matrix_file  = joinpath(dir, "alevin/quants_mat.gz")
    tier_file    = joinpath(dir, "alevin/quants_tier_mat.gz")
    mean_file    = joinpath(dir, "alevin/quants_mean_mat.gz")
    var_file     = joinpath(dir, "alevin/quants_var_mat.gz")
    boot_file    = joinpath(dir, "alevin/quants_boot_mat.gz")
    boot_rows    = joinpath(dir, "alevin/quants_boot_rows.txt")
    whitelist    = joinpath(dir, "alevin/whitelist.txt")

    for f in (barcode_file, gene_file, matrix_file)
        isfile(f) || throw(ArgumentError("alevin file not found: $f (pass the path to quants_mat.gz; alevin_dir is its grandparent)"))
    end

    cell_names = _read_lines(barcode_file)
    gene_names = _read_lines(gene_file)

    cmd_info  = let j = joinpath(dir, "cmd_info.json")
        isfile(j) ? JSON.parsefile(j) : Dict{String,Any}()
    end
    num_boot  = Int(get(cmd_info, "numCellBootstraps", 0))

    keep = trues(length(cell_names))
    if filter_barcodes
        if isfile(whitelist)
            wl = Set(_read_lines(whitelist))
            keep = [c in wl for c in cell_names]
            @info "alevin: keeping $(sum(keep)) of $(length(cell_names)) cell barcodes from whitelist"
        else
            @warn "filter_barcodes=true but no whitelist at $whitelist; keeping all cells"
        end
    end

    counts = only(_read_eds_file(matrix_file, gene_names, cell_names))[:, keep]

    tier = nothing
    if tier_import
        isfile(tier_file) || throw(ArgumentError("tier_import=true but tier file not found: $tier_file"))
        tier = only(_read_eds_file(tier_file, gene_names, cell_names))[:, keep]
    end

    mean_mat = nothing
    var_mat  = nothing
    if num_boot > 0 && !drop_mean_var
        if isfile(mean_file) && isfile(var_file)
            mean_mat = only(_read_eds_file(mean_file, gene_names, cell_names))[:, keep]
            var_mat  = only(_read_eds_file(var_file,  gene_names, cell_names))[:, keep]
        else
            @warn "alevin: mean/variance matrices not found; skipping"
        end
    end

    inf_reps = nothing
    if num_boot > 0 && !drop_inf_reps && isfile(boot_file)
        # CRITICAL: the bootstrap file holds `num_boot` consecutive EDS matrices in
        # one stream. They must be read sequentially — re-opening the file per
        # replicate yields num_boot identical copies.
        boot_cells = isfile(boot_rows) ? _read_lines(boot_rows) : cell_names
        mats = _read_eds_file(boot_file, gene_names, boot_cells; matrices=num_boot)
        inf_reps = [m[:, keep] for m in mats]
    end

    return AlevinData(counts, tier, mean_mat, var_mat, inf_reps,
                      gene_names, cell_names[keep])
end

# ==============================================================================
# countsFromAbundance
# ==============================================================================

"""
    make_counts_from_abundance(counts, abundance, length; method=:scaledTPM) -> Matrix{Float64}

Recompute counts from abundance. Never mutates its inputs.

* `:scaledTPM`       — new counts are TPM, rescaled so each column sums to that
                       sample's library size (column sums of `counts`).
* `:lengthScaledTPM` — TPM × row-mean length, rescaled to library size.

`:dtuScaledTPM` is defined at the transcript level and is applied inside
`summarize_to_gene` (or via `_make_counts_from_abundance` with an explicit
`gene_lengths` per-row scale).
"""
function make_counts_from_abundance(counts::AbstractMatrix{<:Real},
                                    abundance::AbstractMatrix{<:Real},
                                    length_mat::AbstractMatrix{<:Real};
                                    method::Symbol=:scaledTPM)
    method in (:scaledTPM, :lengthScaledTPM) ||
        throw(ArgumentError("make_counts_from_abundance supports :scaledTPM and :lengthScaledTPM (got :$method)"))
    size(counts) == size(abundance) == size(length_mat) ||
        throw(DimensionMismatch("counts, abundance and length must have identical dimensions"))
    new_counts = method == :scaledTPM ? Float64.(abundance) :      # copy — no aliasing
                 abundance .* vec(mean(length_mat, dims=2))
    return _rescale_to_library!(new_counts, counts)
end

# Backwards-compatible wrapper (also the dtuScaledTPM entry point).
function _make_counts_from_abundance(counts_mat::AbstractMatrix{<:Real},
                                     abundance_mat::AbstractMatrix{<:Real},
                                     length_mat::AbstractMatrix{<:Real},
                                     method::Symbol; gene_lengths=nothing)
    if method == :dtuScaledTPM
        gene_lengths === nothing &&
            throw(ArgumentError("dtuScaledTPM requires gene_lengths (per-row median isoform length)"))
        size(abundance_mat, 1) == length(gene_lengths) ||
            throw(DimensionMismatch("gene_lengths has $(length(gene_lengths)) entries, abundance has $(size(abundance_mat,1)) rows"))
        new_counts = abundance_mat .* gene_lengths
        return _rescale_to_library!(new_counts, counts_mat)
    end
    return make_counts_from_abundance(counts_mat, abundance_mat, length_mat; method=method)
end

"""
    _median_length_over_isoform(length_mat, tx_ids, gene_ids) -> Vector{Float64}

Per transcript: mean length across samples, then the median of that value over
all isoforms of the same gene. Rows of `length_mat`, `tx_ids` and `gene_ids`
must correspond (all already cleaned/subset consistently).
"""
function _median_length_over_isoform(length_mat::AbstractMatrix{<:Real},
                                     tx_ids::AbstractVector{String},
                                     gene_ids::AbstractVector{String})
    (size(length_mat, 1) == length(tx_ids) == length(gene_ids)) ||
        throw(DimensionMismatch("length_mat rows, tx_ids and gene_ids must have equal length"))
    ave = vec(mean(length_mat, dims=2))
    per_gene = Dict{String,Vector{Float64}}()
    for (g, l) in zip(gene_ids, ave)
        push!(get!(per_gene, g, Float64[]), l)
    end
    med = Dict(g => median(v) for (g, v) in per_gene)
    return [med[g] for g in gene_ids]
end

# Mutates in place (matches R's .replaceMissingLength). Fully-missing rows get
# `ave_length`; partially-missing rows get the geometric mean of observed entries.
function _replace_missing_length(length_mat::Matrix{Float64}, ave_length::Real)
    for i in axes(length_mat, 1)
        row = @view length_mat[i, :]
        n_missing = count(isnan, row)
        n_missing == 0 && continue
        if n_missing == length(row)
            fill!(row, ave_length)
        else
            s = 0.0; n = 0
            for v in row
                (!isnan(v) && v > 0) && (s += log(v); n += 1)
            end
            fill_v = n > 0 ? exp(s / n) : Float64(ave_length)
            for j in eachindex(row)
                isnan(row[j]) && (row[j] = fill_v)
            end
        end
    end
    return length_mat
end

convert_stringtie_counts(counts::Matrix{Float64}, length::Matrix{Float64},
                         read_length::Float64=75.0) = counts .* length ./ read_length

# ==============================================================================
# Gene-Level Summarization
# ==============================================================================

"""
    summarize_to_gene(tx_data, tx2gene, tx_ids; ...) -> Dict{String,Any}

Sum transcript-level matrices to gene level. All gene-level sums are computed
via a single sparse incidence-matrix product (no per-transcript allocations).
`gene_data["unique_genes"]` gives the row order of every returned matrix —
always use it for rownames; do not re-derive gene ids from `tx2gene`.

Inferential replicates must be `Vector{Matrix}` — one (n_tx × n_samples) matrix
per replicate. `var_reduce=true` replaces them with a per-(gene, sample)
variance (sample variance, mean-subtracted, n-1 denominator).
"""
function summarize_to_gene(tx_data::Dict{String,<:Any}, tx2gene::Tx2GeneMap,
                           tx_ids::Vector{String};
                           var_reduce::Bool=false, ignore_tx_version::Bool=false,
                           ignore_after_bar::Bool=false,
                           counts_from_abundance::Symbol=:no,
                           inf_reps::Union{Vector{Matrix{Float64}},Nothing}=nothing)
    abundance_tx = tx_data["abundance"]
    counts_tx    = tx_data["counts"]
    length_tx    = tx_data["length"]

    counts_from_abundance in CFA_METHODS ||
        throw(ArgumentError("unsupported countsFromAbundance: $counts_from_abundance"))

    n_tx = length(tx_ids)
    size(abundance_tx, 1) == n_tx || throw(DimensionMismatch("abundance has $(size(abundance_tx,1)) rows, expected $n_tx"))
    size(counts_tx, 1)    == n_tx || throw(DimensionMismatch("counts has $(size(counts_tx,1)) rows, expected $n_tx"))
    size(length_tx, 1)    == n_tx || throw(DimensionMismatch("length has $(size(length_tx,1)) rows, expected $n_tx"))
    n_samples = size(abundance_tx, 2)
    (size(counts_tx, 2) == n_samples && size(length_tx, 2) == n_samples) ||
        throw(DimensionMismatch("abundance/counts/length sample counts disagree"))

    clean_ids = _clean_tx_id.(tx_ids, ignore_tx_version, ignore_after_bar)
    clean_map = _clean_tx2gene_dict(tx2gene, ignore_tx_version, ignore_after_bar)

    valid = [haskey(clean_map, t) for t in clean_ids]
    n_drop = count(!, valid)
    n_drop > 0 && @warn "summarize_to_gene: $n_drop of $n_tx quantified transcript(s) not in tx2gene; dropped"

    vidx = findall(valid)
    isempty(vidx) && throw(ArgumentError(
        "none of the quantified transcripts match tx2gene (check ignore_tx_version / ignore_after_bar / the mapping)"))

    clean_ids    = clean_ids[vidx]
    abundance_tx = abundance_tx[vidx, :]
    counts_tx    = counts_tx[vidx, :]
    length_tx    = length_tx[vidx, :]
    inf_reps     = inf_reps === nothing ? nothing : [m[vidx, :] for m in inf_reps]

    gene_ids     = [clean_map[t] for t in clean_ids]
    unique_genes = unique(gene_ids)                       # deterministic: first-appearance order
    gidx         = Dict(g => i for (i, g) in enumerate(unique_genes))
    n_genes      = length(unique_genes)

    # tx → gene incidence; every isoform sum becomes one sparse×dense product
    St = sparse([gidx[g] for g in gene_ids], collect(1:length(gene_ids)), 1.0,
                n_genes, length(gene_ids))

    abundance_gene = St * abundance_tx
    weighted       = St * (abundance_tx .* length_tx)

    counts_gene = if counts_from_abundance == :dtuScaledTPM
        # dtuScaledTPM is defined at transcript level (TPM × median isoform
        # length, rescaled to the tx-level library size) and THEN summed to
        # gene level — which is why the result does not sum to library size.
        med = _median_length_over_isoform(length_tx, clean_ids, gene_ids)
        new_counts = abundance_tx .* med
        _rescale_to_library!(new_counts, counts_tx)
        St * new_counts
    else
        St * counts_tx
    end

    # gene length = abundance-weighted mean tx length; NaN where abundance is 0
    length_gene = weighted ./ ifelse.(abundance_gene .> 0, abundance_gene, NaN)
    observed    = [v for v in length_gene if !isnan(v)]
    ave         = isempty(observed) ? (isempty(length_tx) ? 1000.0 : mean(length_tx)) : mean(observed)
    length_gene = _replace_missing_length(length_gene, ave)

    if counts_from_abundance in (:scaledTPM, :lengthScaledTPM)
        counts_gene = make_counts_from_abundance(counts_gene, abundance_gene, length_gene;
                                                 method=counts_from_abundance)
    end

    inf_reps_gene  = nothing
    variance_gene  = nothing
    if inf_reps !== nothing && !isempty(inf_reps)
        inf_reps_gene = [St * m for m in inf_reps]
        if var_reduce
            n_boot = length(inf_reps_gene)
            n_boot >= 2 || throw(ArgumentError("var_reduce requires ≥ 2 inferential replicates (got $n_boot)"))
            means = zeros(n_genes, n_samples)
            for m in inf_reps_gene; means .+= m; end
            means ./= n_boot
            variance_gene = zeros(n_genes, n_samples)
            for m in inf_reps_gene
                @. variance_gene += (m - means)^2      # fused: no temp allocations
            end
            variance_gene ./= (n_boot - 1)
            inf_reps_gene = nothing
        end
    end

    return Dict{String,Any}(
        "abundance"            => abundance_gene,
        "counts"               => counts_gene,
        "length"               => length_gene,
        "inf_reps"             => inf_reps_gene,
        "variance"             => variance_gene,
        "counts_from_abundance" => counts_from_abundance,
        "unique_genes"         => unique_genes,   # rownames of all matrices above
    )
end

# ==============================================================================
# tx2gene acquisition
# ==============================================================================

"""
    fetch_tx2gene(path; tx_col=1, gene_col=2, delim='\\t', ignore_tx_version=false) -> Tx2GeneMap

Build a `Tx2GeneMap` from a GTF/GFF file (`gene_id` / `transcript_id` attributes)
or a two-column table (plain or gzipped). Table columns are auto-detected from a
header if one is present; pass 1-based `tx_col`/`gene_col` to override.
"""
function fetch_tx2gene(path::AbstractString; tx_col::Int=1, gene_col::Int=2,
                       delim::Char='\t', ignore_tx_version::Bool=false)
    isfile(path) || throw(ArgumentError("tx2gene source does not exist: $path"))
    b = lowercase(basename(path))
    if occursin(".gtf", b) || occursin(".gff", b)
        return _tx2gene_from_gtf(path; ignore_tx_version=ignore_tx_version)
    end
    return _tx2gene_from_table(path; tx_col=tx_col, gene_col=gene_col,
                               delim=delim, ignore_tx_version=ignore_tx_version)
end

function _tx2gene_from_gtf(path::AbstractString; ignore_tx_version::Bool)
    m = Dict{String,String}()
    _open_maybe_gz(path) do io
        for line in eachline(io)
            (isempty(line) || first(line) == '#') && continue
            occursin("transcript_id", line) || continue
            g = match(r"gene_id\s+\"?([^\";]+)\"?", line)
            t = match(r"transcript_id\s+\"?([^\";]+)\"?", line)
            (g === nothing || t === nothing) && continue
            tx   = _clean_tx_id(String(strip(t[1])), ignore_tx_version, false)
            gene = String(strip(g[1]))
            if haskey(m, tx) && m[tx] != gene
                @warn "fetch_tx2gene: transcript '$tx' maps to multiple genes ('$(m[tx])' and '$gene'); keeping the first"
            else
                m[tx] = gene
            end
        end
    end
    isempty(m) && throw(ArgumentError("no gene_id/transcript_id pairs found in $path (is it a GTF/GFF?)"))
    return Tx2GeneMap(m)   # assumes Tx2GeneMap(Dict{String,String}) constructor
end

function _tx2gene_from_table(path::AbstractString; tx_col::Int, gene_col::Int,
                             delim::Char, ignore_tx_version::Bool)
    m = Dict{String,String}()
    _open_maybe_gz(path) do io
        first_line = true
        for line in eachline(io)
            line = rstrip(line, '\r')
            (isempty(line) || first(line) == '#') && continue
            parts = split(line, delim)
            if first_line
                first_line = false
                ti = findfirst(p -> occursin("transcript", lowercase(p)), parts)
                gi = findfirst(p -> occursin("gene", lowercase(p)), parts)
                if ti !== nothing && gi !== nothing
                    tx_col, gene_col = ti, gi
                    continue
                end
            end
            length(parts) >= max(tx_col, gene_col) ||
                throw(ArgumentError("malformed tx2gene row (need ≥ $(max(tx_col, gene_col)) fields): $line"))
            tx = _clean_tx_id(String(strip(parts[tx_col])), ignore_tx_version, false)
            m[tx] = String(strip(parts[gene_col]))
        end
    end
    isempty(m) && throw(ArgumentError("no transcript→gene rows found in $path"))
    return Tx2GeneMap(m)
end

# ==============================================================================
# RSEM gene-level import (tx_in = false)
# ==============================================================================

function _compute_rsem_gene_level(files::Vector{String};
                                  gene_id_col=nothing, abundance_col=nothing,
                                  counts_col=nothing, length_col=nothing,
                                  counts_from_abundance::Symbol=:no, ctx=nothing,
                                  sample_ids::Union{Vector{String},Nothing}=nothing,
                                  ignore_tx_version::Bool=false, ignore_after_bar::Bool=false)
    counts_from_abundance != :no &&
        @warn "countsFromAbundance='$counts_from_abundance' requires transcript-level estimates; returning RSEM gene-level expected counts"

    n_files     = length(files)
    sample_list = sample_ids === nothing ? _default_sample_ids(n_files) : sample_ids

    read_one = (f, sid) -> begin
        if gene_id_col === nothing && abundance_col === nothing &&
           counts_col === nothing && length_col === nothing
            read_quant_file(f, :rsem; sample_id=sid, ignore_tx_version=ignore_tx_version,
                            ignore_after_bar=ignore_after_bar)
        else
            _read_quant_by_cols(f; id_col=gene_id_col, tpm_col=abundance_col,
                                cnt_col=counts_col, len_col=length_col,
                                sample_id=sid, tool="rsem",
                                ignore_tx_version=ignore_tx_version,
                                ignore_after_bar=ignore_after_bar)
        end
    end

    first_rec = read_one(files[1], sample_list[1])
    gene_ids  = first_rec.tx_ids
    n_genes   = length(gene_ids)

    abundance = Matrix{Float64}(undef, n_genes, n_files)
    counts    = Matrix{Float64}(undef, n_genes, n_files)
    lengths   = Matrix{Float64}(undef, n_genes, n_files)
    file_gene_ids  = Vector{Vector{String}}(undef, n_files)
    file_checksums = Dict{String,UInt64}()

    for j in 1:n_files
        rec = read_one(files[j], sample_list[j])
        rec.tx_ids == gene_ids ||
            throw(ArgumentError("gene IDs in $(files[j]) do not match $(files[1])"))
        file_checksums[files[j]] = _tx_checksum(files[j])
        file_gene_ids[j] = rec.tx_ids
        abundance[:, j]  = rec.tpm
        counts[:, j]     = rec.counts
        lengths[:, j]    = max.(rec.efflength, 1.0)
    end

    count_matrix  = CountMatrix(sparse(round.(Int, counts)), gene_ids, sample_list)
    length_offset = log.(max.(lengths, eps(Float64)))

    assays = Dict{String,Matrix{Float64}}("counts" => counts, "abundance" => abundance,
                                          "length" => lengths, "length_offset" => length_offset)
    transcript_counts = Dict{String,Int}(sample_list[j] => length(file_gene_ids[j]) for j in 1:n_files)
    diag = TxImportDiagnostics(file_checksums, transcript_counts, String[],
                               Dict{String,Vector{String}}(), length_offset, nothing,
                               provenance_record("TxImportDiagnostics", "TxImport/rsem-gene-level";
                                                 parameters=(sample_count=n_files, gene_count=n_genes)))
    se = SummarizedExperiment(assays, Dict{Symbol,Vector}(:gene_id => gene_ids),
                              Dict{Symbol,Vector}(:sample_id => sample_list),
                              Dict{Symbol,Any}(:tximport_diagnostics => diag,
                                               :provenance => diag.provenance.id))
    result = TxImportResult(count_matrix, nothing, abundance, lengths, :no, se, nothing,
                            String[], nothing, nothing, diag,
                            provenance_record("TxImportResult", "TxImport/rsem-gene-level";
                                              parameters=(type="rsem", tx_in=false, gene_count=n_genes)))
    return provenance_result!(ctx, result, "tximport"; parents=String[];
                              parameters=(type="rsem", tx_in=false))
end

# ==============================================================================
# Main tximport
# ==============================================================================

"""
    tximport(files; type=:salmon, ...) -> TxImportResult or Dict{String,Any}

Import transcript-level quantification files.

Returns:
* gene level (default) and dense `tx_out=true`: a `TxImportResult` (with a
  `SummarizedExperiment` and full diagnostics);
* `sparse=true` or `type=:alevin`: a `Dict{String,Any}` with sparse matrices
  (mirroring R's list return; matrices are NOT densified).

Inferential replicate layout in the result: `Vector{Matrix}` with one
(features × samples) matrix per replicate.
"""
function tximport(files::Vector{String};
                  type::Symbol=:salmon, tx_in::Bool=true, tx_out::Bool=false,
                  counts_from_abundance::Symbol=:no,
                  tx2gene::Union{Tx2GeneMap,Nothing}=nothing,
                  var_reduce::Bool=false, drop_inf_reps::Bool=false,
                  inf_rep_stat::Union{Function,Nothing}=nothing,
                  ignore_tx_version::Bool=false, ignore_after_bar::Bool=false,
                  gene_id_col=nothing, tx_id_col=nothing, abundance_col=nothing,
                  counts_col=nothing, length_col=nothing,
                  importer::Union{Function,Nothing}=nothing,
                  existence_optional::Bool=false, require_all_files::Bool=true,
                  sparse::Bool=false, sparse_threshold::Float64=1.0,
                  read_length::Float64=75.0,
                  alevin_args::Union{Dict{Symbol,Any},Nothing}=nothing,
                  sample_ids::Union{Vector{String},Nothing}=nothing,
                  options::Union{TxImportOptions,Nothing}=nothing)

    ctx = active_provenance_context()

    if options !== nothing
        counts_from_abundance = options.countsFromAbundance
        ignore_tx_version     = options.ignore_tx_version
        ignore_after_bar      = options.ignore_after_bar
        var_reduce            = options.var_reduce
        drop_inf_reps         = options.drop_inf_reps
        sparse                = options.sparse
        sparse_threshold      = options.sparse_threshold
        tx_in                 = options.tx_in
        tx_out                = options.tx_out
        require_all_files     = options.require_all_files
    end
    want_inf_reps = options !== nothing && options.inferential_replicates

    # ---- validation ----------------------------------------------------------
    counts_from_abundance in CFA_METHODS ||
        throw(ArgumentError("unsupported countsFromAbundance: $counts_from_abundance (expected one of $CFA_METHODS)"))
    type in SUPPORTED_TYPES ||
        throw(ArgumentError("unsupported tximport type: $type (expected one of $SUPPORTED_TYPES)"))
    (type != :none || importer !== nothing) ||
        throw(ArgumentError("type=:none requires a custom importer function"))
    isempty(files) && throw(ArgumentError("tximport requires at least one quantification file"))
    if sample_ids !== nothing
        length(sample_ids) == length(files) ||
            throw(DimensionMismatch("sample_ids ($(length(sample_ids))) must match files ($(length(files)))"))
        allunique(sample_ids) || throw(ArgumentError("sample_ids must be unique"))
    end
    if any(!isnothing, (tx_id_col, abundance_col, counts_col, length_col))
        all(!isnothing, (tx_id_col, abundance_col, counts_col, length_col)) ||
            throw(ArgumentError("custom columns require all of tx_id_col, abundance_col, counts_col, length_col"))
    end

    # ---- alevin: single experiment, already gene-level -----------------------
    if type == :alevin
        length(files) == 1 ||
            throw(ArgumentError("alevin import supports a single experiment (pass the path to quants_mat.gz)"))
        alevin_dir = dirname(dirname(abspath(files[1])))
        args = alevin_args === nothing ? Dict{Symbol,Any}() : alevin_args
        data = read_alevin(alevin_dir;
                           filter_barcodes = get(args, :filterBarcodes, false),
                           tier_import    = get(args, :tierImport, false),
                           drop_mean_var  = get(args, :dropMeanVar, false),
                           drop_inf_reps  = drop_inf_reps)
        # alevin quantifies at gene level; tx2gene/tx_out are not applicable (R does the same)
        result = Dict{String,Any}(
            "counts" => data.counts, "abundance" => nothing, "length" => nothing,
            "counts_from_abundance" => "no",
            "gene_names" => data.gene_names, "cell_names" => data.cell_names)
        data.inf_reps !== nothing && (result["inf_reps"] = data.inf_reps)
        data.variance !== nothing && (result["variance"] = data.variance)
        data.mean    !== nothing && (result["mean"] = data.mean)
        data.tier    !== nothing && (result["tier"] = data.tier)
        return provenance_result!(ctx, result, "tximport/alevin"; parents=String[],
                                  parameters=(type="alevin",))
    end

    # ---- file existence / require_all_files ----------------------------------
    if !existence_optional
        missing_files = [f for f in files if !isfile(f)]
        if !isempty(missing_files)
            require_all_files && throw(ArgumentError("quantification file does not exist: $(missing_files[1])"))
            @warn "require_all_files=false: skipping $(length(missing_files)) missing file(s)"
            keep = [isfile(f) for f in files]
            files = files[keep]
            sample_ids === nothing || (sample_ids = sample_ids[keep])
            isempty(files) && throw(ArgumentError("no quantification files remain after dropping missing files"))
        end
    end

    n_files     = length(files)
    sample_list = sample_ids === nothing ? _default_sample_ids(n_files) : sample_ids

    # ---- RSEM gene-level detection --------------------------------------------
    if type == :rsem && tx_in
        n_gene_files = count(f -> occursin("genes", basename(f)), files)
        if n_gene_files == n_files
            @info "RSEM genes.results detected; importing at gene level (tx_in=false)"
            tx_in = false
        elseif n_gene_files > 0
            throw(ArgumentError("mixing RSEM isoforms.results and genes.results files is not supported"))
        end
    end
    if !tx_in
        type == :rsem || throw(ArgumentError("tx_in=false is only supported for type=:rsem"))
        return _compute_rsem_gene_level(files; gene_id_col=gene_id_col,
                                        abundance_col=abundance_col, counts_col=counts_col,
                                        length_col=length_col,
                                        counts_from_abundance=counts_from_abundance,
                                        ctx=ctx, sample_ids=sample_list,
                                        ignore_tx_version=ignore_tx_version,
                                        ignore_after_bar=ignore_after_bar)
    end

    !tx_out && tx2gene === nothing &&
        throw(ArgumentError("tx2gene is required for gene-level summarization (tx_out=false)"))

    # ---- inferential replicate configuration ----------------------------------
    use_inf_reps  = !drop_inf_reps && type in INFREP_TYPES
    read_var_only = use_inf_reps && var_reduce && tx_out   # variance-only is meaningful at tx level

    if sparse
        tx_out || throw(ArgumentError("sparse import requires tx_out=true"))
        use_inf_reps && throw(ArgumentError("sparse import does not support inferential replicates"))
        counts_from_abundance in (:no, :scaledTPM) ||
            throw(ArgumentError("sparse import supports only countsFromAbundance=:no or :scaledTPM"))
        type == :stringtie && throw(ArgumentError(
            "sparse import does not support stringtie (counts derive from coverage × length)"))
    end
    inf_rep_stat !== nothing && sparse &&
        throw(ArgumentError("inf_rep_stat is incompatible with sparse import"))
    inf_rep_stat !== nothing && read_var_only &&
        throw(ArgumentError("inf_rep_stat requires full inferential replicates (incompatible with var_reduce && tx_out)"))

    # ---- record reader (per-file kallisto h5 detection) ------------------------
    if type == :kallisto && !any(_is_kallisto_h5, files)
        @info "kallisto: importing abundance.h5 is typically faster than abundance.tsv"
    end

    _read_record = (f::String, sid::String) -> begin
        if importer !== nothing
            rec = importer(f)
            rec isa TxQuantRecord ||
                throw(ArgumentError("importer($(repr(f))) must return a TxQuantRecord"))
            rec
        elseif any(!isnothing, (tx_id_col, abundance_col, counts_col, length_col))
            _read_quant_by_cols(f; id_col=tx_id_col, tpm_col=abundance_col,
                                cnt_col=counts_col, len_col=length_col,
                                sample_id=sid, tool=string(type),
                                ignore_tx_version=ignore_tx_version,
                                ignore_after_bar=ignore_after_bar)
        elseif type == :kallisto && _is_kallisto_h5(f)
            read_kallisto_h5(f; sample_id=sid, ignore_tx_version=ignore_tx_version,
                             ignore_after_bar=ignore_after_bar)
        else
            read_quant_file(f, type; sample_id=sid, ignore_tx_version=ignore_tx_version,
                            ignore_after_bar=ignore_after_bar)
        end
    end

    # ---- probe the first file for shape + inferential replicates ---------------
    first_rec = _read_record(files[1], sample_list[1])
    tx_ids    = first_rec.tx_ids
    n_tx      = length(tx_ids)

    if use_inf_reps
        first_inf = _read_inf_reps_for(type, files[1])
        if first_inf === nothing
            want_inf_reps && throw(ArgumentError(
                "inferential_replicates=true but none found for $(files[1])"))
            @warn "no inferential replicates found for $(files[1]); proceeding without them"
            use_inf_reps  = false
            read_var_only = false
        end
    end
    inf_rep_stat !== nothing && !use_inf_reps &&
        @warn "inf_rep_stat was provided but inferential replicates are not available; ignoring it"

    # ---- sparse import ----------------------------------------------------------
    if sparse
        @info "sparse import (tx_out, countsFromAbundance=$counts_from_abundance)"
        counts_I = Int[]; counts_J = Int[]; counts_V = Float64[]
        tpm_I = Int[];    tpm_J = Int[];    tpm_V = Float64[]
        file_checksums  = Dict{String,UInt64}()
        for (j, f) in enumerate(files)
            rec = _read_record(f, sample_list[j])
            rec.tx_ids == tx_ids ||
                throw(ArgumentError("transcript IDs in $f do not match $(files[1])"))
            file_checksums[f] = _tx_checksum(f)
            # scaledTPM must actually be applied here (per column), not just stored as TPM
            col = counts_from_abundance == :scaledTPM ?
                  (let s = sum(rec.tpm); s > 0 ? rec.tpm .* (sum(rec.counts) / s) : copy(rec.tpm); end) :
                  rec.counts
            idx = findall(>=(sparse_threshold), col)
            append!(counts_I, idx)
            append!(counts_J, fill(j, length(idx)))
            append!(counts_V, col[idx])
            if counts_from_abundance == :scaledTPM
                append!(tpm_I, idx)
                append!(tpm_J, fill(j, length(idx)))
                append!(tpm_V, rec.tpm[idx])
            end
        end
        result = Dict{String,Any}(
            "counts"     => sparse(counts_I, counts_J, counts_V, n_tx, n_files),
            "abundance"  => counts_from_abundance == :scaledTPM ?
                            sparse(tpm_I, tpm_J, tpm_V, n_tx, n_files) :
                            spzeros(Float64, n_tx, n_files),
            "length"     => spzeros(Float64, n_tx, n_files),   # lengths are not stored in sparse mode
            "tx_ids"     => tx_ids,
            "samples"    => sample_list,
            "counts_from_abundance" => string(counts_from_abundance),
            "file_checksums" => file_checksums)
        return provenance_result!(ctx, result, "tximport"; parents=String[];
                                  parameters=(type=string(type), tx_out=true, sparse=true,
                                              counts_from_abundance=string(counts_from_abundance)))
    end

    # ---- dense import (threaded; all shared state written per-index) ------------
    abundance_tx = Matrix{Float64}(undef, n_tx, n_files)
    counts_tx    = Matrix{Float64}(undef, n_tx, n_files)
    length_tx    = Matrix{Float64}(undef, n_tx, n_files)

    # NOTE: Vector{Bool}, not BitVector — bits share 64-bit words, so concurrent
    # setindex! from @threads is a data race.
    have_inf       = use_inf_reps ? fill(false, n_files) : nothing
    per_sample_inf = (use_inf_reps && !read_var_only) ?
                     Vector{Union{Matrix{Float64},Nothing}}(nothing, n_files) : nothing
    vars_tx        = read_var_only ? Matrix{Float64}(undef, n_tx, n_files) : nothing
    file_tx_ids    = Vector{Vector{String}}(undef, n_files)
    checksums      = Vector{UInt64}(undef, n_files)   # per-index: Dict writes are not thread-safe

    @threads for j in 1:n_files
        rec = _read_record(files[j], sample_list[j])
        rec.tx_ids == tx_ids ||
            throw(ArgumentError("transcript IDs in $(files[j]) do not match $(files[1])"))
        checksums[j]   = _tx_checksum(files[j])
        file_tx_ids[j] = rec.tx_ids
        abundance_tx[:, j] = rec.tpm
        counts_tx[:, j]    = rec.counts
        # RSEM: clamp 0-length transcripts BEFORE any use of the lengths
        # (covers both the main path and the inf_rep_stat TPM recomputation)
        length_tx[:, j] = type == :rsem ? max.(rec.efflength, 1.0) : rec.efflength

        if use_inf_reps
            rd = _read_inf_reps_for(type, files[j])
            if rd !== nothing
                have_inf[j] = true
                if read_var_only
                    vars_tx[:, j] = rd.vars
                else
                    per_sample_inf[j] = rd.reps
                end
                if inf_rep_stat !== nothing
                    stat_counts = inf_rep_stat(rd.reps)
                    (stat_counts isa AbstractVector{<:Real} && length(stat_counts) == n_tx) ||
                        throw(ArgumentError("inf_rep_stat must return a vector of length $n_tx"))
                    counts_tx[:, j]    = Float64.(stat_counts)
                    abundance_tx[:, j] = _tpm_from_counts(stat_counts, @view length_tx[:, j])
                end
            end
        end
    end

    if use_inf_reps && !all(have_inf)
        n_missing = count(!, have_inf)
        want_inf_reps && throw(ArgumentError(
            "inferential replicates missing for $n_missing of $n_files samples"))
        @warn "inferential replicates missing for $n_missing of $n_files samples; dropping them"
        use_inf_reps  = false
        read_var_only = false
    end

    file_checksums = Dict{String,UInt64}(files[j] => checksums[j] for j in 1:n_files)

    # reorganize per-sample (n_tx × n_boot) into per-replicate (n_tx × n_samples)
    inf_reps_tx = nothing
    if use_inf_reps && !read_var_only
        n_boot = size(per_sample_inf[1], 2)
        all(m -> size(m, 2) == n_boot, per_sample_inf) ||
            throw(ArgumentError("inconsistent number of bootstrap replicates across samples"))
        inf_reps_tx = Vector{Matrix{Float64}}(undef, n_boot)
        for b in 1:n_boot
            M = Matrix{Float64}(undef, n_tx, n_files)
            for j in 1:n_files
                m = per_sample_inf[j]
                M[:, j] = @view m[:, b]
            end
            inf_reps_tx[b] = M
        end
    end

    tx_data = Dict{String,Any}("abundance" => abundance_tx, "counts" => counts_tx,
                               "length" => length_tx)
    inf_reps_tx !== nothing && (tx_data["inf_reps"] = inf_reps_tx)
    (read_var_only && use_inf_reps) && (tx_data["variance"] = vars_tx)

    if type == :stringtie
        tx_data["counts"] = convert_stringtie_counts(tx_data["counts"], tx_data["length"], read_length)
    end

    # ---- transcript-level output -------------------------------------------------
    if tx_out
        tx_counts = tx_data["counts"]
        if counts_from_abundance != :no
            if counts_from_abundance == :dtuScaledTPM
                tx2gene === nothing && throw(ArgumentError("dtuScaledTPM requires tx2gene"))
                clean_map = _clean_tx2gene_dict(tx2gene, ignore_tx_version, ignore_after_bar)
                unmapped  = [t for t in tx_ids if !haskey(clean_map, t)]
                isempty(unmapped) || throw(ArgumentError(
                    "dtuScaledTPM: $(length(unmapped)) transcript(s) not present in tx2gene"))
                med = _median_length_over_isoform(tx_data["length"], tx_ids,
                                                   [clean_map[t] for t in tx_ids])
                tx_counts = _make_counts_from_abundance(tx_counts, tx_data["abundance"],
                                                        tx_data["length"], :dtuScaledTPM;
                                                        gene_lengths=med)
            else
                tx_counts = make_counts_from_abundance(tx_counts, tx_data["abundance"],
                                                       tx_data["length"]; method=counts_from_abundance)
            end
        end

        length_offset = log.(max.(tx_data["length"], eps(Float64)))
        clean_map = tx2gene === nothing ? nothing :
                    _clean_tx2gene_dict(tx2gene, ignore_tx_version, ignore_after_bar)
        ignored_tx = clean_map === nothing ? String[] : [t for t in tx_ids if !haskey(clean_map, t)]
        mapped_tx  = clean_map === nothing ? tx_ids : [t for t in tx_ids if haskey(clean_map, t)]

        transcript_counts = Dict{String,Int}(sample_list[j] => length(file_tx_ids[j]) for j in 1:n_files)
        missing_tx_by_sample = Dict{String,Vector{String}}()
        for j in 1:n_files
            s = Set(file_tx_ids[j])
            missing_tx_by_sample[sample_list[j]] = [t for t in mapped_tx if !(t in s)]
        end
        inf_info = haskey(tx_data, "inf_reps") ? Dict("n_inf_reps" => length(tx_data["inf_reps"])) :
                   haskey(tx_data, "variance") ? Dict("variance_only" => true) : nothing

        diag = TxImportDiagnostics(file_checksums, transcript_counts, ignored_tx,
                                   missing_tx_by_sample, length_offset, inf_info,
                                   provenance_record("TxImportDiagnostics", "TxImport/tximport";
                                       parameters=(sample_count=n_files, tx_out=true,
                                                   ignored_tx_count=length(ignored_tx))))

        assays = Dict{String,Matrix{Float64}}("counts" => tx_counts,
                                              "abundance" => tx_data["abundance"],
                                              "length" => tx_data["length"],
                                              "length_offset" => length_offset)
        haskey(tx_data, "variance") && (assays["variance"] = tx_data["variance"])
        if haskey(tx_data, "inf_reps")
            for (i, m) in enumerate(tx_data["inf_reps"])
                assays["inf_rep_$i"] = m
            end
        end
        se = SummarizedExperiment(assays, Dict{Symbol,Vector}(:transcript_id => tx_ids),
                                  Dict{Symbol,Vector}(:sample_id => sample_list),
                                  Dict{Symbol,Any}(:tximport_diagnostics => diag,
                                                   :provenance => diag.provenance.id))
        result = TxImportResult(nothing, tx_counts, tx_data["abundance"], tx_data["length"],
                                counts_from_abundance, se, tx2gene, ignored_tx,
                                get(tx_data, "inf_reps", nothing),
                                get(tx_data, "variance", nothing), diag,
                                provenance_record("TxImportResult", "TxImport/tximport";
                                    parameters=(type=string(type), tx_out=true,
                                                counts_from_abundance=string(counts_from_abundance),
                                                ignore_tx_version=ignore_tx_version,
                                                ignore_after_bar=ignore_after_bar,
                                                var_reduce=var_reduce)))
        return provenance_result!(ctx, result, "tximport"; parents=String[];
                                  parameters=(type=string(type), tx_out=true,
                                              counts_from_abundance=string(counts_from_abundance)))
    end

    # ---- gene-level summarization ---------------------------------------------------
    gene_data = summarize_to_gene(tx_data, tx2gene, tx_ids;
                                  var_reduce=var_reduce,
                                  ignore_tx_version=ignore_tx_version,
                                  ignore_after_bar=ignore_after_bar,
                                  counts_from_abundance=counts_from_abundance,
                                  inf_reps=get(tx_data, "inf_reps", nothing))

    gene_list = gene_data["unique_genes"]      # single source of truth for row order
    count_matrix = CountMatrix(sparse(round.(Int, gene_data["counts"])), gene_list, sample_list)
    # (fractional counts from countsFromAbundance are preserved in assays["counts"];
    #  CountMatrix rounds for integer-count DE front-ends, like DESeq2 does in R)
    length_offset = log.(max.(gene_data["length"], eps(Float64)))

    assays = Dict{String,Matrix{Float64}}("counts" => gene_data["counts"],
                                          "abundance" => gene_data["abundance"],
                                          "length" => gene_data["length"],
                                          "length_offset" => length_offset)
    gene_data["variance"] !== nothing && (assays["variance"] = gene_data["variance"])
    if gene_data["inf_reps"] !== nothing
        for (i, m) in enumerate(gene_data["inf_reps"])
            assays["inf_rep_$i"] = m
        end
    end

    clean_map  = _clean_tx2gene_dict(tx2gene, ignore_tx_version, ignore_after_bar)
    ignored_tx = [t for t in tx_ids if !haskey(clean_map, t)]
    mapped_tx  = [t for t in tx_ids if haskey(clean_map, t)]
    transcript_counts = Dict{String,Int}(sample_list[j] => length(file_tx_ids[j]) for j in 1:n_files)
    missing_tx_by_sample = Dict{String,Vector{String}}()
    for j in 1:n_files
        s = Set(file_tx_ids[j])
        missing_tx_by_sample[sample_list[j]] = [t for t in mapped_tx if !(t in s)]
    end
    inf_info = gene_data["inf_reps"] !== nothing ? Dict("n_inf_reps" => length(gene_data["inf_reps"])) :
               gene_data["variance"] !== nothing ? Dict("variance_reduced" => true) : nothing

    diag = TxImportDiagnostics(file_checksums, transcript_counts, ignored_tx,
                               missing_tx_by_sample, length_offset, inf_info,
                               provenance_record("TxImportDiagnostics", "TxImport/tximport";
                                   parameters=(sample_count=n_files, gene_count=length(gene_list),
                                               ignored_tx_count=length(ignored_tx))))
    se = SummarizedExperiment(assays, Dict{Symbol,Vector}(:gene_id => gene_list),
                              Dict{Symbol,Vector}(:sample_id => sample_list),
                              Dict{Symbol,Any}(:tximport_diagnostics => diag,
                                               :provenance => diag.provenance.id))
    result = TxImportResult(count_matrix, nothing, gene_data["abundance"], gene_data["length"],
                            counts_from_abundance, se, tx2gene, ignored_tx,
                            gene_data["inf_reps"], gene_data["variance"], diag,
                            provenance_record("TxImportResult", "TxImport/tximport";
                                parameters=(type=string(type),
                                            counts_from_abundance=string(counts_from_abundance),
                                            ignore_tx_version=ignore_tx_version,
                                            ignore_after_bar=ignore_after_bar,
                                            var_reduce=var_reduce)))
    return provenance_result!(ctx, result, "tximport"; parents=String[];
                              parameters=(type=string(type),
                                          counts_from_abundance=string(counts_from_abundance)))
end

end # module