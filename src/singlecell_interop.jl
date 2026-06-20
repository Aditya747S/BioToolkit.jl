using DataFrames
using Distributions: Hypergeometric, ccdf
using HDF5
using PooledArrays
using SHA

using ..BioToolkit: AbstractAnalysisResult, ProvenanceContext, ProvenanceParams, ResultProvenance, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_record, provenance_result!, register_provenance!

export read_h5ad, write_h5ad, CellTypeAnnotationResult, annotate_cell_types, AmbientRNARemovalResult, remove_ambient_rna, SingleCellViewer, interactive_singlecell_viewer, lasso_select_cells, recluster_singlecell_viewer!, cell_hover_text, top_expressed_genes

const _H5AD_ENCODING_TYPE = "anndata"
const _H5AD_ENCODING_VERSION = "0.1.0"
const _H5AD_COLUMN_ORDER_SEPARATOR = '\u001f'

const _H5AD_MARKER_DATABASE = Dict{String,Dict{String,Vector{String}}}(
    "HumanPrimaryCellAtlas" => Dict(
        "B cell" => ["MS4A1", "CD79A", "CD74", "CD37", "HLA-DRA"],
        "Plasma cell" => ["MZB1", "XBP1", "JCHAIN", "IGHG1", "SDC1"],
        "T cell" => ["CD3D", "CD3E", "TRAC", "IL7R", "LTB"],
        "NK cell" => ["NKG7", "GNLY", "PRF1", "KLRD1", "FCGR3A"],
        "Monocyte" => ["LYZ", "S100A8", "S100A9", "FCN1", "CTSS"],
        "Dendritic cell" => ["FCER1A", "CLEC10A", "LILRA4", "GZMB", "IRF7"],
        "Macrophage" => ["C1QA", "C1QB", "C1QC", "APOE", "LST1"],
        "Endothelial cell" => ["PECAM1", "VWF", "KDR", "CLDN5", "EMCN"],
        "Fibroblast" => ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"],
        "Epithelial cell" => ["EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"],
        "Cycling cell" => ["MKI67", "TOP2A", "CDK1", "UBE2C", "BIRC5"]),
    "PanglaoDB" => Dict(
        "B cell" => ["MS4A1", "CD79A", "CD74", "CD37"],
        "T cell" => ["CD3D", "CD3E", "TRAC", "IL7R"],
        "Myeloid" => ["LYZ", "S100A8", "S100A9", "FCN1"],
        "Endothelial" => ["PECAM1", "VWF", "KDR", "CLDN5"],
        "Fibroblast" => ["COL1A1", "COL1A2", "DCN", "LUM"],
        "Epithelial" => ["EPCAM", "KRT8", "KRT18", "KRT19"]))

mutable struct CellTypeAnnotationResult <: AbstractAnalysisResult
    group_ids::Vector{String}
    reference_name::String
    method::Symbol
    reference_labels::Vector{String}
    scores::Matrix{Float64}
    assigned_labels::Vector{String}
    cell_labels::Vector{String}
    provenance::ResultProvenance
end

CellTypeAnnotationResult(group_ids, reference_name, method, reference_labels, scores, assigned_labels, cell_labels) =
    CellTypeAnnotationResult(group_ids, reference_name, method, reference_labels, scores, assigned_labels, cell_labels,
        provenance_record("CellTypeAnnotationResult", "SingleCell/annotate_cell_types";
            parameters=(method=method, reference_name=String(reference_name), group_count=length(group_ids), cell_count=length(cell_labels))))

mutable struct AmbientRNARemovalResult <: AbstractAnalysisResult
    corrected_experiment::SingleCellExperiment
    corrected_counts::SparseMatrixCSC{Int,Int}
    ambient_profile::Vector{Float64}
    contamination_fraction::Vector{Float64}
    background_cells::Vector{String}
    provenance::ResultProvenance
end

AmbientRNARemovalResult(corrected_experiment, corrected_counts, ambient_profile, contamination_fraction, background_cells) =
    AmbientRNARemovalResult(corrected_experiment, corrected_counts, ambient_profile, contamination_fraction, background_cells,
        provenance_record("AmbientRNARemovalResult", "SingleCell/remove_ambient_rna";
            parameters=(gene_count=length(corrected_experiment.gene_ids), cell_count=length(corrected_experiment.cell_ids), background_count=length(background_cells))))

mutable struct SingleCellViewer
    experiment::SingleCellExperiment
    embedding::Matrix{Float64}
    cluster_labels::Vector{Int}
    normalized_counts::Matrix{Float64}
    selected_indices::Vector{Int}
    title::String
    k::Int
    metadata_fields::Vector{String}
end

function _h5ad_jsonify(value)
    if value isa AbstractDict
        return Dict(string(key) => _h5ad_jsonify(item) for (key, item) in value)
    elseif value isa AbstractArray && !(value isa String)
        return [_h5ad_jsonify(item) for item in value]
    elseif value isa ResultProvenance
        return Dict(
            "label" => value.label,
            "source" => value.source,
            "status" => string(value.status),
            "warnings" => copy(value.warnings),
            "errors" => copy(value.errors),
            "fallbacks" => copy(value.fallbacks),
            "notes" => copy(value.notes),
            "parameters" => Dict(string(key) => item for (key, item) in pairs(value.parameters)))
    elseif value === nothing || value isa Number || value isa Bool || value isa String
        return value
    else
        return string(value)
    end
end

function _h5ad_column_to_storage(column)
    if column isa AbstractVector{<:Union{Missing,Number,Bool}}
        if any(ismissing, column)
            return [ismissing(value) ? NaN : Float64(value) for value in column]
        end
        return collect(column)
    elseif column isa AbstractVector{<:String}
        return String.(column)
    elseif column isa PooledArray
        return String.(collect(column))
    elseif any(ismissing, column)
        return [ismissing(value) ? "" : String(value) for value in column]
    else
        return [String(value) for value in column]
    end
end

function _h5ad_table_from_value(value, fallback_index::AbstractVector{<:String})
    if value isa DataFrame
        return copy(value)
    elseif value isa AbstractDict
        return DataFrame(value)
    elseif value === nothing
        return DataFrame(_index=collect(fallback_index))
    else
        return DataFrame(_index=collect(fallback_index))
    end
end

_h5ad_encode_column_order(column_names::AbstractVector{<:String}) = join(String.(column_names), string(_H5AD_COLUMN_ORDER_SEPARATOR))

function _h5ad_decode_column_order(value)
    payload = value isa HDF5.Attribute ? read(value) : value
    text = strip(String(payload))
    isempty(text) && return String[]
    return split(text, _H5AD_COLUMN_ORDER_SEPARATOR)
end

function _h5ad_write_table(parent, name::String, table::DataFrame, fallback_index::AbstractVector{<:String})
    group = haskey(parent, name) ? parent[name] : create_group(parent, name)
    HDF5.attributes(group)["encoding-type"] = "dataframe"
    HDF5.attributes(group)["encoding-version"] = "0.2.0"
    column_names = [String(column_name) for column_name in names(table) if String(column_name) != "_index"]
    HDF5.attributes(group)["column-order"] = _h5ad_encode_column_order(column_names)
    HDF5.attributes(group)["index"] = "_index"
    index = "_index" in names(table) ? table[!, :_index] : collect(fallback_index)
    write(group, "_index", String.(index))
    for column_name in names(table)
        column_name == "_index" && continue
        write(group, String(column_name), _h5ad_column_to_storage(table[!, column_name]))
    end
    return group
end

function _h5ad_read_table(group)
    table = DataFrame()
    attrs = HDF5.attributes(group)
    column_order = haskey(attrs, "column-order") ? _h5ad_decode_column_order(attrs["column-order"]) : String[]
    read_order = isempty(column_order) ? collect(keys(group)) : vcat("_index", [name for name in column_order if name != "_index"])
    for key in read_order
        haskey(group, key) || continue
        table[!, Symbol(key)] = collect(read(group, key))
    end
    return table
end

function _h5ad_write_sparse(parent, name::String, matrix::SparseMatrixCSC)
    group = haskey(parent, name) ? parent[name] : create_group(parent, name)
    HDF5.attributes(group)["encoding-type"] = "csr_matrix"
    HDF5.attributes(group)["encoding-version"] = "0.1.0"
    transposed = sparse(transpose(matrix))
    write(group, "data", transposed.nzval)
    write(group, "indices", Int32.(transposed.rowval .- 1))
    write(group, "indptr", Int32.(transposed.colptr .- 1))
    write(group, "shape", Int64[size(matrix, 1), size(matrix, 2)])
    return group
end

function _h5ad_write_matrix(parent, name::String, matrix::AbstractMatrix{<:Real})
    if matrix isa SparseMatrixCSC
        return _h5ad_write_sparse(parent, name, matrix)
    end
    write(parent, name, Matrix{Float64}(matrix))
    return parent[name]
end

function _h5ad_read_sparse(group)
    data = Float64.(collect(read(group, "data")))
    indices = Int.(collect(read(group, "indices"))) .+ 1
    indptr = Int.(collect(read(group, "indptr"))) .+ 1
    shape = Tuple(Int.(collect(read(group, "shape"))))
    nrows, ncols = shape
    row_indices = Int[]
    col_indices = Int[]
    position = 1
    for row in 1:nrows
        stop = indptr[row + 1] - 1
        run = max(0, stop - indptr[row] + 1)
        if run > 0
            append!(row_indices, fill(row, run))
            append!(col_indices, indices[position:position + run - 1])
            position += run
        end
    end
    return sparse(row_indices, col_indices, data, nrows, ncols)
end

function _h5ad_read_matrix(node)
    if node isa HDF5.Group && haskey(node, "data") && haskey(node, "indices") && haskey(node, "indptr")
        return _h5ad_read_sparse(node)
    end
    return read(node)
end

function _h5ad_metadata_dict(experiment::SingleCellExperiment)
    metadata = Dict{String,Any}(experiment.metadata)
    if !haskey(metadata, "obs")
        metadata["obs"] = DataFrame(_index=copy(experiment.cell_ids))
    end
    if !haskey(metadata, "var")
        metadata["var"] = DataFrame(_index=copy(experiment.gene_ids))
    end
    if !haskey(metadata, "layers")
        metadata["layers"] = Dict{String,Any}()
    end
    metadata["obsm"] = Dict{String,Any}(name => copy(matrix) for (name, matrix) in experiment.reductions)
    return metadata
end

function write_h5ad(experiment::SingleCellExperiment, filepath::String)
    metadata = _singlecell_metadata_with_provenance(experiment.metadata; source="write_h5ad", notes=["written to $(filepath)"], parameters=(filepath=filepath, n_genes=length(experiment.gene_ids), n_cells=length(experiment.cell_ids)))
    h5open(filepath, "w") do file
        HDF5.attributes(file)["encoding-type"] = _H5AD_ENCODING_TYPE
        HDF5.attributes(file)["encoding-version"] = _H5AD_ENCODING_VERSION

        _h5ad_write_matrix(file, "X", experiment.counts)

        obs = _h5ad_table_from_value(get(metadata, "obs", nothing), experiment.cell_ids)
        var = _h5ad_table_from_value(get(metadata, "var", nothing), experiment.gene_ids)
        _h5ad_write_table(file, "obs", obs, experiment.cell_ids)
        _h5ad_write_table(file, "var", var, experiment.gene_ids)

        obsm = create_group(file, "obsm")
        HDF5.attributes(obsm)["encoding-type"] = "dict"
        HDF5.attributes(obsm)["encoding-version"] = "0.1.0"
        for (name, matrix) in experiment.reductions
            write(obsm, startswith(name, "X_") ? name : string("X_", name), Matrix{Float64}(matrix))
        end
        if experiment.spatial_coords !== nothing
            write(obsm, "spatial", Matrix{Float64}(experiment.spatial_coords))
        end

        layers = get(metadata, "layers", Dict{String,Any}())
        if layers isa AbstractDict && !isempty(layers)
            layer_group = create_group(file, "layers")
            HDF5.attributes(layer_group)["encoding-type"] = "dict"
            HDF5.attributes(layer_group)["encoding-version"] = "0.1.0"
            for (name, layer) in layers
                if layer isa AbstractMatrix{<:Real}
                    _h5ad_write_matrix(layer_group, String(name), layer)
                else
                    write(layer_group, String(name), layer)
                end
            end
        end

        uns = create_group(file, "uns")
        HDF5.attributes(uns)["encoding-type"] = "dict"
        HDF5.attributes(uns)["encoding-version"] = "0.1.0"
        sanitized = Dict{String,Any}(k => v for (k, v) in metadata if !(k in ("obs", "var", "layers", "obsm")))
        write(uns, "metadata_json", JSON.json(_h5ad_jsonify(sanitized)))
    end
    _ctx = active_provenance_context()
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(filepath)))
        register_provenance!(_ctx, "write_h5ad"; parents=provenance_parent_ids(experiment), parameters=(filepath=filepath, n_genes=length(experiment.gene_ids), n_cells=length(experiment.cell_ids), hash=provenance_hash))
    end
    return filepath
end

function _h5ad_read_uns(group)
    if group === nothing || !haskey(group, "metadata_json")
        return Dict{String,Any}()
    end
    text = String(read(group, "metadata_json"))
    isempty(strip(text)) && return Dict{String,Any}()
    parsed = JSON.parse(text)
    return parsed isa Dict ? Dict{String,Any}(string(k) => v for (k, v) in parsed) : Dict{String,Any}()
end

function _read_h5ad_file(filepath::String)
    h5open(filepath, "r") do file
        counts = _h5ad_read_matrix(file["X"])
        obs = haskey(file, "obs") ? _h5ad_read_table(file["obs"]) : DataFrame()
        var = haskey(file, "var") ? _h5ad_read_table(file["var"]) : DataFrame()
        cell_ids = "_index" in names(obs) ? String.(obs[!, :_index]) : ["cell$(index)" for index in 1:size(counts, 2)]
        gene_ids = "_index" in names(var) ? String.(var[!, :_index]) : ["gene$(index)" for index in 1:size(counts, 1)]

        metadata = _h5ad_read_uns(haskey(file, "uns") ? file["uns"] : nothing)
        metadata["obs"] = obs
        metadata["var"] = var
        update_provenance!(metadata; source="read_h5ad", notes=["loaded from $(filepath)"], parameters=(filepath=filepath, n_genes=length(gene_ids), n_cells=length(cell_ids)))

        reductions = Dict{String,Matrix{Float64}}()
        spatial_coords = nothing
        if haskey(file, "obsm")
            obsm = file["obsm"]
            for key in keys(obsm)
                matrix = Matrix{Float64}(read(obsm, key))
                if key in ("spatial", "X_spatial", "X_spatial_coords") && size(matrix, 2) >= 2
                    spatial_coords = matrix[:, 1:2]
                elseif startswith(key, "X_")
                    reductions[lowercase(replace(key, "X_" => ""))] = matrix
                else
                    reductions[lowercase(String(key))] = matrix
                end
            end
        end

        if haskey(file, "layers")
            layers = Dict{String,Any}()
            for key in keys(file["layers"])
                node = file["layers"][key]
                layers[String(key)] = node isa HDF5.Group ? _h5ad_read_matrix(node) : read(node)
            end
            metadata["layers"] = layers
        end

        counts_matrix = counts isa SparseMatrixCSC ? sparse(Int.(round.(counts))) : sparse(Int.(round.(Matrix{Float64}(counts))))
        experiment = SingleCellExperiment(counts_matrix, gene_ids, cell_ids; metadata=metadata, spatial_coords=spatial_coords)
        experiment.reductions = reductions
        return experiment
    end
end

function read_h5ad(filepath::String)
    experiment = _read_h5ad_file(filepath)
    _ctx = active_provenance_context()
    if _ctx !== nothing
        provenance_hash = bytes2hex(sha256(read(filepath)))
        register_provenance!(_ctx, "read_h5ad"; parents=String[], parameters=(source=filepath, filepath=filepath, hash=provenance_hash, n_genes=length(experiment.gene_ids), n_cells=length(experiment.cell_ids)))
    end
    return experiment
end

function _group_ids(labels::AbstractVector{<:Integer})
    unique_labels = sort!(unique(Int.(labels)))
    return unique_labels, string.(unique_labels)
end

function _group_labels(experiment::SingleCellExperiment, labels)
    if labels !== nothing
        length(labels) == length(experiment.cell_ids) || throw(ArgumentError("labels must match the number of cells"))
        return Int.(labels), "provided"
    elseif !isempty(experiment.clusters)
        key = first(sort(collect(keys(experiment.clusters))))
        group_labels = experiment.clusters[key]
        length(group_labels) == length(experiment.cell_ids) || throw(ArgumentError("stored cluster labels do not match the number of cells"))
        return Int.(group_labels), key
    else
        return collect(1:length(experiment.cell_ids)), "cell"
    end
end

function _group_expression(matrix::AbstractMatrix{<:Real}, labels::AbstractVector{<:Integer})
    group_ids, group_names = _group_ids(labels)
    profiles = zeros(Float64, size(matrix, 1), length(group_ids))
    for (group_position, group_id) in enumerate(group_ids)
        cell_indices = findall(==(group_id), labels)
        profiles[:, group_position] .= isempty(cell_indices) ? zeros(Float64, size(matrix, 1)) : vec(mean(matrix[:, cell_indices], dims=2))
    end
    return group_ids, group_names, profiles
end

function _marker_overlap_score(query_markers::Set{String}, reference_markers::Vector{String})
    reference_set = Set{String}(reference_markers)
    intersection = length(intersect(query_markers, reference_set))
    union_size = length(union(query_markers, reference_set))
    jaccard = union_size == 0 ? 0.0 : intersection / union_size
    universe = max(length(query_markers) + length(reference_set), 1)
    fisher = if intersection == 0
        1.0
    else
        ccdf(Hypergeometric(length(reference_set), max(universe - length(reference_set), 1), length(query_markers)), intersection - 1)
    end
    return jaccard + (-log10(max(fisher, eps(Float64))) / 10)
end

function _top_markers_for_profile(profile::AbstractVector{<:Real}, gene_ids::Vector{String}; top_n::Int=20, min_expression::Real=0.0)
    ranking = sortperm(profile; rev=true)
    genes = String[]
    for index in ranking
        profile[index] <= min_expression && break
        push!(genes, gene_ids[index])
        length(genes) >= top_n && break
    end
    return genes
end

function _annotation_database(reference::String)
    haskey(_H5AD_MARKER_DATABASE, reference) || throw(ArgumentError("unknown cell type reference database '$reference'"))
    return _H5AD_MARKER_DATABASE[reference]
end

function _reference_profiles(reference::SingleCellExperiment; reference_labels=nothing)
    labels, _ = _group_labels(reference, reference_labels)
    matrix = Matrix{Float64}(normalize_counts(reference))
    group_ids, group_names, profiles = _group_expression(matrix, labels)
    label_names = if reference_labels isa AbstractVector{<:String}
        String.(reference_labels)
    elseif haskey(reference.metadata, "cell_type_labels")
        String.(reference.metadata["cell_type_labels"])
    elseif haskey(reference.metadata, "labels")
        String.(reference.metadata["labels"])
    elseif !isempty(reference.clusters)
        key = first(sort(collect(keys(reference.clusters))))
        string.(reference.clusters[key])
    else
        ["reference_$(index)" for index in 1:length(group_ids)]
    end
    return group_ids, label_names, profiles
end

function _score_against_reference(profile::AbstractVector{<:Real}, reference_profiles::AbstractMatrix{<:Real})
    candidate_scores = zeros(Float64, size(reference_profiles, 2))
    normalized_profile = profile ./ max(norm(profile), eps(Float64))
    for index in 1:size(reference_profiles, 2)
        candidate = view(reference_profiles, :, index)
        candidate_norm = max(norm(candidate), eps(Float64))
        candidate_scores[index] = dot(normalized_profile, candidate ./ candidate_norm)
    end
    return candidate_scores
end

function annotate_cell_types(experiment::SingleCellExperiment; labels=nothing, reference="HumanPrimaryCellAtlas", reference_labels=nothing, method::Symbol=:marker, top_n::Int=20, min_expression::Real=0.0, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    method in (:marker, :reference) || throw(ArgumentError("method must be :marker or :reference"))
    group_labels, grouping_name = _group_labels(experiment, labels)
    expression = Matrix{Float64}(normalize_counts(experiment))
    group_ids, group_names, profiles = _group_expression(expression, group_labels)

    if method == :marker && reference isa String
        database = _annotation_database(reference)
        reference_names = sort!(collect(keys(database)))
        scores = zeros(Float64, length(group_names), length(reference_names))
        assigned = Vector{String}(undef, length(group_names))
        for (group_index, group_name) in enumerate(group_names)
            markers = Set(_top_markers_for_profile(profiles[:, group_index], experiment.gene_ids; top_n=top_n, min_expression=min_expression))
            for (reference_index, cell_type) in enumerate(reference_names)
                scores[group_index, reference_index] = _marker_overlap_score(markers, database[cell_type])
            end
            assigned[group_index] = reference_names[argmax(scores[group_index, :])]
        end
        cell_labels = [assigned[findfirst(==(group_labels[cell_index]), group_ids)] for cell_index in eachindex(group_labels)]
        result = CellTypeAnnotationResult(group_names, String(reference), method, reference_names, scores, assigned, cell_labels)
        experiment.metadata["cell_type_annotations"] = result
        experiment.metadata["cell_type_labels"] = cell_labels
        return provenance_result!(_ctx, result, "annotate_cell_types";
            parents=provenance_parent_ids(experiment, reference),
            parameters=(method=method, reference=String(reference), group_count=length(group_names), cell_count=length(cell_labels)))
    elseif method == :marker && reference isa AbstractDict
        reference_names = sort!(collect(keys(reference)))
        scores = zeros(Float64, length(group_names), length(reference_names))
        assigned = Vector{String}(undef, length(group_names))
        for (group_index, group_name) in enumerate(group_names)
            markers = Set(_top_markers_for_profile(profiles[:, group_index], experiment.gene_ids; top_n=top_n, min_expression=min_expression))
            for (reference_index, cell_type) in enumerate(reference_names)
                scores[group_index, reference_index] = _marker_overlap_score(markers, String.(reference[cell_type]))
            end
            assigned[group_index] = reference_names[argmax(scores[group_index, :])]
        end
        cell_labels = [assigned[findfirst(==(group_labels[cell_index]), group_ids)] for cell_index in eachindex(group_labels)]
        result = CellTypeAnnotationResult(group_names, "custom", method, reference_names, scores, assigned, cell_labels)
        experiment.metadata["cell_type_annotations"] = result
        experiment.metadata["cell_type_labels"] = cell_labels
        return provenance_result!(_ctx, result, "annotate_cell_types";
            parents=provenance_parent_ids(experiment),
            parameters=(method=method, reference="custom", group_count=length(group_names), cell_count=length(cell_labels)))
    else
        if !(reference isa SingleCellExperiment)
            throw(ArgumentError("reference mapping currently expects a SingleCellExperiment reference or a marker database"))
        end
        ref_group_ids, reference_names, reference_profiles = _reference_profiles(reference; reference_labels=reference_labels)
        scores = zeros(Float64, length(group_names), length(reference_names))
        assigned = Vector{String}(undef, length(group_names))
        for group_index in 1:length(group_names)
            scores[group_index, :] .= _score_against_reference(profiles[:, group_index], reference_profiles)
            assigned[group_index] = reference_names[argmax(scores[group_index, :])]
        end
        cell_labels = [assigned[findfirst(==(group_labels[cell_index]), group_ids)] for cell_index in eachindex(group_labels)]
        result = CellTypeAnnotationResult(group_names, "reference_mapping", method, reference_names, scores, assigned, cell_labels)
        experiment.metadata["cell_type_annotations"] = result
        experiment.metadata["cell_type_labels"] = cell_labels
        return provenance_result!(_ctx, result, "annotate_cell_types";
            parents=provenance_parent_ids(experiment, reference),
            parameters=(method=method, reference="reference_mapping", group_count=length(group_names), cell_count=length(cell_labels)))
    end
end

function _empty_droplet_indices(library_sizes::AbstractVector{<:Real}, ambient_quantile::Real)
    n = length(library_sizes)
    n == 0 && return Int[]
    sorted = sort(Float64.(library_sizes))
    position = max(1, min(n, round(Int, ambient_quantile * n)))
    threshold = sorted[position]
    return findall(size -> size <= threshold, library_sizes)
end

function _ambient_profile(counts::AbstractMatrix{<:Real}, indices::AbstractVector{<:Integer})
    if isempty(indices)
        return fill(1.0 / size(counts, 1), size(counts, 1))
    end
    profile = vec(sum(Matrix{Float64}(counts[:, indices]), dims=2))
    total = sum(profile)
    return total > 0 ? profile ./ total : fill(1.0 / size(counts, 1), size(counts, 1))
end

function remove_ambient_rna(experiment::SingleCellExperiment; empty_droplets=nothing, ambient_quantile::Real=0.05, max_contamination::Real=0.25, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    counts = Matrix{Float64}(experiment.counts)
    library_sizes = vec(sum(counts, dims=1))
    background_indices = if empty_droplets === nothing
        _empty_droplet_indices(library_sizes, ambient_quantile)
    elseif empty_droplets isa AbstractVector{Bool}
        findall(identity, empty_droplets)
    else
        Int.(empty_droplets)
    end
    ambient = _ambient_profile(counts, background_indices)
    similarity = [dot(counts[:, cell_index] ./ max(library_sizes[cell_index], eps(Float64)), ambient) for cell_index in 1:size(counts, 2)]
    similarity_min = minimum(similarity)
    similarity_range = maximum(similarity) - similarity_min
    contamination = if similarity_range > 0
        clamp.((similarity .- similarity_min) ./ similarity_range .* max_contamination, 0.0, max_contamination)
    else
        fill(min(max_contamination, ambient_quantile * max_contamination), length(similarity))
    end

    corrected = similar(counts)
    for cell_index in 1:size(counts, 2)
        contamination_counts = contamination[cell_index] * library_sizes[cell_index] .* ambient
        corrected[:, cell_index] .= max.(counts[:, cell_index] .- contamination_counts, 0.0)
    end

    corrected_sparse = sparse(Int.(round.(corrected)))
    corrected_experiment = SingleCellExperiment(corrected_sparse, experiment.gene_ids, experiment.cell_ids; metadata=copy(experiment.metadata), spatial_coords=experiment.spatial_coords)
    corrected_experiment.reductions = Dict{String,Matrix{Float64}}(name => copy(matrix) for (name, matrix) in experiment.reductions)
    corrected_experiment.clusters = Dict{String,Vector{Int}}(name => copy(labels) for (name, labels) in experiment.clusters)
    corrected_experiment.variable_features = Dict{String,Vector{String}}(name => copy(features) for (name, features) in experiment.variable_features)
    corrected_experiment.neighbors = Dict{String,Vector{Vector{Int}}}(name => [copy(neighbor_list) for neighbor_list in neighbors] for (name, neighbors) in experiment.neighbors)
    corrected_experiment.metadata["ambient_rna"] = Dict(
        "ambient_profile" => ambient,
        "contamination_fraction" => contamination,
        "background_cells" => experiment.cell_ids[background_indices])

    result = AmbientRNARemovalResult(corrected_experiment, corrected_sparse, ambient, contamination, experiment.cell_ids[background_indices])
    return provenance_result!(_ctx, result, "remove_ambient_rna";
        parents=provenance_parent_ids(experiment),
        parameters=(ambient_quantile=ambient_quantile, max_contamination=max_contamination, background_count=length(background_indices)))
end

function _default_viewer_embedding(experiment::SingleCellExperiment; embedding::Union{Nothing,AbstractMatrix}=nothing, reduction::String="umap", k::Int=15)
    if embedding !== nothing
        return Matrix{Float64}(embedding)
    elseif haskey(experiment.reductions, reduction)
        return Matrix{Float64}(experiment.reductions[reduction])
    elseif reduction == "pca"
        return run_pca(experiment; n_components=max(2, k))
    else
        return run_umap(experiment; n_neighbors=max(k, 5))
    end
end

function _default_viewer_labels(experiment::SingleCellExperiment, embedding::AbstractMatrix, k::Int)
    if !isempty(experiment.clusters)
        key = first(sort(collect(keys(experiment.clusters))))
        labels = experiment.clusters[key]
        length(labels) == size(embedding, 1) && return Int.(labels)
    end
    if size(embedding, 1) == length(experiment.cell_ids)
        return find_clusters(experiment; embedding=embedding, k=k, method=:leiden)
    end
    return collect(1:size(embedding, 1))
end

function interactive_singlecell_viewer(experiment::SingleCellExperiment; embedding::Union{Nothing,AbstractMatrix}=nothing, reduction::String="umap", k::Int=15, labels=nothing, title::String="Single-cell viewer", metadata_fields::AbstractVector{<:String}=String[])
    coords = _default_viewer_embedding(experiment; embedding=embedding, reduction=reduction, k=k)
    cluster_labels = labels === nothing ? _default_viewer_labels(experiment, coords, k) : Int.(labels)
    length(cluster_labels) == size(coords, 1) || throw(DimensionMismatch("labels must match the number of embedded cells"))
    normalized = normalize_counts(experiment)
    viewer = SingleCellViewer(experiment, coords, cluster_labels, Matrix{Float64}(normalized), Int[], String(title), k, String.(metadata_fields))
    _ctx = active_provenance_context()
    _ctx !== nothing && register_provenance!(_ctx, "interactive_singlecell_viewer"; parents=provenance_parent_ids(experiment), parameters=(reduction=reduction, k=k, label_count=length(cluster_labels), title=String(title)))
    return viewer
end

function _point_in_polygon(point::NTuple{2,Float64}, polygon::AbstractVector{<:NTuple{2,Float64}})
    inside = false
    n = length(polygon)
    n < 3 && return false
    x, y = point
    j = n
    for i in 1:n
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / max(yj - yi, eps(Float64)) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

function lasso_select_cells(viewer::SingleCellViewer, polygon::AbstractVector)
    polygon_points = NTuple{2,Float64}[]
    for point in polygon
        if point isa NTuple{2,<:Real}
            push!(polygon_points, (Float64(point[1]), Float64(point[2])))
        elseif point isa AbstractVector && length(point) >= 2
            push!(polygon_points, (Float64(point[1]), Float64(point[2])))
        else
            throw(ArgumentError("polygon points must be 2D coordinates"))
        end
    end
    selected = Int[]
    for cell_index in 1:size(viewer.embedding, 1)
        candidate = (viewer.embedding[cell_index, 1], viewer.embedding[cell_index, 2])
        _point_in_polygon(candidate, polygon_points) && push!(selected, cell_index)
    end
    viewer.selected_indices = selected
    _ctx = active_provenance_context()
    _ctx !== nothing && register_provenance!(_ctx, "lasso_select_cells"; parents=provenance_parent_ids(viewer), parameters=(selected_count=length(selected), polygon_points=length(polygon_points)))
    return selected
end

function recluster_singlecell_viewer!(viewer::SingleCellViewer; k::Int=viewer.k)
    viewer.k = k
    viewer.cluster_labels = find_clusters(viewer.experiment; embedding=viewer.embedding, k=k, method=:leiden)
    viewer.experiment.clusters["interactive_viewer"] = copy(viewer.cluster_labels)
    _ctx = active_provenance_context()
    _ctx !== nothing && register_provenance!(_ctx, "recluster_singlecell_viewer!"; parents=provenance_parent_ids(viewer), parameters=(k=k, cluster_count=length(unique(viewer.cluster_labels))))
    return viewer.cluster_labels
end

function top_expressed_genes(experiment::SingleCellExperiment, cell_index::Integer; top_n::Int=5, normalize::Bool=true)
    1 <= cell_index <= length(experiment.cell_ids) || throw(ArgumentError("cell_index out of bounds"))
    expression = normalize ? Matrix{Float64}(normalize_counts(experiment)) : Matrix{Float64}(experiment.counts)
    ranked = sortperm(viewer_column(expression, cell_index); rev=true)
    genes = Pair{String,Float64}[]
    for gene_index in ranked[1:min(top_n, length(ranked))]
        push!(genes, experiment.gene_ids[gene_index] => Float64(expression[gene_index, cell_index]))
    end
    _ctx = active_provenance_context()
    _ctx !== nothing && register_provenance!(_ctx, "top_expressed_genes"; parents=provenance_parent_ids(experiment), parameters=(cell_index=cell_index, top_n=top_n, normalize=normalize, gene_count=length(genes)))
    return genes
end

function viewer_column(matrix::AbstractMatrix{<:Real}, column::Integer)
    return view(matrix, :, column)
end

function cell_hover_text(viewer::SingleCellViewer, index; top_n::Int=5)
    cell_index = index isa Tuple ? Int(first(index)) : Int(index)
    1 <= cell_index <= length(viewer.experiment.cell_ids) || return "cell out of bounds"
    coords = viewer.embedding[cell_index, 1:min(2, size(viewer.embedding, 2))]
    lines = String[]
    push!(lines, "cell: $(viewer.experiment.cell_ids[cell_index])")
    push!(lines, "cluster: $(viewer.cluster_labels[cell_index])")
    push!(lines, "embedding: $(round.(coords; digits=3))")
    obs = get(viewer.experiment.metadata, "obs", nothing)
    if obs isa DataFrame && nrow(obs) >= cell_index
        row = obs[cell_index, :]
        for field in viewer.metadata_fields
            symbol = Symbol(field)
            hasproperty(row, symbol) && push!(lines, string(field, ": ", getproperty(row, symbol)))
        end
    end
    push!(lines, "top genes: " * join([string(gene, "=", round(value; digits=3)) for (gene, value) in top_expressed_genes(viewer.experiment, cell_index; top_n=top_n)], ", "))
    text = join(lines, "\n")
    _ctx = active_provenance_context()
    _ctx !== nothing && register_provenance!(_ctx, "cell_hover_text"; parents=provenance_parent_ids(viewer), parameters=(index=cell_index, top_n=top_n))
    return text
end

function cell_hover_text(experiment::SingleCellExperiment, index; top_n::Int=5)
    viewer = interactive_singlecell_viewer(experiment)
    return cell_hover_text(viewer, index; top_n=top_n)
end
