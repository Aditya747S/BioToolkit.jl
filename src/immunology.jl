# ==============================================================================
# immunology.jl — Adaptive immune repertoire and epitope analysis
#
# References:
#   - Bolotin et al. (2015) Nat Methods 12:380-381 (MiXCR concepts)
#   - Reynolds et al. (2020) NAR 48:W449-W454 (B-cell epitope context)
#   - Glanville et al. (2017) Nature 547:94-98 (GLIPH)
#   - Vander Heiden et al. (2020) Immunity 53:323-338 (Alakazam)
#   - Kovaltsuk et al. (2018) J Immunol 201:2502-2509 (OAS / SHM)
# ==============================================================================

module Immunology

using DataFrames
using Distributions
using Statistics
using LinearAlgebra
using ..BioToolkit: AminoAcidAlphabet, BioSequence, DNAAlphabet, threaded_foreach
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!, with_provenance

@inline function _register_immunology_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

export extract_cdr3, clonotype_table, assign_vdj_segments, isotype_switching_summary
export germline_usage_bias, bepipred_like_scores, mhcflurry_like_scores
export repertoire_diversity_metrics, clonal_expansion_test, cdr3_length_spectrum
export tcrdist_like_matrix
export igblast_assign, imgt_highvquest_manifest, paired_clonotype_table, differential_vj_usage

# New exports
export somatic_hypermutation_rate
export bcr_affinity_maturation_score
export antigen_specificity_clustering
export clonotype_trajectory
export v_gene_family_usage
export hla_type_enrichment
export cdr3_physicochemical_properties
export network_centrality_clonotypes
export convergent_clonotype_detection
export lineage_tree_from_clones
export epitope_binding_profile

function extract_cdr3(sequence::AbstractString)
    seq = uppercase(String(sequence))
    m = match(r"C[A-Z]{5,28}(FG|WG|QYF|HF|TF)", seq)

    return m === nothing ? missing : m.match
end

function clonotype_table(sequences::AbstractVector{<:AbstractString}; sample_ids=nothing)
    cdr3 = [extract_cdr3(seq) for seq in sequences]
    key = [ismissing(x) ? "unresolved" : String(x) for x in cdr3]
    tab = combine(groupby(DataFrame(clonotype=key), :clonotype), nrow => :n_cells)
    tab[!, :frequency] = tab.n_cells ./ max(sum(tab.n_cells), 1)

    if sample_ids !== nothing
        length(sample_ids) == length(sequences) || throw(DimensionMismatch("sample_ids must match sequence count"))
        detail = DataFrame(clonotype=key, sample=String.(sample_ids))
        per_sample = combine(groupby(detail, [:sample, :clonotype]), nrow => :n)
        summary_tab = with_provenance(sort(tab, :n_cells, rev=true), "ImmunologyTable", "Immunology/clonotype_table"; parameters=(sequence_count=length(sequences), sample_split=true))
        per_sample_tab = with_provenance(per_sample, "ImmunologyTable", "Immunology/clonotype_table"; notes=["per-sample clonotype counts"], parameters=(sequence_count=length(sequences), sample_split=true))
        result = (summary=summary_tab, per_sample=per_sample_tab)
        _register_immunology_result!(_ctx, summary_tab, "clonotype_table"; parameters=(n_sequences=length(sequences), n_clonotypes=nrow(tab), sample_split=true))
        return result
    end
    result = with_provenance(sort(tab, :n_cells, rev=true), "ImmunologyTable", "Immunology/clonotype_table"; parameters=(sequence_count=length(sequences), sample_split=false))
    _ctx = active_provenance_context()


    return _register_immunology_result!(_ctx, result, "clonotype_table"; parameters=(n_sequences=length(sequences), n_clonotypes=nrow(tab), sample_split=false))
end

function _kmer_overlap_score(query::AbstractString, reference::AbstractString; k::Int=5)
    ncodeunits(query) < k && return 0
    qset = Set{String}(query[i:(i + k - 1)] for i in 1:(ncodeunits(query) - k + 1))
    rset = ncodeunits(reference) < k ? Set{String}() : Set{String}(reference[i:(i + k - 1)] for i in 1:(ncodeunits(reference) - k + 1))
    return length(intersect(qset, rset))
end

function assign_vdj_segments(sequences::AbstractVector{<:AbstractString}, v_db::AbstractDict, j_db::AbstractDict; k::Int=5)
    out = DataFrame(sequence=String[], v_segment=String[], j_segment=String[], v_score=Int[], j_score=Int[])
    for seq in sequences
        s = uppercase(String(seq))

        best_v = "unknown"
        best_v_score = -1
        for (name, ref) in v_db
            sc = _kmer_overlap_score(s, uppercase(String(ref)); k=k)
            if sc > best_v_score
                best_v_score = sc
                best_v = String(name)
            end
        end

        best_j = "unknown"
        best_j_score = -1
        for (name, ref) in j_db
            sc = _kmer_overlap_score(s, uppercase(String(ref)); k=k)
            if sc > best_j_score
                best_j_score = sc
                best_j = String(name)
            end
        end

        push!(out, (String(seq), best_v, best_j, best_v_score, best_j_score))
    end

    return with_provenance(out, "ImmunologyTable", "Immunology/assign_vdj_segments"; parameters=(sequence_count=length(sequences), k=Int(k)))
end

function isotype_switching_summary(clone_ids::AbstractVector, isotypes::AbstractVector)
    length(clone_ids) == length(isotypes) || throw(DimensionMismatch("clone_ids and isotypes must have same length"))
    df = DataFrame(clone=String.(clone_ids), isotype=String.(isotypes))
    grouped = groupby(df, :clone)

    out = DataFrame(clone=String[], n_cells=Int[], n_unique_isotypes=Int[], switches=Int[], dominant_isotype=String[])
    for g in grouped
        counts = combine(groupby(g, :isotype), nrow => :n)
        sort!(counts, :n, rev=true)
        nuniq = nrow(counts)
        push!(out, (g.clone[1], nrow(g), nuniq, max(nuniq - 1, 0), counts.isotype[1]))
    end

    return with_provenance(out, "ImmunologyTable", "Immunology/isotype_switching_summary"; parameters=(clone_count=length(unique(String.(clone_ids))), observation_count=length(clone_ids)))
end

function germline_usage_bias(v_genes::AbstractVector{<:AbstractString})
    counts = combine(groupby(DataFrame(v_gene=String.(v_genes)), :v_gene), nrow => :count)
    sort!(counts, :count, rev=true)
    total = sum(counts.count)
    expected = fill(total / nrow(counts), nrow(counts))
    chi2 = sum((counts.count .- expected) .^ 2 ./ max.(expected, eps(Float64)))
    p = ccdf(Chisq(max(nrow(counts) - 1, 1)), chi2)
    counts[!, :frequency] = counts.count ./ max(total, 1)
    counts[!, :chi_square] = fill(chi2, nrow(counts))
    counts[!, :pvalue] = fill(p, nrow(counts))

    return with_provenance(counts, "ImmunologyTable", "Immunology/germline_usage_bias"; parameters=(gene_count=length(v_genes), category_count=nrow(counts)))
end

function bepipred_like_scores(sequence::AbstractString; window::Int=9)
    scale = Dict(
        'A' => 0.62, 'R' => -2.53, 'N' => -0.78, 'D' => -0.90, 'C' => 0.29,
        'Q' => -0.85, 'E' => -0.74, 'G' => 0.48, 'H' => -0.40, 'I' => 1.38,
        'L' => 1.06, 'K' => -1.50, 'M' => 0.64, 'F' => 1.19, 'P' => 0.12,
        'S' => -0.18, 'T' => -0.05, 'W' => 0.81, 'Y' => 0.26, 'V' => 1.08)
    seq = collect(uppercase(String(sequence)))
    n = length(seq)
    scores = zeros(Float64, n)
    half = max(div(window, 2), 0)
    for i in 1:n
        lo = max(1, i - half)
        hi = min(n, i + half)
        vals = Float64[get(scale, seq[j], 0.0) for j in lo:hi]
        scores[i] = mean(vals)
    end
    prob = 1.0 ./ (1.0 .+ exp.(-2 .* scores))

    return with_provenance(DataFrame(position=1:n, residue=string.(seq), score=scores, epitope_probability=prob), "ImmunologyTable", "Immunology/bepipred_like_scores"; parameters=(sequence_length=n, window=Int(window)))
end

function mhcflurry_like_scores(peptides::AbstractVector{<:AbstractString}; allele::String="HLA-A*02:01")
    anchor = Dict(2 => Set(['L', 'V', 'I', 'M']), 9 => Set(['L', 'V', 'I']))
    out = DataFrame(peptide=String[], allele=String[], score=Float64[], binder=Bool[])
    for pep in peptides
        p = uppercase(String(pep))
        score = 0.0
        for (pos, allowed) in anchor
            if pos <= ncodeunits(p)
                score += p[pos] in allowed ? 1.0 : -0.5
            end
        end
        score += ncodeunits(p) == 9 ? 0.5 : -0.25 * abs(ncodeunits(p) - 9)
        binding = score >= 1.0
        push!(out, (String(pep), String(allele), score, binding))
    end
    sort!(out, :score, rev=true)

    return with_provenance(out, "ImmunologyTable", "Immunology/mhcflurry_like_scores"; parameters=(peptide_count=length(peptides), allele=allele))
end

"""
    repertoire_diversity_metrics(clonotypes)

Compute Shannon, Simpson, inverse Simpson, and observed richness.
"""
function repertoire_diversity_metrics(clonotypes::AbstractVector{<:AbstractString})
    counts = combine(groupby(DataFrame(clonotype=String.(clonotypes)), :clonotype), nrow => :count)
    n = sum(counts.count)
    n > 0 || return (n_clones=0, shannon=0.0, simpson=0.0, inverse_simpson=0.0, richness=0)

    p = counts.count ./ n
    shannon = -sum(pi > 0 ? pi * log(pi) : 0.0 for pi in p)
    simpson = sum(abs2, p)
    inv_simpson = simpson > 0 ? 1 / simpson : 0.0
    result = (n_clones=n, shannon=shannon, simpson=simpson, inverse_simpson=inv_simpson, richness=nrow(counts))
    _ctx = active_provenance_context()


    return _register_immunology_result!(_ctx, result, "repertoire_diversity_metrics"; parameters=(n_clonotypes=length(clonotypes), n_unique=nrow(counts)))
end

function _fisher_right_tail(a::Int, b::Int, c::Int, d::Int)
    (a < 0 || b < 0 || c < 0 || d < 0) && return 1.0
    N = a + b + c + d
    K = a + c
    n = a + b
    if N == 0 || K == 0 || n == 0
        return 1.0
    end
    return ccdf(Hypergeometric(N, K, n), a - 1)
end

"""
    clonal_expansion_test(clonotypes, groups; case_group=nothing)

Fisher exact enrichment test of clonotypes between case and background groups.
"""
function clonal_expansion_test(clonotypes::AbstractVector{<:AbstractString}, groups::AbstractVector{<:AbstractString}; case_group=nothing)
    length(clonotypes) == length(groups) || throw(DimensionMismatch("clonotypes and groups must have equal length"))
    g = String.(groups)
    cg = case_group === nothing ? first(sort!(unique(g))) : String(case_group)
    case_idx = findall(==(cg), g)
    ctrl_idx = setdiff(collect(eachindex(g)), case_idx)
    isempty(case_idx) && throw(ArgumentError("case_group has no observations"))
    isempty(ctrl_idx) && throw(ArgumentError("need at least one non-case observation"))

    clones = sort!(unique(String.(clonotypes)))
    out = DataFrame(clonotype=String[], case_count=Int[], control_count=Int[], odds_ratio=Float64[], pvalue=Float64[])
    for c in clones
        a = count(i -> clonotypes[i] == c, case_idx)
        b = length(case_idx) - a
        cc = count(i -> clonotypes[i] == c, ctrl_idx)
        d = length(ctrl_idx) - cc
        orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (cc + 0.5))
        p = _fisher_right_tail(a, b, cc, d)
        push!(out, (c, a, cc, orr, p))
    end
    out[!, :padj] = sortperm(out.pvalue) |> (perm -> begin
        p = out.pvalue
        q = similar(p)
        m = length(p)
        running = 1.0
        for rank in m:-1:1
            idx = perm[rank]
            running = min(running, p[idx] * m / rank)
            q[idx] = running
        end
        clamp.(q, 0.0, 1.0)
    end)
    sort!(out, [:padj, :pvalue])

    return with_provenance(out, "ImmunologyTable", "Immunology/clonal_expansion_test"; parameters=(clonotype_count=length(clonotypes), group_count=length(unique(groups))))

end

extract_cdr3(sequence::BioSequence{AminoAcidAlphabet}) = extract_cdr3(String(sequence))
clonotype_table(sequences::AbstractVector{<:BioSequence{AminoAcidAlphabet}}; sample_ids=nothing) = clonotype_table(String.(sequences); sample_ids=sample_ids)
assign_vdj_segments(sequences::AbstractVector{<:BioSequence{AminoAcidAlphabet}}, v_db::AbstractDict, j_db::AbstractDict; k::Int=5) = assign_vdj_segments(String.(sequences), v_db, j_db; k=k)
bepipred_like_scores(sequence::BioSequence{AminoAcidAlphabet}; window::Int=9) = bepipred_like_scores(String(sequence); window=window)
mhcflurry_like_scores(peptides::AbstractVector{<:BioSequence{AminoAcidAlphabet}}; allele::String="HLA-A*02:01") = mhcflurry_like_scores(String.(peptides); allele=allele)
cdr3_length_spectrum(sequences::AbstractVector{<:BioSequence{AminoAcidAlphabet}}) = cdr3_length_spectrum(String.(sequences))
cdr3_physicochemical_properties(sequences::AbstractVector{<:BioSequence{AminoAcidAlphabet}}) = cdr3_physicochemical_properties(String.(sequences))
somatic_hypermutation_rate(observed_sequences::AbstractVector{<:BioSequence{AminoAcidAlphabet}}, germline_sequences::AbstractVector{<:BioSequence{AminoAcidAlphabet}}; clone_ids=nothing) = somatic_hypermutation_rate(String.(observed_sequences), String.(germline_sequences); clone_ids=clone_ids)
bcr_affinity_maturation_score(sequences::AbstractVector{<:BioSequence{DNAAlphabet}}) = bcr_affinity_maturation_score(String.(sequences))
antigen_specificity_clustering(cdr3s::AbstractVector{<:BioSequence{AminoAcidAlphabet}}; similarity_threshold::Real=2.0, min_cluster_size::Int=2) = antigen_specificity_clustering(String.(cdr3s); similarity_threshold=similarity_threshold, min_cluster_size=min_cluster_size)
network_centrality_clonotypes(cdr3s::AbstractVector{<:BioSequence{AminoAcidAlphabet}}; similarity_threshold::Real=3.0) = network_centrality_clonotypes(String.(cdr3s); similarity_threshold=similarity_threshold)
convergent_clonotype_detection(sequences::AbstractVector{<:BioSequence{AminoAcidAlphabet}}, vgenes::AbstractVector, jgenes::AbstractVector; max_hamming::Int=1, min_group_size::Int=2) = convergent_clonotype_detection(String.(sequences), vgenes, jgenes; max_hamming=max_hamming, min_group_size=min_group_size)
tcrdist_like_matrix(cdr3_sequences::AbstractVector{<:BioSequence{AminoAcidAlphabet}}; mismatch_penalty::Real=1.0, gap_penalty::Real=1.5, threaded::Bool=true) = tcrdist_like_matrix(String.(cdr3_sequences); mismatch_penalty=mismatch_penalty, gap_penalty=gap_penalty, threaded=threaded)
paired_clonotype_table(alpha_cdr3::AbstractVector{<:BioSequence{AminoAcidAlphabet}}, beta_cdr3::AbstractVector{<:BioSequence{AminoAcidAlphabet}}; sample_ids=nothing) = paired_clonotype_table(String.(alpha_cdr3), String.(beta_cdr3); sample_ids=sample_ids)
lineage_tree_from_clones(sequences::AbstractVector{<:BioSequence{AminoAcidAlphabet}}, germline::BioSequence{AminoAcidAlphabet}; clone_ids::Union{Nothing,AbstractVector}=nothing) = lineage_tree_from_clones(String.(sequences), String(germline); clone_ids=clone_ids)
epitope_binding_profile(peptides::AbstractVector{<:BioSequence{AminoAcidAlphabet}}, known_epitopes::AbstractVector{<:BioSequence{AminoAcidAlphabet}}; allele::String="HLA-A*02:01", max_hamming::Int=1) = epitope_binding_profile(String.(peptides), String.(known_epitopes); allele=allele, max_hamming=max_hamming)

"""
    cdr3_length_spectrum(sequences)

Histogram-like table of inferred CDR3 amino-acid lengths.
"""
function cdr3_length_spectrum(sequences::AbstractVector{<:AbstractString})
    lengths = Int[]
    for seq in sequences
        cdr3 = extract_cdr3(seq)
        ismissing(cdr3) && continue
        push!(lengths, ncodeunits(String(cdr3)))
    end
    if isempty(lengths)
        result = with_provenance(DataFrame(length=Int[], count=Int[], frequency=Float64[]), "ImmunologyTable", "Immunology/cdr3_length_spectrum"; parameters=(sequence_count=length(sequences),))
        return _register_immunology_result!(_ctx, result, "cdr3_length_spectrum"; parameters=(n_sequences=length(sequences), n_unique_lengths=0))
    end
    tab = combine(groupby(DataFrame(length=lengths), :length), nrow => :count)
    sort!(tab, :length)
    tab[!, :frequency] = tab.count ./ sum(tab.count)
    result = with_provenance(tab, "ImmunologyTable", "Immunology/cdr3_length_spectrum"; parameters=(sequence_count=length(sequences),))
    _ctx = active_provenance_context()


    return _register_immunology_result!(_ctx, result, "cdr3_length_spectrum"; parameters=(n_sequences=length(sequences), n_unique_lengths=nrow(tab)))
end

"""
    tcrdist_like_matrix(cdr3_sequences)

Approximate TCRdist-style pairwise distance matrix with mismatch and length
penalties.
"""
function tcrdist_like_matrix(cdr3_sequences::AbstractVector{<:AbstractString}; mismatch_penalty::Real=1.0, gap_penalty::Real=1.5, threaded::Bool=true)
    seqs = uppercase.(String.(cdr3_sequences))
    n = length(seqs)
    D = zeros(Float64, n, n)

    function pair_distance(a::AbstractString, b::AbstractString)
        la = ncodeunits(a)
        lb = ncodeunits(b)
        lmin = min(la, lb)
        mismatch = 0
        @inbounds for i in 1:lmin
            mismatch += a[i] == b[i] ? 0 : 1
        end
        gap = abs(la - lb)
        return mismatch_penalty * mismatch + gap_penalty * gap
    end

    threaded_foreach(n, i -> begin
        @inbounds D[i, i] = 0.0
        for j in (i + 1):n
            d = pair_distance(seqs[i], seqs[j])
            D[i, j] = d
            D[j, i] = d
        end
    end; threaded=threaded)
    _ctx = active_provenance_context()


    return _register_immunology_result!(_ctx, D, "tcrdist_like_matrix"; parameters=(n_sequences=n, mismatch_penalty=Float64(mismatch_penalty), gap_penalty=Float64(gap_penalty)))
end

_tool_cmd(program::String, args::Vector{String}) = `$(program) $(args...)`

"""
    igblast_assign(fasta_path; out_path="igblast.out", dry_run=true)

Plan or run an IgBLAST V(D)J assignment command.
"""
function igblast_assign(fasta_path::AbstractString; out_path::AbstractString="igblast.out", igblast_bin::AbstractString="igblastn", germline_v::AbstractString="", germline_d::AbstractString="", germline_j::AbstractString="", organism::AbstractString="human", outfmt::Int=19, threads::Int=4, dry_run::Bool=true)
    args = String["-query", String(fasta_path), "-organism", String(organism), "-outfmt", string(outfmt), "-num_threads", string(threads), "-out", String(out_path)]
    isempty(germline_v) || append!(args, ["-germline_db_V", String(germline_v)])
    isempty(germline_d) || append!(args, ["-germline_db_D", String(germline_d)])
    isempty(germline_j) || append!(args, ["-germline_db_J", String(germline_j)])
    cmd = _tool_cmd(String(igblast_bin), args)
    if dry_run
        return (cmd=string(cmd), output=String(out_path), status=:planned)
    end
    run(cmd)

    return (cmd=string(cmd), output=String(out_path), status=:ok)
end

"""
    imgt_highvquest_manifest(sample_to_fasta)

Build a submission manifest for IMGT/HighV-QUEST batch uploads.
"""
function imgt_highvquest_manifest(sample_to_fasta::AbstractDict; species::AbstractString="Homo sapiens", receptor_type::AbstractString="TR")
    sample = String[]
    fasta_path = String[]
    for (k, v) in sample_to_fasta
        push!(sample, String(k))
        push!(fasta_path, String(v))
    end
    out = DataFrame(sample=sample, fasta_path=fasta_path)
    out[!, :species] = fill(String(species), nrow(out))
    out[!, :receptor_type] = fill(String(receptor_type), nrow(out))
    out[!, :exists_on_disk] = isfile.(out.fasta_path)

    return with_provenance(out, "ImmunologyTable", "Immunology/imgt_highvquest_manifest"; parameters=(sample_count=length(sample_to_fasta), species=species, receptor_type=receptor_type))
end

"""
    paired_clonotype_table(alpha_cdr3, beta_cdr3)

Create paired-chain clonotype frequencies for TCR/BCR analyses.
"""
function paired_clonotype_table(alpha_cdr3::AbstractVector{<:AbstractString}, beta_cdr3::AbstractVector{<:AbstractString}; sample_ids=nothing)
    length(alpha_cdr3) == length(beta_cdr3) || throw(DimensionMismatch("alpha and beta chains must have equal length"))
    pair = ["$(uppercase(String(a)))|$(uppercase(String(b)))" for (a, b) in zip(alpha_cdr3, beta_cdr3)]
    tab = combine(groupby(DataFrame(paired_clonotype=pair), :paired_clonotype), nrow => :n_cells)
    tab[!, :frequency] = tab.n_cells ./ max(sum(tab.n_cells), 1)
    sort!(tab, :n_cells, rev=true)

    if sample_ids === nothing
        return tab
    end

    length(sample_ids) == length(pair) || throw(DimensionMismatch("sample_ids must match chain counts"))
    detail = DataFrame(sample=String.(sample_ids), paired_clonotype=pair)
    per_sample = combine(groupby(detail, [:sample, :paired_clonotype]), nrow => :n_cells)

    return (
        summary=with_provenance(tab, "ImmunologyTable", "Immunology/paired_clonotype_table"; parameters=(pair_count=length(alpha_cdr3), sample_split=true)),
        per_sample=with_provenance(per_sample, "ImmunologyTable", "Immunology/paired_clonotype_table"; notes=["per-sample paired clonotype counts"], parameters=(pair_count=length(alpha_cdr3), sample_split=true)))
end

"""
    differential_vj_usage(v_genes, j_genes, groups; case_group=nothing)

Differential VJ segment usage by Fisher exact testing.
"""
function differential_vj_usage(v_genes::AbstractVector{<:AbstractString}, j_genes::AbstractVector{<:AbstractString}, groups::AbstractVector{<:AbstractString}; case_group=nothing)
    n = length(v_genes)
    n == length(j_genes) == length(groups) || throw(DimensionMismatch("v_genes, j_genes, and groups must have equal length"))
    pair = ["$(String(v_genes[i]))|$(String(j_genes[i]))" for i in 1:n]
    g = String.(groups)
    case = case_group === nothing ? first(sort!(unique(g))) : String(case_group)
    case_idx = findall(==(case), g)
    ctrl_idx = setdiff(collect(1:n), case_idx)
    isempty(case_idx) && throw(ArgumentError("case_group has no observations"))
    isempty(ctrl_idx) && throw(ArgumentError("need at least one control observation"))

    uniq_pairs = sort!(unique(pair))
    out = DataFrame(vj_pair=String[], case_count=Int[], control_count=Int[], odds_ratio=Float64[], pvalue=Float64[])
    for p in uniq_pairs
        a = count(i -> pair[i] == p, case_idx)
        b = length(case_idx) - a
        c = count(i -> pair[i] == p, ctrl_idx)
        d = length(ctrl_idx) - c
        pval = _fisher_right_tail(a, b, c, d)
        orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        push!(out, (p, a, c, orr, pval))
    end

    out[!, :padj] = sortperm(out.pvalue) |> (perm -> begin
        p = out.pvalue
        q = similar(p)
        m = length(p)
        running = 1.0
        for rank in m:-1:1
            idx = perm[rank]
            running = min(running, p[idx] * m / rank)
            q[idx] = running
        end
        clamp.(q, 0.0, 1.0)
    end)
    sort!(out, [:padj, :pvalue])

    return with_provenance(out, "ImmunologyTable", "Immunology/differential_vj_usage"; parameters=(observation_count=length(v_genes), group_count=length(unique(groups))))
end

# ---------------------------------------------------------------------------
# Somatic Hypermutation Rate (BCR)
# ---------------------------------------------------------------------------

"""
    somatic_hypermutation_rate(sequences, germline_sequences; clone_ids=nothing)

Estimate per-cell and per-clone somatic hypermutation (SHM) rates by
computing the pairwise mismatch fraction against provided germline sequences,
analogous to `alakazam::calcMutationLoad` in Bioconductor.

`sequences`         : observed BCR sequences (same length as germline_sequences)
`germline_sequences`: matched germline V-region sequences

Returns a `DataFrame` with `sequence_index`, `n_mutations`, `seq_length`,
`mutation_rate`, and (if `clone_ids` provided) per-clone summary.
"""
function somatic_hypermutation_rate(
    sequences::AbstractVector{<:AbstractString},
    germline_sequences::AbstractVector{<:AbstractString};
    clone_ids::Union{Nothing,AbstractVector}=nothing)
    n = length(sequences)
    n == length(germline_sequences) || throw(DimensionMismatch("sequences and germline_sequences must have equal length"))

    n_mut = Int[]
    seq_len = Int[]
    for i in 1:n
        obs = uppercase(String(sequences[i]))
        gl  = uppercase(String(germline_sequences[i]))
        L = min(ncodeunits(obs), ncodeunits(gl))
        muts = count(pos -> obs[pos] != gl[pos], 1:L)
        push!(n_mut, muts)
        push!(seq_len, L)
    end
    mut_rate = Float64.(n_mut) ./ max.(Float64.(seq_len), 1.0)

    locus_df = DataFrame(
        sequence_index = 1:n,
        n_mutations    = n_mut,
        seq_length     = seq_len,
        mutation_rate  = mut_rate)
    _ctx = active_provenance_context()

    if clone_ids === nothing
        return _register_immunology_result!(_ctx, locus_df, "somatic_hypermutation_rate"; parameters=(n_sequences=n, with_clones=false))
    end

    n == length(clone_ids) || throw(DimensionMismatch("clone_ids must match sequences length"))
    locus_df[!, :clone_id] = String.(clone_ids)
    clone_df = combine(
        groupby(locus_df, :clone_id),
        :mutation_rate => mean => :mean_mutation_rate,
        :mutation_rate => maximum => :max_mutation_rate,
        :n_mutations   => sum  => :total_mutations,
        nrow => :n_cells)
    sort!(clone_df, :mean_mutation_rate, rev=true)
    result = (cells=locus_df, clones=clone_df)
    _ctx = active_provenance_context()


    return _register_immunology_result!(_ctx, result, "somatic_hypermutation_rate"; parameters=(n_sequences=n, with_clones=true, n_clones=nrow(clone_df)))
end

# ---------------------------------------------------------------------------
# BCR Affinity Maturation Score (AID hotspot density)
# ---------------------------------------------------------------------------

"""
    bcr_affinity_maturation_score(sequences; hotspot_motifs=["WRC","GYW","WA","TW"])

Estimate BCR affinity maturation by counting AID (activation-induced
cytidine deaminase) hotspot motifs in each sequence, analogous to
`SHazaM::findHotspot` in Bioconductor.

Default motifs are IUPAC-encoded AID hotspots (WRC = [A/T]G[A/G]C, etc.).

Returns a `DataFrame` with `sequence`, `seq_length`,
`n_hotspots`, `hotspot_density`, `maturation_score`.
"""
function bcr_affinity_maturation_score(
    sequences::AbstractVector{<:AbstractString};
    hotspot_motifs::AbstractVector{<:AbstractString}=["WRC", "GYW", "WA", "TW"])
    # Expand IUPAC codes to regex character classes
    iupac = Dict(
        'W' => "[AT]", 'R' => "[AG]", 'Y' => "[CT]", 'S' => "[GC]",
        'K' => "[GT]", 'M' => "[AC]", 'B' => "[CGT]", 'D' => "[AGT]",
        'H' => "[ACT]", 'V' => "[ACG]", 'N' => "[ACGT]")
    function iupac_to_regex(motif::String)
        join([haskey(iupac, c) ? iupac[c] : string(c) for c in uppercase(motif)])
    end
    patterns = [Regex(iupac_to_regex(String(m))) for m in hotspot_motifs]

    out_seqs = String[]
    out_len  = Int[]
    out_hot  = Int[]
    out_dens = Float64[]
    out_mat  = Float64[]

    for seq in sequences
        s = uppercase(String(seq))
        L = ncodeunits(s)
        count_hot = 0
        for pat in patterns
            # Count overlapping occurrences
            pos = 1
            while pos <= L
                m = match(pat, s, pos)
                m === nothing && break
                count_hot += 1
                pos = m.offset + 1
            end
        end
        density = L > 0 ? count_hot / L : 0.0
        # Maturation score: sigmoid-scaled density
        mat_score = 1.0 / (1.0 + exp(-10.0 * (density - 0.05)))
        push!(out_seqs, s)
        push!(out_len, L)
        push!(out_hot, count_hot)
        push!(out_dens, density)
        push!(out_mat, mat_score)
    end

    return with_provenance(DataFrame(
        sequence          = out_seqs,
        seq_length        = out_len,
        n_hotspots        = out_hot,
        hotspot_density   = out_dens,
        maturation_score  = out_mat), "ImmunologyTable", "Immunology/bcr_affinity_maturation_score"; parameters=(sequence_count=length(sequences), motif_count=length(hotspot_motifs)))
end

# ---------------------------------------------------------------------------
# Antigen Specificity Clustering (GLIPH2 / CoNGA style)
# ---------------------------------------------------------------------------

"""
    antigen_specificity_clustering(cdr3_sequences; similarity_threshold=3.5, min_cluster_size=2)

Cluster CDR3 sequences into putative antigen-specific groups via single-linkage
hierarchical clustering on TCRdist-like distances, analogous to GLIPH2 and
CoNGA in Bioconductor repertoire workflows.

Returns a `DataFrame` with `sequence`, `cluster_id`, `cluster_size`,
and a `DataFrame` of `cluster_summary`.
"""
function antigen_specificity_clustering(
    cdr3_sequences::AbstractVector{<:AbstractString};
    similarity_threshold::Real=3.5,
    min_cluster_size::Int=2,
    mismatch_penalty::Real=1.0,
    gap_penalty::Real=1.5)
    seqs = uppercase.(String.(cdr3_sequences))
    n = length(seqs)
    n == 0 && return (
        cells=with_provenance(DataFrame(), "ImmunologyTable", "Immunology/antigen_specificity_clustering"; parameters=(sequence_count=0,)),
        cluster_summary=with_provenance(DataFrame(), "ImmunologyTable", "Immunology/antigen_specificity_clustering"; notes=["empty cluster summary"], parameters=(sequence_count=0,)))

    D = tcrdist_like_matrix(seqs; mismatch_penalty=mismatch_penalty, gap_penalty=gap_penalty, threaded=true)
    thresh = Float64(similarity_threshold)

    # Union-Find for single-linkage
    parent = collect(1:n)
    function find!(p, x)
        while p[x] != x
            p[x] = p[p[x]]
            x = p[x]
        end
        return x
    end
    function union!(p, x, y)
        rx = find!(p, x)
        ry = find!(p, y)
        rx != ry && (p[ry] = rx)
    end

    for i in 1:n, j in (i+1):n
        D[i, j] <= thresh && union!(parent, i, j)
    end

    labels = [find!(parent, i) for i in 1:n]
    # Renumber clusters
    unique_roots = sort!(unique(labels))
    root_map = Dict(r => k for (k, r) in enumerate(unique_roots))
    cluster_ids = [root_map[l] for l in labels]

    cell_df = DataFrame(sequence=seqs, cluster_id=cluster_ids)

    # Cluster summary
    summary = combine(
        groupby(cell_df, :cluster_id),
        nrow => :cluster_size,
        :sequence => (seqs -> join(seqs[1:min(3, length(seqs))], "; ")) => :example_sequences)
    filter!(r -> r.cluster_size >= min_cluster_size, summary)
    sort!(summary, :cluster_size, rev=true)

    # Mark cells not in qualifying clusters
    valid_clusters = Set(summary.cluster_id)
    cell_df[!, :in_cluster] = [c in valid_clusters for c in cell_df.cluster_id]

    return (
        cells=with_provenance(cell_df, "ImmunologyTable", "Immunology/antigen_specificity_clustering"; notes=["cell-to-cluster assignments"], parameters=(sequence_count=n, similarity_threshold=Float64(similarity_threshold), min_cluster_size=Int(min_cluster_size))),
        cluster_summary=with_provenance(summary, "ImmunologyTable", "Immunology/antigen_specificity_clustering"; notes=["cluster-level summary"], parameters=(sequence_count=n, similarity_threshold=Float64(similarity_threshold), min_cluster_size=Int(min_cluster_size))))
end

# ---------------------------------------------------------------------------
# Clonotype Trajectory (longitudinal tracking)
# ---------------------------------------------------------------------------

"""
    clonotype_trajectory(clonotypes, timepoints; sample_ids=nothing)

Track clonotype frequencies across ordered timepoints, analogous to
`immunarch::trackClonotypes` in Bioconductor workflows.

`clonotypes` : clonotype assignments for each cell.
`timepoints` : timepoint label for each cell (must be sortable).

Returns a wide-format `DataFrame` with one row per unique clonotype and
columns for each timepoint's frequency, plus `delta_first_last` and
`expansion_class`.
"""
function clonotype_trajectory(
    clonotypes::AbstractVector{<:AbstractString},
    timepoints::AbstractVector;
    sample_ids::Union{Nothing,AbstractVector}=nothing,
    min_max_frequency::Real=0.001)
    n = length(clonotypes)
    n == length(timepoints) || throw(DimensionMismatch("clonotypes and timepoints must have equal length"))

    tp_strings = string.(timepoints)
    sorted_tps = sort!(unique(tp_strings))

    df = DataFrame(clonotype=String.(clonotypes), timepoint=tp_strings)
    sample_ids !== nothing && (df[!, :sample] = String.(sample_ids))

    # Count per timepoint
    freq_df = combine(
        groupby(df, [:timepoint, :clonotype]),
        nrow => :n_cells)

    # Pivot: clonotype × timepoint
    all_clones = sort!(unique(freq_df.clonotype))
    wide = DataFrame(clonotype=all_clones)
    for tp in sorted_tps
        sub = filter(r -> r.timepoint == tp, freq_df)
        total = sum(sub.n_cells)
        freqs = Dict(String(r.clonotype) => r.n_cells / max(total, 1) for r in eachrow(sub))
        wide[!, Symbol("freq_$(tp)")] = [get(freqs, c, 0.0) for c in all_clones]
    end

    # Delta first → last
    first_col = Symbol("freq_$(sorted_tps[1])")
    last_col  = Symbol("freq_$(sorted_tps[end])")
    wide[!, :delta_first_last] = wide[!, last_col] .- wide[!, first_col]

    # Filter clonotypes that never reach min_max_frequency
    max_freq = [maximum(wide[i, Symbol("freq_$(tp)")] for tp in sorted_tps) for i in 1:nrow(wide)]
    filter_mask = max_freq .>= Float64(min_max_frequency)
    wide = wide[filter_mask, :]

    wide[!, :expansion_class] = [
        d > 0.01 ? "expanding" :
        d < -0.01 ? "contracting" : "stable"
        for d in wide.delta_first_last
    ]

    sort!(wide, :delta_first_last, rev=true)
    result = with_provenance(wide, "ImmunologyTable", "Immunology/clonotype_trajectory"; parameters=(clonotype_count=length(clonotypes), timepoint_count=length(sorted_tps)))
    _ctx = active_provenance_context()


    return _register_immunology_result!(_ctx, result, "clonotype_trajectory"; parameters=(clonotype_count=length(clonotypes), timepoint_count=length(sorted_tps), n_tracked=nrow(wide)))
end

# ---------------------------------------------------------------------------
# V-Gene Family Usage
# ---------------------------------------------------------------------------

"""
    v_gene_family_usage(v_genes; family_pattern=r"^(TRB|TRA|IGH|IGL|IGK)V([0-9]+)")

Summarise V-gene usage at the gene-family level (e.g. TRBV12 family),
analogous to `immunarch::geneUsage` with family grouping.

Returns a `DataFrame` with `v_family`, `n_cells`, `frequency`.
"""
function v_gene_family_usage(
    v_genes::AbstractVector{<:AbstractString};
    family_pattern::Regex=r"^([A-Z]+V[0-9]+)",
    group_col::Union{Nothing,AbstractVector}=nothing)
    families = String[]
    for vg in v_genes
        m = match(family_pattern, uppercase(String(vg)))
        push!(families, m === nothing ? "unknown" : m.match)
    end

    if group_col !== nothing
        length(group_col) == length(v_genes) || throw(DimensionMismatch("group_col must match v_genes length"))
        df = DataFrame(v_family=families, group=String.(group_col))
        counts = combine(groupby(df, [:group, :v_family]), nrow => :n_cells)
        counts = transform(groupby(counts, :group), :n_cells => (c -> c ./ max(sum(c), 1)) => :frequency)
        sort!(counts, [:group, :frequency], rev=[false, true])
        return with_provenance(counts, "ImmunologyTable", "Immunology/v_gene_family_usage"; parameters=(gene_count=length(v_genes), grouped=true))
    end

    counts = combine(groupby(DataFrame(v_family=families), :v_family), nrow => :n_cells)
    total = sum(counts.n_cells)
    counts[!, :frequency] = counts.n_cells ./ max(total, 1)
    sort!(counts, :frequency, rev=true)
    return with_provenance(counts, "ImmunologyTable", "Immunology/v_gene_family_usage"; parameters=(gene_count=length(v_genes), grouped=false))
end

# ---------------------------------------------------------------------------
# HLA Type Enrichment
# ---------------------------------------------------------------------------

"""
    hla_type_enrichment(clonotypes, hla_alleles; case_clonotypes=nothing)

Test for enrichment of HLA alleles within expanded clonotypes using Fisher's
exact test, analogous to HLA-EMMA and immunarch HLA enrichment modules.

`clonotypes`    : clonotype label per cell.
`hla_alleles`   : HLA allele string per cell (e.g. "HLA-A*02:01").
`case_clonotypes`: Set of expanded clonotype labels (or nothing to auto-detect top 10%).

Returns a `DataFrame` sorted by `pvalue` with `hla_allele`, `case_count`,
`control_count`, `odds_ratio`, `pvalue`, `padj`.
"""
function hla_type_enrichment(
    clonotypes::AbstractVector{<:AbstractString},
    hla_alleles::AbstractVector{<:AbstractString};
    case_clonotypes::Union{Nothing,AbstractVector}=nothing,
    top_fraction::Real=0.1)
    n = length(clonotypes)
    n == length(hla_alleles) || throw(DimensionMismatch("clonotypes and hla_alleles must have equal length"))

    if case_clonotypes === nothing
        counts = combine(groupby(DataFrame(c=String.(clonotypes)), :c), nrow => :n)
        sort!(counts, :n, rev=true)
        k = max(1, round(Int, nrow(counts) * Float64(top_fraction)))
        case_set = Set(counts.c[1:k])
    else
        case_set = Set(String.(case_clonotypes))
    end

    is_case = [String(c) in case_set for c in clonotypes]
    alleles  = sort!(unique(String.(hla_alleles)))

    out = DataFrame(hla_allele=String[], case_count=Int[], control_count=Int[], odds_ratio=Float64[], pvalue=Float64[])
    for allele in alleles
        a = count(i -> is_case[i] && String(hla_alleles[i]) == allele, 1:n)
        b = count(i -> is_case[i] && String(hla_alleles[i]) != allele, 1:n)
        c = count(i -> !is_case[i] && String(hla_alleles[i]) == allele, 1:n)
        d = count(i -> !is_case[i] && String(hla_alleles[i]) != allele, 1:n)
        or = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        pval = _fisher_right_tail(a, b, c, d)
        push!(out, (allele, a, c, or, pval))
    end
    sort!(out, :pvalue)
    m = nrow(out)
    padj = similar(out.pvalue)
    running = 1.0
    for rank in m:-1:1
        running = min(running, out.pvalue[rank] * m / rank)
        padj[rank] = running
    end
    out[!, :padj] = clamp.(padj, 0.0, 1.0)
    return with_provenance(out, "ImmunologyTable", "Immunology/hla_type_enrichment"; parameters=(observation_count=n, allele_count=length(alleles)))
end

# ---------------------------------------------------------------------------
# CDR3 Physicochemical Properties
# ---------------------------------------------------------------------------

"""
    cdr3_physicochemical_properties(sequences)

Compute amino-acid physicochemical properties for CDR3 regions:
length, charge, hydrophobicity (Kyte-Doolittle), aromaticity,
and instability index (Guruprasad), analogous to `Peptides.jl` /
`alakazam::aminoAcidProperties`.

Returns a `DataFrame` with one row per sequence.
"""
function cdr3_physicochemical_properties(sequences::AbstractVector{<:AbstractString})
    # Amino acid properties
    kd_hydrophobicity = Dict(
        'I' => 4.5, 'V' => 4.2, 'L' => 3.8, 'F' => 2.8, 'C' => 2.5,
        'M' => 1.9, 'A' => 1.8, 'G' => -0.4, 'T' => -0.7, 'S' => -0.8,
        'W' => -0.9, 'Y' => -1.3, 'P' => -1.6, 'H' => -3.2, 'E' => -3.5,
        'Q' => -3.5, 'D' => -3.5, 'N' => -3.5, 'K' => -3.9, 'R' => -4.5)
    charge_map = Dict('R' => 1, 'K' => 1, 'H' => 0, 'D' => -1, 'E' => -1)
    aromatic    = Set(['F', 'W', 'Y', 'H'])
    # Simplified DIWV instability (Guruprasad pairs)
    instab_pairs = Dict("WW"=>1.0,"WC"=>1.0,"WM"=>24.68,"WH"=>24.68,
                        "CW"=>24.68,"CH"=>24.68,"CC"=>1.0,"CF"=>54.96,
                        "KK"=>1.0,"QQ"=>1.0,"RR"=>1.0)

    out_seq  = String[]
    out_cdr3 = String[]
    out_len  = Int[]
    out_chg  = Float64[]
    out_hyd  = Float64[]
    out_aro  = Float64[]
    out_inst = Float64[]

    for seq in sequences
        cdr3_m = extract_cdr3(seq)
        cdr3 = ismissing(cdr3_m) ? uppercase(String(seq)) : String(cdr3_m)
        aa = collect(uppercase(cdr3))
        L = length(aa)

        charge = sum(get(charge_map, c, 0) for c in aa; init=0)
        hydro  = L > 0 ? mean(get(kd_hydrophobicity, c, 0.0) for c in aa) : 0.0
        aro    = L > 0 ? count(c -> c in aromatic, aa) / L : 0.0
        inst   = 0.0
        for i in 1:(L - 1)
            key = string(aa[i]) * string(aa[i+1])
            inst += get(instab_pairs, key, 1.0)
        end
        inst_idx = L >= 2 ? 10.0 * inst / L : 0.0

        push!(out_seq, String(seq))
        push!(out_cdr3, cdr3)
        push!(out_len, L)
        push!(out_chg, Float64(charge))
        push!(out_hyd, hydro)
        push!(out_aro, aro)
        push!(out_inst, inst_idx)
    end

    return with_provenance(DataFrame(
        sequence        = out_seq,
        cdr3            = out_cdr3,
        cdr3_length     = out_len,
        net_charge      = out_chg,
        hydrophobicity  = out_hyd,
        aromaticity     = out_aro,
        instability_idx = out_inst), "ImmunologyTable", "Immunology/cdr3_physicochemical_properties"; parameters=(sequence_count=length(sequences),))
end

# ---------------------------------------------------------------------------
# Network Centrality for Clonotypes
# ---------------------------------------------------------------------------

"""
    network_centrality_clonotypes(cdr3_sequences; similarity_threshold=3.5)

Compute degree and closeness centrality metrics for each clonotype in the
TCRdist/BCRdist similarity network, analogous to network analyses in
CoNGA and `scirpy` / `immunarch`.

Returns a `DataFrame` with `sequence`, `degree`, `closeness_centrality`,
`hub_score`.
"""
function network_centrality_clonotypes(
    cdr3_sequences::AbstractVector{<:AbstractString};
    similarity_threshold::Real=3.5,
    mismatch_penalty::Real=1.0,
    gap_penalty::Real=1.5)
    seqs = uppercase.(String.(cdr3_sequences))
    n = length(seqs)
    n == 0 && return with_provenance(DataFrame(), "ImmunologyTable", "Immunology/network_centrality_clonotypes"; parameters=(sequence_count=0,))

    D = tcrdist_like_matrix(seqs; mismatch_penalty=mismatch_penalty, gap_penalty=gap_penalty, threaded=true)
    thresh = Float64(similarity_threshold)

    # Adjacency (binary)
    A = Float64.(D .<= thresh)
    A[diagind(A)] .= 0.0

    degree = vec(sum(A, dims=2))
    # Closeness: sum of inverse distances (harmonic centrality)
    harmonic = zeros(Float64, n)
    for i in 1:n
        for j in 1:n
            i == j && continue
            harmonic[i] += D[i, j] > 0 ? 1.0 / D[i, j] : 0.0
        end
        harmonic[i] /= max(n - 1, 1)
    end

    # Hub score: approx. HITS — proportional to degree^2 contribution
    hub = degree ./ max(maximum(degree), 1.0)

    return with_provenance(DataFrame(
        sequence            = seqs,
        degree              = Int.(round.(degree)),
        closeness_centrality = harmonic,
        hub_score           = hub), "ImmunologyTable", "Immunology/network_centrality_clonotypes"; parameters=(sequence_count=n, similarity_threshold=Float64(similarity_threshold)))
end

# ---------------------------------------------------------------------------
# Convergent Clonotype Detection
# ---------------------------------------------------------------------------

"""
    convergent_clonotype_detection(sequences, v_genes, j_genes; max_hamming=1, cdr3_col=nothing)

Identify convergent (public) clonotypes — identical or near-identical CDR3
sequences arising from different V(D)J recombinations, analogous to
`scirpy.tl.define_clonotypes` with convergence mode.

Returns a `DataFrame` with `group_id`, `n_sequences`, `sequences`,
`v_genes`, `j_genes`, `is_convergent`.
"""
function convergent_clonotype_detection(
    sequences::AbstractVector{<:AbstractString},
    v_genes::AbstractVector{<:AbstractString},
    j_genes::AbstractVector{<:AbstractString};
    max_hamming::Int=1,
    min_group_size::Int=2)
    n = length(sequences)
    n == length(v_genes) == length(j_genes) ||
        throw(DimensionMismatch("sequences, v_genes, j_genes must have equal length"))

    seqs = uppercase.(String.(sequences))
    v    = String.(v_genes)
    j    = String.(j_genes)

    # Group by CDR3 similarity (within max_hamming)
    parent = collect(1:n)
    function find!(p, x)
        while p[x] != x; p[x] = p[p[x]]; x = p[x]; end; x
    end
    function union!(p, x, y)
        rx = find!(p, x); ry = find!(p, y)
        rx != ry && (p[ry] = rx)
    end

    for i in 1:n, j_idx in (i+1):n
        si = seqs[i]; sj = seqs[j_idx]
        li = ncodeunits(si); lj = ncodeunits(sj)
        lmin = min(li, lj)
        gap  = abs(li - lj)
        gap > max_hamming && continue
        mm = count(pos -> si[pos] != sj[pos], 1:lmin)
        (mm + gap) <= max_hamming && union!(parent, i, j_idx)
    end

    labels = [find!(parent, i) for i in 1:n]
    root_map = Dict(r => k for (k, r) in enumerate(sort!(unique(labels))))
    gids = [root_map[l] for l in labels]

    df = DataFrame(sequence=seqs, v_gene=v, j_gene=j, group_id=gids)
    summary = combine(
        groupby(df, :group_id),
        nrow => :n_sequences,
        :sequence => (x -> join(unique(x), "; ")) => :sequences,
        :v_gene   => (x -> join(unique(x), "; ")) => :v_genes,
        :j_gene   => (x -> join(unique(x), "; ")) => :j_genes)
    summary[!, :n_unique_vj] = [
        length(unique(split(row.v_genes, "; ") .* "|" .* split(row.j_genes, "; ")))
        for row in eachrow(summary)
    ]
    summary[!, :is_convergent] = summary.n_sequences .>= min_group_size .&&
                                  summary.n_unique_vj .> 1
    filter!(r -> r.n_sequences >= min_group_size, summary)
    sort!(summary, :n_sequences, rev=true)
    _ctx = active_provenance_context()


    return _register_immunology_result!(_ctx, summary, "convergent_clonotype_detection"; parameters=(n_sequences=n, max_hamming=max_hamming, min_group_size=min_group_size, n_groups=nrow(summary)))
end

# ---------------------------------------------------------------------------
# Lineage Tree from Clone Sequences (simple parsimony)
# ---------------------------------------------------------------------------

"""
    lineage_tree_from_clones(sequences, germline; clone_ids=nothing)

Build a minimum-spanning-tree-based B-cell lineage tree from a set of
BCR sequences and their inferred germline, analogous to `alakazam::buildPhylipLineage`
in Bioconductor.

Uses pairwise Hamming distance to build an approximate MST via Prim's algorithm.

Returns `(edges=DataFrame, mst_distance_sum)`.
"""
function lineage_tree_from_clones(
    sequences::AbstractVector{<:AbstractString},
    germline::AbstractString;
    clone_ids::Union{Nothing,AbstractVector}=nothing)
    n = length(sequences)
    n >= 2 || throw(ArgumentError("need at least 2 sequences to build a lineage tree"))

    all_seqs = vcat([uppercase(String(germline))], uppercase.(String.(sequences)))
    ids = clone_ids !== nothing ? vcat(["germline"], String.(clone_ids)) :
          vcat(["germline"], ["seq_$(i)" for i in 1:n])
    m = length(all_seqs)

    # Pairwise Hamming (with gap penalty for length differences)
    function hamming(a, b)
        la = ncodeunits(a); lb = ncodeunits(b)
        lmin = min(la, lb)
        mm = count(p -> a[p] != b[p], 1:lmin)
        mm + abs(la - lb)
    end

    dist = [Float64(hamming(all_seqs[i], all_seqs[j])) for i in 1:m, j in 1:m]

    # Prim's MST from germline (node 1)
    in_tree = falses(m)
    in_tree[1] = true
    parent_node = zeros(Int, m)
    edge_weight = fill(Inf, m)
    edge_weight[1] = 0.0
    for v in 2:m
        edge_weight[v] = dist[1, v]
        parent_node[v] = 1
    end
    mst_edges = Tuple{Int,Int,Float64}[]

    for _ in 1:(m - 1)
        # Cheapest edge crossing the cut
        best_v = 0; best_w = Inf
        for v in 1:m
            !in_tree[v] && edge_weight[v] < best_w && (best_w = edge_weight[v]; best_v = v)
        end
        best_v == 0 && break
        in_tree[best_v] = true
        push!(mst_edges, (parent_node[best_v], best_v, best_w))
        for v in 1:m
            !in_tree[v] && dist[best_v, v] < edge_weight[v] &&
                (edge_weight[v] = dist[best_v, v]; parent_node[v] = best_v)
        end
    end

    edge_df = DataFrame(
        parent     = [ids[e[1] == 0 ? 1 : e[1]] for e in mst_edges],
        child      = [ids[e[2]] for e in mst_edges],
        hamming_dist = [e[3] for e in mst_edges])
    total_dist = sum(e[3] for e in mst_edges)
    result = (edges=edge_df, mst_distance_sum=total_dist)
    _ctx = active_provenance_context()


    return _register_immunology_result!(_ctx, result, "lineage_tree_from_clones"; parameters=(n_sequences=n, total_dist=total_dist, n_edges=length(mst_edges)))
end

# ---------------------------------------------------------------------------
# Epitope Binding Profile
# ---------------------------------------------------------------------------

"""
    epitope_binding_profile(peptides, known_epitopes; allele="HLA-A*02:01", max_hamming=1)

Match peptides against a known epitope database and compute binding probability
via anchor-position scoring + sequence similarity, analogous to
`mhcflurry` + `TCRex` output integration.

`known_epitopes`: `AbstractVector` of known binding peptide strings.

Returns a `DataFrame` with `peptide`, `allele`, `mhc_score`, `nearest_epitope`,
`epitope_hamming`, `combined_binding_score`, `predicted_binder`.
"""
function epitope_binding_profile(
    peptides::AbstractVector{<:AbstractString},
    known_epitopes::AbstractVector{<:AbstractString};
    allele::String="HLA-A*02:01",
    max_hamming::Int=1)
    mhc_scores = mhcflurry_like_scores(peptides; allele=allele)
    epi_vec = uppercase.(String.(known_epitopes))

    function min_hamming(pep, epi_list)
        L = ncodeunits(pep)
        best_d = typemax(Int)
        best_e = ""
        for e in epi_list
            le = ncodeunits(e)
            le == L || continue
            d = count(p -> pep[p] != e[p], 1:L)
            d < best_d && (best_d = d; best_e = e)
        end
        return best_e, best_d
    end

    out = copy(mhc_scores)
    nearest_epi  = String[]
    epi_hamming  = Int[]
    combined     = Float64[]
    binder       = Bool[]

    for row in eachrow(mhc_scores)
        p = uppercase(String(row.peptide))
        ne, nhd = min_hamming(p, epi_vec)
        push!(nearest_epi, ne)
        push!(epi_hamming, nhd == typemax(Int) ? -1 : nhd)
        epi_sim = nhd <= max_hamming ? 1.0 : 0.0
        comb_score = 0.6 * row.score + 0.4 * epi_sim
        push!(combined, comb_score)
        push!(binder, row.binder || epi_sim > 0)
    end

    out[!, :nearest_epitope]       = nearest_epi
    out[!, :epitope_hamming]       = epi_hamming
    out[!, :combined_binding_score] = combined
    out[!, :predicted_binder]      = binder
    sort!(out, :combined_binding_score, rev=true)
    _ctx = active_provenance_context()


    return _register_immunology_result!(_ctx, out, "epitope_binding_profile"; parents=provenance_parent_ids(peptides), parameters=(peptide_count=length(peptides), epitope_count=length(known_epitopes), allele=allele, max_hamming=max_hamming))
end

end
