# ==============================================================================
# crispr.jl — CRISPR-Cas guide design & off-target analysis
#
# Fully self-contained module for:
#   - Guide RNA design for SpCas9, Cas12a, CasX, Cas13
#   - On-target efficiency scoring (Rule Set 2 / Doench 2016-like)
#   - Off-target enumeration and scoring (CFD, MIT score)
#   - PAM detection (NGG, TTTV, TTTN, etc.)
#   - Genome-wide guide library design
#   - Base editor window analysis
#   - Prime editing guide design (pegRNA)
#   - CRISPR screen MAGeCK-like analysis
#   - HDR template design
#   - Indel prediction (CRISPR-ML Lindel-like)
#
# References:
#   - Doench et al. (2016) Nature Biotechnology 34:184-191 (Rule Set 2)
#   - Hsu et al. (2013) Nature Biotechnology 31:827-832 (MIT score)
#   - Doench et al. (2014) Nature Biotechnology 32:1262-1267 (CFD)
#   - Anzalone et al. (2019) Nature 576:149-157 (prime editing)
#   - Li et al. (2014) Genome Biology 15:554 (MAGeCK)
#   - Shen et al. (2018) Nature Methods 15:523-525 (Lindel)
# ==============================================================================

module CRISPR

using DataFrames
using Statistics
using LinearAlgebra
using Random

# Import biotypes for type-safe sequence handling
using ..BioToolkit: AASeq, AminoAcidAlphabet, BioSequence, DNAAlphabet, DNASeq, SummarizedExperiment, assay, colData
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!

export GuideRNA, CRISPRSystem, OffTarget, EditingWindow
export design_guides, score_on_target, find_pam_sites, enumerate_off_targets
export cfd_score, mit_score, guide_gc_content
export design_base_editor_guides, analyze_editing_window
export design_pegrna, prime_editing_guide_score
export crispr_screen_analysis, mageck_like_test
export design_library, library_coverage_stats
export design_hdr_template
export predict_indels, indel_distribution
export filter_guides, rank_guides
export guide_specificity_score, genome_wide_offtarget_summary

@inline function _register_crispr_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

"""
    CRISPRSystem

Specification of a CRISPR-Cas system including PAM, protospacer length, and cut position.
"""
struct CRISPRSystem
    name::String
    pam::String             # PAM sequence (e.g. "NGG", "TTTV")
    pam_location::Symbol    # :three_prime or :five_prime
    protospacer_length::Int
    cut_position::Int       # bp upstream of PAM (for SpCas9 = 3)
    strand_specific::Bool
end

# Predefined CRISPR systems
const SpCas9  = CRISPRSystem("SpCas9",  "NGG",  :three_prime, 20, 3, false)
const SaCas9  = CRISPRSystem("SaCas9",  "NNGRRT", :three_prime, 21, 3, false)
const Cas12a  = CRISPRSystem("Cas12a",  "TTTV", :five_prime,  23, 18, false)
const CasX    = CRISPRSystem("CasX",   "TTCN",  :five_prime,  20, 10, false)
const SpRY    = CRISPRSystem("SpRY",    "NRN",  :three_prime, 20, 3, false)  # near-PAMless

export SpCas9, SaCas9, Cas12a, CasX, SpRY

"""
    GuideRNA

A designed guide RNA with target information and predicted scores.
"""
struct GuideRNA
    spacer::String
    chromosome::String
    position::Int
    strand::Int8
    pam::String
    on_target_score::Float64    # Rule Set 2 / Doench-like
    gc_content::Float64
    off_target_count::Int
    specificity_score::Float64  # 0–1 (higher = more specific)
    passes_filters::Bool
    system::String
end

"""
    OffTarget

A predicted off-target site for a guide RNA.
"""
struct OffTarget
    guide::String
    chromosome::String
    position::Int
    strand::Int8
    sequence::String
    mismatches::Int
    bulges::Int
    cfd_score::Float64
    mit_score::Float64
    pam::String
    gene_context::String
end

"""
    EditingWindow

Defines the window of nucleotides amenable to base editing.
"""
struct EditingWindow
    editor_name::String
    window_start::Int   # 1-indexed from PAM-distal end
    window_end::Int
    editable_base::Char
    edited_base::Char
    bystander_positions::Vector{Int}
end

# Base editor definitions
const BE3    = EditingWindow("BE3",    4, 8, 'C', 'T', [3,4,5,6,7,8,9])
const ABE8e  = EditingWindow("ABE8e",  4, 8, 'A', 'G', [3,4,5,6,7,8,9])
const CBE4max = EditingWindow("CBE4max",3,9,'C','T', [2,3,4,5,6,7,8,9,10])

export BE3, ABE8e, CBE4max

# ---------------------------------------------------------------------------
# PAM matching
# ---------------------------------------------------------------------------

const _IUPAC = Dict(
    'N'=>"ACGT",'R'=>"AG",'Y'=>"CT",'S'=>"GC",'W'=>"AT",
    'K'=>"GT",'M'=>"AC",'B'=>"CGT",'D'=>"AGT",'H'=>"ACT",'V'=>"ACG",
    'A'=>"A",'C'=>"C",'G'=>"G",'T'=>"T"
)

"""
    pam_matches(pam_pattern, sequence_segment) → Bool

Check if a DNA sequence segment matches a IUPAC-encoded PAM pattern.
"""
function pam_matches(pam_pattern::AbstractString, seq_seg::AbstractString)
    length(pam_pattern) == length(seq_seg) || return false
    for (p, s) in zip(uppercase(pam_pattern), uppercase(seq_seg))
        s in get(_IUPAC, p, "") || return false
    end

    return true
end

"""
    find_pam_sites(sequence, system; chromosome="chr1") → DataFrame

Find all PAM sites for a CRISPR system in both strands of a sequence.
Accepts `AbstractString` or `BioSequence{DNAAlphabet}` (type-safe).
Returns a DataFrame with protospacer + PAM positions.
"""
function find_pam_sites(
    sequence::BioSequence{DNAAlphabet},
    system::CRISPRSystem;
    chromosome::String="chr1")
    seq  = uppercase(String(sequence))
    n    = length(seq)
    plen = length(system.pam)
    glen = system.protospacer_length
    rows = NamedTuple[]

    # Forward strand
    if system.pam_location == :three_prime
        for i in 1:(n - glen - plen + 1)
            pam_seq = seq[i+glen : i+glen+plen-1]
            if pam_matches(system.pam, pam_seq)
                spacer = seq[i : i+glen-1]
                push!(rows, (chrom=chromosome, position=i, strand=Int8(1),
                    spacer=spacer, pam=pam_seq))
            end
        end
    else  # five_prime PAM (Cas12a)
        for i in (plen+1):(n - glen + 1)
            pam_seq = seq[i-plen : i-1]
            if pam_matches(system.pam, pam_seq)
                spacer = seq[i : i+glen-1]
                push!(rows, (chrom=chromosome, position=i, strand=Int8(1),
                    spacer=spacer, pam=pam_seq))
            end
        end
    end

    # Reverse complement strand
    rc = _reverse_complement_str(seq)
    rc_offset = n
    if system.pam_location == :three_prime
        for i in 1:(length(rc) - glen - plen + 1)
            pam_seq = rc[i+glen : i+glen+plen-1]
            if pam_matches(system.pam, pam_seq)
                spacer = rc[i : i+glen-1]
                pos = rc_offset - (i + glen - 1) + 1  # convert back to fwd coords
                push!(rows, (chrom=chromosome, position=pos, strand=Int8(-1),
                    spacer=spacer, pam=pam_seq))
            end
        end
    else
        for i in (plen+1):(length(rc) - glen + 1)
            pam_seq = rc[i-plen : i-1]
            if pam_matches(system.pam, pam_seq)
                spacer = rc[i : i+glen-1]
                pos = rc_offset - (i + glen - 1) + 1
                push!(rows, (chrom=chromosome, position=pos, strand=Int8(-1),
                    spacer=spacer, pam=pam_seq))
            end
        end
    end

    result = DataFrame(rows)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "find_pam_sites"; parents=provenance_parent_ids(sequence), parameters=(chromosome=chromosome, system=system.name, row_count=nrow(result)))
end

function _reverse_complement_str(s::AbstractString)
    comp = Dict('A'=>'T','T'=>'A','G'=>'C','C'=>'G','N'=>'N',
                'R'=>'Y','Y'=>'R','S'=>'S','W'=>'W','K'=>'M','M'=>'K')
    return String(reverse([get(comp, uppercase(c), 'N') for c in s]))
end

# ---------------------------------------------------------------------------
# GC content
# ---------------------------------------------------------------------------

"""
    guide_gc_content(spacer) → Float64

Accepts `AbstractString` or `BioSequence{DNAAlphabet}`.
"""
function guide_gc_content(seq::BioSequence{DNAAlphabet})
    gc = count(c -> c in ('G','C','g','c'), String(seq))
    result = gc / max(length(seq), 1)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "guide_gc_content"; parents=provenance_parent_ids(seq), parameters=(length=length(seq), gc=result))
end

guide_gc_content(seq::AbstractString) = guide_gc_content(DNASeq(seq))

# ---------------------------------------------------------------------------
# On-target efficiency scoring (Doench 2016 Rule Set 2 simplified)
# ---------------------------------------------------------------------------

# Position-specific single-nucleotide weights derived from Doench 2016
# These are signed weights; sum gives a logit offset from intercept
const _RS2_SINGLE_NUC_WEIGHTS = Dict(
    (1,'G')=>-0.2753771,(2,'A')=>-0.3238875,(2,'C')=>0.17212887,(3,'C')=>-0.1006662,
    (4,'C')=>-0.2018029,(4,'T')=>-0.1747400,(5,'A')=>0.20932776,(5,'C')=>-0.17166690,
    (6,'C')=>0.11385978,(7,'C')=>-0.0671806,(8,'A')=>0.0694723,(8,'T')=>0.18765571,
    (9,'C')=>0.07781066,(9,'G')=>-0.4453400,(10,'C')=>0.27529108,(10,'G')=>-0.22163640,
    (11,'A')=>0.2093,(11,'G')=>-0.09849019,(12,'C')=>-0.29010478,(12,'T')=>0.1564111,
    (13,'G')=>0.07606142,(13,'T')=>-0.2130060736,(14,'C')=>0.1228724,(14,'T')=>-0.10466540,
    (15,'G')=>0.06421542,(15,'T')=>0.0855514,(16,'G')=>0.0498791,(16,'T')=>-0.05312809,
    (17,'C')=>-0.13640294,(17,'G')=>0.1379505,(18,'A')=>0.16827566,(18,'C')=>-0.09963975,
    (19,'A')=>0.28722890,(19,'G')=>-0.21012358,(20,'A')=>-0.06779897,(20,'G')=>0.11098105)
"""
    score_on_target(spacer; pam_context="NGG") → Float64

Estimate on-target cutting efficiency (0–1) using a simplified Rule Set 2
position-specific weight matrix.

Reference: Doench et al. (2016) Nature Biotechnology 34:184-191.
"""
function score_on_target(spacer::BioSequence{DNAAlphabet}; pam_context::String="NGG")
    s = uppercase(String(spacer))
    length(s) < 20 && return 0.0

    intercept = 0.5977
    score = intercept
    for i in 1:min(20, length(s))
        nt = s[i]
        score += get(_RS2_SINGLE_NUC_WEIGHTS, (i, nt), 0.0)
    end

    # GC content penalty
    gc = guide_gc_content(s)
    if gc < 0.25 || gc > 0.75
        score -= 0.15
    end

    # Poly-T penalty (reads off Pol III)
    occursin("TTTT", s) && (score -= 0.20)

    # First position G preference
    s[1] == 'G' && (score += 0.05)

    result = clamp(1.0 / (1.0 + exp(-score)), 0.0, 1.0)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "score_on_target"; parents=provenance_parent_ids(spacer), parameters=(pam_context=pam_context, gc_content=gc, score=result))
end

score_on_target(spacer::AbstractString; pam_context::String="NGG") = score_on_target(DNASeq(spacer); pam_context=pam_context)

# ---------------------------------------------------------------------------
# Off-target scoring: CFD (Cutting Frequency Determination)
# ---------------------------------------------------------------------------

# Simplified mismatch position weights for CFD (Doench 2014)
# Higher weight = more tolerant of mismatch (worse specificity)
const _CFD_MISMATCH_WEIGHTS = [
    0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0,
    0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828,
    0.615, 0.804, 0.685, 0.583
]  # positions 1-20 from PAM-distal end

"""
    cfd_score(guide, off_target) → Float64

Compute the Cutting Frequency Determination (CFD) score between a guide and
an off-target site. Returns a value in [0,1] (1 = perfect match).

Reference: Doench et al. (2014) Nature Biotechnology 32:1262-1267.
"""
function cfd_score(guide::BioSequence{DNAAlphabet}, off_target::BioSequence{DNAAlphabet})
    g = uppercase(String(guide))
    o = uppercase(String(off_target))
    length(g) == length(o) || return 0.0
    n = length(g)
    score = 1.0
    for (i, (gc, oc)) in enumerate(zip(g, o))
        gc == oc && continue
        # Position weight (1-indexed from PAM-distal end, 20 = PAM-proximal)
        pos = min(i, length(_CFD_MISMATCH_WEIGHTS))
        w = _CFD_MISMATCH_WEIGHTS[pos]
        score *= max(1.0 - w, 0.0)
    end
    result = clamp(score, 0.0, 1.0)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "cfd_score"; parents=provenance_parent_ids(guide, off_target), parameters=(length=length(g), score=result))
end

# ---------------------------------------------------------------------------
# Off-target scoring: MIT score (Hsu 2013)
# ---------------------------------------------------------------------------

const _MIT_WEIGHTS = [0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0,
                      0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828,
                      0.615, 0.804, 0.685, 0.583]

"""
    mit_score(guide, off_targets) → Float64

Compute the MIT specificity score for a guide against a list of off-target sequences.
Each off-target reduces the score. Returns a value in [0,100].

Reference: Hsu et al. (2013) Nature Biotechnology 31:827-832.
"""
function mit_score(guide::BioSequence{DNAAlphabet}, off_targets::AbstractVector{<:BioSequence{DNAAlphabet}})
    g = uppercase(String(guide))
    isempty(off_targets) && return 100.0
    total_score = 0.0
    for ot in off_targets
        o = uppercase(String(ot))
        length(g) == length(o) || continue
        mm_positions = findall(i -> g[i] != o[i], 1:length(g))
        n_mm = length(mm_positions)
        n_mm == 0 && continue   # exact match (ignore self)
        # Product of individual weights
        s = prod(1.0 - _MIT_WEIGHTS[min(p, 20)] for p in mm_positions; init=1.0)
        # Distance penalty for clustered mismatches
        if n_mm > 1
            span = maximum(mm_positions) - minimum(mm_positions) + 1
            d_penalty = 1.0 - (n_mm - 1.0) / (span * 1.5)
        else
            d_penalty = 1.0
        end
        total_score += clamp(s * d_penalty, 0.0, 1.0)
    end
    result = clamp(100.0 * (1.0 - total_score / max(length(off_targets), 1)), 0.0, 100.0)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "mit_score"; parents=provenance_parent_ids(guide, off_targets), parameters=(off_target_count=length(off_targets), score=result))
end

# ---------------------------------------------------------------------------
# Off-target enumeration (in-silico, from a target sequence)
# ---------------------------------------------------------------------------

"""
    enumerate_off_targets(guide, genome_seq; max_mismatches=3, chromosome="chr1") → DataFrame

Enumerate off-target sites. Accepts `AbstractString` or `BioSequence{DNAAlphabet}`.
"""
function enumerate_off_targets(
    guide::BioSequence{DNAAlphabet},
    genome_seq::BioSequence{DNAAlphabet};
    max_mismatches::Int=3,
    chromosome::String="chr1",
    system::CRISPRSystem=SpCas9)
    g    = uppercase(String(guide))
    glen = length(g)
    seq  = uppercase(String(genome_seq))
    n    = length(seq)
    plen = length(system.pam)
    ots  = NamedTuple[]

    for strand_sign in (1, -1)
        working = strand_sign == 1 ? seq : _reverse_complement_str(seq)
        nw = length(working)
        for i in 1:(nw - glen - plen + 1)
            site = working[i : i+glen-1]
            mismatches = count(k -> g[k] != site[k], 1:glen)
            mismatches > max_mismatches && continue
            # Check PAM
            pam_seg = working[i+glen : min(i+glen+plen-1, nw)]
            pam_ok = length(pam_seg) == plen && pam_matches(system.pam, pam_seg)
            pam_ok || continue

            push!(ots, (
                guide         = g,
                chromosome    = chromosome,
                position      = strand_sign == 1 ? i : n - (i + glen - 1) + 1,
                strand        = Int8(strand_sign),
                sequence      = site,
                mismatches    = mismatches,
                bulges        = 0,
                cfd           = cfd_score(g, site),
                mit           = mit_score(g, [site]),
                pam           = pam_seg))
        end
    end

    sort!(ots, by = r -> r.mismatches)
    result = DataFrame(ots)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "enumerate_off_targets"; parents=provenance_parent_ids(guide, genome_seq), parameters=(max_mismatches=max_mismatches, chromosome=chromosome, system=system.name, row_count=nrow(result)))
end

# Removed enumerate_off_targets BioSequence wrapper functions - main implementation now uses BioSequence

# ---------------------------------------------------------------------------
# Guide design pipeline
# ---------------------------------------------------------------------------

"""
    design_guides(sequence, system; kwargs...) → DataFrame

Design all possible guide RNAs for a target sequence, score them for
on-target efficiency and filter by quality criteria.

Returns a DataFrame of `GuideRNA`-like rows sorted by on-target score.
"""
function design_guides(
    sequence::BioSequence{DNAAlphabet},
    system::CRISPRSystem=SpCas9;
    chromosome::String="chr1",
    min_gc::Real=0.30,
    max_gc::Real=0.80,
    exclude_poly_t::Bool=true,
    exclude_poly_g::Bool=true,
    min_efficiency::Real=0.0,
    top_n::Union{Int,Nothing}=nothing)
    _ctx = active_provenance_context()
    sites = find_pam_sites(sequence, system; chromosome=chromosome)
    if nrow(sites) == 0
        result = DataFrame()
        return _register_crispr_result!(_ctx, result, "design_guides"; parents=provenance_parent_ids(sequence), parameters=(chromosome=chromosome, system=system.name, min_gc=Float64(min_gc), max_gc=Float64(max_gc), min_efficiency=Float64(min_efficiency), guide_count=0))
    end

    results = NamedTuple[]
    for row in eachrow(sites)
        spacer = row.spacer
        gc = guide_gc_content(spacer)
        (gc < Float64(min_gc) || gc > Float64(max_gc)) && continue
        exclude_poly_t && occursin("TTTT", spacer) && continue
        exclude_poly_g && occursin("GGGG", spacer) && continue
        eff = score_on_target(spacer)
        eff < Float64(min_efficiency) && continue

        push!(results, (
            spacer            = spacer,
            chromosome        = row.chrom,
            position          = row.position,
            strand            = row.strand,
            pam               = row.pam,
            on_target_score   = eff,
            gc_content        = gc,
            off_target_count  = -1,   # requires genome-wide search
            specificity_score = missing,
            system            = system.name))
    end

    sort!(results, by = r -> -r.on_target_score)
    df = DataFrame(results)
    top_n !== nothing && nrow(df) > top_n && (df = df[1:top_n, :])
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, df, "design_guides"; parents=provenance_parent_ids(sequence), parameters=(chromosome=chromosome, system=system.name, min_gc=Float64(min_gc), max_gc=Float64(max_gc), min_efficiency=Float64(min_efficiency), guide_count=nrow(df)))
end

# ---------------------------------------------------------------------------
# Guide filtering and ranking
# ---------------------------------------------------------------------------

"""
    filter_guides(guides; min_efficiency=0.4, max_off_targets=10, min_specificity=0.3) → DataFrame

Apply standard quality filters to a guide DataFrame.
"""
function filter_guides(
    guides::DataFrame;
    min_efficiency::Real=0.4,
    max_off_targets::Int=10,
    min_specificity::Real=0.3,
    exclude_restriction_sites::Vector{String}=String[])
    df = copy(guides)
    mask = trues(nrow(df))

    hasproperty(df, :on_target_score)  && (mask .&= df.on_target_score  .>= Float64(min_efficiency))
    if hasproperty(df, :specificity_score)
        mask .&= [ismissing(value) ? true : Float64(value) >= Float64(min_specificity) for value in df.specificity_score]
    end
    if hasproperty(df, :off_target_count)
        valid_ot = (df.off_target_count .== -1) .| (df.off_target_count .<= max_off_targets)
        mask .&= valid_ot
    end

    for rs in exclude_restriction_sites
        hasproperty(df, :spacer) && (mask .&= .!occursin.(rs, df.spacer))
    end

    result = df[mask, :]
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "filter_guides"; parents=provenance_parent_ids(guides), parameters=(min_efficiency=Float64(min_efficiency), max_off_targets=max_off_targets, min_specificity=Float64(min_specificity), row_count=nrow(result)))
end

"""
    rank_guides(guides; weights=(efficiency=0.6, specificity=0.3, gc=0.1)) → DataFrame

Rank guides by a composite score.
"""
function rank_guides(
    guides::DataFrame;
    eff_weight::Real=0.6,
    spec_weight::Real=0.3,
    gc_weight::Real=0.1)
    df = copy(guides)
    n  = nrow(df)
    n == 0 && return df

        eff  = hasproperty(df, :on_target_score)   ? Float64.(coalesce.(df.on_target_score, 0.5))   : fill(0.5, n)
        spec = hasproperty(df, :specificity_score) ? Float64.(coalesce.(df.specificity_score, 0.5)) : fill(0.5, n)
    gc   = hasproperty(df, :gc_content) ?
            1.0 .- abs.(Float64.(coalesce.(df.gc_content, 0.5)) .- 0.5) ./ 0.5 : fill(0.5, n)

    composite = Float64(eff_weight)*eff .+ Float64(spec_weight)*spec .+ Float64(gc_weight)*gc
    df[!, :composite_score] = composite
    sort!(df, :composite_score, rev=true)
    df[!, :rank] = 1:nrow(df)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, df, "rank_guides"; parents=provenance_parent_ids(guides), parameters=(eff_weight=Float64(eff_weight), spec_weight=Float64(spec_weight), gc_weight=Float64(gc_weight), row_count=nrow(df)))
end

"""
    guide_specificity_score(guide, off_targets_df) → Float64

Compute a combined specificity score from enumerated off-targets.
"""
function guide_specificity_score(guide::BioSequence{DNAAlphabet}, off_targets_df::DataFrame)
    nrow(off_targets_df) == 0 && return 1.0
    cfd_vals = hasproperty(off_targets_df, :cfd) ? Float64.(off_targets_df.cfd) : fill(0.5, nrow(off_targets_df))
    # Aggregate: penalise by sum of CFD scores across all off-targets

    return clamp(1.0 - sum(cfd_vals) / max(length(cfd_vals) * 10, 1), 0.0, 1.0)
end

# ---------------------------------------------------------------------------
# Base editor guide design
# ---------------------------------------------------------------------------

"""
    design_base_editor_guides(sequence, editor; system=SpCas9, kwargs...) → DataFrame

Design guides for base editing: identifies spacers where the target base
falls within the editing window, reporting bystander edits.

Returns a DataFrame with `editable_positions` and `bystander_positions`.
"""
function design_base_editor_guides(
    sequence::BioSequence{DNAAlphabet},
    editor::EditingWindow;
    system::CRISPRSystem=SpCas9,
    kwargs...
)
    _ctx = active_provenance_context()
    guides = design_guides(sequence, system; kwargs...)
    nrow(guides) == 0 && return _register_crispr_result!(_ctx, guides, "design_base_editor_guides"; parents=provenance_parent_ids(sequence), parameters=(editor=editor.editor_name, system=system.name, guide_count=0))

    editable_pos     = Vector{Vector{Int}}(undef, nrow(guides))
    bystander_pos    = Vector{Vector{Int}}(undef, nrow(guides))
    has_target_base  = BitVector(undef, nrow(guides))

    for (i, row) in enumerate(eachrow(guides))
        sp   = row.spacer
        glen = length(sp)
        win_s = max(1, editor.window_start)
        win_e = min(glen, editor.window_end)
        target = editor.editable_base

        edit_pos    = [p for p in win_s:win_e if p <= glen && sp[p] == target]
        bystander   = [p for p in editor.bystander_positions if p <= glen && p ∉ editor.window_start:editor.window_end && sp[p] == target]

        editable_pos[i]   = edit_pos
        bystander_pos[i]  = bystander
        has_target_base[i] = !isempty(edit_pos)
    end

    guides[!, :editable_positions]  = editable_pos
    guides[!, :bystander_positions] = bystander_pos
    guides[!, :has_target_base]     = has_target_base
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, guides, "design_base_editor_guides"; parents=provenance_parent_ids(sequence), parameters=(editor=editor.editor_name, system=system.name, guide_count=nrow(guides)))
end

"""
    analyze_editing_window(spacer, editor) → NamedTuple

Report which positions within the editing window contain the target base.
"""
function analyze_editing_window(spacer::BioSequence{DNAAlphabet}, editor::EditingWindow)
    sp   = uppercase(String(spacer))
    glen = length(sp)
    ws   = max(1, editor.window_start)
    we   = min(glen, editor.window_end)
    target = editor.editable_base

    window_seq  = ws <= we <= glen ? sp[ws:we] : ""
    edit_sites  = [p for p in ws:we if p <= glen && sp[p] == target]
    bystanders  = [p for p in editor.bystander_positions if p <= glen && p ∉ ws:we && sp[p] == target]

    result = (
        spacer            = sp,
        editor            = editor.editor_name,
        window_sequence   = window_seq,
        editable_positions = edit_sites,
        bystander_positions = bystanders,
        n_editable        = length(edit_sites),
        n_bystander       = length(bystanders),
        target_base       = editor.editable_base,
        product_base      = editor.edited_base)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "analyze_editing_window"; parents=provenance_parent_ids(spacer), parameters=(editor=editor.editor_name, editable_count=length(edit_sites), bystander_count=length(bystanders)))
end

# ---------------------------------------------------------------------------
# Prime editing (pegRNA design)
# ---------------------------------------------------------------------------

"""
    design_pegrna(target_sequence, edit; nick_position=17, pbs_length=13, rt_template_length=15) → NamedTuple

Design a prime editing guide RNA (pegRNA) for a desired edit.

Returns spacer, PBS (primer binding site), RT template, and scaffold linker.
Based on Anzalone et al. (2019) Nature 576:149-157.
"""
function design_pegrna(
    target_sequence::AbstractString,
    edit::AbstractString;
    nick_position::Int=17,
    pbs_length::Int=13,
    rt_template_length::Int=15,
    system::CRISPRSystem=SpCas9)
    seq  = uppercase(String(target_sequence))
    edit_seq = uppercase(String(edit))
    n    = length(seq)

    # Spacer: protospacer_length bp ending at nick position
    glen = system.protospacer_length
    spacer_end = min(nick_position + glen - 1, n)
    spacer = seq[max(1, spacer_end - glen + 1) : spacer_end]

    # PBS: reverse complement of sequence downstream of nick
    pbs_start = nick_position + 1
    pbs_end   = min(pbs_start + pbs_length - 1, n)
    pbs_raw   = seq[pbs_start : pbs_end]
    pbs        = _reverse_complement_str(pbs_raw)

    # RT template: edit + flanking seq
    rt_start = max(1, nick_position - rt_template_length + 1)
    rt_raw   = seq[rt_start : nick_position]
    # Incorporate edit in template
    rt_template = _reverse_complement_str(rt_raw * edit_seq)

    efficiency = prime_editing_guide_score(spacer, pbs, rt_template)

    result = (
        spacer       = spacer,
        pbs          = pbs,
        rt_template  = rt_template,
        full_pegrna  = spacer * "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC" * rt_template * pbs,
        nick_position = nick_position,
        efficiency_estimate = efficiency)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "design_pegrna"; parents=String[], parameters=(nick_position=nick_position, pbs_length=pbs_length, rt_template_length=rt_template_length, system=system.name, efficiency_estimate=efficiency))
end

"""
    prime_editing_guide_score(spacer, pbs, rt_template) → Float64

Estimate pegRNA efficiency based on:
- On-target score of spacer
- PBS GC content and length (optimal 10-16 nt)
- RT template secondary structure proxy (GC content)
"""
function prime_editing_guide_score(spacer::BioSequence{DNAAlphabet}, pbs::BioSequence{DNAAlphabet}, rt_template::BioSequence{DNAAlphabet})
    spacer_score = score_on_target(spacer)
    pbs_gc       = guide_gc_content(pbs)
    rt_gc        = guide_gc_content(rt_template)
    pbs_len_score = 1.0 - abs(length(pbs) - 13) / 13  # penalty from optimum

    # Simple linear combination
    score = 0.4 * spacer_score + 0.3 * pbs_gc + 0.2 * rt_gc + 0.1 * pbs_len_score

    return clamp(score, 0.0, 1.0)
end

prime_editing_guide_score(spacer::AbstractString, pbs::AbstractString, rt_template::AbstractString) = prime_editing_guide_score(DNASeq(spacer), DNASeq(pbs), DNASeq(rt_template))

# ---------------------------------------------------------------------------
# HDR template design
# ---------------------------------------------------------------------------

"""
    design_hdr_template(sequence, edit, cut_position; homology_arm_length=80) → NamedTuple

Design an HDR (Homology-Directed Repair) donor template for precise genome editing.
Returns left arm, edit insert, right arm, and full template.
"""
function design_hdr_template(
    sequence::AbstractString,
    edit::AbstractString,
    cut_position::Int;
    homology_arm_length::Int=80,
    silent_pam_mutation::Bool=true,
    system::CRISPRSystem=SpCas9)
    seq  = uppercase(String(sequence))
    n    = length(seq)
    cut  = clamp(cut_position, 1, n)

    left_start = max(1, cut - homology_arm_length)
    left_arm   = seq[left_start : cut]
    right_arm  = seq[min(cut+1, n) : min(cut + homology_arm_length, n)]

    # Optionally mutate PAM to prevent re-cutting
    if silent_pam_mutation && system.pam == "NGG"
        right_arm = _mutate_pam_silent(right_arm, system)
    end

    full_template = left_arm * uppercase(String(edit)) * right_arm
    result = (
        left_arm      = left_arm,
        insert        = uppercase(String(edit)),
        right_arm     = right_arm,
        full_template = full_template,
        template_length = length(full_template),
        cut_position  = cut)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "design_hdr_template"; parents=String[], parameters=(cut_position=cut, homology_arm_length=homology_arm_length, silent_pam_mutation=silent_pam_mutation, system=system.name, template_length=length(full_template)))
end

function _mutate_pam_silent(arm::AbstractString, system::CRISPRSystem)
    # Mutate first NGG to NGA/NGC (silent in most codons)
    i = findfirst("GG", arm)
    i === nothing && return arm
    chars = collect(arm)
    chars[last(i)] = 'A'   # GG → GA: reduces re-cutting, often synonymous
    return String(chars)
end

# ---------------------------------------------------------------------------
# Indel prediction (Lindel-like)
# ---------------------------------------------------------------------------

"""
    predict_indels(spacer; n_samples=1000, seed=1) → DataFrame

Predict the indel distribution resulting from NHEJ repair after Cas9 cutting.
Uses a simplified Lindel-inspired model based on sequence context.

Returns a DataFrame of (indel_type, size, frequency).
"""
function predict_indels(
    spacer::BioSequence{DNAAlphabet};
    n_samples::Int=1000,
    seed::Int=1)
    rng = MersenneTwister(seed)
    s   = uppercase(String(spacer))
    n   = length(s)

    # Simplified model: propensity for insertions vs deletions based on
    # microhomology content and GC content
    gc   = guide_gc_content(s)
    # Score microhomology in the cut-proximal region
    cut_region = s[max(1, n-5):n]
    mh_score = _microhomology_score(cut_region)

    # Insertion probability ≈ 0.2–0.5 depending on sequence
    p_ins = clamp(0.2 + 0.3 * (1 - gc), 0.1, 0.5)
    p_del = 1.0 - p_ins

    rows = NamedTuple[]
    # 1bp insertions (most common NHEJ outcome)
    push!(rows, (indel_type="insertion", size=1, frequency=p_ins * 0.7))
    push!(rows, (indel_type="insertion", size=2, frequency=p_ins * 0.2))
    push!(rows, (indel_type="insertion", size=3, frequency=p_ins * 0.1))

    # Deletions — microhomology mediates larger deletions
    del_sizes  = [1, 2, 3, 5, 7, 10, 15, 20]
    del_probs  = [0.3, 0.2, 0.15, 0.12, 0.1, 0.07, 0.03, 0.03] .* (1.0 + mh_score)
    del_probs ./= sum(del_probs)
    for (sz, prob) in zip(del_sizes, del_probs)
        push!(rows, (indel_type="deletion", size=sz, frequency=p_del * prob))
    end

    df = DataFrame(rows)
    df[!, :frequency] ./= sum(df.frequency)     # normalise
    sort!(df, :frequency, rev=true)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, df, "predict_indels"; parents=provenance_parent_ids(spacer), parameters=(n_samples=n_samples, seed=seed, row_count=nrow(df)))
end

function _microhomology_score(seq::AbstractString)
    # Detect simple di/trinucleotide repeats as microhomology proxy
    n = length(seq)
    n < 4 && return 0.0
    score = 0.0
    for k in 2:min(4, n÷2)
        for i in 1:(n-2k+1)
            score += seq[i:i+k-1] == seq[i+k:i+2k-1] ? 0.1 * k : 0.0
        end
    end
    return clamp(score, 0.0, 1.0)
end

"""
    indel_distribution(spacer) → NamedTuple

Return summary statistics of predicted indel outcomes.
"""
function indel_distribution(spacer::BioSequence{DNAAlphabet})
    _ctx = active_provenance_context()
    df = predict_indels(spacer)
    ins_df = filter(r -> r.indel_type == "insertion", df)
    del_df = filter(r -> r.indel_type == "deletion",  df)
    result = (
        insertion_fraction  = isempty(ins_df) ? 0.0 : sum(ins_df.frequency),
        deletion_fraction   = isempty(del_df) ? 0.0 : sum(del_df.frequency),
        frameshift_fraction = sum(r.frequency for r in eachrow(df) if r.size % 3 != 0; init=0.0),
        most_common_indel   = nrow(df) > 0 ? df[1, :indel_type] * string(df[1, :size]) : "none",
        distribution        = df)


    return _register_crispr_result!(_ctx, result, "indel_distribution"; parents=provenance_parent_ids(spacer), parameters=(distribution_size=nrow(df)))
end

# ---------------------------------------------------------------------------
# Library design
# ---------------------------------------------------------------------------

"""
    design_library(gene_sequences, system; guides_per_gene=6, kwargs...) → DataFrame

Design a genome-wide CRISPR guide library targeting multiple genes.
Selects the top `guides_per_gene` guides per gene after scoring and filtering.

Analogous to the output of Brunello / GeckoV2 library design pipelines.
"""
function design_library(
    gene_sequences::Dict{String,String},
    system::CRISPRSystem=SpCas9;
    guides_per_gene::Int=6,
    include_controls::Int=100,
    seed::Int=1,
    kwargs...
)
    rng = MersenneTwister(seed)
    all_guides = DataFrame[]

    for (gene_name, seq) in gene_sequences
        guides = design_guides(DNASeq(seq), system; kwargs...)
        isempty(guides) && continue
        guides = rank_guides(guides)
        n_sel = min(guides_per_gene, nrow(guides))
        sel = guides[1:n_sel, :]
        sel[!, :gene] .= gene_name
        push!(all_guides, sel)
    end

    # Add non-targeting controls
    if include_controls > 0
        ctrl_seqs = [randstring(rng, "ACGT", system.protospacer_length) for _ in 1:include_controls]
        ctrl_df = DataFrame(
            spacer  = ctrl_seqs,
            gene    = fill("non_targeting_control", include_controls),
            on_target_score = fill(0.0, include_controls),
            gc_content = guide_gc_content.(ctrl_seqs),
            composite_score = fill(0.0, include_controls))
        push!(all_guides, ctrl_df)
    end

    isempty(all_guides) && return DataFrame()
    lib = vcat(all_guides...; cols=:union)
    lib[!, :library_index] = 1:nrow(lib)
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, lib, "design_library"; parents=String[], parameters=(gene_count=length(gene_sequences), guides_per_gene=guides_per_gene, include_controls=include_controls, row_count=nrow(lib)))
end

"""
    library_coverage_stats(library, n_genes) → NamedTuple

Report coverage statistics for a designed library.
"""
function library_coverage_stats(library::DataFrame, n_genes::Int)
    total = nrow(library)
    n_ctrl = hasproperty(library, :gene) ?
             count(g -> g == "non_targeting_control", library.gene) : 0
    n_targeting = total - n_ctrl
    targeting_genes = hasproperty(library, :gene) ? length(unique(filter(g -> g != "non_targeting_control", library.gene))) : 0
    result = (
        total_guides       = total,
        targeting_guides   = n_targeting,
        control_guides     = n_ctrl,
        genes_covered      = targeting_genes,
        gene_coverage_frac = targeting_genes / max(n_genes, 1),
        mean_guides_per_gene = n_targeting / max(targeting_genes, 1))
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "library_coverage_stats"; parents=provenance_parent_ids(library), parameters=(n_genes=n_genes, total_guides=total, genes_covered=targeting_genes))
end

# ---------------------------------------------------------------------------
# CRISPR screen analysis (MAGeCK-like)
# ---------------------------------------------------------------------------

"""
    crispr_screen_analysis(counts_treatment, counts_control, guide_gene_map; method=:rra) → DataFrame

Analyze a CRISPR screen for essentiality or enrichment using a simplified
Robust Rank Aggregation (RRA) approach, analogous to MAGeCK.

`counts_treatment`, `counts_control`: matrices (n_guides × n_replicates).
`guide_gene_map`: DataFrame with `guide` and `gene` columns.

Returns a gene-level result DataFrame sorted by score.

Reference: Li et al. (2014) Genome Biology 15:554.
"""
function crispr_screen_analysis(
    counts_treatment::AbstractMatrix{<:Real},
    counts_control::AbstractMatrix{<:Real},
    guide_gene_map::DataFrame;
    method::Symbol=:rra,
    min_reads::Int=10,
    pseudocount::Real=1.0)
    size(counts_treatment) == size(counts_control) ||
        throw(DimensionMismatch("treatment and control must have same dimensions"))
    n_guides, n_reps = size(counts_treatment)
    nrow(guide_gene_map) == n_guides ||
        throw(DimensionMismatch("guide_gene_map must have one row per guide"))

    pc = Float64(pseudocount)
    # Normalise each sample to total reads × 1e6 (RPM)
    treat_rpm = (Float64.(counts_treatment) .+ pc) ./
                max.(vec(sum(counts_treatment, dims=1))', 1) .* 1e6
    ctrl_rpm  = (Float64.(counts_control)   .+ pc) ./
                max.(vec(sum(counts_control,  dims=1))', 1) .* 1e6

    # Log2 fold change per guide (mean across replicates)
    lfc = vec(mean(log2.(treat_rpm ./ ctrl_rpm), dims=2))

    # Per-guide p-value: simple t-test across replicates
    guide_pvals = zeros(Float64, n_guides)
    for i in 1:n_guides
        t = Float64.(counts_treatment[i,:]) .+ pc
        c = Float64.(counts_control[i,:])   .+ pc
        n_t, n_c = length(t), length(c)
        if n_t >= 2 && n_c >= 2
            d = log2.(t) .- log2.(c)
            guide_pvals[i] = _one_sample_ttest_pvalue(d)
        else
            guide_pvals[i] = 1.0
        end
    end

    # Gene-level aggregation (RRA-like: rank-based)
    genes = unique(guide_gene_map.gene)
    gene_results = NamedTuple[]

    for gene in genes
        idx = findall(g -> g == gene, guide_gene_map.gene)
        g_lfc    = lfc[idx]
        g_pvals  = guide_pvals[idx]
        # Pool p-values: geometric mean (simpler than RRA; same direction)
        pooled_p = prod(g_pvals) ^ (1 / max(length(g_pvals), 1))
        push!(gene_results, (
            gene         = gene,
            n_guides     = length(idx),
            mean_lfc     = mean(g_lfc),
            median_lfc   = median(g_lfc),
            gene_pvalue  = pooled_p,
            direction    = mean(g_lfc) > 0 ? "enriched" : "depleted"))
    end

    sort!(gene_results, by = r -> r.gene_pvalue)

    # BH FDR
    n = length(gene_results)
    padj = [min(gene_results[i].gene_pvalue * n / i, 1.0) for i in 1:n]
    for i in (n-1):-1:1; padj[i] = min(padj[i], padj[i+1]); end

    df = DataFrame(gene_results)
    df[!, :padj] = padj
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, df, "crispr_screen_analysis"; parents=provenance_parent_ids(counts_treatment, counts_control, guide_gene_map), parameters=(method=String(method), min_reads=min_reads, pseudocount=Float64(pseudocount), row_count=nrow(df)))
end

"""
    mageck_like_test(treatment_counts, control_counts, guide_df) → DataFrame

Convenience alias for `crispr_screen_analysis` with MAGeCK-style column names.
"""
mageck_like_test(t, c, g; kwargs...) = crispr_screen_analysis(t, c, g; kwargs...)

"""
    genome_wide_offtarget_summary(guide, ot_df) → NamedTuple

Summary statistics for genome-wide off-target predictions.
"""
function genome_wide_offtarget_summary(guide::BioSequence{DNAAlphabet}, ot_df::DataFrame)
    isempty(ot_df) || nrow(ot_df) == 0 && return (
        guide=guide, n_off_targets=0, n_0mm=0, n_1mm=0, n_2mm=0, n_3mm=0,
        mean_cfd=0.0, specificity_score=1.0
    )
    hasproperty(ot_df, :mismatches) || return (
        guide=guide, n_off_targets=nrow(ot_df), n_0mm=0, n_1mm=0, n_2mm=0, n_3mm=0,
        mean_cfd=0.0, specificity_score=1.0
    )
    mm = Int.(ot_df.mismatches)
    cfd_vals = hasproperty(ot_df, :cfd) ? Float64.(ot_df.cfd) : fill(0.0, nrow(ot_df))
    result = (
        guide             = String(guide),
        n_off_targets     = nrow(ot_df),
        n_0mm             = count(==(0), mm),
        n_1mm             = count(==(1), mm),
        n_2mm             = count(==(2), mm),
        n_3mm             = count(==(3), mm),
        mean_cfd          = mean(cfd_vals),
        specificity_score = guide_specificity_score(guide, ot_df))
    _ctx = active_provenance_context()


    return _register_crispr_result!(_ctx, result, "genome_wide_offtarget_summary"; parents=provenance_parent_ids(guide, ot_df), parameters=(off_target_count=nrow(ot_df), mean_cfd=mean(cfd_vals)))
end

# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

function _one_sample_ttest_pvalue(x::AbstractVector{<:Real})
    n = length(x)
    n < 2 && return 1.0
    m, s = mean(x), std(x)
    s <= 0 && return (m != 0 ? 0.0 : 1.0)
    t = m / (s / sqrt(n))
    df_v = n - 1
    # Normal approximation for df >= 10
    z = abs(t)
    return 2 * (1 - _ncdf(z))
end

function _ncdf(z::Float64)
    0.5 * (1 + erf(z / sqrt(2.0)))
end

end  # module CRISPR
