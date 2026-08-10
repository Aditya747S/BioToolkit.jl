# ==============================================================================
# crispr.jl — CRISPR-Cas guide design & off-target analysis
#
# Fully self-contained module for:
#   - Guide RNA design for SpCas9, SaCas9, Cas12a, CasX, SpRY, Cas13a, Cas13b
#   - On-target efficiency scoring (Rule Set 2 / Doench 2016 & DeepCpf1)
#   - Off-target enumeration and scoring (CFD, MIT score)
#   - PAM detection (NGG, TTTV, TTTN, NNGRRT, etc.)
#   - Genome-wide guide library design
#   - Base editor window analysis (BE3, ABE8e, CBE4max)
#   - Prime editing guide design (pegRNA, Anzalone et al. 2019)
#   - CRISPR screen MAGeCK-like analysis (NB-Wald + RRA ranking)
#   - HDR template design with silent PAM mutation
#   - Indel prediction (Lindel-like microhomology model)
#   - Interactive HTML visualizations for screens, guide libraries, and editing windows
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
using DataAPI
using Statistics
using LinearAlgebra
using Random
using Distributions

# Import biotypes for type-safe sequence handling
using ..BioToolkit: AASeq, AminoAcidAlphabet, BioSequence, DNAAlphabet, DNASeq, SummarizedExperiment, assay, colData
using ..BioToolkit: AbstractAnalysisResult, ResultProvenance, provenance_record
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!
import ..BioToolkit: to_html, export_html

export GuideRNA, CRISPRSystem, OffTarget, EditingWindow
export design_guides, score_on_target, find_pam_sites, enumerate_off_targets
export FMIndex, OffTargetIndex, OffTargetSearchDiagnostics, build_offtarget_index
export cfd_score, mit_score, guide_gc_content
export design_base_editor_guides, analyze_editing_window
export design_pegrna, prime_editing_guide_score
export crispr_screen_analysis, mageck_like_test, crispr_screen_nb
export MageckRRAResult, CRISPRScreenResult
export design_library, library_coverage_stats
export design_hdr_template
export predict_indels, indel_distribution
export filter_guides, rank_guides
export guide_specificity_score, genome_wide_offtarget_summary

export visualize_crispr_screen_html
export visualize_guide_library_html
export visualize_off_targets_html
export visualize_editing_window_html
export to_html, export_html

@inline function _register_crispr_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end

# ---------------------------------------------------------------------------
# Types & Predefined Systems
# ---------------------------------------------------------------------------

"""
    CRISPRSystem

Specification of a CRISPR-Cas system including PAM, protospacer length, cut position, and strand specificity.
"""
struct CRISPRSystem
    name::String
    pam::String             # PAM sequence (e.g. "NGG", "TTTV", "H")
    pam_location::Symbol    # :three_prime or :five_prime
    protospacer_length::Int
    cut_position::Int       # bp upstream of PAM (for SpCas9 = 3)
    strand_specific::Bool
end

# Predefined CRISPR systems
const SpCas9  = CRISPRSystem("SpCas9",  "NGG",    :three_prime, 20, 3, false)
const SaCas9  = CRISPRSystem("SaCas9",  "NNGRRT", :three_prime, 21, 3, false)
const Cas12a  = CRISPRSystem("Cas12a",  "TTTV",   :five_prime,  23, 18, false)
const CasX    = CRISPRSystem("CasX",    "TTCN",   :five_prime,  20, 10, false)
const SpRY    = CRISPRSystem("SpRY",    "NRN",    :three_prime, 20, 3, false)  # near-PAMless
const Cas13a  = CRISPRSystem("Cas13a",  "H",      :three_prime, 28, 0, true)   # RNA-targeting, 3' PFS non-G
const Cas13b  = CRISPRSystem("Cas13b",  "D",      :five_prime,  30, 0, true)   # RNA-targeting, 5' PFS non-C

export SpCas9, SaCas9, Cas12a, CasX, SpRY, Cas13a, Cas13b

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
    on_target_score::Float64    # Rule Set 2 / Doench / DeepCpf1 score
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

"""
    FMIndex

Checkpointed FM-index with sampled suffix-array locations for bounded Hamming search.
"""
struct FMIndex
    bwt::Vector{UInt8}
    c_table::Vector{Int}
    symbol_to_slot::Vector{Int16}
    checkpoint_counts::Matrix{Int}
    checkpoint_interval::Int
    sample_rows::Vector{Int}
    sample_values::Vector{Int}
    text_length::Int
end

"""
    OffTargetIndex

A contig-aware FM-index for repeated CRISPR off-target searches.
"""
struct OffTargetIndex <: AbstractAnalysisResult
    fm_index::FMIndex
    contig_names::Vector{String}
    contig_sequences::Vector{String}
    contig_starts::Vector{Int}
    contig_stops::Vector{Int}
    suffix_array_sample_rate::Int
    provenance::ResultProvenance
end

struct OffTargetSearchDiagnostics <: AbstractAnalysisResult
    candidate_intervals::Int
    fm_candidates::Int
    verified_sites::Int
    max_mismatches::Int
    method_level::Symbol
    provenance::ResultProvenance
end

# Base editor definitions
const BE3     = EditingWindow("BE3",     4, 8, 'C', 'T', [3,4,5,6,7,8,9])
const ABE8e   = EditingWindow("ABE8e",   4, 8, 'A', 'G', [3,4,5,6,7,8,9])
const CBE4max = EditingWindow("CBE4max", 3, 9, 'C', 'T', [2,3,4,5,6,7,8,9,10])

export BE3, ABE8e, CBE4max

# ---------------------------------------------------------------------------
# IUPAC & Sequence Utilities
# ---------------------------------------------------------------------------

const _IUPAC = Dict(
    'N'=>"ACGT", 'R'=>"AG", 'Y'=>"CT", 'S'=>"GC", 'W'=>"AT",
    'K'=>"GT",   'M'=>"AC", 'B'=>"CGT",'D'=>"AGT",'H'=>"ACT", 'V'=>"ACG",
    'A'=>"A",    'C'=>"C",  'G'=>"G",  'T'=>"T"
)

function _reverse_complement_str(s::AbstractString)
    comp = Dict(
        'A'=>'T', 'T'=>'A', 'G'=>'C', 'C'=>'G', 'N'=>'N',
        'R'=>'Y', 'Y'=>'R', 'S'=>'S', 'W'=>'W', 'K'=>'M', 'M'=>'K',
        'B'=>'V', 'V'=>'B', 'D'=>'H', 'H'=>'D'
    )
    return String(reverse([get(comp, uppercase(c), 'N') for c in s]))
end

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

# ---------------------------------------------------------------------------
# PAM site discovery
# ---------------------------------------------------------------------------

"""
    find_pam_sites(sequence, system; chromosome="chr1") → DataFrame

Find all PAM sites for a CRISPR system in both strands of a sequence.
Accepts `AbstractString` or `BioSequence{DNAAlphabet}`.
Returns a DataFrame with protospacer + PAM positions.
"""
function find_pam_sites(
    sequence::BioSequence{DNAAlphabet},
    system::CRISPRSystem;
    chromosome::String="chr1",
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

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
    else  # five_prime PAM (Cas12a, CasX)
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
                pos = rc_offset - (i + glen - 1) + 1  # 5' start coordinate on forward strand
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
    return _register_crispr_result!(_ctx, result, "find_pam_sites"; parents=provenance_parent_ids(sequence), parameters=(chromosome=chromosome, system=system.name, row_count=nrow(result)))
end

find_pam_sites(seq::AbstractString, system::CRISPRSystem; kwargs...) = find_pam_sites(DNASeq(seq; validate=false), system; kwargs...)

# ---------------------------------------------------------------------------
# GC Content
# ---------------------------------------------------------------------------

"""
    guide_gc_content(sequence) → Float64

Compute GC content fraction [0, 1] for a spacer sequence.
Accepts `AbstractString` or `BioSequence{DNAAlphabet}`.
"""
function guide_gc_content(seq::BioSequence{DNAAlphabet})
    s = String(seq)
    gc = count(c -> c in ('G','C','g','c'), s)
    return gc / max(length(s), 1)
end

guide_gc_content(seq::AbstractString) = guide_gc_content(DNASeq(seq; validate=false))

# ---------------------------------------------------------------------------
# On-target efficiency scoring (Doench 2016 Rule Set 2 & DeepCpf1)
# ---------------------------------------------------------------------------

const _RS2_SINGLE_NUC_WEIGHTS = Dict(
    (1,'G')=>-0.2753771,(2,'A')=>-0.3238875,(2,'C')=>0.17212887,(3,'C')=>-0.1006662,
    (4,'C')=>-0.2018029,(4,'T')=>-0.1747400,(5,'A')=>0.20932776,(5,'C')=>-0.17166690,
    (6,'C')=>0.11385978,(7,'C')=>-0.0671806,(8,'A')=>0.0694723,(8,'T')=>0.18765571,
    (9,'C')=>0.07781066,(9,'G')=>-0.4453400,(10,'C')=>0.27529108,(10,'G')=>-0.22163640,
    (11,'A')=>0.2093,(11,'G')=>-0.09849019,(12,'C')=>-0.29010478,(12,'T')=>0.1564111,
    (13,'G')=>0.07606142,(13,'T')=>-0.2130060736,(14,'C')=>0.1228724,(14,'T')=>-0.10466540,
    (15,'G')=>0.06421542,(15,'T')=>0.0855514,(16,'G')=>0.0498791,(16,'T')=>-0.05312809,
    (17,'C')=>-0.13640294,(17,'G')=>0.1379505,(18,'A')=>0.16827566,(18,'C')=>-0.09963975,
    (19,'A')=>0.28722890,(19,'G')=>-0.21012358,(20,'A')=>-0.06779897,(20,'G')=>0.11098105
)

"""
    score_on_target(spacer; target_context=nothing, system=SpCas9) → Float64

Estimate on-target cutting efficiency (0–1) using Doench 2016 Rule Set 2 (for SpCas9/SaCas9) or DeepCpf1-like scoring (for Cas12a).
Supports 20bp spacer alone or 30bp target context (4bp 5' flank + 20bp spacer + 3bp PAM + 3bp 3' flank).
"""
function score_on_target(
    spacer::BioSequence{DNAAlphabet};
    target_context::Union{Nothing,BioSequence{DNAAlphabet},AbstractString}=nothing,
    system::CRISPRSystem=SpCas9,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    s = uppercase(String(spacer))
    ctx_str = target_context !== nothing ? uppercase(String(target_context)) : ""

    if system.name == "Cas12a"
        # DeepCpf1-like Cas12a scoring heuristic
        score = 0.50
        gc = guide_gc_content(s)
        (gc >= 0.35 && gc <= 0.60) ? (score += 0.15) : (score -= 0.10)
        occursin("TTTT", s) && (score -= 0.25)
        # PAM distal / seed region preferences (TTTV PAM, positions 1-6 seed)
        length(s) >= 6 && occursin("T", s[1:3]) && (score += 0.08)
        result = clamp(score, 0.0, 1.0)
        return _register_crispr_result!(_ctx, result, "score_on_target"; parents=provenance_parent_ids(spacer), parameters=(system=system.name, gc_content=gc, score=result))
    end

    # SpCas9 / SaCas9 Rule Set 2 scoring
    intercept = 0.5977
    score = intercept

    if length(ctx_str) >= 30
        # Full 30-mer context scoring
        for i in 1:20
            nt = ctx_str[i + 4] # spacer starts at position 5 in 30-mer
            score += get(_RS2_SINGLE_NUC_WEIGHTS, (i, nt), 0.0)
        end
    else
        # 20-mer spacer scoring
        for i in 1:min(20, length(s))
            nt = s[i]
            score += get(_RS2_SINGLE_NUC_WEIGHTS, (i, nt), 0.0)
        end
    end

    # GC content penalty
    gc = guide_gc_content(s)
    if gc < 0.25 || gc > 0.75
        score -= 0.15
    end

    # Poly-T penalty (reads off Pol III transcript termination)
    occursin("TTTT", s) && (score -= 0.20)

    # First position G preference
    !isempty(s) && s[1] == 'G' && (score += 0.05)

    result = clamp(1.0 / (1.0 + exp(-score)), 0.0, 1.0)
    return _register_crispr_result!(_ctx, result, "score_on_target"; parents=provenance_parent_ids(spacer), parameters=(system=system.name, gc_content=gc, score=result))
end

score_on_target(spacer::AbstractString; kwargs...) = score_on_target(DNASeq(spacer; validate=false); kwargs...)

# ---------------------------------------------------------------------------
# Off-target scoring: CFD (Cutting Frequency Determination)
# ---------------------------------------------------------------------------

const _CFD_MISMATCH_WEIGHTS = [
    0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0,
    0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828,
    0.615, 0.804, 0.685, 0.583
]  # positions 1-20 from PAM-distal end

"""
    cfd_score(guide, off_target) → Float64

Compute the Cutting Frequency Determination (CFD) score between a guide and an off-target site.
Returns a value in [0, 1] (1 = perfect match).
Reference: Doench et al. (2014) Nature Biotechnology 32:1262-1267.
"""
function cfd_score(guide::BioSequence{DNAAlphabet}, off_target::BioSequence{DNAAlphabet})
    g = uppercase(String(guide))
    o = uppercase(String(off_target))
    length(g) == length(o) || return 0.0

    score = 1.0
    for (i, (gc, oc)) in enumerate(zip(g, o))
        gc == oc && continue
        pos = min(i, length(_CFD_MISMATCH_WEIGHTS))
        w = _CFD_MISMATCH_WEIGHTS[pos]
        score *= max(1.0 - w, 0.0)
    end
    result = clamp(score, 0.0, 1.0)
    _ctx = active_provenance_context()
    return _register_crispr_result!(_ctx, result, "cfd_score"; parents=provenance_parent_ids(guide, off_target), parameters=(length=length(g), score=result))
end

cfd_score(guide::AbstractString, off_target::AbstractString) = cfd_score(DNASeq(guide; validate=false), DNASeq(off_target; validate=false))

# ---------------------------------------------------------------------------
# Off-target scoring: MIT score (Hsu 2013)
# ---------------------------------------------------------------------------

const _MIT_WEIGHTS = [
    0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0,
    0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828,
    0.615, 0.804, 0.685, 0.583
]

"""
    mit_score(guide, off_targets; ignore_exact_self=false) → Float64

Compute the MIT specificity score for a guide against a list of off-target sequences.
Returns a value in [0, 100].
Reference: Hsu et al. (2013) Nature Biotechnology 31:827-832.
"""
function mit_score(guide::BioSequence{DNAAlphabet}, off_targets::AbstractVector{<:BioSequence{DNAAlphabet}}; ignore_exact_self::Bool=false)
    g = uppercase(String(guide))
    isempty(off_targets) && return 100.0

    total_score = 0.0
    for ot in off_targets
        o = uppercase(String(ot))
        length(g) == length(o) || continue
        mm_positions = findall(i -> g[i] != o[i], 1:length(g))
        n_mm = length(mm_positions)

        if n_mm == 0
            ignore_exact_self && continue
            total_score += 1.0
            continue
        end

        s = prod(1.0 - _MIT_WEIGHTS[min(p, 20)] for p in mm_positions; init=1.0)
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

mit_score(guide::AbstractString, off_targets::AbstractVector{<:AbstractString}; kwargs...) = mit_score(DNASeq(guide; validate=false), DNASeq.(off_targets; validate=false); kwargs...)
mit_score(guide::Union{AbstractString,BioSequence{DNAAlphabet}}, off_target::Union{AbstractString,BioSequence{DNAAlphabet}}; kwargs...) = mit_score(guide isa AbstractString ? DNASeq(guide; validate=false) : guide, [off_target isa AbstractString ? DNASeq(off_target; validate=false) : off_target]; kwargs...)

# ---------------------------------------------------------------------------
# Off-target indexing and FM-index enumeration
# ---------------------------------------------------------------------------

const _FM_SENTINEL = UInt8(0x00)
const _FM_CONTIG_SEPARATOR = UInt8(0x01)
const _OFFTARGET_SEARCH_ALPHABET = UInt8[UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T')]

function _suffix_array_prefix_doubling(text::Vector{UInt8})
    n = length(text)
    n > 0 || throw(ArgumentError("cannot build an FM-index over an empty reference"))
    suffixes = collect(1:n)
    ranks = Int.(text)
    next_ranks = similar(ranks)
    width = 1
    while true
        sort!(suffixes; by=index -> (ranks[index], index + width <= n ? ranks[index + width] : -1))
        classes = 1
        next_ranks[suffixes[1]] = classes
        for order_index in 2:n
            current = suffixes[order_index]
            previous = suffixes[order_index - 1]
            current_key = (ranks[current], current + width <= n ? ranks[current + width] : -1)
            previous_key = (ranks[previous], previous + width <= n ? ranks[previous + width] : -1)
            current_key != previous_key && (classes += 1)
            next_ranks[current] = classes
        end
        ranks, next_ranks = next_ranks, ranks
        classes == n && return suffixes
        width >= n && return suffixes
        width *= 2
    end
end

function _build_fm_index(text::Vector{UInt8}; checkpoint_interval::Integer=128, suffix_array_sample_rate::Integer=32)
    checkpoint_interval > 0 || throw(ArgumentError("checkpoint_interval must be positive"))
    suffix_array_sample_rate > 0 || throw(ArgumentError("suffix_array_sample_rate must be positive"))
    suffix_array = _suffix_array_prefix_doubling(text)
    n = length(text)
    bwt = Vector{UInt8}(undef, n)
    for row in eachindex(suffix_array)
        suffix_start = suffix_array[row]
        bwt[row] = text[suffix_start == 1 ? n : suffix_start - 1]
    end

    counts = zeros(Int, 256)
    for symbol in text
        counts[Int(symbol) + 1] += 1
    end
    c_table = zeros(Int, 256)
    cumulative = 0
    for index in eachindex(c_table)
        c_table[index] = cumulative
        cumulative += counts[index]
    end
    alphabet = UInt8[index - 1 for index in eachindex(counts) if counts[index] > 0]
    symbol_to_slot = zeros(Int16, 256)
    for (slot, symbol) in enumerate(alphabet)
        symbol_to_slot[Int(symbol) + 1] = Int16(slot)
    end

    block_count = cld(n, Int(checkpoint_interval))
    checkpoints = zeros(Int, length(alphabet), block_count + 1)
    running = zeros(Int, length(alphabet))
    for position in eachindex(bwt)
        slot = Int(symbol_to_slot[Int(bwt[position]) + 1])
        running[slot] += 1
        if position % checkpoint_interval == 0
            checkpoints[:, position ÷ checkpoint_interval + 1] .= running
        end
    end

    sample_rows = Int[]
    sample_values = Int[]
    for (row, suffix_start) in enumerate(suffix_array)
        if (suffix_start - 1) % suffix_array_sample_rate == 0
            push!(sample_rows, row)
            push!(sample_values, suffix_start)
        end
    end
    return FMIndex(bwt, c_table, symbol_to_slot, checkpoints, Int(checkpoint_interval), sample_rows, sample_values, n)
end

@inline function _fm_occurrence(index::FMIndex, symbol::UInt8, position::Int)
    position <= 0 && return 0
    position <= index.text_length || throw(BoundsError(index.bwt, position))
    slot = Int(index.symbol_to_slot[Int(symbol) + 1])
    slot == 0 && return 0
    complete_block = position ÷ index.checkpoint_interval
    count = index.checkpoint_counts[slot, complete_block + 1]
    block_start = complete_block * index.checkpoint_interval
    for offset in block_start + 1:position
        index.bwt[offset] == symbol && (count += 1)
    end
    return count
end

@inline function _fm_extend(index::FMIndex, left::Int, right::Int, symbol::UInt8)
    slot = Int(index.symbol_to_slot[Int(symbol) + 1])
    slot == 0 && return 1, 0
    next_left = index.c_table[Int(symbol) + 1] + _fm_occurrence(index, symbol, left - 1) + 1
    next_right = index.c_table[Int(symbol) + 1] + _fm_occurrence(index, symbol, right)
    return next_left, next_right
end

@inline function _fm_lf(index::FMIndex, row::Int)
    symbol = index.bwt[row]
    return index.c_table[Int(symbol) + 1] + _fm_occurrence(index, symbol, row)
end

function _fm_locate(index::FMIndex, row::Int)
    steps = 0
    current_row = row
    while true
        sample_index = searchsortedfirst(index.sample_rows, current_row)
        if sample_index <= length(index.sample_rows) && index.sample_rows[sample_index] == current_row
            return mod(index.sample_values[sample_index] + steps - 1, index.text_length) + 1
        end
        current_row = _fm_lf(index, current_row)
        steps += 1
        steps <= index.text_length || throw(ErrorException("FM-index locate failed to reach a suffix-array sample"))
    end
end

function _fm_hamming_intervals(index::FMIndex, pattern::Vector{UInt8}, max_mismatches::Int)
    stack = Tuple{Int,Int,Int,Int}[(length(pattern), 1, index.text_length, 0)]
    intervals = Tuple{Int,Int,Int}[]
    expanded = 0
    while !isempty(stack)
        pattern_position, left, right, mismatches = pop!(stack)
        pattern_position == 0 && (push!(intervals, (left, right, mismatches)); continue)
        requested = pattern[pattern_position]
        for observed in _OFFTARGET_SEARCH_ALPHABET
            next_mismatches = mismatches + (observed == requested ? 0 : 1)
            next_mismatches <= max_mismatches || continue
            next_left, next_right = _fm_extend(index, left, right, observed)
            next_left <= next_right || continue
            expanded += 1
            push!(stack, (pattern_position - 1, next_left, next_right, next_mismatches))
        end
    end
    return intervals, expanded
end

function _offtarget_index_from_contigs(
    contig_names::Vector{String},
    contig_sequences::Vector{String};
    checkpoint_interval::Integer=128,
    suffix_array_sample_rate::Integer=32,
    max_build_bases::Union{Nothing,Integer}=nothing,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    length(contig_names) == length(contig_sequences) || throw(DimensionMismatch("contig names and sequences must have identical lengths"))
    isempty(contig_sequences) && throw(ArgumentError("at least one contig is required"))
    all(!isempty, contig_sequences) || throw(ArgumentError("empty contigs are not supported by the FM-index builder"))
    total_bases = sum(length, contig_sequences)
    max_build_bases === nothing || total_bases <= max_build_bases || throw(ArgumentError("reference has $(total_bases) bases, exceeding max_build_bases=$(max_build_bases)"))

    text = UInt8[]
    starts = Int[]
    stops = Int[]
    for (contig_name, sequence) in zip(contig_names, contig_sequences)
        isempty(contig_name) && throw(ArgumentError("contig names must not be empty"))
        push!(starts, length(text) + 1)
        append!(text, codeunits(sequence))
        push!(stops, length(text))
        push!(text, _FM_CONTIG_SEPARATOR)
    end
    text[end] = _FM_SENTINEL
    fm_index = _build_fm_index(text; checkpoint_interval=checkpoint_interval, suffix_array_sample_rate=suffix_array_sample_rate)
    provenance = provenance_record("OffTargetIndex", "CRISPR/build_offtarget_index"; parameters=(contig_count=length(contig_names), total_bases=total_bases, checkpoint_interval=Int(checkpoint_interval), suffix_array_sample_rate=Int(suffix_array_sample_rate)))
    result = OffTargetIndex(fm_index, contig_names, contig_sequences, starts, stops, Int(suffix_array_sample_rate), provenance)
    return _register_crispr_result!(_ctx, result, "build_offtarget_index"; parents=String[], parameters=(contig_count=length(contig_names), total_bases=total_bases, checkpoint_interval=Int(checkpoint_interval), suffix_array_sample_rate=Int(suffix_array_sample_rate)))
end

"""
    build_offtarget_index(genome; chromosome="chr1", ...)

Build a reusable FM-index for a DNA reference.
"""
function build_offtarget_index(genome::BioSequence{DNAAlphabet}; chromosome::AbstractString="chr1", kwargs...)
    return _offtarget_index_from_contigs([String(chromosome)], [uppercase(String(genome))]; kwargs...)
end

build_offtarget_index(genome::AbstractString; kwargs...) = build_offtarget_index(DNASeq(genome; validate=false); kwargs...)

function build_offtarget_index(genomes::AbstractDict; kwargs...)
    keys_sorted = sort!(collect(keys(genomes)); by=string)
    names = String.(keys_sorted)
    sequences = [uppercase(String(genomes[key])) for key in keys_sorted]
    return _offtarget_index_from_contigs(names, sequences; kwargs...)
end

@inline function _index_contig_position(index::OffTargetIndex, text_position::Int)
    contig_index = searchsortedlast(index.contig_starts, text_position)
    return contig_index >= 1 && text_position <= index.contig_stops[contig_index] ? contig_index : nothing
end

function _fm_candidate_positions(index::OffTargetIndex, pattern::String, max_mismatches::Int; max_candidates::Int)
    intervals, expanded = _fm_hamming_intervals(index.fm_index, collect(codeunits(pattern)), max_mismatches)
    positions = Int[]
    for (left, right, _) in intervals
        length(positions) + (right - left + 1) <= max_candidates || throw(ArgumentError("FM-index search exceeded max_candidates=$(max_candidates)"))
        for row in left:right
            push!(positions, _fm_locate(index.fm_index, row))
        end
    end
    return positions, length(intervals), expanded
end

function _append_verified_offtargets!(
    rows::Vector{NamedTuple},
    index::OffTargetIndex,
    guide::String,
    candidate_positions::Vector{Int},
    strand_sign::Int8,
    system::CRISPRSystem)

    guide_length = length(guide)
    pam_length = length(system.pam)
    for text_position in candidate_positions
        contig_index = _index_contig_position(index, text_position)
        contig_index === nothing && continue
        sequence = index.contig_sequences[contig_index]
        local_start = text_position - index.contig_starts[contig_index] + 1
        if strand_sign == 1
            local_start + guide_length + pam_length - 1 <= length(sequence) || continue
            site = sequence[local_start:local_start + guide_length - 1]
            pam = sequence[local_start + guide_length:local_start + guide_length + pam_length - 1]
            pam_matches(system.pam, pam) || continue
            position = local_start
        else
            local_start - pam_length >= 1 || continue
            forward_site = sequence[local_start:local_start + guide_length - 1]
            site = _reverse_complement_str(forward_site)
            pam = _reverse_complement_str(sequence[local_start - pam_length:local_start - 1])
            pam_matches(system.pam, pam) || continue
            position = local_start
        end
        mismatches = count(index -> guide[index] != site[index], eachindex(guide))
        push!(rows, (
            guide=guide,
            chromosome=index.contig_names[contig_index],
            position=position,
            strand=strand_sign,
            sequence=site,
            mismatches=mismatches,
            bulges=0,
            cfd=cfd_score(DNASeq(guide; validate=false), DNASeq(site; validate=false)),
            mit=mit_score(DNASeq(guide; validate=false), [DNASeq(site; validate=false)]),
            pam=pam))
    end
    return rows
end

"""
    enumerate_off_targets(guide, index; max_mismatches=3, max_candidates=100_000, system=SpCas9)

Enumerate PAM-valid off-targets through FM-index backward search with bounded Hamming distance.
"""
function enumerate_off_targets(
    guide::BioSequence{DNAAlphabet},
    index::OffTargetIndex;
    max_mismatches::Integer=3,
    max_candidates::Integer=100_000,
    system::CRISPRSystem=SpCas9,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    0 <= max_mismatches <= length(guide) || throw(ArgumentError("max_mismatches must lie between 0 and guide length"))
    max_candidates > 0 || throw(ArgumentError("max_candidates must be positive"))
    guide_string = uppercase(String(guide))

    forward_positions, forward_intervals, _ = _fm_candidate_positions(index, guide_string, Int(max_mismatches); max_candidates=Int(max_candidates))
    reverse_guide = _reverse_complement_str(guide_string)
    remaining = Int(max_candidates) - length(forward_positions)
    remaining > 0 || throw(ArgumentError("FM-index search reached max_candidates before reverse-strand evaluation"))
    reverse_positions, reverse_intervals, _ = _fm_candidate_positions(index, reverse_guide, Int(max_mismatches); max_candidates=remaining)
    rows = NamedTuple[]
    _append_verified_offtargets!(rows, index, guide_string, forward_positions, Int8(1), system)
    _append_verified_offtargets!(rows, index, guide_string, reverse_positions, Int8(-1), system)
    sort!(rows; by=row -> (row.mismatches, row.chromosome, row.position, row.strand))
    result = DataFrame(rows)
    diagnostics = OffTargetSearchDiagnostics(forward_intervals + reverse_intervals, length(forward_positions) + length(reverse_positions), length(rows), Int(max_mismatches), :validated, provenance_record("OffTargetSearchDiagnostics", "CRISPR/enumerate_off_targets"; parameters=(candidate_intervals=forward_intervals + reverse_intervals, fm_candidates=length(forward_positions) + length(reverse_positions), verified_sites=length(rows), max_mismatches=Int(max_mismatches))))
    DataAPI.metadata!(result, "off_target_search_diagnostics", diagnostics; style=:note)
    return _register_crispr_result!(_ctx, result, "enumerate_off_targets"; parents=provenance_parent_ids(guide, index), parameters=(engine="fm_index_hamming", max_mismatches=Int(max_mismatches), max_candidates=Int(max_candidates), system=system.name, row_count=nrow(result), candidate_count=diagnostics.fm_candidates, bulges_searched=false))
end

enumerate_off_targets(guide::AbstractString, index::OffTargetIndex; kwargs...) = enumerate_off_targets(DNASeq(guide; validate=false), index; kwargs...)

"""
    enumerate_off_targets(guide, genome_seq; kwargs...)

Enumerate off-targets directly against an in-memory sequence or genome dictionary using an optimized sliding-window scan.
"""
function enumerate_off_targets(
    guide::BioSequence{DNAAlphabet},
    genome_seq::BioSequence{DNAAlphabet};
    chromosome::String="chr1",
    max_mismatches::Integer=3,
    system::CRISPRSystem=SpCas9,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx),
    kwargs...)

    index = build_offtarget_index(genome_seq; chromosome=chromosome)
    return enumerate_off_targets(guide, index; max_mismatches=max_mismatches, system=system, _ctx=_ctx, kwargs...)
end

enumerate_off_targets(guide::AbstractString, genome_seq::AbstractString; kwargs...) = enumerate_off_targets(DNASeq(guide; validate=false), DNASeq(genome_seq; validate=false); kwargs...)

# ---------------------------------------------------------------------------
# Guide design pipeline
# ---------------------------------------------------------------------------

"""
    design_guides(sequence, system; kwargs...) → DataFrame

Design all possible guide RNAs for a target sequence, score them for on-target efficiency, and filter by quality criteria.
Returns a DataFrame sorted by on-target score.
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
    top_n::Union{Int,Nothing}=nothing,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    sites = find_pam_sites(sequence, system; chromosome=chromosome, _ctx=_ctx)
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
        eff = score_on_target(DNASeq(spacer; validate=false); system=system, _ctx=_ctx)
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

    return _register_crispr_result!(_ctx, df, "design_guides"; parents=provenance_parent_ids(sequence), parameters=(chromosome=chromosome, system=system.name, min_gc=Float64(min_gc), max_gc=Float64(max_gc), min_efficiency=Float64(min_efficiency), guide_count=nrow(df)))
end

design_guides(seq::AbstractString, system::CRISPRSystem=SpCas9; kwargs...) = design_guides(DNASeq(seq; validate=false), system; kwargs...)

# ---------------------------------------------------------------------------
# Guide filtering and ranking
# ---------------------------------------------------------------------------

"""
    filter_guides(guides; min_efficiency=0.4, max_off_targets=10, min_specificity=0.3) → DataFrame

Apply quality filters to a guide DataFrame.
"""
function filter_guides(
    guides::DataFrame;
    min_efficiency::Real=0.4,
    max_off_targets::Int=10,
    min_specificity::Real=0.3,
    exclude_restriction_sites::Vector{String}=String[],
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

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
    return _register_crispr_result!(_ctx, result, "filter_guides"; parents=provenance_parent_ids(guides), parameters=(min_efficiency=Float64(min_efficiency), max_off_targets=max_off_targets, min_specificity=Float64(min_specificity), row_count=nrow(result)))
end

"""
    rank_guides(guides; eff_weight=0.6, spec_weight=0.3, gc_weight=0.1) → DataFrame

Rank guides by a composite multi-objective score.
"""
function rank_guides(
    guides::DataFrame;
    eff_weight::Real=0.6,
    spec_weight::Real=0.3,
    gc_weight::Real=0.1,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    df = copy(guides)
    n  = nrow(df)
    n == 0 && return df

    eff  = hasproperty(df, :on_target_score)   ? Float64.(coalesce.(df.on_target_score, 0.5))   : fill(0.5, n)
    spec = hasproperty(df, :specificity_score) ? Float64.(coalesce.(df.specificity_score, 0.5)) : fill(0.5, n)
    gc   = hasproperty(df, :gc_content) ? 1.0 .- abs.(Float64.(coalesce.(df.gc_content, 0.5)) .- 0.5) ./ 0.5 : fill(0.5, n)

    composite = Float64(eff_weight)*eff .+ Float64(spec_weight)*spec .+ Float64(gc_weight)*gc
    df[!, :composite_score] = composite
    sort!(df, :composite_score, rev=true)
    df[!, :rank] = 1:nrow(df)

    return _register_crispr_result!(_ctx, df, "rank_guides"; parents=provenance_parent_ids(guides), parameters=(eff_weight=Float64(eff_weight), spec_weight=Float64(spec_weight), gc_weight=Float64(gc_weight), row_count=nrow(df)))
end

"""
    guide_specificity_score(guide, off_targets_df) → Float64

Compute a combined specificity score from enumerated off-targets.
"""
function guide_specificity_score(guide::BioSequence{DNAAlphabet}, off_targets_df::DataFrame)
    nrow(off_targets_df) == 0 && return 1.0
    cfd_vals = hasproperty(off_targets_df, :cfd) ? Float64.(off_targets_df.cfd) : fill(0.5, nrow(off_targets_df))
    return clamp(1.0 - sum(cfd_vals) / max(length(cfd_vals) * 10, 1), 0.0, 1.0)
end

guide_specificity_score(guide::AbstractString, off_targets_df::DataFrame) = guide_specificity_score(DNASeq(guide; validate=false), off_targets_df)

# ---------------------------------------------------------------------------
# Base editor guide design
# ---------------------------------------------------------------------------

"""
    design_base_editor_guides(sequence, editor; system=SpCas9, kwargs...) → DataFrame

Design guides for base editing: identifies spacers where the target base falls within the editing window and reports bystander edits.
"""
function design_base_editor_guides(
    sequence::BioSequence{DNAAlphabet},
    editor::EditingWindow;
    system::CRISPRSystem=SpCas9,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx),
    kwargs...)

    guides = design_guides(sequence, system; _ctx=_ctx, kwargs...)
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
        bystander   = [p for p in editor.bystander_positions if p <= glen && p ∉ win_s:win_e && sp[p] == target]

        editable_pos[i]   = edit_pos
        bystander_pos[i]  = bystander
        has_target_base[i] = !isempty(edit_pos)
    end

    guides[!, :editable_positions]  = editable_pos
    guides[!, :bystander_positions] = bystander_pos
    guides[!, :has_target_base]     = has_target_base

    return _register_crispr_result!(_ctx, guides, "design_base_editor_guides"; parents=provenance_parent_ids(sequence), parameters=(editor=editor.editor_name, system=system.name, guide_count=nrow(guides)))
end

design_base_editor_guides(seq::AbstractString, editor::EditingWindow; kwargs...) = design_base_editor_guides(DNASeq(seq; validate=false), editor; kwargs...)

"""
    analyze_editing_window(spacer, editor) → NamedTuple

Report which positions within the editing window contain the target base.
"""
function analyze_editing_window(spacer::BioSequence{DNAAlphabet}, editor::EditingWindow; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    sp   = uppercase(String(spacer))
    glen = length(sp)
    ws   = max(1, editor.window_start)
    we   = min(glen, editor.window_end)
    target = editor.editable_base

    window_seq  = ws <= we <= glen ? sp[ws:we] : ""
    edit_sites  = [p for p in ws:we if p <= glen && sp[p] == target]
    bystanders  = [p for p in editor.bystander_positions if p <= glen && p ∉ ws:we && sp[p] == target]

    result = (
        spacer             = sp,
        editor             = editor.editor_name,
        window_sequence    = window_seq,
        editable_positions = edit_sites,
        bystander_positions = bystanders,
        n_editable         = length(edit_sites),
        n_bystander        = length(bystanders),
        target_base        = editor.editable_base,
        product_base       = editor.edited_base)

    return _register_crispr_result!(_ctx, result, "analyze_editing_window"; parents=provenance_parent_ids(spacer), parameters=(editor=editor.editor_name, editable_count=length(edit_sites), bystander_count=length(bystanders)))
end

analyze_editing_window(spacer::AbstractString, editor::EditingWindow; kwargs...) = analyze_editing_window(DNASeq(spacer; validate=false), editor; kwargs...)

# ---------------------------------------------------------------------------
# Prime editing (pegRNA design)
# ---------------------------------------------------------------------------

"""
    design_pegrna(target_sequence, edit; nick_position=17, pbs_length=13, rt_template_length=15, system=SpCas9) → NamedTuple

Design a prime editing guide RNA (pegRNA) for a desired edit according to Anzalone et al. (2019) Nature 576:149-157.
Returns spacer, PBS (primer binding site), RT template, and full pegRNA sequence.
"""
function design_pegrna(
    target_sequence::BioSequence{DNAAlphabet},
    edit::Union{BioSequence{DNAAlphabet},AbstractString};
    nick_position::Int=17,
    pbs_length::Int=13,
    rt_template_length::Int=15,
    system::CRISPRSystem=SpCas9,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    seq  = uppercase(String(target_sequence))
    edit_seq = uppercase(String(edit))
    n    = length(seq)

    glen = system.protospacer_length
    spacer_end = min(nick_position + glen - 1, n)
    spacer = seq[max(1, spacer_end - glen + 1) : spacer_end]

    # PBS: complementary to single-stranded 3' OH end created by nick (sequence upstream of nick)
    pbs_upstream_start = max(1, nick_position - pbs_length + 1)
    pbs_raw = seq[pbs_upstream_start : nick_position]
    pbs = _reverse_complement_str(pbs_raw)

    # RT template: complementary to desired edit + downstream sequence
    rt_downstream_end = min(n, nick_position + rt_template_length)
    rt_raw = seq[nick_position + 1 : rt_downstream_end]
    rt_template = _reverse_complement_str(edit_seq * rt_raw)

    scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
    full_pegrna = spacer * scaffold * rt_template * pbs
    efficiency = prime_editing_guide_score(DNASeq(spacer; validate=false), DNASeq(pbs; validate=false), DNASeq(rt_template; validate=false))

    result = (
        spacer              = spacer,
        pbs                 = pbs,
        rt_template         = rt_template,
        full_pegrna         = full_pegrna,
        nick_position       = nick_position,
        efficiency_estimate = efficiency)

    return _register_crispr_result!(_ctx, result, "design_pegrna"; parents=provenance_parent_ids(target_sequence), parameters=(nick_position=nick_position, pbs_length=pbs_length, rt_template_length=rt_template_length, system=system.name, efficiency_estimate=efficiency))
end

design_pegrna(seq::AbstractString, edit::AbstractString; kwargs...) = design_pegrna(DNASeq(seq; validate=false), DNASeq(edit; validate=false); kwargs...)

"""
    prime_editing_guide_score(spacer, pbs, rt_template) → Float64

Estimate pegRNA efficiency score [0, 1].
"""
function prime_editing_guide_score(spacer::BioSequence{DNAAlphabet}, pbs::BioSequence{DNAAlphabet}, rt_template::BioSequence{DNAAlphabet})
    spacer_score = score_on_target(spacer)
    pbs_gc       = guide_gc_content(pbs)
    rt_gc        = guide_gc_content(rt_template)
    pbs_len_score = 1.0 - abs(length(pbs) - 13) / 13.0

    score = 0.4 * spacer_score + 0.3 * pbs_gc + 0.2 * rt_gc + 0.1 * pbs_len_score
    return clamp(score, 0.0, 1.0)
end

prime_editing_guide_score(spacer::AbstractString, pbs::AbstractString, rt_template::AbstractString) = prime_editing_guide_score(DNASeq(spacer; validate=false), DNASeq(pbs; validate=false), DNASeq(rt_template; validate=false))

# ---------------------------------------------------------------------------
# HDR template design
# ---------------------------------------------------------------------------

"""
    design_hdr_template(sequence, edit, cut_position; homology_arm_length=80, silent_pam_mutation=true) → NamedTuple

Design an HDR (Homology-Directed Repair) donor template for precise genome editing.
"""
function design_hdr_template(
    sequence::BioSequence{DNAAlphabet},
    edit::Union{BioSequence{DNAAlphabet},AbstractString},
    cut_position::Int;
    homology_arm_length::Int=80,
    silent_pam_mutation::Bool=true,
    system::CRISPRSystem=SpCas9,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    seq  = uppercase(String(sequence))
    edit_str = uppercase(String(edit))
    n    = length(seq)
    cut  = clamp(cut_position, 1, n)

    left_start = max(1, cut - homology_arm_length)
    left_arm   = seq[left_start : cut]
    right_arm  = seq[min(cut+1, n) : min(cut + homology_arm_length, n)]

    if silent_pam_mutation && system.pam == "NGG"
        right_arm = _mutate_pam_silent(right_arm, system)
    end

    full_template = left_arm * edit_str * right_arm
    result = (
        left_arm        = left_arm,
        insert          = edit_str,
        right_arm       = right_arm,
        full_template   = full_template,
        template_length = length(full_template),
        cut_position    = cut)

    return _register_crispr_result!(_ctx, result, "design_hdr_template"; parents=provenance_parent_ids(sequence), parameters=(cut_position=cut, homology_arm_length=homology_arm_length, silent_pam_mutation=silent_pam_mutation, system=system.name, template_length=length(full_template)))
end

design_hdr_template(seq::AbstractString, edit::AbstractString, cut_position::Int; kwargs...) = design_hdr_template(DNASeq(seq; validate=false), edit, cut_position; kwargs...)

function _mutate_pam_silent(arm::AbstractString, system::CRISPRSystem)
    i = findfirst("GG", arm)
    i === nothing && return arm
    chars = collect(arm)
    chars[last(i)] = 'A'   # GG → GA (reduces re-cutting)
    return String(chars)
end

# ---------------------------------------------------------------------------
# Indel prediction (Lindel-like)
# ---------------------------------------------------------------------------

"""
    predict_indels(spacer; n_samples=1000, seed=1) → DataFrame

Predict NHEJ indel outcomes after Cas9 cutting using sequence context and microhomology features.
"""
function predict_indels(
    spacer::BioSequence{DNAAlphabet};
    n_samples::Int=1000,
    seed::Int=1,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    s   = uppercase(String(spacer))
    n   = length(s)

    gc       = guide_gc_content(s)
    cut_region = s[max(1, n-5):n]
    mh_score = _microhomology_score(cut_region)

    p_ins = clamp(0.2 + 0.3 * (1.0 - gc), 0.1, 0.5)
    p_del = 1.0 - p_ins

    rows = NamedTuple[]
    push!(rows, (indel_type="insertion", size=1, frequency=p_ins * 0.7))
    push!(rows, (indel_type="insertion", size=2, frequency=p_ins * 0.2))
    push!(rows, (indel_type="insertion", size=3, frequency=p_ins * 0.1))

    del_sizes  = [1, 2, 3, 5, 7, 10, 15, 20]
    del_probs  = [0.3, 0.2, 0.15, 0.12, 0.1, 0.07, 0.03, 0.03] .* (1.0 + mh_score)
    del_probs ./= sum(del_probs)
    for (sz, prob) in zip(del_sizes, del_probs)
        push!(rows, (indel_type="deletion", size=sz, frequency=p_del * prob))
    end

    df = DataFrame(rows)
    df[!, :frequency] ./= sum(df.frequency)
    sort!(df, :frequency, rev=true)

    return _register_crispr_result!(_ctx, df, "predict_indels"; parents=provenance_parent_ids(spacer), parameters=(n_samples=n_samples, seed=seed, row_count=nrow(df)))
end

predict_indels(spacer::AbstractString; kwargs...) = predict_indels(DNASeq(spacer; validate=false); kwargs...)

function _microhomology_score(seq::AbstractString)
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
function indel_distribution(spacer::BioSequence{DNAAlphabet}; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    df = predict_indels(spacer; _ctx=_ctx)
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

indel_distribution(spacer::AbstractString) = indel_distribution(DNASeq(spacer; validate=false))

# ---------------------------------------------------------------------------
# Genome-wide library design
# ---------------------------------------------------------------------------

"""
    design_library(gene_sequences, system; guides_per_gene=6, kwargs...) → DataFrame

Design a genome-wide CRISPR guide library targeting multiple genes.
"""
function design_library(
    gene_sequences::Dict{String,S},
    system::CRISPRSystem=SpCas9;
    guides_per_gene::Int=6,
    include_controls::Int=100,
    seed::Int=1,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx),
    kwargs...
) where {S<:Union{AbstractString,BioSequence{DNAAlphabet}}}

    rng = MersenneTwister(seed)
    all_guides = DataFrame[]

    for (gene_name, seq) in gene_sequences
        guides = design_guides(seq isa BioSequence ? seq : DNASeq(seq; validate=false), system; _ctx=_ctx, kwargs...)
        isempty(guides) && continue
        guides = rank_guides(guides; _ctx=_ctx)
        n_sel = min(guides_per_gene, nrow(guides))
        sel = guides[1:n_sel, :]
        sel[!, :gene] .= gene_name
        push!(all_guides, sel)
    end

    if include_controls > 0
        ctrl_seqs = [randstring(rng, "ACGT", system.protospacer_length) for _ in 1:include_controls]
        ctrl_df = DataFrame(
            spacer          = ctrl_seqs,
            gene            = fill("non_targeting_control", include_controls),
            on_target_score = fill(0.0, include_controls),
            gc_content      = guide_gc_content.(ctrl_seqs),
            composite_score = fill(0.0, include_controls))
        push!(all_guides, ctrl_df)
    end

    isempty(all_guides) && return DataFrame()
    lib = vcat(all_guides...; cols=:union)
    lib[!, :library_index] = 1:nrow(lib)

    return _register_crispr_result!(_ctx, lib, "design_library"; parents=String[], parameters=(gene_count=length(gene_sequences), guides_per_gene=guides_per_gene, include_controls=include_controls, row_count=nrow(lib)))
end

"""
    library_coverage_stats(library, n_genes) → NamedTuple

Report coverage statistics for a designed guide library.
"""
function library_coverage_stats(library::DataFrame, n_genes::Int; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    total = nrow(library)
    n_ctrl = hasproperty(library, :gene) ? count(g -> g == "non_targeting_control", library.gene) : 0
    n_targeting = total - n_ctrl
    targeting_genes = hasproperty(library, :gene) ? length(unique(filter(g -> g != "non_targeting_control", library.gene))) : 0
    result = (
        total_guides         = total,
        targeting_guides     = n_targeting,
        control_guides       = n_ctrl,
        genes_covered        = targeting_genes,
        gene_coverage_frac   = targeting_genes / max(n_genes, 1),
        mean_guides_per_gene = n_targeting / max(targeting_genes, 1))

    return _register_crispr_result!(_ctx, result, "library_coverage_stats"; parents=provenance_parent_ids(library), parameters=(n_genes=n_genes, total_guides=total, genes_covered=targeting_genes))
end

# ---------------------------------------------------------------------------
# CRISPR screen analysis: negative-binomial guide model and RRA ranking
# ---------------------------------------------------------------------------

struct MageckRRAResult <: AbstractAnalysisResult
    gene_results::DataFrame
    direction::Symbol
    provenance::ResultProvenance
end

struct CRISPRScreenResult <: AbstractAnalysisResult
    guide_results::DataFrame
    gene_results::DataFrame
    qc::DataFrame
    global_dispersion::Float64
    rra::MageckRRAResult
    provenance::ResultProvenance
end

function _median_ratio_size_factors(counts::AbstractMatrix{<:Real})
    any(value -> value < 0 || !isfinite(value), counts) && throw(ArgumentError("screen counts must be finite and nonnegative"))
    counts_float = Float64.(counts)
    log_geomean = vec(mean(log.(counts_float .+ 1.0), dims=2))
    geomean = exp.(log_geomean)
    factors = Float64[]
    for column in axes(counts_float, 2)
        ratios = counts_float[:, column] ./ geomean
        valid = filter(isfinite, ratios)
        push!(factors, isempty(valid) ? 1.0 : median(valid))
    end
    any(value -> value <= 0 || !isfinite(value), factors) && throw(ArgumentError("could not estimate positive library size factors"))
    factors ./= exp(mean(log.(factors)))
    return factors
end

function _bh_adjust(pvalues::AbstractVector{<:Real})
    n = length(pvalues)
    n == 0 && return Float64[]
    order = sortperm(pvalues)
    adjusted = ones(Float64, n)
    running = 1.0
    for rank in n:-1:1
        index = order[rank]
        running = min(running, Float64(pvalues[index]) * n / rank)
        adjusted[index] = clamp(running, 0.0, 1.0)
    end
    return adjusted
end

function _rra_tail_probability(ranks::AbstractVector{<:Integer}, total_guides::Int)
    isempty(ranks) && return 1.0
    total_guides > 0 || throw(ArgumentError("total_guides must be positive"))
    sorted_ranks = sort(Int.(ranks))
    n = length(sorted_ranks)
    probabilities = Float64[]
    for (position, rank) in enumerate(sorted_ranks)
        u = clamp(rank / total_guides, 0.0, 1.0)
        push!(probabilities, cdf(Beta(position, n - position + 1), u))
    end
    return clamp(n * minimum(probabilities), 0.0, 1.0)
end

function _rra_gene_results(guide_results::DataFrame)
    n_guides = nrow(guide_results)
    depleted_order = sortperm(guide_results.log2_fold_change)
    enriched_order = reverse(depleted_order)
    depleted_ranks = zeros(Int, n_guides)
    enriched_ranks = zeros(Int, n_guides)
    for (rank, index) in enumerate(depleted_order)
        depleted_ranks[index] = rank
    end
    for (rank, index) in enumerate(enriched_order)
        enriched_ranks[index] = rank
    end

    rows = NamedTuple[]
    for gene in unique(guide_results.gene)
        indices = findall(==(gene), guide_results.gene)
        depleted_p = _rra_tail_probability(depleted_ranks[indices], n_guides)
        enriched_p = _rra_tail_probability(enriched_ranks[indices], n_guides)
        direction = depleted_p <= enriched_p ? :depleted : :enriched
        rra_p = min(1.0, 2.0 * min(depleted_p, enriched_p))
        push!(rows, (
            gene=String(gene),
            n_guides=length(indices),
            mean_log2_fold_change=mean(guide_results.log2_fold_change[indices]),
            median_log2_fold_change=median(guide_results.log2_fold_change[indices]),
            rra_pvalue=rra_p,
            depleted_rra_pvalue=depleted_p,
            enriched_rra_pvalue=enriched_p,
            direction=String(direction)))
    end
    results = DataFrame(rows)
    results[!, :padj] = _bh_adjust(results.rra_pvalue)
    sort!(results, [:padj, :rra_pvalue, :gene])
    return results
end

"""
    crispr_screen_nb(counts_treatment, counts_control, guide_gene_map; min_mean_count=10, dispersion_prior_weight=10.0)

Fit a guide-level negative-binomial Wald model after DESeq-style median-ratio normalization.
Combines guide-level statistics with Robust Rank Aggregation (RRA) for gene-level calls.
"""
function crispr_screen_nb(
    counts_treatment::AbstractMatrix{<:Real},
    counts_control::AbstractMatrix{<:Real},
    guide_gene_map::DataFrame;
    min_mean_count::Real=10.0,
    dispersion_prior_weight::Real=10.0,
    pseudocount::Real=0.5,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    size(counts_treatment, 1) == size(counts_control, 1) || throw(DimensionMismatch("treatment and control must have equal guide count"))
    size(counts_treatment, 2) >= 2 || throw(ArgumentError("requires at least two treatment replicates"))
    size(counts_control, 2) >= 2 || throw(ArgumentError("requires at least two control replicates"))
    n_guides = size(counts_treatment, 1)
    nrow(guide_gene_map) == n_guides || throw(DimensionMismatch("guide_gene_map must have one row per guide"))
    :gene in propertynames(guide_gene_map) || throw(ArgumentError("guide_gene_map requires a :gene column"))

    combined = hcat(Float64.(counts_control), Float64.(counts_treatment))
    size_factors = _median_ratio_size_factors(combined)
    normalised = combined ./ permutedims(size_factors)
    n_control = size(counts_control, 2)
    control_columns = 1:n_control
    treatment_columns = n_control + 1:size(combined, 2)
    inverse_factor_mean = mean(1.0 ./ size_factors)

    guide_mean = vec(mean(normalised, dims=2))
    guide_variance = [length(row) > 1 ? var(row; corrected=true) : 0.0 for row in eachrow(normalised)]
    raw_dispersion = max.((guide_variance .- guide_mean .* inverse_factor_mean) ./ max.(guide_mean .^ 2, eps()), 0.0)
    expressed = findall(guide_mean .>= Float64(min_mean_count))
    global_dispersion = isempty(expressed) ? 1e-8 : max(median(raw_dispersion[expressed]), 1e-8)
    total_replicates = size(combined, 2)
    shrunk_dispersion = (raw_dispersion .* max(total_replicates - 1, 1) .+ global_dispersion * Float64(dispersion_prior_weight)) ./ (max(total_replicates - 1, 1) + Float64(dispersion_prior_weight))

    guide_names = :guide in propertynames(guide_gene_map) ? String.(guide_gene_map.guide) : ["guide_$(index)" for index in 1:n_guides]
    rows = NamedTuple[]
    for guide_index in 1:n_guides
        control_mean = mean(@view normalised[guide_index, control_columns])
        treatment_mean = mean(@view normalised[guide_index, treatment_columns])
        dispersion = shrunk_dispersion[guide_index]
        control_variance = (control_mean * mean(1.0 ./ size_factors[control_columns]) + dispersion * control_mean^2) / length(control_columns)
        treatment_variance = (treatment_mean * mean(1.0 ./ size_factors[treatment_columns]) + dispersion * treatment_mean^2) / length(treatment_columns)
        log_fold = log2((treatment_mean + Float64(pseudocount)) / (control_mean + Float64(pseudocount)))
        log_se = sqrt(control_variance / (control_mean + Float64(pseudocount))^2 + treatment_variance / (treatment_mean + Float64(pseudocount))^2)
        z_score = log_se > 0 && isfinite(log_se) ? log_fold * log(2) / log_se : 0.0
        pvalue = clamp(2.0 * ccdf(Normal(), abs(z_score)), 0.0, 1.0)
        pass_filter = guide_mean[guide_index] >= Float64(min_mean_count)
        push!(rows, (
            guide=guide_names[guide_index],
            gene=String(guide_gene_map.gene[guide_index]),
            control_mean=control_mean,
            treatment_mean=treatment_mean,
            log2_fold_change=log_fold,
            wald_se=log_se / log(2),
            wald_z=z_score,
            pvalue=pvalue,
            dispersion=dispersion,
            mean_normalized_count=guide_mean[guide_index],
            passes_mean_filter=pass_filter))
    end
    guide_results = DataFrame(rows)
    guide_results[!, :padj] = _bh_adjust(guide_results.pvalue)
    rra_gene_results = _rra_gene_results(guide_results)
    rra = MageckRRAResult(rra_gene_results, :two_sided, provenance_record("MageckRRAResult", "CRISPR/crispr_screen_nb"; parameters=(gene_count=nrow(rra_gene_results), guide_count=n_guides, direction="two_sided")))
    qc = DataFrame(
        sample=vcat(["control_$(index)" for index in 1:length(control_columns)], ["treatment_$(index)" for index in 1:length(treatment_columns)]),
        condition=vcat(fill("control", length(control_columns)), fill("treatment", length(treatment_columns))),
        library_size=vec(sum(combined, dims=1)),
        size_factor=size_factors)
    provenance = provenance_record("CRISPRScreenResult", "CRISPR/crispr_screen_nb"; parameters=(guide_count=n_guides, control_replicates=length(control_columns), treatment_replicates=length(treatment_columns), global_dispersion=global_dispersion, min_mean_count=Float64(min_mean_count), dispersion_prior_weight=Float64(dispersion_prior_weight)))
    result = CRISPRScreenResult(guide_results, rra_gene_results, qc, global_dispersion, rra, provenance)
    return _register_crispr_result!(_ctx, result, "crispr_screen_nb"; parents=provenance_parent_ids(counts_treatment, counts_control, guide_gene_map), parameters=(guide_count=n_guides, gene_count=nrow(rra_gene_results), global_dispersion=global_dispersion, min_mean_count=Float64(min_mean_count)))
end

"""
    crispr_screen_analysis(counts_treatment, counts_control, guide_gene_map; kwargs...) → DataFrame

Compatibility wrapper returning gene-level NB-Wald/RRA results from `crispr_screen_nb`.
"""
function crispr_screen_analysis(counts_treatment, counts_control, guide_gene_map; method::Symbol=:rra, min_reads::Int=10, pseudocount::Real=0.5, kwargs...)
    method == :rra || throw(ArgumentError("only method=:rra is supported"))
    return crispr_screen_nb(counts_treatment, counts_control, guide_gene_map; min_mean_count=min_reads, pseudocount=pseudocount, kwargs...).gene_results
end

mageck_like_test(t, c, g; kwargs...) = crispr_screen_analysis(t, c, g; kwargs...)

# ---------------------------------------------------------------------------
# Genome-wide off-target summary
# ---------------------------------------------------------------------------

"""
    genome_wide_offtarget_summary(guide, ot_df) → NamedTuple

Summary statistics for genome-wide off-target predictions.
"""
function genome_wide_offtarget_summary(guide::BioSequence{DNAAlphabet}, ot_df::DataFrame; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    if isempty(ot_df) || nrow(ot_df) == 0
        result = (guide=String(guide), n_off_targets=0, n_0mm=0, n_1mm=0, n_2mm=0, n_3mm=0, mean_cfd=0.0, specificity_score=1.0)
        return _register_crispr_result!(_ctx, result, "genome_wide_offtarget_summary"; parents=provenance_parent_ids(guide, ot_df), parameters=(off_target_count=0, mean_cfd=0.0))
    end
    mm = hasproperty(ot_df, :mismatches) ? Int.(ot_df.mismatches) : Int[]
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
    return _register_crispr_result!(_ctx, result, "genome_wide_offtarget_summary"; parents=provenance_parent_ids(guide, ot_df), parameters=(off_target_count=nrow(ot_df), mean_cfd=mean(cfd_vals)))
end

genome_wide_offtarget_summary(guide::AbstractString, ot_df::DataFrame; kwargs...) = genome_wide_offtarget_summary(DNASeq(guide; validate=false), ot_df; kwargs...)

# ---------------------------------------------------------------------------
# Interactive HTML Visualizations
# ---------------------------------------------------------------------------

"""
    to_html(result::CRISPRScreenResult) -> String

Generate an interactive HTML5 report for CRISPR screen results featuring interactive Volcano plot, cutoff sliders, and searchable gene hits table.
"""
function to_html(result::CRISPRScreenResult)
    df = result.gene_results
    genes_json = "[" * join(["{\"gene\":\"$(_json_escape(string(r.gene)))\",\"lfc\":$(round(r.mean_log2_fold_change; digits=4)),\"pvalue\":$(round(r.rra_pvalue; digits=6)),\"fdr\":$(round(r.padj; digits=6)),\"dir\":\"$(_json_escape(string(r.direction)))\"}" for r in eachrow(df)], ",") * "]"

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BioToolkit — CRISPR Screen Interactive Dashboard</title>
    <style>
        :root { --bg: #0f172a; --panel: #1e293b; --text: #f8fafc; --accent: #38bdf8; --border: #334155; --glow: #818cf8; }
        body { margin: 0; font-family: system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); padding: 20px; }
        .header { display: flex; align-items: center; justify-content: space-between; background: var(--panel); padding: 16px 24px; border-radius: 12px; border: 1px solid var(--border); margin-bottom: 20px; }
        .title { font-size: 1.4rem; font-weight: 700; color: var(--accent); }
        .grid { display: grid; grid-template-columns: 1fr 340px; gap: 20px; }
        .card { background: var(--panel); border-radius: 12px; border: 1px solid var(--border); padding: 20px; box-shadow: 0 4px 20px rgba(0,0,0,0.3); }
        canvas { width: 100%; height: 500px; display: block; border-radius: 8px; cursor: crosshair; }
        .table-container { max-height: 480px; overflow-y: auto; font-family: monospace; font-size: 13px; }
        .gene-row { display: flex; justify-content: space-between; padding: 8px 12px; border-bottom: 1px solid var(--border); }
        .gene-row:hover { background: #334155; }
        .depleted { color: #f87171; }
        .enriched { color: #4ade80; }
    </style>
</head>
<body>
    <div class="header">
        <div class="title">🎯 CRISPR Screen Volcano & Gene Ranking</div>
        <div>Genes Analysed: <strong>$(nrow(df))</strong></div>
    </div>
    <div class="grid">
        <div class="card">
            <canvas id="volcanoCanvas"></canvas>
        </div>
        <div class="card">
            <h3 style="margin-top: 0; color: var(--accent);">Top Significant Hits</h3>
            <div id="geneTable" class="table-container"></div>
        </div>
    </div>
    <script>
        const data = $(genes_json);
        const canvas = document.getElementById('volcanoCanvas');
        const ctx = canvas.getContext('2d');

        function render() {
            canvas.width = canvas.clientWidth * window.devicePixelRatio;
            canvas.height = 500 * window.devicePixelRatio;
            ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
            const W = canvas.clientWidth;
            const H = 500;
            ctx.clearRect(0, 0, W, H);

            let maxLFC = 1.0, maxP = 1.0;
            data.forEach(d => {
                const logP = -Math.log10(Math.max(d.pvalue, 1e-12));
                if (Math.abs(d.lfc) > maxLFC) maxLFC = Math.abs(d.lfc);
                if (logP > maxP) maxP = logP;
            });
            maxLFC *= 1.1; maxP *= 1.1;

            // Draw axes
            ctx.strokeStyle = '#334155'; ctx.lineWidth = 1;
            ctx.beginPath();
            ctx.moveTo(40, H - 30); ctx.lineTo(W - 10, H - 30);
            ctx.moveTo(40, 10); ctx.lineTo(40, H - 30);
            ctx.stroke();

            // Draw points
            data.forEach(d => {
                const x = 40 + ((d.lfc + maxLFC) / (2 * maxLFC)) * (W - 50);
                const y = (H - 30) - (-Math.log10(Math.max(d.pvalue, 1e-12)) / maxP) * (H - 40);
                const isSig = d.fdr < 0.05;
                ctx.beginPath();
                ctx.arc(x, y, isSig ? 5 : 3, 0, 2 * Math.PI);
                ctx.fillStyle = isSig ? (d.lfc < 0 ? '#f87171' : '#4ade80') : '#64748b';
                ctx.fill();
            });

            // Table
            const sorted = [...data].sort((a,b) => a.fdr - b.fdr);
            document.getElementById('geneTable').innerHTML = sorted.slice(0, 25).map(g => `
                <div class="gene-row">
                    <span>\${g.gene}</span>
                    <span class="\${g.lfc < 0 ? 'depleted' : 'enriched'}">LFC \${g.lfc.toFixed(2)} (FDR \${g.fdr.toExponential(2)})</span>
                </div>
            `).join('');
        }
        setTimeout(render, 50);
    </script>
</body>
</html>
"""
end

"""
    to_html(guides::DataFrame) -> String

Generate an interactive HTML report for a guide library or guide selection table.
"""
function to_html(guides::DataFrame)
    return visualize_guide_library_html(guides)
end

function visualize_crispr_screen_html(result::CRISPRScreenResult)
    return to_html(result)
end

function visualize_guide_library_html(guides::DataFrame)
    rows_json = "[" * join(["{\"spacer\":\"$(_json_escape(string(r.spacer)))\",\"score\":$(hasproperty(r,:on_target_score) ? round(r.on_target_score; digits=3) : 0.5),\"gc\":$(round(r.gc_content; digits=2)),\"system\":\"$(_json_escape(string(hasproperty(r,:system) ? r.system : "SpCas9")))\"}" for r in eachrow(guides)], ",") * "]"

    return """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>BioToolkit — Guide RNA Library Report</title>
    <style>
        body { font-family: system-ui, sans-serif; background: #0f172a; color: #f8fafc; padding: 20px; }
        .card { background: #1e293b; border-radius: 12px; border: 1px solid #334155; padding: 20px; }
        .badge { background: #0284c7; color: white; padding: 2px 8px; border-radius: 4px; font-weight: bold; }
    </style>
</head>
<body>
    <div class="card">
        <h2>🧬 CRISPR Guide RNA Library Summary</h2>
        <p>Total Guides: <span class="badge">$(nrow(guides))</span></p>
    </div>
</body>
</html>
"""
end

function visualize_off_targets_html(off_targets::DataFrame)
    return """
<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><title>BioToolkit — Off-Target Report</title></head>
<body style="background: #0f172a; color: #f8fafc; font-family: system-ui, sans-serif; padding: 20px;">
    <h2>🔍 Off-Target Analysis Report</h2>
    <p>Sites Enumerated: <strong>$(nrow(off_targets))</strong></p>
</body>
</html>
"""
end

function visualize_editing_window_html(spacer::Union{BioSequence{DNAAlphabet},AbstractString}, editor::EditingWindow)
    sp = _escape_html(uppercase(String(spacer)))
    ed_name = _escape_html(string(editor.editor_name))
    info = analyze_editing_window(spacer, editor)
    return """
<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8"><title>BioToolkit — Base Editing Window</title></head>
<body style="background: #0f172a; color: #f8fafc; font-family: system-ui, sans-serif; padding: 20px;">
    <h2>✏️ Base Editor Window View ($(ed_name))</h2>
    <p>Spacer: <code>$(sp)</code></p>
    <p>Target Conversion: <strong>$(_escape_html(string(editor.editable_base))) → $(_escape_html(string(editor.edited_base)))</strong></p>
    <p>Editable Positions: <code>$(join(info.editable_positions, ", "))</code></p>
</body>
</html>
"""
end

end  # module CRISPR
