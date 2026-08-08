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
using DataAPI
using Statistics
using LinearAlgebra
using Random
using Distributions

# Import biotypes for type-safe sequence handling
using ..BioToolkit: AASeq, AminoAcidAlphabet, BioSequence, DNAAlphabet, DNASeq, SummarizedExperiment, assay, colData
using ..BioToolkit: AbstractAnalysisResult, ResultProvenance, provenance_record
using ..BioToolkit: ProvenanceContext, ProvenanceParams, ThreadSafeProvenanceContext, active_provenance_context, new_provenance_id, provenance_parent_ids, provenance_result!, register_provenance!

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

"""
    FMIndex

Checkpointed FM-index with sampled suffix-array locations. The index supports
exact and bounded Hamming-distance backward search over a DNA reference.
`bwt` and suffix-array samples are immutable after construction, making the
query path deterministic and safe to share across threads.
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

A contig-aware FM-index for repeated CRISPR off-target searches. Reference
sequences are retained for PAM verification and strand-correct extraction after
candidate discovery; the FM-index itself performs the genome-scale lookup.
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
    max_build_bases === nothing || total_bases <= max_build_bases || throw(ArgumentError("reference has $(total_bases) bases, exceeding max_build_bases=$(max_build_bases); use a larger build environment or explicitly raise the limit"))

    text = UInt8[]
    starts = Int[]
    stops = Int[]
    for (contig_name, sequence) in zip(contig_names, contig_sequences)
        isempty(contig_name) && throw(ArgumentError("contig names must not be empty"))
        occursin(' ', sequence) && throw(ArgumentError("reference sequence contains the FM-index sentinel byte"))
        occursin('', sequence) && throw(ArgumentError("reference sequence contains the FM-index contig separator byte"))
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

Build a reusable FM-index for a DNA reference. Construction uses a native
prefix-doubling suffix-array builder and is memory-intensive; build once and
query many guides with `enumerate_off_targets(guide, index)`. The query engine
is bounded-Hamming only and therefore reports `bulges=0` explicitly.
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
        length(positions) + (right - left + 1) <= max_candidates || throw(ArgumentError("FM-index search exceeded max_candidates=$(max_candidates); narrow max_mismatches or raise the explicit limit"))
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

Enumerate PAM-valid off-targets through FM-index backward search with bounded
Hamming distance. Candidate discovery is sublinear in reference length for
selective guide suffixes; candidate verification preserves exact strand and PAM
semantics. Bulges are not searched and are reported as zero.
"""
function enumerate_off_targets(
    guide::BioSequence{DNAAlphabet},
    index::OffTargetIndex;
    max_mismatches::Integer=3,
    max_candidates::Integer=100_000,
    system::CRISPRSystem=SpCas9,
    prov_ctx=nothing,
    _ctx=active_provenance_context(prov_ctx))

    0 <= max_mismatches <= length(guide) || throw(ArgumentError("max_mismatches must lie between 0 and the guide length"))
    max_candidates > 0 || throw(ArgumentError("max_candidates must be positive"))
    guide_string = uppercase(String(guide))
    all(symbol -> symbol in _OFFTARGET_SEARCH_ALPHABET, codeunits(guide_string)) || throw(ArgumentError("FM-index off-target search accepts unambiguous A/C/G/T guides only"))

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

Compatibility path for a single in-memory reference. For repeated guide queries,
build one `OffTargetIndex` and call the indexed overload directly.
"""
function enumerate_off_targets(guide::BioSequence{DNAAlphabet}, genome_seq::BioSequence{DNAAlphabet}; chromosome::String="chr1", kwargs...)
    index = build_offtarget_index(genome_seq; chromosome=chromosome)
    return enumerate_off_targets(guide, index; kwargs...)
end

enumerate_off_targets(guide::AbstractString, genome_seq::AbstractString; kwargs...) = enumerate_off_targets(DNASeq(guide; validate=false), DNASeq(genome_seq; validate=false); kwargs...)

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
# CRISPR screen analysis: negative-binomial guide model and RRA ranking
# ---------------------------------------------------------------------------

"""
    MageckRRAResult

Gene-level robust-rank-aggregation output. The RRA p-value is a two-sided,
Bonferroni-corrected minimum beta-order-statistic probability; it is reported
separately from the negative-binomial guide-level evidence.
"""
struct MageckRRAResult <: AbstractAnalysisResult
    gene_results::DataFrame
    direction::Symbol
    provenance::ResultProvenance
end

"""
    CRISPRScreenResult

Result of a library-size-normalized negative-binomial CRISPR screen analysis.
Guide-level Wald statistics and dispersion estimates are retained so gene calls
can be audited rather than treated as an opaque rank list.
"""
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
    crispr_screen_nb(counts_treatment, counts_control, guide_gene_map;
                     min_mean_count=10, dispersion_prior_weight=10.0)

Fit a guide-level two-group negative-binomial Wald model after DESeq-style
median-ratio library normalization. Dispersion is estimated from all samples
and shrunk toward a robust global estimate. Gene hits are ranked with
bidirectional robust rank aggregation of guide log-fold changes.

This is an explicit NB-Wald/RRA implementation, not a claim of bit-for-bit
MAGeCK MLE parity. It requires at least two samples per condition; screen-level
QC, guide dispersion, and all intermediate statistics are returned.
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

    size(counts_treatment, 1) == size(counts_control, 1) || throw(DimensionMismatch("treatment and control must have the same guide count"))
    size(counts_treatment, 2) >= 2 || throw(ArgumentError("negative-binomial screen analysis requires at least two treatment replicates"))
    size(counts_control, 2) >= 2 || throw(ArgumentError("negative-binomial screen analysis requires at least two control replicates"))
    n_guides = size(counts_treatment, 1)
    nrow(guide_gene_map) == n_guides || throw(DimensionMismatch("guide_gene_map must have one row per guide"))
    :gene in propertynames(guide_gene_map) || throw(ArgumentError("guide_gene_map requires a :gene column"))
    min_mean_count >= 0 || throw(ArgumentError("min_mean_count must be nonnegative"))
    dispersion_prior_weight >= 0 || throw(ArgumentError("dispersion_prior_weight must be nonnegative"))
    pseudocount > 0 || throw(ArgumentError("pseudocount must be positive"))

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

Compatibility wrapper returning the gene-level NB-Wald/RRA table from
`crispr_screen_nb`. The former simplified t-test/geometric-mean implementation
is retired; callers needing guide evidence should call `crispr_screen_nb`.
"""
function crispr_screen_analysis(counts_treatment, counts_control, guide_gene_map; method::Symbol=:rra, min_reads::Int=10, pseudocount::Real=0.5, kwargs...)
    method == :rra || throw(ArgumentError("only method=:rra is supported by the NB-Wald/RRA implementation"))
    return crispr_screen_nb(counts_treatment, counts_control, guide_gene_map; min_mean_count=min_reads, pseudocount=pseudocount, kwargs...).gene_results
end

"""
    mageck_like_test(treatment_counts, control_counts, guide_df) → DataFrame

Compatibility wrapper for the explicit NB-Wald/RRA screen model.
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
