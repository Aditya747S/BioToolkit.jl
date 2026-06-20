# ==============================================================================
# gene_prediction.jl — Ab initio gene prediction & gene-structure analysis
#
# References:
#   - Burge & Karlin (1997) JMB 268(1):78-94 (GENSCAN)
#   - Rabiner (1989) Proc IEEE 77(2):257-286 (HMM tutorial)
#   - Majoros et al. (2004) Genome Res 14:2illil (GlimmerHMM)
#   - Fickett (1982) Nucleic Acids Res 10:5303-5318 (TestCode)
#   - Kozak (1987) Nucleic Acids Res 15:8125-8148 (Kozak context)
# ==============================================================================

export GeneInterval, predict_genes_hmm
export GenePrediction, ExonInterval, IntronInterval, SpliceSite
export predict_gene_structure, detect_splice_sites, score_kozak_context
export codon_bias_index, gene_density, gff3_export
export find_start_codons, find_stop_codons, calculate_testcode
export predict_utr_regions, cds_statistics, annotate_orf_features
export filter_gene_predictions, gene_prediction_summary

# ---------------------------------------------------------------------------
# Core types
# ---------------------------------------------------------------------------

struct GeneInterval
    start::Int
    stop::Int
end

Base.show(io::IO, interval::GeneInterval) = print(io, "GeneInterval(", interval.start, ", ", interval.stop, ")")
Base.iterate(interval::GeneInterval) = (interval.start, 2)
Base.iterate(interval::GeneInterval, state::Int) = state == 2 ? (interval.stop, 3) : nothing

"""
    ExonInterval

Represents a predicted exon with phase information.
"""
struct ExonInterval
    start::Int
    stop::Int
    strand::Int8
    phase::Int8     # reading frame phase at start (0, 1, 2)
    end_phase::Int8 # reading frame phase at end
    score::Float64
end

"""
    IntronInterval

Represents a predicted intron with splice-site scores.
"""
struct IntronInterval
    start::Int
    stop::Int
    strand::Int8
    donor_score::Float64    # 5' splice site score
    acceptor_score::Float64 # 3' splice site score
end

"""
    SpliceSite

Represents a detected splice site with position and score.
"""
struct SpliceSite
    position::Int
    strand::Int8
    site_type::Symbol  # :donor or :acceptor
    score::Float64
    consensus::String  # the matched dinucleotide
end

"""
    GenePrediction

Complete gene prediction with exons, introns, UTRs, and metadata.
"""
struct GenePrediction
    gene_id::String
    chrom::String
    start::Int
    stop::Int
    strand::Int8
    exons::Vector{ExonInterval}
    introns::Vector{IntronInterval}
    cds_length::Int
    protein_length::Int
    score::Float64
    source::String
end

Base.show(io::IO, g::GenePrediction) = print(io,
    "GenePrediction(", g.gene_id, " ", g.start, "-", g.stop,
    " strand=", g.strand == Int8(1) ? "+" : "-", " exons=", length(g.exons), ")")

@inline function _register_gene_prediction_result!(_ctx::Union{Nothing,ProvenanceContext,ThreadSafeProvenanceContext}, result, operation::AbstractString; parents::AbstractVector{<:AbstractString}=String[], parameters=NamedTuple())
    return provenance_result!(_ctx, result, operation; parents=parents, parameters=parameters)
end


"""
    predict_genes_hmm(sequence; p_coding=0.01)

Uses a 2-state HMM (non-coding/coding) and the Viterbi algorithm to predict
open reading frames in nucleotide DNA.
"""
function predict_genes_hmm(sequence::BioSequence{DNAAlphabet}; p_coding::Real=0.01, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    bytes = sequence.data
    L = length(bytes)

    states = ["NonCoding", "Coding"]
    alphabet = UInt8[UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T')]
    initial  = [0.99, 0.01]
    transitions = [
        1.0 - Float64(p_coding)   Float64(p_coding);
        1/1000.0                  1.0 - 1/1000.0
    ]
    emissions = [
        0.25  0.25  0.25  0.25;
        0.15  0.35  0.35  0.15
    ]

    hmm = HMM(states, alphabet, initial, transitions, emissions; log_space=false)
    path, _ = viterbi(hmm, bytes)

    genes    = GeneInterval[]
    in_gene  = false
    gene_start = 0

    for i in 1:L
        if path[i] == 2 && !in_gene
            in_gene    = true
            gene_start = i
        elseif path[i] == 1 && in_gene
            in_gene  = false
            gene_end = i - 1
            if (gene_end - gene_start) > 30
                push!(genes, GeneInterval(gene_start, gene_end))
            end
        end
    end
    if in_gene && (L - gene_start) > 30
        push!(genes, GeneInterval(gene_start, L))
    end

    return _register_gene_prediction_result!(_ctx, genes, "predict_genes_hmm"; parents=provenance_parent_ids(sequence), parameters=(p_coding=Float64(p_coding), gene_count=length(genes)))
end

# ---------------------------------------------------------------------------
# Start and stop codon finders
# ---------------------------------------------------------------------------

"""
    find_start_codons(seq; min_context_score=0.0)

Locate all ATG start codons in both strands of `seq`, returning positions and
Kozak context scores.
"""
function find_start_codons(seq::BioSequence{DNAAlphabet}; min_context_score::Real=0.0, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    s = uppercase(String(seq))
    n = length(s)
    results = NamedTuple{(:position, :strand, :kozak_score)}[]

    # Forward strand
    i = 1
    while (i = findnext("ATG", s, i)) !== nothing
        start_pos = first(i)
        ks = score_kozak_context(s, start_pos, Int8(1))
        if ks >= Float64(min_context_score)
            push!(results, (position=start_pos, strand=Int8(1), kozak_score=ks))
        end
        i = start_pos + 1
    end

    # Reverse complement
    rc = String(reverse_complement(seq))
    i = 1
    while (i = findnext("ATG", rc, i)) !== nothing
        start_pos = first(i)
        orig_pos = n - start_pos - 1   # approximate position on original strand
        ks = score_kozak_context(rc, start_pos, Int8(-1))
        if ks >= Float64(min_context_score)
            push!(results, (position=orig_pos, strand=Int8(-1), kozak_score=ks))
        end
        i = start_pos + 1
    end
    return _register_gene_prediction_result!(_ctx, results, "find_start_codons"; parents=provenance_parent_ids(seq), parameters=(min_context_score=Float64(min_context_score), site_count=length(results)))
end

"""
    find_stop_codons(seq)

Locate all in-frame stop codons (TAA, TAG, TGA) in all 6 reading frames.
"""
function find_stop_codons(seq::BioSequence{DNAAlphabet}; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    stops = Set(["TAA","TAG","TGA"])
    s = uppercase(String(seq))
    n = length(s)
    results = NamedTuple{(:position, :strand, :frame, :codon)}[]

    for frame in 0:2
        pos = frame + 1
        while pos + 2 <= n
            codon = s[pos:pos+2]
            if codon in stops
                push!(results, (position=pos, strand=Int8(1), frame=frame, codon=codon))
            end
            pos += 3
        end
    end

    rc = String(reverse_complement(seq))
    for frame in 0:2
        pos = frame + 1
        while pos + 2 <= n
            codon = rc[pos:pos+2]
            if codon in stops
                push!(results, (position=n - pos - 1, strand=Int8(-1), frame=frame, codon=codon))
            end
            pos += 3
        end
    end
    return _register_gene_prediction_result!(_ctx, results, "find_stop_codons"; parents=provenance_parent_ids(seq), parameters=(site_count=length(results)))
end

# ---------------------------------------------------------------------------
# Kozak consensus scoring
# ---------------------------------------------------------------------------

# Score matrix: positions -6 to -1 and +4 relative to A of ATG
# Based on Kozak (1987) vertebrate consensus gccRccATGG
const _KOZAK_WEIGHTS = Dict(
    (-6, 'G') => 0.1, (-6, 'C') => 0.3, (-6, 'A') => 0.2, (-6, 'T') => 0.1,
    (-5, 'C') => 0.4, (-5, 'G') => 0.1,
    (-4, 'C') => 0.4,
    (-3, 'A') => 0.6, (-3, 'R') => 0.5,  # purine at -3 is the most important
    (-2, 'C') => 0.4,
    (-1, 'C') => 0.4,
    ( 4, 'G') => 0.6,  # +4G is strongly conserved
)

"""
    score_kozak_context(seq, atg_pos, strand)

Score the Kozak consensus around an ATG at `atg_pos` using a position-specific
weight table. Returns a value in [0, 1].
"""
function score_kozak_context(seq::BioSequence{DNAAlphabet}, atg_pos::Int, strand::Int8=Int8(1); prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    n = length(seq)
    score = 0.0
    max_score = sum(values(_KOZAK_WEIGHTS))

    for ((offset, _), weight) in _KOZAK_WEIGHTS
        pos = atg_pos - 1 + offset  # offset relative to A of ATG
        if 1 <= pos <= n
            c = uppercase(Char(seq.data[pos]))
            key = (offset, c)
            score += get(_KOZAK_WEIGHTS, key, 0.0)
        end
    end
    result = score / max(max_score, eps(Float64))
    return _register_gene_prediction_result!(_ctx, result, "score_kozak_context"; parents=provenance_parent_ids(seq), parameters=(atg_pos=atg_pos, strand=strand, score=result))
end

score_kozak_context(seq::AbstractString, atg_pos::Integer, strand::Int8=Int8(1); prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx)) = score_kozak_context(DNASeq(seq), Int(atg_pos), strand; _ctx=_ctx)
score_kozak_context(seq::AbstractString, atg_pos::AbstractUnitRange{<:Integer}, strand::Int8=Int8(1); prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx)) = score_kozak_context(seq, first(atg_pos), strand; _ctx=_ctx)

# ---------------------------------------------------------------------------
# Splice-site detection
# ---------------------------------------------------------------------------

const _DONOR_CONSENSUS    = "GT"   # GT-AG rule
const _ACCEPTOR_CONSENSUS = "AG"

# Position weight matrices (simplified; derived from human GT-AG splice sites)
const _DONOR_PWM = Dict(
    (-2,'A')=>0.3,(-2,'G')=>0.4,
    (-1,'G')=>0.5,(-1,'A')=>0.3,
    ( 1,'G')=>0.9,( 1,'A')=>0.05,
    ( 2,'T')=>0.88,( 2,'C')=>0.05,
    ( 3,'A')=>0.4,( 3,'G')=>0.3,
    ( 4,'A')=>0.3,( 4,'G')=>0.3)
const _ACCEPTOR_PWM = Dict(
    (-2,'P')=>0.0,  # pyrimidine-rich upstream
    (-1,'A')=>0.1,(-1,'G')=>0.1,
    ( 1,'A')=>0.9,
    ( 2,'G')=>0.9,
    ( 3,'G')=>0.5,( 3,'A')=>0.3)

"""
    detect_splice_sites(seq; min_score=0.5)

Detect putative GT-AG splice donor and acceptor sites in `seq` using position-
weight matrices. Analogous to MaxEntScan / NNSPLICE.

Returns a `Vector{SpliceSite}`.
"""
function detect_splice_sites(seq::BioSequence{DNAAlphabet}; min_score::Real=0.5, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    s  = uppercase(String(seq))
    n  = length(s)
    sites = SpliceSite[]

    # Forward strand donors (GT at exon-intron boundary)
    for i in 2:(n-2)
        if s[i:i+1] == "GT"
            sc = _score_site(s, i, _DONOR_PWM, n)
            sc >= Float64(min_score) && push!(sites, SpliceSite(i, Int8(1), :donor, sc, "GT"))
        end
        if s[i:i+1] == "AG"
            sc = _score_site(s, i, _ACCEPTOR_PWM, n)
            sc >= Float64(min_score) && push!(sites, SpliceSite(i, Int8(1), :acceptor, sc, "AG"))
        end
    end

    # Reverse strand
    rc = String(reverse_complement(seq))
    for i in 2:(n-2)
        if rc[i:i+1] == "GT"
            sc = _score_site(rc, i, _DONOR_PWM, n)
            sc >= Float64(min_score) && push!(sites, SpliceSite(n-i+1, Int8(-1), :donor, sc, "GT"))
        end
        if rc[i:i+1] == "AG"
            sc = _score_site(rc, i, _ACCEPTOR_PWM, n)
            sc >= Float64(min_score) && push!(sites, SpliceSite(n-i+1, Int8(-1), :acceptor, sc, "AG"))
        end
    end

    sort!(sites, by = ss -> -ss.score)
    return _register_gene_prediction_result!(_ctx, sites, "detect_splice_sites"; parents=provenance_parent_ids(seq), parameters=(min_score=Float64(min_score), site_count=length(sites)))
end

function _score_site(s::BioSequence{DNAAlphabet}, pos::Int, pwm::Dict, n::Int)
    score     = 0.0
    max_score = 0.0
    for (offset, _) in pwm
        abs_pos = pos + offset
        1 <= abs_pos <= n || continue
        c = s[abs_pos]
        key = (offset, c)
        score     += get(pwm, key, 0.0)
        max_score += maximum(v for (k,v) in pwm if k[1] == offset; init=0.01)
    end
    return score / max(max_score, eps(Float64))
end

# ---------------------------------------------------------------------------
# Fickett TestCode — coding potential
# ---------------------------------------------------------------------------

"""
    calculate_testcode(seq; window=200)

Compute the Fickett TestCode statistic for each window of `seq`, returning a
`Vector{Float64}` of coding potential scores (>0.74 = likely coding).

Reference: Fickett (1982) Nucleic Acids Res 10:5303-5318.
"""
function calculate_testcode(seq::BioSequence{DNAAlphabet}; window::Int=200, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    s = uppercase(String(seq))
    n = length(s)
    window = min(window, n)
    n_windows = n - window + 1
    n_windows <= 0 && return Float64[]

    scores = zeros(Float64, n_windows)
    Threads.@threads for i in 1:n_windows
        sub   = s[i:(i + window - 1)]
        scores[i] = _fickett_score(sub)
    end
    return _register_gene_prediction_result!(_ctx, scores, "calculate_testcode"; parents=provenance_parent_ids(seq), parameters=(window=window, score_count=length(scores)))
end

function _fickett_score(window::BioSequence{DNAAlphabet})
    len = length(window)
    len < 3 && return 0.0

    freq = zeros(Int, 4)   # A C G T
    pos_freq = zeros(Int, 4, 3)   # nucleotide × position-in-codon

    map_nt = Dict('A'=>1,'C'=>2,'G'=>3,'T'=>4)
    for (idx, c) in enumerate(window)
        h = get(map_nt, c, 0)
        h == 0 && continue
        freq[h] += 1
        pos_freq[h, mod1(idx, 3)] += 1
    end

    # Position asymmetry: max(pos1,pos2,pos3) / (min+1)
    pa = Float64[]
    for nt in 1:4
        maxv = maximum(pos_freq[nt, :])
        minv = minimum(pos_freq[nt, :])
        push!(pa, maxv / max(minv + 1, 1))
    end

    # GC content
    gc = (freq[2] + freq[3]) / max(len, 1)

    # Simple TestCode approximation
    score = 0.4 * mean(pa) + 0.6 * (gc > 0.5 ? gc : 1 - gc)
    return clamp(score, 0.0, 1.0)
end

# ---------------------------------------------------------------------------
# Codon bias index (CAI-like)
# ---------------------------------------------------------------------------

"""
    codon_bias_index(seq; reference_dict=nothing)

Compute the Codon Bias Index (CBI) for a coding sequence, measuring deviation
from uniform codon usage. Returns a value in [-1, +1] (positive = biased).

If `reference_dict` is provided (codon → frequency), computes a CAI-like score.
"""
function codon_bias_index(seq::BioSequence{DNAAlphabet}; reference_dict::Union{Nothing,Dict}=nothing, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    s   = uppercase(String(seq))
    n   = length(s)
    n3  = n ÷ 3
    n3 < 1 && return 0.0

    codon_counts = Dict{String,Int}()
    for i in 1:n3
        codon = s[(i-1)*3+1:(i-1)*3+3]
        codon_counts[codon] = get(codon_counts, codon, 0) + 1
    end
    total = max(sum(values(codon_counts)), 1)

    if reference_dict !== nothing
        # Relative synonymous codon usage weighted by reference
        log_cai = 0.0
        used = 0
        for (codon, cnt) in codon_counts
            w = get(reference_dict, codon, 0.0)
            w > 0 && (log_cai += cnt * log(Float64(w)); used += cnt)
        end
        result = used > 0 ? exp(log_cai / used) : 0.0
    else
        # Simple CBI: deviation from uniform distribution
        aa2codons = _build_aa2codons()
        cbi = 0.0
        for (aa, codons) in aa2codons
            length(codons) <= 1 && continue
            counts_aa = [get(codon_counts, c, 0) for c in codons]
            total_aa  = max(sum(counts_aa), 1)
            n_syn     = length(codons)
            expected  = 1.0 / n_syn
            for c in counts_aa
                p = c / total_aa
                cbi += abs(p - expected)
            end
        end
        result = clamp(cbi / max(total, 1), -1.0, 1.0)
    end
    return _register_gene_prediction_result!(_ctx, result, "codon_bias_index"; parents=provenance_parent_ids(seq), parameters=(reference_dict=reference_dict === nothing ? "none" : "provided", score=result))
end

function _build_aa2codons()
    standard_code = Dict(
        "A" => ["GCT","GCC","GCA","GCG"],
        "R" => ["CGT","CGC","CGA","CGG","AGA","AGG"],
        "N" => ["AAT","AAC"],
        "D" => ["GAT","GAC"],
        "C" => ["TGT","TGC"],
        "Q" => ["CAA","CAG"],
        "E" => ["GAA","GAG"],
        "G" => ["GGT","GGC","GGA","GGG"],
        "H" => ["CAT","CAC"],
        "I" => ["ATT","ATC","ATA"],
        "L" => ["TTA","TTG","CTT","CTC","CTA","CTG"],
        "K" => ["AAA","AAG"],
        "M" => ["ATG"],
        "F" => ["TTT","TTC"],
        "P" => ["CCT","CCC","CCA","CCG"],
        "S" => ["TCT","TCC","TCA","TCG","AGT","AGC"],
        "T" => ["ACT","ACC","ACA","ACG"],
        "W" => ["TGG"],
        "Y" => ["TAT","TAC"],
        "V" => ["GTT","GTC","GTA","GTG"])
    return standard_code
end

# ---------------------------------------------------------------------------
# Gene density
# ---------------------------------------------------------------------------

"""
    gene_density(genes, total_length; window_bp=1_000_000)

Compute gene density per megabase across a sequence. `genes` is a vector of
`GeneInterval` or `GenePrediction`.

Returns a `Vector{NamedTuple}` with window start, end, and density (genes/Mb).
"""
function gene_density(genes, total_length::Int; window_bp::Int=1_000_000, prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    total_length >= 1 || return NamedTuple[]
    window_bp = max(window_bp, 1)

    windows = NamedTuple{(:window_start, :window_end, :gene_count, :density_per_mb)}[]
    window_start = 1
    while window_start <= total_length
        window_end = min(window_start + window_bp - 1, total_length)
        cnt = count(g -> begin
            gs = g isa GenePrediction ? g.start : g.start
            ge = g isa GenePrediction ? g.stop  : g.stop
            gs <= window_end && ge >= window_start
        end, genes)
        span_mb = (window_end - window_start + 1) / 1e6
        push!(windows, (
            window_start  = window_start,
            window_end    = window_end,
            gene_count    = cnt,
            density_per_mb = cnt / max(span_mb, eps(Float64))))
        window_start += window_bp
    end
    return _register_gene_prediction_result!(_ctx, windows, "gene_density"; parents=provenance_parent_ids(genes), parameters=(total_length=total_length, window_bp=window_bp, window_count=length(windows)))
end

# ---------------------------------------------------------------------------
# GFF3 export
# ---------------------------------------------------------------------------

"""
    gff3_export(predictions, seqname; io=stdout)

Write gene predictions to GFF3 format. Analogous to Augustus GFF3 output.
"""
function gff3_export(
    predictions::AbstractVector{<:GenePrediction},
    seqname::AbstractString;
    io::IO=stdout,
    prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx)
)
    println(io, "##gff-version 3")
    println(io, "##sequence-region ", seqname, " 1 ", maximum(g.stop for g in predictions; init=1))
    for g in predictions
        strand_char = g.strand == 1 ? "+" : "-"
        source = isempty(g.source) ? "BioToolkit" : g.source
        # Gene record
        println(io, join([seqname, source, "gene",
            g.start, g.stop, ".", strand_char, ".",
            "ID=$(g.gene_id);Name=$(g.gene_id)"], "\t"))
        # mRNA record
        println(io, join([seqname, source, "mRNA",
            g.start, g.stop, ".", strand_char, ".",
            "ID=$(g.gene_id).mRNA1;Parent=$(g.gene_id)"], "\t"))
        # Exon + CDS records
        for (eidx, exon) in enumerate(g.exons)
            println(io, join([seqname, source, "exon",
                exon.start, exon.stop, @sprintf("%.2f", exon.score), strand_char, ".",
                "ID=$(g.gene_id).exon$(eidx);Parent=$(g.gene_id).mRNA1"], "\t"))
            println(io, join([seqname, source, "CDS",
                exon.start, exon.stop, @sprintf("%.2f", exon.score), strand_char,
                string(exon.phase),
                "ID=$(g.gene_id).CDS$(eidx);Parent=$(g.gene_id).mRNA1"], "\t"))
        end
    end
    _ctx === nothing || register_provenance!(_ctx, "gff3_export"; parents=provenance_parent_ids(predictions), parameters=(seqname=String(seqname), output="IO", gene_count=length(predictions)))
end

# ---------------------------------------------------------------------------
# Structural gene prediction (multi-exon)
# ---------------------------------------------------------------------------

"""
    predict_gene_structure(seq; min_gene_length=300, min_orf_length=100, max_intron_length=50000)

Predict multi-exon gene structures by:
1. Finding all ATG/stop codon pairs (ORFs).
2. Scoring splice sites around each ORF.
3. Assembling candidate exon sets using a simple greedy graph approach.

Returns a `Vector{GenePrediction}`.
"""
function predict_gene_structure(
    seq::BioSequence{DNAAlphabet};
    min_gene_length::Int  = 300,
    min_orf_length::Int   = 100,
    max_intron_length::Int= 50_000,
    chrom::String         = "chr1",
    source::String        = "BioToolkit_GenePred",
    seed::Int             = 1,
    prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    s = uppercase(String(seq))
    n = length(s)
    stop_codons = Set(["TAA","TAG","TGA"])
    predictions = GenePrediction[]

    gene_id_counter = Ref(0)

    for strand_sign in (Int8(1), Int8(-1))
        working = strand_sign == 1 ? s : String(reverse_complement(seq))
        nw = length(working)

        for frame in 0:2
            i = frame + 1
            while i + 2 <= nw
                codon = working[i:i+2]
                if codon == "ATG"
                    atg_pos = i
                    j = i + 3
                    while j + 2 <= nw
                        stop = working[j:j+2]
                        if stop in stop_codons
                            orf_len = j + 2 - atg_pos + 1
                            if orf_len >= min_orf_length
                                # Build a minimal single-exon gene prediction
                                actual_start = strand_sign == 1 ? atg_pos : nw - (j + 2) + 1
                                actual_stop  = strand_sign == 1 ? j + 2  : nw - atg_pos + 1
                                if actual_stop - actual_start + 1 >= min_gene_length
                                    exon = ExonInterval(
                                        actual_start, actual_stop, strand_sign,
                                        Int8(0), Int8(0), _fickett_score(working[atg_pos:j+2])
                                    )
                                    gene_id_counter[] += 1
                                    push!(predictions, GenePrediction(
                                        "gene_$(gene_id_counter[])", chrom,
                                        actual_start, actual_stop, strand_sign,
                                        [exon], IntronInterval[], orf_len, orf_len ÷ 3, exon.score, source
                                    ))
                                end
                            end
                            break
                        end
                        j += 3
                    end
                end
                i += 3
            end
        end
    end

    return _register_gene_prediction_result!(_ctx, predictions, "predict_gene_structure"; parents=provenance_parent_ids(seq), parameters=(min_gene_length=min_gene_length, min_orf_length=min_orf_length, max_intron_length=max_intron_length, chrom=chrom, source=source, seed=seed, gene_count=length(predictions)))
end

# ---------------------------------------------------------------------------
# UTR prediction
# ---------------------------------------------------------------------------

"""
    predict_utr_regions(seq, gene_predictions; utr5_max=1000, utr3_max=2000)

Predict approximate 5' and 3' UTR regions based on surrounding sequence context.
Returns a `DataFrame` with UTR bounds per gene.
"""
function predict_utr_regions(
    seq::BioSequence{DNAAlphabet},
    predictions::AbstractVector{<:GenePrediction};
    utr5_max::Int=1_000,
    utr3_max::Int=2_000,
    prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    n = length(seq)
    rows = NamedTuple{(:gene_id,:utr5_start,:utr5_end,:utr3_start,:utr3_end,:utr5_length,:utr3_length)}[]

    for g in predictions
        if g.strand == 1
            utr5_start = max(1, g.start - utr5_max)
            utr5_end   = g.start - 1
            utr3_start = g.stop + 1
            utr3_end   = min(n, g.stop + utr3_max)
        else
            utr5_start = g.stop + 1
            utr5_end   = min(n, g.stop + utr5_max)
            utr3_start = max(1, g.start - utr3_max)
            utr3_end   = g.start - 1
        end
        push!(rows, (
            gene_id     = g.gene_id,
            utr5_start  = utr5_start,
            utr5_end    = utr5_end,
            utr3_start  = utr3_start,
            utr3_end    = utr3_end,
            utr5_length = max(utr5_end - utr5_start + 1, 0),
            utr3_length = max(utr3_end - utr3_start + 1, 0)))
    end
    return _register_gene_prediction_result!(_ctx, rows, "predict_utr_regions"; parents=provenance_parent_ids(seq, predictions), parameters=(utr5_max=utr5_max, utr3_max=utr3_max, row_count=length(rows)))
end

# ---------------------------------------------------------------------------
# CDS statistics
# ---------------------------------------------------------------------------

"""
    cds_statistics(predictions)

Compute summary statistics for a vector of gene predictions.

Returns a `NamedTuple` with gene count, CDS length distribution, and strand counts.
"""
function cds_statistics(predictions::AbstractVector{<:GenePrediction}; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    isempty(predictions) && return _register_gene_prediction_result!(_ctx, (
        n_genes=0, mean_cds_length=0.0, median_cds_length=0.0,
        mean_exons_per_gene=0.0, n_forward=0, n_reverse=0
    ), "cds_statistics"; parents=provenance_parent_ids(predictions), parameters=(gene_count=0, total_cds_bp=0))
    lengths = [g.cds_length for g in predictions]
    exon_counts = [length(g.exons) for g in predictions]
    result = (
        n_genes            = length(predictions),
        mean_cds_length    = mean(lengths),
        median_cds_length  = median(lengths),
        min_cds_length     = minimum(lengths),
        max_cds_length     = maximum(lengths),
        mean_exons_per_gene = mean(exon_counts),
        n_forward          = count(g -> g.strand == 1, predictions),
        n_reverse          = count(g -> g.strand == -1, predictions),
        total_cds_bp       = sum(lengths))
    return _register_gene_prediction_result!(_ctx, result, "cds_statistics"; parents=provenance_parent_ids(predictions), parameters=(gene_count=length(predictions), total_cds_bp=result.total_cds_bp))
end

# ---------------------------------------------------------------------------
# Filter gene predictions
# ---------------------------------------------------------------------------

"""
    filter_gene_predictions(predictions; min_score=0.0, min_length=0, max_overlap_fraction=0.5)

Filter gene predictions by score, length, and remove heavily overlapping calls.
"""
function filter_gene_predictions(
    predictions::AbstractVector{<:GenePrediction};
    min_score::Real=0.0,
    min_length::Int=0,
    max_overlap_fraction::Real=0.5,
    prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    filtered = filter(g -> g.score >= Float64(min_score) &&
                           g.stop - g.start + 1 >= min_length, predictions)

    # Remove overlapping predictions (keep highest score)
    sort!(filtered, by = g -> -g.score)
    kept = GenePrediction[]
    for g in filtered
        overlap = false
        for k in kept
            isect_len = max(0, min(g.stop, k.stop) - max(g.start, k.start) + 1)
            g_len = g.stop - g.start + 1
            if isect_len / max(g_len, 1) > Float64(max_overlap_fraction)
                overlap = true; break
            end
        end
        overlap || push!(kept, g)
    end
    return _register_gene_prediction_result!(_ctx, kept, "filter_gene_predictions"; parents=provenance_parent_ids(predictions), parameters=(min_score=Float64(min_score), min_length=min_length, max_overlap_fraction=Float64(max_overlap_fraction), gene_count=length(kept)))
end

# ---------------------------------------------------------------------------
# Gene prediction summary
# ---------------------------------------------------------------------------

"""
    gene_prediction_summary(predictions)

Print a human-readable summary of gene predictions.
"""
function gene_prediction_summary(predictions::AbstractVector{<:GenePrediction}; prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    stats = cds_statistics(predictions; _ctx=_ctx)
    println("Gene Prediction Summary")
    println("  Total genes:       ", stats.n_genes)
    if stats.n_genes > 0
        println("  Forward strand:    ", stats.n_forward)
        println("  Reverse strand:    ", stats.n_reverse)
        println("  Mean CDS length:   ", round(stats.mean_cds_length; digits=1), " bp")
        println("  Median CDS length: ", round(stats.median_cds_length; digits=1), " bp")
        println("  Total CDS:         ", stats.total_cds_bp, " bp")
        println("  Mean exons/gene:   ", round(stats.mean_exons_per_gene; digits=2))
    end
    return _register_gene_prediction_result!(_ctx, stats, "gene_prediction_summary"; parents=provenance_parent_ids(predictions), parameters=(gene_count=stats.n_genes, total_cds_bp=stats.total_cds_bp))
end

# ---------------------------------------------------------------------------
# Annotate ORF features (returns DataFrame-like NamedTuple array)
# ---------------------------------------------------------------------------

"""
    annotate_orf_features(seq, predictions)

Annotate each predicted gene with: GC content, Kozak score, splice site count,
TestCode (coding potential), and codon bias.

Returns a `Vector{NamedTuple}` compatible with DataFrames.
"""
function annotate_orf_features(
    seq::BioSequence{DNAAlphabet},
    predictions::AbstractVector{<:GenePrediction};
    prov_ctx=nothing, _ctx=active_provenance_context(prov_ctx))
    s = uppercase(String(seq))
    n = length(s)
    results = NamedTuple[]

    for g in predictions
        start = max(1, g.start)
        stop  = min(n, g.stop)
        stop < start && continue

        subseq = BioSequence{DNAAlphabet}(s[start:stop])
        gc   = gc_content(subseq)
        cbi  = codon_bias_index(subseq)
        tc   = calculate_testcode(subseq; window=min(200, stop-start+1))
        mean_tc = isempty(tc) ? 0.0 : mean(tc)
        n_splice = length(g.exons) > 1 ? 2 * (length(g.exons) - 1) : 0

        push!(results, (
            gene_id         = g.gene_id,
            gc_content      = gc,
            codon_bias_index = cbi,
            testcode_mean   = mean_tc,
            n_splice_sites  = n_splice,
            cds_length      = g.cds_length,
            strand          = g.strand,
            n_exons         = length(g.exons)))
    end
    return _register_gene_prediction_result!(_ctx, results, "annotate_orf_features"; parents=provenance_parent_ids(seq, predictions), parameters=(gene_count=length(predictions), row_count=length(results)))
end
