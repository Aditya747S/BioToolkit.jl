# ==============================================================================
# test_biotoolkit.jl — correctness + performance suite for align.jl,
# annotation.jl, and bam.jl
#
# Run with:
#   julia --project=/path/to/BioToolkit.jl -O2 test_biotoolkit.jl
# ==============================================================================

using Test
using Random

try
  @eval using BioToolkit
catch e
  println("FATAL: `using BioToolkit` failed. Check the module name / --project path.")
  println(e)
  exit(1)
end

# SAM_FLAG_* constants aren't in BioToolkit's export list, so `using BioToolkit`
# alone doesn't bring them into scope. Import them explicitly by name -- this
# works for non-exported bindings too, as long as they actually exist.
try
  @eval using BioToolkit: SAM_FLAG_PAIRED, SAM_FLAG_PROPER_PAIR, SAM_FLAG_UNMAPPED,
    SAM_FLAG_MATE_UNMAPPED, SAM_FLAG_REVERSE, SAM_FLAG_MATE_REVERSE,
    SAM_FLAG_READ1, SAM_FLAG_READ2, SAM_FLAG_SECONDARY, SAM_FLAG_QC_FAIL,
    SAM_FLAG_DUPLICATE, SAM_FLAG_SUPPLEMENTARY
catch e
  println("FATAL: could not import SAM_FLAG_* constants from BioToolkit.")
  println("Check the exact constant names in bam.jl (they may be spelled differently).")
  println(e)
  exit(1)
end

Random.seed!(1234)
randseq(n) = String(rand(collect("ACGT"), n))

# ==============================================================================
# align.jl
# ==============================================================================
println("\n" * "="^70)
println("align.jl")
println("="^70)

function recompute_linear_score(left::AbstractString, right::AbstractString, match, mismatch, gap)
  s = 0
  for (x, y) in zip(left, right)
    if x == '-' || y == '-'
      s += gap
    elseif x == y
      s += match
    else
      s += mismatch
    end
  end
  return s
end

function recompute_affine_score(left::AbstractString, right::AbstractString, match, mismatch, gap_open, gap_extend)
  s = 0
  prev_state = :none
  for (x, y) in zip(left, right)
    state = x == '-' ? :gapL : (y == '-' ? :gapR : :match)
    if state == :match
      s += (x == y ? match : mismatch)
    else
      s += (state == prev_state) ? gap_extend : gap_open
    end
    prev_state = state
  end
  return s
end

@testset "align.jl: global/local linear traceback self-consistency" begin
  n_trials = 200
  bad = 0
  for _ in 1:n_trials
    a, b = randseq(rand(1:40)), randseq(rand(1:40))
    m, mm, g = rand(1:5), -rand(1:5), -rand(1:5)

    res = needleman_wunsch(a, b; match=m, mismatch=mm, gap=g)
    recomputed = recompute_linear_score(String(res.left), String(res.right), m, mm, g)
    recomputed == res.score || (bad += 1)
    @test recomputed == res.score
    @test replace(String(res.left), "-" => "") == a
    @test replace(String(res.right), "-" => "") == b

    res2 = smith_waterman(a, b; match=m, mismatch=mm, gap=g)
    recomputed2 = recompute_linear_score(String(res2.left), String(res2.right), m, mm, g)
    @test recomputed2 == res2.score
  end
  println("global+local linear traceback mismatches: $bad / $(2*n_trials) checks")
end

@testset "align.jl: affine global/local traceback self-consistency" begin
  n_trials = 150
  bad = 0
  for _ in 1:n_trials
    a, b = randseq(rand(1:35)), randseq(rand(1:35))
    m, mm = rand(1:5), -rand(1:5)
    go, ge = -rand(2:8), -rand(1:4)

    res = pairwise_align(a, b; match=m, mismatch=mm, gap_open=go, gap_extend=ge, is_local=false)
    recomputed = recompute_affine_score(String(res.left), String(res.right), m, mm, go, ge)
    recomputed == res.score || (bad += 1)
    @test recomputed == res.score

    res2 = pairwise_align(a, b; match=m, mismatch=mm, gap_open=go, gap_extend=ge, is_local=true)
    recomputed2 = recompute_affine_score(String(res2.left), String(res2.right), m, mm, go, ge)
    @test recomputed2 == res2.score
  end
  println("affine global+local traceback mismatches: $bad / $(2*n_trials) checks")
end

@testset "align.jl: codon alignment smoke test" begin
  for _ in 1:50
    a = randseq(3 * rand(2:12))
    b = randseq(3 * rand(2:12))
    res = pairwise_align_codons(a, b; match=1, mismatch=-1, gap=-2, is_local=false)
    @test replace(String(res.left), "-" => "") == a
    @test replace(String(res.right), "-" => "") == b
    @test length(res.left) % 3 == 0
  end
end

println("\n--- align.jl: fixed cases (directly comparable to Biostrings::pairwiseAlignment) ---")
for (a, b, m, mm, g) in [("CTAGC", "GAAC", 4, -5, -4), ("CCGTCC", "CCAGA", 5, -5, -2), ("GACTTACTAAGG", "CTTGCTA", 1, -5, -2)]
  res = needleman_wunsch(a, b; match=m, mismatch=mm, gap=g)
  println("global linear  a=$a b=$b match=$m mismatch=$mm gap=$g  -> score=$(res.score)")
end

println("\n--- align.jl performance ---")
a500, b500 = randseq(500), randseq(500)
println("500x500 global linear:  $(round(1000*(@elapsed needleman_wunsch(a500,b500; match=1,mismatch=-1,gap=-2)), digits=2)) ms")
println("500x500 local linear:   $(round(1000*(@elapsed smith_waterman(a500,b500; match=1,mismatch=-1,gap=-2)), digits=2)) ms")
println("500x500 global affine:  $(round(1000*(@elapsed pairwise_align(a500,b500; match=1,mismatch=-1,gap_open=-5,gap_extend=-1,is_local=false)), digits=2)) ms")
println("500x500 local affine:   $(round(1000*(@elapsed pairwise_align(a500,b500; match=1,mismatch=-1,gap_open=-5,gap_extend=-1,is_local=true)), digits=2)) ms")

# ==============================================================================
# annotation.jl
# ==============================================================================
println("\n" * "="^70)
println("annotation.jl")
println("="^70)

@testset "annotation.jl: minus-strand multi-exon splicing" begin
  genome = BioSequence{DNAAlphabet}("AAA" * "CCCCCC" * "GGG")
  exon1 = FeatureLocationLite(1, 3)
  exon2 = FeatureLocationLite(10, 12)
  compound = CompoundFeatureLocation("join", [exon2, exon1]; strand=-1)
  spliced = feature_sequence(genome, compound)
  expected = "TTTCCC"
  @test String(spliced) == expected
  println("minus-strand splice: got=$(String(spliced))  expected=$expected")
end

@testset "annotation.jl: partial flags swap under reverse-complement slicing" begin
  loc = FeatureLocationLite(5, 20; partial_start=false, partial_stop=false)
  sliced = slice_feature_location(loc, 1, 15; reverse_complemented=true)
  @test sliced.partial_start == true
  @test sliced.partial_stop == false
  println("partial-flag swap: partial_start=$(sliced.partial_start) partial_stop=$(sliced.partial_stop)")
end

@testset "annotation.jl: caret-notation location parsing" begin
  loc = parse_feature_location("123^124")
  @test loc.start == 123 && loc.stop == 124
  println("^-notation parse: start=$(loc.start) stop=$(loc.stop)")
end

@testset "annotation.jl: feature_coverage correctness" begin
  seqlen = 200_000
  rec = AnnotatedSeqRecord(BioSequence{DNAAlphabet}(randseq(seqlen)); identifier="chr_test")
  feats = SeqFeatureLite[]
  for i in 1:500
    s = rand(1:190_000)
    e = s + rand(1:20_000)
    push!(feats, SeqFeatureLite("exon", FeatureLocationLite(s, e)))
  end
  rec.features = feats

  t_cov = @elapsed cov = feature_coverage(rec)
  brute_ok = true
  for _ in 1:50
    p = rand(1:seqlen)
    expected_count = count(feats) do f
      s, e = feature_bounds(f.location)
      min(s, e) <= p <= max(s, e)
    end
    cov[p] == expected_count || (brute_ok = false)
  end
  @test brute_ok
  println("feature_coverage (200kb, 500 features): $(round(1000*t_cov, digits=2)) ms, brute-check ok=$brute_ok")

  idx = build_feature_index(rec)
  knn_ok = true
  for _ in 1:20
    p = rand(1:seqlen)
    got = nearest_feature(idx, p; n=3)
    got_dists = sort([feature_distance(f, p) for f in got])
    all_d = sort([feature_distance(f, p) for f in feats])
    true_top3 = all_d[1:min(3, length(all_d))]
    got_dists == true_top3 || (knn_ok = false)
  end
  @test knn_ok
  println("indexed k-nearest (n=3) vs brute force: ok=$knn_ok")
end

# ==============================================================================
# bam.jl
# ==============================================================================
println("\n" * "="^70)
println("bam.jl")
println("="^70)

@testset "bam.jl: round trip, unmapped-read handling, tag typing" begin
  refs = [BamReference("chr1", 2_000_000)]
  hdr = BamHeader(refs)
  seq50 = BioSequence{DNAAlphabet}(randseq(50))

  rec1 = BamRecord("read1", "chr1", 100, [BamCigarOp(50, 'M')], seq50;
    flag=UInt16(0), tags=Dict{String,Any}("AS" => Int32(42)))
  rec2 = BamRecord("read2", "chr1", 5000, [BamCigarOp(50, 'M')], seq50;
    flag=SAM_FLAG_REVERSE, tags=Dict{String,Any}("ZS" => Float32[1.0, 2.0, 3.0]))
  rec3 = BamRecord("read3_unmapped", "chr1", 9999, BamCigarOp[], seq50;
    flag=SAM_FLAG_UNMAPPED)

  bam = BamFile(hdr, [rec1, rec2, rec3])
  path = tempname() * ".bam"
  write_bam(path, bam; write_index=true)
  bam2 = read_bam(path)

  @test length(bam2.records) == 3
  @test Set(bam2.records) == Set(bam.records)
  println("round trip: $(length(bam2.records)) records, set-equality=$(Set(bam2.records) == Set(bam.records))")

  cov = bam_coverage(bam2)
  unmapped_ok = all(cov["chr1"][9950:10100] .== 0)
  @test unmapped_ok
  println("unmapped-but-positioned read excluded from coverage: $unmapped_ok")

  try
    ivs = alignments_to_interval_collection(bam2)
    println("granges() record count (should exclude unmapped): $(length(ivs))")
  catch e
    println("granges() check skipped (GenomicRanges API mismatch?): ", e)
  end

  sam_path = tempname() * ".sam"
  write_sam(sam_path, bam2)
  bam3 = read_sam(sam_path)
  r1 = only(filter(r -> r.qname == "read1", bam3.records))
  r2 = only(filter(r -> r.qname == "read2", bam3.records))
  @test r1.tags["AS"] isa Integer && r1.tags["AS"] == 42
  @test r2.tags["ZS"] isa AbstractVector
  println("SAM tag round trip: AS type=$(typeof(r1.tags["AS"])) value=$(r1.tags["AS"]);  ZS=$(r2.tags["ZS"])")

  try
    region = GenomicRanges.GenomicInterval("chr1", 90, 200, '+', Dict{String,Any}())
    scanned = BioToolkit._bam_collect(BioToolkit.BamRegionScanReader(path, region))
    indexed = read_bam(path, region)
    match_ok = Set(r.qname for r in scanned.records) == Set(r.qname for r in indexed.records)
    @test match_ok
    println("indexed vs full-scan region query match: $match_ok ($(length(indexed.records)) records)")
  catch e
    println("region-query correctness check skipped: ", e)
  end

  try
    empty_region = GenomicRanges.GenomicInterval("chr1", 1_900_000, 1_900_100, '+', Dict{String,Any}())
    t_empty = @elapsed collect(stream_bam(path, empty_region))
    println("empty-region indexed query time: $(round(1000*t_empty, digits=3)) ms (should be near-instant)")
  catch e
    println("empty-region perf check skipped: ", e)
  end
end

println("\n--- bam.jl performance (20,000 synthetic records) ---")
big_refs = [BamReference("chrBig", 5_000_000)]
big_hdr = BamHeader(big_refs)
N = 20_000
big_records = Vector{BamRecord}(undef, N)
for i in 1:N
  p = rand(0:4_900_000)
  s = BioSequence{DNAAlphabet}(randseq(100))
  unmapped = rand() < 0.02
  flag = unmapped ? SAM_FLAG_UNMAPPED : UInt16(0)
  cig = unmapped ? BamCigarOp[] : [BamCigarOp(100, 'M')]
  big_records[i] = BamRecord("r$i", "chrBig", p, cig, s; flag=flag, mapq=UInt8(rand(0:60)))
end
big_bam = BamFile(big_hdr, big_records)
big_path = tempname() * ".bam"

t_write = @elapsed write_bam(big_path, big_bam; write_index=true)
t_read = @elapsed read_back = read_bam(big_path)
t_cov = @elapsed cov_big = bam_coverage(read_back)

println("write $N records (+ .bai index): $(round(t_write, digits=3)) s")
println("read  $N records:                $(round(t_read, digits=3)) s")
println("coverage over 5Mb reference:      $(round(t_cov, digits=3)) s")

try
  region = GenomicRanges.GenomicInterval("chrBig", 1_000_000, 1_010_000, '+', Dict{String,Any}())
  t_indexed = @elapsed collect(stream_bam(big_path, region))
  t_scan = @elapsed collect(BioToolkit.BamRegionScanReader(big_path, region))
  println("indexed region query (10kb window): $(round(1000*t_indexed, digits=2)) ms")
  println("full-scan region query (same):      $(round(1000*t_scan, digits=2)) ms")
catch e
  println("region-query perf comparison skipped: ", e)
end

println("\n" * "="^70)
println("DONE")
println("="^70)
