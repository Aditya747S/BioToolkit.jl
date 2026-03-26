#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Random
using BioToolkit

const DNA_ALPHABET = collect("ACGTN")
const GFF_STRANDS = ['+', '-', '.']
const VCF_BASES = ['A', 'C', 'G', 'T', 'N', 'R', 'Y', 'S', 'W', 'K', 'M']
const TOKEN_ALPHABET = vcat(collect('A':'Z'), collect('a':'z'), collect('0':'9'), ['_', '.', '-', ';', '=', ':', '/', '|'])

random_dna(rng::AbstractRNG, min_len::Integer=0, max_len::Integer=120) = join(rand(rng, DNA_ALPHABET, rand(rng, min_len:max_len)))

function random_token(rng::AbstractRNG; min_len::Integer=0, max_len::Integer=18)
    length = rand(rng, min_len:max_len)
    length == 0 && return ""
    return join(rand(rng, TOKEN_ALPHABET, length))
end

function maybe_drop_field(rng::AbstractRNG, fields::Vector{String})
    if length(fields) > 1 && rand(rng) < 0.2
        deleteat!(fields, rand(rng, 1:length(fields)))
    end
    return fields
end

function random_vcf_line(rng::AbstractRNG)
    chrom = "chr$(rand(rng, 1:22))"
    pos = rand(rng) < 0.8 ? string(rand(rng, -5:250)) : random_token(rng; min_len=1, max_len=6)
    id = rand(rng) < 0.5 ? random_token(rng; min_len=0, max_len=12) : "."
    ref = string(rand(rng, VCF_BASES))
    alt = rand(rng) < 0.8 ? string(rand(rng, VCF_BASES)) : random_token(rng; min_len=1, max_len=4)
    qual = rand(rng) < 0.5 ? "." : (rand(rng) < 0.7 ? string(rand(rng) * 100) : random_token(rng; min_len=1, max_len=8))
    fields = maybe_drop_field(rng, String[chrom, pos, id, ref, alt, qual])
    rand(rng) < 0.2 && push!(fields, random_token(rng; min_len=0, max_len=10))
    return join(fields, '\t')
end

function random_gff_attributes(rng::AbstractRNG)
    rand(rng) < 0.3 && return "."

    pairs = String[]
    for _ in 1:rand(rng, 0:4)
        key = random_token(rng; min_len=1, max_len=8)
        values = [random_token(rng; min_len=0, max_len=10) for _ in 1:rand(rng, 1:3)]
        push!(pairs, "$(key)=$(join(values, ','))")
    end

    return isempty(pairs) ? "." : join(pairs, ';')
end

function random_gff_line(rng::AbstractRNG)
    chrom = "chr$(rand(rng, 1:22))"
    source = random_token(rng; min_len=1, max_len=10)
    feature = rand(rng) < 0.7 ? rand(rng, ["gene", "exon", "CDS", "mRNA"]) : random_token(rng; min_len=1, max_len=8)
    start = rand(rng) < 0.8 ? string(rand(rng, -10:500)) : random_token(rng; min_len=1, max_len=6)
    stop = rand(rng) < 0.8 ? string(rand(rng, -10:500)) : random_token(rng; min_len=1, max_len=6)
    score = rand(rng) < 0.4 ? "." : (rand(rng) < 0.7 ? string(rand(rng) * 100) : random_token(rng; min_len=1, max_len=7))
    strand = string(rand(rng, GFF_STRANDS))
    phase = rand(rng) < 0.5 ? "." : (rand(rng) < 0.7 ? string(rand(rng, -1:4)) : random_token(rng; min_len=1, max_len=4))
    fields = maybe_drop_field(rng, String[chrom, source, feature, start, stop, score, strand, phase, random_gff_attributes(rng)])
    return join(fields, '\t')
end

function random_bed_line(rng::AbstractRNG)
    chrom = "chr$(rand(rng, 1:22))"
    start = rand(rng) < 0.8 ? string(rand(rng, -10:500)) : random_token(rng; min_len=1, max_len=6)
    stop = rand(rng) < 0.8 ? string(rand(rng, -10:500)) : random_token(rng; min_len=1, max_len=6)
    fields = maybe_drop_field(rng, String[chrom, start, stop])
    rand(rng) < 0.2 && push!(fields, random_token(rng; min_len=0, max_len=10))
    return join(fields, '\t')
end

function assert_parser_behavior(line::AbstractString, parser)
    try
        parser(line)
        return :ok
    catch err
        err isa ArgumentError || rethrow()
        return :argument_error
    end
end

function fuzz_parser(generator, parser; iterations::Integer=5_000, rng::AbstractRNG=Random.default_rng())
    stats = Dict(:ok => 0, :argument_error => 0)
    for _ in 1:iterations
        outcome = assert_parser_behavior(generator(rng), parser)
        stats[outcome] += 1
    end
    return stats
end

function fuzz_reader(file_name::AbstractString, generator, reader; iterations::Integer=500, rng::AbstractRNG=Random.default_rng())
    ok = 0
    argument_error = 0

    mktempdir() do dir
        for _ in 1:iterations
            path = joinpath(dir, file_name)
            open(path, "w") do io
                for _ in 1:rand(rng, 1:6)
                    rand(rng) < 0.2 && write(io, "# comment\n")
                    write(io, generator(rng), '\n')
                end
            end

            try
                reader(path)
                ok += 1
            catch err
                err isa ArgumentError || rethrow()
                argument_error += 1
            end
        end
    end

    return (ok = ok, argument_error = argument_error)
end

function main()
    seed = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : rand(1:1_000_000_000)
    iterations = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5_000
    rng = MersenneTwister(seed)

    println("seed=", seed)
    println("iterations=", iterations)
    println("example_dna=", random_dna(rng, 24, 24))

    vcf_stats = fuzz_parser(random_vcf_line, BioToolkit.parse_vcf_record; iterations=iterations, rng=rng)
    gff_stats = fuzz_parser(random_gff_line, BioToolkit.parse_gff_record; iterations=iterations, rng=rng)
    bed_stats = fuzz_parser(random_bed_line, BioToolkit.parse_bed_record; iterations=iterations, rng=rng)

    reader_iterations = max(100, iterations ÷ 10)
    vcf_file_stats = fuzz_reader("fuzz.vcf", random_vcf_line, BioToolkit.read_vcf; iterations=reader_iterations, rng=rng)
    gff_file_stats = fuzz_reader("fuzz.gff", random_gff_line, BioToolkit.read_gff; iterations=reader_iterations, rng=rng)
    bed_file_stats = fuzz_reader("fuzz.bed", random_bed_line, BioToolkit.read_bed; iterations=reader_iterations, rng=rng)

    println("parse_vcf_record => ", vcf_stats)
    println("parse_gff_record => ", gff_stats)
    println("parse_bed_record => ", bed_stats)
    println("read_vcf => ", vcf_file_stats)
    println("read_gff => ", gff_file_stats)
    println("read_bed => ", bed_file_stats)

    return nothing
end

main()