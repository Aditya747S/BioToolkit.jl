using Downloads
using Pkg

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const BED_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
const BENCHMARK_ENVDIR = joinpath(homedir(), ".cache", "BioToolkit", "bed_compare_env")

function setup_benchmark_environment(envdir::AbstractString)
    mkpath(envdir)
    for filename in ("Project.toml", "Manifest.toml")
        env_file = joinpath(envdir, filename)
        isfile(env_file) && rm(env_file; force=true)
    end
    Pkg.activate(envdir)
    Pkg.develop(path=REPO_ROOT)
    Pkg.add("Arrow")
    Pkg.add("BED")
    Pkg.resolve()
    Pkg.instantiate()
end

function elapsed_ms(f::Function)
    GC.gc()
    start_time = time_ns()
    result = f()
    return (time_ns() - start_time) / 1e6, result
end

function prepare_bed_file(bed_gz_path::AbstractString, bed_path::AbstractString)
    open(pipeline(`gzip -dc $bed_gz_path`), "r") do input_io
        open(bed_path, "w") do out
            for raw_line in eachline(input_io)
                fields = split(raw_line, '\t')
                length(fields) < 8 && continue
                write(out, fields[6], '\t', fields[7], '\t', fields[8], '\n')
            end
        end
    end
    return bed_path
end

function run_benchmark(envdir::AbstractString, bed_path::AbstractString)
    benchmark_code = """
using BioToolkit, BED

function elapsed_ms(f)
    GC.gc()
    start_time = time_ns()
    result = f()
    return (time_ns() - start_time) / 1e6, result
end

function benchmark_parser(parse_function::Function, repetitions::Int)
    timings = Float64[]
    result = nothing
    for _ in 1:repetitions
        elapsed, result = elapsed_ms(parse_function)
        push!(timings, elapsed)
    end
    return timings, result
end

function parse_with_bed(path::AbstractString)
    records = BED.Record[]
    open(path, "r") do io
        for raw_line in eachline(io)
            isempty(raw_line) && continue
            startswith(raw_line, '#') && continue
            startswith(raw_line, "track") && continue
            startswith(raw_line, "browser") && continue
            push!(records, BED.Record(raw_line))
        end
    end
    return records
end

    bio_times, bio_records = benchmark_parser(() -> BioToolkit.read_bed($(repr(bed_path))), 2)
    bed_times, bed_records = benchmark_parser(() -> parse_with_bed($(repr(bed_path))), 2)

bio_count = length(bio_records)
bed_count = length(bed_records)

@assert bio_count == bed_count

println("full BED parse benchmark")
println("  file=rmsk_hg38.bed")
println("  url=", $(repr(BED_URL)))
println("  record_count=", bio_count)
println("  bio_min_ms=", round(minimum(bio_times), digits=4))
println("  bedjl_min_ms=", round(minimum(bed_times), digits=4))
println("  bio_mean_ms=", round(sum(bio_times) / length(bio_times), digits=4))
println("  bedjl_mean_ms=", round(sum(bed_times) / length(bed_times), digits=4))
println("  ratio=", round(minimum(bio_times) / minimum(bed_times), digits=3))
"""

    benchmark_cmd = `env JULIA_PKG_PRECOMPILE_AUTO=0 $(Base.julia_cmd()) --startup-file=no --project=$envdir -e $benchmark_code`
    run(benchmark_cmd)
end

function main()
    mkpath(dirname(BENCHMARK_ENVDIR))
    mktempdir() do dir
        bed_gz_path = joinpath(dir, "hg38.trf.bed.gz")
        bed_path = joinpath(dir, "hg38.trf.bed")

        Downloads.download(BED_URL, bed_gz_path)
        prepare_bed_file(bed_gz_path, bed_path)
        setup_benchmark_environment(BENCHMARK_ENVDIR)
        run_benchmark(BENCHMARK_ENVDIR, bed_path)
    end
end

main()