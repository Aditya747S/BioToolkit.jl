using Downloads
using Pkg

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const MMCIF_URL = "https://files.rcsb.org/download/7KCR.cif"
const BENCHMARK_ENVDIR = joinpath(homedir(), ".cache", "BioToolkit", "mmcif_compare_env")

function setup_benchmark_environment(envdir::AbstractString)
    mkpath(envdir)
    for filename in ("Project.toml", "Manifest.toml")
        env_file = joinpath(envdir, filename)
        isfile(env_file) && rm(env_file; force=true)
    end
    Pkg.activate(envdir)
    Pkg.develop(path=REPO_ROOT)
    Pkg.add("Arrow")
    Pkg.add("PDBTools")
    Pkg.resolve()
    Pkg.instantiate()
end

function run_benchmark(envdir::AbstractString, cif_path::AbstractString)
    benchmark_code = """
using BioToolkit, PDBTools

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

bio_times, bio_structure = benchmark_parser(() -> BioToolkit.read_mmcif($(repr(cif_path))), 5)
pdb_times, pdb_atoms = benchmark_parser(() -> PDBTools.read_mmcif($(repr(cif_path))), 5)

bio_atom_count = length(BioToolkit.structure_atoms(bio_structure))
pdb_atom_count = length(pdb_atoms)

@assert bio_atom_count == pdb_atom_count

println("mmCIF parse benchmark")
println("  file=7KCR.cif")
println("  url=", $(repr(MMCIF_URL)))
println("  atom_count=", bio_atom_count)
println("  bio_min_ms=", round(minimum(bio_times), digits=4))
println("  pdbtools_min_ms=", round(minimum(pdb_times), digits=4))
println("  bio_mean_ms=", round(sum(bio_times) / length(bio_times), digits=4))
println("  pdbtools_mean_ms=", round(sum(pdb_times) / length(pdb_times), digits=4))
println("  ratio=", round(minimum(bio_times) / minimum(pdb_times), digits=3))
"""

    benchmark_cmd = `env JULIA_PKG_PRECOMPILE_AUTO=0 $(Base.julia_cmd()) --startup-file=no --project=$envdir -e $benchmark_code`
    run(benchmark_cmd)
end

function main()
    mkpath(dirname(BENCHMARK_ENVDIR))
    mktempdir() do dir
        cif_path = joinpath(dir, "7KCR.cif")

        Downloads.download(MMCIF_URL, cif_path)
        setup_benchmark_environment(BENCHMARK_ENVDIR)

        run_benchmark(BENCHMARK_ENVDIR, cif_path)
    end
end

main()