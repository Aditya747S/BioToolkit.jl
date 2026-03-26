pushfirst!(LOAD_PATH, abspath(joinpath(@__DIR__, "..", "..")))

using BioToolkit
include(joinpath(@__DIR__, "..", "..", "scripts", "benchmark_helpers.jl"))

using Random

function build_population(population_index::Int, individual_count::Int, locus_count::Int)
    individuals = PopGenIndividual{Int}[]
    for individual_index in 1:individual_count
        loci = Locus{Int}[]
        for locus_index in 1:locus_count
            base_value = isodd(population_index + individual_index + locus_index) ? 1 : 2
            alternate_value = base_value == 1 ? 2 : 1
            genotype = if (population_index + individual_index + locus_index) % 7 == 0
                (base_value, alternate_value)
            elseif (population_index + locus_index) % 5 == 0
                (alternate_value, alternate_value)
            else
                (base_value, base_value)
            end
            push!(loci, Locus{Int}(genotype))
        end
        push!(individuals, PopGenIndividual{Int}("P$(population_index)_I$(individual_index)", loci))
    end
    return Population{Int}("Pop$(population_index)", individuals)
end

function build_populations(population_count::Int, individual_count::Int, locus_count::Int)
    [build_population(population_index, individual_count, locus_count) for population_index in 1:population_count]
end

function write_genepop_fixture(path::AbstractString, populations)
    write_genepop(populations, path; title="BioToolkit popgen benchmark")
    return path
end

function measure_ms(label::AbstractString, reps::Int, func)
    total_ms, _ = repeat_elapsed_ms(func, reps)
    ms = total_ms / reps
    println("  ", label, "=", round(ms; digits=4))
    return ms
end

function main()
    Random.seed!(2024)

    population_count = 3
    individual_count = 240
    locus_count = 24
    populations = build_populations(population_count, individual_count, locus_count)

    fixture_path = tempname() * ".gen"
    parsed_path = tempname() * ".gen"
    write_genepop_fixture(fixture_path, populations)

    parsed_populations = read_genepop(fixture_path)
    parsed_record = read_genepop_record(fixture_path)
    @assert !isempty(parsed_populations)

    read_ms = measure_ms("julia_genepop_read_ms", 10, () -> read_genepop(fixture_path))
    record_read_ms = measure_ms("julia_genepop_record_read_ms", 10, () -> read_genepop_record(fixture_path))
    write_ms = measure_ms("julia_genepop_write_ms", 10, () -> write_genepop(parsed_populations, parsed_path; title="BioToolkit popgen benchmark"))
    allele_ms = measure_ms("julia_allele_frequency_ms", 200, () -> allele_frequencies(parsed_populations[1], 1))
    het_ms = measure_ms("julia_heterozygosity_ms", 200, () -> heterozygosity_expected(parsed_populations[1], 1))
    hwe_ms = measure_ms("julia_hwe_ms", 200, () -> hardy_weinberg_test(parsed_populations[1], 1))
    fst_ms = measure_ms("julia_f_statistics_ms", 100, () -> f_statistics(parsed_populations, 1))
    ld_ms = measure_ms("julia_ld_ms", 100, () -> linkage_disequilibrium(parsed_populations[1], 1, 2))
    pca_ms = measure_ms("julia_population_pca_ms", 20, () -> population_pca(parsed_populations, collect(1:locus_count)))
    split_pops_ms = measure_ms("julia_genepop_split_pops_ms", 25, () -> split_in_pops(parsed_record, ["A", "B", "C"]))
    split_loci_ms = measure_ms("julia_genepop_split_loci_ms", 25, () -> split_in_loci(parsed_record))
    remove_pop_ms = measure_ms("julia_genepop_remove_pop_ms", 25, () -> remove_population!(deepcopy(parsed_record), 0))
    remove_locus_ms = measure_ms("julia_genepop_remove_locus_ms", 25, () -> remove_locus_by_position!(deepcopy(parsed_record), 0))

    println("Julia PopGen benchmark")
    println("  populations=$population_count")
    println("  individuals_per_pop=$individual_count")
    println("  loci=$locus_count")
    println("  julia_genepop_read_ms=$(round(read_ms; digits=4))")
    println("  julia_genepop_record_read_ms=$(round(record_read_ms; digits=4))")
    println("  julia_genepop_write_ms=$(round(write_ms; digits=4))")
    println("  julia_allele_frequency_ms=$(round(allele_ms; digits=4))")
    println("  julia_heterozygosity_ms=$(round(het_ms; digits=4))")
    println("  julia_hwe_ms=$(round(hwe_ms; digits=4))")
    println("  julia_f_statistics_ms=$(round(fst_ms; digits=4))")
    println("  julia_ld_ms=$(round(ld_ms; digits=4))")
    println("  julia_population_pca_ms=$(round(pca_ms; digits=4))")
    println("  julia_genepop_split_pops_ms=$(round(split_pops_ms; digits=4))")
    println("  julia_genepop_split_loci_ms=$(round(split_loci_ms; digits=4))")
    println("  julia_genepop_remove_pop_ms=$(round(remove_pop_ms; digits=4))")
    println("  julia_genepop_remove_locus_ms=$(round(remove_locus_ms; digits=4))")
    println("  julia_stats_note=Biopython provides GenePop parsing/manipulation, not native equivalents for these analytical kernels")

    run(`conda run -n general python $(joinpath(@__DIR__, "popgen_compare.py"))`)
end

main()