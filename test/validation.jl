using Test
using Random
using Statistics
using LinearAlgebra
using Distributions
using BioToolkit

function synthetic_count_matrix(; genes::Int=240, samples::Int=12, spike_genes::Int=6, base_rate::Int=120, fold_change::Int=4, dispersion_size::Int=25, seed::Int=20260325)
    rng = MersenneTwister(seed)
    control = fill(:control, samples ÷ 2)
    treated = fill(:treated, samples ÷ 2)
    design = vcat(control, treated)

    counts = Matrix{Int}(undef, genes, samples)
    for gene in 1:genes
        for sample in 1:samples
            rate = gene <= spike_genes && sample > samples ÷ 2 ? base_rate * fold_change : base_rate
            probability = dispersion_size / (dispersion_size + rate)
            counts[gene, sample] = rand(rng, NegativeBinomial(dispersion_size, probability))
        end
    end

    gene_ids = ["gene_$(index)" for index in 1:genes]
    sample_ids = ["sample_$(index)" for index in 1:samples]
    return BioToolkit.CountMatrix(counts, gene_ids, sample_ids), design
end

function ks_distance(values::AbstractVector{<:Real})
    sorted = sort(Float64.(values))
    n = length(sorted)
    n == 0 && return 0.0

    d_plus = maximum((index / n) - sorted[index] for index in 1:n)
    d_minus = maximum(sorted[index] - ((index - 1) / n) for index in 1:n)
    return max(d_plus, d_minus)
end

@testset "Synthetic validation" begin
    count_matrix, design = synthetic_count_matrix()
    results = BioToolkit.differential_expression(count_matrix, design; min_total=0, shrink=false)

    spiked = results[1:6]
    nulls = results[7:end]

    @test all(abs(result.log2_fold_change - 2.0) < 0.4 for result in spiked)
    @test abs(mean(result.log2_fold_change for result in spiked) - 2.0) < 0.3

    null_pvalues = [result.pvalue for result in nulls]
    @test ks_distance(null_pvalues) < 0.2
    @test 0.35 < mean(null_pvalues) < 0.55
    @test 0.08 < count(<(0.1), null_pvalues) / length(null_pvalues) < 0.24
end

@testset "Alignment baseline" begin
    exact = BioToolkit.pairwise_align("ACGT", "ACGT"; match=2, mismatch=-1, gap=-2)
    @test exact.score == 8
    @test exact.left == "ACGT"
    @test exact.right == "ACGT"

    simple_indel = BioToolkit.pairwise_align("ACGT", "ACGGT"; match=1, mismatch=-1, gap=-1)
    @test simple_indel.score == 3
    @test length(simple_indel.left) == length(simple_indel.right)
end

@testset "Kabsch baseline" begin
    reference = [
        0.0 0.0 0.0;
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0;
    ]
    rotation = [
        0.0 -1.0 0.0;
        1.0 0.0 0.0;
        0.0 0.0 1.0;
    ]
    translation = [3.0, -2.0, 5.0]
    mobile = reference * rotation .+ translation'

    result = BioToolkit.superpose(reference, mobile)
    aligned = mobile * result.rotation .+ result.translation'

    @test result.rmsd < 1e-10
    @test maximum(abs.(aligned .- reference)) < 1e-10
    @test BioToolkit.rmsd(reference, mobile) < 1e-10
end