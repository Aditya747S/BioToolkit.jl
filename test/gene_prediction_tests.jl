using Test
using BioToolkit

@testset "Applied HMM Gene Prediction" begin
    @testset "predict_genes_hmm pipeline" begin
        # 42bp Random noise (GC 25%), then a 60bp "exon" (GC > 60%), then 42bp noise.
        # Background NonCoding states should dominate the edges.
        
        noise_left  = "ATATAAATTTTAAATATATATAAATTTTAAATATATATAAT"
        exon_block  = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" # High GC coding signature
        noise_right = "TAATTTTAATAATATATATAAATTTTAATAATATATATAAT"
        
        seq = noise_left * exon_block * noise_right
        
        genes = predict_genes_hmm(seq; p_coding=0.05)
        
        # It should cleanly isolate the single `exon_block` as a Tuple bounds
        @test length(genes) == 1
        start, stop = genes[1]
        
        # Since Viterbi boundaries aren't an absolute hard knife cut, 
        # mathematically they'll be very close to the exon_block region (~ idx 42)
        @test start > 20
        @test stop < length(seq) - 20
        @test (stop - start) > 20
    end
    
    @testset "Pure Noise Filtering" begin
        # A sequence completely devoid of GC chunks
        seq = "ATATAAATTTTAAATATATATAAATTTTAAATATATATAAT"
        
        genes = predict_genes_hmm(seq; p_coding=1e-5)
        # Should not find any genes as it stays in State 1
        @test isempty(genes)
    end
    
    @testset "Edge Cases" begin
        # Blank sequence
        genes = predict_genes_hmm("")
        @test isempty(genes)
    end
end
