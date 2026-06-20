include(joinpath(@__DIR__, "..", "Examples", "Comparison", "projection_consistency.jl"))

@testset "Projection consistency benchmark" begin
    temp_plot = mktempdir() do tempdir
        plot_path = joinpath(tempdir, "projection_consistency.png")
        report = projection_consistency_report(save_plot=true, plot_path=plot_path)
        @test isfile(plot_path)
        @test report.pc1_correlation > 0.99
        @test report.pc2_correlation > 0.99
        @test report.rmse < 1.0
        report
    end

    @test size(temp_plot.projected_query) == size(temp_plot.true_query)
end