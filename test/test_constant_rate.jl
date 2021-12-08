@testset "ConstantRate" begin
    cell = Atlans.ConstantRate(
        1.0,
        1.0e-3,
    )
    
    @testset "Initialization" begin
        @test typeof(cell) == Atlans.ConstantRate
    end

    @testset "oxidate" begin
        Δt = 1.0
        actual, new_cell = Atlans.oxidate(cell, Δt)
        expected = 0.000999500166624978
        # @test actual == expected
        @test actual ≈ expected
        @test typeof(new_cell) == Atlans.ConstantRate
    end
end