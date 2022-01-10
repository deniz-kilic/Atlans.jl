@testset "ConstantRate" begin
    Δz = 1.0
    α = 1.0e-3
    cell = Atlans.ConstantRate(Δz, α, 0.0)

    @testset "Initialization" begin
        @test typeof(cell) == Atlans.ConstantRate
    end

    @testset "oxidate" begin
        Δt = 1.0
        new_cell = Atlans.oxidate(cell, Δt)
        expected = 0.000999500166624978
        @test new_cell.oxidation ≈ expected
        @test typeof(new_cell) == Atlans.ConstantRate
    end
end
