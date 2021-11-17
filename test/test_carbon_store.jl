@testset "CarbonStore" begin
    cell = Atlans.CarbonStore(
        1.0,
        0.2,
        Atlans.mass_organic(0.2, 1000, 1),
        Atlans.mass_mineral(0.2, 1000, 1),
        Atlans.mass_organic_minimal(Atlans.mass_mineral(0.2, 1000, 1), 0.05),
        1.0e-3,
        1000,
    )
    
    @testset "Initialization" begin
        @test typeof(cell) == Atlans.CarbonStore
    end

    @testset "oxidate" begin
        Δt = 1.0
        actual, new_cell = Atlans.oxidate(cell, Δt)
        expected = 1.0e-7
        # @test actual == expected
        @test actual ≈ expected
        @test typeof(new_cell) == Atlans.CarbonStore
    end
end