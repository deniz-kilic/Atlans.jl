@testset "CarbonStore" begin
    cell = CarbonStore(
        1.0,
        0.2,
        mass_organic(0.2, 1000, 1),
        mass_mineral(0.2, 1000, 1),
        mass_organic_minimal(mass_mineral(0.2, 1000, 1), 0.05),
        1.0e-3,
        1000,
    )
    
    @testset "Initialization" begin
        @test typeof(cell) == CarbonStore
    end

    @testset "oxidate" begin
        Δt = 1.0
        actual, new_cell = oxidate(cell, Δt)
        expected = 1.0e-7
        # @test actual == expected
        @test actual ≈ expected
        @test typeof(new_cell) == CarbonStore
    end


end