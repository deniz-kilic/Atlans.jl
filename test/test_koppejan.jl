@testset "DrainingKoppejan" begin
    cell = Atlans.DrainingKoppejan(
        1.0,
        1.0,
        10_000.0,
        15_000.0,
        15_000.0,
        2,
        1.0e-5,
        1.753846950483396e-17,
        20.0,
        80.0,
        5.0,
        20.0,
        10_000.0,
    )
    
    @testset "Initialization" begin
        @test typeof(cell) == Atlans.DrainingKoppejan
    end

    @testset "compress_γ_wet" begin
        consolidation = 0.1
        actual = Atlans.compress_γ_wet(cell, consolidation)
        expected = 8720.0
        @test actual ≈ expected
    end

    @testset "compress_γ_dry" begin
        consolidation = 0.01
        actual = Atlans.compress_γ_dry(cell, consolidation)
        expected = 15151.515151515152
        @test actual ≈ expected
    end

    @testset "degree of consolidation" begin
        t = 0.01
        actual = Atlans.U(cell, t)
        expected = 0
        @test actual ≈ expected atol = 0.001
    end

    @testset "consolidate" begin
        σ′ = 10_000.0
        Δt = 1.0
        actual, new_cell = Atlans.consolidate(cell, σ′, Δt)
        expected = 0
        @test actual ≈ expected atol=1e-6
        @test typeof(new_cell) == Atlans.DrainingKoppejan
    end

end