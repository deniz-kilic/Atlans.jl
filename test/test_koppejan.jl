@testset "DrainingKoppejan" begin

    Δz = 1.0
    Δz_0 = 0.0
    t = 1.0
    σ′ = 10_000.0
    γ_wet = 15_000.0
    γ_dry = 15_000.0
    c_d = 2
    c_v = 1.0e-5
    U = 1.753846950483396e-17
    Cp = 20.0
    Cs = 80.0
    Cp′ = 5.0
    Cs′ = 20.0
    σ′pre = 10_000.0
    consolidation = 0.0

    cell = Atlans.DrainingKoppejan(
        Δz,
        Δz_0,
        t,
        σ′,
        γ_wet,
        γ_dry,
        c_d,
        c_v,
        U,
        Cp,
        Cs,
        Cp′,
        Cs′,
        σ′pre,
        consolidation,
    )

    @testset "Initialization" begin
        @test typeof(cell) == Atlans.DrainingKoppejan
    end

    @testset "compress_γ_wet" begin
        consolidation = 0.1
        actual = Atlans.compress_γ_wet(cell, consolidation)
        expected = 13333.33333333333
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
        new_cell = Atlans.consolidate(cell, σ′, Δt)
        expected = 0
        @test new_cell.consolidation ≈ expected atol = 1e-6
        @test typeof(new_cell) == Atlans.DrainingKoppejan
    end

end
