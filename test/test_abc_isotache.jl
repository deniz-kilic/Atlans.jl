@testset "DrainingAbcIsotache" begin
    cell = Atlans.DrainingAbcIsotache(
        1.0,
        1.0,
        10000.0,
        15000.0,
        15000.0,
        2,
        0.006912,
        1.0,
        0.01737,
        0.1303,
        0.008686,
        1,#\tau
        0.0,
    )

    @testset "Initialization" begin
        @test typeof(cell) == Atlans.DrainingAbcIsotache
    end

    @testset "τ_intrinsic" begin
        loadstep = 1.0
        actual = Atlans.τ_intermediate(cell, loadstep)
        expected = 0.9987006417333234
        @test actual ≈ expected
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
        expected = 5.7916576227323684e-15
        @test actual ≈ expected
    end

    @testset "consolidate" begin
        σ′ = 10000.0
        Δt = 1.0
        actual, new_cell = Atlans.consolidate(cell, σ′, Δt)
        expected = 0
        @test actual ≈ expected
        @test typeof(new_cell) == Atlans.DrainingAbcIsotache
    end

end
