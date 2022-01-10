@testset "DrainingAbcIsotache" begin

    Δz = 1.0
    t = 1.0
    σ′ = 10000.0
    γ_wet = 15000.0
    γ_dry = 15000.0
    c_d = 2
    c_v = 0.006912
    U = 1.0
    a = 0.01737
    b = 0.1303
    c = 0.008686
    τ = 1.0
    consolidation = 0.0

    cell = Atlans.DrainingAbcIsotache(
        Δz,
        Δz,
        t,
        σ′,
        γ_wet,
        γ_dry,
        c_d,
        c_v,
        U,
        a,
        b,
        c,
        τ,
        consolidation,
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


    @testset "consolidate" begin
        σ′ = 10000.0
        Δt = 1.0
        new_cell = Atlans.consolidate(cell, σ′, Δt)
        expected = 0.0
        @test new_cell.consolidation ≈ expected
        @test typeof(new_cell) == Atlans.DrainingAbcIsotache
    end

end
