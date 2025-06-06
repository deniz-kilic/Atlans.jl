@testset "DrainingAbcIsotache" begin

    Δz = 1.0
    t = 1.0
    σ′ = 10000.0
    γ_wet = 15000.0
    γ_dry = 15000.0
    c_d = 2.0
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

    @testset "consolidate" begin
        σ′ = 10000.0
        Δt = 1.0
        new_cell = Atlans.consolidate(cell, σ′, Δt)
        @test new_cell.consolidation > 0.0
        @test typeof(new_cell) == Atlans.DrainingAbcIsotache
        @test isapprox(new_cell.consolidation, 0.006, atol = 0.0001)
    end

end
