
@testset "DrainingAbc_InterpolatedGroundwater" begin
    # eight cells of 0.25 m

    nlayer = 8
    Δz = fill(0.25, nlayer)
    zbase = 0.0
    ztop = zbase .+ cumsum(Δz)
    z = ztop .- (0.5 .* Δz)

    boundary = Int[1, 8]
    boundary_ϕ = Float64[1.0, 0.0]
    γ_water = 9.81 * 1000.0

    σ′ = fill(0, nlayer)
    γ_wet = fill(18_000.0, nlayer)
    γ_dry = fill(13_000.0, nlayer)
    c_d = fill(2, nlayer)
    c_v = fill(1.0e-5, nlayer)
    a = fill(1.0e-6, nlayer)
    b = fill(1.0e-6, nlayer)
    c = fill(1.0e-6, nlayer)

    groundwater = Atlans.InterpolatedGroundwater(z, Δz, boundary, boundary_ϕ, γ_water)
    consolidation = Atlans.draining_abc_isotache_column(Δz, γ_wet, γ_dry, c_d, c_v, a, b, c)
    oxidation = Vector{Atlans.NullOxidation}(undef, nlayer)

    column = Atlans.SoilColumn(
        0.0,  # x
        0.0,  # y
        z,
        Δz,
        consolidation,
        oxidation,
        groundwater,
    )

    @testset "Constructor" begin
        @test typeof(column) == Atlans.SoilColumn
    end
end
