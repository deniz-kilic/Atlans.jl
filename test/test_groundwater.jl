@testset "InterpolatedGroundwater" begin
    nlayer = 8
    Δz = fill(0.25, nlayer)
    zbase = 0.0
    ztop = zbase .+ cumsum(Δz)
    z = ztop .- (0.5 .* Δz)
    boundary = Int[1, 8]
    boundary_ϕ = Float64[1.0, 0.0]
    confined = fill(true, nlayer)
    dry = fill(false, nlayer)
    ϕ = fill(NaN, nlayer)
    p = fill(NaN, nlayer)

    ig = Atlans.InterpolatedGroundwater(z, Δz, boundary, boundary_ϕ, confined, dry, ϕ, p)

    expected1 = Float64[
        1.0,
        0.8571428571428572,
        0.7142857142857143,
        0.5714285714285714,
        0.4285714285714286,
        0.2857142857142857,
        0.1428571428571429,
        0.0,
    ]

    @testset "Constructor" begin
        @test typeof(ig) == Atlans.InterpolatedGroundwater
    end

    @testset "interpolate_head" begin
        expected2 = reverse(expected1)

        Atlans.interpolate_head!(ig)
        @test all(ig.ϕ .≈ expected1)

        ig.boundary_ϕ[1] = 0.0
        ig.boundary_ϕ[2] = 1.0
        Atlans.interpolate_head!(ig)
        @test all(ig.ϕ .≈ expected2)

        ig.boundary_ϕ[1] = 1.0
        ig.boundary_ϕ[2] = 0.0
        ig.boundary[2] = 4
        expected3 = Float64[
            1.0
            0.6666666666666667
            0.33333333333333337
            0.0
            0.0
            0.0
            0.0
            0.0
        ]
        Atlans.interpolate_head!(ig)
        @test all(ig.ϕ .≈ expected3)

        ig.boundary[1] = 5
        ig.boundary[2] = 8
        expected3 = Float64[
            1.0
            1.0
            1.0
            1.0
            1.0
            0.6666666666666667
            0.33333333333333337
            0.0
        ]
        Atlans.interpolate_head!(ig)
        @test all(ig.ϕ .≈ expected3)
    end

    @testset "flow" begin
        ig.boundary[1] = 1
        ig.boundary[2] = 8
        ig.boundary_ϕ[1] = 1.0
        ig.boundary_ϕ[2] = 0.0
        Atlans.flow!(ig, 0.0)
        @test all(ig.ϕ .≈ expected1)
    end

    @testset "pore_pressure" begin
        ig.ϕ .= 1.0
        Atlans.pore_pressure!(ig)
        @test all(ig.dry[1:4] .== false)
        @test all(ig.dry[5:8] .== true)
    end

    @testset "plot" begin
        Atlans.plot(ig)
    end
end

@testset "Darcy" begin
    using LinearAlgebra

    nlayer = 8
    Δz = fill(0.25, nlayer)
    k = fill(1.0, nlayer)
    SS = fill(0.01, nlayer)
    zbase = 0.0
    ztop = zbase .+ cumsum(Δz)
    z = ztop .- (0.5 .* Δz)
    boundary = Int[1, 8]
    boundary_ϕ = Float64[1.0, 0.0]
    confined = fill(true, nlayer)
    dry = fill(false, nlayer)
    ϕ = fill(NaN, nlayer)
    ϕ_old = fill(NaN, nlayer)
    p = fill(NaN, nlayer)
    conductance = fill(NaN, nlayer - 1)
    dv = fill(NaN, nlayer)
    ev = fill(NaN, nlayer - 1)
    A = SymTridiagonal(dv, ev)
    rhs = fill(0.0, nlayer)

    dc = Atlans.DarcyColumn(
        k,
        z,
        Δz,
        SS,
        boundary,
        boundary_ϕ,
        confined,
        dry,
        conductance,
        ϕ,
        ϕ_old,
        p,
        A,
        rhs,
    )

    boundary = [1, 8]
    ϕ = [2.0, 1.0]

    expected1 = Float64[
        2.0,
        1.8571428571428572,
        1.7142857142857143,
        1.5714285714285714,
        1.4285714285714286,
        1.2857142857142857,
        1.1428571428571429,
        1.0,
    ]

    @testset "Constructor" begin
        @test typeof(dc) == Atlans.DarcyColumn
    end

    @testset "harmonicmean" begin
        @test Atlans.harmonicmean_conductance(1.0, 1.0, 1.0, 1.0) ≈ 1.0
        @test Atlans.harmonicmean_conductance(1.0, 1.0, 0.5, 1.0) ≈ 4.0 / 3.0
        @test Atlans.harmonicmean_conductance(1.0, 0.5, 1.0, 1.0) ≈ 2.0 / 3.0
    end

    @testset "conductance" begin
        Atlans.conductance!(dc)
        @test all(dc.conductance .≈ 4.0)
    end

    @testset "update boundaries" begin
        Atlans.update_boundaries!(dc, boundary, ϕ)
        @test dc.boundary_ϕ == [2.0, 1.0]
        @test dc.boundary == [1, 8]
    end

    @testset "formulate" begin
        # Compute a steady-state: should linearly interpolate
        Δt = 0.0
        Atlans.formulate!(dc, Δt)
        @test dc.A.dv ≈ [1.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, 1.0]
        @test dc.A.ev ≈ [0.0, 4.0, 4.0, 4.0, 4.0, 4.0, 0.0]
        @test dc.rhs ≈ [2.0, -8.0, 0.0, 0.0, 0.0, 0.0, -4.0, 1.0]
    end

    @testset "linearsolve" begin
        Atlans.linearsolve!(dc)
        @test all(dc.ϕ .≈ expected1)
    end

    @testset "steady flow" begin
        dc.A.dv .= 0.0
        dc.A.ev .= 0.0
        dc.rhs .= 0.0
        Atlans.flow!(dc, 0.0)
        @test all(dc.ϕ .≈ expected1)
    end

    @testset "transient flow" begin
        dc.ϕ .= 2.0
        Atlans.update_boundaries!(dc, [1, 8], [0.0, 0.0])

        Atlans.flow!(dc, 0.01)
        ϕ1 = copy(dc.ϕ)
        # Test for symmetry
        @test ϕ1[1:4] ≈ reverse(ϕ1[5:8])

        Atlans.flow!(dc, 0.01)
        # Head should decrease over time
        ϕ2 = copy(dc.ϕ)
        @test all(ϕ2 .<= ϕ1)
    end
end
