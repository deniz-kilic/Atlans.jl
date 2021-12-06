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
