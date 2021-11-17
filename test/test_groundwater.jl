
@testset "SimpleGroundwater" begin
    nlayer = 8
    Δz = fill(0.25, nlayer)
    zbase = 0.0
    ztop = zbase .+ cumsum(Δz)
    z = ztop .- (0.5 .* Δz)
    boundary = Int[1, 8]
    boundary_ϕ = Float64[1.0, 0.0]
    dry = fill(false, nlayer)
    ϕ = fill(NaN, nlayer)
    p = fill(NaN, nlayer)
    γ_water = 9.81 * 1_000.0

    sg = Atlans.SimpleGroundwater(z, Δz, boundary, boundary_ϕ, dry, ϕ, p, γ_water)

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
        @test typeof(sg) == Atlans.SimpleGroundwater
    end

    @testset "interpolate_head" begin
        expected2 = reverse(expected1)

        Atlans.interpolate_head!(sg)
        @test all(sg.ϕ .≈ expected1)

        sg.boundary_ϕ[1] = 0.0
        sg.boundary_ϕ[2] = 1.0
        Atlans.interpolate_head!(sg)
        @test all(sg.ϕ .≈ expected2)

        sg.boundary_ϕ[1] = 1.0
        sg.boundary_ϕ[2] = 0.0
        sg.boundary[2] = 4
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
        Atlans.interpolate_head!(sg)
        @test all(sg.ϕ .≈ expected3)

        sg.boundary[1] = 5
        sg.boundary[2] = 8
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
        Atlans.interpolate_head!(sg)
        @test all(sg.ϕ .≈ expected3)
    end

    @testset "solve" begin
        sg.boundary[1] = 1
        sg.boundary[2] = 8
        sg.boundary_ϕ[1] = 1.0
        sg.boundary_ϕ[2] = 0.0
        Atlans.solve!(sg)
        @test all(sg.ϕ .≈ expected1)
    end

    @testset "pore_pressure" begin
        sg.ϕ .= 1.0
        Atlans.pore_pressure!(sg)
        @test all(sg.dry[1:4] .== false)
        @test all(sg.dry[5:8] .== true)
    end

    @testset "plot" begin
        Atlans.plot(sg)
    end
end
