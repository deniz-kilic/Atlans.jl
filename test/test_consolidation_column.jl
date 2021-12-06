@testset "ConsolidationColumn" begin

    Δz = 1.0
    t = 1.0
    σ′ = 10000.0
    γ_wet = 15000.0
    γ_dry = 10000.0
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

    cells = fill(cell, 4)
    z = collect(0.5:1.0:4.0)
    Δz = fill(1.0, 4)
    σ = fill(0.0, 4)
    σ′ = fill(0.0, 4)
    u = fill(0.0, 4)
    column = Atlans.ConsolidationColumn(cells, z, Δz, σ, σ′, u)

    @testset "ConsolidationColumn constructor" begin
        @test typeof(column) == Atlans.ConsolidationColumn
    end

    @testset "compress_γ_wet" begin
        cell = @set cell.consolidation = 0.1
        actual = Atlans.compress_γ_wet(cell)
        expected = 15576.666666666666
        @test actual ≈ expected
    end

    @testset "compress_γ_dry" begin
        cell = @set cell.consolidation = 0.01
        actual = Atlans.compress_γ_dry(cell)
        expected = 10101.0101010101
        @test actual ≈ expected
    end

    @testset "degree of consolidation" begin
        t = 0.01
        actual = Atlans.U(cell, t)
        expected = 5.7916576227323684e-15
        @test actual ≈ expected
    end

    @testset "cellweight" begin
        zbot = 0.0
        Δz = 1.0
        γ_wet = 10_000.0
        γ_dry = 5000.0

        @test Atlans.weight(2.0, zbot, Δz, γ_wet, γ_dry) ≈ 10_000.0
        @test Atlans.weight(1.0, zbot, Δz, γ_wet, γ_dry) ≈ 10_000.0
        @test Atlans.weight(0.5, zbot, Δz, γ_wet, γ_dry) ≈ 7500.0
        @test Atlans.weight(0.25, zbot, Δz, γ_wet, γ_dry) ≈ 6250.0
        @test Atlans.weight(0.0, zbot, Δz, γ_wet, γ_dry) ≈ 5000.0
        @test Atlans.weight(-1.0, zbot, Δz, γ_wet, γ_dry) ≈ 5000.0
    end

    @testset "submerged_weight" begin
        # Water column of one meter on top
        phreatic_level = 5.0
        actual = Atlans.submerged_weight(column, phreatic_level)
        @test actual ≈ 9810.0

        # Below column top: no water column weight on top
        phreatic_level = 2.5
        actual = Atlans.submerged_weight(column, phreatic_level)
        @test actual ≈ 0.0
    end

    @testset "totalstress" begin
        phreatic_level = 0.0
        Atlans.total_stress!(column, phreatic_level)
        expected = Float64[35000.0, 25000.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        # No effect on midweight!
        phreatic_level = 0.5
        Atlans.total_stress!(column, phreatic_level)
        expected = Float64[35000.0, 25000.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        phreatic_level = 1.5
        Atlans.total_stress!(column, phreatic_level)
        expected = Float64[42500.0, 25000.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        phreatic_level = 2.0
        Atlans.total_stress!(column, phreatic_level)
        expected = Float64[45000.0, 30000.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        phreatic_level = 5.0
        Atlans.total_stress!(column, phreatic_level)
        expected = Float64[62310.0, 47310.0, 32310.0, 17310.0]
        @test all(column.σ .≈ expected)
    end

    @testset "effective_stress" begin
        column.σ .= 10_000.0
        column.p .= 5_000.0
        Atlans.effective_stress!(column)
        @test all(column.σ′ .≈ 5000.0)

        column.σ[4] = 2_000.0
        Atlans.effective_stress!(column)
        @test column.σ′[4] == 0.0
    end

    @testset "plot" begin
        Atlans.plot(column)
    end
end
