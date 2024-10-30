using Setfield

@testset "ConsolidationColumn" begin

    cell = AtlansFixtures.draining_abc_isotache_cell()
    column = AtlansFixtures.draining_abc_isotache_column()

    @testset "constructor" begin
        @test typeof(column) == Atlans.ConsolidationColumn{
            Atlans.DrainingAbcIsotache,
            Atlans.OverConsolidationRatio,
        }
    end

    @testset "surface_level" begin
        @test Atlans.surface_level(column) ≈ 4.0
    end

    @testset "compress_γ_wet" begin
        actual = Atlans.compress_γ_wet(cell, 0.0)
        expected = 15000.0
        @test actual ≈ expected

        actual = Atlans.compress_γ_wet(cell, 0.1)
        expected = 15519.0
        @test actual ≈ expected
    end

    @testset "compress_γ_dry" begin
        actual = Atlans.compress_γ_dry(cell, 0.0)
        expected = 10000.0
        @test actual ≈ expected

        actual = Atlans.compress_γ_dry(cell, 0.01)
        expected = 10100.0
        @test actual ≈ expected
    end

    @testset "degree of consolidation" begin
        t = 0.01
        actual = Atlans.U(cell, t)
        expected = 0.004665987113375191
        @test actual ≈ expected

        t = 1.0e6
        actual = Atlans.U(cell, t)
        @test actual ≈ 1.0
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
        expected = [35000.0, 25000.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        # No effect on midweight!
        phreatic_level = 0.5
        Atlans.total_stress!(column, phreatic_level)
        expected = [35000.0, 25000.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        phreatic_level = 1.2
        Atlans.total_stress!(column, phreatic_level)
        expected = [38500.0, 25000.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        phreatic_level = 1.5
        Atlans.total_stress!(column, phreatic_level)
        expected = [40000.0, 25000.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        phreatic_level = 1.7
        Atlans.total_stress!(column, phreatic_level)
        expected = [41000.0, 26000.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        phreatic_level = 2.0
        Atlans.total_stress!(column, phreatic_level)
        expected = [42500.0, 27500.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        phreatic_level = 5.0
        Atlans.total_stress!(column, phreatic_level)
        expected = [62310.0, 47310.0, 32310.0, 17310.0]
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

    @testset "set_stress" begin
        # phreatic_level at 3.0 m; hydrostatic head.
        phreatic_level = 3.0
        column.p .= [2.5, 1.5, 0.5, 0.0] .* Atlans.γ_water

        Atlans.total_stress!(column, phreatic_level)
        Atlans.effective_stress!(column)
        # Transfer stress to cells
        Atlans.transfer_stress!(column)

        @test all([cell.σ′ for cell in column.cells] .== column.σ′)
    end

    @testset "consolidate" begin
        # Set the phreatic level at 3.0 m
        phreatic_level = 3.0
        Atlans.total_stress!(column, phreatic_level)
        column.p .= [2.5, 1.5, 0.5, 0.0] .* Atlans.γ_water
        Atlans.apply_preconsolidation!(column)

        # Now lower watertable by 0.2
        column.p .= [2.3, 1.3, 0.3, 0.0] .* Atlans.γ_water
        phreatic_level = 2.8
        Δt = 3650.0
        Atlans.consolidate!(column, phreatic_level, Δt)

        # the top [4] cell also experiences creep
        @test all(cell.consolidation > 0.0 for cell in column.cells)
        @test all(cell.Δz < 1.0 for cell in column.cells)
    end

    @testset "prepare_forcingperiod" begin
        Atlans.prepare_forcingperiod!(column)
        for cell in column.cells
            @test cell.t == 0.0
            @test cell.U == 0.0
            @test cell.Δz_0 == cell.Δz
        end
    end
end
