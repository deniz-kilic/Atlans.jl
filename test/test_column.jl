@testset "SoilColumn" begin

    function testing_column()
        column = AtlansFixtures.soil_column_hg_abc_cs()
        # Set initial τ
        Atlans.apply_preconsolidation!(column)
        # Prepare forcing period: store pre-load stress
        Atlans.prepare_forcingperiod!(column, 10.0)
        return column
    end

    function advance(column)
        # Generate exponentially increasing timesteps
        ts = Atlans.ExponentialTimeStepper(1.0, 2)
        timesteps = Atlans.create_timesteps(ts, 3650.0)
        # Run through the timesteps
        return Atlans.advance_forcingperiod!(column, timesteps)
    end

    function advance_set_phreatic()
        column = testing_column()
        Atlans.set_phreatic!(column, 2.8)
        return advance(column)
    end

    function advance_diff_phreatic()
        column = testing_column()
        Atlans.set_phreatic_difference!(column, -0.2)
        return advance(column)
    end

    function advance_deep_subsidence()
        column = testing_column()
        # Note: subsidence is downward, raising requires a negative value.
        Atlans.set_deep_subsidence!(column, -0.2)
        return advance(column)
    end

    @testset "constructor" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        @test typeof(column) == Atlans.SoilColumn{
            Atlans.HydrostaticGroundwater,
            Atlans.DrainingAbcIsotache,
            Atlans.OverConsolidationRatio,
            Atlans.CarbonStore,
        }
    end

    @testset "surface_level" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        @test Atlans.surface_level(column) ≈ 4.0
    end

    function test_initial_stress(column)
        # flow! (changes nothing in hydrostatic column)
        @test Atlans.phreatic_level(column.groundwater) == 3.0
        # pore_pressure!
        @test all(column.groundwater.p .≈ [2.5, 1.5, 0.5, 0.0])
        cc = column.consolidation
        @test all(cc.p .≈ (column.groundwater.p * Atlans.γ_water))
        @test all(cc.σ .≈ [47500.0, 32500.0, 17500.0, 5000.0])
        @test all(cc.σ′ .≈ (cc.σ - cc.p))
        return
    end


    @testset "initial_stress" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        cc = column.consolidation

        Atlans.initial_stress!(column)
        test_initial_stress(column)
        @test all(cc.σ′ .≈ collect(cell.σ′ for cell in cc.cells))
    end

    @testset "apply_preconsolidation" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        cc = column.consolidation

        @test all((cell.τ for cell in cc.cells) .== 1.0)
        Atlans.apply_preconsolidation!(column)
        test_initial_stress(column)
        @test all(cc.σ′ .≈ collect(cell.σ′ for cell in cc.cells))
        @test all((cell.τ for cell in cc.cells) .== 20997.0969628065)
    end

    @testset "prepare_forcingperiod" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        cc = column.consolidation

        Atlans.prepare_forcingperiod!(column, 10.0)
        test_initial_stress(column)
        @test all(cc.σ′ .≈ collect(cell.σ′ for cell in cc.cells))
        @test all((cell.t for cell in cc.cells) .== 0.0)
        @test all((cell.U for cell in cc.cells) .== 0.0)
    end

    @testset "Compute: no load" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        # Set initial τ
        Atlans.apply_preconsolidation!(column)
        # Prepare forcing period: store pre-load stress
        Atlans.prepare_forcingperiod!(column, 10.0)
        # Set a single timestep of 1 day.
        subsidence, consolidation, oxidation = Atlans.advance_timestep!(column, 1.0)

        @test consolidation > 0.0
        @test oxidation > 0.0
        @test subsidence > 0.0
        @test subsidence ≈ (consolidation + oxidation)
    end

    @testset "Compute: lower phreatic" begin
        column = testing_column()
        Atlans.set_phreatic!(column, 2.8)
        subsidence, consolidation, oxidation = advance(column)

        @test consolidation > 0.0
        @test oxidation > 0.0
        @test subsidence > 0.0
        @test subsidence ≈ (consolidation + oxidation)

        # Test for object identity (should all point to same array in memory)
        @test column.z ===
              column.consolidation.z ===
              column.oxidation.z ===
              column.groundwater.z
        @test column.Δz === column.consolidation.Δz === column.oxidation.Δz
    end

    @testset "Compute: different forcing, same output" begin
        s1, c1, o1 = advance_set_phreatic()
        s2, c2, o2 = advance_diff_phreatic()
        s3, c3, o3 = advance_deep_subsidence()

        @test c1 ≈ c2 ≈ c3  # same consolidation
        @test o1 ≈ o2 ≈ o3  # same oxidation
        @test s1 ≈ s2 ≈ s3  # same subsidence
    end

    @testset "NullOxidation" begin
        column = AtlansFixtures.soil_column_hg_abc_null()
        # Set initial τ
        Atlans.apply_preconsolidation!(column)
        # Prepare forcing period: store pre-load stress
        Atlans.prepare_forcingperiod!(column, 10.0)

        Atlans.set_phreatic!(column, 2.8)
        subsidence, consolidation, oxidation = advance(column)

        @test subsidence ≈ consolidation
        @test oxidation == 0.0
    end

    @testset "NullConsolidation" begin
        column = AtlansFixtures.soil_column_hg_null_cs()
        Atlans.prepare_forcingperiod!(column, 10.0)

        Atlans.set_phreatic!(column, 2.8)
        subsidence, consolidation, oxidation = advance(column)

        @test subsidence ≈ oxidation
        @test consolidation == 0.0
    end
end
