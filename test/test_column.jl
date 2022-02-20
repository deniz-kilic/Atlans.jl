@testset "SoilColumn" begin
    
    @testset "constructor" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        @test typeof(column) == Atlans.SoilColumn{
            Atlans.HydrostaticGroundwater,
            Atlans.DrainingAbcIsotache,
            Atlans.OverConsolidationRatio,
            Atlans.CarbonStore,
        }
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
 
        Atlans.prepare_forcingperiod!(column)
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
        Atlans.prepare_forcingperiod!(column)
        # Set a single timestep of 1 day.
        subsidence, consolidation, oxidation = Atlans.advance_timestep!(column, 1.0)

        @test consolidation > 0.0
        @test oxidation > 0.0
        @test subsidence > 0.0
        @test subsidence ≈ (consolidation + oxidation)
    end

    @testset "Compute: lower phreatic" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        # Set initial τ
        Atlans.apply_preconsolidation!(column)
        # Prepare forcing period: store pre-load stress
        Atlans.prepare_forcingperiod!(column)
        # Set a forcing: lower the water table by 0.2 m
        Atlans.set_phreatic!(column, 2.8)
        # Generate exponentially increasing timesteps
        ts = Atlans.ExponentialTimeStepper(1.0, 2)
        timesteps = Atlans.create_timesteps(ts, 3650.0)
        # Run through the timesteps
        subsidence, consolidation, oxidation = Atlans.advance_forcingperiod!(column, timesteps)
        
        @test consolidation > 0.0
        @test oxidation > 0.0
        @test subsidence > 0.0
        @test subsidence ≈ (consolidation + oxidation)
        
        # Test for object identity (should all point to same array in memory)
        @test column.z === column.consolidation.z === column.oxidation.z === column.groundwater.z
        @test column.Δz === column.consolidation.Δz === column.oxidation.Δz
        
        @show subsidence
    end
    
    @testset "Compute: raise base" begin
        # This should give the same result as lowering phreatic level by 0.2 m.
        column = AtlansFixtures.soil_column_hg_abc_cs()
        # Set initial τ
        Atlans.apply_preconsolidation!(column)
        # Prepare forcing period: store pre-load stress
        Atlans.prepare_forcingperiod!(column)
        # Set a forcing: raise the column by 0.2 m.
        Atlans.set_deep_subsidence!(column, -0.2)

        # Generate exponentially increasing timesteps
        ts = Atlans.ExponentialTimeStepper(1.0, 2)
        timesteps = Atlans.create_timesteps(ts, 3650.0)
        # Run through the timesteps
        subsidence, consolidation, oxidation = Atlans.advance_forcingperiod!(column, timesteps)
    end

end