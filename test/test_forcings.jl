@testset "Forcings" begin
    function test_forcings()
        si = Atlans.StageIndexation(AtlansFixtures.stage_indexation_netcdf(), 50)
        sc = Atlans.StageChange(AtlansFixtures.stage_change_netcdf())
        ds = Atlans.DeepSubsidence(AtlansFixtures.deep_subsidence_netcdf())
        t = Atlans.Temperature(AtlansFixtures.temperature_table())
        sur = Atlans.Surcharge(
            AtlansFixtures.simple_surcharge_netcdf(),
            AtlansFixtures.params_table()
        )
        return Atlans.Forcings(
            stage_change=sc,
            stage_indexation=si,
            deep_subsidence=ds,
            temperature=t,
            surcharge=sur
        )
    end

    function testing_model()
        path_csv = AtlansFixtures.params_table()
        path_nc = AtlansFixtures.subsoil_netcdf()

        return Atlans.Model(
            Atlans.HydrostaticGroundwater,
            Atlans.DrainingAbcIsotache,
            Atlans.CarbonStore,
            Atlans.OverConsolidationRatio,
            Atlans.SimpleShrinkage,
            Atlans.AdaptiveCellsize(0.25, 0.01),
            Atlans.ExponentialTimeStepper(1.0, 2),
            path_nc,
            path_csv,
        )
    end

    test_subsidence = [0.05 0.07 0.02; 0.04 0.03 0.05]

    @testset "stage indexation" begin
        model = testing_model()
        forcings = test_forcings()

        model.output.subsidence .= test_subsidence

        f = Atlans.load_forcing!(forcings, :stage_indexation, DateTime("2020-01-01"), model)
        @test all(f.change .== -0.045)

        idx = model.index[1]
        col = model.columns[1]
        col_change = Atlans.get_elevation_shift(f, col, idx)

        @test col_change == -0.045

        @test Atlans.phreatic_level(col.groundwater) == 0.5
        Atlans.apply_forcing!(f, col, idx)
        @test Atlans.phreatic_level(col.groundwater) ≈ 0.4775
    end

    @testset "temperature" begin
        model = testing_model()
        forcings = test_forcings()

        f = Atlans.load_forcing!(forcings, :temperature, DateTime("2020-01-01"), model)

        @test f.temp == 14.0

        col = model.columns[1]
        idx = model.index[1]
        Atlans.apply_forcing!(f, col, idx)
        new_alphas = [c.α for c in col.oxidation.cells]

        @test all(new_alphas .≈ 0.00124611669)
    end

    @testset "deep subsidence" begin
        model = testing_model()
        forcings = test_forcings()

        f = Atlans.load_forcing!(forcings, :deep_subsidence, DateTime("2020-01-01"), model)
        @test all(f.subsidence .≈ 0.05)

        col = model.columns[1]
        idx = model.index[1]

        Atlans.apply_forcing!(f, col, idx)
        @test col.base.z ≈ -0.05
    end

    @testset "stage change" begin
        model = testing_model()
        forcings = test_forcings()

        f = Atlans.load_forcing!(forcings, :stage_change, DateTime("2020-01-01"), model)
        @test all(f.change .≈ -0.1)

        col = model.columns[1]
        idx = model.index[1]

        @test Atlans.phreatic_level(col.groundwater) == 0.5
        Atlans.apply_forcing!(f, col, idx)
        @test Atlans.phreatic_level(col.groundwater) ≈ 0.4
    end

    @testset "Surcharge" begin
        model = testing_model()
        forcings = test_forcings()

        f = Atlans.load_forcing!(forcings, :surcharge, DateTime("2020-01-01"), model)
        @test all(f.lithology .== 2)
        @test all(f.thickness .== 0.5)

        for i in eachindex(model.index)
            idx = model.index[i]
            col = model.columns[i]
            ncells = length(col.z)

            Atlans.apply_forcing!(f, col, idx)
            new_ncells = length(col.z)
            @test new_ncells == ncells + 2
        end
    end
end