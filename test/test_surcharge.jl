@testset "Surcharge" begin
    Δzmax = 0.25
    I = CartesianIndex(1, 1)
    
    function get_surcharge_forcing()
        model = AtlansFixtures.testing_model()
        f = AtlansFixtures.test_forcings()
        f = Atlans.load_forcing!(f, :surcharge, DateTime("2020-01-01"), model)
    end

    function assert_eltype(testcol, col)
        testcol_type = eltype(testcol.cells)
        col_type = eltype(col.cells)
        testcol_type == col_type
    end

    @testset "prepare_domain" begin    
        domain = Atlans.prepare_domain(0.5, 2, Δzmax)
        @test domain.n == 2
        @test all(domain.Δz .== Δzmax)
        @test all(isnan.(domain.z))
        @test all(domain.lithology .== 2)

        # Test when more complex profile is added for Surcharge
        domain = Atlans.prepare_domain([0.1, 0.2, 0.4, 0.1], [2, 2, 3, 2], Δzmax)
        @test domain.n == 5
        @test all(domain.Δz .== [0.1, 0.2, 0.2, 0.2, 0.1])
        @test all(isnan.(domain.z))
        @test all(domain.lithology .== [2, 2, 3, 3, 2])

    end

    @testset "Test corresponding depths" begin
        col = AtlansFixtures.soil_column_hg_abc_cs_shr()
        f = get_surcharge_forcing()

        surcol = Atlans.prepare_surcharge_column(f, col, I)

        z = surcol.z
        Δz = surcol.Δz

        @test all(z .== surcol.groundwater.z)

        @test all(z .== surcol.consolidation.z)
        @test all(Δz .== surcol.consolidation.Δz)

        @test all(z .== surcol.oxidation.z)
        @test all(Δz .== surcol.oxidation.Δz)

        @test all(z .== surcol.shrinkage.z)
        @test all(Δz .== surcol.shrinkage.Δz)
    end
    
    @testset "prepare_surcharge_column_all_processes" begin
        col = AtlansFixtures.soil_column_hg_abc_cs_shr()
        f = get_surcharge_forcing()

        surcol = Atlans.prepare_surcharge_column(f, col, I)
        
        surcol_preconsolidation_type = typeof(surcol.consolidation.preconsolidation)
        col_preconsolidation_type = typeof(col.consolidation.preconsolidation)

        @test supertype(typeof(surcol.groundwater)) == supertype(typeof(col.groundwater))
        @test surcol_preconsolidation_type == col_preconsolidation_type
        @test assert_eltype(surcol.consolidation, col.consolidation)
        @test assert_eltype(surcol.oxidation, col.oxidation)
        @test assert_eltype(surcol.shrinkage, col.shrinkage)
    end

    @testset "prepare_surcharge_column_nullox_nullshr" begin
        col = AtlansFixtures.soil_column_hg_abc_nullox_nullshr()
        f = get_surcharge_forcing()

        surcol = Atlans.prepare_surcharge_column(f, col, I)

        surcol_preconsolidation_type = typeof(surcol.consolidation.preconsolidation)
        col_preconsolidation_type = typeof(col.consolidation.preconsolidation)

        @test supertype(typeof(surcol.groundwater)) == supertype(typeof(col.groundwater))
        @test surcol_preconsolidation_type == col_preconsolidation_type
        @test assert_eltype(surcol.consolidation, col.consolidation)
        @test assert_eltype(surcol.oxidation, col.oxidation)
        @test assert_eltype(surcol.shrinkage, col.shrinkage)
    end

    @testset "prepare_surcharge_column_nullcons" begin
        col = AtlansFixtures.soil_column_hg_nullcons_cs_nullshr()
        f = get_surcharge_forcing()

        surcol = Atlans.prepare_surcharge_column(f, col, I)

        surcol_preconsolidation_type = typeof(surcol.consolidation.preconsolidation)
        col_preconsolidation_type = typeof(col.consolidation.preconsolidation)

        @test supertype(typeof(surcol.groundwater)) == supertype(typeof(col.groundwater))
        @test surcol_preconsolidation_type == col_preconsolidation_type
        @test assert_eltype(surcol.consolidation, col.consolidation)
        @test assert_eltype(surcol.oxidation, col.oxidation)
        @test assert_eltype(surcol.shrinkage, col.shrinkage)
    end
end