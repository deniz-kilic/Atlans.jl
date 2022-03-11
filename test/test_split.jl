@testset "Split column" begin
    @testset "find_split_index" begin
        Δz = fill(0.25, 4)
        z = fill(0.125, 4)
        z[2:end] .+= cumsum(Δz[2:end])

        i, lower, upper = Atlans.find_split_index(z, Δz, 0.25, 0.0)
        @test isnothing(i)
        @test isnan(lower)
        @test isnan(upper)

        i, lower, upper = Atlans.find_split_index(z, Δz, -0.1, 0.0)
        @test isnothing(i)
        @test isnan(lower)
        @test isnan(upper)

        i, lower, upper = Atlans.find_split_index(z, Δz, 1.1, 0.0)
        @test isnothing(i)
        @test isnan(lower)
        @test isnan(upper)

        i, lower, upper = Atlans.find_split_index(z, Δz, 0.125, 0.0)
        @test i == 1
        @test lower ≈ 0.125
        @test upper ≈ 0.125

        i, lower, upper = Atlans.find_split_index(z, Δz, 0.8, 0.0)
        @test i == 4
        @test lower ≈ 0.05
        @test upper ≈ 0.2

        i, lower, upper = Atlans.find_split_index(z, Δz, 0.8, 0.1)
        @test isnothing(i)
        @test isnan(lower)
        @test isnan(upper)
    end

    @testset "shouldsplit" begin
        v = [1, 2, 3, 4]
        @test !(Atlans.shouldsplit(v, 4))
        @test Atlans.shouldsplit(v, 5)
        @test_throws ErrorException("vector is not equal to newlength or newlength - 1") Atlans.shouldsplit(
            v,
            3,
        )
        @test_throws ErrorException("vector is not equal to newlength or newlength - 1") Atlans.shouldsplit(
            v,
            6,
        )
    end

    @testset "cellsplit" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        consolidation = column.consolidation

        Atlans.cellsplit!(consolidation, 1, 4, NaN, NaN)
        @test length(consolidation.cells) == 4

        Atlans.cellsplit!(consolidation, 1, 5, 0.25, 0.75)
        @test length(consolidation.cells) == 5
        @test consolidation.cells[1].Δz == 0.25
        @test consolidation.cells[2].Δz == 0.75
        @test consolidation.cells[3].Δz == 1.0

        column = AtlansFixtures.soil_column_hg_abc_cs()
        consolidation = column.consolidation
        Atlans.cellsplit!(consolidation, 3, 5, 0.25, 0.75)
        @test length(consolidation.cells) == 5
        @test consolidation.cells[1].Δz == 1.0
        @test consolidation.cells[2].Δz == 1.0
        @test consolidation.cells[3].Δz == 0.25
        @test consolidation.cells[4].Δz == 0.75
        @test consolidation.cells[5].Δz == 1.0

        column = AtlansFixtures.soil_column_hg_abc_cs()
        oxidation = column.oxidation
        Atlans.cellsplit!(oxidation, 3, 5, 0.25, 0.75)
        @test length(oxidation.cells) == 5
        @test oxidation.cells[1].Δz == 1.0
        @test oxidation.cells[2].Δz == 1.0
        @test oxidation.cells[3].Δz == 0.25
        @test oxidation.cells[4].Δz == 0.75
        @test oxidation.cells[5].Δz == 1.0
    end

    @testset "zsplit!" begin
        Δz = fill(1.0, 4)
        z = fill(0.5, 4)
        z[2:end] .+= cumsum(Δz[2:end])

        Atlans.zsplit!(Δz, 1, 5, 0.25, 0.75)
        @test all(Δz .== [0.25, 0.75, 1.0, 1.0, 1.0])
    end

    @testset "columnsplit consolidation" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        consolidation = column.consolidation

        Atlans.columnsplit!(consolidation, 1, 5, 0.25, 0.75)
        @test length(consolidation.cells) == 5
        @test length(consolidation.σ) == 5
        @test length(consolidation.σ′) == 5
        @test length(consolidation.p) == 5
        @test length(consolidation.result) == 5
    end

    @testset "columnsplit oxidation" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        oxidation = column.oxidation

        Atlans.columnsplit!(oxidation, 1, 5, 0.25, 0.75)
        @test length(oxidation.cells) == 5
        @test length(oxidation.result) == 5
    end

    @testset "split column" begin
        column = AtlansFixtures.soil_column_hg_abc_cs()
        # First test for no effect
        Atlans.split!(column, 0.0, 0.0)
        @test length(column.z) == 4
        # Current boundary
        Atlans.split!(column, 1.0, 0.0)
        @test length(column.z) == 4
        # Below column
        Atlans.split!(column, -0.5, 0.0)
        @test length(column.z) == 4
        # Above column
        Atlans.split!(column, 10.0, 0.0)
        @test length(column.z) == 4
        # Within tolerance
        Atlans.split!(column, 0.05, 0.1)
        @test length(column.z) == 4
        Atlans.split!(column, 0.95, 0.1)
        @test length(column.z) == 4

        Atlans.split!(column, 0.25, 0.01)
        @test all(column.Δz .== [0.25, 0.75, 1.0, 1.0, 1.0])
        @test all(column.z .== [0.125, 0.625, 1.5, 2.5, 3.5])

        @test length(column.consolidation.z) == 5
        @test length(column.consolidation.Δz) == 5
        @test length(column.consolidation.cells) == 5
        @test length(column.consolidation.σ) == 5
        @test length(column.consolidation.σ′) == 5
        @test length(column.consolidation.p) == 5
        @test length(column.consolidation.result) == 5
        @test length(column.oxidation.z) == 5
        @test length(column.oxidation.Δz) == 5
        @test length(column.oxidation.cells) == 5
        @test length(column.oxidation.result) == 5

    end
end
