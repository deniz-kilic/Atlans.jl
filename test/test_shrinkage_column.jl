@testset "ShrinkageColumn" begin

    column = AtlansFixtures.shrinkage_column()

    @testset "constructor" begin
        @test typeof(column) == Atlans.ShrinkageColumn{
            Atlans.SimpleShrinkage
        }
    end

    @testset "surface_level" begin
        @test Atlans.surface_level(column) ≈ 4.0
    end

    @testset "shrinkage_depth" begin
        phreatic = 2.5
        sz = max(
            phreatic,
            Atlans.surface_level(column) - column.max_shrinkage_depth
        )
        @test sz ≈ 3.0
    end

    @testset "shrink" begin
        phreatic = 2.5
        Δt = 3650.0

        Atlans.shrink!(column, phreatic, Δt)

        expected_shr = [0.0, 0.0, 0.0, 0.12372923070259101]
        expected_Δz = [1.0, 1.0, 1.0, 0.876270769297409]

        @test all(column.result .≈ expected_shr)
        @test all(c.Δz ≈ Δz for (c, Δz) in zip(column.cells, expected_Δz))
    end
end