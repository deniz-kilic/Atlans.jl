@testset "SimpleShrinkage" begin
    Δz = 1.0
    n = 1.7
    m_clay = 0.8
    m_organic = 0.1

    cell = Atlans.SimpleShrinkage(Δz, n, m_clay, m_organic)

    @testset "Initialization" begin
        @test typeof(cell) == Atlans.SimpleShrinkage
        @test cell.sf ≈ 2.45736988
    end

    @testset "Shrink" begin
        Δt = 1.0
        new_cell = Atlans.shrink(cell, Δt)
        @test typeof(new_cell) == Atlans.SimpleShrinkage
        @test new_cell.shrinkage ≈ 3.73750395e-5
    end
end