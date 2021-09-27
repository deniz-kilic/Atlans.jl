
@testset "DrainingAbcIsotache" begin
    cell = DrainingAbcIsotache(
        1.0,
        0.0,
        0.0,
        18000.0,
        13000.0,
        c_d,
        c_v,
        U,
        a,
        b,
        c,
        τ,
    )
    
    @testset "Initialization" begin
        @test typeof(cell) == DrainingAbcIsotache
    end

    @testset "τ_intrinsic" begin
        loadstep = 1.0
        actual = τ_intermediate(cell, loadstep)
        expected = 100.0
        @test actual == expected
        @test actual ≈ expected
    end


end