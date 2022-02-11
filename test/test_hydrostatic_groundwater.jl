@testset "Hydrostatic Groundwater" begin
    z = collect(0.5:1.0:5.0)
    phreatic = Atlans.Phreatic(4.0)
    dry = fill(false, 5)
    p = fill(NaN, 5)
    
    hg1 = Atlans.HydrostaticGroundwater(
        z, phreatic, dry, p
    )

    @testset "Constructor" begin
        @test typeof(hg1) == Atlans.HydrostaticGroundwater
        
        hg2 = Atlans.HydrostaticGroundwater(z, 4.0)
        @test typeof(hg2) == Atlans.HydrostaticGroundwater
        
        @test isequal(hg1.z, hg2.z)
        @test isequal(hg1.dry, hg2.dry)
        @test isequal(hg1.p, hg2.p)
        @test isequal(hg1.phreatic.ϕ, hg2.phreatic.ϕ)
    end
    
    @testset "set_phreatic_difference" begin
        Atlans.set_phreatic_difference!(hg1, -0.5) 
        @test hg1.phreatic.ϕ ≈ 3.5
        Atlans.set_phreatic_difference!(hg1, 0.5) 
        @test hg1.phreatic.ϕ ≈ 4.0
    end

    @testset "set_phreatic" begin
        Atlans.set_phreatic!(hg1, 2.0) 
        @test hg1.phreatic.ϕ ≈ 2.0
        Atlans.set_phreatic!(hg1, 4.0) 
        @test hg1.phreatic.ϕ ≈ 4.0
    end

    @testset "flow" begin
        Atlans.flow!(hg1, 0.0)
    end
    
    @testset "phreatic_level" begin
        @test Atlans.phreatic_level(hg1) == hg1.phreatic.ϕ
    end
    
    @testset "pore_pressure" begin
        Atlans.pore_pressure!(hg1)
        @test all(hg1.dry == [false, false, false, false, true])
        @test all(hg1.p ≈ [3.5, 2.5, 1.5, 0.5, 0.0])
    end
end