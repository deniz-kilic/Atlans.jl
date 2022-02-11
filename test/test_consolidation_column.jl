using Setfield

@testset "ConsolidationColumn" begin

    Δz = 1.0
    t = 0.0
    σ′ = 10000.0
    γ_wet = 15000.0
    γ_dry = 10000.0
    c_d = 2.0
    c_v = 0.006912
    U = 0.0
    a = 0.01737
    b = 0.1303
    c = 0.008686
    τ = 1.0
    consolidation = 0.0
    
    cell = Atlans.DrainingAbcIsotache(
        Δz,
        Δz,
        t,
        σ′,
        γ_wet,
        γ_dry,
        c_d,
        c_v,
        U,
        a,
        b,
        c,
        τ,
        consolidation,
    )

    cells = fill(cell, 4)
    z = collect(0.5:1.0:4.0)
    Δz = fill(1.0, 4)
    σ = fill(NaN, 4)
    σ′ = fill(NaN, 4)
    p = fill(NaN, 4)
    result = fill(NaN, 4)
    preconsolidation = Atlans.OverConsolidationRatio(fill(2.15, 4))
    column = Atlans.ConsolidationColumn(cells, z, Δz, σ, σ′, p, preconsolidation, result)

    @testset "ConsolidationColumn constructor" begin
        @test typeof(column) == Atlans.ConsolidationColumn{Atlans.DrainingAbcIsotache, Atlans.OverConsolidationRatio}
        
        column2 = Atlans.ConsolidationColumn(cells, z, Δz, preconsolidation)
        @test typeof(column2) == Atlans.ConsolidationColumn{Atlans.DrainingAbcIsotache, Atlans.OverConsolidationRatio}
    end

    @testset "compress_γ_wet" begin
        cell = @set cell.consolidation = 0.1
        actual = Atlans.compress_γ_wet(cell)
        expected = 15576.666666666666
        @test actual ≈ expected
    end

    @testset "compress_γ_dry" begin
        cell = @set cell.consolidation = 0.01
        actual = Atlans.compress_γ_dry(cell)
        expected = 10101.0101010101
        @test actual ≈ expected
    end

    @testset "degree of consolidation" begin
        t = 0.01
        actual = Atlans.U(cell, t)
        expected = 5.7916576227323684e-15
        @test actual ≈ expected
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

        phreatic_level = 1.5
        Atlans.total_stress!(column, phreatic_level)
        expected = [42500.0, 25000.0, 15000.0, 5000.0]
        @test all(column.σ .≈ expected)

        phreatic_level = 2.0
        Atlans.total_stress!(column, phreatic_level)
        expected = [45000.0, 30000.0, 15000.0, 5000.0]
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
        
        @test all([cell.σ′ for cell in column.cells]  .== column.σ′)
    end
        
    @testset "consolidate" begin
        Atlans.apply_preconsolidation!(column)
        # Now lower watertable by 0.2
        column.p .= [2.3, 1.3, 0.3, 0.0] .* Atlans.γ_water
        phreatic_level = 2.8
        Δt = 700.0
        Atlans.consolidate!(column, phreatic_level, Δt)
        @show column.cells[1]
        @show column.cells[2]
        @show column.cells[3]
        @show column.cells[4]
    end
    
    @testset "U -> 1.0" begin
        
    end

end