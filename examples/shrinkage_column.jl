using Atlans


function draining_abc_isotache(Δz)
    t = 0.0
    σ′ = NaN
    γ_wet = 12_500.0
    γ_dry = 10_500.0
    c_d = 1.0
    c_v = 0.006912
    U = 0.0
    a = 0.01737
    b = 0.1303
    c = 0.008686
    τ = NaN
    consolidation = NaN

    return Atlans.DrainingAbcIsotache(
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
end


function consolidation_column(z, Δz)
    cells = [draining_abc_isotache(i) for i in Δz]
    σ = fill(NaN, length(z))
    σ′ = fill(NaN, length(z))
    p = fill(NaN, length(z))
    preconsolidation = Atlans.OverConsolidationRatio(fill(2.15, length(z)))
    result = fill(NaN, length(z))

    cc = Atlans.ConsolidationColumn(
        cells, z, Δz, σ, σ′, p, preconsolidation, result
    )
end


function oxidation_column(z, Δz)
    cells = fill(Atlans.NullOxidation(), length(z))
    result = fill(0.0, length(z))
    max_oxidation_depth = NaN

    oc = Atlans.OxidationColumn(cells, z, Δz, result, max_oxidation_depth)
end


function groundwater_column(z)
    phreatic = Atlans.Phreatic(-1.5)
    dry = fill(false, length(z))
    p = fill(NaN, length(z))

    gw = Atlans.HydrostaticGroundwater(z, phreatic, dry, p)
end


function shrinkage_column(z, Δz)
    τ_years = 60.0
    n_vals = [0.7, 0.7, 0.7, 0.7, 0.7, 1.1, 1.2, 1.6, 1.7, 1.7]
    cells = [Atlans.SimpleShrinkage(i, n, τ_years, 3.0, NaN) for (i, n) in zip(Δz, n_vals)]
    result = fill(NaN, length(z))
    max_shrinkage_depth = 1.3

    sc = Atlans.ShrinkageColumn(cells, z, Δz, result, max_shrinkage_depth)
end


function create_soilcolumn(ncells, thickness, zbase)
    x = 0.0
    y = 0.0

    Δz = fill(thickness, ncells)
    z = (zbase .+ cumsum(Δz)) .- (thickness .* Δz)

    consolidation = consolidation_column(z, Δz)
    oxidation = oxidation_column(z, Δz)
    groundwater = groundwater_column(z)
    shrinkage = shrinkage_column(z, Δz)

    soilcolumn = Atlans.SoilColumn(
        zbase,
        x,
        y,
        z,
        Δz,
        groundwater,
        consolidation,
        oxidation,
        shrinkage
    )

end


ncells = 10
thickness = 0.5
zbase = -5.0
τ_years = 30.0

cell = Atlans.SimpleShrinkage(0.5, 1.8, τ_years, 3.0, NaN) # single cell to shrink
soilcolumn = create_soilcolumn(ncells, thickness, zbase) # soilcolumn with all Atlans attributes

phreatic = Atlans.phreatic_level(soilcolumn.groundwater)
time = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 3650]

println(soilcolumn.shrinkage.result)
for t in diff(Float64.(time))
    Atlans.shrink!(soilcolumn.shrinkage, phreatic, t)
end
println(soilcolumn.shrinkage.result)
println(sum(soilcolumn.shrinkage.result))