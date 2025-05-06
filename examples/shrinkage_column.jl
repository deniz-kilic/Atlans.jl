using Atlans


function carbon_store()
    f_organic = 0.2
    ρb = 1000.0
    Δz = 0.5
    f_minimum_organic = 0.05
    α = 1.0e-3
    return Atlans.CarbonStore(Δz, f_organic, f_minimum_organic, ρb, α)
end


function consolidation_column(z, Δz)
    cells = fill(Atlans.NullConsolidation(), length(z))
    σ = fill(NaN, length(z))
    σ′ = fill(NaN, length(z))
    p = fill(NaN, length(z))
    preconsolidation = Atlans.OverConsolidationRatio(fill(2.15, length(z)))
    result = fill(0.0, length(z))

    return Atlans.ConsolidationColumn(cells, z, Δz, σ, σ′, p, preconsolidation, result)
end


function oxidation_column(z, Δz)
    cells = fill(Atlans.NullOxidation(), length(z))
    result = fill(0.0, length(z))

    return Atlans.OxidationColumn(cells, z, Δz, result, NaN, NaN)
end


function groundwater_column(z)
    phreatic = Atlans.Phreatic(-1.5)
    dry = fill(false, length(z))
    p = fill(NaN, length(z))

    return Atlans.HydrostaticGroundwater(z, phreatic, dry, p)
end


function shrinkage_column(z, Δz)
    τ_years = 60.0
    n_vals = [0.7, 0.7, 0.7, 0.7, 0.7, 1.1, 1.2, 1.6, 1.7, 1.7]

    m_clay = 0.8
    m_organic = 0.1

    cells = [Atlans.SimpleShrinkage(i, n, m_clay, m_organic) for (i, n) in zip(Δz, n_vals)]
    result = fill(NaN, length(z))
    Hv0 = 0.3

    return Atlans.ShrinkageColumn(cells, z, Δz, result, Hv0)
end


function create_soilcolumn(ncells, thickness, zbase)
    x = 0.0
    y = 0.0

    Δz = fill(thickness, ncells)
    z = (zbase .+ cumsum(Δz)) .- (thickness .* Δz)

    consolidation = consolidation_column(z, Δz) # NullConsolidation
    oxidation = oxidation_column(z, Δz) # NullOxidation
    groundwater = groundwater_column(z)
    shrinkage = shrinkage_column(z, Δz)

    # oxidation = Atlans.OxidationColumn(
    #     fill(carbon_store(), ncells), z, Δz, fill(0.0, ncells), 1.2
    # )

    return Atlans.SoilColumn(
        zbase,
        x,
        y,
        z,
        Δz,
        groundwater,
        consolidation,
        oxidation,
        shrinkage,
    )

end


ncells = 10
thickness = 0.5
zbase = -5.0

#%%
ad = Atlans.AdaptiveCellsize(0.25, 0.01)
timestepper = Atlans.ExponentialTimeStepper(1.0, 2)
timesteps = Atlans.create_timesteps(timestepper, 3650.0)

soilcolumn = create_soilcolumn(ncells, thickness, zbase) # soilcolumn with all Atlans attributes

Atlans.apply_preconsolidation!(soilcolumn)
Atlans.prepare_forcingperiod!(soilcolumn, 0.01, 0.0, 0.0)
Atlans.set_phreatic_difference!(soilcolumn, -1.0)

#%%
println(soilcolumn.z)
s, c, o, shr = Atlans.advance_forcingperiod!(soilcolumn, timesteps)
println(soilcolumn.z)
