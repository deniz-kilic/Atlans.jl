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

    cc = Atlans.ConsolidationColumn(cells, z, Δz, σ, σ′, p, preconsolidation, result)
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

    m_clay = 0.8
    m_organic = 0.1

    cells = [Atlans.SimpleShrinkage(i, n, m_clay, m_organic) for (i, n) in zip(Δz, n_vals)]
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
        shrinkage,
    )

end


workdir = raw"n:\Projects\11205500\11205981\B. Measurements and calculations\WP3\Fase 3 - prognose case\kem-atlantis\data"
path_nc = joinpath(workdir, "subsoil-model-fase2.nc")
path_csv = joinpath(workdir, "parameters.csv")
subsoil = Atlans.prepare_subsoil_data(path_nc, path_csv)

ncells = 10
thickness = 0.5
zbase = -5.0

soilcolumn = create_soilcolumn(ncells, thickness, zbase) # soilcolumn with all Atlans attributes

timestepper = Atlans.ExponentialTimeStepper(1.0, 2)
timesteps = Atlans.create_timesteps(timestepper, 3650.0)

## From here the model is 'running'
Atlans.apply_preconsolidation!(soilcolumn)
Atlans.prepare_forcingperiod!(soilcolumn, 0.01, 0.0, -1.0)
Atlans.set_phreatic_difference!(soilcolumn, -1.0)
s, c, o, shr = Atlans.advance_forcingperiod!(soilcolumn, timesteps)

# phreatic = Atlans.phreatic_level(soilcolumn.groundwater)
# for t in timesteps
#     Atlans.shrink!(soilcolumn.shrinkage, phreatic, t)
# end
# println(soilcolumn.shrinkage.result)

# @show sum(0.5 * 20 - sum(soilcolumn.Δz))
