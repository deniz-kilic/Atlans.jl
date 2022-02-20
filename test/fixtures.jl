module AtlansFixtures 

using Atlans

function draining_abc_isotache()
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


function carbon_store()
    return Atlans.CarbonStore(
        1.0,
        0.2,
        Atlans.mass_organic(0.2, 1000, 1),
        Atlans.mass_mineral(0.2, 1000, 1),
        Atlans.mass_organic_minimal(Atlans.mass_mineral(0.2, 1000, 1), 0.05),
        1.0e-3,
        1000,
        0.0,
    )
end

function draining_abc_isotache_column()
    cell = draining_abc_isotache()

    cells = fill(cell, 4)
    z = collect(0.5:1.0:4.0)
    Δz = fill(1.0, 4)
    σ = fill(NaN, 4)
    σ′ = fill(NaN, 4)
    p = fill(NaN, 4)
    result = fill(NaN, 4)
    preconsolidation = Atlans.OverConsolidationRatio(fill(2.15, 4))

    return Atlans.ConsolidationColumn(cells, z, Δz, σ, σ′, p, preconsolidation, result)
end

function carbon_store_column()
    cell = carbon_store()

    cells = fill(cell, 4)
    z = collect(0.5:1.0:4.0)
    Δz = fill(1.0, 4)
    result = fill(NaN, 4)
    
    return Atlans.OxidationColumn(cells, z, Δz, result)
end

function hydrostatic_groundwater()
    z = collect(0.5:1.0:4.0)
    phreatic = Atlans.Phreatic(3.0)
    dry = fill(false, 4)
    p = fill(NaN, 4)
    
    return Atlans.HydrostaticGroundwater(z, phreatic, dry, p)
end

function soil_column_hg_abc_cs()
    z = collect(0.5:1.0:4.0)
    Δz = fill(1.0, 4)

    cc = Atlans.ConsolidationColumn(
        fill(draining_abc_isotache(), 4),
        z,
        Δz,
        fill(NaN, 4), # σ
        fill(NaN, 4), # σ′
        fill(NaN, 4), # p
        Atlans.OverConsolidationRatio(fill(2.15, 4)),
        fill(NaN, 4), # result
    )

    oc = Atlans.OxidationColumn(
        fill(carbon_store(), 4),
        z, 
        Δz,
        fill(NaN, 4),
    )
    
    gw = Atlans.HydrostaticGroundwater(
        z,
        Atlans.Phreatic(3.0),
        fill(false, 4),
        fill(NaN, 4),
    )
    
    return Atlans.SoilColumn(
        0.0,
        4.0,
        0.0,
        0.0,
        z,
        Δz,
        gw,
        cc,
        oc,
    )
end

end # module
