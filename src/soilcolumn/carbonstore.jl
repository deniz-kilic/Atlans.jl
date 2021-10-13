using SpecialFunctions: erf

struct CarbonStore <: OxidationProcess
    Δz::Float64  # cell thickness
    f_organic::Float64  # mass fraction organic material
    m_organic::Float64  # mass of organic material
    m_mineral::Float64  # mass of mineral material
    m_minimum_organic::Float64  # no breakdown beyond this level
    α::Float64  # oxidation rate
    ρb::Float64  # bulk density
end

function mass_organic(f_organic, ρb, Δz)
    return f_organic * ρb * Δz
end

function mass_mineral(f_organic, ρb, Δz)
    return (1.0 - f_organic) * ρb * Δz
end

function mass_organic_minimal(m_mineral, f_minimum_organic)
    return (m_mineral * f_minimum_organic / (1.0 - f_minimum_organic))
end

function ρ_bulk(m_organic, m_mineral, Δz)
    return (m_organic + m_mineral) / Δz
end

function fraction_organic(m_organic, m_mineral)
    return m_organic / (m_organic + m_mineral)
end

"""
Empirical equation to compute specific volume of organic material.

As the organic matter in a soil breaks down, density increases. The "airest"
parts are the first to go.
"""
function volume_organic(f_organic, ρb)
    if f_organic == 0
        return 0
    else
        return 0.5 * (f_organic / ρb) * (1.0 + erf((f_organic - 0.2) / 0.1))
    end
end

function oxidate(cs::CarbonStore, Δt::Float64)::Tuple{Float64,CarbonStore}
    if cs.α == 0
        return 0.0, cs
    end
    Δm = min(cs.m_organic - cs.m_minimum_organic, cs.α * cs.Δz * Δt)
    Δm = max(0, Δm)
    m_organic = cs.m_organic - Δm
    oxidation = volume_organic(cs.f_organic, cs.ρb) * Δm
    f_organic = fraction_organic(m_organic, cs.m_mineral)
    Δz = cs.Δz - oxidation
    return oxidation,
    CarbonStore(
        Δz,  # new
        f_organic,  # new
        cs.m_mineral,
        m_organic, # new
        cs.m_minimum_organic,
        cs.α,
        ρ_bulk(m_organic, cs.m_mineral, Δz),  # new
    )
end

