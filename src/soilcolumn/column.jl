"""
x, y, z are all midpoints.
"""
struct SoilColumn
    x::Float
    y::Float
    z::Vector{Float}
    Δz::Vector{Float}
    consolidation::ConsolidationColumn
    oxidation::OxidationColumn
    groundwater::GroundwaterColumn
end

function deep_subsidence(
    column::Union{SoilColumn,ConsolidationColumn,OxidationColumn},
    subsidence::Float,
)
    column.z .-= subsidence
end

function advance_timestep!(c::SoilColumn, Δt::Float)
    total_oxidation = 0.0
    total_consolidation = 0.0

    flow!(c.groundwater, Δt)
    # Exchange pore pressure from groundwater column to consolidation column
    @. c.consolidation.u = c.groundwater.p * γ_water
    consolidate!(c.consolidation, Δt)
    oxidate!(c.oxidation, Δt)

    # Now collect oxidation and consolidation and update Δz and z
    zbot = c.z[1] - 0.5 * c.Δz[1]
    for index in eachindex(c.Δz)
        ccel = c.consolidation[index]
        ocel = c.oxidation[index]
        # Δz should not become negative
        ΔΔz = min(ccel.consolidation + ocel.oxidation, c.Δz[index])
        total_consolidation += ccel.oxidation
        total_oxidation += ocel.oxidation
        # Compute new Δz
        Δz = c.Δz[index] - ΔΔz
        c.Δz[index] = Δz
        # Compute new midpoint
        c.z[index] = zbot + 0.5 * Δz
        zbot += Δz
    end
    return total_consolidation, total_oxidation
end

function advance_stressperiod!(c::SoilColumn, timesteps::Vector{Float})
    for Δt in timesteps
        advance_timestep!(c, Δt)
    end
end

function consolidate!(cc::ConsolidationColumn, Δt::Float)
    # Compute the new effective stress.
    total_stress!(cc, phreatic_level)
    # Pore pressure has been computed by the groundwater column.
    effective_stress!(cc)
    # Using the change in effective stress, compute consolidation.
    for (index, cell) in enumerate(cc.cells)
        σ′ = cc.σ′[index]
        cc.cells[index] = consolidate(cell, σ′, Δt)
    end
end

function oxidate!(oc::OxidationColumn, Δt::Float)
    for (index, cell) in enumerate(oc.cells)
        oc.oxidation[index] = oxidate(cell, Δt)
    end
end


# 2022
function set_surcharge!(c::SoilColumn, material) end
function remove_surcharge!(c::SoilColumn, depth::Float) end
