mutable struct Base
    z::Float
end

mutable struct Surface
    z::Float
end

"""
x, y, z are all midpoints.
"""
struct SoilColumn{G,C,P,O}
    base::Base
    surface::Surface
    x::Float
    y::Float
    z::Vector{Float}
    Δz::Vector{Float}
    groundwater::G
    consolidation::ConsolidationColumn{C,P}
    oxidation::OxidationColumn{O}
    subsidence::Vector{Float}
end

function SoilColumn(
    base::Float,
    surface::Float,
    x,
    y,
    z,
    Δz,
    groundwater,
    consolidation,
    oxidation,
)
    return SoilColumn(
        Base(base),
        Surface(surface),
        x,
        y,
        z,
        Δz,
        groundwater,
        consolidation,
        oxidation,
        fill(NaN, length(z)),
    )
end

function initial_stress!(column)
    flow!(column.groundwater, 0.0)
    pore_pressure!(column.groundwater)
    exchange_pore_pressure!(column)
    total_stress!(column.consolidation, phreatic_level(column.groundwater))
    effective_stress!(column.consolidation)
    transfer_stress!(column.consolidation)
    return
end

function apply_preconsolidation!(column::SoilColumn)
    initial_stress!(column)
    apply_preconsolidation!(column.consolidation)
end

# Forcing

#struct Surcharge{G,C,O}
#    Δz::Vector{Float}
#    groundwater::G
#    consolidation::ConsolidationColumn{C}
#    oxidation::GroundwaterColumn{O}
#end
#
#function set_surcharge!(
#    column::SoilColumn,
#    surcharge,
#)
#    set_surcharge!(column.groundwater, surcharge.groundwater)
#    set_surcharge!(column.consolidation, surcharge.consolidation)
#    set_surcharge!(column.oxidation, surcharge.oxidation)
#    return
#end

function set_deep_subsidence!(
    column::SoilColumn,
    subsidence::Float,
)
    column.base.z -= subsidence
    update_z!(column)
end

function set_aquifer!(
    column::SoilColumn,
    ϕ,
)
    set_aquifer!(column.groundwater, ϕ)
    return
end
 
function set_aquifer_difference!(
    column::SoilColumn,
    Δϕ,
)
    set_aquifer_difference!(column.groundwater, Δϕ)
    return
end

function set_phreatic!(
    column::SoilColumn,
    ϕ,
)
    set_phreatic!(column.groundwater, ϕ)
    return
end

function set_phreatic_difference!(
    column::SoilColumn,
    Δϕ,
)
    set_phreatic_difference!(column.groundwater, Δϕ)
end

# Computation

function exchange_pore_pressure!(column::SoilColumn)
    column.consolidation.p .= column.groundwater.p * γ_water
end


"""
    prepare_timestep!(column)
    
Prepare a single timestep. This updates stresses in the column and accounts for
e.g. drowning of the column.

This computes:

    * Pore pressure
    * Total stress
    * Effective stress
"""
function prepare_timestep!(column::SoilColumn)
    flow!(column.groundwater, 0.0)
    exchange_pore_pressure!(column)
    total_stress!(column.consolidation, phreatic_level(column.groundwater))
    effective_stress!(column.consolidation)
end

"""
    prepare_forcingperiod!(column)
    
Prepare a single forcing period.

This:

    * splits the soil column at maximum oxidation depth (if necessary)
    * computes pore pressure, total stress, effective stress prior to this stress periods loading
    * sets the effective stress in every consolidation cell
    * Reset U and t for DrainingConsolidation processes.
"""
function prepare_forcingperiod!(column::SoilColumn)
    #level = oxidation_depth(column.oxidation)
    #split!(column, level)
    
    initial_stress!(column)

    prepare_forcingperiod!(column.consolidation)
    return
end


"""
    function update_z!(column)
    
Compute new midpoints and surface level.
"""
function update_z!(column::SoilColumn)
    column.z .= column.base.z .+ cumsum(column.Δz) .- 0.5 .* column.Δz
    column.surface.z = column.base.z + sum(column.Δz)
end

"""
    subside!(column)
    
Apply consolidation and oxidation to thickness
"""
function subside!(column::SoilColumn)
    # Δz should not become negative
    column.subsidence .= min.((column.consolidation.result .+ column.oxidation.result), column.Δz)
    column.Δz .-= column.subsidence
    synchronize!(column.groundwater, column.Δz)
    synchronize!(column.consolidation, column.Δz)
    synchronize!(column.oxidation, column.Δz)
    update_z!(column)
end

"""
    advance_timestep!(column, Δt)

Advance a single timestep.

During a timestep the following states are computed:

    * head
    * pore pressure
    * total stress
    * effective stress

Finally, thickness and elevation are updated.
"""
function advance_timestep!(c::SoilColumn, Δt::Float)
    flow!(c.groundwater, Δt)
    exchange_pore_pressure!(c)
    consolidate!(c.consolidation, phreatic_level(c.groundwater), Δt)
    oxidate!(c.oxidation, Δt)
    subside!(c)
    return sum(c.subsidence), sum(c.consolidation.result), sum(c.oxidation.result)
end

"""
    advance_forcingperiod!(column, timesteps)
    
Advances a prepared column by a number of timesteps.
"""
function advance_forcingperiod!(c::SoilColumn, timesteps::Vector{Float})
    subsidence = 0.0
    consolidation = 0.0
    oxidation = 0.0
    for Δt in timesteps
        Δs, Δc, Δo = advance_timestep!(c, Δt)
        subsidence += Δs
        consolidation += Δo
        oxidation += Δc
    end
    return subsidence, consolidation, oxidation
end

# Output

output(oc::OxidationColumn) = sum(cell.oxidation for cell in oc.cells)
output(cc::ConsolidationColumn) = sum(cell.oxidation for cell in cc.cells)
output(gw::GW where {GW<:GroundwaterColumn}) = phreatic_level(gw)
