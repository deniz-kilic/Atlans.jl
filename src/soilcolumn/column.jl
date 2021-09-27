struct SoilColumn
    x::Float
    y::Float
    z::Vector{Float}
    Δz::Vector{Float}
    consolidation::C where {C<:ConsolidationProcess}
    oxidation::O where {O<:OxidationProcess}
    groundwater::G where {G<:GroundwaterProcess}
end

function advance_timestep!(c::SoilColumn, Δt::Float) end

function advance_stressperiod!(c::SoilColumn, timesteps::Vector{Float}) end

function flow!(c::SoilColumn, Δt::Float) end

function consolidate!(c::SoilColumn, Δt::Float) end

function oxidate!(c::SoilColumn, Δt::Float) end

function update!(c::SoilColumn, Δt::Float) end

function set_surcharge!(c::SoilColumn, material) end

function remove_surcharge!(c::SoilColumn, depth::Float) end

function set_groundwatertable!(c::SoilColumn, value::Float) end

function split!(c::SoilColumn, level::Float)
    # split at depth
    # insert new split cell
    # push existing cells "upward"
    # cell = pop!(cells, index)
    #insert!(cells, index, newcell1)
    #insert!(cells, index, newcell2)
end

function oxidation_depth(c::SoilColumn) end
