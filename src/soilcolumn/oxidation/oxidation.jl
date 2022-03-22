struct OxidationColumn{O}
    cells::Vector{O}
    z::Vector{Float}
    Δz::Vector{Float}
    result::Vector{Float}
    max_oxidation_depth::Float
end

oxidation_depth(column::OxidationColumn{NullOxidation}, _, _, _, _) = 0.0
oxidate!(column::OxidationColumn{NullOxidation}, _, _) = nothing
synchronize!(column::OxidationColumn{NullOxidation}, _) = nothing

function oxidation_depth(
    column::OxidationColumn{O},
    surface_level,
    phreatic_level,
    deep_subsidence,
    phreatic_change,
) where {O<:OxidationProcess}
    new_surface = surface_level - deep_subsidence
    new_phreatic = phreatic_level + phreatic_change
    depth = max(0.0, new_surface - new_phreatic)
    depth = min(column.max_oxidation_depth, depth)
    return surface_level - depth
end

function oxidate!(
    column::OxidationColumn{O},
    phreatic_level::Float,
    Δt::Float,
) where {O<:OxidationProcess}
    oxidation_z = max(phreatic_level, surface_level(column) - column.max_oxidation_depth)
    column.result .= 0.0
    for index in reverse(1:length(column.cells))
        if column.z[index] > oxidation_z
            cell = column.cells[index]
            newcell = oxidate(cell, Δt)
            column.cells[index] = newcell
            column.result[index] = newcell.oxidation
        else
            break
        end
    end
    return
end

function synchronize!(column::OxidationColumn{O}, Δz) where {O<:OxidationProcess}
    for i = 1:length(column.cells)
        @set column.cells[i].Δz = Δz[i]
    end
    return
end
