struct OxidationColumn{O}
    cells::Vector{O}
    z::Vector{Float}
    Δz::Vector{Float}
    result::Vector{Float}
    max_oxidation_depth::Float
end

oxidation_depth(column::OxidationColumn{NullOxidation}, _) = 0.0
oxidate!(column::OxidationColumn{NullOxidation}, _, _) = nothing
synchronize!(column::OxidationColumn{NullOxidation}, _) = nothing

function oxidation_depth(
    column::OxidationColumn{O},
    phreatic_level,
) where {O<:OxidationProcess}
    return min(column.max_oxidation_depth, phreatic_level)
end


function oxidate!(
    column::OxidationColumn{O},
    phreatic_level::Float,
    Δt::Float,
) where {O<:OxidationProcess}
    oxidation_z = max(phreatic_level, surface_level(column) - column.max_oxidation_depth)
    column.result .= 0.0
    for (index, cell) in enumerate(column.cells)
        if column.z[index] > oxidation_z
            newcell = oxidate(cell, Δt)
            column.cells[index] = newcell
            column.result[index] = newcell.oxidation
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
