struct OxidationColumn{O}
    cells::Vector{O}
    z::Vector{Float}
    Δz::Vector{Float}
    result::Vector{Float}
    max_oxidation_depth::Float
end

function oxidate!(column::OxidationColumn, phreatic_level::Float, Δt::Float)
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

function synchronize!(column::OxidationColumn, Δz)
    for i = 1:length(column.cells)
        @set column.cells[i].Δz = Δz[i]
    end
    return
end
