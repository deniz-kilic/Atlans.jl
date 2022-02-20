struct OxidationColumn{O}
    cells::Vector{O}
    z::Vector{Float}
    Δz::Vector{Float}
    result::Vector{Float}
end

function oxidate!(column::OxidationColumn, Δt::Float)
    for (index, cell) in enumerate(column.cells)
        newcell = oxidate(cell, Δt)
        column.cells[index] = newcell
        column.result[index] = newcell.oxidation
    end
    return
end

function synchronize!(column::OxidationColumn, Δz)
    for i=1:length(column.cells)
        @set column.cells[i].Δz = Δz[i]
    end
    return
end
