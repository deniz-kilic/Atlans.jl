struct OxidationColumn{O}
    cells::Vector{O}
    z::Vector{Float}
    Δz::Vector{Float}
    result::Vector{Float}
end

function oxidate!(oc::OxidationColumn, Δt::Float)
    for (index, cell) in enumerate(oc.cells)
        oc.oxidation[index] = oxidate(cell, Δt)
    end
end

function result!(column::OxidationColumn)
    for (i, cell) in enumerate(column.cells)
        result[i] = cell.oxidation
    end
    return result
end

function synchronize!(column::OxidationColumn, Δz)
    for i=1:length(column.cells)
        @set column.cells[i].Δz = Δz[i]
    end
    return
end
