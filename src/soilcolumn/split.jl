"""
Logic for splitting a column to accomodate for a moving phreatic level in
combination with organic matter stores.
"""

function find_split_index(z, Δz, level)
    for index in eachindex(z)
        zbot = z[index] - 0.5 * Δz[index]
        ztop = z[index] + 0.5 * Δz[index]
        if level < ztop && level > zbot
            return index, level - zbot, ztop - level
        end
    end
    return 0, NaN, NaN
end

function shouldsplit(vector, newlength)
    if length(vector) == newlength
        return false
    elseif length(vector) != (newlength - 1)
        return true
    else
        error("vector is equal to newlength or newlength - 1")
    end
end

function cellsplit!(
    column::Union{ConsolidationColumn,OxidationColumn},
    index,
    newlength,
    lowerΔz,
    upperΔz,
)
    if shouldsplit(column.cells, newlength)
        # @set also creates a copy
        cell = column.cells[index]
        @set cell.Δz = lowerΔz
        lower = copy(cell)
        @set lower.Δz = upperΔz
        insert!(column.cells, index, lower)
    end
    return
end

function insert_duplicate!(vector::Vector, index, newlength)
    if shouldsplit(vector, newlength)
        insert!(item, index, copy(item[index]))
    end
    return
end

function zsplit!(z, Δz, index, newlength, lowerΔz, upperΔz)
    if shouldsplit(vector, newlength)
        pop!(Δz, index)
        insert!(Δz, index, upperΔz)
        insert!(Δz, index, lowerΔz)

        push!(z, NaN)
        zbot = z[1] - 0.5 * Δz[1]
        for index in eachindex(c.Δz)
            z[index] = zbot + 0.5 * Δz[index]
            zbot += Δz[index]
        end
    end
    return
end

function columnsplit!(ig::InterpolatedGroundwater, index, newlength, lowerΔz, upperΔz)
    zsplit!(ig.z, ig.Δz, index, newlength, lowerΔz, upperΔz)
    # If split, another cell is introduced at the the index.
    # This means that any index equal or higher should be incremented by one.
    for (i, boundary_index) in enumerate(ig.boundary)
        if boundary_index >= index
            ig.boundary[i] = boundary_index + 1
        end
    end
    push!(ig.dry, true)
    push!(ig.p, NaN)
    return
end

function columnsplit!(cc::ConsolidationColumn, index, newlength, lowerΔz, upperΔz)
    zsplit!(cc.z, cc.Δz, index, newlength, lowerΔz, upperΔz)
    cellsplit!(cc, index, newlength, lowerΔz, upperΔz)
    push!(cc.σ, NaN)
    push!(cc.σ′, NaN)
    push!(cc.u, NaN)
end

function columnsplit!(oc::OxidationColumn, index, newlength, lowerΔz, upperΔz)
    zsplit!(oc.z, oc.Δz, index, newlength, lowerΔz, upperΔz)
    cellsplit!(cc, index, newlength, lowerΔz, upperΔz)
end

function split!(c::SoilColumn, level::Float)
    index, lowerΔz, upperΔz = find_split_index(c.z, c.Δz, level)
    newlength = length(c.z) + 1
    if index != 0
        # insert new split cell
        # push existing cells "upward"
        columnsplit!(c.groundwater, index, newlength, lowerΔz, upperΔz)
        columnsplit!(c.consolidation, index, newlength, lowerΔz, upperΔz)
        columnsplit!(c.oxidation, index, newlength, lowerΔz, upperΔz)
    end
    return
end

function split!(
    c::Union{GroundwaterColumn,ConsolidationColumn,OxidationColumn},
    level::Float,
)
    index, lowerΔz, upperΔz = find_split_index(c.z, c.Δz, level)
    newlength = length(c.z) + 1
    if index != 0
        columnsplit!(column, index, newlength, lowerΔz, upperΔz)
    end
    return
end

# 2022
function columnsplit!(column::DarcyColumn, index, newlength) end
