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

function columnsplit!(hg::HydrostaticGroundwater, index, newlength, lowerΔz, upperΔz)
    push!(ig.dry, true)
    push!(ig.ϕ, NaN)
    push!(ig.p, NaN)
    return
end

function columnsplit!(cc::ConsolidationColumn, index, newlength, lowerΔz, upperΔz)
    zsplit!(cc.z, cc.Δz, index, newlength, lowerΔz, upperΔz)
    cellsplit!(cc, index, newlength, lowerΔz, upperΔz)
    push!(cc.σ, NaN)
    push!(cc.σ′, NaN)
    push!(cc.p, NaN)
end

function columnsplit!(oc::OxidationColumn, index, newlength, lowerΔz, upperΔz)
    zsplit!(oc.z, oc.Δz, index, newlength, lowerΔz, upperΔz)
    cellsplit!(cc, index, newlength, lowerΔz, upperΔz)
end


function split!(
    sc::SoilColumn,
    level,
)
    index, lowerΔz, upperΔz = find_split_index(sc.z, sc.Δz, level)
    newlength = length(c.z) + 1
    if index != 0
        columnsplit!(c.groundwater, index, newlength, lowerΔz, upperΔz)
        columnsplit!(c.consolidation, index, newlength, lowerΔz, upperΔz)
        columnsplit!(c.oxidation, index, newlength, lowerΔz, upperΔz)
    end
    return
end
