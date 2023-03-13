"""
Logic for splitting a column to accomodate for a moving phreatic level in
combination with organic matter stores.
"""
struct AdaptiveCellsize
    Δzmax::Float
    split_tolerance::Float
end

function find_split_index(z, Δz, level, tolerance)
    for index in eachindex(z)
        zbot = z[index] - 0.5 * Δz[index]
        ztop = z[index] + 0.5 * Δz[index]
        if level < ztop && level > zbot
            lowerΔz = level - zbot
            upperΔz = ztop - level
            if (lowerΔz < tolerance) || (upperΔz < tolerance)
                return nothing, NaN, NaN
            else
                return index, lowerΔz, upperΔz
            end
        end
    end
    return nothing, NaN, NaN
end

function shouldsplit(vector, newlength)
    if length(vector) == newlength
        return false
    elseif (length(vector) + 1) == newlength
        return true
    else
        error("vector is not equal to newlength or newlength - 1")
    end
end

function cellsplit!(
    column::Union{ConsolidationColumn,OxidationColumn,ShrinkageColumn},
    index,
    newlength,
    lowerΔz,
    upperΔz,
)
    if shouldsplit(column.cells, newlength)
        # @set also creates a copy
        cell = column.cells[index]
        upper = @set cell.Δz = upperΔz
        lower = @set cell.Δz = lowerΔz
        insert!(column.cells, index, lower)
        column.cells[index+1] = upper
    end
    return
end

"""
For CarbonStore, organic and mineral mass should be split according to cell
height (Δz).
"""
function cellsplit!(
    column::OxidationColumn{CarbonStore},
    index,
    newlength,
    lowerΔz,
    upperΔz,
)
    if shouldsplit(column.cells, newlength)
        # @set also creates a copy
        cell = column.cells[index]
        lower_fraction = lowerΔz / cell.Δz
        upper_fraction = upperΔz / cell.Δz

        lower = CarbonStore(
            lowerΔz,
            cell.f_organic,
            cell.f_minimum_organic,
            lower_fraction * cell.m_organic,
            lower_fraction * cell.m_mineral,
            cell.α,
            NaN,
        )
        upper = CarbonStore(
            upperΔz,
            cell.f_organic,
            cell.f_minimum_organic,
            upper_fraction * cell.m_organic,
            upper_fraction * cell.m_mineral,
            cell.α,
            NaN,
        )

        insert!(column.cells, index, lower)
        column.cells[index+1] = upper
    end
    return
end

function zsplit!(Δz, index, newlength, lowerΔz, upperΔz)
    if shouldsplit(Δz, newlength)
        insert!(Δz, index, lowerΔz)
        Δz[index+1] = upperΔz
    end
    return
end

function columnsplit!(hg::HydrostaticGroundwater, index, newlength, lowerΔz, upperΔz)
    push!(hg.dry, false)
    push!(hg.p, NaN)
end

function columnsplit!(cc::ConsolidationColumn, index, newlength, lowerΔz, upperΔz)
    cellsplit!(cc, index, newlength, lowerΔz, upperΔz)
    push!(cc.σ, NaN)
    push!(cc.σ′, NaN)
    push!(cc.p, NaN)
    push!(cc.result, NaN)
end

function columnsplit!(oc::OxidationColumn, index, newlength, lowerΔz, upperΔz)
    cellsplit!(oc, index, newlength, lowerΔz, upperΔz)
    push!(oc.result, NaN)
end

function columnsplit!(sc::ShrinkageColumn, index, newlength, lowerΔz, upperΔz)
    cellsplit!(sc, index, newlength, lowerΔz, upperΔz)
    push!(sc.result, NaN)
end

function split!(sc::SoilColumn, level, tolerance)
    index, lowerΔz, upperΔz = find_split_index(sc.z, sc.Δz, level, tolerance)
    isnothing(index) && return
    newlength = length(sc.z) + 1
    if index != 0
        zsplit!(sc.Δz, index, newlength, lowerΔz, upperΔz)
        push!(sc.z, NaN)
        push!(sc.subsidence, NaN)
        columnsplit!(sc.groundwater, index, newlength, lowerΔz, upperΔz)
        columnsplit!(sc.consolidation, index, newlength, lowerΔz, upperΔz)
        columnsplit!(sc.oxidation, index, newlength, lowerΔz, upperΔz)
        columnsplit!(sc.shrinkage, index, newlength, lowerΔz, upperΔz)
    end
    update_z!(sc)
    return
end
