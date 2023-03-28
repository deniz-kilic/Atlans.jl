"""
    ShrinkageColumn{S}(cells, z, Δz, result, max_shrinkage_depth)

Collection of SimpleShrinkage cells to compute shrinkage for. 

# Arguments:
- `cells::Vector{S}`: Collection of SimpleShrinkage structs.
- `z::Vector{Float}`: Depth of the cells.
- `Δz::Vector{Float}`: Thickness of the cells.
- `result::Vector{Float}`: Computed shrinkage of each cell. All NaNs at t=0.
- `max_shrinkage_depth::Float`: Maximum depth to compute shrinkage for in cells.
"""
struct ShrinkageColumn{S}
    cells::Vector{S}
    z::Vector{Float}  # cell center
    Δz::Vector{Float}  # cell height
    result::Vector{Float}
    max_shrinkage_depth::Float
end


shrinkage_depth(column::ShrinkageColumn{NullShrinkage}, _, _, _, _) = nothing
shrink!(column::ShrinkageColumn{NullShrinkage}, _) = nothing
synchronize_z!(column::ShrinkageColumn{NullShrinkage}, __) = nothing


function shrink!(
    column::ShrinkageColumn{S},
    phreatic_level::Float,
    Δt::Float,
) where {S<:ShrinkageProcess}
    shrinkage_z = max(
        phreatic_level,
        surface_level(column) - column.max_shrinkage_depth
    )
    column.result .= 0.0
    for index in reverse(1:length(column.cells))
        if column.z[index] > shrinkage_z
            cell = column.cells[index]
            newcell = shrink(cell, Δt)
            column.cells[index] = newcell
            column.result[index] = newcell.shrinkage
        else
            break
        end
    end
    return
end


function synchronize_z!(column::ShrinkageColumn{S}, Δz) where {S<:ShrinkageProcess}
    for i = 1:length(column.cells)
        cell = column.cells[i]
        newcell = @set cell.Δz = Δz[i]
        column.cells[i] = newcell
    end
    return
end


function initialize(::Type{SimpleShrinkage}, domain, subsoil, I)
    n = fetch_field(shrinkage_degree, :, I, domain) #TODO: what is fieldname for shrinkage factor in parameter table?
    max_shrinkage_depth = subsoil.data[:phreatic_level][I] + 0.3 #TODO: is prheatic level in nc correct to use for GLG?
    τ = 60.0 # Time dependent factor for shrinkage process
    r = 3.0 # Direction of shrinkage is assumed isoptropic (r=3)

    cells = Vector{SimpleShrinkage}()

    for (i, Δz) in enumerate(domain.Δz)
        cell = SimpleShrinkage(Δz, n[i], τ, r, NaN)
        push!(cells, cell)
    end

    column = ShrinkageColumn(
        cells,
        domain.z,
        domain.Δz,
        fill(NaN, domain.n),
        max_shrinkage_depth
    )
    return column

end