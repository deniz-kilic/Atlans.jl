"""
    ShrinkageColumn{S}(cells, z, Δz, result, max_shrinkage_depth)

Collection of SimpleShrinkage cells to compute shrinkage for. 

# Arguments:
- `cells::Vector{S}`: Collection of structs for SimpleShrinkage cell attributes.
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
shrink!(column::ShrinkageColumn{NullShrinkage}, _, _) = nothing
synchronize_z!(column::ShrinkageColumn{NullShrinkage}, _) = nothing


function shrinkage_depth(
    column::ShrinkageColumn{S},
    surface_level,
    phreatic_level,
    deep_subsidence,
    phreatic_change,
) where {S<:ShrinkageProcess}
    new_surface = surface_level - deep_subsidence
    new_phreatic = phreatic_level + phreatic_change
    depth = max(0.0, new_surface - new_phreatic)
    depth = min(column.max_shrinkage_depth, depth)
    return surface_level - depth
end


function shrink!(
    column::ShrinkageColumn{S},
    phreatic_level::Float,
    Δt::Float,
) where {S<:ShrinkageProcess}
    shrinkage_z = max(phreatic_level, surface_level(column) - column.max_shrinkage_depth)
    column.result .= 0.0
    for index in reverse(1:length(column.cells))
        if column.z[index] > shrinkage_z
            cell = column.cells[index] # SimpleShrinkage cell
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
    n = fetch_field(subsoil, :shrinkage_degree, I, domain)
    L = fetch_field(subsoil, :mass_fraction_lutum, I, domain)
    H = fetch_field(subsoil, :mass_fraction_organic, I, domain)

    max_shrinkage_depth = subsoil.data[:max_shrinkage_depth][I]
    cells = Vector{SimpleShrinkage}()

    for (i, Δz) in enumerate(domain.Δz)
        cell = SimpleShrinkage(Δz, n[i], L[i], H[i])
        push!(cells, cell)
    end

    column = ShrinkageColumn(
        cells,
        domain.z,
        domain.Δz,
        fill(NaN, domain.n),
        max_shrinkage_depth,
    )
    return column
end
