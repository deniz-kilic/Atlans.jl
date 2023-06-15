"""
    ShrinkageColumn{S}(cells, z, Δz, result, Hv0)

Collection of SimpleShrinkage cells to compute shrinkage for. 

# Arguments:
- `cells::Vector{S}`: Collection of structs for SimpleShrinkage cell attributes.
- `z::Vector{Float}`: Depth of the cells.
- `Δz::Vector{Float}`: Thickness of the cells.
- `result::Vector{Float}`: Computed shrinkage of each cell. All NaNs at t=0.
- `Hv0::Float`: Absolute depth above phreatic level to compute shrinkage for in cells.
"""
struct ShrinkageColumn{S}
    cells::Vector{S}
    z::Vector{Float}  # cell center
    Δz::Vector{Float}  # cell height
    result::Vector{Float}
    no_shrinkage_Δz::Float
end


shrinkage_level(column::ShrinkageColumn{NullShrinkage}, _, _) = nothing
shrink!(column::ShrinkageColumn{NullShrinkage}, _, _) = nothing
synchronize_z!(column::ShrinkageColumn{NullShrinkage}, _) = nothing


function shrinkage_level(
    column::ShrinkageColumn{S},
    phreatic_level,
    phreatic_change,
) where {S<:ShrinkageProcess}
    shrinkage_z = phreatic_level + phreatic_change + column.no_shrinkage_Δz
    return shrinkage_z
end


function shrink!(
    column::ShrinkageColumn{S},
    phreatic_level::Float,
    Δt::Float,
) where {S<:ShrinkageProcess}
    shrinkage_z = phreatic_level + column.no_shrinkage_Δz # Doesn't matter is this is above surface level
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


"""
    initialize(::Type{SimpleShrinkage}, domain, subsoil, I)

Initialize a ShrinkageColumn for a domain at location I based subsurface input.
"""
function initialize(::Type{SimpleShrinkage}, domain, subsoil, I)
    n = fetch_field(subsoil, :shrinkage_degree, I, domain)
    L = fetch_field(subsoil, :mass_fraction_lutum, I, domain)
    H = fetch_field(subsoil, :mass_fraction_organic, I, domain)

    no_shrinkage_Δz = subsoil.data[:no_shr_thickness][I]
    cells = Vector{SimpleShrinkage}()

    for (i, Δz) in enumerate(domain.Δz)
        cell = SimpleShrinkage(Δz, n[i], L[i], H[i])
        push!(cells, cell)
    end

    column = ShrinkageColumn(
        cells,
        domain.z,
        domain.Δz,
        fill(0.0, domain.n),
        no_shrinkage_Δz,
    )
    return column
end


"""
    initialize(::Type{SimpleShrinkage}, domain, subsoil, I)

Initialize an empty ShrinkageColumn (i.e. shrinkage is ignored) at location I.
"""
function initialize(::Type{NullShrinkage}, domain, _, _)
    cells = Vector{NullShrinkage}()
    for i in 1:length(domain.Δz)
        push!(cells, Atlans.NullShrinkage())
    end
    return ShrinkageColumn(
        cells,
        domain.z,
        domain.Δz,
        fill(0.0, domain.n),
        NaN
    )
end
