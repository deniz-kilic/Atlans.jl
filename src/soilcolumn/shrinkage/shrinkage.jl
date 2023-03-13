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
    shrinkage_z = max(phreatic_level, surface_level(column) - column.max_shrinkage_depth)
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


function synchronize!(column::ShrinkageColumn{S}, Δz) where {S<:ShrinkageProcess}
    for i = 1:length(column.cells)
        cell = column.cells[i]
        newcell = @set cell.Δz = Δz[i]
        column.cells[i] = newcell
    end
    return
end


function initialize(::Type{SimpleShrinkage}, domain, subsoil, I)
    n = fetch_field(shrinkage_degree, :, I, domain)
    n_residual = fetch_field(residual_shrinkage_degree, :, I, domain)

    cells = Vector{SimpleShrinkage}()
    for (i, Δz) in enumerate(domain.Δz)
        cell = SimpleShrinkage(Δz, n[i], n_residual[i], NaN)
        push!(cells, cell)
    end

    column = ShrinkageColumn(
        cells,
        domain.z,
        domain.Δz,
        fill(NaN, domain.n),
        2 # TODO: find out how to determine max shrinkage depth 
    )
    return column

end