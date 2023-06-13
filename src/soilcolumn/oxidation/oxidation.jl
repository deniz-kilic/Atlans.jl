struct OxidationColumn{O}
    cells::Vector{O}
    z::Vector{Float}
    Δz::Vector{Float}
    result::Vector{Float}
    max_oxidation_depth::Float
    Hv0::Float
end


oxidation_depth(column::OxidationColumn{NullOxidation}, _, _, _, _) = nothing
oxidate!(column::OxidationColumn{NullOxidation}, _, _) = nothing
synchronize_z!(column::OxidationColumn{NullOxidation}, _) = nothing


function oxidation_depth(
    column::OxidationColumn{O},
    surface_level,
    phreatic_level,
    deep_subsidence,
    phreatic_change,
) where {O<:OxidationProcess}
    new_surface = surface_level - deep_subsidence
    new_phreatic = phreatic_level + phreatic_change
    oxidation_z = new_phreatic + column.Hv0 # solution
    depth = max(0.0, new_surface - oxidation_z)
    depth = min(column.max_oxidation_depth, depth)
    return surface_level - depth
end


function oxidate!(
    column::OxidationColumn{O},
    phreatic_level::Float,
    Δt::Float,
) where {O<:OxidationProcess}
    oxidation_z = max( # Doesn't matter if oxidation_z is above surface_level of column
        phreatic_level + column.Hv0,
        surface_level(column) - column.max_oxidation_depth
    )
    column.result .= 0.0
    for index in reverse(1:length(column.cells))
        if column.z[index] > oxidation_z
            cell = column.cells[index]
            newcell = oxidate(cell, Δt)
            column.cells[index] = newcell
            column.result[index] = newcell.oxidation
        else
            break
        end
    end
    return
end


function synchronize_z!(column::OxidationColumn{O}, Δz) where {O<:OxidationProcess}
    for i = 1:length(column.cells)
        cell = column.cells[i]
        newcell = @set cell.Δz = Δz[i]
        column.cells[i] = newcell
    end
    return
end


"""
    initialize(::Type{CarbonStore}, domain, subsoil, I)

Initialize a OxidationColumn for a domain at location I based subsurface input.
"""
function initialize(::Type{CarbonStore}, domain, subsoil, I)
    f_organic = fetch_field(subsoil, :mass_fraction_organic, I, domain)
    f_minimum_organic = fetch_field(subsoil, :minimal_mass_fraction_organic, I, domain)
    α = fetch_field(subsoil, :oxidation_rate, I, domain)
    ρb = fetch_field(subsoil, :rho_bulk, I, domain)

    Hv0 = subsoil.data[:Hv0_ox]
    max_oxidation_depth = subsoil.data[:max_oxidation_depth][I]

    cells = Vector{CarbonStore}()
    for (i, Δz) in enumerate(domain.Δz)
        cell = CarbonStore(Δz, f_organic[i], f_minimum_organic[i], ρb[i], α[i])
        push!(cells, cell)
    end

    column = OxidationColumn(
        cells,
        domain.z,
        domain.Δz,
        fill(0.0, domain.n),
        max_oxidation_depth,
        Hv0,
    )
    return column
end


"""
    initialize(::Type{SimpleShrinkage}, domain, subsoil, I)

Initialize an empty OxidationColumn (i.e. oxidation is ignored) at location I.
"""
function initialize(::Type{NullOxidation}, domain, _, _)
    cells = Vector{NullOxidation}
    for i in 1:length(domain.Δz)
        push!(cells, NullOxidation())
    end
    return OxidationColumn(
        cells,
        domain.z,
        domain.Δz,
        fill(0.0, domain.n),
        NaN,
        NaN
    )
end
