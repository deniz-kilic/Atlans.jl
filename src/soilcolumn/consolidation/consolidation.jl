abstract type AbstractConsolidationColumn end


struct PreOverburdenPressure <: Preconsolidation
    pressure::Vector{Float}
end

struct OverConsolidationRatio <: Preconsolidation
    ratio::Vector{Float}
end

struct ConsolidationColumn{C, P} <: AbstractConsolidationColumn
    cells::Vector{C}
    z::Vector{Float}  # cell center
    Δz::Vector{Float}  # cell height
    σ::Vector{Float}  # Total stress in Pa
    σ′::Vector{Float}  # Effective stress in Pa
    p::Vector{Float}  # Pore pressure in Pa
    preconsolidation::P
    result::Vector{Float}
end

function ConsolidationColumn(cells, z, Δz, preconsolidation)
    n = length(cells)
    return ConsolidationColumn(
        cells,
        z,
        Δz,
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        preconsolidation,
        fill(NaN, n),
    )
end

NullConsolidationColumn = ConsolidationColumn{NullConsolidation, OverConsolidationRatio}
apply_preconsolidation!(::NullConsolidationColumn) = nothing
total_stress!(::NullConsolidationColumn, _) = nothing
effective_stress!(::NullConsolidationColumn) = nothing
transfer_stress!(::NullConsolidationColumn) = nothing
prepare_forcingperiod!(::NullConsolidationColumn) = nothing
consolidate!(::NullConsolidationColumn, _, _) = nothing
synchronize_z!(::NullConsolidationColumn, _) = nothing


function apply_preconsolidation!(
    column::ConsolidationColumn{
        ABC,
        OverConsolidationRatio,
    } where {ABC <: AbstractAbcIsotache},
)
    for (i, cell) in enumerate(column.cells)
        column.cells[i] = set_τ0_ocr(cell, column.preconsolidation.ratio[i])
    end
    return
end


function apply_preconsolidation!(
    column::ConsolidationColumn{
        ABC,
        PreOverburdenPressure,
    } where {ABC <: AbstractAbcIsotache},
)
    for (i, cell) in enumerate(column.cells)
        column.cells[i] = set_τ0_pop(cell, column.preconsolidation.pressure[i])
    end
    return
end

#function apply_preconsolidation!(column::ConsolidationColumn{DrainingKoppejan, OverConsolidationRatio})
#    for (i, cell) in enumerate column.cells 
#        column.cells[i] = set_σ′pre_ocr(cell, column.preconsolidation.ratio[i])
#    end
#    return
#end
#
#function apply_preconsolidation!(column::ConsolidationColumn{DrainingKoppejan, PreOverburdenPressure})
#    for (i, cell) in enumerate column.cells 
#        column.cells[i] = set_σ′pre_pop(cell, column.preconsolidation.pressure[i])
#    end
#    return
#end
#
#function plot(column::ConsolidationColumn)
#    plot(column.σ, column.z, color = :gray, label = "σ")
#    plot!(column.p, column.z, color = :blue, label = "p")
#    plot!(column.σ′, column.z, color = :red, label = "σ′")
#    ybot = column.z .- 0.5 .* column.Δz
#    ytop = column.z[end] + 0.5 * column.Δz[end]
#    hline!(ybot, color = :black, label = "")  # bottom
#    hline!([ytop], color = :black, label = "") # top
#end

"""
Consolidation reduces pore space, pushes out the water.
"""
function compress_γ_wet(cell::C where {C <: ConsolidationProcess}, consolidation::Float)
    old_Δz = cell.Δz + consolidation
    return (cell.γ_wet * old_Δz - consolidation * γ_water) / cell.Δz
end

"""
Consolidation reduces pore space, pushes out the air.
"""
function compress_γ_dry(cell::C where {C <: ConsolidationProcess}, consolidation::Float)
    old_Δz = cell.Δz + consolidation
    return (cell.γ_dry * old_Δz) / cell.Δz
end

"""
Terzaghi, degree of consolidation
"""
function U(cell::C where {C <: ConsolidationProcess}, t::Float)
    t_factor = (cell.c_v * t) / (cell.Δz_0 * cell.c_d)^2
    t_pow3 = t_factor^3
    return (t_pow3 / (t_pow3 + 0.5))^(1.0 / 6.0)
end

"""
Weight of (part of) a single cell
"""
function weight(phreatic_level::Float, zbot::Float, Δz::Float, γ_wet::Float, γ_dry::Float)
    ztop = zbot + Δz
    if phreatic_level >= ztop
        wet_Δz = Δz
        dry_Δz = 0.0
    elseif phreatic_level <= zbot
        wet_Δz = 0.0
        dry_Δz = Δz
    else
        dry_Δz = ztop - phreatic_level
        wet_Δz = phreatic_level - zbot
    end
    weight = wet_Δz * γ_wet + dry_Δz * γ_dry
    return weight
end

function submerged_weight(column::AbstractConsolidationColumn, phreatic)
    surface = column.z[end] + 0.5 * column.Δz[end]
    return (max(0.0, phreatic - surface) * γ_water)
end

"""
Compute total stress for entire column
"""
function total_stress!(column::AbstractConsolidationColumn, phreatic_level)
    cumulative_weight = submerged_weight(column, phreatic_level)

    # Iterate from top to bottom
    for index in reverse(eachindex(column.cells))
        cell = column.cells[index]
        if cell.Δz == 0.0
            continue
        end
        zmid = column.z[index]
        zbot = zmid - 0.5 * cell.Δz
        midweight = weight(phreatic_level, zmid, 0.5 * cell.Δz, cell.γ_wet, cell.γ_dry)

        column.σ[index] = cumulative_weight + midweight

        cumulative_weight += weight(phreatic_level, zbot, cell.Δz, cell.γ_wet, cell.γ_dry)
    end
end

"""
Compute effective stress for entire column
"""
function effective_stress!(column::AbstractConsolidationColumn)
    @. column.σ′ = column.σ - column.p
    column.σ′[column.σ′ .< 0.0] .= 0
    return
end

"""
Transfer computed stress to the cells of the ConsolidationColumn.
"""
function transfer_stress!(column::AbstractConsolidationColumn)
    for (i, cell) in enumerate(column.cells)
        newcell = @set cell.σ′ = column.σ′[i]
        column.cells[i] = newcell
    end
end



function consolidate!(column::ConsolidationColumn, phreatic_level, Δt)
    # Compute the new effective stress.
    total_stress!(column, phreatic_level)
    # Pore pressure has been computed by the groundwater column.
    effective_stress!(column)
    # Using the change in effective stress, compute consolidation.
    for (index, cell) in enumerate(column.cells)
        σ′ = column.σ′[index]
        newcell = consolidate(cell, σ′, Δt)
        column.cells[index] = newcell
        column.result[index] = newcell.consolidation
    end
end

function synchronize_z!(
    column::ConsolidationColumn{C},
    Δz,
) where {C <: ConsolidationProcess}
    for i in 1:length(column.cells)
        cell = column.cells[i]
        newcell = @set cell.Δz = Δz[i]
        column.cells[i] = newcell
    end
    return
end

function update_γ!(
    column::ConsolidationColumn{C},
    shrinkage,
) where {C <: ConsolidationProcess}
    for i in 1:length(column.cells)
        cell = column.cells[i]
        γ_wet = compress_γ_wet(cell, shrinkage[i] + cell.consolidation)
        γ_dry = compress_γ_dry(cell, shrinkage[i] + cell.consolidation)
        newcell = @set cell.γ_dry = γ_dry
        newcell = @set newcell.γ_wet = γ_wet
        column.cells[i] = newcell
    end
    return
end


update_γ!(column::ConsolidationColumn{NullConsolidation}, shrinkage) = nothing
