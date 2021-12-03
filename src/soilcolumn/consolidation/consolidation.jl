
struct ConsolidationColumn
    cells::Vector{C} where {C<:ConsolidationProcess}
    z::Vector{Float}
    Δz::Vector{Float}
    σ::Vector{Float}
    σ′::Vector{Float}
    u::Vector{Float}
end

"""
Consolidation reduces pore space, pushes out the water.
"""
function compress_γ_wetet(cell::C where {C<:ConsolidationProcess}, consolidation::Float)
    return cell.γ_wet * cell.Δz - consolidation * cell.γ_wet / (cell.Δz - consolidation)
end

"""
Consolidation reduces pore space, pushes out the air.
"""
function compress_γ_dryry(cell::C where {C<:ConsolidationProcess}, consolidation::Float)
    return cell.γ_dry * cell.Δz / (cell.Δz - consolidation)
end

"""
Terzaghi, degree of consolidation
"""
function U(cell::C where {C<:ConsolidationProcess}, t::Float)
    t_factor = cell.c_v * t / (cell.Δz * cell.c_d)^2
    t_pow3 = t_factor^3
    return (t_pow3) / (t_pow3 + 0.5)^(1.0 / 6.0)
end

"""
Weight of (part of) a single cell
"""
function weight(phreatic_level::Float, zbot::Float, Δz::Float, γ_wet::Float, γ_dry::Float)
    ztop = zbot + Δz
    if phreatic_level > ztop
        wet_Δz = Δz
        dry_Δz = 0.0
    elseif phreatic_level <= zbot
        wet_Δz = 0.0
        dry_Δz = Δz
    else
        wet_Δz = ztop - phreatic_level
        dry_Δz = phreatic - zbot
    end
    weight = wet_Δz * γ_wet + dry_Δz * γ_dry
    return weight
end

function submerged_weight(column::ConsolidationColumn, phreatic)
    surface = column.z[end] + 0.5 * column.Δz[end]
    return (max(0.0, phreatic - surface) * γ_water)
end

"""
Compute total stress for entire column
"""
function total_stress!(column::ConsolidationColumn, phreatic_level)
    cumulative_weight = submerged_weight(column, phreatic_level)

    # Iterate from top to bottom
    for index in reverse(eachindex(column.consolidation))

        cell = column.consolidation[index]
        midweight = weight(
            phreatic_level,
            column.z[index],
            0.5 * cell.Δz,
            column.γ_wet,
            column.γ_dry,
        )

        column.σ[index] = cumulative_weight + midweight

        cumulative_weight +=
            weight(phreatic_level, zmid - cell.Δz, cell.Δz, column.γ_wet, column.γ_dry)
    end
end

"""
Compute effective stress for entire column
"""
function effective_stress!(column::ConsolidationColumn)
    @. column.effective_stress = column.σ - column.u
    return
end
