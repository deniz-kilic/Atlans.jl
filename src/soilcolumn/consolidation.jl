"""
Consolidation reduces pore space, pushes out the water.
"""
function compress_γ_wet(cell::C where {C<:ConsolidationProcess}, consolidation::Float64)
    return cell.γ_w * cell.Δz - consolidation * cell.γ_w / (cell.Δz - consolidation)
end

"""
Consolidation reduces pore space, pushes out the air.
"""
function compress_γ_dry(cell::C where {C<:ConsolidationProcess}, consolidation::Float64)
    return cell.γ_d * cell.Δz / (cell.Δz - consolidation)
end

"""
Terzaghi, degree of consolidation
"""
function U(cell::C where {C<:ConsolidationProcess}, t::Float64)
    t_factor = cell.c_v * t / (cell.Δz * cell.c_d)^2
    t_pow3 = t_factor^3
    return (t_pow3) / (t_pow3 + 0.5)^(1.0 / 6.0)
end
