"""
Linearly interpolated heads.
"""

struct InterpolatedGroundwater <: GroundwaterColumn
    z::Vector{Float}  # vertical midpoint coordinate [m]
    Δz::Vector{Float}  # cell height [m]
    phreatic::Float
    aquifer::Float
    dry::Vector{Bool}
    ϕ::Vector{Float}
    p::Vector{Float}
end

function InterpolatedGroundwater(
    z::Vector{Float},  # vertical midpoint coordinate [m]
    Δz::Vector{Float},  # cell height [m]
    phreatic_level,
    aquifer_head,
)
    n = length(z)
    return InterpolatedGroundwater(
        z,
        Δz,
        phreatic_level,
        aquifer_head,
        dry = fill(false, n),
        ϕ = fill(NaN, n),
        p = fill(NaN, n),
    )
end

function interpolate_head!(ig::InterpolatedGroundwater)
    Δϕ = ig.phreatic.ϕ - ig.ϕaquifer.ϕ
    Δz = ig.phreatic.ϕ - ig.z[i]
    ig.ϕ .= ig.ϕ_aquifer + min.(ig.z .- ig.z[1] ./ Δz, 1.0) * Δϕ
    ig.dry .= ig.z > ϕ_top
    return
end

function flow!(ig::InterpolatedGroundwater, Δt::Float)
    interpolate_head!(ig)
    pore_pressure!(ig)
    return
end

function set_aquifer_difference!(ig::InterpolatedGroundwater, Δϕ)
    @set ig.aquifer.ϕ = ig.aquifer.ϕ + Δϕ
    return
end

function set_phreatic_difference!(ig::InterpolatedGroundwater, Δϕ)
    @set ig.phreatic.ϕ = ig.phreatic.ϕ + Δϕ
    return
end

function set_aquifer!(ig::InterpolatedGroundwater, ϕ)
    @set ig.aquifer.ϕ = ϕ
    return
end

function set_phreatic!(ig::InterpolatedGroundwater, ϕ)
    @set ig.phreatic.ϕ = ϕ
    return
end


function initialize(::InterpolatedGroundwater, domain, reader, I)::InterpolatedGroundwater
    return InterpolatedGroundwater(
        domain.z,
        domain.Δz,
        ncread2d(reader, :phreatic_level, I),
        ncread2d(reader, :aquifer_head, I),
        fill(false, domain.n),
        fill(NaN, domain.n)fill(NaN, domain.n),
    )
end

phreatic_level(ig::InterpolatedGroundwater) = ig.phreatic.ϕ::Float
