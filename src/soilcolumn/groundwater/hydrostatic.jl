abstract type HydrostaticColumn <: GroundwaterColumn end


mutable struct Phreatic
    ϕ::Float
end

struct HydrostaticGroundwater <: HydrostaticColumn
    z::Vector{Float}
    phreatic::Phreatic
    dry::Vector{Bool}
    p::Vector{Float}
end

function HydrostaticGroundwater(z, phreatic_level)
    n = length(z)
    return HydrostaticGroundwater(z, Phreatic(phreatic_level), fill(false, n), fill(NaN, n))
end

function set_phreatic_difference!(hg::HydrostaticGroundwater, Δϕ)
    hg.phreatic.ϕ = hg.phreatic.ϕ + Δϕ
    return
end

function set_phreatic!(hg::HydrostaticGroundwater, ϕ)
    hg.phreatic.ϕ = ϕ
    return
end


function pore_pressure!(hg::HydrostaticColumn)
    @. hg.p = hg.phreatic.ϕ - hg.z
    @. hg.dry = hg.p < 0.0
    hg.p[hg.dry] .= 0.0
    return
end


flow!(hg::HydrostaticColumn, _) = pore_pressure!(hg)
phreatic_level(hg::HydrostaticColumn) = hg.phreatic.ϕ


function initialize(
    ::Type{HydrostaticGroundwater},
    domain,
    subsoil,
    I,
)::HydrostaticGroundwater
    return HydrostaticGroundwater(domain.z, subsoil.data[:phreatic_level][I])
end
