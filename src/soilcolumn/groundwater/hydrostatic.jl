mutable struct Phreatic
    ϕ::Float
end

struct HydrostaticGroundwater <: GroundwaterColumn
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


function pore_pressure!(hg::HydrostaticGroundwater)
    @. hg.p = hg.phreatic.ϕ - hg.z
    @. hg.dry = hg.p < 0.0
    hg.p[hg.dry] .= 0.0
    return
end


flow!(hg::HydrostaticGroundwater, _) = pore_pressure!(hg)
phreatic_level(hg::HydrostaticGroundwater) = hg.phreatic.ϕ


function initialize(::HydrostaticGroundwater, domain, reader, I)::HydrostaticGroundwater
    return HydrostaticGroundwater(domain.z, ncread2d(reader, :phreatic_level, I))
end
