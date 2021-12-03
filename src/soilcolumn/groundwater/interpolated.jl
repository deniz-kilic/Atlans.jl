"""
Linearly interpolated heads.
"""
struct InterpolatedGroundwater <: GroundwaterColumn
    z::Vector{Float}  # vertical midpoint coordinate [m]
    Δz::Vector{Float}  # cell height [m]
    boundary::Vector{Int}  # location of head boundaries
    boundary_ϕ::Vector{Float}  # head of boundaries
    dry::Vector{Bool}
    ϕ::Vector{Float}
    p::Vector{Float}
end

function InterpolatedGroundwater(
    z::Vector{Float},  # vertical midpoint coordinate [m]
    Δz::Vector{Float},  # cell height [m]
    boundary::Vector{Int},  # location of head boundaries
    boundary_ϕ::Vector{Float},  # head of boundaries
)
    # TODO: check length of other arrays?
    nlayer = length(z)
    dry = fill(false, nlayer)
    ϕ = fill(NaN, nlayer)
    p = fill(NaN, nlayer)
    return InterpolatedGroundwater(z, Δz, boundary, boundary_ϕ, dry, ϕ, p)
end

function interpolate_head!(ig::InterpolatedGroundwater)
    nlayer = length(ig.Δz)
    nbound = length(ig.boundary)
    z1 = ig.z[ig.boundary[1]]
    z2 = ig.z[ig.boundary[2]]
    ϕ1 = ig.boundary_ϕ[1]
    ϕ2 = ig.boundary_ϕ[2]
    Δz = z2 - z1
    Δϕ = ϕ2 - ϕ1
    j = 2
    z2 = ig.z[ig.boundary[j]]
    ϕ2 = ig.boundary_ϕ[j]
    i = 1
    while i < nlayer
        zi = ig.z[i]
        if zi <= z1
            ig.ϕ[i] = ϕ1
            i += 1
        elseif zi <= z2
            w = (zi - z1) / Δz
            ig.ϕ[i] = ϕ1 + w * Δϕ
            i += 1
        elseif j < nbound
            j += 1
            z1 = z2
            ϕ1 = ϕ2
            z2 = ig.z[ig.boundary[j]]
            ϕ2 = ig.boundary_ϕ[j]
        else
            break
        end
    end
    ig.ϕ[i:end] .= ϕ2
    return
end

function set_phreatic_level!(ig::InterpolatedGroundwater, level::Float)
    for index in eachindex(ig.z)
        zbot = ig.z[index] - 0.5 * ig.Δz[index]
        if level > zbot
            ig.boundary[end] = index
            ig.ϕ[end] = level
        end
        break
    end
end

function flow!(ig::InterpolatedGroundwater, Δt::Float)
    if all(ig.boundary_ϕ == first(ig.boundary_ϕ))
        ig.ϕ .= first(ig.boundary_ϕ)
        return
    else
        interpolate_head!(ig)
    end
    return
end
