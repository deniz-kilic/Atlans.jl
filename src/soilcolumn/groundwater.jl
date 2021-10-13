function pore_pressure!(gw::GW where {GW<:GroundwaterProcess})
    @. gw.p = gw.γ_w * (gw.ϕ - (gw.z + 0.5 * gw.Δz))
    gw.p[gw.dry] .= 0.0
end

function plot(gw::GW where {GW<:GroundwaterProcess})
    y = gw.z .+ 0.5 .* Δz
    plot(gw.ϕ, y, color = :blue)
    plot!(gw.p, y, color = :red)
    vline!(gw.z, color = :black)
    vline!(gw.z[end] + gw.Δz[end], color = :black)
end

"""
Linearly interpolated heads.
"""
struct SimpleGroundwater <: GroundwaterProcess
    z::Vector{Float}  # vertical coordinate [m]
    Δz::Vector{Float}  # cell height [m]
    boundary::Vector{Int}  # location of head boundaries
    boundary_ϕ::Vector{Float}  # head of boundaries
    dry::Vector{Bool}
    ϕ::Vector{Float}
    p::Vector{Float}
end

function interpolate_head!(sg::SimpleGroundwater)
    nbound = length(sg.boundary)
    z1 = sg.z[sg.boundary[1]]
    z2 = sg.z[sg.boundary[2]]
    ϕ1 = sg.boundary[1]
    ϕ2 = sg.boundary[2]
    Δz = z2 - z1
    Δϕ = ϕ2 - ϕ1
    j = 2
    z2 = sg.z[sg.boundary[j]]
    ϕ2 = sg.boundary[j]
    i = 1
    while i < nbound
        zi = sg.z[i]
        if zi <= z1
            sg.ϕ[i] = ϕ1
            i += 1
        elseif zi <= z2
            w = (zi - z1) / Δz
            sg.ϕ[i] = ϕ1 + w * Δϕ
            i += 1
        elseif j < nbound
            j += 1
            z1 = z2
            ϕ1 = ϕ2
            z2 = sg.z[sg.boundary[j]]
            ϕ2 = sg.boundary[j]
        else
            break
        end
    end
    sg.ϕ[i:end] .= ϕ2
    return
end

function solve!(sg::SimpleGroundwater)
    if all(sg.boundary_ϕ == first(sg.boundary_ϕ))
        sg.ϕ .= first(sg.boundary_ϕ)
        return
    else
        interpolate_head!(sg)
    end
    return
end


"""
Finite difference vertical 1D groundwater flow.
"""
struct DarcyColumn <: GroundwaterProcess
    k::Vector  # vertical conductivity [m/d]
    z::Vector  # vertical coordinate [m]
    Δz::Vector  # cell height [m]
    dry::Vector{Bool}  # unsaturated flag
    SS::Vector  # specific storage(?) [m/m]
    S_ske::Vector  # skeletal storage [m/m]
    boundary::Vector{Int}  # location of head boundaries
    boundary_ϕ::Vector  # head of boundaries
    γ_w::Float64  # specific weight [kg/m^3]
    # Intermediate
    conductance::Vector{Float}  # [m/d]
    ϕ::Vector{Float}  # average head [m]
    p::Vector{Float}  # pore pressure [m]
    # Numerics
    A::SymTridiagonal
    constant::Vector{Bool}
    rhs::Vector
end

function push!(column::DarcyColumn, values)
    for (field, v) in zip(fieldnames(column), values)
        push!(column.(field), v)
    end
    # update conductance
    return
end

function harmonicmean_conductance(k1, k2, Δz1, Δz2)
    csum = Δz1 / k1 + Δz2 / k2
    return 1.0 / csum
end

function conductance!(column::DarcyColumn)
    k1 = column.kv[1]
    Δz1 = column.Δz[1]
    for i = 2:size(column.ϕ)
        k2 = column.kv[i]
        Δz2 = column.Δz[i]
        column.conductance[i-1] = harmonicmean_conductance(k1, k2, Δz1, Δz2)
    end
    return
end

function update_boundaries!(column::DarcyColumn, boundary::Vector{Int}, ϕ::Vector{Int})
    column.boundary = boundary
    column.boundary_ϕ = ϕ
    column.constant .= false
    for i in boundary
        column.constant[i] = true
    end
    return
end

function formulate!(column::DarcyColumn)
    # A is a symmetric tridiagonal matrix, with attributes
    # d: diagonal
    # ev: super-diagonal
    for i = 1:length(column.ϕ)-1
        if column.constant[i]
            column.A.d[i] = 1.0
            column.rhs[i] = column.ϕ[i]
        end
        j = i + 1
        i_constant = column.constant[i]
        j_constant = column.constant[j]
        conductance = column.conductance[i]
        if i_constant & j_constant
            continue
        elseif j_constant
            column.A.d[i] -= conductance
            column.rhs[i] -= conductance * column.ϕ[j+1]
        elseif i_constant
            column.A.d[j] -= conductance
            column.rhs[j] -= conductance * column.ϕ[i+1]
        else
            column.A.d[i] -= conductance
            column.A.d[j] += conductance
        end
    end
    return
end

function solve!(column::DarcyColumn)
    for (i, ϕ) in zip(column.boundary, column.boundary_ϕ)
        column.ϕ[i] = ϕ
    end
    ldiv!(column.ϕ, column.A, column.rhs)
    return
end
