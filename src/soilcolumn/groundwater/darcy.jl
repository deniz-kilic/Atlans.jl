"""
Finite difference vertical 1D groundwater flow.
"""
struct DarcyColumn <: GroundwaterColumn
    k::Vector  # vertical conductivity [m/d]
    z::Vector  # vertical coordinate [m]
    Δz::Vector  # cell height [m]
    dry::Vector{Bool}  # unsaturated flag
    SS::Vector  # specific storage(?) [m/m]
    S_ske::Vector  # skeletal storage [m/m]
    boundary::Vector{Int}  # location of head boundaries
    boundary_ϕ::Vector  # head of boundaries
    # Intermediate
    conductance::Vector{Float}  # [m/d]
    ϕ::Vector{Float}  # average head [m]
    p::Vector{Float}  # pore pressure [m]
    # Numerics
    A::SymTridiagonal
    constant::Vector{Bool}
    rhs::Vector
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

function flow!(column::DarcyColumn, Δt::Float)
    formulate!(column)
    solve!(column)
end
