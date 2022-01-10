"""
Finite difference vertical 1D groundwater flow.
"""
struct DarcyColumn <: GroundwaterColumn
    k::Vector  # vertical conductivity [m/d]
    z::Vector  # vertical coordinate [m]
    Δz::Vector  # cell height [m]
    SS::Vector  # specific storage(?) [m/m]
    boundary::Vector{Int}  # location of head boundaries
    boundary_ϕ::Vector  # head of boundaries
    confined::Vector{Bool}
    dry::Vector{Bool}
    # Intermediate
    conductance::Vector{Float}  # [m/d]
    ϕ::Vector{Float}  # head [m]
    ϕ_old::Vector{Float}  # head of previous timestep [m]
    p::Vector{Float}  # pore pressure [m]
    # Numerics
    A::SymTridiagonal
    rhs::Vector
end

"""Internodal conductance"""
function harmonicmean_conductance(k1, k2, Δz1, Δz2)
    csum = (0.5 * Δz1) / k1 + (0.5 * Δz2) / k2
    return 1.0 / csum
end

function conductance!(column::DarcyColumn)
    for i = 1:(length(column.ϕ)-1)
        column.conductance[i] = harmonicmean_conductance(
            column.k[i],
            column.k[i+1],
            column.Δz[i],
            column.Δz[i+1],
        )
    end
    return
end

function update_boundaries!(column::DarcyColumn, boundary::Vector{Int}, ϕ::Vector{Float})
    column.boundary .= boundary
    column.boundary_ϕ .= ϕ
    return
end


"""
Central in space, backward in time
"""
function formulate!(column::DarcyColumn, Δt::Float)
    column.rhs .= 0.0
    # dv is diagonal value
    column.A.dv[1] = -column.conductance[1]
    column.A.dv[2:end-1] .= -2 * column.conductance[2:end]
    column.A.dv[end] = -column.conductance[1]
    # ev is off-diagonal value
    column.A.ev .= column.conductance

    # Add storage component
    # SS: storage coefficient for entire cell
    if Δt > 0.0
        column.A.dv .-= (column.SS / Δt)
        column.rhs .-= (column.SS .* column.ϕ_old / Δt)
    end

    # Set boundary conditions 
    n = length(column.ϕ)
    for (i, ϕ) in zip(column.boundary, column.boundary_ϕ)

        # Fix head in cell
        column.A.dv[i] = 1.0
        column.rhs[i] = ϕ

        # Disconnect neighbor and add to rhs; conductance has already included
        # in diagonal above.
        if i > 1
            column.A.ev[i-1] = 0.0
            column.rhs[i-1] -= column.conductance[i-1] * ϕ
        end
        if i < n
            column.A.ev[i] = 0.0
            column.rhs[i+1] -= column.conductance[i] * ϕ
        end

    end
    return
end


function linearsolve!(column::DarcyColumn)
    ldiv!(column.ϕ, column.A, column.rhs)
    return
end


function flow!(column::DarcyColumn, Δt::Float)
    column.ϕ_old .= column.ϕ
    formulate!(column, Δt)
    linearsolve!(column)
end
