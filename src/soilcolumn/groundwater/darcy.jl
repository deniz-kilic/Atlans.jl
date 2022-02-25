"""
Finite difference vertical 1D groundwater flow.
"""

struct NumericalSolution
    ϕ::Vector{Float}  # head [m]
    ϕ_old::Vector{Float}  # head of previous timestep [m]
    ϕ_previous::Vector{Float}  # previous head-guess
    A::SymTridiagonal
    rhs::Vector
    maxiter::Int
    close::Float
end

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
Formulate system of equations for groundwater flow.

Central in space, backward in time
"""
function formulate!(column::DarcyColumn, solution::NumericalSolution, Δt::Float)
    solution.rhs .= 0.0
    # dv is diagonal value
    solution.A.dv[1] = -column.conductance[1]
    solution.A.dv[2:end-1] .= -2 * column.conductance[2:end]
    solution.A.dv[end] = -column.conductance[1]
    # ev is off-diagonal value
    solution.A.ev .= column.conductance

    # Add storage component
    # SS: storage coefficient for entire cell
    if Δt > 0.0
        solution.A.dv .-= (column.SS / Δt)
        solution.rhs .-= (column.SS .* column.ϕ_old / Δt)
    end

    # Set boundary conditions 
    n = length(column.ϕ)
    for (i, ϕ) in zip(column.boundary, column.boundary_ϕ)

        # Fix head in cell
        solution.A.dv[i] = 1.0
        solution.rhs[i] = ϕ

        # Disconnect neighbor and add to rhs; conductance is already included
        # in diagonal above.
        if i > 1
            solution.A.ev[i-1] = 0.0
            solution.rhs[i-1] -= column.conductance[i-1] * ϕ
        end
        if i < n
            solution.A.ev[i] = 0.0
            solution.rhs[i+1] -= column.conductance[i] * ϕ
        end

    end
    return
end



"""
Formulate system of equations for abc-isotache.

Contains a (linear) storage term, and a (linearized) creep flux.
"""
function formulate!(column::ConsolidationColumn, solution::NumericalSolution, Δt::Float)
    for i = 1:length(column.cells)
        cell = column.cells[i]
        hcof, rhs = formulate(cell, solution.ϕ[i], column.z[i], column.σ[i], Δt)
        solution.A.dv[i] += hcof
        solution.rhs[i] -= rhs
    end
end


"""
Formulate system of equations for abc-isotache.
"""
function formulate__!(column::ConsolidationColumn, solution::NumericalSolution, Δt::Float)
    for i = 1:length(column.cells)
        cell = column.cells[i]
        σ′ = column.σ[i] - γ_water * (solution.ϕ[i] - column.z[i])
        solution.rhs[i] -= cell.Δz * Qcreep(cell, σ′, Δt)
    end
end


function linearsolve!(solution::NumericalSolution)
    ldiv!(solution.ϕ, solution.A, solution.rhs)
    return
end


function solve!(column::SoilColumn, solution::NumericalSolution, Δt::Float)
    column.ϕ_old .= column.ϕ

    iter = 1
    Δϕmax = Inf
    while iter <= solution.maxiter && Δϕmax > close
        solution.ϕ_previous .= solution.ϕ
        formulate!(column.groundwater, solution, Δt)
        formulate__!(column.consolidation, solution, Δt)
        linearsolve!(solution)
        Δϕmax = max(abs.(solution.ϕ - solution.ϕ_previous))
        iter += 1
    end
    return
end


function flow!(column::DarcyColumn, Δt::Float)
    column.ϕ_old .= column.ϕ
    formulate!(column, Δt)
    linearsolve!(column)
    pore_pressure!(column)
    return
end


function initialize(::DarcyColumn, reader, I)
    use = domain.use
    lithology = ncread3d(reader, :lithology, I)
    geology = ncread3d(reader, :geology, I)

    k = @view fetch_field(reader, :conductivity, I, lithology, geology)[use]
    SS = @view fetch_field(reader, :specific_storage, I, lithology, geology)[use]
    confined = @view fetch_field(reader, :confined, I, lithology, geology)[use]

    column = DarcyColumn(k[domain.index], domain.z, domain.Δz, SS[domain.index])
end
