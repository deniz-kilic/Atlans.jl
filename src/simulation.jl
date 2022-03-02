struct Reader
    dataset::NCDataset
    params::Dict{Symbol,String}
    times::Vector{DateTime}
end

struct SubsoilReader
    dataset::NCDataset
    tables::Dict{Symbol,Dict{Tuple{Int,Int},Float}}
    params::Dict{Symbol,String}
end

struct Writer
    dataset::NCDataset
    params::Dict{Symbol,String}
end

struct Output
    x::Vector{Float}
    y::Vector{Float}
    phreatic_level::Array{Float}
    consolidation::Array{Float}
    oxidation::Array{Float}
    subsidence::Array{Float}
end

struct Model{G,C,P,O,T}
    columns::Vector{SoilColumn{G,C,P,O}}
    index::Vector{CartesianIndex}
    timestepper::T
    output::Output
end

struct StageIndexation
    percentile::Int
    weir_area::Array{OptionalInt}
    change::Array{OptionalFloat}
    reader::Reader
end

struct DeepSubsidence
    subsidence::Array{OptionalFloat}
    reader::Reader
end

struct StageChange
    change::Array{OptionalFloat}
    reader::Reader
end

struct AquiferHead
    head::Array{OptionalFloat}
    reader::Reader
end

struct Simulation{F}
    model::Model
    clock::Clock
    writer::Writer
    forcing::F
end

# Initialization
"""
Repeat function not in base / stdlib
"""
function repeat_elements(x, v)
    z = similar(x, sum(v))
    i = 0
    for (x, n) in zip(x, v), _ = 1:n
        @inbounds z[i+=1] = x
    end
    return z
end

"""
In how many equal parts should a thick cell be divided?
"""
function discretize(Δz, maxΔz::Float)
    ncells = convert.(Int, ceil.(Δz ./ maxΔz))
    index = repeat_elements(1:length(Δz), ncells)
    return index, ncells
end

"""
Temporary structure used to create SoilColumns.
"""
struct VerticalDomain
    z::Vector{Float} # cell midpoint
    Δz::Vector{Float}  # cell thickness
    geology::Vector{Int}  # cell geology id
    lithology::Vector{Int}  # cell lithology id
    use::UnitRange{Int}  # which part of the column to include
    index::Vector{Int}  # which indices of included column (may contain repeats)
    n::Int  # Total number of cells
end

"""
Temporary structure used to create SoilColumns.
"""
function VerticalDomain(domainbase, modelbase, thickness, Δzmax, geology, lithology)
    ztop = modelbase .+ cumsum(thickness)
    zbot = ztop .- thickness
    zmid = ztop .- 0.5 .* thickness

    base_index = findfirst(zmid .> domainbase)
    top_index = findlast(.!ismissing.(geology))
    use = base_index:top_index

    Δz = @view thickness[use]
    lithology = @view lithology[use]
    geology = @view geology[use]

    index, ncells = discretize(Δz, Δzmax)
    Δz = (Δz./ncells)[index]
    z = zbot[base_index] .+ cumsum(Δz) .- 0.5 .* Δz
    n = sum(ncells)
    return VerticalDomain(z, Δz, geology[index], lithology[index], use, index, n)
end

"""
Initialize a model with specified groundwater, consolidation, and oxidation
processes from a netCDF file and CSV lookup table describing the subsoil
parameters, appropriate for the chosen processes.
"""
function Model(
    groundwater::Type,
    consolidation::Type,
    oxidation::Type,
    preconsolidation::Type,
    timestepper,
    path_subsoil,
    path_lookup,
    Δzmax,
)
    reader = prepare_subsoil_reader(path_subsoil, path_lookup)

    x = ncread(reader, :x)
    y = ncread(reader, :y)
    base = ncread(reader, :base)
    domainbase = ncread(reader, :domainbase)

    columntype = SoilColumn{groundwater,consolidation,preconsolidation,oxidation}
    columns = Vector{columntype}()
    index = Vector{CartesianIndex}()

    for I in CartesianIndices(domainbase)
        ismissing(domainbase[I]) && continue

        geology = ncread3d(reader, :geology, I)
        lithology = ncread3d(reader, :lithology, I)
        thickness = ncread3d(reader, :thickness, I)

        domain =
            VerticalDomain(domainbase[I], base[I], thickness, Δzmax, geology, lithology)
        g_column = initialize(groundwater, domain, reader, I)
        c_column = initialize(consolidation, preconsolidation, domain, reader, I)
        o_column = initialize(oxidation, domain, reader, I)

        column = SoilColumn(
            domainbase[I],
            x[I[1]],
            y[I[2]],
            domain.z,
            domain.Δz,
            g_column,
            c_column,
            o_column,
        )
        # Set values such as preconsolidation stress, τ0, etc.        
        # This requires groundwater: pore pressure, etc.
        apply_preconsolidation!(column)
        push!(columns, column)
        push!(index, I)
    end

    shape = size(domainbase)
    output =
        Output(x, y, fill(NaN, shape), fill(NaN, shape), fill(NaN, shape), fill(NaN, shape))

    return Model(columns, index, timestepper, output)
end

"""
    simulation(model, path_output, stop_time)
    
Setup a simulation from an initialized model.
"""
function Simulation(model::Model, path_output::String, stop_time::DateTime, forcing)
    writer = prepare_writer(path_output, model.output.x, model.output.y)
    clock = Clock(DateTime[], 1, stop_time)
    simulation = Simulation(model, clock, writer, forcing)
    set_periods!(simulation)
    return simulation
end

"""
Advance a single stress period for all columns.

Timesteps are determined by the total duration and the chosen timestepper.
"""
function advance_forcingperiod!(
    model,
    duration;
    deep_subsidence = nothing,
    stage_indexation = nothing,
    stage_change = nothing,
    aquifer_head = nothing,
)
    timesteps = create_timesteps(model.timestepper, duration)
    for (I, column) in zip(model.index, model.columns)
        # Compute pre-loading stresses, set t to 0, etc.
        prepare_forcingperiod!(column)
        # Apply changes
        for forcing in (stage_indexation, deep_subsidence, stage_change, aquifer_head)
            isnothing(forcing) && continue
            apply_forcing!(forcing, column, I)
        end
        # Compute
        subsidence, consolidation, oxidation = advance_forcingperiod!(column, timesteps)
        # Store result
        model.output.phreatic_level[I] = phreatic_level(column.groundwater)
        model.output.subsidence[I] = subsidence
        model.output.consolidation[I] = consolidation
        model.output.oxidation[I] = oxidation
    end
    return
end

function load_forcing!(forcing, key, time, model)
    !haskey(forcing, key) && return nothing
    input = getfield(forcing, key)
    active = read_forcing!(input, time)
    if active
        prepare_forcingperiod!(input, model)
        return input
    else
        return nothing
    end
end


"""
Advance the simulation by a single forcing period.
Reads new input, computes, writes output.
"""
function advance_forcingperiod!(simulation)
    clock = simulation.clock
    time = currenttime(clock)
    duration = periodduration(clock)
    forcing = simulation.forcing
    model = simulation.model

    advance_forcingperiod!(
        model,
        duration;
        deep_subsidence = load_forcing!(forcing, :deep_subsidence, time, model),
        stage_indexation = load_forcing!(forcing, :stage_indexation, time, model),
        stage_change = load_forcing!(forcing, :stage_change, time, model),
        aquifer_head = load_forcing!(forcing, :aquifer_head, time, model),
    )

    write(simulation.writer, clock, simulation.model.output)
    advance!(clock)
    return
end

"""
Collect the period boundaries from the forcing input.
"""
function set_periods!(simulation)
    clock = simulation.clock
    stop_time = clock.stop_time

    alltimes = []
    for forcing in simulation.forcing
        times = forcing.reader.times
        append!(alltimes, times[times.<stop_time])
    end
    push!(alltimes, stop_time)

    clock.times = unique(alltimes)
    return
end

"""
Run all forcing periods of the simulation.
"""
function run!(simulation)
    clock = simulation.clock
    while currenttime(clock) < clock.stop_time
        advance_forcingperiod!(simulation)
    end
    return
end
