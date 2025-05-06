abstract type AbstractTimeStepper end


"""
	ExponentialTimestepper(start::Float, multiplier::T)

Struct to discretize time steps (in days) within each stress period.
"""
struct ExponentialTimeStepper{T} <: AbstractTimeStepper
    start::Float
    multiplier::T
end


"""
Based on the duration and the timestepper, create the required timesteps.
"""
function create_timesteps(timestepper::ExponentialTimeStepper, duration)
    Δt = timestepper.start
    Δt > duration && error("base timestep exceeds duration")
    steps = Vector{Float}()
    remainder = duration

    steps = Vector{Float}()
    while remainder > Δt
        push!(steps, Δt)
        remainder -= Δt
        Δt = Δt * timestepper.multiplier
    end
    push!(steps, remainder)

    return steps
end


"""
	Clock(time::Vector{DateTime}, iteration::int, stop_time::DateTime)

Object to keep track of the stress periods, number of iterations and stop time of
an Atlantis Simulation.
"""
mutable struct Clock
    times::Vector{DateTime}
    iteration::Int
    stop_time::DateTime
end


"""
Advances the clock by one iteration.
"""
function advance!(clock)
    clock.iteration += 1
    return
end


"""
Return current time.
"""
function currenttime(clock)
    return clock.times[clock.iteration]
end


"""
Compute duration of forcing period from current time.
"""
function periodduration(clock)
    dt_milliseconds = clock.times[clock.iteration + 1] - clock.times[clock.iteration]
    return Dates.value(dt_milliseconds) / (24 * 3600 * 1000)
end
