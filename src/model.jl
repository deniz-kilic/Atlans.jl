struct Reader
    dataset::NCDataset
    params::Dict{Symbol,String}
    times::Vector{DateTime}
end

struct Writer
    dataset::NCDataset
    params::Dict{Symbol,String}
end

mutable struct Clock{T}
    time::T
    iteration::Int
    Δt::Day
end

struct Model
    units::Vector{D} where {D<:DrainageUnit}
    clock::Clock
    reader::Reader
    writer::Writer
end

function advance!(clock, time)
    clock.iteration += 1
    clock.time = time
    return clock
end
