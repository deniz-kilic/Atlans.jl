struct StageIndexation <: Forcing
    percentile::Int
    weir_area::Array{OptionalInt}
    change::Array{OptionalFloat}
    reader::Reader
end

struct DeepSubsidence <: Forcing
    subsidence::Array{OptionalFloat}
    reader::Reader
end

struct StageChange <: Forcing
    change::Array{OptionalFloat}
    reader::Reader
end

struct AquiferHead <: Forcing
    head::Array{OptionalFloat}
    reader::Reader
end

#struct Surcharge{G,C,O} <: Forcing
#    surcharge_index::Array{OptionalInt}
#    surcharge::Vector{ColumnSurcharge{G,C,O}}
#    reader::Reader
#end

# Initializing

function StageIndexation(path, percentile)
    reader = prepare_reader(path)
    size = xy_size(reader)
    return StageIndexation(
        percentile,
        Array{OptionalInt}(missing, size),
        Array{OptionalFloat}(missing, size),
        reader,
    )
end

function DeepSubsidence(path)
    reader = prepare_reader(path)
    size = xy_size(reader)
    return DeepSubsidence(Array{OptionalFloat}(missing, size), reader)
end

function StageChange(path)
    reader = prepare_reader(path)
    size = xy_size(reader)
    return StageChange(Array{OptionalFloat}(missing, size), reader)
end

function AquiferHead(path)
    reader = prepare_reader(path)
    size = xy_size(reader)
    return AquiferHead(Array{OptionalFloat}(missing, size), reader)
end

function Surcharge(path, groundwater::Type, consolidation::Type, oxidation::Type)
    error("Not implemented")
end

# Reading

function read_forcing!(si::StageIndexation, time)
    if time in forcing.reader.times
        si.weir_area .= ncread(si.reader, :weir_area, time)
        si.change .= ncread(si.reader, :change, time)
        return true
    end
    return false
end

function read_forcing!(ds::DeepSubsidence, time)
    if time in ds.reader.times
        ds.subsidence .= ncread(ds.reader, :subsidence, time)
        return true
    end
    return false
end

function read_forcing!(sc::StageChange, time)
    if time in sc.reader.times
        sc.change .= ncread(sc.reader, :stage_change, time)
        return true
    end
    return false
end

function read_forcing!(ah::AquiferHead, time)
    if time in ah.reader.times
        ah.weir_area .= ncread(ah.reader, :aquifer_change, time)
        return true
    end
    return false
end

# Prepare forcing period

prepare_forcingperiod!(_, _) = nothing

function prepare_forcingperiod!(si::StageIndexation, model)
    si.change .= 0.0
    weir_areas = si.weir_area
    isarea = Array{Bool}(undef, size(weir_areas))
    for area in unique(weir_areas)
        isarea .= (weir_areas .== area)
        si.change[isarea] = percentile(model.output.subsidence[isarea], si.percentile)
    end
    return
end

# Apply forcing to a column

function apply_forcing!(si::StageIndexation, column, I)
    change = si.change[I]
    (ismissing(change) || change == 0.0) && return

    set_phreatic_difference!(column, change)
    return
end

function apply_forcing!(ds::DeepSubsidence, column, I)
    subsidence = ds.subsidence[I]
    ismissing(subsidence) && return

    set_deep_subsidence!(column, subsidence)
    return
end

function apply_forcing!(si::StageChange, column, I)
    change = si.change[I]
    (ismissing(change) || change == 0.0) && return

    set_phreatic_difference!(column, change)
    return
end

function apply_forcing!(ah::AquiferHead, column, I)
    head = ah.head[I]
    ismissing(head) && return

    set_aquifer(column, head)
    return
end


#function apply_forcing!(sc::Surcharge, column, I)
#    error("Not implemented")
#    return
#end
#
