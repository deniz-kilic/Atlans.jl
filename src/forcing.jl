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
        fill(1.0, size),
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


function Temperature(path)
    temperature = CSV.read(path, DataFrame)
    return Temperature(missing, temperature)
end

# Reading

function read_forcing!(si::StageIndexation, time)
    if time in si.reader.times
        si.weir_area .= ncread(si.reader, :weir_area, time)
        si.factor .= ncread(si.reader, :factor, time)
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


function read_forcing!(t::Temperature, time)
    time_idx = findfirst(t.table.times .== time)
    if !isnothing(time_idx)
        t.temp = t.table.temperature[time_idx]
        return true
    end
    return false
end

# Prepare forcing period

prepare_forcingperiod!(_::Forcing, _::Model) = nothing

function prepare_forcingperiod!(si::StageIndexation, model::Model)
    change_to_negative = -1
    si.change .= 0.0
    weir_areas = si.weir_area
    replace!(weir_areas, missing => typemin(Int64))
    isarea = fill(false, size(weir_areas))

    for area in unique(weir_areas)
        area == typemin(Int64) && continue

        isarea .= (weir_areas .== area)

        if area == -1 # Change is the subsidence per column
            si.change[isarea] .= model.output.subsidence[isarea]
        else
            try
                si.change[isarea] .= nanpercentile(
                    model.output.subsidence[isarea],
                    si.percentile
                ) * change_to_negative
            catch ArgumentError
                continue
            end
        end
    end
    return
end

# Apply forcing to a column


function get_elevation_shift(ds::DeepSubsidence, column, I)
    subsidence = ds.subsidence[I]
    ismissing(subsidence) && return 0.0
    return subsidence
end


function get_elevation_shift(si::StageIndexation, column, I)
    change = si.change[I]
    (ismissing(change) || change == 0.0) && return 0.0
    return change
end


function get_elevation_shift(si::StageChange, column, I)
    change = si.change[I]
    (ismissing(change) || change == 0.0) && return 0.0
    return change
end


function apply_forcing!(si::StageIndexation, column, I)
    change = si.change[I]
    factor = si.factor[I]
    (ismissing(change) || change == 0.0 || ismissing(factor)) && return

    set_phreatic_difference!(column, change * factor)
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


function apply_forcing!(t::Temperature, column, _)
    oc = column.oxidation
    for i in 1:length(oc.cells)
        cell = oc.cells[i]
        newcell = update_alpha(cell, t.temp)
        oc.cells[i] = newcell
    end
end


#function apply_forcing!(sc::Surcharge, column, I)
#    error("Not implemented")
#    return
#end
#


