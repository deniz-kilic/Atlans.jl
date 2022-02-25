# Reading
# -------

function ncread(R::Reader, param::Symbol)
    varname = R.params[param]
    return R.dataset[varname][:]
end

function ncread(R::Reader, param::Symbol, t::TimeType)
    varname = R.params[param]
    times = R.times
    # this behaves like a forward fill interpolation
    i = searchsortedfirst(==(t), times)
    i === nothing && throw(DomainError("time $t not in dataset times $(times)"))
    return R.dataset[varname][:, i]
end

function ncread2d(R::Reader, param::Symbol, index::CartesianIndex)
    varname = R.params[param]
    return R.dataset[varname][index]
end

function ncread3d(R::Reader, param::Symbol, index::CartesianIndex)
    varname = R.params[param]
    return R.dataset[varname][:, index]
end

function column_type(_, name)
    if name == :geology || name == :lithology
        return Int
    elseif name == :geology_name || name == :lithology_name
        return String
    else
        return Float
    end
end

function read_params_table(path)
    df = CSV.read(path, DataFrame; delim = ",", stringtype = String, types = column_type)
    return df
end

function build_lookup_tables(df)
    geology = df[:, :geology]
    lithology = df[:, :lithology]
    tables = Dict{:symbol,Dict{Tuple{Int,Int},Float}}
    for name in names(df, Float)
        tables[name] = Dict(
            (geo, lit) => val for (geo, lit, val) in zip(geology, lithology, df[:, name])
        )
    end
    return tables
end

function lookup(table, geology, lithology)
    n = length(geology)
    found = Vector{Float}(undef, n)
    for i = 1:n
        found[i] = table[(geology[i], lithology[i])]
    end
    return found
end

function prepare_subsoil_reader(path_nc, path_csv)
    ds = Dataset(path_nc, "r")
    df = read_params_table(path_csv)
    tables = build_lookup_tables(df)
    return SubsoilReader(ds, tables, Dict(key => Symbol(key) for key in keys(ds)))
end

function prepare_reader(path)
    ds = Dataset(path, "r")
    times = ds["time"][:]
    return Reader(ds, Dict(key => Symbol(key) for key in keys(ds)), times)
end

function fetch_field(reader, field, I, lithology, geology)
    if haskey(reader.params, field)
        values = ncread3d(reader, field, I)
    elseif haskey(reader.tables, field)
        values = lookup(tables[field], geology, lithology)
    else
        error("Mandatory field not in netCDF or table")
    end
end

fetch_field(reader, ::PreOverburdenPressure, I, lithology, geology) =
    fetch_field(reader, :pop, I, lithology, geology)
fetch_field(reader, ::OverConsolidationRatio, I, lithology, geology) =
    fetch_field(reader, :ocr, I, lithology, geology)

function xy_size(reader)
    xsize = size(reader.ds["x"])
    ysize = size(reader.ds["y"])
    return (xsize, ysize)
end

# Writing
# -------
"Add a new time to the unlimited time dimension, and return the index"
function add_time(ds, time)
    i = length(ds["time"]) + 1
    ds["time"][i] = time
    return i
end

function setup_output_netcdf(path)::NCDataset
    ds = NCDatasets.NCDataset(path, "c")
    defDim(ds, "time", Inf)  # unlimited dimension
    defVar(
        ds,
        "time",
        Float,
        ("time",),
        #attrib = ["units" => time_units, "calendar" => calendar],
    )
    defVar(ds, "x", Float, ("x",))
    defVar(ds, "y", Float, ("y",))

    defVar(ds, "phreatic_level", Float, ("x", "y", "time"))
    defVar(ds, "consolidation", Float, ("x", "y", "time"))
    defVar(ds, "oxidation", Float, ("x", "y", "time"))
    defVar(ds, "subsidence", Float, ("x", "y", "time"))
    return ds
end

function prepare_writer(path)::Writer
    ds = setup_output_netcdf(path)
    return Writer(
        ds,
        Dict(
            :phreatic_level => "phreatic_level",
            :consolidation => "consolidation",
            :oxidation => "oxidation",
            :subsidence => "subsidence",
        ),
    )
end

function ncwrite(writer::Writer, param, values, time_index)
    varname = R.params[param]
    writer.dataset[varname][:, :, time_index] .= values
    return
end

function write(writer, clock, output)
    add_time(writer.ds, clock.time)
    ncwrite(writer, :phreatic_level, output.phreatic_level, clock.iteration)
    ncwrite(writer, :subsidence, output.subsidence, clock.iteration)
    ncwrite(writer, :consolidation, output.consolidation, clock.iteration)
    ncwrite(writer, :oxidation, output.oxidation, clock.iteration)
    return
end
