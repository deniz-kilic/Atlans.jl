# Reading
# -------
function ncread(R::Reader, param::Symbol)
    varname = R.params[param]
    return R.dataset[varname][:]
end

function ncread(R::Reader, param::Symbol, time_index::Int)
    varname = R.params[param]
    return R.dataset[varname][:, time_index]
end

function prepare_reader(path)
    ds = Dataset(path, "r")
    times = ds["time"][:]
    reader = Reader(
        ds,
        Dict(
            :inflow => "inflow",
            :fraction => "fraction",
            :priority => "priority",
            :flow => "flow",
        ),
        times,
    )
    return reader
end

# Writing
# -------
"Add a new time to the unlimited time dimension, and return the index"
function add_time(ds, time)
    i = length(ds["time"]) + 1
    ds["time"][i] = time
    return i
end

function setup_output_netcdf(path, time_units, calendar)::NCDataset
    ds = NCDatasets.NCDataset(path, "c")
    defDim(ds, "time", Inf)  # unlimited dimension
    defVar(
        ds,
        "time",
        Float,
        ("time",),
        attrib = ["units" => time_units, "calendar" => calendar],
    )
    defVar(ds, "shortage", Float, ("node", "time"))
    defVar(ds, "flow", Float, ("edge", "time"))
    return ds
end

function prepare_writer(path, time_units, calendar)::Writer
    ds = setup_output_netcdf(path, time_units, calendar)
    return Writer(ds, topology, Dict(:shortage => "shortage", :flow => "flow"))
end
