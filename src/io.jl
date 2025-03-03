# Reading
# -------
function prepare_subsoil_data(path_nc, path_csv)
    ds = Dataset(path_nc, "r")
    df = read_params_table(path_csv)
    tables = build_lookup_tables(df)
    data = Dict(Symbol(key) => Array(ds[key]) for key in keys(ds))
    return SubsoilData(data, tables)
end


function prepare_reader(path)
    ds = Dataset(path, "r")
    times = ds["time"][:]
    return Reader(ds, Dict(Symbol(key) => key for key in keys(ds)), times)
end


function ncread(reader, param::Symbol)
    varname = reader.params[param]
    return Array(reader.dataset[varname])
end

function ncread(reader, param::Symbol, t::TimeType)
    varname = reader.params[param]
    times = reader.times
    i = searchsortedfirst(times, t)
    i === nothing && throw(DomainError("time $t not in dataset times $(times)"))
    return reader.dataset[varname][:, :, i]
end


function ncread2d(reader, param::Symbol, index)
    varname = reader.params[param]
    return reader.dataset[varname][index]
end


function ncread3d(reader, param::Symbol, index)
    varname = reader.params[param]
    return reader.dataset[varname][:, index]
end


function ncread4d(reader, param::Symbol, t::TimeType)
    varname = reader.params[param]
    times = reader.times
    i = searchsortedfirst(times, t)
    i === nothing && throw(DomainError("time $t not in dataset times $(times)"))
    return reader.dataset[varname][:, :, :, i]
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
    df = CSV.read(path, DataFrame; delim=",", stringtype=String, types=column_type)
    return df
end


function build_lookup_tables(df)
    geology = df[:, :geology]
    lithology = df[:, :lithology]
    tables = Dict{Symbol,Dict{Tuple{Int,Int},Float}}()
    for name in names(df, Float)
        tables[Symbol(name)] = Dict(
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


function fetch_field(subsoil, field, I, domain)
    if haskey(subsoil.data, field)
        values = subsoil.data[field][domain.index, I]
    elseif haskey(subsoil.tables, field)
        values = lookup(subsoil.tables[field], domain.geology, domain.lithology)
    else
        error("Mandatory field not in netCDF or table: $field")
    end
    return values
end


function fetch_field(table, field, domain)
    if haskey(table, field)
        values = lookup(table[field], domain.geology, domain.lithology)
    else
        error("Mandatory field not in Surcharge lookup table: $field")
    end
    return values
end


fetch_field(reader, ::Type{PreOverburdenPressure}, I, domain) =
    fetch_field(reader, :pop, I, domain)

fetch_field(table, ::Type{PreOverburdenPressure}, domain) =
    fetch_field(table, :pop, domain)

fetch_field(reader, ::Type{OverConsolidationRatio}, I, domain) =
    fetch_field(reader, :ocr, I, domain)

fetch_field(table, ::Type{OverConsolidationRatio}, domain) =
    fetch_field(table, :ocr, domain)


function xy_size(reader)
    size_x = length(reader.dataset["x"])
    size_y = length(reader.dataset["y"])
    return (size_x, size_y)
end


function xyz_size(reader)
    size_x = length(reader.dataset["x"])
    size_y = length(reader.dataset["y"])
    size_z = length(reader.dataset["layer"])
    return (size_x, size_y, size_z)
end

# Writing
# -------
function setup_output_netcdf(path, x, y)::NCDataset
    time_units = "days since 1900-01-01 00:00:00"
    calendar = "proleptic_gregorian"

    ds = NCDatasets.NCDataset(path, "c")
    defDim(ds, "time", Inf)  # unlimited dimension
    defVar(
        ds,
        "time",
        Float,
        ("time",),
        attrib=["units" => time_units, "calendar" => calendar,
        "axis" => "T", "long_name" => "time", "standard_name" => "time"],
    )
    defVar(
        ds,
        "x",
        x,
        ("x",),
        attrib=["long_name" => "x-coordinate in cartesian system",
        "standard_name" => "projection_x_coordinate", 
        "axis" => "X", "units" => "Meter",
        "cell_size" => 100.0],
            
    )
    defVar(
        ds,
        "y",
        y,
        ("y",),
        attrib=[ "long_name" => "y-coordinate in cartesian system", 
        "standard_name" => "projection_y_coordinate", "axis" => "Y", 
        "grid_mapping" => "oblique_stereographic", "units" => "Meter",
        "cell_size" => 100.0],
    )

    ds.attrib["Conventions"] = "CF-1.6"
    ds.attrib["title"] = "Land subsidence output"
    ds.attrib["pixel_size_x"] = 100.0
    ds.attrib["pixel_size_y"] = 100.0                    
                    
    defVar(ds, "oblique_stereographic",
            Int32,[],
            attrib=["grid_mapping_name" => "oblique_stereographic",
                    "longitude_of_central_meridian" => 5.38763888888889,
                    "false_easting" => 155000.,
                    "false_northing" => 463000.,
                    "latitude_of_projection_origin" => 52.1561605555556,
                    "scale_factor_at_central_meridian" => 0.9999079,
                    "scale_factor_at_projection_origin" => 0.9999079,
                    "long_name" => "CRS definition",
                    "longitude_prime_meridian" => 0.0,
                    "semi_major_axis" => 6378137.155,
                    "inverse_flattening" => 299.1528128,                                        
                    "spatial_ref" => "PROJCS[\"Amersfoort / RD New\",GEOGCS[\"Amersfoort\",DATUM[\"Amersfoort\",SPHEROID[\"Bessel 1841\",6377397.155,299.1528128,AUTHORITY[\"EPSG\",\"7004\"]],AUTHORITY[\"EPSG\",\"6289\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4289\"]],PROJECTION[\"Oblique_Stereographic\"],PARAMETER[\"latitude_of_origin\",52.15616055555555],PARAMETER[\"central_meridian\",5.38763888888889],PARAMETER[\"scale_factor\",0.9999079],PARAMETER[\"false_easting\",155000],PARAMETER[\"false_northing\",463000],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"X\",EAST],AXIS[\"Y\",NORTH],AUTHORITY[\"EPSG\",\"28992\"]]",                    
                    "EPSG" => 28992]
                    )


    defVar(ds, "phreatic_level", Float, ("x", "y", "time"), 
            attrib=["units" => "m",
            "grid_mapping" => "oblique_stereographic"])
    defVar(ds, "consolidation", Float, ("x", "y", "time"),
            attrib=["units" => "m",
            "coordinates" => "x y",
            "grid_mapping" => "oblique_stereographic"])
    defVar(ds, "oxidation", Float, ("x", "y", "time"), 
            attrib=["units" => "m",
            "coordinates" => "x y",
            "grid_mapping" => "oblique_stereographic"])
    defVar(ds, "shrinkage", Float, ("x", "y", "time"), 
            attrib=["units" => "m",
            "coordinates" => "x y",
            "grid_mapping" => "oblique_stereographic"])
    defVar(ds, "subsidence", Float, ("x", "y", "time"), 
            attrib=["units" => "m",
            "coordinates" => "x y",
            "positive" => "down",
            "grid_mapping" => "oblique_stereographic"])
    return ds
end


"Add a new time to the unlimited time dimension, and return the index"
function add_time(ds, time)
    i = length(ds["time"]) + 1
    ds["time"][i] = time
    return i
end


function prepare_writer(path, x, y)
    ds = setup_output_netcdf(path, x, y)
    return Writer(
        ds,
        Dict(
            :phreatic_level => "phreatic_level",
            :consolidation => "consolidation",
            :oxidation => "oxidation",
            :subsidence => "subsidence",
            :shrinkage => "shrinkage",
        ),
    )
end


function ncwrite(writer::Writer, param, values, time_index)
    varname = writer.params[param]
    writer.dataset[varname][:, :, time_index] .= values
    return
end


function write(writer, clock, output)
    add_time(writer.dataset, currenttime(clock))
    ncwrite(writer, :phreatic_level, output.phreatic_level, clock.iteration)
    ncwrite(writer, :subsidence, output.subsidence, clock.iteration)
    ncwrite(writer, :consolidation, output.consolidation, clock.iteration)
    ncwrite(writer, :oxidation, output.oxidation, clock.iteration)
    ncwrite(writer, :shrinkage, output.shrinkage, clock.iteration)
    return
end
