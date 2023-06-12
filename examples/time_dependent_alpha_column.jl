using CSV
using Dates
using DataFrames
using NCDatasets
using Atlans


function lookup_table(path)
    df = CSV.read(
        path,
        DataFrame,
        delim="\t",
        stringtype=String,
        types=Atlans.column_type
    )
    newnames = [
        :geology,
        :lithology,
        :gamma_wet,
        :gamma_dry,
        :c_v,
        :drainage_factor,
        :ocr,
        :a,
        :b,
        :c,
        :oxidation_rate,
        :minimal_mass_fraction_organic
    ]
    rename!(df, newnames)
    insertcols!(df, :shrinkage_degree => 1.0)
    insertcols!(df, :mass_fraction_organic => 0.7)
    insertcols!(df, :rho_bulk => 1000.0)
    insertcols!(df, :mass_fraction_lutum => 0.1)
    return Atlans.build_lookup_tables(df)
end


function create_xcoord!(ds, x)
    defVar(
        ds,
        "x",
        x,
        ("x",),
        attrib=["standard_name" => "projection_x_coordinate", "axis" => "X"],
    )
end


function create_ycoord!(ds, y)
    defVar(
        ds,
        "y",
        y,
        ("y",),
        attrib=["standard_name" => "projection_y_coordinate", "axis" => "Y"],
    )
end


function create_var!(ds, name, values)
    defVar(
        ds,
        name,
        values,
        ("layer",),
    )
end


function subsoil_netcdf()
    filename = tempname()
    ds = NCDatasets.Dataset(filename, "c")
    defDim(ds, "x", 1)
    defDim(ds, "y", 1)
    defDim(ds, "layer", 11)
    create_xcoord!(ds, [1.0])
    create_ycoord!(ds, [1.0])

    create_var!(ds, "geology", reverse([1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2]))
    create_var!(ds, "lithology", reverse([2, 2, 2, 1, 1, 1, 1, 1, 1, 5, 5]))
    create_var!(ds, "thickness", reverse([0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]))

    phreatic = defVar(ds, "phreatic_level", Float64, ("x", "y"))
    base = defVar(ds, "zbase", Float64, ("x", "y"))
    domainbase = defVar(ds, "domainbase", Float64, ("x", "y"))
    surface_level = defVar(ds, "surface_level", Float64, ("x", "y"))
    max_oxidation_depth = defVar(ds, "max_oxidation_depth", Float64, ("x", "y"))
    max_shrinkage_depth = defVar(ds, "max_shrinkage_depth", Float64, ("x", "y"))

    phreatic .= -0.5
    base .= -4.3
    domainbase .= -4.3
    surface_level .= 0.0
    max_oxidation_depth .= 1.2
    max_shrinkage_depth .= 1.2
    return filename
end


function get_subsoil_data(path)
    ds = Dataset(subsoil_netcdf(), "r")
    tables = lookup_table(path)
    data = Dict(Symbol(key) => ds[key][:] for key in keys(ds))
    return Atlans.SubsoilData(data, tables)
end


function create_column(
    groundwater::Type,
    consolidation::Type,
    oxidation::Type,
    preconsolidation::Type,
    shrinkage::Type,
    path,
)
    sd = get_subsoil_data(path)
    ad = Atlans.AdaptiveCellsize(0.25, 0.01)

    x = sd.data[:x]
    y = sd.data[:y]
    domainbase = sd.data[:domainbase]

    I = CartesianIndex(1, 1)
    d = Atlans.prepare_domain(
        domainbase[I],
        sd.data[:zbase][I],
        sd.data[:surface_level][I],
        sd.data[:thickness][:, I],
        ad.Δzmax,
        sd.data[:geology][:, I],
        sd.data[:lithology][:, I],
    )

    g_column = Atlans.initialize(groundwater, d, sd, I)
    c_column = Atlans.initialize(consolidation, preconsolidation, d, sd, I)
    o_column = Atlans.initialize(oxidation, d, sd, I)
    s_column = Atlans.initialize(shrinkage, d, sd, I)

    column = Atlans.SoilColumn(
        domainbase[I],
        x[I[1]],
        y[I[2]],
        d.z,
        d.Δz,
        g_column,
        c_column,
        o_column,
        s_column,
    )
    return column, d
end


function duration(dates, idx)
    dt_milliseconds = dates[idx+1] - dates[idx]
    return Dates.value(dt_milliseconds) / (24 * 3600 * 1000)
end


function temperature_range()
    outputname = tempname()

    time = DateTime.(DateTime(2020):Year(1):DateTime(2040))

    df = DataFrame(
        times=time,
        temperature=LinRange(15, 16, length(time)),
    )
    CSV.write(outputname, df)
    return outputname
end


function calculate_surface_in_time(path, dates, timestepper, temperature=nothing)
    column = create_column(
        Atlans.HydrostaticGroundwater,
        Atlans.NullConsolidation,
        Atlans.CarbonStore,
        Atlans.OverConsolidationRatio,
        Atlans.NullShrinkage,
        path
    )

    surfaces = Vector{Float64}()
    push!(surfaces, Atlans.surface_level(column))
    for ii in 1:length(dates)-1
        t = dates[ii]
        dt = duration(dates, ii)

        timesteps = Atlans.create_timesteps(timestepper, dt)

        Atlans.prepare_forcingperiod!(column, 0.01, 0.0, 0.0)

        if !isnothing(temperature)
            Atlans.read_forcing!(temperature, t)
            Atlans.apply_forcing!(temperature, column, 0.0)
        end

        sub, cons, ox, shr = Atlans.advance_forcingperiod!(column, timesteps)
        push!(surfaces, Atlans.surface_level(column))

        ii += 1
    end
    return surfaces
end


workdir = raw"c:\Users\knaake\OneDrive - Stichting Deltares\Documents\bodemdaling\atlans_test"
table_path = joinpath(workdir, "parameters_ophoging.tsv")

temperature = Atlans.Temperature(temperature_range())

timestepper = Atlans.ExponentialTimeStepper(1.0, 2)
dates = DateTime.(DateTime(2020):Year(1):DateTime(2040))

s_with_forcing = calculate_surface_in_time(
    table_path,
    dates,
    timestepper,
    temperature
)

s_without_forcing = calculate_surface_in_time(
    table_path,
    dates,
    timestepper
)


df = DataFrame(
    time=dates,
    s_with_forcing=s_with_forcing,
    s_without_forcing=s_without_forcing
)

CSV.write(joinpath(workdir, "result_surface.csv"), df)
