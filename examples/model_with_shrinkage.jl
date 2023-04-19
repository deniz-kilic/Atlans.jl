using Atlans
using CSV
using NCDatasets
using DataFrames

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


function subsoil_netcdf()
    filename = tempname()
    ds = NCDatasets.Dataset(filename, "c")
    defDim(ds, "x", 2)
    defDim(ds, "y", 3)
    defDim(ds, "layer", 4)
    create_xcoord!(ds, [12.5, 37.5])
    create_ycoord!(ds, [87.5, 62.5, 37.5])

    geology = defVar(ds, "geology", Int, ("layer", "x", "y"))
    lithology = defVar(ds, "lithology", Int, ("layer", "x", "y"))
    thickness = defVar(ds, "thickness", Float64, ("layer", "x", "y"))
    phreatic = defVar(ds, "phreatic_level", Float64, ("x", "y"))
    base = defVar(ds, "zbase", Float64, ("x", "y"))
    domainbase = defVar(ds, "domainbase", Float64, ("x", "y"))
    surface_level = defVar(ds, "surface_level", Float64, ("x", "y"))
    max_oxidation_depth = defVar(ds, "max_oxidation_depth", Float64, ("x", "y"))
    max_shrinkage_depth = defVar(ds, "max_shrinkage_depth", Float64, ("x", "y"))

    geology .= 1
    lithology .= 2
    phreatic .= 0.5
    thickness .= 0.25
    base .= 0.0
    domainbase .= 0.0
    surface_level .= 1.0
    max_oxidation_depth .= 1.2
    max_shrinkage_depth .= 1.2
    return filename
end


function params_table()
    filename = tempname()
    df = DataFrame(
        geology_name=["NAWA"],
        lithology_name=["sand"],
        geology=[1],
        lithology=[2],
        gamma_wet=[15000.0],
        gamma_dry=[10000.0],
        drainage_factor=[2.0],
        c_v=[0.006912],
        a=[0.01737],
        b=[0.1303],
        c=[0.008686],
        ocr=[2.15],
        mass_fraction_organic=[0.2],
        minimal_mass_fraction_organic=[0.05],
        oxidation_rate=[0.001],
        rho_bulk=[1000.0],
        mass_fraction_lutum=[0.7],
        shrinkage_degree=[1.4]
    )
    CSV.write(filename, df)
    return filename
end


path_nc = subsoil_netcdf()
path_csv = params_table()


model = Atlans.Model(
    Atlans.HydrostaticGroundwater,
    Atlans.DrainingAbcIsotache,
    Atlans.CarbonStore,
    Atlans.OverConsolidationRatio,
    Atlans.Shrinkage,
    Atlans.AdaptiveCellsize(0.25, 0.01),
    Atlans.ExponentialTimeStepper(1.0, 2),
    path_nc,
    path_csv,
)

