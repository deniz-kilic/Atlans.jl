using Atlans

Δzmax = 0.25

basedir = raw"N:\Projects\11209000\11209259\B. Measurements and calculations\009 effectmodule bodemdaling\data"

model = Atlans.Model(
    Atlans.HydrostaticGroundwater,
    Atlans.DrainingAbcIsotache,
    Atlans.CarbonStore,
    Atlans.OverConsolidationRatio,
    Atlans.AdaptiveCellsize(0.25, 0.01),
    Atlans.ExponentialTimeStepper(1.0, 2),
    joinpath(basedir, "3-input/subsurface_model.nc"),
    joinpath(basedir, "3-input/param_table_test.csv"),
    Δzmax,
)
forcing = (
    stage_change = Atlans.StageChange(joinpath(basedir, "3-input/stage_change.nc")),
)

additional_times = map(DateTime, ["1920-01-01", "1930-01-01", "1940-01-01", "1950-01-01", "1960-01-01", "1970-01-01", "1980-01-01", "1990-01-01", "2000-01-01"])

simulation = Atlans.Simulation(
    model,
    joinpath(basedir, "4-output/test.nc"),
    DateTime("2010-01-01"),
    forcing,
    additional_times,
);
Atlans.run!(simulation);