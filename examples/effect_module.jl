using Atlans
using Dates

Δzmax = 0.25

basedir = raw"c:/Users/knaake/OneDrive - Stichting Deltares/Documents/bodemdaling/for_deniz"

path_nc = joinpath(basedir, "base_subsurface_model.nc")
path_table = joinpath(basedir, "parameters.csv")

model = Atlans.Model(
    Atlans.HydrostaticGroundwater,
    Atlans.DrainingAbcIsotache,
    Atlans.CarbonStore,
    Atlans.OverConsolidationRatio,
    Atlans.NullShrinkage,
    Atlans.AdaptiveCellsize(Δzmax, 0.01),
    Atlans.ExponentialTimeStepper(1.0, 2),
    joinpath(basedir, "base_subsurface_model.nc"),
    joinpath(basedir, "parameters.csv"),
)
forcing = (
    stage_change=Atlans.StageChange(joinpath(basedir, "3-input/REF2017BP18_stage_change.nc")),
)

# additional_times = map(DateTime, ["1920-01-01", "1930-01-01", "1940-01-01", "1950-01-01", "1960-01-01", "1970-01-01", "1980-01-01", "1990-01-01", "2000-01-01"])

simulation = Atlans.Simulation(
    model,
    joinpath(basedir, "4-output/test.nc"),
    DateTime("1940-01-01"),
    forcing,
    additional_times,
);
Atlans.run!(simulation);