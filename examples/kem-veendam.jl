using Atlans
using Dates


path_nc = "examples/subsoil-model.nc"
path_csv = "examples/parameters.csv"
forcing = (
    deep_subsidence = Atlans.DeepSubsidence("examples/subsidence.nc"),
    stage_change = Atlans.StageChange("examples/change.nc"),
)

model = Atlans.Model(
    Atlans.HydrostaticGroundwater,
    Atlans.DrainingAbcIsotache,
    Atlans.CarbonStore,
    Atlans.OverConsolidationRatio,
    Atlans.AdaptiveCellsize(0.25, 0.01),
    Atlans.ExponentialTimeStepper(1.0, 2),
    path_nc,
    path_csv,
);

additional_times = map(DateTime, ["2020-01-01", "2025-01-01", "2030-01-01", "2035-01-01"])
simulation = Atlans.Simulation(
    model,
    "examples/output.nc",
    DateTime("2040-01-01"),
    forcing,
    additional_times,
);
Atlans.run!(simulation);
