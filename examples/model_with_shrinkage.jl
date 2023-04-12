using Atlans


workdir = raw"c:\Users\knaake\OneDrive - Stichting Deltares\Documents\atlans_test\implementatie_rijping"
path_nc = joinpath(workdir, "subsoil-model-fase2.nc")
path_csv = joinpath(workdir, "parameters.csv")


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

additional_times = map(DateTime, ["2020-01-01", "2025-01-01", "2030-01-01", "2035-01-01"])

simulation = Atlans.Simulation(
    model,
    "examples/output.nc",
    DateTime("2040-01-01"),
    forcing,
    additional_times,
)

Atlans.run!(simulation)
