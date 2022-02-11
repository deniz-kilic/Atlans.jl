# # Example of running a simulation
#
# This example demonstrates
#
# * Choosing concepts for groundwater, consolidation, oxidation
# * Initializing the simulation
# * Running a number of stress periods, while
# * Writing to an output dataset.

using Atlans

Δzmax = 0.25

model = atlans.Model(
    atlans.HydrostaticGroundwater,
    atlans.DrainingAbcIsotache,
    atlans.CarbonStore,
    atlans.OverConsolidationRatio,
    atlans.ExponentialTimestepper(1, 2),
    "subsoil-model.nc",
    "parameters.csv",
    Δzmax,
)

simulation = atlans.Simulation(
    model,
    "output.nc",
    DateTime("2020-01-01"),
)

deep_subsidence = atlans.DeepSubsidence("deep-subsidence.nc")
stage_change = atlans.StageChange("stage-change.nc")
set_forcing!(simulation, deep_subsidenc)
set_forcing!(simulation, stage_change)
set_periods!(simulation)

run!(simulation)
