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
forcing = (
    deep_subsidence = Atlans.DeepSubsidence("deep-subsidence.nc"),
    stage_change = Atlans.StageChange("stage-change.nc"),
)
simulation = Atlans.Simulation(model, "output.nc", DateTime("2020-01-01"), forcing)
run!(simulation)
