using Atlans
using Dates
using Profile

workdir = raw"c:\Users\knaake\OneDrive - Stichting Deltares\Documents\projects\bugfix_atlantis"

## input for AdaptiveCellsize
max_voxel_thickness = 0.25 # m
split_tolerance = 0.01 # m

## input for Timestepper
start_day = 1.0
multiplier = 2

@show "Build model"
model = Atlans.Model(
	Atlans.HydrostaticGroundwater,
	Atlans.DrainingAbcIsotache,
	Atlans.CarbonStore,
	Atlans.OverConsolidationRatio,
	Atlans.SimpleShrinkage,
	Atlans.AdaptiveCellsize(max_voxel_thickness, split_tolerance),
	Atlans.ExponentialTimeStepper(start_day, multiplier),
	joinpath(workdir, "data", "ychunk.nc"),
	joinpath(workdir, "data", "parameters.csv"),
)

@show "Build forcings"
forcings = Atlans.Forcings(
	temperature = Atlans.Temperature(
		joinpath(workdir, "data", "temperature_forcing_REF2017VP.csv"),
	),
)

@show "Build simulation"
simulation = Atlans.Simulation(
	model,
	joinpath(workdir, "results", "effectmodule_test_chunk.nc"),
	DateTime("2028-01-01"),
	forcings = forcings,
)
Atlans.run!(simulation)
