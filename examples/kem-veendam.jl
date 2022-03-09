using Revise
using Atlans
using Dates


function run()
    path_nc = "examples/subsoil-model-tiny.nc"
    path_csv = "examples/parameters.csv"
    forcing = (
        deep_subsidence = Atlans.DeepSubsidence("examples/subsidence-tiny.nc"),
        stage_change = Atlans.StageChange("examples/change-tiny.nc"),
    )

    timestepper = Atlans.ExponentialTimeStepper(1.0, 2)

    Δzmax = 0.25
    model = Atlans.Model(
        Atlans.HydrostaticGroundwater,
        Atlans.DrainingAbcIsotache,
        Atlans.CarbonStore,
        Atlans.OverConsolidationRatio,
        timestepper,
        path_nc,
        path_csv,
        Δzmax,
    )

    simulation =
        Atlans.Simulation(model, "examples/output-tiny.nc", DateTime("2018-01-01"), forcing)
    Atlans.run!(simulation)
    return
end

sim = run();
