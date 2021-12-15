struct ConstantRate <: OxidationProcess
    Δz::Float64  # cell thickness
    α::Float64  # oxidation rate
    oxidation::Float  # computed oxidation
end

function oxidate(cr::ConstantRate, Δt::Float64)
    change = cr.Δz * (1.0 - exp(-cr.α * Δt))
    oxidation = min(cr.Δz, change)
    Δz = cr.Δz - oxidation
    return ConstantRate(
        Δz, #new
        cr.α
        oxidation,  #new
    )
end
