struct ConstantRate <: OxidationProcess
    Δz::Float64  # cell thickness
    α::Float64  # oxidation rate
end

function oxidate(cr::ConstantRate, Δt::Float64)::Tuple{Float64,ConstantRate}
    change = cr.Δz * (1.0 - exp(-cr.α * Δt))
    oxidation = min(cr.Δz, change)
    Δz = cr.Δz - oxidation
    return oxidation, 
    ConstantRate(
        Δz, #new
        cr.α
    )
end
