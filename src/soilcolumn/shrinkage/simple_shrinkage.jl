
struct SimpleShrinkage <: ShrinkageProcess
    Δz::Float  # cell thickness
    n::Float # shrinkage degree
    τ::Float # time constant for shrinkage related to soiltype
    r::Float # direction of shrinkage, 3 = isotrope
    shrinkage::Float # computed shrinkage
end


function shrink(ss::SimpleShrinkage, Δt::Float64)
    n_residual = 0.7
    n_next = ss.n + (ss.n - n_residual) * (exp(-Δt / ss.τ) - 1)

    Δn = n_next - ss.n

    factor = 0.001 #TODO: how to determine factor? From subsurface model?
    shr = (factor * Δn)^(1 / ss.r)

    shrinkage = min(ss.Δz, shr)

    return SimpleShrinkage(Δz, n_next, ss.τ, ss.r, shrinkage)
end
