
struct SimpleShrinkage <: ShrinkageProcess
    Δz::Float  # cell thickness
    n::Float
    n_residual::Float
    shrinkage::Float
end
 
function shrink(ss:SimpleShrinkage, Δt::Float64)
    Δn = ss.n_residual - ss.n
    n = ss.n + Δn
    shrinkage = 1.0
    # TODO: volume and γ update?
    return SimpleShrinkage(Δz, n, ss.n_residual, shrinkage)
end
