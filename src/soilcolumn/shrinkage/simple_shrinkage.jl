
struct SimpleShrinkage <: ShrinkageProcess
    Δz::Float  # cell thickness
    n::Float # shrinkage degree
    shrinkage::Float # computed shrinkage
end

function shrink(ss::SimpleShrinkage, Δt::Float64, r::Int64)
    n_residual = 0.7
    Δn = n_residual - ss.n
    n = ss.n + Δn * (exp(-Δt) - 1)

    shrinkage = 1.0
    Δz =
    # TODO: volume and γ update?
        return SimpleShrinkage(Δz, n, shrinkage)
end
