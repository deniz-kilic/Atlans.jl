"""
    SimpleShrinkage(Δz, n, τ, r, shrinkage)

Voxel with attributes to compute shrinkage for.

#Arguments
- `Δz::Float`: Thickness of the Voxel [m].
- `n::Float`: Shrinkage factor of the voxel [-].
- `τ::Float`: Time constant for shrinkage [-].
- `r::Float`: Direction of shrinkage, r is 3 indicates isoptropic [-].
- `shrinkage::Float`: Computed shrinkage or elevation change over time [m].
"""
struct SimpleShrinkage <: ShrinkageProcess
    Δz::Float
    n::Float
    τ::Float
    r::Float
    shrinkage::Float
end


function SimpleShrinkage(Δz::Float, n::Float)
    τ_years = 60.0
    seconds_per_year = 31556926.0
    τ_days = τ_years * seconds_per_year
    r = 3.0
    return SimpleShrinkage(Δz, n, τ_days, r, NaN)
end


"""
    shrink(sc, Δt)

Shrink a voxel for given time interval.

#Arguments
- `sc::SimpleShrinkage`: Voxel to shrink.
- `Δt::Float`: Time interval [days].
"""
function shrink(sc::SimpleShrinkage, Δt::Float64)
    n_residual = 0.7
    n_next = sc.n + (sc.n - n_residual) * (exp(-Δt / sc.τ) - 1)

    Δn = n_next - sc.n

    factor = 0.33 #TODO: get correct factor value from WENR?
    relative_change = (1 + (factor * Δn))^(1.0 / sc.r) - 1

    shrinkage = sc.Δz * relative_change
    Δz = sc.Δz + shrinkage
    Δz = min(sc.Δz, Δz)

    return SimpleShrinkage(Δz, n_next, sc.τ, sc.r, shrinkage)
end
