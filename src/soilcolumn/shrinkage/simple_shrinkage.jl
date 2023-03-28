"""
    SimpleShrinkage(Δz, n, τ, r, shrinkage)

Simple voxel with attributes to compute shrinkage for.

#Arguments
- `Δz::Float`: Thickness of the voxel. [m]
- `n::Float`: Shrinkage factor of the voxel. [-]
- `τ::Float`: Time dependent factor for shrinkage process. [years]
- `r::Float`: Direction of shrinkage, r is 3 indicates isoptropic. [-]
- `shrinkage::Float`: Computed shrinkage or elevation change over time. [m]
"""
struct SimpleShrinkage <: ShrinkageProcess
    Δz::Float
    n::Float
    τ::Float
    r::Float
    shrinkage::Float
end


"""
    shrink(voxel, Δt)

Shrink a voxel for given time interval.

#Arguments
- `voxel::SimpleShrinkage`: Voxel to shrink.
- `Δt::Float`: Time interval. [days]
"""
function shrink(voxel::SimpleShrinkage, Δt::Float64)
    n_residual = 0.7
    n_next = voxel.n + (voxel.n - n_residual) * (exp(-Δt / (voxel.τ * 365.25)) - 1)

    Δn = n_next - voxel.n

    factor = 0.33 #TODO: get correct factor value from WENR?
    relative_change = (1 + (factor * Δn))^(1.0 / voxel.r) - 1

    shrinkage = voxel.Δz * relative_change
    Δz = voxel.Δz + shrinkage
    Δz = min(voxel.Δz, Δz)

    return SimpleShrinkage(Δz, n_next, voxel.τ, voxel.r, shrinkage)
end
