"""
    SimpleShrinkage(Δz, n, τ, r, shrinkage)

Simple voxel with attributes to compute shrinkage for.

#Arguments
- `Δz::Float`: Thickness of the voxel. [m]
- `n::Float`: Shrinkage factor of the voxel. [-]
- `L::Float`: Mass fraction of lutum.
- `H::Float`: Mass fraction of organic.
- `τ::Float`: Time dependent factor for shrinkage process. [days]
- `r::Float`: Direction of shrinkage, r is 3 indicates isoptropic. [-]
- `sf::Float`: TODO: look-up in document [-]
- `shrinkage::Float`: Computed shrinkage or elevation change over time. [m]
"""
struct SimpleShrinkage <: ShrinkageProcess
    Δz::Float
    n::Float
    L::Float
    H::Float
    τ::Float
    r::Float
    sf::Float
    shrinkage::Float
end


function SimpleShrinkage(Δz, n, L, H)
    sf = shrinkage_factor(n, L, H)
    τ = 60.0 * 365.25 # Time dependent factor for shrinkage process
    r = 3.0 # Direction of shrinkage is assumed isoptropic (r=3)
    return SimpleShrinkage(Δz, n, L, H, τ, r, sf, 0.0)
end


function mass_colloidal(L, H)
    return 1 - L - H
end


function mass_solids(R, L, H)
    C = (R / ρR) + (H / ρH) + (L / ρL)
    return C
end


function shrinkage_factor(n, L, H)
    b = 3.0 # constant
    R = mass_colloidal(L, H)
    C = mass_solids(R, L, H)

    pore_volume = (L + b * H) / ρw
    sf = pore_volume / ((n * pore_volume) + (0.2R / ρw) + C)
    return sf
end


"""
    shrink(voxel, Δt)

Shrink a voxel for given time interval.

#Arguments
- `voxel::SimpleShrinkage`: Voxel to shrink.
- `Δt::Float`: Time interval. [days]
"""
function shrink(voxel::SimpleShrinkage, Δt::Float64)
    n_residual = 0.5

    if voxel.n <= n_residual
        return voxel
    end

    n_next = voxel.n + (voxel.n - n_residual) * (exp(-Δt / (voxel.τ)) - 1)

    Δn = voxel.n - n_next

    relative_change = (1 + (voxel.sf * Δn))^(1.0 / voxel.r) - 1

    shrinkage = voxel.Δz * relative_change
    Δz = min(voxel.Δz, voxel.Δz - shrinkage)

    new_sf = shrinkage_factor(n_next, voxel.L, voxel.H)

    return SimpleShrinkage(
        Δz,
        n_next,
        voxel.L,
        voxel.H,
        voxel.τ,
        voxel.r,
        new_sf,
        shrinkage,
    )
end
