const τ_ref = 1.0
const γ_water = 9810.0

abstract type AbcIsotache <: ConsolidationProcess end

struct DrainingAbcIsotache <: AbcIsotache
    Δz::Float64
    t::Float64
    σ′::Float64  # effective stress
    γ_w::Float64  # wet specific mass
    γ_d::Float64  # dry specific mass
    # Degree of consolidation
    c_d::Int  # drainage coefficient
    c_v::Float64  # drainage coefficient
    U::Float64
    # Isotache parameters
    a::Float64
    b::Float64
    c::Float64
    τ::Float64
end

# compute intrinsic time (τ)
function τ_intrinsic(abc::ABC where {ABC<:AbcIsotache}, ocr::Float64)
    if abc.c < 1.0e-4
        return 1.0e-9
    else
        return τ_ref * ocr^((abc.b - abc.a) / abc.c)
    end
end

function τ_intermediate(abc::ABC where {ABC<:AbcIsotache}, loadstep::Float64)
    σ_term = (σ′ - loadstep) / σ′
    abc_term = (abc.b - abc.a) / abc.c
    return abc.τ * σ_term^abc_term
end

function U(abc::ABC where {ABC<:AbcIsotache}, t)
    T = abc.c_v * t / abc.Δz
    U = (T^3 / (T^3 + 0.5))^(1 / 6)
    return U
end

function compress_γ_wet(abc::ABC where {ABC<:AbcIsotache}, consolidation)
    return (abc.γ_w * abc.Δz - consolidation * γ_water) / (abc.Δz - consolidation)
end

function compress_γ_dry(abc::ABC where {ABC<:AbcIsotache}, consolidation)
    return (abc.γ_d * abc.Δz) / (abc.Δz - consolidation)
end

function consolidate(
    abc::DrainingAbcIsotache,
    σ′::Float64,
    Δt::Float64,
)::Tuple{Float64,DrainingAbcIsotache}
    t = abc.t + Δt
    # Degree of consolidation changes
    U = U(abc, t)
    ΔU = U - abc.U #??abc.U??
    # Effective stress changes
    Δσ′ = σ′ - abc.σ′
    σ′ = abc.σ′ + U * Δσ′
    loadstep = ΔU * Δσ′
    # τ changes
    τ_intm = τ_intermediate(abc, loadstep)
    τ = τ_intm + Δt
    # consolidation
    strain = abc.c * log(abc.τ / τ_intm) + log(σ′ / (σ′ - loadstep))
    consolidation = min(Δz, strain * abc.Δz)
    γ_w = compress_γ_wet(abc, consolidation)
    γ_d = compress_γ_dry(abc, consolidation)
    # return new state
    return consolidation,
    DrainingAbcIsotache(
        abc.Δz - consolidation,  # new
        t,  # new
        σ′, # new
        γ_w,  # new
        γ_d,  # new
        abc.c_d,
        abc.c_v,
        U,  # new
        abc.a,
        abc.b,
        abc.c,
        τ,  # new
    )
end


"""
Turn a collection of vectors into a collection of DrainingAbcIsotache cells.
"""
function draining_abc_isotache_column(
    Δz,
    γ_w,
    γ_d,
    c_d,
    c_v,
    a,
    b,
    c,
)::Vector{DrainingAbcIsotache}
    nlayer = length(Δz)
    consolidation = Vector{DrainingAbcIsotache}(undef, nlayer)
    for i = 1:nlayer
        cell = DrainingAbcIsotache(
            Δz[i],
            0.0,  # t
            0.0,  # σ′
            γ_w[i],
            γ_d[i],
            c_d[i],
            c_v[i],
            0.0,  # U
            a[i],
            b[i],
            c[i],
            0.0,  # τ
        )
        consolidation[i] = cell
    end
    return consolidation
end
