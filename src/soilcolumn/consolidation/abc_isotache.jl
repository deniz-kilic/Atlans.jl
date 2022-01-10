abstract type AbcIsotache <: ConsolidationProcess end

struct DrainingAbcIsotache <: AbcIsotache
    Δz::Float
    Δz_0::Float
    t::Float
    σ′::Float  # effective stress
    γ_wet::Float  # wet specific mass
    γ_dry::Float  # dry specific mass
    # Degree of consolidation
    c_d::Int  # drainage coefficient
    c_v::Float  # drainage coefficient
    U::Float
    # Isotache parameters
    a::Float
    b::Float
    c::Float
    τ::Float
    consolidation::Float  # Computed consolidation
end

# compute intrinsic time (τ)
function τ_intrinsic(abc::ABC where {ABC<:AbcIsotache}, ocr::Float)
    if abc.c < 1.0e-4
        return 1.0e-9
    else
        return τ_ref * ocr^((abc.b - abc.a) / abc.c)
    end
end

function τ_intermediate(abc::ABC where {ABC<:AbcIsotache}, loadstep::Float64)
    σ_term = (abc.σ′ - loadstep) / abc.σ′
    abc_term = (abc.b - abc.a) / abc.c
    return abc.τ * σ_term^abc_term
end

"""
Water-filled porespace is reduced as the water is pushed out and the soil
becomes denser.
"""
function compress_γ_wet(abc::ABC where {ABC<:AbcIsotache})
    return (abc.γ_wet * abc.Δz - abc.consolidation * γ_water) / (abc.Δz - abc.consolidation)
end

function compress_γ_dry(abc::ABC where {ABC<:AbcIsotache})
    return (abc.γ_dry * abc.Δz) / (abc.Δz - abc.consolidation)
end

function consolidate(abc::DrainingAbcIsotache, σ′::Float, Δt::Float)
    t = abc.t + Δt
    # Degree of consolidation changes
    Unew = U(abc, t)
    ΔU = Unew - abc.U #??abc.U??
    # Effective stress changes
    Δσ′ = σ′ - abc.σ′
    σ′ = abc.σ′ + Unew * Δσ′
    loadstep = ΔU * Δσ′
    # τ changes
    τ_intm = τ_intermediate(abc, loadstep)
    τ = τ_intm + Δt
    # consolidation
    strain = abc.c * log(abc.τ / τ_intm) + log(σ′ / (σ′ - loadstep))
    consolidation = min(abc.Δz, strain * abc.Δz)
    γ_wet = compress_γ_wet(abc)
    γ_dry = compress_γ_dry(abc)
    # return new state
    return DrainingAbcIsotache(
        abc.Δz - consolidation,  # new
        abc.Δz_0,
        t,  # new
        σ′, # new
        γ_wet,  # new
        γ_dry,  # new
        abc.c_d,
        abc.c_v,
        Unew,  # new
        abc.a,
        abc.b,
        abc.c,
        τ,  # new
        consolidation,  # new
    )
end


"""
Turn a collection of vectors into a collection of DrainingAbcIsotache cells.
"""
function draining_abc_isotache_column(
    Δz,
    γ_wet,
    γ_dry,
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
            Δz_0[i],
            0.0,  # t
            0.0,  # σ′
            γ_wet[i],
            γ_dry[i],
            c_d[i],
            c_v[i],
            0.0,  # U
            a[i],
            b[i],
            c[i],
            0.0,  # τ
            0.0,  # consolidation
        )
        consolidation[i] = cell
    end
    return consolidation
end
