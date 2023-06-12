struct DrainingAbcIsotache <: AbstractAbcIsotache
    Δz::Float
    Δz_0::Float
    t::Float
    σ′::Float  # effective stress
    γ_wet::Float  # wet specific mass
    γ_dry::Float  # dry specific mass
    # Degree of consolidation
    c_d::Float  # drainage factor
    c_v::Float  # drainage coefficient
    U::Float
    # Isotache parameters
    a::Float
    b::Float
    c::Float
    τ::Float
    consolidation::Float  # Computed consolidation
end


function DrainingAbcIsotache(Δz, γ_wet, γ_dry, c_d, c_v, a, b, c)
    return DrainingAbcIsotache(
        Δz,
        Δz,
        0.0,
        NaN,
        γ_wet,
        γ_dry,
        c_d,
        c_v,
        NaN,
        a,
        b,
        c,
        NaN,
        0.0,
    )
end


"""
    consolidate(abc::DrainingAbcIsotache, σ′, Δt)
    
Compute consolidation for a single cell.

The cell contains a state U for Terzaghi's degree of consolidation. This state
is updated every consolidate step. During every Δt, a new U is computed.  The
increase in U can be directly related to the pore pressure, and so the increase
in effective stress is equal to the effective stress prior to loading (abc.σ′)
and the new effective stress (σ′) reached when U becomes 1.0.

The degree of consolidation plays only one role: it distributes the load (final
σ′ - initial σ′) over time and as the column might be submerging, the increase
in σ′ may become lower, providing a negative feedback mechanism.
"""
function consolidate(abc::DrainingAbcIsotache, σ′, Δt)
    if abc.Δz == 0
        new = @set abc.consolidation = 0.0
        return new
    end
    t = abc.t + Δt

    # Degree of consolidation changes
    Unew = U(abc, t)
    ΔU = Unew - abc.U

    # Effective stress changes
    load = σ′ - abc.σ′
    σ′ = abc.σ′ + Unew * load
    loadstep = ΔU * load

    # τ changes
    # This also catches cases where c == 0 for non-creeping soils
    if abc.c < 1.0e-4
        τ⃰ = abc.τ
        τ = τ⃰ + Δt
    else
        τ⃰ = abc.τ * ((σ′ - loadstep) / σ′)^((abc.b - abc.a) / abc.c)
        τ = τ⃰ + Δt
    end

    # consolidation
    elastoplastic = abc.a * log(σ′ / (σ′ - loadstep))
    creep = abc.c * log(τ / τ⃰)
    strain = elastoplastic + creep

    # Thickness should not go below 0
    consolidation = min(abc.Δz, strain * abc.Δz)


    # return new state
    return DrainingAbcIsotache(
        abc.Δz - consolidation,  # new
        abc.Δz_0,
        t,  # new
        abc.σ′,
        abc.γ_wet,
        abc.γ_dry,
        abc.c_d,
        abc.c_v,
        Unew,  # new
        abc.a,
        abc.b,
        abc.c,
        τ,  # new,
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


"""
    initialize(::Type{DrainingAbcIsotache}, domain, subsoil, I)

Initialize a ConsolidationColumn for a domain at location I based subsurface input.
"""
function initialize(
    ::Type{DrainingAbcIsotache},
    preconsolidation::Type,
    domain,
    subsoil,
    I,
)::ConsolidationColumn{DrainingAbcIsotache}
    γ_wet = fetch_field(subsoil, :gamma_wet, I, domain)
    γ_dry = fetch_field(subsoil, :gamma_dry, I, domain)
    c_d = fetch_field(subsoil, :drainage_factor, I, domain)
    c_v = fetch_field(subsoil, :c_v, I, domain)
    a = fetch_field(subsoil, :a, I, domain)
    b = fetch_field(subsoil, :b, I, domain)
    c = fetch_field(subsoil, :c, I, domain)
    precon_values = fetch_field(subsoil, preconsolidation, I, domain)
    precon = preconsolidation(precon_values)

    cells = Vector{DrainingAbcIsotache}()
    for (i, Δz) in enumerate(domain.Δz)
        cell = DrainingAbcIsotache(Δz, γ_wet[i], γ_dry[i], c_d[i], c_v[i], a[i], b[i], c[i])
        push!(cells, cell)
    end

    z = domain.z
    Δz = domain.Δz
    σ = fill(NaN, length(z))
    σ′ = fill(NaN, length(z))
    p = fill(NaN, length(z))
    result = fill(0.0, length(z))

    return ConsolidationColumn(cells, z, Δz, σ, σ′, p, precon, result)
end


"""
    initialize(::Type{NullConsolidation}, domain, subsoil, I)

Initialize an empty ConsolidationColumn (i.e. consolidation is ignored) at location I.
"""
function initialize(::Type{NullConsolidation}, _, domain, _, _)
    cells = Vector{NullConsolidation}()
    for i in 1:length(domain.Δz)
        push!(cells, Atlans.NullConsolidation())
    end
    preconsolidation = OverConsolidationRatio(fill(NaN, length(cells)))
    z = domain.z
    return ConsolidationColumn(
        cells,
        z,
        domain.Δz,
        fill(NaN, length(z)),
        fill(NaN, length(z)),
        fill(NaN, length(z)),
        preconsolidation,
        fill(0.0, length(z))
    )
end


"""
Reset degree of consolidation and time.
"""
function prepare_forcingperiod!(
    column::ConsolidationColumn{DrainingAbcIsotache,P} where {P<:Preconsolidation},
)
    for (i, cell) in enumerate(column.cells)
        cell = @set cell.t = 0.0
        cell = @set cell.U = 0.0
        cell = @set cell.Δz_0 = cell.Δz
        column.cells[i] = cell
    end
    return
end
