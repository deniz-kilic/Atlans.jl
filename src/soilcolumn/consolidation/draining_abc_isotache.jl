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

function DrainingAbcIsotache(
    Δz,
    γ_wet,
    γ_dry,
    c_d,
    c_v,
    a,
    b,
    c,
)
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

function consolidate(abc::DrainingAbcIsotache, σ′::Float, Δt::Float)
    t = abc.t + Δt

    # Degree of consolidation changes
    Unew = U(abc, t)
    ΔU = Unew - abc.U

    # Effective stress changes
    load = σ′ - abc.σ′
    σ′ = abc.σ′ + Unew * load
    loadstep = ΔU * load

    # τ changes
    τ⃰ = abc.τ * ((abc.σ′ - loadstep) / abc.σ′) ^ ((abc.b - abc.a) / abc.c)
    τ = τ⃰ + Δt

    # consolidation
    elastoplastic = abc.a * log(σ′ / (σ′ - loadstep))
    creep = abc.c * log(τ / τ⃰)
    strain = elastoplastic + creep

    # Thickness should not go below 0
    consolidation = min(abc.Δz, strain * abc.Δz)
    γ_wet = compress_γ_wet(abc, consolidation)
    γ_dry = compress_γ_dry(abc, consolidation)
    
    # return new state
    return DrainingAbcIsotache(
        abc.Δz - consolidation,  # new
        abc.Δz_0,
        t,  # new
        abc.σ′,
        γ_wet,  # new
        γ_dry,  # new
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

function initialize(::DrainingAbcIsotache, preconsolidation::Type, domain, reader, I)::ConsolidationColumn{DrainingAbcIsotache}
    use = domain.use
    lithology = domain.lithology
    geology = domain.geology

    γ_wet = @view fetch_field(reader, :γ_wet, I, lithology, geology)[use]
    γ_dry = @view fetch_field(reader, :γ_dry, I, lithology, geology)[use]
    c_d = @view fetch_field(reader, :c_d, I, lithology, geology)[use]
    c_v = @view fetch_field(reader, :c_v, I, lithology, geology)[use]
    a = @view fetch_field(reader, :a, I, lithology, geology)[use]
    b = @view fetch_field(reader, :b, I, lithology, geology)[use]
    c = @view fetch_field(reader, :c, I, lithology, geology)[use]
    precon = @view fetch_field(reader, preconsolidation, I, lithology, geology)[use] 
    
    cells = Vector{DrainingAbcIsotache}()
    for (i, Δz) in zip(domain.index, domain.Δz)
        cell = DrainingAbcIsotache(
            Δz,
            γ_wet[i],
            γ_dry[i],
            c_d[i],
            c_v[i],
            a[i],
            b[i],
            c[i],
        )
        push!(cells, cell)
    end
    
    return ConsolidationColumn(
        cells,
        z,
        Δz,
        σ,
        σ′,
        p,
        precon,
        result,
    )
end

"""
Reset degree of consolidation and time.
"""
function prepare_forcingperiod!(column::ConsolidationColumn{DrainingAbcIsotache,P} where {P <: Preconsolidation})
    for (i, cell) in enumerate(column.cells)
        cell = @set cell.t = 0.0  
        cell = @set cell.U = 0.0 
        cell = @set cell.Δz_0 = cell.Δz   
        column.cells[i] = cell
    end
    return
end
