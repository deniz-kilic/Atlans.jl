
abstract type Koppejan <: ConsolidationProcess end

struct DrainingKoppejan <: Koppejan
    Δz::Float
    Δz_0::Float
    t::Float
    σ′::Float  # effective stress
    γ_wet::Float  # wet specific mass
    γ_dry::Float  # dry specific mass
    c_d::Float  # drainage coefficient
    c_v::Float  # drainage coefficient
    U::Float
    Cp::Float
    Cs::Float
    Cp′::Float
    Cs′::Float
    σ′pre::Float
    consolidation::Float
end

function consolidate(kpj::DrainingKoppejan, σ′::Float, Δt::Float)::DrainingKoppejan
    t = kpj.t + Δt
    # Degree of consolidation changes
    U = Atlans.U(kpj, t)
    ΔU = U - kpj.U
    # Effective stress changes
    Δσ′ = σ′ - kpj.σ′
    # Exit early for a number of cases Koppejan cannot deal with.
    if (kpj.Δz == 0) | (kpj.Cp == 0) | (kpj.σ′ == 0) | (Δσ′ < 0)
        consolidation = 0.0
        return DrainingKoppejan(
            kpj.Δz,  # no change
            kpj.Δz_0,  # no change
            t,  # new
            σ′, # new
            kpj.γ_wet,
            kpj.γ_dry,
            kpj.c_d,
            kpj.c_v,
            U,  # new (?)
            kpj.Cp,
            kpj.Cs,
            kpj.Cp′,
            kpj.Cs′,
            kpj.σ′pre,
            consolidation,
        )
    end

    log10time = log10((1.0 + t) / (1.0 + kpj.t))
    if σ′ <= kpj.σ′pre
        strain = ((ΔU / kpj.Cp) + (1.0 / kpj.Cs) * log10time) * log(σ′ / kpj.σ′)
    else
        strain = (
            ((ΔU / kpj.Cp) + (1.0 / Cs) * log10time) * log(σ′pre / kpj.σ′) +
            ((ΔU / kpj.Cp′) + (1.0 / Cs′) * log10time) * log(σ′ / kpj.σ′pre)
        )
    end

    # consolidation changes
    consolidation = strain * kpj.Δz  # natural strain
    γ_wet = Atlans.compress_γ_wet(kpj, consolidation)
    γ_dry = Atlans.compress_γ_dry(kpj, consolidation)
    # return new state
    return DrainingKoppejan(
        kpj.Δz - consolidation,  # new
        kpj.Δz_0,  # no change
        t,  # new
        σ′, # new
        γ_wet,  # new
        γ_dry,  # new
        kpj.c_d,
        kpj.c_v,
        U,  # new
        kpj.Cp,
        kpj.Cs,
        kpj.Cp′,
        kpj.Cs′,
        kpj.σ′pre,
        consolidation,  # new
    )
end
