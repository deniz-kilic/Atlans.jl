
abstract type Koppejan <: ConsolidationProcess end

struct DrainingKoppejan
    Δz::Float64
    t::Float64
    σ′::Float64  # effective stress
    γ_w::Float64  # wet specific mass
    γ_d::Float64  # dry specific mass
    c_d::Float64  # drainage coefficient
    c_v::Float64  # drainage coefficient
    U::Float64
    Cp::Float64
    Cs::Float64
    Cp′::Float64
    Cs′::Float64
    σ′pre::Float64
end

function consolidate(
    kpj::DrainingKoppejan,
    σ′::Float64,
    Δt::Float64,
)::Tuple{Float64,DrainingKoppejan}
    t = kpj.t + Δt
    # Degree of consolidation changes
    U = U(kpj, t)
    ΔU = U - kpj.U
    # Effective stress changes
    Δσ′ = σ′ - kpj.σ′
    # Exit early for a number of cases Koppejan cannot deal with.
    if (kpj.Δz == 0) | (kpj.Cp == 0) | (kpj.σ′ == 0) | (Δσ′ < 0)
        consolidation = 0.0
        return consolidation,
        DrainingKoppejan(
            kpj.Δz,  # no change
            t,  # new
            σ′, # new
            kpj.γ_w,
            kpj.γ_d,
            kpj.c_d,
            kpj.c_v,
            U,  # new (?)
            kpj.Cp,
            kpj.Cs,
            kpj.Cp′,
            kpj.Cs′,
            kpj.σ′pre,
        )
    end

    # pre-consolidation
    strain = ΔU / kpj.Cp
    if kpj.Cs > 0
        strain += (1.0 / kpj.Cs) * log10(t / kpj.t) * log(σ′ / kpj.σ′)
    end
    # post-consolidation
    if σ′ > kpj.σ′pre
        strain += ΔU / kpj.Cp′
        if kpj.Cs′ > 0
            strain += (1.0 / kpj.Cs′) * log10(t / kpj.t) * log(σ′ / kpj.σ′pre)
        end
    end
    # consolidation changes
    consolidation = min(Δz, strain * kpj.Δz)
    γ_w = Atlans.compress_γ_wet(kpj, consolidation)
    γ_d = Atlans.compress_γ_dry(kpj, consolidation)
    # return new state
    return consolidation,
    DrainingKoppejan(
        kpj.Δz - consolidation,  # new
        t,  # new
        σ′, # new
        γ_w,  # new
        γ_d,  # new
        kpj.c_d,
        kpj.c_v,
        U,  # new
        kpj.Cp,
        kpj.Cs,
        kpj.Cp′,
        kpj.Cs′,
        kpj.σ′pre,
    )
end
