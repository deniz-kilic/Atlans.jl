struct ConsolidationSurcharge{C,P} <: AbstractConsolidationColumn
    cells::Vector{C}
    z::Vector{Float}
    Δz::Vector{Float}
    σ::Vector{Float}  # Total stress in Pa
    σ′::Vector{Float}  # Effective stress in Pa
    p::Vector{Float}  # Pore pressure in Pa
    preconsolidation::P
end


struct OxidationSurcharge{O}
    cells::Vector{O}
    z::Vector{Float}
    Δz::Vector{Float}
end


struct ShrinkageSurcharge{S}
    cells::Vector{S}
    z::Vector{Float}
    Δz::Vector{Float}
end


struct GroundwaterSurcharge <: HydrostaticColumn
    z::Vector{Float}
    phreatic::Phreatic
    dry::Vector{Bool} # Portion of column above the groundwater table
    p::Vector{Float} # pore pressure, is always dry above gw table
end


struct SurchargeColumn
    z::Vector{Float}
    Δz::Vector{Float}
    groundwater::GroundwaterSurcharge
    consolidation::ConsolidationSurcharge
    oxidation::OxidationSurcharge
    shrinkage::ShrinkageSurcharge
end


function initialize(
    ::Type{HydrostaticGroundwater},
    phreatic::Phreatic,
    domain::VerticalDomain
)
    GroundwaterSurcharge(
        domain.z,
        phreatic,
        fill(true, length(domain.z)),
        fill(0.0, length(domain.z))
    )
end


function initialize(
    ::Type{DrainingAbcIsotache},
    preconsolidation::Type,
    domain::VerticalDomain,
    lookup_table::Dict,
)
    γ_wet = fetch_field(lookup_table, :gamma_wet, domain)
    γ_dry = fetch_field(lookup_table, :gamma_dry, domain)
    c_d = fetch_field(lookup_table, :drainage_factor, domain)
    c_v = fetch_field(lookup_table, :c_v, domain)
    a = fetch_field(lookup_table, :a, domain)
    b = fetch_field(lookup_table, :b, domain)
    c = fetch_field(lookup_table, :c, domain)
    precon_values = fetch_field(lookup_table, preconsolidation, domain)
    precon = preconsolidation(precon_values)

    cells = Vector{DrainingAbcIsotache}()
    for (i, Δz) in enumerate(domain.Δz)
        cell = DrainingAbcIsotache(Δz, γ_wet[i], γ_dry[i], c_d[i], c_v[i], a[i], b[i], c[i])
        push!(cells, cell)
    end

    σ = fill(NaN, length(cells))
    σ′ = fill(NaN, length(cells))
    p = fill(NaN, length(cells))

    ConsolidationSurcharge(
        cells, domain.z, domain.Δz, σ, σ′, p, precon)
end


function initialize(
    ::Type{CarbonStore},
    domain::VerticalDomain,
    lookup_table::Dict,
)
    f_minimum_organic = fetch_field(lookup_table, :minimal_mass_fraction_organic, domain)
    α = fetch_field(lookup_table, :oxidation_rate, domain)

    cells = Vector{CarbonStore}()
    for (i, Δz) in enumerate(domain.Δz)
        cell = CarbonStore(Δz, 0.0, f_minimum_organic[i], 0.0, α[i])
        push!(cells, cell)
    end

    OxidationSurcharge(cells, domain.z, domain.Δz)
end


function initialize(
    ::Type{SimpleShrinkage},
    domain::VerticalDomain,
    lookup_table::Dict,
)
    n = 0.5 # minimum shrinkage degree, no shrinkage occurs
    mass_lutum = 0.0
    mass_organic = 0.0
    cells = [SimpleShrinkage(Δz, n, mass_lutum, mass_organic) for Δz in domain.Δz]
    ShrinkageSurcharge(cells, domain.z, domain.Δz)
end


function initialize(::Type{NullConsolidation}, domain, _)
    n = length(domain.z)
    cells = fill(NullConsolidation(), length(n))
    preconsolidation = OverConsolidationRatio(fill(NaN, length(cells)))

    ConsolidationSurcharge(
        cells,
        domain.z,
        domain.Δz,
        fill(NaN, length(n)),
        fill(NaN, length(n)),
        fill(NaN, length(n)),
        preconsolidation
    )
end


function initialize(::Type{NullOxidation}, domain, _)
    cells = fill(NullOxidation(), length(domain.z))
    OxidationSurcharge(cells, domain.z, domain.Δz)
end


function initialize(::Type{NullShrinkage}, domain, _)
    cells = fill(NullShrinkage(), length(domain.z))
    ShrinkageSurcharge(cells, domain.z, domain.Δz)
end


function prepare_timestep!(column::SurchargeColumn, Δt)
    flow!(column.groundwater, Δt)
    exchange_pore_pressure!(column)
    total_stress!(column.consolidation, phreatic_level(column.groundwater))
    effective_stress!(column.consolidation)
end


function exchange_pore_pressure!(column::SurchargeColumn)
    column.consolidation.p .= column.groundwater.p * γ_water
end


function apply_preconsolidation!(
    column::ConsolidationSurcharge{
        ABC,
        OverConsolidationRatio,
    } where {ABC<:AbstractAbcIsotache},
)
    for (i, cell) in enumerate(column.cells)
        column.cells[i] = set_τ0_ocr(cell, column.preconsolidation.ratio[i])
    end
    return
end


function apply_preconsolidation!(
    column::ConsolidationSurcharge{
        ABC,
        PreOverburdenPressure
    } where {ABC<:AbstractAbcIsotache},
)
    for (i, cell) in enumerate(column.cells)
        column.cells[i] = set_τ0_pop(cell, column.preconsolidation.pressure[i])
    end
    return
end