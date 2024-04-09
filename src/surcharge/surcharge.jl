function set_surcharge!(
    column::SoilColumn,
    surcharge::SurchargeColumn,
)
    append!(column.z, surcharge.z)
    append!(column.Δz, surcharge.Δz)

    set_surcharge!(column.groundwater, surcharge.groundwater)
    set_surcharge!(column.consolidation, surcharge.consolidation)
    set_surcharge!(column.oxidation, surcharge.oxidation)
    set_surcharge!(column.shrinkage, surcharge.shrinkage)
end


function prepare_lookup_table(path_table)
    df = read_params_table(path_table)
    build_lookup_tables(df)
end


function prepare_domain(thickness::Float, lithology::Int, Δzmax::Float)
    index, ncells = discretize(thickness, Δzmax)
    Δz = fill(thickness/ncells, ncells)
    z = fill(NaN, ncells)
    geology = fill(1, ncells)
    lithology = fill(lithology, ncells)
    return VerticalDomain(z, Δz, geology, lithology, index, ncells)
end


function prepare_domain(thickness::Vector, lithology::Vector, Δzmax::Float)   
    is_values = .!ismissing.(thickness)
    thickness = thickness[is_values]
    lithology = lithology[is_values]

    index, ncells = discretize(thickness, Δzmax)
    Δz = (thickness./ncells)[index]
    z = fill(NaN, sum(ncells))
    geology = fill(1, sum(ncells))
    n = sum(ncells)
    return VerticalDomain(z, Δz, geology, lithology[index], index, n)
end


function prepare_surcharge_column(sur::Surcharge, column::SoilColumn, I::CartesianIndex)
    domain = prepare_domain(sur.thickness[I], sur.lithology[I], 0.25)
    domain.z .= surface_level(column) .+ cumsum(domain.Δz) .- 0.5 .* domain.Δz
    
    groundwater = initialize(typeof(column.groundwater), column.groundwater.phreatic, domain)
    oxidation = initialize(eltype(column.oxidation.cells), domain, sur.lookup)
    shrinkage = initialize(eltype(column.shrinkage.cells), domain, sur.lookup)

    consolidation = initialize(
        eltype(column.consolidation.cells),
        typeof(column.consolidation.preconsolidation),
        domain,
        sur.lookup
    )

    surcol = SurchargeColumn(
        domain.z,
        domain.Δz,
        groundwater,
        consolidation,
        oxidation,
        shrinkage
    )
    apply_preconsolidation!(surcol)
    return surcol
end


function set_surcharge!(gw::HydrostaticGroundwater, sugw::GroundwaterSurcharge)
    # append!(gw.z, sugw.z)
    append!(gw.dry, sugw.dry)
    append!(gw.p, sugw.p)
end


function set_surcharge!(con::ConsolidationColumn, sucon::ConsolidationSurcharge)
    append!(con.cells, sucon.cells)
    # append!(con.z, sucon.z)
    # append!(con.Δz, sucon.Δz)
    append!(con.σ, sucon.σ)
    append!(con.σ′, sucon.σ′)
    append!(con.p, sucon.p)
    append_preconsolidation!(con.preconsolidation, sucon.preconsolidation)
end


function set_surcharge!(ox::OxidationColumn, suox::OxidationSurcharge)
    append!(ox.cells, suox.cells)
    # append!(ox.z, suox.z)
    # append!(ox.Δz, suox.Δz)
end


function set_surcharge!(shr::ShrinkageColumn, sushr::ShrinkageSurcharge)
    append!(shr.cells, sushr.cells)
    # append!(shr.z, sushr.z)
    # append!(shr.Δz, sushr.Δz)
end


function append_preconsolidation!(
    ocr1::OverConsolidationRatio,
    ocr2::OverConsolidationRatio
)
    append!(ocr1.ratio, ocr2.ratio)
end


function append_preconsolidation!(
    ocr1::PreOverburdenPressure,
    ocr2::PreOverburdenPressure
)
    append!(ocr1.pressure, ocr2.pressure)
end
