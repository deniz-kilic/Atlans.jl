function set_surcharge!(
   column::SoilColumn,
   surcharge::SurchargeColumn,
)
   set_surcharge!(column.groundwater, surcharge.groundwater)
   set_surcharge!(column.consolidation, surcharge.consolidation)
   set_surcharge!(column.oxidation, surcharge.oxidation)
   set_surcharge!(column.shrinkage, surcharge.shrinkage)
   return
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

    surcol = SurchargeColumn(groundwater, consolidation, oxidation, shrinkage)
    apply_preconsolidation!(surcol)
    return surcol
end
