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
    z = 0 .+ cumsum(Δz) .- 0.5 .* Δz
    geology = fill(1, ncells)
    lithology = fill(lithology, ncells)
    return VerticalDomain(z, Δz, geology, lithology, index, ncells)
end


function prepare_domain(
    thickness::Vector{OptionalFloat},
    lithology::Vector{OptionalInt},
    Δzmax::Float
)   
    is_values = .!ismissing.(thickness)
    thickness = thickness[is_values]
    lithology = lithology[is_values]

    index, ncells = discretize(thickness, Δzmax)
    Δz = (thickness./ncells)[index]
    z = 0 .+ cumsum(Δz) .- 0.5 .* Δz
    geology = fill(1, sum(ncells))
    n = sum(ncells)
    return VerticalDomain(z, Δz, geology, lithology[index], index, n)
end
