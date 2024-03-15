struct SurchargeColumn{G,C,O,S}
   Δz::Vector{Float}
   groundwater::G
   consolidation::ConsolidationColumn{C}
   oxidation::OxidationColumn{O}
   shrinkage::ShrinkageColumn{S}
end


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


function prepare_loookup_table(path_table)
    df = read_params_table(path_table)
    build_lookup_tables(df)
end


function prepare_domain(thickness::Float, lithology::Int, Δzmax)
    index, ncells = discretize(thickness, Δzmax)
    Δz = fill(thickness/ncells, ncells)
    z = 0 .+ cumsum(Δz) .- 0.5 .* Δz
    geology = fill(1, ncells)
    lithology = fill(lithology, ncells)
    return VerticalDomain(z, Δz, geology, lithology, index, ncells)
end


function prepare_domain(thickness, lithology, Δzmax)
    index, ncells = discretize(thickness, Δzmax)
    Δz = (thickness./ncells)[index]
    z = 0 .+ cumsum(Δz) .- 0.5 .* Δz
    geology = fill(1, sum(ncells))
    n = sum(ncells)
    return VerticalDomain(z, Δz, geology, lithology[index], index, n)
end
