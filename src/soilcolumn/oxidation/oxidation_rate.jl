"""
    relative_oxidation_rate(T::Float)

Empirical relation between the decay rate of organic soils and air temperature in
Celsius from Hendriks and Vermeulen (1997). Relation is valid between 0° and 19.5°
Celsius.
"""
function relative_oxidation_rate(T::Float)
    krel = 0.264ℯ^(0.127T) / 1.25
    return krel
end


"""
    update_alpha(cell::CarbonStore, T::Float)

Return a new CarbonStore cell with an oxidation rate (α) that is corrected for
air temperature.
"""
function update_alpha(cell::CarbonStore, T::Float)
    krel = relative_oxidation_rate(T)
    new_alpha = cell.α0 * krel
    newcell = @set cell.α = new_alpha
    return newcell
end
