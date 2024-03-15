struct SurchargeColumn{G,C,O,S}
    Δz::Vector{Float}
    groundwater::G
    consolidation::ConsolidationColumn{C}
    oxidation::OxidationColumn{O}
    shrinkage::ShrinkageColumn{S}
end


struct SurchargeGroundwater
    z::Vector{Float}
    dry::Vector{Bool}
    p::Vector{Float}
end


struct SurchargeConsolidation
end


struct OxidationSurcharge
end


struct ShrinkageSurcharge
end


function initialize(
    ::Type{DrainingAbcIsotache},
    domain::VerticalDomain,
    lookup::Dict,
)
end


function initialize(
    ::Type{CarbonStore},
    domain::VerticalDomain,
    lookup::Dict,
)
end


function initialize(
    ::Type{SimpleShrinkage},
    domain::VerticalDomain,
    lookup::Dict,
)
end
 