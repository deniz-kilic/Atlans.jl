struct OxidationColumn
    cells::Vector{O} where {O<:OxidationProcess}
    z::Vector{Float}
    Δz::Vector{Float}
end
