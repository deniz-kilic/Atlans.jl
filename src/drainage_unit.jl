
struct Weir
    x::Float
    y::Float
    level::Float
end

struct SimpleDrainageUnit
    surfacelevel::Vector{Float}
    columns::Vector{SoilColumn}
end

struct WeirDrainageUnit
    columns::Vector{SoilColumn}
    weirs::Vector{Weir}
end
