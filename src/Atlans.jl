module Atlans

using CFTime
using Dates
using LinearAlgebra
using NCDatasets
using Plots

import RecipesBase.plot

const Float = Float64

abstract type ConsolidationProcess end
abstract type OxidationProcess end
abstract type GroundwaterProcess end
abstract type SoilCell end
abstract type DrainageUnit end

struct NullConsolidation <: ConsolidationProcess end
struct NullOxidation <: OxidationProcess end
struct NullGroundwater <: GroundwaterProcess end

include("soilcolumn/column.jl")
include("soilcolumn/groundwater.jl")
include("soilcolumn/consolidation.jl")
include("soilcolumn/koppejan.jl")
include("soilcolumn/abc_isotache.jl")
include("soilcolumn/carbonstore.jl")

include("drainage_unit.jl")

include("model.jl")
include("io.jl")

end # module
