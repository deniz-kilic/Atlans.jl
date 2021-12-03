module Atlans

using Accessors
using CFTime
using Dates
using LinearAlgebra
using NCDatasets
using Plots

import RecipesBase.plot

const Float = Float64
const γ_water = 9810.0

abstract type ConsolidationProcess end
abstract type OxidationProcess end
abstract type GroundwaterColumn end
abstract type SoilCell end
abstract type DrainageUnit end

include("soilcolumn/groundwater/groundwater.jl")
include("soilcolumn/groundwater/interpolated.jl")
include("soilcolumn/groundwater/darcy.jl")

include("soilcolumn/consolidation/consolidation.jl")
include("soilcolumn/consolidation/koppejan.jl")
include("soilcolumn/consolidation/abc_isotache.jl")

include("soilcolumn/oxidation/oxidation.jl")
include("soilcolumn/oxidation/carbonstore.jl")

include("soilcolumn/column.jl")
include("soilcolumn/split.jl")

include("drainage_unit.jl")

include("model.jl")
include("io.jl")

end # module
