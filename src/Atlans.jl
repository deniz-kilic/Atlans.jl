module Atlans

using CSV
using Setfield
using CFTime
using Dates
using LinearAlgebra
using DataFrames
using NCDatasets
using StatsBase
using ProgressLogging
using LoggingExtras
using TerminalLoggers

const Float = Float64
const OptionalFloat = Union{Float,Missing}
const OptionalInt = Union{Int,Missing}
const γ_water = 9810.0
const τ_ref = 1.0

abstract type ConsolidationProcess end
abstract type OxidationProcess end
abstract type GroundwaterColumn end
abstract type ShrinkageProcess end

struct NullConsolidation end
struct NullOxidation end
struct NullShrinkage end

abstract type Preconsolidation end
abstract type AbstractAbcIsotache <: ConsolidationProcess end


include("soilcolumn/groundwater/groundwater.jl")
include("soilcolumn/groundwater/hydrostatic.jl")
#include("soilcolumn/groundwater/interpolated.jl")
#include("soilcolumn/groundwater/darcy.jl")

include("soilcolumn/consolidation/consolidation.jl")
#include("soilcolumn/consolidation/koppejan.jl")
include("soilcolumn/consolidation/abc_isotache.jl")
include("soilcolumn/consolidation/draining_abc_isotache.jl")

include("soilcolumn/oxidation/oxidation.jl")
include("soilcolumn/oxidation/carbonstore.jl")
#include("soilcolumn/oxidation/constant_rate.jl")

include("soilcolumn/column.jl")
include("soilcolumn/split.jl")

include("timesteppers.jl")
include("simulation.jl")
include("forcing.jl")
include("io.jl")
include("logging.jl")

end # module
