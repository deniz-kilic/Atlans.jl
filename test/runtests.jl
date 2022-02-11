using Revise
using Test
using Atlans


@testset "Groundwater" begin
    include("test_hydrostatic_groundwater.jl")
end

@testset "CarbonStore" begin
    include("test_carbon_store.jl")
end

#@testset "ConstantRate" begin
#    include("test_constant_rate.jl")
#end

@testset "DrainingAbcIsotache" begin
    include("test_draining_abc_isotache.jl")
end

##@testset "Koppejan" begin
##    include("test_koppejan.jl")
##end

@testset "ConsolidationColumn" begin
    include("test_consolidation_column.jl")
end

@testset "OxidationColumn" begin
    include("test_oxidation_column.jl")
end
