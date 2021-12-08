using Test
using Atlans

@testset "Groundwater" begin
    include("test_groundwater.jl")
end

@testset "CarbonStore" begin
    include("test_carbon_store.jl")
end

@testset "ConstantRate" begin
    include("test_constant_rate.jl")
end

@testset "AbcIsotache" begin
    include("test_abc_isotache.jl")
end

@testset "Koppejan" begin
    include("test_koppejan.jl")
end

@testset "ConsolidationColumn" begin
    include("test_consolidation_column.jl")
end
                                               
