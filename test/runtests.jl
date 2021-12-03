using Test
using Atlans

@testset "CarbonStore" begin
    include("test_carbon_store.jl")
end

@testset "AbcIsotache" begin
    include("test_abc_isotache.jl")
end

