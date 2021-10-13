using Test
using Atlans


@testset "Atlans" begin
    include("test_carbon_store.jl")
end

@testset "Atlans" begin
    include("test_abc_isotache.jl")
end

