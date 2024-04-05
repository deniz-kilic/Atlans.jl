@testset "Surcharge" begin
    @testset "prepare_domain" begin
        Δzmax = 0.25
        
        domain = Atlans.prepare_domain(1., 2, Δzmax)
        @test domain.n == 4
        @test all(domain.Δz .== 0.25)
        @test all(isnan.(domain.z))
        @test all(domain.lithology .== 2)

        domain = Atlans.prepare_domain([0.1, 0.2, 0.4, 0.1], [2, 2, 3, 2], Δzmax)
        @test domain.n == 5
        @test all(domain.Δz .== [0.1, 0.2, 0.2, 0.2, 0.1])
        @test all(isnan.(domain.z))
        @test all(domain.lithology .== [2, 2, 3, 3, 2])

    end
end