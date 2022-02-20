@testset "TimeSteppers" begin
    @testset "constructor" begin
        timestepper = Atlans.ExponentialTimeStepper(1.0, 2)
        @test typeof(timestepper) == Atlans.ExponentialTimeStepper{Int}
    end
    
    @testset "create_timesteps" begin
        timestepper = Atlans.ExponentialTimeStepper(1.0, 2)

        steps = Atlans.create_timesteps(timestepper, 3.0)
        @test all(steps .== [1.0, 2.0])
        
        steps = Atlans.create_timesteps(timestepper, 4.0)
        @test all(steps .== [1.0, 2.0, 1.0])

        steps = Atlans.create_timesteps(timestepper, 10.0)
        @test all(steps == [1.0, 2.0, 4.0, 3.0])
    end
end