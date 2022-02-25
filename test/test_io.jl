@testset "Input-Output" begin
    using DataFrames

    @testset "read_params_table" begin
        filename = AtlansFixtures.params_table()
        df = Atlans.read_params_table(filename)
        @test typeof(df) == DataFrames.DataFrame
        types = eltype.(eachcol(df))
        @test all(types[1:2] .== String)
        @test all(types[3:4] .== Int)
        @test all(types[5:end] .== Float64)
    end

    @testset "lookup_table" begin
        filename = AtlansFixtures.params_table()
        df = Atlans.read_params_table(filename)
        tables = Atlans.build_lookup_tables(df)
        @test typeof(tables) == Dict{Symbol,Dict{Tuple{Int,Int},Float64}}
        @test collect(keys(tables)) == [:a, :b, :c, :γ_dry, :c_v, :γ_wet, :c_d]
        @test tables[:a][(1,2)] == 0.01737
        
        @test Atlans.lookup(tables[:a], [1, 1], [2, 2]) == [0.01737, 0.01737]
    end
    
    @testset "prepare_subsoil_reader" begin
        path_csv = AtlansFixtures.params_table()
        path_nc = AtlansFixtures.subsoil_netcdf()
        @show path_nc
        #reader = Atlans.prepare_subsoil_reader(path_nc, path_csv)
    end

end
