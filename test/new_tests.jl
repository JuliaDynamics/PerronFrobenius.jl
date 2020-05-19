using CausalityToolsBase
using PerronFrobenius
using DelayEmbeddings
using Simplices 
using Test 

@testset "Grid estimators" begin 
    n = 200
    x = [cos(2π*(t + sin(t)/30)) for t in 1:n] .+ 0.2rand(n)
    push!(x, 5)
    τs = (0, -2, -3)
    js = (1, 1, 1)
    D = genembed(x, τs, js)

    @testset "SingleGrid" begin 
        @test SingleGrid(RectangularBinning(5)) isa SingleGrid
        @test SingleGrid(RectangularBinning(5), "circular") isa SingleGrid

        tog = transferoperatorgenerator(D, SingleGrid(RectangularBinning(5)))
        M = tog()
        all([M[i, :] |> sum .≈ 1.0 for i = 1:size(M, 2)])
    end
end
@testset "Triangulation estimators" begin 
    n = 20
    x = [cos(2π*(t + sin(t)/30)) for t in 1:n] .+ 0.2rand(n)
    push!(x, 5)
    τs = (0, -2, -3)
    js = (1, 1, 1)
    D = genembed(x, τs, js)

    @testset "SimplexExact" begin
        tog = transferoperatorgenerator(D, SimplexExact())
        M = tog()
        @test tog isa PerronFrobenius.TransferOperatorGenerator
        @test all(sum(M, dims = 2) .≈ 1.0)
    end

    @testset "SimplexApprox" begin
        tog = transferoperatorgenerator(D, SimplexApprox())
        M = tog(n = 60)
        @test tog isa PerronFrobenius.TransferOperatorGenerator
        @test all(sum(M, dims = 2) .≈ 1.0)
    end
end