using DelayEmbeddings
using PerronFrobenius
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

        method = SingleGrid(RectangularBinning(5))
        tog = transopergenerator(D, method)
        TO = tog()
        E = TransferOperatorGenerator{SingleGrid{RectangularBinning}}
        @test TO isa TransferOperatorApproximation{<:E}
        M = TO.M
        @test all([M[i, :] |> sum .≈ 1.0 for i = 1:size(M, 2)])

        @test transferoperator(D, method) isa TransferOperatorApproximation{<:E}
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
        tog = transopergenerator(D, SimplexExact())
        TO = tog()
        E = TransferOperatorGenerator{SimplexExact}
        @test TO isa TransferOperatorApproximation{<:E}
        M = TO.M
        @test all(sum(M, dims = 2) .≈ 1.0)
    end

    @testset "SimplexApprox" begin
        tog = transopergenerator(D, SimplexPoint())
        TO = tog(n = 50)
        E = TransferOperatorGenerator{SimplexPoint}
        @test TO isa TransferOperatorApproximation{<:E}
        M = TO.M
        @test all(sum(M, dims = 2) .≈ 1.0)
    end
end

@testset "Invariant measures" begin 

n = 200
x = [cos(2π*(t + sin(t)/30)) for t in 1:n] .+ 0.2rand(n)
push!(x, 5)
τs = (0, -2, -3)
js = (1, 1, 1)
D = genembed(x, τs, js)

@testset "Grid estimators" begin 
    method = SingleGrid(RectangularBinning(5))
    to = transferoperator(D, method)
    iv = invariantmeasure(to)
    @test iv isa InvariantDistribution
end

@testset "Triangulation estimators" begin 
n = 20
x = [cos(2π*(t + sin(t)/30)) for t in 1:n] .+ 0.2rand(n)
push!(x, 5)
τs = (0, -2, -3)
js = (1, 1, 1)
D = genembed(x, τs, js)

    @testset "SimplexPoint" begin 
        method = SimplexPoint()
        to = transferoperator(D, method)
        iv = invariantmeasure(to)
        @test invariantmeasure(to) isa InvariantDistribution
        @test invariantmeasure(D, method) isa InvariantDistribution
    end

    @testset "SimplexExact" begin 
        method = SimplexExact()
        to = transferoperator(D, method)
        iv = invariantmeasure(to)
        @test invariantmeasure(to) isa InvariantDistribution
        @test invariantmeasure(D, method) isa InvariantDistribution
    end

end

end