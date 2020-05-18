using CausalityToolsBase
using PerronFrobenius

# Grid estimator
@test Grid(RectangularBinning(5)) isa Grid
@test Grid(RectangularBinning(5), "circular") isa Grid

# Simplex estimators
n = 20
x = [cos(2π*(t + sin(t)/30)) for t in 1:n] .+ 0.2rand(n)
push!(x, 5)
τs = (0, -2, -3)
js = (1, 1, 1)
D = genembed(x, τs, js)

tog = transferoperatorgenerator(D, SimplexExact())
M = tog()

@test tog isa TransferOperatorGenerator
@test all(sum(M, dims = 2) .≈ 1.0)
