pts = invariantize([rand(3) for i = 1:15])
triang = triangulate(pts)

to_approx =  transferoperator(pts, TriangulationBinning(), ExactIntersection())
to_exact = transferoperator(pts, TriangulationBinning(), ApproximateIntersection())

@test to_approx isa TransferOperatorTriangulationApprox
@test to_exact isa TransferOperatorTriangulationExact

@test all(sum(to_approx.transfermatrix, dims = 2) .≈ 1.0)
@test all(sum(to_exact.transfermatrix, dims = 2) .≈ 1.0)