pts = invariantize([rand(3) for i = 1:15])
triang = triangulate(pts)

μ_exact = invariantmeasure(pts, TriangulationBinning(), ExactIntersection())
μ_approx = invariantmeasure(pts, TriangulationBinning(), ApproximateIntersection())

@test μ_exact isa TransferOperatorTriangulationExact
@test μ_approx isa TransferOperatorTriangulationApprox
