import StateSpaceReconstruction

@testset "Triangulation estimators" begin


	pts = invariantize([rand(3) for i = 1:20])

	TOa = transferoperator_triangulation_approx(pts, n_sample_pts = 50)
	TOe = transferoperator_triangulation_exact(pts)

	@test TOa isa TransferOperatorTriangulationApprox
	@test TOe isa TransferOperatorTriangulationExact
    @test minimum(sum(TOa.transfermatrix, dims = 2)) > 0.99
    @test minimum(sum(TOe.transfermatrix, dims = 2)) > 0.99

    # @testset "Triangulation approximate estimator" begin
    #     E = StateSpaceReconstruction.cembed([diff(rand(15)) for i = 1:3])
    #     E_invariant = invariantize(E)
    #
    #     # Triangulations
    #     triang = triangulate(E)
    #     triang_inv = triangulate(E_invariant)
    #
    #     # Transfer operators from *invariant* triangulations
    #     TO = transferoperator_triang(triang_inv)
    #     TO_approx = transferoperator_triang(triang_inv, exact = false, parallel = false)
    #     TO_approx_rand = transferoperator_triang(triang_inv, exact = false, parallel = false, sample_randomly = true)
    #     invm_approx = invariantmeasure(TO_approx)
    #     invm_approx_rand = invariantmeasure(TO_approx_rand)
    #
    #     @test all(invm_approx.dist .>=  -tol)
    #     @test all(invm_approx_rand.dist .>= -tol)
    #
    #     @test sum(invm_approx.dist) ≈ 1
    #     @test sum(invm_approx_rand.dist) ≈ 1
    # end


end
