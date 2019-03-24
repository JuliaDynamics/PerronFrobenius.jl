tol = 1e-10

@testset "Grid estimator" begin
    points_2D = rand(2, 200)
    points_3D = rand(3, 400)
    E_2D = invariantize(StateSpaceReconstruction.cembed(points_2D))
    E_3D = invariantize(StateSpaceReconstruction.cembed(points_3D))
    ϵ = 3
    bins_visited_by_orbit_2D = assign_bin_labels(E_2D, ϵ)
    bins_visited_by_orbit_3D = assign_bin_labels(E_3D, ϵ)

    bininfo_2D = get_binvisits(bins_visited_by_orbit_2D)
    bininfo_3D = get_binvisits(bins_visited_by_orbit_3D)

    TO_2D = estimate_transferoperator_from_binvisits(bininfo_2D)
    TO_3D = estimate_transferoperator_from_binvisits(bininfo_3D)

    invm_2D = invariantmeasure(TO_2D)
    invm_3D = invariantmeasure(TO_3D)

    @test all(invm_2D.dist .>= -tol)
    @test sum(invm_2D.dist) <= 1 + tol || sum(invm_2D.dist) ≈ 1
end
