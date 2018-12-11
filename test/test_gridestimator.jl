tol = 1e-10

@testset "Grid estimator" begin
    points_2D = rand(2, 200)
    points_3D = rand(3, 400)
    E_2D = invariantize(StateSpaceReconstruction.embed(points_2D))
    E_3D = invariantize(StateSpaceReconstruction.embed(points_3D))
    ϵ = 3
    bins_visited_by_orbit_2D = assign_bin_labels(E_2D, ϵ)
    bins_visited_by_orbit_3D = assign_bin_labels(E_3D, ϵ)

    bininfo_2D = organize_bin_labels(bins_visited_by_orbit_2D)
    bininfo_3D = organize_bin_labels(bins_visited_by_orbit_3D)

    TO_2D = transferoperator_binvisits(bininfo_2D)
    TO_3D = transferoperator_binvisits(bininfo_3D)

    invm_2D = left_eigenvector(TO_2D)
    invm_3D = left_eigenvector(TO_3D)

    @test all(invm_2D.dist .>= -tol)
    @test sum(invm_2D.dist) <= 1 + tol || sum(invm_2D.dist) ≈ 1
end
