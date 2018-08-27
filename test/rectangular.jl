
@testset "Transfer operator from rectangular binning" begin
    @testset "Equidistant" begin
        E = embed([diff(rand(15)) for i = 1:3])
        E_invariant = invariantize(E)
        n_bins = 10

        bins_visited_by_orbit = assign_bin_labels(E, n_bins)
        bins_visited_by_orbit_inv = assign_bin_labels(E_invariant, n_bins)
        bininfo = organize_bin_labels(bins_visited_by_orbit)
        bininfo_inv = organize_bin_labels(bins_visited_by_orbit_inv)
        TO2 = transferoperator(bininfo)
        TO3 = transferoperator(bininfo_inv)

        @test typeof(TO2) <: RectangularBinningTransferOperator
        @test typeof(TO3) <: RectangularBinningTransferOperator

		# Last row might sum to zero, because the last point does not need to
		# be contained in the last bin. However, the remaining row sums must
		# be one.

    end
end
