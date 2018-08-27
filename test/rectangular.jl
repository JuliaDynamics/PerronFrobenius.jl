
@testset "Transfer operator from rectangular binning" begin
    @testset "Equidistant" begin
        E = embed([diff(rand(15)) for i = 1:3])
        E_invariant = invariantize(E)
        n_bins = 10

        equibin = bin_equidistant(E, n_bins)
        equibin_inv = bin_equidistant(E_invariant, n_bins)

        TO2 = transferoperator(equibin)
        TO3 = transferoperator(equibin_inv)

        @test typeof(TO2) <: RectangularBinningTransferOperator
        @test typeof(TO3) <: RectangularBinningTransferOperator

		# Last row might sum to zero, because the last point does not need to
		# be contained in the last bin. However, the remaining row sums must
		# be one.

    end
end
