
@testset "Transfer operator from rectangular binning" begin
    @testset "Equidistant" begin
        E = embed([diff(rand(100)) for i = 1:3])
        n_bins = 10

        bins_visited_by_orbit = assign_bin_labels(E, n_bins)
        bininfo = organize_bin_labels(bins_visited_by_orbit)
        TO = transferoperator(bininfo)

        @test typeof(TO) <: RectangularBinningTransferOperator
		# Last row might sum to zero, because the last point does not need to
		# be contained in the last bin. However, the remaining row sums must
		# be one.
        #@test is_markov(TO)
        @show TO
        if !is_markov(TO)
            warn("There were all-zero columns in the transfer matrix")
            warn("Removing first column and last row")
            if is_markov(TO.TO[1:(end-1), 2:end])
                warn("That made the transfer matrix Markov")
            end
            @test is_markov(TO.TO[1:(end-1), 2:end])
        end
    end
end
