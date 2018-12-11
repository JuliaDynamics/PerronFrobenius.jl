pts = rand(3, 10000)
μF = [4, 4, 4]

# Measure induced at resolution μF by the same resolution

# Accepted discrepancy between induced and original measures is sometimes
# nonzero, because the random starting distribution is *random* every
# time, and will yield slightly different measures every time.

μ_induced1a = PerronFrobenius.μϵF_induced_by_ϵj(pts, μF , μF)
μ_induced1b = PerronFrobenius.induced_measure(pts, μF , μF)

μ_induced2a = PerronFrobenius.μϵF_induced_by_ϵj(pts, μF , [2, 6, 4])
μ_induced2b = PerronFrobenius.induced_measure(pts, μF , [2, 6, 4])

@testset "Induced measure" begin
	@testset "Induced measures sum to 1" begin
		# Test both the precise-name function and the alias.
		@test abs(sum(μ_induced1a) - 1) < 1e-7
		@test abs(sum(μ_induced1b) - 1) < 1e-7
		@test abs(sum(μ_induced2a) - 1) < 1e-7
		@test abs(sum(μ_induced2b) - 1) < 1e-7
	end

	@testset "Induced measure at resultions differ from regular measure" begin
		# If we obtain the same measure by inducing the measure to a final
		# parition from a different partition, then something is very wrong.
		# Check that this is not the case.
		@test !(all(μ_induced1a .≈ μ_induced2a))
		@test !(all(μ_induced1b .≈ μ_induced2b))
	end

	@testset "Comparing regular to induced measure" begin
		# Measure computed regularly, either manually or using the
		# constructor that also stored info.
		TO_orig1 = PerronFrobenius.transferoperator_grid(pts, μF)
		μ_orig1 = PerronFrobenius.invariantmeasure(TO_orig1);
		μ_orig2 = PerronFrobenius.RectangularInvariantMeasure(pts, μF).measure;


		@testset "The same" begin
			@test sum(μ_orig1.dist) ≈ 1
			@test sum(μ_orig2.dist) ≈ 1

			# We need to subset the induced measure at nonzero entries, because
			# the regularly computed measures are computed only on the bins that
			# are visited by the orbit, while the induced measure considers all
			# potential bins within the hypercube enclosing the points.
			# Next, we need to sort the points, because (for the same reason) the
			# order in which the bins are considered is different.
			@test all(sort(μ_orig1.dist) .≈ sort(μ_induced1a[μ_induced1a .> 1e-7]))
			@test all(sort(μ_orig1.dist) .≈ sort(μ_induced1b[μ_induced1b .> 1e-7]))
		end





		@testset "Different" begin
			# If induced from different resolution, the distributions will not be
			# identical.
			#sort(μ_orig1.dist) .≈ sort(μ_induced2a[μ_induced2a .> 0])
			@test !all(sort(μ_orig1.dist) .≈ sort(μ_induced2a[μ_induced2a .> 1e-7]))
			@test !all(sort(μ_orig1.dist) .≈ sort(μ_induced2b[μ_induced2b .> 1e-7]))
		end


	end
end
