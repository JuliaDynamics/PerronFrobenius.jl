using Test
using PerronFrobenius



pts = rand(3, 10000)
μF = [4, 4, 4]

μ_avg = PerronFrobenius.average_measure(pts, μF , [μF, μF])
μ_induced = PerronFrobenius.induced_measure(pts, μF, μF)

@testset "Average measure" begin
	@testset "Average measure sums to 1" begin
		@test sum(μ_avg) ≈ 1
	end

	@testset "Measure induced from multiple same resolutions matches original" begin
		μ_induced = induced_measure(pts, μF, μF)
		@test all(μ_avg .≈ μ_induced)
	end
end
