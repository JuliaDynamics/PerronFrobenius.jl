using StateSpaceReconstruction
using PerronFrobenius
using Test 
E = cembed([diff(rand(100)) for i = 1:4])

@testset "RectangularInvariantMeasure constructor" begin

	# On raw points
	@test typeof(rectangularinvariantmeasure(E.points, RectangularBinning(3))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E.points, RectangularBinning(0.4))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E.points, RectangularBinning([3, 2, 2, 3]))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E.points, RectangularBinning([0.6, 0.5, 0.6, 0.5]))) <: AbstractRectangularInvariantMeasure

	# On embeddings
	@test typeof(rectangularinvariantmeasure(E, RectangularBinning(3))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E, RectangularBinning(0.4))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E, RectangularBinning([3, 2, 2, 3]))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E, RectangularBinning([0.6, 0.5, 0.6, 0.5]))) <: AbstractRectangularInvariantMeasure
end
