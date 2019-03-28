using StateSpaceReconstruction
using PerronFrobenius
using Test 
E = cembed([diff(rand(100)) for i = 1:4])

@testset "RectangularInvariantMeasure constructor" begin

	# On raw points
	@test typeof(rectangularinvariantmeasure(E.points, 3)) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E.points, 0.4)) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E.points, [3, 2, 2, 3])) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E.points, [0.6, 0.5, 0.6, 0.5])) <: AbstractRectangularInvariantMeasure

	# On embeddings
	@test typeof(rectangularinvariantmeasure(E, 3)) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E, 0.4)) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E, [3, 2, 2, 3])) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(E, [0.6, 0.5, 0.6, 0.5])) <: AbstractRectangularInvariantMeasure
end
