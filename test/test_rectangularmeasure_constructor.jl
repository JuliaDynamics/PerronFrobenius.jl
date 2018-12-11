import StateSpaceReconstruction
import PerronFrobenius

E = StateSpaceReconstruction.embed([diff(rand(100)) for i = 1:4])

@testset "RectangularInvariantMeasure constructor" begin

	# On raw points
	@test typeof(RectangularInvariantMeasure(E.points, 3)) <: AbstractRectangularInvariantMeasure
	@test typeof(RectangularInvariantMeasure(E.points, 0.4)) <: AbstractRectangularInvariantMeasure
	@test typeof(RectangularInvariantMeasure(E.points, [3, 2, 2, 3])) <: AbstractRectangularInvariantMeasure
	@test typeof(RectangularInvariantMeasure(E.points, [0.6, 0.5, 0.6, 0.5])) <: AbstractRectangularInvariantMeasure

	# On embeddings
	@test typeof(RectangularInvariantMeasure(E, 3)) <: AbstractRectangularInvariantMeasure
	@test typeof(RectangularInvariantMeasure(E, 0.4)) <: AbstractRectangularInvariantMeasure
	@test typeof(RectangularInvariantMeasure(E, [3, 2, 2, 3])) <: AbstractRectangularInvariantMeasure
	@test typeof(RectangularInvariantMeasure(E, [0.6, 0.5, 0.6, 0.5])) <: AbstractRectangularInvariantMeasure
end
