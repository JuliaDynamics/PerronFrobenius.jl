using StateSpaceReconstruction
using PerronFrobenius
using Test 
using DelayEmbeddings
using StaticArrays
Vpts = [rand(4) for i = 1:50]
Dpts = Dataset([rand(4) for i = 1:50])
Spts = Dataset([rand(4) for i = 1:50])
Mpts = Dataset([rand(4) for i = 1:50])
pts = rand(50, 4)
Tpts = rand(4, 50)

@testset "RectangularInvariantMeasure, constructor" begin

	# On raw points
	@test typeof(rectangularinvariantmeasure(pts, RectangularBinning(3))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(pts, RectangularBinning(0.4))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(pts, RectangularBinning([3, 2, 2, 3]))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(pts, RectangularBinning([0.6, 0.5, 0.6, 0.5]))) <: AbstractRectangularInvariantMeasure

	# On embeddings
	@test typeof(rectangularinvariantmeasure(pts, RectangularBinning(3))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(pts, RectangularBinning(0.4))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(pts, RectangularBinning([3, 2, 2, 3]))) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(pts, RectangularBinning([0.6, 0.5, 0.6, 0.5]))) <: AbstractRectangularInvariantMeasure
end

@testset "RectangularInvariantMeasure, different inputs" begin 
	b = RectangularBinning(3)
	@test typeof(rectangularinvariantmeasure(Vpts, b)) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(Dpts, b)) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(Spts, b)) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(Mpts, b)) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(pts, b)) <: AbstractRectangularInvariantMeasure
	@test typeof(rectangularinvariantmeasure(Tpts, b)) <: AbstractRectangularInvariantMeasure
end
