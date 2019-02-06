@testset "Induced measure constructor" begin
	pts = rand(3,1000)
	ϵF, ϵs = [4, 4, 4], [ [3, 5, 2], [5, 5, 5]]
	iv1 = inducedrectangularinvariantmeasure(pts, ϵF, ϵs[1])
	iv2 = inducedrectangularinvariantmeasure(pts, ϵF, ϵs[2])
	@test iv1 isa InducedRectangularInvariantMeasure
	@test iv2 isa InducedRectangularInvariantMeasure
	@test sum(iv1.measure_induced.dist) ≈ 1
	@test sum(iv2.measure_induced.dist) ≈ 1
end
