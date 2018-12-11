import PerronFrobenius

pts = rand(3, 10000)
ϵF = [0.4, 0.4, 0.4]
ϵs = [ϵF, ϵF]
induced_measures = [InducedRectangularInvariantMeasure(pts, ϵF, ϵᵢ) for ϵᵢ in [ϵF, ϵF]]
avg_measure = PerronFrobenius.average_measure(induced_measures)

@testset "Average measure constructor" begin

    # Construct a AverageRectangularInvariantMeasure for precomputed induced
    # measures and average measure
    @testset "Precomputed induced and average measures" begin
        a1 = AverageRectangularInvariantMeasure(pts, ϵF, ϵs, induced_measures, avg_measure)
        @test a1 isa AverageRectangularInvariantMeasure
    end

    @testset "Computing average measure directly from points" begin
        # Construct a AverageRectangularInvariantMeasure directly from data, final
        # resolution and resolutions from which we induce the measure.
        a2 = AverageRectangularInvariantMeasure(pts, ϵF, ϵs)
        @test a2 isa AverageRectangularInvariantMeasure
    end
end

# PerronFrobenius.AverageRectangularInvariantMeasure(pts, ϵF,ϵs)
