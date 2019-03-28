import PerronFrobenius

pts = rand(3, 10000)
ϵF = [0.4, 0.4, 0.4]
ϵs = [ϵF, ϵF]

induced_measures = [inducedrectangularinvariantmeasure(pts, ϵF, ϵᵢ) for ϵᵢ in [ϵF, ϵF]]
avg_measure = average_measure(induced_measures)

measure_ϵF = invariantmeasure(pts, RectangularBinning(ϵF))

@testset "Average measure constructor" begin
    ϵF = [0.4, 0.4, 0.4]
    ϵs = [ϵF, ϵF]
    pts = rand(3, 5000)

    ϵF = RectangularBinning(ϵF)
    ϵs = [RectangularBinning(ϵi) for (i, ϵi) in enumerate(ϵs)]

    @testset "Computing average measure directly from points" begin
        # Construct a AverageRectangularInvariantMeasure directly from data, final
        # resolution and resolutions from which we induce the measure.
        a2 = averagerectangularinvariantmeasure(pts, ϵF, ϵs)
        @test a2 isa AverageRectangularInvariantMeasure
    end
end

# PerronFrobenius.AverageRectangularInvariantMeasure(pts, ϵF,ϵs)
