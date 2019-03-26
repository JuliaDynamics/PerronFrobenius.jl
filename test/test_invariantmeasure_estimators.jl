

pts = [rand(3) for i = 1:1000]

@test invariantmeasure(pts, RectangularBinning(3)) isa RectangularInvariantMeasure 

#@test invariantmeasure(pts, InducedRectangularBinning(3, 4)) isa InducedRectangularInvariantMeasure 
#@test invariantmeasure(pts, RectangularBinning(3), RectangularBinning(4)) isa InducedRectangularInvariantMeasure 

#@test invariantmeasure(pts, AverageRectangularBinning(3, [4, 5])) isa AverageRectangularInvariantMeasure 
#@test invariantmeasure(pts, RectangularBinning(3), [RectangularBinning(4), RectangularBinning(5)]) isa AverageRectangularInvariantMeasure 
