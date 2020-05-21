const SGA = TransferOperatorGenerator{SingleGrid{RectangularBinning}}
@recipe function f(iv::InvariantMeasure{TransferOperatorApproximation{<:SGA}})
    heatmap(iv.to.M)
end