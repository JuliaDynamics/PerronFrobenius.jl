"""
	TransferOperatorRectangularBinning

A transfer operator estimated over rectangular binning of a set of points, using the
`TransferOperatorEstimatorRectangularBinning` estimator.

The return type for the `TransferOperatorEstimatorRectangularBinning` estimator.

## Fields
-  **`transfermatrix`**. The estimated transfer matrix.

"""
struct TransferOperatorRectangularBinning <: AbstractTransferOperator
	transfermatrix::AbstractArray{Float64, 2}
end

export TransferOperatorRectangularBinning
