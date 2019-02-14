"""
    TransferOperatorTriangulationExact

A transfer operator estimated over a triangulation of a set of points. Transition
probabilities are computed from exact simplex intersections.

The return type for the `TransferOperatorEstimatorTriangulationExact` estimator.

## Fields
- **`transfermatrix`**. The estimated transfer matrix.

"""
struct TransferOperatorTriangulationExact <: AbstractTransferOperator
    transfermatrix::AbstractArray{Float64, 2}
end

export TransferOperatorTriangulationExact
