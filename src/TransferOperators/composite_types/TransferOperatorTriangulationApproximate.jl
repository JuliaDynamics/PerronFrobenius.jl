"""
    TransferOperatorTriangulationApprox

A transfer operator estimated over a triangulation of a set of points. Transition
probabilities are computed from approximate simplex intersections.

The return type for the `transferoperator_triangulation_approx` estimator.

## Fields
-  **`transfermatrix`**. The estimated transfer matrix.

"""
struct TransferOperatorTriangulationApprox <: AbstractTransferOperator
    transfermatrix::AbstractArray{Float64, 2}
end


export TransferOperatorTriangulationApprox
