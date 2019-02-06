import ..TransferOperators:
    transferoperator_triangulation_approx

import StateSpaceReconstruction:
    invariantize,
    triangulate

export TriangulationApproxInvariantMeasure, triangulationapproxinvariantmeasure

"""
TriangulationApproxInvariantMeasure

An invariant measure estimated from a triangulation using the approximate simplex 
intersection estimator `transferoperator_triangulation_approx` for the 
approximation of the transfer matrix.

## Fields
- **`points`**: The points from which the measure is estimated.
- **`invariantized_points`**: The points from which the measure is estimated, 
possibly with the last point moved toward the center of the convex hull of 
the points to ensure no information is lost when computing the transfer matrix.
- **`triangulation`**. A triangulation of `invariantized_points[1:(end - 1)]` which is used to estimate the transfer matrix.

"""
struct TriangulationApproxInvariantMeasure{PT, TT, TOT} <: AbstractTriangulationInvariantMeasure
    points::PT
    invariantized_points::PT
    triangulation::TT
    transferoperator::TOT
    measure::InvariantDistribution
end 



function triangulationapproxinvariantmeasure(points::Vector{Vector{T}}; n_sample_pts = 200, 
    sample_randomly::Bool = false, kwargs...) where {T}
    
    inv_points = invariantize(points)
    triang = triangulate(inv_points[1:(end - 1)])

    TO = transferoperator_triangulation_approx(points, 
        triang.simplexindices, 
        n_sample_pts = n_sample_pts, sample_randomly = sample_randomly)

    measure = invariantmeasure(TO, kwargs...)

    TriangulationInvariantMeasureApprox(
        points, 
        inv_points, 
        triang, 
        TO, 
        measure)
end

