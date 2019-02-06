
import StateSpaceReconstruction:
    invariantize,
    triangulate

import ..TransferOperators:
    transferoperator_triangulation_exact

"""
    TriangulationExactInvariantMeasure

An invariant measure estimated from a triangulation using the exact simplex 
intersection estimator `transferoperator_triangulation_exact` for the 
approximation of the transfer matrix.

## Fields
- **`points`**: The points from which the measure is estimated.
- **`invariantized_points`**: The points from which the measure is estimated, 
possibly with the last point moved toward the center of the convex hull of 
the points to ensure no information is lost when computing the transfer matrix.
- **`triangulation`**. A triangulation of `invariantized_points[1:(end - 1)]` which is used to estimate the transfer matrix.

"""
struct TriangulationExactInvariantMeasure{PT, TT, TOT} <: AbstractTriangulationInvariantMeasure
    points::PT
    invariantized_points::PT
    triangulation::TT
    transferoperator::TOT
    measure::InvariantDistribution
end 


function TriangulationExactInvariantMeasure(points::PT; kwargs...) where {PT}
    inv_points = invariantize(points)
    triang = triangulate(inv_points[1:(end - 1)])
    TO = transferoperator_triangulation_exact(points, triang.simplexindices)
    measure = invariantmeasure(TO, kwargs...)
    
    TriangulationInvariantMeasureExact(
        points, 
        inv_points, 
        triang, 
        TO, 
        measure)
end


export TriangulationExactInvariantMeasure