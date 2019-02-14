import CausalityToolsBase: 
    TriangulationBinningScheme,
    TriangulationBinning,
    SimplexIntersectionType, 
    ExactIntersection

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


function triangulationexactinvariantmeasure(points::PT; kwargs...) where {PT}
    inv_points = invariantize(points)
    triang = triangulate(inv_points[1:(end - 1)])
    TO = transferoperator_triangulation_exact(points, triang.simplexindices)
    measure = invariantmeasure(TO, kwargs...)
    
    TriangulationExactInvariantMeasure(
        points, 
        inv_points, 
        triang, 
        TO, 
        measure)
end

function description(μ::TriangulationExactInvariantMeasure)
    T_pts = typeof(μ.points)
    T_invpts = typeof(μ.invariantized_points)
    triang = μ.triangulation
    to = μ.transferoperator 
    measure = μ.measure
    
    npts = 

    invmeasure_type = typeof(μ)
    return join(["TriangulationExactInvariantMeasure", "\n",
                "  μ.points:\t\t", T_pts,   "\n", 
                 "  μ.invariantized_points:\t", T_invpts, "\n",
                 "  μ.triangulation:\t", triang,
                "  μ.transferoperator:\t", to, 
                "  μ.measure:\t\t", measure])
end

Base.show(io::IO, μ::TriangulationExactInvariantMeasure) = println(io, description(μ))


export TriangulationExactInvariantMeasure, triangulationexactinvariantmeasure