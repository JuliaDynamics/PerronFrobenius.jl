import StaticArrays:
    MVector
import ..AbstractTransferOperator
import ..TransferOperatorTriangulationExact

import StateSpaceReconstruction:
    MutableSimplex,
    orientation



export transferoperator_triangulation_exact

"""
    transferoperator_triangulation_exact(invariant_pts) -> TransferOperatorTriangulationExact

Compute an approximation of the transfer operator over a triangulated partition of the
provided points (`invariant_pts`). The estimator assumes that the points are invariant
under the forward projection of the simplices in the triangulation, which means that
the last point must be contained within the convex hull of all previous points.

*Tip: The `invariantize` function from `StateSpaceReconstruction` can be used to move the
last point toward the center of the embedding until invariance is achieved.*

This function computes *exact* simplex intersections and is thus *slow*. It should
only be used on small datasets.


## Example

```julia
# Create some random points and make sure they are invariant
pts = [rand(3) for i = 1:30]
invariant_pts = invariantize(pts)

# Compute the transfer operator
to = transferoperator_triangulation_exact(invariant_pts)
```

## References
The idea of using a triangulation estimator to compute the transfer operator was
introduced in [1], where it is used to compute physical invariant measures. To compute
a physical invariant measure from a transfer operator `to`, simply call
`invariantmeasure(to)`.

1. Froyland, Gary. "Estimating physical invariant measures and space averages of dynamical systems indicators." Bulletin of the Australian Mathematical Society 56.1 (1997): 157-159.
"""
function transferoperator_triangulation_exact(invariant_pts)
    DT = DelaunayTriangulation(invariant_pts[1:end-1])
    transferoperator_triangulation_exact(invariant_pts, DT)
end

"""
    transferoperator_triangulation_exact(invariant_pts, DT::DelaunayTriangulation)

Approximate the transfer operator over a triangulation of `invariant_pts`. Does the same
as the single-argument version of the function, but requires that the triangulation
is precomputed on `invariant_pts[1:(end - 1)]` (excluding the last point).

## Example
```julia
# Create some random points and make sure they are invariant
pts = [rand(3) for i = 1:30]
invariant_pts = invariantize(pts)

# Triangulate all but the last point (exclude the last, so that an image simplex
# always exists)
DT = DelaunayTriangulation(invariant_pts[1:(end - 1)])

# Compute the transfer operator
to = transferoperator_triangulation_exact(invariant_pts, DT)
```

## References
The idea of using a triangulation estimator to compute the transfer operator was
introduced in [1], where it is used to compute physical invariant measures. To compute
a physical invariant measure from a transfer operator `to`, simply call
`invariantmeasure(to)`.

1. Froyland, Gary. "Estimating physical invariant measures and space averages of dynamical systems indicators." Bulletin of the Australian Mathematical Society 56.1 (1997): 157-159.
"""
function transferoperator_triangulation_exact(invariant_pts,
        DT::DelaunayTriangulation)

    dim = length(invariant_pts[1])
    n_simplices = length(DT)

    #=
    Tolerance for simplex size when filling the transfer matrix
    =#
    ϵ::Float64 = 1e-8 / n_simplices

    # Pre-allocate the simplices
    image_simplex = MutableSimplex(zeros(Float64, dim, dim + 1))
    simplex = MutableSimplex(zeros(Float64, dim, dim + 1))

    # Pre-allocate transfer matrix
    TO = zeros(Float64, n_simplices, n_simplices)

    for j in 1:n_simplices
        # Fill the image simplex
        for k = 1:(dim + 1)
            image_simplex[k] .= 0.0
            image_simplex[k] = invariant_pts[DT[j][k] + 1]
        end

        imvol = abs(orientation(image_simplex))

        for i = 1:n_simplices
            # Fill which we're testing if arrives in the image simplex,
            # and if it does, how much of it overlaps with the image
            # simplex.
            for k = 1:(dim + 1)
                simplex[k] .= 0.0
                simplex[k] = invariant_pts[DT[i][k]]
            end

            # Only compute the entry of the transfer matrix
            # if simplices are of sufficient size.
            vol = abs(orientation(simplex))

            if vol * imvol > 0 && (vol/imvol) > ϵ
                TO[j, i] = (simplex ∩ image_simplex) / imvol
            end
        end
    end
    TransferOperatorTriangulationExact(TO)
end
