import StaticArrays:
    Size

import StateSpaceReconstruction:
    DelaunayTriangulation

"""
    reshape_simplices(pts::Vector{Vector{T}}, DT::DelaunayTriangulation) where T

Creates alternative representations of the simplices that allow for efficient
(mostly non-allocating) checking if a point lies inside a simplex.

This function adds some information needed when calculating the orientations
of the simplices *before* we start the computation (appends a row of ones on top of the
matrix representation of the simplex). Otherwise, the fallback is machinery in Simplices.jl,
which is much slower because the additional information must be added
`(n_simplices^2)*n_sample_pts` times (creating that many new matrices), which takes forever.

Reshaping the simplices and appending that information beforehand avoids excessive
memory allocation.
"""
function reshape_simplices(pts, DT::DelaunayTriangulation) where T
    n_simplices = length(DT.indices)
    n_vertices = length(DT[1])
    dim = n_vertices - 1

    # The [:, j, i] th entry of these two arrays holds the jth vertex of the
    # ith simplex, but instead of having just `dim` vertices, we append a `1`
    # to the end of the vectors. This allows for efficent (non-allocating)
    # computation within the `contains_point_lessalloc!` function. If we instead
    # would have appended the 1's inside that function, we would be performing
    # memory-allocating operations, which are very expensive. Doing this instead
    # gives orders of magnitude speed-ups for sufficiently large triangulations.
    S1 = Array{Float64}(undef, dim + 1, dim + 1, n_simplices)

    # Collect simplices in the form of (dim+1)^2-length column vectors. This
    # also helps with the
    simplices = Size((dim+1)^2, n_simplices)(zeros((dim+1)^2, n_simplices))

    @inbounds for i in 1:n_simplices
        for j in 1:n_vertices
            S1[:, j, i] = vcat(pts[DT[i][j]], 1.0)
        end

        simplices[:, i] = reshape(S1[:, :, i], (dim+1)^2)
    end

    return simplices
end
