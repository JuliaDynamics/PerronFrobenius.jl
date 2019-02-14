include("potentially_intersecting_simplices.jl")
include("reshape_simplices.jl")
include("function_barriers.jl")

import ..AbstractTransferOperator
import ..TransferOperatorTriangulationApprox

import StateSpaceReconstruction:
    DelaunayTriangulation,
    Simplex,
	point_contained
import StateSpaceReconstruction: subsample_coeffs
import StaticArrays: Size, MMatrix
import Statistics: det
import InplaceOps: @!


#####################################
# Exports
#####################################
export transferoperator_triangulation_approx


function mdet_fm_in_x_view!(M, v, n, starts, stops)
    @inbounds for i = 0:(n - 1)
        M[starts[i+1]:stops[i+1]] = view(v, starts[i+1]:stops[i+1])
    end
    @fastmath det(M)
end


function contained_optim!(signs, s_arr, sx, point, dim, temp_arr, starts, stops)
    # Redefine the temporary simplex. This is in-place, so we don't allocate
    # memory. We could also have re-initialised `signs`, but since we're never
    # comparing more than two consecutive signs, this is not necessary.
    rezero!(s_arr)
    rezero!(signs)
    fill_into!(s_arr, sx)

    #Replace first vertex with the point
    fill_at_inds!(s_arr, point, 1:dim)

    # Signed volume
    signs[1] = sign(mdet_fm_in_x_view!(temp_arr, s_arr, dim + 1, starts, stops))

    rezero!(s_arr)
    fill_into!(s_arr, sx) #reset

    for κ = 2:dim # Check remaining signs and stop if sign changes
        # Replace the ith vertex with the point we're cheking (leaving the
        # 1 appended to Vi intact.)
        idxs = ((dim + 1)*(κ - 1)+1):((dim + 1)*(κ - 1)+ 1 + dim - 1)
        fill_at_inds!(s_arr, point, idxs) # ith change

        signs[κ] = sign(mdet_fm_in_x_view!(temp_arr, s_arr, dim + 1, starts, stops))

        if !(signs[κ-1] == signs[κ])
            return false
        end

        rezero!(s_arr)
        fill_into!(s_arr, sx)
    end

    # Last the last vertex with the point in question
    idxs = ((dim + 1)*(dim)+1):((dim+1)^2-1)
    fill_at_inds!(s_arr, point, idxs)

    signs[end] = sign(mdet_fm_in_x_view!(temp_arr, s_arr, dim + 1, starts, stops))

    if !(signs[end-1] == signs[end])
       return false
    else
        return true
    end
end


function innerloop_optim!(inds::Vector{Int}, signs, s_arr, Sj, pt, dim::Int, M, i::Int, temp_arr, starts, stops)
    for j in 1:length(inds)
        if contained_optim!(signs, s_arr, Sj[j], pt, dim, temp_arr, starts, stops)
            M[inds[j], i] += 1.0
        end
    end
end


"""
    transferoperator_triangulation_approx(invariant_pts,
        DT::DelaunayTriangulation, n_sample_pts::Int = 200, sample_randomly::Bool = false)

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
to = transferoperator_triangulation_approx(invariant_pts, DT)
```

## References
The idea of using a triangulation estimator to compute the transfer operator was
introduced in [1], where it is used to compute physical invariant measures. To compute
a physical invariant measure from a transfer operator `to`, simply call
`invariantmeasure(to)`.

1. Froyland, Gary. "Estimating physical invariant measures and space averages of dynamical systems indicators." Bulletin of the Australian Mathematical Society 56.1 (1997): 157-159.

"""
function transferoperator_triangulation_approx(inv_pts,
        DT::DelaunayTriangulation; n_sample_pts::Int = 200, sample_randomly::Bool = false)

    # Some constants used throughout the function
    n_simplices = length(DT)
    dim = length(DT[1]) - 1

    #=
    # Prepare memory-efficient representations of the simplices
    =#
    simplices = reshape_simplices(inv_pts[1:(end- 1)], DT)
    image_simplices = reshape_simplices(inv_pts[2:end], DT)

    #=
    Compute the convex expansions coefficients needed to generate points
    within a simplex. Also, update number of points if a shape-preserving
    simplex subdivision algorithm was used (in that case, we get a
    *minimum* of `n_sample_pts` set of coefficients.
    =#
    convex_coeffs = subsample_coeffs(dim, n_sample_pts, sample_randomly)
    n_coeffs::Int = size(convex_coeffs, 2)
    convex_coeffs = [convex_coeffs[:, i] for i = 1:n_coeffs]

    # Pre-allocated arrays (SizedArrays, for efficiency)
    pt          = Size(1, dim)(zeros(Float64, 1, dim))
    s_arr       = Size((dim+1)^2)(zeros(Float64, (dim+1)^2))
    signs       = Size(dim + 1)(zeros(Float64, dim + 1))

    # Re-arrange simplices so that look-up is a bit more efficient
    simplex_arrs = Vector{Array{Float64, 2}}(undef, n_simplices)
    imsimplex_arrs = Vector{Array{Float64, 2}}(undef, n_simplices)

    for i in 1:n_simplices
        simplex_arrs[i] = hcat(inv_pts[DT[i]]...,)
        imsimplex_arrs[i] = hcat(inv_pts[DT[i] .+ 1]...,)
    end

    # Pre-allocate some information to avoid memory allocation when
    # computing the fractional overlaps. These are a temporary matrix
    # which we use to compute determinants (to check if a point lies
    # inside a simplex), and the indices needed to access its values
    # at different stages of the checking stage (of whether a point
    # lies inside a simplex). Pre-allocating this information and
    # passing it on reduces runtime and memory allocation
    # by half an order of magnitude.

    temp_arr = rand(MMatrix{dim + 1, dim + 1})
    starts = [((dim + 1)*i)+1 for i = 0:dim]
    stops = starts .+ dim

    # The Markov matrix
    M = zeros(Float64, n_simplices, n_simplices)

    for i in 1:n_simplices
        # Find the indices of the simplices who potentially intersect
        # with the image of simplex #i
        idxs = idxs_potentially_intersecting_simplices(inv_pts, DT, i)
        Sj = Vector{AbstractArray}(undef, length(idxs))
        get_simplices_at_inds!(Sj, idxs, simplices)

        # Generate points within the image of simplex #i and check
        # what fraction of those points fall into the simplices.
        # We approximate the intersecting volume between the
        # image simplex and the individual simplices by those
        # fractions.
        @views is = imsimplex_arrs[i]
        for k in 1:n_coeffs
            @! pt = transpose(convex_coeffs[k]) * transpose(is)
            innerloop_optim!(idxs, signs, s_arr, Sj, pt, dim, M, i, temp_arr, starts, stops)
        end
    end
    # Need to normalise, because all we have up until now is counts
    # of how many points inside the image simplex falls into
    # the simplices.
    return TransferOperatorTriangulationApprox(transpose(M) ./ n_coeffs)
end


"""
    transferoperator_triangulation_approx(invariant_pts; n_sample_pts::Int = 200,
        sample_randomly::Bool = false) -> TransferOperatorTriangulationApprox

Compute an approximation of the transfer operator over a triangulated partition of the
provided points (`invariant_pts`). The estimator assumes that the points are invariant
under the forward projection of the simplices in the triangulation, which means that
the last point must be contained within the convex hull of all previous points.

*Tip: The `invariantize` function from `StateSpaceReconstruction` can be used to move the
last point toward the center of the embedding until invariance is achieved.*

This function computes approximate simplex intersections and is slow for very large
triangulations. It should only be used on small to moderately sized datasets (i.e.
in the order of hundreds of points). For larger datasets, use a rectangular estimator.

Returns a  instance.

## Example
```julia
# Create some random points and make sure they are invariant
pts = [rand(3) for i = 1:30]
invariant_pts = invariantize(pts)

# Compute the transfer operator
to = transferoperator_triangulation_approx(invariant_pts)
```

## References
The idea of using a triangulation estimator to compute the transfer operator was
introduced in [1], where it is used to compute physical invariant measures. To compute
a physical invariant measure from a transfer operator `to`, simply call
`invariantmeasure(to)`.

1. Froyland, Gary. "Estimating physical invariant measures and space averages of dynamical systems indicators." Bulletin of the Australian Mathematical Society 56.1 (1997): 157-159.
"""
function transferoperator_triangulation_approx(invariant_pts;
        n_sample_pts::Int = 200, sample_randomly::Bool = false)
    # Triangulate all points but the last (we need the last point
    # for the forward projection of the simplices)
    DT = DelaunayTriangulation(invariant_pts[1:(end - 1)])
    transferoperator_triangulation_approx(invariant_pts, DT)
end
