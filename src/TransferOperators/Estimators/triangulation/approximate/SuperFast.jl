
import StateSpaceReconstruction:
    DelaunayTriangulation,
    Simplex,
    centroid,
    radius,
    MutableSSimplex

import StaticArrays:
    MMatrix,
    MVector

import Statistics:
    transpose!

import InplaceOps:
    @!

export SuperFast



"""
    simplices_within_circumsphere!(inds::Vector{Int}, idx::Int,
        simplices, image_simplices,
        radii_simplices, radii_image_simplices,
        centroids_simplices, centroids_image_simplices)

Find the indices of the simplices whose centroid lies within the circumsphere of the image
of simplex #`idx`. Store the indices in the pre-allocated vector `inds` (the latter part
of the vector might contain zeros if some simplices do not lie within the circumsphere,
so keep track of this outside the function).
"""
function simplices_within_circumsphere!(inds::Vector{Int}, N, idx::Int,
        simplices, image_simplices,
        radii_simplices, radii_image_simplices,
        centroids_simplices, centroids_image_simplices)

    # Because we're re-using `inds`, re-initialize it with zeros
    inds .= 0

    # We've found zero simplices within the circumsphere to begin with.
    N[] = 0

	n_simplices = length(simplices)

	@inbounds for i = 1:n_simplices
        x = transpose(centroids_image_simplices[idx] .- centroids_simplices[i]) * (centroids_image_simplices[idx] .- centroids_simplices[i])
	    δ = x .- ((radii_image_simplices[idx] + radii_simplices[i])^2)

	    if δ[1] < 0
            N[] += 1
            inds[N[]] = i
	    end
	end

end


"""
Get all simplices in the triangulation.
"""
function get_simplices(pts, DT::DelaunayTriangulation)
    dim = length(DT[1]) - 1

    [Simplex(pts[DT[i]]) for i = 1:length(DT)]
end


function mdet_fm_in_x_view!(M, v, n, startstops)
    @inbounds for i = 0:(n - 1)
        M[startstops[i+1]] = view(v, startstops[i+1])
    end
    @fastmath det(M)
end


"""
Fill the 2nd to last row in the i-th column of `prealloc_matrix` with the
2nd to last row in the `i`-th column of `toprow_simplex` (replacing the `i`th vertex).
"""
function fill_vertex(prealloc_matrix::Array{Float64, 2}, pt,
        startstops::Vector{UnitRange{Int64}}, i::Int)
    @show pt, startstops[i]
    prealloc_matrix[startstops[i]] .= pt
end

"""
Fill the 2nd to last row in the i-th column of `prealloc_matrix` with the
2nd to last row in the `i`-th column of `toprow_simplex` (replacing the `i`th vertex).
"""
function refill_vertex(prealloc_matrix::Array{Float64, 2}, toprow_simplex::Array{Float64, 2},
        startstops::Vector{UnitRange{Int64}}, i::Int)

    prealloc_matrix[startstops[i]] .= toprow_simplex[startstops[i]]
end


"""

## Arguments
- **`pt::Vector{Float64}`**: The point
- **`toprow_simplex::Array{Float64, 2}`**: A matrix representation of the simplex where each
    column represents a vertex, but with a row of ones appended on top of it, so that the
    total dimension is `(dim + 1)-by-(dim + 1)`. This is needed when computing determinants.
- **`dim`**: The dimension of the space.
"""
function iscontained_superoptim2!(pt, toprow_simplex, prealloc_matrix, dim, signs, startstops, found::Array{Float64,0})

    # Replace first vertex with the point
    fill_vertex(prealloc_matrix, pt, startstops, 1)

    signs[1] = sign(mdet_fm_in_x_view!(prealloc_matrix, toprow_simplex, dim + 1, startstops))

    refill_vertex(prealloc_matrix, toprow_simplex, startstops, 1)

    for k = 2:dim
        idxs = ((dim + 1)*(k - 1)+1):((dim + 1)*(k - 1)+ 1 + dim - 1)
        fill_vertex(prealloc_matrix, toprow_simplex, startstops, 1)
    end
    found[] = true
end



function SuperFast(inv_pts; n_sample_pts::Int = 200, sample_randomly::Bool = false)

    # Triangulate all points but the last (we need the last point
    # for the forward projection of the simplices)
    DT = DelaunayTriangulation(inv_pts[1:(end - 1)])

    # Some constants used throughout the function
    n_simplices = length(DT)
    dim = length(DT[1]) - 1

    #=
    # Prepare memory-efficient representations of the simplices
    =#
    simplices = get_simplices(inv_pts[1:(end- 1)], DT)
    image_simplices = get_simplices(inv_pts[2:end], DT)

    radii_simplices = [radius(simplex) for simplex in simplices]
    radii_image_simplices = [radius(simplex) for simplex in image_simplices]
    centroids_simplices = [centroid(simplex) for simplex in simplices]
    centroids_image_simplices = [centroid(simplex) for simplex in image_simplices]

    # Matrix representations of the simplices
    matrix_image_simplices = [hcat(image_simplices[i]...,) for i = 1:n_simplices]


    #=
    Compute the convex expansions coefficients needed to generate points
    within a simplex. Also, update number of points if a shape-preserving
    simplex subdivision algorithm was used (in that case, we get a
    *minimum* of `n_sample_pts` set of coefficients.
    =#
    convex_coeffs = subsample_coeffs(dim, n_sample_pts, sample_randomly)
    n_coeffs::Int = size(convex_coeffs, 2)
    convex_coeffs = [reshape(convex_coeffs[:, i], dim + 1, 1) for i = 1:n_coeffs]


    # The Markov matrix
    M = zeros(Float64, n_simplices, n_simplices)

    inds = zeros(Int, n_simplices)
    n_simplices_within_circumsphere = Array{Int, 0}(undef)



    # Indices to fill when computing whether a point lies inside a simplex or not.
    # The top row will always be a row of ones, so we need to know which other
    # indices to fill.
    starts = [((dim + 1)*i)+2 for i = 0:dim]
    stops = starts .+ (dim - 1)
    startstops = [starts[i]:stops[i] for i = 1:length(starts)]

    # Matrix to re-use when computing whether a point lies inside a simplex.
    prealloc_matrix = zeros(Float64, dim + 1, dim + 1)#zeros(MMatrix{dim + 1, dim + 1})
    toprow_inds = [((dim + 1)*i) + 1 for i = 0:dim]
    [prealloc_matrix[i] = 1.0 for i in toprow_inds]

    # Signs vector when checking convex coefficients
    signs = zeros(Int, dim + 1)


    # Matrix representations of the simplices
    toprow_simplices = [[transpose(ones(dim + 1)); hcat(simplices[i]...,)] for i = 1:n_simplices]

    pt = reshape(zeros(Float64, dim), dim, 1)

    transfermatrix = zeros(Float64, n_simplices, n_simplices)

    found = Array{Float64, 0}(undef)

    for im_idx in 1:n_simplices
        # Find the indices of the simplices who potentially intersect
        # with the image of simplex #i
        simplices_within_circumsphere!(inds, n_simplices_within_circumsphere, im_idx,
               simplices, image_simplices,
               radii_simplices, radii_image_simplices,
               centroids_simplices, centroids_image_simplices)


        # Generate points within the image of simplex #i and check
        # what fraction of those points fall into the simplices.
        # We approximate the intersecting volume between the
        # image simplex and the individual simplices by those
        # fractions.

        # View into the image simplex at index `im_idx`
        im_simplex = view(matrix_image_simplices[im_idx], :, :)

        inds_potentially_intersecting = [inds[i] for i = 1:n_simplices_within_circumsphere[]]

        @show inds_potentially_intersecting
        for k in inds_potentially_intersecting
            # Should be a (dim-by(dim+1) * (dim + 1) multiplication -> (dim, 1) vector
            @! pt = im_simplex * convex_coeffs[k]

            # Keep track of whether the simplex is found
            found[] = false

            # Now, loop over the simplices and see if the point we generated inside the
            # image simplex lies inside any of them, if so, update the counts.
            for idx in n_simplices_within_circumsphere
                iscontained_superoptim2!(pt, toprow_simplices[idx], prealloc_matrix, dim, signs, startstops, found)

                #if found[] == true
                #    break
                #end
            end
        end
    end
    # Need to normalise, because all we have up until now is counts
    # of how many points inside the image simplex falls into
    # the simplices.
    #return transpose(transfermatrix) ./ n_coeffs
end
