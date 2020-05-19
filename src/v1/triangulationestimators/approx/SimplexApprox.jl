export SimplexApprox
import StaticArrays: SizedVector, SizedMatrix, MMatrix


include("helper_functions.jl")

"""
    SimplexApprox

A transfer operator estimator using a triangulation partition and approximate simplex intersections [^Diego2019]. 

[^Diego2019]: Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""
struct SimplexApprox <: TransferOperator
    bc::String
    
    function SimplexApprox(bc::String = "circular")
        isboundarycondition(bc, "triangulation")  || error("Boundary condition '$bc' not valid.")
        new(bc)
    end
end
Base.show(io::IO, se::SimplexApprox) = print(io, "SimplexApprox{$(se.bc)}")


""" Generate a TransferOperatorGenerator for an exact simplex estimator."""
function transferoperatorgenerator(pts, method::SimplexApprox)
    # modified points, where the image of each point is guaranteed to lie within the convex hull of the previous points
    invariant_pts = invariantize(pts)
    
    # triangulation of the invariant points. The last point is excluded, so that 
    # the last vertex also can be mapped forward one step in time.
    triang = DelaunayTriangulation(invariant_pts[1:end-1])

    init = (invariant_pts = invariant_pts, triang = triang)
    
    TransferOperatorGenerator(method, pts, init)
end


"""
n: The number of points sampled inside each simplex.
randomsampling: Sample randomly (`random == true`), or using an even simplex splitting routine (`random == false`)
"""
function (tog::TransferOperatorGenerator{<:SimplexApprox})(; tol = 1e-8, randomsampling::Bool = false, n::Int = 100)
    pts, inds = getfield.(Ref(tog.init), (:invariant_pts, :triang))
    D = length(pts[1])
    N = length(inds)
    ϵ = tol / N
    
    # Prepare memory-efficient representations of the simplices
    simplices = reshape_simplices(pts[1:(end- 1)], inds)
    image_simplices = reshape_simplices(pts[2:end], inds)
    
    #=
    Compute the convex expansions coefficients needed to generate points
    within a simplex. Also, update number of points if a shape-preserving
    simplex subdivision algorithm was used (in that case, we get a
    *minimum* of `n_sample_pts` set of coefficients.
    =#
    convex_coeffs = SimplexTypes.subsample_coeffs(D, n, randomsampling)
    n_coeffs::Int = size(convex_coeffs, 2)
    convex_coeffs = [convex_coeffs[:, i] for i = 1:n_coeffs]
    
    # Pre-allocated arrays (SizedArrays, for efficiency)
    L = (D+1)^2
    NVERTICES = D+1
    pt          = SizedMatrix{1, D}(zeros(1, D))
    s_arr       = SizedVector{L}(zeros(L))
    signs       = SizedVector{NVERTICES}(zeros(NVERTICES))
    
    # Re-arrange simplices so that look-up is a bit more efficient
    simplex_arrs = Vector{Array{Float64, 2}}(undef, N)
    imsimplex_arrs = Vector{Array{Float64, 2}}(undef, N)
    for i in 1:N
        simplex_arrs[i] = hcat(pts[inds[i]]...,)
        imsimplex_arrs[i] = hcat(pts[inds[i] .+ 1]...,)
    end
    
    # Pre-allocate some information to avoid memory allocation when
    # computing the fractional overlaps. These are a temporary matrix
    # which we use to compute determinants (to check if a point lies
    # inside a simplex), and the indices needed to access its values
    # at different stages of the checking stage (of whether a point
    # lies inside a simplex). Pre-allocating this information and
    # passing it on reduces runtime and memory allocation
    # by half an order of magnitude.
    temp_arr = rand(MMatrix{D + 1, D + 1})
    starts = [((D + 1)*i)+1 for i = 0:D]
    stops = starts .+ D
    
    
    # The Markov matrix
    M = zeros(Float64, N, N)

    for i in 1:N
        # Find the indices of the simplices who potentially intersect
        # with the image of simplex #i
        idxs = idxs_potentially_intersecting_simplices(pts, inds, i)
        Sj = Vector{AbstractArray}(undef, length(idxs))
        get_simplices_at_inds!(Sj, idxs, simplices)

        # Generate points within the image of simplex #i and check
        # what fraction of those points fall into the simplices.
        # We approximate the intersecting volume between the
        # image simplex and the individual simplices by those
        # fractions.
        @views is = imsimplex_arrs[i]
        for k in 1:n_coeffs
            pt = transpose(convex_coeffs[k]) * transpose(is)
            innerloop_optim!(idxs, signs, s_arr, Sj, pt, D, M, i, temp_arr, starts, stops)
        end
    end
    # Need to normalise, because all we have up until now is counts
    # of how many points inside the image simplex falls into
    # the simplices.
    return transpose(M) ./ n_coeffs

end