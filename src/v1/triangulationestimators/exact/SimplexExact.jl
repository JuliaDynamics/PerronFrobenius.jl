export SimplexExact, transferoperator

"""
    SimplexExact() <: TransferOperator

A transfer operator estimator using a triangulation partition and exact 
simplex intersections[^Diego2019]. 

*Note: due to the exact simplex intersection computations, this estimator is slow compared to [`SimplexApprox`](@ref).*

[^Diego2019]: Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
"""
struct SimplexExact <: TransferOperator
    bc::String
    
    function SimplexExact(bc::String = "circular")
        isboundarycondition(bc, "triangulation")  || error("Boundary condition '$bc' not valid.")
        new(bc)
    end
end
Base.show(io::IO, se::SimplexExact) = print(io, "SimplexExact{$(se.bc)}")


""" Generate a transopergenerator for an exact simplex estimator."""
function transopergenerator(pts, method::SimplexExact)
    # modified points, where the image of each point is guaranteed to lie within the convex hull of the previous points
    invariant_pts = invariantize(pts)
        
    # triangulation of the invariant points. The last point is excluded, so that 
    # the last vertex also can be mapped forward one step in time.
    triang = DelaunayTriangulation(invariant_pts[1:end-1])

    init = (invariant_pts = invariant_pts, triang = triang)

    TransferOperatorGenerator(method, pts, init)
end

function (tog::TransferOperatorGenerator{<:SimplexExact})(; tol = 1e-8)
    invariant_pts, triang = getfield.(Ref(tog.init), (:invariant_pts, :triang))
    
    D = length(invariant_pts[1])
    N = length(triang)
    ϵ = tol / N
    
    # Pre-allocate simplex and its image
    image_simplex = SimplexTypes.MutableSimplex(zeros(Float64, D, D+1))
    simplex = SimplexTypes.MutableSimplex(zeros(Float64, D, D+1))
    M = zeros(Float64, N, N)
    
    for j in 1:N
        for k = 1:(D + 1)
            # Get the vertices of the image simplex
            image_simplex[k] .= invariant_pts[triang[j][k] + 1]
        end
        
        imvol = abs(SimplexTypes.orientation(image_simplex))

        for i = 1:N
            for k = 1:(D + 1)
                # Get the vertices of the source simplex
                simplex[k] .= invariant_pts[triang[i][k]]
            end

            # Only compute the entry of the transfer matrix
            # if simplices are of sufficient size.
            vol = abs(SimplexTypes.orientation(simplex))

            if vol * imvol > 0 && (vol/imvol) > ϵ
                M[j, i] = SimplexTypes.intersect(simplex, image_simplex) / imvol
            end
        end
    end
    params = (; tol = tol)
    TransferOperatorApproximation(tog, M, params)
end

function transferoperator(pts, method::SimplexExact; tol::Real = 1e-8)
    tog = transopergenerator(pts, method)
    tog(tol = tol)
end