export DelaunayTriangulation

import StaticArrays: SVector, MVector 
import Simplices.Delaunay: delaunay

include("addnoise.jl")

struct DelaunayTriangulation <: AbstractDelaunayTriangulation
    indices::Vector{Vector{Int32}}

    function DelaunayTriangulation(vertices::Vector{T}; joggle = 1e-8) where T <: Union{AbstractVector, SVector, MVector}
        # Slightly joggle points to avoid problems with QHull
        if joggle > 0
            addnoise!(vertices; joggle_factor = joggle)
        end
        simplex_inds = delaunay(vertices)
        new(simplex_inds)
    end

    DelaunayTriangulation(vertices::Dataset; joggle = 1e-8) = DelaunayTriangulation(vertices.data)
end