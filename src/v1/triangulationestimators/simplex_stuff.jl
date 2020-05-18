export MutableSimplex, orientation, AbstractSimplex, ∩, intersect
import LinearAlgebra: det 
import Simplices: orientation, radius, volume, centroid, simplexintersection
import Base: intersect

abstract type AbstractSimplex{D, T} end
abstract type AbstractImmutableSimplex{D, T} <: AbstractSimplex{D, T} end
abstract type AbstractMutableSimplex{D, T} <: AbstractSimplex{D, T} end


#######################
# Indexing
#######################
# Return the i-th point as column vector
Base.getindex(s::AbstractSimplex, i) = s.vertices[i]
Base.getindex(s::AbstractSimplex, i, j) = s.vertices[i][j]
Base.getindex(s::AbstractSimplex, i::Colon, j) = hcat(s.vertices[j]...,)
Base.getindex(s::AbstractSimplex, i::Colon, j::Colon) = hcat(s.vertices...,)

Base.length(s::AbstractSimplex) = length(s.vertices)
Base.size(s::AbstractSimplex) = length(s[1]), length(s)
Base.size(s::AbstractSimplex, i) = size(s)[i]
Base.IteratorSize(s::AbstractSimplex) = Base.HasLength()

Base.firstindex(s::AbstractSimplex) = 1
Base.lastindex(s::AbstractSimplex) = length(s)
Base.eachindex(s::AbstractSimplex) = Base.OneTo(length(s))
Base.iterate(s::AbstractSimplex, state = 1) = iterate(s.vertices, state)

###########################################################
# Vector{vertex} representation => Array representation
###########################################################
Base.Array(s::AbstractSimplex) = hcat(s.vertices...,)

#######################
# Sizes
#######################
dimension(s::AbstractSimplex) = length(s[1])
npoints(s::AbstractSimplex) = length(s)
nvertices(s::AbstractSimplex) = length(s)

export dimension, npoints, nvertices

#########################
# Orientation and volumes
#########################
# Orientation (convention: append rows of ones at top)
function orientation(s::AbstractSimplex)
    det([ones(1, dimension(s) + 1); s[:, :]])
end

""" Get the unscaled volume of the simplex. Divide by factorial(dim) to get the true volume."""
volume(s::AbstractSimplex) = abs(orientation(s))

#######################
# Centroid
#######################
function centroid(s::AbstractSimplex)
    D = dimension(s)
    # Results in a dim-by-1 matrix, but we just want a vector, so drop the last dimension
    centroid = dropdims(s[:, :] * (ones(D + 1, 1) / (D + 1)), dims = 2)
end

#######################
# Radius
#######################
function radius(s::AbstractSimplex)
    D = dimension(s)
    centroidmatrix = repeat(centroid(s), 1, D + 1)
    radius = sqrt(maximum(ones(1, D) * ((s[:, :] .- centroidmatrix) .^ 2)))
end


#function orientation(s::MutableSimplex)
#    det([ones(1, dimension(s) + 1); s[:, :]])
#end

function Base.intersect(s1::AbstractSimplex, s2::AbstractSimplex)
    simplexintersection(s1[:, :], s2[:, :])
end

"""
    ∩(s1::AbstractSimplex, s2::AbstractSimplex)

Compute the volume intersection between two simplices.
"""
∩(s1::AbstractSimplex, s2::AbstractSimplex) = intersect(s1, s2)


"""
    MutableSimplex{D, T}

Mutable simplex where vertices are represented by some type of abstract vector.
"""
mutable struct MutableSimplex{D, T} <: AbstractMutableSimplex{D, T}
    vertices::Vector{Vector{T}}

    
    """
    MutableSimplex(pts::Vector{<:AbstractVector{T}})

    Construct a simplex from a vector of vectors.
    """
    function MutableSimplex(pts::Vector{<:AbstractVector{T}}) where {T}
        if !(length(pts) == length(pts[1]) + 1)
            err = """ The input cannot be converted to a simplex.
            Vertices need to have `dim` elements, and there needs to be `dim + 1` vertices.
            """
            throw(DomainError(pts, err))    end
        dim = length(pts[1])
        new{dim, T}([pts[i] for i = 1:length(pts)])
    end


    """
    MutableSimplex(pts::AbstractArray{T, 2}) where T

    Construct a simplex from an array of size `(dim + 1)-by-(dim)` or
    `(dim)-by-(dim + 1)` (faster).
    """
    function MutableSimplex(pts::AbstractArray{T, 2}) where {T}
        s = size(pts)

        if (maximum(s) - minimum(s)) != 1
            err = """ The input cannot be converted to a simplex.
                    size(pts) must be (dim, dim + 1) or (dim + 1, dim).
                """
                throw(DomainError(pts, err))
        end

        if s[1] > s[2] # Rows are points
            dim = s[2]
            return new{dim, T}([pts[i, :] for i = 1:maximum(s)])
        end

        if s[2] > s[1] # Columns are points
            dim = s[1]
            return new{dim, T}([pts[:, i] for i = 1:maximum(s)])
        end
    end
end

# Overwriting the i-th vertex
function Base.setindex!(simplex::MutableSimplex, v, i)
    simplex[i] .= v
end

# Overwriting elements of the i-th vertex
function Base.setindex!(simplex::MutableSimplex, v, i::Int, j)
    if length(v) != length(j)
        err = """
            Trying to overwrite elements $j of vertex $i with $v, which does not
            have the same number of elements as the target.
        """
        throw(ArgumentError(err))
    end
    simplex[i][j] = v
end
