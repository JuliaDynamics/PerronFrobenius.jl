abstract type AbstractTransferOperator end

#########################
# Indexing
#########################
Base.size(to::AbstractTransferOperator) = size(to.TO)
Base.length(to::AbstractTransferOperator) = prod(size(to))
Base.getindex(to::AbstractTransferOperator, i::Vector) = to.TO[:, i]
Base.getindex(to::AbstractTransferOperator, i::Int, j::Int) = to.TO[i, j]
Base.getindex(to::AbstractTransferOperator, i::Colon, j::Int) = getindex(to, i)
Base.getindex(to::AbstractTransferOperator, i::Int, j::Colon) = to[i, :]

#########################
# Pretty printing
#########################
function Base.summary(to::T) where T<:AbstractTransferOperator
    s = size(to)
    l = length(to)
    percent_nonzero = @sprintf("%.4f", count(to.TO .> 0.0)/length(to) * 100)

    transferoperatortype = typeof(to)
    return "$transferoperatortype of size $s with $percent_nonzero% nonzero entries"
end

function matstring(to::T) where T<:AbstractTransferOperator
    return summary(to)
end

Base.show(io::IO, to::T) where {T<:AbstractTransferOperator} = println(io, matstring(to))

##################################################
# Subtype for transfer operators estimated using
# exact simplex volume intersections.
##################################################
type ExactSimplexTransferOperator <: AbstractTransferOperator
    TO::AbstractArray{Float64, 2}
end

##################################################
# Subtype for transfer operators estimated using
# approximate simplex volume intersections.
##################################################
type ApproxSimplexTransferOperator <: AbstractTransferOperator
    TO::AbstractArray{Float64, 2}
end

##################################################
# Subtype for transfer operators estimated from
# a rectangular binning.
##################################################
type RectangularBinningTransferOperator  <: AbstractTransferOperator
    TO::AbstractArray{Float64, 2}
end

"""
    is_markov(TO::AbstractTransferOperator) -> Bool

Is the transfer operator Markov? """
is_markov(TO::T) where T <: AbstractTransferOperator = all(sum(TO.TO, 2) .≈ 1)

"""
    is_almostmarkov(TO::AbstractTransferOperator; tol = 0.01) -> Bool

Is the transfer operator almost Markov?
"""
is_almostmarkov(TO::AbstractTransferOperator; tol = 0.01) = all(sum(TO.TO, 2) .> (1 - tol))

# Load estimator functions.
include("rectangularbinestimators/equidistant.jl")
include("simplexestimators/exact.jl")
include("simplexestimators/pointapprox.jl")

"""
    transferoperator(triang::T; exact = false, parallel = true,
                        n_pts = 200, sample_randomly = false)
                        where T<:AbstractTriangulation

Estimate the transfer operator from a triangulation.
"""
function transferoperator(triang::T; exact = false, parallel = true,
                        n_pts = 200, sample_randomly = false) where T<:AbstractTriangulation
    if exact
        if parallel
            transferoperator_exact_p(triang)
        else
            transferoperator_exact(triang)
        end
    else
        transferoperator_approx(triang, n_pts = n_pts, sample_randomly = sample_randomly)
    end
end


"""
    max_discrep(to1::ExactSimplexTransferOperator,
                to2::ApproxSimplexTransferOperator)

Measures the maximum discrepancy between an exact and an approximate estimate
of the transfer operator.

For this function to work, the transfer operators must have been generated from
the same triangulation (otherwise, the number of simplices don't match up).
"""
function max_discrep(exact::ExactSimplexTransferOperator,
                     approx::ApproxSimplexTransferOperator)

end
