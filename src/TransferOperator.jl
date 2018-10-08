
abstract type AbstractTransferOperator end

#########################
# Indexing
#########################
Base.size(to::AbstractTransferOperator) = size(to.transfermatrix)
Base.length(to::AbstractTransferOperator) = prod(size(to))
Base.sum(to::AbstractTransferOperator, i::Int) = sum(to.transfermatrix, i)

#########################
# Pretty printing
#########################
function Base.summary(to::T) where T<:AbstractTransferOperator
    s = size(to)
    l = length(to)
    percent_nonzero = @sprintf("%.4f", count(to.transfermatrix .> 0.0)/length(to) * 100)

    transferoperatortype = typeof(to)
    return "$transferoperatortype of size $s with $percent_nonzero% nonzero entries"
end

function matstring(to::T) where T<:AbstractTransferOperator
    return summary(to)
end

Base.show(io::IO, to::T) where {T<:AbstractTransferOperator} =
    println(io, matstring(to))


##################################################
# Abstract type for transfer operators estimated
# from triangulations.
##################################################
type TriangulationTransferOperator <: AbstractTransferOperator end

##################################################
# Subtype for transfer operators estimated using
# exact simplex volume intersections.
##################################################
type ExactSimplexTransferOperator <: AbstractTransferOperator
    transfermatrix::AbstractArray{Float64, 2}
end

##################################################
# Subtype for transfer operators estimated using
# approximate simplex volume intersections.
##################################################
type ApproxSimplexTransferOperator <: AbstractTransferOperator
    transfermatrix::AbstractArray{Float64, 2}
end



##################################################
# Subtype for transfer operators estimated from
# a rectangular binning.
##################################################
type RectangularBinningTransferOperator  <: AbstractTransferOperator
    transfermatrix::AbstractArray{Float64, 2}
end


##################################################
# Implementations of the different estimators
##################################################
include("rectangularbinestimators/organize_binvisits.jl")
include("rectangularbinestimators/binvisits_estimator.jl")
include("rectangularbinestimators/gridestimator.jl")
include("simplexestimators/exact.jl")
include("simplexestimators/pointapprox.jl")


##################################################
# Check if the obtained transfer matrices are
# markov.
##################################################
"""
    is_markov(TO::AbstractTransferOperator) -> Bool

Is the transfer operator Markov? """
function is_markov(TO::AbstractTransferOperator)
    all(sum(TO.transfermatrix, 2) .≈ 1)
end

function is_markov(M::AbstractArray{T, 2}) where T
    all(sum(M, 2) .≈ 1)
end

"""
    is_almostmarkov(TO::AbstractTransferOperator; tol = 0.01) -> Bool

Is the transfer operator almost Markov?
"""
function is_almost_markov(TO::AbstractTransferOperator; tol = 1e-3)
    all(sum(TO.transfermatrix, 2) .> (1 - tol))
end

function is_almost_markov(M::AbstractArray{T, 2}; tol = 1e-3) where T
    all(sum(M, 2) .> (1 - tol))
end

function zerocols(M::AbstractArray{T, 2}) where T
    find(sum(M, 1) .== 0)
end

function zerorows(M::AbstractArray{T, 2}) where T
    find(sum(M, 2) .== 0)
end

function zerocols(TO::AbstractTransferOperator)
    find(sum(TO.transfermatrix, 1) .== 0)
end

function zerorows(TO::AbstractTransferOperator)
    find(sum(TO.transfermatrix, 2) .== 0)
end

export zerocols, zerorows
##################################################
# Wrapper combining the different estimators for
# triangulations.
##################################################

"""
    transferoperator(triang::AbstractTriangulation;
                exact = false, parallel = true,
                n_pts = 200, sample_randomly = false)

Estimate the transfer operator from a triangulation.
"""
function transferoperator(triang::AbstractTriangulation;
                            exact = false, parallel = true,
                            n_pts = 200, sample_randomly = false)
    if exact
        if parallel
            transferoperator_exact_p(triang)
        else
            transferoperator_exact(triang)
        end
    else
        transferoperator_approx(triang, n_pts = n_pts,
                                sample_randomly = sample_randomly)
    end
end


#####################################################
# How much do the exact and approximate triangulation
# estimators differ?
#####################################################
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
