
##################################################
# Check if the obtained transfer matrices are
# markov.
##################################################

"""
    is_markov(M::AbstractArray{T, 2}) -> Bool

Is the matrix Markov? """
function is_markov(M::AbstractArray{T, 2}) where T
    all(sum(M, dims = 2) .≈ 1)
end

"""
    is_markov(TO::AbstractTransferOperator) -> Bool

Is the transfer operator Markov? """
function is_markov(TO::AbstractTransferOperator)
    is_markov(TO.transfermatrix)
end

"""
    is_almostmarkov(M::AbstractArray{T, 2}; tol = 0.01) -> Bool

Is the matrix almost Markov?
"""
function is_almost_markov(M::AbstractArray{T, 2}; tol = 1e-3) where T
    all(sum(M, dims = 2) .> (1 - tol))
end

"""
    is_almostmarkov(TO::AbstractTransferOperator; tol = 0.01) -> Bool

Is the transfer operator almost Markov?
"""
function is_almost_markov(TO::AbstractTransferOperator; tol = 1e-3)
    is_almost_markov(TO.transfermatrix)
end


function zerocols(M::AbstractArray{T, 2}) where T
    findall(sum(M, dims = 1) .== 0)
end

function zerorows(M::AbstractArray{T, 2}) where T
    findall(sum(M, dims = 2) .== 0)
end

function zerocols(TO::AbstractTransferOperator)
    findall(sum(TO.transfermatrix, dims = 1) .== 0)
end

function zerorows(TO::AbstractTransferOperator)
    findall(sum(TO.transfermatrix, dims = 2) .== 0)
end


export
zerocols,
zerorows,
is_markov,
is_almost_markov



# #####################################################
# # How much do the exact and approximate triangulation
# # estimators differ?
# #####################################################
# """
#     max_discrep(to1::ExactSimplexTransferOperator,
#                 to2::ApproxSimplexTransferOperator)
#
# Measures the maximum discrepancy between an exact and an approximate estimate
# of the transfer operator.
#
# For this function to work, the transfer operators must have been generated from
# the same triangulation (otherwise, the number of simplices don't match up).
# """
# function max_discrep(exact::ExactSimplexTransferOperator,
#                      approx::ApproxSimplexTransferOperator)
#
# end
