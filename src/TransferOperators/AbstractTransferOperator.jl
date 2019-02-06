using Printf

abstract type AbstractTransferOperator end

#########################
# Indexing
#########################

ATO = AbstractTransferOperator

Base.size(to::AbstractTransferOperator) = size(to.transfermatrix)
Base.length(to::AbstractTransferOperator) = prod(size(to))
Base.sum(to::AbstractTransferOperator, i::Int) = sum(to.transfermatrix, i)

Base.getindex(to::ATO, i) = to.transfermatrix[i]
Base.getindex(to::ATO, i, j) = to.transfermatrix[i, j]

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

export AbstractTransferOperator
