export TransferOperatorApproximation
import Printf: @sprintf

"""
    TransferOperatorApproximation(g, M, params)

An approximation to the transfer operator over some partition.

The transfer operator generator `g` contains the input data, 
initialization steps, partition specification, and the transition
probability estimation method.
Along with the additional `params` (a named tuple), this information
is used to compute the transfer matrix `M`. 

`M`  is a row-stochastic (and partition-dependent)
Markov matrix containing the approximated transition probabilities 
between partition elements. If the input data were an orbit of 
some dynamical system, then every partition element corresponds to 
a unique (discretized) state of the system. Hence, the element ``M_{ij}`` 
gives the probability of jumping between state ``i`` to state ``j`` in 
any particular time step. 
"""
struct TransferOperatorApproximation{G <: TransferOperatorGenerator, T <: Real, P <: NamedTuple}
    g::G
    M::AbstractArray{T, 2}
    params::P
end

function Base.summary(to::T) where T<:TransferOperatorApproximation
    method = to.g.method
    data = to.g.pts 
    s = size(to.M)
    l = length(to.M)
    percent_nonzero = @sprintf("%.4f", count(to.M .> 0.0)/length(to.M) * 100)

    sM = "  Transfer matrix: size $s with $percent_nonzero% nonzero entries"
    sdata = "  Input data: $data"
    smethod = "  Partition and estimation method: $method"
    return "TransferOperatorApproximation\n$(sM)\n$(sdata)\n$(smethod)"
end

Base.show(io::IO, to::T) where {T<:TransferOperatorApproximation} = println(io, summary(to))
