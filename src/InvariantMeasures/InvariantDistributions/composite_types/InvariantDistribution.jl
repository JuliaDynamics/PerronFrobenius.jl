"""
        InvariantDistribution(dist, nonzero_inds)

    Invariant visitiation frequencies over a partitioned state space. Fields are
    `dist` (the probability distribution) and `nonzero_inds` (the indices of `dist`
    with nonzero entries).
"""
struct InvariantDistribution <: AbstractInvariantDistribution
    dist::Vector{Float64} # Distribution over the simplices
    nonzero_inds::Vector{Int} # indices of nonzero entries
end

Base.length(iv::InvariantDistribution) = length(iv.dist)
Base.size(iv::InvariantDistribution) = size(iv.dist)

get_distribution(i::InvariantDistribution) = i.dist

export
InvariantDistribution,
get_distribution