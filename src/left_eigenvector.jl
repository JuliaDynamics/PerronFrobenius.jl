using Parameters

abstract type AbstractInvariantDistribution end
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

"""
Compute the invariant probability distribution from a
[`TransferOperator`](@ref). The distribution is obtained by repeated
application of the transfer operator on an randomly initialised distribution
	left_eigenvector(TO::EquidistantBinningTransferOperator; N=100,
		tolerance=1/10^5, delta=1/10^5)

until the distribution converges.

Optional arguments are `N` (maximum number of iterations), `tolerance` and
`delta`, with the two latter deciding when convergence is achieved.

"""
function left_eigenvector(to::AbstractTransferOperator;
			N::Int = 100, tolerance::Float64 = 1/10^5, delta::Float64 = 1/10^5)
    #=
    # Start with a random distribution `Ρ` (big rho). Normalise it so that it
    # sums to 1 and forms a true probability distribution over the simplices.
    =#
    Ρ = rand(Float64, 1, size(to.TO, 1))
    Ρ = Ρ ./ sum(Ρ, 2)

    #=
    # Start estimating the invariant distribution. We could either do this by
    # finding the left-eigenvector of M, or by repeated application of M on Ρ
    # until the distribution converges. Here, we use the latter approach,
    # meaning that we iterate until Ρ doesn't change substantially between
    # iterations.
    =#
    distribution = Ρ * to.TO

    distance = norm(distribution - Ρ) / norm(Ρ)

    check = floor(Int, 1 / delta)
    check_pts = floor.(Int, collect(1:N).' ./ check) .* collect(1:N).'
    check_pts = check_pts[check_pts .> 0]
    num_checkpts = size(check_pts, 1)
    check_pts_counter = 1

    counter = 1
    while counter <= N && distance >= tolerance
        counter += 1
        Ρ = distribution

        # Apply the Markov matrix to the current state of the distribution
        distribution = Ρ * to.TO

        if (check_pts_counter <= num_checkpts &&
           counter == check_pts[check_pts_counter])

            check_pts_counter += 1
            colsum_distribution = sum(distribution, 2)[1]
            if abs(colsum_distribution - 1) > delta
                distribution = distribution ./ colsum_distribution
            end
        end

        distance = norm(distribution - Ρ) / norm(Ρ)
    end

    # Do the last normalisation and check
    colsum_distribution = sum(distribution, 2)[1]

    if abs(colsum_distribution - 1) > delta
        distribution = distribution ./ colsum_distribution
    end

    # Find partition elements with strictly positive measure.
    simplex_inds_nonzero = find(distribution .> (tolerance/size(to.TO, 1)))

    # Extract the elements of the invariant measure corresponding to these indices
    return PerronFrobenius.InvariantDistribution(vec(distribution),simplex_inds_nonzero)
end
