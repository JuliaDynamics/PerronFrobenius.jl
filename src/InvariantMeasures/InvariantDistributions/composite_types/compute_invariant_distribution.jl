import DelayEmbeddings.Dataset

import ..InvariantMeasures:
    TransferOperatorEstimatorRectangularBinning,
    TransferOperatorRectangularBinning
import CausalityToolsBase: encode, get_minima_and_edgelengths, RectangularBinning

function invariantmeasure(data, ϵ::RectangularBinning,
        estimator = :TransferOperatorEstimatorRectangularBinning;
        kwargs...) where T
    
    if estimator == :TransferOperatorEstimatorRectangularBinning
        # Encode the points (identify which bin, given the 
        # rectangular binning scheme ϵ, each point falls in)
        mini, edgelengths = get_minima_and_edgelengths(data, ϵ)
        encoded_pts = encode(data, mini, edgelengths)

        # Which are the visited bins, which points
        # visits which bin, repetitions, etc...
        binvisits = get_binvisits(encoded_pts)

        # Use that information to estimate transfer operator
        TO = estimate_transferoperator_from_binvisits(binvisits)

        return invariantmeasure(TO; kwargs...)
    end
end

function invariantmeasure(data::StateSpaceReconstruction.AbstractEmbedding,
        ϵ::RectangularBinning,
        estimator = :TransferOperatorEstimatorRectangularBinning;
        kwargs...)

    invariantmeasure(data.points, ϵ, estimator, kwargs...)
end

function invariantmeasure(data::Dataset,
        ϵ::RectangularBinning,
        estimator = :TransferOperatorEstimatorRectangularBinning;
        kwargs...)

    invariantmeasure(transpose(Matrix(data)), ϵ, estimator, kwargs...)
end

"""
    invariantmeasure(data,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator = TransferOperatorEstimatorRectangularBinning;
        kwargs...)

Estimate the invariant measure from a rectangular partition of `data` using
the estimator provided with the `estimator` argument. The default
is `estimator = :TransferOperatorEstimatorRectangularBinning`.

This is done by discretizing the state space into rectangular bins with
edge lengths dictated by the binning scheme `ϵ`. We then approximate the
transfer operator over the discretized state space, and compute the
invariant measure over the bins from the transfer operator.

Returns an [`InvariantDistribution`](@ref) instance.

## Arguments
- **`data`**: The data from which to estimate the invariant measure. The following
    data types are currently accepted:
    - `Dataset` instances from `DynamicalSystems.jl`/`DelayEmbeddings.jl`.
    - `AbstractEmbedding` subtypes from `StateSpaceReconstruction.jl`.
    - `AbstractArray{T, 2}` instances where each column represents a point.
- **`ϵ`**: A valid `RectangularBinning` binning scheme.
"""
function invariantmeasure(data,
        ϵ::RectangularBinning,
        estimator::Symbol; kwargs...)
    invariantmeasure(data, ϵ, estimator, kwargs...)
end

export invariantmeasure


"""
    invariantmeasure(to::AbstractTransferOperator;
            N::Int = 200,
            tolerance::Float64 = 1e-8,
            delta::Float64 = 1e-8)

Compute the invariant probability distribution from a transfer operator.

## Computing an invariant probability distribution
The distribution is taken as a left eigenvector of the transfer matrix,
obtained by repeated application of the transfer operator on a randomly
initialised distribution until the probability distribution converges.

## Keyword arguments
- `N`: the maximum number of iterations.
- `tolerance` and `delta`: decides when convergence is achieved.
"""
function invariantmeasure(to::AbstractTransferOperator;
			N::Int = 200,
            tolerance::Float64 = 1e-8,
            delta::Float64 = 1e-8)
    #=
    # Start with a random distribution `Ρ` (big rho). Normalise it so that it
    # sums to 1 and forms a true probability distribution over the simplices.
    =#
    Ρ = rand(Float64, 1, size(to.transfermatrix, 1))
    Ρ = Ρ ./ sum(Ρ, dims = 2)

    #=
    # Start estimating the invariant distribution. We could either do this by
    # finding the left-eigenvector of M, or by repeated application of M on Ρ
    # until the distribution converges. Here, we use the latter approach,
    # meaning that we iterate until Ρ doesn't change substantially between
    # iterations.
    =#
    distribution = Ρ * to.transfermatrix

    distance = norm(distribution - Ρ) / norm(Ρ)

    check = floor(Int, 1 / delta)
    check_pts = floor.(Int, transpose(collect(1:N)) ./ check) .* transpose(collect(1:N))
    check_pts = check_pts[check_pts .> 0]
    num_checkpts = size(check_pts, 1)
    check_pts_counter = 1

    counter = 1
    while counter <= N && distance >= tolerance
        counter += 1
        Ρ = distribution

        # Apply the Markov matrix to the current state of the distribution
        distribution = Ρ * to.transfermatrix

        if (check_pts_counter <= num_checkpts &&
           counter == check_pts[check_pts_counter])

            check_pts_counter += 1
            colsum_distribution = sum(distribution, dims = 2)[1]
            if abs(colsum_distribution - 1) > delta
                distribution = distribution ./ colsum_distribution
            end
        end

        distance = norm(distribution - Ρ) / norm(Ρ)
    end
    distribution = dropdims(distribution, dims = 1)

    # Do the last normalisation and check
    colsum_distribution = sum(distribution)

    if abs(colsum_distribution - 1) > delta
        distribution = distribution ./ colsum_distribution
    end

    # Find partition elements with strictly positive measure.
    δ = tolerance/size(to.transfermatrix, 1)
    simplex_inds_nonzero = findall(distribution .> δ)

    # Extract the elements of the invariant measure corresponding to these indices
    return PerronFrobenius.InvariantDistribution(distribution,simplex_inds_nonzero)
end

export invariantmeasure


"""
    invariantmeasure(to::AbstractTransferOperator;
            N::Int = 200,
            tolerance::Float64 = 1e-8,
            delta::Float64 = 1e-8)

Compute the invariant probability distribution from a transfer operator.

## Computing an invariant probability distribution
The distribution is taken as a left eigenvector of the transfer matrix,
obtained by repeated application of the transfer operator on a randomly
initialised distribution until the probability distribution converges.

## Keyword arguments
- `N`: the maximum number of iterations.
- `tolerance` and `delta`: decides when convergence is achieved.
"""
function invariantmeasure(to::AbstractTransferOperator;
			N::Int = 200,
            tolerance::Float64 = 1e-8,
            delta::Float64 = 1e-8)
    #=
    # Start with a random distribution `Ρ` (big rho). Normalise it so that it
    # sums to 1 and forms a true probability distribution over the simplices.
    =#
    Ρ = rand(Float64, 1, size(to.transfermatrix, 1))
    Ρ = Ρ ./ sum(Ρ, dims = 2)

    #=
    # Start estimating the invariant distribution. We could either do this by
    # finding the left-eigenvector of M, or by repeated application of M on Ρ
    # until the distribution converges. Here, we use the latter approach,
    # meaning that we iterate until Ρ doesn't change substantially between
    # iterations.
    =#
    distribution = Ρ * to.transfermatrix

    distance = norm(distribution - Ρ) / norm(Ρ)

    check = floor(Int, 1 / delta)
    check_pts = floor.(Int, transpose(collect(1:N)) ./ check) .* transpose(collect(1:N))
    check_pts = check_pts[check_pts .> 0]
    num_checkpts = size(check_pts, 1)
    check_pts_counter = 1

    counter = 1
    while counter <= N && distance >= tolerance
        counter += 1
        Ρ = distribution

        # Apply the Markov matrix to the current state of the distribution
        distribution = Ρ * to.transfermatrix

        if (check_pts_counter <= num_checkpts &&
           counter == check_pts[check_pts_counter])

            check_pts_counter += 1
            colsum_distribution = sum(distribution, dims = 2)[1]
            if abs(colsum_distribution - 1) > delta
                distribution = distribution ./ colsum_distribution
            end
        end

        distance = norm(distribution - Ρ) / norm(Ρ)
    end
    distribution = dropdims(distribution, dims = 1)

    # Do the last normalisation and check
    colsum_distribution = sum(distribution)

    if abs(colsum_distribution - 1) > delta
        distribution = distribution ./ colsum_distribution
    end

    # Find partition elements with strictly positive measure.
    δ = tolerance/size(to.transfermatrix, 1)
    simplex_inds_nonzero = findall(distribution .> δ)

    # Extract the elements of the invariant measure corresponding to these indices
    return PerronFrobenius.InvariantDistribution(distribution, simplex_inds_nonzero)
end
