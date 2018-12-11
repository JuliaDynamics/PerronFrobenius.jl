import DelayEmbeddings.Dataset

include("InvariantDistribution.jl")
include("../TransferOperator/AbstractTransferOperator.jl")
include("../TransferOperator/RectangularPartitionTransferOperators/RectangularBinningTransferOperator.jl")

function invariantmeasure(data::AbstractArray{T, 2},
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator = :transferoperator_grid;
        kwargs...) where T
    if estimator == :transferoperator_grid
        # Identify which bins of the partition resulting from using ϵ each
        # point of the embedding visits.

        # The indices, counting from the start of each coordinate axis
        # in steps given by ϵ
        visited_bins_inds = assign_bin_labels(data, ϵ)

        # The coordinate of the bin origins
        visited_bins_coordinates = assign_coordinate_labels(data, ϵ)

        # Which are the visited bins, which points
        # visits which bin, repetitions, etc...
        binvisits = organize_bin_labels(visited_bins_inds)

        # Use that information to estimate transfer operator
        TO = transferoperator_binvisits(binvisits)

        # Compute invariant measure
        ivm = invariantmeasure(TO; kwargs...)

        return ivm
    end
end

function invariantmeasure(data::StateSpaceReconstruction.AbstractEmbedding,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator = :transferoperator_grid;
        kwargs...)

    invariantmeasure(data.points, ϵ, estimator, kwargs...)
end

function invariantmeasure(data::Dataset,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator = :transferoperator_grid;
        kwargs...)

    invariantmeasure(transpose(Matrix(data)), ϵ, estimator, kwargs...)
end

"""
    invariantmeasure(data,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator = transferoperator_grid;
        kwargs...)

Estimate the invariant measure from a rectangular partition of `data` using
the estimator provided with the `estimator` argument. The default
is `estimator = :transferoperator_grid`.

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
- **`ϵ`**: The binning scheme. The following `ϵ` are valid:
    - `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
    - `ϵ::Float64` divides each axis into intervals of size `ϵ`.
    - `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
    - `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.
- **`estimator`**: A transfer operator estimator yielding a
    `PerronFrobenius.RectangularBinningTransferOperator`. Defaults to
     `:transferoperator_grid`.
- **`kwargs`**: Keyword arguments when calling `invariantmeasure` on the
    transfer operator.
"""
function invariantmeasure(data,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator::Symbol; kwargs...)
    invariantmeasure(data, ϵ, estimator, kwargs...)
end

export invariantmeasure


"""
    left_eigenvector(to::AbstractTransferOperator;
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
function left_eigenvector(to::AbstractTransferOperator;
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

export left_eigenvector


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
