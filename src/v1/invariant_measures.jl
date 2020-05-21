export invariantmeasure, InvariantDistribution

import LinearAlgebra: norm

# Invariant measures (measures over the partition do not have to sum to 1)
abstract type InvariantMeasure end
abstract type AbstractRectangularInvariantMeasure <: InvariantMeasure end
abstract type AbstractTriangulationInvariantMeasure <: InvariantMeasure end

# Invariant distributions (measures over the partition must sum to 1)
abstract type AbstractInvariantDistribution{TO <: TransferOperatorApproximation, T} end

"""
    invariantmeasure(to::TransferOperatorApproximation) → iv::InvariantDistribution
    invariantmeasure(pts, method::SingleGrid) → iv::InvariantDistribution
    invariantmeasure(pts, method::SimplexPoint) → iv::InvariantDistribution
    invariantmeasure(pts, method::SimplexExact) → iv::InvariantDistribution

Compute an invariant probability distribution `iv` from an approximation to the 
transfer operator over a partitioning of the input `pts` using one of the 
available `method`s. 

Alternatively, compute the distribution from a precomputed transfer operator approximation.

## Examples 

### `SingleGrid` estimator 

Use the [`SingleGrid`](@ref) estimator for long time series. 

```julia
using DelayEmbeddings, PerronFrobenius
n = 20000
x = [cos(2π*(t + sin(t)/30)) for t in 1:n] .+ 0.2rand(n)
τs, js = (0, -2, -3), (1, 1, 1)
pts = genembed(x, τs, js);

# Direct approach
b = RectangularBinning(5)
invariantmeasure(pts, SingleGrid(b))

# Indirect approach
to = transferoperator(pts, SingleGrid(b))
invariantmeasure(to)
```

### Triangulation estimators

The [`SimplexApprox`](@ref) estimator can be useful for very short time series. 
The [`SimplexExact`](@ref) estimator can also be useful for very short time series.
It computes transition probabilities from exact simplex intersections, so 
is much slower than the `SimplexApprox` estimator.

```julia
using DelayEmbeddings, PerronFrobenius
n = 30
x = [cos(2π*(t + sin(t)/30)) for t in 1:n] .+ 0.2rand(n)
τs, js = (0, -2, -3), (1, 1, 1)
pts = genembed(x, τs, js);

# Use approximate simplex intersections. Sample each simplex using 
# at least 300 points
invariantmeasure(pts, SimplexPoint(), n = 300)

# Use exact simplex intersections.
invariantmeasure(pts, SimplexExact())
```
"""
function invariantmeasure end 

"""
    InvariantDistribution(to::TransferOperatorApproximation, dist, nonzero_inds)

An invariant distribution `dist` computed from the transfer operator approximation `to`.
`nonzero_inds` contains the indices of `dist` with nonzero entries.
"""
struct InvariantDistribution{TO <: TransferOperatorApproximation, T} <: AbstractInvariantDistribution{TO, T}
    to::TO
    dist::Vector{T} # Distribution over the simplices
    nonzero_inds::Vector{Int} # indices of nonzero entries
end

function Base.summary(iv::InvariantDistribution{TO, T}) where {TO, T}
    method = iv.to.g.method
    data = iv.to.g.pts 
    s = size(iv.to.M, 1)
    l = length(iv.to.M)
    percent_nonzero = @sprintf("%.4f", count(iv.to.M .> 0.0)/length(iv.to.M) * 100)
    f = length(iv.nonzero_inds)/length(iv.dist)
    return "InvariantDistribution over $s partition elements ($(f*100)% with positive measure)"
end

function Base.show(io::IO, iv::InvariantDistribution{TO, T}) where {TO, T}
    s = summary(iv)
    print(io, s)
end

function invariantmeasure(to::TransferOperatorApproximation; 
        N::Int = 200, tol::Float64 = 1e-8, δ::Float64 = 1e-8)
    M = to.M
    #=
    # Start with a random distribution `Ρ` (big rho). Normalise it so that it
    # sums to 1 and forms a true probability distribution over the simplices.
    =#
    Ρ = rand(Float64, 1, size(M, 1))
    Ρ .= Ρ ./ sum(Ρ, dims = 2)

    #=
    # Start estimating the invariant distribution. We could either do this by
    # finding the left-eigenvector of M, or by repeated application of M on Ρ
    # until the distribution converges. Here, we use the latter approach,
    # meaning that we iterate until Ρ doesn't change substantially between
    # iterations.
    =#
    distribution = Ρ * M

    distance = norm(distribution - Ρ) / norm(Ρ)

    check = floor(Int, 1 / δ)
    check_pts = floor.(Int, transpose(collect(1:N)) ./ check) .* transpose(collect(1:N))
    check_pts = check_pts[check_pts .> 0]
    num_checkpts = size(check_pts, 1)
    check_pts_counter = 1

    counter = 1
    while counter <= N && distance >= tol
        counter += 1
        Ρ = distribution

        # Apply the Markov matrix to the current state of the distribution
        distribution = Ρ * M

        if (check_pts_counter <= num_checkpts &&
           counter == check_pts[check_pts_counter])

            check_pts_counter += 1
            colsum_distribution = sum(distribution, dims = 2)[1]
            if abs(colsum_distribution - 1) > δ
                distribution = distribution ./ colsum_distribution
            end
        end

        distance = norm(distribution - Ρ) / norm(Ρ)
    end
    distribution = dropdims(distribution, dims = 1)

    # Do the last normalisation and check
    colsum_distribution = sum(distribution)

    if abs(colsum_distribution - 1) > δ
        distribution = distribution ./ colsum_distribution
    end

    # Find partition elements with strictly positive measure.
    Δ = δ/size(M, 1)
    positive_measure_inds = findall(distribution .> Δ)

    return InvariantDistribution(to, distribution, positive_measure_inds)
end

