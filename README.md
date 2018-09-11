# PerronFrobenius.jl

Implementations of Ulam's approximation of the transfer operator (Perron-Frobenius operator).

`PerronFrobenius.jl` is part of the  [`CausalityTools.jl`](https://github.com/kahaaga/CausalityTools.jl package.
Given an embedding (state space reconstruction) of a set of observations,
it can compute the transfer (Perron-Frobenius) operator and invariant measures
of partitioned state spaces.

## Basic usage
The transfer operator can be computed from subtypes of `Partition` from `StateSpaceReconstruction.jl`.

The basis for computing the transfer operator using this package is a set of time series.

```julia
# Create three red noise time series.
ts = [cumsum(randn(100)) for i = 1:3]
```

Now, create a state space reconstruction (embedding) of the red noise time series.

```
using PerronFrobenius
E = embed(ts, [1, 2, 3], [0, 0, 0]) # embed all time series with no relative lag
```

### Coarse graining

Then the transfer operator can be computed from variuos coarse grainings (partitions) of the state space.

```julia
# Partition the reconstruction by subdividing into disjoint simplices
triang = triangulate(E) #

# A simplex partitioning is not always invariant under the forward linear map, which may bias the estimate. In this case, we can
# adjust the last point of the embedding so that it falls inside the convex hull of the preceding points. This way, information
# is not lost when taking powers of the transfer operator.
triang_lininv = triangulate(invariantize(E))

# Partition the embedding into a rectangular box grid.
equibinning = bin_equidistant(E, 5) # bin into rectangular boxes, five boxes along each dimension
```


## Transfer operator estimation

The `transferoperator` function computes the Perron Frobenius operator. There are different estimators for each type of state space partition, but these are all available from `transferoperator`. Currently, the available partitionings are `Triangulation`, `LinearlyInvariantTriangulation` and `EquidistantBinning` (type is inferred when partitioning from different embeddings).

### 1. From invariant simplex partitions
The simplex intersection approximation of the transfer operator was used in [1]. This package handles embedding of arbitrary dimension by utilising [`SimplexIntersection.jl`](https://github.com/kahaaga/StateSpaceReconstruction.jl).

The transfer operator is guaranteed to be Markov if computed from a `LinearlyInvariantTriangulation`. You have some options as to how
the operator is estimated. The approximate approach is fastest for simplex partitions. It produces (with default settings) estimates that usually deviate less than 10% from exact estimates.

```julia
to_triang = transferoperator(triang_lininv, exact = false) # transition probabilities computed by approximate simplex intersection
to_triang = transferoperator(triang_lininv, exact = true)  # transition probabilities computed by exact simplex intersection
to_triang = transferoperator(triang_lininv, exact = true, parallel = true)  # exact intersection, run in parallel (this is the default).
to_triang = transferoperator(triang_lininv, exact = true, parallel = false)  # exact intersection, don't run in parallel.
```
### 2. From simplex partitions
The same works for the regular `Triangulation`, but the resulting operator is not guaranteed to be fully Markov (but will likely be
very close).

```julia
to_triang = transferoperator(triang, exact = false) # transition probabilities computed by approximate simplex intersection
to_triang = transferoperator(triang, exact = true)  # transition probabilities computed by exact simplex intersection
to_triang = transferoperator(triang, exact = true, parallel = true)  # exact intersection, run in parallel (this is the default).
to_triang = transferoperator(triang, exact = true, parallel = false)  # exact intersection, don't run in parallel.
```
### 2. From rectangular bin partitions
The fastest estimator dispatches on `EquidistantBinning`s, which is the recommended partition method for all but the shortest time series (< ~300 pts).

```julia
to_equibin = transferoperator(equibinning)
```

## References

[^1]: Gary Froyland, Bulletin of the Australian Mathematical Society 56, 157 (1997).
