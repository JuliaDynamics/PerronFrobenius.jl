# PerronFrobenius.jl

Package for computing the transfer operator (Perron-Frobenius operator) from time series data. 

## Basic usage
The transfer operator can be computed from subtypes of `Partition` from the [`StateSpaceReconstruction.jl`](https://github.com/kahaaga/StateSpaceReconstruction.jl) package. This means that you have to make an state space reconstruction (embedding) of your data before computing the transfer operator. 

Say you've done a state space reconstruction on some time series data, for example as follows.

```julia
# Create a set of red noise time series and do a state space reconstruction of them.
ts = [cumsum(randn(100)) for i = 1:3] 
```

`PerronFrobenius` re-exports `StateSpaceReconstruction.jl`, so you can use everything available in that package directly by just doing `using PerronFrobenius`. First, do a state space reconstruction.

```
using PerronFrobenius 
E = embed(ts, [1, 2, 3], [0, 0, 0]) # embed all time series with no relative lag
```

### Coarse graining

Then the transfer operator can be computed from variuos coarse grainings (partitions) of the state space.  Currently, the available partitionings in `StateSpaceReconstruction.jl` are `Triangulation`, `LinearlyInvariantTriangulation` and `EquidistantBinning`. 

```julia 
# Partition the reconstruction by subdividing into disjoint simplices
triang = triangulate(E) 

# A simplex partitioning is not always invariant under the forward linear map, which may bias the estimate. In this case, we can 
# adjust the last point of the embedding so that it falls inside the convex hull of the preceding points. This way, information
# is not lost when taking powers of the transfer operator. 
triang_lininv = triangulate(invariantize(E)) 

# Partition the embedding into a rectangular box grid.
equibinning = bin_equidistant(E, 5) # bin into rectangular boxes, five boxes along each dimension
```


### Estimate the transfer operator from the partitions

The `transferoperator` function computes the Perron Frobenius operator. There are different estimators for each type of state space partition, but these are all available from `transferoperator`. 

The transfer operator is guaranteed to be Markov if computed from a `LinearlyInvariantTriangulation`. You have some options as to how 
the operator is estimated. The approximate approach is fastest for simplex partitions. It produces (with default settings) estimates that usually deviate less than 10% from exact estimates.

```julia
to_triang = transferoperator(triang_lininv, exact = false) # transition probabilities computed by approximate simplex intersection
to_triang = transferoperator(triang_lininv, exact = true)  # transition probabilities computed by exact simplex intersection
to_triang = transferoperator(triang_lininv, exact = true, parallel = true)  # exact intersection, run in parallel (this is the default).
to_triang = transferoperator(triang_lininv, exact = true, parallel = false)  # exact intersection, don't run in parallel.
```

The same works for the regular `Triangulation`, but the resulting operator is not guaranteed to be fully Markov (but will likely be 
very close).

```julia
to_triang = transferoperator(triang, exact = false) # transition probabilities computed by approximate simplex intersection
to_triang = transferoperator(triang, exact = true)  # transition probabilities computed by exact simplex intersection
to_triang = transferoperator(triang, exact = true, parallel = true)  # exact intersection, run in parallel (this is the default).
to_triang = transferoperator(triang, exact = true, parallel = false)  # exact intersection, don't run in parallel.
```

The fast estimator is the one that dispatches on `EquidistantBinning`s. 


to_equibin = transferoperator(equibinning)
```

