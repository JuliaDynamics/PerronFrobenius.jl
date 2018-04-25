# PerronFrobenius.jl

Package for computing the transfer operator (Perron-Frobenius operator) from time series data. 

## Basic usage
The `julia transferentropy` function works on `Partition`s structs from the [`StateSpaceReconstruction.jl`](https://github.com/kahaaga/StateSpaceReconstruction.jl) package. Say you've done a state space reconstruction on some time series data, for example as follows.

```julia
# Create a set of red noise time series and do a state space reconstruction of them.
using PerronFrobenius

ts = [cumsum(randn(100)) for i = 1:3] 
E = embed(ts, [1, 2, 3], [0, 0, 0]) # embed all time series with no relative lag
```

Then the transfer operator can be computed on variuos coarse grainings (partitions) of the state space. 

```julia 
# Partition the reconstruction in various ways
triang = triangulate(E) 
triang_lininv = triangulate(invariantize(E)) # ensure embedding (and, hence, the triangulation) is invariant under forward linear map
equibinning = bin_equidistant(E, 5) # bin into rectangular boxes, five boxes along each dimension

# Estimate the transfer operator from the partitions
to_triang = transferoperator(triang)
to_triang = transferoperator(triang_lininv)
to_equibin = transferoperator(equibinning)
```

