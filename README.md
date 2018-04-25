# PerronFrobenius.jl

Package for computing the transfer operator (Perron-Frobenius operator) from time series data. 

## Structure

```julia
E = embed([cumsum(randn(100)) for i = 1:3])
```
