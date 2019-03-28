# Overview

For short time series, the most reliable estimates of the transfer operator are obtained by using a triangulation of the state space as our partition. This approach is computationally costly because it has to compute N-dimensional simplex intersections. However, it gives robust estimates of ergodic transition probabilities down to as little as a few hundred points.

For longer time series, we use a rectangular grid to discretize the embedding.
This approach is must faster, because no volume intersections have to be explicitly computed.

```@setup s
using CausalityTools
using Plots; pyplot()
```

## Short time series (triangulation estimator)
Embed some random short time series (n = 30), make sure that the embedding is
invariant, triangulate the embedding, and estimate the transfer operator on the
states (partition elements) of the discretized state space. We then compute the
invariant measure over the states from the transfer matrix.

```@repl s
ts = [diff(rand(30)) for i = 1:3]
E = invariantize(cembed(ts))
triang = triangulate(E)
TO = transferoperator(triang)
```

```@repl s
# Plot the transfer operator and the invariant distribution
maxprob = min(1, maximum(TO.TO)*1.1) # for plotting
heatmap(TO.TO, clims=(0, maxprob));
xlabel!("Target simplex # (j)");
ylabel!("Source simplex # (i)");
savefig("transferoperator-short-ex.svg"); nothing #hide
```

![](transferoperator-short-ex.svg)

## Long time series (gridding estimator)
The gridding approach is fast and can be used on a large number of points
in higher dimensions. Here, we compute the transfer operator for 5000
4-dimensional data points, subdividing each coordinate axis into
intervals of length 0.3.

```@repl s
data = rand(4, 5000)
TO = TransferOperatorEstimatorRectangularBinning(data, 0.7)
maxprob = min(1, maximum(TO.TO)*1.1);
heatmap(TO.TO, clims=(0, maxprob));
xlabel!("Target box # (j)");
ylabel!("Source box # (i)");
savefig("transferoperator-long-ex.svg"); nothing #hide
```

![](transferoperator-long-ex.svg)


This also works on embeddings from `StateSpaceReconstruction.jl`:

```@repl s
data = rand(3, 3000)
E = cembed(data)
TO = TransferOperatorEstimatorRectangularBinning(E, 0.3)
maxprob = min(1, maximum(TO.TO)*1.1);
heatmap(TO.TO, clims=(0, maxprob));
xlabel!("Target box # (j)");
ylabel!("Source box # (i)");
savefig("transferoperator-long-ex2.svg"); nothing #hide
```

![](transferoperator-long-ex2.svg)
