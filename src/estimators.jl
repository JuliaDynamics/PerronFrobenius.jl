
# Different types of estimators based either on simplex partitioning or rectangular binning.
abstract type SimplexEstimator <: TransferOperatorEstimator end
abstract type RectangularBinningEstimator <: TransferOperatorEstimator end

# Load estimator functions.
include("simplexestimators/exact.jl")
include("simplexestimators/pointcloud.jl")
include("rectangularbinestimators/singlepoint.jl")

"""
    @with_kw struct ExactSimplexEstimator <: SimplexEstimator
        estimate::Function = TO_exact
        estimate_parallel::Function = TO_exact_parallel
    end
Estimates the Perron Frobenius operator (Transfer operator) from a simplex-based
partitioning of the state space. Calculates transition probabilities between volumes based
on *exact* intersections between the simplices.
"""
@with_kw struct ExactSimplexEstimator <: SimplexEstimator
    estimate::Function = TO_exact
    estimate_parallel::Function = TO_exact_parallel
end

"""
    @with_kw struct PointCloudSimplexEstimator <: SimplexEstimator
        estimate::Function = TO_pointcloud
    end

Estimates the Perron Frobenius operator (Transfer operator) from a simplex-based
partitioning of the state space. Calculates transition probabilities between volumes based
on approximate intersections between the simplices. This is done by approximating simplices
as a distribution of points contained within the simplex.

Uses `TO_pointcloud(t::Triangulation; n_randpts::Int = 100, sample_randomly::Bool = true)`
under the hood. To check the documentation run `?TO_pointcloud`.

## Examples
Say we want to estimate the transfer operator from some random set of points. We could do

```julia
# Define some triangulation
pts = rand(20, 3)

# Embed the points and triangulate them.
t = triang_from_embedding(Embedding(pts))
cloudestimator = PointCloudSimplexEstimator()

# Estimate the transfer operator
TO = cloudestimator.estimate(t::Triangulation, n_randpots )
```
"""
@with_kw struct PointCloudSimplexEstimator <: SimplexEstimator
    estimate::Function = TO_pointcloud
end


"""
    @with_kw struct RectBinEstimator <: RectangularBinningEstimator
        estimate::Function = TO_from_binning
    end

Estimates the Perron Frobenius operator (Transfer operator) by a regular binning of the
embedding (into rectangular bins with regular size along each axis.)
"""
@with_kw struct RectBinEstimator <: RectangularBinningEstimator
    estimate::Function = TO_from_binning
end
