module PerronFrobenius

using Reexport
using StaticArrays
using InplaceOps
using RecipesBase

import Parameters:
            @with_kw,
            @unpack

using Simplices:
            subsample_coeffs,
            simplexintersection
@reexport using StateSpaceReconstruction

# Abstract estimator type
abstract type TransferOperatorEstimator end

include("TransferOperator.jl")
include("estimators.jl")
include("transferoperator.jl")
include("left_eigenvector.jl")

export
    # Transfer operator types
    TransferOperator,
    ApproxSimplexTransferOperator,
    ExactSimplexTransferOperator,
    EquidistantBinningTransferOperator,

    # Methods on transfer operator types
    is_markov,
    is_almostmarkov,
    InvariantDistribution, left_eigenvector,

    # Transfer operator estimators
    transferoperator,
    transferoperator_exact,
    transferoperator_exact_p,
    transferoperator_approx

end
