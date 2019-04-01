__precompile__(true)

module PerronFrobenius

using Reexport
@reexport using CausalityToolsBase
using StateSpaceReconstruction

import StateSpaceReconstruction:
    GroupSlices,
    Embeddings.AbstractEmbedding,
    DelaunayTriangulation,
    cembed,
    invariantize


import Simplices:
    subsample_coeffs,
    simplexintersection

using Distributed
using SharedArrays
using StaticArrays
using SparseArrays
using InplaceOps
using RecipesBase
using Printf
using LinearAlgebra
using DelayEmbeddings

# Abstract estimator type
abstract type TransferOperatorEstimator end

include("BinningSchemes/BinningSchemes.jl")
include("TransferOperators/TransferOperators.jl")
include("InvariantMeasures/InvariantMeasures.jl")

include("is_markov.jl")

end
