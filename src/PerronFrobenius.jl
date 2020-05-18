#__precompile__(true)

module PerronFrobenius

include("v1/binvisits.jl")
include("v1/api.jl")
include("v1/Grid.jl")

#using Reexport
#@reexport using CausalityToolsBase

#import StateSpaceReconstruction:
#    GroupSlices,
#    Embeddings.AbstractEmbedding,
#    DelaunayTriangulation,
#    cembed,
#    invariantize


#import Simplices:
#    subsample_coeffs,
#    simplexintersection

#using Distributed
#using SharedArrays
#using StaticArrays
#using SparseArrays
#using InplaceOps
#using RecipesBase
#using Printf
#using LinearAlgebra
#using DelayEmbeddings

# Abstract estimator type
#abstract type TransferOperatorEstimator end

#include("BinningSchemes/BinningSchemes.jl")
#include("TransferOperators/TransferOperators.jl")
#include("InvariantMeasures/InvariantMeasures.jl")

#include("is_markov.jl")

end
