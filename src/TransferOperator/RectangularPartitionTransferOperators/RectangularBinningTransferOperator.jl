using Reexport

using StateSpaceReconstruction
using LinearAlgebra
using SparseArrays

include("../AbstractTransferOperator.jl")

##################################################
# Subtype for transfer operators estimated from
# a rectangular binning.
##################################################
struct RectangularBinningTransferOperator <: AbstractTransferOperator
    transfermatrix::AbstractArray{Float64, 2}
end


include("rectangularbinning_estimators/estimator_binvisits.jl")
include("rectangularbinning_estimators/estimator_grid.jl")

export RectangularBinningTransferOperator
