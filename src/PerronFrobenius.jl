module PerronFrobenius

using Reexport
using StaticArrays, InplaceOps
import GroupSlices: groupslices, firstinds, groupinds
import Simplices: subsample_coeffs, simplexintersection
@reexport using SimplexSplitting
import SimplexSplitting: maybeintersecting_simplices
import Parameters: @with_kw

# Internal module imports
using TimeSeries: SingleTimeSeries



# Abstract estimator type
abstract type TransferOperatorEstimator end

include("estimators.jl")
include("transferoperator.jl")


export
    TransferOperator,
    TransferOperatorEstimator,

    SimplexEstimator,
    ExactSimplexEstimator,
    PointCloudSimplexEstimator,

    RectangularBinningEstimator,
    RectBinEstimator
end
