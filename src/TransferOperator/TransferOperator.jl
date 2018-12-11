include("AbstractTransferOperator.jl")

include("RectangularPartitionTransferOperators/RectangularBinningTransferOperator.jl")
include("TriangulationTransferOperators/TriangulationTransferOperators.jl")


##################################################
# Wrapper combining the different estimators for
# triangulations.
##################################################
#include("transferoperator_triang.jl")

#include("simplexestimators/exact.jl")
#include("simplexestimators/pointapprox.jl")
