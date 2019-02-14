using Reexport

@reexport module TransferOperators

    # Abstract types and composite types for different estimators.
    include("AbstractTransferOperator.jl")

    include("composite_types/TransferOperatorRectangularBinning.jl")
    include("composite_types/TransferOperatorTriangulationExact.jl")
    include("composite_types/TransferOperatorTriangulationApproximate.jl")
    include("composite_types/plot_recipes.jl")

    # Constructors for the estimators.
    include("Estimators/Estimators.jl")


end # module

"""
    TransferOperators

A module containing types and functions to approximate the transfer operator
from numerical data. The estimators in this module are the basis for computing
invariant measures, which is done in the `InvariantMeasures` module.
"""
TransferOperators
