using Reexport

@reexport module InvariantMeasures

    include("AbstractInvariantMeasure.jl")
	include("invariant_measure_from_transferoperator.jl")

    # Estimators and types for rectangular (box) partitions
    include("composite_types/RectangularInvariantMeasure/RectangularInvariantMeasure.jl")
    include("composite_types/InducedRectangularInvariantMeasure/InducedRectangularInvariantMeasure.jl")
    include("composite_types/AverageRectangularInvariantMeasure/AverageRectangularInvariantMeasure.jl")

    # Estimators and types for triangulation partitions
    include("composite_types/TriangulationExactInvariantMeasure/TriangulationExactInvariantMeasure.jl")
    include("composite_types/TriangulationApproxInvariantMeasure/TriangulationApproxInvariantMeasure.jl")

    # Plot recipes
    include("plot_recipes/plot_recipes.jl")
    
    include("invariantmeasure_estimators.jl")

end # module

"""
    InvariantMeasures

A module containing types and functions to estimate invariant measures from data by
using transfer operator estimators from the `TransferOperators` module.
"""
InvariantMeasures
