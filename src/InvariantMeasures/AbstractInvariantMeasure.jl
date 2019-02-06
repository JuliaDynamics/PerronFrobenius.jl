
include("InvariantDistributions/InvariantDistributions.jl")

abstract type InvariantMeasure end
abstract type AbstractRectangularInvariantMeasure <: InvariantMeasure end
abstract type AbstractTriangulationInvariantMeasure <: InvariantMeasure end

abstract type InvariantMeasureEstimator end


export
    InvariantMeasure,
    InvariantMeasureEstimator,
    AbstractRectangularInvariantMeasure,
    AbstractTriangulationInvariantMeasure