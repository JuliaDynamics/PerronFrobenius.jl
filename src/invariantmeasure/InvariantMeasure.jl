abstract type InvariantMeasure end
abstract type AbstractRectangularInvariantMeasure <: InvariantMeasure end

include("InvariantDistribution.jl")
include("compute_invariant_distribution.jl")
include("RectangularInvariantMeasure.jl")
include("average_invariant_measure.jl")
include("InducedRectangularInvariantMeasure.jl")
include("AverageRectangularInvariantMeasure.jl")

include("prettyprinting.jl")

export
InvariantMeasure,
AbstractRectangularInvariantMeasure
