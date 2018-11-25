abstract type InvariantMeasure end

include("compute_invariant_distribution.jl")
include("RectangularInvariantMeasure.jl")

export
InvariantMeasure,
RectangularInvariantMeasure
