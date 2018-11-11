if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end

using Test
using StateSpaceReconstruction
using PerronFrobenius

include("triangulations.jl")
include("rectangular.jl")
include("invariantmeasure.jl")
