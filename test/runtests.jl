if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end

using Test
using CausalityToolsBase
using StateSpaceReconstruction
using PerronFrobenius


# Individual constructors (i.e. rectangularinvariantmeasure(pts, ϵ))
include("test_rectangularmeasure_constructor.jl")
#include("test_inducedmeasure.jl")
#include("test_inducedmeasure_constructor.jl")
#include("test_averagemeasure.jl")
#include("test_averagemeasure_constructor.jl")

# Test the estimators (i.e. invariantmeasure(pts, binningscheme))
include("test_invariantmeasure_estimators.jl")

include("test_gridestimator.jl")
include("test_transferoperator_rectangular.jl")

#include("triangulations.jl")
include("test_triangulationestimators.jl")