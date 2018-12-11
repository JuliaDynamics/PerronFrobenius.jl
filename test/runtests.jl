if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end

using Test
using StateSpaceReconstruction
using PerronFrobenius


include("test_rectangularmeasure_constructor.jl")
include("test_inducedmeasure.jl")
include("test_inducedmeasure_constructor.jl")
include("test_averagemeasure.jl")
include("test_averagemeasure_constructor.jl")


include("test_gridestimator.jl")
include("test_transferoperator_rectangular.jl")
