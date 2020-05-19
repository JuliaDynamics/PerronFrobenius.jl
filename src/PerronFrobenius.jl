
module PerronFrobenius
    using Requires 

    include("v1/api.jl")
    include("v1/boundary_conditions.jl")
    include("v1/gridestimators/Grid.jl")

    function __init__()
        @require Simplices="d5428e67-3037-59ba-9ab1-57a04f0a3b6a" begin
            include("v1/triangulationestimators/TriangulationEstimators.jl")
        end
    end
end
