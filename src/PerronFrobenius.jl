
module PerronFrobenius
    import CausalityToolsBase 
    export RectangularBinning
    
    include("v1/api.jl")
    include("v1/boundary_conditions.jl")
    include("v1/transfer_operators.jl")
    include("v1/invariant_measures.jl")
    include("v1/gridestimators/SingleGrid.jl")
    include("v1/invariant_measures_grid.jl")
    
    using Requires 
    function __init__()
        @require Simplices="d5428e67-3037-59ba-9ab1-57a04f0a3b6a" begin
            include("v1/triangulationestimators/TriangulationEstimators.jl")
        end

        @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
            #include("v1/plot_recipes/invariant_measures.jl")
        end
    end
    #import .TriangulationEstimators: SimplexPoint, SimplexExact 
    #export SimplexPoint, SimplexExact
end
