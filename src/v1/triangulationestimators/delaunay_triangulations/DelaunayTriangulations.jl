using Reexport

@reexport module DelaunayTriangulations
    using Requires
    function __init__()
        @require Simplices="d5428e67-3037-59ba-9ab1-57a04f0a3b6a" begin
            include("AbstractDelaunayTriangulation.jl")
            include("composite_types/DelaunayTriangulation.jl")
        end 
    end
    #include("plot_recipes/plot_recipes.jl")
end

"""
    DelaunayTriangulations

A module handling the creation of Delaunay triangulations from data.
"""
DelaunayTriangulations
