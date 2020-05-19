using Reexport

@reexport module DelaunayTriangulations
    include("AbstractDelaunayTriangulation.jl")
    include("composite_types/DelaunayTriangulation.jl")

    #include("plot_recipes/plot_recipes.jl")
end

"""
    DelaunayTriangulations

A module handling the creation of Delaunay triangulations from data.
"""
DelaunayTriangulations
