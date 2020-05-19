using Reexport 

@reexport module TriangulationEstimators
    export transferoperatorgenerator
    import ..transferoperatorgenerator
    import ..TransferOperatorGenerator
    import ..TransferOperator
    import ..isboundarycondition

    # Requirements and useful types
    include("simplex_types/simplex_types.jl")
    include("delaunay_triangulations/DelaunayTriangulations.jl")

    # Triangulation estimators
    include("common.jl")
    include("approx/SimplexApprox.jl")
    include("exact/SimplexExact.jl")
end
