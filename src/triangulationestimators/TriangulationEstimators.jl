using Reexport 

@reexport module TriangulationEstimators
    
    import ..TransferOperator
    import ..TransferOperatorGenerator
    import ..TransferOperatorApproximation
    import ..transopergenerator
    import ..transferoperator
    import ..isboundarycondition
    import ..invariantmeasure
    import ..InvariantDistribution
    
    using Requires
    function __init__()
        @require Simplices="d5428e67-3037-59ba-9ab1-57a04f0a3b6a" begin
            export TriangulationBasedTransferOperator
            
            """ Supertype of all triangulation-based transfer operator estimators """
            abstract type TriangulationBasedTransferOperator <: TransferOperator end

            # Requirements and useful types
            include("simplex_types/simplex_types.jl")
            include("delaunay_triangulations/DelaunayTriangulations.jl")

            # Triangulation estimators
            include("common.jl")
            include("point/SimplexPoint.jl")
            include("exact/SimplexExact.jl")
            include("invariant_measure_simplex.jl")
        end
    end
   

 
end
