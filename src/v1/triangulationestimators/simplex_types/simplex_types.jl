
module SimplexTypes
    using Requires
    function __init__()
        @require Simplices="d5428e67-3037-59ba-9ab1-57a04f0a3b6a" begin
            import .Simplices 
            
            include("AbstractSimplex.jl")
            include("ImmutableSimplex.jl")
            include("MutableSimplex.jl")
            include("simplex_subsampling.jl")
        end
    end
end
