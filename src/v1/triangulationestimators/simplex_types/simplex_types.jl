using Reexport 

@reexport module SimplexTypes 
    include("AbstractSimplex.jl")
    include("ImmutableSimplex.jl")
    include("MutableSimplex.jl")
    include("simplex_subsampling.jl")
end
