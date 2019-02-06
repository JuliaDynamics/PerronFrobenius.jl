using Reexport

@reexport module InvariantDistributions
    include("AbstractInvariantDistribution.jl")
    include("composite_types/InvariantDistribution.jl")
	include("plot_recipes.jl")
end # module

"""
    InvariantDistributions

A module containing abstract and composite types for invariant distributions. These
composite types hold the distributions obtained by the estimators in the parent
`InvariantMeasures` module.
"""
InvariantDistributions
