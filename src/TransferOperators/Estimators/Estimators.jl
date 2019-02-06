using Reexport

@reexport module Estimators
	# Rectangular estimators
	using StateSpaceReconstruction
	using LinearAlgebra
	using SparseArrays

	############################
	# Rectangular estimators
	############################
    include("rectangularbinning/binvisits/estimator_binvisits.jl")
    include("rectangularbinning/grid/estimator_grid.jl")

	############################
	# Triangulation estimators
	############################
	include("triangulation/exact/transferoperator_triangulation_exact.jl")
	include("triangulation/approximate/transferoperator_triangulation_approximate.jl")
	#include("triangulation/approximate/SuperFast.jl")


	include("transferoperator_estimators.jl")

end # module
