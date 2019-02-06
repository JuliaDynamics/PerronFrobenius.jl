using Reexport 

@reexport module BinningSchemes 

import CausalityToolsBase: RectangularBinningScheme

VALIDRECTBINS = Union{Int, Float64, Vector{Int}, Vector{Float64}, Tuple{Vector{Tuple{Float64,Float64}},Int64}}

"""
InducedRectangularBinning

Instructions for creating a rectangular target partition and a rectangular source partition.

## Fields 

- **`ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}}`**: The instructions for deciding 
the edge lengths of the rectangular boxes. The following `ϵ` are valid:
    1. `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
    2. `ϵ::Float` divides each axis into intervals of size `ϵ`.
    3. `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
    4. `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.
"""
struct InducedRectangularBinning <: RectangularBinningScheme
    ϵ_target_partition::VALIDRECTBINS
    ϵ_source_partition::VALIDRECTBINS
end

"""
AverageRectangularBinning

Instructions for constructing a rectangular target partition and multiple source 
partitions.

## Fields 

- **`ϵ_target_partition`::Union{Int, Float64, Vector{Int}, Vector{Float64}}`**: The 
instructions for deciding the edge lengths of the rectangular boxes of the target 
partition. The following `ϵ` are valid:
    1. `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
    2. `ϵ::Float` divides each axis into intervals of size `ϵ`.
    3. `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
    4. `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.

- **`ϵs_source_partitions::Vector{Union{Int, Float64, Vector{Int}, Vector{Float64}}}`**: The 
instructions for deciding the edge lengths of the rectangular boxes of the source 
partitions. The following `ϵ` are valid:
    1. `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
    2. `ϵ::Float` divides each axis into intervals of size `ϵ`.
    3. `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
    4. `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.
"""
struct AverageRectangularBinning <: RectangularBinningScheme
    ϵ_target_partition::VALIDRECTBINS
    ϵs_source_partitions::VALIDRECTBINS
end



export
InducedRectangularBinning,
AverageRectangularBinning



end # module 


""" 
    BinningSchemes

Defines rectangular binning schemes used in the PerronFrobenius package not present 
in CausalityToolsBase.
""" 
BinningSchemes