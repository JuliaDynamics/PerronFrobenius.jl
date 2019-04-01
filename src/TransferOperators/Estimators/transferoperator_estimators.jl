import CausalityToolsBase: 
    BinningScheme,
    RectangularBinningScheme,
    RectangularBinning,
    TriangulationBinningScheme,
    TriangulationBinning,
    SimplexIntersectionType, 
    ExactIntersection, 
    ApproximateIntersection

import StaticArrays: SVector, MVector 
import DelayEmbeddings: Dataset
import CausalityToolsBase: get_minima_and_edgelengths, encode

""" 
    transferoperator(points, binning_scheme::RectangularBinning; kwargs...)

Discretize `points` using the provided `binning_scheme` and compute the transfer operator 
over the partition elements. 


# Example

Assume we have enough points that a rectangular partition yields a good estimate of the 
invariant measure. Then the transfer operator over the partition can be computed as follows.

```julia 
pts = [rand(3) for i = 1:2000]

# Use rectangular boxes constructed by subdividing each coordinate 
# axis into 10 subintervals of equal length.
transferoperator(pts, RectangularBinning(10))
```
""" 
function transferoperator(points, binning_scheme::RectangularBinning; kwargs...) end

######################## 
# Triangulation binnings
######################## 
"""
    transferoperator(pts, ϵ::TriangulationBinning, 
        simplex_intersection_type::ApproximateIntersection; 
        n::Int = 200, sample_randomly::Bool = false) -> TransferOperatorTriangulationApprox

Estimate the invariant measure over the state space defined by `pts` using a triangulation 
of the phase space as the partition, using approximate simplex intersections to compute transition
probabilities between the states (simplices). 

`n` is the number of points that each simplex is 
sampled with, and `sample_randomly` indicates whether points used be sampled randomly 
within each simpled (`sample_randomly = true`) or by a regular simplex splitting routine 
(`sample_randomly = false`, which is default).

# Example

Assume we have sufficiently few points that a triangulation approach involving simplex 
intersections is computationally feasible to compute the invariant measure. Then 
a transfer operator can be computed as follows.

```julia 
pts = [rand(3) for i = 1:30]
transferoperator(pts, TriangulationBinning(), ApproximateIntersection())
```
"""
function transferoperator(pts, ϵ::TriangulationBinning, 
        simplex_intersection_type::ApproximateIntersection; 
        n::Int = 200, sample_randomly::Bool = false)
    
    transferoperator_triangulation_approx(pts, 
        n_sample_pts = n, 
        sample_randomly = sample_randomly)
end


"""
    transferoperator(pts, ϵ::TriangulationBinning, 
        simplex_intersection_type::ExactIntersection) -> TransferOperatorTriangulationExact

Estimate the invariant measure over the state space defined by `pts` using a triangulation 
of the phase space as the partition, using exact simplex intersections to compute transition
probabilities between the states (simplices).

# Example

Assume we have sufficiently few points that a triangulation approach involving simplex 
intersections is computationally feasible to compute the invariant measure. Then 
a transfer operator can be computed as follows.

```julia 
pts = [rand(3) for i = 1:30]
transferoperator(pts, TriangulationBinning(), ExactIntersection())
```
"""
function transferoperator(pts, ϵ::TriangulationBinning, 
        simplex_intersection_type::ExactIntersection)
    
    transferoperator_triangulation_exact(pts)
end

######################## 
# Rectangular binnings
######################## 

function transferoperator(points::Vector{T}, binning_scheme::RectangularBinning;  
    allocate_frac::Float64 = 1.0, boundary_condition = :none) where {T <:Union{SVector, MVector, Vector}}

    # Identify which bins of the partition resulting from using ϵ each
    # point of the embedding visits.
    mini, edgelengths = get_minima_and_edgelengths(points, binning_scheme)
    encoded_pts = encode(points, mini, edgelengths)

    #visited_bins = assign_bin_labels(points, binning_scheme.ϵ)

    # Which are the visited bins, which points
    # visits which bin, repetitions, etc...
    binvisits = get_binvisits(encoded_pts)

    # Use that information to estimate transfer operator
    estimate_transferoperator_from_binvisits(binvisits,
                        allocate_frac = allocate_frac,
                        boundary_condition = boundary_condition)
end


function transferoperator(points::AbstractArray{T, 2}, binning_scheme::RectangularBinning;  
        allocate_frac::Float64 = 1.0, boundary_condition = :none) where {T}

    if size(points, 1) > size(points, 2)
        points = transpose(points)
    end

    pts = [points[:, i] for i = 1:maximum(size(points))]

    # Identify which bins of the partition resulting from using ϵ each
    # point of the embedding visits.
    mini, edgelengths = get_minima_and_edgelengths(pts, binning_scheme)
    encoded_pts = encode(pts, mini, edgelengths)

    # Which are the visited bins, which points
    # visits which bin, repetitions, etc...
    binvisits = get_binvisits(encoded_pts)

    # Use that information to estimate transfer operator
    estimate_transferoperator_from_binvisits(binvisits,
                        allocate_frac = allocate_frac,
                        boundary_condition = boundary_condition)
end

function transferoperator(points::Vector{Vector{T}}, binning_scheme::RectangularBinning;  
    allocate_frac::Float64 = 1.0, boundary_condition = :none) where {T}

    transferoperator(hcat(points...,), binning_scheme, allocate_frac = allocate_frac,
        boundary_condition = boundary_condition)
end

function transferoperator(points::Vector{SVector{T}}, binning_scheme::RectangularBinning;  
    allocate_frac::Float64 = 1.0, boundary_condition = :none) where {T}

    transferoperator(Array(hcat(points...,)), binning_scheme, allocate_frac = allocate_frac,
        boundary_condition = boundary_condition)
end

function transferoperator(points::Vector{MVector{T}}, binning_scheme::RectangularBinning;  
    allocate_frac::Float64 = 1.0, boundary_condition = :none) where {T}

    transferoperator(Array(hcat(points...,)), binning_scheme, allocate_frac = allocate_frac,
        boundary_condition = boundary_condition)
end

function transferoperator(points::Dataset, binning_scheme::RectangularBinning;  
    allocate_frac::Float64 = 1.0, boundary_condition = :none) where {T}

    transferoperator(transpose(Matrix(points)), binning_scheme, allocate_frac = allocate_frac,
        boundary_condition = boundary_condition)
end


######################## 
# Triangulation binnings
########################

""" 
    transferoperator(points, binning_scheme::TriangulationBinningScheme,
        simplex_intersection::SimplexIntersectionType; kwargs...)

Discretize `points` using the provided `binning_scheme` and compute the transfer operator 
over the partition elements using the type of simplex intersections indicated by 
`simplex_intersection_type`.
""" 
transferoperator(invariant_pts, binning_scheme::TriangulationBinning, 
    simplex_intersection_type::SimplexIntersectionType; kwargs...)

export transferoperator