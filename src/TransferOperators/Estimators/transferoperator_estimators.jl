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


""" 
    transferoperator(points, binning_scheme::RectangularBinningScheme; kwargs...)

Discretize `points` using the provided `binning_scheme` and compute the transfer operator 
over the partition elements. 
""" 
function transferoperator(points, binning_scheme::RectangularBinningScheme; kwargs...) end

######################## 
# Triangulation binnings
######################## 
"""
    transferoperator(pts, ϵ::TriangulationBinning, simplex_intersection_type::ApproximateIntersection; 
        n::Int = 200, sample_randomly::Bool = false) -> TransferOperatorTriangulationApprox

Estimate the invariant measure over the state space defined by `pts` using a triangulation 
of the phase space as the partition, using approximate simplex intersections to compute transition
probabilities between the states (simplices). 

`n` is the number of points that each simplex is 
sampled with, and `sample_randomly` indicates whether points used be sampled randomly 
within each simpled (`sample_randomly = true`) or by a regular simplex splitting routine 
(`sample_randomly = false`, which is default).

Example: If `pts = [rand(3) for i = 1:30`], then run
`transferoperator(pts, TriangulationBinning(), ApproximateIntersection())`. 
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

Example: If `pts = [rand(3) for i = 1:30`], then run
`transferoperator(pts, TriangulationBinning(), ExactIntersection())`. 
"""
function transferoperator(pts, ϵ::TriangulationBinning, 
        simplex_intersection_type::ExactIntersection)
    
    transferoperator_triangulation_exact(pts)
end

######################## 
# Rectangular binnings
######################## 

function transferoperator(points::AbstractArray{T, 2}, binning_scheme::RectangularBinning;  
        allocate_frac::Float64 = 1.0, boundary_condition = :none) where {T}

    if size(points, 1) > size(points, 2)
        points = transpose(points)
    end

    # Identify which bins of the partition resulting from using ϵ each
    # point of the embedding visits.
    visited_bins = assign_bin_labels(points, binning_scheme.ϵ)

    # Which are the visited bins, which points
    # visits which bin, repetitions, etc...
    binvisits = organize_bin_labels(visited_bins)

    # Use that information to estimate transfer operator
    TransferOperatorEstimatorRectangularBinVisits(binvisits,
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