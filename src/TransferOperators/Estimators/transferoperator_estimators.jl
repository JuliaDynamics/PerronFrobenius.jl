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

function transferoperator(invariant_pts, binning_scheme::TriangulationBinning, 
        simplex_intersection_type::ApproximateIntersection; 
        n_sample_pts::Int = 200, sample_randomly::Bool = false)

    transferoperator_triangulation_approx(invariant_pts, n_sample_pts = n_sample_pts,
        sample_randomly = sample_randomly)
end


""" 
    transferoperator(points, binning_scheme::TriangulationBinningScheme,
        simplex_intersection::SimplexIntersectionType; kwargs...)

Discretize `points` using the provided `binning_scheme` and compute the transfer operator 
over the partition elements using the type of simplex intersections indicated by 
`simplex_intersection_type`.
""" 
transferoperator(invariant_pts, binning_scheme::TriangulationBinning, 
    simplex_intersection_type::SimplexIntersectionType; kwargs...)


function transferoperator(invariant_pts, binning_scheme::TriangulationBinning, 
        simplex_intersection_type::ApproximateIntersection; 
        n_sample_pts::Int = 200, sample_randomly::Bool = false)

    transferoperator_triangulation_approx(invariant_pts, n_sample_pts = n_sample_pts,
        sample_randomly = sample_randomly)
end


function transferoperator(invariant_pts, binning_scheme::TriangulationBinning, 
    simplex_intersection_type::ExactIntersection; 
    n_sample_pts::Int = 200, sample_randomly::Bool = false)

    transferoperator_triangulation_exact(invariant_pts)
end

export transferoperator