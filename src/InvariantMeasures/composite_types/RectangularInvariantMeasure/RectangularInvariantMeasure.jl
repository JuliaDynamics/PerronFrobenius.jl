import StateSpaceReconstruction:
    AbstractEmbedding,
    assign_bin_labels,
    assign_coordinate_labels,
    TransferOperatorEstimatorRectangularBinVisits

import ..TransferOperators:
    BinVisits,
    TransferOperatorRectangularBinning,
    organize_bin_labels,
    TransferOperatorEstimatorRectangularBinVisits

import DelayEmbeddings:
    Dataset

import StaticArrays:
    SVector, MVector


export RectangularInvariantMeasure, rectangularinvariantmeasure




"""
    RectangularInvariantMeasure{T} <: AbstractRectangularInvariantMeasure where {T}

A RectangularInvariantMeasure created from a set of points representing a
state space by discretizing the state space into rectangular bins with
edge lengths dictated by the binning scheme `ϵ`. The invariant measure is
then computed from an approximation to the transfer operator over the
discretized state space.

## Fields
- **`points::AbstractArray{T, 2}`**: The points for which to estimate the invariant measure.
    Each column is a point.

- **`ϵ`**: The binning scheme. The following `ϵ` are valid:
        1. `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
        2. `ϵ::Float` divides each axis into intervals of size `ϵ`.
        3. `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
        4. `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.

- **`ϵ_absolute`**: `ϵ` converted to absolute edge lengths.

- **`visited_bins_inds::AbstractArray{Int,2}`**: Counting from the start of each
    coordinate axis in steps dictated by `ϵ`, which bins do each
    of the points lie in? One column vector of indices for each point in `points`
    (one index per coordinate axis).

- **`visited_bins_coordinates::AbstractArray{Float64,2}`**: The coordinates of the origin
    of each bin that is visited. Each column vector is a bin origin. The bins (rectangular
    boxes) have edge lengths given by `ϵ`.

- **`binvisits::BinVisits`**: A `BinVisits` instance, indicating which points
    visits which bin. The indices in this `BinVisits` instance refers to the unique columns 
    of `visited_bins_inds` (one column for each point in `points`, of which some may 
    be identical), or, equivalently, 
    `visited_bins_coordinates` (also one column for each point in `points`).

- **`transfermatrix::TransferOperatorRectangularBinning`**: The transfer matrix from which
    the invariant measure is obtained.

- **`measure::InvariantDistribution`**: The invariant measure over the visited bins.
"""
struct RectangularInvariantMeasure{T} <: AbstractRectangularInvariantMeasure where {T}
    points::AbstractArray{T, 2}
    ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}}
    ϵ_absolute::Vector{Float64}
    visited_bins_inds::AbstractArray{Int,2}
    visited_bins_coordinates::AbstractArray{Float64,2}
    binvisits::BinVisits
    transfermatrix::TransferOperatorRectangularBinning
    measure::InvariantDistribution
end


"""
    rectangularinvariantmeasure(data,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator = TransferOperatorEstimatorRectangularBinning;
        kwargs...)

Estimate the invariant measure from a rectangular partition of `data`.

This is done by discretizing the state space into rectangular bins with
edge lengths dictated by the binning scheme `ϵ`. We then approximate the
transfer operator over the discretized state space, and compute the
invariant measure over the bins from the transfer operator.

Returns a `RectangularInvariantMeasure` instance.

## Arguments
- **`data`**: The data from which to estimate the invariant measure. The following data
    types are currently accepted:
    - `Dataset` instances from `DynamicalSystems.jl`/`DelayEmbeddings.jl`.
    - `AbstractEmbedding` subtypes from `StateSpaceReconstruction.jl`.
    - `AbstractArray{T, 2}` instances where each column represents a point.
    - `Vector{Vector{T}}`
    - `Vector{SVector{D, T}}`
    - `Vector{MVector{D, T}}`

- **`ϵ`**: The binning scheme. The following `ϵ` are valid:
    - `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
    - `ϵ::Float64` divides each axis into intervals of size `ϵ`.
    - `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
    - `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.

- **`estimator`**: A transfer operator estimator yielding a
    `TransferOperatorRectangularBinning`. Defaults to
    `:TransferOperatorEstimatorRectangularBinning`.

- **`kwargs`**: Keyword arguments when calling `invariantmeasure` on the
    transfer operator.
"""
function rectangularinvariantmeasure end

function rectangularinvariantmeasure(data::AbstractArray{T, 2},
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator::Symbol = :TransferOperatorEstimatorRectangularBinning;
        kwargs...) where {T}

    if estimator == :TransferOperatorEstimatorRectangularBinning
        # Identify which bins of the partition resulting from using ϵ each
        # point of the embedding visits.

        # Get absolute bin sizes 
        ϵ_absolute = minima_and_stepsizes(data, ϵ)[2]

        # The indices, counting from the start of each coordinate axis
        # in steps given by ϵ
        visited_bins_inds = assign_bin_labels(data, ϵ)

        # The coordinate of the bin origins
        visited_bins_coordinates = assign_coordinate_labels(data, ϵ)

        # Which are the visited bins, which points
        # visits which bin, repetitions, etc...
        binvisits = organize_bin_labels(visited_bins_inds)

        # Use that information to estimate transfer operator
        TO = TransferOperatorEstimatorRectangularBinVisits(binvisits)

        # Compute invariant measure
        ivm = invariantmeasure(TO; kwargs...)

        RectangularInvariantMeasure(
            data,
            ϵ,
            ϵ_absolute,
            visited_bins_inds,
            visited_bins_coordinates,
            binvisits,
            TO,
            ivm
        )
    end
end

function rectangularinvariantmeasure(data::Dataset,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator::Symbol = :TransferOperatorEstimatorRectangularBinning;
        kwargs...)

    rectangularinvariantmeasure(transpose(Matrix(data)), ϵ, estimator, kwargs...)
end

function rectangularinvariantmeasure(data::AbstractEmbedding,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator = :TransferOperatorEstimatorRectangularBinning;
        kwargs...)

    rectangularinvariantmeasure(data.points, ϵ, estimator, kwargs...)
end

function rectangularinvariantmeasure(data::Vector{Vector{T}},
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator = :TransferOperatorEstimatorRectangularBinning;
        kwargs...) where {T}

    rectangularinvariantmeasure(hcat(data...,), ϵ, estimator, kwargs...)
end

function rectangularinvariantmeasure(data::Vector{SVector{D, T}},
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator = :TransferOperatorEstimatorRectangularBinning;
        kwargs...) where {D, T}

    rectangularinvariantmeasure(Array(hcat(data...,)), ϵ, estimator, kwargs...)
end


function rectangularinvariantmeasure(data::Vector{MVector{D, T}},
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        estimator = :TransferOperatorEstimatorRectangularBinning;
        kwargs...) where {D, T}

    rectangularinvariantmeasure(Array(hcat(data...,)), ϵ, estimator, kwargs...)
end



function summarise(invm::RectangularInvariantMeasure)
    D = size(invm.points, 1)
    npoints = size(invm.points, 2)

    points_str = "  points: $npoints $D-dimensional points\n"

    # Discretization
    ϵ = invm.ϵ
    ϵ_abs = invm.ϵ_absolute

    ϵ_str = "  ϵ: $ϵ\n"
    ϵ_abs_str = "  ϵ_absolute: $ϵ_abs\n"

    n_visited_bins = size(unique(invm.visited_bins_inds, dims = 2), 2)
    coord_minima = tuple(minimum(invm.visited_bins_coordinates, dims = 2)...,)
    coord_maxima = tuple(maximum(invm.visited_bins_coordinates, dims = 2)...,)

    inds_str = "  visited_bins_inds: $n_visited_bins unique bins (rectangular boxes) are visited by the points\n"
    coords_str = "  visited_bins_coords: Bins are distributed within the hypercube enclosing \n\tx_{min} =$coord_minima to \n\tx_{max} = $coord_maxima\n"
    bv = invm.binvisits
    binvisits_str = "  binvisits: $bv"

    TO = invm.transfermatrix
    iv = invm.measure
    transfermatrix_str = "  transfermatrix: $TO"
    measure_str = "  measure: $iv"
    return join(["RectangularInvariantMeasure\n", points_str, ϵ_str, ϵ_abs_str, inds_str, coords_str,
                binvisits_str, transfermatrix_str, measure_str])
end

Base.show(io::IO, invm::RectangularInvariantMeasure) = println(io, summarise(invm))
