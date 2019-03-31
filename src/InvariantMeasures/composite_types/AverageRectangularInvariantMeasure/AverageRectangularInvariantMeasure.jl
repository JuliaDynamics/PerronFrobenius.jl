import DelayEmbeddings:
    Dataset

import ..TransferOperators: get_binvisits
import CausalityToolsBase

import StaticArrays:
    SVector,
    MVector

include("average_invariant_measure.jl")

export AverageRectangularInvariantMeasure, averagerectangularinvariantmeasure

"""
    average_measure(induced_measures::Vector{InducedRectangularInvariantMeasure})

Compute the average measure from a set of precomputed induced measures.
Assumes that all induced measures were computed at the same final resolution
`ϵF`.
"""
function average_measure(induced_measures::Vector{InducedRectangularInvariantMeasure{T}}) where T
    n_measures = length(induced_measures)
    n_states = length(induced_measures[1].measure_induced.dist)
    measures = zeros(Float64, n_measures, n_states)

    for i in 1:n_measures
        measures[i, :] = induced_measures[i].measure_induced.dist
    end

    avg_measure = zeros(Float64, n_states)
    for i in 1:n_states
        avg_measure[i] = (1/n_measures) * sum(measures[:, i])
    end

    return InvariantDistribution(avg_measure, findall(avg_measure .> 0))
end

export average_measure


"""
    struct AverageRectangularInvariantMeasure{T} <: AbstractRectangularInvariantMeasure where {T}
        points::AbstractArray{T, 2}
        ϵF::Union{Int, Float64, Vector{Int}, Vector{Float64}}
        ϵs::Union{Int, Float64, Vector{Int}, Vector{Float64}}
        induced_measures::Vector{InducedRectangularInvariantMeasure}
        averaged_measure::RectangularInvariantMeasure
    end

An average rectangular invariant measure induced at resolution
`ϵF` by the invariant measure computed at multiple resolutions `ϵs`.

## Fields

- **`points`**: The points for which to estimate the invariant measure.

- **`ϵs`**: The binning schemes. Must be an iterable of `ϵ`s of the following
    types:
        1. `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
        2. `ϵ::Float` divides each axis into intervals of size `ϵ`.
        3. `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
        4. `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.

- **`induced_measures`**: The measures induced by each of the resolutions in `ϵ`.

- **`measure_ϵF`**: The invariant measure computed at resolution `ϵF`.

- **`averaged_measure_ϵF`**: The average rectangular invariant measure at resolution `ϵF`
    induced by the measures computed at resolutions `ϵs[1], ϵs[2], …, ϵs[end]`.
"""
struct AverageRectangularInvariantMeasure{T} <: AbstractRectangularInvariantMeasure where {T}
    points::AbstractArray{T, 2}
    ϵF::RectangularBinning
    ϵs::AbstractArray{RectangularBinning, 1}
    induced_measures::Vector{InducedRectangularInvariantMeasure{T}}
    measure_ϵF::RectangularInvariantMeasure
    measure_ϵF_averaged::RectangularInvariantMeasure
end

"""
    AverageRectangularInvariantMeasure(points, ϵF, ϵⱼ)

Compute the average invariant measure the boxes of a hyperrectangular subdivision
of the hypercube spanned by `points` at resolution
`ϵF` induced by the invariant measure computed for the same points at multiple
resolutions provided by `ϵs`.

This constructor returns an `AverageRectangularInvariantMeasure` instance.

## Arguments

- **`points`**: The points for which to estimate the invariant measure.

- **`ϵF`**: The final resolution at which to induce the measure. Must be a valid 
    `RectangularBinning` instances. See docs for `RectangularBinning` for details.

- **`ϵs`**: The multiple resolutions from which the measure at resolution `ϵF`
    is induced. Must be an iterable of valid `RectangularBinning` instances.
    See docs for `RectangularBinning` for details.
"""
function averagerectangularinvariantmeasure(points::AbstractArray{T, 2}, 
        ϵF::RectangularBinning, 
        ϵs::AbstractArray{RectangularBinning, 1};
        estimator::Symbol = :TransferOperatorEstimatorRectangularBinning,
        kwargs...) where {T}
    
    ptsvec = [points[:, i] for i = 1:size(points, 2)]
    edgelengths_ϵF = CausalityToolsBase.get_edgelengths(ptsvec, ϵF)
    edgelengths_ϵs = [CausalityToolsBase.get_edgelengths(ptsvec, ϵi) for (k, ϵi) in enumerate(ϵs)]

    induced_measures = [inducedrectangularinvariantmeasure(points, edgelengths_ϵF, edgelengths_ϵs) for ϵᵢ in ϵs]
    avg_measure = average_measure(induced_measures)

    # Check that measure sums to 1! With too few input points and too fine partitions, 
    # this might not be the case
    if sum(avg_measure.dist[isfinite.(avg_measure.dist)]) < 0.99 # some tolerance allowed
        @warn "The number of input points is too low for these resolutions. The average measure will not sum to 1."
    end

    measure_ϵF::RectangularInvariantMeasure = rectangularinvariantmeasure(points, ϵF)
    avg_measure_ϵF::RectangularInvariantMeasure = RectangularInvariantMeasure(
        measure_ϵF.points,
        measure_ϵF.binning_scheme,
        measure_ϵF.axisminima,
        measure_ϵF.edgelengths,
        hcat(measure_ϵF.encoded_points...,),
        measure_ϵF.visited_bins_coordinates,
        measure_ϵF.binvisits,
        measure_ϵF.transfermatrix,
        avg_measure
    )

    AverageRectangularInvariantMeasure(
        points,
        ϵF,#get_edgelengths(points, ϵF.ϵ),#
        ϵs,#[get_edgelengths(points, ϵs[i]) for i = 1:length(ϵs)],
        induced_measures,
        measure_ϵF,
        avg_measure_ϵF
    )
end

averagerectangularinvariantmeasure(points::AbstractArray{T, 2}, ϵF, ϵs) where T = 
    averagerectangularinvariantmeasure(
        points, 
        RectangularBinning(ϵF), 
        [RectangularBinning(ϵi) for (i, ϵi) in enumerate(ϵs)]
        )

function averagerectangularinvariantmeasure(points::Vector{Vector{T}}, ϵF, ϵs) where {T}
    averagerectangularinvariantmeasure(hcat(points...,), ϵF, ϵs)
end

function averagerectangularinvariantmeasure(points::Vector{SVector{T}}, ϵF, ϵs) where {T}
    averagerectangularinvariantmeasure(Array(hcat(points...,)), ϵF, ϵs)
end

function averagerectangularinvariantmeasure(points::Vector{MVector{T}}, ϵF, ϵs) where {T}
    averagerectangularinvariantmeasure(Array(hcat(points...,)), ϵF, ϵs)
end

function averagerectangularinvariantmeasure(points::Dataset, ϵF, ϵs) where {T}
    averagerectangularinvariantmeasure(transpose(Matrix(points)), ϵF, ϵs)
end



function summarise(invm::AverageRectangularInvariantMeasure)
    ϵF = invm.ϵF
    ϵs = invm.ϵs
    ϵF_str = " partition resulting from binning scheme  ϵF = $ϵF"
    ϵs_str = " partitions resulting from binning schemes ϵs = $ϵs"
    D = size(invm.points, 1)
    npts = size(invm.points, 2)

    measure_type = typeof(invm)

    pts_str = "$npts $D-dimensional points\n"
    return join([measure_type, " from ", pts_str, " AVERAGED at", ϵF_str, 
                "\n INDUCED by", ϵs_str])

end

Base.show(io::IO, invm::AverageRectangularInvariantMeasure) = println(io, summarise(invm))
