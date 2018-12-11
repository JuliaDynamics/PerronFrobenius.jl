
include("RectangularInvariantMeasure.jl")
include("InducedRectangularInvariantMeasure.jl")
include("average_invariant_measure.jl")


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
        avg_measure::InvariantDistribution
    end

An average rectangular invariant measure induced at resolution
`ϵF` by the invariant measure computed at multiple resolutions `ϵs`.

## Fields

- **`points`**: The points for which to estimate the invariant measure.
- **`ϵs`**: The binning schemes. Must be an iterable of `ϵ`s of the following
    types:
    - `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
    - `ϵ::Float` divides each axis into intervals of size `ϵ`.
    - `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
    - `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.
- **`induced_measures`**: The measures induced by each of the resolutions in `ϵ`.
- **`avg_measure`**: The average invariant measure over the visited bins
    induced by the measures at resolutions `ϵs[1], ϵ[2], …, ϵ[end]`.
"""
struct AverageRectangularInvariantMeasure{T} <: AbstractRectangularInvariantMeasure where {T}
    points::AbstractArray{T, 2}
    ϵF::Union{Int, Float64, Vector{Int}, Vector{Float64}}
    ϵs
    induced_measures::Vector{PerronFrobenius.InducedRectangularInvariantMeasure{T}}
    avg_measure::PerronFrobenius.InvariantDistribution
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
- **`ϵF`**: The final resolution at which to induce the measure. Must
    - `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
    - `ϵ::Float` divides each axis into intervals of size `ϵ`.
    - `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
    - `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.
- **`ϵs`**: The multiple resolutions from which the measure at resolution `ϵF`
    is induced. Must be an iterable of valid `ϵ` (see the `ϵF` argument for
    details).
"""
function AverageRectangularInvariantMeasure(points, ϵF, ϵs)
    induced_measures = [InducedRectangularInvariantMeasure(points, ϵF, ϵᵢ) for ϵᵢ in ϵs]
    avg_measure = average_measure(induced_measures)
    AverageRectangularInvariantMeasure(
        points,
        ϵF,
        ϵs,
        induced_measures,
        avg_measure
    )
end

export AverageRectangularInvariantMeasure
