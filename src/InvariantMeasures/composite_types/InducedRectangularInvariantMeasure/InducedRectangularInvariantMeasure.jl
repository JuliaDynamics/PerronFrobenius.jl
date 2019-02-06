
import StaticArrays: SVector, MVector
import DelayEmbeddings: Dataset

include("induced_invariant_measure.jl")

export InducedRectangularInvariantMeasure, inducedrectangularinvariantmeasure

"""
    struct InducedRectangularInvariantMeasure{T} <: AbstractRectangularInvariantMeasure where {T}
        points::AbstractArray{T, 2}
        ϵF::Union{Int, Float64, Vector{Int}, Vector{Float64}}
        ϵj::Union{Int, Float64, Vector{Int}, Vector{Float64}}
        bins_ϵF::AbstractArray{Float64, 2}
        bins_visited_ϵj::AbstractArray{Float64, 2}
        transferoperator_ϵj::TransferOperatorRectangularBinning
        measure_ϵj::InvariantDistribution
        measure_induced::InvariantDistribution
    end

An induced RectangularInvariantMeasure. Created from a set of points by
discretizing the state space into rectangular bins with edge lengths dictated
by the binning scheme `ϵ`.

The invariant measure is first computed at resolution `ϵj`, then induced at the
final resolution `ϵF` by considering the overlap of the bins between resolutions.


## Fields
- **`points`**: The points for which to estimate the invariant measure.

- **`ϵF`**: The binning scheme at the final resolution, expressed as absolute
edge lengths. The following `ϵ` are valid, and will all be converted into
    the `ϵ::Vector{Float64` because knowing edge lengths at both resolutions
    is necessary.
        1. `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
        2. `ϵ::Float` divides each axis into intervals of size `ϵ`.
        3. `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
        4. `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.
- **`ϵF_absolute`**: `ϵF` converted to absolute edge lengths.

- **`ϵj`**: The binning scheme at the resolution from which we induce the measure.

- **`ϵj_absolute`**: `ϵj` converted to absolute edge lengths.

- **`bins_ϵF`**: The origins of all bins covering the hyperrectangle spanned
    by a box covering of `points` at resolution `ϵF`.

- **`bins_visited_ϵj`**: As for `bins_ϵF`, but at resolution `ϵj` and only
    for nonempty bins.

- **`transferoperator_ϵj`**: The transfer operator from from which the invariant
    measure is obtained at resolution `ϵj`.

- **`measure_ϵj`**: The invariant measure obtained from at resolution `ϵj`,
    which is obtained from `transferoperator_ϵj`.

- **`measure_induced`**: The invariant measure induced at resolution `ϵF`
    from resolution `ϵj`.
"""
struct InducedRectangularInvariantMeasure{T} <: AbstractRectangularInvariantMeasure where {T}
    points::AbstractArray{T, 2}
    ϵF::Union{Int, Float64, Vector{Int}, Vector{Float64}}
    ϵF_absolute::Vector{Float64}
    ϵj::Union{Int, Float64, Vector{Int}, Vector{Float64}}
    ϵj_absolute::Vector{Float64}
    bins_ϵF::AbstractArray{Float64, 2}
    bins_visited_ϵj::AbstractArray{Float64, 2}
    transferoperator_ϵj::TransferOperatorRectangularBinning
    measure_ϵj::InvariantDistribution
    measure_induced::InvariantDistribution
end


"""
    inducedrectangularinvariantmeasure(points, ϵF, ϵⱼ)

Compute the invariant measure the boxes of a hyperrectangular subdivision
of the hypercube spanned by `points` at resolution
`ϵF` induced by the invariant measure computed for the same points at resolution
`ϵⱼ`. Returns a `InducedRectangularInvariantMeasure` instance.
"""
function inducedrectangularinvariantmeasure(points::AbstractArray{T, 2}, ϵF, ϵⱼ) where {T}
    ϵF_absolute = minima_and_stepsizes(points, ϵF)[2]
    ϵⱼ_absolute = minima_and_stepsizes(points, ϵⱼ)[2]

    # All possible bins in the final partition given by ϵF
    allbins_ϵF = boxorigins(points, ϵF_absolute)

    # The unique bins visited by the orbit in the partition given by ϵⱼ
    visitedbins_ϵⱼ = unique(assign_coordinate_labels(points, ϵⱼ_absolute), dims = 2)

    # The transfer operator in the partition given by ϵⱼ
    TOϵⱼ = TransferOperatorEstimatorRectangularBinning(points, ϵⱼ_absolute)

    # The measure of the of the visited bins in the partition given by ϵⱼ
    μϵⱼ = invariantmeasure(TOϵⱼ)

    # The induced measure at resolution ϵF. This does not filter
    # zero entries, so be careful when comparing with the distribution
    # obtained using invariantmeasure(::TransferOperatorRectangularBinning)
    μϵF = μ_allbins_ϵF_induced_by_binningscheme_ϵⱼ(
        points,
        ϵF,
        allbins_ϵF,
        ϵⱼ,
        visitedbins_ϵⱼ,
        μϵⱼ.dist)

    InducedRectangularInvariantMeasure(
        points,
        ϵF,
        ϵF_absolute,
        ϵⱼ,
        ϵⱼ_absolute,
        allbins_ϵF,
        visitedbins_ϵⱼ,
        TOϵⱼ,
        μϵⱼ,
        InvariantDistribution(μϵF, findall(μϵF .> 0))
    )
end

function inducedrectangularinvariantmeasure(points::Vector{Vector{T}}, ϵF, ϵⱼ) where {T}
    inducedrectangularinvariantmeasure(hcat(points...,), ϵF, ϵⱼ)
end


function inducedrectangularinvariantmeasure(points::Vector{SVector{T}}, ϵF, ϵⱼ) where {T}
    inducedrectangularinvariantmeasure(Array(hcat(points...,)), ϵF, ϵⱼ)
end

function inducedrectangularinvariantmeasure(points::Vector{MVector{T}}, ϵF, ϵⱼ) where {T}
    inducedrectangularinvariantmeasure(Array(hcat(points...,)), ϵF, ϵⱼ)
end

function inducedrectangularinvariantmeasure(points::Dataset, ϵF, ϵⱼ) where {T}
    inducedrectangularinvariantmeasure(transpose(Matrix(points)), ϵF, ϵⱼ)
end


function summarise(invm::InducedRectangularInvariantMeasure)
    ϵF = invm.ϵF
    ϵⱼ = invm.ϵj
    ϵF_str = " partition resulting from binning scheme ϵF = $ϵF"
    ϵⱼ_str = " partition resulting from binning scheme ϵⱼ = $ϵⱼ"
    D = size(invm.points, 1)
    npts = size(invm.points, 2)

    measure_type = typeof(invm)
    measure_type_str = "$measure_type"
    pts_str = "$npts $D-dimensional points"
    measure_type_str*" from "*pts_str*" at"*ϵF_str*" induced by"*ϵⱼ_str

end

Base.show(io::IO, invm::InducedRectangularInvariantMeasure) = println(io, summarise(invm))
