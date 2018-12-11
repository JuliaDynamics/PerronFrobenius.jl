include("induced_invariant_measure.jl")
include("InvariantDistribution.jl")

"""
    struct InducedRectangularInvariantMeasure{T} <: AbstractRectangularInvariantMeasure where {T}
        points::AbstractArray{T, 2}
        ϵF::Union{Int, Float64, Vector{Int}, Vector{Float64}}
        ϵj::Union{Int, Float64, Vector{Int}, Vector{Float64}}
        bins_ϵF::AbstractArray{Float64, 2}
        bins_visited_ϵj::AbstractArray{Float64, 2}
        transferoperator_ϵj::RectangularBinningTransferOperator
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
    - `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
    - `ϵ::Float` divides each axis into intervals of size `ϵ`.
    - `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
    - `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.
- **`ϵj`**: The binning scheme at the resolution from which we induce the measure.
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
    ϵj::Union{Int, Float64, Vector{Int}, Vector{Float64}}
    bins_ϵF::AbstractArray{Float64, 2}
    bins_visited_ϵj::AbstractArray{Float64, 2}
    transferoperator_ϵj::RectangularBinningTransferOperator
    measure_ϵj::InvariantDistribution
    measure_induced::InvariantDistribution
end


"""
    InducedRectangularInvariantMeasure(points, ϵF, ϵⱼ)

Compute the invariant measure the boxes of a hyperrectangular subdivision
of the hypercube spanned by `points` at resolution
`ϵF` induced by the invariant measure computed for the same points at resolution
`ϵⱼ`. This constructor returns an `InducedRectangularInvariantMeasure` instance.
"""
function InducedRectangularInvariantMeasure(points, ϵF, ϵⱼ)
    ϵF = minima_and_stepsizes(points, ϵF)[2]
    ϵⱼ = minima_and_stepsizes(points, ϵⱼ)[2]

    # All possible bins in the final partition given by ϵF
    allbins_ϵF = coords_bin_origins(points, ϵF)

    # The unique bins visited by the orbit in the partition given by ϵⱼ
    visitedbins_ϵⱼ = unique(assign_coordinate_labels(points, ϵⱼ), dims = 2)

    # The transfer operator in the partition given by ϵⱼ
    TOϵⱼ = transferoperator_grid(points, ϵⱼ)

    # The measure of the of the visited bins in the partition given by ϵⱼ
    μϵⱼ = left_eigenvector(TOϵⱼ)

    # The induced measure at resolution ϵF. This does not filter
    # zero entries, so be careful when comparing with the distribution
    # obtained using left_eigenvector(::RectangularBinningTransferOperator)
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
        ϵⱼ,
        allbins_ϵF,
        visitedbins_ϵⱼ,
        TOϵⱼ,
        μϵⱼ,
        InvariantDistribution(μϵF, findall(μϵF .> 0))
    )
end

export InducedRectangularInvariantMeasure
