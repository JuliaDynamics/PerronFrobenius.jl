import DelayEmbeddings:
    Dataset

import StateSpaceReconstruction: 
    Embeddings,
    minima_and_stepsizes
    
import StaticArrays:
    SVector,
    MVector

import ..TransferOperators:
    TransferOperatorEstimatorRectangularBinning

import Base.Cartesian

"""
    bin_origins_subdivided_boxcovering(axisminima, stepsizes,
            n_intervals_eachaxis) -> Array{Float64, 2}

Consider a (single) hyperrectangular box covering a point cloud. We're given the
minimal along each coordinate axis (`axisminima`), a set of `stepsizes` along
each axis (which might be different for each axis), and a set of
intervals ()`n_intervals_per_axis`) to which divide each axis into.

Together, this information gives a complete box covering of the points. This
function returns the origins of all the boxes that the original (single) big
hyperrectangle gets subdivided into, given the `stepsizes` and
`n_intervals_per_axis`.

The origins are returned as an array where each column is the origin of one
box.
"""
function bin_origins_subdivided_boxcovering(axisminima, stepsizes, n_intervals_eachaxis)
    D = length(axisminima)
    coords = zeros(Float64, D, prod(n_intervals_eachaxis))

    n_origins_generated = 0
    for I in Iterators.product([1:i for i in n_intervals_eachaxis]...,)
        n_origins_generated += 1
        for dim in 1:D
            coords[dim, n_origins_generated] = axisminima[dim] + stepsizes[dim]*(I[dim] - 1)
        end
    end
    coords
end

export bin_origins_subdivided_boxcovering



"""
    boxorigins(::StateSpaceReconstruction.AbstractEmbedding, ϵ) -> Array{Float64, 2}

Consider the box covering of the embedding induced by the binning scheme `ϵ`.
This functions returns the origin of each of those bins (always counting
from minima to maxima along each axis, with boxes extending 1% past the extrema 
in both directions along each axis).
"""
function boxorigins end 

function boxorigins(points::AbstractArray{T, 2}, ϵ) where T
    # Make sure that the array contains points as columns.
    if size(points, 1) > size(points, 2)
        throw(DomainError("The dimension of the dataset exceeds the number of points."))
    end

    D = size(points, 1)
    n_pts = size(points, 2)

    axisminima = zeros(Float64, D)
    top = zeros(Float64, D)

    for i = 1:D
        axisminima[i] = minimum(points[i, :])
        top[i] = maximum(points[i, :])
    end
    axisminima = axisminima - (top - axisminima) / 100
    top = top + (top - axisminima) / 100

    stepsizes = Vector{Float64}(undef, D)
    if typeof(ϵ) <: Float64
        stepsizes = [ϵ for i in 1:D]
    elseif typeof(ϵ) == Vector{Float64}
        stepsizes = ϵ
    elseif typeof(ϵ) <: Int
        stepsizes = (top - axisminima) / ϵ
    elseif typeof(ϵ) == Vector{Int}
        stepsizes = (top - axisminima) ./ ϵ
    end
    n_intervals_eachaxis = floor.(Int, (top .- axisminima) ./ stepsizes) .+ 1

    bin_origins_subdivided_boxcovering(axisminima,
        stepsizes,
        n_intervals_eachaxis)
end


function boxorigins(E::Embeddings.AbstractEmbedding, ϵ)
    boxorigins(E.points, ϵ)
end

function boxorigins(D::Dataset, ϵ)
    boxorigins(transpose(Matrix(D)), ϵ)
end

function boxorigins(pts::Vector{Vector{T}}, ϵ) where T
    c = boxorigins(hcat(pts...,), ϵ)
    [c[:, i] for i = 1:maximum(size(c))]
end

function boxorigins(pts::Vector{SVector{D, T}}, ϵ) where {D, T}
    c = boxorigins(Array(hcat(pts...,)), ϵ)
    [c[:, i] for i = 1:maximum(size(c))]
end

function boxorigins(pts::Vector{MVector{D, T}}, ϵ) where {D, T}
    c = boxorigins(Array(hcat(pts...,)), ϵ)
    [c[:, i] for i = 1:maximum(size(c))]
end

export boxorigins



######################################################
# Some function barriers to reduce memory allocations
######################################################
function fill_edges!(edges, a::Float64, b::Float64, c::Float64, d::Float64)
    edges[1] = a
    edges[2] = b
    edges[3] = c
    edges[4] = d
end

function find_lengths_intersections!(sorted_edges, edges, permutation,
        temp, tmp_reset, T, intersecting_lengths, c)

    # partialsort reduces allocations by half
    partialsortperm!(permutation, edges, 1:2, rev = true)
    sort!(edges, rev = true)

    T .= temp[permutation]

    if T[1] != T[2]
        intersecting_lengths[c] = edges[2] - edges[3]
    else
        intersecting_lengths[c] = 0.0
    end

    temp .= tmp_reset
end


######################################################
"""
    μ_of_bₐ_induced_by_bins_from_binningscheme_ϵⱼ(
        ϵF::Vector{Float64},
        bin_origin_bₐ::Vector{Float64},
        ϵⱼ::Vector{Float64},
        binorigins_visitedbins_ϵⱼ::Array{Float64, 2},
        μϵⱼ::Vector{Float64}) -> Float64


## Returns

The measure ``\\mu(b_a)`` of the bin ``b_a`` in the partition
given by the binning scheme `ϵF` induced by the elements of another
partition given by the binning scheme `ϵⱼ` overlapping with ``b_a``.

## Explanation

Consider a rectangular binning ``B_{\\epsilon_{F}}`` of some points, given by the
partition scheme `ϵF`. Next, consider a separate rectanglar binning
``B_{\\epsilon_{j}}`` of the same points, but formed from the partition scheme `ϵⱼ`.

Both partitions are now made up of rectangular boxes, but differ in that the
boxes have different sizes. Thus, an invariant measure computed over
``B_{\\epsilon_{F}}`` will in general be different from an invariant measure
computed over ``B_{\\epsilon_{j}}`` (the first partition may be coarse than
the second, for example). However, for a given box ``b_a \\in B_{\\epsilon_{F}}``,
there will always be some boxes ``b_k \\in B_{\\epsilon_{j}}`` overlapping with it.
Each of the ``b_k \\in B_{\\epsilon_{j}}`` have a measure ``\\mu_k`` associated
with it. We want to know what the measure of ``b_a`` would be from the perspective
of ``B_{\\epsilon_{j}}``.

This function finds the measure of the target bin ``b_a \\in B_{\\epsilon_{F}}``
induced by the bins ``\\{b_k \\in B_{\\epsilon_{j}} : b_k \\cap b_a \\neq \\emptyset \\}``
(that is, by the ``b_k`` having have some intersection with ``b_a``).


## Arguments

- `ϵF`: Edge lengths of boxes in the ``B_{\\epsilon_{F}}`` rectangular partition.
- `bin_origin_bₐ`: Origin of the a-th bin in ``B_{\\epsilon_{F}}``.
- `ϵⱼ`: Edge lengths of boxes in the ``B_{\\epsilon_{j}}`` rectangular partition.
- `binorigins_visitedbins_ϵⱼ`: The origin of each visited bin in
    ``P_{\\epsilon_{j}}``. Array where columns are the coordinates of the
    bin origins.
- `μϵⱼ`: the measures of the visited bins in ``B_{\\epsilon_{j}}``.
"""
function μ_of_bₐ_induced_by_bins_from_binningscheme_ϵⱼ(
            ϵF::Vector{Float64},
            bin_origin_bₐ::Vector{Float64},
            ϵⱼ::Vector{Float64},
            binorigins_visitedbins_ϵⱼ::Array{Float64, 2},
            μϵⱼ::Vector{Float64})

    dim = length(ϵF)
    volume_ϵⱼ = prod(ϵⱼ)
    n_visited_bins_ϵⱼ = length(μϵⱼ)
    binning_target_bin1 = falses(n_visited_bins_ϵⱼ)
    binning_target_bin2 = zeros(n_visited_bins_ϵⱼ)
    edges = zeros(Float64, 4)
    sorted_edges = zeros(Float64, 4)
    permutation = zeros(Int, 4)
    temp::Vector{Int} = [1, 1, 2, 2]
    tmp_reset::Vector{Int} = [1, 1, 2, 2]
    T::Vector{Number} = zeros(Number, size(temp))
    intersecting_lengths = zeros(Float64, dim)

    for b = 1:n_visited_bins_ϵⱼ
        for c = 1:dim
            fill_edges!(edges,
               bin_origin_bₐ[c] + ϵF[c],
               bin_origin_bₐ[c],
               binorigins_visitedbins_ϵⱼ[c, b] + ϵⱼ[c],
               binorigins_visitedbins_ϵⱼ[c, b])

            find_lengths_intersections!(
                sorted_edges, edges, permutation,
                temp, tmp_reset, T,
                intersecting_lengths, c)
        end

        intersecting_volume = prod(intersecting_lengths)
        if intersecting_volume > 0
            binning_target_bin1[b] = true
            binning_target_bin2[b] = intersecting_volume
        end
    end
    ρₐ = transpose(μϵⱼ[binning_target_bin1]) *
        binning_target_bin2[binning_target_bin1]/volume_ϵⱼ
end


export μ_of_bₐ_induced_by_bins_from_binningscheme_ϵⱼ

"""
    μ_allbins_ϵF_induced_by_binningscheme_ϵⱼ(ϵF,
        binorigins_ϵF::Array{Float64, 2},
        ϵⱼ,
        binorigins_visitedbins_ϵⱼ::Array{Float64, 2},
        μϵⱼ::Vector{Float64})  -> Vector{Float64}

## Returns
The measure ``\\mu(b_a)`` of the bin ``b_a`` in the partition
given by the binning scheme `ϵF` induced by the elements of another
partition given by the binning scheme `ϵⱼ` overlapping with ``b_a``.


## Explanation
Consider a rectangular binning ``B_{\\epsilon_{F}}`` of some points, given by the
partition scheme `ϵF`. Next, consider a separate rectanglar binning
``B_{\\epsilon_{j}}`` of the same points, but formed from the partition scheme `ϵⱼ`.

Both partitions are now made up of rectangular boxes, but differ in that the
boxes have different sizes. Thus, an invariant measure computed over
``B_{\\epsilon_{F}}`` will in general be different from an invariant measure
computed over ``B_{\\epsilon_{j}}`` (the first partition may be coarse than
the second, for example). However, for a given box ``b_a \\in B_{\\epsilon_{F}}``,
there will always be some boxes ``b_k \\in B_{\\epsilon_{j}}`` overlapping with it.
Each of the ``b_k \\in B_{\\epsilon_{j}}`` have a measure ``\\mu_k`` associated
with it. We want to know what the measure of ``b_a`` would be from the perspective
of ``B_{\\epsilon_{j}}``.

This function finds the measure of the target bin ``b_a \\in B_{\\epsilon_{F}}``
induced by the bins ``\\{b_k \\in B_{\\epsilon_{j}} : b_k \\cap b_a \\neq \\emptyset \\}``
(that is, by the ``b_k`` having have some intersection with ``b_a``).


## Arguments

- `ϵF`: Edge lengths along each coordinate axis in ``B_{\\epsilon_{F}}``.
- `binorigins_ϵF`: Coordinates of all bins ``b_i \\in B_{\\epsilon_{F}}``.
- `ϵⱼ`: Edge lengths along each coordinate axis in ``B_{\\epsilon_{j}}``.
- `binorigins_visitedbins_ϵⱼ`: Coordinates of the origin of each visited bin
    in ``P_{\\epsilon_{j}}``. Each column of this array is a bin origin.
- `μϵⱼ`: the measures of the visited bins in ``B_{\\epsilon_{j}}``.
"""
function μ_allbins_ϵF_induced_by_binningscheme_ϵⱼ(
            points,
            ϵF,
            binorigins_ϵF::Array{Float64, 2},
            ϵⱼ,
            binorigins_visitedbins_ϵⱼ::Array{Float64, 2},
            μϵⱼ::Vector{Float64})

    ϵF = minima_and_stepsizes(points, ϵF)[2]
    ϵⱼ = minima_and_stepsizes(points, ϵⱼ)[2]
    n = size(binorigins_ϵF, 2)
    ρ = zeros(Float64, n)

    # Loop over the bins bᵢ in the partition induced by the binning scheme ϵF,
    # and find the measure of those bins in terms of the measure of the
    # bins of the other partition.
    for i = 1:n
        bₐ = binorigins_ϵF[:, i]
        ρ[i] = μ_of_bₐ_induced_by_bins_from_binningscheme_ϵⱼ(
                    ϵF,
                    bₐ,
                    ϵⱼ,
                    binorigins_visitedbins_ϵⱼ, # only need to check the visited bins in the other partition
                    μϵⱼ)
    end

    return ρ
end

export μ_allbins_ϵF_induced_by_binningscheme_ϵⱼ


"""
    μϵF_induced_by_ϵj(points, ϵF, ϵⱼ)

## Returns
The measures ``\\{\\mu(b_f)\\}`` of the bins ``b_f \\in B_{\\epsilon_{F}}``
(``B_{\\epsilon_{F}}`` is a partition given by the binning scheme `ϵF`)
induced by the bins ``b_j \\in B_{\\epsilon_{j}}`` (``B_{\\epsilon_{j}}``
is a partition given by the binning scheme `ϵⱼ`) which have different
measures ``\\{\\mu(b_j)\\}`` associated with them. The partitions may
have a different number of elements, and
``\\{\\mu(b_f)\\} = \\{\\mu(b_j)\\}`` only if
``B_{\\epsilon_{F}} = B_{\\epsilon_{j}}``.

This gives us a way of comparing the measures obtained at different
resolutions of a system.


## Explanation

Consider a cloud of `points` and two separate binning schemes, `ϵF` and `ϵⱼ`.
The binning schemes dictate how to create a rectangular partition covering the
points.

Create the partition ``B_{\\epsilon_{F}}`` from `ϵF` and a separate partition
``B_{\\epsilon_{j}}`` from `ϵⱼ`. Then, compute an invariant measure over each
partition separately.

The boxes in ``B_{\\epsilon_{F}}`` are have sizes different from
the boxes in ``B_{\\epsilon_{j}}``. Thus, the invariant measure associated
with each of the partitions convey something about the dynamics at different
scales. We may want to know how the dynamics exposed by
``B_{\\epsilon_{j}}`` is expressed in terms of another coarse-graining
``B_{\\epsilon_{F}}``.


This function computes the measure of the bins ``b_f \\in B_{\\epsilon_{F}}``
induced by the bins ``b_j \\in B_{\\epsilon_{j}}``. That is, for each ``b_f``,
find the set of bins
``\\{b_k \\in B_{\\epsilon_{j}} : b_k \\cap b_f \\neq \\emptyset \\}``.

Taking this further, we can decide on a final rectangular grid (resolution)
and compute the measure over that grid as an average over multiple separate
partitions (finer or coarser, depending on the specific application).

## Arguments
- **`points`**: The points representing observations from the system from
    which to compute an invariant measure.
- **`ϵF`**: Edge lengths along each coordinate axis in the partition ``B_{\\epsilon_{F}}``.
- **`ϵⱼ`**: Edge lengths along each coordinate axis in the partition ``B_{\\epsilon_{j}}``.
"""
function μϵF_induced_by_ϵj(points, ϵF, ϵⱼ)
    ϵF = minima_and_stepsizes(points, ϵF)[2]
    ϵⱼ = minima_and_stepsizes(points, ϵⱼ)[2]

    # All possible bins in the final partition given by ϵF
    allbins_ϵF = boxorigins(points, ϵF)

    # The unique bins visited by the orbit in the partition given by ϵⱼ
    visitedbins_ϵⱼ = unique(assign_coordinate_labels(points, ϵⱼ), dims = 2)

    # The transfer operator in the partition given by ϵⱼ
    TOϵⱼ = TransferOperatorEstimatorRectangularBinning(points, ϵⱼ)

    # The measure of the of the visited bins in the partition given by ϵⱼ
    μϵⱼ = invariantmeasure(TOϵⱼ).dist

    # The induced measure at resolution ϵF. This does not filter
    # zero entries, so be careful when comparing with the distribution
    # obtained using invariantmeasure(::TransferOperatorRectangularBinning)
    μϵF = μ_allbins_ϵF_induced_by_binningscheme_ϵⱼ(
        points,
        ϵF,
        allbins_ϵF,
        ϵⱼ,
        visitedbins_ϵⱼ,
        μϵⱼ)
end


μϵF_induced_by_ϵj(E::Embeddings.AbstractEmbedding, ϵF, ϵⱼ) =
    μϵF_induced_by_ϵj(E.points, ϵF, ϵⱼ)

export μϵF_induced_by_ϵj


"""
    induced_measure(points, ϵF, ϵⱼ)

Express the invariant measure of a rectangular box covering of `points`, at
the resolution given by the binning scheme `ϵⱼ` as the invariant measure at
a separate resolution given by the binning scheme `ϵF`.

See the documentation for [`μϵF_induced_by_ϵⱼ`](@ref) for more details.

## Arguments
- **`points`**: The points representing observations from the system from
    which to compute an invariant measure.
- **`ϵF`**: Edge lengths along each coordinate axis in the partition ``B_{\\epsilon_{F}}``.
- **`ϵⱼ`**: Edge lengths along each coordinate axis in the partition ``B_{\\epsilon_{j}}``.
"""
function induced_measure(points, ϵF, ϵⱼ)
    μϵF_induced_by_ϵj(points,  ϵF, ϵⱼ)
end

export induced_measure
