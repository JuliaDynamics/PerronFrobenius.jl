
"""
    transferoperator_grid(E::AbstractEmbedding,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}};
        allocate_frac::Float64 = 1.0,
        boundary_condition = :exclude) ->
        RectangularBinningTransferOperator

Estimates the transfer operator for an embedding using the binning scheme
specified by `ϵ`:

- `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
- `ϵ::Float` divides each axis into intervals of size `ϵ`.
- `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
- `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.

`allocate_frac` controls what fraction of the total number of
possible transitions (``n_states^2``) we pre-allocate for. Typically,
we don't need more than 10 percent of the total number of possible
transitions.

`boundary_condition` controls what to do with the forward
map of the last point of the embedding. The default, `:exclude`,
simply ignores the last point.
"""
function transferoperator_grid(
        E::AbstractEmbedding,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}};
        allocate_frac::Float64 = 1.0,
        boundary_condition::Symbol = :none)


    valid_boundary_conditions = [:none, :exclude, :circular,
                                :invariantize]
    if !(boundary_condition ∈ valid_boundary_conditions)
        error("Boundary condition $boundary_condition not valid.")
    end

    if boundary_condition == :invariantize
        warn("Invariantizing embedding")
        E = invariantize(E)
    end

    # Identify which bins of the partition resulting from using ϵ each
    # point of the embedding visits.
    visited_bins = assign_bin_labels(E, ϵ)

    # Which are the visited bins, which points
    # visits which bin, repetitions, etc...
    binvisits = organize_bin_labels(visited_bins)

    # Use that information to estimate transfer operator
    transferoperator(binvisits, allocate_frac, boundary_condition)
end

"""
    transferoperator_grid(points::AbstractArray{T, 2},
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}},
        allocate_frac::Float64 = 1.0,
        boundary_condition = :exclude) ->
        RectangularBinningTransferOperator

Estimates the transfer operator for a set of points (assumed to be
provided as an array where each point is a column vector) with the
binning scheme specified by `ϵ`:

- `ϵ::Int` divides each axis into `ϵ` intervals of the same size.
- `ϵ::Float` divides each axis into intervals of size `ϵ`.
- `ϵ::Vector{Int}` divides the i-th axis into `ϵᵢ` intervals of the same size.
- `ϵ::Vector{Float64}` divides the i-th axis into intervals of size `ϵᵢ`.

`allocate_frac` controls what fraction of the total number of
possible transitions (``n_states^2``) we pre-allocate for. Typically,
we don't need more than 10 percent of the total number of possible
transitions.

`boundary_condition` controls what to do with the forward
map of the last point of the embedding. The default, `:exclude`,
simply ignores the last point.
"""
function transferoperator_grid(
        points::AbstractArray{T, 2},
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}};
        allocate_frac::Float64 = 1.0,
        boundary_condition = :exclude) where T

    # Identify which bins of the partition resulting from using ϵ each
    # point of the embedding visits.
    visited_bins = assign_bin_labels(points, ϵ)

    # Which are the visited bins, which points
    # visits which bin, repetitions, etc...
    binvisits = organize_bin_labels(visited_bins)

    # Use that information to estimate transfer operator
    transferoperator(binvisits, allocate_frac, boundary_condition)
end
