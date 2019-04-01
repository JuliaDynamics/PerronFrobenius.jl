import CausalityToolsBase: 
    BinningScheme,
    RectangularBinning

import ..BinningSchemes:
    InducedRectangularBinning,
    AverageRectangularBinning

"""
    invariantmeasure(data, binning_scheme::BinningScheme; kwargs...)

Partition `data` according to the the given `binning_scheme` and compute the invariant
measure over the partition elements.
"""
function invariantmeasure(data, binning_scheme::BinningScheme; kwargs...)

end


######################################
# Triangulation binnign schemes
######################################

"""
    invariantmeasure(pts, ϵ::TriangulationBinning, 
        simplex_intersection_type::ExactIntersection) -> TriangulationExactInvariantMeasure

Estimate the invariant measure over the state space defined by `pts` using a triangulation 
of the phase space as the partition, using exact simplex intersections to compute transition
probabilities between the states (simplices).

# Example

Assume we have sufficiently few points that a triangulation approach involving simplex 
intersections is computationally feasible to compute the invariant measure. Then 
a transfer operator can be computed as follows.

```julia 
pts = [rand(3) for i = 1:30]
invariantmeasure(pts, TriangulationBinning(), ExactIntersection())
```
"""
function invariantmeasure(pts, ϵ::TriangulationBinning, 
        simplex_intersection_type::ExactIntersection)
    triangulationexactinvariantmeasure(pts)
end


"""
    invariantmeasure(pts, ϵ::TriangulationBinning, 
        simplex_intersection_type::ApproximateIntersection;
        n::Int = 100, sample_randomly::Bool = false) -> TriangulationApproxInvariantMeasure

Estimate the invariant measure over the state space defined by `pts` using a triangulation 
of the phase space as the partition, using exact simplex intersections to compute transition
probabilities between the states (simplices).


# Example

Assume we have sufficiently few points that a triangulation approach involving simplex 
intersections is computationally feasible to compute the invariant measure. Then 
a transfer operator can be computed as follows.

```julia
pts = [rand(3) for i = 1:30]
invariantmeasure(pts, TriangulationBinning(), ApproximateIntersection())
```
"""
function invariantmeasure(pts, ϵ::TriangulationBinning, 
        simplex_intersection_type::ApproximateIntersection;
        n::Int = 100, sample_randomly::Bool = false)
    triangulationapproxinvariantmeasure(pts, n_sample_pts = n, sample_randomly = sample_randomly)
end

######################################
# Providing customized binning schemes
######################################
# Regular rectnagular invariant measure

""" 
    invariantmeasure(data, binning_scheme::RectangularBinning; kwargs...)

Partition `data` according to the the given `binning_scheme` and compute the invariant
measure over the partition elements. 


# Example

Assume we have enough points that a rectangular partition yields a good estimate of the 
invariant measure. Then the measure over the partition can be computed as follows.

```julia 
pts = [rand(3) for i = 1:2000]

# Use rectangular boxes constructed by subdividing each coordinate 
# axis into 10 subintervals of equal length.
binning_scheme = RectangularBinning(10)
invariantmeasure(pts, binning_scheme)
```
""" 
function invariantmeasure(data, binning_scheme::RectangularBinning; kwargs...)
    rectangularinvariantmeasure(
        data, 
        binning_scheme, 
        kwargs...)
end

# Induced measures
function invariantmeasure(data, binning_scheme::InducedRectangularBinning; kwargs...)
    inducedrectangularinvariantmeasure(
        data, 
        binning_scheme.ϵ_target_partition, 
        binning_scheme.ϵ_source_partition, 
        kwargs...)
end

# If you don't want to construct an induced invariant measure.
function invariantmeasure(data, 
    ϵ_target_partition::RectangularBinning,
    ϵ_source_partition::RectangularBinning; 
    kwargs...)

inducedrectangularinvariantmeasure(data, ϵ_target_partition.ϵ, ϵ_source_partition.ϵ, 
    kwargs...)
end

function invariantmeasure(data, binning_scheme::AverageRectangularBinning; kwargs...)
    averagerectangularinvariantmeasure(
        data, 
        binning_scheme.ϵ_target_partition, 
        binning_scheme.ϵs_source_partitions, 
        kwargs...)
end

function invariantmeasure(data, 
    ϵ_target_partition::RectangularBinning,
    ϵs_source_partitions::Vector{RectangularBinning}; 
    kwargs...)

    averagerectangularinvariantmeasure(
        data, 
        ϵ_target_partition.ϵ, 
        [ϵs_source_partitions[i].ϵ for i in 1:length(ϵs_source_partitions)], 
        kwargs...)
end




export invariantmeasure