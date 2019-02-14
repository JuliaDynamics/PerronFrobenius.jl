"""
    μϵF_induced_by_ϵjs(pts, ϵF, ϵjs)

## Returns
The measures ``\\{\\mu(b_f)\\}`` of the bins ``b_f \\in B_{\\epsilon_{F}}``
(``B_{\\epsilon_{F}}`` is a partition given by the binning scheme `ϵF`)
induced by the measures of the bins ``b_j \\in B_{\\epsilon_{j}}``
(``B_{\\epsilon_{j}}`` is a partition given by the binning scheme `ϵⱼ`)
for several different ``B_{\\epsilon_{j}}`` constructed from the
binning schemes `ϵⱼs`

This is essentially averaging the invariant measures obtained at different
resolutions, expressing the average at the final resolution given by
the binning scheme `ϵF`.


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

Extending this to multiple partitions (many different resolutions),
we can decide on a final rectangular grid (resolution) and compute the measure
over that grid as an average over multiple separate partitions (finer, coarser, or both,
depending on the specific application).

## Arguments
- **`points`**: The points representing observations from the system from
    which to compute an invariant measure.
- **`ϵF`**: Edge lengths along each coordinate axis in the partition ``B_{\\epsilon_{F}}``.
- **`ϵⱼ`**: Edge lengths along each coordinate axis in the partition ``B_{\\epsilon_{j}}``.
"""
function μϵF_induced_by_ϵjs(points, ϵF, ϵjs)
    n_partitions = length(ϵjs)

    μ_across_binsizes = Vector{Vector{Float64}}(undef, n_partitions)

    for i=1:n_partitions
        μ_across_binsizes[i] = μϵF_induced_by_ϵj(points, ϵF, ϵjs[i])
    end

    # Gather in a matrix and transpose. Now the i-th column
    # and a-th row give the measure indiced for the i-th
    # state of the target partition induced by the a-th
    # bin size ϵₐ.
    μ = transpose(hcat(μ_across_binsizes...,))
end


export μϵF_induced_by_ϵjs

"""
    average_measure(points, ϵF, ϵjs)

Express the invariant measure of a rectangular box covering of `points`, at
the resolution given by the binning scheme `ϵⱼ` as the invariant measure at
a separate resolution given by the binning scheme `ϵF`.

See the documentation for [`μϵF_induced_by_ϵjs`](@ref) for more details.

## Arguments
- **`points`**: The points representing observations from the system from
    which to compute an invariant measure.
- **`ϵF`**: Edge lengths along each coordinate axis in the partition ``B_{\\epsilon_{F}}``.
- **`ϵjs`**: Edge lengths along each coordinate axis in the (multiple) partitions ``B_{\\epsilon_{j}}``.
"""
function average_measure(points, ϵF, ϵjs)
    n_partitions = length(ϵjs)

    measures = μϵF_induced_by_ϵjs(points, ϵF, ϵjs)

    # Compute the average measure for the final partition
    n_states_μF = size(measures, 2)
    avg_measure = zeros(Float64, n_states_μF)
    for i in 1:n_states_μF
        avg_measure[i] = (1/n_partitions) * sum(measures[:, i])
    end

    return avg_measure
end

export average_measure
