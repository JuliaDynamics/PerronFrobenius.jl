"""
    transferoperator(bv::BinVisits)

Estimate transfer operator from information about which bin
gets visited by which points of the orbit.

## Example usage

```
# Random orbit.
orbit = rand(1000, 3)
x, y, z = orbit[:, 1], orbit[:, 2], orbit[:, 3]

# Embedding E = {(y(t+1), y(t), x(t))}
E = embed([x, y], [2, 2, 1], [1, 0, 0]);

# Which bin sizes to use along each dimension?
ϵ = [.4, .2, .4]

# Identify which bins of the partition resulting
# from using ϵ each point of the embedding visits.
visited_bins = assign_bin_labels(E, ϵ)

# Which are the visited bins, which points
# visits which bin, repetitions, etc...
binvisits = organize_bin_labels(visited_bins)

# Use that information to estimate transfer operator
TO = transferoperator2(binvisits)

# Verify that TO is Markov (NB. last row might be zero!,
# then this test might fail)
all(sum(TO, 2) .≈ 1)
```
"""
function transferoperator(bv::BinVisits)
    first_visited_by = bv.first_visited_by
    repeated_visitors = bv.repeated_visitors
    visits_whichbin = bv.visits_whichbin

    # Initialise transfer (Perron-Frobenius) operator as a sparse matrix
    # (keep track of indices and values in separate columns for now)
    n_possible_nonzero_entries = length(first_visited_by)^2
    I = zeros(Int, n_possible_nonzero_entries)
    J = zeros(Int, n_possible_nonzero_entries)
    V = zeros(Float64, n_possible_nonzero_entries)

    max_row = 0

    # Preallocate target index for the case where there is only
    # one point of the orbit visiting a bin.
    target_bin_j::Int = 0
    n_visitsᵢ::Int = 0
    for i in 1:length(first_visited_by)
        # How many times is this bin visited?
        n_visitsᵢ = length(repeated_visitors[i])

        # If there is only one point and it happens to be the last,
        # we skip it, because we don't know its image. Otherwise,
        # we
        if n_visitsᵢ == 1
            if !(i == visits_whichbin[end])

                # Keep track of how many transitions we've considered.
                max_row += 1

                I[max_row] = i
                # If this part of the code gets visited, there is
                # just one point visiting the i-th bin. To which bin
                # does this point jump if we shift the point one
                # time step ahead in time along its orbit?
                target_bin_j = visits_whichbin[repeated_visitors[i] + 1][1]

                # We now know that exactly one point, the i-th, does the
                # transition from i to the target j.
                J[max_row] = target_bin_j
                V[max_row] = 1
            end

        # If more than one point visits the bin,
        elseif n_visitsᵢ > 1
            # identify find which points of the orbit visits the i-th bin.
            visiting_pts = repeated_visitors[i]

            # If the visited bin under consideration is the last,
            # we  exclude the last point that visits it - that point
            # will have nowhere to go along its orbit in the next
            # time step, because it is the last point in the orbit.
            if i == visits_whichbin[end]
                n_visitsᵢ = length(visiting_pts) - 1
                visiting_pts = visiting_pts[1:(end - 1)]
            end

            # To which boxes do each of the points visiting the i-th bin
            # jump in the next time step?
            target_bins_j = visits_whichbin[visiting_pts .+ 1]


            # How many unique target bins are there?
            unique_target_bins = unique(target_bins_j)

            # For each of the j unique target bins count how many
            # points jump from the i-th bin to that target bin.
            for j in 1:length(unique_target_bins)
                # How many points jump from the i-th bin to this j-th bin?
                n_transitions_i_to_j = sum(target_bins_j .== unique_target_bins[j])

                # Keep track of how many transitions we've considered.
                max_row += 1

                # We now know how many points jumps from the i-th bin
                # to the j-th bin. Normalize by the number of points
                # visiting the i-th bin
                I[max_row] = i
                J[max_row] = unique_target_bins[j]
                V[max_row] = n_transitions_i_to_j / n_visitsᵢ


            end
        end
    end

    # Combine indices and values into a sparse representation of the Perron-
    # Frobenius operator (transfer operator). Filter out the nonzero entries.
    # The number of rows in the matrix is given by the number of unique points
    # we have in the embedding.
    n_rows = length(first_visited_by)
    TO = sparse(I[1:max_row], J[1:max_row], V[1:max_row], n_rows, n_rows)
    RectangularBinningTransferOperator(Array(TO))
end

"""
 transferoperator(
        E::AbstractEmbedding,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}}
) -> RectangularBinningTransferOperator

Estimates the transfer operator given an embedding and ϵ,
the latter controlling the binning procedure. See docs
for `StateSpaceReconstruction.assign_bin_labels`.
"""
function transferoperator(
        E::AbstractEmbedding,
        ϵ::Union{Int, Float64, Vector{Int}, Vector{Float64}})

    # Identify which bins of the partition resulting from using ϵ each
    # point of the embedding visits.
    visited_bins = assign_bin_labels(E, ϵ)

    # Which are the visited bins, which points
    # visits which bin, repetitions, etc...
    binvisits = organize_bin_labels(visited_bins)

    # Use that information to estimate transfer operator
    transferoperator(binvisits)
end
