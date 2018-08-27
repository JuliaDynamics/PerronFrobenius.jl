"""
    transferoperator(eb::RectangularBinning)

Estimate transfer operator from an equidistant, rectangular binning.
"""
function transferoperator(rb::RectangularBinning)
    @unpack dim, n_pts, bottom, top,
        stepsizes, inds_nonempty_bins, first_inds, group_inds, all_inds = rb

    # Initialise transfer (Perron-Frobenius) operator as a sparse matrix
    # (keep track of indices and values in separate columns for now)
    n_possible_nonzero_entries = length(first_inds)^2
    I = zeros(Int, n_possible_nonzero_entries)
    J = zeros(Int, n_possible_nonzero_entries)
    V = zeros(Float64, n_possible_nonzero_entries)

    # For all unique points
    max_row = 0
    for i in 1:length(first_inds)
        # How many times does this point appear?
        n_appearancesᵢ = length(group_inds[i])

        # If there is only one point and it happens to be the last, we skip it,
        # because we don't know its image.
        if n_appearancesᵢ == 1
            if !(i == all_inds[end])
                max_row += 1

                I[max_row] = i
                J[max_row] = all_inds[group_inds[i] + 1][1]
                V[max_row] = 1
            end

        # Otherwise, we have at least two points
        else
            # Find the indices of these (at least two) points
            indices = group_inds[i]

            # If one of the points is the last one, exclude it
            if i == all_inds[end]
                n_appearancesᵢ = length(indices) - 1
                indices = indices[1:(end - 1)]
            end

            # Find the (row) indices of the image points
            bin_image_indices = all_inds[indices .+ 1]

            # Go over the unique image points, and find the number of times
            # each of these unique image points are repeated.
            unique_indices = unique(bin_image_indices)
            for j in 1:length(unique_indices)
                num_of_i_in_j = sum(bin_image_indices .== unique_indices[j])
                max_row += 1

                I[max_row] = i
                J[max_row] = unique_indices[j]
                V[max_row] = num_of_i_in_j / n_appearancesᵢ
            end
        end
    end

    # Combine indices and values into a sparse representation of the Perron-
    # Frobenius operator (transfer operator). Filter out the nonzero entries.
    # The number of rows in the matrix is given by the number of unique points
    # we have in the embedding.
    n_rows = length(first_inds)
    TO = sparse(I[1:max_row], J[1:max_row], V[1:max_row], n_rows, n_rows)
    RectangularBinningTransferOperator(Array(TO))
end
