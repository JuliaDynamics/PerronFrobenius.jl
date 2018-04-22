function indexin_rows(A1::Array{Int, 2}, A2::Array{Int, 2})
    inds = []
    for j = 1:size(A1, 1)
        for i = 1:size(A2, 1)
            if all(A1[j, :] .== A2[i, :])
                push!(inds, i)
            end
        end
    end
    return inds
end

function indexin_rows(A1::Array{Float64, 2}, A2::Array{Float64, 2})
    inds = []
    for j = 1:size(A1, 1)
        for i = 1:size(A2, 1)
            if all(A1[j, :] .== A2[i, :])
                push!(inds, i)
            end
        end
    end
    return inds
end

function unique_rows_info(embedding)
    first_inds = GroupSlices.firstinds(GroupSlices.groupslices(embedding, 1))
    group_inds = GroupSlices.groupinds(GroupSlices.groupslices(embedding, 1))
    all_inds = indexin_rows(embedding, unique(embedding, 1))
    first_inds, group_inds, Int.(all_inds)
end



"""
    `TO_from_binning(embedding::AbstractArray{Float64, 2})`

Estimate transfer operator from a binning of the embedding, which is a
`dim`-by-`n_points` array of points forming the embedding.

"""
function TO_from_binning(embedding::Embedding, n_bins::Int)
    emb = embedding.embedding
    dim, n_pts = size(emb, 2), size(emb, 1)

    bottom = [minimum(emb[:, i]) for i in 1:dim]
    top = [maximum(emb[:, i]) for i in 1:dim]


    bottom = bottom - (top - bottom) / 100
    top = top + (top - bottom) / 100

    stepsizes = (top - bottom) / n_bins

    # Indices of the bins. The coordinates of each point of the original
    # embedding are assigned an integer number indicating which bin along
    # the respective dimension it falls into.
    binning = zeros(Int, n_pts, dim)

    for i = 1:n_pts
        binning[i, :] = ceil((emb[i, :] - bottom) ./ stepsizes)
    end

    first_inds, group_inds, all_inds = unique_rows_info(binning)

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
    Array(TO)
end

TO_from_binning(embedding::AbstractArray{Float64, 2}, n_bins::Int) = TO_from_binning(Embedding(embedding), n_bins)

"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
function coords_bottom(b::RectangularBinning, i::Int)
    b.bottom .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end

"""
Find the coordinates of the bottom/origin of the i-th nonempty bin.
"""
coords_origin(b::RectangularBinning, i::Int) = coords_bottom(b, i)

"""
Find the coordinates of the top of the i-th nonempty bin.
"""
function coords_top(b::RectangularBinning, i::Int)
    b.top .+ b.inds_nonempty_bins[i, :] .* b.stepsizes
end
