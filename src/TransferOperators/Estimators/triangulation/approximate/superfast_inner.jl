
"""

## Arguments
- **`pt::Vector{Float64}`**: The point
- **`toprow_simplex::Array{Float64, 2}`**: A matrix representation of the simplex where each
    column represents a vertex, but with a row of ones appended on top of it, so that the
    total dimension is `(dim + 1)-by-(dim + 1)`. This is needed when computing determinants.
- **`dim`**: The dimension of the space.
"""
function iscontained_superoptim!(pt, toprow_simplex, dim, signs, startstops)

    # Replace first vertex with the point
    #fill_vert

    found = true
end

"""
Fill the 2nd to last row in the i-th column of `prealloc_matrix` with the
2nd to last row in the `i`-th column of `toprow_simplex` (replacing the `i`th vertex).
"""
function fill_vertex(prealloc_matrix::Array{Float64, 2}, toprow_simplex::Array{Float64, 2},
        startstops::Vector{UnitRange{Int64}}, i::Int)

    prealloc_matrix[startstops[i]] .= toprow_simplex[startstops[i]]
end


function contained_superoptim!(signs, s_arr, sx, point, dim, temp_arr, starts, stops)
    # Redefine the temporary simplex. This is in-place, so we don't allocate
    # memory. We could also have re-initialised `signs`, but since we're never
    # comparing more than two consecutive signs, this is not necessary.
    rezero!(s_arr)
    rezero!(signs)
    fill_into!(s_arr, sx)

    #Replace first vertex with the point
    fill_at_inds!(s_arr, point, 1:dim)

    # Signed volume
    signs[1] = sign(mdet_fm_in_x_view!(temp_arr, s_arr, dim + 1, starts, stops))

    rezero!(s_arr)
    fill_into!(s_arr, sx) #reset

    for κ = 2:dim # Check remaining signs and stop if sign changes
        # Replace the ith vertex with the point we're cheking (leaving the
        # 1 appended to Vi intact.)
        idxs = ((dim + 1)*(κ - 1)+1):((dim + 1)*(κ - 1)+ 1 + dim - 1)
        fill_at_inds!(s_arr, point, idxs) # ith change

        signs[κ] = sign(mdet_fm_in_x_view!(temp_arr, s_arr, dim + 1, starts, stops))

        if !(signs[κ-1] == signs[κ])
            return false
        end

        rezero!(s_arr)
        fill_into!(s_arr, sx)
    end

    # Last the last vertex with the point in question
    idxs = ((dim + 1)*(dim)+1):((dim+1)^2-1)
    fill_at_inds!(s_arr, point, idxs)

    signs[end] = sign(mdet_fm_in_x_view!(temp_arr, s_arr, dim + 1, starts, stops))

    if !(signs[end-1] == signs[end])
       return false
    else
        return true
    end
end


function innerloop_superoptim!(inds::Vector{Int}, signs, s_arr, Sj, pt, dim::Int, M, i::Int, temp_arr, starts, stops)
    for j in 1:length(inds)
        if contained_optim!(signs, s_arr, Sj[j], pt, dim, temp_arr, starts, stops)
            M[inds[j], i] += 1.0
        end
    end
end
