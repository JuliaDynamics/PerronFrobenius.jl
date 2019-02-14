
""" Re-zeros the array `a`. """
function rezero!(a)
    @inbounds for i in eachindex(a)
        a[i] = 0.0
    end
end

""" Fill the elements of vector `from` into vector `into`. """
function fill_into!(into, from)
    @inbounds for i in eachindex(into)
        into[i] = from[i]
    end
end

""" Fill vector `into` at indices `inds` with the elements at `from[inds]`. """
function fill_at_inds!(into, from, inds)
    @inbounds for i in 1:length(inds)
        into[inds[i]] = from[i]
    end
end

""" Check whether a point is contained within the simplex currently being checked """
function contained!(signs, s_arr, sx, point, dim)
    # Redefine the temporary simplex. This is in-place, so we don't allocate
    # memory. We could also have re-initialised `signs`, but since we're never
    # comparing more than two consecutive signs, this is not necessary.
    rezero!(s_arr)
    rezero!(signs)
    fill_into!(s_arr, sx)

    #Replace first vertex with the point
    fill_at_inds!(s_arr, point, 1:dim)

    # Signed volume
    signs[1] = sign(det(reshape(s_arr, dim + 1, dim + 1)))

    rezero!(s_arr)
    fill_into!(s_arr, sx) #reset

    for κ = 2:dim # Check remaining signs and stop if sign changes
        # Replace the ith vertex with the point we're cheking (leaving the
        # 1 appended to Vi intact.)
        idxs = ((dim + 1)*(κ - 1)+1):((dim + 1)*(κ - 1)+ 1 + dim - 1)
        fill_at_inds!(s_arr, point, idxs) # ith change

        signs[κ] = sign(det(reshape(s_arr, dim + 1, dim + 1)))

        if !(signs[κ-1] == signs[κ])
            return false
        end

        rezero!(s_arr)
        fill_into!(s_arr, sx)
    end

    # Last the last vertex with the point in question
    idxs = ((dim + 1)*(dim)+1):((dim+1)^2-1)
    fill_at_inds!(s_arr, point, idxs)

    signs[end] = sign(det(reshape(s_arr, dim + 1, dim + 1)))

    if !(signs[end-1] == signs[end])
       return false
    else
        return true
    end
end

""" The inner loop in the computation of the transfer matrix. """
function innerloop!(inds::Vector{Int}, signs, s_arr, Sj, pt, dim::Int, M, i::Int)
    for j in 1:length(inds)
        if contained!(signs, s_arr, Sj[j], pt, dim)
            M[inds[j], i] += 1.0
        end
    end
end

""" Get the simplices at the given indices. """
function get_simplices_at_inds!(simps, inds::Vector{Int}, simplices)
    for i in 1:length(inds)
        simps[i] = simplices[:, inds[i]]
    end
end
