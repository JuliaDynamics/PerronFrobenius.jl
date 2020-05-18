export Grid
using CausalityToolsBase
using SparseArrays

include("v1/binvisits.jl")

# TODO: add projected boundary condition?
function isboundarycondition(bc)
    bc ∈ ["circular", "random"]
end

"""
    Grid(b::BinningScheme, bc::String = "none", f::Real = 1.0) <: TransferOperatorEstimator

A simple grid estimator for approximating the transfer operator. 
`bc` is the boundary condition, and `f` is the allocation factor.
"""
struct Grid{B <: BinningScheme} <: TransferOperator
    b::B
    bc::String
    f::Real
     
    function Grid(b::B, bc::String = "circular", f::Real = 1.0) where B
        isboundarycondition(bc)  || error("Boundary condition '$bc' not valid.")
        new{B}(b, bc, f)
    end
end
Base.show(io::IO, g::Grid) = print(io, "Grid{$(g.b)}")


function transferoperatorgenerator(pts, method::Grid)
    b = method.b
    
    # Calculate minima and edgelengths given the chosen grid
    mini, el = get_minima_and_edgelengths(pts, b)
    
    # Tie points to grid cells using integer tuples
    encoded_pts = encode(pts, mini, el)
    
    # Coordinates of the visited bins relative to grid.
    vb = [el .* pt .+ mini for pt in encoded_pts]
    
    # TODO: this is the expensive operation. Speed ups possible by 
    # making versions of the groupslice functions that operate on 
    # integers tuples instead of this matrix of integers?
    bv = get_binvisits(hcat(encoded_pts...,))
    
    init = (vb = vb, bv = bv)
    TransferOperatorGenerator(method, pts, init)
end


Base.show(io::IO, g::TransferOperatorGenerator) = print(io, "TransferOperatorGenerator{$(g.method)}")

# Stricly speaking, for the single-grid estimator, making a generator is not necessary.
# However, follow the conventions for the rest of the estimators and use it.
function (tog::TransferOperatorGenerator{<:Grid})(bc)
    boundary_condition = tog.method.bc

    alloc_frac = tog.method.f
    if !(0 < alloc_frac <= 1)
        error("allocation fraction needs to be on the interval (0, 1]")
    end

    bv, vb = getfield.(Ref(tog.init), (:bv, :vb))
    first_visited_by = bv.first_visited_by
    visitors = bv.visitors
    visits_whichbin = bv.visits_whichbin
    n_visited = length(first_visited_by)
    valid_boundary_conditions = [:none, :exclude, :circular, :invariantize]


    # Initialise transfer (Perron-Frobenius) operator as a sparse matrix
    # (keep track of indices and values in separate columns for now)
    n_possible_nonzero_entries = n_visited^2
    N = ceil(Int, n_possible_nonzero_entries / alloc_frac)

    I = zeros(Int, N)
    J = zeros(Int, N)
    P = zeros(N)

    # Preallocate target index for the case where there is only
    # one point of the orbit visiting a bin.
    target_bin_j::Int = 0
    n_visitsᵢ::Int = 0

    # Keep track of how many transitions we've considered.
    k = 0

    if boundary_condition == "circular"
        append!(visits_whichbin, [1])
    elseif boundary_condition == "random"
        append!(visits_whichbin, [rand(1:length(visits_whichbin))])
    end

    # Loop over the visited bins bᵢ
    for i in 1:n_visited
        # How many times is this bin visited?
        n_visitsᵢ = length(visitors[i])

        # If both conditions below are true, then there is just one
        # point visiting the i-th bin. If there is only one visiting point and
        # it happens to be the last, we skip it, because we don't know its
        # image.
        if n_visitsᵢ == 1 && !(i == visits_whichbin[end])
            k += 1
            # To which bin does the single point visiting bᵢ jump if we shift it one time step ahead along its orbit?
            target_bin_j = visits_whichbin[visitors[i][1] + 1][1]

            # We now know that exactly one point (the i-th) does the
            # transition from i to the target j.
            I[k] = i
            J[k] = target_bin_j
            P[k] = 1
        end
        # If more than one point of the orbit visits the i-th bin, we
        # identify the visiting points and track which bins bⱼ they end up
        # in after the forward linear map of the points.
        if n_visitsᵢ > 1
            timeindices_visiting_pts = visitors[i]

            # TODO: Introduce circular boundary condition. Simply excluding
            # might lead to a cascade of loosing points.

            # If bᵢ is the bin visited by the last point in the orbit, then
            # the last entry of `visiting_pts` will be the time index of the
            # last point of the orbit. In the next time step, that point will
            # have nowhere to go along its orbit (precisely because it is the
            # last data point). Thus, we exclude it.
            if i == visits_whichbin[end]
                #warn("Removing last point")
                n_visitsᵢ = length(timeindices_visiting_pts) - 1
                timeindices_visiting_pts = timeindices_visiting_pts[1:(end - 1)]
            end

            # To which boxes do each of the visitors to bᵢ jump in the next
            # time step?
            target_bins = visits_whichbin[timeindices_visiting_pts .+ 1]
            unique_target_bins = unique(target_bins)

            # Count how many points jump from the i-th bin to each of
            # the unique target bins, and use that to calculate the transition
            # probability from bᵢ to bⱼ.
            for j in 1:length(unique_target_bins)
                n_transitions_i_to_j = sum(target_bins .== unique_target_bins[j])

                k += 1
                I[k] = i
                J[k] = unique_target_bins[j]
                P[k] = n_transitions_i_to_j / n_visitsᵢ
            end
        end
    end

    # Combine indices and values into a sparse representation of the Perron-
    # Frobenius operator (transfer operator). Filter out the nonzero entries.
    # The number of rows in the matrix is given by the number of unique points
    # we have in the embedding.
    TO = Array(sparse(I[1:k], J[1:k], P[1:k], n_visited, n_visited))
    
    # There may be boxes which are visited by points of the orbit, but not by
    # any image points.
    # zc = zerocols(TO)
    # l = length(zc)
    # if l > 0
    #     col = zc[1]
    #     warn("There were $l all-zero columns. Column $col is the culprit -> normalizing ...")
    #     TO = TO ./ sum(TO, 2)
    # end
    # If the first bin is not visited, make every point starting in that bin
    # jump to the same bin with probability one.
    # zc = zerocols(TO)
    # if length(zc) > 0
    #     warn("There were $l all-zero columns: column $col")
    #     @show zc
    #     i = zc[1]
    #     TO[i, 1] = 1.0
    # end
end

"""
    transferoperator(pts, method::Grid, bc::String)

Compute the transfer operator from the phase/state space `points`,
using boundary condition `bc` ("circular" or "random"), over a
rectangular partition of the state space.
"""
function transferoperator(pts, method::Grid, bc::String = "circular")
    tog = transferoperatorgenerator(pts, method)
    tog(bc)
end
