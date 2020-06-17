import DelayEmbeddings: AbstractDataset 

function invariantize(pts::AbstractDataset, bc = "circular")
    invariant_pts = copy(pts.data)
    if bc == "circular"
        return push!(invariant_pts, pts.data[1])
    elseif bc == "random"
        return push!(invariant_pts, pts.data[rand(1:length(pts))])
    end
    return invariant_pts
end

