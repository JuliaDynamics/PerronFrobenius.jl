import StateSpaceReconstruction:
    DelaunayTriangulation,
    Simplex,
    centroid,
    radius

"""
idxs_potentially_intersecting_simplices(all_pts::Vector{Vector{T}},
        DT::DelaunayTriangulation, idx::Int) where T

Given a delaunay triangulation `DT` of `all_pts[1:(end - 1)]`, find the indices of
the simplices that potentially intersect with the image simplex with index `idx`.
"""
function idxs_potentially_intersecting_simplices(all_pts,
	        DT::DelaunayTriangulation, idx::Int) where T
	# Vector that will hold the indices of the simplices potentially
	# intersecting with the image of simplex #idx
	inds_potential_simplices = Int[]

	n_simplices = length(DT)

	original_pts = all_pts[1:(end - 1)]
	forward_pts = all_pts[2:end]

	simplices = [ Simplex(original_pts[DT[i]]) for i = 1:n_simplices]
	image_simplices = [Simplex(forward_pts[DT[i]]) for i = 1:n_simplices]

	cs = [centroid(simplices[i]) for i = 1:n_simplices]
	cs_im = [centroid(image_simplices[i]) for i = 1:n_simplices]

	rs = [radius(simplices[i]) for i = 1:n_simplices]
	rs_im = [radius(image_simplices[i]) for i = 1:n_simplices]

	@inbounds for i = 1:n_simplices
	    δ = transpose(cs_im[idx] .- cs[i]) * (cs_im[idx] .- cs[i]) .- ((rs_im[idx] + rs[i])^2)
	    if δ[1] < 0
	        push!(inds_potential_simplices, i)
	    end
	end

	return inds_potential_simplices
end
