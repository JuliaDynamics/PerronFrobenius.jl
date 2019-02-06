"""

##
```julia
# Some random points
pts = [rand(3) for i = 1:32]

# Invariantize the points
pts = invariantize(pts)

# Create a forwardmap
ForwardmapTriangulation(pts, triangulate_f(pts[1:end-1]), triangulate_f(pts[2:end]))
```
"""
struct ForwardmapTriangulation{DT, TT} <: AbstractForwardMapTriangulation{DT, TT}
    data::DT
    triang::TT
    forwardmap_triang::TT
end

function ForwardMapTriangulation(pts::DT)
    ForwardmapTriangulation(pts, triangulate_f(pts[1:end-1]), triangulate_f(pts[2:end]))
end

####################
# Pretty printing
####################
function summarise(fmt::ForwardmapTriangulation)
    _type = typeof(fmt)
    n_simplices_triang = nsimplices(fmt.triang)
    n_simplices_forwardmap_triang = nsimplices(fmt.forwardmap_triang)
    n_pts = length(fmt.data)
    summary = "$_type"
end

Base.show(io::IO, fmt::ForwardmapTriangulation) = println(io, summarise(fmt))


export ForwardMapTriangulation
