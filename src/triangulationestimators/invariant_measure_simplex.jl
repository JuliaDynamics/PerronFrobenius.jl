
function invariantmeasure(pts, method::TriangulationEstimators.SimplexExact; tol = 1e-8)
    tog = transopergenerator(pts, method)
    to = tog(tol = tol)
    invariantmeasure(to)
end

function invariantmeasure(pts, method::TriangulationEstimators.SimplexPoint; tol = 1e-8, randomsampling::Bool = false, n::Int = 100)
    tog = transopergenerator(pts, method)
    to = tog(tol = tol, randomsampling = randomsampling, n = n)
    invariantmeasure(to)
end