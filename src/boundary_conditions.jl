
# TODO: add projected boundary condition for grid estimator?
function isboundarycondition(bc, method::String)
    if method == "grid"
        bc ∈ ["circular", "random"]
    elseif method ∈ ["triangulation"]
        bc ∈ ["circular", "random"]
    else
        error("method $method not defined")
    end
end