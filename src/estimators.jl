# Load estimator functions.
include("rectangularbinestimators/equidistant.jl")
include("simplexestimators/exact.jl")
include("simplexestimators/pointapprox.jl")

function transferoperator(t::Triangulation;
                            exact = false, parallel = true,
                            n_pts = 200, sample_randomly = false)
    if exact
        if parallel
            transferoperator_exact_p(t)
        else
            transferoperator_exact(t)
        end
    else
        transferoperator_approx(t,
                                n_pts = n_pts,
                                sample_randomly = sample_randomly)
    end
end
