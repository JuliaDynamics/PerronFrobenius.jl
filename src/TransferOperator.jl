abstract type TransferOperator <: AbstractArray{Float64, 2} end

type ExactSimplexTransferOperator <: TransferOperator
    TO::AbstractArray{Float64, 2}
end

type ApproxSimplexTransferOperator <: TransferOperator
    TO::AbstractArray{Float64, 2}
end

type EquidistantBinningTransferOperator  <: TransferOperator
    TO::AbstractArray{Float64, 2}
end

is_markov(TO::TransferOperator) = all(sum(TO.TO, 2) .≈ 1)
is_almostmarkov(TO::TransferOperator; tol = 0.01) = all(sum(TO.TO, 2) .> (1 - tol))

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


"""
    max_discrep(to1::ExactSimplexTransferOperator,
                to2::ApproxSimplexTransferOperator)

Measures the maximum discrepancy between an exact and an approximate estimate
of the transfer operator.

For this function to work, the transfer operators must have been generated from
the same triangulation (otherwise, the number of simplices don't match up).
"""
function max_discrep(exact::ExactSimplexTransferOperator,
                     approx::ApproxSimplexTransferOperator)

end
