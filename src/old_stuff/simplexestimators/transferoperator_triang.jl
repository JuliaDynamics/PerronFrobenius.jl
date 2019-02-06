
"""
    transferoperator(triang::AbstractTriangulation;
                exact = false, parallel = true,
                n_pts = 200, sample_randomly = false)

Estimate the transfer operator from a triangulation.
"""
function transferoperator_triang(triang::AbstractTriangulation;
                            exact = false, parallel = true,
                            n_pts = 200, sample_randomly = false)
    if exact
        if parallel
            transferoperator_exact_p(triang)
        else
            transferoperator_exact(triang)
        end
    else
        transferoperator_approx(triang, n_pts = n_pts,
                                sample_randomly = sample_randomly)
    end
end
