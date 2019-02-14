using RecipesBase
using LaTeXStrings


@recipe function f(to::TransferOperatorTriangulationExact)
    seriestype := :heatmap

    xlabel --> "i"
    ylabel --> "j"
    colorbar_title --> L"p(S_i \to im(S_j))"
    xticks --> false
    yticks --> false
    size --> (600, 500)
    to.transfermatrix
end

@recipe function f(to::TransferOperatorTriangulationApprox)
    seriestype := :heatmap


    xlabel --> "i"
    ylabel --> "j"
    colorbar_title --> L"p(S_i \to im(S_j))"
    xticks --> false
    yticks --> false
    size --> (600, 500)

    to.transfermatrix
end


@recipe function f(to::TransferOperatorRectangularBinning)
    seriestype := :heatmap


    xlabel --> "i"
    ylabel --> "j"
    colorbar_title --> L"p(b_i \to im(b_j))"
    xticks --> false
    yticks --> false
    size --> (600, 500)

    to.transfermatrix
end