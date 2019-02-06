"""
    ForwardMapTriangulation{DT, TT}

A triangulation containing a triangulation of the

## Fields
- **`data`**. The data
- **`triang`**. A triangulation of points ``p_1, p_2, \\ldots, p_{N-1}``.
- **`forwardmap_triang`**. A triangulation of points ``p_2, p_2, \\ldots, p_{N}``.

"""
abstract type AbstractForwardMapTriangulation{DT, TT} end
