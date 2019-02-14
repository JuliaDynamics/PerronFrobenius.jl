## Release 0.4.0

### New functionality

- **Invariant measures** may be computed on triangulations with the syntax `invariantmeasure(pts, ϵ::TriangulationBinning, 
    simplex_intersection_type::ApproxIntersection)` for approximate simplex intersections, and `invariantmeasure(pts, ϵ::TriangulationBinning, simplex_intersection_type::ExactIntersection)` for example intersections. If for example `pts = [rand(3) for i = 1:50]`, then you can do `invariantmeasure(pts, TriangulationBinning(), ApproxIntersection())`. This will return `TransferOperatorTriangulationExact` or `TransferOperatorTriangulationApprox` instances, which holds all the information used to compute the measure.
- **Transfer operators** may be computed without storing any other information by using the same syntax, e.g. `invariantmeasure(pts, TriangulationBinning(), ApproxIntersection())`.

### Bug fixes
- Fixed redundant method definitions.