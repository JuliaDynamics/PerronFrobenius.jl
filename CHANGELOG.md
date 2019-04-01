## Release 0.6.0

### Breaking changes
- Rectangular binning schemes must now be explicitly specified (not inferred from numeric type). 

### Improvements
- Documentation fixes.
- Added more examples to doc strings.

## Release 0.5.2

- Require newest StateSpaceReconstruction.jl version, so that invariantizing works properly after PyCall syntax update.

## Release 0.5.1

### Breaking changes
- `rectangularinvariantemeasure` now always needs an instance of `RectangularBinning` to 
    specify the binning scheme.

## Release 0.4.0

### New functionality

- **Invariant measures** may be computed on triangulations with the syntax `invariantmeasure(pts, ϵ::TriangulationBinning, 
    simplex_intersection_type::ApproxIntersection)` for approximate simplex intersections, and `invariantmeasure(pts, ϵ::TriangulationBinning, simplex_intersection_type::ExactIntersection)` for example intersections. If for example `pts = [rand(3) for i = 1:50]`, then you can do `invariantmeasure(pts, TriangulationBinning(), ApproxIntersection())`. This will return `TransferOperatorTriangulationExact` or `TransferOperatorTriangulationApprox` instances, which holds all the information used to compute the measure.
- **Transfer operators** may be computed without storing any other information by using the same syntax, e.g. `invariantmeasure(pts, TriangulationBinning(), ApproxIntersection())`.

### Bug fixes
- Fixed redundant method definitions.
