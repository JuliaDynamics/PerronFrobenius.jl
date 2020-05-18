using CausalityToolsBase
using PerronFrobenius

@test Grid(RectangularBinning(5)) isa Grid
@test Grid(RectangularBinning(5), "circular") isa Grid