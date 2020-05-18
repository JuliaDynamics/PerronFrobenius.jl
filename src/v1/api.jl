export TransferOperator, transferoperator, transferoperatorgenerator

import CausalityToolsBase.BinningScheme

""" Supertype of all transfer operator estimators. """
abstract type TransferOperator end 

struct TransferOperatorGenerator{E <: TransferOperator, X, A}
    method::E # estimator with its input parameters
    pts::X    # the phase space / reconstruted state space points
    init::A   # pre-initialized thins that speed up estimation process
end

"""
    transferoperatorgenerator(pts, method::TransferOperator) → to::TransferOperatorGenerator

Initialize a generator that creates transfer operators on demand, based on the given `method`.
This is efficient, because some things can be initialized and reused.

To generate a transfer operator, call `to` as a function with the , e.g.:
```julia
to = transferoperatorgenerator(x, method)
for i in 1:1000
    s = to()
    # do stuff with s and or x
    result[i] = stuff
end
```
"""
function transferoperatorgenerator end

"""
    transferoperator(pts, method::TransferOperator) → to

Create a transfer operator `to` from the phase/state space points `x` based on the given `method`.
If you want to generate more than one transfer operator from `x`, you should use [`transferoperatorgenerator`](@ref).
"""
function transferoperator(pts, method::TransferOperator)
    to = transferoperatorgenerator(pts, method)
    to()
end
