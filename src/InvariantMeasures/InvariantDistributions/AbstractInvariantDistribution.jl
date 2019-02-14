export AbstractInvariantDistribution



abstract type AbstractInvariantDistribution end

####################
# Pretty printing
####################
function summarise(iv::AbstractInvariantDistribution)
    invdist_type = typeof(iv)
    n_bins = length(iv.dist)
    n_nonzero = length(iv.nonzero_inds)
    invdist_str = "$invdist_type where $n_nonzero out of $n_bins bins have positive measure\n"
    return invdist_str
end


Base.show(io::IO, iv::AbstractInvariantDistribution) = println(io, summarise(iv))
