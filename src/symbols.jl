function appendsymbols(symbols::Vector{Symbol}, N::Int)
    nSymbols = length(symbols)
    @assert nSymbols <= N
    @assert nSymbols > 0
    return vcat(symbols, [ symbol(:var, n) for n in 1:(N-nSymbols)) ])
end
