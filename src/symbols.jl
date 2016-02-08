function appendsymbols(symbols::Vector{Symbol}, N::Int)
    return vcat(symbols, [ symbol(:var, n) for n in 1:(N-length(symbols)) ])
end
