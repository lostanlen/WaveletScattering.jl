function appendsymbols!(symbols::Vector{Symbol}, N::Int)
    symbols = vcat(symbols, [ symbol(:var, n) for n in 1:(length(symbols)-N) ])
end
