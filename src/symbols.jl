function appendsymbols!(symbols::Vector{Symbol}, N::Int)
    length(symbols) == N && return
    symbols = vcat(symbols, [ symbol(:var, n) for n in 1:(N-length(symbols)) ])
end
