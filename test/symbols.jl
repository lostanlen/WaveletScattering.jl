using Base.Test
# symbols.jl
import WaveletScattering: appendsymbols

@test appendsymbols(Symbol[], 0) == []
@test appendsymbols(Symbol[], 1) == [:var1]
@test appendsymbols(Symbol[], 2) == [:var1, :var2]
@test appendsymbols([:time], 0) == [:time]
@test appendsymbols([:time], 1) == [:time, :var1]
@test appendsymbols([:time], 2) == [:time, :var1, :var2]
@test appendsymbols([:time, :space], 0) == [:time, :space]
@test appendsymbols([:time, :space], 1) == [:time, :space, :var1]
@test appendsymbols([:time, :space], 2) == [:time, :space, :var1, :var2]
