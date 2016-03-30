using Base.Test
# symbols.jl
import WaveletScattering: appendsymbols

@test_throws AssertionError appendsymbols(Symbol[], 1)
@test_throws AssertionError appendsymbols([:time], 0)
@test appendsymbols([:time], 1) == [:time]
@test appendsymbols([:time], 2) == [:time, :var1]
@test appendsymbols([:time], 3) == [:time, :var1, :var2]
@test appendsymbols([:time, :space], 2) == [:time, :space]
@test appendsymbols([:time, :space], 3) == [:time, :space, :var1]
