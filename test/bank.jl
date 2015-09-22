using Base.Test
# bank.jl
import WaveletScattering: Behavior, FourierNonOriented1DBank,
    FourierOriented1DBank
# morlet1d.jl
import WaveletScattering: Morlet1DSpec

# Behavior
js = Int8[0,0,1,1,2,2,3,3]
b = Behavior(js)
@test b.γ_range == 0:7
@test b.log2_oversampling == 0
@test b.min_log2_resolution == -2

# FourierNonOriented1DBank
spec = Morlet1DSpec()
bank = FourierNonOriented1DBank(spec)
