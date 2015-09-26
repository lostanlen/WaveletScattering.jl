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
# exception handling
bank = FourierNonOriented1DBank(Morlet1DSpec())
@test_throws ErrorException FourierNonOriented1DBank{Float32}(
    Morlet1DSpec(Float64))
@test_throws ErrorException FourierNonOriented1DBank{Float64}(
    Morlet1DSpec(Float32))
# single-core mode
spec = Morlet1DSpec(nFilters_per_octave = 24, max_scale = 4410,
    nOctaves = 8, log2_size = 16)
bank = FourierNonOriented1DBank(spec)
@test length(bank.metas) == length(bank.ψs)
@test typeof(bank.ψs) ==
    Array{WaveletScattering.AbstractFourier1DFilter{Float32},1}
@test typeof(bank.ϕ) == WaveletScattering.Symmetric1DFilter{Float32}
# multi-core mode

# FourierOriented1DBank
# exception handling
bank = FourierOriented1DBank(Morlet1DSpec())
@test_throws ErrorException FourierOriented1DBank{Float32}(
    Morlet1DSpec(Float64))
@test_throws ErrorException FourierOriented1DBank{Float64}(
    Morlet1DSpec(Float32))
# single-core mode
spec = Morlet1DSpec(nFilters_per_octave = 24, max_scale = 4410,
    nOctaves = 8, log2_size = 16)
bank = FourierOriented1DBank(spec)
@test length(bank.metas) == length(bank.ψs)
@test typeof(bank.ψs) ==
    Array{WaveletScattering.AbstractFourier1DFilter{Float32},2}
@test typeof(bank.ϕ) == WaveletScattering.Symmetric1DFilter{Float32}
