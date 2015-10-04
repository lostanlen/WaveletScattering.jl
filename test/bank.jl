using Base.Test
# bank.jl
import WaveletScattering: Behavior, FourierNonOriented1DBank,
    FourierOriented1DBank
# fourierfilter.jl
import WaveletScattering: AbstractFourier1DFilter, Symmetric1DFilter
# morlet1d.jl
import WaveletScattering: Morlet1DSpec

# exception handling
# FourierNonOriented1DBank
@test_throws ErrorException FourierNonOriented1DBank{Float32}(
    Morlet1DSpec(Float64))
@test_throws ErrorException FourierNonOriented1DBank{Float64}(
    Morlet1DSpec(Float32))
# FourierOriented1DBank
@test_throws ErrorException FourierOriented1DBank{Float32}(
    Morlet1DSpec(Float64))
@test_throws ErrorException FourierOriented1DBank{Float64}(
    Morlet1DSpec(Float32))

# single-core mode
# FourierNonOriented1DBank
spec = Morlet1DSpec(nFilters_per_octave = 24, max_scale = 4410,
    nOctaves = 8, log2_size = 16)
bank = FourierNonOriented1DBank(spec)
@test length(bank.metas) == length(bank.ψs)
@test typeof(bank.ψs) ==
    Array{AbstractFourier1DFilter{Float32},1}
@test typeof(bank.ϕ) == Symmetric1DFilter{Float32}
# FourierOriented1DBank
spec = Morlet1DSpec(nFilters_per_octave = 24, max_scale = 4410,
    nOctaves = 8, log2_size = 16)
bank = FourierOriented1DBank(spec)
@test size(bank.metas) == size(bank.ψs)
@test typeof(bank.ψs) ==
    Array{AbstractFourier1DFilter{Float32},2}
@test typeof(bank.ϕ) == Symmetric1DFilter{Float32}
