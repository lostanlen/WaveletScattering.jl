using Base.Test
# filter.jl
import WaveletScattering: AbstractFilter, spin!
# domain.jl
import WaveletScattering: FourierDomain
# fourier1dfilter.jl
import WaveletScattering: FullResolution1DFilter, spin

ψ = FullResolution1DFilter(Float32[0.01, 0.1, 0.2, 0.3])
ψ2 = Float32(2.0) * ψ
@test isa(ψ2, FullResolution1DFilter{Float32})
@test_approx_eq ψ2.coeff Float32[0.02, 0.2, 0.4, 0.6]


ψ = FullResolution1DFilter(Float32[0.01, 0.1, 0.2, 0.3])
ψ2 = ψ .* Float32(2.0)
@test isa(ψ2, FullResolution1DFilter{Float32})
@test_approx_eq ψ2.coeff Float32[0.02, 0.2, 0.4, 0.6]

ψs = Array{AbstractFilter{Float32,FourierDomain{1}}}(2, 1, 1)
ψs[1, 1, 1] = ψ .* 2.0
spin!(ψs)
@test ψs[1, 1, 1].coeff == (ψ .* 2.0).coeff
@test ψs[2, 1, 1].coeff == spin(ψ .* 2.0).coeff
