using Base.Test
# bank1d.jl
import WaveletScattering: Bank1D
# spec1d.jl
import WaveletScattering: Spec1D

W = Bank1D(Spec1D())
@test ndims(W) == 1

x = map(Float32, randn(1 << W.spec.log2_size[1]))
@test W(x)
