using Base.Test
# bank1d.jl
import WaveletScattering: Bank1D
# spec1d.jl
import WaveletScattering: Spec1D

bank = Bank1D(Spec1D())
@test ndims(bank) == 1
