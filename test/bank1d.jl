using Base.Test
# bank1d.jl
import WaveletScattering: Bank1D
# spec1d.jl
import WaveletScattering: Spec1D

W = Bank1D(Spec1D(n_filters_per_octave=4, n_octaves=8))
@test ndims(W) == 1
