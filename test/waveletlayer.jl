using Base.Test
# waveletlayer.jl
import WaveletScattering: WaveletLayer
# bank1d.jl
import WaveletScattering: Bank1D
# spec1d.jl
import WaveletScattering: Spec1D

bank = Bank1D(Spec1D())
bottoms = [:fourier]
name = "wavelets"
tops = [:wavelets]
waveletlayer = WaveletLayer(bank, bottoms, name, tops)

@test Mocha.can_do_bp(waveletlayer)
@test !Mocha.has_neuron(waveletlayer)
@test !Mocha.has_param(waveletlayer)
