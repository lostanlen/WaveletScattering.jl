using Base.Test
# waveletlayer.jl
import WaveletScattering: WaveletLayer
# bank1d.jl
import WaveletScattering: NullBank, Bank1D
# spec1d.jl
import WaveletScattering: Spec1D

bank = Bank1D(Spec1D())
bottoms = [:fourier]
name = "wavelets"
pathkeys = [PathKey(:time)]
tops = [:wavelets]
waveletlayer = WaveletLayer(bank=bank,
    bottoms=bottoms, name=name, pathkeys=pathkeys, tops=tops)
@test Mocha.can_do_bp(waveletlayer)
@test !Mocha.has_neuron(waveletlayer)
@test !Mocha.has_param(waveletlayer)

@test_throws AssertionError WaveletLayer(bank=NullBank(),
    bottoms=bottoms, name=name, pathkeys=pathkeys, tops=tops)
@test_throws AssertionError WaveletLayer(bank=bank,
    bottoms=Symbol[], name=name, pathkeys=pathkeys, tops=tops)
@test_throws AssertionError WaveletLayer(bank=bank,
    bottoms=bottoms, name=name, pathkeys=pathkeys, tops=[:top1, :top2])
@test_throws AssertionError WaveletLayer(bank=bank,
    bottoms=[:bottom1, :bottom2], name=name, pathkeys=pathkeys, tops=tops)
