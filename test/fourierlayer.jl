using Base.Test
# fourierlayer.jl
import WaveletScattering: FourierLayer

fourier = FourierLayer(
        bottoms = [:signal],
        tops = [:fourier],
        pathkeys = [PathKey(:time)])

@test Mocha.can_do_bp(fourier)
@test !Mocha.has_neuron(fourier)
@test !Mocha.has_param(fourier)
