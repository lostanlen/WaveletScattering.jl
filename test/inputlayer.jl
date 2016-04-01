using Base.Test
# fourierlayer.jl
import WaveletScattering: InputLayer

signal = InputLayer(
        tops = [:signal],
        symbols = [:time, :chunk],
        data = rand(Float32, 256, 2))

@test !Mocha.can_do_bp(signal)
@test !Mocha.has_neuron(signal)
@test !Mocha.has_param(signal)
@test Mocha.is_source(signal)
