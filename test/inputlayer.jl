using Base.Test
# fourierlayer.jl
import WaveletScattering: InputLayer

inputlayer = InputLayer(
        tops = [:signal],
        symbols = [:time, :chunk],
        data = rand(Float32, 256, 2))

@test !Mocha.can_do_bp(inputlayer)
@test !Mocha.has_neuron(inputlayer)
@test !Mocha.has_param(inputlayer)
@test Mocha.is_source(inputlayer)
