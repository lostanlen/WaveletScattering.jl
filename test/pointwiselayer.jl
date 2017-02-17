using Base.Test
import WaveletScattering: PointwiseLayer

modulus = PointwiseLayer(bottoms = [:bottom], tops = [:top])

@test !Mocha.can_do_bp(modulus)
@test !Mocha.has_neuron(modulus)
@test !Mocha.has_param(modulus)
