using Base.Test
import WaveletScattering: PointwiseLayer, Modulus

modulus = PointwiseLayer(
        bottoms = [:fourier],
        tops = [:modulus],
        ρ = Modulus())


@test Mocha.can_do_bp(modulus)
@test !Mocha.has_neuron(modulus)
@test !Mocha.has_param(modulus)
