using Base.Test
# import modulus.jl
import WaveletScattering: Modulus, SquaredModulus

@test_approx_eq Modulus()([-2.0, 3.0 + 4.0im, -0.0]) [2.0, 5.0, 0.0]

@test_approx_eq SquaredModulus()([-2.0, 3.0 + 4.0im, -0.0]) [4.0, 25.0, 0.0]
