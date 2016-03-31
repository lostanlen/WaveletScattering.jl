using Base.Test
# weighting.jl
import WaveletScattering: EqualWeighting, weight_frequencies

# EqualWeighting
weighting = EqualWeighting()
@test_approx_eq weight_frequencies(weighting, [0.4, 0.2, 0.1]) [1.0, 1.0, 1.0]
