using Base.Test
# weighting.jl
import WaveletScattering: EqualWeighting, LoudnessWeighting, weight_frequencies

# EqualWeighting
weighting = EqualWeighting()
@test_approx_eq weight_frequencies(weighting, [0.4, 0.2, 0.1]) [1.0, 1.0, 1.0]

# LoudnessWeighting
weighting = LoudnessWeighting(44100)
@test_approx_eq weight_frequencies(weighting, 0.0) 0.0
@test_approx_eq weight_frequencies(LoudnessWeighting(12000.0), [1e48]) 0.0
