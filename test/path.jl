import WaveletScattering: Literal, PathKey

# Literal
@test isimmutable(Literal(:time))
@test Literal(:time).depth == 1
@test isimmutable(Literal((:γ, 2)))
@test Literal((:γ, 2)).depth == 2
