using Base.Test
# pointwise.jl
import WaveletScattering: Identity, Log1P

ρ = Identity()
@test_approx_eq ρ([1.0, 2.0, 3.0]) [1.0, 2.0, 3.0]

threshold = 0.1
ρ = Log1P(threshold)
@test_approx_eq ρ([0.0]) [0.0]
@test_approx_eq expm1(ρ([1.0, 2.0, 3.0])) [1.0, 2.0, 3.0]*threshold
@test_throws DomainError ρ([- 1.001 / threshold])
