using Base.Test
# spec2d.jl
import WaveletScattering: Spec2D

spec = Spec2D()
@test isa(spec, Spec2D)
