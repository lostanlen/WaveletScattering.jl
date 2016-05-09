using Base.Test
# spec2d.jl
import WaveletScattering: Spec2D, default_max_aspectratio, default_nOrientations
# waveletclass.jl
import WaveletScattering: MexicanHat, Morlet

# Spec2D constructor
spec = Spec2D()
@test isa(spec, Spec2D)

# default_max_aspectratio
@test_approx_eq default_max_aspectratio(Morlet(), nothing) 2.0
@test_approx_eq default_max_aspectratio(MexicanHat(), nothing) 1.0
@test_approx_eq default_max_aspectratio(Morlet(), 3.0) 3.0

# default_nOrientations
@test default_nOrientations(Morlet(), nothing) == 4
@test default_nOrientations(MexicanHat(), nothing) == 1
@test default_nOrientations(Morlet(), 3) == 3
