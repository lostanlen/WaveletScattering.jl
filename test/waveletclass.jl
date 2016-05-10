using Base.Test
import WaveletScattering: default_motherfrequency, isdyadic, issteerable,
    Gammatone, MexicanHat, Morlet

@test_approx_eq default_motherfrequency(Morlet(), 1) 0.39

@test Gammatone() == Gammatone(4)
@test !isdyadic(Gammatone())
@test_throws MethodError issteerable(Gammatone())

@test isdyadic(MexicanHat())
@test !issteerable(MexicanHat())

@test !isdyadic(Morlet())
@test issteerable(Morlet())
