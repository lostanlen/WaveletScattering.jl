import WaveletScattering: isdyadic, issteerable, Gammatone, MexicanHat, Morlet

@test Gammatone() == Gammatone(4)
@test !isdyadic(Gammatone())
@test_throws MethodError issteerable(Gammatone())

@test isdyadic(MexicanHat())
@test !issteerable(MexicanHat())

@test !isdyadic(Morlet())
@test issteerable(Morlet())
