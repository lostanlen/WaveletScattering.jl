import WaveletScattering: issteerable, MexicanHat, Morlet

@test !issteerable(MexicanHat())
@test issteerable(Morlet())
