abstract RedundantWaveletClass

immutable Gammatone <: RedundantWaveletClass end

immutable MexicanHat <: RedundantWaveletClass end
issteerable(::MexicanHat) = false

immutable Morlet <: RedundantWaveletClass end
issteerable(::Morlet) = true
