abstract RedundantWaveletClass

immutable Gammatone <: RedundantWaveletClass
    order::Int
end
Gammatone() = Gammatone(4)
isdyadic(::Gammatone) = false

immutable MexicanHat <: RedundantWaveletClass end
isdyadic(::MexicanHat) = true
issteerable(::MexicanHat) = false

immutable Morlet <: RedundantWaveletClass end
isdyadic(::Morlet) = false
issteerable(::Morlet) = true

"""
By neglecting the low-frequency corrective term, we write the Morlet wavelet as
a Gaussian of variance σ in the Fourier domain. Its 3 dB bandwidth, defined as
the full width at half maximum (FWHM) of the squared magnitude in the Fourier
domain, is then equal to b = 2σ*sqrt(log(2)).

Therefore, for a given center frequency ω and a quality factor Q, the variance
σ of the Gaussian is equal to σ = b / (2*sqrt(log(2))). In the spatial domain,
this amounts to a Gabor wavelet (a Gaussian modulated by a sine wave, without
any low-frequency corrective term) of variance 1/σ. Its spatial scale (FWTM)
is equal to s = 2*sqrt(log(10))/σ. We conclude that the uncertainty
constant of the Morlet wavelet, defined as the (scale*bandwidth) constant,
is equal to
    h = b*s = sqrt(log(10)/log(2)) = 1.8226...
"""
uncertainty(class::Morlet) = sqrt(log(10.0) / log(2.0))
