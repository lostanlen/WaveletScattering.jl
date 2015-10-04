immutable Morlet1DSpec{T<:Number} <: Abstract1DSpec{T}
    ɛ::Float64
    log2_size::Tuple{Int}
    max_qualityfactor::Float64
    max_scale::Float64
    motherfrequency::Float64
    nFilters_per_octave::Int
    nOctaves::Int
    signaltype::Type{T}
    function call{T<:Number}(::Type{Morlet1DSpec{T}}, ::Type{T};
                             ɛ = default_ɛ(T), log2_size = 15,
                             max_qualityfactor = nothing, max_scale = Inf,
                             nFilters_per_octave = nothing, nOctaves = nothing,
                             tuningfrequency = nothing)
        "Integer log2_size is automatically converted to one-element tuple"
        isa(log2_size, Int) && (log2_size = tuple(log2_size))
        max_qualityfactor, nFilters_per_octave =
             default_max_qualityfactor(max_qualityfactor, nFilters_per_octave),
             default_nFilters_per_octave(nFilters_per_octave, max_qualityfactor)
        motherfrequency = tune_motherfrequency(tuningfrequency, Morlet1DSpec,
                                               nFilters_per_octave)
        nOctaves = default_nOctaves(nOctaves, Morlet1DSpec, log2_size,
                                    max_qualityfactor, max_scale,
                                    motherfrequency, nFilters_per_octave)
        spec = new{T}(ɛ, log2_size, max_qualityfactor, max_scale,
                      motherfrequency, nFilters_per_octave, nOctaves, T)
        checkspec(spec) && return spec
    end
end

"By default, `Morlet1DSpec operates on single-precision real input (Float32)."
Morlet1DSpec(T=Float32; args...) = Morlet1DSpec{T}(T; args...)

"""In the special case `nFilters_per_octave=1`, we manually set `ξ=0.39`, which
is more accurate with the Littlewood-Paley energy conservation criterion than
the generic fallback `ξ=0.4`, which is only valid when the wavelet has a
symmetric profile in the Fourier domain. This is no longer the case for
nFilters_per_octave==max_qualityfactor==1 as the Morlet low-frequency corrective
term is no longer negligible."""
default_motherfrequency(::Type{Morlet1DSpec}, nFilters_per_octave) =
    nFilters_per_octave==1 ? 0.39 : inv(3.0 - exp2(-1.0/nFilters_per_octave))

"""Computes gauss{T<:Real}(ω, den::T) = exp(- ω*ω/den).
The Gaussian bell curve is defined as gauss(ω) = exp(- ω² / 2σ²).
For performance reasons, we memoize the denominator 2σ², which is computed only
once in the caller morlet1d.
Also note that the exponentiation ω^2 is replaced by the explicit product ω*ω.
"""
gauss{T<:Real}(ω, den::T) = @fastmath convert(T, exp(- ω*ω/den))::T

"""Computes a one-dimensional Morlet wavelet in the Fourier domain.
A Morlet wavelet of center frequency ξ and of variance σ looks almost like
a Gaussian bell curve. To ensure that the wavelet has a vanishing moment, we
substract a corrective term around the zeroth frequency. Since we operate over
signals of finite length N, the corrective term must also be applied around the
frequencies -N, +N, and +2N."""
function fourierwavelet{T<:Real}(meta::AbstractMeta, spec::Morlet1DSpec{T})
    """1. **Gaussian denominator `den = 2σ²`**
    The FWHM (full width at half maximum) bw of a Gaussian bell curve of
    variance `σ` is defined by the equation
        `g(±bw/2) = exp(- (bw/2)²/(2σ²)) = 1/sqrt(2)`
    which leads to
        `bw² = 2 log(2) * 2σ²`.
    The denominator `den = 2σ²` of the Gaussian is thus equal to
        `den = 2σ² = bw² / (2 log(2))."""
    @inbounds log2_length = spec.log2_size[1]
    halfN = 1 << (log2_length - 1)
    N = halfN << 1
    center = N * T(meta.centerfrequency)
    bw = N * T(meta.bandwidth)
    den = @fastmath bw * bw / T(2.0 * log(2.0))

    """5. **Periodization**"""
    morlet = reshape(morlet, (div(length(morlet), nPeriods), nPeriods))
    morlet = squeeze(sum(morlet, 2), (2,))
    """6. **Trimming to true support boundaries**"""
    return AbstractFourier1DFilter(morlet, spec)
end

function morlet{T<:Number}(center::T, den::T, N::Int, nPeriods::Int)
    halfN = N >> 1
    pstart = - (nPeriods-1)>>1
    pstop = (nPeriods-1)>>1 + iseven(nPeriods)
    ωstart = - halfN + pstart * N
    ωstop = halfN + pstop * N - 1
    @inbounds begin
        gauss_center = T[ gauss(ω-center, den) for ω in ωstart:ωstop]
        gauss_0 = T[ gauss(ω, den)
            for ω in (ωstart-(pstart*N)):(ωstop+(pstop*N)) ]
        corrective_gaussians = [ gauss_0[1 + ω + p*N]
            for ω in 0:(N*nPeriods-1), p in 1:nPeriods ]
    end
    b = [ gauss(p*N - center) for p in pstart:pstop ]
    A = [ gauss((q-p)*N - center) for p in 0:(nPeriods-1), q in 0:(nPeriods-1) ]
    corrective_factors = A \ b
    return gauss_center - sum(corrective_factors .* corrective_gaussians, 2)
end

function scalingfunction{T<:Number}(spec::Morlet1DSpec{T})
    halfN = 1 << (spec.log2_size[1] - 1)
    bw = T( (1 << (spec.log2_size[1] - spec.nOctaves)) *
         uncertainty(spec) * spec.motherfrequency )
    den = @fastmath bw * bw / T(2.0 * log(2.0))
    lastω = round(Int, min(bw * sqrt(2.0 / spec.ɛ), halfN - 1))
    leg = T[ gauss(ω, den) for ω in 1:lastω ]
    leg = leg[1:findlast(leg .> spec.ɛ)]
    return Symmetric1DFilter(leg, one(T))
end

"""
By neglecting the low-frequency corective term, we write the Morlet wavelet as a
Gaussian of variance σ in the Fourier domain. Its 3 dB bandwidth, defined as the
full width at half maximum (FWHM) of the squared magnitude in the Fourier
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
uncertainty(::Type{Morlet1DSpec}) = sqrt(log(10.0) / log(2.0))
uncertainty{T<:Number}(::Type{Morlet1DSpec{T}}) = uncertainty(Morlet1DSpec)
