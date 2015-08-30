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
                             ɛ=default_ɛ(T), log2_size=15,
                             max_qualityfactor=nothing, max_scale=Inf,
                             nFilters_per_octave=nothing, nOctaves=nothing,
                             tuningfrequency=nothing)
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

"""Computes gauss{T<:Real}(ω::T, den::T) = exp(- ω*ω/den).
The Gaussian bell curve is defined as gauss(ω) = exp(- ω² / 2σ²).
For performance reasons, we memoize the denominator 2σ², which is computed only
once in the caller morlet1d.
Also note that the exponentiation ω^2 is replaced by the explicit product ω*ω.
"""
gauss{T<:Real}(ω::T, den::T) = @fastmath exp(- ω*ω/den)

"""Computes a one-dimensional Morlet wavelet in the Fourier domain.
A Morlet wavelet of center frequency ξ and of variance σ looks almost like
a Gaussian bell curve. To ensure that the wavelet has a vanishing moment, we
substract a corrective term around the zeroth frequency. Since we operate over
signals of finite length N, the corrective term must also be applied around the
frequency N."""
function fourierwavelet{T<:Real}(meta::AbstractMeta, spec::Morlet1DSpec{T})
    """1. **Gaussian denominator `den = 2σ²`**
    The FWHM (full width at half maximum) bw of a Gaussian bell curve of
    variance `σ` is defined by the equation
        `g(±bw/2) = exp(- (bw/2)²/(2σ²)) = 1/2`
    which leads to
        `bw² = 4 * log(2) * 2σ²`.
    The denominator `den = 2σ²` of the Gaussian is thus equal to
        `den = 2σ² = bw² / (4 log(2))."""
    log2_length = spec.log2_size[1]
    half_length = 1 << (log2_length - 1)
    N = T(half_length << 1)
    center = N * T(meta.centerfrequency)
    bw = N * T(meta.bandwidth)
    den = @fastmath bw * bw / T(4.0 * log(2.0))
    """2. **Morlet low-frequency corrective terms**
    The one-dimensional Morlet wavelet of center frequency c is defined in the
    Fourier domain under the form
        `ψ(ω) = g(ω-c) - corr0 * g(ω) - corrN * g(ω-N)`
    where `corr0` and `corrN` are corrective terms to the Gaussian bell curve
    g(ω) of variance σ to ensure one vanishing moment.
    These terms satisfy the equations `ψ(0) = 0` and `ψ(N) = 0`, which leads to
    the system
        `ψ(0) = g(c)   + corr0 * g(0) - corrN * g(N) = 0`
        `ψ(N) = g(N-c) + corr0 * g(N) - corrN * g(0) = 0`
    which, recalling that `g(0) = 1` and that `g` is symmetric, is inverted as
        `corr0 = (g(N-c) - g(N)*g(c)) / (1 - g(N)²)`
        `corrN = g(c) - g(N) * (g(N-c) - g(N)*g(c)) / (1 - g(N)²)`."""
    gauss_N = gauss(N, den)
    gauss_center = gauss(center, den)
    gauss_N_minus_center = gauss(N-center, den)
    corrN = (gauss_N_minus_center - gauss_N*gauss_center)/(1 - gauss_N*gauss_N)
    corr0 = gauss_center - gauss_N * corrN
    """3. **Conservative support boundaries**
    Since the Morlet wavelet has a fast (Gaussian-like) decay in the frequency
    domain, we may spare unnecessary computations by specifying analytically
    the support over which it will be non-negligible.
    Given a floating-point number `ɛ` (defaulting to machine precision
    `eps(T)`), the support above `ɛ` of the Gaussian `g(ω) of variance `σ` is
    defined by
        `g(ω) = exp(- ω² / (2σ²)) = 1/ɛ`,
    which leads to
        `ω² = log(ɛ) * 2σ²`.
    Recalling that `bw² = 4 * log(2) * 2σ²` (see paragraph 1. above), we have
        `ω² = log(ɛ) * bw² / (4 * log(2))`, and then
        `ω  = ± bw/2 * sqrt(log(ɛ)/log(2)) = ± bw/2 * sqrt(log2(ɛ))
    Let ρ be the constant bw/2 * sqrt(log2(ɛ)).
    The ɛ bounds of each of the three terms are:
        1. main Gaussian: c ± ρ
        2. corr0 Gaussian: ± ρ*(1+sqrt(-log2(corr0)))
        3. corrN Gaussian: N ± ρ*(1+sqrt(-log2(corrN)))
    By the triangular inequality, we take the union of those three intervals
    to get a conservative superset of the ɛ support of the Morlet wavelet.
    Finally, we bound the obtained set by -N/2 and 3N/2, in order to compute
    at most 2N Fourier-domain Morlet coefficients."""
    if spec.ɛ == 0.0
        first = -half_length + 1
        last = 3half_length
    else
        ρ = @fastmath 0.5 * bw * sqrt(-log2(spec.ɛ))
        first = center - ρ
        last = center + ρ
        if corr0 != 0
            corr0_first = - ρ * (1 + sqrt(-log2(corr0)))
            corr0_last = - corr0_first
            first = min(first, corr0_first)
            last = max(last, corr0_last)
        end
        if corrN != 0
            corrN_first = N - ρ * (1 + sqrt(-log2(corrN)))
            corrN_last = 2N - corrN_first
            first = min(first, corrN_first)
            last = max(last, corrN_last)
        end
        first = floor(Int, min(first, -half_length + 1))
        last = ceil(Int, max(last, 3half_length))
    end
    "4. **Computational comprehension of the Morlet 1D wavelet**"
    @inbounds ωs = convert(Vector{T}, collect(first:last))
    @fastmath @inbounds morlet =
        T[ morletfourier1d(ω, center, den, N, corr0, corrN) for ω in ωs ]
    """5. **Trimming to true support boundaries**
    We look for the true ɛ boundaries of the vector above by looking
    at the first (resp. last) coefficient for which `|ψ|²(ω) > ɛ²`."""
    ɛ2 = T(spec.ɛ * spec.ɛ)
    morlet2 = abs2(morlet)
    sub_first = findfirst(morlet2 .> ɛ2)
    sub_last = findlast(morlet2 .> ɛ2)
    morlet = morlet[sub_first:sub_last]
    first = first + (sub_first-1)
    last = last - (length(morlet)-sub_last)
    "6. **Construction of AbstractFourier1DFilter object**"
    AbstractFourier1DFilter(morlet, first, last, log2_length)
end

"""Computes one coefficient of the Morlet wavelet of center frequency `center`
and denominator `den = 2σ²` in the frequency domain, at the frequency ω.
`N` is the signal length. `corr0` and `corrN` are the Morlet low-frequency
corrective terms (computed by the `fourierwavelet` function above)."""
morletfourier1d(ω, center, den, N, corr0, corrN) =
    gauss(ω-center, den) - corr0 * gauss(ω, den) - corrN * gauss(ω-N, den)

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
