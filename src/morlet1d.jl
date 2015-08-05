immutable Morlet1DSpec{T<:Number} <: Abstract1DSpec{T}
    ɛ::Float64
    log2_length::Int
    max_qualityfactor::Float64
    max_scale::Float64
    motherfrequency::Float64
    nFilters_per_octave::Int
    nOctaves::Int
    signaltype::Type{T}
    function call{T<:Number}(::Type{Morlet1DSpec{T}}, signaltype::Type{T};
                             ɛ=default_ɛ(T), log2_length=15,
                             max_qualityfactor=nothing, max_scale=Inf,
                             mother_centerfrequency=nothing,
                             nFilters_per_octave=nothing, nOctaves=nothing)
        max_qualityfactor, nFilters_per_octave =
            default_max_qualityfactor(max_qualityfactor, nFilters_per_octave),
            default_nFilters_per_octave(max_qualityfactor, nFilters_per_octave)
        motherfrequency =
            default_motherfrequency(Morlet1DSpec, nFilters_per_octave)
        nOctaves =
            default_nOctaves(Morlet1DSpec, nOctaves, log2_length,
                             max_qualityfactor, max_scale, nFilters_per_octave)
        checkspec(ɛ, log2_length, max_qualityfactor,
                  max_scale, nFilters_per_octave, nOctaves)
        new{T}(ɛ, log2_length, max_qualityfactor,
               max_scale, nFilters_per_octave, nOctaves, signaltype)
    end
end

"By default, `Morlet1DSpec operates on single-precision real input (Float32)."
Morlet1DSpec(T=Float32; args...) = Morlet1DSpec{T}(T; args...)

"""Returns the maximal number octaves in a Morlet 1d filter bank such that all
wavelet scales are below 2^(log2_length)."""
function default_nOctaves(::Type{Morlet1DSpec}, nOctaves::Void, log2_length,
                          max_qualityfactor, max_scale, nFilters_per_octave)
    log2_nFilters_per_octave = ceil(Int, log2(nFilters_per_octave))
    log2_max_qualityfactor = ceil(Int, log2(max_qualityfactor))
    if max_scale < (exp2(log2_length)+default_ɛ(T))
        gap = max(1+log2_nFilters_per_octave, 2)
    else
        gap = max(1+log2_nFilters_per_octave, 2+log2_max_qualityfactor)
    end
    nOctaves = log2_length - gap
end

"""The dimensionless mother center frequency ξ (corresponding to a log-period
γ=0) is computed as the midpoint between the center frequency of the second
center frequency ξ*2^(-1/nFilters_per_octave) (corresponding to γ=1) and the
negative mother center frequency (1-ξ).

Hence the equation 2ξ = ξ*2^(-1/nFilters_per_octave) + (1-ξ), of which we
derive ξ = 1 / (3.0 - 2^(1/nFilters_per_octave)).

In the special case `nFilters_per_octave=1`, we manually set `ξ=0.39`, which is
more accurate with the Littlewood-Paley energy conservation criterion than the
formula `ξ=0.4`."""
function default_motherfrequency(::Type{Morlet1DSpec}, nFilters_per_octave)
    nFilters_per_octave==1 && return 0.39
    return inv(3.0 - exp2(-1.0/nFilters_per_octave))
end

"""Computes in closed form the bandwidths, center frequencies, quality factors,
and scales of all wavelets in a given Morlet 1d filter bank.

The dimensionless center frequencies (between 0.0 and 0.5) are of the form
`ω = ξ*2^(-γs/spec.nFilters_per_octave)`, where `ξ` is the dimensionless mother
frequency and the `γs` are natural integers.

There is a classical tradeoff between spatial and frequential localizations
in a filter bank. We address it by supporting two user specifications
* spatial localization: `max_scale` sets the maximal wavelet scale, in the sense
  of squared-magnitude full width at tenth maximum (FWTM).
* frequential localization: `max_qualityfactor` set the quality factor (ratio
  between center frequency and 3dB bandwidth) in absence of spatial localization
  constraints.

For each center frequency, the quality factor and the scale are governed by the
following criteria, in decreasing priority order:
1. quality factor is equal or greater than 1.0
2. scale is equal or smaller than max_scale
3. quality factor is equal to max_qualityfactor

The bandwidth of a wavelet is defined as the full width at half maximm (FWHM)
of the squared magnitude in the Fourier domain, i.e. the 3dB bandwidth.
By neglecting the low-frequency terms, we write the Morlet wavelet as a
Gaussian of variance σ in the Fourier domain. Its bandwidth (FWHM) is then
equal to 2σ*sqrt(log(2)).

Therefore, for a given center frequency ω and a quality factor Q, the variance
σ of the Gaussian is equal to σ = ω / (2*Q*sqrt(log(2))). In the spatial domain,
this amounts to a Gabor wavelet (a Gaussian modulated by a sine wave, without
any low-frequency corrective term) of variance 1/σ. Its spatial scale (FWTM)
is equal to s = 2*sqrt(log(10))/σ. We conclude that the quality factor Q and the
scale s are linked by the equation:
    s = sqrt(log(10)/log(2)) * 1/Q.

To localize Morlet wavelets according to user-defined max_qualityfactor and
max_scale, we proceed with the following steps:
1. compute ""unbounded scales"" sqrt(log(10)/log(2)) * 1/max_qualityfactor.
2. bound scales from above s = min(unbounded_scales, max_scale)
3. compute ""unbounded quality factors"" sqrt(log(2)/log(10)) * 1/s
4. bound quality factors from below q = max(unbounded_qualityfactors, 1.0)
5. compute corresponding scales."""
function localize{T<:Number}(spec::Morlet1DSpec{T})
    RealT = realtype(T)
    γs = gammas(spec)
    resolutions = exp2(γs/spec.nFilters_per_octave))
    centerfrequencies = spec.mother_centerfrequency * resolutions
    scale_multiplier = sqrt(log(10.0)/log(2.0))
    unbounded_scales =
        scale_multiplier * (spec.max_qualityfactor./centerfrequencies)
    scales = min(unbounded_scales, spec.max_scale)
    unbounded_qualityfactors = scales .* centerfrequencies / scale_multiplier
    # we also bound qualityfactors from above for better numerical accuracy
    qualityfactors = clamp(unbounded_qualityfactors, 1.0, spec.max_qualityfactor)
    bandwidths = resolutions ./ qualityfactors
    scales = scale_multiplier * qualityfactors./centerfrequencies
    return (bandwidths, centerfrequencies, qualityfactors, scales)
end
