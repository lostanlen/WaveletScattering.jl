immutable Morlet1DSpec{T<:Number} <: Abstract1DSpec{T}
    ɛ::Float64
    log2_size::Tuple{Int}
    max_qualityfactor::Float64
    max_scale::Float64
    motherfrequency::Float64
    nFilters_per_octave::Int
    nOctaves::Int
    signaltype::Type{T}
    function call{T<:Number}(::Type{Morlet1DSpec{T}}, signaltype::Type{T};
                             ɛ=default_ɛ(T), log2_length=15,
                             max_qualityfactor=nothing, max_scale=Inf,
                             nFilters_per_octave=nothing, nOctaves=nothing,
                             tuningfrequency=nothing)
        isa(log2_size, Int) && log2_size = tuple(log2_size)
        max_qualityfactor, nFilters_per_octave =
            default_max_qualityfactor(max_qualityfactor, nFilters_per_octave),
            default_nFilters_per_octave(max_qualityfactor, nFilters_per_octave)
        motherfrequency = tune_motherfrequency(tuningfrequency, Morlet1DSpec,
                                               nFilters_per_octave)
        nOctaves =
            default_nOctaves(nOctaves, Morlet1DSpec, log2_size,
                             max_qualityfactor, max_scale, nFilters_per_octave)
        new{T}(ɛ, log2_length, max_qualityfactor, max_scale, motherfrequency,
               nFilters_per_octave, nOctaves, signaltype)
        checkspec(spec)
    end
end

"By default, `Morlet1DSpec operates on single-precision real input (Float32)."
Morlet1DSpec(T=Float32; args...) = Morlet1DSpec{T}(T; args...)

"""Returns the maximal number octaves in a Morlet 1d filter bank such that all
wavelet scales are below 2^(log2_length)."""
function default_nOctaves(nOctaves::Void, ::Type{Morlet1DSpec}, log2_length,
    max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave)

end


"""In the special case `nFilters_per_octave=1`, we manually set `ξ=0.39`, which
is more accurate with the Littlewood-Paley energy conservation criterion than
the generic fallback `ξ=0.4`, which is only valid when the wavelet has a
symmetric profile in the Fourier domain. This is no longer the case for
nFilters_per_octave==max_qualityfactor==1 as the Morlet low-frequency corrective
term is no longer negligible."""
default_motherfrequency(::Type{Morlet1DSpec}, nFilters_per_octave) =
    nFilters_per_octave==1 ? 0.39 : inv(3.0 - exp2(-1.0/nFilters_per_octave))

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

function qualityfactors(spec::AbstractSpec)
    γs = gammas(spec)
    centerfrequencies = centerfrequencies(spec)
    heisenberg = heisenberg(spec)
    bandwidths = centerfrequencies/spec.max_qualityfactor
    scales = heisenberg * bandwidths
    scales = min(scales, spec.max_scale)
    qualityfactors = scales .* centerfrequencies / heisenberg
    # we also bound qualityfactors from above for better numerical accuracy
    qualityfactors = clamp(qualityfactors, 1.0, spec.max_qualityfactor)
end
bandwidths(spec::AbstractSpec) =
    bandwidths = centerfrequencies ./ qualityfactors
    scales = heisenberg / bandwidths
    return (bandwidths, centerfrequencies, qualityfactors, scales)
end

centerfrequencies(spec::AbstractSpec) =
    exp2(gammas(spec)/spec.nFilters_per_octave)

heisenberg(spec::Morlet1DSpec) = sqrt(log(10.0) / log(2.0))
