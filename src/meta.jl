"""A `Meta` object contains all the meta-information to identify a wavelet
within a filter bank."""
abstract AbstractMeta

"""A `NonOrientedMeta` object contains all the meta-information to identify a
non-oriented wavelet within a filter bank. Fields:
* γ log-scale. 2^(-γ) is proportional to center frequency
* θ orientation. θ is proportional to angle (2d) or spin (1d)
* χ chroma. χ is mod(γ, nFilters_per_octave)
* j octave. j is div(γ, nFilters_per_octave)
* bandwidth ∈]0,1] is the width at -3dB, expressed in fraction of signal length
* centerfrequency ∈]0,1] is expressed in fraction of signal length
* qualityfactor ∈[1,max_qualityfactor] is equal to centerfrequency/bandwidth
* scale is the FWTM (full width at tenth maximum) in spatial domain"""
immutable NonOrientedMeta <: AbstractMeta
    γ::Int16
    χ::Int8
    bandwidth::Float64
    centerfrequency::Float64
    j::Int8
    qualityfactor::Float64
    scale::Float64
end

"""Returns the 3dB bandwidths, i.e. the full widths at half maximum (FWHM) of
the squared magnitude in the Fourier domain, of a given spec.
Bandwidths are decreasing because they are indexed by `γ`"""
bandwidths(spec::AbstractSpec) = centerfrequencies(spec) ./ qualityfactors(spec)

"""Returns the center frequencies of a given spec. They are exponentially
decreasing because they are indexed by `γ`. The first coefficient corresponds
to the so-called ""mother"" frequency, i.e. γ=0."""
centerfrequencies(spec::AbstractSpec) =
    spec.motherfrequency * exp2(-gammas(spec)/spec.nFilters_per_octave)

"""Returns the chroma indices `χs`, i.e. locations within the octave, of a
wavelet spec. Chroma indices range from `0` to `nFilters_per_octave-1`. The
convention is that higher chroma indices `χs` mean *lower* center frequencies.
Log-periods `γs`, chromas ``χs`, and octaves `js` are linked by
    γ = j + nFilters_per_octave * χ"""
function chromas(spec::AbstractSpec)
    repmat(collect(0:(spec.nFilters_per_octave-1)), spec.nOctaves)
end

"""Returns the wavelet log-period integer indices `γs`. Center frequencies are
proportional to 2^(-γ). γ ranges from 0 to nFilters_per_octave*nOctaves, where
γ=0 corresponds to the mother frequency. The convention is that higher indices
`γs` mean *lower* center frequencies. Log-periods γs`, chromas ``χs`, and
octaves `js` are linked by
    γ = j + nFilters_per_octave * χ"""
gammas(spec::AbstractSpec) =
    collect(0:(spec.nFilters_per_octave * spec.nOctaves-1))

"""Returns the octave indices js of a wavelet spec.
Octave indices range from `0` to `nOctaves-1`. The convention is that
higher octave indices `js` mean *lower* center frequencies. Log-periods
``γs`, chromas ``χs`, and octaves `js` are linked by
    γ = j + nFilters_per_octave * χ"""
function octaves(spec::AbstractSpec)
    vec(repmat(transpose(collect(0:(spec.nOctaves-1))), spec.nFilters_per_octave))
end

"""Returns the quality factors (ratios of center frequencies over bandwidths).

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

To localize Morlet wavelets according to user-defined max_qualityfactor and
max_scale, we proceed with the following steps:
1. compute ""unbounded scales""  h/max_qualityfactor.
2. bound scales from above s = min(unbounded_scales, max_scale)
3. compute ""unbounded quality factors"" 1/(h*s)
4. bound quality factors from below q = max(unbounded_qualityfactors, 1.0)
5. compute corresponding scales.
"""
function qualityfactors(spec::AbstractSpec)
    bandwidths = centerfrequencies(spec)/spec.max_qualityfactor
    scales = min(uncertainty(spec) ./ bandwidths, spec.max_scale)
    qualityfactors = scales .* centerfrequencies(spec) / uncertainty(spec)
    # we also bound qualityfactors from above for better numerical accuracy
    qualityfactors = clamp(qualityfactors, 1.0, spec.max_qualityfactor)
end

"""Returns the scales of a wavelet spec, defined as the full width at tenth
maximum (FWTM) of the squared-magnitude spatial support."""
scales(spec::AbstractSpec) = uncertainty(spec) ./ bandwidths(spec)

"""Fallback of the uncertainty constant from the spec to its type. The RHS
method must be specifically implemented by AbstractSpec concrete subtypes."""
uncertainty(spec::AbstractSpec, args...) = uncertainty(typeof(spec), args...)
