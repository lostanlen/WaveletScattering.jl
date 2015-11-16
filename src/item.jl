"""A `Meta` object contains all the meta-information to identify a wavelet
within a filter bank."""
abstract AbstractMeta{G<:AbstractPointGroup}

immutable MetaVector <: AbstractMetaVector{TrivialGroup}
    metas::Vector{NonOrientedMeta}
end



"""A `NonOrientedMeta` object contains all the meta-information to identify a
non-oriented wavelet within a filter bank. Fields:
* `γ` log-scale. `2^(-γ)` is proportional to center frequency
* `χ` chroma. `χ` is `mod(γ, nFilters_per_octave)`
* `j` octave. `j` is `div(γ, nFilters_per_octave)`
* `bandwidth` ∈]0,1] the width at -3dB, expressed in fraction of signal length
* `centerfrequency` ∈]0,1] is expressed in fraction of signal length
* `qualityfactor` ∈[1,max_qualityfactor] is equal to `centerfrequency/bandwidth`
* `scale` is the FWTM (full width at tenth maximum) in spatial domain"""
immutable NonOrientedMeta
    γ::Int16
    χ::Int8
    bandwidth::Float64
    centerfrequency::Float64
    j::Int8
    qualityfactor::Float64
    scale::Float64
end

"""An `OrientedMeta` object contains all the meta-information to identify
an oriented wavelet within a filter bank. Fields:
* `γ` log-scale. `2^(-γ)` is proportional to center frequency
* `θ` orientation, i.e. angle (in 2d) or sign of center frequency (in 1d)
* `χ` chroma. `χ` is `mod(γ, nFilters_per_octave)`
* `j` octave. `j` is `div(γ, nFilters_per_octave)`
* `bandwidth` ∈]0,1] the width at -3dB, expressed in fraction of signal length
* `centerfrequency` ∈]0,1] expressed in fraction of signal length
* `qualityfactor` ∈[1,max_qualityfactor] is equal to `centerfrequency/bandwidth`
* `scale` is the FWTM (full width at tenth maximum) in spatial domain"""
immutable OrientedMeta{G<:Union{ReflectionGroup,RotationGroup}}
    γ::Int16
    θ::Int8
    χ::Int8
    bandwidth::Float64
    centerfrequency::Float64
    j::Int8
    pointgroup::G
    qualityfactor::Float64
    scale::Float64
end

"""Returns the 3dB bandwidths, i.e. the full widths at half maximum (FWHM) of
the squared magnitude in the Fourier domain, of a given spec.
Bandwidths are decreasing because they are indexed by `γ`"""
bandwidths(spec::AbstractSpec) = centerfrequencies(spec) ./ qualityfactors(spec)

"""Returns the center frequencies of a given spec. They are exponentially
decreasing because they are indexed by `γ`. The first coefficient corresponds
to the so-called ""mother"" frequency, i.e. `γ=0`."""
centerfrequencies(spec::AbstractSpec) =
    spec.motherfrequency * exp2(-gammas(spec)/spec.nFilters_per_octave)

"""Returns the chroma indices `χs`, i.e. locations within the octave, of a
wavelet spec. Chroma indices range from `0` to `nFilters_per_octave-1`. The
convention is that higher chroma indices `χs` mean *lower* center frequencies.
Log-periods `γs`, chromas `χs`, and octaves `js` are linked by
    γ = j + nFilters_per_octave * χ"""
function chromas(spec::AbstractSpec)
    repmat(collect(0:(spec.nFilters_per_octave-1)), spec.nOctaves)
end

"""Returns the wavelet log-period integer indices `γs`. Center frequencies are
proportional to `2^(-γ)`. `γ` ranges from 0 to `nFilters_per_octave*nOctaves`,
where `γ = 0` corresponds to the mother frequency. The convention is that higher
indices `γs` mean *lower* center frequencies. Log-periods γs`, chromas `χs`, and
octaves `js` are linked by
    `γ = j + nFilters_per_octave * χ`"""
gammas(spec::AbstractSpec) =
    collect(0:(spec.nFilters_per_octave * spec.nOctaves-1))

"""Returns the octave indices `js` of a wavelet spec.
Octave indices range from `0` to `nOctaves-1`. The convention is that
higher octave indices `js` mean *lower* center frequencies. Log-periods
`γs`, chromas `χs`, and octaves `js` are linked by
    `γ = j + nFilters_per_octave * χ`"""
function octaves(spec::AbstractSpec)
    vec(repmat(transpose(collect(0:(spec.nOctaves-1))),
        spec.nFilters_per_octave))
end

"""Returns the quality factors (ratios of center frequencies over bandwidths).

There is a classical tradeoff between spatial and frequential localizations
in a filter bank. We address it by supporting two user specifications:
* spatial localization: `max_scale` sets the maximal wavelet scale, in the sense
  of squared-magnitude full width at tenth maximum (FWTM).
* frequential localization: `max_qualityfactor` sets the quality factor (ratio
  between center frequency and 3dB bandwidth) in absence of spatial localization
  constraints.

For each center frequency, the quality factor and the scale are governed by the
following criteria, in decreasing priority order:
1. quality factor is equal or greater than `1.0`
2. scale is equal or smaller than `max_scale`
3. quality factor is equal to `max_qualityfactor`

To localize Morlet wavelets according to user-defined max_qualityfactor and
max_scale, we proceed with the following steps:
1. compute center frequencies `ξs` and uncertainty `h`,
2. compute unbounded scales `max_qualityfactor/(h*ξs)`,
3. bound scales from above by `max_scale`,
4. compute unbounded quality factors `scales.*ξs/h`, and
5. bound quality factors from below by `1.0`.
"""
function qualityfactors(spec::AbstractSpec)
    h = uncertainty(spec)
    ξs = centerfrequencies(spec)
    unbounded_scales = h * spec.max_qualityfactor ./ ξs
    scales = min(unbounded_scales, spec.max_scale)
    unbounded_qualityfactors = scales .* ξs / h
    # we also clamp qualityfactors from above for better numerical accuracy
    return clamp(unbounded_qualityfactors, 1.0, spec.max_qualityfactor)
end

"""Returns the scales of a wavelet spec, defined as the full width at tenth
maximum (FWTM) of the squared-magnitude spatial support."""
scales(spec::AbstractSpec) = uncertainty(spec) ./ bandwidths(spec)

"""Fallback of the uncertainty constant from the spec to its type. The RHS
method must be specifically implemented by AbstractSpec concrete subtypes."""
uncertainty(spec::AbstractSpec, args...) = uncertainty(typeof(spec), args...)
