"""A `ψMeta` object contains all the item-information to identify
an oriented wavelet within a filter bank. Fields:
* `γ` log-scale. `2^(-γ)` is proportional to center frequency
* `θ` orientation, i.e. angle (in 2d) or sign of center frequency (in 1d)
* `χ` chroma. `χ` is `mod(γ, nFilters_per_octave)`
* `j` octave. `j` is `div(γ, nFilters_per_octave)`
* `bandwidth` ∈]0,1] the width at -3dB, expressed in fraction of signal length
* `centerfrequency` ∈]0,1] expressed in fraction of signal length
* `qualityfactor` ∈[1,max_qualityfactor] is equal to `centerfrequency/bandwidth`
* `scale` is the FWTM (full width at tenth maximum) in spatial domain"""
immutable ΨMeta
    γ::Int16
    θ::Int8
    χ::Int8
    bandwidth::Float64
    centerfrequency::Float64
    j::Int8
    qualityfactor::Float64
    scale::Float64
end

immutable ΦMeta
    bandwidth::Float64
    scale::Float64
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
function waveletmetas(spec::Spec1D)
    nΘs = get_nΘs(spec.pointgroup)
    ħ = uncertainty(spec)
    ψmetas = Matrix{ΨMeta}(spec.nOctaves, spec.nFilters_per_octave, nΘs)
    for j in 0:(spec.nOctaves-1)
        for χ in 0:(spec.nFilters_per_octave-1)
            γ = j * nFilters_per_octave + χ
            resolution = exp2(-γ / spec.nFilters_per_octave)
            centerfrequency = spec.motherfrequency * resolution
            unbounded_scale = ħ * spec.max_qualityfactor / centerfrequency
            scale = min(unbounded_scale, spec.max_scale)
            unbounded_qualityfactor = scale * centerfrequency / ħ
            qualityfactor =
                clamp(unbounded_qualityfactor, 1.0, spec.max_qualityfactor)
            bandwidth = centerfrequency / qualityfactor
            for θ in 0:(spec.nΘs-1)
                ψmetas[1+j, 1+χ, 1+θ] = ΨMeta(γ, θ, χ, bandwidth,
                    centerfrequency, j, qualityfactor, scale)
            end
        end
    end
end

"""Fallback of the uncertainty constant from the spec to its class. The RHS
method must be specifically implemented by AbstractSpec concrete subtypes."""
uncertainty(spec::AbstractSpec) = uncertainty(spec.class)
