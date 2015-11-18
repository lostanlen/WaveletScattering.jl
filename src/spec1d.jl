immutable Spec1D{T<:FFTW.fftwReal,D<:LineDomains,
        G<:LineGroups,W<:RedundantWaveletClass} <: AbstractSpec{T,D,G}
    ɛ::Float64
    ϕmeta::ΦMeta
    ψmetas::Array{ΨMeta,3}
    class::W
    domain::D
    log2_size::Tuple{Int}
    max_qualityfactor::Float64
    max_scale::Float64
    motherfrequency::Float64
    nFilters_per_octave::Int
    nOctaves::Int
    pointgroup::G
    signaltype::Type{T}
    function call{T,D,G,W}(
            ::Type{Spec1D} ;
            class::W = Morlet(),
            pointgroup::G = TrivialGroup(),
            signaltype::Type{T} = Float32,
            domain::D = FourierDomain(1),
            ɛ = default_ɛ(signaltype),
            log2_size = 15,
            max_qualityfactor = nothing,
            max_scale = Inf,
            nFilters_per_octave = nothing,
            nOctaves = nothing,
            tuningfrequency = nothing)
        "Integer log2_size is automatically converted to one-element tuple"
        isa(log2_size, Int) && (log2_size = tuple(log2_size))
        max_qualityfactor, nFilters_per_octave =
             default_max_qualityfactor(max_qualityfactor, nFilters_per_octave),
             default_nFilters_per_octave(nFilters_per_octave, max_qualityfactor)
        motherfrequency =
            tune_motherfrequency(tuningfrequency, class, nFilters_per_octave)
        nOctaves = default_nOctaves(nOctaves, class, log2_size,
            max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave)
        ψmetas = waveletmetas(spec)
        ϕmeta = lowpassmeta(spec)
        spec = new{T,D,G,W}(ɛ, class, domain, log2_size, max_qualityfactor,
            max_scale, motherfrequency, nFilters_per_octave, nOctaves,
            pointgroup, signaltype)
        checkspec(spec) && return spec
    end
end

"""
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
