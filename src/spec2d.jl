immutable Spec2D{T<:Real,D<:PlaneDomains,
        G<:PlaneGroups,W<:RedundantWaveletClass} <: AbstractSpec{T,D,G,W}
    ɛ::Float64
    ϕmeta::ΦMeta
    ψmetas::Array{ΨMeta,3}
    class::W
    domain::D
    log2_size::Int
    max_aspectratio::Float64
    max_qualityfactor::Float64
    max_scale::Float64
    motherfrequency::Float64
    nFilters_per_octave::Int
    nOctaves::Int
    nOrientations::Int
    pointgroup::G
    signaltype::Type{T}

    function call{T,D,W}(
            ::Type{Spec2D} ;
            class::W = Morlet(),
            signaltype::Type{T} = Float32,
            domain::D = FourierDomain(2),
            ɛ = default_ɛ(signaltype),
            log2_size = 5,
            max_aspectratio = nothing,
            max_qualityfactor = nothing,
            max_scale = Inf,
            nFilters_per_octave = nothing,
            nOctaves = nothing,
            nOrientations = 4)
        """Infer point group from wavelet"""
        G = issteerable(class) ? RotationGroup() : TrivialGroup()
        """Maximum aspect ratio (length-to-width) of the wavelets"""
        max_aspectratio = default_max_aspectratio(class, max_aspectratio)
        """Maximum quality factor and number of filters per octave are
        set simultaneously, so that they default to each other if either
        is present."""
        max_qualityfactor, nFilters_per_octave =
             default_max_qualityfactor(max_qualityfactor, nFilters_per_octave),
             default_nFilters_per_octave(nFilters_per_octave, max_qualityfactor)
        """Center frequency of the mother wavelet (`motherfrequency`)."""
        motherfrequency = default_motherfrequency(class, nFilters_per_octave)
        """By default, the number of octaves in the filter bank is such that
        all wavelet scales do not exceed the signal length, nor the
        user-defined maximum scale `max_scale` if present."""
        nOctaves = default_nOctaves(nOctaves, class, log2_size,
            max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave)
        """The number of orientations of a one-dimensional filter bank is
        equal to 1 for non-steerable real wavelets (`pointgroup::TrivialGroup`), equal to  (`pointgroup::ReflectionGroup`)."""
        nΘs = get_nOrientations(pointgroup)
        ħ = uncertainty(class)
        """The meta-information of the wavelets is computed as a 3-d array,
        whose dimensions respectively correspond to
        * the orientation variable `θ`,
        * the chroma variable `χ`, and
        * the octave variable `j`."""
        ψmetas = Array{ΨMeta}(nΘs, nFilters_per_octave, nOctaves)
        for j in 0:(nOctaves-1), χ in 0:(nFilters_per_octave-1)
            γ = j * nFilters_per_octave + χ
            resolution = exp2(-γ / nFilters_per_octave)
            centerfrequency = motherfrequency * resolution
            unbounded_scale =
                ħ * max_aspectratio * max_qualityfactor / centerfrequency
            scale = min(unbounded_scale, max_scale)
            unbounded_q = scale * centerfrequency / (ħ * max_aspectratio)
            qualityfactor = clamp(unbounded_q, 1.0, max_qualityfactor)
            unbounded_scale =
                ħ * max_aspectratio * qualityfactor / centerfrequency
            scale = min(unbounded_scale, max_scale)
            aspectratio = scale * centerfrequency / (ħ * scale)
            bandwidth = centerfrequency / qualityfactor
            for θ in 0:(nΘs-1)
                ψmetas[1+θ, 1+χ, 1+j] = ΨMeta(γ, θ, χ, bandwidth,
                    centerfrequency, j, qualityfactor, scale)
            end
        end
        """The bandwidth of the lowpass filter `ϕ` is such that the Fourier
        transform of `ϕ` will cover the low-frequency part of the spectrum
        that is not covered by the wavelets `ψs`."""
        ϕbandwidth = motherfrequency * exp2(-nOctaves)
        ϕscale = ħ / ϕbandwidth
        ϕmeta = ΦMeta(ϕbandwidth, ϕscale)
        spec = new{T,D,G,W}(ɛ, ϕmeta, ψmetas, class, domain, log2_size,
            max_aspectratio, max_qualityfactor, max_scale, motherfrequency,
            nFilters_per_octave, nOctaves, pointgroup, signaltype)
        """Before returning the `spec`, we call the function `checkspec`
        which enforces properties of the wavelet filter bank to satisfy null
        mean, limited spatial support, and Littlewood-Paley inequality."""
        # checkspec(spec) && return spec
        return spec
    end
end
