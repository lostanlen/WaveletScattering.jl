immutable Spec2D{T<:Real,D<:PlaneDomains,
        G<:PlaneGroups,W<:RedundantWaveletClass} <: AbstractSpec{T,D,G,W}
    ɛ::Float64
    ϕmeta::ΦMeta
    ψmetas::Array{ΨMeta2D,3}
    class::W
    domain::D
    log2_size::Int
    max_aspectratio::Float64
    max_qualityfactor::Float64
    max_scale::Float64
    motherfrequency::Float64
    n_filters_per_octave::Int
    n_octaves::Int
    pointgroup::G
    signaltype::Type{T}

    function (::Type{Spec2D}){T,D,W}( ;
            class::W = Morlet(),
            signaltype::Type{T} = Float32,
            domain::D = FourierDomain(2),
            ɛ = default_ɛ(signaltype),
            log2_size = 5,
            max_aspectratio = nothing,
            max_qualityfactor = nothing,
            max_scale = Inf,
            n_filters_per_octave = nothing,
            n_octaves = nothing,
            n_orientations = nothing)
        n_orientations = default_n_orientations(class, n_orientations)
        """Infer point group from wavelet."""
        pointgroup = issteerable(class) ?
            RotationGroup(n_orientations) : TrivialGroup()
        """Maximum aspect ratio (length-to-width) of the wavelets."""
        max_aspectratio = default_max_aspectratio(class, max_aspectratio)
        """Maximum quality factor and number of filters per octave are
        set simultaneously, so that they default to each other if either
        is present."""
        max_qualityfactor, n_filters_per_octave =
            default_max_qualityfactor(max_qualityfactor, n_filters_per_octave),
            default_n_filters_per_octave(n_filters_per_octave, max_qualityfactor)
        """Center frequency of the mother wavelet (`motherfrequency`)."""
        motherfrequency = default_motherfrequency(class, n_filters_per_octave)
        """By default, the number of octaves in the filter bank is such that
        all wavelet scales do not exceed the signal length, nor the
        user-defined maximum scale `max_scale` if present."""
        n_octaves = default_n_octaves(n_octaves, class, log2_size,
            max_qualityfactor, max_scale, motherfrequency, n_filters_per_octave)
        """The number of orientations of a one-dimensional filter bank is
        equal to 1 for non-steerable real wavelets (`pointgroup::TrivialGroup`),
        equal to n_orientations otherwise (`pointgroup::RotationGroup`)."""
        nΘs = get_n_orientations(pointgroup)
        ħ = uncertainty(class)
        """The meta-information of the wavelets is computed as a 3-d array,
        whose dimensions respectively correspond to
        * the orientation variable `θ`,
        * the chroma variable `χ`, and
        * the octave variable `j`."""
        ψmetas = Array{ΨMeta2D}(nΘs, n_filters_per_octave, n_octaves)
        for j in 0:(n_octaves-1), χ in 0:(n_filters_per_octave-1)
            γ = j * n_filters_per_octave + χ
            resolution = exp2(-γ / n_filters_per_octave)
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
                ψmetas[1+θ, 1+χ, 1+j] = ΨMeta2D(γ, θ, χ,
                    aspectratio, bandwidth, centerfrequency, j,
                    qualityfactor, scale)
            end
        end
        """The bandwidth of the lowpass filter `ϕ` is such that the Fourier
        transform of `ϕ` will cover the low-frequency part of the spectrum
        that is not covered by the wavelets `ψs`."""
        ϕbandwidth = motherfrequency * exp2(-n_octaves)
        ϕscale = ħ / ϕbandwidth
        ϕmeta = ΦMeta(ϕbandwidth, ϕscale)
        G = typeof(pointgroup)
        spec = new{T,D,G,W}(ɛ, ϕmeta, ψmetas, class, domain, log2_size,
            max_aspectratio, max_qualityfactor, max_scale, motherfrequency,
            n_filters_per_octave, n_octaves, pointgroup, signaltype)
        """Before returning the `spec`, we call the function `checkspec`
        which enforces properties of the wavelet filter bank to satisfy null
        mean, limited spatial support, and Littlewood-Paley inequality."""
        checkspec(spec) && return spec
    end
end

"""Enforces properties of two-dimensional wavelets.

* For steerable wavelets, such as `Morlet`, `max_aspectratio` must be `≧1.0` and
`n_orientations` must be `≧2`.

* For non-steerable wavelets, such as `MexicanHat`, `max_aspectratio` must be
equal to `1.0` and `n_orientations` must be equal to `1`."""
function checkspec(spec::Spec2D)
    if issteerable(spec.class)
        if spec.max_aspectratio < 1.0
            error("`max_aspectratio` must be `≧1.0`",
                "for steerable wavelets.")
        end
        if get_n_orientations(spec.pointgroup) < 2
            error("`n_orientations` must be `≧2` for steerable wavelets")
        end
    else
        if spec.max_aspectratio != 1.0
            error("`max_aspectratio` must be equal to `1.0`",
                "for non-steerable wavelets.")
        end
        if get_n_orientations(spec.pointgroup) != 1
            error("`n_orientations` must be equal to `1`",
                "for non-steerable wavelets")
        end
    end
    checkspec_super(spec) && return true
end

"""By default, the maximum aspect ratio is set to `2.0` for steerable
wavelets, such as `Morlet`. For non-steerable wavelets, such as `MexicanHat`,
the maximum aspect ratio must be equal to `1.0`."""
default_max_aspectratio(class::RedundantWaveletClass, max_aspectratio::Void) =
    issteerable(class) ? 2.0 : 1.0
default_max_aspectratio(class::RedundantWaveletClass, max_aspectratio::Any) =
    Float64(max_aspectratio)

"""By default, the number of orientations is set to `4` for steerable
wavelets, such as `Morlet`. For non-steerable wavelets, such as `MexicanHat`,
the number of orientations must be equal to `1`."""
default_n_orientations(class::RedundantWaveletClass, n_orientations::Void) =
    issteerable(class) ? 4 : 1

function default_n_orientations(class::RedundantWaveletClass, n_orientations::Any)
    if issteerable(class)
        return Int(n_orientations)
    else
        error("`n_orientations` must be equal to `1` for steerable wavelets")
    end
end
