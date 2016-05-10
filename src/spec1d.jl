immutable Spec1D{T<:Real,D<:LineDomains,
        G<:LineGroups,W<:RedundantWaveletClass} <: AbstractSpec{T,D,G,W}
    ɛ::Float64
    ϕmeta::ΦMeta
    ψmetas::Array{ΨMeta1D,3}
    class::W
    domain::D
    log2_size::Int
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
        """Maximum quality factor and number of filters per octave are
        set simultaneously, so that they default to each other if either
        is present."""
        max_qualityfactor, nFilters_per_octave =
             default_max_qualityfactor(max_qualityfactor, nFilters_per_octave),
             default_nFilters_per_octave(nFilters_per_octave, max_qualityfactor)
        """The center frequency of the mother wavelet (`motherfrequency`)
        may be tuned to user-defined frequency (e.g. 440 Hz in music) or left
        as a default."""
        motherfrequency =
            tune_motherfrequency(tuningfrequency, class, nFilters_per_octave)
        """By default, the number of octaves in the filter bank is such that
        all wavelet scales do not exceed the signal length, nor the
        user-defined maximum scale `max_scale` if present."""
        nOctaves = default_nOctaves(nOctaves, class, log2_size,
            max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave)
        """The number of orientations of a one-dimensional filter bank is
        equal to 1 for real input (`pointgroup::TrivialGroup`), equal to 2 for
        complex input (`pointgroup::ReflectionGroup`)."""
        nΘs = get_nOrientations(pointgroup)
        ħ = uncertainty(class)
        """The meta-information of the wavelets is computed as a 3-d array,
        whose dimensions respectively correspond to
        * the orientation variable `θ`,
        * the chroma variable `χ`, and
        * the octave variable `j`."""
        ψmetas = Array{ΨMeta1D}(nΘs, nFilters_per_octave, nOctaves)
        for j in 0:(nOctaves-1), χ in 0:(nFilters_per_octave-1)
            γ = j * nFilters_per_octave + χ
            resolution = exp2(-γ / nFilters_per_octave)
            centerfrequency = motherfrequency * resolution
            unbounded_scale = ħ * max_qualityfactor / centerfrequency
            scale = min(unbounded_scale, max_scale)
            unbounded_q = scale * centerfrequency / ħ
            qualityfactor = clamp(unbounded_q, 1.0, max_qualityfactor)
            bandwidth = centerfrequency / qualityfactor
            for θ in 0:(nΘs-1)
                ψmetas[1+θ, 1+χ, 1+j] = ΨMeta1D(γ, θ, χ, bandwidth,
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
            max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave,
            nOctaves, pointgroup, signaltype)
        """Before returning the `spec`, we call the function `checkspec`
        which enforces properties of the wavelet filter bank to satisfy null
        mean, limited spatial support, and Littlewood-Paley inequality."""
        checkspec(spec) && return spec
    end
end

function checkspec(spec::Spec1D)
    checkspec_super(spec) && return true
end

"""Given a dimensionless tuning frequency, returns the maximal admissible
mother frequency such that the subsequent wavelets will be in tune with the
tuning frequency.

For example, to tune a 12-chroma Morlet filter bank to a concert pitch of
440 Hz at a sample rate of 44,1 kHz:
    Morlet1DSpec(nFilters_per_octave = 12, tuning_frequency = 440.0/44100.0)"""
function tune_motherfrequency(tuningfrequency, spectype, nFilters_per_octave)
    max_centerfrequency =
        default_motherfrequency(spectype, nFilters_per_octave)
    tuning_ratio = max_centerfrequency / tuningfrequency
    tuning_γ = floor(nFilters_per_octave*log2(tuning_ratio))
    return tuningfrequency * exp2(tuning_γ/nFilters_per_octave)
end
tune_motherfrequency(tuningfrequency::Void, spectype, nFilters_per_octave) =
    default_motherfrequency(spectype, nFilters_per_octave)
