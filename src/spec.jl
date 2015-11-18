"""An `AbstractSpec` object contains all the immutable specifications in a
wavelet filter bank"""
abstract AbstractSpec{
    T<:Number,
    D<:AbstractDomain,
    G<:AbstractPointGroup,
    W<:RedundantWaveletClass}

"""Enforces properties of the wavelets to satisfy null mean, limited spatial
support, and Littlewood-Paley inequality.

* The truncation threshold `ɛ` must be `in [0.0, 1.0[`.

* The signal length must be a power of two above 4, i.e. `log2_size > 1`.

* The maximum required quality factor `max_qualityfactor` must be between `1.0`
and the number of filters per octaves.

* The maximum scale must `max_scale` must be at least about 5 times above the
maximum quality factor.

* The mother frequency `motherfrequency` must be in `]0.0, 0.5]`.

* The lowest center frequency must be greater or equal than the number of
per octaves, i.e. `(log2_size-nOctaves) >= 1 + log2(nFilters_per_octave)`."""
function checkspec(spec::AbstractSpec)
    if (spec.ɛ >= 1.0) || (spec.ɛ < 0.0) || (spec.ɛ === -0.0)
        error("`ɛ` must be in `[0.0, 1.0[`. A typical value is `1e-4`.")
    end
    if any(collect(spec.log2_size) .< 2)
        error("Too short signal length.\n",
        "`log2_size =` ", spec.log2_size, "\n",
        "All elements in `log2_size` must be `≧2`")
    end
    if spec.max_qualityfactor < 1.0
        error("Too small maximum quality factor.\n",
        "`max_qualityfactor =` ", spec.max_qualityfactor, "\n",
        "`max_qualityfactor` must be `≧1.0`.")
    end
    if spec.motherfrequency <= 0.0 || spec.motherfrequency > 0.5
        error("`motherfrequency` must be in `]0.0, 0.5]`.")
    end
    if spec.nFilters_per_octave < 1
        error("Too few filters per octave.\n",
        "`nFilters_per_octave = `", spec.nFilters_per_octave, "\n",
        "`nFilters_per_octave` must be `≧1`.")
    end
    if spec.nOctaves < 1
        error("Too few octaves.\n",
        "`nOctaves = `", spec.nOctaves, "\n",
        "`nOctaves` must be `≧1`")
    end
    if spec.max_qualityfactor > spec.nFilters_per_octave
        error("Too few filters per octave for the given quality factor.\n",
        "`max_qualityfactor = `", spec.max_qualityfactor, "\n",
        "`nFilters_per_octave = `", spec.nFilters_per_octave, "\n",
        """The inequality `nFilters_per_octave ≧ max_qualityfactor` must be
        satisfied.""")
    end
    log2_sizes = collect(spec.log2_size)
    if any(spec.nOctaves .>= log2_sizes)
        error("Too many octaves.\n",
        "`log2_size = `", spec.log2_size, "\n",
        "`nOctaves = `", spec.nOctaves, "\n",
        """The inequality `minimum(log2_size) > nOctaves` must be satisfied.
        Either increase `log2_size` or decrease `nOctaves`.""")
    end
    if any(log2_sizes-spec.nOctaves .< 1+log2(spec.nFilters_per_octave))
        error("Too many filters per octave for the given length.\n",
        "`log2_size = `", spec.log2_size, "\n",
        "`log2(nFilters_per_octave) = `", log2(spec.nFilters_per_octave), "\n",
        "`nOctaves = `", spec.nOctaves, "\n",
        """The inequality
        `minimum(log2_size)-nOctaves ≧ 1 + log2(nFilters_per_octave)`
        must be satisfied. Either increase `log2_size`, decrease `nOctaves`,
        or decrease `nFilters_per_octave`.""")
    end
    maximumscale = mapreduce(m -> m.scale, max, spec.metas)
    if maximumscale > (spec.max_scale + 1e-3)
        error("Required time-frequency localization is too tight.\n",
        "`max_qualityfactor = `", spec.max_qualityfactor, "\n",
        "`max_scale = `", spec.max_scale, "\n",
        "`motherfrequency = `", spec.motherfrequency, "\n",
        "The wavelet ", typeof(spec), "cannot have both a bandwidth `< ``",
        spec.motherfrequency / spec.max_qualityfactor,
        "and a scale `<` ", spec.max_scale, ".\n",
        "Either decrease `max_qualityfactor` or decrease `max_scale`.")
    end
    if maximumscale > (2.0^minimum(spec.log2_size) + 1e-3)
        min_resolution = 2.0^(-spec.nOctaves/spec.nFilters_per_octave)
        min_centerfrequency = spec.motherfrequency * min_resolution
        max_bandwidth = min_centerfrequency / spec.max_qualityfactor
        size = tuple(1 .<< collect(spec.log2_size)...)
        error("Spatial localization is coarser than signal length.\n",
        "`log2_size = ``", spec.log2_size,
        "`max_qualityfactor = ``", spec.max_qualityfactor,
        "`max_scale = ``", spec.max_scale,
        "`motherfrequency = ``", spec.motherfrequency, "\n",
        "`nOctaves = ``", spec.nOctaves, "\n",
        "The wavelet ", typeof(spec), "cannot have both a bandwidth ``< ",
        "motherfrequency*2^(-nOctaves)/qualityfactor = ``", max_bandwidth,
        "and a scale `< 2^(log2_size) = `", size, ".\n",
        """Either increase log2_size, decrease max_qualityfactor,
        set `max_scale <= log2_size`, or decrease `nOctaves`.""")
    end
    return true
end

"""Given an signal input type `T`, returns a conservative default value for
the truncation threshold `ɛ of the wavelets.
The chosen values are powers of ten within the same order of magnitude as
floating-point machine precision `eps(T)`."""
default_ɛ(T::Type{Float16}) = 1e-3
default_ɛ(T::Type{Float32}) = 1e-7
default_ɛ(T::Type{Float64}) = 1e-15
default_ɛ{RealT}(T::Type{Complex{RealT}}) = default_ɛ(RealT)

"""Given a maximum quality factor and a number of filter per octaves (both of
which may be `Void`), returns the maximum quality factor in a wavelet filter
bank."""
default_max_qualityfactor(max_q::Real, nfo::Any) = Float64(max_q)
default_max_qualityfactor(max_q::Void, nfo::Integer) = Float64(nfo)
default_max_qualityfactor(max_q::Void, nfo::Void) = 1.0

"""The dimensionless mother center frequency `ξ` (corresponding to a log-period
`γ=0`) is computed as the midpoint between the center frequency of the second
center frequency `ξ*2^(-1/nFilters_per_octave)` (corresponding to `γ=1`) and the
negative mother center frequency `(1-ξ)`. Hence the equation
`2ξ = ξ*2^(-1/nFilters_per_octave) + (1-ξ)`, of which we
derive `ξ = 1 / (3 - 2^(1/nFilters_per_octave))`. This formula is valid
only when the wavelet is a symmetric bump in the Fourier domain."""
default_motherfrequency(class::RedundantWaveletClass, nFilters_per_octave) =
    inv(3.0 - exp2(-inv(nFilters_per_octave)))

"""Given a maximum quality factor and a number of filter per octaves (both of
which may be `Void`), returns the default number of filters per octave in a
wavelet filter bank."""
default_nFilters_per_octave(nfo::Integer, max_q::Any) = Int(nfo)
default_nFilters_per_octave(nfo::Void, max_q::Real) = ceil(Int, max_q)
default_nFilters_per_octave(nfo::Void, max_q::Void) = 1

"""Returns the maximal number octaves in a filter bank such that all scales are
below `2^(log2_size)`."""
default_nOctaves(nOctaves::Int, args...) = nOctaves
function default_nOctaves(nOctaves::Void, class::RedundantWaveletClass,
        log2_size::Tuple, max_qualityfactor::Float64, max_scale::Real,
        motherfrequency::Float64, nFilters_per_octave::Int, args...)
    siglength = 1 << minimum(log2_size)
    if max_scale > siglength
        min_centerfrequency = uncertainty(class) / siglength * max_qualityfactor
    else
        min_centerfrequency = uncertainty(class) / max_scale * 1.0
    end
    nOctaves_a = floor(Int, log2(motherfrequency / min_centerfrequency))
    nOctaves_b = minimum(log2_size) - 1 - ceil(Int, log2(nFilters_per_octave))
    return min(nOctaves_a, nOctaves_b)
end

"""Given a dimensionless tuning frequency, returns the maximal admissible
mother frequency such that the subsequent wavelets will be in tune with the
tuning frequency.

For example, to tune a 12-chroma Morlet filter bank to a concert pitch of
440 Hz at a sample rate of 44,1 kHz:
    Morlet1DSpec(nFilters_per_octave = 12, tuning_frequency = 440.0/44100.0)"""
function tune_motherfrequency(tuningfrequency, spectype, nFilters_per_octave)
    max_centerfrequency = default_motherfrequency(spectype, nFilters_per_octave)
    tuning_ratio = max_centerfrequency / tuningfrequency
    tuning_γ = floor(nFilters_per_octave*log2(tuning_ratio))
    return tuningfrequency * exp2(tuning_γ/nFilters_per_octave)
end
tune_motherfrequency(tuningfrequency::Void, spectype, nFilters_per_octave) =
    default_motherfrequency(spectype, nFilters_per_octave)

"""Fallback of the uncertainty constant from the spec to its class. The RHS
method must be specifically implemented by AbstractSpec concrete subtypes."""
uncertainty(spec::AbstractSpec) = uncertainty(spec.class)
