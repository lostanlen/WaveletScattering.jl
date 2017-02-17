"""An `AbstractSpec` object contains all the immutable specifications in a
wavelet filter bank"""
abstract AbstractSpec{
    T<:Number,
    D<:AbstractDomain,
    G<:AbstractPointGroup,
    W<:RedundantWaveletClass}

"""Enforces properties of the wavelets to satisfy null mean, limited spatial
support, and Littlewood-Paley inequality.

* If the wavelet is dyadic, the maximum required quality factor
`max_qualityfactor` must be equal to `1.0`.

* The truncation threshold `ɛ` must be `in [0.0, 1.0[`.

* The signal length must be a power of two above 4, i.e. `log2_size > 1`.

* The maximum required quality factor `max_qualityfactor` must be between `1.0`
and the number of filters per octaves.

* The maximum scale must `max_scale` must be at least about 5 times above the
maximum quality factor.

* The mother frequency `motherfrequency` must be in `]0.0, 0.5]`.

* The lowest center frequency must be greater or equal than the number of
per octaves, i.e. `(log2_size-n_octaves) >= 1 + log2(n_filters_per_octave)`."""
function checkspec_super(spec::AbstractSpec)
    if isdyadic(spec.class) && (spec.max_qualityfactor != 1.0)
        error("Maximum quality factor must be equal to `1.0` ",
        "for dyadic wavelets.\n",
        "`class` = ", spec.class, "\n",
        "`max_qualityfactor` = ", spec.max_qualityfactor, ".")
    end
    if (spec.ɛ >= 1.0) || (spec.ɛ < 0.0) || (spec.ɛ === -0.0)
        error("`ɛ` must be in `[0.0, 1.0[`. A typical value is `1e-4`.")
    end
    if spec.log2_size .< 2
        error("Too short signal length.\n",
        "`log2_size =` ", spec.log2_size, "\n",
        "`log2_size` must be `≧2`")
    end
    if spec.max_qualityfactor < 1.0
        error("Too small maximum quality factor.\n",
        "`max_qualityfactor =` ", spec.max_qualityfactor, "\n",
        "`max_qualityfactor` must be `≧1.0`.")
    end
    if spec.motherfrequency <= 0.0 || spec.motherfrequency > 0.5
        error("`motherfrequency` must be in `]0.0, 0.5]`.")
    end
    if spec.n_filters_per_octave < 1
        error("Too few filters per octave.\n",
        "`n_filters_per_octave = `", spec.n_filters_per_octave, "\n",
        "`n_filters_per_octave` must be `≧1`.")
    end
    if spec.n_octaves < 1
        error("Too few octaves.\n",
        "`n_octaves = `", spec.n_octaves, "\n",
        "`n_octaves` must be `≧1`")
    end
    if spec.max_qualityfactor > spec.n_filters_per_octave
        error("Too few filters per octave for the given quality factor.\n",
        "`max_qualityfactor = `", spec.max_qualityfactor, "\n",
        "`n_filters_per_octave = `", spec.n_filters_per_octave, "\n",
        """The inequality `n_filters_per_octave ≧ max_qualityfactor` must be
        satisfied.""")
    end
    if spec.n_octaves >= spec.log2_size
        error("Too many octaves.\n",
        "`log2_size = `", spec.log2_size, "\n",
        "`n_octaves = `", spec.n_octaves, "\n",
        """The inequality `minimum(log2_size) > n_octaves` must be satisfied.
        Either increase `log2_size` or decrease `n_octaves`.""")
    end
    if spec.log2_size-spec.n_octaves < 1+log2(spec.n_filters_per_octave)
        error("Too many filters per octave for the given length.\n",
        "`log2_size = `", spec.log2_size, "\n",
        "`log2(n_filters_per_octave) = `", log2(spec.n_filters_per_octave), "\n",
        "`n_octaves = `", spec.n_octaves, "\n",
        """The inequality
        `minimum(log2_size)-n_octaves ≧ 1 + log2(n_filters_per_octave)`
        must be satisfied. Either increase `log2_size`, decrease `n_octaves`,
        or decrease `n_filters_per_octave`.""")
    end
    empirical_max_ψscale = mapreduce(get_scale, max, spec.ψmetas)
    empirical_max_scale = max(empirical_max_ψscale, spec.ϕmeta.scale)
    if empirical_max_scale > (spec.max_scale + 1e-3)
        error("Required time-frequency localization is too tight.\n",
        "`max_qualityfactor = `", spec.max_qualityfactor, "\n",
        "`max_scale = `", spec.max_scale, "\n",
        "`motherfrequency = `", spec.motherfrequency, "\n",
        "The wavelet ", typeof(spec), "cannot have both a bandwidth `< ``",
        spec.motherfrequency / spec.max_qualityfactor,
        "and a scale `<` ", spec.max_scale, ".\n",
        "Either decrease `max_qualityfactor` or decrease `max_scale`.")
    end
    if empirical_max_scale > (2.0^spec.log2_size + 1e-3)
        min_resolution = 2.0^(-spec.n_octaves/spec.n_filters_per_octave)
        min_centerfrequency = spec.motherfrequency * min_resolution
        max_bandwidth = min_centerfrequency / spec.max_qualityfactor
        size = 1 .<< spec.log2_size
        error("Spatial localization is coarser than signal length.\n",
        "`log2_size = ``", spec.log2_size,
        "`max_qualityfactor = ``", spec.max_qualityfactor,
        "`max_scale = ``", spec.max_scale,
        "`motherfrequency = ``", spec.motherfrequency, "\n",
        "`n_octaves = ``", spec.n_octaves, "\n",
        "The wavelet ", typeof(spec), "cannot have both a bandwidth ``< ",
        "motherfrequency*2^(-n_octaves)/qualityfactor = ``", max_bandwidth,
        "and a scale `< 2^(log2_size) = `", size, ".\n",
        """Either increase log2_size, decrease max_qualityfactor,
        set `max_scale <= log2_size`, or decrease `n_octaves`.""")
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
wavelet `ξ*2^(-1/n_filters_per_octave)` (corresponding to `γ=1`) and the
negative mother center frequency `(1-ξ)`. Hence the equation
`2ξ = ξ*2^(-1/n_filters_per_octave) + (1-ξ)`, of which we
derive `ξ = 1 / (3 - 2^(-1/n_filters_per_octave))`. This formula is valid
only when the wavelet is a symmetric bump in the Fourier domain."""
default_motherfrequency(class::RedundantWaveletClass, n_filters_per_octave) =
    inv(3.0 - exp2(-inv(n_filters_per_octave)))

"""Given a maximum quality factor and a number of filter per octaves (both of
which may be `Void`), returns the default number of filters per octave in a
wavelet filter bank."""
default_n_filters_per_octave(nfo::Integer, max_q::Any) = Int(nfo)
default_n_filters_per_octave(nfo::Void, max_q::Real) = ceil(Int, max_q)
default_n_filters_per_octave(nfo::Void, max_q::Void) = 1

"""Returns the maximal number octaves in a filter bank such that all scales are
below `2^(log2_size)`."""
default_n_octaves(n_octaves::Int, args...) = n_octaves
function default_n_octaves(n_octaves::Void, class::RedundantWaveletClass,
        log2_size::Int, max_qualityfactor::Float64, max_scale::Real,
        motherfrequency::Float64, n_filters_per_octave::Int, args...)
    siglength = 1 << minimum(log2_size)
    if max_scale > siglength
        min_centerfrequency = uncertainty(class) / siglength * max_qualityfactor
    else
        min_centerfrequency = uncertainty(class) / max_scale * 1.0
    end
    n_octaves_a = floor(Int, log2(motherfrequency / min_centerfrequency))
    n_octaves_b = minimum(log2_size) - 1 - ceil(Int, log2(n_filters_per_octave))
    return min(n_octaves_a, n_octaves_b)
end

"""Fallback of the uncertainty constant from the spec to its class. The RHS
method must be specifically implemented by `AbstractSpec` concrete subtypes."""
uncertainty(spec::AbstractSpec) = uncertainty(spec.class)
