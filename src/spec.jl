"""An AbstractSpec object contains all the immutable specifications in a
wavelet filter bank"""
abstract AbstractSpec
abstract Abstract1DSpec{T<:Number} <: AbstractSpec
abstract Abstract2DSpec{T<:Number} <: AbstractSpec

"""Enforces properties of the wavelets to satisfy null mean, limited spatial
support, and Littlewood-Paley inequality.

* The truncation threshold ɛ must be in [0.0, 1.0[.

* The signal length must be a finite power of two above 4.

* The maximum required quality factor must be between 1.0 and the number of
filters per octaves.

* The maximum scale must be at least 5 times above the maximum quality factor.

* The lowest center frequency must be higher than 1, i.e. log2_length>nOctaves.

* The lowest center frequency must be higher than half the number of
filters per octaves, i.e. (log2_length-nOctaves) >= log2(nFilters_per_octave).
"""
function checkspec(ɛ, log2_length, max_qualityfactor, max_scale,
                  nFilters_per_octave, nOctaves)
    if ɛ>=one(ɛ) || ɛ<zero(ɛ)
        error("ɛ must be in [0.0, 1.0[. A typical value is 1e-4.")
    end
    if log2_length < 2
        error("Too short signal length.\n",
        "log2_length = ", log2_length, "\n",
        "log2_length must be ≧2")
    end
    if max_qualityfactor < 1.0
        error("Too small maximum quality factor.\n",
        "max_qualityfactor = ", max_qualityfactor, "\n",
        "max_qualityfactor must be ≧1.0.")
    end
    if nFilters_per_octave < 1
        error("Too few filters per octave.\n",
        "nFilters_per_octave = ", nFilters_per_octave, "\n",
        "nFilters_per_octave must be ≧1.")
    end
    if nOctaves < 1
        error("Too few octaves.\n",
        "nOctaves = ", nOctaves, "\n",
        "nOctaves must be ≧1")
    end
    if max_qualityfactor > nFilters_per_octave
        error("Too few filters per octave for the given quality factor.\n",
        "max_qualityfactor = ", max_qualityfactor, "\n",
        "nFilters_per_octave = ", nFilters_per_octave, "\n",
        """The inequality nFilters_per_octave ≧ max_qualityfactor must be
        satisfied.""")
    end
    if nOctaves >= log2_length
        error("To many octaves.\n",
        "log2_length = ", log2_length, "\n",
        "nOctaves = ", nOctaves, "\n",
        """The inequality log2_length > nOctaves must be satisfied. Either
        increase log2_length or decrease nOctaves.""")
    end
    if (log2_length-nOctaves) < log2(nFilters_per_octave)
        error("Too many filters per octave for the given length.\n",
        "log2_length = ", log2_length, "\n",
        "log2(nFilters_per_octave) = ", log2(nFilters_per_octave), "\n",
        "nOctaves = ", nOctaves, "\n",
        """The inequality log2_length-nOctaves ≧ log2(nFilters_per_octave)
        must be satisfied. Either increase log2_length, decrease nOctaves,
        or decrease nFilters_per_octave.""")
    end
    if max_scale < (5.0*max_qualityfactor)
        error("Too small maximum scale for the given quality factor.\n",
        "max_qualityfactor = ", max_qualityfactor, "\n",
        "max_scale = ", max_scale, "\n",
        "The ratio (max_qualityfactor/max_scale) must be ≧5.0.")
    end
    return true
end

"""Given an signal input type T, returns a conservative default value for
the truncation threshold ɛ of the wavelets.
The chosen values are powers of ten within the same order of magnitude as
floating-point machine precision eps(T)."""
default_ɛ(T::Type{Float16}) = 1e-3
default_ɛ(T::Type{Float32}) = 1e-7
default_ɛ(T::Type{Float64}) = 1e-15
default_ɛ{RealT}(T::Type{Complex{RealT}}) = default_ɛ(RealT)

"""Given a maximum quality factor and a number of filter per octaves (both of
which may be Void), returns the maximum quality factor in a wavelet filter
bank."""
default_max_qualityfactor(max_q::Float64, nfo::Any) = max_q
default_max_qualityfactor(max_q::Void, nfo::Integer) = Float64(nfo)
default_max_qualityfactor(max_q::Void, nfo::Void) = 1.0

"""Given a maximum quality factor and a number of filter per octaves (both of
which may be Void), returns the default number of filters per octave in a
wavelet filter bank."""
default_nFilters_per_octave(max_q::Any, nfo::Int) = nfo
default_nFilters_per_octave(max_q::Float64, nfo::Void) = ceil(Int, max_q)
default_nFilters_per_octave(max_q::Void, nfo::Void) = 1

"Generic fallback for default_nOctaves when nOctaves is not a Void input."
default_nOctaves{T<:AbstractSpec}(::Type{T}, nOctaves::Int, args...) = nOctaves

"""Returns the type parameter of a complex type.
For example, realtype(Complex{Float32}) returns Float32.
For numeric real types, e.g. Float32, it is a no-op."""
realtype{T<:Real}(::Type{T}) = T
realtype{T<:Real}(::Type{Complex{T}}) = T

"""Returns a tuple of log-scales `γs`, chromas `χs`, and octaves `js` in a
given spec. The input spec must have fields nFilters_per_octave and nOctaves."""
function specgammas(spec::AbstractSpec)
    nΓs = spec.nFilters_per_octave * spec.nOctaves
    γs = collect(0:(nΓs-1))
    chromas = collect(0:(spec.nFilters_per_octave-1))
    χs = repmat(chromas, spec.nOctaves)
    octaves = collect(0:(spec.nOctaves-1))
    js = vec(repmat(transpose(octaves), spec.nFilters_per_octave))
    return (γs, χs, js)
end

"""Given a dimensionless tuning frequency, returns the maximal admissible
mother frequency such that the subsequent wavelets will be in tune with the
tuning frequency.

```ξ = tune(Morlet1DSpec, 12, 440.0/44100.0)
Morlet1DSpec(nFilters_per_octave=12, mother_frequency=ξ)``
tunes a 12-chroma filter bank to a concert pitch of 440Hz."""
function tune{T<:AbstractSpec}(::Type{T}, nFilters_per_octave, tuning_frequency)
    max_centerfrequency = default_motherfrequency(T, nFilters_per_octave)
    tuning_ratio = max_centerfrequency / tuning_frequency
    tuning_γ = floor(log2(spec.nFilters_per_octave*tuning_ratio))
    return tuning_frequency * exp2(tuning_γ/spec.nFilters_per_octave)
end
