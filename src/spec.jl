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
    if nOctaves < 1
        error("Too few octaves.\n",
        "nOctaves = ", nOctaves, "\n",
        "nOctaves must be ≧1")
    end
    if (log2_length-nOctaves) < 2
        error("Wavelet support is too large in low frequencies.\n",
        "Either increase log2_length or decrease nOctaves.\n",
        "log2_length = ", log2_length, "\n",
        "nOctaves = ", nOctaves, "\n",
        "The difference (log2_length-nOctaves) must be ≧2.")
    end
    if (log2_length-nOctaves) < 1+log2(nFilters_per_octave)
        error("Too many filters per octave for the given length.\n",
        "log2_length = ", log2_length, "\n",
        "log2(nFilters_per_octave) = ", log2(nFilters_per_octave), "\n",
        "nOctaves = ", nOctaves, "\n",
        """The inequality log2_length-nOctaves > log2(nFilters_per_octave)
        must be satisfied. Either increase log2_length, decrease nOctaves,
        or decrease nFilters_per_octave.""")
    end
    if max_qualityfactor < 1.0
        error("Too small maximum quality factor.\n",
        "max_qualityfactor = ", max_qualityfactor, "\n",
        "max_qualityfactor must be ≧1.0.")
    end
    if max_scale < (5.0*max_qualityfactor)
        error("Too small maximum scale for the given quality factor.\n",
        "max_qualityfactor = ", max_qualityfactor, "\n",
        "max_scale = ", max_scale, "\n",
        "The ratio (max_qualityfactor/max_scale) must be ≧5.0.")
    end
    if nFilters_per_octave < 1
        error("Too few filters per octave.\n",
        "nFilters_per_octave = ", nFilters_per_octave, "\n",
        "nFilters_per_octave must be ≧1.")
    end
    if max_qualityfactor > nFilters_per_octave
        error("Too few filters per octave for the given quality factor.\n",
        "max_qualityfactor = ", max_qualityfactor, "\n",
        "nFilters_per_octave = ", nFilters_per_octave, "\n",
        """The inequality nFilters_per_octave≧max_qualityfactor must be
        satisfied.""")
    end
    return true
end

default_epsilon(T::Type{Float16}) = 1e-3
default_epsilon(T::Type{Float32}) = 1e-7
default_epsilon(T::Type{Float64}) = 1e-15
default_epsilon{RealT}(T::Type{Complex{RealT}}) = default_epsilon(RealT)


"""Provides the type parameter of a complex type.

For example, realtype(Complex{Float32}) yields Float32.
For numeric real types, e.g. Float32, it is a no-op."""
realtype{T<:Real}(::Type{T}) = T
realtype{T<:Real}(::Type{Complex{T}}) = T


"""Yields log-scales γs, chromas χs, and octaves js of a given spec.

The input spec should have fields nFilters_per_octave and nOctaves."""
function specgammas(spec::AbstractSpec)
    nΓs = spec.nFilters_per_octave * spec.nOctaves
    γs = collect(0:(nΓs-1))
    chromas = collect(0:(spec.nFilters_per_octave-1))
    χs = repmat(chromas, spec.nOctaves)
    octaves = collect(0:(spec.nOctaves-1))
    js = vec(repmat(transpose(octaves), spec.nFilters_per_octave))
    return (γs, χs, js)
end
