# A Spec object contains the immutable specifications of a filter bank
abstract AbstractSpec
abstract Abstract1DSpec{T<:Number} <: AbstractSpec
abstract Abstract2DSpec{T<:Number} <: AbstractSpec

# checkspec enforces properties of the wavelets to satisfy null mean and
# Littlewood-Paley inequality. See documentation for details.
function checkspec(ɛ, log2_length, max_qualityfactor, max_scale,
  nFilters_per_octave, nOctaves)
    for (k,v) in Dict("ɛ"=>ɛ, "log2_length"=>log2_length,
      "max_qualityfactor"=>max_qualityfactor,
      "nFilters_per_octave"=>nFilters_per_octave, "nOctaves"=>nOctaves)
        isfinite(v) || errror(k, "cannot be infinite")
    end
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
    if (log2_length-nOctaves) <= log2(nFilters_per_octave)
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
    if max_qualityfactor > oftype(max_qualityfactor, nFilters_per_octave)
        error("Too few filters per octave for the given quality factor.\n",
        "max_qualityfactor = ", max_qualityfactor, "\n",
        "nFilters_per_octave = ", nFilters_per_octave, "\n",
        """The inequality nFilters_per_octave≧max_qualityfactor must be
        satisfied.""")
    end
end

# realtype provides the type parameter of a complex type e.g. Complex{Float32}.
# For numeric real types, e.g. Float32, it is a no-op.
realtype{T<:Real}(::Type{T}) = T
realtype{T<:Real}(::Type{Complex{T}}) = T

# specgammas yields log-scales γs, chromas χs, and octaves js of a given spec.
function specgammas(spec::AbstractSpec)
  nΓs = spec.nFilters_per_octave * spec.nOctaves
  γs = collect(0:(nΓs-1))
  chromas = collect(0:(spec.nFilters_per_octave-1))
  χs = repmat(chromas, spec.nOctaves)
  octaves = collect(0:(spec.nOctaves-1))
  js = vec(repmat(transpose(octaves), spec.nFilters_per_octave))
  return (γs, χs, js)
end
