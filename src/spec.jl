# A Spec object contains the immutable specifications of a filter bank
abstract AbstractSpec
abstract Abstract1DSpec{T<:Number} <: AbstractSpec
abstract Abstract2DSpec{T<:Number} <: AbstractSpec

# checkspec enforces properties of the wavelets to satisfy null mean and
# Littlewood-Paley inequality. See documentation for details.
function checkspec(max_qualityfactor, nFilters_per_octave)
  max_qualityfactor<1.0 || error("max_qualityfactor cannot be lower than 1.0.")
  nFilters_per_octave<1 || error("nFilters_per_octave cannot be lower than 1.")
  qualityfactor>nFilters_per_octave ||
      error("max_qualityfactor cannot be lower than nFilters_per_octave")
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
