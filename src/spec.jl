# A Spec object contains the immutable specifications of a filter bank
abstract AbstractSpec
abstract Abstract1DSpec{T<:Number} <: AbstractSpec
abstract Abstract2DSpec{T<:Number} <: AbstractSpec

function specgammas(spec::AbstractSpec)
    nΓs = spec.nFilters_per_octave * spec.nOctaves
    γs = collect(0:(nΓs-1))
    chromas = collect(0:(spec.nFilters_per_octave-1))
    χs = repmat(chromas, spec.nOctaves)
    octaves = collect(0:(spec.nOctaves-1))
    js = vec(repmat(transpose(octaves), spec.nFilters_per_octave))
    return (γs, χs, js)
end

realtype{T<:Real}(::Type{T}) = T
realtype{T<:Real}(::Type{Complex{T}}) = T
