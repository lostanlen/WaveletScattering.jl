abstract AbstractSpec
abstract AbstractSpec1D{T<:Number} <: AbstractSpec
abstract AbstractSpec2D{T<:Number} <: AbstractSpec

function specgammas(spec::AbstractSpec)
    nΓs = spec.nFilters_per_octave * spec.nOctaves
    γs = collect(0:(nΓs-1))
    chromas = collect(0:(spec.nFilters_per_octave-1))
    χs = repmat(chromas, spec.nOctaves)
    octaves = collect(0:(spec.nOctaves-1))
    js = vec(repmat(transpose(octaves), spec.nFilters_per_octave))
    return (γs, χs, js)
end
