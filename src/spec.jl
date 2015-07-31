# A Spec object contains the immutable specifications of a filter bank
abstract AbstractSpec
abstract Abstract1DSpec{T<:Number} <: AbstractSpec
abstract Abstract2DSpec{T<:Number} <: AbstractSpec

function checkspec(opts::Options)
    if opts[:max_qualityfactor]<1.0
        error("max_qualityfactor cannot be lower than 1.0.")
    end
    if opts[:max_qualityfactor]>opts[:nFilters_per_octave]
        error("max_qualityfactor cannot be greater than nFilters_per_octave.")
    end
end

realtype{T<:Real}(::Type{T}) = T
realtype{T<:Real}(::Type{Complex{T}}) = T

function specgammas(spec::AbstractSpec)
    nΓs = spec.nFilters_per_octave * spec.nOctaves
    γs = collect(0:(nΓs-1))
    chromas = collect(0:(spec.nFilters_per_octave-1))
    χs = repmat(chromas, spec.nOctaves)
    octaves = collect(0:(spec.nOctaves-1))
    js = vec(repmat(transpose(octaves), spec.nFilters_per_octave))
    return (γs, χs, js)
end
