immutable Morlet1DSpec{T<:Number} <: AbstractSpec1D{T}
    ɛ # ::realtype(T) (imposed by inner constructor)
    signaltype::Type{T}
    log2_length::Int
    max_qualityfactor # ::realtype(T) (imposed by inner constructor)
    max_scale # ::realtype(T) (imposed by inner constructor)
    nFilters_per_octave::Int
    nOctaves::Int
    function Morlet1DSpec{T<:Number}(ɛ, log2_length, max_qualityfactor,
                                     max_scale, nFilters_per_octave, nOctaves)
        @assert isa(ɛ, realtype(T))
        @assert isa(max_qualityfactor, realtype(T))
        @assert isa(max_scale, realtype(T))
        new{T}(ɛ, log2_length, max_qualityfactor, max_scale,
            nFilters_per_octave, nOctaves)
    end
end
