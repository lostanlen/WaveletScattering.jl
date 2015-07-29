immutable MorletSpec1D{T<:Number} <: AbstractSpect1D{T}
    ɛ
    log2_length::Int
    max_qualityfactor
    max_scale
    nFilters_per_octave::Int
    nOctaves::Int
    function Morlet1DSpec{T<:Number}(ɛ, log2_length, max_qualityfactor,
                                     max_scale, nFilters_per_octave, nOctaves)
        @assert isa(ɛ, realtype(T))
        @assert isa(max_qualityfactor, realtype(T))
        @assert isa(max_scale, realtype(T))
        new(ɛ, log2_length, max_qualityfactor, max_scale,
            nFilters_per_octave, nOctaves)
    end
end
