immutable Spec1D{T<:FFTW.fftwReal,D<:LineDomains,
        G<:LineGroups,W<:RedundantWaveletClass} <: AbstractSpec{T,D,G}
    ɛ::Float64
    class::Type{W}
    domain::Type{D}
    log2_size::Tuple{Int}
    max_qualityfactor::Float64
    max_scale::Float64
    motherfrequency::Float64
    nFilters_per_octave::Int
    nOctaves::Int
    pointgroup::Type{G}
    signaltype::Type{T}
    function call{T<:FFTW.fftwReal,D<:LineDomains,G<:LineGroups}(
            ::Type{Spec1D}, ::Type{W} = Morlet, ::Type{T} = Float32,
            ::Type{D} = FourierDomain{1}, ::Type{G} = TrivialGroup;
            ɛ = default_ɛ(T), log2_size = 15, max_qualityfactor = nothing,
            max_scale = Inf, nFilters_per_octave = nothing, nOctaves = nothing,
            tuningfrequency = nothing)
        "Integer log2_size is automatically converted to one-element tuple"
        isa(log2_size, Int) && (log2_size = tuple(log2_size))
        max_qualityfactor, nFilters_per_octave =
             default_max_qualityfactor(max_qualityfactor, nFilters_per_octave),
             default_nFilters_per_octave(nFilters_per_octave, max_qualityfactor)
        motherfrequency =
            tune_motherfrequency(tuningfrequency, W, nFilters_per_octave)
        nOctaves = default_nOctaves(nOctaves, W, log2_size,
            max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave)
        spec = new{T,D,G,W}(ɛ, W, D, log2_size, max_qualityfactor, max_scale,
            motherfrequency, nFilters_per_octave, nOctaves, G, T)
        checkspec(spec) && return spec
    end
end
