immutable Spec1D{T<:FFTW.fftwReal,D<:LineDomains,
        G<:LineGroups,W<:RedundantWaveletClass} <: AbstractSpec{T,D,G}
    ɛ::Float64
    ϕmeta::ΦMeta
    ψmetas::Array{ΨMeta,3}
    class::W
    domain::D
    log2_size::Tuple{Int}
    max_qualityfactor::Float64
    max_scale::Float64
    motherfrequency::Float64
    nFilters_per_octave::Int
    nOctaves::Int
    pointgroup::G
    signaltype::Type{T}
    function call{T,D,G,W}(
            ::Type{Spec1D} ;
            class::W = Morlet(),
            pointgroup::G = TrivialGroup(),
            signaltype::Type{T} = Float32,
            domain::D = FourierDomain(1),
            ɛ = default_ɛ(signaltype),
            log2_size = 15,
            max_qualityfactor = nothing,
            max_scale = Inf,
            nFilters_per_octave = nothing,
            nOctaves = nothing,
            tuningfrequency = nothing)
        "Integer log2_size is automatically converted to one-element tuple"
        isa(log2_size, Int) && (log2_size = tuple(log2_size))
        max_qualityfactor, nFilters_per_octave =
             default_max_qualityfactor(max_qualityfactor, nFilters_per_octave),
             default_nFilters_per_octave(nFilters_per_octave, max_qualityfactor)
        motherfrequency =
            tune_motherfrequency(tuningfrequency, class, nFilters_per_octave)
        nOctaves = default_nOctaves(nOctaves, class, log2_size,
            max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave)
        ψmetas = waveletmetas(spec)
        ϕmeta = lowpassmeta(spec)
        spec = new{T,D,G,W}(ɛ, class, domain, log2_size, max_qualityfactor,
            max_scale, motherfrequency, nFilters_per_octave, nOctaves,
            pointgroup, signaltype)
        checkspec(spec) && return spec
    end
end
