immutable Spec1D{T<:Real,D<:LineDomains,
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
        nΘs = get_nOrientations(pointgroup)
        ħ = uncertainty(class)
        ψmetas = Array{ΨMeta}(nΘs, nFilters_per_octave, nOctaves)
        for j in 0:(nOctaves-1)
            for χ in 0:(nFilters_per_octave-1)
                γ = j * nFilters_per_octave + χ
                resolution = exp2(-γ / nFilters_per_octave)
                centerfrequency = motherfrequency * resolution
                unbounded_scale = ħ * max_qualityfactor / centerfrequency
                scale = min(unbounded_scale, max_scale)
                unbounded_q = scale * centerfrequency / ħ
                qualityfactor = clamp(unbounded_q, 1.0, max_qualityfactor)
                bandwidth = centerfrequency / qualityfactor
                for θ in 0:(nΘs-1)
                    ψmetas[1+θ, 1+χ, 1+j] = ΨMeta(γ, θ, χ, bandwidth,
                        centerfrequency, j, qualityfactor, scale)
                end
            end
        end
        ϕbandwidth = motherfrequency * exp2(-nOctaves)
        ϕscale = ħ / ϕbandwidth
        ϕmeta = ΦMeta(ϕbandwidth, ϕscale)
        spec = new{T,D,G,W}(ɛ, ϕmeta, ψmetas, class, domain, log2_size,
            max_qualityfactor, max_scale, motherfrequency, nFilters_per_octave,
            nOctaves, pointgroup, signaltype)
        checkspec(spec) && return spec
    end
end
